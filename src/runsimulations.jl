
function makeODEProblem(
    simulation::GlymphaticSimulation,
    parameters,
)#::Tuple{ODEProblem, NamedTuple}

    timedifferential! = simulation.encloser(parameters)
    #@unpack D, Dm, v, T_max = parameters
    #p = [D, Dm, v]
    #p = ComponentArray(D=D, Dm=Dm, v=v)
    p = ComponentArray{Union{Bool,Int64,String,Float64,ForwardDiff.Dual{Nothing, Float64, 5}}}(parameters)
    p = collect(p)
    simulation_problem = ODEProblem(
        timedifferential!,
        simulation.initial_condition,
        (0.0, parameters.T_max),
        ComponentArray(parameters),
    )

    return simulation_problem, p
end

function run_simulation!(
    simulation::GlymphaticSimulation,
    prior_or_posterior_parameters::Union{NamedTuple, Vector{Float64}};
    saveat::Union{Nothing, Vector{Float64}}=nothing,
)
    simulation_problem, _ = makeODEProblem(
        simulation,
        prior_or_posterior_parameters,
    )

    if !isnothing(saveat)
        simulation.ode_solution = solve(
            simulation_problem,
            alg_hints=[:stiff],
            saveat=saveat,
        )
    else
        simulation.ode_solution = solve(simulation_problem, alg_hints=[:stiff])
    end

    solution_matrix = transpose(hcat(simulation.ode_solution.u...))

    (bm, m_v, m_csf, m_b) = state_mass_over_time(
        simulation.ode_solution,
        prior_or_posterior_parameters
    )

    depth = prior_or_posterior_parameters.brain_conc_depth
    conc_loc_index = argmin((prior_or_posterior_parameters.r_space .- depth).^2)
    brain_conc = round.(solution_matrix[:, conc_loc_index], digits=5)

    σ_brain_conc = get(prior_or_posterior_parameters, :σ_brain_conc, 1.0)
    brain_conc_obs = rand.(Normal.(brain_conc, σ_brain_conc))

    σ_brain_mass = get(prior_or_posterior_parameters, :σ_brain_mass, 1.0)
    bm_obs = rand.(Normal.(bm, σ_brain_mass))

    σ_csf_mass = get(prior_or_posterior_parameters, :σ_csf_mass, 1.0)
    m_csf_obs = rand.(Normal.(m_csf, σ_csf_mass))

    c_v = solution_matrix[:, end-2]
    σ_vent = get(prior_or_posterior_parameters, :σ_vent, 1.0)
    c_v_obs = rand.(Normal.(c_v, σ_vent))

    c_csf = solution_matrix[:, end-1]
    σ_csf_conc = get(prior_or_posterior_parameters, :σ_csf_conc, 1.0)
    c_csf_obs = rand.(Normal.(c_csf, σ_csf_conc))

    c_blood = solution_matrix[:, end]
    σ_blood = get(prior_or_posterior_parameters, :σ_blood, 1.0)
    c_blood_obs = rand.(Normal.(c_blood, σ_blood))

    solution_df = DataFrame(
        time=simulation.ode_solution.t,
        brain_mass=bm,
        brain_mass_obs=bm_obs,
        ventricle_mass=m_v,
        CM_mass=m_csf,
        CM_mass_obs=m_csf_obs,
        blood_mass=m_b,
        brain_conc=brain_conc,
        brain_conc_obs=brain_conc_obs,
        ventricle_conc=c_v,
        ventricle_conc_obs=c_v_obs,
        CM_conc=c_csf,
        CM_conc_obs=c_csf_obs,
        blood_conc=c_blood,
        blood_conc_obs=c_blood_obs,
    )

    return solution_df, simulation_problem
end

function run_prior_simulation!(
    simulation::GlymphaticSimulation;
    endtransform::Function=identity,
    savedir::String,
    savefile::String,
)

    simulation.prior_predictions_df, simulation_problem = run_simulation!(
        simulation,
        simulation.prior_parameters
    )

    num_obs = size(simulation.data, 1)
    num_cols = size(simulation.data, 2)
    missing_data = fill(missing, (num_obs, num_cols - 1))
    missing_model = simulation.turingmodelgenerator(
        simulation.data[1:end, 1],
        missing_data,
        simulation_problem,
        simulation.prior_parameters
    )
    n_samples = 500
    simulation.prior_chains = sample(missing_model, Prior(), n_samples)
    simulation.prior_predictions_df = summarize_prediction_chains(
        simulation.prior_chains,
        names(simulation.data)[2:end],
        endtransform,
    )
    simulation.prior_predictions_df[!, "time (h)"] = simulation.data[!, "time (h)"]

    CSV.write(savedir * savefile, simulation.prior_predictions_df)

    return simulation.prior_predictions_df
end

function chain_model_params(chains::Chains)
    params = [k for k in keys(get_sections(chains, :parameters))]
              #if !contains(string(k), "σ")]
    return params
end

function sample_posterior!(
    simulation::GlymphaticSimulation;
    savedir::String,
    predictions_file::String,
    chains_file::String = "tmp_chains.jls",
    loadfile = "",
    num_samples = 500,
    num_chains = 3,
    datapretransform::Function=identity,
    predictionsendtransform::Function=identity,
    WAICsfile::String="",
)
    simulation_prob, p = makeODEProblem(
        simulation,
        simulation.prior_parameters,
    )

    #turing_data = Turing.TArray(Array(convert.(Float32, datapretransform.(simulation.data[1:end, 2:end]))))
    turing_data = Turing.TArray(Array(datapretransform.(simulation.data[1:end, 2:end])))
    simulation.turing_model = simulation.turingmodelgenerator(
        simulation.data[1:end, 1],
        turing_data,
        simulation_prob,
        simulation.prior_parameters
    )

    if length(loadfile) > 0
        simulation.chains = deserialize(loadfile)
    elseif isnothing(simulation.chains) && Threads.nthreads() == num_chains + 1
        #sample(model, NUTS(), MCMCThreads(), 1000, 4)
        @info "sampling with " * string(Threads.nthreads() - 1) * " threads."
        simulation.chains = sample(
            simulation.turing_model,
            NUTS(.65),
            MCMCThreads(),
            num_samples,
            Threads.nthreads() - 1,
            θ=p,
        )
        serialize(savedir * chains_file, simulation.chains)
    elseif isnothing(simulation.chains)
        simulation.chains = mapreduce(
            c -> sample(simulation.turing_model,
                NUTS(.65),
                num_samples,
                θ=p,
            ),
            chainscat,
            1:num_chains,
        )
        serialize(savedir * chains_file, simulation.chains)
    else
        println("Keeping chains already sampled.")
    end

    params = chain_model_params(simulation.chains)
    vals = [mean(simulation.chains[p]) for p in params]

    parameter_update = NamedTuple(zip(params, vals))

    simulation.posterior_parameters = update_namedtuple(
        simulation.prior_parameters;
        parameter_update...
    )

    missing_model = simulation.turingmodelgenerator(
        simulation.data[1:end, 1],
        fill(missing, size(simulation.data[1:end, 2:end])),
        simulation_prob,
        simulation.posterior_parameters,
    )
    post_pred = predict(missing_model, simulation.chains)
    simulation.posterior_predictions_df = summarize_prediction_chains(
        post_pred,
        names(simulation.data)[2:end],
        predictionsendtransform,
    )
    simulation.posterior_predictions_df[!, "time (h)"] = simulation.data[1:end, 1]
    CSV.write(savedir * predictions_file, simulation.posterior_predictions_df)

    simulation_prob, _ = makeODEProblem(
        simulation,
        simulation.posterior_parameters,
    )
    simulation.ode_solution = solve(simulation_prob, alg_hints=[:stiff])

    if length(WAICsfile) != 0
        calc_and_save_WAIC!(simulation)
        open(savedir * WAICsfile, "w") do io
            write(io, "WAIC,std\n$(simulation.WAIC[1]),$(simulation.WAIC[2])\n")
        end
    end

    return simulation.chains
end

function makeparamcombinationstring(fittingparameters::NamedTuple)::String
    paramcombinationstring = ""
    if fittingparameters.estimate_D
        paramcombinationstring *= "_D"
    end
    if fittingparameters.estimate_Dm
        paramcombinationstring *= "_Dm"
    end
    if fittingparameters.estimate_v
        paramcombinationstring *= "_v"
    end

    return paramcombinationstring
end

export run_simulation!, run_prior_simulation!, sample_posterior!, makeparamcombinationstring
