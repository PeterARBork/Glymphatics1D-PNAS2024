
abstract type GlymphaticModel end

mutable struct GlymphaticSimulation <: GlymphaticModel
    prior_parameters::NamedTuple
    initial_condition::Vector{Float64}
    encloser::Function
    fitting_parameters::Union{NamedTuple, Nothing}
    turingmodelgenerator::Function
    prior_predictions_df::Union{DataFrame, Nothing}
    posterior_predictions_df::Union{DataFrame, Nothing}
    ode_solution::Union{ODESolution, Nothing}
    data::Union{DataFrame, Nothing}
    turing_model::Union{Turing.Model, Nothing}
    chains::Union{Chains, Nothing}
    prior_chains::Union{Chains, Nothing}
    posterior_parameters::Union{NamedTuple, Nothing}
    WAIC::Union{Tuple{Float64, Float64}, Nothing}

    GlymphaticSimulation(
        prior_params::NamedTuple,
        initial_condition::Vector{Float64},
        fitting_parameters::NamedTuple,
        turingmodelgenerator::Function,
    ) = new(
        prior_params,
        initial_condition,
        enclosethetimedifferential,
        fitting_parameters,
        turingmodelgenerator,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )

    GlymphaticSimulation(
        prior_params::NamedTuple,
        initial_condition::Vector{Float64},
    ) = new(
        prior_params,
        initial_condition,
        enclosethetimedifferential,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )

    #=
    GlymphaticSimulation(
        prior_params::NamedTuple,
        initial_condition::Vector{Float64},
        fitting_parameters::NamedTuple,
        turingmodelgenerator::Function,
        encloser::Function,
    ) = new(
        prior_params,
        initial_condition,
        encloser,
        fitting_parameters,
        turingmodelgenerator,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )
    =#

end

function discretize_line(Δr, R)
    r_space = range(Δr, step=Δr, stop=R)
    J = length(r_space)

    ord_deriv = 2
    ord_approx = 3
    ∇ = CenteredDifference(1, ord_approx, Δr, J)
    Δ = CenteredDifference(ord_deriv, ord_approx, Δr, J)

    return r_space, ∇, Δ
end

function init_slab_with_mass(r_space, S, mass, center, stddev, nc=2)
    J = length(r_space)
    init = zeros(J + nc)
    #init[1:J] = (mass / S / α) .* exp.(-(1/2) * ((r_space .- center) / stddev).^2) / (stddev * sqrt(2 * pi))
    init[1:J] = (mass / S) .* exp.(-(1/2) * ((r_space .- center) / stddev).^2) / (stddev * sqrt(2 * pi))
    return init
end

function initialize_striatal_infusion_normaldist(
    T_max,
    injected_mass;
    paramsfile="data/inputs/ratparametersfromtables.csv",
    infusion_spread_stdev=0.4,
    Turing_space_steps=1,
    fitting_parameters::NamedTuple=NamedTuple(),
    turingmodelgenerator::Function,
    k_p::Real=0.0,
    )

    tmp_ps = loadparams(paramsfile)
    Δr_bayes = tmp_ps.L / Turing_space_steps
    ps_bayes = loadparams(paramsfile; T_max, k_p, injected_mass, Δr=Δr_bayes, fitting_parameters...)

    order_app = 2
    ∇ = CenteredDifference(1, order_app, ps_bayes.Δr, length(ps_bayes.r_space))
    Δ = CenteredDifference(2, order_app, ps_bayes.Δr, length(ps_bayes.r_space))
    bc = NeumannBC((3.14, 1.0), ps_bayes.Δr)
    Δ, _ = Array(Δ * bc)
    ∇, _ = Array(∇ * bc)

    prior_parameters = update_namedtuple(
        ps_bayes;
        ∇,
        Δ,
    )

    initial_condition = init_slab_with_mass(
        prior_parameters.r_space,
        prior_parameters.S,
        injected_mass,
        prior_parameters.L / 2,
        infusion_spread_stdev,
        3)

    if length(fitting_parameters) > 0
        return GlymphaticSimulation(prior_parameters, initial_condition, fitting_parameters, turingmodelgenerator)
    else
        return GlymphaticSimulation(prior_parameters, initial_condition)
    end
end

function initialize_striatal_infusion_with_initial_condition(
    T_max,
    ic_brain_measurements_concentrations,
    ic_brain_measurements_locations,
    c_v_ic,
    c_c_ic,
    c_b_ic;
    paramsfile="data/inputs/ratparametersfromtables.csv",
    Turing_space_steps=1,
    fitting_parameters::NamedTuple=NamedTuple(),
    turingmodelgenerator::Function
)

    tmp_ps = loadparams(paramsfile)
    Δr_bayes = tmp_ps.L / Turing_space_steps
    ps_bayes = loadparams(paramsfile; T_max, Δr=Δr_bayes, fitting_parameters...)

    order_app = 2
    ∇ = CenteredDifference(1, order_app, ps_bayes.Δr, length(ps_bayes.r_space))
    Δ = CenteredDifference(2, order_app, ps_bayes.Δr, length(ps_bayes.r_space))
    bc = NeumannBC((3.14, 1.0), ps_bayes.Δr)
    Δ, _ = Array(Δ * bc)
    ∇, _ = Array(∇ * bc)

    prior_parameters = update_namedtuple(
        ps_bayes;
        ∇,
        Δ,
    )

    icbrain = icfromdf(
        ic_brain_measurements_concentrations,
        ic_brain_measurements_locations,
        prior_parameters.r_space,
    )
    initial_condition = zeros(length(prior_parameters.r_space) + 3)
    initial_condition[1:end-3] = icbrain
    initial_condition[end-2] = c_v_ic
    initial_condition[end-1] = c_c_ic
    initial_condition[end] = c_b_ic

    injected_mass = sum(prior_parameters.Δr * prior_parameters.S * initial_condition[1:end-3])
    injected_mass += c_v_ic * prior_parameters.V_v
    injected_mass += c_c_ic * prior_parameters.V_c
    injected_mass += c_b_ic * prior_parameters.V_b
    prior_parameters = update_namedtuple(
        prior_parameters;
        injected_mass=injected_mass,
    )

    prior_parameters = addbrainconcindices(prior_parameters, ic_brain_measurements_locations)

    if length(fitting_parameters) > 0
        return GlymphaticSimulation(prior_parameters, initial_condition, fitting_parameters, turingmodelgenerator)
    else
        return GlymphaticSimulation(prior_parameters, initial_condition)
    end
end

function addbrainconcindices(params::NamedTuple, mlocations::Vector{Float64})
    findindex(locationmm) = argmin((params.r_space .- locationmm).^2)

    for (i, mloc) in enumerate(mlocations)
        params = update_namedtuple(
            params;
            Dict(Symbol("brain_conc_index" * string(i)) => findindex(mloc))...
        )
    end
    #=
    params = update_namedtuple(
        params;
        brain_conc_index1=findindex(mlocations[1]),
        brain_conc_index2=findindex(mlocations[2]),
        brain_conc_index3=findindex(mlocations[3]),
        brain_conc_index4=findindex(mlocations[4]),
        brain_conc_index5=findindex(mlocations[5]),
        brain_conc_index6=findindex(mlocations[6]),
        brain_conc_index7=findindex(mlocations[7]),
    )
    =#
    return params
end

function summarize_prediction_chains(
    chains::Chains,
    variable_names::Vector{String},
    endtransform::Function=identity
)

    function get_indices_from_chain_index(chain_index::Symbol)
        collect(parse(Int, s) for s in split(split(split(string(chain_index), "[")[2], "]")[1], ","))
    end

    function get_row_from_chain_index(chain_index::Symbol)
        get_indices_from_chain_index(chain_index::Symbol)[1]
    end

    function get_col_from_chain_index(chain_index::Symbol)
        get_indices_from_chain_index(chain_index::Symbol)[2]
    end

    chain_pred_keys = sort([k for k in keys(chains) if contains(string(k), "data")])
    num_rows = maximum([get_row_from_chain_index(k) for k in chain_pred_keys])
    num_columns = maximum([get_col_from_chain_index(k) for k in chain_pred_keys])
    num_chains = size(chains, 3)
    num_samples_per_chain = size(chains, 1)
    total_num_samples = num_chains * num_samples_per_chain
    #solution_matrix = Array{Float32}(undef, num_rows, num_columns, total_num_samples)
    means = Array{Union{Missing, Float64}, 2}(missing, num_rows, num_columns)
    low = similar(means)
    med = similar(means)
    high = similar(means)

    quantiles_chain_df = quantile(chains; q=[0.32, 0.5, 0.68])

    for chain_index in chain_pred_keys
        row = get_row_from_chain_index(chain_index)
        column = get_col_from_chain_index(chain_index)
        means[row, column] = endtransform.(mean(getindex(chains, chain_index)[:]))
        low[row, column] = endtransform.(quantiles_chain_df[chain_index, Symbol("32.0%")])
        med[row, column] = endtransform.(quantiles_chain_df[chain_index, Symbol("50.0%")])
        high[row, column] = endtransform.(quantiles_chain_df[chain_index, Symbol("68.0%")])
    end

    col_name_suffixes = ["mean", "low-cred", "median", "high-cred"]
    var_name_combos = collect(Base.Iterators.product(variable_names, col_name_suffixes))
    colnames = [a * "_" * b for (a, b) in var_name_combos[:]]

    return DataFrame(hcat(means, low, med, high)[:, :, 1], colnames)
end

export GlymphaticSimulation, discretize_line, summarize_prediction_chains
export initialize_striatal_infusion_normaldist, initialize_striatal_infusion_with_initial_condition
export summarize_prediction_chains, addbrainconcindices
