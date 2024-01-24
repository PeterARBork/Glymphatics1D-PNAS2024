function parseandanalyzeCserrdataset(
    fitting_parameters::NamedTuple;
    n_samples::Int64=1000,
    Turing_space_steps::Int64=15,
    turingmodelgenerator=posteriorsamplergeneratorCserr,
    raw_datafile="data/inputs/scraped_data.csv",
    parsed_datafile="data/parsed/Cserr_parsed.csv",
    savedir_results="data/results/Cserr/",
    savedir_figures="visualizations/Cserr/",
    data_label="Cserr",
    tracer::String="PEG900",
)
    @assert(:estimate_D in keys(fitting_parameters))
    @assert(:estimate_Dm in keys(fitting_parameters))
    @assert(:estimate_v in keys(fitting_parameters))

    paramcombinationstring = makeparamcombinationstring(fitting_parameters)

    k_p_raw_Cserr = 0.2 * 10^(-5) # per second
    k_p = k_p_raw_Cserr * 3600
    standard_D = 0.46

    if tracer == "PEG900"
        D = standard_D
        μ_D = D
        println("Using Cserrs reference for k_p = " * string(k_p))
    elseif tracer == "PEG4000"
        D = (900/4000)^(1/3) * standard_D
        μ_D = D
        parsed_datafile = replace(parsed_datafile, "_parsed" => "_PEG4000_parsed")
        savedir_results = replace(savedir_results, "Cserr/" => "Cserr/PEG4000/")
        savedir_figures = replace(savedir_figures, "Cserr/" => "Cserr/PEG4000/")
    elseif tracer == "albumin"
        D = (900/69000)^(1/3) * standard_D
        μ_D = D
        parsed_datafile = replace(parsed_datafile, "_parsed" => "_albumin_parsed")
        savedir_results = replace(savedir_results, "Cserr/" => "Cserr/albumin/")
        savedir_figures = replace(savedir_figures, "Cserr/" => "Cserr/albumin/")
        k_p_raw_Cserr = 0.2 * 10^(-5)
        k_p = k_p_raw_Cserr * 3600
        println("Using Cserrs reference for k_p = " * string(k_p))
    end

    # parse the data
    println("Parsing data.")
    Cserr_df = parse_Cserr_data(raw_datafile, parsed_datafile, tracer)

    T_max = 30.0
    injected_mass = Cserr_df[1, 2]

    # setup the parameter-fitting combination for simulations
    println("Setting up simulations to fit each parameter-combination.")
    sim = initialize_striatal_infusion_normaldist(
        T_max,
        injected_mass;
        paramsfile="data/inputs/ratparametersfromtables.csv",
        infusion_spread_stdev=0.4,
        k_p=k_p,
        Turing_space_steps=Turing_space_steps,
        fitting_parameters=fitting_parameters,
        turingmodelgenerator=turingmodelgenerator,
    )
    sim.prior_parameters = update_namedtuple(
        sim.prior_parameters;
        D = D,
        μ_D = μ_D,
    )

    # load data onto each
    load_data!(sim, parsed_datafile)

    # run each prior simulation
    run_prior_simulation!(
        sim;
        savedir=savedir_results,
        savefile="prior_predictions" * paramcombinationstring * ".csv",
    )
    plot_ode_solution_mass_distribution(
        sim,
        "% infused mass";
        savedir=savedir_figures,
        savefile="prior_mass_distribution" * paramcombinationstring * ".svg",
    )

    # plot priors
    p1 = plot(
        sim.prior_parameters.r_space,
        sim.initial_condition[1:end-3],
        xlabel="mm",
        ylabel="μg / μL",
        label="",
        thickness_scaling=2,
    )
    savefig(
        p1,
        savedir_figures * "Cserr_initial_condition" * paramcombinationstring * ".svg",
    )

    plot_simulation_with_data(
        sim,
        "% infused mass";
        posterior_or_prior="prior",
        data_label=data_label,
        simulation_column="m_br_mean",
        data_column="m_br",
        cred_column_base="m_br",
        rescale_hours_to_mins=false,
        savefile="prior_brain_mass" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    plot_simulation_with_data(
        sim,
        "% inj. mass in CSF";
        posterior_or_prior="prior",
        data_label=data_label,
        simulation_column="m_c_mean",
        data_column="m_c",
        cred_column_base="m_c",
        rescale_hours_to_mins=false,
        savefile="prior_CM_conc" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    # Estimate the posterior with NUTS sampling
    println("Sample posterior")
    sample_posterior!(
        sim;
        savedir=savedir_results,
        predictions_file="posterior_predictions" * paramcombinationstring * ".csv",
        chains_file="chains" * paramcombinationstring * ".jls",
        num_samples=n_samples,
        WAICsfile="WAICs" * paramcombinationstring * ".csv",
        #loadfile=savedir_results * "chains" * paramcombinationstring * ".jls",
    )

    # plot the chains
    p = plot(sim.chains)
    savefig(p, savedir_figures * "chains_plot" * paramcombinationstring * ".svg")

    # plot the posteriors of the main models
    plot_simulation_with_data(
        sim,
        "% infused mass";
        posterior_or_prior="posterior",
        data_label=data_label,
        simulation_column="m_br_mean",
        data_column="m_br",
        cred_column_base="m_br",
        rescale_hours_to_mins=false,
        savefile="posterior_brain_mass" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    plot_simulation_with_data(
        sim,
        "% inj. mass in CSF";
        posterior_or_prior="posterior",
        data_label=data_label,
        simulation_column="m_c_mean",
        data_column="m_c",
        cred_column_base="m_c",
        rescale_hours_to_mins=false,
        savefile="posterior_csf_mass" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    # animate the main model
    println("animate solution")
    simani = animate_solution_mass(
        sim.ode_solution,
        sim.posterior_parameters,
        savedir_figures * "animation" * paramcombinationstring * ".gif",
    )

end # parseandanalyzeCserrdataset

function parseandanalyzeCserrdatasetrobustness(
    fitting_parameters::NamedTuple;
    n_samples::Int64=1000,
    Turing_space_steps::Int64=15,
    raw_datafile="data/inputs/scraped_data.csv",
    parsed_datafile="data/parsed/Cserr_parsed.csv",
    savedir_results="data/results/Cserr/",
    savedir_figures="visualizations/Cserr/",
    data_label="Cserr",
    tracer::String="PEG900",
)
    @assert(:estimate_D in keys(fitting_parameters))
    @assert(:estimate_Dm in keys(fitting_parameters))
    @assert(:estimate_v in keys(fitting_parameters))

    paramcombinationstring = makeparamcombinationstring(fitting_parameters)

    k_p_raw_Cserr = 0.2 * 10^(-5) # per second
    k_p = k_p_raw_Cserr * 3600
    standard_D = 0.46

    if tracer == "PEG900"
        D = standard_D
        μ_D = D
        println("Using Cserrs reference for k_p = " * string(k_p))
    elseif tracer == "PEG4000"
        D = (900/4000)^(1/3) * standard_D
        μ_D = D
        parsed_datafile = replace(parsed_datafile, "_parsed" => "_PEG4000_parsed")
        savedir_results = replace(savedir_results, "Cserr/" => "Cserr/PEG4000/")
        savedir_figures = replace(savedir_figures, "Cserr/" => "Cserr/PEG4000/")
    elseif tracer == "albumin"
        D = (900/69000)^(1/3) * standard_D
        μ_D = D
        parsed_datafile = replace(parsed_datafile, "_parsed" => "_albumin_parsed")
        savedir_results = replace(savedir_results, "Cserr/" => "Cserr/albumin/")
        savedir_figures = replace(savedir_figures, "Cserr/" => "Cserr/albumin/")
        k_p_raw_Cserr = 0.2 * 10^(-5)
        k_p = k_p_raw_Cserr * 3600
        println("Using Cserrs reference for k_p = " * string(k_p))
    end

    # parse the data
    println("Parsing data.")
    Cserr_df = parse_Cserr_data(raw_datafile, parsed_datafile, tracer)

    T_max = 30.0
    injected_mass = Cserr_df[1, 2]

    # setup the parameter-fitting combination for simulations
    println("Setting up simulations to fit each parameter-combination.")
    sim = initialize_striatal_infusion_normaldist(
        T_max,
        injected_mass;
        animal="rat",
        infusion_spread_stdev=0.4,
        k_p=k_p,
        Turing_space_steps=Turing_space_steps,
        fitting_parameters=fitting_parameters,
        turingmodelgenerator=posteriorsamplergeneratorCserrallphysprior,
    )
    sim.prior_parameters = update_namedtuple(
        sim.prior_parameters;
        D = D,
        μ_D = μ_D,
    )
    sim.encloser = enclosethetimedifferentialnophyisoparams

    # load data onto each
    load_data!(sim, parsed_datafile)

    # run each prior simulation
    run_prior_simulation!(
        sim;
        savedir=savedir_results,
        savefile="prior_predictions" * paramcombinationstring * ".csv",
    )
    plot_ode_solution_mass_distribution(
        sim,
        "% infused mass";
        savedir=savedir_figures,
        savefile="prior_mass_distribution" * paramcombinationstring * ".svg",
    )

    # plot priors
    p1 = plot(
        sim.prior_parameters.r_space,
        sim.initial_condition[1:end-3],
        xlabel="mm",
        ylabel="μg / μL",
        label="",
        thickness_scaling=2,
    )
    savefig(
        p1,
        savedir_figures * "Cserr_initial_condition" * paramcombinationstring * ".svg",
    )

    plot_simulation_with_data(
        sim,
        "% infused mass";
        posterior_or_prior="prior",
        data_label=data_label,
        simulation_column="m_br_mean",
        data_column="m_br",
        cred_column_base="m_br",
        rescale_hours_to_mins=false,
        savefile="prior_brain_mass" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    plot_simulation_with_data(
        sim,
        "% inj. mass in CSF";
        posterior_or_prior="prior",
        data_label=data_label,
        simulation_column="m_c_mean",
        data_column="m_c",
        cred_column_base="m_c",
        rescale_hours_to_mins=false,
        savefile="prior_CM_conc" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    # Estimate the posterior with NUTS sampling
    println("Sample posterior")
    sample_posterior!(
        sim;
        savedir=savedir_results,
        predictions_file="posterior_predictions" * paramcombinationstring * ".csv",
        chains_file="chains" * paramcombinationstring * ".jls",
        num_samples=n_samples,
        WAICsfile="WAICs" * paramcombinationstring * ".csv",
        #loadfile=savedir_results * "chains" * paramcombinationstring * ".jls",
    )

    # plot the chains
    p = plot(sim.chains)
    savefig(p, savedir_figures * "chains_plot" * paramcombinationstring * ".svg")

    # plot the posteriors of the main models
    plot_simulation_with_data(
        sim,
        "% infused mass";
        posterior_or_prior="posterior",
        data_label=data_label,
        simulation_column="m_br_mean",
        data_column="m_br",
        cred_column_base="m_br",
        rescale_hours_to_mins=false,
        savefile="posterior_brain_mass" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    plot_simulation_with_data(
        sim,
        "% inj. mass in CSF";
        posterior_or_prior="posterior",
        data_label=data_label,
        simulation_column="m_c_mean",
        data_column="m_c",
        cred_column_base="m_c",
        rescale_hours_to_mins=false,
        savefile="posterior_csf_mass" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    # animate the main model
    println("animate solution")
    simani = animate_solution_mass(
        sim.ode_solution,
        sim.posterior_parameters,
        savedir_figures * "animation" * paramcombinationstring * ".gif",
    )

end # parseandanalyzeCserrdataset

export parseandanalyzeCserrdataset, parseandanalyzeCserrdatasetrobustness
