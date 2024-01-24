
function parseandanalyzePladataset(
    fitting_parameters::NamedTuple;
    n_samples::Int64=1000,
    Turing_space_steps::Int64=15,
    raw_datafile="data/inputs/scraped_data.csv",
    parsed_datafile="data/parsed/Pla_parsed.csv",
    savedir_results="data/results/Pla/",
    savedir_figures="visualizations/Pla/",
    data_label="Pla",
    turingmodelgenerator=posteriorsamplergeneratorPla,
)
    @assert(:estimate_D in keys(fitting_parameters))
    @assert(:estimate_Dm in keys(fitting_parameters))
    @assert(:estimate_v in keys(fitting_parameters))

    paramcombinationstring = makeparamcombinationstring(fitting_parameters)

    # parse the data
    @info "Parsing data."
    Pla_df = parse_Pla_data(raw_datafile, parsed_datafile)

    T_max = 2.0
    injected_mass = Pla_df[1, 2]

    # setup the parameter-fitting combination for simulations
    @info "Setting up simulations to fit each parameter-combination."
    sim = initialize_striatal_infusion_normaldist(
        T_max,
        injected_mass;
        paramsfile="data/inputs/mouseparametersfromtables.csv",
        infusion_spread_stdev=0.4,
        Turing_space_steps=Turing_space_steps,
        fitting_parameters=fitting_parameters,
        turingmodelgenerator=turingmodelgenerator,
    )

    # load data onto each
    load_data!(sim, parsed_datafile)

    # run each prior simulation
    @info "Running prior simulation"
    run_prior_simulation!(
        sim;
        savedir=savedir_results,
        savefile="prior_predictions" * paramcombinationstring * ".csv",
    )
    plot_ode_solution_mass_distribution(
        sim,
        "μg";
        savedir=savedir_figures,
        savefile="prior_mass_distribution" * paramcombinationstring * ".svg",
    )

    # plot priors
    plot_simulation_with_data(
        sim,
        "μg";
        posterior_or_prior="prior",
        data_label=data_label,
        simulation_column="m_br_mean",
        data_column="m_br",
        cred_column_base="m_br",
        rescale_hours_to_mins=true,
        savefile="prior_brain_mass" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    plot_simulation_with_data(
        sim,
        "blood [μg/mL]";
        posterior_or_prior="prior",
        data_label=data_label,
        simulation_column="c_b_mean",
        data_column="c_b",
        cred_column_base="c_b",
        rescale_hours_to_mins=true,
        savefile="prior_blood_conc" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topleft,
        scale_ys=1000,
    )

    # Estimate the posterior with NUTS sampling
    @info "Sample posterior"
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
        "μg";
        posterior_or_prior="posterior",
        data_label=data_label,
        simulation_column="m_br_mean",
        data_column="m_br",
        cred_column_base="m_br",
        rescale_hours_to_mins=true,
        savefile="posterior_brain_mass" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topright,
    )

    plot_simulation_with_data(
        sim,
        "blood [μg/mL]";
        posterior_or_prior="posterior",
        data_label=data_label,
        simulation_column="c_b_mean",
        data_column="c_b",
        cred_column_base="c_b",
        rescale_hours_to_mins=true,
        savefile="posterior_blood_conc" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:topleft,
        scale_ys=1000,
    )

    # animate the main model
    @info "animate solution"
    simani = animate_solution_mass(
        sim.ode_solution,
        sim.posterior_parameters,
        savedir_figures * "animation" * paramcombinationstring * ".gif",
    )

    return sim

end # parseandanalyzePladataset

export parseandanalyzePladataset
