function parseandanalyzeGadobutrolinfusion(
    fitting_parameters::NamedTuple;
    n_samples::Int64=1000,
    Turing_space_steps::Int64=15,
    raw_datafile="data/inputs/gadobutrol intrastriatal/coronal_May21_TACs.csv",
    gadobutrol_parenchyma_measurement_locations="data/inputs/gadobutrol intrastriatal/coronal horizontal line segment distances.csv",
    parsed_datafile="data/parsed/gadobutrol_infusion_parsed.csv",
    savedir_results="data/results/gadobutrol_infusion/",
    savedir_figures="visualizations/gadobutrol_infusion/",
    data_label="gadobutrol_infusion",
)
    animalID = join(split(split(raw_datafile, "/")[end], "_")[1:end-1], "_")
    parsed_datafile = replace(parsed_datafile, "gadobutrol_infusion" => animalID)
    if fitting_parameters.estimate_all_phys
        savedir_results = savedir_results * animalID * "/robustness/"
        savedir_figures = savedir_figures * animalID * "/robustness/"
    else
        savedir_results = savedir_results * animalID * "/"
        savedir_figures = savedir_figures * animalID * "/"
    end
    @info "animalID = $animalID"
    @info "parsed_datafile = $parsed_datafile"
    @info "savedir_results = $savedir_results"
    @info "savedir_figures = $savedir_figures"

    @assert(:estimate_D in keys(fitting_parameters))
    @assert(:estimate_Dm in keys(fitting_parameters))
    @assert(:estimate_v in keys(fitting_parameters))

    paramcombinationstring = makeparamcombinationstring(fitting_parameters)

    # parse the data
    println("Parsing data.")
    gadobutrol_infusion_df = parsegadobutrolinfusion(
        raw_datafile,
        parsed_datafile;
        signalcolbasename="DCE_",
        framesperhour=10,
        baselinesubtract=false,
        normalizedf=true,
        normcolumn="c_c",
        croppump_finish_time=0.0,
        brainconcstring="c(0.25mm)",
        numparenchymaroi=7,
    )
    gadobutrol_parenchyma_measurement_locations = DataFrame(CSV.File(gadobutrol_parenchyma_measurement_locations))
    gadobutrol_parenchyma_measurement_locations = gadobutrol_parenchyma_measurement_locations[1:7, 2]

    T_max = 1.1 * maximum(gadobutrol_infusion_df[:, "time (h)"])

    # setup the parameter-fitting combination for simulations
    println("Setting up simulations to fit each parameter-combination.")
    ic_measurements = Vector(gadobutrol_infusion_df[1, Cols(contains("parenchyma"))])
    sim = initialize_striatal_infusion_with_initial_condition(
        T_max,
        ic_measurements,
        gadobutrol_parenchyma_measurement_locations,
        gadobutrol_infusion_df[1, :c_v],
        gadobutrol_infusion_df[1, :c_c],
        0.0;
        paramsfile="data/inputs/ratparametersfromtables.csv",
        Turing_space_steps=Turing_space_steps,
        fitting_parameters=fitting_parameters,
        turingmodelgenerator=posteriorsamplergeneratorGadobutrolinfusion,
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
        "infused mass";
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
    plot!(
        gadobutrol_parenchyma_measurement_locations,
        ic_measurements,
        color="green",
        marker=:circle,
        label="",
    )
    scatter!(
        [sim.prior_parameters.r_space[1], sim.prior_parameters.r_space[end]],
        [gadobutrol_infusion_df[1, :c_v], gadobutrol_infusion_df[1, :c_c]],
        color="green",
        marker=:circle,
        label="",
    )
    scatter!(
        [sim.prior_parameters.r_space[1], sim.prior_parameters.r_space[end]],
        [sim.initial_condition[end-2], sim.initial_condition[end-1]],
        color="blue",
        marker=:cross,
        label="",
    )
    savefig(
        p1,
        savedir_figures * "gadobutrol_infusion_initial_condition" * paramcombinationstring * ".svg",
    )

    plot_simulation_with_data(
        sim,
        "conc. CM";
        posterior_or_prior="prior",
        data_label=data_label,
        simulation_column="c_c_mean",
        data_column="c_c",
        cred_column_base="c_c",
        rescale_hours_to_mins=false,
        savefile="prior_CM_conc" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:none,
    )

    plot_simulation_with_data(
        sim,
        "conc. ventricle";
        posterior_or_prior="prior",
        data_label=data_label,
        simulation_column="c_v_mean",
        data_column="c_v",
        cred_column_base="c_v",
        rescale_hours_to_mins=false,
        savefile="prior_ventricle_conc" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:none,
    )

    plot_simulation_with_data(
        sim,
        "paren. 1";
        posterior_or_prior="prior",
        data_label=data_label,
        simulation_column="parenchyma 1_mean",
        data_column="parenchyma 1",
        cred_column_base="parenchyma 1",
        rescale_hours_to_mins=false,
        savefile="prior_parenchyma 1" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:none,
    )

    # Estimate the posterior with NUTS sampling
    println("Sample posterior")
    sample_posterior!(
        sim;
        savedir=savedir_results,
        predictions_file="posterior_predictions" * paramcombinationstring * ".csv",
        chains_file="chains" * paramcombinationstring * ".jls",
        num_samples=n_samples,
        #loadfile=savedir_results * "chains" * paramcombinationstring * ".jls",
    )

    # plot the chains
    p = plot(sim.chains)
    savefig(p, savedir_figures * "chains_plot" * paramcombinationstring * ".svg")

    # plot the posteriors of the main models
    plot_simulation_with_data(
        sim,
        "conc. CM";
        posterior_or_prior="posterior",
        data_label=data_label,
        simulation_column="c_c_mean",
        data_column="c_c",
        cred_column_base="c_c",
        rescale_hours_to_mins=false,
        savefile="posterior_CM_conc" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:none,
    )

    plot_simulation_with_data(
        sim,
        "conc. ventricle";
        posterior_or_prior="posterior",
        data_label=data_label,
        simulation_column="c_v_mean",
        data_column="c_v",
        cred_column_base="c_v",
        rescale_hours_to_mins=false,
        savefile="posterior_ventricle_conc" * paramcombinationstring * ".svg",
        savedir=savedir_figures,
        legend=:none,
    )

    for i ∈ 1:7
        plot_simulation_with_data(
            sim,
            "paren. " * string(i);
            posterior_or_prior="posterior",
            data_label=data_label,
            simulation_column="parenchyma " * string(i) * "_mean",
            data_column="parenchyma " * string(i),
            cred_column_base="parenchyma " * string(i),
            rescale_hours_to_mins=false,
            savefile="posterior_parenchyma " * string(i) * paramcombinationstring * ".svg",
            savedir=savedir_figures,
            legend=:none,
        )
    end

    # animate the main model
    println("animate solution")
    simani = animate_solution_mass(
        sim.ode_solution,
        sim.posterior_parameters,
        savedir_figures * "animation" * paramcombinationstring * ".gif",
    )

end # parseandanalyzeGadobutrolinfusion

export parseandanalyzeGadobutrolinfusion
