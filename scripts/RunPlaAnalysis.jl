
using Glymphatics1D

sim = parseandanalyzePladataset(
    (estimate_D=true, estimate_Dm=false, estimate_v=false),
);

sim = parseandanalyzePladataset(
    (estimate_D=true, estimate_Dm=true, estimate_v=false),
);

sim = parseandanalyzePladataset(
    (estimate_D=true, estimate_Dm=false, estimate_v=true),
);

parseandanalyzePladataset(
    (estimate_D=true, estimate_Dm=true, estimate_v=true, estimate_all_phys=false);
    #n_samples=100,
);

@time parseandanalyzePladataset(
    (estimate_D=true, estimate_Dm=true, estimate_v=true, estimate_all_phys=true);
    n_samples=1000,
    savedir_results="data/results/Pla/robustness/",
    savedir_figures="visualizations/Pla/robustness/",
    data_label="Pla",
    #turingmodelgenerator=posteriorsamplergeneratorPlaAllparams,
);

savedir_figures="visualizations/Pla/"
savedir_results="data/results/Pla/"
paramcombinationstring = "_D_Dm_v"
T_max = 30
injected_mass = 100
n_samples = 1000

extrapolparseddatafile = "data/parsed/Pla_PEG1kDa_parsed.csv"
FITC_or_Rhodamine = "FITC"
# for extrapolation to FITC-PEG 1 kDa
extrapolsim = initialize_striatal_infusion_normaldist(
    T_max,
    injected_mass;
    animal="mouse",
    infusion_spread_stdev=0.4,
    Turing_space_steps=15,
    fitting_parameters=(estimate_D=true, estimate_Dm=true, estimate_v=true),
    turingmodelgenerator=posteriorsamplergeneratorPla,
)
extrapolsim.posterior_parameters = update_namedtuple(
    sim.posterior_parameters;
    T_max,
    injected_mass,
)
extrapolsim.chains = sim.chains
load_data!(extrapolsim, extrapolparseddatafile)
timedifferential! = enclosethetimedifferential(extrapolsim.prior_parameters)
@unpack D, Dm, v, T_max = extrapolsim.prior_parameters
p = [D, Dm, v]
simulation_prob = ODEProblem(
    timedifferential!,
    extrapolsim.initial_condition,
    (0.0, T_max),
    p,
)

#turing_data = Turing.TArray(Array(convert.(Float32, datapretransform.(simulation.data[1:end, 2:end]))))
turing_data = Turing.TArray(Array(datapretransform.(simulation.data[1:end, 2:end])))
simulation.turing_model = simulation.turingmodelgenerator(
    simulation.data[1:end, 1],
    turing_data,
    simulation_prob,
    simulation.prior_parameters
)
missing_model = extrapolsim.turingmodelgenerator(
    extrapolsim.data[1:end, 1],
    fill(missing, size(extrapolsim.data[1:end, 2:end])),
    simulation_prob,
    extrapolsim.posterior_parameters,
)

sample_posterior!(
    extrapolsim;
    savedir=savedir_results,
    predictions_file=FITC_or_Rhodamine * "_PEG_extrapolation_posterior_predictions" * paramcombinationstring * ".csv",
    chains_file=FITC_or_Rhodamine * "_PEG_extrapolation_chains" * paramcombinationstring * ".jls",
    num_samples=n_samples,
    WAICsfile=FITC_or_Rhodamine * "_PEG_extrapolation_WAICs" * paramcombinationstring * ".csv",
    #loadfile=savedir_results * "chains" * paramcombinationstring * ".jls",
)

# plot priors
plot_simulation_with_data(
    extrapolsim,
    "remaining tracer [%]";
    posterior_or_prior="posterior",
    data_label="Pla $FITC_or_Rhodamine-PEG",
    simulation_column="m_br_mean",
    data_column="m_br",
    cred_column_base="m_br",
    rescale_hours_to_mins=false,
    savefile="$(FITC_or_Rhodamine)_PEG_extrapolation_brain_mass" * paramcombinationstring * ".svg",
    savedir=savedir_figures,
    legend=:topright,
)

extrapolparseddatafile = "data/parsed/Pla_Rhodamine_parsed.csv"
FITC_or_Rhodamine = "Rhodamine"
# for extrapolation to FITC-PEG 1 kDa
extrapolsim = initialize_striatal_infusion_normaldist(
    T_max,
    injected_mass;
    animal="mouse",
    infusion_spread_stdev=0.4,
    Turing_space_steps=15,
    fitting_parameters=(estimate_D=true, estimate_Dm=true, estimate_v=true),
    turingmodelgenerator=posteriorsamplergeneratorPla,
)
extrapolsim.prior_parameters = update_namedtuple(
    sim.posterior_parameters;
    T_max,
    injected_mass,
    D = extrapolsim.prior_parameters.D / 4,
    Dm = extrapolsim.prior_parameters.Dm / 4,
)
extrapolsim.chains = sim.chains
extrapolsim.chains[:, :D, :] = extrapolsim.chains[:, :D, :] ./ 4
extrapolsim.chains[:, :Dm, :] = extrapolsim.chains[:, :Dm, :] ./ 4
load_data!(extrapolsim, extrapolparseddatafile)
sample_posterior!(
    extrapolsim;
    savedir=savedir_results,
    predictions_file=FITC_or_Rhodamine * "_PEG_extrapolation_posterior_predictions" * paramcombinationstring * ".csv",
    chains_file=FITC_or_Rhodamine * "_PEG_extrapolation_chains" * paramcombinationstring * ".jls",
    num_samples=n_samples,
    WAICsfile=FITC_or_Rhodamine * "_PEG_extrapolation_WAICs" * paramcombinationstring * ".csv",
    #loadfile=savedir_results * "chains" * paramcombinationstring * ".jls",
)

# plot priors
plot_simulation_with_data(
    extrapolsim,
    "brain tracer mass [%]";
    posterior_or_prior="posterior",
    data_label="Pla $FITC_or_Rhodamine-PEG",
    simulation_column="m_br_mean",
    data_column="m_br",
    cred_column_base="m_br",
    rescale_hours_to_mins=false,
    savefile="$(FITC_or_Rhodamine)_PEG_extrapolation_brain_mass" * paramcombinationstring * ".svg",
    savedir=savedir_figures,
    legend=:topright,
)
