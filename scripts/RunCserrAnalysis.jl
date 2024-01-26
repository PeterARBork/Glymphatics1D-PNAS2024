
using Glymphatics1D

parseandanalyzeCserrdataset(
    (estimate_D=true, estimate_Dm=false, estimate_v=false),
)

parseandanalyzeCserrdataset(
    (estimate_D=true, estimate_Dm=true, estimate_v=false),
)

@time parseandanalyzeCserrdataset(
    (estimate_D=true, estimate_Dm=false, estimate_v=true, estimate_all_phys=false),
)

@time parseandanalyzeCserrdataset(
    (estimate_D=true, estimate_Dm=true, estimate_v=true, estimate_all_phys=false),
)

@time parseandanalyzeCserrdataset(
    (estimate_D=true, estimate_Dm=true, estimate_v=true, estimate_all_phys=true),
    n_samples=1000,
    savedir_results="data/results/Cserr/robustness/",
    savedir_figures="visualizations/Cserr/robustness/",
    data_label="Cserr",
    tracer="PEG900",
)

parseandanalyzeCserrdataset(
    (estimate_D=true, estimate_Dm=true, estimate_v=true, estimate_all_phys=false),
    n_samples=1000,
    savedir_results="data/results/Cserr/",
    savedir_figures="visualizations/Cserr/",
    data_label="Cserr",
    tracer="PEG4000",
)

parseandanalyzeCserrdataset(
    (estimate_D=true, estimate_Dm=true, estimate_v=true, estimate_all_phys=false),
    n_samples=1000,
    savedir_results="data/results/Cserr/",
    savedir_figures="visualizations/Cserr/",
    data_label="Cserr",
    tracer="albumin",
)
