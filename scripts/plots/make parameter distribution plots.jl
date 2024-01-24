using Glymphatics1D

using Distributions
using Plots, StatsPlots
using Plots.PlotMeasures
using Turing

results_dir = "data/results/"
visualizations_dir = "visualizations/"

D = truncated(LogNormal(log(0.4), 0.5^2), 0.1, 20)

Cserr_chains = deserialize(results_dir * "Cserr/chains_D_Dm_v.jls")
Pla_chains = deserialize(results_dir * "Pla/chains_D_Dm_v.jls")

p1 = makemarginalkdeplot(
    Pla_chains[:, :v, :],
    Pla_chains[:, :Dm, :],
    visualizations_dir;
    clip=((-3.1, 3), (-2.5, 2.5)),
    scale_b=100,
    saveas="Pla/joint posterior v Dm.pdf",
)

p4 = makemarginalkdeplot(
    Cserr_chains[:, :v, :],
    Cserr_chains[:, :Dm, :],
    visualizations_dir;
    scale_b=100,
    #left_margin=0px,
    #bottom_margin=0px,
    clip=((-1.5, 3), (-3.5, 3)),
    saveas="Cserr/joint posterior v Dm.pdf",
)
