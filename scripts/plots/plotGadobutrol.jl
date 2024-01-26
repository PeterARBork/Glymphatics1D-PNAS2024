using Glymphatics1D
using DataFrames, CSV
using Statistics, Distributions
using Plots
using MCMCChains#, Turing
using Plots, Plots.PlotMeasures
using LaTeXStrings
using Serialization


recordings = [
    "210215_mouse03",
    "210314_mouse02",
    "210314_mouse03",
    "210322_mouse04",
    "210330_mouse05",
    "210330_mouse06",
    "210330_mouse07",
]

#for rob in ["", "robustness/"]
rob = ""
for recID in recordings
    p = gadobutrol_overviewplot(
        (estimate_D=true, estimate_Dm=true, estimate_v=true),
        recID,
        "data/parsed/gadobutrol_infusion_parsed.csv",
        "data/results/gadobutrol_infusion/$recID/" * rob,
        "visualizations/gadobutrol_infusion/$recID/" * rob,
        "visualizations/gadobutrol_infusion/$recID/" * rob,
        "data/inputs/gadobutrol intrastriatal/coronal horizontal line segment distances.csv",
        xtickfontvalign=:bottom,
        ytickfonthalign=:left,
        bottom_margin=20px,
        left_margin=10px,
    )
end
#end
display(p)

fitting_parameters = (estimate_D=true, estimate_Dm=true, estimate_v=true)
jointps = plot(layout=(3, 3))
recID = recordings[1]
savedir_results = "data/results/gadobutrol_infusion/$recID/"
paramcombinationstring = makeparamcombinationstring(fitting_parameters)
loadfile = savedir_results * "/chains" * paramcombinationstring * ".jls"
chains = read(loadfile, Chains)

lpp = pointwise_loglikelihoods(model, get_sections(chains, :parameters))


function gadobutrol_overviewplot(
    fitting_parameters=(estimate_D=true, estimate_Dm=true, estimate_v=true),
    recordingID="210215_mouse03",
    parsed_datafile="data/parsed/gadobutrol_infusion_parsed.csv",
    savedir_results="data/results/gadobutrol_infusion/",
    savedir_figures="visualizations/gadobutrol_infusion/",
    savedir_figure_overview = "visualizations/gadobutrol_infusion/",
    gadobutrol_parenchyma_measurement_locations="data/inputs/gadobutrol intrastriatal/coronal horizontal line segment distances.csv";
    varargs...,
)
    paramcombinationstring = makeparamcombinationstring(fitting_parameters)
    parsed_datafile = replace(parsed_datafile, "gadobutrol_infusion"=>recordingID)
    savedir_figure_overview = savedir_figure_overview * "/"
    loadfile = savedir_results * "/chains" * paramcombinationstring * ".jls"

    predsdf = DataFrame(CSV.File(savedir_results * "/posterior_predictions" * paramcombinationstring * ".csv"))
    datadf = DataFrame(CSV.File(parsed_datafile))
    chains = deserialize(loadfile)

    predsdf[!, Not("time (h)")] = predsdf[:, Not("time (h)")] ./ 100
    datadf[!, Not("time (h)")] = datadf[:, Not("time (h)")] ./ 100

    p = plotGadobutrol(
        datadf,
        predsdf,
        chains,
        gadobutrol_parenchyma_measurement_locations;
        varargs...,
    )
    annotate!((0.05, 0.9), text("Parenchyma 1", 9, :left), subplot=6,)
    for pari in range(2, 7)
        annotate!((0.95, 0.9), text("Parenchyma $pari", 9, :right), subplot=pari + 5,)
    end

    D = truncated(LogNormal(log(0.46), 0.5^2), 0.01, 20)
    plot!(
        D,
        subplot=1,
        label="prior",
        color=:black,
        xlims=(0.2, 0.7),
        yticks=nothing,
        xticks=[0.2, 0.4, 0.6],
    )
    plot!(
        truncated(Normal(0.0, 1.0), -2, 2),
        subplot=3,
        label="prior",
        color=:black,
    )
    _, xl = xlims(p[2])
    plot!(
        range(0, xl, length=5), ones(5), label="prior", color=:black,
    )

    savefig(p, savedir_figure_overview * "$recordingID overview.pdf")

    jointp = makemarginalkdeplot(
        chains[:, :v, :],
        chains[:, :D_m, :],
        savedir_figures;
        scale_b=100,
        #left_margin=0px,
        #bottom_margin=0px,
        clip=((-3.0, 2.5), (-3.5, 3)),
        saveas="joint posterior v Dm.pdf",
    )
    display(jointp)

    return p
end
