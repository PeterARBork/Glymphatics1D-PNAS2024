using Glymphatics1D
using DataFrames, CSV
using Statistics, Distributions
using Plots
using MCMCChains#, Turing
using Plots, Plots.PlotMeasures
using LaTeXStrings
using Serialization

fitting_parameters=(estimate_D=true, estimate_Dm=true, estimate_v=true)
paramcombinationstring = makeparamcombinationstring(fitting_parameters)
savedir_figure_overview = "visualizations/"
parseddir = "data/parsed/"
resultsdir = "data/results/"


Cserr_chains = deserialize(resultsdir * "Cserr/chains$paramcombinationstring.jls")
Pla_chains = deserialize(resultsdir * "Pla/chains$paramcombinationstring.jls")

Pladatadf = DataFrame(CSV.File(parseddir * "Pla_parsed.csv"))
Cserrdatadf = DataFrame(CSV.File(parseddir * "Cserr_parsed.csv"))

Plapredsdf = DataFrame(CSV.File(resultsdir * "Pla/posterior_predictions$paramcombinationstring.csv"))
Cserrpredsdf = DataFrame(CSV.File(resultsdir * "Cserr/posterior_predictions$paramcombinationstring.csv"))

scbdarkorange = palette(:seaborn_colorblind)[4];
darkorange = RGBA(scbdarkorange.r, scbdarkorange.g, scbdarkorange.b, 1.0);
darkorangealphahalf = RGBA(scbdarkorange.r, scbdarkorange.g, scbdarkorange.b, 0.250);

scblightorange = palette(:seaborn_colorblind)[2];
lightorange = RGBA(scblightorange.r, scblightorange.g, scblightorange.b, 1.0);
lightorangealphahalf = RGBA(scblightorange.r, scblightorange.g, scblightorange.b, 0.250);
background_color_Cserr = RGB(235/255,242/255,232/255);
background_color_Pla = RGB((249, 245, 225) ./ 255...);

begin
    gr()
    p = plot(
        thickness_scalign=2,
        layout=(2, 4),
        #link=:y,
        size=(1200, 600),
        tickfontsize=10,
        labelfontsize=12,
        legendfontsize=10,
        bottom_margin=20px,
        left_margin=25px,
    )

    sp = 1
    plot!(
        subplot=sp,
        grid=nothing,
        showaxis=false,
        xticks=nothing,
        yticks=nothing,
        #background_color_subplot=darkorange,
    )
    annotate!(0.40, 0.50, "experimental\ncartoon", subplot=sp)
    annotate!(-0.3 * xlims(p[sp])[2], 0.95 * ylims(p[sp])[2], "a", subplot=sp)

    sp = 7
    pPladatadf = vcat(Pladatadf, 0.99999 .* Pladatadf[:, :])
    pPlapredsdf = vcat(Plapredsdf, 0.99999 .* Plapredsdf[:, :])
    plotglymphaxis!(
        pPladatadf,
        pPlapredsdf,
        ylabel="brain mass [μg]";
        simulation_column="m_br_mean",
        data_column="m_br",
        cred_column_base="m_br",
        subplot=sp,
        legend=false,
        background_color_subplot=background_color_Pla,
        ylims=(-1, 21),
    )
    xticks!(p[sp], [0, 1, 2], ["0", "1", "2"])
    #yticks!(p[sp], [15, 20], ["15", "20"])
    annotate!(-0.25 * xlims(p[sp])[2], 1.05 * ylims(p[sp])[2], "d", subplot=sp)

    sp = 8
    sPladatadf = copy(Pladatadf)
    sPlapredsdf = copy(Plapredsdf)
    sPladatadf[!, Not("time (h)")] .*= 1000
    sPlapredsdf[!, Not("time (h)")] .*= 1000
    plotglymphaxis!(
        sPladatadf,
        sPlapredsdf,
        ylabel="blood conc. [μg/mL]";
        simulation_column="c_b_mean",
        data_column="c_b",
        cred_column_base="c_b",
        subplot=sp,
        legend=false,
        background_color_subplot=background_color_Pla,
        #ylims=(minrow1, maxrow1)
    )
    xticks!(p[sp], [0, 1, 2], ["0", "1", "2"])
    ylims!(p[sp], -0.10, 1.90)

    sp = 5
    plotglymphaxis!(
        Cserrdatadf,
        Cserrpredsdf;
        simulation_column="m_br_mean",
        data_column="m_br",
        cred_column_base="m_br",
        subplot=sp,
        legend=true,
        ylabel="brain tracer mass [%]",
        background_color_subplot=background_color_Cserr,
        ylims=(-5, 105),
        #legend=:bottomright,
    )
    xticks!(p[sp], [0, 10, 20], ["0", "10", "20"])
    annotate!(-0.35 * xlims(p[sp])[2], 1.05 * ylims(p[sp])[2], "c", subplot=sp)

    sp = 6
    plotglymphaxis!(
        Cserrdatadf,
        Cserrpredsdf;
        simulation_column="m_c_mean",
        data_column="m_c",
        cred_column_base="m_c",
        subplot=sp,
        legend=false,
        ylabel="CSF tracer mass [%]",
        background_color_subplot=background_color_Cserr,
        #ylims=(minrow1, maxrow1),
        #legend=:bottomright,
    )
    xticks!(p[sp], [0, 10, 20], ["0", "10", "20"])
    yticks!(p[sp], [0.0, 0.5, 1.0, 1.5], ["0.0", "0.5", "1.0", "1.5"])
    #yticks!(p[4], [0.0, 2.0, 4.0], ["0.0", "2.0", "4.0"])
    ylims!(p[sp], -0.10, 1.90)

    for i ∈ 5:8
        xlabel!("hours", subplot=i)
    end

    sp=2
    glymph_density!("D", Cserr_chains, "mm² / h", "PEG-900"; subplot=sp, colorindex=3, xlabeldistance=0.16)
    glymph_density!("D", Pla_chains, "mm² / h", "Pla"; subplot=sp, colorindex=lightorange, xlabeldistance=0.16)
    D = truncated(LogNormal(log(0.46), 0.5^2), 0.1, 20)
    p = plot!(
        D,
        label="Prior",
        color="black",
        subplot=sp,
        xlims=(0.0, 2.0),
        yticks=nothing,
    )
    annotate!(-0.2 * xlims(p[sp])[2], 0.95 * ylims(p[sp])[2], "b", subplot=sp)

    sp = 3
    glymph_density!("D_m", Cserr_chains, "10⁻² mm² / h", "PEG-900"; subplot=sp, scaleby=100.0, colorindex=3, xlabeldistance=0.16)
    glymph_density!("D_m", Pla_chains, "10⁻² mm² / h", "Pla"; subplot=sp, scaleby=100.0, colorindex=lightorange, xlabeldistance=0.16)
    plot!(yticks=nothing, subplot=sp)
    _, xl = xlims(p[sp])
    plot!(range(0, xl, length=30), ones(30), subplot=sp, label="Prior", color="black")
    #annotate!(-0.25 * xlims(p[sp])[2], 0.95 * ylims(p[sp])[2], "g", subplot=sp)

    sp = 4
    if :v in names(Cserr_chains)
        glymph_density!("v", Cserr_chains, "mm / h", "PEG-900"; subplot=sp, colorindex=3, xlabeldistance=0.16)
        glymph_density!("v", Pla_chains, "mm / h", "Pla"; subplot=sp, colorindex=lightorange, xlabeldistance=0.16)
        v = truncated(Normal(0.0, 1.0), -3, 3)
        plot!(v, label="Prior", color="black", subplot=sp, yticks=nothing)
        plot!(yticks=nothing, subplot=sp, legend=:topleft)
        #annotate!(-0.55 * (xlims(p[sp])[2] - xlims(p[sp])[1]), 0.95 * ylims(p[sp])[2], "h", subplot=sp)
    else
        plot!(
            subplot=sp,
            grid=nothing,
            showaxis=false,
            xticks=nothing,
            yticks=nothing,
            #background_color_subplot=darkorange,
        )
        #annotate!(0.40, 0.50, "experimental\ncartoon", :color, subplot=sp)
        #annotate!(-0.2 * xlims(p[sp])[2], 0.95 * ylims(p[sp])[2], "a", subplot=sp)
    end
end

savefig(p, savedir_figure_overview * "efflux post pred and params $paramcombinationstring.pdf")
savefig(p, savedir_figure_overview * "efflux post pred and params $paramcombinationstring.svg")

p = plot(
    60 * Pladatadf."time (h)",
    1000 * Pladatadf."c_b",
    color=palette(:Blues)[6],
    marker=:circle,
    markerstrokecolor=:navyblue,
    label="",
    xlab="minutes since injection",
    ylab="DB53 blood signal",
    #left_margin=-20px,
    #bottom_margin=-20px,
    thickness_scaling=2,
    labelfontsize=8,
    tickfontsize=6,
    grid=false,
    background_color_outside=RGBA(1.0, 1.0, 1.0, 0.0),
)

savefig("visualizations/Pla/Pla_raw.svg")
savefig("visualizations/Pla/Pla_raw.pdf")

#pyplot()
p = plot(
    Cserrdatadf[:, "time (h)"],
    Cserrdatadf[:, "m_br"],
    color=RGB((202, 136, 148) ./ 255...),
    marker=:circle,
    markerstrokecolor=RGB((202, 136, 148) ./ 255...),
    markercolor=RGB((231, 193, 200) ./ 255...),
    label="Brain",
    xlab="hours since injection",
    ylab="remaining tracer [%]",
    ylims=(0.1, 150),
    xlims=(-1, 32),
    left_margin=0px,
    bottom_margin=0px,
    yscale=:log10,
    thickness_scaling=2,
    labelfontsize=8,
    tickfontsize=6,
    grid=false,
    legendforegroundcolor=:white,
    legendfontsize=6,
    background_color_outside=RGBA(1.0, 1.0, 1.0, 0.0),
)
plot!(
    Cserrdatadf[2:end, "time (h)"],
    Cserrdatadf[2:end, "m_c"],
    color=RGB((122, 196, 193) ./ 255...),
    markerstrokecolor=RGB((122, 196, 193) ./ 255...),
    markercolor=RGB((180, 219, 219) ./ 255...),
    marker=:circle,
    label="CSF",
)
yt = [1, 10, 100]
yticks!(yt, [string(t) for t in yt])
xt = [0, 10, 20, 30]
xticks!(xt, [string(t) for t in xt])

savefig("visualizations/Cserr/Cserr_raw.svg")
savefig("visualizations/Cserr/Cserr_raw.pdf")


Cserr_PEG4000_chains = deserialize(resultsdir * "Cserr/PEG4000/chains$paramcombinationstring.jls")
Cserr_albumin_chains = deserialize(resultsdir * "Cserr/albumin/chains$paramcombinationstring.jls")
CserrdatadfPEG4000 = DataFrame(CSV.File(parseddir * "Cserr_PEG4000_parsed.csv"))
Cserrdatadfalbumin = DataFrame(CSV.File(parseddir * "Cserr_albumin_parsed.csv"))

CserrpredsdfPEG4000 = DataFrame(CSV.File(resultsdir * "Cserr/PEG4000/posterior_predictions$paramcombinationstring.csv"))
Cserrpredsdfalbumin = DataFrame(CSV.File(resultsdir * "Cserr/albumin/posterior_predictions$paramcombinationstring.csv"))

D_PEG900 = truncated(LogNormal(log(0.46), 0.5^2), 0.1, 20)
D_PEG4000 = truncated(LogNormal(log(0.46 * (900/4000)^(1/3)), 0.5^2), 0.1, 20)
D_albumin = truncated(LogNormal(log(0.46 * (900/69000)^(1/3)), 0.5^2), 0.1, 20)
v_prior = truncated(Normal(0.0, 1.0), -2.0, 2.0)

PEG900_color = palette(:Blues_9)[end-4];
PEG4000_color = palette(:Blues_9)[end-2];
albumin_color = palette(:Blues_9)[end];
background_color_Cserr = RGB(232/255,244/255,230/255);

p = plot(
    layout=(5, 3),
    size=(600, 800),
    xticks=[0, 10, 20],
    titlefontsize=10,
    background_color_subplot=background_color_Cserr,
    tickfontsize=9,
    labelfontsize=9,
    left_margin=15px,
    legendfontsize=8,foreground_color_border=background_color_Cserr,
    legend_foreground_color=:white,
    legend_background_color=RGBA(1, 1, 1, 0.75),
)
plotglymphaxis!(
    Cserrdatadf,
    Cserrpredsdf;
    simulation_column="m_br_mean",
    data_column="m_br",
    cred_column_base="m_br",
    ylabel="Brain mass [%]",
    title="PEG900 (0.9 kDa)",
    subplot=1,
    simcolor=PEG900_color,
    titlefontcolor=PEG900_color,
    legend=:topright,
)
plotglymphaxis!(
    Cserrdatadf,
    Cserrpredsdf;
    simulation_column="m_c_mean",
    data_column="m_c",
    cred_column_base="m_c",
    ylabel="CSF mass [%]",
    xlabel="hours",
    subplot=4,
    legend=false,
    simcolor=PEG900_color,
    titlefontcolor=PEG900_color,
)
plotglymphaxis!(
    CserrdatadfPEG4000,
    CserrpredsdfPEG4000;
    simulation_column="m_br_mean",
    data_column="m_br",
    cred_column_base="m_br",
    ylabel="",
    title="PEG4000 (4 kDa)",
    subplot=2,
    simcolor=PEG4000_color,
    titlefontcolor=PEG4000_color,
)
plotglymphaxis!(
    CserrdatadfPEG4000,
    CserrpredsdfPEG4000;
    simulation_column="m_c_mean",
    data_column="m_c",
    cred_column_base="m_c",
    ylabel="",
    xlabel="hours",
    legend=false,
    subplot=5,
    simcolor=PEG4000_color,
    titlefontcolor=PEG4000_color,
)
plotglymphaxis!(
    Cserrdatadfalbumin,
    Cserrpredsdfalbumin;
    simulation_column="m_br_mean",
    data_column="m_br",
    cred_column_base="m_br",
    ylabel="",
    title="Albumin (69 kDa)",
    subplot=3,
    simcolor=albumin_color,
    titlefontcolor=albumin_color,
)
plotglymphaxis!(
    Cserrdatadfalbumin,
    Cserrpredsdfalbumin;
    simulation_column="m_c_mean",
    data_column="m_c",
    cred_column_base="m_c",
    ylabel="",
    xlabel="hours",
    subplot=6,
    legend=false,
    simcolor=albumin_color,
    titlefontcolor=albumin_color,
)

for i ∈ 1:3
    plot!(ylims=(10, 105), subplot=i)#yticks=[0.0, 0.5, 1.0, 1.5])
end
for i ∈ 4:6
    plot!(ylims=(-0.2, 2.0), subplot=i)
end
#display(p)

#=p = plot(
    layout=(3, 3),
)=#

glymph_density!("D", Cserr_chains, "mm² / h", ""; xlabeldistance=0.25, subplot=7, color=PEG900_color, legend=false,)
glymph_density!("D", Cserr_PEG4000_chains, "mm² / h", ""; xlabeldistance=0.25, subplot=8, color=PEG4000_color, legend=false,)
glymph_density!("D", Cserr_albumin_chains, "mm² / h", ""; xlabeldistance=0.25, subplot=9, color=albumin_color, legend=false,)
p = plot!(
    D_PEG900,
    label="",
    color="gray",
    subplot=7,
    xlims=(0.0, 1.0),
    yticks=nothing,
)
plot!(
    D_PEG4000,
    label="",
    color="gray",
    subplot=8,
    ylabel="",
    yticks=nothing,
    xlims=(0.0, 1.0),
)
plot!(
    D_albumin,
    label="",
    ylabel="",
    color="gray",
    subplot=9,
    yticks=nothing,
    xlims=(0.0, 1.0),
)
plot!(xticks=[0.0, 0.3, 0.6, 0.9], subplot=7)
plot!(xticks=[0.0, 0.3, 0.6, 0.9], subplot=8)
plot!(xticks=[0.0, 0.3, 0.6, 0.9], subplot=9)

glymph_density!("D_m", Cserr_chains, "10⁻² mm² / h", ""; xlabeldistance=0.25, subplot=10, scaleby=100.0, xlims=(0.0, 0.1), yticks=nothing, color=PEG900_color)
glymph_density!("D_m", Cserr_PEG4000_chains, "10⁻² mm² / h", ""; xlabeldistance=0.25, subplot=11, scaleby=100.0, xlims=(0.0, 0.1), yticks=nothing, color=PEG4000_color)
glymph_density!("D_m", Cserr_albumin_chains, "10⁻² mm² / h", ""; xlabeldistance=0.25, subplot=12, scaleby=100.0, xlims=(0.0, 1), yticks=nothing, color=albumin_color)
for i in 10:12
    plot!(xticks=range(xlims(p[i])..., length=3), subplot=i)
end


glymph_density!("v", Cserr_chains, "mm / h", ""; xlabeldistance=0.25, subplot=13, color=PEG900_color)
glymph_density!("v", Cserr_PEG4000_chains, "mm / h", ""; xlabeldistance=0.25, subplot=14, color=PEG4000_color)
glymph_density!("v", Cserr_albumin_chains, "mm / h", ""; xlabeldistance=0.25, subplot=15, color=albumin_color)
for i ∈ [13, 14, 15]
    plot!(
        v_prior,
        label="",
        color="gray",
        subplot=i,
        xlims=(-1, 1),
        xticks=[-0.5, 0.0, 0.5],
        yticks=nothing,
    )
end
display(p)

savefig(p, savedir_figure_overview * "Cserr/postparams_preds_allweights.svg")
savefig(p, savedir_figure_overview * "Cserr/postparams_preds_allweights.pdf")

p1 = makemarginalkdeplot(
    Pla_chains[:, :v, :],
    Pla_chains[:, :D_m, :],
    visualizations_dir;
    clip=((-3.1, 3), (-2.5, 2.5)),
    scale_b=100,
    saveas="Pla/joint posterior v Dm.pdf",
)

p4 = makemarginalkdeplot(
    Cserr_chains[:, :v, :],
    Cserr_chains[:, :D_m, :],
    visualizations_dir;
    scale_b=100,
    clip=((-1.5, 3), (-3.5, 3)),
    saveas="Cserr/joint posterior v Dm.pdf",
)
