
function makemarginalkdeplot(
    chain_a,
    chain_b,
    visualizations_dir;
    clip=((-3, 3), (-3, 3)),
    scale_b=1,
    saveas="",
    subplot=1,
)
    p1 = StatsPlots.marginalkde(
        chain_a,
        scale_b * chain_b,
        linewidth=2,
        grid=false,
        clip=clip,#((-1.5, 3), (-3.5, 3))
        subplot=subplot,
    )
    #xlims!(-1.0, 0.6)
    #ylims!(0.0, 0.003)
    plot!(
        size=(600, 500),
        tickfontsize=12,
        labelfontsize=14,
        yformatter=y->string(round(y, digits=2)),
        #left_margin=left_margin,
        #bottom_margin=bottom_margin,
    )
    ticks = xticks(p1)
    locs, labels = ticks[2]
    xticks!(p1[2], locs[1:2:end-1], labels[1:2:end-1])
    ticks = yticks(p1)
    locs, labels = ticks[2]
    yticks!(p1[2], locs[1:2:end], labels[1:2:end])

    orderofmagnitude = log10(scale_b)
    if orderofmagnitude == 0.0
        scalestring = ""
    elseif orderofmagnitude == 1.0
        scalestring = "10¹ "
    elseif orderofmagnitude == 2.0
        scalestring = "10² "
    else
        scalestring = "10^" * string(orderofmagnitude) * " "
    end

    xlabel!(p1[2], "v [mm / h]")
    ylabel!(p1[2], "Dm [" * scalestring * "mm² / h]")
    ylabel!(p1[1], "Pr(v)")
    xlabel!(p1[3], "Pr(Dm)")

    if length(saveas) > 1
        savefig(visualizations_dir * saveas)
    end

    return p1
end

function glymph_density!(
    param::String,
    chains::Chains,
    units::String,
    label::String;
    subplot::Int=1,
    scaleby::Float64=1.0,
    colorindex::Union{Missing, Integer, RGBA}=missing,
    xlabeldistance=0.22,
    varargs...,
)
    if !ismissing(colorindex)
        pd = Dict(:color=>colorindex)
    else
        pd = Dict()
    end

    density!(
        scaleby .* vec(chains[Symbol(param)]),
        label=label,
        xlabel="",
        ylabel="",
        subplot=subplot,
        linewidth=2,
        #left_margin=-10px,
        legend_foreground_color=:white;
        pd...,
        varargs...
    )
    param = replace(param, "_m" => "ₘ")
    xlabel=param * " [" * units * "]"
    annotate!((0.5, -xlabeldistance), text(xlabel, 9, :center), subplot=subplot)
    ylabel="Pr(" * param * ")"
    annotate!((-0.1, 0.5), text(ylabel, 9, :center, rotation=90), subplot=subplot)
end

function plotDposteriorfourdatasets!(
    Lee_chains::Chains,
    Stanton_chains::Chains,
    Cserr_chains::Chains,
    Pla_chains::Chains,
    subplot::Integer=1,
)

    D = truncated(LogNormal(log(0.4), 0.5^2), 0.1, 20)

    glymph_density!("D", Lee_chains, "mm² / h", "Lee"; subplot)
    glymph_density!("D", Stanton_chains, "mm² / h", "Stanton"; subplot)
    glymph_density!("D", Cserr_chains, "mm² / h", "Cserr"; subplot)
    glymph_density!("D", Pla_chains, "mm² / h", "Pla"; subplot)

    p = plot!(
        D,
        label="Prior",
        color="black",
        subplot=subplot,
        xlims=(0.0, 2.0),
    )
    #xlims!(0.0, 2)

    plot!(
        legend=:topright,
        #bottom_margin=[0mm 0mm],
        #left_margin=[0mm 0mm],
        #margin=0px,
        #bbox_inches="tight",
        yticks=false,
        subplot=subplot,
    )
    #savefig(visualizations_dir * "posterior D.pdf")
    return p

end

function plotposteriorparametersforfourtrainingsets(
    Lee_chains::Chains,
    Stanton_chains::Chains,
    Cserr_chains::Chains,
    Pla_chains::Chains,
    visualizations_dir,
)
    p = plot(
        layout=(4, 1),
        size=(600, 1200),
        yticks=false,
        legendfontsize=8,
        legend_foreground_color=:white,
        thickness_scaling=2,
        #bottom_margin=-30px,
        #left_margin=-5px,
        labelfontsize=10,
        tickfontsize=8,
    )
    plotDposteriorfourdatasets!(
        Lee_chains,
        Stanton_chains,
        Cserr_chains,
        Pla_chains,
        1,
    )
    #savefig(visualizations_dir * "posterior D.pdf")

    glymph_density!("Dm", Lee_chains, "mm²/h", ""; subplot=2)
    glymph_density!("Dm", Stanton_chains, "mm²/h", ""; subplot=2)
    #plot!(legend=:topleft, subplot=2)
    #savefig(visualizations_dir * "posterior Dm Lee and Stanton.pdf")
    glymph_density!("Dm", Cserr_chains, "10⁻³ mm²/h", ""; subplot=3, scaleby=1000.0, colorindex=3)
    glymph_density!("Dm", Pla_chains, "10⁻³ mm²/h", ""; subplot=3, scaleby=1000.0, colorindex=4)
    #savefig(visualizations_dir * "posterior Dm Cserr and Pla.pdf")

    glymph_density!("v", Lee_chains, "mm / h", ""; subplot=4)
    glymph_density!("v", Stanton_chains, "mm / h", ""; subplot=4)
    glymph_density!("v", Cserr_chains, "mm / h", ""; subplot=4)
    glymph_density!("v", Pla_chains, "mm / h", ""; subplot=4)
    plot!(legend=:top, subplot=4)
    #savefig(visualizations_dir * "posterior v all.pdf")
    savefig(visualizations_dir * "posterior parameter distributions for four training sets.pdf")

    return p
end

export makemarginalkdeplot, glymph_density!
export plotDposteriorfourdatasets!
export plotposteriorparametersforfourtrainingsets
