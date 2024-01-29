
function plot_ode_solution_mass_distribution(
    simulation::GlymphaticSimulation,
    ylabel::String;
    savefile::String="tmp_brain_mass_model.svg",
    savedir::String="../Visualizations/",
    legendloc::Symbol=:right,
)
    t = simulation.ode_solution.t
    (bm, m_v, m_csf, m_b) = state_mass_over_time(
        simulation.ode_solution,
        simulation.prior_parameters,
    )
    total_mass = bm + m_v + m_csf + m_b

    p1 = plot(
        t,
        bm,
        label="brain",
        thickness_scaling=2,
        legend=legendloc,
    )
    plot!(t, m_v, label="ventricle")
    plot!(t, m_csf, label="CM")
    plot!(t, m_b, label="blood")
    plot!(t, total_mass, label="sum")
    xlabel!("hours")
    ylabel!(ylabel)

    display(p1)
    savefig(p1, savedir * savefile)

    return p1
end

function plot_data_and_prediction_dfs(
    data_df::DataFrame,
    solution_df::DataFrame,
    ylabel::String;
    simulation_column::String,
    data_column::String,
    cred_column_base::String="",
    data_label::String="",
    rescale_hours_to_mins::Bool=false,
    savefile::String="tmp_sim_data.svg",
    savedir::String="../Visualizations/",
    legend::Symbol=:bottomright,
    scale_ys::Real=1.0,
    model_label::String="model",
)
    sim_time = solution_df."time (h)"
    data_time = data_df."time (h)"

    if rescale_hours_to_mins
        sim_time = sim_time * 60
        data_time = data_time * 60
        xlabel = "minutes"
    else
        xlabel = "hours"
    end

    yformatter(number) = string(round(number, digits=1))
    p1 = plot(
        thickness_scaling=2,
        data_time,
        data_df[!, data_column] * scale_ys,
        marker=:circle,
        markersize=2,
        label=data_label,
        color="green",
        legend=legend,
        left_margin=-25px,
        bottom_margin=-20px,
        yformatter = yformatter,
    )
    plot!(
        sim_time,
        solution_df[!, simulation_column] * scale_ys,
        label=model_label,
        color="darkblue"
    )

    if any(cred_column_base * "_low-cred" .== names(solution_df))
        lower_col = cred_column_base * "_low-cred"
        upper_col = cred_column_base * "_high-cred"
        keep_rows = [!ismissing(row[lower_col]) for row in eachrow(solution_df)]

        plot!(
            sim_time[keep_rows],
            solution_df[keep_rows, lower_col] * scale_ys,
            fillrange=solution_df[keep_rows, upper_col] * scale_ys,
            alpha=0.1,
            color=:blue,
            label="",
            linewidth=0,
        )
    end

    xlabel!(xlabel)
    ylabel!(ylabel)
    #display(p1)
    savefig(p1, savedir * savefile)

    return p1
end

function plot_simulation_with_data(
    simulation::GlymphaticSimulation,
    ylabel::String;
    posterior_or_prior::String,
    simulation_column::String,
    data_column::String,
    cred_column_base::String="",
    data_label::String="",
    rescale_hours_to_mins::Bool=false,
    savefile::String="tmp_sim_data.svg",
    savedir::String="../Visualizations/",
    legend::Symbol=:bottomright,
    scale_ys::Real=1.0,
    model_label::String="model",
)
    solution_df = select_post_or_prior_predictions(simulation, posterior_or_prior)

    p = plot_data_and_prediction_dfs(
        simulation.data,
        solution_df,
        ylabel;
        simulation_column,
        data_column,
        cred_column_base,
        data_label,
        rescale_hours_to_mins,
        savefile,
        savedir,
        legend,
        scale_ys,
        model_label,
    )

    return p
end

function select_post_or_prior_predictions(
    simulation::GlymphaticSimulation,
    posterior_or_prior::String
)
    if posterior_or_prior == "prior"
        solution_df = simulation.prior_predictions_df
    elseif posterior_or_prior == "posterior"
        solution_df = simulation.posterior_predictions_df
    else
        error("parameter posterior_or_prior must be either posterior or prior")
    end

    return solution_df
end

function add_dev!(WAIC, std, y_index::Int, subplot::Int=1)
    plot!(
        [WAIC - std, WAIC + std],
        y_index * ones(2),
        label="",
        color=:black,
        subplot=subplot,
    )
    scatter!([WAIC], [y_index], label="", color=:black, subplot=subplot)
end

function plot_compare_WAICs!(
    WAICs::Vector{Float64},
    stds::Vector{Float64},
    labels::Array{String, 1};
    subplot::Int=1,
)
    y_labels_indices = Dict()
    for (i, label) in enumerate(labels)
        add_dev!(WAICs[i], stds[i], i, subplot)
        y_labels_indices[i] = label
    end

    ylims!(0.5, length(labels) + .50)
    yticks!(range(1, length=length(labels)), labels)

    min_vals = WAICs .- stds
    max_vals = WAICs .+ stds
    min_val = round(minimum(min_vals), digits=-1)
    max_val = round(maximum(max_vals), digits=-1)
    step_length = round(abs(max_val - min_val) / 2, digits=-1)
    xticks!(min_val:step_length:max_val, subplot=subplot)
    #=if log(maximum(abs.(min_vals))) > 2
        #xticks!(round(minimum(min_vals), digits=-2):100:round(maximum(max_vals), digits=-2))
        xticks!([-400, -200, 0, 200])
    end=#
    xlabel!("deviance", subplot=subplot)
end

function animate_solution_mass(sol, params, saveas)
    T = size(sol)[2]
    J = length(params.r_space)

    if params.T_max > 5
        timepoints = 0.1:0.2:params.T_max
    else
        timepoints = 0.1:0.05:params.T_max
    end

    function label_func(t)
        if t > 2
            return "t = " * string(round(t, digits=1)) * " h"
        end
        return "t = " * string(round(60 * t)) * " min"
    end

    max_tissue_conc = maximum(hcat(sol.u...)[1:end-3, :])

    anim = @animate for ti ∈ timepoints
        p1 = plot(params.r_space, sol(ti)[1:J], title="Brain",
                    xlabel="depth [mm]", ylabel="μg / mm³",
                    titlefontsize=11, tickfontsize=11, guidefontsize=11,
                    label="", ylim=(0, max_tissue_conc), color=:darkgreen,
                    linewidth=2)
        annotate!(([0], 1.05 * max_tissue_conc, text(label_func(ti), 10)))

        p2 = barplot_timepoint(sol, params, ti)
        #annotate!(([2], -0.1 * params.injected_mass, text(label_func(ti), 10)))

        plot(p1, p2, layout = (1, 2))

    end every 1
    gif(anim, saveas, fps=5)
end

function barplot_timepoint(sol, params, ti)
    J = length(params.r_space)
    num_comps = size(sol)[1] - J
    vols = Diagonal([params.V_v, params.V_c, params.V_b])

    if num_comps == 3
        vent_csf_blood = vols * sol(ti)[J + 1: J + 3]
        p2 = bar([1], [vent_csf_blood[1]], title="Vent, CSF & blood",
                    ylabel="μg", labelfontsize=11,
                    titlefontsize=11, tickfontsize=11,
                    label="", ylim=(0, 1.1 * params.injected_mass),
                 c=:darkblue)
        bar!([2], [vent_csf_blood[2]], label="", c=:blue)
        bar!([3], [vent_csf_blood[3]], label="", c=:red)
        xticks!([1.0, 2.0, 3.0], ["vent.", "CSF", "blood"])
    elseif num_comps == 2
        vols = Diagonal([params.V_c, params.V_b])
        csf_blood = vols * sol(ti)[J + 1: J + 2]
        brain_m = slab_int(sol(ti)[1:J], params)

        p2 = bar([2], [csf_blood[1]], title="CSF & blood",
                    ylabel="μg", titlefontsize=10,
                    label="", ylim=(0, 1.1 * params.injected_mass),
                 c=:darkblue)
        bar!([3], [csf_blood[2]], label="", c=:red)
        bar!([1], [brain_m], label="", c=:orange)
        xticks!([1.0, 2.0, 3.0], ["Brain", "CSF", "blood"])
    end
    return p2
end

function plotglymphaxis!(
    datadf,
    predsdf,
    ylabel=nothing;
    simulation_column,
    data_column,
    cred_column_base,
    subplot,
    simcolor=:blue,
    varargs...,
)
    ylabel = isnothing(ylabel) ? data_column : ylabel

    plot!(
        datadf."time (h)",
        datadf[:, data_column],
        subplot=subplot,
        marker=:circle,
        markersize=3,
        markerstrokecolor=:green,
        linewidth=2,
        color=:green,
        label="data",
        legend_background_color=:white,
        legend_foreground_color=:white;
        varargs...
    )
    plot!(
        predsdf."time (h)",
        predsdf[:, simulation_column],
        subplot=subplot,
        label="model",
        color=simcolor,
    )
    lower_col = cred_column_base * "_low-cred"
    upper_col = cred_column_base * "_high-cred"
    keep_rows = [!ismissing(row[lower_col]) for row in eachrow(predsdf)]

    plot!(
        predsdf[keep_rows, "time (h)"],
        Vector{Float64}(predsdf[keep_rows, lower_col]),
        fillrange=Vector{Float64}(predsdf[keep_rows, upper_col]),
        alpha=0.1,
        color=simcolor,
        label="",
        linewidth=0,
        subplot=subplot,
    )
end

function plotGadobutrol(
    datadf,
    predsdf,
    chains,
    gadobutrol_parenchyma_measurement_locations;
    varargs...,
)
    axesgreen = RGB((135, 182, 151) ./ 255...);
    lightbluebackground = RGB((238, 248, 250) ./ 255...);
    par1col = RGB((232, 213, 211) ./ 255...);
    par2col = RGB((233, 221, 213) ./ 255...);
    par3col = RGB((244, 227, 213) ./ 255...);
    par4col = RGB((245, 216, 224) ./ 255...);
    par5col = RGB((237, 217, 233) ./ 255...);
    par6col = RGB((232, 218, 233) ./ 255...);
    par7col = RGB((223, 213, 226) ./ 255...);
    parcolors = [par1col, par2col, par3col, par4col, par5col, par6col, par7col]

    minrow1 = 0.95 * minimum(skipmissing(hcat(
        predsdf."parenchyma 1_low-cred",
        datadf."c_v",
        predsdf."c_c_low-cred"
    )))
    maxrow1 = 1.05 * maximum(skipmissing(hcat(
        predsdf."parenchyma 1_high-cred",
        datadf."c_v",
        predsdf."c_v_high-cred",
        predsdf."c_c_high-cred",
    )))
    parcat = hcat(
        datadf."parenchyma 1",
        datadf."parenchyma 2",
        datadf."parenchyma 3",
        predsdf."parenchyma 1_high-cred",
        predsdf."parenchyma 2_high-cred",
        predsdf."parenchyma 3_high-cred",
        predsdf."parenchyma 1_low-cred",
        predsdf."parenchyma 2_low-cred",
        predsdf."parenchyma 3_low-cred",
    )
    minrow2, maxrow2 = 0.95 * minimum(skipmissing(parcat)), 1.05 * maximum(skipmissing(parcat))
    minrow2, maxrow2 = round.([minrow2, maxrow2], digits=1)
    parcat = hcat(
        datadf."parenchyma 4",
        datadf."parenchyma 5",
        datadf."parenchyma 6",
        datadf."parenchyma 7",
        predsdf."parenchyma 4_high-cred",
        predsdf."parenchyma 5_high-cred",
        predsdf."parenchyma 6_high-cred",
        predsdf."parenchyma 7_high-cred",
        predsdf."parenchyma 4_low-cred",
        predsdf."parenchyma 5_low-cred",
        predsdf."parenchyma 6_low-cred",
        predsdf."parenchyma 7_low-cred",
    )
    minrow3, maxrow3 = 0.95 * minimum(skipmissing(parcat)), 1.05 * maximum(skipmissing(parcat))
    minrow3, maxrow3 = round.([minrow3, maxrow3], digits=1)
    minrow4, maxrow4 = minimum([minrow2, minrow2]), maximum([maxrow2, maxrow2])

    p = plot(
        thickness_scalign=2,
        layout=(4, 4),
        #link=:y,
        size=(1200, (4/3) * 500),
        foreground_color_border=axesgreen,
        legend_foreground_color=:white,
        legend_background_color=RGBA(1, 1, 1, 0.75),
        tickfontsize=9,
        labelfontsize=9,
        legendfontsize=9;
        varargs...,
    )
    plot!(
        subplot=5,
        grid=nothing,
        showaxis=false,
        xticks=nothing,
        yticks=nothing,
    )
    vthirdy = (maxrow1 - minrow1) / 3
    plotglymphaxis!(
        datadf,
        predsdf;
        simulation_column="c_v_mean",
        data_column="c_v",
        cred_column_base="c_v",
        subplot=4,
        ylabel="",
        ylims=(minrow1, maxrow1),
        yticks=round.([minrow1 + vthirdy, maxrow1 - vthirdy], digits=1),
        xlabel="",
        background_color_subplot=lightbluebackground,
        legend=:bottomleft,
    )
    annotate!((0.5, -0.3), text("hours", 9, :center), subplot=4)

    subplotidx = 6:12
    paridx = 1:7
    for (subploti, pari) ∈ zip(subplotidx, paridx)
        simulation_column="parenchyma " * string(pari) * "_mean"
        data_column="parenchyma " * string(pari)
        cred_column_base="parenchyma " * string(pari)
        ylims = subploti > 8 ? (minrow3, maxrow3) : (minrow2, maxrow2)
        thirddiff3 = round((maxrow3 - minrow3) / 3, digits=1)
        thirddiff2 = round((maxrow2 - minrow2) / 3, digits=1)
        yticks = subploti > 8 ? [minrow3 + thirddiff3, maxrow3 - thirddiff3] : [minrow2 + thirddiff2, maxrow2 - thirddiff2]
        legend = subploti == 6 ? :topright : nothing

        plotglymphaxis!(
            datadf,
            predsdf;
            simulation_column,
            data_column,
            cred_column_base,
            subplot=subploti,
            ylims=ylims,
            background_color_subplot=parcolors[pari],
            xlabel="",
            ylabel="",
            yticks=yticks,
            legend=legend,
        )
        annotate!((0.5, -0.3), text("hours", 9, :center), subplot=subploti)
    end

    for i in [4, 6, 9, 13]
        annotate!((-0.22, 0.5), text("MR signal", 9, :center, rotation=90), subplot=i)
    end
    #plot!(subplot=6, ylabel="MR signal")
    #plot!(subplot=9, ylabel="MR signal")

    glymph_density!(
        "D",
        chains,
        "mm² / h",
        "posterior";
        subplot=1,
        background_color_subplot=lightbluebackground,
        yticks=nothing,
        legend=:topleft,
    )
    glymph_density!(
        "D_m",
        chains,
        "10⁻² mm² / h",
        "";
        subplot=2,
        scaleby=100.0,
        background_color_subplot=lightbluebackground,
        yticks=nothing,
    )
    glymph_density!(
        "v",
        chains,
        "mm / h",
        "";
        subplot=3,
        background_color_subplot=lightbluebackground,
        yticks=nothing,
        #xticks=[0.22, 0.24, 0.26],
    )

    for i ∈ 1:12
        if i == 5
            continue
        end
        locs, labs = xticks(p[i])
        if length(locs) % 2 == 1
            mididx = convert(Int, ceil(length(locs) / 2))
            locs = [locs[1], locs[mididx], locs[end]]
            labs = [labs[1], labs[mididx], labs[end]]
            xticks!(locs, labs, subplot=i)
        elseif length(locs) > 4
            locs = locs[1:2:end]
            labs = labs[1:2:end]
            xticks!(locs, labs, subplot=i)
        end
    end

    gadobutrol_parenchyma_measurement_locations = DataFrame(CSV.File(gadobutrol_parenchyma_measurement_locations))
    gadobutrol_parenchyma_measurement_locations = gadobutrol_parenchyma_measurement_locations[1:7, 2]
    t_is = convert.(Int, round.(range(1, stop=size(datadf, 1), length=4), digits=0))
    for (spi, t_i) in enumerate(t_is)
        plotspaceaxis!(
            datadf,
            predsdf,
            datadf[t_i, "time (h)"],
            gadobutrol_parenchyma_measurement_locations;
            subplot=12 + spi,
            #link=:y,
            ylims=(minrow4, maxrow4),
        )
    end

    return p
end

function plotspaceaxis!(
    datadf,
    predsdf,
    timepoint,
    gadobutrol_parenchyma_measurement_locations;
    varargs...,
)
    datatimepointindex = argmax(timepoint .== datadf[1:end, "time (h)"])
    predstimepointindex = argmax(timepoint .== predsdf[1:end, "time (h)"])
    pdf = predsdf[1:end, Cols(r"time", r"par.*mean")]
    pldf = predsdf[1:end, Cols(r"time", r"par.*low-cred")]
    phdf = predsdf[1:end, Cols(r"time", r"par.*high-cred")]
    ddf = datadf[1:end, Cols(r"time", r"par")]
    plot!(
        gadobutrol_parenchyma_measurement_locations,
        Vector(ddf[datatimepointindex, 2:end]),
        marker=:circle,
        markersize=2,
        color=:green,
        label="",
        foreground_color_border=:black,
        xticks=[0, 1, 2],
        xlabel="depth (mm.)";
        varargs...,
    )
    plot!(
        gadobutrol_parenchyma_measurement_locations,
        Vector(pdf[predstimepointindex, 2:end]),
        label="",
        color=1;
        varargs...,
    )

    #lower_col = roiname * "_low-cred"
    #upper_col = roiname * "_high-cred"
    #keep_rows = [!ismissing(row[lower_col]) for row in eachrow(predsdf)]
    plot!(
        gadobutrol_parenchyma_measurement_locations,
        Vector(pldf[predstimepointindex, 2:end]),
        fillrange=Vector(phdf[predstimepointindex, 2:end]),
        #Vector{Float64}(predsdf[keep_rows, lower_col]),
        #fillrange=Vector{Float64}(predsdf[keep_rows, upper_col]),
        alpha=0.1,
        color=:blue,
        label="",
        linewidth=0;
        varargs...,
    )

    timepoint = round(timepoint, digits=1)
    annotate!((0.95, 0.9), text(L"t=%$timepoint h", 9, :right); varargs...,)

end

export plot_ode_solution_mass_distribution, plot_simulation_with_data, plot_compare_WAICs!
export animate_solution_mass
export plot_data_and_prediction_dfs
export plotglymphaxis!
export plotGadobutrol
