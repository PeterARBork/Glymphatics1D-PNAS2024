
function load_data!(simulation::GlymphaticSimulation, source_file::String)
    simulation.data = DataFrame(CSV.File(source_file))

    return simulation.data
end

function parse_Lee_data(source_file::String, parsed_datafile::String)
    rawdf = DataFrame(CSV.File(source_file))

    GdDOTA_g_per_mol = 558.65
    GdDOTA_ug_per_umol = GdDOTA_g_per_mol

    Lee_df = rawdf[rawdf.source .== "Lee", ["time (h)", "m_br"]]
    m_br_umol = Lee_df.m_br
    m_br_nmol = 1000 * m_br_umol
    model_space_to_whole_brain_volume = (1.0 - 0.125) / 2
    Lee_df.m_br = m_br_nmol * model_space_to_whole_brain_volume
    Lee_df = Lee_df[2:end, :]
    Lee_df[!, "time (h)"] .-= Lee_df[1, "time (h)"]
    Lee_df[!, "m_br"] .-= Lee_df[1, "m_br"]

    CSV.write(parsed_datafile, Lee_df)

    return Lee_df
end

function parse_PlaPEG_data(source_file::String, parsed_datafile::String)
    rawdf = DataFrame(CSV.File(source_file))

    PlaPEG_df = rawdf[rawdf.source .== "PlaPEG", ["time (h)", "m_br"]]
    m_br = PlaPEG_df.m_br
    PlaPEG_df.m_br = m_br

    CSV.write(parsed_datafile, PlaPEG_df)

    return PlaPEG_df
end

function parse_Cserr_data(
    raw_datafile::String,
    parsed_datafile::String,
    tracer::String="PEG900",
)
    rawdf = DataFrame(CSV.File(raw_datafile))

    Cserr_df = rawdf[rawdf.source .== "Cserr 1981", ["time (h)", "m_br", "m_c", "tracer"]]
    PEG900_df = Cserr_df[Cserr_df.tracer .== "PEG-900", :]
    PEG4000_df = Cserr_df[Cserr_df.tracer .== "PEG-4000", :]
    albumin_df = Cserr_df[Cserr_df.tracer .== "Albumin", :]

    if tracer == "PEG900"
        CSV.write(parsed_datafile, PEG900_df[1:end, 1:end-1])
        return PEG900_df
    elseif tracer == "PEG4000"
        CSV.write(parsed_datafile, PEG4000_df[1:end, 1:end-1])
        return PEG4000_df
    elseif tracer == "albumin"
        CSV.write(parsed_datafile, albumin_df[1:end, 1:end-1])
        return albumin_df
    else
        error("tracer-argument needs to be either PEG900, PEG4000 or albumin.")
    end
end

function parse_Pla_data(raw_datafile::String, parsed_datafile::String)
    rawdf = DataFrame(CSV.File(raw_datafile))
    Pla_df = rawdf[rawdf.source .== "Requena 2022", ["time (h)", "m_br", "c_b", "tracer"]]
    Pla_df = Pla_df[Pla_df.tracer .== "DB53 median", ["time (h)", "m_br", "c_b"]]

    # rawdf has the ratios of tracer mass in brain, translate here via known
    # infused tracer mass at 20 μg
    infused_tracer_mass = 20 # μg
    Pla_df[!, 2] = infused_tracer_mass * Pla_df[!, 2]

    # rawdf has the blood concentration in μg / mL, but since mm^3 = μL, we'll scale
    # this data to match with the units of the simulation
    mL_per_μL = 10^(-3)
    Pla_df[!, 3] = mL_per_μL * Pla_df[!, 3]

    CSV.write(parsed_datafile, Pla_df)

    return Pla_df
end

function parse_Stanton_data(
    source_file::String,
    parsed_datafile::String,
    α::Float64,
    downsample_rate::Int;
    treatment_group::String="kx1-6",
    pump_finish_time::Float64=0.2,
)
    rawdf = DataFrame(CSV.File(source_file))

    keep_cols_setup = ["time (h)", "c(1.5mm)", "c_v", "c_c", "setup"]
    keep_cols = ["time (h)", "c(1.5mm)", "c_v", "c_c"]

    Stanton_df = rawdf[rawdf[!, "time (h)"] .> pump_finish_time, keep_cols_setup]
    Stanton_df = Stanton_df[Stanton_df[!, "setup"] .== treatment_group, keep_cols]
    Stanton_df[!, "time (h)"] = Stanton_df[!, "time (h)"] .- pump_finish_time
    #Stanton_df[!, "c(1.5mm)"] = Stanton_df[!, "c(1.5mm)"] ./ α

    Stanton_df = Stanton_df[2:downsample_rate:end, :]

    CSV.write(parsed_datafile, Stanton_df)

    return Stanton_df
end # function parse_Stanton_data

function parseNatalieIntrastriatal(
    sourcefile::String,
    parsedfile::String;
    signalcolbasename="DCE_",
    framesperhour::Int=10,
    baselinesubtract::Bool=true,
    normalizedf::Bool=true,
    cropfrommax::Bool=true,
    normcolumn::String="c_c",
    croppump_finish_time::Float64=0.0,
    brainconcstring::String="c(0.25mm)",
    numparenchymaroi=7,
)
    df = parsegadobutrolinfusion(
        sourcefile,
        parsedfile;
        signalcolbasename,
        framesperhour,
        baselinesubtract,
        normalizedf,
        cropfrommax,
        normcolumn,
        croppump_finish_time,
        brainconcstring,
        numparenchymaroi,
    )

    df = allowmissing(df)
    df[2:end, :c_c] .= NaN

    CSV.write(parsedfile, df)

    return df
end

function parsegadobutrolinfusion(
    sourcefile::String,
    parsedfile::String;
    signalcolbasename="DCE_",
    framesperhour::Int=10,
    baselinesubtract::Bool=true,
    normalizedf::Bool=true,
    cropfrommax::Bool=true,
    normcolumn::String="c_c",
    croppump_finish_time::Float64=0.0,
    brainconcstring::String="c(0.25mm)",
    numparenchymaroi=7,
)

    rawdf = DataFrame(CSV.File(sourcefile))
    parseddf = DataFrame()

    VOInames = rawdf."Label Name"
    for VOIname in VOInames
        if VOIname == "Clear Label"
            continue
        end
        rawdataVOIdf = filter(:"Label Name" => n -> VOIname == n, rawdf)

        meanFITCdf = rawdataVOIdf[:, Cols(x -> all(contains.(x, ["mean", signalcolbasename])))]

        parseddf[!, VOIname] = Matrix(meanFITCdf)'[:, 1]
    end

    if baselinesubtract
        baseline = mean.(eachcol(parseddf[1:2, :]))
        parseddf .-= repeat(baseline', size(parseddf, 1))
        #parseddf[2, Not(:"time (h)")] = zeros(size(parseddf, 2) - 1)
        parseddf = parseddf[2:end, :]
    end

    lasttime = (size(parseddf, 1) - 1) / framesperhour
    parseddf[!, "time (h)"] = collect(0:1/framesperhour:lasttime)
    parseddf = parseddf[parseddf[!, "time (h)"] .> croppump_finish_time, :]
    parseddf[!, "time (h)"] = parseddf[!, "time (h)"] .- croppump_finish_time

    if size(parseddf[:, Cols(contains("SAS"))], 1) == size(parseddf, 1)
        parseddf[!, "c_c"] = mean.(eachrow(parseddf[:, Cols(contains("SAS"))]))
    end
    if size(parseddf[:, Cols(contains("ventricle"))], 1) == size(parseddf, 1)
        parseddf[!, "c_v"] = mean.(eachrow(parseddf[:, Cols(contains("ventricle"))]))
    end
    if size(parseddf[:, Cols(contains("parench"))], 1) == size(parseddf, 1)
        parseddf[!, brainconcstring] = mean.(eachrow(parseddf[:, Cols(contains("parench"))]))
    end

    if normalizedf
        normto = maximum(parseddf[:, normcolumn])
        parseddf[!, Not("time (h)")] = 100 .* parseddf[:, Not("time (h)")] ./ normto
    end

    if cropfrommax
        maxindex = argmax(sum(eachcol(parseddf[:, Cols(contains("parenchyma"))])))
        parseddf = parseddf[maxindex:end, :]
        parseddf[!, "time (h)"] .-= parseddf[1, "time (h)"]
    end

    parenchyma_cols_sorted = sort([n for n in names(parseddf) if contains(n, "parenchyma")])
    select!(
        parseddf,
        Symbol("time (h)"),
        :c_v,
        :c_c,
        #Not([Symbol("time (h)"), :c_v, :c_c]),
        Cols(parenchyma_cols_sorted)
    )
    #parseddf = parseddf[:, sort(names(parseddf))]

    parseddf = parseddf[:, 1:numparenchymaroi+3]

    CSV.write(parsedfile, parseddf)

    return parseddf
end

"""
    initial_condition(
        concentrations_df,
        outputlocations,
        concentration_measurement_locations,
    )

Creates an initial condition vector for the brian on locations outputlocations with the
concentrations interpolated from the concentrations_df and concentration_measurement_locations.
"""
function icfromdf(
    ic_measurements,
    concentration_measurement_locations,
    outputlocations,
)
    @assert(size(ic_measurements)==size(concentration_measurement_locations))

    itp = LinearInterpolation(
        concentration_measurement_locations,
        ic_measurements,
        extrapolation_bc=Line(),
    )
    output = itp(outputlocations)

    return output
end

export load_data!, parse_Lee_data, parse_Cserr_data, parse_Pla_data, parse_Stanton_data
export parsegadobutrolinfusion, icfromdf
export parseNatalieIntrastriatal
