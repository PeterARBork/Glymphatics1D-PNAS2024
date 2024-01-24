module Glymphatics1D

using DataFrames, CSV, Parameters, ComponentArrays
using DifferentialEquations, DiffEqOperators
using ForwardDiff, Turing, FillArrays
import MCMCChains: get_sections
import LinearAlgebra: Diagonal
import Statistics: mean, std
using Interpolations
using StatsPlots, Plots, Plots.PlotMeasures, LaTeXStrings
using Logging, Serialization


Turing.setadbackend(:forwarddiff)

include("structures.jl")
include("glymphaticmodel.jl")
include("parameters.jl")
include("parsing.jl")
include("runsimulations.jl")
include("glymphaticplots.jl")
include("CserrAnalyses.jl")
include("GadobutrolinfusionAnalyses.jl")
include("PlaAnalyses.jl")
include("distributionplots.jl")

end # module
