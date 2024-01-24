using Plots, Plots.PlotMeasures
using LaTeXStrings
pyplot()

savedir = "visualizations/"
# Exported from maple
cg = (Pe,phi__1,phi__3) -> (-16 * Pe * phi__1 * phi__3 + 16 * exp(Pe) * Pe * phi__1 * phi__3 + 4 * exp(Pe) * Pe ^ 2 * phi__3 ^ 2 + 16 * exp(Pe) * phi__1 * phi__3 ^ 2 + 4 * exp(Pe) * Pe ^ 3 * phi__3 - 4 * Pe ^ 3 * phi__3 + 4 * Pe ^ 2 * phi__3 ^ 2 + 16 * phi__1 * phi__3 ^ 2 - (Pe ^ 2 + 4 * phi__1) ^ (3//2) * phi__3 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) - (Pe ^ 2 + 4 * phi__1) ^ (3//2) * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 + Pe ^ 3 * exp(Pe) * sqrt(Pe ^ 2 + 4 * phi__1) - (Pe ^ 2 + 4 * phi__1) ^ (3//2) * exp(Pe) * Pe + (Pe ^ 2 + 4 * phi__1) ^ (3//2) * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 - 4 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 ^ 2 + 4 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 ^ 2 - 8 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 * phi__3 ^ 3 + 8 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 * phi__3 ^ 3 + 8 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 ^ 2 * phi__3 + 24 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 ^ 2 * phi__3 + exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * (Pe ^ 2 + 4 * phi__1) ^ (3//2) * phi__3 ^ 2 + 3 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * (Pe ^ 2 + 4 * phi__1) ^ (3//2) * phi__3 ^ 2 - exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 3 * sqrt(Pe ^ 2 + 4 * phi__1) - 8 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * phi__3 ^ 3 - 12 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 2 * phi__3 ^ 2 - 5 * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 2 * phi__3 - 4 * sqrt(Pe ^ 2 + 4 * phi__1) * Pe * phi__3 ^ 2 - 4 * sqrt(Pe ^ 2 + 4 * phi__1) * phi__1 * phi__3 - 4 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 3 * phi__3 - 32 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 * phi__3 ^ 2 + (Pe ^ 2 + 4 * phi__1) ^ (3//2) * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__3 - 8 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * phi__3 ^ 3 + 4 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 2 * phi__3 ^ 2 - exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 3 * sqrt(Pe ^ 2 + 4 * phi__1) + 4 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 3 * phi__3 + (Pe ^ 2 + 4 * phi__1) ^ (3//2) * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe - (Pe ^ 2 + 4 * phi__1) ^ (3//2) * exp(Pe) * phi__3 + 4 * sqrt(Pe ^ 2 + 4 * phi__1) * Pe * phi__1 + (Pe ^ 2 + 4 * phi__1) ^ (3//2) * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe + 8 * exp(Pe) * sqrt(Pe ^ 2 + 4 * phi__1) * phi__3 ^ 3 + (Pe ^ 2 + 4 * phi__1) ^ (3//2) * phi__3 + 8 * sqrt(Pe ^ 2 + 4 * phi__1) * phi__3 ^ 3 + Pe ^ 3 * sqrt(Pe ^ 2 + 4 * phi__1) - (Pe ^ 2 + 4 * phi__1) ^ (3//2) * Pe + exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 2 * phi__1 - exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 2 * phi__1 + 4 * exp(Pe) * sqrt(Pe ^ 2 + 4 * phi__1) * phi__1 * phi__3 + 4 * Pe * phi__1 * exp(Pe) * sqrt(Pe ^ 2 + 4 * phi__1) - 4 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 * sqrt(Pe ^ 2 + 4 * phi__1) + exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * (Pe ^ 2 + 4 * phi__1) ^ (3//2) * Pe * phi__3 - exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 3 * phi__3 - exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 2 * phi__3 ^ 2 + exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * (Pe ^ 2 + 4 * phi__1) ^ (3//2) * Pe * phi__3 - exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 3 * phi__3 - 3 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 2 * phi__3 ^ 2 + 4 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 2 * phi__1 * phi__3 + 4 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 2 * phi__1 * phi__3 - 4 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 * phi__3 ^ 2 + 4 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 * phi__3 ^ 2 + 3 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 2 * phi__3 - 4 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * sqrt(Pe ^ 2 + 4 * phi__1) * phi__3 ^ 2 - 3 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 2 * phi__3 - 4 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe * phi__3 ^ 2 - 4 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * sqrt(Pe ^ 2 + 4 * phi__1) * phi__1 * phi__3 - 8 * exp(Pe / 2 + sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 * phi__3 + 8 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 * phi__3 - 4 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 * sqrt(Pe ^ 2 + 4 * phi__1) + 5 * exp(Pe) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe ^ 2 * phi__3 + 4 * exp(Pe / 2 - sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 * phi__3 * sqrt(Pe ^ 2 + 4 * phi__1) + 12 * exp(Pe) * sqrt(Pe ^ 2 + 4 * phi__1) * Pe * phi__3 ^ 2) * exp(-Pe / 2) / phi__1 / (sqrt(Pe ^ 2 + 4 * phi__1) * exp(-sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__3 + sqrt(Pe ^ 2 + 4 * phi__1) * exp(-sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__3 ^ 2 + sqrt(Pe ^ 2 + 4 * phi__1) * exp(sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__3 + 3 * sqrt(Pe ^ 2 + 4 * phi__1) * exp(sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__3 ^ 2 + exp(-sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 2 * phi__3 - exp(-sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__3 ^ 2 - 2 * exp(-sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__3 ^ 3 + exp(sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe ^ 2 * phi__3 + exp(sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__3 ^ 2 + 2 * exp(sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__3 ^ 3 - sqrt(Pe ^ 2 + 4 * phi__1) * exp(-sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 + sqrt(Pe ^ 2 + 4 * phi__1) * exp(sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 - exp(-sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 + 2 * exp(-sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 * phi__3 + exp(sqrt(Pe ^ 2 + 4 * phi__1) / 2) * Pe * phi__1 + 6 * exp(sqrt(Pe ^ 2 + 4 * phi__1) / 2) * phi__1 * phi__3) / 4

cg(4.8, 0.1, 1)
cg(9.4, 0.1, 1)
cg(9.4, 0.1, 1) / cg(4.8, 0.1, 1)
cg(25, 0.1, 1) / cg(5, 0.1, 1)

Peline = 0.01:0.01:5.0
phi__3_line = 0.01:0.01:1.0

Peline = exp.(collect(log(0.1):0.1:log(10.0)))
phi__3_line = exp.(collect(log(0.01):0.1:log(1.0)))

estimatedparams = [
    (v = 0.16, D = 0.55, D_m = 5e-3),
    (v = 0.09, D=0.45, D_m = 4e-3),
    (v = 0.19, D = 0.58, D_m = 5.5 * 1e-3),
    (v = 0.18, D = 0.57, D_m = 4.1 * 1e-3),
    (v = 0.08, D = 0.38, D_m = 4.5e-3),
    (v = 0.15, D = 0.57, D_m = 4.9e-3),
    (v = 0.09, D = 0.28, D_m = 3.4e-3),
]
mouse02_params = (v=0.16, D=0.55, D_m = 5e-3)

L = 2.0
L_m = 0.01
α = 0.2
calcγ(p) = α * (p[:D_m] / L_m) / (p[:D] / L)
calcPe(p) = p[:v] * L / p[:D]
estimated_γs = calcγ.(estimatedparams)
estimated_Pes = calcPe.(estimatedparams)
mouse02_γ = calcγ(mouse02_params)
mouse02_Pe = calcPe(mouse02_params)

#gr()
pyplot()
p = contour(
    Peline,
    phi__3_line,
    (x, y) -> cg(x, 0.1, y),
    fill=true,
    #c=cgrad(:lajolla, rev=true),
    #c=cgrad(:algae, rev=false),
    c=cgrad(:tempo, rev=false),
    #linecolor=:black,
    #linewidth=0.5,
    xscale=:log10,
    yscale=:log10,
    xlabel="Péclet",
    ylabel="relative surface transport",#"γ",
    thickness_scaling=2,
    #bottom_margin=-20px,
    #left_margin=-30px,
    right_margin=-30px,
    #labelfontsize=8,
    #tickfontsize=6,
    #colorbar_title_location=:top,
    colorbar_title="endogenous waste",
    colorbar_ticks=[0.2, 0.5, 0.8],
    foreground_color_border=:white,
    #yguidefontrotation=-90,
    labelfontsize=9,
    tickfontsize=7,
    colorbar_tickfontsize=7,
)
xticks!([0.1, 1.0, 10.0], ["0.1", "1.0", "10.0"])
yticks!([0.012, 0.10, 1.0], ["0.01", "0.1", "1.0"])

for (Pe, γ) in zip(estimated_Pes, estimated_γs)
    annotate!(Pe, γ, text("+", :gray, :center, 10))
end

annotate!(mouse02_Pe, mouse02_γ, text("+", :black, :center, 10))
#scatter!([Pe, ], [γ, ], label="", color=:black)

savefig(p, savedir * "nondimensionalwastemass.pdf")
savefig(p, savedir * "nondimensionalwastemass.svg")
