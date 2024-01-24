
function interpretparametersymbol(s::AbstractString)
    s = replace(
        s,
        "alpha" => "α",
        "tau" => "τ",
        "mu" => "μ",
        "Delta" => "Δ",
    )
    return Symbol(s)
end

function loadparams(paramsfile; kwargs...)
    ks = keys(kwargs)

    csvfilecontent = CSV.File(paramsfile)
    paramsymbols = tuple(interpretparametersymbol.((fi.symbol for fi in csvfilecontent))...)
    paramvalues = tuple((fi.value for fi in csvfilecontent)...)

    baseparameters = merge(NamedTuple{paramsymbols}(paramvalues), kwargs)

    r_space, ∇, Δ = discretize_line(baseparameters.Δr, baseparameters.L)

    brain_conc_depth = baseparameters.L / 2
    brain_conc_index = argmin((r_space .- brain_conc_depth).^2)

    V_v = baseparameters.V_c / 2
    V_c = baseparameters.V_c / 2

    estimate_D = :estimate_D in ks ? kwargs[:estimate_D] : true
    estimate_Dm = :estimate_Dm in ks ? kwargs[:estimate_Dm] : true
    estimate_v = :estimate_v in ks ? kwargs[:estimate_v] : true
    estimate_all_phys = :estimate_all_phys in ks ? kwargs[:estimate_all_phys] : false

    computedparameters = (;
        V_v, V_c, r_space,
        ∇, Δ,
        brain_conc_depth, brain_conc_index,
        estimate_D, estimate_Dm, estimate_v, estimate_all_phys,
    )

    allparameters = merge(
        baseparameters,
        computedparameters,
    )

    return allparameters
end


function get_mouse_params(geometry; kwargs...)
    ks = keys(kwargs)
    D_DB53_tort1_7 = 1.28e-06 # cm^2 / s = 1.28 * 10^(-6) * (10 mm)^2 / s
                          #         = 1.28 * 10^(-4) mm^2 / s = 36 * 1.28 10^(-2) mm / h
    D = 36 * 1.28 * 10^(-2) # mm^2 / h
    v = 0.0 # mm / h

    #println("Using slab_ventricles geometry")
    R = :R in ks ? kwargs[:R] : 2.0 # mm
    S = :S in ks ? kwargs[:S] : 52 # R^2
    V_t = S * R

    # Pardridge 2016 references book by Davson from '87 that mouse V_c = 35 microlitre = 35 mm^3.
    V_c = 36 # mm^3 = microliter
    V_v = V_c / 2
    V_c = V_v

    α = :α in ks ? kwargs[:α] : 0.2
    V = α * V_t

    V_b = 1490.0 # mm^3, so 1500 = 1.5 mL

    k_e = :k_e in ks ? kwargs[:k_e] : 0.0 # per hour
    k_c = :k_c in ks ? clamp(kwargs[:k_c], 0.0, Inf) : 0.04 # per hour
    Lm = :Lm in ks ? kwargs[:Lm] : 0.01 # How deep is the ependymal layer?
    Dm = 0.005 # Lm * (k_e + k_c)

    k_p = 0.0 # per hour
    Q_l = 20 # mm^3 per hour, 100 nL / min, Guojun's estimate
    Q_r = 3.6 # mm^3 per hour, scaled from human data in Tuura 2021

    σ_brain_conc = 1.0 # from Stanton fig 4
    σ_brain_mass = 2.0 # nmol from Lee et al 2018 fig 6c
    σ_vent = 1.0 # from Stanton
    σ_csf_conc = 1.0 # 10 from Stanton fig 2
    σ_csf_mass = 0.1 # from Cserr fig 3
    σ_blood = 0.00020 # from Pla

    Δr = :Δr in ks ? kwargs[:Δr] : 0.01
    r_space, ∇, Δ = discretize_line(Δr, R)
    T_max = :T_max in ks ? kwargs[:T_max] : 4.5

    bc_x=zeros(Real, length(r_space))
    bc_xx=zeros(Real, length(r_space))

    brain_conc_depth = :brain_conc_depth in ks ? kwargs[:brain_conc_depth] : R / 2
    brain_conc_index = argmin((r_space .- brain_conc_depth).^2)

    estimate_D = :estimate_D in ks ? kwargs[:estimate_D] : true
    estimate_Dm = :estimate_Dm in ks ? kwargs[:estimate_Dm] : true
    estimate_v = :estimate_v in ks ? kwargs[:estimate_v] : true

    params = (D = D, μ_D = D, v = v, Dm=Dm, Lm=Lm,
    R = R, V_t = V_t, α = α,
    V = V, V_c = V_c, V_v = V_v, V_b = V_b,
    k_p = k_p, Q_l = Q_l, Q_r = Q_r, S = S,
    Δr = Δr, r_space = r_space, ∇ = ∇, Δ = Δ,  bc_x=bc_x, bc_xx=bc_xx,
    J = length(r_space), T_max = T_max, brain_conc_depth=brain_conc_depth,
    brain_conc_index,
    σ_brain_conc, σ_brain_mass, σ_vent, σ_csf_conc, σ_csf_mass, σ_blood,
    estimate_D=estimate_D, estimate_Dm=estimate_Dm, estimate_v=estimate_v,)
    update_namedtuple(params; kwargs...)
end

function get_rat_params(geometry; kwargs...)
    ks = keys(kwargs)
    D_DB53 = 36 * 1.28 * 10^(-2) # mm^2 / h
    D = D_DB53

    v = 0.00 # mm / h

    R = :R in ks ? kwargs[:R] : 3.0 # mm
    S = :S in ks ? kwargs[:S] : 145 # R^2
    V_t = R * S

    # Pardridge 2016 references book by Davson from '87 that mouse V_c = 35 microlitre = 35 mm^3.
    V_c = 370 # mm^3 = microliter
    V_v = V_c / 2
    V_c = V_v

    α = :α in ks ? kwargs[:α] : 0.2
    V = α * V_t

    V_b = 10.2 * 10^3 # mm^3, so 1500 = 1.5 mL

    k_e = :k_e in ks ? kwargs[:k_e] : 0.0 # per hour
    k_c = :k_c in ks ? kwargs[:k_c] : 0.05 # per hour
    Lm = :Lm in ks ? kwargs[:Lm] : 0.01 # How deep is the ependymal layer?
    Dm = 0.005 # Lm * (k_e + k_c)

    k_p = 0.0 # per hour
    Q_l = 200 # mm^3 per hour, 100 nL / min, Guojun's estimate
    Q_r = 36

    Δr = :Δr in ks ? kwargs[:Δr] : 0.01
    r_space, ∇, Δ = discretize_line(Δr, R)
    T_max = :T_max in ks ? kwargs[:T_max] : 4.5

    bc_x=zeros(Real, length(r_space))
    bc_xx=zeros(Real, length(r_space))

    brain_conc_depth = :brain_conc_depth in ks ? kwargs[:brain_conc_depth] : R / 2
    brain_conc_index = argmin((r_space .- brain_conc_depth).^2)

    mu_D = D
    tau_D_square = 0.5^2
    D_lower_trunc = 0.01
    D_upper_trunc = 20

    tau_v_square = 1.0
    v_lower_trunc = -10
    v_upper_trunc = 10

    Dm_overline = 1.0

    σ_brain_conc = 1.0 # from Stanton fig 4
    σ_brain_mass = 2.0 # nmol from Lee et al 2018 fig 6c
    σ_vent = 1.0 # from Stanton
    σ_csf_conc = 1.0 # 10 from Stanton fig 2
    σ_csf_mass = 0.1 # from Cserr fig 3
    σ_blood = 0.00020 # from Pla

    estimate_D = :estimate_D in ks ? kwargs[:estimate_D] : true
    estimate_Dm = :estimate_Dm in ks ? kwargs[:estimate_Dm] : true
    estimate_v = :estimate_v in ks ? kwargs[:estimate_v] : true

    params = (D = D, v = v, Dm=Dm, Lm=Lm,
    R = R, V_t = V_t, α = α,
    V = V, V_c = V_c, V_v = V_v, V_b = V_b, k_e = k_e,
    k_c = k_c,
    k_p = k_p, Q_l = Q_l, Q_r = Q_r, S = S,
    Δr = Δr, r_space = r_space, ∇ = ∇, Δ = Δ,
    J = length(r_space), T_max = T_max, brain_conc_depth=brain_conc_depth,
    σ_brain_conc, σ_brain_mass, σ_vent, σ_csf_conc, σ_csf_mass, σ_blood,
    estimate_D=estimate_D, estimate_Dm=estimate_Dm, estimate_v=estimate_v,
    geometry = geometry, bc_x=bc_x, bc_xx=bc_xx,
    tau_D_square = tau_D_square, D_lower_trunc = D_lower_trunc,
    D_upper_trunc = D_upper_trunc, μ_D = mu_D,
    tau_v_square = tau_v_square, v_lower_trunc = v_lower_trunc,
    v_upper_trunc, brain_conc_index,
    Dm_overline = Dm_overline,)

    update_namedtuple(params; kwargs...)

end

function update_namedtuple(nt::NamedTuple; kwargs...)
    #println("WARNING: You're responsible for computations required on e.g. r_p if setting r_e, V and S if setting R or α and r_space if setting R!")
    tmp_d = Dict(
        string(k) => v
        for (k, v) in zip(keys(nt), values(nt))
    )

    kwarg_d = Dict(kwargs)
    for (k, v) in zip(keys(kwarg_d), values(kwarg_d))
        tmp_d[string(k)] = v
    end

    return NamedTuple((Symbol(k), v) for (k, v) in zip(keys(tmp_d), values(tmp_d)))
end

export loadparams, update_namedtuple
