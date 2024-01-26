
"""
    enclosethetimedifferential(fixedparameters)

Provides the closed derivative function, which knows all the
fixed parameters but takes the D, Dm and v parameters in its p argument.
"""
function enclosethetimedifferential(fixedparameters::NamedTuple)::Function

    @unpack (r_space, Δr, Δ, ∇) = fixedparameters
    bc_x = zeros(Real, length(r_space))
    bc_xx = zeros(Real, length(r_space))

    function timedifferentialclosure!(du, u, p, t)

        @unpack (α, D, v, k_p, V_c, Q_l, Q_r, V_b,
        S, L_m, D_m, V_v, ) = p
        c = u[1:end - 3]
        c_v = u[end - 2]
        c_c = u[end - 1]
        c_b = u[end]

        J_B0 = (D_m/L_m) * (α * c_v - c[1])
        J_BL = (D_m/L_m) * (c[end] - α * c_c)

        grad_0 = (v/D) * c[1] - J_B0 / D
        grad_L = (v/D) * c[end] - J_BL / D

        bc_x[1] = grad_0 / 2
        bc_x[end] = grad_L / 2
        grad_c = ∇ * c + bc_x

        bc_xx[1] = -grad_0 / Δr
        bc_xx[end] = grad_L / Δr
        Lap_c = Δ * c + bc_xx

        # Important here to use the ECS concentration c (not the average tisssue concentration).
        #C = sum(Δr .* S * α * (k_p * (c .- c_b)))
        C = sum(Δr .* S * (k_p * (c .- c_b)))

        dc_dt = D * Lap_c - v * grad_c .- k_p * (c .- c_b)
        du[1:end-3] = dc_dt[1:end]

        #dcv_dt = -S * α * J_B0 / V_v - (Q_l / V_v) * c_v + (Q_r / V_v) * c_c
        dcv_dt = -S * J_B0 / V_v - (Q_l / V_v) * c_v + (Q_r / V_v) * c_c
        du[end-2] = dcv_dt

        #dcc_dt = S * α * J_BL / V_c + (Q_l / V_c) * c_v - (Q_l/V_c) * c_c - (Q_r / V_v) * c_c
        dcc_dt = S * J_BL / V_c + (Q_l / V_c) * c_v - (Q_l/V_c) * c_c - (Q_r / V_v) * c_c
        du[end-1] = dcc_dt

        dcb_dt = (Q_l/V_b) * c_c + C / V_b
        du[end] = dcb_dt
    end # closure

    return timedifferentialclosure!

end # enclosethetimedifferential

function slab_int(c, p)
    # Using α factor because this integral is used for total brain mass, where
    # the relevant volume to sum over is the ECS volume.
    sum(c .* p.S .* p.Δr)# * p.α)
end

function brain_mass_over_time(c, p)
    r_space, Δr = p.r_space, p.Δr
    J = length(r_space)
    mt = [slab_int(c[1:J, ti], p) for ti ∈ 1:size(c)[2]]
    return mt
end

function state_mass_over_time(states, p)
    if size(states)[1] - length(p.r_space) == 2
        bm_t = brain_mass_over_time(states[1:end-2, :], p)
        m_csf = p.V_c * states[end - 1, :]
        m_b = p.V_b * states[end, :]
        return bm_t, m_csf, m_b
    elseif size(states)[1] - length(p.r_space) == 3
        bm_t = brain_mass_over_time(states[1:end-3, :], p)
        m_v = p.V_v * states[end - 2, :]
        m_csf = p.V_c * states[end - 1, :]
        m_b = p.V_b * states[end, :]
        return bm_t, m_v, m_csf, m_b
    else
        error("Could not determine state space shape.")
    end
end

brain_int(sol, p) = Array(sol)[1:end-3, :]' * p.Δr * p.S * ones(length(p.r_space)) #* p.α

@model function posteriorsamplergeneratorCserr(timepoints, data, prob, parameters)

    σ_brain_mass ~ Uniform(0.0, parameters.overline_o_bm_Cserr)
    σ_csf_mass ~ Uniform(0.0, parameters.overline_o_c_Cserr)

    if parameters.estimate_D
        D ~ truncated(LogNormal(log(parameters.μ_D), parameters.τ_D^2), parameters.D_lower, parameters.D_upper)
    else
        D = parameters.μ_D
    end
    if parameters.estimate_Dm
        D_m ~ Uniform(0.0, parameters.overline_D_m)
    else
        D_m = parameters.D_m
    end
    if parameters.estimate_v
        v ~ truncated(Normal(0.0, parameters.τ_v), -parameters.v_truncation, parameters.v_truncation)
    else
        v = parameters.v
    end
    if parameters.estimate_all_phys
        α ~ truncated(Normal(0.2, 0.025), 0.1, 0.3)
        k_p ~ truncated(Normal(parameters.k_p, parameters.k_p/10), 0.0, 1000.0)
        Q_l ~ truncated(Normal(parameters.Q_l, parameters.Q_l/10), 0.0, 10 * parameters.Q_l)
        Q_r ~ truncated(Normal(parameters.Q_r, parameters.Q_r/10), 0.0, 10 * parameters.Q_r)
        V_c ~ truncated(Normal(parameters.V_c, parameters.V_c/10), 0.0, 10 * parameters.V_c)
        V_b ~ truncated(Normal(parameters.V_b, parameters.V_b/10), 0.0, 10 * parameters.V_b)
        V_v ~ truncated(Normal(parameters.V_v, parameters.V_v/10), 0.0, 10 * parameters.V_v)
        S ~ truncated(Normal(parameters.S, parameters.S/10), 0.0, 10*parameters.S)
        L_m ~ truncated(Normal(parameters.L_m, parameters.L_m/10), 0.0, 10*parameters.L_m)
    else
        α = parameters.α
        k_p = parameters.k_p
        V_c = parameters.V_c
        Q_l = parameters.Q_l
        Q_r = parameters.Q_r
        V_b = parameters.V_b
        S = parameters.S
        L_m = parameters.L_m
        V_v = parameters.V_v
    end

    p = ComponentArray(
        D=D,
        D_m=D_m,
        v=v,
        α=α,
        k_p=k_p,
        V_c=V_c,
        Q_l=Q_l,
        Q_r=Q_r,
        V_b=V_b,
        S=S,
        L_m=L_m,
        V_v=V_v,
    )

    re_prob = remake(prob, p=p)
    sol = solve(re_prob, alg_hints=[:stiff], saveat=timepoints, reltol=1e-2)
    if sol.retcode != :Success
        println(sol.retcode)
        println("(D, D_m, v) = ", string(D), ", ", string(D_m), ", ", string(v), ")")
        Turing.@addlogprob! -Inf
        return
    end

    brain_pred = brain_int(sol, parameters)
    m_csf_pred = parameters.V_c * sol[end - 1, :] # + params.V_v * sol[end - 2, :]

    for i = 1:size(brain_pred)[1]
        data[i, 1] ~ Normal(brain_pred[i], σ_brain_mass)
        data[i, 2] ~ Normal(m_csf_pred[i], σ_csf_mass)
    end

    #= Alternative 2 with MvNormal
    data[:, 2] ~ MvNormal(brain_pred, σ_brain)
    data[:, 3] ~ MvNormal(m_csf_pred, σ_csf_mass)
    =#
end

@model function posteriorsamplergeneratorPla(timepoints, data, prob, parameters)

    σ_brain_mass ~ Uniform(0.0, parameters.overline_o_bm_Pla)
    σ_blood ~ Uniform(0.0, parameters.overline_o_bl_Pla)

    if parameters.estimate_D
        D ~ truncated(LogNormal(log(parameters.μ_D), parameters.τ_D^2), parameters.D_lower, parameters.D_upper)
    else
        D = parameters.μ_D
    end
    if parameters.estimate_Dm
        D_m ~ Uniform(0.0, parameters.overline_D_m)
    else
        D_m = parameters.D_m
    end
    if parameters.estimate_v
        v ~ truncated(Normal(0.0, parameters.τ_v), -parameters.v_truncation, parameters.v_truncation)
    else
        v = parameters.v
    end
    if parameters.estimate_all_phys
        Q_l ~ truncated(Normal(parameters.Q_l, parameters.Q_l / 5), parameters.Q_l / 10, 10 * parameters.Q_l)
        Q_r ~ truncated(Normal(parameters.Q_r, parameters.Q_r / 5), parameters.Q_r / 10, 10 * parameters.Q_r)
        V_c ~ truncated(Normal(parameters.V_c, parameters.V_c / 5), parameters.V_c / 10, 10 * parameters.V_c)
        V_b ~ truncated(Normal(parameters.V_b, parameters.V_b / 5), parameters.V_b / 10, 10 * parameters.V_b)
        V_v ~ truncated(Normal(parameters.V_v, parameters.V_v / 5), parameters.V_v / 10, 10 * parameters.V_v)
        S ~ truncated(Normal(parameters.S, parameters.S / 5), parameters.S / 10, 10 * parameters.S)
        L_m ~ truncated(Normal(parameters.L_m, parameters.L_m / 5), parameters.L_m / 5, 5 * parameters.L_m)
        α ~ truncated(Normal(0.2, 0.01), 0.15, 0.25)
    else
        α = parameters.α
        V_c = parameters.V_c
        Q_l = parameters.Q_l
        Q_r = parameters.Q_r
        V_b = parameters.V_b
        S = parameters.S
        L_m = parameters.L_m
        V_v = parameters.V_v
    end

    p = ComponentArray(
        D=D,
        D_m=D_m,
        v=v,
        α=α,
        k_p=parameters.k_p,
        V_c=V_c,
        Q_l=Q_l,
        Q_r=Q_r,
        V_b=V_b,
        S=S,
        L_m=L_m,
        V_v=V_v,
    )

    re_prob = remake(prob, p=p)
    sol = solve(re_prob, alg_hints=[:stiff], saveat=timepoints, reltol=1e-2)
    if sol.retcode != :Success
        println(sol.retcode)
        println("(D, D_m, v) = ", string(D), ", ", string(D_m), ", ", string(v), ")")
        Turing.@addlogprob! -Inf
        return
    end

    brain_pred = brain_int(sol, parameters)
    c_b_pred = sol[end, :]

    # Measurements were made at the start and end of the experiment, ...
    data[1, 1] ~ Normal(brain_pred[1], σ_brain_mass)
    data[end, 1] ~ Normal(brain_pred[end], σ_brain_mass)
    # ... so to make predictions on brain mass in the intervening period as part of this
    # MCMC run, uncomment the section below.
    #=for i = 1:size(brain_pred)[1]
        data[i, 1] ~ Normal(brain_pred[i], σ_brain_mass)
    end=#

    for i = 1:size(brain_pred)[1]
        data[i, 2] ~ Normal(c_b_pred[i], σ_blood)
    end

end

@model function posteriorsamplergeneratorGadobutrolinfusion(timepoints, data, prob, parameters)

    σ_brain_conc ~ Uniform(0.0, parameters.overline_o_b_Bork)
    σ_vent ~ Uniform(0.0, parameters.overline_o_v_Bork)
    σ_csf_conc ~ Uniform(0.0, parameters.overline_o_c_Bork)

    if parameters.estimate_D
        D ~ truncated(LogNormal(log(parameters.μ_D), parameters.τ_D^2), parameters.D_lower, parameters.D_upper)
    else
        D = parameters.μ_D
    end
    if parameters.estimate_Dm
        D_m ~ Uniform(0.0, parameters.overline_D_m)
    else
        D_m = parameters.D_m
    end
    if parameters.estimate_v
        v ~ truncated(Normal(0.0, parameters.τ_v), -parameters.v_truncation, parameters.v_truncation)
    else
        v = parameters.v
    end
    if parameters.estimate_all_phys
        Q_l ~ truncated(Normal(parameters.Q_l, parameters.Q_l / 5), parameters.Q_l / 10, 10 * parameters.Q_l)
        Q_r ~ truncated(Normal(parameters.Q_r, parameters.Q_r / 5), parameters.Q_r / 10, 10 * parameters.Q_r)
        V_c ~ truncated(Normal(parameters.V_c, parameters.V_c / 5), parameters.V_c / 10, 10 * parameters.V_c)
        V_b ~ truncated(Normal(parameters.V_b, parameters.V_b / 5), parameters.V_b / 10, 10 * parameters.V_b)
        V_v ~ truncated(Normal(parameters.V_v, parameters.V_v / 5), parameters.V_v / 10, 10 * parameters.V_v)
        S ~ truncated(Normal(parameters.S, parameters.S / 5), parameters.S / 10, 10 * parameters.S)
        L_m ~ truncated(Normal(parameters.L_m, parameters.L_m / 5), parameters.L_m / 5, 5 * parameters.L_m)
        α ~ truncated(Normal(0.2, 0.01), 0.15, 0.25)
    else
        α = parameters.α
        V_c = parameters.V_c
        Q_l = parameters.Q_l
        Q_r = parameters.Q_r
        V_b = parameters.V_b
        S = parameters.S
        L_m = parameters.L_m
        V_v = parameters.V_v
    end

    p = ComponentArray(
        D=D,
        D_m=D_m,
        v=v,
        α=α,
        k_p=parameters.k_p,
        V_c=V_c,
        Q_l=Q_l,
        Q_r=Q_r,
        V_b=V_b,
        S=S,
        L_m=L_m,
        V_v=V_v,
    )

    re_prob = remake(prob, p=p)
    sol = solve(re_prob, alg_hints=[:stiff], saveat=timepoints, reltol=1e-2)
    if sol.retcode != :Success
        println(sol.retcode)
        println("(D, D_m, v) = ", string(D), ", ", string(D_m), ", ", string(v), ")")
        Turing.@addlogprob! -Inf
        return
    end

    c_parenchyma1 = sol[parameters.brain_conc_index1, :]
    c_parenchyma2 = sol[parameters.brain_conc_index2, :]
    c_parenchyma3 = sol[parameters.brain_conc_index3, :]
    c_parenchyma4 = sol[parameters.brain_conc_index4, :]
    c_parenchyma5 = sol[parameters.brain_conc_index5, :]
    c_parenchyma6 = sol[parameters.brain_conc_index6, :]
    c_parenchyma7 = sol[parameters.brain_conc_index7, :]
    c_v_pred = sol[end - 2, :]
    c_c_pred = sol[end - 1, :]

    for i = 1:size(c_parenchyma1)[1]
        data[i, 1] ~ Normal(c_v_pred[i], σ_vent)
        data[i, 2] ~ Normal(c_c_pred[i], σ_csf_conc)
        data[i, 3] ~ Normal(c_parenchyma1[i], σ_brain_conc)
        data[i, 4] ~ Normal(c_parenchyma2[i], σ_brain_conc)
        data[i, 5] ~ Normal(c_parenchyma3[i], σ_brain_conc)
        data[i, 6] ~ Normal(c_parenchyma4[i], σ_brain_conc)
        data[i, 7] ~ Normal(c_parenchyma5[i], σ_brain_conc)
        data[i, 8] ~ Normal(c_parenchyma6[i], σ_brain_conc)
        data[i, 9] ~ Normal(c_parenchyma7[i], σ_brain_conc)
    end
end

lppd(lpp) = sum(log(mean(exp.(lpp[k]))) for k in keys(lpp))
pWAIC(lpp) = sum(var(lpp[k]) for k in keys(lpp))
WAIC(lpp) = -2 * (lppd(lpp) - pWAIC(lpp))

lppd_obs(lpp) = [log(mean(exp.(lpp[k]))) for k in keys(lpp)]
pWAIC_obs(lpp) = (var(lpp[k]) for k in keys(lpp))
WAIC_obs(lpp) = -2 * (lppd_obs(lpp) .- pWAIC_obs(lpp))

function WAIC(model::Turing.Model, chains::Chains)

    num_cases = size(chains, 2)

    lpp = pointwise_loglikelihoods(model, get_sections(chains, :parameters))

    WAIC_each_observation = WAIC_obs(lpp)

    WAIC_std = sqrt(num_cases * var(WAIC_each_observation))
    WAIC = sum(WAIC_each_observation)

    return WAIC, WAIC_std
end

function calc_and_save_WAIC!(simulation::GlymphaticSimulation)
    simulation.WAIC = WAIC(simulation.turing_model, simulation.chains)

    return simulation.WAIC
end

export enclosethetimedifferential, state_mass_over_time, calc_and_save_WAIC!
export posteriorsamplergeneratorCserr
export posteriorsamplergeneratorPla
export posteriorsamplergeneratorGadobutrolinfusion
