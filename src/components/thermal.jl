"""
   get_thermal(p=nothing)

Plugs the different models for the thermal machines together. The thermal nodes consist of:
    - Fifth Order Machine Model + Shaft
    - Exciter
    - Power System Stabilizer
    - PQ Load
"""
function get_thermal(p=nothing)
    shaft = get_shaft_thermal()
    machine_fifth_t = get_machine_fifth_thermal() # we ended up using the fifth order model for the thermal nodes as well
    exciter = get_exciter()
    rfc = get_rfc()
    PQ = get_internal_current()

    thermal_no_rfc = IOSystem([
        machine_fifth_t.P_e => shaft.P_e,
        machine_fifth_t.v_d => exciter.v_d,
        machine_fifth_t.v_q => exciter.v_q,
        exciter.E_f       => machine_fifth_t.E_f,
        machine_fifth_t.P_e => exciter.P_e], # pss is part of the exciter 
        [machine_fifth_t, shaft, exciter],
        namespace_map = [
            # Machine Parameter names
            machine_fifth_t.X_d   => :X_d,
            machine_fifth_t.X_q′  => :X_q′,
            machine_fifth_t.X_q′′ => :X_q′′,
            machine_fifth_t.X_d′  => :X_d′,
            machine_fifth_t.X_q′′ => :X_q′′,
            machine_fifth_t.T_q′′ => :T_q′′,
            machine_fifth_t.T_d′  => :T_d′,
            machine_fifth_t.T_d′′ => :T_d′′,
            exciter.E_f         => :E_f,
            # Shaft parameters
            shaft.D     => :D, # damping
            shaft.H     => :H, # Inertia Constant
            shaft.Ω_H   => :Ω_H,
            shaft.P_m => :P_m_shaft,
            exciter.P_m => :P_m_exciter],
        outputs = [machine_fifth_t.v_d, machine_fifth_t.v_q, shaft.ω, shaft.δ], 
        name = :thermal_no_rfc) 

    thermal = IOSystem([
        rfc.i_d => thermal_no_rfc.i_d,
        rfc.i_q => thermal_no_rfc.i_q,
        thermal_no_rfc.v_d => rfc.v_d,
        thermal_no_rfc.v_q => rfc.v_q,
        thermal_no_rfc.δ   => rfc.δ, 
        rfc.u_r          => PQ.u_r,
        rfc.u_i          => PQ.u_i,
        PQ.i_int_r       => rfc.i_r,
        PQ.i_int_i       => rfc.i_i],
        [rfc, thermal_no_rfc, PQ],
        outputs = [rfc.u_r, rfc.u_i, thermal_no_rfc.ω, thermal_no_rfc.δ], 
        name = :thermal)
    thermal =  connect_system(thermal)
    isnothing(p) ? thermal :  set_p(thermal, p)
end


"""
    get_machine_fifth_thermal()

Fifth order machine model for the thermal nodes. The original publication uses sixth order machines.
Voltage equations of the salient rotor using the fifth order model from J. Machowski 'Power System Dynamics' (section 11.103 on page 455).
"""
function get_machine_fifth_thermal()
    @parameters t
    dt = Differential(t)

    @parameters T_q′′ T_d′ T_d′′ X_d X_d′ X_q′ X_d′′ X_q′′ E_f(t) i_q(t) i_d(t) #X_q T_q′
    @variables v_q(t) v_d′(t) v_q′(t) P_e(t) v_d(t)

    machine_fifth = IOBlock([
        dt(v_q)  ~ (1/T_d′) * (E_f - v_q + i_d * (X_d - X_d′)),
        v_d      ~ 0,
        dt(v_q′) ~ (1/T_d′′) * (v_q - v_q′ + (X_d′ - X_d′′) * i_d),
        dt(v_d′) ~ (1/T_q′′) * (v_d - v_d′ - (X_q′ - X_q′′) * i_q),
        P_e      ~ (v_d′ * i_d + v_q′ * i_q) + (X_d′′ - X_q′′) * i_d * i_q],
        [i_d, i_q, E_f],
        [v_d, v_q, P_e],
        name = :machine_fifth)
    return machine_fifth
end


"""
    get_exciter()
Exciter connected to the thermal machines.
"""
function get_exciter()
    @parameters t
    dt = Differential(t)
    @parameters v_ref v_filt(t) v_stab(t) T1 T2 i_1(t) i_2(t) v_d(t) v_q(t) dti(t) i(t)
    @variables Δv(t) o(t) Δ(t) V(t) E_f(t) dto(t)
    pss = get_pss()

    vol_measurement = IOBlock([V ~ (v_d^2 + v_q^2)^(1/2)], [v_d, v_q], [V], name = :vol_measurement)
    vol_filter = IOBlock([dto ~ (1 / T1) * (i - o), dt(o) ~ dto], [i], [o, dto], name = :vol_filter)
    vol_diff = IOBlock([Δv ~ v_ref + v_filt - v_stab], [v_filt, v_stab], [Δv], name = :vol_diff)
    damping = IOBlock([dt(o) ~ 1/T2 * (T1 * dti - o)], [dti], [o], name = :damping) # H(s) = T1 * s / (T2 * s + 1) 
    main_reg = IOBlock([dt(o) ~ 1/T2 * (i * T1 - o)], [i], [o], name = :main_reg) # H(s) = T1/(T2 * s + 1)
    delta = IOBlock([Δ ~ i_1 - i_2], [i_1, i_2], [Δ], name = :delta)

    @parameters i_1(t) i_2(t)
    @variables o(t)
    max_block = IOBlock([o ~ max(i_1, i_2)], [i_1, i_2], [o])
    saturation = IOBlock([E_f ~ 5 / (1 + exp(2.5 - v_filt))], [v_filt], [E_f], name = :saturation)

    # Connecting the blocks into the excitation system
    exciter = IOSystem([
        vol_measurement.V => vol_filter.i,
        vol_filter.o      => vol_diff.v_filt,
        pss.v_stab => vol_diff.v_stab,
        vol_filter.dto    => damping.dti,
        vol_filter.o => max_block.i_1,
        main_reg.o => max_block.i_2,
        max_block.o => saturation.v_filt,  
        vol_diff.Δv => delta.i_1,
        damping.o         => delta.i_2,
        delta.Δ           => main_reg.i],
        [vol_measurement, main_reg, delta, vol_filter, vol_diff, max_block, damping, saturation, pss],
        namespace_map = [
            vol_filter.o => :E_f,# Electrical Field Voltage 
            main_reg.T1  => :main_reg_T1,
            main_reg.T2  => :main_reg_T2,
            damping.T1   => :damp_T1,
            damping.T2   => :damp_T2,
            vol_measurement.v_d => :v_d,
            vol_measurement.v_q => :v_q,
            vol_filter.T1 => :vol_filter_T1], 
            outputs=[:E_f], 
            name = :exciter)
    return connect_system(exciter)
end

"""
   get_pss()

Power System Stabilizer connected to the thermal nodes.
"""
function get_pss()
    @parameters t
    dt = Differential(t)
    # Power System Stabilizer (PSS) Blocks
    @parameters T1 T2 K i(t) dti(t) P_e(t) P_m
    @variables o(t) dto(t) ΔP(t) v_stab(t)

    # Creating the different blocks for the pss
    delta_P = IOBlock([ΔP ~ P_m - P_e], [P_e], [ΔP], name = :delta_P)
    filter_1 = IOBlock([dto ~ (1 / T1) * (i - o), dt(o) ~ dto], [i], [o, dto], name = :filter_1)
    washout = IOBlock([dto ~ dti - 2/9 * o, dt(o) ~ dto], [dti], [o, dto], name = :washout) # H(s) = 4.5s / (4.5s + 1)
    lead_lag_1 = IOBlock([dto ~ (1 / T2) * (dti * T1 + i - o), dt(o) ~ dto], [i, dti], [o, dto], name = :lead_lag_1)
    lead_lag_2 = IOBlock([dto ~ (1 / T2) * (dti * T1 + i - o), dt(o) ~ dto], [i, dti], [o, dto], name = :lead_lag_2)
    gain = IOBlock([o ~ K * i], [i], [o], name = :gain)

    # A sigmoid function is used to approximate the saturation to avoid jumps in the ODEFunction
    saturation = IOBlock([v_stab ~ (0.1 / (1 + exp(-50*i))) - 0.05], [i], [v_stab], name = :saturation) 

    # Connecting the blocks into the pss system
    pss = IOSystem([
        delta_P.ΔP     => filter_1.i,
        filter_1.dto   => washout.dti,
        washout.o      => lead_lag_1.i,
        washout.dto    => lead_lag_1.dti,
        lead_lag_1.o   => lead_lag_2.i,
        lead_lag_1.dto => lead_lag_2.dti,
        lead_lag_2.o   => gain.i,
        gain.o         => saturation.i],
        [delta_P, filter_1, washout, lead_lag_1, lead_lag_2, gain, saturation],
        namespace_map = [
            saturation.v_stab => :v_stab,
            delta_P.P_m       => :P_m,
            delta_P.P_e       => :P_e,
            # parameter names as in the paper
            filter_1.T1    => :filter_T1,
            lead_lag_1.T1  => :T1, 
            lead_lag_1.T2  => :T2, 
            lead_lag_2.T1  => :T3, 
            lead_lag_2.T2  => :T4, 
            gain.K         => :K_pss], 
            outputs = [:v_stab], 
            name = :pss)
    return connect_system(pss)
end