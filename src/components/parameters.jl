function get_hydro_machine_p()
    machine_p = [
        :X_d   => 1.1,
        :X_q′  => 0.7,
        :X_d′  => 0.25,
        :X_d′′ => 0.2,
        :X_q′′ => 0.2,
        :T_d′  => 5.0,
        :T_d′′ => 0.05,
        :T_q′′ => 0.1,
        :E_f   => 1.0]
    return machine_p
end

function get_exciter_p()
    exciter_p = [
    :main_reg_T1 => 300,   
    :main_reg_T2 => 0.001,
    :damp_T1 => 0.001, 
    :damp_T2 => 0.1,   
    :v_ref => 1.0,
    :vol_filter_T1 => 0.02]

    return exciter_p
end

function get_thermal_machine_p()
    thermal_p = [
        :T_d′  => 7.0,
        :T_q′′ => 0.05,
        :T_d′′ => 0.05,
        :X_d   => 2.2,
        :X_d′  => 0.3,
        :X_q′  => 0.4,
        :X_d′′ => 0.2,
        :X_q′′ => 0.2,
        :P_m_exciter   => 0.9,
        :P_m_shaft => 0.9]

    return thermal_p
end

function get_thermal_4_p()
    machine_p = get_thermal_machine_p()
    pss_p = [:filter_T1 => 0.01, :T1 => 0.1323, :T2 => 0.6743, :T3 => 0.1323, :T4 => 0.6743, :K_pss => 0.15]
    exciter_p = get_exciter_p()
    shaft_p = [:H => 6.0, :Ω_H => 2π * 50.0]

    return vcat(machine_p, shaft_p, pss_p, exciter_p)
end

function get_thermal_5_p()
    machine_p = get_thermal_machine_p()
    pss_p = [:filter_T1 => 0.01, :T1 => 0.1323, :T2 => 0.6743, :T3 => 0.1323, :T4 => 0.6743, :K_pss => 0.15]
    exciter_p = get_exciter_p()
    shaft_p = [:H => 6.0, :Ω_H => 2π * 50.0]
    
    return vcat(machine_p, shaft_p, pss_p, exciter_p, )
end

function get_hydro_1_p()
    machine_p = get_hydro_machine_p()
    hydro_gov_p = [
        :T_y => 0.2,
        :T_w => 0.7] # waterways time constant

    shaft_p = [:H => 3.0, :Ω_H => 2π * 50.0]

    vcat(machine_p, hydro_gov_p, shaft_p)
end

function get_hydro_2_p()
    machine_p = get_hydro_machine_p()
    
    hydro_gov_p = [
        :T_y => 0.2,
        :T_w => 1.4] # waterways time constant

    shaft_p = [:H => 3.0, :Ω_H => 2π * 50.0]

    vcat(machine_p, hydro_gov_p, shaft_p)
end

function get_hydro_3_p()
    machine_p = get_hydro_machine_p()
    hydro_gov_p = [
        :T_y => 0.2,
        :T_w => 1.4]

    shaft_p = [:H => 3.0, :Ω_H => 2π * 50.0]

    vcat(machine_p, hydro_gov_p, shaft_p, )
end

function get_wt_simple_p()
    return [:z => 5.8 * 8 * 0.001]
end
