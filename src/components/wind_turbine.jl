"""
    get_simple_wind_gov()

Adapted Wind Turbine Design from [1]. Allows for meaningful reactions after over frequencies.

[1] J. Björk et. al., Variable-Speed Wind Turbine Control Designed for Coordinated Fast Frequency Reserves, 2021
"""
function get_simple_wind_gov()
    @parameters t P_ref(t) z P_mpp
    dt = Differential(t)
    @variables P_e(t) ΔP_e(t) ΔP(t)
    IOBlock([dt(P_e) ~ max(0.0, dt(P_ref)) - max(P_e - P_ref, (P_ref + P_e - 2 * P_mpp) * z)], [P_ref], [P_e], name = :wt_dc)
end


"""
    get_simple_wind_turbine(inverter::IOBlock, ffr_controller::IOBlock)

Connects the wind turbine the inverter and the ffr controller in front of the wind turbine.
"""
function get_simple_wind_turbine(inverter::IOBlock, ffr_controller::IOBlock)
    gov = get_simple_wind_gov()

    gov_ffr = IOSystem([ffr_controller.P_ref => gov.P_ref], [ffr_controller, gov]) |> connect_system

    wt_simple = IOSystem([gov_ffr.P_e => inverter.P],
                         [gov_ffr, inverter],
                         outputs = [inverter.i_r, inverter.i_i], name = :wt_simple)
    return connect_system(wt_simple)
end

"""
    get_easy_inverter()

Assuming a grid following inverter: i = P / u* which is simply a voltage following current injector
"""
function get_easy_inverter()
    @parameters t

    @parameters u_r(t) u_i(t) P(t)
    @variables v(t) i_r(t) i_i(t)
    return IOBlock([v ~ u_r^2 + u_i^2, i_r ~ (P * u_r) / v , i_i ~ (P * u_i) / v], [u_r, u_i, P], [i_r, i_i], name = :easy_inverter)
end
