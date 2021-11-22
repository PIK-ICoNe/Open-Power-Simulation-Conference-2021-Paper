"""
    get_hydro_wind_simple(inverter::IOBlock, ffr_controller::IOBlock, fcr_controller::IOBlock)

Plugs the hydro model and the simple wind turbine together.
"""
function get_hydro_wind_simple(inverter::IOBlock, ffr_controller::IOBlock, fcr_controller::IOBlock)
    hydro = get_hydro(fcr_controller)
    wind = get_simple_wind_turbine(inverter, ffr_controller)
    current_sum = get_current_sum()

    hydro_wind = IOSystem([
        current_sum.i_r_in => hydro.i_r,
        current_sum.i_i_in => hydro.i_i,
        wind.i_r           => current_sum.i_r_wind,
        wind.i_i           => current_sum.i_i_wind,
        hydro.u_r          => wind.u_r,
        hydro.u_i          => wind.u_i,
        hydro.ω            => wind.ω],
        [wind, hydro, current_sum],
        namespace_map = [current_sum.i_r => :i_r, current_sum.i_i => :i_i],
        outputs = [hydro.u_r, hydro.u_i, hydro.δ, hydro.ω],
        name = :hydro_wind)

    return connect_system(hydro_wind)
end

"""
    get_thermal_wind_simple(inverter::IOBlock, ffr_controller::IOBlock)

Plugs the thermal model and the simple wind turbine together.
"""
function get_thermal_wind_simple(inverter::IOBlock, ffr_controller::IOBlock)
    thermal = get_thermal()
    wind = get_simple_wind_turbine(inverter, ffr_controller)
    current_sum = get_current_sum()

    thermal_wind = IOSystem([
        current_sum.i_r_in   => thermal.i_r,
        current_sum.i_i_in   => thermal.i_i,
        wind.i_r             => current_sum.i_r_wind,
        wind.i_i             => current_sum.i_i_wind,
        thermal.u_r          => wind.u_r,
        thermal.u_i          => wind.u_i,
        thermal.ω            => wind.ω],
        [wind, thermal, current_sum],
        namespace_map = [current_sum.i_r => :i_r, current_sum.i_i => :i_i],
        outputs = [thermal.u_r, thermal.u_i, thermal.δ, thermal.ω],
        name = :thermal_wind)

    return connect_system(thermal_wind)
end
