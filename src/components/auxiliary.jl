export get_constant_Pref
"""
    get_shaft()

Shaft of the machine with an additional proportional control term.
"""
function get_shaft()
    @parameters t
    dt = Differential(t)

    @parameters P_m(t) P_e(t) D H Ω_H
    @variables δ(t) ω(t)
    shaft = IOBlock([
        dt(δ) ~ ω * Ω_H,
        dt(ω) ~ 1/(2H) * (P_m - P_e - D * ω)],
        [P_m, P_e],
        [δ, ω], name=:shaft)

    return shaft
end

"""
    get_shaft_thermal()

Shaft of the thermal machine with an additional proportional control term.
The mechanical power P_m is a parameter as the thermal nodes don't inherit governors. 
"""
function get_shaft_thermal()
    @parameters t
    dt = Differential(t)

    @parameters P_m P_e(t) D H Ω_H
    @variables δ(t) ω(t)
    shaft_thermal = IOBlock([
        dt(δ) ~ ω * Ω_H,
        dt(ω) ~ 1/(2H) * (P_m - P_e - D * ω)],
        [P_e],
        [δ, ω], name= :shaft_thermal)

    return shaft_thermal
end

"""
    get_rfc()

Reference Frame Conversion
"""
function get_rfc()
    @parameters t

    @parameters δ(t) v_d(t) v_q(t) i_r(t) i_i(t)
    @variables i_d(t) i_q(t) u_r(t) u_i(t) v_h(t) 

    rfc = IOBlock([
        u_r ~ v_q * cos(δ)  + v_d * sin(δ),
        u_i ~ v_q * sin(δ)  - v_d * cos(δ),
        i_q ~ i_r * cos(-δ) - i_i * sin(-δ),
        i_d ~ i_r * sin(-δ) + i_i * cos(-δ)],
        [δ, v_d, v_q, i_r, i_i], 
        [i_d, i_q, u_r, u_i], 
        name=:rfc)
    return rfc
end

function get_internal_current()
    @parameters t u_r(t) u_i(t) i_r(t) i_i(t) P_inj Q_inj
    @variables i_int_r(t) i_int_i(t)

    @named inj = IOBlock([i_int_r ~ i_r - real((P_inj + im*Q_inj)/(u_r + im*u_i)),
                          i_int_i ~ i_i + imag((P_inj + im*Q_inj)/(u_r + im*u_i))],
                         [u_r, u_i, i_r, i_i], [i_int_r, i_int_i])
end

function get_current_sum()
    @parameters t i_r(t) i_i(t) i_r_wind(t) i_i_wind(t)
    @variables i_r_in(t) i_i_in(t)

    IOBlock([i_r_in ~ i_r - i_r_wind, i_i_in ~ i_i - i_i_wind], [i_r, i_i, i_r_wind, i_i_wind], [i_r_in, i_i_in], name=:current_sum)
end

function get_constant_Pref(p...)
    @parameters t P_fix ω(t)
    @variables P_ref(t)
    blk = IOBlock([P_ref ~ P_fix], [ω], [P_ref]; name=:Pfix)
    set_p(blk, Dict(p))
end