"""
    get_controlled_swing(@nospecialize(controller), params...)

Connects the controller and the swing equation.
"""
function get_controlled_swing(@nospecialize(controller), params...)
    spec = BlockSpec([:ω], [:P_ref])
    @assert spec(controller) "Controller block has to be ω ↦ P_ref"

    swing = _get_swing()

    @named controlled_swing = IOSystem([controller.P_ref => swing.P_ref,
                                        swing.ω => controller.ω],
                                       [swing, controller],
                                       outputs=[swing.u_r, swing.u_i])
    con = connect_system(controlled_swing)
    return set_p(con, Dict(params...))
end

"""
    _get_swing()
"""
function _get_swing()
    @parameters t P_ref(t) γ H i_r(t) i_i(t) P(t)
    @variables ω(t) u_r(t) u_i(t)
    dt = Differential(t)

    @named shaft = IOBlock([dt(ω) ~ 1/(2H) * (P_ref - γ * ω - P),
                            dt(u_r) ~ - u_i * ω,
                            dt(u_i) ~   u_r * ω],
                           [P, P_ref], [u_r, u_i, ω])

    power = _get_power()
    internal_current = get_internal_current()

    @named swing = IOSystem([power.P => shaft.P,
                             shaft.u_r => power.u_r,
                             shaft.u_i => power.u_i,
                             shaft.u_r => internal_current.u_r,
                             shaft.u_i => internal_current.u_i,
                             internal_current.i_int_r => power.i_r,
                             internal_current.i_int_i => power.i_i,
                             ],
                            [power, shaft, internal_current];
                            outputs=[shaft.u_r, shaft.u_i, shaft.ω])

    connect_system(swing)
end

"""
    _get_power()
"""
function _get_power()
    @parameters t u_r(t) u_i(t) i_r(t) i_i(t)
    @variables P(t)

    @named power = IOBlock([P ~ i_r * u_r + i_i * u_i], [u_r, u_i, i_r, i_i], [P])
end