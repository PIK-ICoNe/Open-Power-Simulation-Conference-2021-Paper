"""
   get_hydro(fcr_controller::IOBlock, p=nothing)

Plugs the different models for the hydro nodes together. The hydro nodes consist of:
    - Fifth Order Machine Model + Shaft
    - Governor + FCR Controller
    - PQ Load
"""
function get_hydro(fcr_controller::IOBlock, p=nothing)
   # Loading the different components
   machine_fifth = get_machine_fifth()
   shaft = get_shaft()
   hydro_gov = get_hydro_gov(fcr_controller)
   rfc = get_rfc()
   PQ = get_internal_current()
   
   hydro_no_rfc = IOSystem([
      shaft.ω           => hydro_gov.ω,
      machine_fifth.P_e => shaft.P_e,
      hydro_gov.P_m     => shaft.P_m],
      [machine_fifth, hydro_gov, shaft],
      namespace_map = [
            shaft.ω => :ω,
            shaft.δ => :δ,
            # Shaft parameters
            shaft.D => :D,
            shaft.H => :H,
            shaft.Ω_H => :Ω_H,
            # Machine Parameter names
            machine_fifth.X_d   => :X_d,
            machine_fifth.X_q′  => :X_q′,
            machine_fifth.X_q′′ => :X_q′′,
            machine_fifth.X_d′  => :X_d′,
            machine_fifth.X_q′′ => :X_q′′,
            machine_fifth.T_q′′ => :T_q′′,
            machine_fifth.T_d′  => :T_d′,
            machine_fifth.T_d′′ => :T_d′′,
            machine_fifth.E_f => :E_f,],
      outputs = [machine_fifth.v_d, machine_fifth.v_q, shaft.ω, shaft.δ],
      name=:hydro_no_rfc)

   hydro = IOSystem([
      rfc.i_d          => hydro_no_rfc.i_d,
      rfc.i_q          => hydro_no_rfc.i_q,
      hydro_no_rfc.v_d => rfc.v_d,
      hydro_no_rfc.v_q => rfc.v_q,
      hydro_no_rfc.δ   => rfc.δ,
      rfc.u_r          => PQ.u_r,
      rfc.u_i          => PQ.u_i,
      PQ.i_int_r       => rfc.i_r,
      PQ.i_int_i       => rfc.i_i],
      [rfc, hydro_no_rfc, PQ],
      outputs = [rfc.u_r, rfc.u_i, hydro_no_rfc.ω, hydro_no_rfc.δ], 
      name = :hydro)
   hydro = connect_system(hydro)
   isnothing(p) ? hydro :  set_p(hydro, p)
end

"""
   get_hydro_gov(fcr_controller::IOBlock)

Governor connected to the hydro nodes. Depends on the fcr_controller which gives the reference power P_ref.
"""
function get_hydro_gov(fcr_controller::IOBlock)
    @parameters t
    dt = Differential(t)

    @parameters P_ref(t) T_y T_w i(t) i1(t) i2(t)
    @variables gate(t) P_m(t) o(t) flow(t)

    servo = IOBlock([dt(gate) ~ (1 / T_y) * (P_ref - gate)], [P_ref], [gate], name = :servo) # H(s) = 1 / (T_y * s + 1)
    square = IOBlock([o ~ i^2], [i], [o], name = :square)
    water_lag = IOBlock([dt(flow) ~ (1 / T_w) * i], [i], [flow], name = :water_lag) # H(s) = 1 / (Tw * s)
    mult = IOBlock([P_m ~ i1 * i2], [i1, i2], [P_m], name = :mult)

    hydro_gov = IOSystem([fcr_controller.P_ref        => servo.P_ref,
                          water_lag.flow / servo.gate => square.i,
                          square.o                    => mult.i1,
                          1.0 - square.o              => water_lag.i,
                          water_lag.flow              => mult.i2],
                         [fcr_controller, servo, square, water_lag, mult],
                         outputs = [mult.P_m],
                         name = :hydro_gov)
    return connect_system(hydro_gov)
end

"""
   get_machine_fifth()

Fifth order machine model for the hydro nodes. 
Voltage equations of the salient rotor using the fifth order model from J. Machowski 'Power System Dynamics' (section 11.103 on page 455).
"""
function get_machine_fifth()
   @parameters t
   dt = Differential(t)

   @parameters T_q′′ T_d′ T_d′′ X_d X_d′ X_q′ X_d′′ X_q′′ E_f i_q(t) i_d(t)
   @variables v_d(t) v_q(t) v_d′(t) v_q′(t) P_e(t)

   machine_fifth = IOBlock([
                     dt(v_q)  ~ (1/T_d′)  * (E_f - v_q  + (X_d - X_d′) * i_d ),
                     dt(v_q′) ~ (1/T_d′′) * (v_q - v_q′ + (X_d′ - X_d′′) * i_d),
                     dt(v_d′) ~ (1/T_q′′) * (v_d - v_d′ - (X_q′ - X_q′′) * i_q),
                     v_d ~ 0,
                     P_e ~ (v_d′ * i_d + v_q′ * i_q) + (X_d′′ - X_q′′) * i_d * i_q ],
                     [i_d, i_q], [v_d, v_q, P_e], name = :machine_fifth)

   return machine_fifth
end
