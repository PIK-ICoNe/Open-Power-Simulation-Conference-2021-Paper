using Nordic5
import Nordic5.wrap_node_p
using Graphs
using OrdinaryDiffEq
using Plots
using ProBeTune
using BlockSystems
using ModelingToolkit
using NetworkDynamics
using DelimitedFiles

"""
    get_g()
Network Structure of the Nordic5 Bus System, without a connection to central europe.
"""
function get_g()
    g = SimpleGraph(5)
    add_edge!(g, 1, 4)
    add_edge!(g, 2, 4)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 5)
    add_edge!(g, 1, 2)
    return g
end

"""
    generate_nordic5(P_ref_vec)
Generates the ODEFunction of the Nordic5 system.
The controllers in front of the hydro and wind nodes are fixed for these simulations. 
They return a fixed P_ref_0 which is simply the power flow in the synchronous state.
All parameters, except the proportional gains D_i at the shafts, are fixed.
"""
function generate_nordic5(P_ref_vec)
    g = get_g()
    edge = StaticLine(-100im)

    # node 1, only hydro
    hydro1_c = get_constant_Pref(:P_fix => P_ref_vec[1])
    hydro1 = get_hydro(hydro1_c, get_hydro_1_p())
    hydro1 = set_p(hydro1, :Q_inj => 0.0)
    node1 = ODEVertex(hydro1, [:P_inj, :D])

    # node 2, hydro and wind
    hydro2_c = get_constant_Pref(:P_fix => P_ref_vec[2])
    wind2_c  = get_constant_Pref(:P_fix => P_ref_vec[3])
    wind_hydro_2_para = Dict(vcat(get_hydro_2_p(), get_wt_simple_p(),
                                :Q_inj => 0.0,
                                :P_mpp => P_ref_vec[3]))
    hydro_wind = get_hydro_wind_simple(get_easy_inverter(), wind2_c, hydro2_c)
    hydro_wind = set_p(hydro_wind, wind_hydro_2_para)
    node2 = ODEVertex(hydro_wind, [:P_inj, :D])

    # node 3, only hydro
    hydro3_c = get_constant_Pref(:P_fix => P_ref_vec[4])
    hydro3 = get_hydro(hydro3_c, get_hydro_3_p())
    hydro3 = set_p(hydro3, :Q_inj => 0.0)
    node3 = ODEVertex(hydro3, [:P_inj, :D])

    # node 4, wind and thermal
    wind4_c  = get_constant_Pref(:P_fix => P_ref_vec[5])
    wind_thermal_4_para = Dict(vcat(get_thermal_4_p(), get_wt_simple_p(),
                                    :Q_inj => -0.0,
                                    :P_mpp => P_ref_vec[5],
                                    :P_m_exciter => P_ref_vec[6],
                                    :P_m_shaft => P_ref_vec[6]))
    thermal_wind = get_thermal_wind_simple(get_easy_inverter(), wind4_c)
    thermal_wind = set_p(thermal_wind, wind_thermal_4_para)
    node4 = ODEVertex(thermal_wind, [:P_inj, :D])

    # node 5, only thermal
    thermal5 = get_thermal()
    thermal5_para = Dict(vcat(get_thermal_5_p(),
                            :Q_inj => -0.0,
                            :P_m_exciter => P_ref_vec[7],
                            :P_m_shaft => P_ref_vec[7]))
    thermal5 = set_p(thermal5, thermal5_para)
    node5 = ODEVertex(thermal5, [:P_inj, :D]);

    return network_dynamics([node1, node2, node3, node4, node5], edge, g)
end

"""
   generate_spec(P_ref_vec)
ODEFunction consisting of 5 swing equations with additional proportional control at each node. 
All parameters, except the proportional gains D_i, are fixed.
"""
function generate_spec(P_ref_vec)
    g = get_g()

    @parameters t ω(t) D P_ref_0
    @variables P_ref(t)
    controller = IOBlock([P_ref ~ P_ref_0 - D * ω], [ω], [P_ref])

    edge = StaticLine(-100im)

    node1 = ODEVertex(get_controlled_swing(controller, :Q_inj=>0, :γ=>0, P_ref_0 => P_ref_vec[1], :H => 3.0), [:P_inj, :D])
    node2 = ODEVertex(get_controlled_swing(controller, :Q_inj=>0, :γ=>0, P_ref_0 => P_ref_vec[2] + P_ref_vec[3], :H => 3.0), [:P_inj, :D])
    node3 = ODEVertex(get_controlled_swing(controller, :Q_inj=>0, :γ=>0, P_ref_0 => P_ref_vec[4], :H => 3.0), [:P_inj, :D])
    node4 = ODEVertex(get_controlled_swing(controller, :Q_inj=>0, :γ=>0, P_ref_0 => P_ref_vec[5] + P_ref_vec[6], :H => 6.0), [:P_inj, :D])
    node5 = ODEVertex(get_controlled_swing(controller, :Q_inj=>0, :γ=>0, P_ref_0 => P_ref_vec[7], :H => 6.0), [:P_inj, :D])

    network_dynamics([node1, node2, node3, node4, node5], edge, g)
end

"""
    solve_sys_spec(pbt::NDProblem, in, p_sys, p_spec)
Solve the system and specification for given parameters
Parameters:
- `pbt::PBTProblem`
- `i`: input to the system
- `p_sys`: parameters of the system
- `p_spec`: parameters of specification
"""
function ProBeTune.solve_sys_spec(pbt::NDProblem, in, p_sys, p_spec)
    # injected power is base injected power + input
    P_inj = pbt.P_inj .+ in
       
    # wrap the parameters
    p_sys = wrap_node_p(P_inj, p_sys)
    
    y0_sys = convert(eltype(p_sys), pbt.y0_sys)

    prob_sys = ODEProblem(pbt.nd_sys, y0_sys, pbt.tspan, (p_sys, nothing))
    sol_sys = solve(prob_sys, Rodas4())

    p_spec = wrap_node_p(P_inj, p_spec)
    y0_spec = convert(eltype(p_spec), pbt.y0_spec)
    
    prob_spec = ODEProblem(pbt.nd_spec, y0_spec, pbt.tspan, (p_spec, nothing))
    sol_spec = solve(prob_spec, Rodas4())

    return sol_sys, sol_spec
end

"""
    save_loss_callback(p, loss)
Writes the the loss into a .txt-file after each iteration of the optimizer.
"""
function save_loss_callback(p, loss)
    output = joinpath(@__DIR__, "../data/loss.txt")
    if isfile(output)
        open(output, "a") do io
            writedlm(io, [loss])
        end  
    else
        writedlm(output, [loss])
    end
    return false
end

# Reference Power for the different nodes in [p.u.], adapted from the power flow calculations in table 1.5: https://github.com/joakimbjork/Nordic5/blob/main/Nordic5_doc.pdf 
P_ref_1_h = 1.8
P_ref_2_h = 0.9
P_ref_2_w = 0.3
P_ref_3_h = 0.3
P_ref_4_t = 0.7
P_ref_4_w = 0.3
P_ref_5_t = 0.7

P_refs = vcat(P_ref_1_h, P_ref_2_h, P_ref_2_w, P_ref_3_h, P_ref_4_w, P_ref_4_t, P_ref_5_t)
pinj = [P_ref_1_h, (P_ref_2_h, P_ref_2_w), P_ref_3_h, (P_ref_4_t, P_ref_4_w), P_ref_5_t]

# Generates the ODEFunctions for the system and specification, as well as a initial condition
nordic5 = generate_nordic5(P_refs)
x0_sys = u0_guess(nordic5.syms, pinj)

spec = generate_spec(P_refs)
x0_spec = u0_guess(spec.syms, pinj)

# Consumed Power at the nodes in [p.u.], adapted from the power flow calculations in table 1.5: https://github.com/joakimbjork/Nordic5/blob/main/Nordic5_doc.pdf 
P_load_1 = -1.8
P_load_2 = -0.65
P_load_3 = -0.2
P_load_4 = -1.6
P_load_5 = -0.75