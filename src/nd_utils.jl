using NetworkDynamics
using ModelingToolkit: getname, value
using BlockSystems

function NetworkDynamics.ODEVertex(iob::IOBlock, p_order=[])
    if isempty(iob.inputs)
        return _slack_vertex(iob)
    end

    spec = BlockSpec([:i_r, :i_i], [:u_r, :u_i])
    @assert spec(iob) "Block has to follow PowerDynamics i/o conventions!"
    @assert length(p_order) == length(iob.iparams) "Provide order of all iparams!"

    # parameters may be given as Num oder Symbol types
    gen = generate_io_function(iob,
                               f_states=[iob.u_r, iob.u_i],
                               f_inputs=[iob.i_r, iob.i_i],
                               f_params=p_order, warn=false,
                               type=:ode);

    f = (du, u, edges, p, t) -> gen.f_ip(du, u, flowsum(edges), p, t)
    vars = Symbol.(gen.states)
    ODEVertex(f = f, dim = length(vars), sym = vars, mass_matrix = gen.massm)
end

function _slack_vertex(iob::IOBlock)
    spec = BlockSpec([], [:u_r, :u_i])
    @assert spec(iob) "Block has to have :u_r and :u_i as outputs!"

    # parameters may be given as Num oder Symbol types
    gen = generate_io_function(iob,
                               f_states=[iob.u_r, iob.u_i],
                               f_inputs=[],
                               f_params=[], warn=false,
                               type=:ode);

    f = (du, u, edges, p, t) -> gen.f_ip(du, u, nothing, p, t)
    vars = Symbol.(gen.states)
    ODEVertex(f = f, dim = length(vars), sym = vars, mass_matrix = gen.massm)
end

# allocation free oop aggregator. might be more difficult for more-dimensional systems
# unfortunately there are no vector variables in MDK and we can't model the aggregators
# as an IOSystem.
function flowsum(edges)
    i_r, i_i = 0.0, 0.0
    for e in edges
        i_r += e[1]
        i_i += e[2]
    end
    return (i_r, i_i)
end

"""
    StaticLine(Y)

A static model that represents a line with an admittance Y, defined like in PowerDynamics.jl.
"""
function StaticLine(Y)
    f = (i, v_s, v_d, p, t) -> begin
        u_s = v_s[1] + im*v_s[2]
        u_d = v_d[1] + im*v_d[2]
        icomplex = Y * (u_d - u_s)
        i[1] = real(icomplex)
        i[2] = imag(icomplex)
        i[3] = -real(icomplex)
        i[4] = -imag(icomplex)
    end
    StaticEdge(f = f, dim=4, sym=[:i_r, :i_i, :i_r, :i_i], coupling=:fiducial)
end

function get_current(nd, state, node)
    gd = nd(state, nothing, 0.0, GetGD)
    total_current(get_dst_edges(gd, node))
end

function total_current(edges)
    # Keeping with the convention of negative sign for outgoing currents
    current = 0.0im
    for e in edges
        current -= e[1] + e[2]*im
    end
    current
end

get_u0(nd, dict) = getindex.(Ref(dict), nd.syms)