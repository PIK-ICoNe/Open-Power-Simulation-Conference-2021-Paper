using DiffEqFlux
using ProBeTune
using Graphs
using OrderedCollections

Base.@kwdef struct NDProblem{IS,SPEC,SYS} <: AbstractPBTProblem
    p::PBTProblemParameters
    input_sampler::IS
    nd_spec::SPEC
    nd_sys::SYS
    y0_spec::Vector{Float64}
    y0_sys::Vector{Float64}
    P_inj::Vector{Float64}
    tspan::NTuple{2, Float64}
end

"""
    make_p_dict(paras)
Creates a dict which gives the position of a parameter, originally given in an array of arrays, one for each node, in a flat array.
"""
function make_p_dict(paras)
    p_dict = OrderedDict{Int, Vector{Int}}()
    counter = 1 # gives the position in the flat array
    for node in 1:length(paras) 
        if paras[node] == [] # accounting for nodes which don't have tuneable parameters
            p_dict[node] = (counter, 0.0) # we will fill empty parameter slots with zeros 
            counter += 1
        elseif length(paras[node]) > 1 # nodes with multiple tuneable parameters
            dict_entry = []
            for j in 1:length(paras[node])
                push!(dict_entry, counter)
                counter += 1
            end
            p_dict[node] = dict_entry
        else # nodes with only one tuneable parameter
            p_dict[node] = [counter]
            counter += 1
        end
    end
    return p_dict
end


"""
    wrap_node_p(p_base, p_additional)
Wraps parameters for the simulation. p_base contains the perturbed parameters while p_additional accounts for the tuneable parameters.
"""
function wrap_node_p(p_base, p_additional)
    N = length(p_base)

    # each vertex gets a tuple of parameters, P_inj first and than the others
    return [(vcat(p_base[i], p_additional[i]...)) for i in 1:N]
end

"""
    get_ordered(dict, unordered)
Orders the flat array "unordered" and orders it for the ODEFunction using dict.
"""
function get_ordered(dict, unordered)
    ordered = Vector{Vector{Any}}(undef, 0)
    for key in keys(dict)
        if typeof(dict[key]) == Tuple{Int64, Float64}
            push!(ordered, [0.0]) # adding a zero for nodes which don't have a tuneable parameter
        else
            push!(ordered, [unordered[val] for val in dict[key]])
        end
    end
    return ordered
end

"""
    solve_sys_spec_n(pbt::NDProblem, n::Int, p, p_dict; solver_options...)
Solve system and spec for the nth input sample given a stacked parameter array p and the parameter dict p_dict.
Parameters:
- `pbt`: existing PBTProblem instance
- `n`: number of input in the array
- `p`: combined array of sys and spec parameters
- `p_dict`: Dict giving the positions of the parameters in a flat array and with respect to an ordered array for NetworkDynamics.jl
"""
function ProBeTune.solve_sys_spec_n(pbt::NDProblem, n::Int, p, p_dict; solver_options...)
    i = pbt.input_sampler[n]

    p = get_ordered(p_dict, p) # orders the flat parameter array

    num_node_sys = nv(pbt.nd_sys.f.graph)

    p_sys = p[1:num_node_sys] # each node has its own array in the parameter array
    p_spec = p[num_node_sys + 1:end]

    solve_sys_spec(pbt, i, p_sys, p_spec; solver_options...)
end


"""
    pbt_loss(pbt::NDProblem, p, p_dict; solver_options...)
pbt_loss provides the loss function for the Probabilistic Tuning Problem.
"""
function ProBeTune.pbt_loss(pbt::NDProblem, p, p_dict; solver_options...)
    loss = 0.0
    for n in 1:pbt.p.N_samples
        sol_sys, sol_spec = solve_sys_spec_n(pbt, n, p, p_dict; solver_options...)
        loss += pbt.p.output_metric(sol_sys, sol_spec)
    end

    loss / pbt.p.N_samples
end


"""
    pbt_tuning(pbt::NDProblem, p, p_dict; optimizer = DiffEqFlux.ADAM(0.01), optimizer_options = (:maxiters => 100,), solver_options...)
Tune the system to the specification.
# Arguments:
- `pbt`: PBT problem
- `p`: stacked array of initial system and specification parameters
- `p_dict`: Dict giving the positions of the parameters in a flat array and with respect to an ordered array for NetworkDynamics.jl
- `optimizer`: choose optimization algorithm (default `DiffEqFlux.ADAM(0.01)`)
- `optimizer_options`: choose optimization options (default `(:maxiters => 100,)`)
- `solver_options...` all further options are passed through to the differential equations solver
"""
function ProBeTune.pbt_tuning(pbt::NDProblem, p, p_dict; optimizer = DiffEqFlux.ADAM(0.01), optimizer_options = (:maxiters => 100,), solver_options...)
    DiffEqFlux.sciml_train(
      x -> pbt_loss(pbt, x, p_dict; solver_options...),
      p,
      optimizer;
      optimizer_options...
      )
end
