using DiffEqFlux
using ProBeTune
using Graphs
using OrderedCollections
using GalacticOptim
export flatten_p, unflatten_p, pbt_loss, pbt_tuning

Base.@kwdef struct NDProblem{IS,SPEC,SYS} <: AbstractPBTProblem
    p::PBTProblemParameters
    input_sampler::IS
    nd_spec::SPEC
    nd_sys::SYS
    p_spec_structure::Vector{Int}
    p_sys_structure::Vector{Int}
    y0_spec::Vector{Float64}
    y0_sys::Vector{Float64}
    P_inj::Vector{Float64}
    tspan::NTuple{2, Float64}
end

"""
    flatten_p(pbt_prob, p_sys_unflat, p_specs_unflat)

Transform two nested parameter objects `p_sys_unflat` and a vector of
`p_specs_unflat` to a flat stacked vector.

    p_sys, p_spec1, p_spec2, ...

This method should only be necessary to create the initial flat parameters!
"""
function flatten_p(pbt_prob, p_sys_unflat, p_specs_unflat)
    @assert length(p_specs_unflat) == pbt_prob.p.N_samples

    # start with flat sys parameters
    p_flat  = _flatten_p(pbt_prob.p_sys_structure,  p_sys_unflat)
    @assert pbt_prob.p.size_p_sys == length(p_flat)

    # no append for each spec
    for p_spec_unflat in p_specs_unflat
        p_spec_flat = _flatten_p(pbt_prob.p_spec_structure, p_spec_unflat)
        @assert pbt_prob.p.size_p_spec == length(p_spec_flat)
        append!(p_flat, p_spec_flat)
    end
    return p_flat
end

"""
    unflatten_p(pbt_prob, p_sys_flat, p_spec_flat)

Unflatten given flat arrays of sys and spec parameters according to
parameter structure.
Returns tuple of `p_sys_unflat` and `p_spec_unflat`.
"""
function unflatten_p(pbt_prob, p_sys_flat, p_spec_flat)
    p_sys_unflat  = _unflatten_p(pbt_prob.p_sys_structure,  p_sys_flat)
    p_spec_unflat = _unflatten_p(pbt_prob.p_spec_structure, p_spec_flat)
    (p_sys_unflat, p_spec_unflat)
end

"""
    _flatten_p(structure, p_unflat)

Check whether given structure matches the unflat parameters
and return vector of flatten Parameters.
"""
function _flatten_p(structure, p_unflat::Matrix)
    @assert all(structure .== size(p_unflat)[2])
    transpose(p_unflat)[:]
end
function _flatten_p(structure, p_unflat::Vector{T}) where T
    @assert all(length.(p_unflat) .== structure)
    vcat((p_unflat...)...)
end

"""
    _unflatten_p(structure, p_flat)

Transform the flat parameter vector to a `Vector{Vector}` according
to the `structure`.
"""
function _unflatten_p(structure, p_flat)
    @assert sum(structure) == length(p_flat)
    p = Vector{Vector{eltype(p_flat)}}(undef, length(structure))
    pos = 1
    for (i, size) in enumerate(structure)
        p[i] = p_flat[pos:pos+(size-1)]
        pos += size
    end
    return p
end

"""
    wrap_node_p(p_base, p_additional)

For each node take the base parameters and append the additional parameters.
"""
function wrap_node_p(p_base, p_additional)
    @assert length(p_base) == length(p_additional)
    N = length(p_base)
    # each vertex gets a vector of parameters, P_base first and than the others
    return [vcat(p_base[i]..., p_additional[i]...) for i in 1:N]
end

"""
    pbt_loss(pbt::NDProblem, p, p_dict; solver_options...)
pbt_loss provides the loss function for the Probabilistic Tuning Problem.
Uses multi threading to speed up the calculations.
"""
function ProBeTune.pbt_loss(pbt::NDProblem, p; solver_options...)
    losses = zeros(eltype(p), pbt.p.N_samples)

    Threads.@threads for n in 1:pbt.p.N_samples
        _pbt = deepcopy(pbt)
        sol_sys, sol_spec = solve_sys_spec_n(_pbt, n, copy(p); solver_options...)
        losses[n] = _pbt.p.output_metric(sol_sys, sol_spec)
    end

    loss = sum(losses)
    loss / pbt.p.N_samples
end

"""
    pbt_tuning(pbt::NDProblem, p, p_dict; optimizer = DiffEqFlux.ADAM(0.01), optimizer_options = (:maxiters => 100,), solver_options...)
Tune the system to the specification.
# Arguments:
- `pbt`: PBT problem
- `p`: stacked array of initial system and specification parameters
- `optimizer`: choose optimization algorithm (default `DiffEqFlux.ADAM(0.01)`)
- `optimizer_options`: choose optimization options (default `(:maxiters => 100,)`)
- `solver_options...` all further options are passed through to the differential equations solver
"""
function ProBeTune.pbt_tuning(pbt::NDProblem, p; optimizer=DiffEqFlux.ADAM(0.01), optimizer_options=(:maxiters => 100,), solver_options...)
    DiffEqFlux.sciml_train(
      x -> pbt_loss(pbt, x; solver_options...),
      p,
      optimizer, GalacticOptim.AutoForwardDiff();
      optimizer_options...
      )
end
