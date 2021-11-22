export u0_guess, plot_sol, plot_flows, plot_angles, plot_mag, plot_sol!, plot_flows!, plot_angles!, plot_mag!

node_idx(s::Symbol) = node_idx(String(s))
node_idx(s::String) = parse(Int, match(r"_([0-9]+)$", s)[1])

u0_guess(sv::Vector{Symbol}, P_ref=nothing) = Float64[u0_guess(String(s), P_ref) for s in sv]

"""
    u0_guess(s::String, P_ref)

Type-specific initial guess for the power grids using the reference power P_ref.
"""
function u0_guess(s::String, P_ref)
    if occursin(r"^u_r", s)
        return 1.0
    elseif occursin(r"^u_i", s)
        return 0.0
    elseif occursin(r"^ω", s)
        return 0.0
    elseif occursin(r"^δ", s)
        return 0.0
    elseif occursin(r"^v_q", s)
        return 1.0
    elseif occursin(r"^v_d", s)
        return 0.0
    elseif occursin(r"^gate", s) || occursin(r"^flow", s)
        idx = node_idx(s)
        if P_ref !== nothing
            return first(P_ref[idx])
        else
            println("Don't know P_ref for vertex $idx, return 1.0 for $s")
            return 1.0
        end
    elseif occursin(r"P_e_[0-9]+$", s)
        idx = node_idx(s)
        if P_ref !== nothing
            return last(P_ref[idx])
        else
            println("Don't know P_ref for vertex $idx, return 1.0 for $s")
            return 1.0
        end
    else
        println("Don't know what to do about $s zero it is!")
        return 0.0
    end
end


plot_sol(sol; kwargs...) = plot_sol!(plot(), sol; kwargs...)
plot_sol!(sol; kwargs...) = plot_sol!(plot!(), sol; kwargs...)
function plot_sol!(p, sol; sym="", idxs=:all)
    syms = sol.prob.f.syms
    mask = map(syms) do s
        s = String(s)
        if idxs !== :all && node_idx(s) ∉ idxs
            return false
        end
        if !occursin(sym, s)
            return false
        end
        return true
    end
    symidx = (1:length(syms))[mask]
    if isempty(symidx)
        println("Nothing matches!")
    else
        plot!(p, sol; vars=symidx)
    end
end

plot_angles(sol; kwargs...) = plot_angles!(plot(), sol; kwargs...)
plot_angles!(sol; kwargs...) = plot_angles!(plot!(), sol; kwargs...)
function plot_angles!(p, sol; idxs=:all, tsteps=sol.t)
    syms = sol.prob.f.syms
    nd = sol.prob.f.f
    g = nd.graph

    if idxs === :all
        idxs = 1:nv(g)
    end

    for vidx in idxs
        ur_idx = findfirst(s -> occursin(Regex("u_r_$vidx\$"), String(s)), syms)
        ui_idx = findfirst(s -> occursin(Regex("u_i_$vidx\$"), String(s)), syms)
        @assert ur_idx !== nothing
        @assert ui_idx !== nothing
        θ = zeros(length(tsteps))
        for (i, t) in enumerate(tsteps)
            solt = sol(t)
            θ[i] = atan(solt[ui_idx], solt[ur_idx])
        end
        plot!(p, tsteps, θ, label="θ at $vidx")
    end
    return p
end

plot_mag(sol; kwargs...) = plot_angles!(plot(), sol; kwargs...)
plot_mag!(sol; kwargs...) = plot_angles!(plot!(), sol; kwargs...)
function plot_mag!(p, sol; idxs=:all, tsteps=sol.t)
    syms = sol.prob.f.syms
    nd = sol.prob.f.f
    g = nd.graph

    if idxs === :all
        idxs = 1:nv(g)
    end

    for vidx in idxs
        ur_idx = findfirst(s -> occursin(Regex("u_r_$vidx\$"), String(s)), syms)
        ui_idx = findfirst(s -> occursin(Regex("u_i_$vidx\$"), String(s)), syms)
        @assert ur_idx !== nothing
        @assert ui_idx !== nothing
        v = zeros(length(tsteps))
        for (i, t) in enumerate(tsteps)
            solt = sol(t)
            v[i] = sqrt(solt[ui_idx]^2, solt[ur_idx]^2)
        end
        plot!(p, tsteps, v, label="Vm at $vidx")
    end
    return p
end

plot_flows(sol; kwargs...) = plot_flows!(plot(), sol; kwargs...)
plot_flows!(sol; kwargs...) = plot_flows!(plot!(), sol; kwargs...)
function plot_flows!(p, sol; idxs=:all, tsteps=sol.t)
    nd = sol.prob.f.f
    g = nd.graph

    if idxs === :all
        idxs = 1:ne(g)
    end

    for eidx in idxs
        edge = collect(edges(g))[eidx]
        P = zeros(length(tsteps))
        Q = zeros(length(tsteps))

        for (i, t) in enumerate(tsteps)
            gd = nd(sol(t), sol.prob.p, t, GetGD)
            edata = get_edge(gd, eidx)
            j  = Complex(edata[1], edata[2])
            srcdata = get_src_vertex(gd, eidx)
            u1 = Complex(srcdata[1], srcdata[2])
            dstdata = get_dst_vertex(gd, eidx)
            u2 = Complex(dstdata[1], dstdata[2])
            u = u2# - u2
            S = conj(j) * u
            P[i] = real(S)
            Q[i] = imag(S)
        end
        plot!(p, tsteps, P, label="P at $(edge.src) -> $(edge.dst)")
        plot!(p, tsteps, Q, label="Q at $(edge.src) -> $(edge.dst)")
    end
    return p
end
