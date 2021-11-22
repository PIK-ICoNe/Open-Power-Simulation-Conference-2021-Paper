using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using ColorSchemes
using LaTeXStrings
import Nordic5.get_ordered
include(joinpath(@__DIR__, "../scripts/nordic5_frequency_metric.jl"))
default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)



# Loading the simulation data
loss = readdlm(joinpath(@__DIR__, "../data/loss_25.txt"), '\t', Float64, '\n')
p_sys = readdlm(joinpath(@__DIR__, "../data/p_sys_start.txt"), '\t', Float64)
p_spec = readdlm(joinpath(@__DIR__, "../data/p_spec_start.txt"), '\t', Float64)
ptune = readdlm(joinpath(@__DIR__, "../data/ptune_25.txt"), '\t', Float64, '\n')
p_dict = make_p_dict(vcat(p_sys, p_spec))
p_flat = vcat(vcat(p_sys, p_spec)...)

# Plotting the loss over the optimizer iterations
plot(loss, yaxis = "Loss", xaxis = "Iter", lw = 3, legend = false)
png("loss")


# Simulation time span
tspan = (0.0, 50.0)
ω_idx_sys = findall(map(x -> occursin("ω", string(x)), nordic5.syms))
ω_idx_spec = findall(map(x -> occursin("ω", string(x)), spec.syms))

# Initial balanced Power Flow
P_inj = [P_load_1, P_load_2, P_load_3, P_load_4, P_load_5]

# Power Perturbation on node 4
ΔP = [0.0,0.0,0.0,-1,0.0] 

# Simulation of the untuned system after a perturbation on node 4
p = wrap_node_p(P_inj .+  ΔP, p_sys)
prob = ODEProblem(nordic5, x0_sys, tspan, (p, nothing))
sol_sys = solve(prob, Rodas4())

p = wrap_node_p(P_inj .+ ΔP, p_spec)
prob = ODEProblem(spec, x0_spec, tspan, (p, nothing))
sol_spec = solve(prob, Rodas4())

# Simulation of the tuned system after a perturbation on node 4
p_t = get_ordered(p_dict, ptune)
p = wrap_node_p(P_inj .+ ΔP, p_t[1:5])
prob = ODEProblem(nordic5, x0_sys, tspan, (p, nothing))
sol_sys_t = solve(prob, Rodas4())

p = wrap_node_p(P_inj .+ ΔP, p_t[6:10])
prob = ODEProblem(spec, x0_spec, tspan, (p, nothing))
sol_spec_t = solve(prob, Rodas4())

# Plotting the frequency transients of the untuned system
plot(sol_spec, vars = ω_idx_spec,  xaxis = L"t [s]", yaxis = L"\omega [rad/s]", lw = 3, title = "Untuned System", label = ["" "" "" "" "Spec"],  palette = range(ColorSchemes.leonardo[end], ColorSchemes.leonardo[1], length = 5))
plot!(sol_sys, vars = ω_idx_sys,  xaxis = L"t[s]", lw = 3, label = ["" "" "" "" "Sys"], palette = range(ColorSchemes.berlin[1], ColorSchemes.berlin[10], length = 5))
png("untuned")

# Plotting the frequency transients of the tuned system
plot(sol_spec_t, vars = ω_idx_spec, xaxis = L"t", yaxis = L"\omega [rad/s]", lw = 3, title = "Tuned System", label = ["" "" "" "" "Spec"],  palette = range(ColorSchemes.leonardo[end], ColorSchemes.leonardo[1], length = 5))
plt = plot!(sol_sys_t, vars = ω_idx_sys, xaxis = L"t[s]", lw = 3, label = ["" "" "" "" "Sys"], palette = range(ColorSchemes.berlin[1], ColorSchemes.berlin[10], length = 5))
png("tuned")

# Graph Plot of the Nordic5 System
using GraphMakie, CairoMakie, Colors
using GraphMakie.NetworkLayout, FileIO

g = SimpleGraph(9)
add_edge!(g, 1, 4)
add_edge!(g, 2, 4)
add_edge!(g, 2, 3)
add_edge!(g, 3, 5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 6)
add_edge!(g, 2, 7)
add_edge!(g, 4, 8)
add_edge!(g, 4, 9)

node_color = [:mediumturquoise, :black, :mediumturquoise,  :black, :gold, :royalblue4, :mediumturquoise,:royalblue4, :gold]
node_marker = [:circle, :hline, :circle,  :hline, :circle, :circle, :circle,:circle, :circle]

fixed_layout(_) = [(0,0), (1.5, 2), (3, 1),  (1.5, -1), (3, 0), (1.4, 2.1), (1.6, 2.1), (1.4, -1.1), (1.6, -1.1)]

node_size = ones(9) * 25
node_size[2] = 30
node_size[4] = 30

f, ax, p = graphplot(g, edge_width = 3.0, layout=fixed_layout, node_size = node_size, node_marker = node_marker, node_color = node_color)

hidedecorations!(ax)
hidespines!(ax)
ax.aspect = DataAspect()

FileIO.save(joinpath(@__DIR__, "network.svg"), f, resolution = (1000, 1000))