using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using ColorSchemes
using LaTeXStrings
using Statistics
include(joinpath(@__DIR__, "../scripts/nordic5_frequency_metric.jl"))
default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)

# Loading the simulation data
loss = readdlm(joinpath(@__DIR__, "../data/loss.txt"), '\t', Float64, '\n')
ptune = readdlm(joinpath(@__DIR__, "../data/ptune.txt"), '\t', Float64, '\n')
perturbations =  readdlm(joinpath(@__DIR__, "../data/input_samples.txt"), '\t', Float64, '\n') # all perturbations
ΔP = perturbations[1,:] # choosing the first perturbation as an example for the plots
p_sys = readdlm(joinpath(@__DIR__, "../data/p_sys_start.txt"), '\t', Float64)
p_init = readdlm(joinpath(@__DIR__, "../data/p_init.txt"), '\t', Float64) # all parameters
p_spec = p_init[6:10] # parameters of the first sample
d_init = mean(readdlm(joinpath(@__DIR__, "../data/distances_init.txt"), '\t', Float64, '\n')) # Initial behavioral distance

println("The initial behavioral distance is: ", d_init)
println("The loss has been reduced by a factor of: ", loss[1]/loss[end])
#println("The final behavioral distance is: ", d_init)

plot(loss, yaxis = "Loss", xaxis = "Iterations", lw = 3, legend = false)
png("loss")

tspan = (0.0, 50.0)
ω_idx_sys = findall(map(x -> occursin("ω", string(x)), nordic5.syms))
ω_idx_spec = findall(map(x -> occursin("ω", string(x)), spec.syms))

# Simulating the untuned system and specification after a perturbation on node 4
P_inj = [P_load_1, P_load_2, P_load_3, P_load_4, P_load_5]
p = wrap_node_p(P_inj .+ ΔP, p_sys)
prob = ODEProblem(nordic5, x0_sys, tspan, (p, nothing))
sol_sys = solve(prob, Rodas4())

p = wrap_node_p(P_inj .+ ΔP, p_spec)
prob = ODEProblem(spec, x0_spec, tspan, (p, nothing))
sol_spec = solve(prob, Rodas4())

# Plotting the untuned system and specification after a perturbation on node 4
plot(sol_spec, vars = ω_idx_spec,  xaxis = L"t", yaxis = L"\omega", lw = 3, label = ["" "" "" "" "Spec"],  palette = range(ColorSchemes.leonardo[end], ColorSchemes.leonardo[1], length = 5))
Plots.plot!(sol_sys, vars = ω_idx_sys,  xaxis = L"t", lw = 3, label = ["" "" "" "" "Sys"], palette = range(ColorSchemes.berlin[1], ColorSchemes.berlin[10], length = 5))
png("untuned")

# Simulating the tuned system and specification after a perturbation on node 4
p = wrap_node_p(P_inj .+ ΔP, ptune[1:5])
prob = ODEProblem(nordic5, x0_sys, tspan, (p, nothing))
sol_sys_t = solve(prob, Rodas4())

p = wrap_node_p(P_inj .+ ΔP, ptune[6:10])
prob = ODEProblem(spec, x0_spec, tspan, (p, nothing))
sol_spec_t = solve(prob, Rodas4())

# Plotting the tuned system and specification after a perturbation on node 4
plot(sol_spec_t, vars = ω_idx_spec, xaxis = L"t", yaxis = L"\omega", lw = 3, label = ["" "" "" "" "Spec"],  palette = range(ColorSchemes.leonardo[end], ColorSchemes.leonardo[1], length = 5))
Plots.plot!(sol_sys_t, vars = ω_idx_sys, xaxis = L"t", lw = 3, label = ["" "" "" "" "Sys"], palette = range(ColorSchemes.berlin[1], ColorSchemes.berlin[10], length = 5))
png("tuned")