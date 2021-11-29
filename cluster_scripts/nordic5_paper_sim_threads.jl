using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
include(joinpath(@__DIR__,"../scripts/nordic5_frequency_metric.jl"))

@info "Julia runs with $(Threads.nthreads()) threads!"

N_samples = 10

pbt_prob = NDProblem(;p = PBTProblemParameters(N_samples = N_samples,
                                               size_p_spec = 5,
                                               size_p_sys = 5,
                                               output_metric = frequency_metric),
                     input_sampler = ConstantInputSampler(single_node_perturbation),
                     nd_spec = spec,
                     nd_sys = nordic5,
                     p_spec_structure = [1,1,1,1,1],
                     p_sys_structure  = [1,1,1,1,1],
                     y0_spec = x0_spec,
                     y0_sys = x0_sys,
                     P_inj = [P_load_1, P_load_2, P_load_3, P_load_4, P_load_5],
                     tspan = (0.0, 50.0));

p_sys = readdlm(joinpath(@__DIR__, "../data/p_sys_start.txt"), '\t', Float64)
p_spec = readdlm(joinpath(@__DIR__, "../data/p_spec_start.txt"), '\t', Float64)

p_flat = flatten_p(pbt_prob, p_sys, repeat([p_spec], N_samples, 1))

dist_init, p_init, distances_init  = behavioural_distance(pbt_prob, p_flat; verbose=true, optimizer_options=(:maxiters => 10,))
#writedlm(joinpath(@__DIR__, "../data/p_init.txt"), p_init)
#writedlm(joinpath(@__DIR__, "../data/distances.txt"), distances_init)

ptune = pbt_tuning(pbt_prob, p_init; optimizer_options=(:maxiters=>250, :cb=> save_loss_callback))
#writedlm(joinpath(@__DIR__, "../data/input_samples.txt"), pbt_prob.input_sampler._cache)
#writedlm(joinpath(@__DIR__, "../data/p_tune.txt"), ptune)

dist_end, p_end, distances_end = behavioural_distance(pbt_prob, ptune; verbose=true, optimizer_options=(:maxiters => 1000,))

#writedlm(joinpath(@__DIR__, "../data/p_end.txt"), p_end)
#writedlm(joinpath(@__DIR__, "../data/distances_end.txt"), distances_end)
