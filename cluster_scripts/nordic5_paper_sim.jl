using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
include(joinpath(@__DIR__,"../scripts/nordic5_frequency_metric.jl"))

pbt_prob = NDProblem(;p = PBTProblemParameters(N_samples = 25,
                      size_p_spec = 0,
                      size_p_sys = 0,
                      output_metric = frequency_metric),
                     input_sampler = ConstantInputSampler(single_node_perturbation),
                     nd_spec = spec,
                     nd_sys = nordic5,
                     y0_spec = x0_spec,
                     y0_sys = x0_sys,
                     P_inj = [P_load_1, P_load_2, P_load_3, P_load_4, P_load_5],
                     tspan = (0.0, 50.0))

p_sys = readdlm(joinpath(@__DIR__, "../data/p_sys_start.txt"), '\t', Float64)
p_spec = readdlm(joinpath(@__DIR__, "../data/p_spec_start.txt"), '\t', Float64)
p_dict = make_p_dict(vcat(p_sys, p_spec))
p_flat = vcat(vcat(p_sys, p_spec)...)

ptune = pbt_tuning(pbt_prob, p_flat, p_dict; optimizer_options=(:maxiters=>300, :cb=> save_loss_callback))

writedlm(joinpath(@__DIR__, "../data/p_input_sampler_cluster.txt"), pbt_prob.input_sampler._cache)
writedlm(joinpath(@__DIR__, "../data/ptune_25_cluster.txt"), ptune)
