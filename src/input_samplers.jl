"""
    single_node_perturbation(rng, n)

Choses a random node in the network i and samples a normally distributed power perturbation a.
"""
function single_node_perturbation(rng, n)
    N = 5
    i = rand(rng, 1:N)
    a = zeros(N)
    a[i] = randn(rng)
    return a
end
