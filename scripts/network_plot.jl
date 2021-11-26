using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using ColorSchemes
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