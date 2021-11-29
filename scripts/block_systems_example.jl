using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using ModelingToolkit
using BlockSystems
using OrdinaryDiffEq
using Plots

@parameters t M D P_m(t) P_e(t)
@variables ω(t)
dt = Differential(t)

swing = IOBlock([dt(ω) ~ 1/M * (P_m - D*ω - P_e)], # equations
                 [P_m, P_e],                       # inputs
                 [ω])                              # outputs

swing = set_p(swing, :D=>0.5, :M=>1)

# system without pid controller, fixed pref
@variables P_m(t)
pfix = IOBlock([P_m ~ 1], [], [P_m])

wo_pid = IOSystem([pfix.P_m => swing.P_m], [pfix, swing]) |> connect_system

# system with pid controller for p_ref
@parameters input(t)
@variables int(t) out(t) pid(t)
pid = IOBlock([dt(int) ~ input,
               pid ~ input + int + dt(input),
               out ~ 1 - pid],
              [input],
              [out])

w_pid = IOSystem([pid.out => swing.P_m,
                  swing.ω => pid.input],
                 [swing, pid];
                 namespace_map=[pid.out => :P_m]
                 )
w_pid = connect_system(w_pid)

# connect step function
function simulate(sys)
    @variables P_e(t)
    step = IOBlock([P_e ~ 1 - 0.25/(1 + exp(-30*(t)))], [], [P_e])

    sys = IOSystem([step.P_e => sys.P_e], [step, sys]) |> connect_system

    gen = generate_io_function(sys)

    u0 = zeros(length(gen.states))
    odef = ODEFunction((du, u, p, t) -> gen.f_ip(du, u, nothing, p, t); mass_matrix=gen.massm, syms=Symbol.(gen.states))
    prob = ODEProblem(odef, u0, (-0.5, 7))

    sol = solve(prob, Rodas4())
end

sol1 = simulate(wo_pid)
sol2 = simulate(w_pid)

size = (500, 300)
linewidth=2
palette=:tab10
grid=false
xlabel = "Time"

p1 = plot(sol1; linewidth, legend=:none, palette, grid, xlabel, ylabel="Frequency / Power", size, title="fixed P_m")
p2 = plot(sol2; linewidth=2, legend=:right, vars=[1,3,2], palette, grid, xlabel, size, title="PID controller")
p = plot(p1, p2, size=(550,250))
savefig(p, "pid_example.pdf")
