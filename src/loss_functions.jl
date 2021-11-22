"""
    frequency_metric(sol_sys, sol_spec)

Difference between the frequency transients of the system and specification.
"""
function frequency_metric(sol_sys, sol_spec)
    ω_idx_sys = idx_containing(sol_sys.prob.f, "ω")
    ω_idx_spec = idx_containing(sol_spec.prob.f, "ω")

    if length(ω_idx_sys) != length(ω_idx_spec)
        error("System and Spec need the same number of oscillators")
    end

    # Only successful simulations are compared
    if sol_sys.retcode == :Success && sol_spec.retcode == :Success
        s = 0.0

        #  runs over all nodes in the networks
        for i in 1:length(ω_idx_sys)
            # Compares 100 points in the trajectory
            for t in range(sol_sys.prob.tspan...; length=100)
                s += (sol_sys(t)[ω_idx_sys[i]] - sol_spec(t)[ω_idx_spec[i]])^2
            end
        end
        return s
    else
        @warn "The solution of either sys or spec was unsuccessful!"
        return Inf # Solvers failing is bad.
    end
end
