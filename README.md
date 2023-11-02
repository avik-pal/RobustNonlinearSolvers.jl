# RobustNonlinearSolvers

Provides an auto-switching algorithm that switches between a continuous and discrete solver
to balance between robustness and speed for solving nonlinear equations.

## Usage Example

```julia
using RobustNonlinearSolvers, NonlinearSolve, OrdinaryDiffEq

function newton_fails(u, p)
    return 0.010000000000000002 .+
           10.000000000000002 ./ (1 .+
            (0.21640425613334457 .+
             216.40425613334457 ./ (1 .+
              (0.21640425613334457 .+
               216.40425613334457 ./
               (1 .+ 0.0006250000000000001(u .^ 2.0))) .^ 2.0)) .^ 2.0) .-
           0.0011552453009332421u .- p
end

p = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
u0 = [-10.0, -1.0, 1.0, 2.0, 3.0, 4.0, 10.0]
prob = NonlinearProblem{false}(newton_fails, u0, p)

# Fails
sol1 = solve(prob, NewtonRaphson())
sol1.retcode
sol1.resid

# Fails
sol2 = solve(prob, PseudoTransient())
sol2.retcode
sol2.resid

# Tune and Pass
sol3 = solve(prob, PseudoTransient(; alpha_initial = 1.0))
sol3.retcode
sol3.resid
sol3.stats

# Continuous Solver


# Our Solver

```