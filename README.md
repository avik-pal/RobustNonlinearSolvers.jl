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
sol1.stats

# Fails
sol2 = solve(prob, PseudoTransient())
sol2.retcode
sol2.resid
sol2.stats

# Tune and Pass
sol3 = solve(prob, PseudoTransient(; alpha_initial = 1.0))
sol3.retcode
sol3.resid
sol3.stats

# Continuous Solver: Adaptive Time Stepping
sol4 = solve(prob, ContinuousNonlinearSolveAlgorithm(ImplicitEuler()))
sol4.retcode
sol4.resid
sol4.stats

# Continuous Solver: Fixed Time Stepping -- fails
sol5 = solve(prob, ContinuousNonlinearSolveAlgorithm(Euler()); dt = 1.0)
sol5.retcode
sol5.resid
sol5.stats

# Auto Switching
sol6 = solve(prob, CompositeNonlinearSolveAlgorithm())
sol6.retcode
sol6.resid
sol6.stats
```

```julia
function quadratic_f!(du, u, p)
    du .= abs2.(u) .- p
end

u0 = -ones(1000) .* 10.0
prob = NonlinearProblem(quadratic_f!, u0, 2.0)

# PT -- doesn't work
sol1 = solve(prob, PseudoTransient())

# PT -- again doesn't work
sol2 = solve(prob, PseudoTransient(; alpha_initial = 10.0))

# Newton -- works
sol3 = solve(prob, NewtonRaphson())

# Implicit Euler -- doesn't work
sol4 = solve(prob, ContinuousNonlinearSolveAlgorithm(ImplicitEuler()))

# Tsit5 -- doesn't works
sol5 = solve(prob, ContinuousNonlinearSolveAlgorithm(Tsit5()))

# Auto Switching -- works
sol6 = solve(prob, CompositeNonlinearSolveAlgorithm())
```
