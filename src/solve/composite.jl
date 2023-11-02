@concrete struct CompositeNonlinearSolveAlgorithm <:
                 NonlinearSolve.AbstractNonlinearSolveAlgorithm
    nl_alg
    ode_alg
    switching_alg
end
