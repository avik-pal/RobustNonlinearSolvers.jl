module RobustNonlinearSolvers

import PrecompileTools: @recompile_invalidations

@recompile_invalidations begin
    using DiffEqCallbacks, NonlinearSolve, OrdinaryDiffEq, SciMLBase
    import ConcreteStructs: @concrete
end

import NonlinearSolve: perform_step!
import SciMLBase: __init, __solve, solve!

include("solve/continuous.jl")
include("solve/composite.jl")

export ContinuousNonlinearSolveAlgorithm

end
