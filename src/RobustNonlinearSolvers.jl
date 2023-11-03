module RobustNonlinearSolvers

import PrecompileTools: @recompile_invalidations

@recompile_invalidations begin
    using DiffEqCallbacks, NonlinearSolve, OrdinaryDiffEq, SciMLBase
    import ConcreteStructs: @concrete
end

import NonlinearSolve: perform_step!
import SciMLBase: init, __init, solve, __solve, solve!, NLStats

abstract type AbstractSwitchingAlgorithm end

include("solve/continuous.jl")
include("solve/composite.jl")
include("switching.jl")

export ContinuousNonlinearSolveAlgorithm, CompositeNonlinearSolveAlgorithm
export PseudoTransientSwitching

end
