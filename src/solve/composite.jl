@concrete struct CompositeNonlinearSolveAlgorithm <:
                 NonlinearSolve.AbstractNonlinearSolveAlgorithm
    nl_alg
    ode_alg
    switching_alg
end

function CompositeNonlinearSolveAlgorithm(; nl_alg = NewtonRaphson(),
    ode_alg = Tsit5(), switching_alg = PseudoTransientSwitching())
    ode_alg = ode_alg isa ContinuousNonlinearSolveAlgorithm ?
              ode_alg : ContinuousNonlinearSolveAlgorithm(ode_alg)
    return CompositeNonlinearSolveAlgorithm(nl_alg, ode_alg, switching_alg)
end

function __init(prob::NonlinearProblem, alg::CompositeNonlinearSolveAlgorithm, args...;
    reltol = nothing, abstol = nothing, termination_condition = nothing,
    maxiters::Int = 1000, kwargs...)
    ode_alg_cache = __init(prob, alg.ode_alg; kwargs..., maxiters, reltol, abstol,
        termination_condition)
    nlsolve_alg_cache = __init(prob, alg.nl_alg; kwargs..., maxiters, reltol, abstol,
        termination_condition)

    switching_alg_cache = init_switching_algorithm(alg.switching_alg,
        NonlinearSolve.get_fu(ode_alg_cache); kwargs...)

    return CompositeNonlinearSolveCache{isinplace(prob)}(prob, alg, prob.f,
        nlsolve_alg_cache.fu1, prob.u0,
        2, ode_alg_cache, nlsolve_alg_cache, abstol, reltol, false, maxiters,
        ReturnCode.Default, nlsolve_alg_cache.tc_cache, NLStats(1, 0, 0, 0, 0),
        alg.switching_alg, switching_alg_cache)
end

@concrete mutable struct CompositeNonlinearSolveCache{S, iip} <:
                         NonlinearSolve.AbstractNonlinearSolveCache{iip}
    prob
    alg
    f
    fu1
    u
    current::Int
    ode_alg_cache
    nlsolve_alg_cache
    abstol
    reltol
    force_stop::Bool
    maxiters::Int
    retcode::ReturnCode.T
    tc_cache
    stats
    switching_alg::S
    switching_alg_cache
end

# TODO: NonlinearSolve.jl needs a proper interface for interator interface
function perform_step!(cache::CompositeNonlinearSolveCache{<:PseudoTransientSwitching})
    if cache.current == 1
        perform_step!(cache.nlsolve_alg_cache)
        cache.retcode = cache.nlsolve_alg_cache.retcode
        new_current = cache.switching_alg_cache(NonlinearSolve.get_fu(cache.nlsolve_alg_cache))
        @show new_current, cache.stats.nsteps
        cache.force_stop = cache.nlsolve_alg_cache.force_stop
        cache.u = NonlinearSolve.get_u(cache.nlsolve_alg_cache)
        cache.fu1 = NonlinearSolve.get_fu(cache.nlsolve_alg_cache)
    elseif cache.current == 2
        perform_step!(cache.ode_alg_cache)
        cache.retcode = cache.ode_alg_cache.retcode
        cache.force_stop = cache.ode_alg_cache.force_stop
        new_current = cache.switching_alg_cache(NonlinearSolve.get_fu(cache.ode_alg_cache))
        @show new_current, cache.stats.nsteps
        @show cache.switching_alg_cache.alpha
        cache.u = NonlinearSolve.get_u(cache.ode_alg_cache)
        cache.fu1 = NonlinearSolve.get_fu(cache.ode_alg_cache)
    else
        error("Unreachable Reached!")
    end

    if new_current != cache.current
        if cache.current == 1
            NonlinearSolve.set_u!(cache.ode_alg_cache, cache.u)
            NonlinearSolve.set_fu!(cache.ode_alg_cache, cache.fu1)
        elseif cache.current == 2
            NonlinearSolve.set_u!(cache.nlsolve_alg_cache, cache.u)
            NonlinearSolve.set_fu!(cache.nlsolve_alg_cache, cache.fu1)
        else
            error("Unreachable Reached!")
        end
    end

    cache.current = new_current

    return nothing
end
