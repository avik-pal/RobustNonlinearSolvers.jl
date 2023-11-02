@concrete struct ContinuousNonlinearSolveAlgorithm <:
                 NonlinearSolve.AbstractNonlinearSolveAlgorithm
    alg
end

function __init(prob::NonlinearProblem, alg::ContinuousNonlinearSolveAlgorithm, args...;
        reltol = nothing, abstol = nothing, termination_condition = nothing, kwargs...)
    T = eltype(prob.u0)
    f = if isinplace(prob)
        (du, u, p, _) -> prob.f(du, u, p)
    else
        (u, p, _) -> prob.f(u, p)
    end

    fu1 = NonlinearSolve.evaluate_f(prob, prob.u0)

    abstol, reltol, tc_cache = NonlinearSolve.init_termination_cache(abstol, reltol, fu1,
        prob.u0, termination_condition)

    callback = TerminateSteadyState(abstol, reltol, tc_cache)

    if haskey(kwargs, :callback)
        callback = CallbackSet(callback, kwargs[:callback])
    end

    odeprob = ODEProblem{isinplace(prob)}(f, prob.u0, (T(0), T(Inf)), prob.p)
    odecache = __init(odeprob, alg.alg; kwargs..., callback)

    return ContinuousNonlinearSolveCache(prob.f, odecache, abstol, reltol, tc_cache,
        NLStats(1, 0, 0, 0, 0))
end

@concrete mutable struct ContinuousNonlinearSolveCache <:
                         NonlinearSolve.AbstractNonlinearSolveCache
    f
    odecache
    abstol
    reltol
    tc_cache
    stats
end

function perform_step!(cache::ContinuousNonlinearSolveCache)
    SciMLBase.step!(cache.odecache)
    return nothing
end
