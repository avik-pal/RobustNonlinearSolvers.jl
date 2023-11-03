@concrete struct ContinuousNonlinearSolveAlgorithm <:
                 NonlinearSolve.AbstractNonlinearSolveAlgorithm
    alg
end

function ContinuousNonlinearSolveAlgorithm(; alg = Tsit5())
    return ContinuousNonlinearSolveAlgorithm(alg)
end

function __init(prob::NonlinearProblem, alg::ContinuousNonlinearSolveAlgorithm, args...;
    reltol = nothing, abstol = nothing, termination_condition = nothing,
    maxiters::Int = 1000, kwargs...)
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
    odecache = __init(odeprob, alg.alg; kwargs..., callback, maxiters)

    return ContinuousNonlinearSolveCache{isinplace(prob)}(prob, alg, prob.f, fu1, prob.u0,
        odecache, abstol, reltol, false, maxiters, ReturnCode.Default, tc_cache,
        NLStats(1, 0, 0, 0, 0))
end

@concrete mutable struct ContinuousNonlinearSolveCache{iip} <:
                         NonlinearSolve.AbstractNonlinearSolveCache{iip}
    prob
    alg
    f
    fu1
    u
    odecache
    abstol
    reltol
    force_stop::Bool
    maxiters::Int
    retcode::ReturnCode.T
    tc_cache
    stats
end

function perform_step!(cache::ContinuousNonlinearSolveCache{iip}) where {iip}
    SciMLBase.step!(cache.odecache)
    cache.u = cache.odecache.u

    if iip
        cache.f(cache.fu1, cache.u, cache.prob.p)
    else
        cache.fu1 = cache.f(cache.u, cache.prob.p)
    end

    if SciMLBase.successful_retcode(cache.odecache.sol.retcode)
        cache.force_stop = true
        cache.retcode = ReturnCode.Success
    end
    return nothing
end

# TODO: `reinit!`
