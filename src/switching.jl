@concrete struct PseudoTransientSwitching
    alpha_threshold
    alpha_initial
end

function PseudoTransientSwitching(; alpha_initial = 0.001, alpha_threshold = 1.0)
    return PseudoTransientSwitching(alpha_threshold, alpha_initial)
end

function init_switching_algorithm(alg::PseudoTransientSwitching, fu;
    internalnorm = NonlinearSolve.DEFAULT_NORM, kwargs...)
    return PseudoTransientSwitchingCache(internalnorm(fu), alg.alpha_initial,
        alg.alpha_threshold, internalnorm)
end

@concrete mutable struct PseudoTransientSwitchingCache
    res_norm
    alpha
    alpha_threshold
    internalnorm
end

function (cache::PseudoTransientSwitchingCache)(fu)
    new_norm = cache.internalnorm(fu)
    @show new_norm
    cache.alpha *= cache.res_norm / new_norm
    cache.res_norm = new_norm
    # 1 --> Discrete Solver, 2 --> Continuous Solver
    return ifelse(cache.alpha > cache.alpha_threshold, 1, 2)
end
