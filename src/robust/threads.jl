################################################################################
# HarmonicPowerModels.jl                                                       #
# Robust Harmonic Hosting Capacity (RHHC): thread helpers.                     #
#                                                                               #
# The robust method spends its time in independent evaluations of the         #
# moments at parameter points. These helpers run such evaluations as tasks on #
# the threads Julia was started with (`julia -t auto`) and return the results  #
# in input order, so every result is independent of the scheduling. The       #
# evaluations only read the shared network and PCE data; their only shared    #
# writes (caches, the dictionary) go through `_cached` or happen after the    #
# parallel part, in a fixed order.                                             #
#                                                                               #
# The matrices are small (one row per bus), so BLAS threads only compete with #
# these tasks: run with `LinearAlgebra.BLAS.set_num_threads(1)`.              #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

# f applied to every element of xs, as one task per element; results in the order of xs
function _tmap(f, xs)
    Threads.nthreads() == 1 && return map(f, xs)
    tasks = [Threads.@spawn f(x) for x in xs]
    return [fetch(t) for t in tasks]
end

# f applied to every element of xs for its side effects (on disjoint data), one task each
function _tforeach(f, xs)
    if Threads.nthreads() == 1
        foreach(f, xs)
    else
        foreach(wait, [Threads.@spawn f(x) for x in xs])
    end
    return nothing
end

# value of f at x from `cache`, computed outside the lock on a miss; the cached values are
# deterministic, so a concurrent duplicate computation only costs time
function _cached(f, cache::Dict, lk::ReentrantLock, x)
    v = lock(() -> get(cache, x, nothing), lk)
    v === nothing || return v
    v = f(x)
    lock(() -> (cache[x] = v), lk)
    return v
end
