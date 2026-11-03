################################################################################
# HarmonicPowerModels.jl                                                       #
# Robust Harmonic Hosting Capacity (RHHC): thread helpers.                     #
#                                                                               #
# The robust method spends its time in independent evaluations of the         #
# moments at parameter points. These helpers run such evaluations as tasks on #
# the threads Julia was started with (`julia -t auto`) and return the results  #
# in input order, so every result is independent of the scheduling. The       #
# evaluations only read the shared network and PCE data; their only shared    #
# writes (the admittance cache, under a lock, and the dictionary, after the   #
# parallel part, in a fixed order) do not depend on the scheduling.          #
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
