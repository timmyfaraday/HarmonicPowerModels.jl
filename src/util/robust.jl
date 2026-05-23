using JuMP

# Multiplies interval by `-1`
function _flip(i::MOI.Interval)
    return MOI.Interval(-i.upper, -i.lower)
end

# Returns `τ::JuMP.AffExpr` constrained by
# `τ ≤ min_{ci in [cil,ciu], d in [dl,du]} c' x + d`
# by rewriting as
# `τ ≤ -max_{ci in [cil,ciu], d in [dl,du]} (-c)' x - d`
# or equivalently
# `τ ≤ -max_{ci in [-ciu,-dil], d in [-du,-dl]} c' x + d`
function min_box(model, c, x, d)
    return -max_box(model, _flip.(c), x, _flip(d))
end

# Returns `τ::JuMP.AffExpr` constrained by
# `τ ≥ max_{ci in [cil,ciu], d in [dl,du]} c' x + d`
# by rewriting it into
# `τ ≥ du + max_{bi in [-1, 1]} sum xi * (cil+diu)/2 + bi * xi * (dil-diu)/2`
# or equivalently to
# `τ ≥ du + sum xi * (cil+diu)/2 + |xi * (dil-diu)/2|`
# or
# `τ ≥ du + sum xi * (cil+diu)/2 + ||(xi * (dil-diu)/2)_i||_1`
# where `||⋅||_1` is the norm one cone.
# In case `cil = ciu` we get `xi * 0` so we can remove it from the
# norm one cone.
function max_box(model, c, x, d)
    τ = zero(JuMP.AffExpr)
    JuMP.add_to_expression!(τ, d.upper)
    # We don't create it in case all interval
    # have same lower and upper bounds
    norm_one = nothing
    for i in eachindex(x)
        center = (c[i].lower + c[i].upper) / 2
        JuMP.add_to_expression!(τ, center, x[i])
        if c[i].lower != c[i].upper
            if isnothing(norm_one)
                norm_one = JuMP.AffExpr[@variable(model)]
                JuMP.add_to_expression!(τ, norm_one[])
            end
            # Shift and scale the interval to `[-1, 1]`
            radius = (c[i].upper - c[i].lower) / 2
            push!(norm_one, radius * x[i])
        end
    end
    @constraint(model, norm_one in MOI.NormOneCone(length(norm_one)))
    return τ
end
