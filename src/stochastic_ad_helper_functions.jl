"""
    istrue(x)

Check if a variable is unambiguously `true` in a way that works for stochastic triples.

If any of the branches of the stochastic triple are `false`, then the result is `false`.
"""
istrue(x) = istrue(convert(Bool, x))
istrue(x::Bool) =  x
function istrue(x::StochasticAD.StochasticTriple)
    primary = isone(StochasticAD.value(x))
    perts = StochasticAD.alltrue(iszero, x.Δs)
    primary && perts
end

"""
    isfalse(x)

Check if a variable is unambiguously `false` in a way that works for stochastic triples.

If any of the branches of the stochastic triple are `true`, then the result is `false`.
"""
isfalse(x) = isfalse(convert(Bool, x))
isfalse(x::Bool) = !x
function isfalse(x::StochasticAD.StochasticTriple)
    primary = iszero(StochasticAD.value(x))
    perts = StochasticAD.alltrue(iszero, x.Δs)
    primary && perts
end

"""
    make_dual_number(primal::Float64, derivative::Float64)

Create a dual number with the given primal and derivative values.

The dual number is a `StochasticAD.StochasticTriple`, but with only an infinitesimal
perturbation.
"""
make_dual_number(primal::Float64, derivative::Float64) =
    StochasticAD.dual_number(x -> primal + derivative * x, 0)
