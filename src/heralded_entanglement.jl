"""
    HeraldedEntanglement

Recipe for an edge on which heralded entanglement generation takes place.

In heralded entanglement generation, a sequence of attempts at entanglement generation is
performed. Success of these attempts is "heralded", that is, the nodes that are performing
heralded entanglement generation learn whether the attempt was a succes and entanglement
was generated, or wether it was a failure. A series of attempts is performed until a success
is obtained.

It is assumed that all attempts are identical and have the same success probability.
That is, they are modeled as identical and independent Bernoulli trials.
As a consequence, the number of attempts until success is geometrically distributed.
"""
abstract type HeraldedEntanglement{T<:EdgePhysicalRepresentation, S<:EntangledStateType} <:
    EdgeRecipe{T} end

"""
    success_probability(edge::HeraldedEntanglement, nodes::Tuple{Vararg{NodeRecipe}})
    success_probability(x::EdgeAndNodes)
    success_probability(edge::HeraldedEntanglement)

The probability that a single attempt at heralded entanglement generation succeeds.
This may depend not only on the properties of the edge itself, but also on the properties
of the nodes that are connected by the edge.
If an `EdgeAndNodes` is provided, the edge and nodes are automatically unpacked.
If only an `HeraldedEntanglement` is provided, the node recipes are assumed to be a pair of
the type `SimpleNode{PerfectNode}`.
"""
function succcess_probability end
success_probability(x::EdgeAndNodes) = success_probability(x[1], x[2])
success_probability(x::HeraldedEntanglement) =
    success_probability(x, (SimpleNode(PerfectNode()), SimpleNode(PerfectNode())))

"""
    attempt_duration(edge::HeraldedEntanglement, nodes::Tuple{Vararg{NodeRecipe}})
    attempt_duration(x::EdgeAndNodes)
    attempt_duration(edge::HeraldedEntanglement)

The duration of a single attempt at heralded entanglement generation.
This may depend not only on the properties of the edge itself, but also on the properties
of the nodes that are connected by the edge.
If an `EdgeAndNodes` is provided, the edge and nodes are automatically unpacked.
If only an `HeraldedEntanglement` is provided, the node recipes are assumed to be a pair of
the type `SimpleNode{PerfectNode}`.
"""
function attempt_duration end
attempt_duration(x::EdgeAndNodes) = attempt_duration(x[1], x[2])
attempt_duration(x::HeraldedEntanglement) =
    attempt_duration(x, (SimpleNode(PerfectNode()), SimpleNode(PerfectNode())))

"""
    entangled_state_fidelity(edge::HeraldedEntanglement, nodes::Tuple{Vararg{NodeRecipe}})
    entangled_state_fidelity(x::EdgeAndNodes)
    entangled_state_fidelity(edge::HeraldedEntanglement)

The fidelity of the entangled states that are generated when a success is heralded.
This may depend not only on the properties of the edge itself, but also on the properties
of the nodes that are connected by the edge.
If an `EdgeAndNodes` is provided, the edge and nodes are automatically unpacked.
If only an `HeraldedEntanglement` is provided, the node recipes are assumed to be a pair of
the type `SimpleNode{PerfectNode}`.
"""
entangled_state_fidelity(x::EdgeAndNodes) = entangled_state_fidelity(x[1], x[2])
entangled_state_fidelity(x::HeraldedEntanglement) =
    entangled_state_fidelity(x, (SimpleNode(PerfectNode()), SimpleNode(PerfectNode())))

"""
    generation_duration(edge::HeraldedEntanglement, nodes::Tuple{Vararg{NodeRecipe}})
    generation_duration(x::EdgeAndNodes)
    generation_duration(edge::HeraldedEntanglement)

Determine the expected duration until the first success is heralded.
If an `EdgeAndNodes` is provided, the edge and nodes are automatically unpacked.
If only an `HeraldedEntanglement` is provided, the node recipes are assumed to be a pair of
the type `SimpleNode{PerfectNode}`.
"""
generation_duration(x::Union{EdgeAndNodes, HeraldedEntanglement}) = 
    attempt_duration(x) / success_probability(x)
generation_duration(edge::HeraldedEntanglement, nodes::Tuple{Vararg{NodeRecipe}}) =
    attempt_duration(edge, nodes) / success_probability(edge, nodes)

"""
    NumberOfAttempts(success_probability)
    NumberOfAttempts(x::EdgeAndNodes)

The number of identical and independent attempts needed to obtain a success.

Every attempt is an independent Bernoulli process with the same success probability.
The number of attempts until success is a random variable that is geometrically distributed.
However, it deviates slightly from `Geometric` included in `Distributions.jl`,
as it uses a different convention:
that one counts the number of _failed_ attempts until the first success.
`NumberOfAttempts` also counts the succesful attempt.

If an `EdgeAndNodes` is provided, first `success_probability` is called on the input.
"""
struct NumberOfAttempts{T<:Real} <: Distributions.DiscreteUnivariateDistribution
    success_probability::T
    ρ::Distributions.LocationScale{Int, Distributions.Discrete, Distributions.Geometric{T}}
    function NumberOfAttempts(success_probability::T) where T <: Real
        ρ = Distributions.Geometric(success_probability) + 1
        new{T}(success_probability, ρ)
    end
end
NumberOfAttempts(x::EdgeAndNodes) = NumberOfAttempts(success_probability(x))
Distributions.succprob(x::NumberOfAttempts) = x.success_probability
Distributions.failprob(x::NumberOfAttempts) = 1 - x.success_probability
Distributions.params(x::NumberOfAttempts) = (x.success_probability,)
Distributions.partype(::NumberOfAttempts{T}) where T = T
Base.show(io::IO, x::NumberOfAttempts{T}) where T = print(io,
    "NumberOfAttempts{$T}(success_probability=$(x.success_probability))")

"""
    mean_exponential_number_of_attempts(x::NumberOfAttempts, base)
    mean_exponential_number_of_attempts(success_probability, base)

Calculate the expected value of `base` raised to the power of the number of attempts
until success.
"""
function mean_exponential_number_of_attempts(x::NumberOfAttempts, base)
    mean_exponential_number_of_attempts(Distributions.succprob(x), base)
end
function mean_exponential_number_of_attempts(success_probability, base)
    denom = 1 - base * (1 - success_probability)
    nom = base * success_probability
    nom / denom
end

"""
    Duration(success_probability, attempt_duration)
    Duration(x::EdgeAndNodes)

The duration until success for trials with fixed success probability and duration.

This is a random variable that represents the amount of time required until a successful
attempt has taken place, given that each attempt has the same success probability and the
same duration.
This random variable is used to model the duration of heralded entanglement generation.
Includes the time required to finish the successful attempt.

If an `EdgeAndNodes` is provided, first the `success_probability` and `attempt_duration` are
determined from that.
"""
struct Duration{S<:Real, T<:Real} <: Distributions.DiscreteUnivariateDistribution
    success_probability::S
    attempt_duration::T
    ρ::Distributions.DiscreteUnivariateDistribution
    function Duration(success_probability::S, attempt_duration::T) where {S, T}
        ρ = NumberOfAttempts(success_probability) * attempt_duration
        new{S, T}(success_probability, attempt_duration, ρ)
    end
end
Duration(x::EdgeAndNodes) = Duration(success_probability(x), attempt_duration(x))
Distributions.params(x::Duration) = (x.success_probability, x.attempt_duration)
Distributions.partype(::Duration{S, T}) where {S, T}= (S, T)
Base.show(io::IO, x::Duration{S, T}) where {S, T} = print(io,
    "Duration{$S, $T}(success_probability=$(x.success_probability), " *
    "attempt_duration=$(x.attempt_duration))")

_unary_functions_to_extend = [
    :(Statistics.mean), :(Statistics.median), :(Distributions.mode), :(Statistics.var),
    :(Distributions.skewness), :(Distributions.kurtosis), :(Distributions.entropy),
    :(Distributions.sampler), :(Distributions.minimum), :(Distributions.maximum),
]
_binary_functions_to_extend = [
    :(Distributions.pdf), :(Distributions.logpdf), :(Distributions.cdf),
    :(Distributions.logcdf), :(Distributions.ccdf), :(Distributions.logccdf),
    :(Distributions.quantile), :(Distributions.cquantile), :(Distributions.insupport),
    :(Distributions.mgf), :(Distributions.cf),
]
for dist in [:NumberOfAttempts, :Duration]
    @eval function Base.rand(rng::Random.AbstractRNG, x::$dist)
        rand(rng, x.ρ)
    end
    for f in _unary_functions_to_extend
        @eval function $(f)(x::$dist)
            $f(x.ρ)
        end
    end
    for f in _binary_functions_to_extend
        @eval function $(f)(x::$dist, y::Real)
            $f(x.ρ, y)
        end
    end
    @eval Distributions.kldivergence(x::$dist, y) =
        Distributions.kldivergence(x.ρ, y)
    @eval Distributions.kldivergence(x, y::$dist) =
        Distributions.kldivergence(x, y.ρ)
    @eval Distributions.kldivergence(x::$dist, y::$dist) =
        Distributions.kldivergence(x.ρ, y.ρ)
end