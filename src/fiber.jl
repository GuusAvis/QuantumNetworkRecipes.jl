"""Speed of light in fiber in km/s."""
const speed_of_light_in_fiber = 200_000

"""
    fiber_attenuation_efficiency(fiber_length, attenuation_coefficient=0.2)

Probability of photon traveling through fiber without being lost due to attenuation.

Given by ``10^{-αL/10}``, where ``α`` is the attenuation coefficient [dB/km] and ``L`` is
the fiber length [km].
"""
function fiber_attenuation_efficiency(fiber_length, attenuation_coefficient=0.2)
    10^(-attenuation_coefficient * fiber_length / 10)
end

"""
    success_probability_with_multiplexing(success_probability, num_modes)

The probability of success for a multiplexed heralded entanglement generation scheme.

`success_probability` is the probability of success for a single mode, and `num_modes` is
the number of modes supported for multiplexing.
"""
function success_probability_with_multiplexing(success_probability, num_modes)
    isone(num_modes) && return success_probability
    1 - (1 - success_probability)^num_modes
end

"""
    FiberBasedEdge

An abstract type for network edges formed by optical fibers.
"""
abstract type FiberBasedEdge <: EdgePhysicalRepresentation end

"""
    num_multiplexing_modes(x::FiberBasedEdge)
    num_multiplexing_modes(x::EdgeRecipe{FiberBasedEdge})

Number of multiplexing modes supported by the edge. Defaults to 1.
"""
num_multiplexing_modes(x::FiberBasedEdge) = 1
num_multiplexing_modes(x::EdgeRecipe{T}) where T<:FiberBasedEdge =
    num_multiplexing_modes(physical_representation(x))

"""
    num_multiplexing_modes(edge, nodes)
    
The number of multiplexing modes supported over the edge that is connected by the nodes.
"""
# num_multiplexing_modes(edge::Union{E, EdgeRecipe{E}}, nodes::Tuple{Vararg{N}})  where
#         {E<:FiberBasedEdge, N<:AbstractNodeWithPhotonSource} =
#     min(num_multiplexing_modes(edge), minimum(num_multiplexing_modes.(nodes)))

"""
    FiberSegment

A network edge formed by a segment of optical fiber.
"""
struct FiberSegment <: FiberBasedEdge
    """
    Total length [km] of the fiber segment.
    """
    link_length::Float64
    """
    Attenuation coefficient [dB/km] of the optical fiber.
    """
    attenuation_coefficient::Float64
    """
    Efficiency of the fiber segment. This is the probability that a photon sent through the
    fiber is detected on the other side if there were no fiber attenuation.
    """
    efficiency::Float64
    """Number of modes supported for multiplexing."""
    num_multiplexing_modes::Integer
    function FiberSegment(link_length, attenuation_coefficient=0.2,
        efficiency=1.0, num_multiplexing_modes=1)
        link_length ≥ 0 || DomainError(link_length)
        attenuation_coefficient ≥ 0 || DomainError(attenuation_coefficient)
        0 ≤ efficiency ≤ 1 || DomainError(efficiency)
        num_multiplexing_modes ≥ 1 || throw(DomainError(num_multiplexing_modes))
        new(link_length, attenuation_coefficient, efficiency, num_multiplexing_modes)
    end
end
num_multiplexing_modes(x::FiberSegment) = x.num_multiplexing_modes
edge_length(x::FiberSegment) = x.link_length

"""
    FiberSegmentWithMidpoint

A network edge formed by two segments of optical fibers joined at a some "midpoint station".
"""
struct FiberSegmentWithMidpoint <: FiberBasedEdge
    """Fiber length [km] between "A" node and midpoint station."""
    length_a::Real
    """Fiber length [km] between "B" node and midpoint station."""
    length_b::Real
    """Attenuation coffiecient [dB/km] of fiber between "A" node and midpoint station."""
    attenuation_coefficient_a::Real
    """Attenuation coffiecient [dB/km] of fiber between "B" node and midpoint station."""
    attenuation_coefficient_b::Real
    """Efficiency of the fiber segment between "A" node and midpoint station."""
    efficiency_a::Real
    """Efficiency of the fiber segment between "B" node and midpoint station."""
    efficiency_b::Real
    """Number of modes supported for multiplexing."""
    num_multiplexing_modes::Integer
    function FiberSegmentWithMidpoint(length_a, length_b, attenuation_coefficient_a=0.2,
        attenuation_coefficient_b::Real=0.2, efficiency_a=1.0, efficiency_b=1.0,
        num_multiplexing_modes=1)
        for len in [length_a, length_b]
            len ≥ 0 || throw(DomainError(len))
        end
        for att in [attenuation_coefficient_a, attenuation_coefficient_b]
            att ≥ 0 || throw(DomainError(att))
        end
        for eff in [efficiency_a, efficiency_b]
            0 ≤ eff ≤ 1 || throw(DomainError(eff))
        end
        num_multiplexing_modes ≥ 1 || throw(DomainError(num_multiplexing_modes))
        new(length_a, length_b, attenuation_coefficient_a, attenuation_coefficient_b,
            efficiency_a, efficiency_b, num_multiplexing_modes)
    end
end
num_multiplexing_modes(x::FiberSegmentWithMidpoint) = x.num_multiplexing_modes
edge_length(x::FiberSegmentWithMidpoint) = x.length_a + x.length_b

FiberSegmentWithMidpoint(link_length::Real, attenuation_coefficient::Real=0.2,
    efficiency::Real=1.0, num_multiplexing_modes::Integer=1) =
    FiberSegmentWithMidpoint(
        link_length / 2, link_length / 2, attenuation_coefficient, attenuation_coefficient,
        sqrt(efficiency), sqrt(efficiency), num_multiplexing_modes
    )

Base.convert(::Type{FiberSegment}, x::FiberSegmentWithMidpoint) = FiberSegment(
    edge_length(x),
    (x.attenuation_coefficient_a * x.length_a + x.attenuation_coefficient_b * x.length_b) /
    edge_length(x),
    x.efficiency_a * x.efficiency_b,
    x.num_multiplexing_modes
)
Base.convert(::Type{FiberSegmentWithMidpoint}, x::FiberSegment) =
    FiberSegmentWithMidpoint(edge_length(x), x.attenuation_coefficient, x.efficiency,
    x.num_multiplexing_modes)

# TODO not sure if this function really adds much, perhaps should be removed
"""
edge_efficiency(x)

The total efficiency associated with a network edge.
"""
edge_efficiency(x::FiberSegment) =
x.efficiency * fiber_attenuation_efficiency(edge_length(x), x.attenuation_coefficient)
edge_efficiency(x::FiberSegmentWithMidpoint) = 
    x.efficiency_a * x.efficiency_b *
        fiber_attenuation_efficiency(x.length_a, x.attenuation_coefficient_a) *
        fiber_attenuation_efficiency(x.length_b, x.attenuation_coefficient_b)
edge_efficiency(x::EdgeRecipe{N}) where N<:FiberBasedEdge =
    edge_efficiency(physical_representation(x))