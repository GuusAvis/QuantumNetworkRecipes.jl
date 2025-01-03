"""
    SimplifiedSingleClickElementaryLink

Simplified recipe for single-click heralded entanglement generation.

This can be thought of as a simplified type of `SingleClickElementaryLink`.
The simplifications include the following:
1. it is assumed that the midpoint is located exactly at the center of the fiber segment;
2. it is assumed that the bright-state parameter is the same on both sides;
3. simplified formulas are used to calculate the success probability and fidelity.

In order to reflect the fact that midpoint placement is not a degree of freedom,
this recipe is defined for a `FiberSegment` instead of a `FiberSegmentWithMidpoint`,
even though its physical representation is given as a `FiberSegmentWithMidpoint`.
It is also required that the two nodes partaking in entanglement generation have the same
emission efficiency.
"""
struct SimplifiedSingleClickElementaryLink{T<:EntangledStateType} <:
    MidHeraldedDetection{FiberSegment, T}
    physical_edge::FiberSegment
    bright_state_parameter::Real
end
physical_representation(edge::SimplifiedSingleClickElementaryLink) =
    convert(FiberSegmentWithMidpoint, edge.physical_edge)
    
function success_probability(edge::SimplifiedSingleClickElementaryLink,
        nodes::Tuple{SimpleNode{N}, SimpleNode{N}}) where
        N<:Union{PerfectNode, NodeWithMemory}
    bsp = edge.bright_state_parameter
    eff = sqrt(edge_efficiency(edge.physical_edge))
    succ_prob_single_link = 2 * bsp * eff
    success_probability_with_multiplexing(succ_prob_single_link,
        num_multiplexing_modes(edge, nodes))
end
function success_probability(edge::SimplifiedSingleClickElementaryLink,
        nodes::Tuple{SimpleNode{N}, SimpleNode{N}}) where N<:AbstractNodeWithPhotonSource
    eff_nodes = emission_efficiency.(nodes)
    eff_nodes[1] ≈ eff_nodes[2] || throw(DomainError("Nodes should have same efficiency."))
    eff = sqrt(edge_efficiency(edge.physical_edge)) * eff_nodes[1]
    bsp = edge.bright_state_parameter
    succ_prob_single_link = 2 * bsp * eff
    success_probability_with_multiplexing(succ_prob_single_link,
        num_multiplexing_modes(edge, nodes))
end

entangled_state_fidelity(edge::SimplifiedSingleClickElementaryLink,
        ::Tuple{SimpleNode{N}, SimpleNode{N}}) where N<:Union{PerfectNode, NodeWithMemory,
        NodeWithPhotonSource} =
    1 - edge.bright_state_parameter

# TODO implement more-complicated models for full version of single click.
"""
    SingleClickElementaryLink

Recipe for single-click heralded entanglement generation.
"""
struct SingleClickElementaryLink{T<:EntangledStateType} <:
        MidHeraldedDetection{FiberSegmentWithMidpoint, T}
    physical_edge::FiberSegmentWithMidpoint
    bright_state_parameter_a::Real
    bright_state_parameter_b::Real
end

Base.convert(::Type{SingleClickElementaryLink{T}},
        x::SimplifiedSingleClickElementaryLink{T}) where T =
    SingleClickElementaryLink{T}(
        physical_representation(x),
        x.bright_state_parameter,
        x.bright_state_parameter
    )