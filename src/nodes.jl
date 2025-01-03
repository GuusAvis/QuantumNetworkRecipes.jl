"""A node with perfect hardware."""
struct PerfectNode <: NodePhysicalRepresentation end

"""Parameterizes the type of a quantum memory."""
abstract type MemoryType end

"""Quantum memory that is subject to depolarizing noise."""
struct DepolarizingMemory <: MemoryType end

"""Parameterizes the type of a Pauli operator."""
abstract type PauliType end

"""Pauli operator X."""
struct PauliX <: PauliType end

"""Pauli operator Y."""
struct PauliY <: PauliType end

"""Pauli operator Z."""
struct PauliZ <: PauliType end

"""Quantum memory that is subject to Pauli noise."""
struct PauliMemory{T<:PauliType} <: MemoryType end

"""Abstract type for nodes that have quantum memory with a specified type.

Must implement the function `coherence_time(x::AbstractNodeWithMemory)`.
"""
abstract type AbstractNodeWithMemory{T<:MemoryType} <: NodePhysicalRepresentation end

"""
    coherence_time(x::AbstractNodeWithMemory)

Return the coherence time of the quantum memory of the node `x`.
"""
function coherence_time end
coherence_time(x::NodeRecipe{T}) where {T<:AbstractNodeWithMemory} =
    coherence_time(physical_representation(x))

"""A node with a quantum memory of a specific time and coherence time.

This node is taken to be perfect (like `PerfectNode`) but with a finite coherence time.
"""
struct NodeWithMemory{T<:MemoryType,R<:Real} <: AbstractNodeWithMemory{T}
    coherence_time::R
end
coherence_time(x::NodeWithMemory) = x.coherence_time
NodeWithMemory{T}(coherence_time::R) where {T,R} = NodeWithMemory{T,R}(coherence_time)

"""Abstract type for nodes that can emit photons and has a quantum memory.

Must implement the functions `emission_efficiency(x::AbstractNodeWithPhotonSource)`,
`cycle_time(x::AbstractNodeWithPhotonSource)`, and
`num_multiplexing_modes(x::AbstractNodeWithPhotonSource)`.
"""
abstract type AbstractNodeWithPhotonSource{T<:MemoryType} <: AbstractNodeWithMemory{T} end

"""
    emission_efficiency(x::AbstractNodeWithPhotonSource)

Return the probability that a photon is actually emitted and leaves in the correct mode.
"""
function emission_efficiency end
emission_efficiency(x::NodeRecipe{T}) where {T<:AbstractNodeWithPhotonSource} =
    emission_efficiency(physical_representation(x))

"""
    cycle_time(x::AbstractNodeWithPhotonSource)

Return the time required for emitting one photon.
"""
function cycle_time end
cycle_time(x::NodeRecipe{T}) where {T<:AbstractNodeWithPhotonSource} =
    cycle_time(physical_representation(x))

"""
    num_multiplexing_modes(x::AbstractNodeWithPhotonSource)

Return the number of multiplexing modes supported by this node.
"""
function num_multiplexing_modes end
num_multiplexing_modes(x::NodeRecipe{T}) where {T<:AbstractNodeWithPhotonSource} =
    num_multiplexing_modes(physical_representation(x))

"""A node that can emit photons."""
struct NodeWithPhotonSource{T<:MemoryType,R_coh<:Real,R_eff<:Real,R_time<:Real} <:
       AbstractNodeWithPhotonSource{T}
    coherence_time::R_coh
    """Probability that a photon is actually emitted and leaves in the correct mode."""
    emission_efficiency::R_eff
    """Time required for emitting one photon."""
    cycle_time::R_time
    """Number of multiplexing modes supported by this node."""
    num_multiplexing_modes::Int
end
NodeWithPhotonSource{T}(
        coherence_time::R_coh, emission_efficiency::R_eff, cycle_time::R_time,
        num_multiplexing_modes::Int) where {T,R_coh,R_eff,R_time} =
    NodeWithPhotonSource{T,R_coh,R_eff,R_time}(
        coherence_time, emission_efficiency, cycle_time, num_multiplexing_modes)
coherence_time(x::NodeWithPhotonSource) = x.coherence_time
emission_efficiency(x::NodeWithPhotonSource) = x.emission_efficiency
cycle_time(x::NodeWithPhotonSource) = x.cycle_time
num_multiplexing_modes(x::NodeWithPhotonSource) = x.num_multiplexing_modes

"""Simple node recipe that does not come pre-equiped with particular local protocols."""
struct SimpleNode{T<:NodePhysicalRepresentation} <: NodeRecipe{T}
    physical_node::T
end
physical_representation(x::SimpleNode) = x.physical_node

WernerCompatibleNode = SimpleNode{T} where {T<:Union{PerfectNode,
        NodeWithMemory{DepolarizingMemory}}}