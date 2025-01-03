"""
    CutoffProtocol

Protocol that specifies a cutoff policy for entanglement generation.

Cutoff policies can be used to discard entangled states that have been stored in quantum
memory for too long and have become too noisy. Usually, cutoff policies represent a tradeoff
between the fidelity of created entangled states and the rate at which they can be created.
"""
abstract type CutoffProtocol end

"""
    LocalCutoffProtocol

Quantum states are discarded based on how long they have been stored in a single node.

In this cutoff protocol, individual nodes decide whether they should keep or discard their
entangled qubits based solely on local information, namely based on how long the qubits have
been stored in its quantum memory. The protocol is specified in terms of a cutoff time,
which is the time after which a qubit is discarded. (Qubits that are exactly one cutoff
time old are not discarded; only qubits that are older than the cutoff time.)

Note that while decisions are made locally, it is important that nodes inform other nodes
in the network when they discard a qubit. Typically it is desirable to discard any qubits
that were entangled to the qubit that triggered a cutoff as soon as possible in order to
free up resources.

It is assumed that the end nodes of a repeater chain never discard their qubits.
Firstly, if they do discard qubits, this could result in a situation where one end node
believes it has an end-to-end entangled state while the other end node has discarded it
(because the associated classical message has not yet arrived).
Secondly, assuming that end nodes store their qubits until end-to-end entanglement is
confirmed, the nodes will typically store their qubits much longer than the repeaters,
which may trigger frequent undesired cutoffs. On the other hand, if the end nodes directly
measure their qubits, they do not need to store their qubits at all and there is no need for
a cutoff protocol.
"""
struct LocalCutoffProtocol <: CutoffProtocol
    """
    Time [s] after which nodes decide to discard entangled qubits.
    If the nodes in the network are `AbstractNodeWithCutoff`s, each node individually uses
    the cutoff time specified in its recipe. Each node that is not a `NodeWithCutoff` uses
    the value set here instead.
    A value of `Inf` (default) indicates that no cutoff time is used.
    """
    cutoff_time::Real
    LocalCutoffProtocol(cutoff_time::Real=Inf) = new(cutoff_time)
end

"""
    AbstractNodeWithCutoff

Abstract type for nodes that run a `LocalCutoffProtocol` with an individual cutoff time.

Implementations of this type must also implement the function `cutoff_time()`.
"""
abstract type AbstractNodeWithCutoff{T<:NodePhysicalRepresentation} <: NodeRecipe{T} end

"""
    cutoff_time(x::AbstractNodeWithCutoff)

Return the cutoff time in seconds to be used in the `LocalCutoffProtocol` for this node.
"""
function cutoff_time end

"""
    NodeWithCutoff{T}(physical_node::T, cutoff_time::Float64=Inf) where
        T<:NodePhysicalRepresentation

Node recipe that specifies a cutoff time for `LocalCutoffProtocol`.
"""
struct NodeWithCutoff{T<:NodePhysicalRepresentation} <: AbstractNodeWithCutoff{T}
    """
    Cutoff time [s] used by this node when running `LocalCutoffProtocol`.
    A value of `Inf` (default) indicates that no cutoff time is used.
    """
    physical_node::T
    cutoff_time::Real
    NodeWithCutoff(physical_node::T, cutoff_time::Float64=Inf) where
        T<:NodePhysicalRepresentation = new{T}(physical_node, cutoff_time)
end

cutoff_time(x::NodeWithCutoff) = x.cutoff_time
physical_representation(x::NodeWithCutoff) = x.physical_node
Base.convert(::Type{SimpleNode}, x::NodeWithCutoff{T}) where T =
    SimpleNode(physical_representation(x))
success_probability(x::HeraldedEntanglement, nodes::Tuple{NodeWithCutoff, NodeWithCutoff}) =
    success_probability(x, convert.(SimpleNode, nodes))
attempt_duration(x::HeraldedEntanglement, nodes::Tuple{NodeWithCutoff, NodeWithCutoff}) =
    attempt_duration(x, convert.(SimpleNode, nodes))
entangled_state_fidelity(x::HeraldedEntanglement, nodes::Tuple{NodeWithCutoff,
    NodeWithCutoff}) = entangled_state_fidelity(x, convert.(SimpleNode, nodes))