"""
    EdgeRepresentation

Representation of a single edge in a (quantum) network.
"""
abstract type EdgeRepresentation end

"""
    NodeRepresentation

Representation of a single node in a (quantum) network.
"""
abstract type NodeRepresentation end

"""
    NetworkRepresentation

Representation of a (quantum) network.
"""
abstract type NetworkRepresentation end

"""
    EdgePhysicalRepresentation

Representation of the physical properties of a network edge.
"""
abstract type EdgePhysicalRepresentation <: EdgeRepresentation end

"""
    NodePhysicalRepresentation

Representation of the physical properties of a network node.
"""
abstract type NodePhysicalRepresentation <: NodeRepresentation end

"""
    NetworkPhysicalRepresentation

Representation of the physical properties of a (quantum) network.
"""
abstract type NetworkPhysicalRepresentation <: NetworkRepresentation end

"""
    edges(x::NetworkRepresentation)

Return iterator over all the edges (`EdgeRepresentation`) in the network.
"""
function edges end

"""
    nodes(x::NetworkRepresentation)

Return iterator over all nodes (`NodeRepresentation`) in the network.
"""
function nodes end

"""
    vertices(x::NetworkRepresentation)

Return `nodes(x)`, i.e., iterator over all the nodes (`NodeRepresentation`) in the network.
"""
vertices(x::NetworkRepresentation) = nodes(x)

"""
    EdgeAndNodes

Tuple of an edge and the nodes it connects.

# Examples
```jldoctest
julia> struct mock_edge <: EdgeRepresentation end

julia> struct mock_node <: NodeRepresentation end

julia> edge_and_nodes = (mock_edge(), (mock_node(), mock_node()))
(mock_edge(), (mock_node(), mock_node()))

julia> edge_and_nodes isa EdgeAndNodes
true

julia> edge = edge_and_nodes[1]
mock_edge()

julia> nodes = edge_and_nodes[2]
(mock_node(), mock_node())

julia> node1 = nodes[1]
mock_node()

julia> node2 = nodes[2]
mock_node()
```
"""
EdgeAndNodes = Tuple{EdgeRepresentation, Tuple{Vararg{NodeRepresentation}}}

"""
    edges_and_nodes(x::NetworkRepresentation)

Return iterator over `EdgeAndNodes`s of edges and the nodes each connects in the network.
"""
function edges_and_nodes end

"""
    num_edges(x::NetworkRepresentation)

Number of edges in the network.
"""
num_edges(x::NetworkRepresentation) = length(edges(x))

"""
    num_nodes(x::NetworkRepresentation)

Number of nodes in the network.
"""
num_nodes(x::NetworkRepresentation) = length(nodes(x))

"""
    edge_length(x::EdgeRepresentation)
    edge_length(x::NetworkRepresentation)

Length associated with an edge, or the sum of the lengths of all the edges in the network.
"""
function edge_length end
edge_length(x::NetworkRepresentation) = sum(edge_length(e) for e in edges(x))

"""
    physical_representation(x::NetworkRepresentation)
    physical_representation(x::EdgeRepresentation)
    physical_representation(x::NodeRepresentation)

Return a purely physical representation of `x`.
"""
function physical_representation end
physical_representation(x::NetworkPhysicalRepresentation) = x
physical_representation(x::EdgePhysicalRepresentation) = x
physical_representation(x::NodePhysicalRepresentation) = x

# TODO how are interfaces implemented for network representation?
# Are edges the "elements"?
# length(x::QuantumNetworkRecipe) = num_edges(x)


"""
    ChainPhysicalRepresentation

Representation of the physical properties of a repeater chain.
"""
struct ChainPhysicalRepresentation{E<:EdgePhysicalRepresentation,
    N<:NodePhysicalRepresentation} <: NetworkPhysicalRepresentation
    edges::Vector{E}
    nodes::Vector{N}
    function ChainPhysicalRepresentation(edges, nodes)
        length(nodes) == length(edges) + 1 ||
            throw(DimensionMismatch("There should be one more node than edges in a chain."))
        new{eltype(edges), eltype(nodes)}(edges, nodes)
    end
end
edges(x::ChainPhysicalRepresentation) = x.edges
nodes(x::ChainPhysicalRepresentation) = x.nodes
edges_and_nodes(x::ChainPhysicalRepresentation) = [(e, (n1, n2))
    for (e, n1, n2) in zip(edges(x), nodes(x)[1:end-1], nodes(x)[2:end])]

"""
    EdgeRecipe

A recipe for the operation of a single edge in a network.

Must specify the physical characteristics of the edge and may contain additional information
about the method used to operate it.
Typically used as a building block for `NetworkRecipe`s.
In that case, the `EdgeRecipe` may not need to provide a complete discription of what
happens at the edge, as long as the `NetworkRecipe` does.
"""
abstract type EdgeRecipe{T<:EdgePhysicalRepresentation} <: EdgeRepresentation end

"""
    edge_length(x::EdgeRecipe)

Length associated with an edge.

If this function has no method for a particular `EdgeRecipe`, it returns the `edge_length`
of its `physical_representation`.
"""
edge_length(x::EdgeRecipe) = edge_length(physical_representation(x))

"""
    NodeRecipe

A recipe for the operation of a single node in a network.

Must specify the physical characteristics of the node and may contain additional information
about the method used to operate it.
Typically used as a building block for `NetworkRecipe`s.
In that case, the `NodeRecipe` may not need to provide a complete discription of what
happens at the node, as long as the `NetworkRecipe` does.
"""
abstract type NodeRecipe{T<:NodePhysicalRepresentation} <: NodeRepresentation end

"""
    NetworkRecipe

"Recipe" for the operation of a network.

Must specify both the physical characteristics of the network and the protocol used to
operate it. A recipe must be complete enough to allow for assessment of the network
through, e.g., simulations.
"""
abstract type NetworkRecipe <: NetworkRepresentation end

"""
    ChainProtocol

Abstract type for repeater-chain protocols.
"""
abstract type ChainProtocol end

"""
    ChainRecipe(edges, nodes, chain_protocol)

Recipe for the operation of a chain of nodes (i.e., a linear network).
"""
struct ChainRecipe{E<:EdgeRecipe,N<:NodeRecipe,P<:ChainProtocol} <: NetworkRecipe
    """
    Vector of `EdgeRecipe`s. The ith element of the vector describes the ith edge in the
    chain.
    """
    edges::Vector{E}
    """
    Vector of `NodeRecipe`s. The ith element of the vector describes the ith node in the
    chain.
    """
    nodes::Vector{N}
    """
    The `ChainProtocol` used to operate the chain.
    """
    chain_protocol::P
    function ChainRecipe(edges, nodes, chain_protocol)
        length(nodes) == length(edges) + 1 ||
            throw(DimensionMismatch("There should be one more node than edges in a chain."))
        new{eltype(edges), eltype(nodes), typeof(chain_protocol)}(
            edges, nodes, chain_protocol)
    end
end

edges(x::ChainRecipe) = x.edges
nodes(x::ChainRecipe) = x.nodes
edges_and_nodes(x::ChainRecipe) = [(e, (n1, n2))
    for (e, n1, n2) in zip(edges(x), nodes(x)[1:end-1], nodes(x)[2:end])]
physical_representation(x::ChainRecipe) = ChainPhysicalRepresentation(
    physical_representation.(edges(x)), physical_representation.(nodes(x))
)
Base.convert(::Type{ChainPhysicalRepresentation}, x::ChainRecipe) =
    physical_representation(x)