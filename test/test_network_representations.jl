# Test if the abstract types exist and are abstract
@testset "Abstract Types" begin
    @test isdefined(Main, :EdgeRepresentation)
    @test isdefined(Main, :NodeRepresentation)
    @test isdefined(Main, :NetworkRepresentation)
    @test isdefined(Main, :NetworkPhysicalRepresentation)
    @test isdefined(Main, :EdgePhysicalRepresentation)
    @test isdefined(Main, :NodePhysicalRepresentation)
    @test isdefined(Main, :NetworkRecipe)
    @test isdefined(Main, :EdgeRecipe)
    @test isdefined(Main, :NodeRecipe)

    @test isabstracttype(EdgeRepresentation)
    @test isabstracttype(NodeRepresentation)
    @test isabstracttype(NetworkRepresentation)
    @test isabstracttype(NetworkPhysicalRepresentation)
    @test isabstracttype(EdgePhysicalRepresentation)
    @test isabstracttype(NodePhysicalRepresentation)
    @test isabstracttype(NetworkRecipe)
    @test isabstracttype(EdgeRecipe)
    @test isabstracttype(NodeRecipe)
end

# Test if the functions exist end that they are really functions.
@testset "Functions" begin
    @test isdefined(Main, :edges)
    @test isdefined(Main, :nodes)
    @test isdefined(Main, :num_edges)
    @test isdefined(Main, :num_nodes)
    @test isdefined(Main, :edges_and_nodes)
    @test isdefined(Main, :edge_length)
    @test isdefined(Main, :physical_representation)
end

struct MockEdgePhysical <: EdgePhysicalRepresentation
    length
end

struct MockNodePhysical <: NodePhysicalRepresentation end

struct MockEdgeRecipe <: EdgeRecipe{MockEdgePhysical}
    physical_edge::MockEdgePhysical
end

struct MockNodeRecipe <: NodeRecipe{MockNodePhysical}
    physical_node::MockNodePhysical
    identifier
end

@struct_hash_equal struct SingleEdgeNetworkPhysical{E<:EdgePhysicalRepresentation,
    N<:NodePhysicalRepresentation} <: NetworkPhysicalRepresentation
    edge::E
    nodes::Vector{N}
end

@struct_hash_equal struct SingleEdgeNetwork{E<:EdgeRecipe,N<:NodeRecipe} <: NetworkRecipe
    edge::E
    nodes::Vector{N}
end

QuantumNetworkRecipes.edge_length(x::MockEdgePhysical) = x.length
QuantumNetworkRecipes.edge_length(x::MockEdgeRecipe) =
    edge_length(physical_representation(x))
QuantumNetworkRecipes.physical_representation(x::MockEdgeRecipe) = x.physical_edge
QuantumNetworkRecipes.physical_representation(x::MockNodeRecipe) = x.physical_node
QuantumNetworkRecipes.physical_representation(x::SingleEdgeNetwork) =
    SingleEdgeNetworkPhysical(physical_representation(x.edge),
        physical_representation.(x.nodes))
QuantumNetworkRecipes.edges(x::SingleEdgeNetworkPhysical) = [x.edge]
QuantumNetworkRecipes.edges(x::SingleEdgeNetwork) = [x.edge]
QuantumNetworkRecipes.nodes(x::SingleEdgeNetworkPhysical) = x.nodes
QuantumNetworkRecipes.nodes(x::SingleEdgeNetwork) = x.nodes
QuantumNetworkRecipes.edges_and_nodes(x::SingleEdgeNetworkPhysical) =
    [(x.edge, (x.nodes[1], x.nodes[2]))]
QuantumNetworkRecipes.edges_and_nodes(x::SingleEdgeNetwork) =
    [(x.edge, (x.nodes[1], x.nodes[2]))]

@testset "Functions to check if representations are well-defined." begin
    @testset "EdgePhysicalRepresentation" begin
        test_edge_physical_representation(MockEdgePhysical(1))
    end

    @testset "NodePhysicalRepresentation" begin
        test_node_physical_representation(MockNodePhysical())
    end

    @testset "NetworkPhysicalRepresentation" begin
        test_network_physical_representation(SingleEdgeNetworkPhysical(
            MockEdgePhysical(1), [MockNodePhysical(), MockNodePhysical()]
        ))
    end

    @testset "EdgeRecipe" begin
        test_edge_recipe(MockEdgeRecipe(MockEdgePhysical(1)))
    end

    @testset "NodeRecipe" begin
        test_node_recipe(MockNodeRecipe(MockNodePhysical(), 1))
    end

    @testset "NetworkRecipe" begin
        test_network_recipe(SingleEdgeNetwork(MockEdgeRecipe(MockEdgePhysical(1)),
            [MockNodeRecipe(MockNodePhysical(), 1), MockNodeRecipe(MockNodePhysical(), 2)]
        ))
    end
end

@testset "Mock representations of networks, edges and nodes." begin
    physical_edge = MockEdgePhysical(1)
    physical_nodes = [MockNodePhysical(), MockNodePhysical()]
    physical_network = SingleEdgeNetworkPhysical(physical_edge, physical_nodes)

    edge_recipe = MockEdgeRecipe(physical_edge)
    node_recipes = [MockNodeRecipe(physical_nodes[1], 1),
        MockNodeRecipe(physical_nodes[2], 2)]
    network_recipe = SingleEdgeNetwork(edge_recipe, node_recipes)

    @testset "edge_length" begin
        @test edge_length(physical_edge) == 1
        @test edge_length(edge_recipe) == 1
        @test edge_length(physical_network) == 1
    end

    @testset "physical_representation" begin
        @test physical_representation(physical_edge) == physical_edge
        @test physical_representation(physical_network) == physical_network
        @test physical_representation(edge_recipe) == physical_edge
        @test physical_representation(network_recipe) == physical_network
        for (p, r) in zip(physical_nodes, node_recipes)
            @test physical_representation(p) == p
            @test physical_representation(r) == p
        end
    end

    @testset "edges" begin
        @test edges(physical_network) == [physical_edge]
        @test edges(network_recipe) == [edge_recipe]
    end

    @testset "num_edges" begin
        @test num_edges(physical_network) == 1
        @test num_edges(network_recipe) == 1
    end

    @testset "num_nodes" begin
        @test num_nodes(physical_network) == 2
        @test num_nodes(network_recipe) == 2
    end

    @testset "nodes" begin
        @test nodes(physical_network) == physical_nodes
        @test nodes(network_recipe) == node_recipes
    end

    @testset "edges_and_nodes" begin
        physical_edges_and_nodes = edges_and_nodes(physical_network)
        @test length(physical_edges_and_nodes) == 1
        physical_edges_and_nodes, = physical_edges_and_nodes
        @test physical_edges_and_nodes isa EdgeAndNodes
        @test physical_edges_and_nodes == (physical_edge, Tuple(physical_nodes))
        recipe_edges_and_nodes = edges_and_nodes(network_recipe)
        @test length(recipe_edges_and_nodes) == 1
        recipe_edges_and_nodes, = recipe_edges_and_nodes
        @test recipe_edges_and_nodes isa EdgeAndNodes
        @test recipe_edges_and_nodes == (edge_recipe, Tuple(node_recipes))
    end
end

struct MockProtocol <: ChainProtocol end

@struct_hash_equal ChainPhysicalRepresentation
@struct_hash_equal ChainRecipe

@testset "Representations of chain networks." begin
    edges_physical = [MockEdgePhysical(1), MockEdgePhysical(2)]
    nodes_physical = [MockNodePhysical(), MockNodePhysical(), MockNodePhysical()]
    chain_physical = ChainPhysicalRepresentation(edges_physical, nodes_physical)

    @testset "ChainPhysicalRepresentation" begin
        @test chain_physical isa ChainPhysicalRepresentation{MockEdgePhysical,
            MockNodePhysical}
        @test num_edges(chain_physical) == 2
        @test num_nodes(chain_physical) == 3
        @test edges(chain_physical) == edges_physical
        @test nodes(chain_physical) == nodes_physical
        @test edges_and_nodes(chain_physical) == [
            (edges_physical[1], (nodes_physical[1], nodes_physical[2])),
            (edges_physical[2], (nodes_physical[2], nodes_physical[3]))
        ]
        @test physical_representation(chain_physical) == chain_physical
        @test edge_length(chain_physical) == 3
        @test_throws DimensionMismatch ChainPhysicalRepresentation(
            [edges_physical[1]], nodes_physical
        )
    end

    edge_recipes = [MockEdgeRecipe(physical_edge) for physical_edge in edges_physical]
    node_recipes = [MockNodeRecipe(n, i) for (i, n) in enumerate(nodes_physical)]
    protocol = MockProtocol()
    chain_recipe = ChainRecipe(edge_recipes, node_recipes, protocol)

    @testset "ChainRecipe" begin
        @test chain_recipe isa ChainRecipe{MockEdgeRecipe,MockNodeRecipe, MockProtocol}
        @test num_edges(chain_recipe) == 2
        @test num_nodes(chain_recipe) == 3
        @test edges(chain_recipe) == edge_recipes
        @test nodes(chain_recipe) == node_recipes
        @test edges_and_nodes(chain_recipe) == [
            (edge_recipes[1], (node_recipes[1], node_recipes[2])),
            (edge_recipes[2], (node_recipes[2], node_recipes[3]))
        ]
        @test physical_representation(chain_recipe) == chain_physical
        @test convert(ChainPhysicalRepresentation, chain_recipe) == chain_physical
        @test edge_length(chain_recipe) == 3
        @test_throws DimensionMismatch ChainRecipe(
            [edge_recipes[1]], node_recipes, protocol
        )
    end
end