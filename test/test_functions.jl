function test_edge_physical_representation(x::EdgePhysicalRepresentation)
    @testset "test if $(typeof(x)) is a well-defined `EdgePhysicalRepresentation`" begin
        @test edge_length(x) isa Real
        @test physical_representation(x) === x
    end
end

function test_node_physical_representation(x::NodePhysicalRepresentation)
    @testset "test if $(typeof(x)) is a well-defined `NodePhysicalRepresentation`" begin
        @test physical_representation(x) === x
    end
end

function test_network_physical_representation(x::NetworkPhysicalRepresentation)
    @testset "test if $(typeof(x)) is a well-defined `NetworkPhysicalRepresentation`" begin
        @test num_edges(x) isa Integer
        @test num_nodes(x) isa Integer
        for edge in edges(x)
            @test edge isa EdgePhysicalRepresentation
        end
        for node in nodes(x)
            @test node isa NodePhysicalRepresentation
        end
        for edge_and_nodes in edges_and_nodes(x)
            @test edge_and_nodes isa EdgeAndNodes
        end
        @test edge_length(x) isa Real
        @test physical_representation(x) === x
    end
end

function test_edge_recipe(x::EdgeRecipe)
    @testset "test if $(typeof(x)) is a well-defined `EdgeRecipe`" begin
        @test edge_length(x) isa Real
        @test physical_representation(x) isa EdgePhysicalRepresentation
        test_edge_physical_representation(physical_representation(x))
    end
end

function test_node_recipe(x::NodeRecipe)
    @testset "test if $(typeof(x)) is a well-defined `NodeRecipe`" begin
        @test physical_representation(x) isa NodePhysicalRepresentation
        test_node_physical_representation(physical_representation(x))
    end
end

function test_network_recipe(x::NetworkRecipe)
    @testset "test if $(typeof(x)) is a well-defined `NetworkRecipe`" begin
        @test num_edges(x) isa Integer
        @test num_nodes(x) isa Integer
        for edge in edges(x)
            @test edge isa EdgeRecipe
        end
        for node in nodes(x)
            @test node isa NodeRecipe
        end
        for edge_and_nodes in edges_and_nodes(x)
            @test edge_and_nodes isa EdgeAndNodes
        end
        @test edge_length(x) isa Real
        @test physical_representation(x) isa NetworkPhysicalRepresentation
        test_network_physical_representation(physical_representation(x))
    end
end

function test_node_with_memory(x::AbstractNodeWithMemory)
    @testset "test if $(typeof(x)) is a well-defined `AbstractNodeWithMemory`" begin
        @test coherence_time(x) isa Real
        test_node_physical_representation(x)
    end
end

function test_node_with_photon_source(x::AbstractNodeWithPhotonSource)
    @testset "test if $(typeof(x)) is a well-defined `AbstractNodeWithPhotonSource`" begin
        @test emission_efficiency(x) isa Real
        @test cycle_time(x) isa Real
        @test num_multiplexing_modes(x) isa Integer
        test_node_with_memory(x)
    end
end

function test_heralded_entanglement_recipe(x::HeraldedEntanglement)
    @testset "test if $(typeof(x)) is a well-defined `HeraldedEntanglement`" begin
        test_edge_recipe(x)
        @test 0 ≤ success_probability(x) ≤ 1
        @test attempt_duration(x) isa Real
        @test 0 ≤ entangled_state_fidelity(x) ≤ 1
    end
end

function test_fiber_based_edge(x::FiberBasedEdge)
    @testset "test if $(typeof(x)) is a well-defined `FiberBasedEdge`" begin
        @test num_multiplexing_modes(x) isa Integer
        test_edge_physical_representation(x)
    end
end

function stoch_trip(val::Real, inf_pert::Real, fin_pert::Real, prob::Real)
    Δs = StochasticAD.similar_new(StochasticAD.create_Δs(PrunedFIsBackend(), Int),
        fin_pert, prob)
    StochasticAD.StochasticTriple{0}(val, inf_pert, Δs)
end

"""
    test_perturbation(Δ_weight_tuples)

Test if a perturbation has the right values.

The format returned by `StochasticAD.perturbations()` depends on the back end;
sometimes they are tuples, sometimes they are named tuples.
This function tests the values correctly independently of the back end.
"""
function perturbation_is(x, Δ_weight_tuples)
    for (pert, tuple) in Iterators.zip(perturbations(x), Δ_weight_tuples)
        perturbation_is_inner(pert, tuple) || return false
    end
    true
end
perturbation_is(x, Δ_weight_tuples::Tuple{N, M}) where {N<:Number, M<:Number} =
    perturbation_is(x, (Δ_weight_tuples,))

function perturbation_is_inner(pert::Tuple, tuple)
    # if there is no perturbation, we don't check if the epsilon is correct,
    # it is meaningless in that case
    pert[1] ≈ tuple[1] || return false
    iszero(pert[1]) && return true
    pert[2] ≈ tuple[2]
end
function perturbation_is_inner(pert::NamedTuple, tuple)
    pert.Δ ≈ tuple[1] || return false
    iszero(pert.Δ) && return true
    pert.weight ≈ tuple[2]
end

x = stoch_trip(1., 2., 3., 4.)
@test value(x) == 1.
@test perturbation_is(x, (3., 4.))
@test perturbation_is(x, ((3., 4.),))