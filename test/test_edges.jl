@testset "EntangledStateType" begin
    @test isabstracttype(EntangledStateType)
    @test typeof(WernerState()) <: EntangledStateType
    @test typeof(RState()) <: EntangledStateType
end

@testset "WernerState" begin
    @testset "werner_parameter_from_fidelity" begin
        f = QuantumNetworkRecipes.werner_parameter_from_fidelity
        @test f(0.75) ≈ 0.6666666666666666
        @test f(0.25) ≈ 0.0
        @test f(1.0) ≈ 1.0
    end

    @testset "fidelity_from_werner_parameter" begin
        f = QuantumNetworkRecipes.fidelity_from_werner_parameter
        @test f(0.6666666666666666) ≈ 0.75
        @test f(0.0) ≈ 0.25
        @test f(1.0) ≈ 1.0
    end

    @testset "Fidelity to Werner parameter and back again." begin
        f = QuantumNetworkRecipes.werner_parameter_from_fidelity
        g = QuantumNetworkRecipes.fidelity_from_werner_parameter
        # compositions should give identity function
        a = g ∘ f
        b = f ∘ g
        for fidelity in [0.3, 0.5, 0.75, 0.9]
            @test a(fidelity) ≈ fidelity
        end
        for werner_parameter in [0.1, 0.3, 0.5, 0.75, 0.9]
            @test b(werner_parameter) ≈ werner_parameter
        end
    end
end

@testset "heralded_entanglement.jl" begin
    include("test_heralded_entanglement.jl")
end

@testset "simple heralded entanglement" begin
    physical = SimpleHeraldedPhysical(10.0, 0.9, 15.0, 0.95)
    test_edge_physical_representation(physical)
    @test physical.edge_length == 10.0
    @test edge_length(physical) == 10.0
    @test physical.success_probability == 0.9
    @test physical.attempt_duration == 15.0
    @test physical.fidelity == 0.95

    recipe = SimpleHeraldedRecipe{WernerState}(physical)
    test_heralded_entanglement_recipe(recipe)
    @test physical_representation(recipe) === physical
    @test edge_length(recipe) == 10.0

    for phys_node in [PerfectNode(), NodeWithMemory{DepolarizingMemory}(10),
        NodeWithMemory{PauliMemory{PauliZ}}(20)]
        edge_and_nodes = (recipe, (SimpleNode(phys_node), SimpleNode(phys_node)))
        for x in [recipe, edge_and_nodes]
            @test success_probability(x) == 0.9
            @test attempt_duration(x) == 15.0
            @test entangled_state_fidelity(x) == 0.95
        end
    end
end

@testset "fiber.jl" begin
    include("test_fiber.jl")
end

@testset "DirectTransmissionOverFiber" begin
    phys = FiberSegment(100.0, 0.0, 0.9, 3)
    recipe = DirectTransmissionOverFiber{WernerState}(phys, 0.95)
    test_heralded_entanglement_recipe(recipe)
    @test physical_representation(recipe) === phys
    @test edge_length(recipe) == 100.0
    @test success_probability(recipe) ≈ 0.999
    @test attempt_duration(recipe) ≈ 2 * 100.0 /
        QuantumNetworkRecipes.speed_of_light_in_fiber
    @test entangled_state_fidelity(recipe) == 0.95

    sending_node_phys = NodeWithPhotonSource{PauliMemory{PauliZ}}(5, .6, 10, 100)
    # node with largest efficiency is picked as the sender
    sending_node = SimpleNode(sending_node_phys)
    receiving_node_phys = NodeWithPhotonSource{PauliMemory{PauliZ}}(5, .1, 3, 5)
    receiving_node = SimpleNode(receiving_node_phys)
    phys = FiberSegment(1E-5, 0.0, 0.9, 10)
    recipe = DirectTransmissionOverFiber{WernerState}(phys, 0.95)
    nodes = (sending_node, receiving_node)
    # num multiplexing modes is only limited by link and sending node
    @test num_multiplexing_modes(recipe, nodes) == 10
    succ_prob_one_mode = 0.9 * 0.6  # link and sending node
    @test success_probability(recipe, nodes) ≈ 1 - (1 - succ_prob_one_mode) ^ 10
    # link is very short so only sender node limits cycle time
    @test attempt_duration(recipe, nodes) == 10
end

include("test_single_click_heralded_entanglement.jl")