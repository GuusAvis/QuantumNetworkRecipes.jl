test_single_sequential_repeater_is_perfect(phys_nodes) = begin
    recipe_nodes = [SimpleNode(n) for n in phys_nodes]
    phys_edges = [SimpleHeraldedPhysical(20, 1, 100, 1),
        SimpleHeraldedPhysical(30, 1, 10, 1)]
    recipe_edges = [SimpleHeraldedRecipe{WernerState}(e) for e in phys_edges] 
    network = SingleSequentialRepeater(recipe_edges, recipe_nodes)
    test_network_recipe(network)
    phys_network = physical_representation(network)
    @test phys_network isa ChainPhysicalRepresentation
    for net in [network, phys_network]
        @test num_edges(net) == 2
        @test num_nodes(net) == 3
    end
    @test edges(network) == recipe_edges
    @test edges(phys_network) == phys_edges
    @test nodes(network) == recipe_nodes
    @test nodes(phys_network) == phys_nodes
    for method in [Analytical(), Sampling(10)]
        @test generation_duration(network, method) == (110, 0.)
        @test QuantumNetworkRecipes.werner_parameter(network, method) == (1, 0.)
        @test entangled_state_fidelity(network, method) == (1, 0.)
        @test qber_x(network, method) == (0, 0.)
        @test qber_y(network, method) == (0, 0.)
        @test qber_z(network, method) == (0, 0.)
        @test skr_bb84(network, method) == (1 / generation_duration(network, method)[1], 0.)
    end
end
@testset "PerfectNode" begin
    phys_nodes = [PerfectNode(), PerfectNode(), PerfectNode()]
    test_single_sequential_repeater_is_perfect(phys_nodes)
end

@testset "infinite coherence time" begin
    phys_nodes = [NodeWithMemory{DepolarizingMemory}(Inf) for _ in 1:3]
    test_single_sequential_repeater_is_perfect(phys_nodes)
end

@testset "lots of decoherence" begin
    phys_nodes = [NodeWithMemory{DepolarizingMemory}(1E-5) for _ in 1:3]
    recipe_nodes = [SimpleNode(n) for n in phys_nodes]
    phys_edges = [SimpleHeraldedPhysical(20, 1E-8, 100, 1),
        SimpleHeraldedPhysical(30, 1E-8, 10, 1)]
    network = SingleSequentialRepeater([SimpleHeraldedRecipe{WernerState}(e)
        for e in phys_edges], recipe_nodes)
    for method in [Analytical(), Sampling(10)]
        @test generation_duration(network, method)[1] > 1E4
        @test QuantumNetworkRecipes.werner_parameter(network, method)[1] ≈ 0
        @test entangled_state_fidelity(network, method)[1] ≈ 1 / 4
        @test qber_x(network, method)[1] ≈ 1 / 2
        @test qber_y(network, method)[1] ≈ 1 / 2
        @test qber_z(network, method)[1] ≈ 1 / 2
        @test skr_bb84(network, method, false)[1] < 0
        if method isa Analytical  # test errors in estimates
            @test generation_duration(network, method)[2] == 0.
            @test QuantumNetworkRecipes.werner_parameter(network, method)[2] == 0.
            @test entangled_state_fidelity(network, method)[2] == 0.
            @test qber_x(network, method)[2] == 0.
            @test skr_bb84(network, method)[2] == 0.
        else
            @test generation_duration(network, method)[2] > 1.
            @test QuantumNetworkRecipes.werner_parameter(network, method)[2] ≈ 0
            @test entangled_state_fidelity(network, method)[2] ≈ 0.
            @test qber_x(network, method)[2] ≈ 0.
            @test skr_bb84(network, method, false)[2] > 0.
        end
    end
end

@testset "no decoherence but bad fidelity" begin
    recipe_nodes = [SimpleNode(PerfectNode()) for _ in 1:3]
    phys_edges = [SimpleHeraldedPhysical(20, 1, 100, 1 / 4),
        SimpleHeraldedPhysical(30, 1, 10, 1 / 4)]
    network = SingleSequentialRepeater([SimpleHeraldedRecipe{WernerState}(e)
        for e in phys_edges], recipe_nodes)
    for method in [Analytical(), Sampling(10)]
        @test QuantumNetworkRecipes.werner_parameter(network, method)[1] ≈ 0
        @test QuantumNetworkRecipes.werner_parameter(network, method)[2] ≈ 0
        @test entangled_state_fidelity(network, method)[1] ≈ 1 / 4
        @test entangled_state_fidelity(network, method)[2] ≈ 0.
        @test qber_x(network, method)[1] ≈ 1 / 2
        @test qber_x(network, method)[2] ≈ 0.
        @test qber_y(network, method)[1] ≈ 1 / 2
        @test qber_y(network, method)[2] ≈ 0.
        @test qber_z(network, method)[1] ≈ 1 / 2
        @test qber_z(network, method)[2] ≈ 0.
        @test skr_bb84(network, method, false)[1] < 0
    end
end