@testset "LocalCutoffProtocol" begin
    protocol = LocalCutoffProtocol()
    @test typeof(protocol) == LocalCutoffProtocol
    @test protocol isa CutoffProtocol
    @test protocol.cutoff_time == Inf
    protocol = LocalCutoffProtocol(10.)
    @test protocol.cutoff_time == 10.
end
    
@testset "NodeWithCutoff" begin
    physical_node = PerfectNode()
    node = NodeWithCutoff(physical_node)
    @test typeof(node) == NodeWithCutoff{PerfectNode}
    @test node isa AbstractNodeWithCutoff
    @test physical_representation(node) === physical_node
    test_node_recipe(node)
    @test cutoff_time(node) == Inf
    node = NodeWithCutoff(physical_node, 10.)
    @test cutoff_time(node) == 10.
    node = NodeWithCutoff(physical_node, Inf)
    @test cutoff_time(node) == Inf
    simple_node = convert(SimpleNode, node)
    @test typeof(simple_node) == SimpleNode{PerfectNode}

    @testset "heralded entanglement" begin
        succ_prob = 0.45
        attempt_dur = 10.0
        fidelity = 0.95
        edge = SimpleHeraldedRecipe{WernerState}(
            SimpleHeraldedPhysical(10.0, succ_prob, attempt_dur, fidelity)
        )
        edge_and_nodes = (edge,
            (NodeWithCutoff(PerfectNode(), 20.), NodeWithCutoff(PerfectNode(), 30.))
        )
        @test success_probability(edge_and_nodes) == succ_prob
        @test attempt_duration(edge_and_nodes) == attempt_dur
        @test entangled_state_fidelity(edge_and_nodes) == fidelity
    end
end