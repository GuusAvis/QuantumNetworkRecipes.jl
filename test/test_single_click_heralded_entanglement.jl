@testset "simplified single click" begin
    edge = SimplifiedSingleClickElementaryLink{WernerState}(FiberSegment(10, 0.2, 0.1, 3),
        0.2)
    test_heralded_entanglement_recipe(edge)
    edge2 = SimplifiedSingleClickElementaryLink{WernerState}(FiberSegment(10, 0.2, 0.1, 3),
        0.3)
    @test success_probability(edge2) > success_probability(edge)
    @test entangled_state_fidelity(edge2) < entangled_state_fidelity(edge)
    edge3 = SimplifiedSingleClickElementaryLink{WernerState}(FiberSegment(10, 0.2, 0.1, 4),
        0.2)
    @test success_probability(edge3) > success_probability(edge)
    @test entangled_state_fidelity(edge3) ≈ entangled_state_fidelity(edge)
    edge4 = SimplifiedSingleClickElementaryLink{WernerState}(FiberSegment(50, 0.2, 0.1, 3),
        0.2)
    @test success_probability(edge4) < success_probability(edge)
    @test entangled_state_fidelity(edge4) ≈ entangled_state_fidelity(edge)
    nodes = Tuple(SimpleNode(NodeWithPhotonSource{DepolarizingMemory}(10, .5, 5, 1))
        for _ in 1:2)
    @test num_multiplexing_modes(edge, nodes) == 1
    @test entangled_state_fidelity(edge, nodes) ≈ entangled_state_fidelity(edge)
    edge5 = SimplifiedSingleClickElementaryLink{WernerState}(FiberSegment(10, 0.2, 0.1, 1),
        0.2)
    # single click has sqrt(total efficiency) scaling, so the extra 50% loss in the nodes
    # becomes a total 50% smaller success probability
    @test success_probability(edge5, nodes) ≈ success_probability(edge5) / 2
end