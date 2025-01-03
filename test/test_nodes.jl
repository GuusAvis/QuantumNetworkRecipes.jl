# Test PerfectNode
@testset "PerfectNode" begin
    @test typeof(PerfectNode()) == PerfectNode
    test_node_physical_representation(PerfectNode())
end

# Test DepolarizingMemory
@testset "DepolarizingMemory" begin
    @test typeof(DepolarizingMemory()) == DepolarizingMemory
end

@testset "PauliType" begin
    @test PauliX() isa PauliType
    @test PauliY() isa PauliType
    @test PauliZ() isa PauliType
end

# Test PauliMemory
@testset "PauliMemory" begin
    @test typeof(PauliMemory{PauliX}()) == PauliMemory{PauliX}
    @test typeof(PauliMemory{PauliY}()) == PauliMemory{PauliY}
    @test typeof(PauliMemory{PauliZ}()) == PauliMemory{PauliZ}
end

# Test NodeWithMemory
@testset "NodeWithMemory" begin
    coh_time_float::Float64 = 1.0
    node = NodeWithMemory{DepolarizingMemory}(coh_time_float)
    @test typeof(node) == NodeWithMemory{DepolarizingMemory, Float64}
    @test node isa AbstractNodeWithMemory
    test_node_with_memory(node)
    coh_time_int::Int = 1
    node = NodeWithMemory{PauliMemory{PauliX}}(coh_time_int)
    @test typeof(node) == NodeWithMemory{PauliMemory{PauliX}, Int}
    @test node.coherence_time == coh_time_int
    test_node_with_memory(node)
end

@testset "NodeWithPhotonSource" begin
    coh_time = 5.0
    emission_eff = 0.9
    cycle_t = 1.0
    num_modes = 2
    node = NodeWithPhotonSource{DepolarizingMemory}(
        coh_time, emission_eff, cycle_t, num_modes)
    test_node_with_photon_source(node)
    @test coherence_time(node) == coh_time
    @test emission_efficiency(node) == emission_eff
    @test cycle_time(node) == cycle_t
    @test num_multiplexing_modes(node) == num_modes
end

# Test SimpleNode
@testset "SimpleNode" begin
    physical_node = PerfectNode()
    node = SimpleNode(physical_node)
    @test typeof(node) == SimpleNode{PerfectNode}
    @test physical_representation(node) === physical_node
    test_node_recipe(node)
end