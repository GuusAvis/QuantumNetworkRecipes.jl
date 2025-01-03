edge_efficiency = QuantumNetworkRecipes.edge_efficiency

# Test fiber_attenuation_efficiency function
@testset "fiber_attenuation_efficiency" begin
    f = QuantumNetworkRecipes.fiber_attenuation_efficiency
    @test f(0) ≈ 1.0
    @test f(10, 1.) ≈ 0.1
    @test f(40, 0.5) ≈ 0.01
end

@testset "success_probability_with_multiplexing" begin
    f = QuantumNetworkRecipes.success_probability_with_multiplexing
    @test f(1., 1) ≈ 1.
    @test f(1., 100) ≈ 1.
    @test f(0., 1) ≈ 0.
    @test f(0., 1000) ≈ 0.
    @test f(0.4, 2) > f(.4, 1)
    @test f(0.9, 3) ≈ 0.999
end

# Test FiberSegment struct
@testset "FiberSegment" begin
    fiber_segment = FiberSegment(100, 0.1, 0.2)
    test_fiber_based_edge(fiber_segment)
    @test fiber_segment.link_length == edge_length(fiber_segment) == 100
    @test fiber_segment.attenuation_coefficient == 0.1
    @test fiber_segment.efficiency == 0.2
    @test edge_efficiency(fiber_segment) ≈ 0.02  # 0.1 fiber attenuation, 0.2 efficiency
    fiber_segment = FiberSegment(100, 0., 0.9, 3)
    # multiplexing is not included in `edge_efficiency`
    @test edge_efficiency(fiber_segment) ≈ 0.9
end

# Test FiberSegmentWithMidpoint struct
@testset "FiberSegmentWithMidpoint" begin
    fiber_segment_with_midpoint = FiberSegmentWithMidpoint(5, 5, 0.2, 0.3, 0.8, 0.7, 2)
    test_fiber_based_edge(fiber_segment_with_midpoint)
    @test fiber_segment_with_midpoint.length_a == 5
    @test fiber_segment_with_midpoint.length_b == 5
    @test fiber_segment_with_midpoint.attenuation_coefficient_a == 0.2
    @test fiber_segment_with_midpoint.attenuation_coefficient_b == 0.3
    @test fiber_segment_with_midpoint.efficiency_a == 0.8
    @test fiber_segment_with_midpoint.efficiency_b == 0.7
    @test fiber_segment_with_midpoint.num_multiplexing_modes == 2
    fiber_segment_with_midpoint = FiberSegmentWithMidpoint(10, 100, 1., 0.1, 0.5, 0.5, 1)
    @test edge_length(fiber_segment_with_midpoint) == 110
    # 0.1 fiber attenuation on both sides, 0.5 efficiency on both sides
    @test edge_efficiency(fiber_segment_with_midpoint) ≈ 0.0025
    fiber_segment_with_midpoint = FiberSegmentWithMidpoint(10, 0., 0.5, 2)
    # multiplexing is not included in `edge_efficiency`
    @test edge_efficiency(fiber_segment_with_midpoint) ≈ 0.5
end

@testset "success_probability_with_multiplexing" begin
    f = QuantumNetworkRecipes.success_probability_with_multiplexing
    # with success probability of 0.9, there is only a probability of 0.1^N that all N modes
    # fail, resulting in a success probability of 1 - 0.1^N.
    @test f(0.9, 1) ≈ 0.9
    @test f(0.9, 2) ≈ 0.99
    @test f(0.9, 3) ≈ 0.999
end

@testset "convert fiber segments" begin
    fiber_segment = FiberSegment(10, 0.2, 0.25, 10)
    fiber_segment_with_midpoint = convert(FiberSegmentWithMidpoint, fiber_segment)
    @test fiber_segment_with_midpoint.length_a == 5
    @test fiber_segment_with_midpoint.length_b == 5
    @test fiber_segment_with_midpoint.attenuation_coefficient_a == 0.2
    @test fiber_segment_with_midpoint.attenuation_coefficient_b == 0.2
    @test fiber_segment_with_midpoint.efficiency_a ≈ 0.5
    @test fiber_segment_with_midpoint.efficiency_b ≈ 0.5
    @test fiber_segment_with_midpoint.num_multiplexing_modes == 10
    @test edge_efficiency(fiber_segment_with_midpoint) ≈ edge_efficiency(fiber_segment)
    @test edge_length(fiber_segment_with_midpoint) ≈ edge_length(fiber_segment)

    fiber_segment_with_midpoint = FiberSegmentWithMidpoint(100, 10, 0.1, 1., 0.5, 0.5, 3)
    fiber_segment = convert(FiberSegment, fiber_segment_with_midpoint)
    @test fiber_segment.link_length == 110
    @test fiber_segment.attenuation_coefficient ≈ 20 / 110  # total 20 dB on 110 km
    @test fiber_segment.efficiency ≈ 0.25  # 0.5^2 = 0.25
    @test fiber_segment.num_multiplexing_modes == 3
    @test edge_efficiency(fiber_segment_with_midpoint) ≈ edge_efficiency(fiber_segment)
    @test edge_length(fiber_segment_with_midpoint) ≈ edge_length(fiber_segment)
    fiber_segment_with_midpoint = FiberSegmentWithMidpoint(100, .3, .5, 1)
    fiber_segment = convert(FiberSegment, fiber_segment_with_midpoint)
    @test fiber_segment.link_length == 100
    @test fiber_segment.attenuation_coefficient ≈ .3
    @test fiber_segment.efficiency ≈ .5
end