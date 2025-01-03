using Test, StructEquality
using QuantumNetworkRecipes
using Statistics, Distributions, StochasticAD
import StochasticAD.value

include("test_functions.jl")
@testset "Representations of networks, edges and nodes." begin
    include("test_network_representations.jl")
end
@testset "Performance-metric base functions." begin
    include("test_performance_metrics.jl")
end
@testset "Parameterizing network nodes." begin
    include("test_nodes.jl")
end
@testset "Parameterizing network edges." begin
    include("test_edges.jl")
end
@testset "Cutoffs." begin
    include("test_cutoffs.jl")
end
@testset "Stochastic AD helper functions." begin
    include("test_stochastic_ad_helper_functions.jl")
end
@testset "Single sequential repeater." begin
    include("test_single_sequential_repeater.jl")
end
@testset "SwapASAP" begin
    include("test_swap_asap.jl")
end