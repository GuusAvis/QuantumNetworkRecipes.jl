module QuantumNetworkRecipes

import Distributions
using Random, Statistics, StochasticAD
import Graphs.edges, Graphs.vertices

export
    # network_representations.jl
    EdgeRepresentation, NodeRepresentation, NetworkRepresentation,
    EdgePhysicalRepresentation, NodePhysicalRepresentation, NetworkPhysicalRepresentation,
    EdgeRecipe, NodeRecipe, NetworkRecipe, ChainPhysicalRepresentation, ChainRecipe,
    ChainProtocol, EdgeAndNodes, edges, nodes, edges_and_nodes, num_edges, num_nodes,
    edge_length, physical_representation,
    # performance_metrics.jl
    EvaluationMethod, Analytical, Sampling, DiscreteEventSimulation, generation_duration,
    entangled_state_fidelity, skr_bb84, qber_x, qber_y, qber_z, binary_entropy,
    binary_entropy_derivative, skr_bb84_cost_fct,
    # nodes.jl
    PerfectNode, DepolarizingMemory, PauliType, PauliX, PauliY, PauliZ, PauliMemory,
    AbstractNodeWithMemory, coherence_time, NodeWithMemory,
    AbstractNodeWithPhotonSource, emission_efficiency, cycle_time, num_multiplexing_modes,
    NodeWithPhotonSource, SimpleNode,
    # edges.jl
    EntangledStateType, WernerState, RState,
    NumberOfAttempts, Duration,
    HeraldedEntanglement, SimpleHeraldedPhysical, SimpleHeraldedRecipe,
    success_probability, attempt_duration,
    FiberBasedEdge, FiberSegment, FiberSegmentWithMidpoint, DirectTransmissionOverFiber,
    MidHeraldedDetection, SimplifiedSingleClickElementaryLink,
    # cutoffs.jl
    CutoffProtocol, LocalCutoffProtocol, AbstractNodeWithCutoff, cutoff_time,
    NodeWithCutoff,
    # single_sequential_repeater.jl
    SingleSequentialRepeater,
    # swap_asap.jl
    SwapASAP, SwapASAPRestartPerState, SwapASAPWithoutCutoff, SwapASAPWithCutoff

include("network_representations.jl")
include("performance_metrics.jl")
include("nodes.jl")
include("edges.jl")
include("cutoffs.jl")
include("stochastic_ad_helper_functions.jl")
include("single_sequential_repeater.jl")
include("swap_asap.jl")

end