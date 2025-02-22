"""
    SingleSequentialRepeater

Recipe for a quantum network consisting of one sequential repeater and two end nodes.

The sequential repeater generates entanglement first along the least-lossy link, stores it
in memory, generates entanglement along the second link, and then performs entanglement
swapping.
"""
struct SingleSequentialRepeater{E<:EdgeRecipe, N<:NodeRecipe} <: NetworkRecipe
    edges::Vector{E}
    nodes::Vector{N}
    function SingleSequentialRepeater(edges, nodes)
        length(edges) == 2 || throw(DomainError)
        length(nodes) == 3 || throw(DomainError)
        new{eltype(edges), eltype(nodes)}(edges, nodes)
    end
end
edges(x::SingleSequentialRepeater) = x.edges
nodes(x::SingleSequentialRepeater) = x.nodes
edges_and_nodes(x::SingleSequentialRepeater) =
    [(x.edges[1], (x.nodes[1], x.nodes[2])), (x.edges[2], (x.nodes[2], x.nodes[3]))]
physical_representation(x::SingleSequentialRepeater) = ChainPhysicalRepresentation(
    physical_representation.(edges(x)), physical_representation.(nodes(x))
)

function generation_duration(x::SingleSequentialRepeater{E, N}, ::Analytical) where
    {E<:HeraldedEntanglement, N<:SimpleNode}
    duration = sum(generation_duration(e) for e in edges_and_nodes(x))
    error = 0.
    duration, error
end

function sample_link_durations(x::SingleSequentialRepeater{E, N}, number_of_samples) where
        {E<:HeraldedEntanglement, N}
    t_random_vars = [Duration(e) for e in edges_and_nodes(x)]
    samples = Matrix{Real}(undef, length(t_random_vars), number_of_samples)
    for col in eachcol(samples)
        col .= rand.(t_random_vars)
    end
    samples
end

function generation_duration(x::SingleSequentialRepeater{E, N}, method::Sampling) where
        {E<:HeraldedEntanglement, N}
    ts = sample_link_durations(x, method.number_of_samples)
    generation_duration(x, method, ts)
end

function generation_duration(::SingleSequentialRepeater{E, N}, ::Sampling,
        link_duration_samples) where
        {E<:HeraldedEntanglement, N<:SimpleNode}
    durations = sum.(eachcol(link_duration_samples))
    _calculate_mean_and_error(durations)
end

"""
    werner_parameter(x::SingleSequentialRepeater{E, N}[, method::EvaluationMethod]) where
        {E<:HeraldedEntanglement{T, WernerState},
        N<:SimpleNode{NodeWithMemory{DepolarizingMemory}}}
    
Determine the werner parameter of the states created by the network.

For this setup, all the noise is of the depolarizing form, and hence the final state is a
Werner state (see also `WernerState`). The Werner parameter of the state can be determined
using this function, which is helpful to calculate function such as the fidelity and the
quantum-bit error rate (QBER).

Defaults to analytical evaluation method.
"""
function werner_parameter(x::SingleSequentialRepeater{E, SimpleNode{PerfectNode}},
        ::EvaluationMethod) where {E<:HeraldedEntanglement{T, WernerState}} where T
    werner_param = prod(werner_parameter_from_fidelity(entangled_state_fidelity(e)) for
        e in edges_and_nodes(x))
    error = 0.
    werner_param, error
end

werner_parameter(x::SingleSequentialRepeater{E, SimpleNode{N}}) where
        {T, E<:HeraldedEntanglement{T, WernerState},
        N<:AbstractNodeWithMemory{DepolarizingMemory}} =
    werner_parameter(x, Analytical())

function werner_parameter(x::SingleSequentialRepeater{E, SimpleNode{N}}, ::Analytical) where
        {T, E<:HeraldedEntanglement{T, WernerState},
        N<:AbstractNodeWithMemory{DepolarizingMemory}}
    edge_werner_parameters = [werner_parameter_from_fidelity(entangled_state_fidelity(e))
        for e in edges_and_nodes(x)]
    slow_edge_and_nodes = sort(edges_and_nodes(x), by=generation_duration)[2]
    p = success_probability(slow_edge_and_nodes)
    t = attempt_duration(slow_edge_and_nodes)
    coherence_times = [physical_representation(n).coherence_time
        for n in slow_edge_and_nodes[2]]
    depol_param_per_attempt = prod(exp(-t / s) for s in coherence_times)
    avg_depol_param = mean_exponential_number_of_attempts(p, depol_param_per_attempt)
    werner_param = avg_depol_param * prod(edge_werner_parameters)
    error = 0.
    werner_param, error
end

function werner_parameter(x::SingleSequentialRepeater{E, SimpleNode{N}},
        method::Sampling) where
        {T, E<:HeraldedEntanglement{T, WernerState},
        N<:AbstractNodeWithMemory{DepolarizingMemory}}
    slow_edge_and_nodes = sort(edges_and_nodes(x), by = generation_duration)[2]
    ts = [rand(Duration(slow_edge_and_nodes)) for _ in 1:method.number_of_samples]
    werner_parameter(x, method, ts)
end

function werner_parameter(x::SingleSequentialRepeater{E, SimpleNode{N}},
        ::Sampling, slow_link_duration_samples) where
        {T, E<:HeraldedEntanglement{T, WernerState},
        N<:AbstractNodeWithMemory{DepolarizingMemory}}
    edge_werner_parameters = [werner_parameter_from_fidelity(entangled_state_fidelity(e))
        for e in edges_and_nodes(x)]
    slow_edge_and_nodes = sort(edges_and_nodes(x), by = generation_duration)[2]
    coherence_times = [physical_representation(n).coherence_time
        for n in slow_edge_and_nodes[2]]
    effective_coherence_time = 1 / mean(1 ./ coherence_times)
    werner_params = [exp(- t / effective_coherence_time) for t in
        slow_link_duration_samples] * prod(edge_werner_parameters)
    _calculate_mean_and_error(werner_params)
end
    
function entangled_state_fidelity(x::SingleSequentialRepeater{E, N},
        method::EvaluationMethod) where {T, E<:HeraldedEntanglement{T, WernerState},
        N<:WernerCompatibleNode}
    fidelity_from_werner_parameter(werner_parameter(x, method)...)
end

function qber_x(x::SingleSequentialRepeater{E, N},
        method::EvaluationMethod) where {T, E<:HeraldedEntanglement{T, WernerState},
        N<:WernerCompatibleNode}
    werner_param, werner_param_error = werner_parameter(x, method)
    qber = (1 - werner_param) / 2
    qber_error = werner_param_error / 2
    qber, qber_error
end

function qber_y(x::SingleSequentialRepeater{E, N},
        method::EvaluationMethod) where {T, E<:HeraldedEntanglement{T, WernerState},
        N<:WernerCompatibleNode}
    qber_x(x, method)
end

function qber_z(x::SingleSequentialRepeater{E, N},
        method::EvaluationMethod) where {T, E<:HeraldedEntanglement{T, WernerState},
        N<:WernerCompatibleNode}
    qber_x(x, method)
end

function skr_bb84(recipe::SingleSequentialRepeater{E, N}, method::Sampling,
        positive_only::Bool=true) where
        {T, E<:HeraldedEntanglement{T, WernerState},
        N<:AbstractNodeWithMemory{DepolarizingMemory}}
    link_duration_samples = sample_link_durations(recipe, method.number_of_samples)
    duration, duration_error = generation_duration(recipe, method, link_duration_samples)
    index_slow_link = argmax(generation_duration.(edges_and_nodes(recipe)))
    slow_link_duration_samples = link_duration_samples[index_slow_link, :]
    werner_param, werner_param_error =
        werner_parameter(recipe, method, slow_link_duration_samples)
    qber = (1 - werner_param) / 2
    qber_error = werner_param_error / 2
    _skr_bb84(qber, qber_error, qber, qber_error, duration, duration_error, positive_only)
end