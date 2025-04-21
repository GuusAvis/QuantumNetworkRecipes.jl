"""
    SwapASAP

Protocol for a quantum-repeater chain where nodes swap entanglement as soon as possible.

It is assumed that all the links can be used to generate entanglement in parallel.
We also assume that the nodes and edges become idle after swapping entanglement,
and that a new round of the protocol is only started after end-to-end entanglement has been
successfully created.
Decoherence takes place at the end nodes until the final state is ready;
to avoid this, set their coherence times to infinity.
"""
abstract type SwapASAP <: ChainProtocol end

"""
    SwapASAPRestartPerState

Swap-ASAP protocol where the repeater chain restarts for each new entangled state.

Each elementary link produces entanglement once and then becomes idle.
When all links have been produced, there will be an end-to-end entangled state due to
entanglement swapping.
Only then, if more entanglement is desired, do the elementary links all restart
simultaneously to start producing the next end-to-end entangled state.
Entangled states are never discarded, i.e., there is no cutoff policy.
While it may be inefficient that links are idle after producing entanglement, this protocol
is relatively simple and hence easy to understand and study.
It also makes sure that the statistics of entanglement distribution are independent.
"""
struct SwapASAPRestartPerState <: SwapASAP end

"""
    SwapASAPWithoutCutoff

Swap-ASAP protocol without a cutoff policy. Links are regenerated as soon as possible.

In this protocol, elementary links are never needlessly idle; whenever both nodes are free,
they will restart entanglement generation. This means that when end-to-end entanglement
is delivered, there may already be some entanglement in the repeater chain that can be
used to generate the next end-to-end entangled state (sometimes referred to as "residual
entanglement"), making this protocol more efficient than `SwapASAPRestartPerState`
(though it may perform worse in terms of fidelity, as the residual entanglement undergoes
decoherence while waiting for the end-to-end entangled state to be created).
Note that the statistics of different end-to-end entangled states are not independent in
this protocol (for instance, the time it takes to generate the second end-to-end entangled
state is correlated to the time required for the first).
"""
struct SwapASAPWithoutCutoff <: SwapASAP end

"""
    SwapASAPWithCutoff{T<:CutoffProtocol}

Swap-ASAP protocol with a cutoff policy. Links are regenerated as soon as possible.
"""
struct SwapASAPWithCutoff{T<:CutoffProtocol} <: SwapASAP
    cutoff_protocol::T
end

"""
    sample_link_durations(x::ChainRecipe{E, N, P}, number_of_samples) where
        {E<:HeraldedEntanglement, N, P<:SwapASAP}

Draw a number of samples for the durations of generating each of the links in the chain.

The samples are provided as an m-by-n matrix, where m is the number of links in the chain and n is
the number of samples. Each entry represents the duration of generating the corresponding link.
"""
function sample_link_durations end

function sample_link_durations(x::ChainRecipe{E, N, SwapASAPRestartPerState},
        number_of_samples) where {E<:HeraldedEntanglement, N}
    random_vars = [Duration(e) for e in edges_and_nodes(x)]
    samples = Matrix{typeof(rand(random_vars[1]))}(undef, length(random_vars), number_of_samples)
    for col in eachcol(samples)
        col .= rand.(random_vars)
    end
    samples
end

"""
    _find_free_time(sample, link_index)

Determine at what time the connection at `link_index` is free to generate a new link.

A link can be generated only if the communication qubits at both nodes are free.
For most links, this means that swaps need to have taken place at two different repeaters,
which can only be done if both the link itself and its two neighbors have been generated.
At the edges, we assume that end nodes hold on to their qubits until end-to-end entanglement
is confirmed, and hence these links are only free for the next pair after all links are
generated.

# Examples
```jldoctest
julia> sample = [0, 6, 3, 4, 9]
5-element Vector{Int64}:
 0
 6
 3
 4
 9

julia> _find_free_time(sample, 1)
9

julia> _find_free_time(sample, 5)
9

julia> _find_free_time(sample, 2)
6

julia> _find_free_time(sample, 3)
6
```
"""
function _find_free_time(sample, link_index)
    if link_index == 1 || link_index == length(sample)
        return maximum(sample)
    else
        return max(sample[link_index - 1], sample[link_index], sample[link_index + 1])
    end
end

function sample_link_durations(x::ChainRecipe{E, N, SwapASAPWithoutCutoff},
        number_of_samples) where {E<:HeraldedEntanglement, N}
    random_variables = [Duration(e) for e in edges_and_nodes(x)]
    samples = Matrix{typeof(rand(random_variables[1]))}(undef, length(random_variables), number_of_samples)
    samples[:, 1] .= rand.(random_variables)
    for j in 2:number_of_samples
        previous_sample = view(samples, :, j-1)
        end_time_previous_sample = maximum(previous_sample)
        for (i, var) in enumerate(random_variables)
            # as free_times may be earlier than the end of the previous sample, some links may
            # have a head start in entanglement generation
            # as a result, some of the links may be finished already at negative times (i.e., 
            # before the previous end-to-end link was finished)
            samples[i, j] = _find_free_time(previous_sample, i) + rand(var) - end_time_previous_sample
        end
    end
    samples
end

"""
    _find_first_cutoff_link(sample, cutoff_times)

Determine which link triggers its cutoff time first, and when.

Returns both the index of the link and the time at which the cutoff occurs.
If no cutoff is triggered, the index is returned as `0` and the cutoff moment `false`.
This function assumes that cutoffs only occur at repeater nodes, hence `cutoff_times` should
only contain cutoff times for the repeaters (i.e., its length should be one less than that
of `sample`).

# Examples
```jldoctest
julia> sample = [4, 1, 3, 8]
4-element Vector{Int64}:
 4
 1
 3
 8

julia> _find_first_cutoff_link(sample, [10 for _ in 1:3])
(0, false)

julia> _find_first_cutoff_link(sample, [1 for _ in 1:3])
(2, 2)

julia> _find_first_cutoff_link(sample, [3 for _ in 1:3])
(3, 6)
```
"""
function _find_first_cutoff_link(sample, cutoff_times)
    length(cutoff_times) == length(sample) - 1 || throw(DimensionMismatch)

    index = 0
    first_cutoff_moment = Inf
    for i in 1:(length(sample) - 1)
        waiting_time = abs(sample[i + 1] - sample[i])
        waiting_time <= cutoff_times[i] && continue  # cutoff not triggered
        cutoff_moment = min(sample[i + 1], sample[i]) + cutoff_times[i]
        cutoff_moment >= first_cutoff_moment && continue  # not the first cutoff
        # this is the earliest cutoff moment found so far
        first_cutoff_moment = cutoff_moment
        index = sample[i] <= sample[i + 1] ? i : i + 1
    end
    first_cutoff_moment == Inf && (first_cutoff_moment = false)
    index, first_cutoff_moment
end

"""
    _find_swapped_links(sample, link_index, time)

Determine links in `sample` that were swapped with the one at `link_index` at time `time`.

Returns a range containing the indices of the swapped links.
If `sample[link_index] > time`, return an empty range. Otherwise, it always includes `link_index`.
Assumes a swap-asap protocol is performed.

# Examples
```jldoctest
julia> sample = [0, 6, 3, 4, 9]
5-element Vector{Int64}:
 0
 6
 3
 4
 9

julia> _find_swapped_links(sample, 3, 0)  # link cannot be swapped before it is created
3:2

julia> _find_swapped_links(sample, 3, 10)  # after last link generated, all is swapped
1:5

julia> _find_swapped_links(sample, 3, 4)  # swap happened on the right, but not on the left
3:4

julia> _find_swapped_links(sample, 3, 6)  # swapping on the left also gives link 1
1:4
```
"""
function _find_swapped_links(sample, link_index, time)
    sample[link_index] > time && return link_index:link_index-1
    i = link_index
    while i > 1 && sample[i-1] <= time
        i -= 1
    end
    j = link_index
    while j < length(sample) && sample[j+1] <= time
        j += 1
    end
    i:j
end

"""
    _resample_for_earliest_cutoff(sample, cutoff_times, duration_samples)

Create updated sample to account for the earliest cutoff that takes place.

A qubit is discarded when it has been stored longer than the cutoff time of its node,
and with it the entangled state it was part of, which may consist of several swapped
elementary links. For each of these links the generation value is increased to represent
the generation of a new link, starting at the time when both qubits were free.

This function resamples only the earliest cutoff triggered according to the link-generation
times contained in the sample. If no cutoff is triggered, the sample is not changed.
The function returns `(true, updated_sample)` if resampling took place, and
`(false, old_sample)` if there was no cutoff.
"""
function _resample_for_earliest_cutoff(sample, cutoff_times, duration_samples)
    link_index, cutoff_moment = _find_first_cutoff_link(sample, cutoff_times)
    iszero(link_index) && return false, sample
    swapped_links = _find_swapped_links(sample, link_index, cutoff_moment)
    new_sample = copy(sample)
    for index in swapped_links
        free_time = _find_free_time(sample, index)
        restart_time = min(free_time, cutoff_moment)
        new_sample[index] = restart_time + duration_samples[index]
    end
    true, new_sample
end

function _resample_for_earliest_cutoff(sample::AbstractVector{T}, cutoff_times, duration_samples) where
        T<:StochasticAD.StochasticTriple
    StochasticAD.propagate(_resample_for_earliest_cutoff, sample, cutoff_times, duration_samples)
end
function _resample_for_earliest_cutoff(sample, cutoff_times::AbstractVector{T}, duration_samples) where
        T<:StochasticAD.StochasticTriple
    StochasticAD.propagate(_resample_for_earliest_cutoff, sample, cutoff_times, duration_samples)
end
function _resample_for_earliest_cutoff(sample, cutoff_times, duration_samples::AbstractVector{T}) where
        T<:StochasticAD.StochasticTriple
    StochasticAD.propagate(_resample_for_earliest_cutoff, sample, cutoff_times, duration_samples)
end
function _resample_for_earliest_cutoff(sample::AbstractVector{T}, cutoff_times::AbstractVector{S}, duration_samples) where
        {T<:StochasticAD.StochasticTriple, S<:StochasticAD.StochasticTriple}
    StochasticAD.propagate(_resample_for_earliest_cutoff, sample, cutoff_times, duration_samples)
end
function _resample_for_earliest_cutoff(sample::AbstractVector{T}, cutoff_times, duration_samples::AbstractVector{S}) where
        {T<:StochasticAD.StochasticTriple, S<:StochasticAD.StochasticTriple}
    StochasticAD.propagate(_resample_for_earliest_cutoff, sample, cutoff_times, duration_samples)
end
function _resample_for_earliest_cutoff(sample, cutoff_times::AbstractVector{T}, duration_samples::AbstractVector{S}) where
        {T<:StochasticAD.StochasticTriple, S<:StochasticAD.StochasticTriple}
    StochasticAD.propagate(_resample_for_earliest_cutoff, sample, cutoff_times, duration_samples)
end
function _resample_for_earliest_cutoff(sample::AbstractVector{T}, cutoff_times::AbstractVector{S}, duration_samples::AbstractVector{U}) where
        {T<:StochasticAD.StochasticTriple, S<:StochasticAD.StochasticTriple, U<:StochasticAD.StochasticTriple}
    StochasticAD.propagate(_resample_for_earliest_cutoff, sample, cutoff_times, duration_samples)
end

"""
    _resample_until_no_cutoff(sample, random_variables, cutoff_times)

Resample for the earliest cutoff until it does not trigger a cutoff.

This function can be used to simulate a swap-asap protocol with a cutoff policy.
Given an initial sample of link-generation times, it repeatedly updates the sample until it
contains only the generation times of the links used in the final end-to-end state.
"""
function _resample_until_no_cutoff(sample, random_variables, cutoff_times)
    resample = true
    new_sample = promote(sample, rand.(random_variables), cutoff_times)[1]
    # Stop when every branch is done resampling
    while !isfalse(resample)
        # Give every branch the same random numbers
        duration_samples = rand.(random_variables)
        resample, new_sample = _resample_for_earliest_cutoff(new_sample, cutoff_times, duration_samples)
    end
    new_sample
end

# Custom dispatch for efficiency
function _resample_until_no_cutoff(sample::AbstractVector{Float64}, random_variables::Vector{Duration{Float64, Float64}}, cutoff_times::Vector{Float64})
    cur_sample = copy(sample)
    next_sample = copy(sample)
    while true
        link_index, cutoff_moment = _find_first_cutoff_link(cur_sample, cutoff_times)
        iszero(link_index) && return cur_sample
        swapped_links = _find_swapped_links(cur_sample, link_index, cutoff_moment)
        for index in swapped_links
            free_time = _find_free_time(cur_sample, index)
            restart_time = min(free_time, cutoff_moment)
            next_sample[index] = restart_time + rand(random_variables[index])
        end
        cur_sample[swapped_links] = view(next_sample, swapped_links)
    end
end

function sample_link_durations(
        x::ChainRecipe{E, N, SwapASAPWithCutoff{LocalCutoffProtocol{T}}}, number_of_samples
        ) where {E<:HeraldedEntanglement, N, T<:Real}
    cutoff_times = [
        n isa AbstractNodeWithCutoff ? cutoff_time(n) :
            x.chain_protocol.cutoff_protocol.cutoff_time
        for n in nodes(x)[begin + 1:end - 1]
    ]
    random_variables = [Duration(e) for e in edges_and_nodes(x)]
    output_type = typeof(sum(rand.(random_variables)) + sum(cutoff_times))
    samples = Matrix{output_type}(undef, length(random_variables), number_of_samples)
    samples[:, 1] = _resample_until_no_cutoff(rand.(random_variables), random_variables, cutoff_times)
    for j in 2:number_of_samples
        previous_sample = view(samples, :, j-1)
        end_time_previous_sample = maximum(previous_sample)
        for (i, var) in enumerate(random_variables)
            samples[i, j] = _find_free_time(previous_sample, i) + rand(var) - end_time_previous_sample
        end
        samples[:, j] = _resample_until_no_cutoff(samples[:, j], random_variables, cutoff_times)
    end
    samples
end

function generation_duration(x::ChainRecipe{E, N, P}, method::Sampling) where
        {E, N, P<:SwapASAP}
    ts = sample_link_durations(x, method.number_of_samples)
    generation_duration(x, method, ts)
end

function generation_duration(::ChainRecipe{E, N, P}, ::Sampling,
        link_duration_samples) where {E, N, P<:SwapASAP}
    durations = maximum.(eachcol(link_duration_samples))
    _calculate_mean_and_error(durations)
end

function werner_parameter(x::ChainRecipe{E, SimpleNode{PerfectNode}, P},
        ::EvaluationMethod) where {T, E<:HeraldedEntanglement{T, WernerState},
        P<:SwapASAP}
    werner_param = prod(werner_parameter_from_fidelity(entangled_state_fidelity(e))
        for e in edges_and_nodes(x))
    error = 0  # with perfect memory, there is no randomness
    werner_param, error
end

"""
    mem_depolar_param_from_link_durations(x::ChainRecipe, durations)

Obtain a single sample for the total memory depolarizing parameter from single samples for
the durations of generating each of the links.
"""
function mem_depolar_param_from_link_durations(x::ChainRecipe{E, SimpleNode{N}, P},
        durations) where {E, N<:AbstractNodeWithMemory{DepolarizingMemory}, P<:SwapASAP}
    num_edges(x) == length(durations) || throw(DimensionMismatch())
    num_repeaters = num_nodes(x) - 2  # all nodes except the two end nodes are repeaters
    werner_param = 1
    for i in 1:num_repeaters
        repeater_storage_time = abs(durations[i+1] - durations[i])
        repeater = nodes(x)[i+1]
        repeater_depolar_param = exp(-repeater_storage_time / coherence_time(repeater))
        werner_param *= repeater_depolar_param
    end
    completion_time = maximum(durations)
    end_node_1_depolar_param = exp(-(completion_time - durations[begin]) /
        coherence_time(nodes(x)[begin]))
    end_node_2_depolar_param = exp(-(completion_time - durations[end]) /
        coherence_time(nodes(x)[end]))
    werner_param * end_node_1_depolar_param * end_node_2_depolar_param
end

function werner_parameter(x::ChainRecipe{E, SimpleNode{N}, P},
        method::Sampling) where {T, E<:HeraldedEntanglement{T, WernerState},
        N<:AbstractNodeWithMemory{DepolarizingMemory}, P<:SwapASAP}
    ts = sample_link_durations(x, method.number_of_samples)
    werner_parameter(x, method, ts)
end

function werner_parameter(x::ChainRecipe{E, SimpleNode{N}, P}, ::Sampling,
        link_duration_samples) where
        {T, E<:HeraldedEntanglement{T, WernerState},
        N<:AbstractNodeWithMemory{DepolarizingMemory}, P<:SwapASAP}
    werner_param_links = prod(werner_parameter_from_fidelity(entangled_state_fidelity(e))
        for e in edges_and_nodes(x))
    ts = eachcol(link_duration_samples)
    mem_depolar_params = [mem_depolar_param_from_link_durations(x, t) for t in ts]
    werner_params = mem_depolar_params * werner_param_links
    _calculate_mean_and_error(werner_params)
end

function entangled_state_fidelity(x::ChainRecipe{E, N, P},
        method::EvaluationMethod) where {T, E<:HeraldedEntanglement{T, WernerState},
        N<:WernerCompatibleNode, P<:SwapASAP}
    fidelity_from_werner_parameter(werner_parameter(x, method)...)
end

function qber_x(x::ChainRecipe{E, N, P},
        method::EvaluationMethod) where {T, E<:HeraldedEntanglement{T, WernerState},
        N<:WernerCompatibleNode, P<:SwapASAP}
    werner_param, werner_param_error = werner_parameter(x, method)
    qber = (1 - werner_param) / 2
    qber_error = werner_param_error / 2
    qber, qber_error
end

function qber_y(x::ChainRecipe{E, N, P},
        method::EvaluationMethod) where {T, E<:HeraldedEntanglement{T, WernerState},
        N<:WernerCompatibleNode, P<:SwapASAP}
    qber_x(x, method)
end

function qber_z(x::ChainRecipe{E, N, P},
        method::EvaluationMethod) where {T, E<:HeraldedEntanglement{T, WernerState},
        N<:WernerCompatibleNode, P<:SwapASAP}
    qber_x(x, method)
end

function skr_bb84(recipe::ChainRecipe{E, N, P}, method::Sampling,
        positive_only::Bool=true) where
        {T, E<:HeraldedEntanglement{T, WernerState}, N<:WernerCompatibleNode, P<:SwapASAP}
    link_duration_samples = sample_link_durations(recipe, method.number_of_samples)
    duration, duration_error = generation_duration(recipe, method, link_duration_samples)
    werner_param, werner_param_error =
        werner_parameter(recipe, method, link_duration_samples)
    qber = (1 - werner_param) / 2
    qber_error = werner_param_error / 2
    _skr_bb84(qber, qber_error, qber, qber_error, duration, duration_error, positive_only)
end