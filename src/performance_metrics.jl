"""
    EvaluationMethod

Abstract type to specify how to evaluate the performance of a (quantum) network.
"""
abstract type EvaluationMethod end

"""
    Analytical

Calculate network performance directly using analytically-derived functions.

While this method typically has the largest precision and smallest computational overhead,
it may only be available for the simplest of cases.
"""
struct Analytical <: EvaluationMethod end

"""
    Sampling

Estimate network performance by Monte-Carlo simulation through sampling.

All the random elements of the execution of a network protocol are sampled from their
probability distributions, and these samples are used to construct a single execution of the
protocol. Doing so many times allows estimating statistical properties of the protocol.
Note that constructing the execution from the samples and then determining the desired
properties from that construction may require some analytical work still. Hence this method
may not be suited for more complicated cases. It can be thought of as taking the middle
ground between analytical evaluation and discrete-event simulation in terms of how generally
it can be applied, how much analytical work is required to implement the evaluation
method and how computationally expensive it is.
"""
struct Sampling <: EvaluationMethod
    """The number of samples that should be taken per evaluation."""
    number_of_samples::Int
end

"""
    DiscreteEventSimulation

Estimate network performance through discrete-even Monte-Carlo simulation.

Simulate a quantum-network protocol event-by-event while tracking the state of the network,
including the states of quantum systems. This method may be computationally expensive
but can in principle be used to evaluate even the most complicated of network protocols
and hardware models.
"""
struct DiscreteEventSimulation <: EvaluationMethod
    """The number of samples that should be taken per evaluation."""
    number_of_samples::Int
end

"""
    skr_bb84(recipe::NetworkRecipe, method::EvaluationMethod; positive_only::Bool=true)
 
Determine the asymptotic secret-key rate achievable using BB84 quantum-key distribution.

Return both an estimate of the key rate and the standard error in the estimate.

The secret-key rate is the number of bits of secret key that can be produced per second.
It is given by ``skr = max(0, R (1 - h(Q_x) - h(Q_y)))``, where ``R`` is the rate
at which a raw key is produced, ``Q_x`` (``Q_y``) is the Quantum-Bit Error Rate (QBER)
for measurements in the X (Y) basis, and ``h()`` is the binary-entropy function.

Note that the ``max'' in the definition of ``skr'' guarantees that the secret-key rate is
non-negative, as it is impossible to produce a negative amount of secret key.
When the QBER is too large, the key rate just becomes zero.
However, as a result of the ``max'' function, the secret-key rate is not everywhere
differentiable. Moreover, when doing optimizations over the secret-key rate,
the ``max'' can create a large barren plateau that make it hard to find the optimum.
Therefore, it is possible to remove ``max'' by setting the argument `positive_only` to
`false`.

It is here assumed that the sifting rate is unity, i.e., the two parties performing BB84
together always measure in the same basis. In a real QKD session, they pick their bases
randomly so that an eavesdropper cannot predict them, and as a result they would sometimes
measure in different bases. However, by biasing their measurements to one of the bases,
in the asymptotic limit considered here, the fraction of measurement mismatches can be made
negligible.

For more information about the BB84 protocol, see, e.g., 
[arxiv:2002.07305](https://arxiv.org/abs/2002.07305) and references therein.

# Implementation details
While a general version of this function has been defined in terms of the QBER functions
and duration function, specific recipes may define their own implementation of this function
for efficiency gain (for example, one may want to avoid taking separate samples to estimate
each of the QBERs and the duration).

# StochasticAD
When derivatives are extracted using the `StochasticAD` package, by propagating a
`StochasticAD.StochasticTriple` through the function, the first value that is returned is
a dual number of the form "mean value + (mean value of derivative)ϵ".
The second value that is returned is of the form "standard error + (standard error of
derivative)ϵ".
"""
function skr_bb84(recipe::NetworkRecipe, method::EvaluationMethod, positive_only::Bool=true)
    qx, qx_error = qber_x(recipe, method)
    qz, qz_error = qber_z(recipe, method)
    t, t_error = generation_duration(recipe, method)
    skr, skr_error = _skr_bb84(qx, qx_error, qz, qz_error, t, t_error, positive_only)
    return skr, skr_error
end

function _skr_bb84(qber_x, qber_x_error, qber_z, qber_z_error, duration, duration_error,
        positive_only::Bool)
    secret_fraction = 1 - binary_entropy(qber_x) - binary_entropy(qber_z)
    az = iszero(qber_z_error) ? 0. : binary_entropy_derivative(qber_z) * qber_z_error
    ax = iszero(qber_x_error) ? 0. : binary_entropy_derivative(qber_x) * qber_x_error
    secret_fraction_error = sqrt(ax^2 + az^2)
    rate = 1 / duration
    rate_error = duration_error / duration ^ 2
    skr = secret_fraction * rate
    skr_error = sqrt(
        secret_fraction_error ^ 2 * rate ^ 2 + rate_error ^ 2 * secret_fraction ^ 2
    )
    if skr < 0. && positive_only
        skr_error = max(skr + skr_error, 0.)
        skr = 0.
    end
    skr, skr_error
end

"""
    generation_duration(x::NetworkRecipe, method::EvaluationMethod)

Determine the expected value of the duration of entanglement distribution or transmission.

Return both an estimate of the generation duration and the standard error in the estimate.

# StochasticAD
When derivatives are extracted using the `StochasticAD` package, by propagating a
`StochasticAD.StochasticTriple` through the function, the first value that is returned is
a dual number of the form "mean value + (mean value of derivative)ϵ".
The second value that is returned is of the form "standard error + (standard error of
derivative)ϵ".
"""
function generation_duration end
generation_duration(x::EdgeAndNodes) = generation_duration(x[1], x[2])

"""
    qber_x(x::NetworkRecipe, method::EvaluationMethod)

Determine the Quantum Bit Error Rate (QBER) for measurements in the X basis.

Return both an estimate of the QBER and the standard error in the estimate.

This is the probability that if both parties use the X measurement basis, they receive
the same result.

# StochasticAD
When derivatives are extracted using the `StochasticAD` package, by propagating a
`StochasticAD.StochasticTriple` through the function, the first value that is returned is
a dual number of the form "mean value + (mean value of derivative)ϵ".
The second value that is returned is of the form "standard error + (standard error of
derivative)ϵ".
"""
function qber_x end

"""
    qber_y(x::NetworkRecipe, method::EvaluationMethod)

Determine the Quantum Bit Error Rate (QBER) for measurements in the Y basis.

Return both an estimate of the QBER and the standard error in the estimate.

This is the probability that if both parties use the Y measurement basis, they receive
the same result.

# StochasticAD
When derivatives are extracted using the `StochasticAD` package, by propagating a
`StochasticAD.StochasticTriple` through the function, the first value that is returned is
a dual number of the form "mean value + (mean value of derivative)ϵ".
The second value that is returned is of the form "standard error + (standard error of
derivative)ϵ".
"""
function qber_y end

"""
    qber_z(x::NetworkRecipe, method::EvaluationMethod)

Determine the Quantum Bit Error Rate (QBER) for measurements in the Z basis.

Return both an estimate of the QBER and the standard error in the estimate.

This is the probability that if both parties use the Z measurement basis, they receive
the same result.

# StochasticAD
When derivatives are extracted using the `StochasticAD` package, by propagating a
`StochasticAD.StochasticTriple` through the function, the first value that is returned is
a dual number of the form "mean value + (mean value of derivative)ϵ".
The second value that is returned is of the form "standard error + (standard error of
derivative)ϵ".
"""
function qber_z end

"""
    binary_entropy(p)

Calculate the binary entropy of the probability `p`.
"""
function binary_entropy(p)
    0 <= p <= 1 || throw(DomainError("p must be between 0 and 1, not $p."))
    iszero(p) && return zero(p)
    iszero(1 - p) && return zero(p)
    -p * log2(p) - (1 - p)log2(1 - p)
end

"""
    binary_entropy_derivative(p)

Calculate the derivative of the binary entropy of the probability `p`.
"""
binary_entropy_derivative(p) = -log2(p / (1 - p))

"""
    entangled_state_fidelity(x::NetworkRecipe, method::EvaluationMethod)

Fidelity of entangled states that are shared w.r.t. target state.

Return both an estimate of the fidelity and the standard error in the estimate.

We here use the definition for fidelity between density matrices ``ρ`` and ``σ``
``F(ρ, σ) = (tr √[√ρ σ √ρ])^2``. Note that, for pure states, this expression reduces to
``F(ρ, |ψ><ψ|) = <ψ|ρ|ψ>``.

If an `EdgeAndNodes` is provided, the edge and nodes are automatically unpacked.

# StochasticAD
When derivatives are extracted using the `StochasticAD` package, by propagating a
`StochasticAD.StochasticTriple` through the function, the first value that is returned is
a dual number of the form "mean value + (mean value of derivative)ϵ".
The second value that is returned is of the form "standard error + (standard error of
derivative)ϵ".
"""
function entangled_state_fidelity end

"""
    _calculate_mean_and_error(samples)
    _calculate_mean_and_error(samples::Vector{<:StochasticTriple})

Calculate both the mean and the standard error from a list of samples.

When the samples are `StochasticAD.StochasticTriple`s, propagating them (after smoothing)
through the catch-all method would return a mean of the type
"mean value + (mean value of derivative)ϵ" and a standard error of the type
"standard error + (derivative of the standard error)ϵ".
That is, the last value quantifies the rate of change of the size of the standard error.

Instead, this function has a special implementation such that the second returned value
is instead of the type "standard error + (standard error of the derivative)ϵ".
This function should be used when implementing performance metrics based on sampling,
such that performance-metric functions can be used to obtain not only estimates and
standard errors and estimates of derivatives, but also standard errors in the derivative
estimates.
"""
function _calculate_mean_and_error(samples)
    mean(samples), std(samples) / sqrt(length(samples))
end

function _calculate_mean_and_error(samples::Vector{<:StochasticTriple})
    primal_samples = StochasticAD.value.(samples)
    primal_mean, primal_error = _calculate_mean_and_error(primal_samples)
    deriv_samples = StochasticAD.derivative_contribution.(samples)
    deriv_mean, deriv_error = _calculate_mean_and_error(deriv_samples)
    make_dual_number(primal_mean, deriv_mean), make_dual_number(primal_error, deriv_error)
end