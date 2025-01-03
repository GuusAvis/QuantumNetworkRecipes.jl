import QuantumNetworkRecipes: _skr_bb84


# Test EvaluationMethod abstract type
@testset "EvaluationMethod" begin
    @test typeof(Analytical()) <: EvaluationMethod
    @test typeof(Sampling(100)) <: EvaluationMethod
    @test Sampling(100).number_of_samples == 100
    @test typeof(DiscreteEventSimulation(100)) <: EvaluationMethod
    @test DiscreteEventSimulation(100).number_of_samples == 100
end

# define a mock network recipe for testing
struct MockNetwork <: NetworkRecipe
    error_probability
    generation_duration
    generation_duration_error
end
MockNetwork(error_probability) = MockNetwork(error_probability, 1., 0.)
QuantumNetworkRecipes.qber_x(x::MockNetwork, ::Analytical) = (x.error_probability / 2, 0.)
QuantumNetworkRecipes.qber_y(x::MockNetwork, ::Analytical) = (x.error_probability / 2, 0.)
QuantumNetworkRecipes.qber_z(x::MockNetwork, ::Analytical) = (x.error_probability / 2, 0.)
QuantumNetworkRecipes.generation_duration(x::MockNetwork, ::Analytical) =
    (x.generation_duration, x.generation_duration_error)
QuantumNetworkRecipes.entangled_state_fidelity(x::MockNetwork, ::Analytical) =
    (1 - 3 * x.error_probability / 4, 0.)

# Test skr_bb84 function
@testset "skr_bb84" begin
    @test skr_bb84(MockNetwork(0), Analytical()) == (1., 0.)
    @test skr_bb84(MockNetwork(0), Analytical(), true) == (1., 0.)
    @test skr_bb84(MockNetwork(0), Analytical(), false) == (1., 0.)
    @test skr_bb84(MockNetwork(1), Analytical(), false) == (-1., 0.)
    @test skr_bb84(MockNetwork(1), Analytical(), true) == (0., 0.)
    @test skr_bb84(MockNetwork(1), Analytical()) == (0., 0.)
    for only_positive in [true, false]
        @test _skr_bb84(0., 0., 0., 0., 5., 0., only_positive) ==
            skr_bb84(MockNetwork(0., 5., 0.), Analytical(), only_positive) == (1. / 5., 0.)
        @test _skr_bb84(0., 0., 0., 0., 1., 1., only_positive) ==
            skr_bb84(MockNetwork(0., 1., 1.), Analytical(), only_positive) == (1., 1.)
    end
    @test skr_bb84(MockNetwork(0., -1., 0.5), Analytical()) == (0., 0.)
    @test skr_bb84(MockNetwork(0., -1., 1.), Analytical()) == (0., 0.)
    @test skr_bb84(MockNetwork(0., -1., 2.), Analytical()) == (0., 1.)
    for only_positive in [true, false]
        @test iszero(_skr_bb84(0.2, 0., 0.3, 0., 32, 0., only_positive)[2])
        @test _skr_bb84(0.2, 0.04, 0., 0., 3., 0.2, only_positive) ==
            _skr_bb84(0., 0., 0.2, 0.04, 3., 0.2, only_positive)
    end

    interval_zero_to_half = _skr_bb84(0.1, 0.01, 0., 0., 5., 0.2, false)
    interval_half_to_one = _skr_bb84(0.9, 0.01, 0., 0., 5., 0.2, false)
    @test interval_zero_to_half[1] ≈ interval_half_to_one[1]
    @test interval_zero_to_half[2] ≈ interval_half_to_one[2]

    resultx = _skr_bb84(0.1, 0.01, 0., 0., 1., 0., false)
    resultz = _skr_bb84(0., 0., 0.1, 0.01, 1., 0., false)
    expected = (1 - binary_entropy(0.1), binary_entropy_derivative(0.1) * 0.01)
    @test resultx[1] ≈ expected[1] ≈ resultz[1]
    @test resultx[2] ≈ expected[2] ≈ resultz[2]
end

# Test binary_entropy function
@testset "binary_entropy" begin
    @test binary_entropy(0.5) ≈ 1.0
    @test binary_entropy(0.0) ≈ 0.0
    @test isapprox(binary_entropy(1E-7), 0.0; atol = 1E-5)
    @test binary_entropy(1.0) ≈ 0.0
    @test isapprox(binary_entropy(1 - 1E-7), 0.0; atol = 1E-5)
end

# Test binary_entropy_derivative function
@testset "binary_entropy_derivative" begin
    @test binary_entropy_derivative(0.5) ≈ 0.0
    @test binary_entropy_derivative(0.0) ≈ Inf
    @test binary_entropy_derivative(1.0) ≈ -Inf
end

@testset "calculate mean and errors" begin
    # normal numbers
    f = QuantumNetworkRecipes._calculate_mean_and_error
    x = [1., 1., 1.]
    @test f(x) == (1., 0.)
    x = [1., 2., 3.]  # samples with mean 2 and variance 1 (sample variance)
    @test f(x) == (2., 1. / sqrt(3))
    x = [10., 20., 30.]  # mean 20, variance 100
    @test f(x) == (20., 10. / sqrt(3))

    # dual numbers
    g = QuantumNetworkRecipes.make_dual_number
    x = [g(1., 30.), g(2., 20.), g(3., 10.)]
    # Primal: mean 11. Derivative: mean 101. Both standard error 1 / sqrt(3).
    res = f(x)
    @test value.(res) == (2., 1. / sqrt(3))
    @test derivative_contribution.(res) == (20., 10. / sqrt(3))

    # stochastic triples
    x = [
        stoch_trip(1., 30., 0., 0.),
        stoch_trip(2., 0., 40., 0.5),
        stoch_trip(3., 0., 50., 0.2)
        ]
    res = f(x)
    @test value.(res) == (2., 1. / sqrt(3))
    @test derivative_contribution.(res) == (20., 10. / sqrt(3))
end