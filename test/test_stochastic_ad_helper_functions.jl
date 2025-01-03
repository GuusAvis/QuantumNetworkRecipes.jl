@testset "istrue" begin
    f = QuantumNetworkRecipes.istrue
    @test f(true)
    @test !f(false)
    @test f(stoch_trip(1, 0., 0, 1.))
    @test !f(stoch_trip(0, 0., 0, 1.))
    @test !f(stoch_trip(1, 0., -1, 1.))
    @test !f(stoch_trip(0, 0., 1, 1.))
end

@testset "isfalse" begin
    f = QuantumNetworkRecipes.isfalse
    @test !f(true)
    @test f(false)
    @test !f(stoch_trip(1, 0., 0, 1.))
    @test f(stoch_trip(0, 0., 0, 1.))
    @test !f(stoch_trip(1, 0., -1, 1.))
    @test !f(stoch_trip(0, 0., 1, 1.))
end

@testset "dual number" begin
    f = QuantumNetworkRecipes.make_dual_number
    primal = 1.
    deriv = 5.
    x = f(primal, deriv)
    @test x isa StochasticAD.StochasticTriple
    @test value(x) == primal
    @test derivative_contribution(x) == StochasticAD.delta(x) == deriv
end