@testset "NumberOfAttempts" begin

    number_of_attempts = NumberOfAttempts(.1)
    @test number_of_attempts isa NumberOfAttempts{Float64}
    @test mean(number_of_attempts) == 10
    @test median(number_of_attempts) > 1
    @test mode(number_of_attempts) == 1
    @test var(number_of_attempts) > 0
    @test skewness(number_of_attempts) != 0
    @test kurtosis(number_of_attempts) != 0
    @test entropy(number_of_attempts) > 0
    @test pdf(number_of_attempts, 2) < pdf(number_of_attempts, 1)
    @test logpdf(number_of_attempts, 2) < logpdf(number_of_attempts, 1)
    @test cdf(number_of_attempts, 2) > cdf(number_of_attempts, 1)
    @test logcdf(number_of_attempts, 2) > logcdf(number_of_attempts, 1)
    @test logcdf(number_of_attempts, 8.) ≈ log(cdf(number_of_attempts, 8.))
    @test cdf(number_of_attempts, 10) + ccdf(number_of_attempts, 10) ≈ 1
    @test logccdf(number_of_attempts, 10) ≈ log(ccdf(number_of_attempts, 10))
    @test quantile(number_of_attempts, 1.) == Inf
    @test quantile(number_of_attempts, .3) ≈ cquantile(number_of_attempts, .7)
    @test kldivergence(number_of_attempts, number_of_attempts) ≈ 0.

    @test rand(number_of_attempts) >= 1
    @test rand(number_of_attempts) isa Int
    number_of_attempts = NumberOfAttempts(1. - eps(1.))
    @test mean(number_of_attempts) ≈ 1
    samples = [rand(number_of_attempts) for _ in 1:10]
    for sample in samples
        @test sample ≈ 1
    end
    @test rand(NumberOfAttempts(.00001)) > 10
    number_of_attempts = NumberOfAttempts(1)
    @test mean(number_of_attempts) == 1
    samples = [rand(number_of_attempts) for _ in 1:10]
    for sample in samples
        @test sample == 1
    end

end

@testset "Duration" begin
    duration = Duration(.5, 10)
    @test duration isa Duration{Float64, Int}
    @test mean(duration) ≈ 20.
    @test median(duration) > 1
    @test mode(duration) == 10
    @test var(duration) > 0
    @test skewness(duration) != 0
    @test kurtosis(duration) != 0
    @test entropy(duration) > 0
    @test pdf(duration, 20) < pdf(duration, 10)
    @test logpdf(duration, 20) < logpdf(duration, 10)
    @test cdf(duration, 20) > cdf(duration, 10)
    @test logcdf(duration, 20) > logcdf(duration, 10)
    @test logcdf(duration, 8.) ≈ log(cdf(duration, 8.))
    @test cdf(duration, 50) + ccdf(duration, 50) ≈ 1
    @test logccdf(duration, 50) ≈ log(ccdf(duration, 50))
    @test quantile(duration, 1.) == Inf
    @test quantile(duration, .3) ≈ cquantile(duration, .7)
    @test kldivergence(duration, duration) ≈ 0.

    duration = QuantumNetworkRecipes.Duration(1. - eps(1.), 100.)
    samples = [rand(duration) for _ in 1:10]
    for sample in samples
        @test sample ≈ 100.
    end
    duration = QuantumNetworkRecipes.Duration(1, 100.)
    samples = [rand(duration) for _ in 1:10]
    for sample in samples
        @test sample == 100.
    end
end