function test_mem_depolar_param_from_link_durations(swap_asap::SwapASAP)
    @testset "mem_depolar_param_from_link_durations" begin
        f = QuantumNetworkRecipes.mem_depolar_param_from_link_durations
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, 5, 1))
            for _ in 1:5]
        nodes = [SimpleNode(NodeWithMemory{DepolarizingMemory}(300)) for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        # wait 10 at each repeater, 20 at end nodes = 80 total
        durations = [10, 20, 30, 20, 10]
        @test f(chain, durations) ≈ exp(-80 / 300)
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, 5, 1/4))
            for _ in 1:5]
        nodes = [SimpleNode(NodeWithMemory{DepolarizingMemory}(Inf)) for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        # even though the links are fully depolarized,
        # that should not be reflected in this fct
        @test f(chain, durations) ≈ 1
    end
end

function test_perfect(swap_asap::SwapASAP)
    @testset "perfect" begin
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, 5, 1))
            for _ in 1:5]
        nodes = [SimpleNode(PerfectNode()) for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        @test edge_length(chain) == 50
        for method in [Analytical(), Sampling(10)]
            @test QuantumNetworkRecipes.werner_parameter(chain, method)[1] ≈ 1.
            @test QuantumNetworkRecipes.werner_parameter(chain, method)[2] ≈ 0.
            @test entangled_state_fidelity(chain, method)[1] ≈ 1.
            @test entangled_state_fidelity(chain, method)[2] ≈ 0.
            @test qber_x(chain, method)[1] ≈ 0.
            @test qber_x(chain, method)[2] ≈ 0.
            @test qber_y(chain, method)[1] ≈ 0.
            @test qber_y(chain, method)[2] ≈ 0.
            @test qber_z(chain, method)[1] ≈ 0.
            @test qber_z(chain, method)[2] ≈ 0.
        end
        @test generation_duration(chain, Sampling(10))[1] ≈ 5.
        @test generation_duration(chain, Sampling(10))[2] ≈ 0.
        nodes = [SimpleNode(NodeWithMemory{DepolarizingMemory}(Inf)) for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        @test QuantumNetworkRecipes.werner_parameter(chain, Sampling(10))[1] ≈ 1.
        @test QuantumNetworkRecipes.werner_parameter(chain, Sampling(10))[2] ≈ 0.
        @test generation_duration(chain, Sampling(10))[1] ≈ 5.
        @test generation_duration(chain, Sampling(10))[2] ≈ 0.
        edges = [DirectTransmissionOverFiber{WernerState}(FiberSegment(10, 0, 1, 1), 1)
            for _ in 1:5]
        nodes = [SimpleNode(NodeWithPhotonSource{DepolarizingMemory}(Inf, 1, 10, 100))
            for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        @test generation_duration(chain, Sampling(10))[1] ≈ 10
        @test generation_duration(chain, Sampling(10))[2] ≈ 0.
        @test QuantumNetworkRecipes.werner_parameter(chain, Sampling(10))[1] ≈ 1
        @test QuantumNetworkRecipes.werner_parameter(chain, Sampling(10))[2] ≈ 0.
    end
end

function test_max_mixed(swap_asap::SwapASAP)
    @testset "max mixed" begin
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, 1, 1/4))
            for _ in 1:5]
        nodes = [SimpleNode(PerfectNode()) for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        @test entangled_state_fidelity(chain, Analytical())[1] ≈ 1 / 4
        @test entangled_state_fidelity(chain, Analytical())[2] ≈ 0.
        @test entangled_state_fidelity(chain, Sampling(5))[1] ≈ 1 / 4
        @test entangled_state_fidelity(chain, Sampling(5))[2] ≈ 0.
        nodes = [SimpleNode(NodeWithMemory{DepolarizingMemory}(Inf)) for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        @test entangled_state_fidelity(chain, Sampling(5))[1] ≈ 1 / 4
        @test entangled_state_fidelity(chain, Sampling(5))[2] ≈ 0.
    end
end

function test_decoherence(swap_asap::SwapASAP)
    @testset "decoherence" begin
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, .1, 1, 1))
            for _ in 1:5]
        nodes = [SimpleNode(NodeWithMemory{DepolarizingMemory}(50)) for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        @test entangled_state_fidelity(chain, Sampling(10))[1] < .8
        @test entangled_state_fidelity(chain, Sampling(10))[2] > 0.001
        # below is a way of testing variance before error estimates were introduced to
        # this package, it's a bit outdated but doesn't hurt to keep around
        samples = [generation_duration(chain, Sampling(1))[1] for _ in 1:10]
        @test std(samples) > 0
        samples = [entangled_state_fidelity(chain, Sampling(1))[1] for _ in 1:10]
        @test std(samples) > 0
        @test generation_duration(chain, Sampling(10))[1] > 1
        @test generation_duration(chain, Sampling(10))[2] > 0.001
    end
end

function test_stochastic_ad(swap_asap::SwapASAP)
    @testset "compatibility with stochastic AD" begin
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(
            stochastic_triple(10.), stochastic_triple(.1), stochastic_triple(5.),
            stochastic_triple(.96))) for _ in 1:5]
        nodes = [SimpleNode(NodeWithMemory{DepolarizingMemory}(stochastic_triple(50.)))
            for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        dur, dur_error = generation_duration(chain, Sampling(5))
        fid, fid_error = entangled_state_fidelity(chain, Sampling(5))
        for x in [dur, dur_error, fid, fid_error]
            @test x isa StochasticAD.StochasticTriple
        end
    end
end

@testset "restart per state" begin
    swap_asap = SwapASAPRestartPerState()

    test_mem_depolar_param_from_link_durations(swap_asap)
    test_perfect(swap_asap)
    test_max_mixed(swap_asap)
    test_decoherence(swap_asap)
    test_stochastic_ad(swap_asap)

    @testset "sample link durations" begin
        f = QuantumNetworkRecipes.sample_link_durations
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, 8, 1))
            for _ in 1:5]
        nodes = [SimpleNode(PerfectNode()) for _ in 1:6]
        chain = ChainRecipe(edges, nodes, swap_asap)
        samples = f(chain, 3)
        @test size(samples) == (5, 3)
        @test all(samples .== 8)
    end
end

@testset "without cutoff" begin

    @testset "find free time" begin
        f = QuantumNetworkRecipes._find_free_time
        x = [3, 8, 2, 4, 1]
        @test f(x, 1) == f(x, 2) == f(x, 3) == f(x, 5) == 8
        @test f(x, 4) == 4
        x = [4, 1, 2, 3, 1, 3]
        @test f(x, 1) == f(x, 2) == f(x, 6) == 4
        @test f(x, 3) == f(x, 4) == f(x, 5) == 3

        @testset "Stochastic AD compatibility" begin
            # link in the middle
            x = [10., 6., stoch_trip(4., 0., 4., 1.), 6., 10.]
            ready_time = f(x, 3)
            # zeroth order: link 3 swaps both neighbors at time 6
            # first order: link 3 is the limiting factor and only free at time 8
            @test value(ready_time) == 6.
            @test perturbation_is(ready_time, (2., 1.))

            # link at the edge
            x = [stoch_trip(10., 1., 10., 2.), 10., 15., 10.]
            ready_time = f(x, 1)
            @test value(ready_time) == 15.
            @test perturbation_is(ready_time, (5., 2.))

            # stochastic index
            x = [5., 5., 5, 10., 10., 10.]
            ready_time = f(x, stoch_trip(2, 3, 3, 1.))
            # ready at 5 for index 2, ready at 10 for index 5
            @test value(ready_time) == 5.
            @test perturbation_is(ready_time, (5., 1.))

            # stochastic index that makes a difference between center or edge
            x = [5., 10., 5., 5., 5.]
            ready_time = f(x, stoch_trip(4, 1, 1, 1.))
            # 4 is center, 5 is edge and will have to wait until time 10
            @test value(ready_time) == 5.
            @test perturbation_is(ready_time, (5., 1.))
        end
    end

    swap_asap = SwapASAPWithoutCutoff()

    test_mem_depolar_param_from_link_durations(swap_asap)
    test_perfect(swap_asap)
    test_max_mixed(swap_asap)
    test_decoherence(swap_asap)
    test_stochastic_ad(swap_asap)

    @testset "sample link durations" begin
        f = QuantumNetworkRecipes.sample_link_durations
        link_times = [5., 4., 3., 2., 3., 4.]  # deterministic time to generate links
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, i, 1))
            for i in link_times]
        nodes = [SimpleNode(PerfectNode()) for _ in 1:7]
        chain = ChainRecipe(edges, nodes, swap_asap)
        samples = f(chain, 4)
        @test size(samples) == (6, 4)
        @test samples[:, 1] == link_times
        # link 1 and 6 can only restart at time 5, after the end-to-end link is finished,
        # because they are end nodes
        # link 2 can only restart at time 5 because the first repeater is busy until then
        # link 3 and 5 can restart at time 4, link 4 at time 3, giving them a head start
        # of 1 and 2 time units respectively; effectively, link 3 and 5 finish at time
        # 4 + 3 = 7, and link 4 finishes at time 3 + 2 = 5. Relative to the starting time of
        # the second round, i.e. counted after delivering the first entangled state, these
        # give completion times of 2 and 0 respectively
        @test samples[:, 2] == [5., 4., 2., 0., 2., 4.]
        # in the third round, all links except the 4th one become free at the same time as
        # last time, and hence should give the same result again
        # but link 4 can restart at time 2, making it finish at time 4 after completion of
        # the first entangled state, and hence time -1 after the completion of the second
        # entangled state (i.e., one time unit before it was completed)
        @test samples[:, 3] == [5., 4., 2., -1., 2., 4.]
    end
end

@testset "local cutoff" begin
    @testset "find first cutoff link" begin
        f = QuantumNetworkRecipes._find_first_cutoff_link
        sample = [6, 3, 10, 5, 4, 3]
        cutoff_times = [2 for _ in 1:5]
        index, cutoff_moment = f(sample, cutoff_times)
        @test index == 2
        @test cutoff_moment == 5  # 3 + 2
        sample[2] = 10  # this one no longer triggers cutoff
        index, cutoff_moment = f(sample, cutoff_times)
        @test index == 4
        @test cutoff_moment == 7  # 5 + 2
        sample = [2, 7, 5, 3, 3]
        cutoff_times = [3 for _ in 1:4]
        index, cutoff_moment = f(sample, cutoff_times)
        @test index == 1
        @test cutoff_moment == 5  # 2 + 3
        sample[1] = 5  # now no cutoff is triggered anymore
        index, cutoff_moment = f(sample, cutoff_times)
        @test iszero(index)
        @test cutoff_moment == false
        sample = [10, 2, 3, 10]
        cutoff_times = [7, 2, 3]
        index, cutoff_moment = f(sample, cutoff_times)
        @test index == 3  # not the oldest link that is discarded, but it is the first
        @test cutoff_moment == 6  # 3 + 3

        @testset "StochasticAD compatibility" begin
            st = stoch_trip(10., 1., 5., 2.,)
            @test st isa StochasticAD.StochasticTriple
            sample = [0, st, 10]
            cutoff_times = [12., 12.]
            # if st is 10, there is no cutoff and hence the result is (0, false)
            # if st is 15, link 1 is discarded at time 12, and hence the result is (1, 12)
            # hence, we expect (0 + (1 with probability 2.0ϵ),
            # 0.0 + (12. with probability 2.0ϵ))
            # (Note: false = 0.0, allowing treatment by StochasticAD.)
            (index, cutoff_moment) = f(sample, cutoff_times)
            @test value(index) == 0
            @test value(cutoff_moment) == false
            @test perturbation_is(index, (1, 2.))
            @test perturbation_is(cutoff_moment, (12., 2.))

            sample = [100., 12, st, 100]
            cutoff_times = [12., 12., 12.]
            # if st is 10, link 3 is cut off at time 10 + 12 = 22
            # if st is 15, link 2 is cut off at time 12 + 12 = 24
            (index, cutoff_moment) = f(sample, cutoff_times)
            @test value(index) == 3
            @test value(cutoff_moment) == 22
            @test perturbation_is(index, (-1, 2.))
            @test perturbation_is(cutoff_moment, (2., 2.))

            sample = [st, st]
            cutoff_times = [1.]
            # st should always be the same value, and hence there is never a cutoff
            # if they are accidentally treated as independent perturbations though,
            # there will be cutoffs in perturbations (when one is 15 and the other is 10)
            (index, cutoff_moment) = f(sample, cutoff_times)
            @test iszero(index)
            @test iszero(cutoff_moment)

            # if the two instances of the sample are not independent,
            # then there are perturbation terms where there is a cutoff
            counter1 = counter2 = 0
            for _ in 1:10
                st1 = stoch_trip(10., 1., 5., 2.,)
                st2 = stoch_trip(10., 1., 5., 2.,)
                sample = [st1, st2]
                (index, cutoff_moment) = f(sample, cutoff_times)
                @test iszero(value(index))
                @test iszero(value(cutoff_moment))
                # whichever way the perturbation falls, the cutoff is always at time 11
                # both triples contribute 2ϵ probability, so 4ϵ in total
                @test perturbation_is(cutoff_moment, (11., 4.))

                # the index will depend on which way the perturbation falls
                # instead of accounting for all branches, it seems StochasticAD
                # randomly chooses a branch; we should see both branches (cutoff index
                # 1 or 2) with equal probability
                # while each of the perturbations occurs with probability 2ϵ,
                # the "collapsing" of the branch means that we will see probabilities 4ϵ
                # (4ϵ * .5 to account for the fact that each appears only half of the times
                # will give 2ϵ again)
                perturbation_is(index, (1, 4.)) && (counter1 += 1)
                perturbation_is(index, (2, 4.)) && (counter2 += 1)
            end
            @test counter1 + counter2 == 10
            @test counter1 > 0
            @test counter2 > 0

            sample = [12., 0., 10.]
            cutoff_times = [st, st]
            # cutoff at link 2 if cutoff is 10, not if it is 15
            (index, cutoff_moment) = f(sample, cutoff_times)
            @test value(index) == 2
            @test value(cutoff_moment) == 10
            @test perturbation_is(index, (-2, 2.))
            @test perturbation_is(cutoff_moment, (-10., 2.))

            # both cutoff time and sample need to be in perturbation to get cutoff
            # this is a second-order event and should hence not happen
            sample = [10., stoch_trip(9., 0., -4., 1.), 10.]
            cutoff_times = [stoch_trip(6., 0., -3., 1.) for _ in 1:2]
            (index, cutoff_moment) = f(sample, cutoff_times)
            @assert iszero(index)
            @assert iszero(cutoff_moment)
        end
    end

    @testset "find swapped links" begin
        _to_set(x) = Set([i for (i, val) in Iterators.enumerate(x) if val])
        f = _to_set ∘ QuantumNetworkRecipes._find_swapped_links
        sample = [3, 2, 6, 4, 10]
        @test f(sample, 2, 1) == Set()
        @test f(sample, 2, 2) == Set([2])
        @test f(sample, 2, 3) == Set([1, 2])
        @test f(sample, 2, 6) == f(sample, 2, 7) == Set(1:4)
        @test f(sample, 2, 10) == Set(1:5)
        @test f(sample, 4, 3) == Set()
        @test f(sample, 4, 4) == Set([4])
        @test f(sample, 4, 7) == Set(1:4)
        @test f(sample, 5, 10) == Set(1:5)

        @testset "StochasticAD compatibility" begin
            f = QuantumNetworkRecipes._find_swapped_links

            # stochastic sample
            sample = [10., stoch_trip(15., 0., 10., 1.), 10.]
            @test f(sample, 1, 5.) == [false, false, false]
            @test f(sample, 1, 30.) == [true, true, true]
            # at time 20, swap happened only in the perturbation
            f(sample, 1, 20.)

            # stochastic index
            sample = [10, 100, 10, 100]
            index = stoch_trip(1, 0, 2, 1.)
            @test f(sample, index, 5.) == [false, false, false, false]
            swapped_links = f(sample, index, 15.)
            @test value.(swapped_links) == [true, false, false, false]
            @test all(perturbation_is.(swapped_links,
                [(-1, 1.), (0, 1.), (1, 1.), (0, 1.)]))

            # stochastic time
            sample = [10., 20.]
            time = stoch_trip(5., 0., 10., 3.)
            swapped_links = f(sample, 1, time)
            @test value.(swapped_links) == [false, false]
            @test all(perturbation_is.(swapped_links, [(1., 3.), (0., 3.)]))
            time = stoch_trip(15., 0., 10., 3.)
            swapped_links = f(sample, 1, time)
            @test value.(swapped_links) == [true, false]
            @test all(perturbation_is.(swapped_links, [(0., 3.), (1., 3.)]))

            # stochastic sample and index
            sample = [stoch_trip(10., 0., 10., 3.), 15.]
            index = stoch_trip(1, 0, 1, 10.)
            swapped_links = f(sample, index, 12)
            # both perturbations result in no swaps at all, with total 10 + 3 = 13ϵ
            @test value.(swapped_links) == [true, false]
            @test all(perturbation_is.(swapped_links, [(-1., 13.), (0., 13)]))

            # stochastic sample and time
            sample = [stoch_trip(10., 0., -5., 3.), 100.]
            time = stoch_trip(8., 0., 10., 10.)
            swapped_links = f(sample, 1, time)
            # either perturbation makes the first entry from false to true
            @test value.(swapped_links) == [false, false]
            @test all(perturbation_is.(swapped_links, [(1., 13.), (0., 13.)]))

            # stochastic index and time
            sample = [10., 20.]
            index = stoch_trip(1, 0, 1, 7.)
            time = stoch_trip(15., 1., -10., 2.)
            swapped_links = f(sample, index, time)
            # either perturbation results in no swapped links
            @test value.(swapped_links) == [true, false]
            @test all(perturbation_is.(swapped_links, [(-1., 9.), (0., 9.)]))

            # stochastic sample, index, and time
            sample = [stoch_trip(10., 1., 100., 1.), 100.]
            index = stoch_trip(1, 0, 1, 7.)
            time = stoch_trip(15., 1., -10., 4.)
            swapped_links = f(sample, index, time)
            # all perturbations result in no swapped links
            @test value.(swapped_links) == [true, false]
            @test all(perturbation_is.(swapped_links, [(-1., 12.), (0., 12.)]))
        end
    end

    @testset "resample for earliest cutoff" begin
        f = QuantumNetworkRecipes._resample_for_earliest_cutoff
        sample = [10, 1, 10]
        cutoff_times = [2, 2] 
        random_variables = [Dirac(i) for i in [100, 4, 200]]
        resampled, sample = f(sample, random_variables, cutoff_times)
        # at time 3, link 2 will trigger a cutoff and will have to regenerate
        # regeneration takes 4 time units, meaning the new link is finished at time
        # 1 + 2 + 4 = 7
        @test sample == [10, 7, 10]
        resampled, sample = f(sample, random_variables, cutoff_times)
        @test sample == [10, 13, 10]
        sample = [20, 9, 13, 14, 20]
        cutoff_times = [8 for _ in 1:4]
        random_variables = [Dirac(100 * i) for i in 1:5]
        resampled, sample = f(sample, random_variables, cutoff_times)
        # at time 17, link 2 triggers a cutoff; at that time, it will have been swapped
        # already with link 3 and link 4, which means they are all discarded
        # link 2 starts regenerating at time 17, and finishes at time 217
        # link 4 also starts regenerating at the cutoff and finishes at the time 417
        # link 3 on the other hand becomes available after swapping on both sides at time 14
        # and hence it already starts regenerating at that time, finishing at time 314.
        @test sample == [20, 217, 314, 417, 20]

        @testset "StochasticAD compatibility" begin
            # stochastic cutoff time
            sample = [10., 4., 10.]
            random_variables = [Dirac(10.), Dirac(4.), Dirac(10.)]
            st = stoch_trip(2., 0., 10., 2.)
            cutoff_times = [st for _ in 1:2]
            resampled, sample = f(sample, random_variables, cutoff_times)
            @test convert(Bool, value(resampled)) == true
            @test value.(sample) == [10., 10., 10.]
            @test perturbation_is(resampled, (-1, 2.))
            @test all(perturbation_is.(sample, [(0, 2.), (-6., 2.), (0, 2.)]))

            # stochastic sample, cutoff at 0th order, not at 1st order
            sample = [10., st, 10.]
            cutoff_times = [4. for _ in 1:2]
            # if finished at time 2, cutoff at 6 and regenerate at 10
            # if finished at time 12, no cutoff
            resampled, sample = f(sample, random_variables, cutoff_times)
            @test value(resampled) == true
            @test perturbation_is(resampled, (-1, 2.))
            @test value.(sample) == [10., 10., 10.]
            @test all(perturbation_is.(sample, [(0, 2.), (2., 2.), (0, 2.)]))

            # stochastic sample, cutoff at 1st order, not at 0th order
            sample = [10., 12 - st, 10.]
            cutoff_times = [2. for _ in 1:2]
            # 0th order: all finish at time 10
            # 1st order: link 2 ready at time 0, cutoff at time 2, regen at time 6
            resampled, sample = f(sample, random_variables, cutoff_times)
            @test value(resampled) == false
            @test perturbation_is(resampled, (1, 2.))
            @test value.(sample) == [10., 10., 10.]
            @test all(perturbation_is.(sample, [(0, 2.), (-4., 2.), (0, 2.)]))

            # stochastic, coupled sample
            sample = [st, st, st]
            cutoff_times = [1., 1.]
            resampled, sample = f(sample, random_variables, cutoff_times)
            # the samples are coupled so should never be a cutoff
            @test iszero(resampled)
            @test sample == [st, st, st]

            # stochastic sample and cutoff time
            sample = [10., 8 + st, 10.]
            cutoff_times = [st for _ in 1:2]
            # 0th order: finish at same time with small cutoff
            # 1st order: finish late with large cutoff
            # both cases, cutoff is not triggered
            resampled, sample = f(sample, random_variables, cutoff_times)
            @test iszero(resampled)
            @test value.(sample) == [10., 10., 10.]
            @test all(perturbation_is.(sample, [(0, 2.), (10., 2.), (0, 2.)]))

            # resampled value is stochastic triple
            random_variables = [Dirac(10.), Dirac(st), Dirac(10.)]
            sample = [10., 7., 10.]
            cutoff_times = [2. for _ in 1:2]
            # cutoff is triggered at time 9
            # 0th order: regenerates at time 11
            # 1st order: regenerates at time 21
            resampled, sample = f(sample, random_variables, cutoff_times)
            @test resampled == true
            @test value.(sample) == [10., 11., 10.]
            @test all(perturbation_is.(sample, [(0., 0.), (10., 2.), (0., 0.)]))
            # in the first case, resampling again has no effect
            # in the second case, if we resample two more times,
            # the other links are discarded at 12 and regenerated at 22
            resampled, sample = f(sample, random_variables, cutoff_times)
            @test value(resampled) == false
            @test perturbation_is(resampled, (1, 2.))
            resampled, sample = f(sample, random_variables, cutoff_times)
            @test value(resampled) == false
            @test perturbation_is(resampled, (1, 2.))
            @test value.(sample) == [10., 11., 10.]
            @test all(perturbation_is.(sample, [(12., 2.), (10., 2.), (12., 2.)]))

        end
    end

    @testset "resample until no cutoff" begin
        f = QuantumNetworkRecipes._resample_until_no_cutoff
        sample = [10, 1, 10]
        cutoff_times = [2, 2] 
        random_variables = [Dirac(3) for _ in 1:3]
        sample = f(sample, random_variables, cutoff_times)
        # middle link will cutoff at 3, regenerate at 6, cut off at 8 and regenrate at 11
        @test sample == [10, 11, 10]
        sample = [10, 1, 10]
        cutoff_times = [8, 9]
        random_variables = [Dirac(5), Dirac(11), Dirac(5)]
        sample = f(sample, random_variables, cutoff_times)
        # middle link will cutoff at 9 and regenerate at 20
        # then, link 1 will cut off at 18 and regenerate at 23
        # finally, link 3 will cut off at time 19 and then regenerate at 24
        @test sample == [23, 20, 24]
        # the next test is probabilistic as completion times are sampled,
        # but the resampling function should be such that no matter what we always end
        # with a sample where no qubit is stored longer than the cutoff time
        for _ in 1:10
            sample = [4, 30, 18, 3, 20]
            cutoff_time = 8
            cutoff_times = [cutoff_time for _ in 1:4]
            random_variables = [Geometric(0.1) for _ in 1:5]
            sample = f(sample, random_variables, cutoff_times)
            diffs = [abs(sample[i + 1] - sample[i]) for i = 1:4]
            @test all(diffs .<= cutoff_time)
            @test iszero(QuantumNetworkRecipes._find_first_cutoff_link(sample,
                cutoff_times)[1])
        end

        @testset "Stochastic AD compatibility" begin
            st = stoch_trip(2., 0., 10., 2.)

            # only resample once (copied from "resample for earliest cutoff") 
            # cutoff at 0th order, not at 1st
            sample = [10., st, 10.]
            random_variables = [Dirac(10.), Dirac(4.), Dirac(10.)]
            cutoff_times = [4. for _ in 1:2]
            sample = f(sample, random_variables, cutoff_times)
            @test value.(sample) == [10., 10., 10.]
            @test all(perturbation_is.(sample, [(0, 2.), (2., 2.), (0, 2.)]))

            # resample multiple times (copied from "resample for earliest cutoff")
            random_variables = [Dirac(10.), Dirac(st), Dirac(10.)]
            sample = [10., 7., 10.]
            cutoff_times = [2. for _ in 1:2]
            sample = f(sample, random_variables, cutoff_times)
            @test value.(sample) == [10., 11., 10.]
            @test all(perturbation_is.(sample, [(12., 2.), (10., 2.), (12., 2.)]))
        end
    end

    @testset "sample link durations" begin

        f = QuantumNetworkRecipes.sample_link_durations
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, t, 1))
            for t in [10, 4, 10]]
        nodes = [SimpleNode(PerfectNode()) for _ in 1:4]
        chain = ChainRecipe(edges, nodes, SwapASAPWithCutoff(LocalCutoffProtocol(2.)))
        # start with [10, 4, 10], then middle node cutoff at 6, regenerates at 10
        # this gives [10, 10, 10] which means we get a clean slate for the next one
        # (no residual entangelement)
        samples = f(chain, 5)
        @test size(samples) == (3, 5)
        @test all(samples .== 10)

        # now a test with residual entanglement and no cutoffs
        # this test is essentially copied from the one for swap-asap without cutoffs
        link_times = [5., 4., 3., 2., 3., 4.]
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, t, 1))
            for t in link_times]
        nodes = [SimpleNode(PerfectNode()) for _ in 1:7]
        # default cutoff is infinity, i.e. no cutoff
        chain = ChainRecipe(edges, nodes, SwapASAPWithCutoff(LocalCutoffProtocol()))
        samples = f(chain, 4)
        @test size(samples) == (6, 4)
        @test samples[:, 1] == link_times
        @test samples[:, 2] == [5., 4., 2., 0., 2., 4.]
        @test samples[:, 3] == [5., 4., 2., -1., 2., 4.]

        # now we test residual entanglement in the presence of cutoffs
        link_times = [8., 6., 2., 5., 2.]
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, t, 1))
            for t in link_times]
        nodes = [SimpleNode(PerfectNode()) for _ in 1:6]
        chain = ChainRecipe(edges, nodes, SwapASAPWithCutoff(LocalCutoffProtocol(2.)))
        # links 1 and 2 don't trigger cutoffs and swap at time 8
        # link 3 is discarded at 4, regenerated at 6 and swapped at 6
        # link 5 is also discarded at 4 and regenerated at 6
        # link 3, 4 and 5 are swapped at time 6
        # this gives links 3 and 4 a 2 headstart next round
        # note that link 5 does not get such a head start, as end nodes hold on to their
        # qubits until the end-to-end link is finished
        samples = f(chain, 3)
        @test samples[:, 1] == [8., 6., 6., 5., 6.]
        # links 3, 4 are ready at times 0, 3 due to headstarts
        # link 5 is ready at time 2
        # links 3 triggers cutoffs at 2 and then regenerates at 4
        # links 3, 4, 5 are then ready at times 4, 3, 2;
        # link 4 and 5 swap at time 3, links 3 and 4 at time 4
        # then at time 6, link 2 and 3 are swapped and at time 8, links 1 and 2,
        # hence no further cutoffs are triggered
        @test samples[:, 2] == [8., 6., 4., 3., 2.]
        # now, link 3 has a headstart of 2 and link 4 of 4
        # therefore they regenerate their links in the next round at times 0 and 1
        # giving [8, 6, 0, 1, 2]
        # at time 1, links 3 and 4 swap, and at time 2, links 4 and 5 swap
        # however, at time 2, a cutoff is triggered in link 3 (or more accurately,
        # by the second repeater node), which means that links 3, 4 and 5 are all discarded
        # at this time, links 3 and 5 start anew, finishing at time 4
        # link 4 also restarts and finishes at time 7
        # now we have [8, 6, 4, 7, 4]
        # then at time 6, links 2 and 3 swap and are immediately discarded,
        # while at the same time link 5 is discarded, giving [8, 12, 8, 7, 6]
        # now links 1, 3, 4, 5 are all discarded at time 10 (4 and 5 because they swapped
        # with link 3) and have to regenerate, at this time, except for link 4 which was
        # already free at time 8 giving [18, 12, 12, 13, 12]
        # at time 14, links 2-5 are discarded; links 2 and 5 restart at this time,
        # but links 3 and 4 restart at time 13, giving [18, 20, 15, 18, 16]
        # then link 3 is discarded at 17 and regenerated at 19
        # we then, finally, have an end-to-end pair!
        @test samples[:, 3] == [18., 20., 19., 18., 16.]

        # test non-uniform cutoff times
        link_times = [8., 3., 8., 3., 8., 3., 8., 3., 8.]
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, t, 1))
            for t in link_times]
        nodes = [NodeWithCutoff(PerfectNode(), t) for
            t in [Inf, 2., 2., 3., 3., 2., 3., 3., 2., Inf]]
        chain = ChainRecipe(edges, nodes, SwapASAPWithCutoff(LocalCutoffProtocol()))
        samples = f(chain, 1)
        # the first 3 has a cutoff time of 2 at both its nodes,
        # meaning it will be cut off at time 5 and regenerate at time 8
        # the second 3 has a cutoff time of 3 at both its nodes,
        # meaning it will be cut off at time 6 and regenerate at time 9
        # the third 3 has a cutoff time of 2 on one side and 3 on the other;
        # the shorter will be the limiting factor and hence it will be like the first
        # same for the final 3, but with the cutoff sides reversed
        @test samples == [[8., 8., 8., 9., 8., 8., 8., 8., 8.];;]

        # check that end nodes never trigger cutoffs
        link_times = [10., 20., 30., 40., 50., 40., 30., 20., 10.]
        edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, t, 1))
            for t in link_times]
        nodes = [SimpleNode(PerfectNode()) for _ in 1:10]
        chain = ChainRecipe(edges, nodes, SwapASAPWithCutoff(LocalCutoffProtocol(15.)))
        # each of the repeaters can swap after 10 time units and hence never triggers
        # a cutoff, the end nodes need to store for 40 time units
        # they shouldn't trigger a cutoff either though because end nodes are never
        # supposed to do so
        sample = f(chain, 1)
        @test sample == [link_times;;]

        @testset "Stochastic AD compatibility" begin

            # stochastic cutoff time
            edges = [SimpleHeraldedRecipe{WernerState}(SimpleHeraldedPhysical(10, 1, t, 1))
                for t in [10, 4, 10]]
            nodes = [SimpleNode(PerfectNode()) for _ in 1:4]
            chain = ChainRecipe(edges, nodes,
                SwapASAPWithCutoff(LocalCutoffProtocol(stoch_trip(2., 0., 10., 3.))))
            # for cutoff time of 2, link 2 restarts at time 6 and finishes at time 10
            # for cutoff time of 12, link 2 does not need to restart and finishes at 4
            sample = f(chain, 1)
            @test value.(sample) == [[10, 10, 10];;]
            @test all(perturbation_is.(sample, [[(0, 3.), (-6., 3.), (0, 3.)];;]))

            # TODO more extensive testing
        end
    end

    swap_asap = SwapASAPWithCutoff(LocalCutoffProtocol(10.))

    test_mem_depolar_param_from_link_durations(swap_asap)
    test_perfect(swap_asap)
    test_max_mixed(swap_asap)
    test_decoherence(swap_asap)
    # test_stochastic_ad(SwapASAPWithCutoff(LocalCutoffProtocol(stochastic_triple(10.))))

end