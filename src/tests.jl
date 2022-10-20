
using Test
include("structs.jl")
include("simulation.jl")
# include("plotting.jl")
# include("tools.jl")


@testset "Begin smelting..." begin

    @testset "Checking for Perkeo defaults" begin

        dummy = Detector()

        @test dummy.t_meas == 2e-7
        @test dummy.t_dead == 1.5e-6
        # FIXME: make t_meas as input for rum_sim None and use Detector default
    end

    @testset "Pure rate" begin
        # Check if theo and sim rate are equal for low activity, no bs, tiny integration window
        meas_time = 1.
        bi_rate = 5e2
        mono_src = Dict("mono" => Source(bi_rate, 400., 0.))

        for i = 1:2
            res = run_sim(mono_src; n_det=i, meas_time=meas_time, with_bs=false, t_meas=1e-8)
            @test (length(res) / i) ≈ (bi_rate * meas_time) atol= sqrt(2. * bi_rate * meas_time)
        end
    end

    @testset "Single Detector Sanity" begin
        # Check that second detector buffer is kept empty during single detector runs
        meas_time = 0.1
        bi_rate = 5e2
        mono_src = Dict("mono" => Source(bi_rate, 400., 0.))

        res = run_sim(mono_src; n_det=1, meas_time=meas_time, with_bs=false)

        inv_count = 0
        for e in res
            if e.t_trig[2] > 0 || e.e_ind[2] > 0
                inv_count += 1
            end
        end
        @test inv_count == 0
    end

    @testset "Dead time" begin
        # Check if theo and dead time corrected sim rate are equal for high activity, no bs, tiny integration window
        # meas_time = 0.5
        # for i = 1:2
        #     res = run_sim(mono_src; n_det=i, meas_time=meas_time, with_bs=false, t_meas=1e-8)
        #     @test (length(res) / i) ≈ (bi_rate * meas_time) atol= sqrt(1.0 * bi_rate * meas_time)
        # end
        # FIXME: Dead time calculation is broken
    end
end  




# ToDo: add check if bs electron too late for integration

"""
sanity checks

2) expect perfect deltatriggertime (camel) plot for mono energetic source with bs, check time diff with formula

3) expect uniform deltatriggertime plot for mono energetic source without bs

4) single detector, no bs: accidental rate should be the same as analytics
    10s meas time
    meas: ca. 37 acc evs * (4873 dt corr trig / 4541 obs trig ) = 40 [counts]
    anal: 4600 * 4600 * (2e-7) * 10 = 42 [counts]
    30s meas time
    meas: 192 [counts]
    anal: 127 [counts]

    redo dead time corr
    10s meas time
    meas:  41 [counts] (41 / 45731 / 45419)
    30s meas time
    meas: 158 [counts] (157 / 138256 / 137307)

    ERROR: Dead time corr probably wrong with t_meas and n_det


5) single detector, with bs but short integration window: should see two peaks

    ERROR: see only one and second detector gets tof trigger -> it shouldnt.
"""

