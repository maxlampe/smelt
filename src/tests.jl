

using Test
using Statistics
include("structs.jl")
include("simulation.jl")
include("tools.jl")


@info "Starting tests"
@testset "Begin smelting..." begin
    @info "Checking for Perkeo defaults"
    @testset "Checking for Perkeo defaults" begin

        dummy = Detector()
        @test dummy.t_meas == 3e-7
        @test dummy.t_dead == 1.5e-6
    end

    @info "Pure rate"
    @testset "Pure rate" begin
        # Check if theo and sim rate are equal for low activity, no bs, tiny integration window
        meas_time = 0.6
        src_rate = 5e2
        mono_src = [Source("mono", src_rate, 400., 0., 0., 0., 0., 0)]

        res, stat = run_sim(mono_src; meas_time=meas_time, with_bs=false, t_meas=1e-8, sim_step=1e-8)
        rate_exp = src_rate * meas_time
        @test length(res) ≈ rate_exp atol=3.0 * sqrt(rate_exp)

    end

    @info "Dead time"
    @testset "Dead time" begin
        # Check if theo and dead time corrected sim rate are equal for high activity, no bs, tiny integration window
        meas_time = 0.8
        src_rate = 5e4
        mono_src = [Source("mono", src_rate, 400., 0., 0., 0., 0., 0)]

        res, stat = run_sim(mono_src; meas_time=meas_time, with_bs=false, t_meas=1e-8, sim_step=1e-8)
        rate_exp = src_rate * meas_time
        @test stat["ev_dt_corr"] ≈ rate_exp atol=3.0 * sqrt(rate_exp)
    end

    @info "Accidental rate"
    @testset "Accidental rate" begin
        # Check if theo and meas accidental rate agree for high activity, no bs
        meas_time = 0.6
        src_rate = 1e5
        mono_src = [Source("mono", src_rate, 400., 0., 0., 0., 0., 0)]
        # use default, checked in first test if Perkeo value
        dummy = Detector()
        meas_wnd = dummy.t_meas

        res, stat = run_sim(mono_src; meas_time=meas_time, with_bs=false, sim_step=1e-8)
        # get number of accidentals
        count_acc = 0
        for e in res
            count_acc += e.n_sum
        end
        rate_acc = (count_acc * stat["ev_dt_corr"] / stat["ev_count"]) / meas_time
        rate_exp = (src_rate)^2 * meas_wnd
        @test rate_acc ≈ rate_exp atol=3.0 * sqrt(rate_exp)
    end

    @info "Backscattering"
    @testset "Backscattering" begin
        # Check if one of the two detectors with bs gets two peaks if they match the empirical fraction in energy and height

        meas_time = 0.3
        src_rate = 5e3
        src_energy = 1000.
        mono_src = [Source("mono", src_rate, src_energy, 0., 0., 0., 0., 0)]
        res, stat = run_sim(mono_src; meas_time=meas_time, with_bs=true, sim_step=1e-8)

        for i = 1:2
            es_mn = []
            es_bs0 = []
            es_bs1 = []
            for e in res
                if e.e_ind[i] > 0.
                    if e.e_ind[i] > 700.
                        push!(es_mn, e.e_ind[i])
                    elseif e.e_ind[i] > 500.
                        push!(es_bs0, e.e_ind[i])
                    else
                        push!(es_bs1, e.e_ind[i])
                    end
                end
            end
            # check probability of backscattering (should be doubled for two detectors)
            tot_bs = length(es_bs0) + length(es_bs1)
            p_exp = (tot_bs) / length(es_mn)
            p_exp_err = sqrt( (sqrt(tot_bs) / length(es_mn))^2 + (sqrt(length(es_mn)) * tot_bs / length(es_mn)^2 )^2 )
            @test p_exp ≈ 2.0 * prob_bs(src_energy) atol=(3.0 * p_exp_err)

            # check fraction of energy being backscattered - should see both parts of the fraction now
            f_bs_exp = mean(es_bs0) / mean(es_mn)
            @test f_bs_exp ≈ e_bs_frac(src_energy) atol=0.01
            f_bs_exp = mean(es_bs1) / mean(es_mn)
            @test f_bs_exp ≈ (1. - e_bs_frac(src_energy)) atol=0.01
        end
    end

    @info "Trigger"
    @testset "Trigger" begin
        # Check if single detector with bs gets two peaks if they match the empirical fraction in energy and height

        meas_time = 0.2
        src_rate = 1e4
        src_energy = 20.
        mono_src = [Source("mono", src_rate, src_energy, 0., 0., 0., 0., 0)]
        # Check if trigger function deviates from Perkeo defaults
        @test 0.5456 ≈ trigger_func(src_energy) atol=0.004

        res, stat = run_sim(mono_src; meas_time=meas_time, with_bs=false)
        p_sum = stat["rate_dt_corr"] / (src_rate)
        @test p_sum ≈ trigger_func(src_energy) atol=0.04
    end
end 
@info "Done!"


"""
missing checks
1) add test to check with no bs, if both trig is a third of all acc with two detectors and half for one
2) expect perfect deltatriggertime (camel) plot for mono energetic source with bs, check time diff with formula
"""