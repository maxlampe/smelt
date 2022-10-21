
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
        @test dummy.t_meas == 2e-7
        @test dummy.t_dead == 1.5e-6
        # FIXME: make t_meas as input for rum_sim None and use Detector default
    end

    @info "Single Detector Sanity"
    @testset "Single Detector Sanity" begin
        # Check that second detector buffer is kept empty during single detector runs even with bs
        meas_time = 0.1
        src_rate = 5e4
        mono_src = Dict("mono" => Source(src_rate, 400., 0.))

        res, stat = run_sim(mono_src; n_det=1, meas_time=meas_time, with_bs=true)
        inv_count = 0
        for e in res
            if e.t_trig[2] > 0 || e.e_ind[2] > 0 || e.t_trig[1] <= 0 || e.e_ind[1] <= 0
                inv_count += 1
            end
        end
        @test inv_count == 0
    end

    @info "Pure rate"
    @testset "Pure rate" begin
        # Check if theo and sim rate are equal for low activity, no bs, tiny integration window
        meas_time = 0.5
        src_rate = 5e2
        mono_src = Dict("mono" => Source(src_rate, 400., 0.))

        for i = 1:2
            res, stat = run_sim(mono_src; n_det=i, meas_time=meas_time, with_bs=false, t_meas=1e-8)
            # half rate for one detector due to it only covering one hemisphere
            rate_exp = src_rate * meas_time * 0.5 * i
            @test length(res) ≈ rate_exp atol=sqrt(4.0 * rate_exp)
        end
    end

    @info "Dead time"
    @testset "Dead time" begin
        # Check if theo and dead time corrected sim rate are equal for high activity, no bs, tiny integration window
        meas_time = 0.5
        src_rate = 5e4
        mono_src = Dict("mono" => Source(src_rate, 400., 0.))
        for i = 1:2
            res, stat = run_sim(mono_src; n_det=i, meas_time=meas_time, with_bs=false, t_meas=1e-8)
            # half rate for one detector due to it only covering one hemisphere
            rate_exp = src_rate * meas_time * 0.5 * i
            @test stat["ev_dt_corr"] ≈ rate_exp atol=sqrt(4.0 * rate_exp)
        end
    end

    @info "Accidental rate"
    @testset "Accidental rate" begin
        # Check if theo and meas accidental rate agree for high activity, no bs
        meas_time = 0.5
        src_rate = 5e4
        mono_src = Dict("mono" => Source(src_rate, 400., 0.))
        # use default, checked in first test if Perkeo value
        dummy = Detector()
        meas_wnd = dummy.t_meas
        for i = 1:2
            res, stat = run_sim(mono_src; n_det=i, meas_time=meas_time, with_bs=false)
            # get number of accidentals
            count_acc = 0
            for e in res
                if e.n_sum > 0
                    count_acc += 1
                end
            end
            rate_acc = (count_acc * stat["ev_dt_corr"] / stat["ev_count"]) / meas_time
            rate_exp = (src_rate * 0.5 * i)^2 * meas_wnd
            @test rate_acc ≈ rate_exp atol=sqrt(4.0 * rate_exp)
        end
    end

    @info "Backscattering, single detector"
    @testset "Backscattering, single detector" begin
        # Check if single detector with bs gets two peaks if they match the empirical fraction in energy and height

        meas_time = 0.3
        src_rate = 5e3
        src_energy = 1000.
        mono_src = Dict("mono" => Source(src_rate, src_energy, 0.))
        res, stat = run_sim(mono_src; n_det=1, meas_time=meas_time, with_bs=true)

        es_mn = []
        es_bs= []
        for e in res
            if e.e_det > 0.
                if e.e_det > 700.
                    push!(es_mn, e.e_det)
                else
                    push!(es_bs, e.e_det)
                end
            end
        end
        # check probability of backscattering
        p_exp = length(es_bs) / length(es_mn)
        p_exp_err = sqrt( (sqrt(length(es_bs)) / length(es_mn))^2 + (sqrt(length(es_mn)) * length(es_bs) / length(es_mn)^2 )^2 )
        @test p_exp ≈ prob_bs(src_energy) atol=(2.0* p_exp_err)
        # check fraction of energye being backscattered
        f_bs_exp = es_bs[1] / es_mn[1]
        @test f_bs_exp ≈ e_bs_frac(src_energy) atol=0.01
    end

    @info "Backscattering, double detector"
    @testset "Backscattering, double detector" begin
        # Check if one of the two detectors with bs gets two peaks if they match the empirical fraction in energy and height

        meas_time = 0.3
        src_rate = 5e3
        src_energy = 1000.
        mono_src = Dict("mono" => Source(src_rate, src_energy, 0.))
        res, stat = run_sim(mono_src; n_det=2, meas_time=meas_time, with_bs=true)

        es_mn = []
        es_bs0 = []
        es_bs1 = []
        for e in res
            if e.e_ind[1] > 0.
                if e.e_ind[1] > 700.
                    push!(es_mn, e.e_ind[1])
                elseif e.e_ind[1] > 500.
                    push!(es_bs0, e.e_ind[1])
                else
                    push!(es_bs1, e.e_ind[1])
                end
            end
        end
        # check probability of backscattering (should be doubled for two detectors)
        tot_bs = length(es_bs0) + length(es_bs1)
        p_exp = (tot_bs) / length(es_mn)
        p_exp_err = sqrt( (sqrt(tot_bs) / length(es_mn))^2 + (sqrt(length(es_mn)) * tot_bs / length(es_mn)^2 )^2 )
        @test p_exp ≈ 2.0 * prob_bs(src_energy) atol=(1.41 * p_exp_err)

        # check fraction of energy being backscattered - should see both parts of the fraction now
        f_bs_exp = mean(es_bs0) / mean(es_mn)
        @test f_bs_exp ≈ e_bs_frac(src_energy) atol=0.01
        f_bs_exp = mean(es_bs1) / mean(es_mn)
        @test f_bs_exp ≈ (1. - e_bs_frac(src_energy)) atol=0.01
    end
end 
@info "Done!"

# ToDo : discuss probs and rates for 2 detector setup
#       currently assuming perkeo-like setup: one detector sees only half of source activity, to detectors see full activity
#       this model makes the expected, detected rate dependend on the number of detectors

# ToDo: add check if bs electron too late for integration

# ToDo: discuss tolerance size (2 sigma) for unit tests

"""
sanity checks
2) expect perfect deltatriggertime (camel) plot for mono energetic source with bs, check time diff with formula
3) expect uniform deltatriggertime plot for mono energetic source without bs
"""