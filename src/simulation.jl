

using Printf
using ProgressBars
include("structs.jl")
include("tools.jl")


# Only works for dtime * SIM_STEP << 1
function create_energy(srcs::Vector{Vector{Float64}}, dtime::Float64)::Vector{Vector{Float64}}
    energy = []
    for src in srcs
        # rate effectively halved for one detector, as one only covers one hemisphere
        if (src[1] * dtime) >= rand()
            det = Int64(round(rand()) + 1)
            ev = [0., 0.]
            ev[det] = clamp(randn() * src[3] + src[2], 0., Inf)
            append!(energy, [ev])
            
            # create correlated event
            if src[4] > 0.
                if src[4] >= rand()
                    ev = [0., 0.]
                    # assign correlate event
                    if src[7] == -1
                        # opposite detector
                        det_c = (det % 2) + 1
                    elseif src[7] == 1
                        # same detector
                        det_c = det 
                    else
                        # random
                        det_c = Int64(round(rand()) + 1)
                    end
                    ev[det_c] = clamp(randn() * src[6] + src[5], 0., Inf)
                    append!(energy, [ev])
                end
            end
        end
    end
    return energy
end


# Approximate inverse CDF sampling (smaller angles more likely than they actually are)
function create_ang()::Float64
    # return rand() * 47.
    return (1.0 - (1.0 / (1.0 - 0.9 * rand()))) * (47. / 9.) + 47.
end


function check_trigger(en::Float64)

    if en > 60.0
        return true
    elseif en > 8.0
        p_trig = trigger_func(en)
        if p_trig >= rand()
            return true
        end
    end
    return false
end


function run_sim(
    srcs::Vector{Source};
    meas_time::Float64 = 1e0,
    sim_step::Float64 = 5e-9,
    with_bs::Bool = false,
    with_ang::Bool = false,
    with_bs_var::Bool = false,
    t_meas::Union{Float64, Nothing} = nothing,
    verbose::Bool = false,
)::Tuple{Vector{Event}, Dict}

    if !isnothing(t_meas)
        det_main = Detector(t_meas=t_meas)
    else
        det_main = Detector()
    end

    # calculate theoretical, pure rate from source
    theo_rate = 0.
    for src in srcs
        theo_rate += src.rate
    end

    #convert to vector for faster evaluation
    src_vec = Vector{Vector{Float64}}(undef, length(srcs))
    for (i, src) in enumerate(srcs)
        src_vec[i] = [src.rate, src.energy, src.width, src.p_cor, src.e_cor, src.w_cor, src.d_cor]
    end

    n_time_steps = meas_time/sim_step
    # Should be large enough buffer
    data = Vector{Event}(undef, floor(Int, theo_rate * meas_time * 2.0)) # n_time_steps * 0.01
    data_ind = 1
    curr_ev = Event()
    t_last_ev = 0.

    if verbose
        print("Smelting... - ")
    end

    # main loop
    if verbose
        iter = ProgressBar(1:n_time_steps)
    else
        iter = (1:n_time_steps)
    end
    for t in iter

        if verbose
            if t % (0.25 * n_time_steps) == 0
                print(100. * t/n_time_steps, "% - ")
            end
        end

        t_now = t * sim_step
        dt_last_ev = t_now - t_last_ev

        # check if detector is in dead time
        if det_main.is_ready == false
            if dt_last_ev >= det_main.t_dead
                det_main.is_ready = true
            else
                continue
            end
        end

        # Create events in this time step for both detectors [n_ev x 2]
        events = create_energy(src_vec, sim_step)
        n_ev = lastindex(events)
        if n_ev > 0
            
            # assign each ev an angle (at detector) (ANGEL FLAG) [n_ev x 2]
            if with_ang
                ang_now = [[create_ang(), create_ang()] for i=1:n_ev]
            else
                ang_now = [[35.7, 35.7] for i=1:n_ev]
            end

            # calc tof (to detector ~3.2m) for timing [n_ev x 2]
            dist = [3.7, 2.7]
            tofs_ini = [[t_tof(events[j][i], ang_now[j][i], dist[i]) for i=1:2] for j=1:n_ev ]

            # do bs: split energies (or create new events?), calculate new tofs (from detector to detector) (BS FLAG, VAR FLAG)
            if with_bs
                events_bs = []
                tof_bs = []

                # calc for each event, probability of p_bs [n_ev x 2]
                p_bs = [[prob_bs(events[j][i], ang_now[j][i]) for i=1:2] for j=1:n_ev ]
                # p = rand(n_ev, 2)
                p = [rand(2) for i=1:n_ev]

                for n = 1:n_ev
                    for d = 1:2
                        if events[n][d] <= 0.
                            continue
                        end

                        if p_bs[n][d] >= p[n][d]
                            frac_bs = e_bs_frac(events[n][d], ang_now[n][d])
                            if with_bs_var
                                norm = frac_bs * 2.6 + 1.6
                                frac_bs = clamp(frac_bs + randn() * frac_bs / norm, 0., 1.)
                            end
                            
                            energy_bs = (1.0 - frac_bs) * events[n][d]
                            events[n][d] = events[n][d] - energy_bs

                            # create new event for bs electron
                            e_bs = [0., 0.]
                            e_bs[(d % 2) + 1] = energy_bs
                            append!(events_bs, [e_bs])

                            # add original event tof to use as offset for global timing
                            t_bs = [Inf, Inf]
                            t_bs[(d % 2) + 1] = tofs_ini[n][d] + t_tof(energy_bs, create_ang())
                            append!(tof_bs, [t_bs])
                        end
                    end
                end
            end

            # Combine regular and backscatter events
            if with_bs
                append!(events, events_bs)
                append!(tofs_ini, tof_bs)
                n_ev = lastindex(events)
            end

            if det_main.is_measuring == true

                # check if one detector has not triggered (one must have triggered)
                d_untrig = -1.
                t_trig = []
                if minimum(curr_ev.t_trig) < 0.
                    d_untrig = argmin(curr_ev.t_trig)
                end

                # process events
                curr_ev.n_sum += n_ev
                for i = 1:n_ev
                    curr_ev.e_ind += events[i]

                    if d_untrig >= 0
                        # check if triggered, add time of trigger
                        if check_trigger(events[i][d_untrig])
                            append!(t_trig, tofs_ini[i][d_untrig])
                        end
                    end
                end

                if d_untrig >= 0 && lastindex(t_trig) > 0
                    curr_ev.t_trig[d_untrig] = t_now + minimum(t_trig)
                end


            else
                # sort events by time on detector (tof)
                t_perm = sortperm([minimum(tofs_ini[i]) for i=1:n_ev])

                # ToDo: add check if bs electron too late for integration! (Should be negligible)
                for perm in t_perm
                    if det_main.is_measuring == false
                        t1 = check_trigger(events[perm][1])
                        t2 = check_trigger(events[perm][2])

                        if t1 || t2

                            det_main.is_measuring = true
                            timing = [-1., -1.]
                            # if both trigger in the same time step, this order gives a bias towards det = 1
                            if t2
                                timing[2] = t_now + tofs_ini[perm][2]
                                det_trig = 2
                            end
                            if t1
                                timing[1] = t_now + tofs_ini[perm][1]
                                det_trig = 1
                            end

                            
                            curr_ev = Event(timing, events[perm], det_trig, 0)
                            # get index of untriggered detector (-1 if both have triggered)
                            d_untrig = (t1 * t2) ? -1. : (t2 * 1 + t1 * 2)
                        end


                    else
                        curr_ev.n_sum += 1
                        curr_ev.e_ind += events[perm]
                        if d_untrig >= 0
                            # check if triggered, add time of trigger
                            if check_trigger(events[perm][d_untrig])
                                curr_ev.t_trig[d_untrig] = t_now + tofs_ini[perm][d_untrig]
                            end
                        end

                    end

                end
                
            end
        end

        # check if detector still in measurement window
        if det_main.is_measuring == true
            if curr_ev.t_trig[1] < 0 || curr_ev.t_trig[2] < 0
                int_last_ev = max(curr_ev.t_trig...)
            else
                int_last_ev = min(curr_ev.t_trig...)
            end
            
            if  (t_now - int_last_ev) >= det_main.t_meas
                det_main.is_ready = false
                det_main.is_measuring = false
                t_last_ev = int_last_ev
                
                data[data_ind] = curr_ev
                data_ind += 1

                curr_ev = Event()
            end
        end

    end

    # remove rest of event buffer
    data = data[1:data_ind - 1]
    stats = Dict(
        "ev_count" => length(data),
        "ev_dt_corr" => dead_time_corr(length(data), t_meas=meas_time),
        "rate_obs_trig" => length(data) / meas_time,
        "rate_theo_pure" => theo_rate,
        "rate_dt_corr" => dead_time_corr(length(data), t_meas=meas_time) / meas_time,
    )
    if verbose
        println()
        for stat in stats
            @printf "%s\t%.0f \n" stat[1] stat[2]
        end    
    end

    return data, stats
end
