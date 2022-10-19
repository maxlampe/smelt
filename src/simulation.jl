
include("structs.jl")


# ToDo: currently only works for dtime * SIM_STEP << 1
function create_energy(src_dict::Dict, dtime::Float64)::Float64
    energy = 0.
    for src in src_dict
        if src[2].rate * dtime >= rand()
            energy += randn() * src[2].width + src[2].energy
        end
    end
    return energy
end


function run_sim(src_dict::Dict; meas_time::Float64 = 1e0, sim_step::Float64 = 1e-8, n_det::Int64 = 1)
    det_main = Detector()

    n_time_steps = meas_time/sim_step
    # Should be large enough buffer
    data = Vector{Event}(undef, convert(Int, n_time_steps * 0.01))
    data_ind = 1
    curr_ev = Event()
    t_last_ev = 0.

    for t = 1:(n_time_steps)
        if t % (0.1 * n_time_steps) == 0
            println("At ", 100. * t/n_time_steps, "%")
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

        # Independent event for each detector.
        for i = 1:n_det
            # create energy from all sources at this time step
            e_now = create_energy(src_dict, sim_step)

            # detect energy
            if det_main.is_measuring == true
                if e_now >= 0.1
                    if curr_ev.t_trig[i] < 0
                        curr_ev.t_trig[i] = t_now
                    end
                    curr_ev.e_det += e_now
                    curr_ev.n_sum += 1
                    curr_ev.e_ind[i] += e_now
                end
            else
                if e_now >= det_main.thresh
                    # Todo: This adds trigger bias towards first detector. Should not be a problem as dt = 1e-8
                    timing = [-1., -1.]
                    timing[i] = t_now
                    e_ind = [0., 0.]
                    e_ind[i] = e_now
                    curr_ev = Event(timing, e_now, e_ind, 0, 0)
                    det_main.is_measuring = true
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

    data = data[1:data_ind - 1]

    println("Event rate ", length(data))
    println("DT corr: ", dead_time_corr(length(data)))
    tot_rate = 0.
    for src in src_dict
        tot_rate += src[2].rate
    end
    println("Theor. rate of sources", tot_rate * meas_time)

    return data
end
