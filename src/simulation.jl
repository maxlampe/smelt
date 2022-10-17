
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


function run_single_det(src_dict::Dict, meas_time::Float64 = 1e0, sim_step::Float64 = 1e-8)
    det_main = Detector(3e-7, 1.5e-6, 1., true, false)

    n_time_steps = meas_time/sim_step
    data = Vector{Event}(undef, convert(Int, n_time_steps))
    data_ind = 1
    curr_ev = Event(-1., -1., -1)
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

        # create energy from all sources at this time step
        e_now = create_energy(src_dict, sim_step)

        # detect energy
        if det_main.is_measuring == true
            curr_ev.e_det += e_now
        else
            if e_now >= det_main.thresh
                curr_ev = Event(t_now, e_now, 0)
                det_main.is_measuring = true
            end
        end

        # check if detector still in measurement window
        if det_main.is_measuring == true
            if  (t_now - curr_ev.t_trig) >= det_main.t_meas
                det_main.is_ready = false
                det_main.is_measuring = false
                t_last_ev = curr_ev.t_trig
                
                data[data_ind] = curr_ev
                data_ind += 1

                curr_ev = Event(-1., -1., -1)
            end
        end
        
    end

    return data[1:data_ind - 1]
end
