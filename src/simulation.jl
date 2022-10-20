

using Printf
include("structs.jl")
include("tools.jl")


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


function run_sim(
    src_dict::Dict;
    meas_time::Float64 = 1e0,
    sim_step::Float64 = 5e-9,
    n_det::Int64 = 1,
    with_bs::Bool = false,
    t_meas::Float64 = 2e-7,
)
    det_main = Detector(t_meas=t_meas)

    # calculate theoretical, pure rate from source
    theo_rate = 0.
    for src in src_dict
        theo_rate += src[2].rate
    end

    n_time_steps = meas_time/sim_step
    # Should be large enough buffer
    data = Vector{Event}(undef, convert(Int, theo_rate * meas_time * n_det * 2.0)) # n_time_steps * 0.01
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

            e_split = [0.0, 0.0]
            e_split[i] = e_now
            t_tof_now = 0.0

            # do backscattering process
            if e_now > 0.1 && with_bs
                p = prob_bs(e_now)
                if p >= rand()
                    frac_bs = e_bs_frac(e_now)
                    e_split[i] = e_now * frac_bs
                    e_split[(i % 2) + 1] = e_now * (1.0 - frac_bs)
                    t_tof_now = t_tof(e_now)   
                end
            end

            # detect energy - i is main
            for j = 1:n_det
                if det_main.is_measuring == true
                    if e_split[j] >= 0.1
                        # Set trigger time, if not set. Add tof for secondary detector
                        if curr_ev.t_trig[j] < 0
                            curr_ev.t_trig[j] = t_now + abs(i - j) * t_tof_now
                        end
                        curr_ev.e_det += e_now
                        curr_ev.n_sum += 1
                        curr_ev.e_ind[j] += e_split[j]
                    end
                else
                    if e_split[j] >= det_main.thresh
                        # Todo: This adds trigger bias towards first detector. Should not be a problem as dt = 1e-8
                        timing = [-1., -1.]
                        timing[i] = t_now
                        if t_tof_now > 0.0
                            timing[(i % 2) + 1] = t_now + t_tof_now
                        end

                        # e_ind = 
                        # e_ind[j] = e_split[j]

                        curr_ev = Event(timing, e_now, e_split, 0, 0)
                        det_main.is_measuring = true
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

    data = data[1:data_ind - 1]

    @printf "Obs. trigger rate %.0f \n" length(data) / (n_det * meas_time)
    println("Event counts ", length(data))
    @printf "DT corr. rate  %.0f \n" dead_time_corr(length(data)) / (n_det * meas_time)
    @printf "Theor. pure rate of srcs %.0f \n" theo_rate

    return data
end
