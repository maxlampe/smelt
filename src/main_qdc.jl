
include("structs.jl")
include("plotting.jl")


function create_energy(e_max::Float64 = 500.)
    return rand() * e_max
    # return e_max
end

function split_energy(energy::Float64, n_pmts::Int64 = 16)
    # split equally among all 
    split = zeros(n_pmts).+ (energy / n_pmts)
    return split
end

function corr_energy(e::Float64, m::Float64, c::Float64, del_t::Int64=16)
    return e * (1. + del_t * m) - del_t * c
end

function run_qdcsim(
    n_evs::Int64 = 1;
    n_pmts::Int64 = 16,
    with_grad::Bool = false,
    e_max::Float64=1100.,
)

    det = DetectorQDC()
    # data buffer
    data = Vector{Event}(undef, n_evs)

    # main loop
    for ev_c in (1:n_evs)
        curr_ev = Event()
        en_kev = create_energy(e_max)
        en_ch = en_kev * det.gain
        curr_ev.e_raw = en_ch

        # ToDo: Expand with uneven split and differntiate between detectors (2x8)
        en_split = split_energy(en_ch, n_pmts)
        curr_ev.e_ind = [0., 0.]

        # iterate over PMTs and "add" each PMT to event≈
        for n in (1:n_pmts)
            curr_en = en_split[n]
            if with_grad
                curr_en = corr_energy(curr_en, det.grad_m[n], det.grad_c[n])
            end

            if n <= 8
                curr_ev.e_ind[1] += curr_en
            else
                curr_ev.e_ind[2] += curr_en
            end

        end

        # set event trigger and time with dummies?
        curr_ev.n_sum += 1
        data[ev_c] = curr_ev
    end

    return data
end


res = run_qdcsim(20000; with_grad=true)
plot2D_e_raw_e_det_rel(res; k_b=450., k_off=0.133)
# plot2D_e_raw_e_det_abs(res)

# 18: 600, 0.166
# 16: 500, 0.144 / 450, 0.133
# 13: 350, 0.108
