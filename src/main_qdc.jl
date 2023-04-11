
include("structs.jl")
include("plotting.jl")
include("birks.jl")
include("simulation.jl")

using Measures
using Statistics
using LsqFit


function create_energy(e_max::Float64 = 500., fix_e::Bool = false)
    if fix_e
        return e_max
    else
        return rand() * e_max
    end
end


pmt_nph = [
    1.036100e+04, 1.975300e+04, 1.921400e+04, 1.173900e+04,
    1.205700e+04, 1.863100e+04, 1.811600e+04, 1.354500e+04
]
norm = sum(pmt_nph)
pmt_dist = pmt_nph./norm


function split_energy(energy::Float64, n_pmts::Int64 = 16, split_equ::Bool = false)
    charge_delay = false

    # split equally among all 
    if split_equ
        split = zeros(n_pmts).+ (energy / n_pmts)
    else
        # ang = 35.7
        ang = create_ang()

        # calc probability of p_bs
        p_bs = prob_bs(energy, ang)
        p = rand()
        energy_bs = 0.

        # Fixme: This does only do bs without mirroring!
        if p_bs >= p
            charge_delay = true
            frac_bs = e_bs_frac(energy, ang)
            norm = frac_bs * 2.6 + 1.6
            frac_bs = clamp(frac_bs + randn() * frac_bs / norm, 0., 1.)
            
            energy_bs = (1.0 - frac_bs) * energy
            energy = energy - energy_bs
        end
        # Assign detector
        if rand() > 0.5
            split_energy = [energy, energy_bs] 
        else
            split_energy = [energy_bs, energy]
        end

        split = vcat(pmt_dist.* split_energy[1], pmt_dist.* split_energy[2])
    end

    return split, charge_delay
end


function corr_energy(e::Float64, m::Float64, c::Float64, del_t::Int64=18)
    # e_corr = e - (c + m * e) * del_t
    return e * (1. - del_t * m) - del_t * c
end


function run_qdcsim(
    n_evs::Int64 = 1;
    n_pmts::Int64 = 16,
    with_grad::Bool = false,
    e_max::Float64=1100.,
    fix_e::Bool=false,
)

    det = DetectorQDC()
    # data buffer
    data = Vector{Event}(undef, n_evs)

    # main loop
    for ev_c in (1:n_evs)
        curr_ev = Event()
        en_kev = create_energy(e_max, fix_e)
        en_ch = en_kev * det.gain
        curr_ev.e_raw = en_ch

        # ToDo: Expand with uneven split and differntiate between detectors (2x8)
        en_split, delayed = split_energy(en_kev, n_pmts).* det.gain
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

    return data, det
end


res, det = run_qdcsim(40000; with_grad=true)
plot2D_e_raw_e_det_abs(res; k_b=320., k_off=-1.9, gain=det.gain)
plot2D_e_raw_e_det_rel(res; k_b=470., k_off=0., k_mul=1.159, gain=det.gain)

# rel add
# 20: 650, 0.177
# 18: 600, 0.166
# 16: 500, 0.144 / 450, 0.133
# 13: 350, 0.108

# abs add
# 20: 370, -2.15
# 18: 320, -1.90
# 16: 280, -1.70
# 13: 235, -1.50

function gen_qdc_data(x_in; n_smp::Int64 = 500)
    mus = []
    sigs = []
    for en in x_in
        curr_res, det = run_qdcsim(n_smp; with_grad=true, e_max=en, fix_e=true,)
        ev_ens = []
        for ev in curr_res
                push!(ev_ens, (ev.e_ind[1] + ev.e_ind[2]) / ev.e_raw)
        end

        mu = mean(ev_ens)
        sig = std(ev_ens)
        push!(mus, mu)
        push!(sigs, sig)
    end

    return mus, sigs
end

function gen_birks_comp(
    k_b::Float64 = 150., 
    k_off::Float64 = 0.,
    k_mul::Float64 = 1.,
    fit_birks::Bool=true,
)
    # xs = range(50, 1100, length=50)
    xs = [
        # n_bins = [10., 16., 36., 41., 52., 68.]; norm = 11.
        # Cd
        # 60., 
        # Ce
        # 105.,
        # Sn
        330., 345., 360.,
        # Bi mid
        420., 460., 480., 520.,
        # Cs
        570., 600., 640., 680., 740., 
        # Bi high
        870., 910., 940, 980., 1020., 1050., 1100., 
    ]
    # xs = vcat(range(250, 1100, length=50))
    # xs = vcat(range(50, 1100, length=50))
    mus, sigs = gen_qdc_data(xs)

    # sigs = sigs .* 3.
    # sigs[1] *= 6
    # sigs[2] *= 6

    if fit_birks
        lb = [50., 0.5]
        ub = [1500., 1.5] 
        wt = 1 ./ sigs .^2

        fit = curve_fit(f_birks, xs, mus, wt, [400., 1.159], lower=lb, upper=ub)
        println(coef(fit))
        println(stderror(fit))

        p1 = scatter(xs, mus, yerror=sigs, label="Sim Fit Data", ms=2, xaxis=false, xticks = 0:200:1100) #, markerstrokewidth=2)
        x_plot = vcat(range(10, 1100, length=100))
        y_plot = f_birks(x_plot, coef(fit))
        kb_plot = round(coef(fit)[1], digits=1)
        kb_plot_err =  round(stderror(fit)[1], digits=1)
        plot!(x_plot, y_plot, label="Fit (kB = $(kb_plot) +- $(kb_plot_err))")
        mus, sigs = gen_qdc_data(x_plot)
        plot!(x_plot, mus, grid=true, label="Sim Full")
        ylims!(0.7, 1.1)
        ylabel!("rel. Deviation [ ]")

        p2 = plot(x_plot, y_plot .- mus, label="Birks - Sim", xticks=0:200:1100, yticks=-0.01:0.005:0.01)
        ylims!(-0.0095, 0.0095)
        xlabel!("Raw Energy [keV]")

        plot(p1, p2 ,layout=Plots.grid(2,1, heights=[0.75, 0.25]), bottom_margin=[-5mm 0mm])
        savefig("qdc_nonlin_comp")

    else
        plot(xs, mus, grid=false, yerror=sigs)
        x_plot = range(10, 1100, length=100)
        y_plot = [f_birks(x_i, [k_b, k_off, k_mul]) for x_i in x_plot]
        plot!(x_plot, y_plot)
        ylims!(0.7, 1.1)
        xlabel!("Raw Energy [kch]")
        ylabel!("rel. Deviation [ ]")
    end
end

gen_birks_comp(470., 0., 1.159)


x_plot = vcat(range(10, 1100, length=100))
y_plot = f_birks(x_plot, [490., 1.])
y_plot2 = f_birks(x_plot, [490., 1.]).* f_birks(x_plot, [150., 1.]) * 1.05
plot(x_plot, y_plot)
plot!(x_plot, y_plot2)