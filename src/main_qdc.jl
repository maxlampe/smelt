
include("structs.jl")
include("plotting.jl")
include("birks.jl")
include("simulation.jl")

using Measures
using Statistics
using LsqFit


function create_energy(e_max::Float64 = 500., fix_e::Bool = false, with_resolution::Bool = false)
    if fix_e
        if with_resolution
            # resolution = sqrt(alpha * e * PE) with empirical gain factor alpha, as std of poisson is sqrt(mu) 
            resol = sqrt(e_max / 0.6 * 2.4) 

            # resol = -1.178e-4 * e_max ^ 2. + 0.24 * e_max + 24.13
            # # FWHM to Sig
            # resol = resol / 2.355

            return clamp(e_max + randn() * resol, 1., 3. * e_max) 
        else
            return e_max
        end
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
        p_bs_raw, p_mirror = prob_bs(energy, ang, false)
        p_bs = p_bs_raw * p_mirror
        p = rand()
        energy_bs = 0.

        if p_bs_raw >= p
            charge_delay = true
            p2 = rand()
            if p_mirror >= p2
                frac_bs = e_bs_frac(energy, ang)
                norm = frac_bs * 2.6 + 1.6
                frac_bs = clamp(frac_bs + randn() * frac_bs / norm, 0., 1.)
                
                energy_bs = (1.0 - frac_bs) * energy
                energy = energy - energy_bs
            end
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


function corr_energy(e::Float64, m::Float64, c::Float64, del_t::Int64=17) # 17
    # e_corr = e - (c + m * e) * del_t
    # return e * (1. - del_t * m) - del_t * c
    return e - (c + m * e) * del_t
end


function corr_energy_bs(e::Float64, e_corr::Float64, del_m::Float64, del_c::Float64, del_t::Int64=10) # 10
    # e_corr_bs = e_corr - (del_c + del_m * e) * del_t
    return e_corr - (del_c + del_m * e) * del_t
end 


function run_qdcsim(
    n_evs::Int64 = 1;
    n_pmts::Int64 = 16,
    with_grad::Bool = false,
    with_delay::Bool = false,
    with_prebirks::Bool=false,
    with_resolution=false,
    e_max::Float64=1100.,
    fix_e::Bool=false,
)

    # Pure scintilator non-linearity from Diss. Saul 2018: 123. +- 14 [nm/keV]
    # Include dead layer? Says something like 158?
    k_b_theo = 123.
    det = DetectorQDC()
    if with_prebirks
        birks_norm = f_birks([1000.], [k_b_theo, 1.])[1]
    end

    # data buffer
    data = Vector{Event}(undef, n_evs)

    # main loop
    for ev_c in (1:n_evs)
        curr_ev = Event()
        en_kev = create_energy(e_max, fix_e, with_resolution)
        if with_prebirks
            en_kev = en_kev * f_birks([en_kev], [k_b_theo, 1.])[1] / birks_norm
        end
        en_ch = en_kev * det.gain
        curr_ev.e_raw = en_ch

        # ToDo: Expand with uneven split and differntiate between detectors (2x8)
        en_split, delayed = split_energy(en_kev, n_pmts)
        en_split = en_split.* det.gain
        curr_ev.e_ind = [0., 0.]

        # iterate over PMTs and "add" each PMT to event≈
        for n in (1:n_pmts)
            curr_en = en_split[n]
            if with_grad
                curr_en_corr = corr_energy(curr_en, det.grad_m[n], det.grad_c[n])
                if with_delay && delayed
                    delta_m = det.grad_m_bs[n] - det.grad_m[n]
                    delta_c = det.grad_c_bs[n] - det.grad_c[n]
                    curr_en_corr = corr_energy_bs(curr_en, curr_en_corr, delta_m, delta_c)
                end
                curr_en = curr_en_corr
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


res, det = run_qdcsim(40000; with_grad=true, with_delay=true, with_prebirks=false, with_resolution=false)
plot2D_e_raw_e_det_abs(res; k_b=320., k_off=-1.9, gain=det.gain)
# plot2D_e_raw_e_det_rel(res; k_b=470., k_off=0., k_mul=1.159, gain=det.gain)


function gen_qdc_data(x_in; n_smp::Int64 = 2000)
    mus = []
    sigs = []
    for en in x_in
        curr_res, det = run_qdcsim(n_smp; with_grad=true, e_max=en, fix_e=true, with_delay=true, with_prebirks=true, with_resolution=true)
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
    # Use with energy resolution approximation
    xs = [
        # Cd, Ce, Sn, Bi mid, Cs, Bi high
        # 60., 105., 
        344., 460., 640., 990.,
    ]
    # xs = [
    #     # n_bins = [10., 16., 36., 41., 52., 68.]; norm = 11.
    #     # Cd
    #     # 60., 
    #     # Ce
    #     # 105.,
    #     # Sn
    #     330., 345., 360.,
    #     # Bi mid
    #     420., 460., 480., 520.,
    #     # Cs
    #     570., 600., 640., 680., 740., 
    #     # Bi high
    #     870., 910., 940, 980., 1020., 1050., 1100., 
    # ]
    # xs = vcat(range(250, 1100, length=50))
    # xs = vcat(range(50, 1100, length=50))
    mus, sigs = gen_qdc_data(xs)

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
        mus_full, sigs_full = gen_qdc_data(x_plot)
        plot!(x_plot, mus_full, grid=true, label="Sim Full")
        # ylims!(0.7, 1.1)
        ylims!(0.4, 1.1)
        ylabel!("rel. Deviation [ ]")

        p2 = plot(x_plot, y_plot .- mus_full, label="Birks - Sim", xticks=0:200:1100, yticks=-0.01:0.005:0.01)
        ylims!(-0.0095, 0.0095)
        xlabel!("Raw Energy [keV]")

        # return plot(p1, p2 ,layout=Plots.grid(2,1, heights=[0.75, 0.25]), bottom_margin=[-5mm 0mm])
        # savefig("qdc_nonlin_comp")

        # Fit simulation nonlinearity empiricaly

        wt_full = 1 ./ sigs_full .^2
        # fit_exp = curve_fit(f_emp, x_plot, mus_full, wt_full, [1.1, 1., -0.005])
        # fit_exp = curve_fit(f_emp, x_plot, mus_full, wt_full, [1.1, 1., -0.005, 1., -0.01])
        fit_exp = curve_fit(f_emp, x_plot, mus_full, [1.1, 1., -0.005, 1., -0.01, 0.5, -0.05])
        # fit_exp = curve_fit(f_emp, x_plot, mus_full, [1.1, 1., -0.005, 0.5, -0.001])
        # fit_exp = curve_fit(f_emp, x_plot, mus_full, [1.04, 40.])

        
        y_plot_exp = f_emp(x_plot, coef(fit_exp))
        println(coef(fit_exp))
        # println(stderror(fit_exp))
        p4 = plot(x_plot, mus_full, grid=true, label="Smelt Sim. Continuous", lw=2)
        plot!(x_plot, y_plot_exp, label="Fit (empirical)", lw=2)
        plot!(x_plot, f_birks(x_plot, [470., 1.15]), label="Birks", lw=2)
        
        plot!(x_plot, 1.05 .- 35. ./ x_plot)
        ylims!(0.0, 1.1)
        xlabel!("Energy [keV]")
        ylabel!("relative Deviation [a.u.]")

        rchi2 = sum(((y_plot_exp - mus_full) ./ sigs_full).^2)
        println(rchi2, lastindex(mus_full))
        println(rchi2 / (lastindex(mus_full) - 7))

        # savefig("qdc_smeltsim_qdcempirical.pdf")

        return p4

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


function plot_energy_resolution()
    xplot = vcat(range(50, 1100, length=50))
    plot(xplot, (xplot ./ 0.6 * 2.4) .^ 0.5 .* 2.355, label="Model")
    scatter!([340.], [90.], label="Sn")
    scatter!([1000.], [146.], label="Bi")
    scatter!([105.], [48.], label="Ce")
    scatter!([610.], [122.], label="Cs")
    xlabel!("Energy [keV]")
    ylabel!("FWHM [keV]")
    # savefig("fwhm_model")
end
plot_energy_resolution()

function f_expo(x_en::Vector{Float64}, par::Vector{Float64})
    a_0 = par[1]
    a_1 = par[2]
    # a_2 = par[3]

    y = [(x - (a_0 * x + a_1)) / (x) for x in x_en]

    return y
end 

function f_emp(x_en::Vector{Float64}, par::Vector{Float64})
    # a_0 = par[1]
    # a_1 = par[2]
    # a_2 = par[3]
    # y = a_0 .-  exp.(x_en .* a_2) .* a_1

    a_0 = par[1]
    a_1 = par[2]
    a_2 = par[3]
    a_3 = par[4]
    a_4 = par[5]
    a_5 = par[6]
    a_6 = par[7]
    y = a_0 .-  exp.(x_en .* a_2) .* a_1 .-  exp.(x_en .* a_4) .* a_3 .-  exp.(x_en .* a_6) .* a_5
    # [1.02130, 0.32546, -0.00438, 1.99055, -0.02619, 73.5903, -0.14360]
    # rChi2 0.17906381977672683

    # a_0 = par[1]
    # a_1 = par[2]
    # a_2 = par[3]
    # a_3 = par[4]
    # a_4 = par[5]

    # y = a_0 .-  exp.(x_en .* a_2) .* a_1 .-  exp.(x_en .* a_4) .* a_3 
    # rChi2 32.

    # a_0 = par[1]
    # a_1 = par[2]
    # y = a_0 .- a_1 ./ x_en
    # #  a_i = [1.05, -35.]

    return y
end 

x_plot = vcat(range(10, 1100, length=100))
plot(x_plot, f_emp(x_plot, [1.02130, 0.32546, -0.00438, 1.99055, -0.02619, 73.5903, -0.14360]), ylims=[0.3, 1.1])
plot!(x_plot, 1.05 .- 35. ./ x_plot)



y_plot = [(1.1 - exp(-0.005 * x)) for x in x_plot]

y_plot = [(x - (1.54 * x - 32.5)) / (x) for x in x_plot]
plot(x_plot, y_plot)
