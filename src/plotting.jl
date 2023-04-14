

using Plots
include("birks.jl")


function plot_e_hist(
    data;
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 3.36,
    filename::String = "detsum",
)
    es = []
    for e in data
        push!(es, e.e_ind[1] + e.e_ind[2])
    end

    histogram(es, bins=0:e_bin:e_max, title="DetSum")
    ylabel!("Counts [ ]")
    xlabel!("Energy [keV]")
    # savefig(filename)
end


function plot_e_ind_hist(
    data;
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 3.36,
    filename::String = "det_ind",
)
    es0 = []
    es1 = []
    for e in data
        if e.e_ind[1] > 0.0
            push!(es0, e.e_ind[1])
        end
        if e.e_ind[2] > 0.0
            push!(es1, e.e_ind[2])
        end
    end

    plot(
        histogram(es0, bins=0:e_bin:e_max, title="Det0"),
        histogram(es1, bins=0:e_bin:e_max, title="Det1"),
    )
    ylabel!("Counts [ ]")
    xlabel!("Energy [keV]")
    # savefig(filename)
end


function plot_e_hist_tofcut(
    data;
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 3.36,
    filename::String = "detsum_tofcut",
)
    es = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 && abs(e.t_trig[1] - e.t_trig[2]) < 2e-8
            push!(es, e.e_ind[1] + e.e_ind[2])
        end
    end

    histogram(es, bins=0:e_bin:e_max, title="DetSum")
    ylabel!("Counts [ ]")
    xlabel!("Energy [keV]")
    # savefig(filename)
end


function plot_e_acc_hist(
    data;
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 3.36,
    filename::String = "detsum_acc",
)
    es = []
    for e in data
        if e.n_sum > 0
            push!(es, e.e_ind[1] + e.e_ind[2])
        end
    end
    println("No accidentals: ", length(es))
    histogram(es, bins=0:e_bin:_max, title="DetSum - Only Pile-Ups")
    ylabel!("Counts [ ]")
    xlabel!("Energy [keV]")
    # savefig(filename)
end


function plot_e_both_trig(
    data;
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 3.36,
    filename::String = "detsum_bothtrig",
)
    es = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(es, e.e_ind[1] + e.e_ind[2])
        end
    end
    
    histogram(es, bins=0:e_bin:e_max, title="DetSum - Both Triggered")
    ylabel!("Counts [ ]")
    xlabel!("Energy [keV]")
    # savefig(filename)
end


function plot_diff_both_trig(data; filename::String = "dtt")
    dt = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(dt, e.t_trig[1] - e.t_trig[2])
        end
    end

    histogram(dt, bins=(-2.5e-7):(1.0e-8):(2.5e-7), title="DeltaTriggerTime")
    ylabel!("Counts [ ]")
    xlabel!("TrigTime0 - TrigTime1 [10 ns]")
    # savefig(filename)
end


function plot2D_e_both_trig(
    data;
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 3.36,
    filename::String = "dtt_vs_detsum",
)
    es = []
    dt = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(es, e.e_ind[1] + e.e_ind[2])
            push!(dt, (e.t_trig[1] - e.t_trig[2]) * 1e8)
        end
    end
    
    histogram2d(
        dt,
        es,
        bins=((-2.5e1 ):(2):(2.5e1), 0:e_bin:e_max),
        title="DeltaTriggerTime vs DetSum",
        c=:imola,
        colorbar_title = "\nCounts [ ]",
        right_margin = 5Plots.mm,
    )
    xlabel!("TrigTime0 - TrigTime1 [10 ns]")
    ylabel!("Energy [keV]")
    # savefig(filename)
end


function plot2D_e_ind_both_trig(
    data;
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 3.36,
    e_min::Float64 = 0.,
    filename::String = "det0_vs_det1_bothtrig",
)
    es0 = []
    es1 = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0
            push!(es0, e.e_ind[1])
            push!(es1, e.e_ind[2])
        end
    end
    
    histogram2d(
        es0,
        es1,
        bins=(e_min:e_bin:e_max, e_min:e_bin:e_max),
        title="Energy: Det0 vs Det1 (Both Triggered)",
        c=:imola,
        colorbar_title = "\nCounts [ ]",
        right_margin = 5Plots.mm,
    )
    xlabel!("Energy0 [keV]")
    ylabel!("Energy1 [keV]")
    xlims!(0., e_max)
    ylims!(0., e_max)
    # savefig(filename)
end


function plot2D_e_ind(
    data;
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 3.36,
    e_min::Float64 = 0.,
    filename::String = "det0_vs_det1",
)
    es0 = []
    es1 = []
    for e in data
        if e.e_ind[1] > 0 || e.e_ind[2] > 0
            push!(es0, e.e_ind[1])
            push!(es1, e.e_ind[2])
        end
    end
    
    histogram2d(
        es0,
        es1,
        bins=(e_min:e_bin:e_max, e_min:e_bin:e_max),
        title="Energy: Det0 vs Det1",
        # c=:imola,
        c=cgrad(:imola),
        colorbar_title = "\nCounts [ ]",
        right_margin = 5Plots.mm,
    )
    xlabel!("Energy0 [keV]")
    ylabel!("Energy1 [keV]")
    xlims!(0., e_max)
    ylims!(0., e_max)
    # savefig(filename)
end

function plot2D_e_raw_e_det_rel(
    data;
    e_max::Float64 = 35.,
    e_bin::Float64 = 0.104,
    k_b::Float64 = 150., 
    k_off::Float64 = 0.,
    k_mul::Float64 = 1.,
    gain::Float64 = 30.9,
    filename::String = "eraw_vs_detsum",
)
    es = []
    rel_e = []
    for e in data
            push!(rel_e, (e.e_ind[1] + e.e_ind[2]) / e.e_raw)
            push!(es, 0.001 * e.e_raw)
    end

    histogram2d(
        es,
        rel_e,
        bins=(0:e_bin:e_max, 0.70:0.001:1.1),
        title="Raw vs Corrected Energy (rel)",
        c=:imola,
        colorbar_title = "\nCounts [ ]",
        right_margin = 5Plots.mm,
    )
    # x_plot, y_plot = calc_birks(k_b, k_off)
    # plot!(x_plot * gain * 0.001, ((y_plot./x_plot).* k_mul).+ k_off)
    ylims!(0.7, 1.1)
    xlabel!("Raw Energy [kch]")
    ylabel!("rel. Deviation [ ]")
    # savefig(filename)
end

function plot2D_e_raw_e_det_abs(
    data;
    e_max::Float64 = 35.,
    e_bin::Float64 = 0.104,
    k_b::Float64 = 150., 
    k_off::Float64 = 0.,
    k_mul::Float64 = 1.,
    gain::Float64 = 30.9,
    filename::String = "qdc_abs_nonlin",
)
    es = []
    rel_e = []
    for e in data
            push!(rel_e, 0.001 * ((e.e_ind[1] + e.e_ind[2]) - e.e_raw))
            push!(es, 0.001 * e.e_raw)
    end

    histogram2d(
        es,
        rel_e,
        bins=(0:e_bin:e_max, -1.2:0.01:1.2),
        title="Raw vs Corrected Energy (abs)",
        c=:imola,
        colorbar_title = "\nCounts [ ]",
        right_margin = 5Plots.mm,
    )
    # x_plot, y_plot = calc_birks(k_b, k_off)
    # plot!(x_plot * gain * 0.001, (((x_plot.-y_plot) * gain * 0.001).* k_mul).+ k_off)
    xlabel!("Raw Energy [kch]")
    ylabel!("abs. Deviation [kch]")
    # savefig(filename)
end