
using Plots


function plot_e_hist(
    data;
    gain::Float64 = 1.0,
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 20.,
    filename::String = "detsum",
)
    es = []
    for e in data
        push!(es, e.e_ind[1] + e.e_ind[2])
    end

    histogram(es, bins=0:(gain * e_bin):(gain * e_max), title="DetSum")
    ylabel!("Counts [ ]")
    xlabel!("Energy [$(gain) keV]")
    # savefig(filename)
end

function plot_e_ind_hist(
    data;
    gain::Float64 = 1.0,
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 20.,
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
        histogram(es0, bins=0:(gain * e_bin):(gain * e_max), title="Det0"),
        histogram(es1, bins=0:(gain * e_bin):(gain * e_max), title="Det1"),
    )
    ylabel!("Counts [ ]")
    xlabel!("Energy [$(gain) keV]")
    # savefig(filename)
end

function plot_e_acc_hist(
    data;
    gain::Float64 = 1.0,
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 20.,
    filename::String = "detsum_acc",
)
    es = []
    for e in data
        if e.n_sum > 0
            push!(es, e.e_ind[1] + e.e_ind[2])
        end
    end
    println("No accidentals: ", length(es))
    histogram(es, bins=0:(gain * e_bin):(gain * e_max), title="DetSum - Only Pile-Ups")
    ylabel!("Counts [ ]")
    xlabel!("Energy [$(gain) keV]")
    # savefig(filename)
end

function plot_e_both_trig(
    data;
    gain::Float64 = 1.0,
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 20.,
    filename::String = "detsum_bothtrig",
)
    es = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(es, e.e_ind[1] + e.e_ind[2])
        end
    end
    
    histogram(es, bins=0:(gain * e_bin):(gain * e_max), title="DetSum - Both Triggered")
    ylabel!("Counts [ ]")
    xlabel!("Energy [$(gain) keV]")
    # savefig(filename)
end


function plot_diff_both_trig(data; filename::String = "dtt")
    dt = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(dt, e.t_trig[1] - e.t_trig[2])
        end
    end
    
    histogram(dt, bins=(-2.5e-7):(2e-8):(2.5e-7), title="DeltaTriggerTime")
    ylabel!("Counts [ ]")
    xlabel!("TrigTime0 - TrigTime1 [10 ns]")
    # savefig(filename)
end


function plot2D_e_both_trig(
    data;
    gain::Float64 = 1.0,
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 20.,
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
        bins=((-2.5e1 ):(2):(2.5e1), 0:(gain * e_bin):(gain * e_max)),
        title="DeltaTriggerTime vs DetSum",
        c=:imola,
        colorbar_title = "\nCounts [ ]",
        right_margin = 5Plots.mm,
    )
    xlabel!("TrigTime0 - TrigTime1 [10 ns]")
    ylabel!("Energy [$(gain) keV]")
    # savefig(filename)
end

function plot2D_e_ind_both_trig(
    data;
    gain::Float64 = 1.0,
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 20.,
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
        bins=(0:(gain * e_bin):(gain * e_max), 0:(gain * e_bin):(gain * e_max)),
        title="Energey: Det0 vs Det1 (Both Triggered)",
        c=:imola,
        colorbar_title = "\nCounts [ ]",
        right_margin = 5Plots.mm,
    )
    xlabel!("Energy0 [$(gain) keV]")
    ylabel!("Energy1 [$(gain) keV]")
    # savefig(filename)
end

function plot2D_e_ind(
    data;
    gain::Float64 = 1.0,
    e_max::Float64 = 1.8e3,
    e_bin::Float64 = 20.,
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
        bins=(0:(gain * e_bin):(gain * e_max), 0:(gain * e_bin):(gain * e_max)),
        title="Energy: Det0 vs Det1",
        # c=:imola,
        c=cgrad(:imola),
        colorbar_title = "\nCounts [ ]",
        right_margin = 5Plots.mm,
    )
    xlabel!("Energy0 [$(gain) keV]")
    ylabel!("Energy1 [$(gain) keV]")
    # savefig(filename)
end
