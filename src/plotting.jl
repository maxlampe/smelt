
using Plots


function plot_e_hist(data, gain::Float64 = 1.0)
    es = []
    for e in data
        push!(es, e.e_det)
    end

    histogram(es, bins=0:(gain * 33):(gain * 2.2e3), title="DetSum")
    ylabel!("Counts [ ]")
    xlabel!("Energy [$(gain) keV]")
end

function plot_e_ind_hist(data, gain::Float64 = 1.0)
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
        histogram(es0, bins=0:(gain * 33):(gain * 2.2e3), title="Det0"),
        histogram(es1, bins=0:(gain * 33):(gain * 2.2e3), title="Det1"),
    )
    ylabel!("Counts [ ]")
    xlabel!("Energy [$(gain) keV]")
end

function plot_e_acc_hist(data, gain::Float64 = 1.0)
    es = []
    for e in data
        if e.n_sum > 0
            push!(es, e.e_det)
        end
    end
    println("No accidentals: ", length(es))
    histogram(es, bins=0:(gain * 33):(gain * 2.2e3), title="DetSum - Only Pile-Ups")
    ylabel!("Counts [ ]")
    xlabel!("Energy [$(gain) keV]")
end

function plot_e_both_trig(data, gain::Float64 = 1.0)
    es = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(es, e.e_det)
        end
    end
    
    histogram(es, bins=0:(gain * 33):(gain * 2.2e3), title="DetSum - Both Triggered")
    ylabel!("Counts [ ]")
    xlabel!("Energy [$(gain) keV]")
end


function plot_diff_both_trig(data)
    dt = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(dt, e.t_trig[1] - e.t_trig[2])
        end
    end
    
    histogram(dt, bins=(-2.5e-7):(2e-8):(2.5e-7), title="DeltaTriggerTime")
    ylabel!("Counts [ ]")
    xlabel!("TrigTime0 - TrigTime1 [10 ns]")
end


function plot2D_e_both_trig(data, gain::Float64 = 1.0)
    es = []
    dt = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(es, e.e_det)
            push!(dt, (e.t_trig[1] - e.t_trig[2]) * 1e8)
        end
    end
    
    histogram2d(dt, es, nbins=40, title="DeltaTriggerTime vs DetSum", c=:imola)
    xlabel!("TrigTime0 - TrigTime1 [10 ns]")
    ylabel!("Energy [$(gain) keV]")
end

function plot2D_e_ind_both_trig(data, gain::Float64 = 1.0)
    es0 = []
    es1 = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0
            push!(es0, e.e_ind[1])
            push!(es1, e.e_ind[2])
        end
    end
    
    histogram2d(es0, es1, bins=(0:(gain * 50):(gain * 1.7e3), 0:(gain * 50):(gain * 1.7e3)), title="Det0 vs Det1 (Both Triggered)", c=:imola)
    xlabel!("Energy0 [$(gain) keV]")
    ylabel!("Energy1 [$(gain) keV]")
end

function plot2D_e_ind(data, gain::Float64 = 1.0)
    es0 = []
    es1 = []
    for e in data
        if e.e_ind[1] > 0 || e.e_ind[2] > 0
            push!(es0, e.e_ind[1])
            push!(es1, e.e_ind[2])
        end
    end
    
    histogram2d(es0, es1, bins=(0:(gain * 50):(gain * 1.7e3), 0:(gain * 50):(gain * 1.7e3)), title="Det0 vs Det1 (Both Triggered)", c=:imola)
    xlabel!("Energy0 [$(gain) keV]")
    ylabel!("Energy1 [$(gain) keV]")
end
