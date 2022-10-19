
using Plots


function plot_e_hist(data, gain::Float64 = 1.0)
    es = []
    for e in data
        push!(es, e.e_det)
    end

    histogram(es, bins=0:(gain * 50):(gain * 2.2e3))
end

function plot_e_acc_hist(data, gain::Float64 = 1.0)
    es = []
    for e in data
        if e.n_sum > 0
            push!(es, e.e_det)
        end
    end

    histogram(es, bins=0:(gain * 50):(gain * 2.2e3))
end

function plot_e_both_trig(data, gain::Float64 = 1.0)
    es = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(es, e.e_det)
        end
    end
    
    histogram(es, bins=0:(gain * 50):(gain * 2.2e3))
end

function plot2D_e_both_trig(data, gain::Float64 = 1.0)
    es = []
    dt = []
    for e in data
        if e.t_trig[1] > 0 && e.t_trig[2] > 0 
            push!(es, e.e_det)
            push!(dt, e.t_trig[1] - e.t_trig[2])
        end
    end
    
    histogram2d(dt, es, nbins=40)
end
