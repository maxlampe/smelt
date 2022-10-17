using Plots

function plot_e_hist(data)
    es = []
    for e in data
        push!(es, e.e_det)
    end

    histogram(es, bins=0:500:60e3)
end