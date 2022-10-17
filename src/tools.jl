

function dead_time_corr(counts::Int64, t_dead::Float64 = 1.5e-6, t_meas::Float64 = 1.0)
    return counts / (1.0 - (t_dead * counts) / t_meas)
end