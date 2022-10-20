

function dead_time_corr(counts::Int64; t_dead::Float64 = 1.5e-6, t_meas::Float64 = 1.0)
    return counts / (1.0 - (t_dead * counts) / t_meas)
end


function velo_rel(en::Float64, mass::Float64 = 511.)
    frac = 1 / ( (en / mass) + 1)
    return 3e8 * (1 - frac ^ 2 ) ^ 0.5
end


function r_gyr(en::Float64, b_0::Float64 = 0.084, b_1::Float64 = 0.165)
    sq = sqrt(en * (en + 2. * 511.)) * 1.6e-16
    b_ratio = sqrt(b_1 / b_0)
    norm = 1. / (3e8 * 1.6e-19 * b_1) 
    return sq * norm * b_ratio
end


function t_tof(en::Float64, thet::Float64 = 35.7 * pi / 180.)
    return 6.4 / (velo_rel(en) * cos(thet))
end


function exp_dec(x::Float64, a::Float64, b::Float64, c::Float64)
    return c + a * exp(-b * x)
end


function prob_bs(en::Float64)
    p_raw = exp_dec(en, 0.061, 0.00083, 0.033)
    p_mirror = exp_dec(en, 0.032, 0.00084, 0.471) 
    return p_raw * p_mirror
end


function e_bs_frac(en::Float64)
    return exp_dec(en, -0.232, 0.001, 0.660) 
end
