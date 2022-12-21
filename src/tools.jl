using BasicInterpolators

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


function t_tof(en::Float64, thet::Float64 = 35.7)
    return 6.4 / (velo_rel(en) * cos(thet * pi / 180.))
end


function exp_dec(x::Float64, a::Float64, b::Float64, c::Float64)
    return c + a * exp(-b * x)
end


# function prob_bs(en::Float64, ang::Float64 = 35.7)
#     p_raw = exp_dec(en, 0.061, 0.00083, 0.033)
#     p_mirror = exp_dec(en, 0.032, 0.00084, 0.471) 
#     return p_raw * p_mirror
# end
function prob_bs(en::Float64, ang::Float64 = 35.7)
    p_raw = interpol_p_bs_raw(ang, en)
    p_mirror = interpol_p_bs_mirror(ang, en) 
    return p_raw * p_mirror
end

# function e_bs_frac(en::Float64)
#     return exp_dec(en, -0.232, 0.001, 0.660) 
# end
function e_bs_frac(en::Float64, ang::Float64 = 35.7)
    return interpol_e_bs_frac(ang, en) 
end


function trigger_func(en::Float64; a::Float64 = 0.093, p::Float64 = 0.78)

    val = 1.0 - (1.0 - p)^(a * en) * (1.0 + (a * p * en) / (1.0 - p))
    if val < 0.0
        return 0.0
    else
        return val
    end
end


data = [
    [0., 50.], # 0 deg
    [0., 200.],
    [0., 400.],
    [0., 600.],
    [0., 800.],
    [0., 1000.],
    [15., 50.], # 15 deg
    [15., 200.],
    [15., 400.],
    [15., 600.],
    [15., 800.],
    [15., 1000.],
    [35.7, 50.], # 35.7 deg
    [35.7, 200.],
    [35.7, 400.],
    [35.7, 600.],
    [35.7, 800.],
    [35.7, 1000.],
    [42., 50.], # 42 deg
    [42., 200.],
    [42., 400.],
    [42., 600.],
    [42., 800.],
    [42., 1000.],
    [47., 50.], # 47 deg
    [47., 200.],
    [47., 400.],
    [47., 600.],
    [47., 800.],
    [47., 1000.],
]
data = reduce(vcat,transpose.(data))
# raw probability of backscattering
p_bs_raw = [
    0.046, 0.036, 0.03, 0.026, 0.022, 0.019, # 0 deg
    0.052, 0.043, 0.037, 0.032, 0.028, 0.024, # 15 deg
    0.092, 0.084, 0.077, 0.070, 0.065, 0.059, # 35.7 deg
    0.116, 0.108, 0.101, 0.094, 0.088, 0.082, # 42 deg
    0.141, 0.133, 0.126, 0.119, 0.112, 0.105, # 47 deg
]
# probability to be reflected due to magnetic mirror effect
p_bs_mirror = [
    0.460, 0.448, 0.443, 0.44, 0.44, 0.44, # 0 deg
    0.466, 0.456, 0.449, 0.446, 0.446, 0.442, # 15 deg
    0.503, 0.498, 0.494, 0.492, 0.488, 0.486, # 35.7 deg
    0.527, 0.524, 0.522, 0.52, 0.517, 0.514, # 42 deg
    0.552, 0.552, 0.55, 0.548, 0.546, 0.544, # 47 deg
]
# fraction of energy deposited in main detector
frac_e_bs = [
    0.442, 0.492, 0.531, 0.563, 0.585, 0.603, # 0 deg
    0.442, 0.492, 0.531, 0.563, 0.585, 0.603, # 15 deg (= 0 deg)
    0.429, 0.467, 0.506, 0.537, 0.559, 0.581, # 35.7 deg
    0.419, 0.456, 0.490, 0.519, 0.542, 0.564, # 42 deg
    0.407, 0.441, 0.473, 0.502, 0.525, 0.546, # 47 deg
]
# frac, gaussian, 500
# mirr, gaussian, 750
# p_raw, gaussian, 600
interpol_p_bs_raw = RBFInterpolator(data, p_bs_raw, 500, gaussian)
interpol_p_bs_mirror = RBFInterpolator(data, p_bs_mirror, 750, gaussian)
interpol_e_bs_frac = RBFInterpolator(data, frac_e_bs, 500, gaussian)


