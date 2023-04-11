

using Parameters


@with_kw mutable struct Detector
    t_meas::Float64 = 3e-7
    t_dead::Float64 = 1.5e-6
    thresh::Float64 = 1.0
    is_ready::Bool= true
    is_measuring::Bool = false
end

# {
# 0: [[-0.0027726302119042904, 4.570725538613477], [-0.0029869500693498214, 5.665745447865509]], 
# 1: [[-0.0029366127395134534, 3.4537377697297393], [-0.003128027880529408, 4.3098087230480715]], 
# 2: [[-0.0033613507916551142, 4.907752978279222], [-0.0034747499158074486, 5.7618482958693855]], 
# 3: [[-0.0037734805031322586, 2.241773885191816], [-0.003967817605706616, 2.879599787711115]], 
# 4: [[-0.00307167001168234, 2.864614759296568], [-0.0033038002943305824, 3.8253692681421647]], 
# 5: [[-0.0030294520426235568, 4.03886453341743], [-0.0032875460101800913, 5.109554456853292]], 
# 6: [[-0.0028199306577301486, 2.893378071172577], [-0.002999059160725512, 3.7311722430827845]], 
# 7: [[-0.003168238958527181, 1.9486375219125216], [-0.0033585689857879046, 2.896290449097631]], 
# 8: [[-0.0027022291886093473, 2.8794068605780905], [-0.0032894121389988643, 4.511791339287047]], 
# 9: [[-0.0026675477947363285, 2.8794068605780905], [-0.003217170817858253, 4.418720738354981]], 
# 10: [[-0.002515804915587131, 2.912107718286143], [-0.0029987469811330685, 3.954548505472813]], 
# 11: [[-0.0020574023033302982, 6.996926904145113], [-0.0025855172841346046, 8.24046943476849]], 
# 12: [[-0.003092297887553726, 5.326732324708251], [-0.003422874313732227, 6.786966101391472]], 
# 13: [[-0.0026341094286689393, 4.822570156748241], [-0.0030026601843005513, 6.681170740991173]], 
# 14: [[-0.003506345432429798, 3.0528724234756015], [-0.0037990557859252407, 4.36027615636128]], 
# 15: [[-0.003655698586626174, 2.975215554164139], [-0.003963755150201819, 4.579327215778789]]
# }
@with_kw mutable struct DetectorQDC
    gain::Float64 = 30.79
    n_pmts::Int64 = 16
    grad_c::Vector{Float64} = [
        4.570725538613477, 3.4537377697297393, 4.907752978279222, 2.241773885191816, 
        2.864614759296568, 4.03886453341743, 2.893378071172577, 1.9486375219125216, 
        2.8794068605780905, 2.8794068605780905, 2.912107718286143, 6.996926904145113, 
        5.326732324708251, 4.822570156748241, 3.0528724234756015, 2.975215554164139
        ]
    grad_m::Vector{Float64} = [
        -0.0027726302119042904, -0.0029366127395134534, -0.0033613507916551142, -0.0037734805031322586, 
        -0.00307167001168234, -0.0030294520426235568, -0.0028199306577301486, -0.003168238958527181,
        -0.0027022291886093473, -0.0026675477947363285, -0.002515804915587131, -0.0020574023033302982, 
        -0.003092297887553726, -0.0026341094286689393, -0.003506345432429798, -0.003655698586626174,
    ]
    grad_c_bs::Vector{Float64} = [
        5.665745447865509, 4.3098087230480715, 5.7618482958693855, 2.879599787711115,
        3.8253692681421647, 5.109554456853292, 3.7311722430827845, 2.896290449097631,
        4.511791339287047, 4.418720738354981, 3.954548505472813, 8.24046943476849,
        6.786966101391472, 6.681170740991173, 4.36027615636128, 4.579327215778789,
    ]
    grad_m_bs::Vector{Float64} = [
        -0.002986950069349821, -0.003128027880529408, -0.0034747499158074486, -0.003967817605706616,
        -0.0033038002943305824, -0.0032875460101800913, -0.002999059160725512, -0.0033585689857879046,
        -0.0032894121389988643, -0.003217170817858253, -0.0029987469811330685, -0.0025855172841346046,
        -0.003422874313732227, -0.0030026601843005513, -0.0037990557859252407, -0.003963755150201819,
    ]
end


@with_kw mutable struct Event
    t_trig::Vector{Float64} = [-1., -1.]
    e_ind::Vector{Float64} = [-1., -1]
    e_raw::Float64 = -1.
    det::Int64 = -1
    n_sum::Int64 = 0
end


struct Source
    label::String
    rate::Float64
    energy::Float64
    width::Float64
    p_cor::Float64
    e_cor::Float64
    w_cor::Float64
    d_cor::Int64
end
