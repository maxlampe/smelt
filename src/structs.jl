

using Parameters


@with_kw mutable struct Detector
    t_meas::Float64 = 3e-7
    t_dead::Float64 = 1.5e-6
    thresh::Float64 = 1.0
    is_ready::Bool= true
    is_measuring::Bool = false
end


@with_kw mutable struct DetectorQDC
    gain::Float64 = 30.9
    n_pmts::Int64 = 16
    grad_c::Vector{Float64} = [4., 4., 5., 2., 3., 4., 3., 2., 3., 3., 4., 7., 6., 5., 3., 3.]
    grad_m::Vector{Float64} = [
        0.00321, 0.00340, 0.00360, 0.00385, 0.00348, 0.00327, 0.00307, 0.00349,
        0.00291, 0.00291, 0.00327, 0.00267, 0.00358, 0.00281, 0.00364, 0.00390,
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
