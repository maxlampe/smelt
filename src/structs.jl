
using Parameters

@with_kw mutable struct Detector
    t_meas::Float64 = 3e-7
    t_dead::Float64 = 1.5e-6
    thresh::Float64 = 1.0
    is_ready::Bool= true
    is_measuring::Bool = false
end


@with_kw mutable struct Event
    t_trig::Vector{Float64} = [-1., -1.]
    e_ind::Vector{Float64} = [-1., -1]
    det::Int64 = -1
    n_sum::Int64 = 0
end


struct Source
    label::String
    rate::Float64
    energy::Float64
    width::Float64
end
