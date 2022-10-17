mutable struct Detector
    t_meas::Float64
    t_dead::Float64
    thresh::Float64
    is_ready::Bool
    is_measuring::Bool
end

mutable struct Event
    t_trig::Float64
    e_det::Float64
    det::Int8
end

struct Source
    rate::Float64
    energy::Float64
    width::Float64
end
