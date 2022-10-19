
# include("structs.jl")
include("simulation.jl")
include("plotting.jl")
include("tools.jl")

# Bi trigger rate ca. 4500 Hz; Bg trigger rate ca 500 Hz
bi_rate = 4.5e3
gain = 1. # [1/keV] or 30.89 [ch/keV]
srcs = Dict(
    "bi1" => Source(bi_rate * 0.7, gain * 450., gain * 35.),
    "bi1_1" => Source(bi_rate * 0.25, gain * 520., gain * 35.),
    "bi1_2" => Source(bi_rate * 0.1, gain * 570., gain * 30.),
    "bi2" => Source(bi_rate * 4., gain * 1000., gain * 60.),
    "bi0" => Source(bi_rate * 1.2, gain * 70., gain * 15.),
)


res = run_sim(srcs; n_det=2, meas_time=10.)
plot_e_hist(res, gain)

println(dead_time_corr(length(res)))
println(bi_rate * (1. + 4. + 1.2))