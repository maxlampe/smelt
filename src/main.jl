
include("structs.jl")
include("simulation.jl")
include("plotting.jl")
include("tools.jl")


bi_rate = 5e4
gain = 30.89
srcs = Dict(
    "bi1" => Source(bi_rate * 0.7, gain * 450., gain * 35.),
    "bi1_1" => Source(bi_rate * 0.2, gain * 530., gain * 35.),
    "bi1_2" => Source(bi_rate * 0.1, gain * 570., gain * 30.),
    "bi2" => Source(bi_rate * 4., gain * 1000., gain * 60.),
    "bi0" => Source(bi_rate * 1.2, gain * 70., gain * 15.),
)


res = run_single_det(srcs)
plot_e_hist(res)

println(dead_time_corr(length(res)))
println(bi_rate * (1. + 4. + 1.2))