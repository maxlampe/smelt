using BenchmarkTools
# include("structs.jl")
include("simulation.jl")
# include("tools.jl")

bi_rate = 4.6e3
gain = 1.
srcs = [
    Source("bi1", bi_rate * 0.7 / 7.0, gain * 450., gain * 35.),
    Source("bi1_1", bi_rate * 0.25 / 7.0, gain * 520., gain * 35.),
    Source("bi1_2", bi_rate * 0.1 / 7.0, gain * 570., gain * 30.),
    Source("bi2", bi_rate * 4. / 7.0, gain * 1000., gain * 60.),
    Source("bi0/bg", bi_rate * 1.95 / 7.0, gain * 70., gain * 19.),
]

# res, stat = run_sim(srcs; n_det=2, meas_time=0.1, with_bs=true)
@btime run_sim($srcs; n_det=2, meas_time=0.1, with_bs=true, verbose=false, sim_step=1e-8)
@code_warntype run_sim(srcs; n_det=2, meas_time=0.1, with_bs=true)
