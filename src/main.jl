
# include("structs.jl")
include("simulation.jl")
include("plotting.jl")
include("tools.jl")

# Bi trigger rate ca. 4500 Hz; Bg trigger rate ca 500 Hz
bi_rate = 4.6e3
gain = 1. # [1/keV] or 30.89 [ch/keV]
srcs = Dict(
    "bi1" => Source(bi_rate * 0.7 / 7.0, gain * 450., gain * 35.),
    "bi1_1" => Source(bi_rate * 0.25 / 7.0, gain * 520., gain * 35.),
    "bi1_2" => Source(bi_rate * 0.1 / 7.0, gain * 570., gain * 30.),
    "bi2" => Source(bi_rate * 4. / 7.0, gain * 1000., gain * 60.),
    "bi0/bg" => Source(bi_rate * 1.95 / 7.0, gain * 70., gain * 19.),
)

res, stat = run_sim(srcs; n_det=2, meas_time=20.0, with_bs=true)
plot_e_hist(res, gain)
plot_diff_both_trig(res)
plot_e_acc_hist(res, gain)
plot_e_both_trig(res, gain)
plot2D_e_both_trig(res, gain)

res2, stat2 = run_sim(srcs; n_det=2, meas_time=20.0, with_bs=false)
plot_e_hist(res2, gain)
plot_diff_both_trig(res2)
plot_e_acc_hist(res2, gain)
plot_e_both_trig(res2, gain)
plot2D_e_both_trig(res2, gain)

# using BenchmarkTools

# res, stat = run_sim(srcs; n_det=2, meas_time=0.1, with_bs=true)
# @btime run_sim(srcs; n_det=2, meas_time=0.1, with_bs=true)
# @code_warntype run_sim(srcs; n_det=2, meas_time=0.1, with_bs=true)
