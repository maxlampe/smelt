

include("simulation.jl")
include("plotting.jl")


# Bi trigger rate ca. 4500 Hz; Bg trigger rate ca 500 Hz
bi_rate = 4.6e3
gain = 1. # [1/keV] or 30.89 [ch/keV]
srcs = [
    # Source("bi1", bi_rate * 0.7 / 7.0, gain * 450., gain * 35.),
    # Source("bi1_1", bi_rate * 0.25 / 7.0, gain * 520., gain * 35.),
    # Source("bi1_2", bi_rate * 0.1 / 7.0, gain * 570., gain * 30.),
    # Source("bi2", bi_rate * 4. / 7.0, gain * 1000., gain * 60.),
    # Source("bi0/bg", bi_rate * 1.95 / 7.0, gain * 70., gain * 19.),
    Source("test", bi_rate, gain * 1000., gain * 10.),
]

res, stat = run_sim(srcs; n_det=1, meas_time=0.5, with_bs=true, t_meas=1e-8)
print(stat)
plot_e_hist(res, gain)
plot_e_ind_hist(res, gain)
plot_diff_both_trig(res)
plot_e_acc_hist(res, gain)
plot_e_both_trig(res, gain)
plot2D_e_both_trig(res, gain)
plot2D_e_ind_both_trig(res, gain)
plot2D_e_ind(res, gain)

res2, stat2 = run_sim(srcs; n_det=2, meas_time=8.0, with_bs=true)
plot_e_hist(res2, gain)
plot_e_ind_hist(res2, gain)
plot_diff_both_trig(res2)
plot_e_acc_hist(res2, gain)
plot_e_both_trig(res2, gain)
plot2D_e_both_trig(res2, gain)
plot2D_e_ind_both_trig(res2, gain)
plot2D_e_ind(res2, gain)

# TODO check 2k peak likelihood with BS