

include("simulation.jl")
include("plotting.jl")
include("srcs.jl")


res, stat = run_sim(srcs_sn; meas_time=1., sim_step=1e-8, with_bs=true, with_ang=true, with_bs_var=true, verbose=true)
print(stat)
plot_e_hist(res; e_max=800.)
plot_diff_both_trig(res)
plot2D_e_both_trig(res; e_max=0.6e3)
plot2D_e_ind_both_trig(res; e_max=800.)
plot_e_hist_tofcut(res; e_max=800.)

# -----------------------------------------------------------------

res, stat = run_sim(srcs_cd; meas_time=5., sim_step=1e-8, with_bs=true, with_ang=true, with_bs_var=true, verbose=true)
print(stat)
plot_e_hist(res; e_max=200.)
plot_diff_both_trig(res)
plot2D_e_both_trig(res; e_max=0.3e3)
plot2D_e_ind_both_trig(res; e_max=200.)
plot_e_hist_tofcut(res; e_max=200.)

# -----------------------------------------------------------------

res, stat = run_sim(srcs_ce; meas_time=5., sim_step=1e-8, with_bs=true, with_ang=true, with_bs_var=true, verbose=true)
print(stat)
plot_e_hist(res; e_max=300.)
plot_diff_both_trig(res)
plot2D_e_both_trig(res; e_max=0.3e3)
plot2D_e_ind_both_trig(res; e_max=300.)
plot_e_hist_tofcut(res; e_max=300.)

# -----------------------------------------------------------------

# bi_rate = 4.7e3
# srcs = [
#     Source("bi1", bi_rate * 0.7 / 7.0, gain * 450., gain * 35., 0.0, 0., 0., 0),
#     # Source("bi1_1", bi_rate * 0.25 / 7.0, gain * 520., gain * 35., 0.0, 0., 0., 0),
#     # Source("bi1_2", bi_rate * 0.1 / 7.0, gain * 570., gain * 30., 0.0, 0., 0., 0),
#     Source("bi1_1", bi_rate * 0.25 / 7.0, gain * 520., gain * 35., 0.05, 50., 20., 0),
#     Source("bi1_2", bi_rate * 0.1 / 7.0, gain * 570., gain * 30., 0.05, 50., 20., 0),
#     # Source("bi2", bi_rate * 4. / 7.0, gain * 1000., gain * 60., 0.03, 570., 45., 0),
#     Source("bi2", bi_rate * 2.05 / 7.0, gain * 1000., gain * 60., 0.05, 570., 45., 0),
#     Source("bi2", bi_rate * 2.05 / 7.0, gain * 1000., gain * 60., 0.03, 50., 20., 0),
#     Source("bi0/bg", bi_rate * 1.85 / 7.0, gain * 70., gain * 19., 0.0, 0., 0., 0),
#     # Source("test", bi_rate, gain * 1000., gain * 10.),
# ]

# res, stat = run_sim(srcs; n_det=2, meas_time=60., with_bs=true, with_ang=true, with_bs_var=true)
# print(stat)
# plot_e_hist(res; gain=gain)
# plot_e_hist_tofcut(res; gain=gain, e_max=2000., e_bin=20.)
# # plot_e_ind_hist(res; gain=gain)
# plot_diff_both_trig(res)
# # plot_e_acc_hist(res; gain=gain)
# # plot_e_both_trig(res; gain=gain)
# plot2D_e_both_trig(res; gain=gain, e_max=2.e3)
# plot2D_e_ind_both_trig(res; gain=gain, e_max=1500., e_bin=25., e_min=50.)
# plot2D_e_ind(res; gain=gain, e_max=1500., e_bin=25., e_min=0.)
