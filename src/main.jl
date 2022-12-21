

include("simulation.jl")
include("plotting.jl")


source_ind = 0
gain = 1. # [1/keV] or 30.89 [ch/keV]

# Bi trigger rate ca. 4500 Hz; Bg trigger rate ca 500 Hz
if source_ind == 0
    cd_rate = 5.4e3
    srcs = [
        Source("cd1", cd_rate, gain * 60., gain * 17.),
    ]
elseif source_ind == 1
    ce_rate = 3.3e3
    srcs = [
        Source("ce1", ce_rate * 0.8, gain * 116., gain * 17.),
        Source("ce0", ce_rate * 0.2, gain * 30., gain * 17.),
    ]
elseif source_ind == 2
    bi_rate = 4.7e3
    srcs = [
        Source("bi1", bi_rate * 0.7 / 7.0, gain * 450., gain * 35.),
        Source("bi1_1", bi_rate * 0.25 / 7.0, gain * 520., gain * 35.),
        Source("bi1_2", bi_rate * 0.1 / 7.0, gain * 570., gain * 30.),
        Source("bi2", bi_rate * 4. / 7.0, gain * 1000., gain * 60.),
        Source("bi0/bg", bi_rate * 1.95 / 7.0, gain * 70., gain * 19.),
        # Source("test", bi_rate, gain * 1000., gain * 10.),
    ]
elseif source_ind == 3
    sn_rate = 2.63e3
    srcs = [
        Source("sn1", sn_rate * 0.8, gain * 369., gain * 35.),
        Source("sn0/bg", sn_rate * 0.15, gain * 40., gain * 13.),
        Source("sn0/bg", sn_rate * 0.05, gain * 70., gain * 30.),
    ]
elseif source_ind == 4
    cs_rate = 7.9e3
    srcs = [
        Source("cs1", cs_rate * 0.1, gain * 620., gain * 40.),
        Source("cs0_0", cs_rate * 0.35, gain * 40., gain * 50.),
        Source("cs0_1", cs_rate * 0.22, gain * 150., gain * 60.),
        Source("cs0_2", cs_rate * 0.12, gain * 250., gain * 60.),
        Source("cs0_3", cs_rate * 0.11, gain * 380., gain * 60.),
    ]
end

res, stat = run_sim(srcs; n_det=2, meas_time=60.0, with_bs=true, with_ang=true)
print(stat)
plot_e_hist(res, gain)
# plot_e_ind_hist(res, gain)
plot_diff_both_trig(res)
# plot_e_acc_hist(res, gain)
plot_e_both_trig(res, gain)
plot2D_e_both_trig(res, gain)
# plot2D_e_ind_both_trig(res, gain)
plot2D_e_ind(res, gain)

# res2, stat2 = run_sim(srcs; n_det=2, meas_time=5.0, with_bs=true, with_ang=false)
# plot_e_hist(res2, gain)
# plot_diff_both_trig(res2)
# plot2D_e_ind(res2, gain)
# plot2D_e_ind_both_trig(res2, gain)