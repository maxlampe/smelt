

using BenchmarkTools
include("simulation.jl")
include("srcs.jl")


# Runtime checks 
# 3.637s / 5.36 Gb
@btime run_sim($srcs_ce; meas_time=0.1, sim_step=5e-9, with_bs=false, with_ang=false, with_bs_var=false, verbose=false)
# 3.637s / 5.36 Gb
@btime run_sim($srcs_ce; meas_time=0.1, sim_step=5e-9, with_bs=true, with_ang=true, with_bs_var=true, verbose=false)
# 1.831s / 2.68 Gb
@btime run_sim($srcs_ce; meas_time=0.1, sim_step=1e-8, with_bs=true, with_ang=true, with_bs_var=true, verbose=false)
# 3.131s / 3.13 Gb (verbose makes it significantly slower!)
@btime run_sim($srcs_ce; meas_time=0.1, sim_step=1e-8, with_bs=true, with_ang=true, with_bs_var=true, verbose=true)


# Type checks
@code_warntype run_sim(srcs_ce; meas_time=0.1, with_bs=true, with_ang=true, with_bs_var=true)
