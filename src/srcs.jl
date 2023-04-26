

include("structs.jl")

# [keV]

# obs data trigger rate 5200 / obs bg rate 530 [Hz]
# target obs rate 4670 [Hz]
cd_rate = 4.93e3
# srcs_cd = [
#     Source("cd1", cd_rate, 55., 17., 0.22, 14., 4., 0),
# ]


# obs data trigger rate 3200 / obs bg rate 530 [Hz]
# target obs rate 2670 [Hz]
ce_rate = 3.03e3
# srcs_ce = [
#     Source("ce0", ce_rate * 0.83, 106., 25., 0.22, 22., 5., 0),
#     Source("ce1", ce_rate * 0.15, 18., 7., 0., 0., 0., 0),
#     Source("ce2", ce_rate * 0.02, 40., 28., 0., 0., 0., 0),
# ]


# obs data trigger rate 2500 / obs bg rate 530 [Hz]
# target obs rate 1970 [Hz]
sn_rate = 2.e3
# srcs_sn = [
#     Source("sn1", sn_rate * 0.9, 340., 35., 0., 0., 0., 0),
#     Source("sn0", sn_rate * 0.1, 22., 10., 0., 0., 0., 0),
#     Source("snx1", sn_rate * 0.04, 65., 100., 0., 0., 0., 0),
#     Source("snx2", sn_rate * 0.03, 145., 100., 0., 0., 0., 0),
# ]
