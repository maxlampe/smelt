

using Statistics
using Plots
include("tools.jl")
include("simulation.jl")


function n_turns(en::Float64, dist::Float64 = 3.2)
    return dist / (cos(create_ang() * pi / 180.) * r_gyr(en))
end


function e_loss_turn(en::Float64, b::Float64 = 0.15)
    return 2.65e4 * b *  (en * 1e-6)^3 
end


function synch_loss(en::Float64)
    return e_loss_turn(en) * n_turns(en)
end


# mean E loss due to synchroton radiation for 1 MeV electron in Perkeo
es = LinRange(1., 1500., 30)
losses = [mean([synch_loss(e) for i=1:1000]) * 1.1 for e in es]

plot(es, losses, title="Synchroton Loss (PERKEO III)", label="Loss")
ylabel!("Energy Loss [eV]")
xlabel!("Energy [keV]")
savefig("synch_loss")
