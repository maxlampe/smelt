
using HypergeometricFunctions
using Plots


estar_nentries = 41
# dE / dx from ESTAR, (energy [keV], dE/dx [keV/100nm])
estar_data = [
  [1., 0.0126936], [1.223, 0.0111146], [1.495, 0.00968532], [1.827, 0.00841286],
  [2.234, 0.00728489], [2.732, 0.00629107], [3.34, 0.00542006],
  [4.083, 0.00466051], [4.992, 0.00400003], [6.103, 0.0034283],
  [7.462, 0.00293398], [9.123, 0.00250982], [11.15, 0.00214553],
  [13.64, 0.00183283], [16.67, 0.00156658], [20.38, 0.00133954],
  [24.92, 0.00114655], [30.47, 0.000983083], [37.25, 0.000844486],
  [45.54, 0.000727354], [55.68, 0.000628694], [68.07, 0.000545722],
  [83.22, 0.000476165], [101.7, 0.000418166], [124.4, 0.000369766],
  [152.1, 0.000329827], [185.9, 0.000296906], [227.3, 0.000270178],
  [277.9, 0.000248506], [339.8, 0.000231374], [415.5, 0.000217752],
  [507.9, 0.000207019], [621., 0.00019897], [759.2, 0.000193087],
  [928.2, 0.000189166], [1135., 0.000186792], [1387., 0.000185657],
  [1696., 0.000185657], [2074., 0.000186586], [2536., 0.000188237],
  [3100., 0.000190404]
]


function getBirks(energy::Float64, estar_i, estar_k, estar_y)
    tar_i = 1
    # find correct interpolation point
    for i in (2:estar_nentries)
        if estar_data[i][1] > energy
            tar_i = i
            break
        end
    end
    if tar_i == 1
        return EstarIntegral(0, energy, estar_k[1], estar_y[1])
    else
        return estar_i[tar_i - 1] + EstarIntegral(estar_data[tar_i - 1][1], energy, estar_k[tar_i - 1], estar_y[tar_i - 1])
    end
end


function EstarIntegral(e0::Float64, e1::Float64, k::Float64, y::Float64)
    if (e0 >= e1 || e1 <= 0.)
        return 0.
    elseif (y == 0.)
        return (e1 - e0)/(1. + k)
    end
    
    invy = 1. / y
    i0 = 0.
    if e0 <= 0.
        if y < 0.
            i0 = k ^ (-invy) * pi * invy / (sin(pi * invy))
        else
            i0 = 0.
        end
    else
        i0 = e0 * _₂F₁(1., invy, 1. + invy, -k * e0 ^ y)
    end

    i1 = e1 * _₂F₁(1., invy, 1. + invy, -k * e1 ^ y)
    return i1 - i0;
end


function cacheIntegral6(k_B::Float64)
    estar_y = Vector{Float64}(undef, estar_nentries)
    estar_k = Vector{Float64}(undef, estar_nentries)
    estar_i = Vector{Float64}(undef, estar_nentries)

    # update values for log-log linear interpolation
    for i in (2:estar_nentries)
        estar_y[i - 1] = (log(estar_data[i][2]) - log(estar_data[i - 1][2])) / (log(estar_data[i][1]) - log(estar_data[i - 1][1]))
        estar_k[i - 1] = k_B * estar_data[i - 1][2] * estar_data[i - 1][1] ^ (-estar_y[i - 1])
    end

    #update cached integrals
    estar_i[1] = EstarIntegral(0., estar_data[1][1], estar_k[1], estar_y[1])
    for i in (2:estar_nentries)
        estar_i[i] = estar_i[i - 1] + EstarIntegral(estar_data[i - 1][1], estar_data[i][1], estar_k[i - 1], estar_y[i - 1])
    end

    return estar_i, estar_k, estar_y
 end


function calc_birks(kB::Float64 = 150., with_plot::Bool = false)

    estar_i, estar_k, estar_y = cacheIntegral6(kB)
    x = range(1, 1100, length=100)
    y = [(getBirks(x_i, estar_i, estar_k, estar_y)) for x_i in x]

    if with_plot
        plot(x, y)
        ylims!(0.7, 1.1)
    end

    return x, y
end 


function f_birks(x_en::Vector{Float64}, par::Vector{Float64})
    k_B = par[1]
    k_fac = par[2]
    estar_i, estar_k, estar_y = cacheIntegral6(k_B)

    y = [getBirks(xx, estar_i, estar_k, estar_y) for xx in x_en]
    y = ((y.* k_fac)./ x_en)

    return y
end 


function f_birks_150(x_en::Vector{Float64}, par::Vector{Float64})
    k_B = par[1]
    k_fac = par[2]
    estar_i, estar_k, estar_y = cacheIntegral6(k_B)

    y = [getBirks(xx, estar_i, estar_k, estar_y) for xx in x_en]
    y = ((y.* k_fac)./ x_en)

    return y.* f_birks(x_en, [150., 1.])
end 