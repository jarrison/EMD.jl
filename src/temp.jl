include("gEMD.jl")
using Plots,BenchmarkTools

t = LinRange(0,1,10^7)
const signal = 3sin.(2π*8*t) .+ sin.(2π*4*t)
Eav = zeros(length(t))

function sift!(Eav, signal, d_osf, nsifts=3)
    N = length(signal)
    e1 = zeros(N+d_osf-1)
    e2 = zeros(N+d_osf-1)
    avg = zeros(N)
    decomp = copy(signal)
    decomp2 = copy(signal)
    mirror_buffer = copy(signal)
    w = div(d_osf-1,2)

    for j in 1:nsifts
        gEMD.mirror!(mirror_buffer, decomp,d_osf)
        gEMD.stream_minmax!(e1,e2,mirror_buffer,d_osf)
        avg .= (e1[d_osf:end] .+ e2[d_osf:end]) ./ 2.0
        savg = gEMD.moving_average(avg,d_osf)
        Eav .= Eav .+ savg
        decomp = decomp2 - savg
        fill!(resize!(avg,N),0.0)
    end
    return decomp
end

using Profile
@btime sift!($Eav,$signal,111,3)
# Juno.profiler()
