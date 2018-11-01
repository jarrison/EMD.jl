include("gEMD.jl")
using Plots,BenchmarkTools
gr()

function sift(signal, d_osf, nsifts=3)
    N = length(signal)
    e1 = zeros(N+d_osf-1)
    e2 = zeros(N+d_osf-1)
    avg = zeros(N)
    Eav = zeros(N)
    w = div(d_osf-1,2)
    decomp = copy(signal)
    for j in 1:nsifts
        e1,e2 = gEMD.stream_minmax(decomp,d_osf)
        avg = (e1 .+ e2) ./ 2.0
        savg = gEMD.moving_average(avg,d_osf)
        Eav .= Eav .+ savg
        decomp .= decomp .- savg
    end
    return Eav, decomp
end

function main(signal;maximfs=5,nsifts=5)
    N = length(signal)
    a_new = copy(signal)
    imfs = zeros(N,maximfs)

    for i in 1:maximfs
        exs = gEMD.find_extrema_count(a_new)
        d_osf =  div(N,div(exs,2))
        if iseven(d_osf)
            d_osf += 1
        end
        if d_osf >= div(N,3) || exs < 5
            imfs[:,i] .= a_new
            break
        end

        Eav, snew = sift(a_new,d_osf,nsifts)
        imfs[:,i] .= snew
        a_new .= Eav
    end
        return imfs
end
using Interpolations

using Base.Threads,Statistics

function get_spline(s)
    mids = gEMD.midpoints(s)
    N = length(s)
    @show mids[1]
    itp = interpolate(s[mids],BSpline(Cubic(Periodic(OnGrid()))))
    extrapolate(itp,Reflect())
    sitp = scale(itp, range(mids[1],stop=mids[end],length=length(mids)))
    sig = zeros(length(s))

    for i in 1:N
        sig[i] = sitp(i)
    end
    return sig
end

function bEMD(s)
    sig = copy(s)
    m = copy(s)
    sd = 1.0
    imfs = zeros(length(s),5)
    for i in 1:3
        while sd > 0.25
            avg =  get_spline(m)
            h = m .- avg
            sd = sum(abs2,(m.-h))/ sum(m.^2)
            m .-= avg
        end
        imfs[:,i] = m
        m .= sig .- m
    end
    return imfs
end


t = LinRange(0,1,10^3)
s = 3sin.(2π*8*t) .+ sin.(2π*4*t)+ .05*randn(length(t))
mids = gEMD.midpoints(s)

sig = get_spline(s)
plot(s)
plot!(mids,s[mids])
