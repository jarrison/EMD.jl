include("gEMD.jl")
using Interpolations
using Statistics

function midpoints(s)
    nex, exts = gEMD.find_extrema(s)
    extvals = s[exts]
    midpts = exts[2:end] .-div.(diff(exts),2)
    midpts = convert(Array{Int64},midpts)
    push!(midpts,length(s))
    pushfirst!(midpts,1)
    mvals = s[midpts]

    return convert(Array{Int64},midpts), mvals
end

function get_spline(s)
    mids,mvals = midpoints(s)
    N = length(s)
    itp = interpolate(mvals,BSpline(Cubic(Periodic(OnGrid()))))
    extrapolate(itp,Reflect())
    sitp = scale(itp, range(mids[1],stop=mids[end],length=length(mids)))
    sig = zeros(length(s))

    for i in 1:N
        sig[i] = sitp(i)
    end
    return sig
end

function bsift(s)
    m = copy(s)
    sd = 1.0
    while sd > 0.25
        avg =  get_spline(m)
        h = m - avg
        sd = sum(abs2,(m.-h)) / sum(m.^2)
        @show sd
        m = h
    end
    return m
end

function bEMD(s)
    sig = copy(s)
    imfs = zeros(length(s),5)
    for i in 1:2
        h = bsift(sig)
        imfs[:,i] = h
        sig -= h
    end
    return imfs
end
