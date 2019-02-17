module bEMD
using Dierckx, ElasticArrays
include("utils.jl")

function midpoints(s,t)
    maxs, mins, exts = find_extrema(s)
    mpx = Float64[]
    mpy = Float64[]
    for i=1:(length(exts)-1)
        push!(mpx, (t[exts[i]] + t[exts[i+1]]) /2.0)
        push!(mpy, (s[exts[i]] + s[exts[i+1]]) /2.0)
    end
    mpx,mpy
end

function get_spline(s,t)
    mpt, mps = midpoints(s,t)
    spl = Spline1D(mpt,mps,bc="extrapolate")
    return spl(t)
end

function get_spline_env(s,t)
    maxs, mins, = find_extrema(s)
    max_env = Spline1D(t[maxs],s[maxs])
    min_env = Spline1D(t[mins],s[mins])
    return (1/2)*(max_env(t)+min_env(t))
end

function bsift(s,t)
    h = ElasticArray{Float64}(undef,length(s),0)
    append!(h,s)
    count = 1
    sd = 1.0
    while sd > 0.25
        count > 5 && break
        spls = get_spline(h[:,end],t)
        append!(h, h[:,end]-spls)
        sd = sum(abs2.(h[:,end-1]-h[:,end]))/sum(h[:,end-1].^2)
        count +=1
    end
    return h[:,end]
end

function bEMD(s,t;maximfs=5)
    sig = copy(s)
    N=length(s)
    imfs = ElasticArray{Float64}(undef,N,0)
    for i in 1:maximfs
        h = bsift(sig,t)
        append!(imfs, h)
        sig -= h
    end
    return imfs
end
end
