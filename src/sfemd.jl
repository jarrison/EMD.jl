include("utils.jl")

""" Sift with preallocations"""
function sift!(Eav, decomp, d_osf, nsifts=5)
    N = length(decomp)
    e1 = zeros(N)
    e2 = zeros(N)
    avg = zeros(N)
    w = div(d_osf-1,2)
    for j in 1:nsifts
        stream_minmax(e1, e2, decomp, d_osf)
        for i in 1:N
            avg[i] = (e1[i] + e2[i]) / 2
        end
        avg .= moving_average(avg, d_osf)
        Eav .= Eav .+ avg
        decomp .= decomp .- avg
    end
    return nothing
end

function sift(signal, d_osf, nsifts=5)
    N = length(signal)
    e1 = zeros(N)
    e2 = zeros(N)
    avg = zeros(N)
    Eav = zeros(N)
    w = div(d_osf-1,2)
    decomp = copy(signal)

    for j in 1:nsifts
        stream_minmax(e1, e2, decomp, d_osf)
        avg = (e1 .+ e2) ./ 2
        savg = moving_average(avg,d_osf)
        Eav .=  Eav .+ savg
        decomp .= decomp .- savg
    end
    return Eav, decomp
end

""" Succint and Fast Empirical Mode Decomposition"""
function sfEMDallocate(signal;maximfs=10,nsifts=5)
    N = length(signal)
    a_new = copy(signal)
    imfs = zeros(N,maximfs)
    for i in 1:maximfs
        exs = find_extrema_count(a_new)
        d_osf =  div(N, div(exs,2))
        if iseven(d_osf)
            d_osf += 1
        end
        if d_osf >= div(N,3) || exs < 5
            imfs[:,i] .= a_new
            break
        end
        Eav, decomp = sift(a_new, d_osf,nsifts)
        imfs[:,i] .= decomp
        a_new .= Eav
    end
        return imfs
end

function sfEMD(signal;maximfs=10,nsifts=5)
    N = length(signal)
    a_new = copy(signal)
    imfs = zeros(N,maximfs)
    Eav = zeros(N)
    for i in 1:maximfs
        exs = find_extrema_count(a_new)
        d_osf =  div(N, div(exs,2))
        if iseven(d_osf)
            d_osf += 1
        end
        if d_osf >= div(N,3) || exs < 5
            imfs[:,i] .= a_new
            break
        end
        fill!(Eav,0)
        sift!(Eav, a_new, d_osf, nsifts)
        imfs[:,i] .= a_new
        a_new .= Eav
    end
        return imfs
end
