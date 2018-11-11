include("utils.jl")

function sift(signal, d_osf, nsifts=5)
    N = length(signal)
    e1 = zeros(N+d_osf-1)
    e2 = zeros(N+d_osf-1)
    avg = zeros(N)
    Eav = zeros(N)
    w = div(d_osf-1,2)
    decomp = copy(signal)
    
    for j in 1:nsifts
        e1,e2 = stream_minmax(decomp,d_osf)
        avg = (e1 .+ e2) ./ 2.0
        savg = moving_average(avg,d_osf)
        Eav .= Eav .+ savg
        decomp .= decomp .- savg
    end
    return Eav, decomp
end

""" Succint and Fast Empirical Mode Decomposition"""
function sfEMD(signal;maximfs=10,nsifts=5)
    N = length(signal)
    a_new = copy(signal)
    imfs = zeros(N,maximfs)
    for i in 1:maximfs
        exs = find_extrema_count(a_new)
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
