function sift!(Eav, decomp, d_osf, nsifts=5)
    N = length(decomp)
    e1 = zeros(N)
    e2 = zeros(N)
    avg = zeros(N)
    w = max(div(d_osf-1, 2), 3)
    if iseven(w)
        w += 1
    end
    for j in 1:nsifts
        stream_minmax(e1, e2, decomp, d_osf)
        e1 .= moving_average(e1, d_osf)
        e2 .= moving_average(e2, d_osf)
        avg .= (e1 .+ e2)./2
        avg .= moving_average(avg, d_osf)
        Eav .= Eav .+ avg
        decomp .-= avg
    end
    return nothing
end

function sEMD(signal;maximfs=10,nsifts=5)
    N = length(signal)
    a_new = copy(signal)
    imfs = ElasticArray{Float64}(undef,N,0)
    Eav = zeros(N)
    for i in 1:maximfs
        exs = find_extrema_count(a_new)
        d_osf =  div(N, div(exs,2))
        iseven(d_osf) && (d_osf+=1)

        if d_osf >= div(N,3) || exs < 5
            append!(imfs, a_new)
            break
        end
        fill!(Eav,0)
        sift!(Eav, a_new, d_osf, nsifts)
        append!(imfs, a_new)
        a_new .= Eav
    end
        return imfs
end
