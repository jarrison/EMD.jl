module sEMDT
using Base.Threads,Distributed, SharedArrays
function myextrema(ret::AbstractVector{T},X) where T
           max = X[1]
           min = X[1]
           @inbounds for i = 2:lastindex(X)
               if isnan(X[i])
                   return NaN
               end
               if X[i] > max
                   max = X[i]
               end
               if X[i] < min
                   min = X[i]
               end
           end
           ret[1] = max
           ret[2] = min
           nothing
end


function myargextrema(ret::AbstractVector{T},X) where T
           max = X[1]
           min = copy(X[1])
           maxloc = 1
           minloc = 1
           @inbounds for i = 2:lastindex(X)
               if isnan(X[i])
                   return NaN
               end
               if X[i] > max
                   maxloc = i
               end
               if X[i] < min
                   minloc = i
               end
           end
           ret[1] = maxloc
           ret[2] = minloc
           nothing
end


function mirror(A,d::Int=3)
    w = div(d-1,2)

    a = reverse(view(A,2:w+1))
    b = A
    c = copy(A[length(A)-w:length(A)-1])
    reverse!(c)
    ret = [a;b;c]
    return ret
end
export mirror


function find_extrema(A,d::Int=3)

    signal = mirror(A,3)
    n = length(A)
    maxs = Int[]
    mins = Int[]

    @inbounds for i in 1:lastindex(A)
        vv = view(signal,i:i+2)
        amax = argmax(vv)
        amin = argmin(vv)
        if amax == 2
            append!(maxs,i)
        elseif amin == 2
            append!(mins,i)
        else
            continue
        end
    end

    return maxs, mins
end
export find_extrema

function _allrelextrema(A, tol=0.0)
    y = mirror(A)
    # compute difference between successive values (like the slope)
    slope = diff(y)
    slope[abs.(slope) .< tol::Float64] .= 0.0

    # remove all zeros while tracking original indices
    nonzero = (slope .!= 0.0)
    slope = slope[nonzero]
    indices = 1:lastindex(slope)
    indices = indices[nonzero]

    # we just want the sign of the slope
    slope_sign = zeros(Int8,length(slope))
    slope_sign[slope .> 0] .= 1
    slope_sign[slope .< 0] .= -1

    # so that we can find the sign of the curvature at points with differing
    # slope signs to either side
    curve_sign = diff(slope_sign)
    arg_curve_chng = findall(!iszero,curve_sign)
    i0 = indices[arg_curve_chng]
    i1 = indices[arg_curve_chng] .+ 1
    i = div.((i0 .+ i1) , 2)

    return i
end
export _allrelextrema

function get_envelopes(A,d_osf::Int=3)
    mir = mirror(A,d_osf)
    env1 = similar(A)
    env2 = similar(A)
    p = [Vector{Float64}(undef,2) for i in 1:nthreads()]

    @threads for i in 1:lastindex(A)
        # println( threadid())
        win = view(mir,i:d_osf+i-1)

        myextrema(p[threadid()],win)
        @inbounds env1[i] = p[threadid()][1]
        @inbounds env2[i] = p[threadid()][2]
    end
    env1, env2
end
export get_envelopes

function ma(A,d_osf::Int=3)
    w = div(d_osf-1,2)
    B = similar(A)
    cumsum!(B,A)
    B[d_osf:end] .- B[1:(end-d_osf+1)] ./ d_osf
end
export ma

function sift(A,d_osf::Int,n_sifts::Int)
    Eav = similar(A)
    A_new = copy(A)
    fill!(Eav,0.0)

    for i in 1:n_sifts

        env1, env2 = get_envelopes(A_new,d_osf)

        senv1 = moving_average(env1,d_osf)

        senv2 = moving_average(env2,d_osf)
        avg = (senv1.+senv2)./2.0
        savg = moving_average(avg,d_osf)

        Eav += savg
        A_new .-= savg
    end
    return Eav, A_new
end
export sift

function sEMD(A::Array{Float64};max_sifts::Int=5,max_imfs::Int=5)
    imfs = []
    A_new = copy(A)
    for i in 1:max_imfs

        exs = _allrelextrema(A_new)

        #calculate d_osf
        neg= div(length(exs),2)
        if neg <= 5
            pushfirst!(imfs,A_new)
            break
        end
        d_osf = div(length(A_new),neg)

        if d_osf > div(length(A_new),3)
            pushfirst!(imfs,A_new)
            break
        end
        if  d_osf & 1 != true
            d_osf -= 1
        end
        Eav, A_new = sift(A_new,d_osf,max_sifts)
        pushfirst!(imfs,A_new)
        A_new = Eav

    end
    return imfs
end
export sEMD
end
