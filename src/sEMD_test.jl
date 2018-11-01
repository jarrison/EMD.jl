using BenchmarkTools,Base.Threads, Distributed
include("Semd.jl")
const τ = 2π # Glory to the τ overlords
function mirror(A,d::Int=3)
    w = div(d-1,2)

    a = reverse(view(A,2:w+1))
    b = A
    c = copy(A[length(A)-w:length(A)-1])
    reverse!(c)
    ret = [a;b;c]
    return ret
end

"""
Fast Implementation of a minmax order statistics filter.
"""
function stream_minmax!(env1::Array{Float64},env2::Array{Float64},a,w::Int64)
    upper = Int[]
    lower = Int[]

    pushfirst!(upper,1)
    pushfirst!(lower,1)

    for i in 2:lastindex(a)
        if i >= w
            @views env1[i] = a[upper[1]]
            @views env2[i] = a[lower[1]]
        end
        if a[i] > a[i-1]
            pop!(upper)
            while isempty(upper) != true && a[i] > a[upper[end]]
                pop!(upper)
            end
        elseif a[i] <= a[i-1]
            pop!(lower)
            while isempty(lower) != true && a[i] < a[lower[end]]
                pop!(lower)
            end
        end
        push!(upper,i)
        push!(lower,i)
        if i == w + upper[1]
            popfirst!(upper)
        elseif i == w + lower[1]
            popfirst!(lower)
        end
    end
        nothing
end

"""
Moving Average using cumulative sum algorithm (see wikipedia if you care
why this works).
"""
function moving_average2(A,d_osf::Int=3)
    sum=0
    result = similar(A)

    for i in 1:d_osf
        sum+= A[i]
        result[i] = sum / i
    end

    for i in d_osf:lastindex(A)
        sum = sum - A[i-d_osf+1] + A[i]
        result[i] = sum/d_osf
    end
    result
end

function moving_average(A,d_osf::Int=3)
    mirr = mirror(A,d_osf)
    ret = similar(mirr)
    cumsum!(ret,mirr)
    s = (ret[d_osf:end] .- ret[1:end-d_osf+1]) ./ d_osf

end

function fast_sift(A,d_osf::Int,n_sifts::Int)
    Eav = similar(A)
    A_new = copy(A)
    fill!(Eav,0.0)
    mir = mirror(A_new,d_osf)
    env1 = similar(mir)
    env2 = similar(mir)
    for i in 1:n_sifts
        mir = mirror(A_new,d_osf)
        stream_minmax!(env1,env2,mir,d_osf)
    end
    return env1[(d_osf):end], env2[(d_osf):end]
end

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

t = LinRange(0, 1e-3, convert(Int, 1e4))
const a = 1sin.(75e3*t) .+ .01sin.(3*t)

function EMD(signal)

    imfs = []

    new_a = copy(signal)
    for j in 1:3
        exs = _allrelextrema(new_a)
        Eav = zeros(length(signal))
        #calculate d_osf
        neg= div(length(exs),2)
        if neg <= 5
            pushfirst!(imfs,new_a)
            break
        end
        d_osf = div(length(new_a),neg)

        if d_osf > div(length(new_a),3)
            pushfirst!(imfs,new_a)
            break
        end
        if  d_osf & 1 != true
            d_osf -= 1
        end
        env1 = zeros(length(signal)+d_osf-1)
        env2 = zeros(length(signal)+d_osf-1)
        for i in 1:3
            mir = mirror(new_a,d_osf)
            stream_minmax!(env1, env2, mir, d_osf)
            avg = ((env1+env2)/2)[d_osf:end]
            savg = moving_average(avg,div(d_osf,2))
            new_a -= savg
            Eav += savg
        end
        push!(imfs,new_a)
        new_a = Eav
    end
    imfs
end
imfs =  EMD(a)

plot(e)
plot!(imf)
gui()
