function mirror(A,d::Int=3)
    w = div(d-1,2) #width of window, d is an ODD integer
    mA = [A[(w+1):-1:2];A;A[lastindex(A)-1:-1:(length(A)-w)]]
end

function find_extrema(A::Array,d::Int=3)
    mA = mirror(A)
    count = 0
    maxs = Int[]
    mins = Int[]
    exts = Int[]
    for i in 1:lastindex(mA)-d+1
        win = view(mA,i:i+d-1)
        if (win[2] > win[1] && win[2] > win[3])
                count+=1
                push!(maxs,i)
                push!(exts,i)
        elseif (win[2] < win[1] && win[2] < win[3])
                count+=1
                push!(mins,i)
                push!(exts,i)
        end
    end
    maxs, mins,exts
end

function find_extrema_count(A::Array,d::Int=3)
    dd = diff(sign.(diff(A)))
    return length(dd[dd .!= 0]) + 2
end

function find_extrema_minmax(A::Array,d::Int=3)
    dd = diff(sign.(diff(A)))
    return findall(dd .< 0), findall(dd .> 0)
end

function stream_minmax(env1, env2, A, d::Int)
    a = mirror(A, d)
    N = length(a)
    upper = UInt32[] #buffer for indices
    lower = UInt32[]
    push!(upper,1)
    push!(lower,1)
    w = d
    for i in 2:lastindex(a)
        if i >= w
            env1[i-w+1] = a[upper[1]]
            env2[i-w+1] = a[lower[1]]
        end

        if a[i] > a[i-1]
            pop!(upper)
            # remove maxima from buffer that are less than our new one
            while isempty(upper) != true
                if a[i] <= a[upper[end]]
                    break
                end
                pop!(upper)
            end
        else
            pop!(lower)
            # remove minima from buffer that are greater than our new one
            while isempty(lower) != true
                if a[i] >= a[lower[end]]
                    break
                end
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
        env1[end-w] = a[upper[1]]
        env2[end-w] = a[lower[1]]
    end
    return nothing
end

function moving_average(A, d::Int=3)
    A = mirror(A,d)
    T = typeof(one(eltype(A))/1)
    ret = Vector{T}(undef, length(A) - d + 1)
    id = 1 / d
    s = sum(view(A, 1:d))
    ret[1] = s * id
    @inbounds for n = 1:length(ret)-1
        s += A[n+d] - A[n]
        ret[n+1] = s * id
    end
    return ret
end
