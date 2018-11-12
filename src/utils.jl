function mirror!(A::Array,d::Int=3)
    w = div(d-1,2) #width of window, d is an ODD integer
    prepend!(A,zeros(w))
    append!(A,zeros(w))

    @views A[1:w] = reverse(A[w+2:2w+1])
    @views A[end-w+1:end] = reverse(A[end-2w:end-w-1])
    return nothing
end

function mirror(A::AbstractArray,d::Int=3)
    mir = copy(A)
    w = div(d-1,2) #width of window, d is an ODD integer
    prepend!(mir,zeros(w))
    append!(mir,zeros(w))
    mir[1:w] = reverse(mir[w+2:2w+1])
    mir[end-w+1:end] = reverse(mir[end-2w:end-w-1])
    return mir
end

function find_extrema_count(A::Array,d::Int=3)
    count = 0
    mA = mirror(A)
    for i in 1:lastindex(mA)-2
        win = view(mA,i:i+2)

        @views if (win[2] > win[1] && win[2] > win[3]) || (win[2] < win[1] && win[2] < win[3])
                count += 1
        end
    end
    return count
end

function find_extrema(A::Array,d::Int=3)
    count = 0
    mA = mirror(A)
    exts = Int[]

    for i in 1:lastindex(mA)-d+1
        win = view(mA,i:i+d-1)
        if (win[2] > win[1] && win[2] > win[3])
                count+=1
                push!(exts,i)
        elseif (win[2] < win[1] && win[2] < win[3])
                count+=1
                push!(exts,i)
        end
    end
    count, exts
end

function midpoints(s)
    nex, exts = find_extrema(s)
    # @views extvals = s[exts]
    midpts = exts[2:end] .-div.(diff(exts),2)
    # extvals = diff(extvals)./2

    pushfirst!(midpts,1)
    push!(midpts,length(s))

    return convert(Array{Int64},midpts)
end

function stream_minmax(env1, env2, a, d::Int64)
    a = mirror(a,d)
    N = length(a)
    upper = Int[] #buffer for indices
    lower = Int[]

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
