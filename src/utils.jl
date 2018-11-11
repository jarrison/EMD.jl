
function mirror!(A::Array,d::Int=3)
    w = div(d-1,2) #width of window, d is an ODD integer
    prepend!(A,zeros(w))
    append!(A,zeros(w))
    A[1:w] = reverse(A[w+2:2w+1])
    A[end-w+1:end] = reverse(A[end-2w:end-w-1])
    return nothing
end

function mirror(A::Array,d::Int=3)
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
        if (win[2] > win[1] && win[2] > win[3]) || (win[2] < win[1] && win[2] < win[3])
                count+=1
        end
    end
    count
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

function stream_minmax(a,d::Int64)
    a = mirror(a,d)
    N = length(a)
    env1 = zeros(N-d+1)
    env2 = zeros(N-d+1)
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
        # env1[N-w+1] = mir[upper[1]]
        # env2[N-w+1] = mir[lower[1]]
    end
    return env1, env2
end

function moving_average(A,d::Int=3)
    mirror!(A,d)
    T = typeof(one(Float64)/1)
    ret = Vector{T}(undef, length(A) - d + 1)
    invd = inv(d)
    cs = cumsum(A)
    ret[1] = cs[d] * invd
    for n = 1:length(ret)-1
        ret[n+1] = (cs[n+d] - cs[n]) * invd
    end
    return ret
end
