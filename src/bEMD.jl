using Libdl
push!(DL_LOAD_PATH,"/home/jarrison/Dropbox/Code/Rust Projects/rust-spline/target/release")

struct rBSpline
    handle::Ptr{Cvoid}
end

function rBSpline(x,y,itptype=3)
    ccall((:new_spline,:librustspline),rBSpline,
                                        (Ptr{Cfloat},Ptr{Cfloat},Csize_t,Cuchar),
                                        x,y,length(x),itptype)
end

sample(spl::rBSpline,t) = ccall((:sample_spline,:librustspline),Cfloat,
                                            (Ptr{rBSpline},Cfloat),spl.handle,t)



function midpoints(s,t)
    maxs, mins, exts = find_extrema(s)
    mpx = Float64[]
    mpy = Float64[]
    for i=1:(length(exts)-1)
        push!(mpx, (t[exts[i]] + t[exts[i+1]]) /2.0)
        push!(mpy, (s[exts[i]] + s[exts[i+1]]) /2.0)
    end

    mpxx = [-mpx[3];-mpx[2];0.0 ;mpx; t[end];2.0t[end]-mpx[end-3]; 2.0t[end]-mpx[end-2]]
    mpyy = [mpy[3];mpy[2];s[1];mpy;s[end];mpy[end-3];mpy[end-2]]
    mpxx,mpyy
end

function get_spline(s,t)
    mpt, mps = midpoints(s,t)
    mpt=mpt .|> Float32
    mps=mps .|> Float32
    bspl = rBSpline(mpt,mps,2)

    return [sample(bspl,Float32(tt)) for tt in t]
end

function bsift(s,t;sdthresh=0.2)
    count = 1
    sd = 1.0
    hlast = copy(s)
    hnew = copy(s)
    while sd > sdthresh
        spls = get_spline(hlast,t)
        @. hnew = hlast - spls
        sd = sum(abs2.(hlast-hnew))/sum(hlast.^2)
        @. hlast = hnew
        count +=1
    end
    @show count
    return hnew
end

function bEMD(s,t;maximfs=5,sdthresh=0.2)
    sig = copy(s)
    N=length(s)
    imfs = ElasticArray{Float64}(undef,N,0)
    siggs = ElasticArray{Float64}(undef,N,0)
    for i in 1:maximfs
        h = bsift(sig,t)
        append!(imfs, h)
        append!(siggs,sig.-h)
    end
    return convert(Array,imfs)
end
