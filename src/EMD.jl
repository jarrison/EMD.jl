module EMD

using Dierckx, ElasticArrays
include("utils.jl")
include("bEMD.jl")
include("sfemd.jl")
export sEMD, bEMD

end
