module EMD

using Dierckx, ElasticArrays
include("utils.jl")
include("bEMD.jl")
include("sEMD.jl")
export sEMD, bEMD

end
