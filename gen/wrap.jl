using Clang

rustpath = "/home/jarrison/Dropbox/Code/Rust Projects/rust-spline"

headers = [joinpath(rustpath,"include","myheader.h")]

wc = init(; headers=headers,
                output_file = joinpath(@__DIR__,"librust_spline_api.jl"),
                common_file = joinpath(@__DIR__,"librust_spline_common.jl"),
                # clang_args = [],
                clang_includes = ["usr/include/c++/v1",CLANG_INCLUDE],
                header_wrapped=(root,current)->root == current,
                header_library=x->":librust_spline",
                clang_diagnostics=true)

run(wc)
