#=* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                       *
*    This file is part of the HiGHS linear optimization suite           *
*                                                                       *
*    Written and engineered 2008-2023 by Julian Hall, Ivet Galabova,    *
*    Leona Gottwald and Michael Feldmeier                               *
*    Available as open-source under the MIT License                     *
*                                                                       *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *=#

import Clang: Generators

highs_src = joinpath(dirname(dirname(@__DIR__)), "src")
c_api = joinpath(highs_src, "interfaces", "highs_c_api.h")
libhighs_filename = joinpath(@__DIR__, "libhighs.jl")

Generators.build!(
    Generators.create_context(
        [c_api, joinpath(highs_src, "util", "HighsInt.h")],
        [Generators.get_default_args(); "-I$highs_src"; "-I$(@__DIR__)"],
        Dict{String,Any}(
            "general" => Dict{String,Any}(
                "output_file_path" => libhighs_filename,
                "library_name" => "libhighs",
                "print_using_CEnum" => false,
                "extract_c_comment_style" => "doxygen",
            )
        ),
    ),
)

write(
    libhighs_filename,
    replace(
        read(libhighs_filename, String),
        "[`HighsInt`](@ref)" => "`HighsInt`",
    ),
)

open(libhighs_filename, "a") do io
    for line in readlines(c_api)
        m = match(r"const HighsInt kHighs([a-zA-Z]+) = (-?[0-9]+);", line)
        if m === nothing
            continue
        end
        println(io, "const kHighs$(m[1]) = HighsInt($(m[2]))")
    end
end
