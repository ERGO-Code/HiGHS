include("logparser.jl")
include("../../Julia/ipx.jl")

module matwriter

using MAT
using logparser
using ipx

"""
 IPXRuntimes(testset::AbstractString, logdir::AbstractString,
             outfile::AbstractString)

 writes Matlab binary file to @outfile containing arrays time_total,
 time_starting_basis, time_kkt_factorize, time_kkt_solve and time_crossover. The
 runtimes are extracted from the IPX logfiles. If a model was not solved
 (including crossover), then all 5 arrays have entry NaN.


 @testset: file from testsets/ specifying the set of LP models

 @logdir: directory holding subdirectory ipx/ that contains the logfiles
"""
function IPXRuntimes(testset::AbstractString, logdir::AbstractString,
                     outfile::AbstractString)
    models, = readdlm(testset, header=true)
    ids = convert(Array{Int}, models[:,1])
    n = length(ids)
    time_total = Array{Float64}(n)
    time_starting_basis = Array{Float64}(n)
    time_kkt_factorize = Array{Float64}(n)
    time_kkt_solve = Array{Float64}(n)
    time_crossover = Array{Float64}(n)
    for i = 1:n
        logfile = @sprintf("%04d.log", ids[i])
        ipxinfo = ipx.ParseInfo(joinpath(logdir, "ipx", logfile))
        if ipxinfo.status_crossover == 1 || ipxinfo.status_crossover == 2
            time_total[i] = ipxinfo.time_total
            time_starting_basis[i] = ipxinfo.time_starting_basis
            time_kkt_factorize[i] = ipxinfo.time_kkt_factorize
            time_kkt_solve[i] = ipxinfo.time_kkt_solve
            time_crossover[i] = ipxinfo.time_crossover
        else
            time_total[i] = NaN
            time_starting_basis[i] = NaN
            time_kkt_factorize[i] = NaN
            time_kkt_solve[i] = NaN
            time_crossover[i] = NaN
        end
    end
    fout = matopen(outfile, "w")
    write(fout, "time_total", time_total)
    write(fout, "time_starting_basis", time_starting_basis)
    write(fout, "time_kkt_factorize", time_kkt_factorize)
    write(fout, "time_kkt_solve", time_kkt_solve)
    write(fout, "time_crossover", time_crossover)
    close(fout)
end

end
