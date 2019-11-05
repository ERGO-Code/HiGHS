include("logparser.jl")
include("../../Julia/ipx.jl")

module tblwriter

using logparser
using ipx

# Give each group a label that takes less width in the table.
const grouplabel = Dict("castro/cta" => "1a",
                        "castro/irp" => "1b",
                        "castro/L1" => "1c",
                        "castro/Linf" => "1d",
                        "hall" => "2",
                        "kennington" => "3",
                        "meszaros/misc" => "4a",
                        "meszaros/New" => "4b",
                        "meszaros/problematic" => "4c",
                        "meszaros/stochlp" => "4d",
                        "miplib2010" => "5",
                        "mittelmann/fome" => "6a",
                        "mittelmann/misc" => "6b",
                        "mittelmann/nug" => "6c",
                        "mittelmann/pds" => "6d",
                        "mittelmann/rail" => "6e",
                        "netlib" => "7")

"""
 Results(testset::AbstractString, logdir::AbstractString,
         outfile::AbstractString, simplex::Bool=true, latex::Bool=false)

 Writes table with problem dimensions and solution times to @outfile.

 @testset: file from testsets/ specifying the set of LP models

 @logdir: directory holding subdirectories ipx/, ipx_cleanup/, gurobi_barrier/
 and gurobi_simplex/ that contain the logfiles

 @simplex: if true, then dual simplex is included in the output

 @latex: if true, then output is a Latex table
"""
function Results(testset::AbstractString, logdir::AbstractString,
                 outfile::AbstractString, simplex::Bool=true, latex::Bool=false)
    models, = readdlm(testset, header=true)
    ids = convert(Array{Int}, models[:,1])
    names = convert(Array{String}, models[:,2])
    groups = convert(Array{String}, models[:,3])
    rows = convert(Array{Int}, models[:,4])
    cols = convert(Array{Int}, models[:,5])
    entries = convert(Array{Int}, models[:,6])
    fout = open(outfile, "w")
    colsep = "   "
    lineend = "\n"
    if latex
        colsep = " & "
        lineend = "\\\\\n"
    end
    # --------------------------------------------------------------------------
    # Header
    # --------------------------------------------------------------------------
    if latex
        @printf(fout, "\\begin{longtable}")
        @printf(fout, "{llrrr|rr@{\\hspace{1pt}}lrr|rr@{\\hspace{1pt}}lr")
        if simplex @printf(fout, "|r"); end
        @printf(fout, "}\n")
        @printf(fout, "\\hline\n")
        # First header line
        @printf(fout, " ")
        @printf(fout, "%3s", "src")
        @printf(fout, "%s", colsep)
        @printf(fout, "%-25s", "name")
        @printf(fout, "%s", colsep)
        @printf(fout, "%9s", "\$\\bar{m}\$")
        @printf(fout, "%s", colsep)
        @printf(fout, "%9s", "\$\\bar{n}\$")
        @printf(fout, "%s", colsep)
        @printf(fout, "%9s", "\$\\nnz\$")
        @printf(fout, "%s", colsep)
        @printf(fout, "%s", "\\multicolumn{5}{|c|}{\\ipx}")
        @printf(fout, "%s", colsep)
        @printf(fout, "%s", "\\multicolumn{4}{|c}{\\gurobi\\ barrier}")
        if simplex
            @printf(fout, "%s", colsep)
            @printf(fout, "%s", "simplex")
        end
        @printf(fout, "%s", lineend)
        # Second header line
        @printf(fout, " ")
        @printf(fout, "%3s", "")
        @printf(fout, "%s", colsep)
        @printf(fout, "%-25s", "")
        @printf(fout, "%s", colsep)
        @printf(fout, "%9s", "")
        @printf(fout, "%s", colsep)
        @printf(fout, "%9s", "")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "total")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "IPM")
        @printf(fout, "%s", colsep)
        @printf(fout, "%s", " ") # flag for imprecise solution
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "push")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "cleanup")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "total")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "IPM")
        @printf(fout, "%s", colsep)
        @printf(fout, "%s", " ") # flag for imprecise solution
        @printf(fout, "%s", colsep)
        @printf(fout, "%9s", "crossover")
        @printf(fout, "%s", "\\\\\n")
        @printf(fout, "%s", "\\hline \\endhead\n")
        @printf(fout, "%s", "\\hline \\endfoot\n")
        @printf(fout, "%s", "\\hline \\\\")
        @printf(fout, "%s", "\\multicolumn{10}{l}{f: failed, t: time limit, r: IPM solution reported not optimal}\\endlastfoot\n")
    else
        @printf(fout, " ")
        @printf(fout, "%3s", "src")
        @printf(fout, "%s", colsep)
        @printf(fout, "%-25s", "name")
        @printf(fout, "%s", colsep)
        @printf(fout, "%9s", "rows")
        @printf(fout, "%s", colsep)
        @printf(fout, "%9s", "cols")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "nnz")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "IPXtotal")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "IPX_IPM")
        @printf(fout, "%s", colsep)
        @printf(fout, "%s", " ") # flag for imprecise solution
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "IPXpush")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "IPXclean")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "GRBtotal")
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "GRB_IPM")
        @printf(fout, "%s", colsep)
        @printf(fout, "%s", " ") # flag for imprecise solution
        @printf(fout, "%s", colsep)
        @printf(fout, "%8s", "GRB_cros")
        if simplex
            @printf(fout, "%s", colsep)
            @printf(fout, "%8s", "dsimplex")
        end
        @printf(fout, "%s", lineend)
    end
    # --------------------------------------------------------------------------
    # Body
    # --------------------------------------------------------------------------
    for i = 1:length(ids)
        logfile = @sprintf("%04d.log", ids[i])
        simplexlog = logparser.GurobiSimplex(
            joinpath(logdir, "gurobi_dsimplex", logfile))
        barrierlog = logparser.GurobiBarrier(
            joinpath(logdir, "gurobi_barrier", logfile))
        ipxinfo = ipx.ParseInfo(joinpath(logdir, "ipx", logfile))
        # Logfile from simplex clean-up exists only when IPX crossover
        # terminated successfully.
        cleanuplog = logparser.simplexlog()
        if ipxinfo.status_crossover == 1 || ipxinfo.status_crossover == 2
            cleanuplog = logparser.GurobiSimplex(
                joinpath(logdir, "ipx_cleanup", logfile))
        end
        @printf(fout, " ")
        @printf(fout, "%3s", grouplabel[groups[i]])
        @printf(fout, "%s", colsep)
        name = names[i]
        if latex
            name = replace(name, "_" => "\\_")
            name = string("\\ct{", name, "}")
        end
        @printf(fout, "%-25s", name)
        @printf(fout, "%s", colsep)
        @printf(fout, "%9d", rows[i])
        @printf(fout, "%s", colsep)
        @printf(fout, "%9d", cols[i])
        @printf(fout, "%s", colsep)
        @printf(fout, "%8d", entries[i])
        @printf(fout, "%s", colsep)
        #-----------------------------------------------------------------------
        # IPX total
        #-----------------------------------------------------------------------
        if cleanuplog.status == logparser.STATUS_OPTIMAL
            # Model was solved to basic solution.
            time_total = ipxinfo.time_total
            if cleanuplog.iter > 0
                # clean-up was needed
                time_total += cleanuplog.time_total
            end
            @printf(fout, "%8.1f", time_total)
        else
            # Failed for whatever reason. Leave "total" field empty.
            @printf(fout, "%8s", "")
        end
        @printf(fout, "%s", colsep)
        #-----------------------------------------------------------------------
        # IPX: time IPM and flag if imprecise
        #-----------------------------------------------------------------------
        if ipxinfo.status_ipm == 1 || ipxinfo.status_ipm == 2
            @printf(fout, "%8.1f", ipxinfo.time_total-ipxinfo.time_crossover)
        elseif ipxinfo.status_ipm == 5
            @printf(fout, "%8s", "t") # time limit
        else
            @printf(fout, "%8s", "f") # failed for other reason
        end
        @printf(fout, "%s", colsep)
        if ipxinfo.status_ipm == 2
            @printf(fout, "%s", "r") # IPM solution imprecise
        else
            @printf(fout, "%s", " ")
        end
        @printf(fout, "%s", colsep)
        #-----------------------------------------------------------------------
        # IPX push phases
        #-----------------------------------------------------------------------
        if ipxinfo.status_crossover == 0
            @printf(fout, "%8s", "") # crossover not run
        elseif ipxinfo.status_crossover == 1 || ipxinfo.status_crossover == 2
            @printf(fout, "%8.1f", ipxinfo.time_crossover)
        elseif ipxinfo.status_crossover == 5
            @printf(fout, "%8s", "t") # time limit
        else
            @printf(fout, "%8s", "f") # failed for other reason
        end
        @printf(fout, "%s", colsep)
        #-----------------------------------------------------------------------
        # IPX clean-up
        #-----------------------------------------------------------------------
        if cleanuplog.iter > 0
            if cleanuplog.status == logparser.STATUS_OPTIMAL
                @printf(fout, "%8.1f", cleanuplog.time_total)
            elseif cleanuplog.status == logparser.STATUS_TIME_LIMIT
                @printf(fout, "%8s", "t") # time limit
            else
                @printf(fout, "%8s", "f") # failed for other reason
            end
        else
            @printf(fout, "%8s", "") # no simplex clean-up performed
        end
        @printf(fout, "%s", colsep)
        #-----------------------------------------------------------------------
        # Gurobi barrier total
        #-----------------------------------------------------------------------
        if barrierlog.status_crossover == logparser.STATUS_OPTIMAL
            @printf(fout, "%8.1f",
                    barrierlog.time_total-barrierlog.time_presolve)
        else
            # Model not solved for whatever reason. Leave "total" field empty.
            @printf(fout, "%8s", "")
        end
        @printf(fout, "%s", colsep)
        #-----------------------------------------------------------------------
        # Gurobi IPM: time and flag if sub-optimal
        #-----------------------------------------------------------------------
        if barrierlog.status_crossover == logparser.STATUS_OPTIMAL ||
            barrierlog.status_ipm == logparser.STATUS_OPTIMAL ||
            barrierlog.status_ipm == logparser.STATUS_IMPRECISE
            # Model was solved, so IPM either optimal, sub-optimal or numerical
            # trouble. Report IPM time in any case.
            time_ipm = barrierlog.time_total - barrierlog.time_crossover -
                barrierlog.time_presolve
            @printf(fout, "%8.1f", time_ipm)
        elseif barrierlog.status_ipm == logparser.STATUS_TIME_LIMIT
            @printf(fout, "%8s", "t")
        else                    # failed for other reason
            @printf(fout, "%8s", "f")
        end
        @printf(fout, "%s", colsep)
        if barrierlog.status_crossover == logparser.STATUS_OPTIMAL &&
            barrierlog.status_ipm != logparser.STATUS_OPTIMAL
            @printf(fout, "%s", "r") # IPM solution not optimal
        else
            @printf(fout, "%s", " ")
        end
        @printf(fout, "%s", colsep)
        #-----------------------------------------------------------------------
        # Gurobi crossover
        #-----------------------------------------------------------------------
        if barrierlog.status_crossover == logparser.STATUS_UNKNOWN
            @printf(fout, "%8s", "") # not run
        elseif barrierlog.status_crossover == logparser.STATUS_OPTIMAL
                @printf(fout, "%8.1f", barrierlog.time_crossover)
        elseif barrierlog.status_crossover == logparser.STATUS_TIME_LIMIT
            @printf(fout, "%8s", "t")
        else                    # failed for other reason
            @printf(fout, "%8s", "f")
        end
        #-----------------------------------------------------------------------
        # Gurobi simplex
        #-----------------------------------------------------------------------
        if simplex
            @printf(fout, "%s", colsep)
            if simplexlog.status == logparser.STATUS_OPTIMAL
                @printf(fout, "%8.1f",
                        simplexlog.time_total-simplexlog.time_presolve)
            elseif simplexlog.status == logparser.STATUS_TIME_LIMIT
                @printf(fout, "%8s", "t") # time limit
            else
                @printf(fout, "%8s", "f") # failed for other reason
            end
        end
        @printf(fout, "%s", lineend)
    end
    if latex
        @printf(fout, "\\end{longtable}\n")
    end
    close(fout)
end

"""
 Comparison(testset::AbstractString, logdir::AbstractString,
            outfile::AbstractString)

 Writes table comparing the geometric mean of runtime ratios to @outfile.

 @testset: file from testsets/ specifying the set of LP models

 @logdir: directory holding subdirectories ipx/, ipx_cleanup/ and
 gurobi_barrier/ that contain the logfiles
"""
function Comparison(testset::AbstractString, logdir::AbstractString,
                    outfile::AbstractString)
    models, = readdlm(testset, header=true)
    ids = convert(Array{Int}, models[:,1])
    # Build arrays of runtimes of IPX and Gurobi. Put NaN if not solved.
    n = length(ids)
    time_ipx = Array{Float64}(n)
    time_grb = Array{Float64}(n)
    for i = 1:n
        logfile = @sprintf("%04d.log", ids[i])
        barrierlog = logparser.GurobiBarrier(
            joinpath(logdir, "gurobi_barrier", logfile))
        ipxinfo = ipx.ParseInfo(joinpath(logdir, "ipx", logfile))
        # Logfile from simplex clean-up exists only when IPX crossover
        # terminated successfully.
        cleanuplog = logparser.simplexlog()
        if ipxinfo.status_crossover == 1 || ipxinfo.status_crossover == 2
            cleanuplog = logparser.GurobiSimplex(
                joinpath(logdir, "ipx_cleanup", logfile))
        end
        if cleanuplog.status == logparser.STATUS_OPTIMAL
            time_ipx[i] = ipxinfo.time_total + cleanuplog.time_total
        else
            time_ipx[i] = NaN
        end
        if barrierlog.status_crossover == logparser.STATUS_OPTIMAL
            time_grb[i] = barrierlog.time_total - barrierlog.time_presolve
        else
            time_grb[i] = NaN
        end
    end
    # Print table.
    fout = open(outfile, "w")
    @printf(fout, " %6s %9s %10s %10s %13s\n", "subset", "instances",
                    "IPX/Gurobi", "IPX faster", "Gurobi faster")
    print_comparison_line(fout, time_ipx, time_grb, 1)
    print_comparison_line(fout, time_ipx, time_grb, 10)
    print_comparison_line(fout, time_ipx, time_grb, 100)
    close(fout)
end

function print_comparison_line(fout, time_ipx, time_grb, subset)
    i = find(isfinite.(time_ipx) .& isfinite.(time_grb) .&
             (max(time_ipx,time_grb).>subset))
    n = length(i)
    ratio = prod((time_ipx[i]./time_grb[i]).^(1/n))
    ipx_faster = sum(time_ipx[i] .< time_grb[i])
    grb_faster = sum(time_ipx[i] .> time_grb[i])
    @printf(fout, " %6d %9d %10.2f %10d %13d\n",
            subset, n, ratio, ipx_faster, grb_faster)
end

"""
 BasisUpdates(testset::AbstractString, logdir::AbstractString,
              outfile::AbstractString)

 Writes table with number of basis updates performed by IPX (incl. crossover and
 simplex clean-up) and by Gurobi dual simplex to @outfile. The number of basis
 updates is divided by min(m,n), where m-by-n is the dimension of the presolved
 constraint matrix.

 @testset: file from testsets/ specifying the set of LP models

 @logdir: directory holding subdirectories ipx/, ipx_cleanup/ and
 gurobi_dsimplex/ that contain the logfiles
"""
function BasisUpdates(testset::AbstractString, logdir::AbstractString,
                      outfile::AbstractString)
    models, = readdlm(testset, header=true)
    ids = convert(Array{Int}, models[:,1])
    names = convert(Array{String}, models[:,2])
    rows = convert(Array{Int}, models[:,4])
    cols = convert(Array{Int}, models[:,5])
    fout = open(outfile, "w")
    @printf(fout, "%6s %8s %8s\n", "name", "simplex", "IPX")
    for i = 1:length(ids)
        logfile = @sprintf("%04d.log", ids[i])
        simplexlog = logparser.GurobiSimplex(
            joinpath(logdir, "gurobi_dsimplex", logfile))
        ipxinfo = ipx.ParseInfo(joinpath(logdir, "ipx", logfile))
        # Logfile from simplex clean-up exists only when IPX crossover
        # terminated successfully.
        cleanuplog = logparser.simplexlog()
        if ipxinfo.status_crossover == 1 || ipxinfo.status_crossover == 2
            cleanuplog = logparser.GurobiSimplex(
                joinpath(logdir, "ipx_cleanup", logfile))
        end
        updates_simplex = simplexlog.iter
        updates_ipx = ipxinfo.updates_ipm + ipxinfo.updates_crossover +
            cleanuplog.iter
        dim = min(rows[i], cols[i])
        @printf(fout, "%6s %8.2f %8.2f\n", names[i], updates_simplex/dim,
                updates_ipx/dim)
    end
    close(fout)
end

end
