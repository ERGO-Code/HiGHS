module logparser

const STATUS_UNKNOWN = -1
const STATUS_OPTIMAL = 0
const STATUS_IMPRECISE = 1
const STATUS_TIME_LIMIT = 2
const STATUS_ITER_LIMIT = 3
const STATUS_OUT_OF_MEM = 4
const STATUS_FAILED = 5

type simplexlog
    status::Int
    iter::Int
    time_total::Float64
    time_presolve::Float64
    function simplexlog()
        return new(STATUS_UNKNOWN, 0, 0.0, 0.0)
    end
end

type barrierlog
    status_ipm::Int
    status_crossover::Int
    iter_ipm::Int
    time_total::Float64
    time_crossover::Float64
    time_presolve::Float64
    function barrierlog()
        return new(STATUS_UNKNOWN, STATUS_UNKNOWN, 0, 0.0, 0.0, 0.0)
    end
end

"""
 GurobiSimplex(logfile) -> simplexlog

 Parses logfile according to the format of the Gurobi simplex log.
"""
function GurobiSimplex(logfile::AbstractString)
    lines = readlines(logfile)
    log = simplexlog()
    for line in lines
        s = match(r"^Presolve time\: ([\.0-9]+)s$", line)
        if s != nothing
            log.time_presolve = parse(Float64, s.captures[1])
        end
        s = match(r"^Solved in ([0-9]+) iterations and ([\.0-9]+) seconds$",
                  line)
        if s != nothing
            log.iter = parse(Int, s.captures[1])
            log.time_total = parse(Float64, s.captures[2])
            log.status = STATUS_OPTIMAL
        end
        s = match(r"^Stopped in ([0-9]+) iterations and ([\.0-9]+) seconds$",
                  line)
        if s != nothing
            log.iter = parse(Int, s.captures[1])
            log.time_total = parse(Float64, s.captures[2])
        end
        if match(r"^Time limit reached$", line) != nothing
            log.status = STATUS_TIME_LIMIT
        end
        if match(r"^Iteration limit reached$", line) != nothing
            log.status = STATUS_ITER_LIMIT
        end
    end
    return log
end
 
"""
 GurobiBarrier(logfile) -> barrierlog

 Parses logfile according to the format of the Gurobi barrier log.
"""
function GurobiBarrier(logfile::AbstractString)
    lines = readlines(logfile)
    log = barrierlog()
    l = 1                       # line number
    time_barrier = 0.0          # time when barrier stopped

    # Parse barrier log
    while l <= length(lines)
        line = lines[l]
        if match(r"^Crossover log...$", line) != nothing break; end
        s = match(r"^Presolve time\: ([\.0-9]+)s$", line)
        if s != nothing
            log.time_presolve = parse(Float64, s.captures[1])
        end
        s = match(r"^Barrier performed ([0-9]+) iterations in ([\.0-9]+) seconds$",
                  line)
        if s != nothing
            log.iter_ipm = parse(Int, s.captures[1])
            time_barrier = parse(Float64, s.captures[2])
        end
        s = match(r"^Barrier solved model in ([0-9]+) iterations and ([\.0-9]+) seconds$",
                  line)
        if s != nothing
            log.iter_ipm = parse(Int, s.captures[1])
            log.status_ipm = STATUS_OPTIMAL
            time_barrier = parse(Float64, s.captures[2])
        end
        if match(r"^Sub-optimal termination", line) != nothing
            log.status_ipm = STATUS_IMPRECISE
        end
        if match(r"^Numerical trouble encountered$", line) != nothing
            log.status_ipm = STATUS_FAILED
            l = l+1
            break               # continue with crossover/simplex log
        end
        l = l+1
    end

    # Parse crossover log
    while l <= length(lines)
        line = lines[l]
        s = match(r"^Solved in ([0-9]+) iterations and ([\.0-9]+) seconds$",
                  line)
        if s != nothing
            log.time_total = parse(Float64, s.captures[2])
            log.status_crossover = STATUS_OPTIMAL
        end
        s = match(r"^Stopped in ([0-9]+) iterations and ([\.0-9]+) seconds$",
                  line)
        if s != nothing
            log.time_total = parse(Float64, s.captures[2])
        end
        if match(r"^Time limit reached$", line) != nothing
            log.status_crossover = STATUS_TIME_LIMIT
        end
        if match(r"^Iteration limit reached$", line) != nothing
            log.status_crossover = STATUS_ITER_LIMIT
        end
        l = l+1
    end

    log.time_crossover = log.time_total - time_barrier
    return log
end

end
