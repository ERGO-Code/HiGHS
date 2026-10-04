using HPRLP
using Printf

length(ARGS) == 5 || error("usage: neso_hpr_julia.jl MODEL PREFIX CR TOL TIME_LIMIT")
println("NESO Julia input: ", ARGS[1])
flush(stdout)
model = HPRLP.build_from_mps(ARGS[1])
flush(stdout)
params = HPRLP.HPRLP_parameters()
params.use_gpu = true
params.device_number = 0
params.presolve = "NONE"
params.warm_up = true
params.use_Curtis_Reid_scaling = ARGS[3] == "true"
params.stoptol = parse(Float64, ARGS[4])
params.time_limit = parse(Float64, ARGS[5])
result = HPRLP.optimize(model, params)
flush(stdout)
for (suffix, values) in ((".primal.txt", result.x), (".dual.txt", result.y), (".reduced_cost.txt", result.z))
    open(ARGS[2] * suffix, "w") do io
        for v in Array(values)
            @printf(io, "%.17g\n", v)
        end
    end
end
open(ARGS[2] * ".json", "w") do io
    @printf(io, "{\"status\":\"%s\",\"iterations\":%d,\"time\":%.17g,\"objective\":%.17g,\"residual\":%.17g,\"gap\":%.17g}\n",
            result.status, result.iter, result.time, result.primal_obj, result.residuals, result.gap)
end
