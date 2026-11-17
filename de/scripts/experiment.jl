# Experimento completo: todas las conjeturas × clases × algoritmos × corridas.
#
# Uso (desde la carpeta del proyecto):
#   julia -t auto scripts/experiment.jl                         # protocolo completo
#   julia -t auto scripts/experiment.jl --runs 5 --time 60      # prueba corta
#   julia -t auto scripts/experiment.jl --ids 197,289 --algs mde
#   julia -t auto scripts/experiment.jl --only refuted --time 120
#
# Cada corrida se guarda apenas termina en results/runs.csv, así que se puede
# cortar con Ctrl+C y retomar: las corridas ya hechas se saltean.

include(joinpath(@__DIR__, "..", "src", "GraffitiDE.jl"))
using .GraffitiDE
using LinearAlgebra, Printf, Dates

BLAS.set_num_threads(1)   # el paralelismo es entre corridas, no dentro de cada autovalor

function argval(flag, default)
    i = findfirst(==(flag), ARGS)
    return i === nothing ? default : ARGS[i+1]
end

runs     = parse(Int, argval("--runs", "30"))
tlimit   = parse(Float64, argval("--time", "900"))
algs     = Symbol.(split(argval("--algs", "mdes,mde,de,rs"), ","))
only     = argval("--only", "all")              # all | refuted | open
idsarg   = argval("--ids", "")
outdir   = argval("--out", joinpath(@__DIR__, "..", "results"))

ids = isempty(idsarg) ? sort(collect(keys(CONJECTURES))) : parse.(Int, split(idsarg, ","))
only == "refuted" && (ids = filter(i -> !isnan(CONJECTURES[i].ref_time), ids))
only == "open"    && (ids = filter(i -> isnan(CONJECTURES[i].ref_time), ids))

mkpath(joinpath(outdir, "counterexamples"))
csv = joinpath(outdir, "runs.csv")
header = "conj,class,alg,run,seed,success,time_s,evals,n,E,graph6,restarts,ref_time_s,timestamp,verified"
isfile(csv) || open(io -> println(io, header), csv, "w")

done = Set{Tuple{Int,String,String,Int}}()
for line in Iterators.drop(eachline(csv), 1)
    f = split(line, ",")
    length(f) >= 4 && push!(done, (parse(Int, f[1]), f[2], f[3], parse(Int, f[4])))
end

tasks = [(id, cls, alg, r) for id in ids for cls in CONJECTURES[id].classes for alg in algs for r in 1:runs
         if !((id, String(cls), String(alg), r) in done)]

seedof(id, cls, r) = 1_000_000 * id + 1000 * findfirst(==(cls), [:any, :tree, :k3free, :girth5]) + r

@printf("%d corridas pendientes, %d hilos, %.0f s por corrida (estimado: %.1f h)\n",
        length(tasks), Threads.nthreads(), tlimit, length(tasks) * tlimit / Threads.nthreads() / 3600)

# Precompila con una corrida corta para que el tiempo de compilación no cuente.
search(166, :any, :mde; seed = 0, params = Params(time_limit = 2.0))

lk = ReentrantLock()
Threads.@threads :dynamic for (id, cls, alg, r) in tasks
    seed = seedof(id, cls, r)    # misma semilla para los tres algoritmos: comparación pareada
    res = search(id, cls, alg; seed = seed, params = Params(time_limit = tlimit))
    c = CONJECTURES[id]
    ver = res.success && verify(id, res.A)
    lock(lk) do
        open(csv, "a") do io
            @printf(io, "%d,%s,%s,%d,%d,%s,%.3f,%d,%d,%.6g,%s,%d,%s,%s,%s\n",
                    id, cls, alg, r, seed, res.success, res.time, res.evals, res.n, res.E,
                    res.g6, res.restarts, isnan(c.ref_time) ? "" : string(c.ref_time), now(), ver)
        end
        if res.success
            fn = joinpath(outdir, "counterexamples", "G$(id)_$(cls)_$(alg)_run$(r).txt")
            open(fn, "w") do io
                println(io, "Graffiti $(id) ($(c.statement)), clase $(cls), algoritmo $(alg), corrida $(r)")
                println(io, "n = $(res.n), E(G) = $(res.E), tiempo = $(round(res.time, digits = 2)) s, evaluaciones = $(res.evals)")
                println(io, "verificado al reetiquetar vértices: $(ver)")
                println(io, "graph6: $(res.g6)")
                println(io, "aristas: $(res.edges)")
            end
        end
        @printf("[%s] G%-4d %-7s %-4s run %2d  %s  %8.2f s  n=%3d  E=% .4f\n",
                Dates.format(now(), "HH:MM:SS"), id, cls, alg, r,
                res.success ? (ver ? "REFUTADA" : "REF?(no verif.)") : "--------", res.time, res.n, res.E)
        flush(stdout)
    end
end
println("Listo. Resultados en $(csv)")
