# Resume results/runs.csv por conjetura, clase y algoritmo.
# Uso: julia scripts/summary.jl [results/runs.csv]
# Genera results/summary.csv y results/summary.tex (tabla para el artículo).

using Statistics, Printf

csv = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "results", "runs.csv")
outdir = dirname(csv)

rows = [split(l, ",") for l in Iterators.drop(eachline(csv), 1) if !isempty(strip(l))]
# columnas: conj,class,alg,run,seed,success,time_s,evals,n,E,graph6,restarts,ref_time_s,timestamp
groups = Dict{Tuple{Int,String,String},Vector{Vector{SubString{String}}}}()
for f in rows
    push!(get!(groups, (parse(Int, f[1]), f[2], f[3]), []), f)
end

algorder = Dict("mdes" => 1, "mde" => 2, "de" => 3, "rs" => 4)
ks = sort(collect(keys(groups)); by = k -> (k[1], k[2], get(algorder, k[3], 9)))

med(v) = isempty(v) ? NaN : median(v)
fmt(x) = isnan(x) ? "--" : (x < 10 ? @sprintf("%.2f", x) : @sprintf("%.0f", x))

open(joinpath(outdir, "summary.csv"), "w") do io
    println(io, "conj,class,alg,runs,successes,success_rate,median_time_s,median_evals,min_n,ref_time_s")
    for k in ks
        g = groups[k]
        ok = [f for f in g if f[6] == "true" && (length(f) < 15 || f[15] == "true")]
        t = [parse(Float64, f[7]) for f in ok]
        ev = [parse(Float64, f[8]) for f in ok]
        nmin = isempty(ok) ? "" : string(minimum(parse(Int, f[9]) for f in ok))
        @printf(io, "%d,%s,%s,%d,%d,%.3f,%s,%s,%s,%s\n", k[1], k[2], k[3], length(g), length(ok),
                length(ok) / length(g), fmt(med(t)), fmt(med(ev)), nmin, g[1][13])
    end
end

open(joinpath(outdir, "summary.tex"), "w") do io
    println(io, "% Éxitos / corridas, mediana del tiempo (s) y orden mínimo del contraejemplo")
    println(io, "\\begin{tabular}{rl l rrr r}")
    println(io, "\\hline")
    println(io, raw"Graffiti & Clase & Alg. & Éxitos & $\tilde t$ (s) & $n_{\min}$ & Ref. (s)\\")
    println(io, "\\hline")
    for k in ks
        g = groups[k]
        ok = [f for f in g if f[6] == "true" && (length(f) < 15 || f[15] == "true")]
        t = [parse(Float64, f[7]) for f in ok]
        nmin = isempty(ok) ? "--" : string(minimum(parse(Int, f[9]) for f in ok))
        ref = isempty(g[1][13]) ? "--" : g[1][13]
        @printf(io, "%d & %s & %s & %d/%d & %s & %s & %s\\\\\n", k[1], k[2], k[3], length(ok), length(g),
                fmt(med(t)), nmin, ref)
    end
    println(io, "\\hline")
    println(io, "\\end{tabular}")
end
println("Escrito: summary.csv y summary.tex en $(outdir)")
