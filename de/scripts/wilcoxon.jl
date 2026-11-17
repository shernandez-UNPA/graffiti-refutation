# Comparaciones pareadas entre algoritmos (test de Wilcoxon de rangos con signo).
# Uso: julia scripts/wilcoxon.jl [results/runs.csv] [pares]
#   pares por defecto: de-rs,mde-de,mdes-mde,des-de
# Compara evaluaciones de la función objetivo, emparejando por corrida (misma
# semilla). Una corrida fallida cuenta como peor que cualquier éxito.
# Aproximación normal con corrección por empates; corrección de Holm.

using Printf

csv = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "results", "runs.csv")
pairs = [Tuple(Symbol.(split(p, "-"))) for p in split(length(ARGS) >= 2 ? ARGS[2] : "de-rs,mde-de,mdes-mde,des-de", ",")]

R = Dict{Tuple{Int,String,Symbol},Dict{Int,Float64}}()
for line in Iterators.drop(eachline(csv), 1)
    f = split(line, ",")
    length(f) < 8 && continue
    ok = f[6] == "true" && (length(f) < 15 || f[15] == "true")
    d = get!(R, (parse(Int, f[1]), String(f[2]), Symbol(f[3])), Dict{Int,Float64}())
    d[parse(Int, f[4])] = ok ? parse(Float64, f[8]) : Inf
end

normccdf2(z) = erfc_(abs(z) / sqrt(2))
# erfc por fracción continua / serie (suficiente para p-valores)
function erfc_(x)
    t = 1 / (1 + 0.5x)
    return t * exp(-x^2 - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
           t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
           t * (-0.82215223 + t * 0.17087277)))))))))
end

function wilcoxon(a::Dict, b::Dict)
    d = Float64[]
    for k in keys(a)
        haskey(b, k) || continue
        x, y = a[k], b[k]
        x == y && continue
        push!(d, isinf(x) ? 1e18 : (isinf(y) ? -1e18 : x - y))
    end
    n = length(d)
    n == 0 && return (1.0, 0.0, 0.0, 0)
    o = sortperm(abs.(d))
    r = zeros(n); tie = 0.0; i = 1
    while i <= n
        j = i
        while j < n && abs(d[o[j+1]]) == abs(d[o[i]]); j += 1; end
        for k in i:j; r[o[k]] = (i + j) / 2; end
        t = j - i + 1; tie += t^3 - t; i = j + 1
    end
    Wp = sum(r[d .> 0]; init = 0.0); Wm = sum(r[d .< 0]; init = 0.0)
    μ = n * (n + 1) / 4; σ2 = n * (n + 1) * (2n + 1) / 24 - tie / 48
    σ2 <= 0 && return (1.0, Wp, Wm, n)
    return (normccdf2((min(Wp, Wm) - μ) / sqrt(σ2)), Wp, Wm, n)
end

tests = []
for (id, cls) in sort(unique([(k[1], k[2]) for k in keys(R)])), (a, b) in pairs
    haskey(R, (id, cls, a)) && haskey(R, (id, cls, b)) || continue
    p, Wp, Wm, n = wilcoxon(R[(id, cls, a)], R[(id, cls, b)])
    push!(tests, (p = p, id = id, cls = cls, a = a, b = b, win = Wm > Wp ? a : (Wp > Wm ? b : :empate)))
end
sort!(tests; by = t -> t.p)
m = length(tests)
println("Comparaciones: $m (Holm, α = 0.05)")
@printf("%-6s %-8s %-12s %-10s %-9s %s\n", "conj", "clase", "comparación", "p", "Holm", "menos evaluaciones")
stop = false
for (k, t) in enumerate(tests)
    global stop
    sig = !stop && t.p <= 0.05 / (m - k + 1)
    sig || (stop = true)
    @printf("%-6d %-8s %-12s %-10.2e %-9s %s\n", t.id, t.cls, "$(t.a) vs $(t.b)", t.p, sig ? "sí" : "no", t.win)
end
