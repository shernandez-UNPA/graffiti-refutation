# Búsqueda de contraejemplos. Algoritmos:
#   :rs    búsqueda aleatoria (línea base)
#   :de    evolución diferencial DE/rand/1/bin
#   :mde   DE memética: DE + buscador local sobre el mejor individuo
#   :des   DE con función objetivo suavizada, sin buscador local
#   :mdes  DE memética con función objetivo suavizada: selección por
#          (E + E_suave)/2, donde E_suave reemplaza los conteos de autovalores
#          por conteos con sigmoide. El éxito se decide siempre con E exacta. Todas usan el mismo currículum por orden:
# se empieza en nmin, y la población final de orden n, extendida con claves
# nuevas, es la población inicial de orden n+1. Al llegar a nmax sin éxito
# se reinicia desde nmin con una población nueva, hasta agotar el tiempo.

using Random

Base.@kwdef struct Params
    NP::Int = 40              # tamaño de población
    F::Float64 = 0.5          # factor de escala
    CR::Float64 = 0.9         # probabilidad de cruza
    gmax::Int = 150           # generaciones máximas por orden
    stall::Int = 40           # generaciones sin mejora antes de pasar a n+1
    ls_start::Int = 30        # generación (dentro de cada orden) desde la que actúa el buscador local
    ls_tries::Int = 30        # intentos del buscador local por generación
    nmin::Int = 4             # orden inicial
    time_limit::Float64 = 900.0   # segundos por corrida (15 min, como Roucairol y Cazenave)
    max_evals::Int = typemax(Int)
end

mutable struct Result
    conj::Int
    cls::Symbol
    alg::Symbol
    seed::Int
    success::Bool
    time::Float64       # segundos hasta el contraejemplo (o hasta agotar el tiempo)
    evals::Int          # evaluaciones de la función objetivo
    n::Int              # orden del contraejemplo (o del mejor grafo hallado)
    E::Float64          # violación del contraejemplo (o la mejor hallada)
    g6::String
    edges::String
    restarts::Int
    A::Matrix{Float64}
end

"""
Verificación de un contraejemplo: la violación debe superar TOL también al
reetiquetar los vértices al azar (cambia el redondeo de los autovalores).
"""
function verify(id::Int, A; trials = 5, rng = Xoshiro(12345))
    c = CONJECTURES[id]
    n = size(A, 1)
    violation(c, A) > TOL || return false
    for _ in 1:trials
        p = randperm(rng, n)
        violation(c, A[p, p]) > TOL || return false
    end
    return true
end

struct Done <: Exception end

mutable struct State
    c::Conjecture
    cls::Symbol
    dec::Decoder
    t0::Float64
    p::Params
    evals::Int
    bestE::Float64
    bestA::Matrix{Float64}
    found::Bool
    tfound::Float64
    smooth::Bool
end

"Conjeturas en las que un lado es un conteo de autovalores (la suavización cambia algo)."
const SMOOTHABLE = Set([166, 189, 29, 30, 21, 39, 40, 20, 712, 295])

function evaluate!(S::State, x, n)
    (time() - S.t0 > S.p.time_limit || S.evals >= S.p.max_evals) && throw(Done())
    A = decode(S.dec, x, n, S.cls)
    E = violation(S.c, A)
    S.evals += 1
    if E > S.bestE
        S.bestE = E
        S.bestA = A
    end
    if E > TOL
        S.found = true
        S.tfound = time() - S.t0
        throw(Done())
    end
    S.smooth || return E
    Es = violation(S.c, A, true)
    return isfinite(E) ? (E + Es) / 2 : E
end

reflect(v) = v < 0 ? min(-v, 1.0) : (v > 1 ? max(2 - v, 0.0) : v)

"Agrega n claves nuevas (los pares del vértice n+1) a cada individuo."
function extend!(pop::Vector{Vector{Float64}}, n, rng)
    for x in pop
        append!(x, rand(rng, n))
    end
end

"Buscador local: cambia de lado del umbral claves elegidas al azar en el mejor individuo."
function local_search!(S::State, x, Ex, n, rng)
    P = npairs(n)
    y = copy(x)
    for _ in 1:S.p.ls_tries
        copyto!(y, x)
        k = 1 + rand(rng, 1:P)
        t = y[1]
        y[k] = y[k] > t ? rand(rng) * t : t + rand(rng) * (1 - t)
        if rand(rng) < 0.1
            y[1] = clamp(y[1] + 0.1 * (rand(rng) - 0.5), 0.0, 1.0)
        end
        Ey = evaluate!(S, y, n)
        if Ey >= Ex
            copyto!(x, y)
            Ex = Ey
        end
    end
    return Ex
end

function run_order!(S::State, pop, n, alg, rng)
    p = S.p
    NP = length(pop)
    D = genome_length(n)
    fit = [evaluate!(S, x, n) for x in pop]
    best = maximum(fit)
    stall = 0
    trial = zeros(D)
    for g in 1:p.gmax
        for i in 1:NP
            if alg === :rs
                rand!(rng, trial)
            else
                r1, r2, r3 = i, i, i
                while r1 == i; r1 = rand(rng, 1:NP); end
                while r2 == i || r2 == r1; r2 = rand(rng, 1:NP); end
                while r3 == i || r3 == r1 || r3 == r2; r3 = rand(rng, 1:NP); end
                jr = rand(rng, 1:D)
                xi, a, b, c = pop[i], pop[r1], pop[r2], pop[r3]
                @inbounds for j in 1:D
                    trial[j] = (rand(rng) < p.CR || j == jr) ? reflect(a[j] + p.F * (b[j] - c[j])) : xi[j]
                end
            end
            Et = evaluate!(S, trial, n)
            if Et >= fit[i]
                copyto!(pop[i], trial)
                fit[i] = Et
            end
        end
        if (alg === :mde || alg === :mdes) && g >= p.ls_start
            ib = argmax(fit)
            fit[ib] = local_search!(S, pop[ib], fit[ib], n, rng)
        end
        newbest = maximum(fit)
        if newbest > best + 1e-12
            best = newbest; stall = 0
        else
            stall += 1
            stall >= p.stall && break
        end
    end
    return pop
end

"""
    search(id, cls, alg; seed, params) -> Result

Busca un contraejemplo de la conjetura `id` en la clase `cls` con el
algoritmo `alg` (:rs, :de, :des, :mde o :mdes).
"""
function search(id::Int, cls::Symbol, alg::Symbol; seed::Int = 1, params::Params = Params())
    c = CONJECTURES[id]
    rng = Xoshiro(seed)
    S = State(c, cls, Decoder(c.nmax + 1), time(), params, 0, -Inf, zeros(1, 1), false, NaN,
              (alg === :mdes || alg === :des) && id in SMOOTHABLE)
    restarts = -1
    try
        while true
            restarts += 1
            n = params.nmin
            pop = [rand(rng, genome_length(n)) for _ in 1:params.NP]
            while n <= c.nmax
                run_order!(S, pop, n, alg, rng)
                extend!(pop, n, rng)
                n += 1
            end
        end
    catch e
        e isa Done || rethrow()
    end
    A = S.bestA
    nA = size(A, 1)
    return Result(id, cls, alg, seed, S.found,
                  S.found ? S.tfound : time() - S.t0,
                  S.evals, nA, S.bestE,
                  nA <= 62 ? graph6(A) : "", edgelist(A), max(restarts, 0), A)
end
