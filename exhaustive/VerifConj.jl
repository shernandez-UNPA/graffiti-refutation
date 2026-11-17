"""
    VerifConj

Verificación exhaustiva de conjeturas sobre invariantes de grafos.

Lee grafos en formato graph6 (por ejemplo, la salida de `nauty-geng`), calcula
invariantes y evalúa cada conjetura en todos los grafos, en paralelo con hilos.

Solo usa la biblioteca estándar de Julia (LinearAlgebra, Printf).
Representación interna: listas de adyacencia como máscaras de bits (`UInt64`),
por lo que admite grafos de hasta 62 vértices (de sobra para búsqueda exhaustiva).
"""
module VerifConj

using LinearAlgebra
using Printf

export Grafo, leer_graph6!, grafo_desde_graph6, escribir_graph6, grafo_desde_aristas,
       n_vertices, n_aristas, grados, grado_max, grado_min, matriz_adyacencia,
       espectro_adyacencia, λ1, espectro_laplaciano, μ1_laplaciano,
       espectro_laplaciano_sin_signo, q1_sin_signo, espectro_distancias, ∂1_distancias,
       matriz_distancias, conexo, diametro, radio, numero_matching, numero_independencia,
       Conjetura, Resultado, verificar_flujo, imprimir_resultados, escribir_resultados

# ---------------------------------------------------------------------------
# Estructura de grafo con cachés de invariantes
# ---------------------------------------------------------------------------

"""
    Grafo

Grafo simple no dirigido con `n ≤ 62` vértices. Guarda buffers reutilizables y
cachés de invariantes que se invalidan cada vez que se carga un grafo nuevo con
`leer_graph6!`. Un `Grafo` no debe compartirse entre hilos.
"""
mutable struct Grafo
    n::Int
    adj::Vector{UInt64}          # adj[i+1] = vecinos del vértice i (bits 0-based)
    M::Matrix{Float64}           # buffer para matrices reales n×n
    A::Matrix{Float64}           # matriz de adyacencia (caché, no se modifica)
    okAmat::Bool
    D::Matrix{Int}               # matriz de distancias
    memo::Vector{Int8}           # memo para matching / independencia (2^n)
    espA::Vector{Float64}; okA::Bool
    espL::Vector{Float64}; okL::Bool
    espQ::Vector{Float64}; okQ::Bool
    espD::Vector{Float64}; okDesp::Bool
    okD::Bool; conexo::Bool
    matching::Int
    indep::Int
end

Grafo() = Grafo(0, UInt64[], zeros(0, 0), zeros(0, 0), false, zeros(Int, 0, 0), Int8[],
                Float64[], false, Float64[], false, Float64[], false, Float64[], false,
                false, false, -1, -1)

function _preparar!(g::Grafo, n::Int)
    n > 62 && throw(ArgumentError("solo se admiten grafos de hasta 62 vértices (n = $n)"))
    if g.n != n
        g.n = n
        g.adj = zeros(UInt64, n)
        g.M = zeros(n, n)
        g.A = zeros(n, n)
        g.D = zeros(Int, n, n)
        g.memo = n <= 24 ? Vector{Int8}(undef, 1 << n) : Int8[]
    else
        fill!(g.adj, 0)
    end
    g.okA = g.okL = g.okQ = g.okDesp = g.okD = g.okAmat = false
    g.matching = -1
    g.indep = -1
    return g
end

# ---------------------------------------------------------------------------
# Formato graph6
# ---------------------------------------------------------------------------

"""
    leer_graph6!(g, s) -> g

Carga en `g` el grafo codificado en la cadena `s`, en formato **graph6** o
**sparse6** (este último empieza con `:`; es el que produce `nauty-gentreeg`).
No admite el encabezado opcional `>>graph6<<`. Admite n ≤ 62.
"""
function leer_graph6!(g::Grafo, s::AbstractString)
    b = codeunits(s)
    b[1] == UInt8(':') && return _leer_sparse6!(g, b)
    n = Int(b[1]) - 63
    (0 <= n <= 62) || throw(ArgumentError("graph6 inválido o n > 62: $(repr(s))"))
    _preparar!(g, n)
    adj = g.adj
    pos = 2
    bits = 0
    quedan = 0
    @inbounds for j in 1:n-1, i in 0:j-1
        if quedan == 0
            bits = Int(b[pos]) - 63
            pos += 1
            quedan = 6
        end
        quedan -= 1
        if (bits >> quedan) & 1 == 1
            adj[i+1] |= UInt64(1) << j
            adj[j+1] |= UInt64(1) << i
        end
    end
    return g
end

# Decodificador sparse6 (formato descrito en formats.txt de nauty; mismo
# algoritmo que el lector de networkx). Ignora aristas repetidas y lazos.
function _leer_sparse6!(g::Grafo, b)
    n = Int(b[2]) - 63
    (0 <= n <= 62) || throw(ArgumentError("sparse6 inválido o n > 62"))
    _preparar!(g, n)
    adj = g.adj
    k = 1
    while (1 << k) < n
        k += 1
    end
    pos = 3
    L = length(b)
    d = 0
    dlen = 0
    v = 0
    @inbounds while true
        if dlen < 1
            pos > L && break
            d = Int(b[pos]) - 63; pos += 1
            dlen = 6
        end
        dlen -= 1
        bit = (d >> dlen) & 1
        x = d & ((1 << dlen) - 1)
        xlen = dlen
        completo = true
        while xlen < k
            if pos > L
                completo = false
                break
            end
            d = Int(b[pos]) - 63; pos += 1
            x = (x << 6) + d
            xlen += 6
        end
        completo || break
        x >>= (xlen - k)
        dlen = xlen - k
        d &= (1 << dlen) - 1
        bit == 1 && (v += 1)
        if x >= n || v >= n
            break
        elseif x > v
            v = x
        elseif x != v
            adj[x+1] |= UInt64(1) << v
            adj[v+1] |= UInt64(1) << x
        end
    end
    return g
end

grafo_desde_graph6(s::AbstractString) = leer_graph6!(Grafo(), s)

"""
    escribir_graph6(g) -> String

Codifica el grafo `g` en formato graph6.
"""
function escribir_graph6(g::Grafo)
    n = g.n
    out = UInt8[UInt8(n + 63)]
    acc = 0
    k = 0
    for j in 1:n-1, i in 0:j-1
        acc = (acc << 1) | Int((g.adj[i+1] >> j) & 1)
        k += 1
        if k == 6
            push!(out, UInt8(acc + 63)); acc = 0; k = 0
        end
    end
    if k > 0
        push!(out, UInt8((acc << (6 - k)) + 63))
    end
    return String(out)
end

"""
    grafo_desde_aristas(n, aristas) -> Grafo

Construye un grafo con vértices `1:n` a partir de una lista de pares `(u, v)`.
"""
function grafo_desde_aristas(n::Int, aristas)
    g = _preparar!(Grafo(), n)
    for (u, v) in aristas
        u == v && throw(ArgumentError("no se admiten lazos"))
        g.adj[u] |= UInt64(1) << (v - 1)
        g.adj[v] |= UInt64(1) << (u - 1)
    end
    return g
end

# ---------------------------------------------------------------------------
# Invariantes básicos
# ---------------------------------------------------------------------------

n_vertices(g::Grafo) = g.n
n_aristas(g::Grafo) = sum(count_ones, g.adj; init = 0) ÷ 2
grados(g::Grafo) = [count_ones(a) for a in g.adj]
grado_max(g::Grafo) = maximum(count_ones, g.adj; init = 0)
grado_min(g::Grafo) = minimum(count_ones, g.adj; init = 0)

function _llenar_adyacencia!(M, g::Grafo)
    n = g.n
    @inbounds for j in 1:n, i in 1:n
        M[i, j] = Float64((g.adj[i] >> (j - 1)) & 1)
    end
    return M
end

"""
    matriz_adyacencia(g) -> Matrix{Float64}

Matriz de adyacencia 0/1 de `g`. Se calcula una vez por grafo y se reutiliza:
**no debe modificarse** (copiarla si hace falta).
"""
function matriz_adyacencia(g::Grafo)
    if !g.okAmat
        _llenar_adyacencia!(g.A, g)
        g.okAmat = true
    end
    return g.A
end

# ---------------------------------------------------------------------------
# Invariantes espectrales (autovalores en orden creciente)
# ---------------------------------------------------------------------------

"Autovalores de la matriz de adyacencia A, en orden creciente."
function espectro_adyacencia(g::Grafo)
    if !g.okA
        _llenar_adyacencia!(g.M, g)
        g.espA = eigvals!(Symmetric(g.M))
        g.okA = true
    end
    return g.espA
end

"Radio espectral (mayor autovalor) de la matriz de adyacencia."
λ1(g::Grafo) = g.n == 0 ? 0.0 : espectro_adyacencia(g)[end]

function _laplaciano!(M, g::Grafo, signo::Float64)
    _llenar_adyacencia!(M, g)
    @inbounds for i in 1:g.n
        for j in 1:g.n
            M[i, j] *= signo
        end
        M[i, i] = count_ones(g.adj[i])
    end
    return M
end

"Autovalores del laplaciano L = D − A, en orden creciente."
function espectro_laplaciano(g::Grafo)
    if !g.okL
        _laplaciano!(g.M, g, -1.0)
        g.espL = eigvals!(Symmetric(g.M))
        g.okL = true
    end
    return g.espL
end

"Mayor autovalor del laplaciano."
μ1_laplaciano(g::Grafo) = g.n == 0 ? 0.0 : espectro_laplaciano(g)[end]

"Autovalores del laplaciano sin signo Q = D + A, en orden creciente."
function espectro_laplaciano_sin_signo(g::Grafo)
    if !g.okQ
        _laplaciano!(g.M, g, 1.0)
        g.espQ = eigvals!(Symmetric(g.M))
        g.okQ = true
    end
    return g.espQ
end

"Mayor autovalor del laplaciano sin signo."
q1_sin_signo(g::Grafo) = g.n == 0 ? 0.0 : espectro_laplaciano_sin_signo(g)[end]

# ---------------------------------------------------------------------------
# Distancias (BFS sobre máscaras de bits)
# ---------------------------------------------------------------------------

const INF = typemax(Int) ÷ 4

function _distancias!(g::Grafo)
    g.okD && return g.D
    n = g.n
    D = g.D
    fill!(D, INF)
    todos_alcanzados = true
    @inbounds for s in 0:n-1
        D[s+1, s+1] = 0
        visitados = UInt64(1) << s
        frontera = visitados
        d = 0
        while frontera != 0
            d += 1
            siguiente = UInt64(0)
            f = frontera
            while f != 0
                v = trailing_zeros(f)
                f &= f - 1
                siguiente |= g.adj[v+1]
            end
            siguiente &= ~visitados
            visitados |= siguiente
            frontera = siguiente
            t = siguiente
            while t != 0
                v = trailing_zeros(t)
                t &= t - 1
                D[s+1, v+1] = d
            end
        end
        if count_ones(visitados) != n
            todos_alcanzados = false
        end
    end
    g.conexo = n == 0 || todos_alcanzados
    g.okD = true
    return D
end

"Matriz de distancias (entradas `INF` si el grafo no es conexo)."
matriz_distancias(g::Grafo) = _distancias!(g)

conexo(g::Grafo) = (_distancias!(g); g.conexo)

"Diámetro (lanza error si el grafo no es conexo)."
function diametro(g::Grafo)
    conexo(g) || throw(DomainError(g.n, "el diámetro no está definido para grafos no conexos"))
    return maximum(g.D)
end

"Radio (mínima excentricidad; lanza error si el grafo no es conexo)."
function radio(g::Grafo)
    conexo(g) || throw(DomainError(g.n, "el radio no está definido para grafos no conexos"))
    return minimum(maximum(view(g.D, i, :)) for i in 1:g.n)
end

"Autovalores de la matriz de distancias (solo grafos conexos), en orden creciente."
function espectro_distancias(g::Grafo)
    if !g.okDesp
        conexo(g) || throw(DomainError(g.n, "matriz de distancias no definida: grafo no conexo"))
        @inbounds for j in 1:g.n, i in 1:g.n
            g.M[i, j] = g.D[i, j]
        end
        g.espD = eigvals!(Symmetric(g.M))
        g.okDesp = true
    end
    return g.espD
end

"Radio espectral de la matriz de distancias."
∂1_distancias(g::Grafo) = espectro_distancias(g)[end]

# ---------------------------------------------------------------------------
# Invariantes combinatorios exactos (programación dinámica sobre subconjuntos)
# Costo O(2^n · Δ): perfectamente viable para n ≤ 12.
# ---------------------------------------------------------------------------

function _matching(adj, memo, S::UInt64)
    S == 0 && return 0
    @inbounds m = memo[S+1]
    m >= 0 && return Int(m)
    v = trailing_zeros(S)
    S2 = S & ~(UInt64(1) << v)
    mejor = _matching(adj, memo, S2)
    @inbounds N = adj[v+1] & S2
    while N != 0
        u = trailing_zeros(N)
        N &= N - 1
        mejor = max(mejor, 1 + _matching(adj, memo, S2 & ~(UInt64(1) << u)))
    end
    @inbounds memo[S+1] = mejor
    return mejor
end

"Número de matching μ(G): tamaño de un emparejamiento máximo (exacto)."
function numero_matching(g::Grafo)
    if g.matching < 0
        g.n <= 24 || throw(ArgumentError("numero_matching exacto implementado para n ≤ 24"))
        fill!(g.memo, Int8(-1))
        g.matching = _matching(g.adj, g.memo, (UInt64(1) << g.n) - 1)
    end
    return g.matching
end

function _indep(adj, memo, S::UInt64)
    S == 0 && return 0
    @inbounds m = memo[S+1]
    m >= 0 && return Int(m)
    v = trailing_zeros(S)
    sin_v = S & ~(UInt64(1) << v)
    @inbounds con_v = S & ~(adj[v+1] | (UInt64(1) << v))
    mejor = max(_indep(adj, memo, sin_v), 1 + _indep(adj, memo, con_v))
    @inbounds memo[S+1] = mejor
    return mejor
end

"Número de independencia α(G) (exacto)."
function numero_independencia(g::Grafo)
    if g.indep < 0
        g.n <= 24 || throw(ArgumentError("numero_independencia exacto implementado para n ≤ 24"))
        fill!(g.memo, Int8(-1))
        g.indep = _indep(g.adj, g.memo, (UInt64(1) << g.n) - 1)
    end
    return g.indep
end

# ---------------------------------------------------------------------------
# Conjeturas y motor de verificación
# ---------------------------------------------------------------------------

"""
    Conjetura(nombre, descripcion, margen)

`margen(g)` debe devolver un número tal que la conjetura afirma `margen(g) ≥ 0`
para todo grafo `g` de la clase considerada. Por ejemplo, para `λ₁ ≥ 2m/n`
el margen es `λ1(g) - 2n_aristas(g)/n_vertices(g)`.
"""
struct Conjetura
    nombre::String
    descripcion::String
    margen::Function
end

mutable struct Acumulador
    revisados::Int
    contraejemplos::Int
    dudosos::Int
    margen_min::Float64
    g6_min::String
    n_min_contra::Int
    ejemplos::Vector{String}
end

Acumulador() = Acumulador(0, 0, 0, Inf, "", typemax(Int), String[])

"""
Resultado de una conjetura:
- `revisados`: grafos evaluados.
- `contraejemplos`: grafos con margen < −ε.
- `dudosos`: grafos con |margen| ≤ ε (igualdad numérica: revisar con más precisión).
- `margen_min`, `g6_min`: menor margen observado y un grafo que lo alcanza.
- `n_min_contra`: menor orden en que apareció un contraejemplo (`typemax(Int)` si no hubo).
- `ejemplos`: primeros contraejemplos encontrados (graph6).
"""
struct Resultado
    conjetura::Conjetura
    revisados::Int
    contraejemplos::Int
    dudosos::Int
    margen_min::Float64
    g6_min::String
    n_min_contra::Int
    ejemplos::Vector{String}
end

function _procesar_bloque(lineas, rango, conjeturas, ε, max_ejemplos)
    g = Grafo()
    acs = [Acumulador() for _ in conjeturas]
    for idx in rango
        s = lineas[idx]
        leer_graph6!(g, s)
        for (k, c) in enumerate(conjeturas)
            ac = acs[k]
            m = Float64(c.margen(g))
            ac.revisados += 1
            if m < ac.margen_min
                ac.margen_min = m
                ac.g6_min = s
            end
            if m < -ε
                ac.contraejemplos += 1
                ac.n_min_contra = min(ac.n_min_contra, g.n)
                length(ac.ejemplos) < max_ejemplos && push!(ac.ejemplos, s)
            elseif m <= ε
                ac.dudosos += 1
            end
        end
    end
    return acs
end

function _combinar!(a::Acumulador, b::Acumulador, max_ejemplos)
    a.revisados += b.revisados
    a.contraejemplos += b.contraejemplos
    a.dudosos += b.dudosos
    if b.margen_min < a.margen_min
        a.margen_min = b.margen_min
        a.g6_min = b.g6_min
    end
    a.n_min_contra = min(a.n_min_contra, b.n_min_contra)
    for e in b.ejemplos
        length(a.ejemplos) >= max_ejemplos && break
        push!(a.ejemplos, e)
    end
    return a
end

"""
    verificar_flujo(io, conjeturas; ε = 1e-9, lote = 50_000, max_ejemplos = 1000,
                    progreso = true) -> (total, resultados)

Lee grafos graph6 de `io` (una línea por grafo; ignora líneas vacías y las que
empiezan con `>`), los procesa por lotes repartidos entre los hilos disponibles
(`julia -t auto`) y devuelve el total de grafos y un `Resultado` por conjetura.
No guarda los grafos: memoria acotada independientemente de cuántos haya.
"""
function verificar_flujo(io::IO, conjeturas::Vector{Conjetura};
                         ε::Float64 = 1e-9, lote::Int = 50_000,
                         max_ejemplos::Int = 1000, progreso::Bool = true)
    globales = [Acumulador() for _ in conjeturas]
    lineas = String[]
    sizehint!(lineas, lote)
    total = 0
    t0 = time()
    nt = Threads.nthreads()
    while !eof(io)
        empty!(lineas)
        while length(lineas) < lote && !eof(io)
            l = readline(io)
            (isempty(l) || first(l) == '>') && continue
            push!(lineas, l)
        end
        isempty(lineas) && break
        tam = cld(length(lineas), nt)
        tareas = [Threads.@spawn _procesar_bloque(lineas, r, conjeturas, ε, max_ejemplos)
                  for r in Iterators.partition(eachindex(lineas), tam)]
        for t in tareas
            parcial = fetch(t)
            for k in eachindex(globales)
                _combinar!(globales[k], parcial[k], max_ejemplos)
            end
        end
        total += length(lineas)
        if progreso
            dt = time() - t0
            @printf(stderr, "\r  %d grafos procesados (%.0f grafos/s)", total, total / max(dt, 1e-9))
        end
    end
    progreso && println(stderr)
    resultados = [Resultado(c, a.revisados, a.contraejemplos, a.dudosos, a.margen_min,
                            a.g6_min, a.n_min_contra, a.ejemplos)
                  for (c, a) in zip(conjeturas, globales)]
    return total, resultados
end

# ---------------------------------------------------------------------------
# Salida
# ---------------------------------------------------------------------------

function imprimir_resultados(io::IO, total, resultados, segundos)
    @printf(io, "\nGrafos procesados: %d   Tiempo: %.2f s   Hilos: %d\n\n",
            total, segundos, Threads.nthreads())
    @printf(io, "%-18s %14s %14s %12s %14s  %s\n",
            "conjetura", "contraejemplos", "dudosos(|m|≤ε)", "n mínimo", "margen mínimo", "grafo con margen mínimo")
    println(io, "-"^110)
    for r in resultados
        nmin = r.n_min_contra == typemax(Int) ? "—" : string(r.n_min_contra)
        @printf(io, "%-18s %14d %14d %12s %14.6g  %s\n",
                r.conjetura.nombre, r.contraejemplos, r.dudosos, nmin, r.margen_min, r.g6_min)
    end
    println(io)
end

"""
    escribir_resultados(dir, total, resultados, segundos)

Escribe en `dir` un `resumen.md`, un `resumen.csv` y, por cada conjetura con contraejemplos, un
archivo `<nombre>.g6` con los primeros contraejemplos encontrados.
"""
function escribir_resultados(dir::AbstractString, total, resultados, segundos)
    mkpath(dir)
    open(joinpath(dir, "resumen.md"), "w") do io
        println(io, "# Resumen de verificación\n")
        @printf(io, "- Grafos procesados: %d\n- Tiempo: %.2f s\n- Hilos: %d\n\n",
                total, segundos, Threads.nthreads())
        println(io, "| Conjetura | Descripción | Contraejemplos | Dudosos | n mínimo | Margen mínimo | Grafo (graph6) |")
        println(io, "|---|---|---|---|---|---|---|")
        for r in resultados
            nmin = r.n_min_contra == typemax(Int) ? "—" : string(r.n_min_contra)
            @printf(io, "| %s | %s | %d | %d | %s | %.6g | `%s` |\n",
                    r.conjetura.nombre, r.conjetura.descripcion, r.contraejemplos,
                    r.dudosos, nmin, r.margen_min, r.g6_min)
        end
    end
    open(joinpath(dir, "resumen.csv"), "w") do io
        println(io, "conjetura,revisados,contraejemplos,dudosos,n_min,margen_min,g6_min")
        for r in resultados
            nmin = r.n_min_contra == typemax(Int) ? "" : string(r.n_min_contra)
            # graph6/sparse6 pueden contener comas: se citan entre comillas
            g6 = replace(r.g6_min, "\"" => "\"\"")
            @printf(io, "%s,%d,%d,%d,%s,%.10g,\"%s\"\n", r.conjetura.nombre, r.revisados,
                    r.contraejemplos, r.dudosos, nmin, r.margen_min, g6)
        end
    end
    for r in resultados
        isempty(r.ejemplos) && continue
        open(joinpath(dir, r.conjetura.nombre * ".g6"), "w") do io
            foreach(e -> println(io, e), r.ejemplos)
        end
    end
end

end # module
