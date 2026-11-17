# Codificación por claves aleatorias (random keys).
#
# Un individuo para orden n es x ∈ [0,1]^(1+P(n)), con P(n) = n(n-1)/2:
#   x[1]      = umbral t (controla la densidad),
#   x[1+k]    = clave del par k (pares ordenados por su vértice mayor j, así
#               el genoma de orden n es un prefijo del de orden n+1).
#
# Decodificación (tipo Kruskal): se recorren los pares por clave decreciente.
# Un par se agrega si su clave supera t o si une dos componentes, y siempre
# que sea factible para la clase. Resultado: un grafo conexo y factible por
# construcción en las cuatro clases.
#   :any     cualquier grafo conexo
#   :tree    árbol (árbol generador de máxima clave)
#   :k3free  sin triángulos
#   :girth5  cintura >= 5 (sin C3 ni C4)

npairs(n) = n * (n - 1) ÷ 2
genome_length(n) = 1 + npairs(n)

"Tabla de pares (i,j), i<j, ordenada por j y luego i, hasta nmax vértices."
function pair_table(nmax)
    P = Vector{Tuple{Int,Int}}(undef, npairs(nmax))
    k = 0
    for j in 2:nmax, i in 1:j-1
        k += 1
        P[k] = (i, j)
    end
    return P
end

mutable struct Decoder
    pairs::Vector{Tuple{Int,Int}}
    parent::Vector{Int}
    order::Vector{Int}
end
Decoder(nmax) = Decoder(pair_table(nmax), zeros(Int, nmax), Int[])

function uf_find(p, x)
    while p[x] != x
        p[x] = p[p[x]]
        x = p[x]
    end
    return x
end

"¿Agregar la arista ij mantiene la clase?"
function feasible(cls, A, i, j)
    cls === :any && return true
    n = size(A, 1)
    if cls === :k3free
        @inbounds for w in 1:n
            (A[i, w] == 1 && A[j, w] == 1) && return false
        end
        return true
    elseif cls === :girth5
        # Prohibido un camino de longitud 2 o 3 entre i y j.
        @inbounds for a in 1:n
            A[i, a] == 1 || continue
            A[a, j] == 1 && return false
            for b in 1:n
                (A[j, b] == 1 && A[a, b] == 1) && return false
            end
        end
        return true
    end
    return true
end

"Decodifica x en la matriz de adyacencia de un grafo de orden n."
function decode(dec::Decoder, x::AbstractVector, n::Int, cls::Symbol)
    P = npairs(n)
    t = x[1]
    keys = view(x, 2:1+P)
    resize!(dec.order, P)
    sortperm!(dec.order, keys; rev = true)
    p = dec.parent
    for v in 1:n
        p[v] = v
    end
    A = zeros(n, n)
    ncomp = n
    for k in dec.order
        i, j = dec.pairs[k]
        ri, rj = uf_find(p, i), uf_find(p, j)
        joins = ri != rj
        if cls === :tree
            joins || continue
        else
            (keys[k] > t || joins) || continue
            feasible(cls, A, i, j) || continue
        end
        A[i, j] = 1.0; A[j, i] = 1.0
        if joins
            p[ri] = rj
            ncomp -= 1
            cls === :tree && ncomp == 1 && break
        end
    end
    return A
end

"Formato graph6 (estándar de nauty y House of Graphs), para n <= 62."
function graph6(A)
    n = size(A, 1)
    @assert n <= 62
    bits = Bool[]
    for j in 2:n, i in 1:j-1
        push!(bits, A[i, j] == 1)
    end
    while length(bits) % 6 != 0
        push!(bits, false)
    end
    out = Char[Char(63 + n)]
    for k in 1:6:length(bits)
        v = 0
        for b in 0:5
            v = 2v + bits[k+b]
        end
        push!(out, Char(63 + v))
    end
    return String(out)
end

"Lista de aristas en texto (vértices numerados desde 0)."
edgelist(A) = join(["$(i-1)-$(j-1)" for j in 2:size(A, 1) for i in 1:j-1 if A[i, j] == 1], " ")
