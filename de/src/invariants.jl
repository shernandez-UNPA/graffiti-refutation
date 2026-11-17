# Invariantes de grafos usados por las conjeturas de Graffiti.
# Las definiciones siguen el código de referencia de Roucairol y Cazenave
# (github.com/RoucairolMilo/refutationGBR, src/models/conjectures/invariants.rs),
# para que los resultados sean comparables con su Tabla 1.
# El grafo se representa con una matriz de adyacencia A::Matrix{Float64} (0/1, simétrica).

using LinearAlgebra
using Statistics

nv(A) = size(A, 1)
nedges(A) = round(Int, sum(A) / 2)
degrees(A) = vec(sum(A, dims = 1))

"Autovalores (ordenados de menor a mayor) de una matriz simétrica."
spec(M) = eigvals(Symmetric(M))

"Matriz de distancias por BFS (el grafo se supone conexo)."
function dist_matrix(A)
    n = nv(A)
    D = fill(Inf, n, n)
    nbrs = [findall(!iszero, view(A, :, v)) for v in 1:n]
    queue = Vector{Int}(undef, n)
    for s in 1:n
        D[s, s] = 0.0
        head, tail = 1, 1
        queue[1] = s
        while head <= tail
            u = queue[head]; head += 1
            du = D[s, u]
            for w in nbrs[u]
                if D[s, w] == Inf
                    D[s, w] = du + 1
                    tail += 1
                    queue[tail] = w
                end
            end
        end
    end
    return D
end

function laplacian(A)
    L = -copy(A)
    d = degrees(A)
    for i in 1:nv(A)
        L[i, i] = d[i]
    end
    return L
end

"Desviación estándar poblacional (divide por n), como en el código de referencia."
popstd(v) = isempty(v) ? 0.0 : sqrt(sum(abs2, v .- mean(v)) / length(v))

"Índice de Randić: suma sobre aristas uv de 1/sqrt(d(u)d(v))."
function randic(A)
    d = degrees(A); n = nv(A); s = 0.0
    for i in 1:n, j in i+1:n
        A[i, j] == 1 && (s += 1 / sqrt(d[i] * d[j]))
    end
    return s
end

"Índice armónico: suma sobre aristas uv de 2/(d(u)+d(v))."
function harmonic(A)
    d = degrees(A); n = nv(A); s = 0.0
    for i in 1:n, j in i+1:n
        A[i, j] == 1 && (s += 2 / (d[i] + d[j]))
    end
    return s
end

"Temperatura de cada vértice: d(v)/(n - d(v))."
temperature(A) = (d = degrees(A); n = nv(A); d ./ (n .- d))

"Grado dual de v: promedio de los grados de sus vecinos."
function dual_degree(A, v)
    d = degrees(A)
    nb = findall(!iszero, view(A, :, v))
    return sum(d[nb]) / length(nb)
end
dual_degrees(A) = [dual_degree(A, v) for v in 1:nv(A)]

"""
Matriz de gravedad corregida (Brewster, Dinneen y Faber), la que usan
Roucairol y Cazenave: G[u,v] = d(u) d(v) / ((n-1) dist(u,v)) para u ≠ v.
"""
function gravity_matrix(A, D = dist_matrix(A))
    n = nv(A); d = degrees(A)
    G = zeros(n, n)
    for i in 1:n, j in 1:n
        i != j && (G[i, j] = d[i] * d[j] / ((n - 1) * D[i, j]))
    end
    return G
end

"Para cada vértice, cuántos vértices hay a distancia par (positiva)."
even_vec(D) = [count(j -> j != i && iseven(Int(D[i, j])), 1:size(D, 1)) for i in 1:size(D, 1)]
"Para cada vértice, cuántos vértices hay a distancia impar."
odd_vec(D) = [count(j -> j != i && isodd(Int(D[i, j])), 1:size(D, 1)) for i in 1:size(D, 1)]

"Deficiencia de v: pares de vecinos de v que no son adyacentes."
function deficiency(A, v)
    nb = findall(!iszero, view(A, :, v)); s = 0
    for a in 1:length(nb), b in 1:a-1
        A[nb[a], nb[b]] == 0 && (s += 1)
    end
    return s
end

"Diferencias consecutivas del vector ordenado (de mayor a menor)."
function vec_derivative(v)
    s = sort(v)
    return [s[end-i+1] - s[end-i] for i in 1:length(s)-1]
end

"""
Moda con tolerancia 1e-5, replicando el código de referencia: devuelve 0 si
todos los valores son distintos; ante empate gana el primer valor (en orden
creciente) que alcanza la multiplicidad máxima.
"""
function mode_tol(v)
    vals = Float64[]; cnt = Int[]
    best, bestval = 0, 0.0
    for x in v
        k = findfirst(y -> abs(y - x) < 1e-5, vals)
        if k === nothing
            push!(vals, x); push!(cnt, 1)
            if 1 > best
                best, bestval = 1, x
            end
        else
            cnt[k] += 1
            if cnt[k] > best
                best, bestval = cnt[k], vals[k]
            end
        end
    end
    return best == 1 ? 0.0 : bestval
end

"Segundo mayor autovalor (λ ordenados de menor a mayor)."
second_largest(λ) = λ[end-1]
"Segundo menor autovalor."
second_smallest(λ) = λ[2]
