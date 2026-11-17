# Pruebas de los invariantes y de los decodificadores.
# Uso: julia test/runtests.jl

include(joinpath(@__DIR__, "..", "src", "GraffitiDE.jl"))
using .GraffitiDE
using .GraffitiDE: nedges, dist_matrix, laplacian, temperature, even_vec, odd_vec, Decoder, decode
using Test, LinearAlgebra, Random

complete(n) = ones(n, n) - I
function path(n)
    A = zeros(n, n)
    for i in 1:n-1
        A[i, i+1] = A[i+1, i] = 1
    end
    A
end
function cycle(n)
    A = path(n); A[1, n] = A[n, 1] = 1; A
end
function star(n)  # n vértices, centro 1
    A = zeros(n, n)
    for i in 2:n
        A[1, i] = A[i, 1] = 1
    end
    A
end

istree(A) = nedges(A) == size(A, 1) - 1 && all(isfinite, dist_matrix(A))
connected(A) = all(isfinite, dist_matrix(A))
has_triangle(A) = tr(A^3) > 0.5
# cintura >= 5 si no hay triángulos ni 4-ciclos: (A^2)_{ij} <= 1 para i != j
function girth_ge5(A)
    has_triangle(A) && return false
    B = A^2
    n = size(A, 1)
    return all(B[i, j] <= 1.5 for i in 1:n, j in 1:n if i != j)
end

@testset "invariantes" begin
    @test spec(complete(4)) ≈ [-1, -1, -1, 3]
    @test dist_matrix(path(4))[1, 4] == 3
    @test GraffitiDE.harmonic(path(3)) ≈ 4 / 3
    @test GraffitiDE.randic(star(5)) ≈ 4 / 2
    @test temperature(complete(3)) == [2.0, 2.0, 2.0]
    @test even_vec(dist_matrix(cycle(6))) == fill(2, 6)
    @test odd_vec(dist_matrix(cycle(6))) == fill(3, 6)
    @test graph6(complete(3)) == "Bw"
    @test graph6(path(4)) == "Ch"
    @test length(CONJECTURES) == 35
end

@testset "decodificadores" begin
    rng = Xoshiro(7)
    dec = Decoder(30)
    for n in (5, 12, 25), _ in 1:50
        x = rand(rng, 1 + n * (n - 1) ÷ 2)
        @test connected(decode(dec, x, n, :any))
        @test istree(decode(dec, x, n, :tree))
        A = decode(dec, x, n, :k3free); @test connected(A) && !has_triangle(A)
        A = decode(dec, x, n, :girth5); @test connected(A) && girth_ge5(A)
    end
end

@testset "contraejemplo publicado: Graffiti 197 con C17" begin
    # Roucairol y Cazenave (2025) refutan Graffiti 197 con el ciclo de 17 vértices.
    @test violation(CONJECTURES[197], cycle(17)) > TOL
end

@testset "todas las conjeturas se evalúan sin error" begin
    rng = Xoshiro(3)
    dec = Decoder(20)
    for (id, c) in CONJECTURES, cls in c.classes
        x = rand(rng, 1 + 15 * 14 ÷ 2)
        E = violation(c, decode(dec, x, 15, cls))
        @test !isnan(E)
    end
end
