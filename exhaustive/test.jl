# Pruebas de los invariantes con grafos de valores conocidos.
# Uso:  julia test.jl

include(joinpath(@__DIR__, "VerifConj.jl"))
using .VerifConj
using Test

completo(n) = grafo_desde_aristas(n, [(i, j) for i in 1:n for j in i+1:n])
camino(n) = grafo_desde_aristas(n, [(i, i + 1) for i in 1:n-1])
ciclo(n) = grafo_desde_aristas(n, [[(i, i + 1) for i in 1:n-1]; [(n, 1)]])
estrella(n) = grafo_desde_aristas(n, [(1, i) for i in 2:n])
# Petersen: ciclo exterior 1..5, estrella interior 6..10, radios i—i+5
petersen() = grafo_desde_aristas(10,
    [[(i, i % 5 + 1) for i in 1:5];
     [(5 + i, 5 + (i + 1) % 5 + 1) for i in 1:5];
     [(i, i + 5) for i in 1:5]])

@testset "VerifConj" begin
    @testset "graph6 (ida y vuelta)" begin
        for g in (completo(5), camino(7), ciclo(6), estrella(8), petersen())
            h = grafo_desde_graph6(escribir_graph6(g))
            @test h.adj == g.adj
        end
        @test escribir_graph6(completo(4)) == "C~"
    end

    @testset "sparse6 (contra graph6, convertidos con copyg)" begin
        for (s6, g6) in ((":An", "A_"), (":I`ESxol]~", "IhCS?C@?O"),
                         (":Fw@GsV", "F??Fw"), (":FkGM@Gs_QLN", "FCf~o"))
            @test grafo_desde_graph6(s6).adj == grafo_desde_graph6(g6).adj
        end
    end

    @testset "completo K_n" begin
        for n in 2:8
            g = completo(n)
            @test λ1(g) ≈ n - 1
            @test μ1_laplaciano(g) ≈ n
            @test numero_matching(g) == n ÷ 2
            @test numero_independencia(g) == 1
            @test diametro(g) == 1
        end
    end

    @testset "camino, ciclo y estrella" begin
        @test λ1(camino(3)) ≈ sqrt(2)
        @test λ1(ciclo(7)) ≈ 2
        @test diametro(camino(6)) == 5
        @test radio(camino(7)) == 3
        @test numero_matching(camino(7)) == 3
        @test numero_independencia(ciclo(7)) == 3
        for n in 3:9
            @test λ1(estrella(n)) ≈ sqrt(n - 1)
            @test numero_matching(estrella(n)) == 1
            @test numero_independencia(estrella(n)) == n - 1
        end
    end

    @testset "Petersen" begin
        g = petersen()
        @test grados(g) == fill(3, 10)
        esp = espectro_adyacencia(g)
        @test esp ≈ [-2, -2, -2, -2, 1, 1, 1, 1, 1, 3]
        @test numero_independencia(g) == 4
        @test numero_matching(g) == 5
        @test diametro(g) == 2
        @test ∂1_distancias(g) ≈ 15        # regular de distancias: suma de cada fila = 3·1 + 6·2
    end

    @testset "conexidad" begin
        g = grafo_desde_aristas(4, [(1, 2), (3, 4)])
        @test !conexo(g)
        @test conexo(ciclo(5))
    end
end
