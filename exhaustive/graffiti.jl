#!/usr/bin/env julia
#
# Verificación exhaustiva de las 35 conjeturas de Graffiti de la Tabla 1 de
# Roucairol y Cazenave (2025), con las MISMAS definiciones que el experimento
# de evolución diferencial del artículo para TAG (carpeta graffiti/, copiada de
# 2026 - TAG/GraffitiDE/src sin cambios).
#
# Cada conjetura es f(G) ≤ g(G); violation devuelve E(G) = f(G) − g(G) y un grafo
# es contraejemplo si E(G) > TOL = 1e-4 (igual que en GraffitiDE).
#
# La clase de grafos la fija el generador; hay que pasar la misma con --clase:
#
#   any     nauty-geng -c n          (conexos)
#   tree    nauty-gentreeg n         (árboles; sale en sparse6, se lee igual)
#   k3free  nauty-geng -ct n         (conexos sin triángulos)
#   girth5  nauty-geng -ctf n        (conexos sin triángulos ni 4-ciclos: cintura ≥ 5)
#
# Ejemplos:
#   for n in $(seq 3 9); do nauty-geng -cq $n; done | julia -t auto graffiti.jl --clase any
#   for n in $(seq 4 16); do nauty-gentreeg -q $n; done | julia -t auto graffiti.jl --clase tree --salida res_tree
#   nauty-geng -cq 10 | julia -t auto graffiti.jl --clase any --ids 29,30 --salida res_29_30_n10
#
# Opciones:
#   --clase any|tree|k3free|girth5   (obligatoria)
#   --ids 29,30,166                  solo esas conjeturas (por defecto, todas las de la clase)
#   --estado R|O                     solo refutadas (R) o abiertas (O) según la Tabla 1
#   --salida DIR                     resumen.md + contraejemplos .g6
#   --entrada ARCHIVO                graph6 desde archivo en vez de la entrada estándar
#   --silencioso

include(joinpath(@__DIR__, "VerifConj.jl"))
using .VerifConj
include(joinpath(@__DIR__, "..", "de", "src", "invariants.jl"))
include(joinpath(@__DIR__, "..", "de", "src", "conjectures.jl"))
using Printf

function opciones(args)
    op = Dict{String,String}()
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--silencioso"
            op["silencioso"] = "1"; i += 1
        elseif startswith(a, "--") && i < length(args)
            op[a[3:end]] = args[i+1]; i += 2
        else
            error("opción no reconocida: $a")
        end
    end
    haskey(op, "clase") || error("falta --clase (any, tree, k3free o girth5)")
    return op
end

function seleccionar(op)
    clase = Symbol(op["clase"])
    clase in (:any, :tree, :k3free, :girth5) || error("clase desconocida: $clase")
    cs = sort([c for c in values(CONJECTURES) if clase in c.classes]; by = c -> c.id)
    if haskey(op, "ids")
        ids = parse.(Int, split(op["ids"], ","))
        cs = filter(c -> c.id in ids, cs)
    end
    if haskey(op, "estado")
        e = only(op["estado"])
        cs = filter(c -> c.status == e, cs)
    end
    isempty(cs) && error("ninguna conjetura de la Tabla 1 corresponde a esa selección")
    return clase, cs
end

function main(args)
    op = opciones(args)
    clase, cs = seleccionar(op)
    conjs = [Conjetura("G$(c.id)_$(clase)",
                       "[$(c.status)] $(c.statement)",
                       g -> -violation(c, matriz_adyacencia(g)))
             for c in cs]
    println("Clase: $clase   Conjeturas: ", join([c.id for c in cs], ", "))
    io = haskey(op, "entrada") ? open(op["entrada"]) : stdin
    t0 = time()
    # margen = −E(G); contraejemplo si margen < −TOL, es decir E(G) > TOL.
    total, res = verificar_flujo(io, conjs; ε = TOL, progreso = !haskey(op, "silencioso"))
    dt = time() - t0
    imprimir_resultados(stdout, total, res, dt)
    if haskey(op, "salida")
        escribir_resultados(op["salida"], total, res, dt)
        println("Resultados escritos en ", op["salida"])
    end
end

main(ARGS)
