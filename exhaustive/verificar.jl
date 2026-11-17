#!/usr/bin/env julia
#
# Uso (en WSL / Linux):
#
#   nauty-geng -c 9 | julia -t auto verificar.jl
#   nauty-geng -c 10 | julia -t auto verificar.jl --salida resultados_n10
#   for n in 3 4 5 6 7 8 9 10; do nauty-geng -cq $n; done | julia -t auto verificar.jl
#
# Opciones:
#   --conjeturas ARCHIVO   archivo con la lista CONJETURAS (por defecto conjeturas.jl)
#   --salida DIR           carpeta donde escribir resumen.md y los contraejemplos (.g6)
#   --eps X                tolerancia numérica (por defecto 1e-9)
#   --entrada ARCHIVO      leer graph6 de un archivo en vez de la entrada estándar
#   --silencioso           no mostrar el progreso

include(joinpath(@__DIR__, "VerifConj.jl"))
using .VerifConj

function opciones(args)
    op = Dict{String,String}("conjeturas" => joinpath(@__DIR__, "conjeturas.jl"),
                             "eps" => "1e-9")
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--silencioso"
            op["silencioso"] = "1"
            i += 1
        elseif startswith(a, "--") && i < length(args)
            op[a[3:end]] = args[i+1]
            i += 2
        else
            error("opción no reconocida: $a")
        end
    end
    return op
end

function main(args)
    op = opciones(args)
    include(op["conjeturas"])
    conjs = Base.invokelatest(() -> CONJETURAS)
    ε = parse(Float64, op["eps"])
    io = haskey(op, "entrada") ? open(op["entrada"]) : stdin
    t0 = time()
    total, res = Base.invokelatest(verificar_flujo, io, conjs;
                                   ε = ε, progreso = !haskey(op, "silencioso"))
    dt = time() - t0
    imprimir_resultados(stdout, total, res, dt)
    if haskey(op, "salida")
        escribir_resultados(op["salida"], total, res, dt)
        println("Resultados escritos en ", op["salida"])
    end
end

main(ARGS)
