#!/usr/bin/env julia
#
# Combina las porciones de un orden (salida de correr_por_partes.sh) en un único
# resultado con el mismo formato que una corrida entera:
#   resumen.csv, resumen.md y un .g6 por conjetura (primeros 1000 contraejemplos).
#
# Uso:  julia combinar_partes.jl DIR_PARTES DIR_FINAL

using Printf

function campos_csv(linea)
    out = String[]; buf = IOBuffer(); comillas = false
    cs = collect(linea); i = 1
    while i <= length(cs)
        c = cs[i]
        if comillas
            if c == '"' && i < length(cs) && cs[i+1] == '"'
                write(buf, '"'); i += 1
            elseif c == '"'
                comillas = false
            else
                write(buf, c)
            end
        elseif c == '"'
            comillas = true
        elseif c == ','
            push!(out, String(take!(buf)))
        else
            write(buf, c)
        end
        i += 1
    end
    push!(out, String(take!(buf)))
    return out
end

mutable struct Total
    revisados::Int
    contra::Int
    dudosos::Int
    nmin::Int
    margen::Float64
    g6::String
end

function main(args)
    length(args) == 2 || error("uso: julia combinar_partes.jl DIR_PARTES DIR_FINAL")
    base, final = args
    partes = parse(Int, strip(read(joinpath(base, "PARTES"), String)))
    dirs = [joinpath(base, @sprintf("p%02dde%02d", r, partes)) for r in 0:partes-1]
    faltan = [d for d in dirs if !isfile(joinpath(d, "resumen.csv"))]
    isempty(faltan) || error("faltan $(length(faltan)) de $partes porciones, por ejemplo $(faltan[1])")

    orden = String[]
    tot = Dict{String,Total}()
    for d in dirs
        ls = readlines(joinpath(d, "resumen.csv"))
        for l in ls[2:end]
            f = campos_csv(l)
            nombre = f[1]
            if !haskey(tot, nombre)
                push!(orden, nombre)
                tot[nombre] = Total(0, 0, 0, typemax(Int), Inf, "")
            end
            t = tot[nombre]
            t.revisados += parse(Int, f[2])
            t.contra += parse(Int, f[3])
            t.dudosos += parse(Int, f[4])
            isempty(f[5]) || (t.nmin = min(t.nmin, parse(Int, f[5])))
            m = parse(Float64, f[6])
            if m < t.margen
                t.margen = m
                t.g6 = f[7]
            end
        end
    end

    mkpath(final)
    open(joinpath(final, "resumen.csv"), "w") do io
        println(io, "conjetura,revisados,contraejemplos,dudosos,n_min,margen_min,g6_min")
        for nombre in orden
            t = tot[nombre]
            nmin = t.nmin == typemax(Int) ? "" : string(t.nmin)
            @printf(io, "%s,%d,%d,%d,%s,%.10g,\"%s\"\n", nombre, t.revisados, t.contra,
                    t.dudosos, nmin, t.margen, replace(t.g6, "\"" => "\"\""))
        end
    end
    open(joinpath(final, "resumen.md"), "w") do io
        println(io, "# Resumen de verificación (combinado de $partes porciones)\n")
        println(io, "| Conjetura | Revisados | Contraejemplos | Dudosos | n mínimo | Margen mínimo | Grafo |")
        println(io, "|---|---|---|---|---|---|---|")
        for nombre in orden
            t = tot[nombre]
            nmin = t.nmin == typemax(Int) ? "—" : string(t.nmin)
            @printf(io, "| %s | %d | %d | %d | %s | %.6g | `%s` |\n", nombre, t.revisados,
                    t.contra, t.dudosos, nmin, t.margen, t.g6)
        end
    end
    # contraejemplos: concatenar los .g6 de las porciones (hasta 1000 por conjetura)
    for nombre in orden
        tot[nombre].contra == 0 && continue
        lineas = String[]
        for d in dirs
            p = joinpath(d, nombre * ".g6")
            isfile(p) || continue
            for l in eachline(p)
                length(lineas) >= 1000 && break
                push!(lineas, l)
            end
        end
        open(joinpath(final, nombre * ".g6"), "w") do io
            foreach(l -> println(io, l), lineas)
        end
    end

    total = isempty(orden) ? 0 : tot[orden[1]].revisados
    println("Grafos revisados en total: ", total, "  (compararlo con el conteo de geng -u)")
    for nombre in orden
        t = tot[nombre]
        @printf("  %-14s contraejemplos=%d  dudosos=%d  margen_min=%.6g\n",
                nombre, t.contra, t.dudosos, t.margen)
    end
    println("Resultado combinado en ", final)
end

main(ARGS)
