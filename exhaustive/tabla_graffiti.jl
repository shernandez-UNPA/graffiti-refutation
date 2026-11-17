#!/usr/bin/env julia
#
# Junta los resultados de resultados/graffiti/<clase>/nNN/resumen.csv y arma la
# tabla por conjetura: hasta qué orden se verificó, orden mínimo de un
# contraejemplo y cuántos contraejemplos hay de ese orden.
#
# Opcionalmente compara con el orden mínimo hallado por la evolución diferencial
# (summary.csv de GraffitiDE, columna min_n; se toma el mínimo entre algoritmos).
#
# Uso:
#   julia tabla_graffiti.jl
#   julia tabla_graffiti.jl --de "/mnt/c/Users/sebas/OneDrive/Documentos/Publicaciones/2026 - TAG/GraffitiDE/results/summary.csv"
#
# Escribe resultados/graffiti/tabla.md, tabla.csv y tabla.tex (fragmento tabular).

include(joinpath(@__DIR__, "..", "de", "src", "invariants.jl"))
include(joinpath(@__DIR__, "..", "de", "src", "conjectures.jl"))
using Printf

const RAIZ = joinpath(@__DIR__, "resultados", "graffiti")

# Parser mínimo de CSV con campos opcionalmente entre comillas dobles.
function campos_csv(linea)
    out = String[]; buf = IOBuffer(); comillas = false; i = 1
    cs = collect(linea)
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

function leer_de(path)
    d = Dict{Tuple{Int,String},Int}()
    lineas = readlines(path)
    cab = campos_csv(lineas[1])
    ic, icl, imn = findfirst(==("conj"), cab), findfirst(==("class"), cab), findfirst(==("min_n"), cab)
    for l in lineas[2:end]
        f = campos_csv(l)
        isempty(f[imn]) && continue
        k = (parse(Int, f[ic]), f[icl])
        d[k] = min(get(d, k, typemax(Int)), parse(Int, f[imn]))
    end
    return d
end

struct Fila
    id::Int
    clase::String
    estado::Char
    enunciado::String
    nverif::Int            # mayor n verificado (todos los n entre el menor y este)
    ndesde::Int
    nmin::Int              # orden mínimo de contraejemplo (0 si no hay)
    cant_nmin::Int         # contraejemplos de ese orden
    g6::String             # un contraejemplo de orden mínimo (el de margen mínimo en ese n)
    margen::Float64        # margen mínimo observado en todo lo verificado
    revisados::Int
end

function main(args)
    de = Dict{Tuple{Int,String},Int}()
    if length(args) >= 2 && args[1] == "--de"
        de = leer_de(args[2])
    end
    isdir(RAIZ) || error("no hay resultados en $RAIZ (correr antes correr_graffiti.sh)")
    filas = Fila[]
    for clase in sort(readdir(RAIZ))
        dclase = joinpath(RAIZ, clase)
        isdir(dclase) || continue
        ns = sort([parse(Int, d[2:end]) for d in readdir(dclase)
                   if startswith(d, "n") && isfile(joinpath(dclase, d, "resumen.csv"))])
        isempty(ns) && continue
        # tramo contiguo desde el menor n
        hasta = ns[1]
        for k in 2:length(ns)
            ns[k] == hasta + 1 ? (hasta = ns[k]) : break
        end
        por_conj = Dict{String,Vector{Tuple{Int,Vector{String}}}}()
        for n in ns
            n > hasta && break
            ls = readlines(joinpath(dclase, @sprintf("n%02d", n), "resumen.csv"))
            for l in ls[2:end]
                f = campos_csv(l)
                push!(get!(por_conj, f[1], Tuple{Int,Vector{String}}[]), (n, f))
            end
        end
        for (nombre, regs) in por_conj
            id = parse(Int, match(r"^G(\d+)_", nombre).captures[1])
            c = CONJECTURES[id]
            nmin, cant, g6 = 0, 0, ""
            margen = Inf; rev = 0
            for (n, f) in sort(regs; by = first)
                rev += parse(Int, f[2])
                margen = min(margen, parse(Float64, f[6]))
                if nmin == 0 && parse(Int, f[3]) > 0
                    nmin, cant, g6 = n, parse(Int, f[3]), f[7]
                end
            end
            push!(filas, Fila(id, clase, c.status, c.statement, hasta, ns[1], nmin, cant, g6, margen, rev))
        end
    end
    sort!(filas; by = f -> (f.estado == 'R' ? 0 : 1, f.clase, f.id))

    ncol(f) = f.nmin == 0 ? "—" : string(f.nmin)
    decol(f) = (v = get(de, (f.id, f.clase), 0)) == 0 ? "—" : string(v)

    open(joinpath(RAIZ, "tabla.md"), "w") do io
        println(io, "| Conj. | Clase | Estado | Enunciado | Verificado n = | n mínimo | # de ese orden | n mínimo DE | Contraejemplo mínimo | Margen mínimo |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        for f in filas
            @printf(io, "| %d | %s | %s | %s | %d–%d | %s | %s | %s | `%s` | %.4g |\n",
                    f.id, f.clase, f.estado, f.enunciado, f.ndesde, f.nverif, ncol(f),
                    f.nmin == 0 ? "—" : string(f.cant_nmin), decol(f), f.g6, f.margen)
        end
    end
    open(joinpath(RAIZ, "tabla.csv"), "w") do io
        println(io, "conj,clase,estado,n_desde,n_hasta,n_min,cant_n_min,n_min_DE,g6_min,margen_min,revisados")
        for f in filas
            @printf(io, "%d,%s,%s,%d,%d,%s,%s,%s,\"%s\",%.10g,%d\n", f.id, f.clase, f.estado,
                    f.ndesde, f.nverif, f.nmin == 0 ? "" : string(f.nmin),
                    f.nmin == 0 ? "" : string(f.cant_nmin), get(de, (f.id, f.clase), ""),
                    replace(f.g6, "\"" => "\"\""), f.margen, f.revisados)
        end
    end
    open(joinpath(RAIZ, "tabla.tex"), "w") do io
        println(io, "\\begin{tabular}{rllrrrr}\n\\hline")
        println(io, "Conj. & Class & Status & Verified up to \$n\$ & Min.\\ \$n\$ & \\# at min.\\ \$n\$ & Min.\\ \$n\$ (DE) \\\\\n\\hline")
        for f in filas
            @printf(io, "%d & %s & %s & %d & %s & %s & %s \\\\\n", f.id, f.clase, f.estado, f.nverif,
                    ncol(f), f.nmin == 0 ? "--" : string(f.cant_nmin), replace(decol(f), "—" => "--"))
        end
        println(io, "\\hline\n\\end{tabular}")
    end
    for f in filas
        @printf("%5d %-7s %s  n=%2d–%-2d  n_min=%-3s  #=%-6s  DE=%-3s  margen_min=%.4g\n",
                f.id, f.clase, f.estado, f.ndesde, f.nverif, ncol(f),
                f.nmin == 0 ? "—" : string(f.cant_nmin), decol(f), f.margen)
    end
    println("\nTablas escritas en ", RAIZ, " (tabla.md, tabla.csv, tabla.tex)")
end

main(ARGS)
