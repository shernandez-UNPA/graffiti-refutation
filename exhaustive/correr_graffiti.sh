#!/usr/bin/env bash
#
# Corre la verificación exhaustiva de las conjeturas de Graffiti para una clase,
# orden por orden, guardando un resultado por cada n en
#   resultados/graffiti/<clase>/n<NN>/   (resumen.md, resumen.csv, contraejemplos .g6)
#
# Es reanudable: si un orden ya tiene resumen.csv, lo saltea.
#
# Uso:
#   bash correr_graffiti.sh CLASE NMIN NMAX [opciones extra para graffiti.jl]
#
# Ejemplos:
#   bash correr_graffiti.sh any 3 10
#   bash correr_graffiti.sh tree 4 22
#   bash correr_graffiti.sh k3free 4 14
#   bash correr_graffiti.sh girth5 5 18
#   bash correr_graffiti.sh any 11 11 --estado O        # solo las abiertas, n = 11
#
# Hilos: por defecto todos (-t auto). Para fijarlos: HILOS=8 bash correr_graffiti.sh ...

set -euo pipefail
cd "$(dirname "$0")"

clase=${1:?falta la clase (any, tree, k3free, girth5)}
nmin=${2:?falta NMIN}
nmax=${3:?falta NMAX}
shift 3
hilos=${HILOS:-auto}

# Nombres de los generadores: en Ubuntu llevan el prefijo nauty-
GENG=$(command -v nauty-geng || command -v geng) || { echo "no encuentro geng (sudo apt install nauty)"; exit 1; }
GENTREEG=$(command -v nauty-gentreeg || command -v gentreeg) || true

for n in $(seq "$nmin" "$nmax"); do
    dir="resultados/graffiti/$clase/n$(printf %02d "$n")"
    if [[ -f "$dir/resumen.csv" ]]; then
        echo "== $clase n=$n: ya hecho, se saltea"
        continue
    fi
    case "$clase" in
        any)    gen=("$GENG" -cq "$n") ;;
        tree)   gen=("$GENTREEG" -q "$n") ;;
        k3free) gen=("$GENG" -ctq "$n") ;;
        girth5) gen=("$GENG" -ctfq "$n") ;;
        *) echo "clase desconocida: $clase"; exit 1 ;;
    esac
    echo "== $clase n=$n  ($(date '+%F %T'))"
    tmp="$dir.parcial"
    rm -rf "$tmp"
    "${gen[@]}" | julia -t "$hilos" graffiti.jl --clase "$clase" --salida "$tmp" "$@"
    mv "$tmp" "$dir"
done

echo "Listo. Tabla resumen:  julia tabla_graffiti.jl"
