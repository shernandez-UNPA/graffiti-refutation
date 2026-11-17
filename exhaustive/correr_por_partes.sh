#!/usr/bin/env bash
#
# Corre UN orden grande partido en porciones independientes (opción res/mod de
# geng/gentreeg), para poder pausar y retomar perdiendo como mucho una porción.
#
# Cada porción queda en  resultados/graffiti/<clase>/partes_nNN/pRRdeMM/
# Cuando están todas, se combinan en  resultados/graffiti/<clase>/nNN/
# (mismo formato que correr_graffiti.sh, así que tabla_graffiti.jl lo toma igual
# y correr_graffiti.sh lo saltea).
#
# Uso:
#   bash correr_por_partes.sh CLASE N PARTES [opciones extra para graffiti.jl]
#
# Ejemplos:
#   bash correr_por_partes.sh k3free 14 10     # ~446 M grafos en 10 porciones
#   bash correr_por_partes.sh girth5 18 50     # ~1 800 M grafos en 50 porciones
#
# Pausar: Ctrl+C. Retomar: el mismo comando (saltea las porciones terminadas).
# Importante: usar siempre el mismo número de PARTES para un mismo orden.
#
# Hilos: por defecto todos (-t auto). Para fijarlos: HILOS=8 bash correr_por_partes.sh ...

set -euo pipefail
cd "$(dirname "$0")"

clase=${1:?falta la clase (any, tree, k3free, girth5)}
n=${2:?falta N}
partes=${3:?falta PARTES}
shift 3
hilos=${HILOS:-auto}

GENG=$(command -v nauty-geng || command -v geng) || { echo "no encuentro geng (sudo apt install nauty)"; exit 1; }
GENTREEG=$(command -v nauty-gentreeg || command -v gentreeg) || true

nn=$(printf %02d "$n")
final="resultados/graffiti/$clase/n$nn"
base="resultados/graffiti/$clase/partes_n$nn"

if [[ -f "$final/resumen.csv" ]]; then
    echo "== $clase n=$n: ya está completo en $final"
    exit 0
fi
mkdir -p "$base"
# Evita mezclar porciones de particiones distintas
if [[ -f "$base/PARTES" && "$(cat "$base/PARTES")" != "$partes" ]]; then
    echo "Este orden se empezó con $(cat "$base/PARTES") partes; usá ese mismo número."
    exit 1
fi
echo "$partes" > "$base/PARTES"

inicio=$(date +%s)
for r in $(seq 0 $((partes - 1))); do
    dir="$base/p$(printf %02d "$r")de$(printf %02d "$partes")"
    if [[ -f "$dir/resumen.csv" ]]; then
        continue
    fi
    case "$clase" in
        any)    gen=("$GENG" -cq "$n" "$r/$partes") ;;
        tree)   gen=("$GENTREEG" -q "$n" "$r/$partes") ;;
        k3free) gen=("$GENG" -ctq "$n" "$r/$partes") ;;
        girth5) gen=("$GENG" -ctfq "$n" "$r/$partes") ;;
        *) echo "clase desconocida: $clase"; exit 1 ;;
    esac
    hechas=$(find "$base" -mindepth 2 -maxdepth 2 -name resumen.csv | wc -l)
    echo "== $clase n=$n  porción $((r + 1))/$partes  (terminadas: $hechas)  $(date '+%F %T')"
    tmp="$dir.parcial"
    rm -rf "$tmp"
    "${gen[@]}" | julia -t "$hilos" graffiti.jl --clase "$clase" --salida "$tmp" --silencioso "$@"
    mv "$tmp" "$dir"
done

echo "== Combinando $partes porciones en $final"
julia combinar_partes.jl "$base" "$final"
echo "Listo ($(( ($(date +%s) - inicio) / 60 )) min en esta sesión)."
