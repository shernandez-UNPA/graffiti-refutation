# Verificación exhaustiva de conjeturas en grafos — código

Pipeline en Julia que lee grafos en formato graph6 (salida de `nauty-geng`),
calcula invariantes y evalúa conjeturas sobre **todos** los grafos de un orden dado,
en paralelo con hilos. Solo usa la biblioteca estándar de Julia: no hace falta
instalar paquetes.

## Archivos

| Archivo | Contenido |
|---|---|
| `VerifConj.jl` | Módulo: lector/escritor graph6, invariantes, motor de verificación paralelo. |
| `conjeturas.jl` | Lista de conjeturas a verificar (es lo que se edita). |
| `verificar.jl` | Programa principal (línea de comandos). |
| `test.jl` | Pruebas de los invariantes con grafos de valores conocidos (K_n, P_n, C_n, estrellas, Petersen). |

## Requisitos (WSL / Ubuntu)

```bash
sudo apt install nauty                         # en Ubuntu: nauty-geng, nauty-showg, ...
curl -fsSL https://install.julialang.org | sh  # Julia (juliaup)
```

## Uso

Desde WSL, con la carpeta de Windows montada en `/mnt/c/...`:

```bash
cd "/mnt/c/Users/sebas/OneDrive/Documentos/Publicaciones/2027 - Grafos/Verificacion exhaustiva Julia/codigo"

julia test.jl                                              # 1) pruebas

nauty-geng -cq 9 | julia -t auto verificar.jl               # 2) todos los conexos de 9 vértices

# 3) todos los conexos de 3 a 10 vértices, guardando resultados
for n in 3 4 5 6 7 8 9 10; do nauty-geng -cq $n; done \
  | julia -t auto verificar.jl --salida resultados_n3a10
```

`-c` = solo conexos, `-q` = sin mensajes de geng. `-t auto` usa todos los núcleos.

Opciones de `verificar.jl`: `--conjeturas ARCHIVO`, `--salida DIR`, `--eps X`
(tolerancia, por defecto `1e-9`), `--entrada ARCHIVO` (graph6 desde archivo), `--silencioso`.

### Dividir el trabajo (n = 11 y más)

`geng` admite `res/mod` para generar solo una porción del espacio, sin solapamientos.
Cada porción se puede correr por separado (en otra terminal, otra máquina o otro día):

```bash
for r in 0 1 2 3 4 5 6 7 8 9; do
  nauty-geng -cq 11 $r/10 | julia -t auto verificar.jl --silencioso --salida res_n11_parte$r
done
```

## Cómo agregar una conjetura

En `conjeturas.jl`, cada conjetura afirma `margen(g) ≥ 0`. Por ejemplo, para
"λ₁ ≥ 2m/n":

```julia
Conjetura("grado_medio", "λ₁ ≥ 2m/n", g -> λ1(g) - 2 * n_aristas(g) / n_vertices(g))
```

Invariantes disponibles: `n_vertices`, `n_aristas`, `grados`, `grado_max`, `grado_min`,
`λ1`, `espectro_adyacencia`, `μ1_laplaciano`, `espectro_laplaciano`, `q1_sin_signo`,
`espectro_laplaciano_sin_signo`, `∂1_distancias`, `espectro_distancias`, `diametro`,
`radio`, `conexo`, `numero_matching`, `numero_independencia` (estos dos últimos son
**exactos**, por programación dinámica sobre subconjuntos: viables para n ≤ 12).

Los invariantes se calculan una sola vez por grafo aunque los usen varias conjeturas.

## Salida

- **contraejemplos**: grafos con margen < −ε.
- **dudosos**: grafos con |margen| ≤ ε, es decir, casos de igualdad numérica que conviene
  re-verificar con más precisión.
- **n mínimo**: menor orden en que apareció un contraejemplo (el contraejemplo mínimo).
- **margen mínimo** y el grafo (graph6) que lo alcanza: indica qué tan ajustada es la cota.
- Con `--salida DIR`: `resumen.md` y un `.g6` por conjetura con los primeros 1000 contraejemplos.
  Para verlos: `nauty-showg -e archivo.g6` (lista de aristas).

## Validación realizada (conexos, 3 ≤ n ≤ 9: 273 191 grafos)

| Conjetura | Contraejemplos | Dudosos | Comentario |
|---|---|---|---|
| `wagner_2_1` (λ₁ + μ ≥ √(n−1) + 1) | 0 | 7 | Los 7 casos de igualdad son las estrellas K₁,ₙ₋₁ (una por cada n). Coherente con que el contraejemplo conocido tiene 19 vértices. |
| `ctrl_grado_medio` (λ₁ ≥ 2m/n) | 0 | 53 | 53 = cantidad de grafos conexos regulares con 3 a 9 vértices (OEIS A005177). |
| `ctrl_laplaciano` (μ₁(L) ≤ n) | 0 | 15 072 | Igualdad cuando el complemento es no conexo. |
| `ctrl_falsa` (λ₁ ≥ Δ) | 273 138 | 53 | Todos los no regulares la violan; n mínimo = 3 (P₃). 273 138 + 53 = 273 191. |

Velocidad de referencia: ~80 000 grafos/s con 2 hilos y estas 4 conjeturas
(n = 10, en una máquina de prueba). Con más núcleos escala casi linealmente porque
los lotes son independientes.

---

## Conjeturas de Graffiti (Tabla 1 de Roucairol y Cazenave, 2025)

`graffiti.jl` verifica las 35 conjeturas con **las mismas definiciones** que el
experimento de evolución diferencial del artículo para TAG: la carpeta `graffiti/`
es una copia sin cambios de `2026 - TAG/GraffitiDE/src/invariants.jl` y
`conjectures.jl`. Un grafo es contraejemplo si E(G) = f(G) − g(G) > 10⁻⁴, igual que en GraffitiDE.

La clase la fija el generador:

| `--clase` | Generador | Grafos |
|---|---|---|
| `any` | `nauty-geng -c n` | conexos |
| `tree` | `nauty-gentreeg n` (sale en sparse6; el lector lo acepta) | árboles |
| `k3free` | `nauty-geng -ct n` | conexos sin triángulos |
| `girth5` | `nauty-geng -ctf n` | conexos de cintura ≥ 5 |

### Corrida completa (reanudable)

```bash
bash correr_graffiti.sh any 3 10
bash correr_graffiti.sh tree 4 22
bash correr_graffiti.sh k3free 4 14
bash correr_graffiti.sh girth5 5 17
julia tabla_graffiti.jl --de "/mnt/c/Users/sebas/OneDrive/Documentos/Publicaciones/2026 - TAG/GraffitiDE/results/summary.csv"
```

### Órdenes grandes, en porciones (pausables)

Para un orden que tarda horas (`k3free 14`, `girth5 18`), `correr_por_partes.sh`
lo parte en porciones independientes con el `res/mod` de geng. Al cortar con Ctrl+C
se pierde como mucho la porción en curso; el mismo comando retoma desde la siguiente.

```bash
bash correr_por_partes.sh k3free 14 10    # 10 porciones de ~20 min
bash correr_por_partes.sh girth5 18 50    # 50 porciones de ~20 min
```

Hay que usar siempre el mismo número de porciones para un mismo orden (el script lo
controla). Al terminar la última, `combinar_partes.jl` junta todo en
`resultados/graffiti/<clase>/nNN/`, con el mismo formato que una corrida entera, así que
`tabla_graffiti.jl` lo toma igual y `correr_graffiti.sh` lo saltea. Verificado con
`k3free 12` en 4 porciones (interrumpido y retomado): resultado idéntico a la corrida entera.

Cada orden queda en `resultados/graffiti/<clase>/nNN/`. Si se corta (Ctrl+C), al
volver a lanzar el mismo comando se saltean los órdenes terminados.
`tabla_graffiti.jl` arma `tabla.md`, `tabla.csv` y `tabla.tex` con el orden mínimo de
contraejemplo, cuántos hay de ese orden y la comparación con la DE.

Tamaño de cada clase (grafos a procesar):

| n | conexos | sin triángulos | cintura ≥ 5 |
|---|---|---|---|
| 10 | 11 716 571 | 9 832 | — |
| 12 | — | 1 144 061 | — |
| 13 | — | 19 425 052 | — |
| 14 | — | 445 781 050 | 275 480 |
| 16 | — | — | 17 772 647 |
| 17 | — | — | 179 593 823 |
| 18 | — | — | ~1 800 000 000 (estimado, ~10× n = 17) |

(conteos obtenidos con `geng -u`.) Los árboles son pocos: hasta n = 22 son unos 5,6 millones en total.

### Resultados preliminares (máquina de prueba, 2 hilos)

Rangos verificados: conexos 3–9, árboles 4–18, sin triángulos 4–12, cintura ≥ 5 5–15.

| Conj. | Clase | n mínimo (exhaustivo) | # de ese orden | n mínimo DE |
|---|---|---|---|---|
| 189 | any | **3** | 1 (K₃) | 4 |
| 166 | any | 5 | 1 | 5 |
| 711 | any | 5 | 2 | 5 |
| 29 | any | 7 | 1 | 7 |
| 30 | any | 9 | 1 | 9 |
| 321 | k3free | 4 | 3 | 4 |
| 715 | tree | 12 | 2 | 12 |
| 301 | tree | 14 | 1 | 14 |
| 302 | tree | 16 | 2 | 16 |
| 137, 139 | any | > 9 | — | 31, 24 |
| 197 | k3free | > 12 | — | 14 |
| 289 | girth5 | > 15 | — | 19 |

Ninguna conjetura abierta tiene contraejemplos en esos rangos.
