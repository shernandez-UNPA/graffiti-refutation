# GraffitiDE: evolución diferencial para refutar conjeturas de Graffiti

Código del experimento de la sección 6 del manuscrito para *Theory and Applications of Graphs* (TAG). Compara cuatro algoritmos en las 35 conjeturas de Graffiti de la Tabla 1 de Roucairol y Cazenave (2025):

| Código | Algoritmo |
| --- | --- |
| `rs` | Búsqueda aleatoria (línea base) |
| `de` | Evolución diferencial DE/rand/1/bin |
| `mde` | DE memética: DE + buscador local sobre el mejor individuo, activado a partir de la generación `ls_start` |
| `des` | DE con función objetivo suavizada, sin buscador local |
| `mdes` | DE memética con función objetivo suavizada (los conteos de autovalores se reemplazan por conteos con sigmoide para orientar la búsqueda; el éxito se decide siempre con la violación exacta) |

Solo usa la biblioteca estándar de Julia (LinearAlgebra, Random, Statistics, Printf, Dates). No hay que instalar paquetes.

## Requisitos

Julia 1.10 o posterior. Probado con Julia 1.13.

## Estructura

```
src/invariants.jl    invariantes (espectros, distancias, Randić, armónico, gravedad, ...)
src/conjectures.jl   las 35 conjeturas como E(G) = f(G) - g(G)
src/decoders.jl      codificación por claves aleatorias y decodificadores por clase
src/search.jl        rs, de, mde, mdes con currículum por orden
test/runtests.jl     pruebas de invariantes y decodificadores
scripts/experiment.jl  experimento completo (paralelo, reanudable)
scripts/summary.jl     tabla resumen (CSV y LaTeX)
```

## Cómo correrlo (Windows)

Abrir una terminal en la carpeta del proyecto:

```
cd "C:\Users\sebas\OneDrive\Documentos\Publicaciones\2026 - TAG\GraffitiDE"
julia test\runtests.jl
```

Prueba corta (13 conjeturas ya refutadas, 1 corrida de 60 s cada una):

```
julia -t 8 scripts\experiment.jl --only refuted --algs mdes --runs 1 --time 60
```

Usar `-t 8` (un hilo por núcleo físico del Ryzen 7 5700U) y no 16: con hyperthreading cada corrida se hace más lenta y los tiempos dejan de ser comparables.

### Protocolo propuesto

**Etapa A, comparación con la literatura.** Las 13 conjeturas que Roucairol y Cazenave refutaron, 4 algoritmos, 30 corridas, 15 min como límite:

```
julia -t 8 scripts\experiment.jl --only refuted --runs 30 --time 900
```

Como peor caso son unas 58 h. En la práctica es mucho menos, porque `mdes` y `mde` resuelven casi todo en segundos; el tiempo lo consumen `rs` y `de` cuando fallan. Para una primera pasada, `--time 300` la deja en un tercio.

**Etapa B, conjeturas abiertas.** Solo `mdes`, 10 corridas de 15 min por instancia (32 instancias): unas 10 h.

```
julia -t 8 scripts\experiment.jl --only open --algs mdes --runs 10 --time 900
```

Cada corrida se agrega a `results\runs.csv` apenas termina. Se puede cortar con Ctrl+C y volver a lanzar el mismo comando: las corridas hechas se saltean.

Resumen:

```
julia scripts\summary.jl
```

Genera `results\summary.csv` y `results\summary.tex`. Los contraejemplos quedan en `results\counterexamples\`, en formato graph6 (se pueden cargar en House of Graphs o con nauty) y como lista de aristas.

**Etapa C, ablación.** DE con objetivo suavizado y sin buscador local (`des`) en las 4 conjeturas refutadas donde la suavización cambia algo (29, 30, 166, 189). Menos de 1 h.

```
julia -t 8 scripts\experiment.jl --ids 29,30,166,189 --algs des --runs 30 --time 300
```

Tests de Wilcoxon pareados con corrección de Holm:

```
julia scripts\wilcoxon.jl
```

Hay un `.bat` para cada etapa (`etapa_A.bat`, `etapa_B.bat`, `etapa_C.bat`): doble clic y listo.

## Decisiones de diseño

- **Codificación.** Un individuo de orden n es x ∈ [0,1]^(1+n(n−1)/2): un umbral t y una clave por par de vértices. Los pares se recorren por clave decreciente (estilo Kruskal); un par se agrega si su clave supera t o si une dos componentes, y siempre que respete la clase (árbol, sin triángulos, cintura ≥ 5). Todos los grafos decodificados son conexos y factibles por construcción.
- **Currículum.** Los pares están ordenados por su vértice mayor, así que el genoma de orden n es prefijo del de orden n+1. La población final de orden n, extendida con claves nuevas, es la inicial de orden n+1. Al llegar a `nmax` sin éxito se reinicia desde `nmin`.
- **Buscador local.** Sobre el mejor individuo, desde la generación `ls_start` de cada orden: elige claves al azar y las pasa al otro lado del umbral (agrega o quita una arista), y acepta si la violación no empeora.
- **Semillas.** La semilla depende de la conjetura, la clase y la corrida, no del algoritmo: la comparación entre algoritmos es pareada (sirve para Wilcoxon de rangos con signo).
- **Verificación.** Cada contraejemplo se reevalúa con 5 reetiquetados aleatorios de los vértices. La columna `verified` del CSV lo registra, y el resumen solo cuenta los verificados.

## Diferencias con el código de referencia

Las definiciones replican el código de Roucairol y Cazenave (github.com/RoucairolMilo/refutationGBR) para que la comparación sea justa, con estas excepciones:

1. **Umbral de cero.** En las conjeturas 166, 189, 20, 302 y 715 el código de referencia decide el signo de un autovalor con umbrales entre 0 y 1e-13. Con esos umbrales un autovalor nulo calculado como 1e-16 cuenta como positivo y aparecen falsos contraejemplos (pasaba con la 20 en un árbol de 5 vértices). Acá se usa 1e-9; los umbrales mayores (1e-8 a 1e-4) se mantienen.
2. **Cintura ≥ 5.** El código de referencia solo permite agregar una arista entre vértices a distancia ≥ 5, lo que en realidad fuerza cintura ≥ 6. Acá se exige distancia ≥ 4, que es exactamente cintura ≥ 5.

## Cuestiones a revisar contra la fuente original ("Written on the Wall")

- **"Distancia media".** En las conjeturas 39, 40, 143 y 154 se promedian las n² entradas de la matriz de distancias, incluida la diagonal, como en el código de referencia. La definición usual excluye la diagonal.
- **Conjeturas 262 y 698.** Con la formulación del código de referencia no pueden fallar (en la 698, el laplaciano es semidefinido positivo). Probablemente hay un error de transcripción.
- **Condiciones extra.** En las conjeturas 195 y 198, el código de referencia solo acepta un contraejemplo si hay al menos tantos pares a distancia par como impar. Acá se aplica esa hipótesis (función `even_ge_odd`). Sin ella, P4 refuta la 195, lo que confirma que es parte de la conjetura. Hay que verificarla en la fuente original. En la 197 el código de referencia calcula la condición pero no la aplica, y acá tampoco.
- **Gravity, harmonic y range.** Se usan las definiciones corregidas que proponen Roucairol y Cazenave (gravedad de Brewster, Dinneen y Faber; armónico estándar; rango como máximo menos mínimo).

## Resultados de la validación inicial (nube, 2 núcleos, 60 s por corrida)

**Conjeturas ya refutadas (13), algoritmo `mdes`, 3 corridas cada una:** 39 de 39 corridas con contraejemplo verificado.

| Graffiti | Clase | Mediana del tiempo (s) | n mínimo | Mejor tiempo publicado (s) |
| --- | --- | --- | --- | --- |
| 137 | any | 5,0 | 33 | 513 (contraejemplos desde 67 vértices) |
| 139 | any | 0,62 | 24 | 36 (desde 50) |
| 197 | k3free | 0,23 | 14 | < 1 (el contraejemplo publicado es C17) |
| 289 | girth5 | 0,62 | 19 | 6 |
| las otras 9 | | ≤ 0,36 | | < 1 |

Los tiempos publicados son en Rust sobre otra máquina; para el artículo conviene comparar también evaluaciones de la función objetivo (columna `evals`).

**Ablación de la suavización (Graffiti 29):** `mde` 0 de 3 corridas en 60 s; `mdes` 3 de 3, todas en menos de 0,4 s. Por búsqueda exhaustiva, el contraejemplo más chico tiene 7 vértices (840 grafos etiquetados entre unos 2 millones).

**Conjeturas abiertas (22), sonda de 60 s con `mdes`:** ninguna refutada. Observaciones:

- **195:** sin la hipótesis Σeven ≥ Σodd del código de referencia, P4 la refuta, lo que indica que la hipótesis forma parte de la conjetura. Con la hipótesis, la mejor violación hallada es 0 (igualdad).
- **714:** la mejor violación es −0,0101 con 100 vértices, y corresponde a K₁₀₀. Para Kₙ el lado izquierdo vale 1 y el derecho n/(n−1): la cota se ajusta cuando n crece pero no se viola.
- **21, 284, 712:** la mejor violación hallada es 0 (casos de igualdad).

## Resultados de la etapa A (laptop, 8 hilos, 30 corridas de hasta 300 s)

- `de`, `mde` y `mdes` refutaron las 13 conjeturas en las 30 corridas. `rs` falló en todas las corridas de 137, 139 y 301.
- `de` usó significativamente menos evaluaciones que `rs` en 30, 137, 139, 197, 289 y 301 (Wilcoxon con Holm); `rs` fue mejor en 29 y 166, cuyos contraejemplos tienen 7 y 5 vértices.
- El buscador local no ayudó: `de` usó significativamente menos evaluaciones que `mde` en 8 conjeturas y nunca fue peor.
- La suavización fue decisiva en la 29: `mdes` 30 de 30 corridas, `mde` 2 de 30.
- Contraejemplos más chicos que los publicados: 31 vértices en la 137 (publicado: desde 67), 24 en la 139 (desde 50), 14 en la 197 (C17).
