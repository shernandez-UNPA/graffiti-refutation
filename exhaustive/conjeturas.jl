# Lista de conjeturas a verificar.
#
# Cada conjetura afirma  margen(g) ≥ 0  para todo grafo g de la clase.
# El motor cuenta como contraejemplo todo grafo con margen < −ε, y como "dudoso"
# todo grafo con |margen| ≤ ε (casos de igualdad que conviene revisar con más precisión).
#
# La clase de grafos la decide la generación: por ejemplo `nauty-geng -c` solo
# produce grafos conexos, `-t` sin triángulos, `nauty-gentreeg` árboles, etc.
#
# Para agregar una conjetura: copiar una entrada y cambiar nombre, descripción y margen.
# Invariantes disponibles (ver VerifConj.jl):
#   n_vertices, n_aristas, grados, grado_max, grado_min,
#   λ1, espectro_adyacencia, μ1_laplaciano, espectro_laplaciano,
#   q1_sin_signo, espectro_laplaciano_sin_signo, ∂1_distancias, espectro_distancias,
#   diametro, radio, conexo, numero_matching, numero_independencia

CONJETURAS = Conjetura[

    # --- Conjetura real (benchmark de la literatura) -------------------------
    Conjetura("wagner_2_1",
        "λ₁ + μ ≥ √(n−1) + 1 (conexos, n ≥ 3). Refutada por Wagner (2021) con un grafo de 19 vértices; la estrella da igualdad.",
        g -> λ1(g) + numero_matching(g) - sqrt(n_vertices(g) - 1) - 1),

    # --- Controles: sirven para validar el pipeline --------------------------
    Conjetura("ctrl_grado_medio",
        "Teorema: λ₁ ≥ 2m/n. No debe haber contraejemplos; los grafos regulares dan igualdad.",
        g -> λ1(g) - 2 * n_aristas(g) / n_vertices(g)),

    Conjetura("ctrl_laplaciano",
        "Teorema: μ₁(L) ≤ n. No debe haber contraejemplos; igualdad si el complemento es no conexo.",
        g -> n_vertices(g) - μ1_laplaciano(g)),

    Conjetura("ctrl_falsa",
        "Afirmación falsa: λ₁ ≥ Δ. Debe dar contraejemplos desde n = 3 (el camino P₃).",
        g -> λ1(g) - grado_max(g)),
]
