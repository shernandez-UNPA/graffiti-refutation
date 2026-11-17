# Las 35 conjeturas de Graffiti de la Tabla 1 de Roucairol y Cazenave (2025),
# "Refutation of Spectral Graph Theory Conjectures with Search Algorithms".
#
# Cada conjetura es una desigualdad f(G) <= g(G). La función `violation`
# devuelve E(G) = f(G) - g(G); un grafo con E(G) > TOL es un contraejemplo.
# Las fórmulas replican el código de referencia de los autores
# (github.com/RoucairolMilo/refutationGBR), incluidas sus tolerancias,
# para que la comparación sea justa. Las particularidades se anotan en `nota`.

const TOL = 1e-4

# Umbral para decidir si un autovalor calculado es cero. El código de
# referencia usa umbrales de 0 a 1e-13 en algunas conjeturas (166, 189, 20,
# 302, 715); con ellos un autovalor nulo calculado como 1e-16 cuenta como
# positivo y aparecen falsos contraejemplos. Usamos 1e-9 en esos casos y
# mantenemos los umbrales mayores (1e-8 a 1e-4) del código de referencia.
const ZT = 1e-9

struct Conjecture
    id::Int
    classes::Vector{Symbol}   # :any, :tree, :k3free, :girth5
    status::Char              # 'R' refutada antes, 'O' abierta
    nmax::Int                 # orden máximo explorado por los autores
    ref_time::Float64         # mejor tiempo (s) entre sus 8 algoritmos; NaN si no la refutaron
    statement::String
    nota::String
    f::Function               # A -> E(G)
end

# Abreviaturas
λA(A) = spec(A)
λL(A) = spec(laplacian(A))
λD(D) = spec(D)
# Conteos de autovalores. Con s = true se usa un conteo suave (sigmoide de
# ancho τ) que sirve de desempate cuando un lado de la desigualdad es entero.
const τ = 0.05
σ(z) = 1 / (1 + exp(-z))
npos(λ, tol, s = false)    = s ? sum(σ.((λ .- tol) ./ τ)) : count(>(tol), λ)
nneg(λ, tol, s = false)    = s ? sum(σ.((-λ .- tol) ./ τ)) : count(<(-tol), λ)
nnonpos(λ, tol, s = false) = s ? sum(σ.((tol .- λ) ./ τ)) : count(<=(tol), λ)

function c166(A, s = false)  # sqrt(m) <= número de autovalores no positivos
    λ = λA(A)
    return sqrt(nedges(A)) - nnonpos(λ, ZT, s)
end
function c321(A, s = false)  # índice <= promedio del vector "even"
    D = dist_matrix(A)
    return λA(A)[end] - mean(even_vec(D))
end
function c189(A, s = false)  # moda de autovalores laplacianos <= cantidad de autovalores no positivos
    return mode_tol(λL(A)) - nnonpos(λA(A), ZT, s)
end
function c289(A, s = false)  # λ2 <= promedio de grados duales
    return second_largest(λA(A)) - mean(dual_degrees(A))
end
function c302(A, s = false)  # amplitud de autovalores positivos (λmax - mínimo positivo) <= promedio de grados duales
    λ = λA(A)
    pos = filter(>(ZT), λ)
    sp = isempty(pos) ? 0.0 : λ[end] - minimum(pos)
    return sp - mean(dual_degrees(A))
end
function c715(A, s = false)  # amplitud de autovalores no positivos <= promedio de grados mayores que el grado medio
    λ = λA(A)
    np = filter(<=(ZT), λ)
    sp = isempty(np) ? 0.0 : maximum(np) - minimum(np)
    d = degrees(A); m̄ = mean(d)
    hi = d[d .> m̄]
    isempty(hi) && return -Inf      # grafo regular: el código de referencia da NaN
    return sp - mean(hi)
end
function c29(A, s = false)   # índice de Randić <= autovalores negativos de la matriz de distancias
    return randic(A) - nneg(λD(dist_matrix(A)), 1e-4, s)
end
function c30(A, s = false)   # autovalores positivos de la matriz de distancias <= suma de temperaturas
    return npos(λD(dist_matrix(A)), 1e-4, s) - sum(temperature(A))
end
function c301(A, s = false)  # amplitud de autovalores positivos <= índice armónico
    λ = λA(A)
    pos = filter(>(1e-7), λ)
    sp = isempty(pos) ? 0.0 : maximum(pos) - minimum(pos)
    return sp - harmonic(A)
end
function c137(A, s = false)  # λ2 <= índice armónico
    return second_largest(λA(A)) - harmonic(A)
end
function c139(A, s = false)  # -(segundo menor autovalor) <= índice armónico
    return -second_smallest(λA(A)) - harmonic(A)
end
function c711(A, s = false)  # rango de la deficiencia <= rango del espectro
    λ = λA(A)
    def = [deficiency(A, v) for v in 1:nv(A)]
    return (maximum(def) - minimum(def)) - (λ[end] - λ[1])
end
function c197(A, s = false)  # -(segundo menor autovalor) <= rango del espectro de la matriz de gravedad
    μ = spec(gravity_matrix(A))
    return -second_smallest(λA(A)) - (μ[end] - μ[1])
end
function c140(A, s = false)  # desviación estándar del espectro <= índice armónico
    return popstd(λA(A)) - harmonic(A)
end
function c290(A, s = false)  # -(segundo menor autovalor) <= m / promedio de la matriz de gravedad
    return -second_smallest(λA(A)) - nedges(A) / mean(gravity_matrix(A))
end
function c284(A, s = false)  # mínimo grado dual <= -(menor autovalor de distancias)
    nv(A) < 3 && return -Inf
    return minimum(dual_degrees(A)) + λD(dist_matrix(A))[1]
end
"""
Hipótesis que el código de referencia exige antes de aceptar un contraejemplo
de las conjeturas 195 y 198: hay al menos tantos pares a distancia par como a
distancia impar. Si no se cumple, la violación se fuerza a un valor negativo
que mide cuánto falta para cumplirla (así la búsqueda conserva una dirección).
"""
function even_ge_odd(E, D)
    se, so = sum(even_vec(D)), sum(odd_vec(D))
    se >= so && return E
    return min(E, 0.0) - (so - se) / size(D, 1)^2
end

function c195(A, s = false)  # índice <= máximo del vector "even" (si pares pares >= pares impares)
    D = dist_matrix(A)
    return even_ge_odd(λA(A)[end] - maximum(even_vec(D)), D)
end
function c21(A, s = false)   # autovalores negativos <= suma de autovalores positivos
    nv(A) < 3 && return -Inf
    λ = λA(A)
    return nneg(λ, 1e-4, s) - sum(filter(>(1e-4), λ))
end
function c39(A, s = false)   # desviación estándar de las distancias <= autovalores positivos
    D = dist_matrix(A)
    return popstd(vec(D)) - npos(λA(A), 1e-4, s)
end
function c143(A, s = false)  # varianza de autovalores positivos <= m / distancia media
    D = dist_matrix(A)
    pos = filter(>(1e-4), λA(A))
    return popstd(pos)^2 - nedges(A) / mean(D)
end
function c154(A, s = false)  # desviación estándar del espectro <= n / distancia media
    D = dist_matrix(A)
    return popstd(λA(A)) - nv(A) / mean(D)
end
function c20(A, s = false)   # cantidad de autovalores positivos <= su suma
    pos = filter(>(ZT), λA(A))
    return (s ? npos(λA(A), ZT, true) : length(pos)) - sum(pos; init = 0.0)
end
function c40(A, s = false)   # desviación estándar de las distancias <= autovalores negativos
    D = dist_matrix(A)
    return popstd(vec(D)) - nneg(λA(A), 1e-4, s)
end
function c254(A, s = false)  # mínima brecha del espectro laplaciano <= suma de 1/(vértices a distancia impar)
    D = dist_matrix(A)
    gap = minimum(vec_derivative(λL(A)))
    return gap - sum(1 ./ odd_vec(D))
end
function c262(A, s = false)  # -índice <= máximo del vector "even"
    D = dist_matrix(A)
    return -λA(A)[end] - maximum(even_vec(D))
end
function c712(A, s = false)  # mínima temperatura <= autovalores no positivos
    return minimum(temperature(A)) - nnonpos(λA(A), 1e-5, s)
end
function c198(A, s = false)  # menor autovalor <= n / promedio de gravedad (si pares pares >= pares impares)
    D = dist_matrix(A)
    return even_ge_odd(λA(A)[1] - nv(A) / mean(gravity_matrix(A, D)), D)
end
function c714(A, s = false)  # -(promedio de autovalores no positivos) <= suma de 1/temperatura
    np = filter(<=(1e-5), λA(A))
    return -mean(np) - sum(1 ./ temperature(A))
end
function c219(A, s = false)  # segundo mayor autovalor de la matriz de gravedad <= n(n-1)/2 - m
    n = nv(A)
    return second_largest(spec(gravity_matrix(A))) - (n * (n - 1) / 2 - nedges(A))
end
function c322(A, s = false)  # suma de 1/even (entradas no nulas) <= rango del espectro de distancias
    D = dist_matrix(A)
    ev = even_vec(D)
    μ = λD(D)
    return sum(1 ./ ev[ev .!= 0]; init = 0.0) - (μ[end] - μ[1])
end
function c292(A, s = false)  # menor autovalor positivo <= n / promedio de la matriz de gravedad
    pos = filter(>(1e-8), λA(A))
    return minimum(pos) - nv(A) / mean(gravity_matrix(A))
end
function c295(A, s = false)  # autovalores positivos de distancias <= n / promedio de la matriz de gravedad
    D = dist_matrix(A)
    return npos(λD(D), 1e-8, s) - nv(A) / mean(gravity_matrix(A, D))
end
function c129(A, s = false)  # desviación estándar del espectro laplaciano <= índice de Randić
    return popstd(λL(A)) - randic(A)
end
function c698(A, s = false)  # norma de los autovalores laplacianos negativos <= índice de Randić
    neg = filter(<(-1e-8), λL(A))
    return sqrt(sum(abs2, neg; init = 0.0)) - randic(A)
end
function c252(A, s = false)  # mínima brecha del espectro laplaciano <= suma de 1/grado dual
    return minimum(vec_derivative(λL(A))) - sum(1 ./ dual_degrees(A))
end

const AT = [:any, :tree]
const NOTA_MEDIA = "La distancia media promedia las n² entradas de la matriz de distancias, incluida la diagonal (como el código de referencia)."
const NOTA_TRIVIAL = "Con esta formulación la desigualdad no puede fallar; probablemente hay un error de transcripción en la fuente."

const CONJECTURES = Dict{Int,Conjecture}(c.id => c for c in [
    Conjecture(166, [:any],    'R', 50, 0.0,   "sqrt(m) <= #autovalores no positivos", "", c166),
    Conjecture(321, [:k3free], 'R', 50, 0.0,   "índice <= promedio de even", "", c321),
    Conjecture(189, [:any],    'R', 50, 0.0,   "moda(espectro laplaciano) <= #autovalores no positivos", "Empates en la moda se resuelven como en el código de referencia.", c189),
    Conjecture(289, [:girth5], 'R', 50, 6.0,   "λ2 <= promedio de grados duales", "", c289),
    Conjecture(302, [:tree],   'R', 50, 0.0,   "λmax - mínimo autovalor positivo <= promedio de grados duales", "", c302),
    Conjecture(715, [:tree],   'R', 50, 0.0,   "amplitud de autovalores no positivos <= promedio de grados mayores que la media", "", c715),
    Conjecture(29,  [:any],    'R', 50, 0.0,   "Randić <= #autovalores negativos de distancias", "", c29),
    Conjecture(30,  [:any],    'R', 50, 0.0,   "#autovalores positivos de distancias <= suma de temperaturas", "", c30),
    Conjecture(301, [:tree],   'R', 50, 0.0,   "amplitud de autovalores positivos <= armónico", "", c301),
    Conjecture(137, [:any],    'R', 67, 513.0, "λ2 <= armónico", "", c137),
    Conjecture(139, [:any],    'R', 50, 36.0,  "-λ_{n-1} <= armónico", "", c139),
    Conjecture(711, [:any],    'R', 50, 0.0,   "rango de deficiencias <= rango del espectro", "", c711),
    Conjecture(197, [:k3free], 'O', 50, 0.0,   "-λ_{n-1} <= rango del espectro de gravedad", "Abierta en la literatura; Roucairol y Cazenave la refutan con C17.", c197),
    Conjecture(140, AT,        'R', 50, NaN,   "desviación estándar del espectro <= armónico", "Refutación atribuida a Favaron et al., no localizada por los autores.", c140),
    Conjecture(290, [:girth5], 'O', 50, NaN,   "-λ_{n-1} <= m / promedio de gravedad", "", c290),
    Conjecture(284, [:girth5], 'O', 50, NaN,   "mínimo grado dual <= -λmin(distancias)", "", c284),
    Conjecture(195, [:any],    'O', 50, NaN,   "índice <= máximo de even", "Se exige Σeven >= Σodd, como en el código de referencia; sin esa condición P4 la refuta.", c195),
    Conjecture(21,  AT,        'O', 50, NaN,   "#autovalores negativos <= suma de autovalores positivos", "", c21),
    Conjecture(39,  AT,        'O', 50, NaN,   "desviación estándar de distancias <= #autovalores positivos", NOTA_MEDIA, c39),
    Conjecture(143, AT,        'O', 100, NaN,  "varianza de autovalores positivos <= m / distancia media", NOTA_MEDIA, c143),
    Conjecture(154, AT,        'O', 50, NaN,   "desviación estándar del espectro <= n / distancia media", NOTA_MEDIA, c154),
    Conjecture(20,  AT,        'O', 50, NaN,   "#autovalores positivos <= suma de autovalores positivos", "", c20),
    Conjecture(40,  AT,        'O', 50, NaN,   "desviación estándar de distancias <= #autovalores negativos", NOTA_MEDIA, c40),
    Conjecture(254, AT,        'O', 50, NaN,   "mínima brecha laplaciana <= suma de 1/odd", "", c254),
    Conjecture(262, AT,        'O', 50, NaN,   "-índice <= máximo de even", NOTA_TRIVIAL, c262),
    Conjecture(712, AT,        'O', 50, NaN,   "mínima temperatura <= #autovalores no positivos", "", c712),
    Conjecture(198, [:any],    'O', 50, NaN,   "λmin <= n / promedio de gravedad", "Se exige Σeven >= Σodd, como en el código de referencia.", c198),
    Conjecture(714, [:any],    'O', 100, NaN,  "-promedio de autovalores no positivos <= suma de 1/temperatura", "", c714),
    Conjecture(219, [:k3free], 'O', 50, NaN,   "λ2(gravedad) <= n(n-1)/2 - m", "", c219),
    Conjecture(322, [:k3free], 'O', 50, NaN,   "suma de 1/even <= rango del espectro de distancias", "Con la definición de 'range' como cantidad de valores distintos, C4 la refuta.", c322),
    Conjecture(292, [:girth5], 'O', 50, NaN,   "menor autovalor positivo <= n / promedio de gravedad", "", c292),
    Conjecture(295, [:girth5], 'O', 50, NaN,   "#autovalores positivos de distancias <= n / promedio de gravedad", "", c295),
    Conjecture(129, [:any],    'O', 50, NaN,   "desviación estándar laplaciana <= Randić", "", c129),
    Conjecture(698, [:any],    'O', 50, NaN,   "norma de autovalores laplacianos negativos <= Randić", NOTA_TRIVIAL * " (el laplaciano es semidefinido positivo).", c698),
    Conjecture(252, [:any],    'O', 50, NaN,   "mínima brecha laplaciana <= suma de 1/grado dual", "", c252),
])

"Violación E(G) robusta: NaN o error numérico cuentan como -Inf. Con s = true, versión suavizada."
function violation(c::Conjecture, A, s::Bool = false)
    e = try
        c.f(A, s)
    catch
        -Inf
    end
    return isnan(e) ? -Inf : float(e)
end
