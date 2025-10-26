# Código Sage: Wiedemann que llama a Padé vía EEA normalizado (Alg. 3.14)
from sage.all import *

# ---------- utilidades ----------
def rev_d(poly, d):
    """
    rev_d: reverso de grado d de un polinomio 'poly' en R = F[x].
    Si poly = sum_{i=0..deg} a_i x^i, entonces
    rev_d(poly) = a_0 x^d + a_1 x^{d-1} + ... + a_deg x^{d - deg}.
    Devuelve polinomio en el mismo anillo que 'poly'.
    """
    R = poly.parent()
    x = R.gen()
    coeffs = [poly.coefficient(i) for i in range(d+1)]  # si i>deg -> coef = 0
    res = R.zero()
    for i, a in enumerate(coeffs):
        res += a * x**(d - i)
    return res


from sage.all import *

def egcd_pade_canonico(R, pol_f, pol_g, k):
    """
    Extended Euclidean Algorithm (EEA) 'canónico' para Padé sobre R = Zmod(q)[x].
    Sigue el algoritmo 3.14, con normalización t(0)=1.
    """
    # R = PolynomialRing(Zmod(q), 'x'); x = R.gen()
    f = R(pol_f)
    g = R(pol_g)

    # Verificar que f y g no sean nulos
    if f.is_zero() or g.is_zero():
        raise ValueError("f y g deben ser no nulos.")

    # Normalización inicial (leading coefficient = 1)
    inv_lc_f = 1 / f.leading_coefficient()
    inv_lc_g = 1 / g.leading_coefficient()

    r = [f * inv_lc_f, g * inv_lc_g]
    s = [R(inv_lc_f), R(0)]
    t = [R(0),         R(inv_lc_g)]

    i = 1
    while True:
        if r[i].is_zero():
            j = i - 1
            rj, tj = r[j], t[j]

            # Si cumple grados y coprimalidad
            if rj.degree() < k and rj.gcd(tj) == 1:
                # Normalizar para que t(0) = 1 si es posible
                if not tj(0).is_zero():
                    inv_t0 = 1 / tj(0)
                    rj *= inv_t0
                    tj *= inv_t0

                tau = tj.leading_coefficient()
                return ("ok", rj / tau, tj / tau, j)

            return ("unsat", rj, tj, j)

        if r[i].degree() < k:
            j = i
            rj, tj = r[j], t[j]

            # Si cumple grados y coprimalidad
            if rj.gcd(tj) == 1:
                # Normalizar para que t(0) = 1 si es posible
                if not tj(0).is_zero():
                    inv_t0 = 1 / tj(0)
                    rj *= inv_t0
                    tj *= inv_t0

                tau = tj.leading_coefficient()
                return ("ok", rj / tau, tj / tau, j)

            return ("unsat", rj, tj, j)

        # División y actualización
        qi, ri1 = r[i-1].quo_rem(r[i])
        si1 = s[i-1] - qi * s[i]
        ti1 = t[i-1] - qi * t[i]

        # Normalizar el nuevo resto
        if not ri1.is_zero():
            inv_lc = 1 / ri1.leading_coefficient()
            ri1 *= inv_lc
            si1 *= inv_lc
            ti1 *= inv_lc

        r.append(ri1)
        s.append(si1)
        t.append(ti1)
        i += 1



# ---------- Pade approximation (n,n) desde secuencia ----------
def Pade_approximation_from_sequence(seq, n, F):
    """
    seq : lista [a0, a1, ..., a_{2n-1}] (elementos en F)
    n   : parámetro (usamos (n,n)-Padé)
    F   : campo base (GF(p) u otro)
    Retorna:
      ("ok", m) con m el polinomio minimal (rev_d(t)) o
      ("no", reason)
    Algoritmo:
      - construye f(x) = sum_{i=0..2n-1} a_i x^i (en la descripción ese f es el polinomio de la secuencia)
      - g(x) := x^(2n)
      - llama al EEA normalizado y aplica las condiciones del corolario:
         si gcd(rj,tj)=1 y x ∤ tj y cumple la congruencia, entonces
         d = max(1+deg rj, deg tj)
         retornar rev_d(tj)
    """
    R = PolynomialRing(F, 'x'); x = R.gen()
    # f(x) <- polinomio con coef seq a0 + a1 x + ... + a_{2n-1} x^{2n-1}
    f = R(list(seq))   # Sage interpreta lista como coeficientes
    g = x**(2*n)

    # Chequeos
    if f.is_zero():
        return ("no", "f es cero (secuencia nula)")
    # g no será cero, es x^(2n)

    # Usamos EEA normalizado con k = n (buscar r_j con deg < n)
    res = egcd_pade_canonico(R, f, g, k=n)
    if res[0] != "ok":
        return ("no", "eea returned unsat or failed")

    _, rj, tj, j = res

    # Condiciones extra: gcd(rj,tj)==1 (ya chequeado) y x no divide tj
    if tj % x == 0:
        return ("no", "x divide a t_j (tiene factor x)")

    # d = max{1+deg rj, deg tj}
    d = max(1 + Integer(rj.degree()), Integer(tj.degree()))
    m = rev_d(tj, d)   # rev_d(t) según especificación
    return ("ok", m)


# ---------- Wiedemann usando la subrutina Padé ----------
def Wiedemann_via_pade(A, b, max_tries=50):
    """
    Implementa el pseudocódigo:
      - si b==0 -> devuelve 1 (polinomio 1)
      - elige u aleatorio en F^n hasta obtener secuencia no nula
      - construye seq = [u^T A^i b for i=0..2n-1]
      - llama Pade_approximation_from_sequence(seq, n)
      - si m(A)*b == 0 devuelve m, sino reintenta
    """
    if b.is_zero():
        # Devuelvo el polinomio 1 en el anillo polinomial sobre el campo base
        return PolynomialRing(A.base_ring(), 'x')(1)

    n = A.nrows()
    F = A.base_ring()
    V = b.parent()
    print(V)
    tries = 0

    while tries < max_tries:
        tries += 1
        u = V.random_element()
        # Generamos la secuencia c_i = u^T A^i b para i=0..2n-1
        v = b
        seq = []
        for i in range(2*n):
            seq.append(u.dot_product(v))
            v = A * v

        # Si la secuencia es trivial (todo 0), reintenta
        if all(ci == 0 for ci in seq):
            # reintentar con otro u
            continue

        # Intenta Padé (n,n) sobre esta secuencia
        stat, res = Pade_approximation_from_sequence(seq, n, F)
        if stat != "ok":
            # reintentar
            continue

        m = res  # polinomio en R = F[x]
        # verificar m(A)*b == 0
        # Para evaluar m(A) * b: evaluamos polinomio m en la matriz A
        MA = m(A)   # Sage permite evaluar polinomio en matriz
        test = MA * b
        if test == vector(F, [0]*n):
            # éxito: devolvemos el polinomio minimal m
            return m
        # si no anula, reintentar con otro u

    raise ValueError("Wiedemann_via_pade: no convergió tras {} intentos".format(max_tries))


# ---------- Ejemplo de prueba (pequeño, estable) ----------
def test_example():
    F = GF(5)    # campo suficientemente grande para reducir degeneraciones
    A = Matrix(F, [[1,1,4],
                   [4,0,3],
                   [1,2,4]])
    b = vector(F, [3,1,2])

    print("Matrix A:\n", A)
    print("b =", b)

    m = Wiedemann_via_pade(A, b, max_tries=1000000)
    print("polinomio minimal m(x) encontrado =", m)
    print("m(A)*b =", m(A) * b)   # debe dar vector cero

test_example()