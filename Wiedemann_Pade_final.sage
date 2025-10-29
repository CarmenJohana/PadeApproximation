from sage.all import *

# ---------- Utilidades ----------
def rev_d(poly, d):
    """Reverso de grado d de un polinomio."""
    R = poly.parent()
    x = R.gen()
    coeffs = [poly.coefficient(i) for i in range(d+1)]
    res = R.zero()
    for i, a in enumerate(coeffs):
        res += a * x**(d - i)
    return res


def egcd_pade_canonico(R, pol1, pol2, k):
    
    f = R(pol1)
    g = R(pol2)
    inv_lc1= 1/f.leading_coefficient()
    inv_lc2 = 1/g.leading_coefficient()
    q = []
    i = 1
    r = [g*inv_lc2, f*inv_lc1]
    rho = [inv_lc2, inv_lc1]
    s = [R(1/rho[0]), R(0)]
    t = [R(0), R(1/rho[1])]
    
    
    while not r[i].is_zero():
        
        qi, ri1 = r[i-1].quo_rem(r[i])
        q.append(qi) # q_i

        
        if ri1!=0:
            lc = ri1.leading_coefficient()
            inv_lc = 1 / lc 
        else:
            inv_lc = 1
        rho.append(inv_lc)
        
        
        r.append(ri1/rho[i+1]) # r_{i+1} = r_{i-1} - q_i r_i
        s.append((s[i-1] - qi*s[i])/rho[i+1])
        t.append((t[i-1] - qi*t[i])/rho[i+1])
       
        if ri1.degree() < k and t[i+1].degree() <= k and r[i+1].gcd(t[i+1]) == 1:
            tj = t[i+1]/(t[i+1].constant_coefficient())

            return ("ok", ri1, tj, i+1)    
          
        i += 1
        
def Pade_approximation_from_sequence(seq, n, F):
    """
    Calcula el polinomio minimal de una secuencia usando Padé (en lugar de Berlekamp–Massey).
    """
    R = PolynomialRing(F, 'x'); x = R.gen()
    f = R(list(seq))
    g = x**(2*n)
    
   
    if f.is_zero():
        return ("no", "f es cero")
    res = egcd_pade_canonico(R, f, g, k=n)
    if res[0] != "ok":
        return ("no", "eea falló")
    _, rj, tj, j = res
    
    
    d = max(1 + Integer(rj.degree()), Integer(tj.degree()))
    m = rev_d(tj, d)
    return ("ok", m)


# ---------- Wiedemann ----------
def Wiedemann(A, b):
    """
    Implementación del algoritmo de Wiedemann con Padé para polinomio minimal.
    """
    if not A.is_square():
        raise ValueError("La matriz no es cuadrada.")
    if A.is_singular():
        raise ValueError("La matriz es singular (no tiene inversa).")

    F = A.base_ring()
    n = A.nrows()
    V = VectorSpace(F, n)

    y = V.random_element()
    while True:
        u = V.random_element()
        v = b
        seq = [u.dot_product(v)]
        for i in range(1, 2*n):
            seq.append(u.dot_product(A^i*v))
        status, m = Pade_approximation_from_sequence(seq, n, F)
        
        if status != "ok":
            continue  # vuelve a intentar con otro u
            
            
        R = m.parent(); x = R.gen()
        h = -(m - m(0)) / (m(0)*x)
        y = h(A) * b
        
        # if(m(A)*b != vector(F, [0]*n)):
            # continue
        
        if m(A)*b == vector(F, [0]*n):
            return y
     
   
# ---------- Ejemplo ----------

F = Zmod(163)
R = PolynomialRing(F, 'x'); x = R.gen()
V = VectorSpace(F,4)
VV = MatrixSpace(F, 4, 4)

n = 2
A = VV.random_element()
b = V.random_element()
print("----- Matriz A -----")
print(A)
print("----- Vector b -----")
print(b)

y = Wiedemann(A, b)
print("Resultado para y: ")
print("y: ", y)
print("A*y: ", A*y)
