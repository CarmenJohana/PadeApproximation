from sage.all import *

# ---------- utilidades ----------
def rev_d(poly, d):
    """Reverso de grado d de un polinomio."""
    R = poly.parent()
    x = R.gen()
    coeffs = [poly.coefficient(i) for i in range(d+1)]
    res = R.zero()
    for i, a in enumerate(coeffs):
        res += a * x**(d - i)
    return res

'''
def egcd_pade_canonico(R, pol_f, pol_g, k):
    """
    Extended Euclidean Algorithm canónico para Padé sobre R = GF(p)[x].
    Normaliza t(0)=1 cuando es posible.
    """
    f = R(pol_f)
    g = R(pol_g)

    if f.is_zero() or g.is_zero():
        raise ValueError("f y g deben ser no nulos.")

    inv_lc_f = 1 / f.leading_coefficient()
    inv_lc_g = 1 / g.leading_coefficient()

    r = [f * inv_lc_f, g * inv_lc_g]
    s = [R(inv_lc_f), R(0)]
    t = [R(0), R(inv_lc_g)]

    i = 1
    while True:
        if r[i].is_zero():
            j = i - 1
        elif r[i].degree() < k:
            j = i
        else:
            qi, ri1 = r[i-1].quo_rem(r[i])
            si1 = s[i-1] - qi * s[i]
            ti1 = t[i-1] - qi * t[i]
            if not ri1.is_zero():
                inv_lc = 1 / ri1.leading_coefficient()
                ri1 *= inv_lc
                si1 *= inv_lc
                ti1 *= inv_lc
            r.append(ri1)
            s.append(si1)
            t.append(ti1)
            i += 1
            continue

        rj, tj = r[j], t[j]
        if rj.gcd(tj) == 1:
            if not tj(0).is_zero():
                inv_t0 = 1 / tj(0)
                rj *= inv_t0
                tj *= inv_t0
            tau = tj.leading_coefficient()
            return ("ok", rj / tau, tj / tau, j)
        else:
            return ("unsat", rj, tj, j)
'''



def egcd_pade_canonico(R, pol1, pol2, k):
    
    #R = PolynomialRing(Zmod(q), 'x')
    f = R(pol1)
    g = R(pol2)
    inv_lc1= 1/f.leading_coefficient()
    inv_lc2 = 1/g.leading_coefficient()
    
    r = [f,g]
    #r = [f*inv_lc1, g*inv_lc2]
    q = []
    i = 1
    rho = [inv_lc1, inv_lc2]
    #s = [R(1/rho[0]), R(0)]
    t = [R(0), R(R(1)/R(rho[1]))]
    s = [R(1), R(0)]
    #t = [R(0), R(1)]
    while not r[i].is_zero():
        
        
        qi, ri1 = r[i-1].quo_rem(r[i])
        
        q.append(qi) # q_i
        
        
        if ri1!=0:
            lc = ri1.leading_coefficient()
            inv_lc = 1 / lc
            #ri1 *= inv_lc
            
        else:
            inv_lc = 1
        rho.append(inv_lc)
        
        
        r.append(ri1) # r_{i+1} = r_{i-1} - q_i r_i
        #print(i, rho)
        #s.append((s[i-1] - qi*s[i])*rho[i+1])    # s_{i+1} = s_{i-1} - q_i s_i
        #t.append((t[i-1] - qi*t[i])*rho[i+1])    # t_{i+1} = t_{i-1} - q_i t_i
        s.append(s[i-1] - qi*s[i])
        t.append(t[i-1] - qi*t[i])
        
        if ri1.degree()<k and s[i+1].degree()<=k and ri1.gcd(s[i+1])==1:
            print("s: ", s[i+1])
            print("s: ", s[i+1]/s[i+1](0))
            
            return ("ok", ri1, s[i+1]/s[i+1](0), i+1)
        
        i += 1
        print("q : ", q)
        print("r : ", r)
        print("t: ", t)
        print("s: ", s)
        
     
    return r, t


def Pade_approximation_from_sequence(seq, n, F):
    """
    Calcula el polinomio minimal de una secuencia usando Padé (en lugar de Berlekamp–Massey).
    """
    R = PolynomialRing(F, 'x'); x = R.gen()
    f = R(list(seq))
    g = x**(2*n)
    
    print("f: ",f, " and ", g) 
    if f.is_zero():
        return ("no", "f es cero")
    res = egcd_pade_canonico(R, f, g, k=n)
    if res[0] != "ok":
        return ("no", "eea falló")
    _, rj, tj, j = res
    
    print("t_j: ", tj)
    
    if tj % x == 0:
        return ("no", "x divide t_j")
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
        #u = vector(F, [0,6])
        print("Valor de u: ", u)
        
        v = b
        seq = [u.dot_product(v)]
        for i in range(1, 2*n):
            v = A * v
            seq.append(u.dot_product(v))
        print("Secuencia: ", seq)
        status, m = Pade_approximation_from_sequence(seq, n, F)
        
        print("m:", m)
        if status != "ok":
            #continue  # vuelve a intentar con otro u
            return
        # Si m(x) = c0 + c1 x + ... + ck x^k
        # definimos h(x) = (m(x) - m(0)) / (x * m(0))
        R = m.parent(); x = R.gen()
        h = -(m - m(0)) / (m(0)*x)

        print("h: ", h)
        y = h(A) * b
        
        print("El valor de m(A)*b es: ", m(A)*b)
        
        if m(A)*b == vector(F, [0]*n):
            return y


# ---------- Ejemplo ----------

F = Zmod(7)
R = PowerSeriesRing(F, 'x')
V = VectorSpace(F, 2)
VV = MatrixSpace(F, 2, 2)
A = matrix(F, [[0,6], [2,6]])
b = vector(F, [4,1])
print(A)
print(b)
y = Wiedemann(A,b)
print(y)
