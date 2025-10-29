from sage.all import *

# ---------- Utilidades ----------
def rev_d(poly, d):
    """Reverso de grado d de un polinomio."""
    R = poly.parent()
    x = R.gen()
    coeffs = [poly[i] for i in range(d+1)]  # <- aquí
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
    #r = [g*inv_lc2, f*inv_lc1]
    r = [g*inv_lc2, f*inv_lc1]
    rho = [g.leading_coefficient(),f.leading_coefficient()]
    s = [R(1/rho[0]), R(0)]
    t = [R(0), R(1/rho[1])]
    
    
    while not r[i].is_zero():
        
        qi, ri1 = r[i-1].quo_rem(r[i]) # q[i] = r[i-1] quo r[i], ri1 = r[i-1] - q[i]*r[i]
        q.append(qi) # q_i
        lc = 1
        
        if ri1!=0:
            lc = ri1.leading_coefficient()
            inv_lc = 1 / lc 
            
        else:
            inv_lc = 1
        rho.append(lc)
        
        r.append(ri1 * inv_lc)            # r_{i+1} ahora es monico
        s.append((s[i-1] - qi*s[i]) * inv_lc)
        t.append((t[i-1] - qi*t[i]) * inv_lc)
        #r.append(ri1/rho[i+1]) # r_{i+1} = r_{i-1} - q_i r_i
        #s.append((s[i-1] - qi*s[i])/rho[i+1])
        #t.append((t[i-1] - qi*t[i])/rho[i+1])
       
        if ri1.degree() < k and t[i+1].degree() <= k and r[i+1].gcd(t[i+1]) == 1:
            tj = t[i+1]/(t[i+1].constant_coefficient())

            return ("ok", ri1, tj, i+1)    
          
        i += 1
        
    return ("no", None, None, None)
        
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


from sage.all import *
import csv

N = list(range(3,5))  # dimensiones
num_trials = 50       # repeticiones por combinación q,n
Q_vals = list(primes(5, 50))  # más primos

for n in N:
    probs = []

    for q in Q_vals:
        res = []
        F = Zmod(q)
        U = VectorSpace(F, n)
        
        for _ in range(num_trials):
            while True:
                A = random_matrix(F, n, n)
                b = random_vector(F, n)
                f,_ = A.cyclic_subspace(b, var="x", basis="iterates")
                if f.degree() == n:
                    break

            u = U.random_element()
            status, m = Pade_approximation_from_sequence([u*(A^j)*b for j in range(2*n)], n, F)
            if status == "ok":
                res.append(int(f == m))
            else:
                res.append(0)

        probs.append(mean(res))

    with open(f"datos/datos_prob{n}.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["q", "prob"])
        for a, b in zip(Q_vals, probs):
            writer.writerow([a, b])

