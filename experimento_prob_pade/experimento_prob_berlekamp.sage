from sage.matrix.berlekamp_massey import berlekamp_massey
from numpy import mean
import csv

N = list(range(3,5))

for n in N:
    probs = []
    Q_vals = list(primes(5, 20))

    for q in Q_vals:
        res = []
        #x = polygen(F)
        F = Zmod(q)
        U = VectorSpace(F, n)
        for i in range(20):

            while True:
                A = random_matrix(F, n, n)
                b = random_vector(F,n)
                f,_ = A.cyclic_subspace(b, var="x", basis="iterates")
                if f.degree()==n:
                    break

            u = U.random_element()
            res.append(
                int(f == berlekamp_massey(
                    [u*(A^j)*b for j in range(2*n)]
                ))
            )

        probs.append(mean(res))

    with open(f"datos/datos_prob{n}.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["q", "prob"])
        for a,b in zip(Q_vals, probs):
            writer.writerow([a,b])


