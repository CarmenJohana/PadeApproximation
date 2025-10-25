from sage.matrix.berlekamp_massey import berlekamp_massey
def Wiedemann(A,b):
    
    if A.is_singular():
        raise ValueError("la matriz es singular")
    if not A.is_square():
        raise ValueError("La matriz no es cuadrada")

    V = b.parent()  
    F = A.base_ring()
    R = F['x']
    x = R.gen()

    y = V([1,1])
    
    while A*y != b:
        
        u = V.random_element()
        v = b
        c = []
        c.append(u.dot_product(v))

        for i in range(1, 2*n):
            v = A * v
            c.append(u*v)
              
        m = berlekamp_massey(c)
        h =  (m(0) - m)/x
        h = h/m(0)
        y = h(A) * b
        
    return y
F = Zmod(7)
R = PowerSeriesRing(F, 'x')
V = VectorSpace(F, 2)
VV = MatrixSpace(F, 2, 2)
u = V.random_element()
A = VV.random_element()
b = V.random_element()
u = V.random_element()
n = 2
c = [u.dot_product(A**i * b) for i in range(2*n-1, -1, -1)]
c_2 = [ u * A**i * b for i in range(2*n-1, -1, -1)]
m_in = R(c)

A = VV.random_element()
b = V.random_element()
print(A)
print(b)
y = Wiedemann(A,b)
print(y)

print("Finalizado")