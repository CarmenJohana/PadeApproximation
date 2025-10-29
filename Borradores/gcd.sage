def egcd_traditional(q, pol1, pol2):
    
    R = PolynomialRing(Zmod(q), 'x')
    f = R(pol1)
    g = R(pol2)
    inv_lc1= 1/f.leading_coefficient()
    inv_lc2 = 1/g.leading_coefficient()
    
    r = [f*inv_lc1, g*inv_lc2,]

    q = []
    i = 1
    rho = [inv_lc1, inv_lc2]
    s = [R(1/rho[0]), R(0)]
    t = [R(0), R(R(1)/R(rho[1]))]
    
    while not r[i].is_zero():
        
        
        qi, ri1 = r[i-1].quo_rem(r[i])
        
        q.append(qi) # q_i
        
        
        if ri1!=0:
            lc = ri1.leading_coefficient()
            inv_lc = 1 / lc
            ri1 *= inv_lc
            
        else:
            inv_lc = 1
        rho.append(inv_lc)
        
        r.append(ri1) # r_{i+1} = r_{i-1} - q_i r_i
        #print(i, rho)
        s.append((s[i-1] - qi*s[i])*rho[i+1])    # s_{i+1} = s_{i-1} - q_i s_i
        t.append((t[i-1] - qi*t[i])*rho[i+1])    # t_{i+1} = t_{i-1} - q_i t_i
        i += 1
        
     
    return r, t


print(egcd_traditional(5,[0,0,0,0,1], [1,2,3,4]))