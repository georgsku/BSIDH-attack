def linearCombinations(m, n):
    for k1 in range(0 , m):
        for k2 in range(0 , n):
            yield (k1, k2)

def genCRTconstants(order):
    facorder = factor(order)
    Ms = {}
    for (l,e) in facorder:
        for ei in range(1 , e+1 ):
            Ms[l**(e)] = [l**(i) * order//(l**e) for i in range(ei)]
    ys = {}
    for (l,e) in facorder:
        for ei in range(1 , e+1 ):
            ys[l**e] = pow(Ms[l**e][0], -1 , l**e)
    Ms[1 ] = 1 
    ys[1 ] = 1 
    return Ms, ys

def recPart(recoveredSize, ks, Ms, ys):
    facSize = factor(recoveredSize)
    recovered = 0 
    for l, e in facSize:
        for i in range(e):
            try:
                recovered += ks[l**(i + 1)]*Ms[l**(i + 1)][i]*(ys[l**(i + 1)])
            except KeyError:
                for key in ks:
                    if key.quo_rem(l) == 0:
                        recovered += ks[l**e]*Ms[key][i]*(ys[key])
                        

    return recovered