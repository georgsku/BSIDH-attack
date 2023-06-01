from sage.all import parallel
import itertools

def possibly_parallel(num_cores):
    if num_cores == 1:
        def _wrap(fun):
            def _fun(args):
                for a in args:
                    yield ((a,), None), fun(a)
            return _fun
        return _wrap
    return parallel(num_cores)

def supersingular_gens(E):
    """
    Compute generators of E, assuming E is supersingular
    with smooth order (p+1)^2 with factors 2 and 3 only.
    This is faster than the PARI method.
    """
    # Find a random point of order (p+1) (probability 1/3)
    p = E.base_ring().characteristic()
    while True:
        P = E.random_point()
        if ((p+1)//2) * P != 0 and ((p+1)//3) * P != 0:
            break

    while True:
        Q = E.random_point()
        if ((p+1)//2) * Q != 0 and ((p+1)//3) * Q != 0:
            # but is it linearly independent? (probability 1/3)
            w = P.weil_pairing(Q, p+1)
            if w**((p+1)/2) != 1 and w**((p+1)//3) != 1:
                return P, Q


def generate_distortion_map(E):
    if E.a_invariants() != (0,6,0,1,0):
        raise NotImplementedError
    return E.isogeny(E.lift_x(ZZ(1)), codomain=E)

def genCRTconstants(order):
    facorder = factor(order)
    Ms = {}
    for (l,e) in facorder:
        for _ in range(1 , e+1 ):
            Ms[l**(e)] = order//(l**e)
    ys = {}
    for (l,e) in facorder:
        for _ in range(1 , e+1 ):
            ys[l**e] = Ms[l**e]**(-1) % l**e
    Ms[1] = 1 
    ys[1] = 1 
    print("Generated CRT constants, Ms:",Ms, "and ys:", ys)
    return Ms, ys

def recPart(recoveredSize, ks, Ms, ys):
    facSize = factor(recoveredSize)
    recovered = 0 
    for l, e in facSize:
        for i in range(e):
            for key in ks:
                q, r = key.quo_rem(l)
                if r == 0:
                    recovered += l^i * ks[key][i]*Ms[key]*ys[key]
    return recovered

def get_CRT_values(value, Ms, ys, recoveredSize, N):
    M = 0
    Y = 0
    recoveredFactors = {}
    for l, e in factor(N):
        recoveredFactors[l] = 0

    for l, e in factor(recoveredSize):
        recoveredFactors[l] = e

    for key in Ms:
        q,r = key.quo_rem(value)
        if r == 0:
            M = Ms[key] * value**recoveredFactors[value]
            Y = ys[key]
    return M, Y

def find_valid_c(M, N):
    factors_M = Factors(M)
    factors_N = Factors(N)

    M_list = increment_combinations(factors_M.get_list_of_exponents())
    N_list = increment_combinations(factors_N.get_list_of_exponents())

    for M_el in M_list:
        exp_after_subtract_M = factors_M.get_factors_after_subtract(M_el)
        M_tmp = factors_M.list_of_exponents_to_factors(exp_after_subtract_M)
        for b in range(1, len(N_list)):
            N_el = N_list[b]
            exp_after_subtract_N = factors_N.get_factors_after_subtract(N_el)
            N_tmp = factors_N.list_of_exponents_to_factors(exp_after_subtract_N)
            C = M_tmp - N_tmp
            if C % 4 == 1 and C > 1:
                u_v = find_u_v(C)
                if u_v != None:
                    print("Found something! with C: ",C,", u: ",u_v[0], " and v; ",u_v[1],", and list: ", exp_after_subtract_N, "Sum of list", sum(exp_after_subtract_N),"and ", factor(M_tmp), "Removed: ", factor(N/N_tmp))
                    #print("Trying to find a better alpha")
                    #M_tmp, alpha, u_v = find_bigger_alpha(M_tmp, N_tmp, u_v)
                    return C, M_tmp, N_tmp, reduce((lambda x, y: x + y), M_el), reduce((lambda x, y: x + y), N_el),  u_v[0], u_v[1], N / N_tmp
                    #return C, M_tmp, N_tmp, alpha, reduce((lambda x, y: x + y), N_el),  u_v[0], u_v[1], N / N_tmp

def find_bigger_alpha(M, N, u_v):
    M_tmp = M
    best_M = M_tmp
    alpha = 0
    for i in range(1, factor(M_tmp)[0][1]):
        M_tmp = M_tmp/2
        C = M_tmp - N
        
        # No need to check further
        if C <= 1: break

        if C % 4 == 1:
            u_v = find_u_v(C)
            if u_v != None:
                alpha = i
                best_M = M_tmp
                print("Found a better alpha", alpha, factor(best_M), "with u,v", u_v)
    return best_M, alpha, u_v

def increment_combinations(lst):
    result = [[]]
    for i in lst:
        temp = []
        for j in result:
            for k in range(i+1):
                temp.append(j+[k])
        result = temp
    result = sorted(result, key=lambda x: (sum(x), x[-1]))
    return result

def find_u_v(C):
    assert C % 4 == 1
    Q = BinaryQF([1, 0, 4])
    return Q.solve_integer(C)

def create_guesses(order):
    factors = factor(order)
    guesses = []
    for base, exp in factors:
        for _ in range(exp):
            guesses.append([ZZ(i).digits(base, padto=1)[0] for i in range(base)])
    
    if len(guesses) == 1:
        return guesses[0]
    return list(itertools.product(*guesses))

def compute_isogeny_composition_chain(E, S, order):
    factors = list(factor(order))
    ϕs = None
    E_tmp = E
    S_tmp = S

    for l, e in factors:
        ϕs, E_tmp, S_tmp = compute_l_isogeny_chain(E_tmp, S_tmp, l, e, ϕs, order)
        
    return ϕs

def compute_l_isogeny_chain(E, S, l, e, ϕs, order):
    S_tmp = S
    E_tmp = E
    order_tmp = order

    for _ in range(e):
        order_tmp = order_tmp/l
        R_tmp = S_tmp
        R_tmp = order_tmp*R_tmp
        if (R_tmp.order() > 300):
                ϕ_k = EllipticCurveHom_velusqrt(E_tmp, R_tmp)
        else:
            ϕ_k = E_tmp.isogeny(R_tmp)
        S_tmp = ϕ_k(S_tmp)
        E_tmp = ϕ_k.codomain()
        if ϕs is None:
            ϕs = ϕ_k
        else:
            ϕs = ϕ_k * ϕs

    return ϕs, E_tmp, S_tmp

class Factors:

    def __init__(self, num):
        self.value = num
        self.factors = factor(num)
        self.list = list(factor(num))
        self.list_length = len(list(factor(num)))

    def get_number_of_factors(self):
        number_of_factors = 0 
        for i in range (self.list_length):
            number_of_factors += self.factors[i][1]
        return number_of_factors

    def get_number_of_bases(self):
        return len(self.list)
    
    def get_list_of_bases(self):
        return [self.list[i][0] for i in range(self.list_length)]

    def get_list_of_exponents(self):
        return [self.list[i][1] for i in range(self.list_length)]

    def create_all_combinations_off_exponents(self):
        numbers = []
        for i in self.get_list_of_exponents():
            l = [ZZ(j).digits(i + 1, padto=1) for j in range(i + 1)]
            g = [item for sublist in l for item in sublist]
            numbers.append(g)
        return list(itertools.product(*numbers))
    
    def get_factors_after_subtract(self, subtract):
        return [ self.list[i][1] - int(subtract[i]) for i in range(self.get_number_of_bases())]

    def list_of_exponents_to_factors(self, list):
        values = [self.get_list_of_bases()[i]**list[i] for i in range(self.get_number_of_bases())]
        return reduce((lambda x, y: x * y), values)