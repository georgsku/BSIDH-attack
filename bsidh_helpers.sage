import itertools
def find_u_v(c):
    u = 0
    v = 0
    iterations = 0
    assert c % 4 == 1
    while True:
        if iterations > 10*c:
            return False, False
        if u > c:
            print("Sorry, couldnt find anything ")
            return False, False
        if u**2 + 4*v**2 == c:
            return u,v
        if u**2 + 4*v**2 > c:
            v += 1
            u = 0
        u += 1
        iterations += 1

def find_valid_c(M, N):
    print("Finding a suitible C")
    factors_M = Factors(M)
    factors_N = Factors(N)

    M_list = factors_M.create_all_combinations_off_exponents()
    N_list = factors_N.create_all_combinations_off_exponents()
    for a in range(1, len(M_list)):
        exp_after_subtract_M = factors_M.get_factors_after_subtract(M_list[a])
        M_tmp = factors_M.list_of_exponents_to_factors(exp_after_subtract_M)
        for b in range(1, len(N_list)):
            exp_after_subtract_N = factors_N.get_factors_after_subtract(N_list[b])
            N_tmp = factors_N.list_of_exponents_to_factors(exp_after_subtract_N)
            C = M_tmp - N_tmp
            if C % 4 == 1 and C > 1:
                print(C, exp_after_subtract_N)
                u, v = find_u_v(C)
                if u and v:
                    print("Found something! with C: ",C,", u: ",u, " and v; ",v,", and list: ", exp_after_subtract_N, "and ", factor(M_tmp) )
                    return C, M_tmp, N_tmp, reduce((lambda x, y: x + y), M_list[a]), reduce((lambda x, y: x + y), N_list[b]),  u, v, N / N_tmp

def get_inverse(M, n):
    y = 1
    while True:
        if M/n * y % n == 1:
            return y
        y += 1

def computeCompositionIsogeny(E, S, order):
    factors = list(factor(order))
    ϕs = None
    index = 0
    E_tmp = E
    S_tmp = S
    while True: 
        try:   
            factors[index]
        except IndexError:
            """ List out of range """
            break
        ϕs, E_tmp, S_tmp = computeIsogeny(E_tmp, S_tmp, factors[index][0], factors[index][1], ϕs, order)
        index = index + 1
    return ϕs

def computeIsogeny(E, S, l, e, ϕs, order):
    S_tmp = S
    E_tmp = E
    order_tmp = order
    for k in range(e):
        order_tmp = order_tmp/l
        R_tmp = S_tmp
        R_tmp = order_tmp*R_tmp
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