from cmath import log
from distutils.log import Log
from re import A
import random

class BSIDH():

    """ TODO """
    """ Finne en bedre måte å finne bobs basis punkt """
    """ Ikke hard kode k_a hemmelighet """
    """ Ikke hardkode composition av isogeniene """

    def __init__(self, p, M, N):
        self.p = p
        self.M = M
        self.N = N

        assert (self.p + 1) % self.M == 0
        assert (self.p - 1) % self.N == 0
        assert gcd(self.M, self.N) == 1

        self.F = GF(self.p**4)
        i = sqrt(self.F(-1))

        self.E = self.E_a(self.F, 6)
        print("Starting curve: ", self.E)

        assert self.E.is_supersingular()
        assert self.E.order() == (self.p**2-1)**2

        self.P_A, self.Q_A = [sqrt(self.E.order())/(self.M)*K for K in self.E.gens()]
        self.P_B, self.Q_B = [sqrt(self.E.order())/(self.N)*K for K in self.E.gens()]

        assert self.P_A.order() == self.Q_A.order() == self.M
        assert self.P_B.order() == self.Q_B.order() == self.N

        print("Alice Factors: ", factor(self.M))
        print("Bobs Factors: ", factor(self.N))

    def derive_alice(self, k_A, pub_key_bob):
        (E_B, B_P_A, B_Q_A) = pub_key_bob
        B_S_A = B_P_A + k_A*B_Q_A
        ϕs_A= self.compute_isogeny_composition_chain(E_B, B_S_A, self.M)
        return ϕs_A.codomain().j_invariant()

    def key_gen_alice(self):
        while True:
            k_A = random.randint(1, self.M)
            S_A = self.P_A + k_A*self.Q_A
            que, r = self.M.quo_rem(S_A.order())
            if r == 0:
                break
        ϕ_A = self.compute_isogeny_composition_chain(self.E, S_A, self.M)

        pub_key_alice = (ϕ_A.codomain(), ϕ_A(self.P_B), ϕ_A(self.Q_B))
        print("Alice public key", pub_key_alice[0].j_invariant())
        return [k_A, pub_key_alice]

    def key_gen_bob(self):      
        while True:
            k_B = random.randint(1, self.N)
            S_B = self.P_B + k_B*self.Q_B
            que, r = self.N.quo_rem(S_B.order())
            if r == 0:
                break

        ϕ_B = self.compute_isogeny_composition_chain(self.E, S_B, self.N)
        
        pub_key_bob = (ϕ_B.codomain(), ϕ_B(self.P_A), ϕ_B(self.Q_A))
        print("Bobs public key", pub_key_bob[0].j_invariant())
        return [k_B, pub_key_bob]

    def derive_bob(self, k_B, pub_key_alice):
        (E_A, A_P_B, A_Q_B) = pub_key_alice
        A_S_B = A_P_B + k_B*A_Q_B
        ϕs_B = self.compute_isogeny_composition_chain(E_A, A_S_B, self.N)
        return ϕs_B.codomain().j_invariant()

    def E_a(self, F, a):
        return EllipticCurve(F, [0, a, 0, 1, 0])

    def compute_isogeny_composition_chain(self, E, S, order):
        factors = list(factor(order))
        ϕs = None
        E_tmp = E
        S_tmp = S

        for l, e in factors:
            ϕs, E_tmp, S_tmp = self.compute_l_isogeny_chain(E_tmp, S_tmp, l, e, ϕs, order)
            
        return ϕs

    def compute_l_isogeny_chain(self, E, S, l, e, ϕs, order):
        print("Computing ", e, " isogenies of degree", l)
        S_tmp = S
        E_tmp = E
        order_tmp = order

        for _ in range(e):
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