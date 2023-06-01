from cmath import log
from distutils.log import Log
from re import A
import random

class SIDH():

    def __init__(self, e_A, e_B):
        self.e_A = e_A
        self.e_B = e_B
        
        p = (2**e_A) * (3**e_B) - 1
        print("Initing sidh with eA and eB: ", e_A, e_B)
        F.<i> = GF(p**2, modulus=x**2 + 1)
        # Starting curve
        self.E = self.E_a(F, 6)
        assert self.E.is_supersingular()

        #self.P_A, self.Q_A = [sqrt(self.E.order())/(2^e_A)*K for K in self.E.gens()]
        self.P_A = self.E(70*i + 332, 371*i + 191)
        self.Q_A = self.E(269*i + 403, 302*i + 3)
        assert self.P_A.order() == self.Q_A.order() == 2^self.e_A

        #self.P_B, self.Q_B = [sqrt(self.E.order())/(3^e_B)*K for K in self.E.gens()]
        self.P_B = self.E(i + 420, 112*i + 230)
        self.Q_B = self.E(265*i + 329, 377*i + 124)
        assert self.P_B.order() == self.Q_B.order() == 3^self.e_B

        print(self.E.j_invariant())


    def derive_alice(self, k_A, pub_key_bob):
        (E_B, B_P_A, B_Q_A) = pub_key_bob
        B_S_A = B_P_A + k_A*B_Q_A
        print("B_S_A", B_S_A)
        print("j_invariant", E_B.j_invariant())
        ϕs_A = self.compute_isogenies(E_B, B_S_A, 2, self.e_A)
        # shared j_invariant
        return ϕs_A[-1].codomain().j_invariant()

    def key_gen_alice(self):
        #k_A = random.randint(1, (2**self.e_A - 1))
        S_A = self.P_A + k_A*self.Q_A
        # S_A remains in E[2^e_A]
        assert S_A.order() == 2^self.e_A
        # Alice computes secret 2^4-isogeny
        ϕ_A_direct = self.E.isogeny(S_A)
        assert ϕ_A_direct.degree() == 2^self.e_A

        ϕs_A = self.compute_isogenies(self.E, S_A, 2, 4)
        # we can get is a "composite map".
        ϕ_A = ϕs_A[3] * ϕs_A[2] * ϕs_A[1] * ϕs_A[0]
        # Alice's public key.
        pub_key = (ϕ_A.codomain(), ϕ_A(self.P_B), ϕ_A(self.Q_B))
        return [k_A, pub_key]

    def key_gen_bob(self):      
        #k_B = random.randint(1, (3**self.e_B - 1))
        S_B = self.P_B + k_B*self.Q_B
        assert S_B.order() == 3^self.e_B

        ϕs_B = self.compute_isogenies(self.E, S_B, 3, self.e_B)
        ϕ_B = ϕs_B[2] * ϕs_B[1] * ϕs_B[0]

        pub_key = (ϕ_B.codomain(), ϕ_B(self.P_A), ϕ_B(self.Q_A))
        return [k_B, pub_key]

    def derive_bob(self, k_B, pub_key_alice):
        (E_A, A_P_B, A_Q_B) = pub_key_alice
        A_S_B = A_P_B + k_B*A_Q_B
        ϕs_B = self.compute_isogenies(E_A, A_S_B, 3, self.e_B)
        return ϕs_B[-1].codomain().j_invariant()

    def E_a(self, F, a):
        return EllipticCurve(F, [0, a, 0, 1, 0])

    def compute_isogenies(self, E, S, l, e):
        """Returns a list of e l-isogenies making up the target isogeny."""
        S_tmp = S
        E_tmp = E
        ϕs = []
        for k in range(e):
            # Generate a point in the l-torsion via repeated
            # doubling (l=2) or tripling (l=3) starting from S_tmp.
            R_tmp = S_tmp
            for _ in range(e-k-1):
                R_tmp = l*R_tmp
            assert R_tmp.order() == l
            # Generate an l-isogeny to take a step in the graph.
            ϕ_k = E_tmp.isogeny(R_tmp)
            assert ϕ_k.degree() == l
            S_tmp = ϕ_k(S_tmp)
            # Each S_tmp lies above a point in the l-isogeny's kernel,
            # so its order decreases by a factor of l per step.
            assert S_tmp.order() == l^(e-k-1)
            E_tmp = ϕ_k.codomain()
            ϕs.append(ϕ_k)
        return ϕs
