import random
import time

class SIDH():

    def __init__(self, e_A, e_B):
        self.e_A = e_A
        self.e_B = e_B
        
        self.p = (2**e_A) * (3**e_B) - 1
        print("Initing sidh with eA and eB: ", e_A, e_B)
        F.<i> = GF(self.p**2, modulus=x**2 + 1)
        # Starting curve
        self.E = self.E_a(F, 6)
        assert self.E.is_supersingular()

        self.P_A, self.Q_A = [sqrt(self.E.order())/(2^e_A)*K for K in self.E.gens()]
        assert self.P_A.order() == self.Q_A.order() == 2^self.e_A

        self.P_B, self.Q_B = [sqrt(self.E.order())/(3^e_B)*K for K in self.E.gens()]
        assert self.P_B.order() == self.Q_B.order() == 3^self.e_B


    def derive_alice(self, k_A, pub_key_bob):
        (E_B, B_P_A, B_Q_A) = pub_key_bob
        B_S_A = B_P_A + k_A*B_Q_A
        ϕs_A = self.compute_isogenies(E_B, B_S_A, 2, self.e_A)
        # shared j_invariant
        return ϕs_A[-1].codomain().j_invariant()

    def key_gen_alice(self):
        k_A = random.randint(1, (2**self.e_A - 1))
        S_A = self.P_A + k_A*self.Q_A
        # S_A remains in E[2^e_A]
        assert S_A.order() == 2^self.e_A
        # Alice computes secret 2^4-isogeny
        ϕ_A_direct = self.E.isogeny(S_A)
        assert ϕ_A_direct.degree() == 2^self.e_A

        ϕs_A = self.compute_isogenies(self.E, S_A, 2, self.e_A)
        # we can get is a "composite map".
        ϕ_A = ϕs_A[3] * ϕs_A[2] * ϕs_A[1] * ϕs_A[0]
        # Alice's public key.
        pub_key = (ϕ_A.codomain(), ϕ_A(self.P_B), ϕ_A(self.Q_B))
        return [k_A, pub_key]

    def key_gen_bob(self, key):      
        k_B = random.randint(1, (3**self.e_B - 1))
        k_B = key
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
        for k in range(e):
            R_tmp = S_tmp
            for _ in range(e-k-1):
                R_tmp = l*R_tmp
            assert R_tmp.order() == l

            ϕ_k = E_tmp.isogeny(R_tmp)
            assert ϕ_k.degree() == l

            S_tmp = ϕ_k(S_tmp)
            assert S_tmp.order() == l^(e-k-1)

            E_tmp = ϕ_k.codomain()
            if ϕs is None:
                ϕs = ϕ_k
            else:
                ϕs = ϕ_k * ϕs
        return ϕs

for i in range(3):
    tid = time.time()
    sidh = SIDH(216, 137)
    alice_private_key, alice_pub_key = sidh.key_gen_alice()
    bob_private_key, bob_pub_key = sidh.key_gen_bob()

    alice_j_invariant = sidh.derive_alice(alice_private_key, bob_pub_key)
    bob_j_invariant = sidh.derive_bob(bob_private_key, alice_pub_key)
    if (alice_j_invariant == bob_j_invariant):
        print("complete!")
        print("it took: ", time.time() - tid)
