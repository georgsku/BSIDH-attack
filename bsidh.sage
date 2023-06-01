from re import A
import random

class BSIDH():

    def __init__(self, p, M, N):
        self.p = p
        self.M = M # (p + 1)-side
        self.N = N # (p - 1)-side
        
        _, r1 = (self.p + 1).quo_rem(self.M)
        _, r2 = (self.p - 1).quo_rem(self.N)
        assert r1 == 0
        assert r2 == 0
        assert gcd(self.M, self.N) == 1

        F.<z4> = GF(p^4, name="z4")
        self.F = F
        self.Fp2 = GF(p^2, name="i", modulus=x^2 + 1)
        Fp4.<z4> = GF(p^4, name="z4")
        self.Fp4 = Fp4
        self.Fp4.register_coercion(Hom(self.Fp2,self.Fp4).an_element())
        self.F = self.Fp4
        i = sqrt(self.F(-1))

        self.E = EllipticCurve(self.F, [0, 6, 0, 1, 0])

        assert self.E.is_supersingular()
        assert self.E.order() == (self.p**2-1)**2
        self.P_A, self.Q_A = [sqrt(self.E.order())/(self.M)*K for K in self.E.gens()]
        self.P_B, self.Q_B = [sqrt(self.E.order())/(self.N)*K for K in self.E.gens()]
        
        """ self.P_A = self.E.lift_x(136*z4^3 + 17*z4^2 + 204*z4 + 123) 
        self.Q_A = self.E.lift_x(218*z4^3 + 135*z4^2 + 327*z4 + 222) 
        self.P_B = self.E.lift_x(380*z4^3 + 263*z4^2 + 139*z4 + 292) 
        self.Q_B = self.E.lift_x(336*z4^3 + 42*z4^2 + 73*z4 + 120)  """

        print(self.P_A)
        print(self.Q_A)

        print(self.P_B)
        print(self.Q_B)

        assert self.P_A.order() == self.Q_A.order() == self.M
        assert self.P_B.order() == self.Q_B.order() == self.N
        print("Staring curve with j_invariant", self.E.j_invariant())
        #G = self.E.isogeny_ell_graph(5, directed=False, label_by_j=True)
        #print(G.edges())


    def derive_alice(self, k_A, pub_key_bob):
        (E_B, B_P_A, B_Q_A) = pub_key_bob
        B_S_A = B_P_A + k_A*B_Q_A
        ϕs_A= self.compute_isogeny_composition_chain(E_B, B_S_A, self.M)
        return ϕs_A.codomain().j_invariant()

    def key_gen_alice(self):
        while True:
            k_A = random.randint(0, self.M)
            S_A = self.P_A + k_A*self.Q_A
            _, r = self.M.quo_rem(S_A.order())
            if r == 0:
                break
        ϕ_A = self.compute_isogeny_composition_chain(self.E, S_A, self.M)

        pub_key_alice = (ϕ_A.codomain(), ϕ_A(self.P_B), ϕ_A(self.Q_B))
        print("alice secret generator:", S_A)
        print("Alice public key", pub_key_alice[0].j_invariant())
        print("Alice public key", self.Fp2(self.Fp4(pub_key_alice[0].j_invariant())))
        return [k_A, pub_key_alice]

    def key_gen_bob(self):      
        while True:
            k_B = random.randint(0, self.N)
            k_B = 3550
            S_B = self.P_B + k_B*self.Q_B
            _, r = self.N.quo_rem(S_B.order())
            if r == 0:
                break

        ϕ_B = self.compute_isogeny_composition_chain(self.E, S_B, self.N)
        print("boib secret generator:", S_B)
        pub_key_bob = (ϕ_B.codomain(), ϕ_B(self.P_A), ϕ_B(self.Q_A))
        print("Bobs public key", pub_key_bob[0].montgomery_model(), pub_key_bob)
        print("Bobs public key", self.Fp2(self.Fp4(pub_key_bob[0].j_invariant())))
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
        #print("Computing ", e, " isogenies of degree", l, "from node", E.j_invariant())
        S_tmp = S
        E_tmp = E
        order_tmp = order

        for _ in range(e):
            order_tmp = order_tmp/l
            R_tmp = S_tmp
            R_tmp = order_tmp*R_tmp
            #ϕ_k = E_tmp.isogeny(R_tmp)
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
            #print("intermediate node:", E_tmp.j_invariant(), "Translated into:", self.Fp2(self.Fp4(E_tmp.j_invariant())), "curve:", E_tmp.montgomery_model())
        return ϕs, E_tmp, S_tmp