load('castryck_decru_attack_bsidh.sage')
#load('castryck_decru_attack_bsidh_uvtable.sage')
load('bsidh.sage')
load('uvtable_example_2.sage')

debug = False

p = 53380256806900813486663140872422666278845625642838258109389995943373349847039
M = 2^110
N = 3^34 * 11 * 17 * 19^2 * 29 * 37 * 53^2 * 97 * 107

Alices_factors = Factors(M)
bobs_factors = Factors(N)

a = Alices_factors.get_number_of_factors()
b = bobs_factors.get_number_of_factors()

bsidh = BSIDH(p, M, N)

alice_private_key, alice_pub_key = bsidh.key_gen_alice()
bob_private_key, bob_pub_key = bsidh.key_gen_bob()

alice_j_invariant = bsidh.derive_alice(alice_private_key, bob_pub_key)
bob_j_invariant = bsidh.derive_bob(bob_private_key, alice_pub_key)

if alice_j_invariant == bob_j_invariant:
    print("Alice and bob calculated the same shared secret: \n", alice_j_invariant)
else:
    print("Alice and bob calculated different secrets.. ")
    print(alice_j_invariant, bob_j_invariant)

# Generation of the endomorphism 2i
two_i = generate_distortion_map(bsidh.E)

p = bsidh.p
E = bsidh.E
P2 = bsidh.P_A
Q2 = bsidh.Q_A
P3 = bsidh.P_B
Q3 = bsidh.Q_B
EB = bob_pub_key[0]
P_B = bob_pub_key[1]
Q_B = bob_pub_key[2]

print(f"Running the attack against Baby BISDH parameters, which has a prime: {p}")
print(f"If all goes well then the following digits should be found:", bob_private_key)
print("\n --- BEGINNING ATTACK --- \n ")

# ===================================
# =====  ATTACK  ====================
# ===================================


def RunAttackOnBSIDH(num_cores):
    return CastryckDecruAttack(E, P2, Q2, EB, P_B, Q_B, two_i, num_cores=num_cores)

if __name__ == '__main__' and '__file__' in globals():
    if '--parallel' in sys.argv:
        # Set number of cores for parallel computation
        num_cores = 4
    else:
        num_cores = 1
    
    print(f"Performing the attack in parallel using {num_cores} cores")
    recovered_key = RunAttackOnBSIDH(num_cores)
