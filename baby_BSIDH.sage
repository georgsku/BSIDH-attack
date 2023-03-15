load('castryck_decru_attack_bsidh.sage')
load('bsidh.sage')

p = 8191
M = 2**12
N = 3**2 * 5 * 7 * 13

Alices_factors = Factors(M)
bobs_factors = Factors(N)

a = Alices_factors.get_number_of_factors()
b = bobs_factors.get_number_of_factors()

sidh = BSIDH(p, M, N)

alice_private_key, alice_pub_key = sidh.key_gen_alice()
bob_private_key, bob_pub_key = sidh.key_gen_bob()

alice_j_invariant = sidh.derive_alice(alice_private_key, bob_pub_key)
bob_j_invariant = sidh.derive_bob(bob_private_key, alice_pub_key)

if alice_j_invariant == alice_j_invariant:
    print("Alice and bob calculated the same shared secret: \n", alice_j_invariant)
else:
    print("Alice and bob calculated different secrets.. ")
    print(alice_j_invariant, bob_j_invariant)

# Generation of the endomorphism 2i
two_i = generate_distortion_map(sidh.E)

p = sidh.p
E = sidh.E
P2 = sidh.P_A
Q2 = sidh.Q_A
P3 = sidh.P_B
Q3 = sidh.Q_B
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
        num_cores = os.cpu_count()
        print(f"Performing the attack in parallel using {num_cores} cores")
    else:
        num_cores = 1
    recovered_key = RunAttackOnBSIDH(num_cores)
