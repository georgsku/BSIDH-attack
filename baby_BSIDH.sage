load('castryck_decru_attack_bsidh.sage')
#load('castryck_decru_attack_bsidh_uvtable.sage')
load('bsidh.sage')
load('uvtable_example.sage')

debug = False
running_scheme = "BSIDH"

SIKE_parameters = {
    "p7" : (127),
    "p13" : (8191),
    "p17" : (131071),
    "p19" : (524287),
    "p31" : (2147483647),
    "p61" : (2305843009213693951),  
    "p127" : (170141183460469231731687303715884105727),  
}

bsidh_settings = "p13"
p = SIKE_parameters[bsidh_settings]
M = (p + 1)
N = (p - 1)/(2)

Alices_factors = Factors(M)
bobs_factors = Factors(N)

a = Alices_factors.get_number_of_factors()
b = bobs_factors.get_number_of_factors()

sidh = BSIDH(p, M, N)

alice_private_key, alice_pub_key = sidh.key_gen_alice()
bob_private_key, bob_pub_key = sidh.key_gen_bob()

alice_j_invariant = sidh.derive_alice(alice_private_key, bob_pub_key)
bob_j_invariant = sidh.derive_bob(bob_private_key, alice_pub_key)

if alice_j_invariant == bob_j_invariant:
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
        num_cores = 2
    else:
        num_cores = 1
    
    print(f"Performing the attack in parallel using {num_cores} cores")
    recovered_key = RunAttackOnBSIDH(num_cores)
