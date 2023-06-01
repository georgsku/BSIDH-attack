# Python imports
import time
from itertools import product
from multiprocessing import Pool


# Load Sage Files
load('helpers.sage')
load('richelot_aux.sage')
load('speedup.sage')
load('uvtable.sage')

# ===================================
# =====  ATTACK  ====================
# ===================================

def CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i, num_cores=1):
    tim = time.time()

    skB = [] # TERNARY DIGITS IN EXPANSION OF BOB'S SECRET KEY
    digits_left = b
    remaining_order_M = 2^a
    remaining_order_N = 3^b
    ai = a
    bi = b

    while digits_left >= 1:

        _, remaining_order_M, remaining_order_N, alpha_i, beta_i, u, v, order,  = find_valid_c(remaining_order_M, remaining_order_N)

        ai = ai - alpha_i
        bi = bi - beta_i
        alp = a - ai


        print(f"Determination of first {beta_i} ternary digits. We are working with 2^{ai}-torsion.")

        @possibly_parallel(num_cores)
        def CheckGuess(first_digits):
            print(f"Testing digits: {first_digits}")
            digits = skB + first_digits
            scalar = sum(3^k*d for k,d in enumerate(digits))

            tauhatkernel = 3^bi * (P3 + scalar*Q3)
            assert tauhatkernel.order() == 3^(b - bi)

            C, P_c, Q_c, _ = AuxiliaryIsogeny(beta_i, u, v, E_start, P2, Q2, tauhatkernel, two_i)
            assert P_c.order() == 2^a

            return Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai)

        guesses = [ZZ(i).digits(3, padto=beta_i) for i in range(3^beta_i-1)]
        
        for result in CheckGuess(guesses):
            ((first_digits,), _), is_split = result
            if is_split is not None:
                print("Glue-and-split! These are most likely the secret digits.")
                skB += first_digits
                break

        else:
            print("All other guesses failed, so first digits must be all 2!")
            skB += [2]*beta_i

        print(skB)
        if b - len(skB) == 1: 
            break


        
    @possibly_parallel(num_cores)
    def CheckGuess(i):
        bobskey = key + i*3^(b-1)
        bobscurve, _ = Pushing3Chain(E_start, P3 + bobskey*Q3, b)
        return bobscurve.j_invariant() == EB.j_invariant()

    key = sum([skB[i]*3^(i) for i in range(b-1)])
    print(f"Determination of last {1} ternary digits. We are brute-forcing this.")

    for result in CheckGuess([0..2]):
        ((i,), _), found = result
        if found:
            bobskey = key + i*3^(b-1)
            break

    if found:
        print(f"Bob's secret key revealed as: {bobskey}")
        print(f"In ternary, this is: {Integer(bobskey).digits(base=3)}")
        return bobskey
    else:
        print("Something went wrong.")
        return None