import time
from itertools import product
from helpers import possibly_parallel
import itertools

load('richelot_aux.sage')
load('bsidh_helpers.sage')
load('CRTrepresentasjon.sage')

# ===================================
# =====  ATTACK  ====================
# ===================================

def CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i, num_cores=1):
    tim = time.time()

    ks = {}
    factorOrder = {}
    recoveredSizeFactor = {}
    for l, e in factor(N):
        ks[l**e] = -1  #-1 means not yet recovered (no need to initialize this, other than readability)
        factorOrder[l] = e
        recoveredSizeFactor[l] = 0

    recoveredSize = 1
    remainingOrderM = M
    remainingOrderN = N
    ai = a
    bi = b

    Ms, ys = genCRTconstants(N)

    while recoveredSize != N:
        print("Remaining order:", factor(remainingOrderN))
        _, remainingOrderM, remainingOrderN, alpha_i, beta_i, u, v, order = find_valid_c(remainingOrderM, remainingOrderN)
        
        ai = ai - alpha_i
        bi = bi - beta_i
        alp = a - ai
        facOrder = factor(order)

        print(f"Determination of the next {beta_i} digits. We are working with 2^{ai}-torsion.")

        @possibly_parallel(num_cores)
        def CheckGuess(ki):
            print(f"Testing digits: {ki}")

            recoveredPart = recPart(recoveredSize, ks, Ms, ys)
            scalar = recoveredPart + sum(k * Ms[facOrder[i][0]**factorOrder[facOrder[i][0]]][recoveredSizeFactor[facOrder[i][0]]] * ys[facOrder[i][0]**factorOrder[facOrder[i][0]]] for i,k in enumerate(ki))

            tauhatkernel = remainingOrderN * (P3 + scalar*Q3)
            assert tauhatkernel.order() == recoveredSize * order

            C, P_c, Q_c, _ = AuxiliaryIsogeny(beta_i, u, v, E_start, P2, Q2, tauhatkernel, two_i, order)
            assert P_c.order() == M

            return Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai)

        factors = Factors(order).get_list_of_bases()
        guesses = []
        for base in factors:
            guesses.append([ZZ(i).digits(base, padto=1)[0] for i in range(base)])

        guesses = list(itertools.product(*guesses))
        
        print("We will try the following digits", guesses)
        for result in CheckGuess(guesses):
            ((ki,), _), is_split = result
            if is_split is not None:
                print("Glue-and-split! These are most likely the secret digits.", ki)
                print(ki)
                for i,k in enumerate(ki):
                    ks[facOrder[i][0]] = k
                break

        else:
            print("Sorry, something wrong happened. Im out of guesses")
            print("However i found these digits:", ks)
            return None

        recoveredSize *= order
        for o in factor(order):
            recoveredSizeFactor[o[0]] += 1
        print(ks)

        if (len(list(factor(remainingOrderN))) == 1):
            break

    print(f"Determination of last 1 CRT digit. We are brute-forcing this.")
    
    def CheckGuess(i):
        recoveredPart = recPart(recoveredSize, ks, Ms, ys)
        scalar = i * Ms[remainingOrderN**factorOrder[remainingOrderN]][recoveredSizeFactor[remainingOrderN]] * ys[remainingOrderN ** factorOrder[remainingOrderN]]
        bobskey = (recoveredPart + scalar) % N
        print(i, bobskey)

        bobs_S_B = P3 + bobskey*Q3
        ϕ_B = compute_isogeny_composition_chain(E_start, bobs_S_B, N)
        return ϕ_B.codomain().j_invariant() == EB.j_invariant()

    for i in range(9):
        result = CheckGuess(i)
        if result:
            print("Found correct last digit:", i)
            ks[remainingOrderN ** factorOrder[remainingOrderN]] = i
            recoveredSizeFactor[remainingOrderN] += 1
            break

    if ks:
        key = recPart(N, ks, Ms, ys) % N
        print(f"Bob's secret key revealed as: {key}")
        print(f"Altogether this took {time.time() - tim} seconds.")
        return ks
    else:
        print("Something went wrong.")
        print(f"Altogether this took {time.time() - tim} seconds.")
        return None
