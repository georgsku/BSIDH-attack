import time
from avisogenies_sage import *
from sage.rings.all import Integer


load('richelot_aux.sage')
load('helpers.sage')
load('speedup.sage')

# ===================================
# =====  ATTACK  ====================
# ===================================

integer_types = (int, Integer)

def CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i, num_cores):
    tim = time.time()

    finding_u_v_time = 0

    ks = {}
    factorOrder = {}
    recoveredSizeFactor = {}
    for l, e in factor(N):
        ks[l**e] = [-1 for _ in range(e)]
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

        t0 = time.time()

        _, remainingOrderM, remainingOrderN, alpha_i, beta_i, u, v, order,  = find_valid_c(remainingOrderM, remainingOrderN)

        ai = ai - alpha_i
        bi = bi - beta_i
        alp = a - ai
        facOrder = factor(order)

        print(f"Determination of the next {beta_i} digits.")
        
        @possibly_parallel(num_cores)
        def CheckGuess(digits):
            total_guess_time = time.time()
            print(f"Testing digits: {digits}")
            scalar = recPart(recoveredSize, ks, Ms, ys)
            recoveredSizeTmp = recoveredSize
            step = 0
            for l, e in facOrder:
                for i in range(e):
                    m, y = get_CRT_values(l, Ms, ys, recoveredSizeTmp, N)       
                    recoveredSizeTmp *= l
                    if isinstance(digits, integer_types + (Integer,)):
                        scalar += digits * m * y
                    else:
                        scalar += digits[step] * m * y
                    step += 1
                step = e

            tauhatkernel = remainingOrderN * (P3 + scalar*Q3)
            assert tauhatkernel.order() == recoveredSize * order

            AuxiliaryIsogenyTime = time.time()
            C, P_c, Q_c, _ = AuxiliaryIsogeny(beta_i, u, v, E_start, P2, Q2, tauhatkernel, two_i, order * recoveredSize)
            print("\t AuxiliaryIsogenyTime:", time.time() - AuxiliaryIsogenyTime)
            assert P_c.order() == M
            print(C.montgomery_model())
            print(P_c)
            print(Q_c)
            Does22ChainSplitTime = time.time()
            result = Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai)
            print("\t Does22ChainSplitTime:", time.time() - Does22ChainSplitTime)

            print("\t Guessing took:", time.time() - total_guess_time)
            return result

        guesses = create_guesses(order)
        for result in CheckGuess(guesses):
            ((digits,), _), is_split = result
            if is_split and is_split != 'NO DATA':
                print("\n--- ---")
                print("Glue-and-split! These are most likely the secret digits.", digits)
                print("--- ---\n")
                ki_pos = 0
                for l, e in facOrder:
                    for i in range(e):
                        for key in ks:
                            q, r = key.quo_rem(l)
                            if r == 0:
                                if isinstance(digits, integer_types + (Integer,)):
                                    ks[key][recoveredSizeFactor[l] + i] = digits
                                else:
                                    ks[key][recoveredSizeFactor[l] + i] = digits[ki_pos]
                        ki_pos += 1
                break
        else:
            print("Sorry, something wrong happened. Im out of guesses")
            print("However i found these digits:", ks)
            return None
        
        recoveredSize *= order
        for l, e in factor(order):
            recoveredSizeFactor[l] += e

        print("I have found these CRT digits:",ks)
        print("Beginning new round! \n")
        rem = Factors(remainingOrderN)
        if (rem.get_number_of_factors() <= 1):
            break

    print(f"Determination of last #{rem.get_number_of_factors()} CRT digits. We are brute-forcing this.")

    @possibly_parallel(num_cores)
    def CheckGuess(ki):
        last_digit_guess_time = time.time()
        print("Testing final ki:", ki)
        recoveredPart = recPart(recoveredSize, ks, Ms, ys)
        recoveredSizeTmp = recoveredSize
        
        m, y = get_CRT_values(remainingOrderN, Ms, ys, recoveredSizeTmp, N)
        if isinstance(ki, integer_types + (Integer,)):
            scalar = ki * m * y
        else:
            scalar = ki[0] * Ms[remainingOrderN] * ys[remainingOrderN]
        bobskey = (recoveredPart + scalar) % N
        bobscurve, _ = PushingLChain(E_start,  P3 + bobskey*Q3, a, N)
        print("last_digit_guess_time", time.time() - last_digit_guess_time)
        return bobscurve.j_invariant() == EB.j_invariant()
    
    guesses = create_guesses(remainingOrderN)
    for result in CheckGuess(guesses):
        ((digit,), _), is_split = result
        if is_split:
            ki_pos = 0
            for l, e in factor(remainingOrderN):
                for i in range(e):
                    for key in ks:
                        q, r = key.quo_rem(l)
                        if r == 0:
                            ks[key][recoveredSizeFactor[l] + i] = digit
                    ki_pos += 1
            print("Found correct last digit:", digit)
            recoveredSizeFactor[remainingOrderN] += 1
            break
    else:
        print("Sorry, could not find the last digit...")
    
    if ks:
        print(ks)
        key = recPart(N, ks, Ms, ys) % N
        print(f"Bob's secret key revealed as: {key}")
        print(f"Altogether this took {time.time() - tim} seconds.")
        return ks
    else:
        print("Something went wrong.")
        print(f"Altogether this took {time.time() - tim} seconds.")
        return None
