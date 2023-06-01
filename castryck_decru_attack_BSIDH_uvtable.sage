import time
from avisogenies_sage import *
from sage.rings.all import Integer


load('richelot_aux.sage')
load('helpers.sage')
load('speedup.sage')
load('uvtable_example_2.sage')

# ===================================
# =====  ATTACK  ====================
# ===================================

integer_types = (int, Integer)
guessing_times = []

def CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i, num_cores):
    tim = time.time()

    ks = {}
    factorOrder = {}
    recoveredSizeFactor = {}
    for l, e in factor(N):
        ks[l**e] = [-1 for _ in range(e)]
        factorOrder[l] = e
        recoveredSizeFactor[l] = 0

    Ms, ys = genCRTconstants(N)

    recoveredSize = 1
    remainingOrderN = N
    ai = a
    alp = M
    uv_step = 0
    guessing_number_step = 0
    
    while uv_step < len(uvtable):
        print("Remaining order:", factor(remainingOrderN))

        ai, order, u, v = uvtable[uv_step]
        remainingOrderN = remainingOrderN/order
        beta_i = len(factor(order))

        alp = a - ai
        facOrder = factor(order)

        print("From uvtable:")
        print("u:", u, "v:", v)
        print("revovering:", facOrder)

        tuple  = find_valid_c(2^ai, remainingOrderN*order)
        print(find_u_v(2^ai - remainingOrderN))
        print(tuple)

        print(f"Determination of the next {beta_i} digits. We are working with 2^{ai}", "torsion")
        
        @possibly_parallel(num_cores)
        def CheckGuess(digits):
            total_guess_time = time.time()
            print(f"Testing digits: {digits}")
            scalar = recPart(recoveredSize, ks, Ms, ys)
            recoveredSizeTmp = recoveredSize
            step = 0
            for l, e in facOrder:   
                for _ in range(e):
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
            Does22ChainSplitTime = time.time()
            result = Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai)
            print("\t Does22ChainSplitTime:", time.time() - Does22ChainSplitTime)
            tot_guess_time =  time.time() - total_guess_time
            guessing_times.append([tot_guess_time, order])
            print("\t Guessing took:", time.time() - total_guess_time)
            return result

        guesses = create_guesses(order)
        for result in CheckGuess(guesses):
            ((digits,), _), is_split = result
            print(result)
            if is_split and is_split != 'NO DATA':
                print("\n--- ---")
                print("Glue-and-split! These are most likely the secret digits.", digits)
                print("--- ---\n")
                ki_pos = 0
                guessing_number_step += 1
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
        uv_step += 1

    print(f"Determination of last 1 CRT digits. We are brute-forcing this.")

    @possibly_parallel(num_cores)
    def CheckGuess(ki):
        last_digit_guess_time = time.time()
        print("Testing final ki:", ki)
        recoveredPart = recPart(recoveredSize, ks, Ms, ys)
        recoveredSizeTmp = recoveredSize
        
        #m, y = get_CRT_values(remainingOrderN, Ms, ys, recoveredSizeTmp, N)
        scalar = ki * Ms[remainingOrderN] * ys[remainingOrderN]
        bobskey = (recoveredPart + scalar) % N
        bobscurve, _ = PushingLChain(E_start,  P3 + bobskey*Q3, a, N)
        print(bobscurve)
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
        key = recPart(N, ks, Ms, ys) % N
        print(f"Bob's secret key revealed as: {key}")
        print(f"Altogether this took {time.time() - tim} seconds.")
        print(f"Individual guessing times where:", guessing_times)
        return ks
    else:
        print("Something went wrong.")
        print(f"Altogether this took {time.time() - tim} seconds.")
        print(f"Individual guessing times where:", guessing_times)
        return None
