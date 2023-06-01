load("helpers.sage")

def find_lowest_alpha(N, M):
        for i in range(factor(N)[0][1] + 1):
            C = 2^i - M
            if C % 4 and C > 1:
                u_v = find_u_v(C)
                if u_v is not None:
                    return i, u_v[0], u_v[1]
        return None

def prime_combinations(int):
    # Get the prime factors
    prime_factors = [p for p, e in list(factor(int)) for _ in range(e)]

    # Sort in descending order
    prime_factors.sort(reverse=True)
    
    # Copy the prime factors
    prime_factors_combinations = list(prime_factors)

    result = []
    while prime_factors_combinations:
        # Pop the largest prime factor
        largest_prime_factor = prime_factors_combinations.pop(0)
        # Generate all combinations of the largest prime factor and each of the other prime factors, 
        # starting with the smallest, and append them to the result
        result.extend(largest_prime_factor * p for p in reversed(prime_factors_combinations))

    # Concatenate the prime factors and the combinations and return
    return prime_factors + result


# Function to find optimal sequence
# Returns a list list of form [alpha_i, factor to remove, u, v]
def optimal_sequence(N, M):
    # Initialize variables: temporary M, recovered M and results 
    M_tmp = M
    recovered_M = 1
    results = [[-1,-1,-1,-1]] 

    # Iterate while there are at least one prime factor of M_tmp
    while len([p for p, e in list(factor(M_tmp)) for _ in range(e)]) > 1:
        # Iterate over each prime factor combination of M_tmp
        for prime_factor in prime_combinations(M_tmp):
            # Multiply the current prime factor to the recovered_M
            recovered_M *= prime_factor

            # If recovered_M modulo 4 is not equal to 3
            # then divide the recovered_M by the prime factor and continue the loop
            if recovered_M % 4 != 3: 
                recovered_M /= prime_factor
                continue

            # Call the function find_lowest_alpha with parameters N and recovered_M
            # If the result is None, then continue the loop
            result = find_lowest_alpha(N, recovered_M)
            if result is None: continue

            # If there is a result, unpack the result into N_low, u, and v
            N_low, u, v = result
            # Divide the M_tmp by the prime factor
            M_tmp /= prime_factor
            # Append the result to the results list
            results[-1][1] = prime_factor
            results.append([N_low, -1, u, v])
            # Need to update factor of prevous list such that u_v list is correct
            # Break the loop
            break
    # Return the results list
    results[-1][1] = M_tmp
    results.reverse()
    return results