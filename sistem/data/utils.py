import numpy as np
from scipy.stats import beta
from scipy.optimize import newton_krylov, NoConvergence

from sistem.utilities.utilities import get_reg_id

# Estimating alpha/beta from point on lorenz curve
def get_alpha_beta(x0, y0):
    def F(P):
        Aa, Bb = P[0], P[1]
        X = beta.cdf(float(Aa)/(Aa+Bb), Aa, Bb) - x0
        Y = beta.cdf(float(Aa)/(Aa+Bb), Aa+1, Bb) - y0
        return [X, Y]
    guess = [10, 5]
    max_n = 1000
    n = 1
    while n < max_n:
        try:
            sol = newton_krylov(F, guess, method = 'lgmres', verbose = 0, rdiff = 0.1, maxiter=50)
            break
        except NoConvergence as e:
            guess = np.random.rand(2) * 10 + 0.1
        except ValueError as e:
            guess = np.random.rand(2) * 10 + 0.1
        n += 1
    if n == max_n:
        # Error has occurred here
        return [1.38, 1.38]
    else:
        return sol

def bezier_coef(points):
    n = len(points) - 1

    C = 4 * np.identity(n)
    np.fill_diagonal(C[1:], 1)
    np.fill_diagonal(C[:, 1:], 1)
    C[0, 0] = 2
    C[n - 1, n - 1] = 7
    C[n - 1, n - 2] = 2

    P = [2 * (2 * points[i] + points[i + 1]) for i in range(n)]
    P[0] = points[0] + 2 * points[1]
    P[n - 1] = 8 * points[n - 1] + points[n]

    A = np.linalg.solve(C, P)
    B = [0] * n
    for i in range(n - 1):
        B[i] = 2 * points[i + 1] - A[i + 1]
    B[n - 1] = (A[n - 1] + points[n]) / 2

    return A, B

def single_cubic_bezier(a, b, c, d):
    return lambda t: np.power(1 - t, 3) * a + 3 * np.power(1 - t, 2) * t * b + 3 * (1 - t) * np.power(t, 2) * c + np.power(t, 3) * d

def multiple_cubic_bezier(points):
    A, B = bezier_coef(points)
    return [single_cubic_bezier(points[i], A[i], B[i], points[i + 1]) for i in range(len(points) - 1)]
    
# Generates starting coverage points at set intervals 
def gen_start_interval(num_windows, interval, Aa, Bb):
    # Sample point every interval bin
    x = [i for i in range(0, num_windows, interval)]
    if x[-1] != num_windows-1:
        if num_windows - x[-1] >= interval/2:
            x.append(num_windows-1)
        else:
            x[-1] = num_windows-1
    # draw coverage from beta distribution
    y = list(np.random.beta(Aa, Bb, len(x)))
    points = [np.array(p) for p in zip(x, y)]
    return points

# Smooths points across bins with bezier curves
def gen_coverage(num_windows, interval, Aa, Bb):
    points = gen_start_interval(num_windows, interval, Aa, Bb)
    A, B = bezier_coef(points)
    curves = [single_cubic_bezier(points[i], A[i], B[i], points[i + 1]) for i in range(len(points) - 1)]

    new_points = []
    for i in range(len(points) - 1):
        f = curves[i]
        gaps = points[i+1][0] - points[i][0] + 1
        coords = [f(t) for t in np.linspace(0, 1, int(gaps))]
        if i == 0:
            new_points.append((round(coords[0][0]), coords[0][1]))
        new_points += [(round(i[0]), i[1]) for i in coords[1:]]
    
    new_points.sort(key=lambda x: x[0])
    # Multiply by 2 so that values are in range 0-2. Value of 1 indicates mean coverage.
    new_points = [2*max(min(i[1], 1), 0) for i in new_points]
    return new_points

def bin_region_values(vals, r2b):
    """
    Given an array of values and the number of regions in each bin, returns a new array containing the binwise averages. Allows the ratio to be non-integral.

    Arguments:
        vals (array[int]): Value of arrays (e.g. Copy numbers, readcounts.)
        r2b (int, float): The number of regions per bin (e.g. 2, 2.5).
    """
    res = []
    tot_sum = 0
    n = len(vals)
    while tot_sum < n:
        idx = int(tot_sum)

        left_w = np.ceil(tot_sum) - tot_sum
        if left_w == 0:
            left_w = 1
        
        binned_vals = [vals[idx]]
        weights = [left_w]

        if tot_sum + r2b > n:
            for i in range(idx + 1, n):
                binned_vals.append(vals[i])
                weights.append(1)
            
        else:
            full_w = int(r2b - left_w)
            right_w = (r2b - left_w) % 1

            binned_vals.extend([vals[idx + i + 1] for i in range(full_w)])
            
            weights.extend([1 for i in range(full_w)])
            if right_w > 0:
                binned_vals.append(vals[idx + full_w + 1])
                weights.append(right_w)

        res.append(float(np.average(binned_vals, weights=weights)))
        tot_sum += r2b

    return res

def bin_region_readcounts(vals, r2b):
    """
    Given an array of readcounts and the number of regions in each bin, returns a new array containing the binwise sums. Allows the ratio to be non-integral. Normalizez for the number of regions.

    Arguments:
        vals (array[int]): Value of arrays (e.g. Copy numbers, readcounts.)
        r2b (int, float): The number of regions per bin (e.g. 2, 2.5).
    """
    res = []
    tot_sum = 0
    n = len(vals)
    while tot_sum < n:
        idx = int(tot_sum)

        left_w = np.ceil(tot_sum) - tot_sum
        if left_w == 0:
            left_w = 1
        
        binned_vals = [vals[idx]]
        weights = [left_w]

        if tot_sum + r2b > n:
            for i in range(idx + 1, n):
                binned_vals.append(vals[i])
                weights.append(1)
            
        else:
            full_w = int(r2b - left_w)
            right_w = (r2b - left_w) % 1

            binned_vals.extend([vals[idx + i + 1] for i in range(full_w)])
            
            weights.extend([1 for i in range(full_w)])
            if right_w > 0:
                binned_vals.append(vals[idx + full_w + 1])
                weights.append(right_w)

        next_val = 0
        for i,x in enumerate(binned_vals):
            next_val += x*weights[i]
        res.append(round(next_val * (min(r2b, n - tot_sum)) / r2b))
        tot_sum += r2b

    return res

def count_region_CNs(cell, total=False):
    profile = {}
    for chrname, nregions in cell.library.regions.items():
        CNs = [np.zeros(nregions, dtype=int), np.zeros(nregions, dtype=int)]
        for chromosome in cell.genome[chrname]:
            counts = np.bincount([get_reg_id(x) for x in chromosome.seq], minlength=nregions)
            CNs[chromosome.allele] += counts
        
        if total:
            profile[chrname] = CNs[0] + CNs[1]
        else:
            profile[chrname] = CNs
    return profile


def get_mutated_basepairs(observed_cells, chrnames, region_len):
    mutated_bps = {chrname: set() for chrname in chrnames}
    for chrname in chrnames:
        for cell in observed_cells:
            for chromosome in cell.genome[chrname]:
                for r, positions in chromosome.SNVs.items():
                    v = get_reg_id(int(r))
                    for p in positions:
                        loc = int(v*region_len + p)
                        mutated_bps[chrname].add(loc)
        mutated_bps[chrname] = list(mutated_bps[chrname])
        mutated_bps[chrname].sort()
    mutated_bps = [(chrname, pos) for chrname in chrnames for pos in mutated_bps[chrname]]
    return mutated_bps