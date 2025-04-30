import numpy as np
from sklearn.preprocessing import MinMaxScaler

from collections import defaultdict
import os

def get_num_migrations(clone, p):
    nmigrations = np.random.binomial(clone.popsize, p)
    return nmigrations

def draw_migrations(nmigrations, tprobs):
    migration_counts = np.random.multinomial(nmigrations, tprobs / sum(tprobs))
    migrations = (migration_counts > 0).astype(int)
    return migrations

def create_random_distances(nsites):
    points = [np.random.uniform(size=2) for i in range(nsites)]
    dists = []

    for i in range(nsites):
        for j in range(i, nsites):
            dists.append(np.linalg.norm(points[i] - points[j]))
    dists = np.array(dists).reshape(-1, 1)
    norm_dists = MinMaxScaler().fit_transform(dists)

    return points, norm_dists

def compute_diff(A, B):
    '''
    Returns the number of elements unique to A or B. 
    '''
    return len(A.union(B)) - len(A.intersection(B))

def compute_score(A, D):
    score = 0
    for i in range(D.shape[0] - 1):
        for j in range(i, D.shape[0]):
            diff = compute_diff(A[i], A[j])
            score += np.abs(D[i,j] - diff)
    return score

def get_average_dist(A, D, i):
    dists = []
    for j in range(D.shape[0]):
        if i != j:
            diff = compute_diff(A[i], A[j])
            dists.append(np.abs(D[i,j] - diff))
    return np.mean(dists)

def get_total_dist(A, D, i):
    dists = []
    for j in range(D.shape[0]):
        if i != j:
            diff = compute_diff(A[i], A[j])
            dists.append(np.abs(D[i,j] - diff))
    return np.sum(dists)

def find_improvement(A, D, used, unused, score):
    order = [i for i in range(1, D.shape[0])]
    order.sort(key = lambda x: get_average_dist(A, D, x), reverse=True)
    candidates = used[:]
    #print(order)
    if len(unused):
        candidates.append(unused[0])
    for i in order:
        scores, els = [], []
        for c in candidates:
            if c not in A[i]:
                A[i].add(c)
                els.append(c)
                scores.append(compute_score(A, D))
                A[i].remove(c)
                #print(i, c, scores[-1])
        #print(i, scores, els)
        if len(scores):
            idx = np.argmin(scores)
            best = scores[idx]
            if best < score:
                A[i].add(els[idx])
                if len(unused) and els[idx] == unused[0]:
                    unused.pop(0)
                    used.append(els[idx])
                return best
    return score

def assign_n_alt_drivers(D, num_elements):
    nsites = D.shape[0]
    A = {i: set() for i in range(nsites)}
    scores = [compute_score(A, D)]
    used = []
    unused = [e for e in range(num_elements)]
    while len(scores) <= 1 or scores[-1] < scores[-2]:
        next_score = find_improvement(A, D, used, unused, scores[-1])
        #print(next_score, used, unused, A)
        scores.append(next_score)
    return A, scores[-1]

def get_dist_matrix(A):
    if isinstance(A, list):
        alt = {i: set(v) for i,v in enumerate(A)}
    else:
        alt = A
    mat = np.zeros((len(A), len(A)), dtype=int)
    for i in range(len(A)):
        for j in range(i, len(A)):
            if i == j:
                mat[i,i] = 0
            else:
                diff = compute_diff(alt[i], alt[j])
                mat[i,j] = mat[j,i] = diff
    return mat

def generate_migration_graph(tree, site_ids):
    edges = defaultdict(int)
    for node in tree.traverse():
        psite = site_ids[node.site]
        for child in node.children:
            csite = site_ids[child.site]
            if psite != csite:
                edges[f'{psite}-{csite}'] += 1
    return dict(edges)

def save_migration_graph(out_dir, edges):
    with open(os.path.join(out_dir, f'migration_graph.tsv'), 'w+') as f:
        for i,v in edges.items():
            f.write(f'{i}\t{v}\n')