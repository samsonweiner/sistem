import numpy as np
from collections import defaultdict

class Graph:
    def __init__(self):
        self.graph = defaultdict(list)
        self.V = set()

    def add_edge(self, u, v):
        self.V.add(u)
        self.V.add(v)
        self.graph[u].append(v)
    
    def get_nodes(self):
        return list(self.V)

    def is_cyclic_util(self, v, visited, rec_stack):
        visited[v] = True
        rec_stack[v] = True

        for neighbor in self.graph[v]:
            if not visited[neighbor]:
                if self.is_cyclic_util(neighbor, visited, rec_stack):
                    return True
            elif rec_stack[neighbor]:
                return True

        rec_stack[v] = False
        return False

    def is_cyclic(self):
        visited = [False] * len(self.V)
        rec_stack = [False] * len(self.V)
        
        for node in range(len(self.V)):
            if not visited[node]:
                if self.is_cyclic_util(node, visited, rec_stack):
                    return True
        return False

def read_graph(graph_path):
    edges = {}
    with open(graph_path) as f:
        for line in f:
            info = line.strip().split('\t')
            edges[info[0]] = int(info[1])
    return edges

def get_pattern(graph):
    pattern = ''
    counts = np.array(list(graph.values()))
    if np.all(counts == 1):
        pattern += 'm'
    else:
        pattern += 'p'

    edges = [tuple(e.split('-')) for e in graph.keys()]
    if np.all(np.array([e[0] for e in edges]) == 'P'):
        pattern += 'PS'
        return pattern
    
    seeded_sites = defaultdict(list)
    for e in edges:
        seeded_sites[e[1]].append(e[0])
    if np.all(np.array([len(v) for v in seeded_sites.values()]) == 1):
        pattern += 'S'
        return pattern
    
    g = Graph()
    for a,b in edges:
        g.add_edge(a, b)
    if g.is_cyclic():
        return pattern
    
    for a,b in edges:
        if (b,a) in edges:
            pattern += 'R'
            return pattern

    pattern += 'M'
    return pattern