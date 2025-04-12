import multiprocessing as mp
from collections import defaultdict
from functools import partial

# ------------ Multiprocessing Helpers ------------
def dbg_defaultdict():
    return defaultdict(int)

reads_global = None  # Shared data for string graph

# ------------ FASTA File Reader ------------
def read_fasta(filepath):
    reads = []
    with open(filepath, 'r') as f:
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    reads.append(seq)
                    seq = ""
            else:
                seq += line
        if seq:
            reads.append(seq)
    return reads

# ------------ Parallel De Bruijn Graph ------------
def process_de_bruijn_chunk(chunk, k):
    local_graph = defaultdict(dbg_defaultdict)
    for read in chunk:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            local_graph[prefix][suffix] += 1
    return local_graph

def build_de_bruijn_graph(reads, k):
    if not reads:
        return defaultdict(dbg_defaultdict)
    
    num_processes = mp.cpu_count()
    chunk_size = max(1, len(reads) // num_processes)
    chunks = [reads[i:i+chunk_size] for i in range(0, len(reads), chunk_size)]
    
    with mp.Pool(num_processes) as pool:
        results = pool.starmap(process_de_bruijn_chunk, [(chunk, k) for chunk in chunks])
    
    # Merge results
    graph = defaultdict(dbg_defaultdict)
    for local_graph in results:
        for prefix, suffixes in local_graph.items():
            for suffix, count in suffixes.items():
                graph[prefix][suffix] += count
    return graph

# ------------ Parallel String Graph ------------
def init_pool(reads):
    global reads_global
    reads_global = reads.copy()

def process_row(args):
    i, min_overlap = args
    edges = []
    a = reads_global[i]
    for j in range(len(reads_global)):
        if i == j:
            continue
        b = reads_global[j]
        olen = overlap(a, b, min_overlap)
        if olen > 0:
            edges.append((a, b, olen))
    return edges

def overlap(a, b, min_length=3):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def build_string_graph(reads, min_overlap=3):
    if not reads:
        return defaultdict(list)
    
    n = len(reads)
    with mp.Pool(initializer=init_pool, initargs=(reads,)) as pool:
        results = pool.map(process_row, [(i, min_overlap) for i in range(n)])
    
    graph = defaultdict(list)
    for edges in results:
        for a, b, olen in edges:
            graph[a].append((b, olen))
    return graph

# ------------ Transitive Reduction ------------
def transitive_reduction(graph):
    reduced_graph = defaultdict(list)
    for src in graph:
        direct_dests = set(dest for dest, _ in graph[src])
        indirect_dests = set()
        for intermediate, _ in graph[src]:
            for final, _ in graph.get(intermediate, []):
                indirect_dests.add(final)
        for dest, olen in graph[src]:
            if dest not in indirect_dests:
                reduced_graph[src].append((dest, olen))
    return reduced_graph

# ------------ Main Program ------------
def main():
    fasta_file = input("Enter path to FASTA file: ").strip()
    reads = read_fasta(fasta_file)
    print(f"\nLoaded {len(reads)} reads.")

    # De Bruijn Graph
    k = int(input("Enter k-mer size for De Bruijn Graph (e.g. 4): "))
    dbg = build_de_bruijn_graph(reads, k)
    with open("de_bruijn_graph.txt", "w") as f:
        f.write("--- De Bruijn Graph ---\n")
        for src in dbg:
            for dest in dbg[src]:
                count = dbg[src][dest]
                f.write(f"{src} -> {dest} [count={count}]\n")

    # String Graph
    min_ov = int(input("\nEnter minimum overlap length for String Graph (e.g. 3): "))
    sg = build_string_graph(reads, min_ov)
    sg_reduced = transitive_reduction(sg)
    with open("string_graph.txt", "w") as f:
        f.write("--- String Graph (with Transitive Reduction) ---\n")
        for src in sg_reduced:
            for dest, olen in sg_reduced[src]:
                f.write(f"{src} -> {dest} [overlap={olen}]\n")

if __name__ == "__main__":
    main()
