from collections import defaultdict

# ------------ FASTA File Reader (Pure Python) ------------
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

# ------------ De Bruijn Graph Construction ------------
def build_de_bruijn_graph(reads, k):
    graph = defaultdict(list)
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            graph[prefix].append(suffix)
    return graph

# ------------ String Graph Construction ------------
def overlap(a, b, min_length=3):
    """Return length of longest suffix of 'a' matching prefix of 'b'."""
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def build_string_graph(reads, min_overlap=3):
    graph = defaultdict(list)
    for a in reads:
        for b in reads:
            if a != b:
                olen = overlap(a, b, min_overlap)
                if olen > 0:
                    graph[a].append((b, olen))
    return graph

# ------------ Transitive Reduction for String Graph ------------
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

# ------------ Main Runner ------------
def main():
    fasta_file = input("Enter path to FASTA file: ").strip()
    reads = read_fasta(fasta_file)
    print(f"\nLoaded {len(reads)} reads.")

    # --- De Bruijn Graph ---
    k = int(input("Enter k-mer size for De Bruijn Graph (e.g. 4): "))
    dbg = build_de_bruijn_graph(reads, k)
    print("\n--- De Bruijn Graph ---")
    for node in dbg:
        print(f"{node} -> {', '.join(dbg[node])}")

    # --- String Graph (with Transitive Reduction) ---
    min_ov = int(input("\nEnter minimum overlap length for String Graph (e.g. 3): "))
    sg = build_string_graph(reads, min_ov)
    sg_reduced = transitive_reduction(sg)

    print("\n--- String Graph (with Transitive Reduction) ---")
    for src in sg_reduced:
        for dest, olen in sg_reduced[src]:
            print(f"{src} -> {dest} [overlap={olen}]")

if __name__ == "__main__":
    main()
