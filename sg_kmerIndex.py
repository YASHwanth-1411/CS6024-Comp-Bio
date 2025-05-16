from collections import defaultdict, deque, Counter
import heapq
import subprocess
import multiprocessing
import os

MIN_OVERLAP = 8
KMER_SIZE = 8

# Step 1: Parse FASTA Reads
def parse_reads_from_fasta(filename):
    reads = []
    with open(filename) as f:
        read = ""
        for line in f:
            if line.startswith('>'):
                if read:
                    reads.append(read.lstrip('N'))
                    read = ""
            else:
                read += line.strip()
        if read:
            reads.append(read.lstrip('N'))
    return reads

# Step 2: Build k-mer index
def build_kmer_index(reads, k):
    index = defaultdict(set)
    for i, read in enumerate(reads):
        for j in range(len(read) - k + 1):
            index[read[j:j + k]].add(i)
    return index

# Step 3: Find overlaps using k-mer index
def find_overlaps(i, reads, kmer_index):
    a = reads[i]
    candidates = set()
    for j in range(len(a) - KMER_SIZE + 1):
        kmer = a[-(j + KMER_SIZE):-j or None]
        candidates.update(kmer_index.get(kmer, set()))
    candidates.discard(i)

    overlaps = {}
    for j in candidates:
        b = reads[j]
        max_len = min(len(a), len(b))
        for k in range(max_len, MIN_OVERLAP - 1, -1):  # Start from longest
            if a[-k:] == b[:k]:
                overlaps[(i, j)] = b[:k]
                break  # Only take the longest valid overlap for this (i, j)
    return overlaps


def build_string_graph(reads):
    kmer_index = build_kmer_index(reads, KMER_SIZE)
    args = [(i, reads, kmer_index) for i in range(len(reads))]

    graph = {}
    with multiprocessing.Pool() as pool:
        results = pool.starmap(find_overlaps, args)
        for res in results:
            graph.update(res)
    return graph

# Step 3: Remove Transitive Edges using DFS
def remove_transitive_edges(graph):
    adj = defaultdict(list)
    for (src, dst), label in graph.items():
        adj[src].append((dst, label))

    non_transitive = graph.copy()

    def dfs(start, current, visited, target):
        visited.add(current)
        for neighbor, _ in adj[current]:
            if neighbor == target and current != start:
                if (start, target) in non_transitive:
                    non_transitive.pop((start, target), None)
                continue
            if neighbor not in visited:
                dfs(start, neighbor, visited, target)

    for (start, target) in list(graph.keys()):
        visited = set()
        dfs(start, start, visited, target)

    return non_transitive

# Step 4: Collapse Chains
def collapse_chains(graph):
    adj = defaultdict(list)
    rev = defaultdict(list)
    for (u, v), label in graph.items():
        adj[u].append((v, label))
        rev[v].append((u, label))

    in_deg = defaultdict(int)
    out_deg = defaultdict(int)
    for u in adj:
        for v, _ in adj[u]:
            in_deg[v] += 1
            out_deg[u] += 1

    new_graph = {}
    visited_nodes = set()

    for u in adj:
        if (in_deg[u] != 1 or out_deg[u] != 1):
            for v, label in adj[u]:
                path = [(u, v)]
                combined_label = label
                current = v
                while in_deg[current] == 1 and out_deg[current] == 1 and current in adj:
                    next_node, next_label = adj[current][0]
                    path.append((current, next_node))
                    combined_label += next_label
                    visited_nodes.add(current)
                    current = next_node
                new_graph[(u, current)] = combined_label

    for (u, v), label in graph.items():
        if u not in visited_nodes and v not in visited_nodes and (u, v) not in new_graph:
            new_graph[(u, v)] = label

    return new_graph

# Step 5: Estimate Edge Weights
def estimate_edge_weights(graph, reads):
    edge_weights = defaultdict(int)
    label_count = Counter()
    for _, label in graph.items():
        label_count[label] += 1
    for edge, label in graph.items():
        edge_weights[edge] = label_count[label]
    return edge_weights

# Step 6: Flow-Based Classification
def classify_edges_by_flow(graph, edge_weights, epsilon=0.01):
    print("No. of edges:", len(graph))
    classification = {}
    for edge in graph:
        classification[edge] = "required"

    for _ in range(5):
        updated = False
        in_flow = defaultdict(float)
        out_flow = defaultdict(float)

        for (u, v), weight in edge_weights.items():
            in_flow[v] += weight
            out_flow[u] += weight

        for node in set(in_flow.keys()).union(set(out_flow.keys())):
            if abs(in_flow[node] - out_flow[node]) > epsilon:
                updated = True

        if not updated:
            break
    for (u, v), w in edge_weights.items():
        if w == 0:
            classification[(u, v)] = "not_required"
        elif w == 1:
            classification[(u, v)] = "unreliable"
        else:
            classification[(u, v)] = "required"
    print("No. of classified edges:", len(classification))
    return classification
def write_edge_classification(classification, output_path):
    with open(output_path, "w") as f:
        for (u, v), status in classification.items():
            f.write(f"Edge ({u}, {v}): {status}\n")
    print(f"Edge classifications written to {output_path}")
# Step 7: Construct Final Genome
def construct_genome(graph, classification, reads):
    adj = defaultdict(list)
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)

    for (u, v), label in graph.items():
        if classification.get((u, v)) == "required":
            adj[u].append((v, label))
            in_deg[v] += 1
            out_deg[u] += 1

    start_nodes = [node for node in adj if in_deg[node] == 0]
    if not start_nodes:
        start_nodes = [next(iter(adj))]

    genome = ""
    visited = set()

    for start in start_nodes:
        current = start
        path = [current]
        while adj[current]:
            next_node, label = adj[current][0]
            if next_node in visited:
                break
            genome += label
            visited.add(current)
            current = next_node
            path.append(current)

    genome = reads[start] + genome
    return genome

# Step 8: QUAST Verification
def run_quast(assembled_genome, reference_genome):
    cmd = [
        "python3", "/mnt/e/desktop/CompBio/quast-5.2.0/quast.py",
        assembled_genome,
        "-r", reference_genome,
        "-o", "quast_output"
    ]
    subprocess.run(cmd)


# ---- MAIN DRIVER ----
def main():
    fasta_file = input("Enter path to FASTA file: ").strip()
    reference_file = input("Enter path to reference genome file (FASTA): ").strip()

    reads = parse_reads_from_fasta(fasta_file)
    print(f"Parsed {len(reads)} reads.")

    graph = build_string_graph(reads)
    print(f"Initial edges: {len(graph)}")

    graph = remove_transitive_edges(graph)
    print(f"After transitive edge removal: {len(graph)} edges")

    graph = collapse_chains(graph)
    print(f"After chain collapsing: {len(graph)} edges")

    edge_weights = estimate_edge_weights(graph, reads)
    classification = classify_edges_by_flow(graph, edge_weights)
    write_edge_classification(classification, "edge_classification.txt")
    genome = construct_genome(graph, classification, reads)
    with open("assembled_genome.fasta", "w") as f:
        f.write(">assembled\n")
        f.write(genome + "\n")
    print("Final genome written to assembled_genome.fasta")

    run_quast("assembled_genome.fasta", reference_file)
    print("QUAST analysis completed. Check 'quast_output' folder.")

if __name__ == "__main__":
    main()
