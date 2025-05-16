from collections import defaultdict, deque, Counter
import heapq

# Parameters
MIN_OVERLAP = 7

# Step 1: Parse FASTA Reads
def parse_reads_from_fasta(filename):
    reads = []
    with open(filename) as f:
        read = ""
        for line in f:
            if line.startswith('>'):
                if read:
                    reads.append(read.lstrip('N'))  # Remove leading Ns
                    read = ""
            else:
                read += line.strip()
        if read:
            reads.append(read.lstrip('N'))
    return reads

# Step 2: Burrows-Wheeler Transform and FM-index

def build_bwt(text):
    text += "$"
    rotations = [text[i:] + text[:i] for i in range(len(text))]
    rotations.sort()
    bwt = ''.join([row[-1] for row in rotations])
    return bwt

def build_suffix_array(text):
    text += "$"
    return sorted(range(len(text)), key=lambda i: text[i:])

def build_fm_index(bwt):
    ranks = defaultdict(list)
    totals = defaultdict(int)
    for c in "ACGT$":
        ranks[c].append(0)  # index 0 for before any character
    for i, char in enumerate(bwt):
        for c in "ACGT$":
            ranks[c].append(ranks[c][-1] + (1 if char == c else 0))
        totals[char] += 1
    return ranks, totals


def count_pattern(pattern, bwt, ranks, totals):
    top = 0
    bottom = len(bwt) - 1
    while top <= bottom and pattern:
        symbol = pattern[-1]
        pattern = pattern[:-1]
        top = totals[symbol] + ranks[symbol][top]
        bottom = totals[symbol] + ranks[symbol][bottom + 1] - 1
    return bottom - top + 1 if top <= bottom else 0

# Step 3: Find maximum suffix-prefix overlap
def find_overlap(a, b, min_length=MIN_OVERLAP):
    max_len = 0
    for i in range(min_length, min(len(a), len(b)) + 1):
        if a[-i:] == b[:i]:
            max_len = i
    return max_len

# Step 4: Build Overlap Graph using BWT + FM-index
def build_string_graph(reads):
    concat_reads = '#'.join(reads) + '#'
    bwt = build_bwt(concat_reads)
    sa = build_suffix_array(concat_reads)
    ranks, totals = build_fm_index(bwt)

    graph = {}
    for i, a in enumerate(reads):
        suffix = a[-MIN_OVERLAP:]
        if count_pattern(suffix, bwt, ranks, totals) == 0:
            continue
        for j, b in enumerate(reads):
            if i != j:
                olap_len = find_overlap(a, b)
                if olap_len > 0:
                    graph[(i, j)] = reads[j][:olap_len]
    return graph

# Step 5: Remove Transitive Edges using DFS
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

# Step 6: Collapse chains
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

# Step 7: Estimate Edge Weights Based on Read Support
def estimate_edge_weights(graph, reads):
    edge_weights = defaultdict(int)
    label_count = Counter()
    for _, label in graph.items():
        label_count[label] += 1
    for edge, label in graph.items():
        edge_weights[edge] = label_count[label]
    return edge_weights

# Step 8: Simulate Min-Cost Flow and Edge Classification
def classify_edges_by_flow(graph, edge_weights, epsilon=0.01):
    classification = {}
    for edge in graph:
        classification[edge] = "required"  # initialize all required

    for _ in range(5):
        updated = False
        in_flow = defaultdict(float)
        out_flow = defaultdict(float)

        for (u, v), weight in edge_weights.items():
            in_flow[v] += weight
            out_flow[u] += weight

        for node in set(in_flow.keys()).union(set(out_flow.keys())):
            in_wt = in_flow[node]
            out_wt = out_flow[node]
            if abs(in_wt - out_wt) > epsilon:
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

    return classification

# Step 9: Output Graph
def write_graph(graph, reads, output_path):
    with open(output_path, "w") as f:
        for (src, dst), label in graph.items():
            f.write(f"R{src} -> R{dst} [label={label}]\n")
    print(f"Graph written to {output_path}")

# Step 10: Write classified edges
def write_edge_classification(classification, output_path):
    with open(output_path, "w") as f:
        for (u, v), status in classification.items():
            f.write(f"Edge ({u}, {v}): {status}\n")
    print(f"Edge classifications written to {output_path}")

# Step 11: Construct Final Genome
def construct_genome(graph, classification, reads):
    adj = defaultdict(list)
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)

    for (u, v), label in graph.items():
        if classification.get((u, v)) == "required":
            adj[u].append((v, label))
            in_deg[v] += 1
            out_deg[u] += 1

    # Find start node (in_deg = 0)
    start_nodes = [node for node in adj if in_deg[node] == 0]
    if not start_nodes:
        start_nodes = [next(iter(adj))]  # fallback

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

    # Prepend the first read from the start node
    genome = reads[start] + genome
    return genome


# ---- MAIN DRIVER ----
def main():
    fasta_file = input("Enter path to FASTA file: ").strip()
    reads = parse_reads_from_fasta(fasta_file)
    print(f"Parsed {len(reads)} reads.")

    graph = build_string_graph(reads)
    print(f"Initial edges: {len(graph)}")

    graph = remove_transitive_edges(graph)
    print(f"After transitive edge removal: {len(graph)} edges")

    graph = collapse_chains(graph)
    print(f"After chain collapsing: {len(graph)} edges")

    write_graph(graph, reads, "string_graph.txt")

    edge_weights = estimate_edge_weights(graph, reads)
    classification = classify_edges_by_flow(graph, edge_weights)

    write_edge_classification(classification, "edge_classification.txt")

    genome = construct_genome(graph, classification, reads)
    with open("assembled_genome.txt", "w") as f:
        f.write(genome + "\n")
    print("Final genome written to assembled_genome.txt")


if __name__ == "__main__":
    main()
