import os
import re
from collections import defaultdict
from Bio import SeqIO
import subprocess

# ------------ Graph Reading & Genome Reconstruction ------------

def parse_de_bruijn_graph(file_path):
    """Read De Bruijn graph from text file with error handling"""
    graph = defaultdict(list)
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    k = None
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('---'): continue
                match = re.match(r"(\w+)\s->\s(\w+)\s\[count=(\d+)\]", line.strip())
                if match:
                    src, dest, count = match.groups()
                    count = int(count)
                    if k is None: 
                        k = len(src) + 1  # k = node_length + 1
                    
                    # Store edges with multiplicity
                    graph[src].extend([dest]*count)
                    out_degree[src] += count
                    in_degree[dest] += count
                    
        return graph, in_degree, out_degree, k
    
    except FileNotFoundError:
        raise SystemExit(f"Error: De Bruijn graph file not found at {file_path}")
    except Exception as e:
        raise SystemExit(f"Error parsing graph file: {str(e)}")

def find_eulerian_start(in_degree, out_degree):
    """Find optimal start node considering circular genomes"""
    start = None
    for node in out_degree:
        # Check for potential start nodes
        diff = out_degree[node] - in_degree.get(node, 0)
        if diff == 1:
            return node  # Classic Eulerian start
        elif diff == -1:
            continue  # Would be end node
        elif out_degree[node] > 0:
            start = node  # Candidate for circular genome
    
    return start or next(iter(out_degree))  # Fallback to first node

def hierholzer_algorithm(graph, start_node):
    """Iterative Hierholzer's algorithm with edge removal"""
    stack = [start_node]
    path = []
    edge_counts = {src: len(dests) for src, dests in graph.items()}
    
    current = start_node
    while stack:
        if edge_counts.get(current, 0) > 0:
            stack.append(current)
            next_node = graph[current].pop()
            edge_counts[current] -= 1
            current = next_node
        else:
            path.append(current)
            current = stack.pop()
    
    return path[::-1]  # Reverse to get correct order

def reconstruct_genome(path, k):
    """Build genome sequence from Eulerian path"""
    if not path:
        return ""
    
    genome = [path[0]]
    for node in path[1:]:
        genome.append(node[-1])
    
    # Circular genome handling
    if len(path) > 1 and genome[0] == genome[-1][-(k-1):]:
        return ''.join(genome)[:-(k-1)]
    
    return ''.join(genome)

# ------------ Validation & Comparison ------------

def validate_assembly(assembly_path, reference_path):
    """Run comprehensive validation pipeline"""
    os.makedirs("validation_results", exist_ok=True)
    
    # 1. QUAST Quality Assessment
    # subprocess.run([
    #     "quast.py",
    #     "-r", reference_path,
    #     "-o", "validation_results/quast",
    #     assembly_path
    # ], check=True)
    
    # 2. MUMmer Alignment
    subprocess.run([
        "nucmer",
        "--prefix=validation_results/alignment",
        reference_path,
        assembly_path
    ], check=True)
    
    subprocess.run([
        "dnadiff",
        "-d", "validation_results/alignment.delta",
        "-p", "validation_results/dnadiff"
    ], check=True)
    
    # 3. Basic Metrics
    ref_length = sum(len(r) for r in SeqIO.parse(reference_path, "fasta"))
    asm_length = sum(len(r) for r in SeqIO.parse(assembly_path, "fasta"))
    
    print(f"\nAssembly Length: {asm_length:,}bp")
    print(f"Reference Length: {ref_length:,}bp")
    print(f"Length Difference: {abs(asm_length - ref_length):,}bp")

# ------------ Main Workflow ------------

def main():
    # Configuration
    graph_file = "de_bruijn_graph.txt"
    ref_genome = "Datasets/EColi/ref_genome/ecoli_k12_reference.fasta"
    assembly_file = "reconstructed_genome.fasta"
    
    # 1. Reconstruct Genome
    print("Reading De Bruijn graph...")
    graph, in_degree, out_degree, k = parse_de_bruijn_graph(graph_file)
    
    print(f"Graph stats: {len(graph)} nodes, k={k}")
    print(f"Semi-balanced nodes: {sum(1 for n in graph if abs(out_degree[n] - in_degree.get(n,0)) == 1)}")
    
    start_node = find_eulerian_start(in_degree, out_degree)
    print(f"Starting reconstruction from node: {start_node}")
    
    eulerian_path = hierholzer_algorithm(graph, start_node)
    genome = reconstruct_genome(eulerian_path, k)
    
    # Save assembly
    with open(assembly_file, "w") as f:
        f.write(f">Reconstructed_EColi\n")
        for i in range(0, len(genome), 80):
            f.write(genome[i:i+80] + "\n")
    
    # 2. Validate Against Reference
    print("\nValidating assembly...")
    validate_assembly(assembly_file, ref_genome)
    
    print("\nValidation complete! Check:")
    print("- validation_results/quast/report.html")
    print("- validation_results/dnadiff.report")

if __name__ == "__main__":
    main()
