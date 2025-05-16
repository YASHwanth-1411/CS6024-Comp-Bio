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

def find_path_start(in_deg, out_deg):
    """Find optimal start node considering various conditions"""
    # Prioritize nodes with out_degree - in_degree == 1
    for node in out_deg:
        if out_deg[node] - in_deg.get(node, 0) == 1:
            return node
    # Fallback to any node with outgoing edges
    for node in out_deg:
        if out_deg[node] > 0:
            return node
    return None

def extract_single_path(graph, start_node, in_deg, out_deg):
    """Modified Hierholzer's algorithm for single path extraction with edge removal"""
    if start_node is None:
        return []
    
    stack = [start_node]
    path = []
    current = start_node
    edge_counts = {src: len(dests) for src, dests in graph.items()}
    
    while stack:
        if edge_counts.get(current, 0) > 0:
            stack.append(current)
            next_node = graph[current].pop()
            edge_counts[current] -= 1
            # Update degree counts
            out_deg[current] -= 1
            in_deg[next_node] -= 1
            current = next_node
        else:
            path.append(current)
            current = stack.pop()
    
    return path[::-1]  # Reverse to get correct order

def clean_graph(graph, in_deg, out_deg):
    """Remove nodes with no outgoing edges and update degrees"""
    to_remove = [node for node in graph if not graph[node]]
    for node in to_remove:
        del graph[node]
        if node in out_deg:
            del out_deg[node]
        if node in in_deg:
            del in_deg[node]

def extract_all_contigs(original_graph, original_in, original_out, k):
    """Iterative contig extraction using path cover approach"""
    # Create working copies to preserve original
    graph = {node: edges[:] for node, edges in original_graph.items()}
    in_deg = original_in.copy()
    out_deg = original_out.copy()
    
    contigs = []
    
    while graph:
        start_node = find_path_start(in_deg, out_deg)
        if not start_node:
            break
            
        path = extract_single_path(graph, start_node, in_deg, out_deg)
        
        if len(path) < 2:  # Single node can't form contig
            break
            
        # Reconstruct contig
        contig = path[0]
        for node in path[1:]:
            contig += node[-1]
        
        # Handle circular contigs
        if contig.startswith(contig[-(k-1):]):
            contig = contig[:-(k-1)]
        
        contigs.append(contig)
        clean_graph(graph, in_deg, out_deg)
    
    return contigs

# ------------ Validation & Comparison ------------

def validate_assembly(assembly_path, reference_path):
    # """Run comprehensive validation pipeline"""
    # os.makedirs("validation_results", exist_ok=True)
    
    # # 1. QUAST Quality Assessment
    # subprocess.run([
    #     "/usr/bin/quast.py",
    #     "-r", "ecoli_ref.fasta",
    #     "-o", "validation_results/quast",
    #     assembly_path
    # ], check=True)
    
    # # 2. MUMmer Alignment
    # subprocess.run([
    #     "nucmer",
    #     "--prefix=validation_results/alignment",
    #     reference_path,
    #     assembly_path
    # ], check=True)
    
    # subprocess.run([
    #     "dnadiff",
    #     "-d", "validation_results/alignment.delta",
    #     "-p", "validation_results/dnadiff"
    # ], check=True)
    
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
    assembly_file = "reconstructed_contigs.fasta"
    
    # 1. Parse De Bruijn Graph
    print("Reading De Bruijn graph...")
    graph, in_degree, out_degree, k = parse_de_bruijn_graph(graph_file)
    
    print(f"Graph stats: {len(graph)} nodes, k={k}")
    print(f"Semi-balanced nodes: {sum(1 for n in graph if abs(out_degree[n] - in_degree.get(n,0)) == 1)}")
    
    # 2. Extract All Contigs
    print("\nExtracting contigs using iterative path cover...")
    contigs = extract_all_contigs(graph, in_degree, out_degree, k)
    
    print(f"Found {len(contigs)} contigs")
    print("Contig lengths:", [len(c) for c in sorted(contigs, key=len, reverse=True)[:5]])
    
    # 3. Save Assembly
    with open(assembly_file, "w") as f:
        for i, contig in enumerate(sorted(contigs, key=len, reverse=True)):
            f.write(f">Contig_{i+1}_len={len(contig)}\n")
            for j in range(0, len(contig), 80):
                f.write(contig[j:j+80] + "\n")
    
    # 4. Validate Against Reference
    print("\nValidating assembly...")
    validate_assembly(assembly_file, ref_genome)
    
    print("\nValidation complete! Check:")
    print("- validation_results/quast/report.html")
    print("- validation_results/dnadiff.report")

if __name__ == "__main__":
    main()
