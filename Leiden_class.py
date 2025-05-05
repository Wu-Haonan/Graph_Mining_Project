import edlib
from Bio import SeqIO
import itertools
import pandas as pd
import igraph as ig
import leidenalg

def read_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, 'fasta'):
        id = record.id.split("|")[0].split(":")[1]
        sequences[id] = str(record.seq)
    return sequences

def build_graph_with_edlib(sequences, k):
    ids = list(sequences.keys())
    g = ig.Graph()
    g.add_vertices(ids)

    for (id1, seq1), (id2, seq2) in itertools.combinations(sequences.items(), 2):
        result = edlib.align(seq1, seq2, mode="NW", task="distance", k=k)
        dist = result['editDistance']
        # print(id1, id2, dist)
        if dist != -1:
            g.add_edge(id1, id2, weight=1/(1+dist))

    return g

def run_leiden(graph):
    partition = leidenalg.find_partition(graph, leidenalg.ModularityVertexPartition, weights="weight")
    return partition

if __name__ == '__main__':
    fasta_file = './data/PRJNA324093_Dnr4_10k.fasta'  
    max_distance = 40 

    sequences = read_fasta(fasta_file)
    graph = build_graph_with_edlib(sequences, k=max_distance)
    partition = run_leiden(graph)

    cluster_data = []
    for vid, cluster in zip(graph.vs['name'], partition.membership):
        print(f"{vid}\tCluster {cluster}")
        cluster_data.append((vid, cluster))

    df = pd.DataFrame(cluster_data, columns=["Sequence_ID", "Cluster"])
    df.to_csv("leiden_clusters.csv", index=False)

