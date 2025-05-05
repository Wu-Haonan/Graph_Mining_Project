import pandas as pd
from sklearn.metrics import adjusted_rand_score

def evaluate_clustering(predicted_csv, ground_truth_csv):
    pred_df = pd.read_csv(predicted_csv)
    pred_dict = dict(zip(pred_df['Sequence_ID'], pred_df['Cluster']))

    truth_df = pd.read_csv(ground_truth_csv)
    truth_dict = dict(zip(truth_df['Sequence_ID'], truth_df['Cluster']))


    common_ids = list(set(pred_dict.keys()) & set(truth_dict.keys()))
    pred_labels = [pred_dict[seq_id] for seq_id in common_ids]
    truth_labels = [truth_dict[seq_id] for seq_id in common_ids]

    print(f"Ground truth cluster count: {len(set(truth_labels))}")
    print(f"Leiden cluster count: {len(set(pred_labels))}")

    ari = adjusted_rand_score(truth_labels, pred_labels)
    return ari

if __name__ == '__main__':
    predicted_file = "leiden_clusters.csv"  
    ground_truth_file = "ground_truth.csv" 

    ari_score = evaluate_clustering(predicted_file, ground_truth_file)
    print(f"Adjusted Rand Index (ARI): {ari_score:.4f}")
