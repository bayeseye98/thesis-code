import pandas as pd
import glob
import os

# BLAST columns (outfmt 6 format)
columns = ["query", "subject", "percent_identity", "alignment_length", "mismatches", "gap_opens",
           "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"]

# Find all files vsFungi_blastp.out
fungi_files = glob.glob("*vsFungi_blastp.out")
sample_names = [f.replace("vsFungi_blastp.out", "") for f in fungi_files]

for sample in sample_names:
    try:
        # Build corresponding file paths
        fungi_file = f"{sample}vsFungi_blastp.out"
        saccharo_file = f"{sample}vsSaccharomycotina_blastp.out"
        self_file = f"{sample}vsItself_blastp.out"

        # Check that all required files exist
        if not all(os.path.exists(f) for f in [fungi_file, saccharo_file, self_file]):
            print(f"⚠️ Missing files for {sample}, skipping.")
            continue

        # Load files
        fungi_blast = pd.read_csv(fungi_file, sep="\t", header=None, names=columns)
        saccharo_blast = pd.read_csv(saccharo_file, sep="\t", header=None, names=columns)
        self_blast = pd.read_csv(self_file, sep="\t", header=None, names=columns)

        # Get all queries
        all_queries = set(fungi_blast["query"]).union(saccharo_blast["query"], self_blast["query"])

        # Compute best bitscores
        fungi_scores = fungi_blast.groupby("query")["bitscore"].max().reindex(all_queries, fill_value=0)
        saccharo_scores = saccharo_blast.groupby("query")["bitscore"].max().reindex(all_queries, fill_value=0)
        self_scores = self_blast.groupby("query")["bitscore"].max().reindex(all_queries, fill_value=0)

        # Normalization
        self_scores_safe = self_scores.replace(0, 1)  # avoid division by zero
        norm_fungi = fungi_scores / self_scores_safe
        norm_saccharo = saccharo_scores / self_scores_safe

        # Create DataFrame
        results = pd.DataFrame({
            "query": list(all_queries),
            "fungi_bitscore": norm_fungi,
            "saccharomycotina_bitscore": norm_saccharo,
            "self_bitscore": self_scores
        })

        # Add subjects
        results["fungi_subject"] = fungi_blast.groupby("query")["subject"].first().reindex(all_queries, fill_value="No Subject")
        results["saccharomycotina_subject"] = saccharo_blast.groupby("query")["subject"].first().reindex(all_queries, fill_value="No Subject")

        # AI calculation
        results["AI"] = norm_fungi - norm_saccharo

        # Save output
        output_file = f"{sample}AlienIndex_results.csv"
        results.to_csv(output_file, index=False)
        print(f"✔️ {output_file} created.")

    except Exception as e:
        print(f"❌ Error for {sample}: {e}")
