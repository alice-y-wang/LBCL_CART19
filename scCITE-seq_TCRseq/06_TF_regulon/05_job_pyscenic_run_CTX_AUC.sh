#!/bin/bash
#SBATCH -c 16
#SBATCH --job-name=Pyscenic_CTX_AUC_Lymphoma
#SBATCH --mem=164G
#SBATCH --mail-type=END,FAIL,BEGIN          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -t 5-10:00:00
#SBATCH --output=scenic-CTX-%j.out

source ~/mambaforge/etc/profile.d/conda.sh
conda activate pyscenic_new
unset PYTHONPATH


# --- CAR ---
/path/to/env/pyscenic ctx \
  CAR_consensus_adj.csv \
  hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname analysis/pyscenic/metacell_CAR_loom.loom \
  --output CAR_consensus_reg.csv \
  --mask_dropouts \
  --num_workers 16

# --- NonCAR ---
/path/to/env/pyscenic ctx \
  NonCAR_consensus_adj.csv \
  hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname analysis/pyscenic/metacell_NonCAR_loom.loom \
  --output NonCAR_consensus_reg.csv \
  --mask_dropouts \
  --num_workers 16

# --- Mono ---
/path/to/env/pyscenic ctx \
  Mono_consensus_adj.csv \
  hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname analysis/pyscenic/metacell_Mono_loom.loom \
  --output Mono_consensus_reg.csv \
  --mask_dropouts \
  --num_workers 16

# --- AUCell Scoring ---
/path/to/env/pyscenic aucell \
  analysis/pyscenic/metacell_CAR_loom.loom \
  CAR_consensus_reg.csv \
  --output CAR_consensus_output.loom \
  --num_workers 16

/path/to/env/pyscenic aucell \
  analysis/pyscenic/metacell_NonCAR_loom.loom \
  NonCAR_consensus_reg.csv \
  --output NonCAR_consensus_output.loom \
  --num_workers 16

/path/to/env/pyscenic aucell \
  analysis/pyscenic/metacell_Mono_loom.loom \
  Mono_consensus_reg.csv \
  --output Mono_consensus_output.loom \
  --num_workers 16