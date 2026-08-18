#!/bin/bash
#SBATCH -c 20
#SBATCH --job-name=Pyscenic_GRN_Lymphoma
#SBATCH --mem=300G
#SBATCH --mail-type=END,FAIL,BEGIN          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -t 5-10:00:00
#SBATCH --output=scenic-grn-%j.out

source ~/mambaforge/etc/profile.d/conda.sh
conda activate pyscenic_new
unset PYTHONPATH


# Run PySCENIC GRN for each dataset
/path/to/env/pyscenic grn \
  analysis/pyscenic/metacell_CAR_loom.loom \
  analysis/pyscenic/human_tfs.txt \
  -o adj_CAR_ix.csv --num_workers 20 --seed ix

/path/to/env/pyscenic grn \
  analysis/pyscenic/metacell_NonCAR_loom.loom \
  analysis/pyscenic/human_tfs.txt \
  -o adj_NonCAR_ix.csv --num_workers 20 --seed ix

/path/to/env/pyscenic grn \
  analysis/pyscenic/metacell_Mono_loom.loom \
  analysis/pyscenic/human_tfs.txt \
  -o adj_mono_ix.csv --num_workers 20 --seed ix