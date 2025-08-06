#!/bin/bash
#BSUB -q oversubscribed
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000

set -euo pipefail
 
source source_me.sh
export CONFIG="${PROJECT_DIR}/commands/germline_variants.config"
export REVISION="0.3.2"
 
 
# Load module dependencies
module load nextflow-23.10.0
module load /software/modules/ISG/singularity/3.11.4
  
nextflow pull 'https://github.com/team113sanger/dermatlas_germlinepost_nf'
 
nextflow run 'https://github.com/team113sanger/dermatlas_germlinepost_nf' \
-resume \
-r "${REVISION}" \
-c "${CONFIG}" \
-profile farm22