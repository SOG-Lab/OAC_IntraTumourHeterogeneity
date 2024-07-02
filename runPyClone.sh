#!/bin/bash

#SBATCH --job-name=PyClone
#SBATCH --output=PyCloneVI.output
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --array=0-28

projectPath=RESULTS/
cd $projectPath/PyCloneVI/
dataPath=/DATA/ForPyClone/
pyOutDir=PyCloneResults/

## Activate the conda environment
source activate pyclone-vi

donors=(OESO_0001 OESO_0003 OESO_0008 OESO_0009 OESO_0014 OESO_0015 OESO_0023 OESO_0025 OESO_0031 OESO_0040 OESO_0051 OESO_0053 OESO_P0056 OESO_0064 OESO_0070 OESO_0073 OESO_0093 OESO_0096 OESO_0098 OESO_0109 OESO_0113 OESO_0117 OESO_0119 OESO_0120 OESO_0121 OESO_0125 OESO_0118 OESO_6129 OESO_0227)
donor=${donors[$SLURM_ARRAY_TASK_ID]}

clones=40 
restarts=100

echo --- run Pyclone for donor "$donor" ---
## >&2 echo --- run Pyclone for donor $donor ---

pyOut=$pyOutDir/${donor}_${restarts}rS_ascatNGS.h5
  
pyclone-vi fit -i $dataPath/${donor}_pileup_forPyClone_ascatNGS.tsv -o $pyOut -c $clones -d beta-binomial -r $restarts
pyclone-vi write-results-file -i $pyOut -o $pyOutDir/${donor}_clones_${restarts}rS_ascatNGS.tsv

## If pre only file exists run PyClone for pre only
file=$dataPath/${donor}_pileup_forPyClone_preOnly_ascatNGS.tsv
if test -f "$file"; then
  pyOut=$pyOutDir/${donor}_${restarts}rS_preOnly_ascatNGS.h5
  pyclone-vi fit -i $file -o $pyOut -c $clones -d beta-binomial -r $restarts
  pyclone-vi write-results-file -i $pyOut -o $pyOutDir/${donor}_clones_${restarts}rS_preOnly_ascatNGS.tsv
fi
  
echo --- all done for donor "$donor" ---
>&2 echo --- all done for donor "$donor" ---

##done 
exit

