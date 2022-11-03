# TCRanno
Quantitative and Qualitative Annotations of TCR Repertoire Specificity

## Requirements
Python>=3.8 \
tensorflow>=2.4.1, levenshtein, pandas, matplotlib, comut, palettable

## Installation
pip3 install tensorflow levenshtein pandas comut palettable \
pip3 install tcranno

## Usage
Step 1: run tcr2tcr (qualitative annotations) for input repertoire
```
from tcranno import *
model = model_predict.load_encoder()
DB, DB_VDJ, AO_key, AO_map = core_analysis.load_DB(ref_DB='IEDB')
infile = 'example_input_repertoire.tsv'
tcr2tcr_output = 'example_tcr2tcr_output.tsv'
core_analysis.tcr2tcr(infile, tcr2tcr_output, model, DB, DB_VDJ, AO_key, AO_map, header=True, cdr3_aa_col=0, frequency=True, frequency_col=1, sep='\t', k=10)

# single-column input (no header, only a list of CDR3 sequences, one sequence per row)
core_analysis.tcr2tcr(infile, tcr2tcr_output, model, DB, DB_VDJ, AO_key, AO_map, header=False, k=10)

# speed up using multi-processing
#core_analysis.tcr2tcr(infile, outfile, model, DB, DB_VDJ, AO_key, AO_map, header=True, cdr3_aa_col=0, frequency=True, frequency_col=1, sep='\t', k=10, thread = -1)
```
Step 2: run tcr2ept, tcr2ag, tcr2org (quantitative annotations) based on tcr2tcr output
```
outprefix = 'example'
repertoire_analysis.tcr2ept(tcr2tcr_output, outprefix, AO_key=AO_key, AO_map=AO_map, intermediate=True, k=30, limit=1e-4)
repertoire_analysis.tcr2ag(tcr2tcr_output, outprefix, AO_key=AO_key, AO_map=AO_map, intermediate=True, k=20, limit=1e-4)
repertoire_analysis.tcr2org(tcr2tcr_output, outprefix, AO_key=AO_key, AO_map=AO_map, intermediate=True, k=10, limit=1e-4)
```
