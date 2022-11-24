# TCRanno
Qualitative and Quantitative Annotations of TCR Repertoire Specificity

## Requirements
Python>=3.8 \
tensorflow>=2.4.1, levenshtein, pandas, matplotlib, comut, palettable

## Installation
pip3 install tensorflow levenshtein pandas comut palettable \
pip3 install tcranno

## Quick Start
Step 1: run tcr2tcr (qualitative annotations) for input repertoire.
```
from tcranno import *
model = model_predict.load_encoder()
DB, DB_VDJ, AO_map = core_analysis.load_DB(ref_DB='IEDB')
infile = 'example/example_input_repertoire.tsv'
outprefix = 'example'
core_analysis.tcr2tcr(infile=infile, outprefix=outprefix, encoder=model, DB=DB, DB_VDJ=DB_VDJ, AO_map=AO_map, header=True, cdr3_aa_col=0, frequency=True, frequency_col=1, sep='\t', k=10)

# speed up using multi-processing
core_analysis.tcr2tcr(infile=infile, outprefix=outprefix, encoder=model, DB=DB, DB_VDJ=DB_VDJ, AO_map=AO_map, header=True, cdr3_aa_col=0, frequency=True, frequency_col=1, sep='\t', k=10, thread = -1)

# If the input file has a header and cdr3_aa column, but no frequency/count column
aa_col_name = 'name-of-the-cdr3-aa-column'
core_analysis.tcr2tcr(infile=infile, outprefix=outprefix, encoder=model, DB=DB, DB_VDJ=DB_VDJ, AO_map=AO_map, header=True, cdr3_aa_col=aa_col_name, sep='\t', k=10)

# If single-column input with no header, only a list of CDR3 sequences, one sequence per row
core_analysis.tcr2tcr(infile=infile, outprefix=outprefix, encoder=model, DB=DB, DB_VDJ=DB_VDJ, AO_map=AO_map, header=False)
```
Step 2: run tcr2ept, tcr2ag, tcr2org (quantitative annotations) based on tcr2tcr output generated in Step 1.
```
tcr2tcr_output = 'example_tcr2tcr_output.tsv'
repertoire_analysis.tcr2ept(tcr2tcr_output, outprefix, AO_map=AO_map, is_tcr2tcr=True, k=30, limit=1e-4)
repertoire_analysis.tcr2ag(tcr2tcr_output, outprefix, AO_map=AO_map, is_tcr2tcr=True, k=20, limit=1e-4)
repertoire_analysis.tcr2org(tcr2tcr_output, outprefix, AO_map=AO_map, is_tcr2tcr=True, k=10, limit=1e-4)
```
Alternatively, if you don't want to run the above commands within python environment, you can run the wrappers below that do the same thing.
```
python3 run_tcr2tcr.py --infile example_input_repertoire.tsv --outprefix example --cdr3_aa_col 0 --frequency_col 1 --k 10
python3 run_tcr2eao.py --infile example_tcr2tcr_output.tsv --is_tcr2tcr True --outprefix example --anno_type tcr2ept --k 30
python3 run_tcr2eao.py --infile example_tcr2tcr_output.tsv --is_tcr2tcr True --outprefix example --anno_type tcr2ag --k 20
python3 run_tcr2eao.py --infile example_tcr2tcr_output.tsv --is_tcr2tcr True --outprefix example --anno_type tcr2org --k 10
```
Optional: repertoire specificity landscape visualization (require Comut).
```
#choosing 'all' for anno_type will produce three plots (tcr2ept, tcr2ag, tcr2org)
python3 plot_landscape.py --tcr2tcr example_tcr2tcr_output.tsv --outprefix example --tcr2ept example_tcr2ept.tsv --tcr2ag example_tcr2ag.tsv --tcr2org example_tcr2org.tsv --anno_type all
```

## Parameters
1. core_analysis.tcr2tcr
```
--infile: 
    input file name (path), can be either a repertoire file (including at least the cdr3_aa column) or a single-column (cdr3_aa) file

--outprefix: 
    the prefix of the tcr2tcr output file to be generated

--encoder: 
    the path to the model that is used to generate embeddings (i.e., discriminative features) as a component of the tcr2tcr algorithm. Default is the pretrained model.

--DB: 
    the reference database of known TCRs, can be 'IEDB' (default) or 'DB_FULL'. 'DB_FULL' contains ~200k known TCRs; 'IEDB' is a subset of 'DB_FULL' and contains ~120k known TCRs. 'IEDB' was used throughout the original study becuase it was relatively easy for benchmarking and it has a relatively balanced sample sizes for different epitopes. Customized DB can be made and used (need to specify the path to the customized DB), as long as the DB follows the {tcr1:[ept1,ept2,ept3,...], tcr2:[eptx,epty,eptz,...],...} dictionary format.

--DB_VDJ: 
    the precomputed VDJ segments of each TCR sequences in the reference database. Since 'IEDB' is a subset of 'DB_FULL', the default DB_VDJ is the VDJ segments of the DB_FULL sequences.
    
--AO_map: 
    the epitope-to-antigen-to-organism map for all the epitopes in the reference database. The default AO_map contains mapping information for all epitopes (total: 1290) in the 'DB_FULL' database.
    
--header: 
    whether the input file has a header, defualt=False.

--cdr3_aa_col: 
    the column name (a string) or the column index (a number starting from 0, i.e., 0 for the first column) for the cdr3_aa column. If cdr3_aa_col is not given, the programme automatically take it as a single-column file without header.

--frequency: 
    whether the input file has a frequency column, defualt=False.

--frequency_col: 
    the column name (a string) or the column index (a number starting from 0, i.e., 0 for the first column) for the freuquency column.
    
--count: 
    whether the input file has a count column, defualt=False. If frequency column is given, the count column will be ignored.

--count_col: 
    the column name (a string) or the column index (a number starting from 0, i.e., 0 for the first column) for the count column.

--sep: 
    the delimitor of the input file, defualt='\s+'.

--k: 
    the number of closest TCRs within the reference database to be output for each query if not completely matching. Defualt=10. Note that the quantitative annotations (tcr2ept/tcr2ag/tcr2org) only take into account the best choice (k=1) for each input TCR sequence, so k=1 will be enough for the subsequent steps. Choose k=1 if you don't need multiple matched TCR sequences for each input, which can save space and time.

--thread: 
    the number of threads (cpus) to be used, default=1. Specify thread=-1 for all cpus.
```

2. repertoire_analysis.tcr2ept, repertoire_analysis.tcr2ag, repertoire_analysis.tcr2org
```
--infile: 
    input file name (path), a tcr2tcr output file

--outprefix: 
    the prefix of the tcr2ept/tcr2ag/tcr2org output file to be generated
    
--AO_map: 
    the epitope-to-antigen-to-organism map for all the epitopes in the reference database. The default AO_map contains mapping information for all epitopes (total: 1290) in the 'DB_FULL' database.

--k: 
    the number of top epitopes/antigens/organisms (with largest TCR-specific fractions) to be output. Default=30 for tcr2ept, default=20 for tcr2ag, default=10 for tcr2org.

--hit-rate: 
    the expected hit rate for the predicted match. Used to calculate a weighted fraction. Default=0.2.

--limit: 
    the filter which discard epitopes having a TCR-sepcific fraction lower than the specified limit. Defualt=1e-4. Setting limit=0 will turn off the filter.

```
