from tcranno import *
# run tcr2tcr with defualt single cpu (thread=1); to speed up with mutiple cpus, use thread > 1 or thread = -1 (all cpus).

model = model_predict.load_encoder()
DB, DB_VDJ, AO_map = core_analysis.load_DB(ref_DB='IEDB')
infile = 'example_input_repertoire.tsv'
outprefix = 'example'
core_analysis.tcr2tcr(infile=infile, outprefix=outprefix, encoder=model, DB=DB, DB_VDJ=DB_VDJ, AO_map=AO_map, header=True, cdr3_aa_col=0, frequency=True, frequency_col=1, sep='\t', k=10, perform_stats=True, limit=1e-4)
tcr2tcr_output = 'example_tcr2tcr_output.tsv'
repertoire_analysis.tcr2ept(tcr2tcr_output, outprefix, AO_map=AO_map, k=30)
repertoire_analysis.tcr2ag(tcr2tcr_output, outprefix, AO_map=AO_map, k=20)
repertoire_analysis.tcr2org(tcr2tcr_output, outprefix, AO_map=AO_map, k=10)
