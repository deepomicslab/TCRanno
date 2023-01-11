#!/bin/bash
# run tcr2tcr with single cpu; to speed up with mutiple cpus, use the --t paramerter.
python3 ../run_tcr2tcr.py --infile example_input_repertoire.tsv --outprefix example --cdr3_aa_col 0 --frequency_col 1 --limit 1e-4
# run tcr2ept after tcr2tcr
python3 ../run_tcr2eao.py --infile example_tcr2tcr_output.tsv --is_tcr2tcr True --outprefix example --anno_type tcr2ept --k 30
# run tcr2ag after tcr2tcr
python3 ../run_tcr2eao.py --infile example_tcr2tcr_output.tsv --is_tcr2tcr True --outprefix example --anno_type tcr2ag --k 20
# run tcr2org after tcr2tcr
python3 ../run_tcr2eao.py --infile example_tcr2tcr_output.tsv --is_tcr2tcr True --outprefix example --anno_type tcr2org --k 10

# run tcr2ept, tcr2ag, tcr2org all in one command
python3 ../run_tcr2eao.py --infile example_tcr2tcr_output.tsv --is_tcr2tcr True --outprefix example --anno_type all --k 20
# run tcr2ept, tcr2ag, tcr2org all in one command from input repertoire
python3 ../run_tcr2eao.py --infile example_input_repertoire.tsv --is_tcr2tcr False --outprefix example --cdr3_aa_col 0 --frequency_col 1 --anno_type all --k 20
# landscape plot
#python3 ../plot_landscape.py --tcr2tcr example_tcr2tcr_output.tsv --outprefix example --tcr2ept example_tcr2ept.tsv --anno_type tcr2ept
#python3 ../plot_landscape.py --tcr2tcr example_tcr2tcr_output.tsv --outprefix example --tcr2ag example_tcr2ag.tsv --anno_type tcr2ag
#python3 ../plot_landscape.py --tcr2tcr example_tcr2tcr_output.tsv --outprefix example --tcr2org example_tcr2org.tsv --anno_type tcr2org
python3 ../plot_landscape.py --tcr2tcr example_tcr2tcr_output.tsv --outprefix example --tcr2ept example_tcr2ept.tsv --tcr2ag example_tcr2ag.tsv --tcr2org example_tcr2org.tsv --anno_type all
