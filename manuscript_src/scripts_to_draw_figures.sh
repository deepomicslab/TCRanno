#!/bin/bash

# Suppl figure S4
python3 all_tools_performance.py benchmark/TOTAL.cdr3s.txt benchmark/TOTAL_benchmark_set_filtered.pkl Total
python3 all_tools_performance.py benchmark/SARS2.cdr3s.txt benchmark/filtered_SARS2_for_test.pkl SARS-CoV2
python3 all_tools_performance.py benchmark/FLU.cdr3s.txt benchmark/Influenza-M1-filtered.pkl Flu
python3 all_tools_performance.py benchmark/EBV.cdr3s.txt benchmark/EBV_2ept_set.pkl EBV
python3 all_tools_performance.py benchmark/YF.cdr3s.txt benchmark/Yellow_fever_set.pkl YFV

# figure 2
python3 violin_and_stacked_barplot.py

# figure 3
python3 keyword_fraction_boxplots_merged_cmv.py CMV_786_1e-4_tcranno_output/CMV_pos CMV+ CMV_786_1e-4_tcranno_output/CMV_neg CMV- tcr2org:CMV,tcr2ag:IE1,tcr2ag:pp65 CM 
python3 keyword_fraction_boxplots_merged_cmv.py CMV_786_1e-4_tcranno_output/CMV_pos CMV+ CMV_786_1e-4_tcranno_output/CMV_neg CMV- tcr2org:CMV,tcr2ag:IE1,tcr2ag:pp65 PM
python3 keyword_shared_intersect_seq_merged_venn_boxplot.py

# figure 4
python3 keyword_fraction_boxplots_merged_covid_org_6pop.py
python3 keyword_fraction_boxplots_merged_covid.py

#figure 5
python3 SLE_vs_ctrl_butterfly_plots.py

# figure 6
python3 all_cancers_violin.py
