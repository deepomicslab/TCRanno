python3 run_tcranno_8cpus.py ImmuneAccess_TCR_data/cancer_1e-4 ImmuneAccess_TCR_data/cancer_1e-4_tcranno_output
python3 run_tcranno_8cpus.py ImmuneAccess_TCR_data/CMV_786_1e-4 ImmuneAccess_TCR_data/CMV_786_1e-4_tcranno_output
python3 run_tcranno_8cpus.py ImmuneAccess_TCR_data/covid1485_1e-4 ImmuneAccess_TCR_data/covid1485_1e-4_tcranno_output
python3 run_tcranno_8cpus.py ImmuneAccess_TCR_data/covid_vaccine_1e-4 ImmuneAccess_TCR_data/covid_vaccine_1e-4_tcranno_output
python3 run_tcranno_8cpus.py ImmuneAccess_TCR_data/SLE_1e-4 ImmuneAccess_TCR_data/sle_1e-4_tcranno_output

python3 population_analysis.py CMV_786_1e-4_tcranno_output/CMV_pos tcr2ept 30 CMV_pos.pop.tcr2ept.tsv
python3 population_analysis.py CMV_786_1e-4_tcranno_output/CMV_neg tcr2ept 30 CMV_neg.pop.tcr2ept.tsv
python3 population_analysis.py CMV_786_1e-4_tcranno_output/CMV_pos,CMV_786_1e-4_tcranno_output/CMV_neg,CMV_786_1e-4_tcranno_output/CMV_unknown tcr2ept 30 CMV_all.pop.tcr2ept.tsv

python3 population_analysis.py CMV_786_1e-4_tcranno_output/CMV_pos tcr2ag 30 CMV_pos.pop.tcr2ag.tsv
python3 population_analysis.py CMV_786_1e-4_tcranno_output/CMV_neg tcr2ag 30 CMV_neg.pop.tcr2ag.tsv
python3 population_analysis.py CMV_786_1e-4_tcranno_output/CMV_pos,CMV_786_1e-4_tcranno_output/CMV_neg,CMV_786_1e-4_tcranno_output/CMV_unknown tcr2ag 30 CMV_all.pop.tcr2ag.tsv

python3 population_analysis.py CMV_786_1e-4_tcranno_output/CMV_pos tcr2org 30 CMV_pos.pop.tcr2org.tsv
python3 population_analysis.py CMV_786_1e-4_tcranno_output/CMV_neg tcr2org 30 CMV_neg.pop.tcr2org.tsv
python3 population_analysis.py CMV_786_1e-4_tcranno_output/CMV_pos,CMV_786_1e-4_tcranno_output/CMV_neg,CMV_786_1e-4_tcranno_output/CMV_unknown tcr2org 30 CMV_all.pop.tcr2org.tsv

python3 population_analysis.py covid_vaccine_1e-4_tcranno_output/pre-vaccination tcr2ept 30 covid_pre_vaccination.pop.tcr2ept.tsv
python3 population_analysis.py covid_vaccine_1e-4_tcranno_output/post-vaccination tcr2ept 30 covid_post_vaccination.pop.tcr2ept.tsv
python3 population_analysis.py covid_vaccine_1e-4_tcranno_output/pre-vaccination tcr2ag 30 covid_pre_vaccination.pop.tcr2ag.tsv
python3 population_analysis.py covid_vaccine_1e-4_tcranno_output/post-vaccination tcr2ag 30 covid_post_vaccination.pop.tcr2ag.tsv
python3 population_analysis.py covid_vaccine_1e-4_tcranno_output/pre-vaccination tcr2org 30 covid_pre_vaccination.pop.tcr2org.tsv
python3 population_analysis.py covid_vaccine_1e-4_tcranno_output/post-vaccination tcr2org 30 covid_post_vaccination.pop.tcr2org.tsv

python3 population_analysis.py covid1485_1e-4_tcranno_output tcr2ept 30 covid1485.pop.tcr2ept.tsv
python3 population_analysis.py covid_conval_tcranno_output tcr2ept 30 covid_conval.pop.tcr2ept.tsv
python3 population_analysis.py covid_deceased_tcranno_output tcr2ept 30 covid_deceased.pop.tcr2ept.tsv
python3 population_analysis.py covid1485_1e-4_tcranno_output tcr2ag 30 covid1485.pop.tcr2ag.tsv
python3 population_analysis.py covid_conval_tcranno_output tcr2ag 30 covid_conval.pop.tcr2ag.tsv
python3 population_analysis.py covid_deceased_tcranno_output tcr2ag 30 covid_deceased.pop.tcr2ag.tsv
python3 population_analysis.py covid1485_1e-4_tcranno_output tcr2org 30 covid1485.pop.tcr2org.tsv
python3 population_analysis.py covid_conval_tcranno_output tcr2org 30 covid_conval.pop.tcr2org.tsv
python3 population_analysis.py covid_deceased_tcranno_output tcr2org 30 covid_deceased.pop.tcr2org.tsv

python3 population_analysis.py SLE_1e-4_tcranno_output tcr2ept 30 SLE.pop.tcr2ept.tsv
python3 population_analysis.py SLE_1e-4_tcranno_output tcr2ag 30 SLE.pop.tcr2ag.tsv
python3 population_analysis.py SLE_1e-4_tcranno_output tcr2org 30 SLE.pop.tcr2org.tsv
