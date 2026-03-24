# BOSSrussia
Analysis of ice-associated seal abundance in the western Bering Sea (2012 and 2013) and in the Sea of Okhotsk (2013). This repository is associated with the paper

"Abundance and distribution of ice-associated seals in the western
Bering Sea and Sea of Okhotsk, 2012-2013"  by P.B. Conn, I.S. Trukhanova, P.L. Boveng, A.N. Vasiliev, and V.I. Chernook (to be submitted to Marine Mammal Science sometime in 2026)

Analysis of western Bering Sea data in 2012 and 2013 is conducted with the R scripts
"run_wBS_2012.R" and "run_wBS_2013.R", with bootstrap-based variance requiring additional runs of 
"run_wBS_2012_boot.R" and "run_wBS_2013_boot.R."  Sensitivity analysis scripts follow a similar naming convention.
All of these scripts are located in the "\inst" directory and require compiling Template Model Builder (TMB) C++ functions
within the "\src" directory.  Prospective users will likely need to modify paths to various files within the R scripts
to reflect their local installation directories before running.

Similarly, the Okhotsk analysis for 2013 can be replicated using the scripts "run_Okhotsk_Jan2022.R", "run_Okhotsk_boot_Apr2021.R", and "run_Okhotsk_sens.R"


