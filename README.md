********************************************************************************************************************************************************************************************************************************
**
**                                                                                          U S A G E :
**
**      default on lxplus: make changes in include/configuration.h according to your signal samples. Run setup.sh; make; ./bin/combine_AmupiTrack_1D.exe v1-simple_nsim_dg_wf_cheb3_PT_1D
**      This folds your signal samples with the raw tracking asymmetries in 1D according to the following rules of the "systematics" string:
**                 v1 : version 1 of files in the dump directory
**                 simple: choice of variable binning; look into Binning directory for more options 
**                 nsim: sequential fits, 'sim' for simultaneous fits (they don't give good results atm)
**                 dg: signal shape. other options: 'cb','sg' for crystal ball and single gaussian signal shapes
**                 wf: fitting specification. 'wf' fits for withfactor between wide and narrow gaussian; other options: 'sgf', 'def' for single gaussian signal in the fail category, and fixed widthfactor respectively
**                 cheb3: background shape. 'cheb{1,2,3}' for chebychev polinomials of Nth order; 'exp' for exponential background
**                  PT_1D: kinematic dimensions; other options: 'P','ETA'
**                 
********************************************************************************************************************************************************************************************************************************
