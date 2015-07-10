********************************************************************************************************************************************************************************************************************************
**
**                                                                                          U S A G E :
**
**      default on lxplus: make changes in include/configuration.h according to your signal samples. Run setup.sh; make; ./bin/combine_ATrack_1Dfor2Particles.exe v1-coarse_nsim_dg_wf_cheb3_PT_1D
**      This folds your signal samples with the raw tracking asymmetries in 1D according to the following rules of the "systematics" string:
**                 v1 : version 1 of files in the dump directory
**                 coarse: choice of variable binning; look into Binning directory for more options and check in ATrackRaw if corresponding asymmetries are given. Otherwise run fits for them
**                 nsim: sequential fits, 'sim' for simultaneous fits (they don't give good results atm)
**                 dg: signal shape. other options: 'cb','sg' for crystal ball and single gaussian signal shapes
**                 wf: fitting specification. 'wf' fits for withfactor between wide and narrow gaussian; other options: 'sgf', 'def' for single gaussian signal in the fail category, and fixed widthfactor respectively
**                 cheb3: background shape. 'cheb{1,2,3}' for chebychev polinomials of Nth order; 'exp' for exponential background
**                 PT_1D: kinematic dimensions; other options: 'P','ETA'
**                 
********************************************************************************************************************************************************************************************************************************
