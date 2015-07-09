#ifndef FIT_JPSI_H
#define FIT_JPSI_H

/*!

  \date July, 1st 2015
  \author Marian Stahl
  \brief J/psi fitter for tracking asymmetries
  This code is a complete mess and needs huge effort to be cleaned up properly

*/

//C(++)
#include <vector>
//ROOT
#include "TObject.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLine.h"
//RooFit
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooArgList.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooGlobalFunc.h"
//base-class
#include "configuration.h"

using namespace RooFit;
using namespace std;

//global variables
int cpus_for_fit = 1;
const bool erroroption = false;
const bool simple_fits = false;
extern TString temp;

struct parameter_values{
  double start;
  double lo;
  double hi;
};

//functions to check if the fit converged and to re-fit the samples
double* Fit_Jpsi_nsim(RooDataHist *plusdatapass, RooDataHist *plusdatafail, RooDataHist *minusdatapass, RooDataHist *minusdatafail, configuration* myconfig);
double* Fit_Jpsi_sim(RooDataHist *plusdatapass, RooDataHist *plusdatafail, RooDataHist *minusdatapass, RooDataHist *minusdatafail, configuration* myconfig);

bool convcheck(RooFitResult *fitres, vector<RooRealVar*> vars,vector<RooRealVar*> posis, RooAbsPdf *model, RooDataHist *data, configuration* myconfig);
bool errcheck(vector<RooRealVar*> posis);
bool paratlimit(vector<RooRealVar*> vars);
void setnewparranges(vector<RooRealVar*> vars);

double* Fit_Jpsi(RooDataHist *plusdatapass, RooDataHist *plusdatafail, RooDataHist *minusdatapass, RooDataHist *minusdatafail, configuration* myconfig){
  if(myconfig->is_simultaneous()){
    cout << "Welcome to the fitter. You chose to use a simultaneous fit" << endl;
    if(myconfig->get_verbosity() > 1)cout << "signalfitmodel: " <<  myconfig->get_sigmodel() << "\t backgroundfitmodel: " << myconfig->get_bkgmodel() << endl;
    double* asymmetries;
    asymmetries = Fit_Jpsi_sim(plusdatapass, plusdatafail, minusdatapass, minusdatafail, myconfig);
    return asymmetries;
  }
  else{
    cout << "Welcome to the fitter. You chose NOT to use a simultaneous fit" << endl;
    double* asymmetries;
    asymmetries = Fit_Jpsi_nsim(plusdatapass, plusdatafail, minusdatapass, minusdatafail, myconfig);
    return asymmetries;
  }
}

double* Fit_Jpsi_nsim(RooDataHist *plusdatapass, RooDataHist *plusdatafail, RooDataHist *minusdatapass, RooDataHist *minusdatafail, configuration* myconfig){

  double binsize = (myconfig->get_Jpsi_M_max() - myconfig->get_Jpsi_M_min())/myconfig->get_fitbinning();
  double pull_threshold = 4.5;

  //cout << "building model.." << endl;

  // --- Observable ---
  RooRealVar Jpsi_M("J_psi_1S_M","M_{inv}(#mu^{+}#mu^{-})",myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max(),"MeV") ;

  // --- Signal J/psis ---
  // Define RooRealVars as pointers to get the correct values after convergence checks (C++: "arguments passed by reference")
  RooRealVar *Jpsi_sigmean = new RooRealVar("#mu_{J/#psi}","#mu_{J/#psi}",myconfig->get_fitparams().MJpsistart,myconfig->get_fitparams().MJpsilo,myconfig->get_fitparams().MJpsihi) ;
  RooRealVar *Jpsi_sigma1 = new RooRealVar("#sigma_{J/#psi,1}","#sigma_{J/#psi,1}",myconfig->get_fitparams().sigma1start,myconfig->get_fitparams().sigma1lo,myconfig->get_fitparams().sigma1hi) ;//,myconfig->get_fitparams().sigma1lo,myconfig->get_fitparams().sigma1hi
  RooGaussian Jpsi_gauss1("Jpsi_gauss1","Gaussian PDF 1",Jpsi_M,*Jpsi_sigmean,*Jpsi_sigma1) ;

  RooRealVar *CB_n = new RooRealVar("CB_n","CB parameter n",myconfig->get_fitparams().fixCBn) ;
  RooRealVar *CB_alpha = new RooRealVar("CB_alpha","CB parameter alpha",myconfig->get_fitparams().CBalphastart,myconfig->get_fitparams().CBalphalo,myconfig->get_fitparams().CBalphahi) ;
  RooCBShape Jpsi_CB("Jpsi_CB","Crystal Ball shape",Jpsi_M,*Jpsi_sigmean,*Jpsi_sigma1,*CB_alpha,*CB_n);

  RooRealVar *widthfactor;
  if(myconfig->get_shared_params() == 1)
    widthfactor = new RooRealVar("widthfactor","widthfactor",myconfig->get_fitparams().wfstart,myconfig->get_fitparams().wflo,myconfig->get_fitparams().wfhi) ;
  else
    widthfactor = new RooRealVar("widthfactor","widthfactor",myconfig->get_fitparams().fixwf);
  RooFormulaVar *Jpsi_sigma2 = new RooFormulaVar("gen","@0*@1",RooArgList(*Jpsi_sigma1,*widthfactor));
  //RooRealVar *Jpsi_sigma2 = new RooRealVar("#sigma_{J/#psi,2}","#sigma_{J/#psi,2}",sigma2start,sigma2lo,sigma2hi) ;
  RooGaussian Jpsi_gauss2("Jpsi_gauss2","Gaussian PDF 2",Jpsi_M,*Jpsi_sigmean,*Jpsi_sigma2) ;

  RooRealVar *signalfraction = new RooRealVar("signalfraction", "signalfraction",myconfig->get_fitparams().sfracstart,myconfig->get_fitparams().sfraclo,myconfig->get_fitparams().sfrachi);
  RooAddPdf *signal;

  if(myconfig->get_sigmodel() == 1){
    signal = new RooAddPdf("signal","s",RooArgList(Jpsi_gauss1,Jpsi_gauss2),*signalfraction) ;
  }
  if(myconfig->get_sigmodel() == 2){
    signal  = new RooAddPdf("signal","s",RooArgList(Jpsi_CB,Jpsi_gauss2),*signalfraction) ;//
  }

  // --- Background ---
  RooRealVar *cheba = new RooRealVar("cheba","chebBackground parameter a",myconfig->get_fitparams().chebastart,myconfig->get_fitparams().chebalo,myconfig->get_fitparams().chebahi ) ;
  RooRealVar *chebb = new RooRealVar("chebb","chebBackground parameter b",myconfig->get_fitparams().chebbstart,myconfig->get_fitparams().chebblo,myconfig->get_fitparams().chebbhi) ;
  RooRealVar *chebc = new RooRealVar("chebc","chebBackground parameter c",myconfig->get_fitparams().chebcstart,myconfig->get_fitparams().chebclo,myconfig->get_fitparams().chebchi) ;
  RooRealVar *tau = new RooRealVar("tau","#tau",myconfig->get_fitparams().taustart,myconfig->get_fitparams().taulo,myconfig->get_fitparams().tauhi);

  RooArgList bkg_par_list(*cheba);
  if(myconfig->get_bkgmodel() == 2)
    bkg_par_list.add(*chebb);
  if(myconfig->get_bkgmodel() == 3){
    bkg_par_list.add(*chebb);
    bkg_par_list.add(*chebc);
  }

  RooChebychev chebbkg("chebbkg","cheb Background PDF",Jpsi_M,bkg_par_list);
  RooExponential expbkg("expbkg","exp Background PDF",Jpsi_M,*tau);

  // --- Signal + Background ---
  RooRealVar *signal_yield = new RooRealVar("N_{J/#psi}","Jpsi signal yield",myconfig->get_fitparams().Jpsisig1start,myconfig->get_fitparams().Jpsisiglo,myconfig->get_fitparams().Jpsisighi);
  RooRealVar *background_yield = new RooRealVar("N_{bkg}","background yield",myconfig->get_fitparams().bkgstart,myconfig->get_fitparams().bkglo,myconfig->get_fitparams().bkghi);
  RooExtendPdf* sig;
  if(myconfig->get_sigmodel() == 0)sig = new RooExtendPdf("sig","Extendend signal PDF",Jpsi_gauss1,*signal_yield);
  else sig = new RooExtendPdf("sig","Extendend signal PDF",*signal,*signal_yield);
  RooExtendPdf* bkg;
  if(myconfig->get_bkgmodel() == 0)bkg = new RooExtendPdf("bkg","Extendend background PDF",expbkg,*background_yield);
  else bkg = new RooExtendPdf("bkg","Extendend background PDF",chebbkg,*background_yield);
  RooAddPdf* sum = new RooAddPdf("sum","Full Model pdf", RooArgList(*sig, *bkg));

  RooRealVar *effz = new RooRealVar("effz","effz",myconfig->get_fitparams().effstart,myconfig->get_fitparams().efflo,myconfig->get_fitparams().effhi);
  //if(large_asymmetry)effz->setRange(0.25,1.0);
  RooFormulaVar *signal_yield_fail = new RooFormulaVar("signal_yield_fail","@0*((1/@1)-1)",RooArgList(*signal_yield,*effz));
  RooExtendPdf* sig_fail;
  if(myconfig->get_sigmodel() == 0 || myconfig->get_shared_params() == 2)sig_fail = new RooExtendPdf("sig_fail","Extendend signal PDF",Jpsi_gauss1,*signal_yield_fail);
  else sig_fail = new RooExtendPdf("sig_fail","Extendend signal PDF",*signal,*signal_yield_fail);
  RooAddPdf* sum_fail = new RooAddPdf("sum_fail","Full Model pdf", RooArgList(*sig_fail,*bkg));

  //Some things for convergence checks
  //the vector vars_for_convcheck contains the parameters which are checked for convergence. Took signal and background yield away, to not spoil the efficiencies
  vector<RooRealVar*> vars_for_convcheck;
  vars_for_convcheck.push_back(Jpsi_sigmean);vars_for_convcheck.push_back(Jpsi_sigma1);
  if(myconfig->get_sigmodel() > 0){vars_for_convcheck.push_back(widthfactor);vars_for_convcheck.push_back(signalfraction);}
  if(myconfig->get_sigmodel() > 1){vars_for_convcheck.push_back(CB_alpha);}//vars_for_convcheck.push_back(CB_n);
  if(myconfig->get_bkgmodel() == 0)vars_for_convcheck.push_back(tau);
  if(myconfig->get_bkgmodel() == 1)vars_for_convcheck.push_back(cheba);
  if(myconfig->get_bkgmodel() == 2){vars_for_convcheck.push_back(cheba);vars_for_convcheck.push_back(chebb);}
  if(myconfig->get_bkgmodel() == 3){vars_for_convcheck.push_back(cheba);vars_for_convcheck.push_back(chebb);vars_for_convcheck.push_back(chebc);}

  //varibales to be reset to their initial values after fit to mu+ probe spectra
  vector<RooRealVar*> vars;
  vars.push_back(Jpsi_sigmean);
  vars.push_back(Jpsi_sigma1);
  vars.push_back(signal_yield);
  vars.push_back(background_yield);
  if(myconfig->get_sigmodel() > 0){
    if(myconfig->get_shared_params() == 1)
      vars.push_back(widthfactor);
    vars.push_back(signalfraction);
  }
  if(myconfig->get_sigmodel() > 1){vars.push_back(CB_alpha);}//vars.push_back(CB_n);
  if(myconfig->get_bkgmodel() == 0)vars.push_back(tau);
  if(myconfig->get_bkgmodel() == 1)vars.push_back(cheba);
  if(myconfig->get_bkgmodel() == 2){vars.push_back(cheba);vars.push_back(chebb);}
  if(myconfig->get_bkgmodel() == 3){vars.push_back(cheba);vars.push_back(chebb);vars.push_back(chebc);}

  vector<RooRealVar*> shared_vars;
  shared_vars.push_back(Jpsi_sigmean);
  if(myconfig->get_shared_params() != 2)shared_vars.push_back(Jpsi_sigma1);
  shared_vars.push_back(signal_yield);
  if(myconfig->get_sigmodel() > 0){
    if(myconfig->get_shared_params() == 1)
      shared_vars.push_back(widthfactor);
    shared_vars.push_back(signalfraction);
  }
  if(myconfig->get_sigmodel() > 1){shared_vars.push_back(CB_alpha);}//shared_vars.push_back(CB_n);

  vector<RooRealVar*> fail_floating_vars;
  fail_floating_vars.push_back(background_yield);
  fail_floating_vars.push_back(effz);
  if(myconfig->get_shared_params() == 2)fail_floating_vars.push_back(Jpsi_sigma1);
  if(myconfig->get_bkgmodel() == 0)fail_floating_vars.push_back(tau);
  if(myconfig->get_bkgmodel() == 1)fail_floating_vars.push_back(cheba);
  if(myconfig->get_bkgmodel() == 2){fail_floating_vars.push_back(cheba);fail_floating_vars.push_back(chebb);}
  if(myconfig->get_bkgmodel() == 3){fail_floating_vars.push_back(cheba);fail_floating_vars.push_back(chebb);fail_floating_vars.push_back(chebc);}

  vector<RooRealVar*> posis;
  posis.push_back(signal_yield);posis.push_back(background_yield);
  vector<RooRealVar*> reco_posis;
  reco_posis.push_back(background_yield);

  vector<parameter_values> initial_range;
  for (vector<RooRealVar*>::const_iterator siter = vars.begin() ;  vars.end() != siter; ++siter) {
    parameter_values the_range;
    RooRealVar *iter = *siter;
    the_range.start = iter->getValV();
    the_range.lo = iter->getMin();
    the_range.hi = iter->getMax();
    initial_range.push_back(the_range);
  }

  vector<parameter_values> initial_range_for_fail;
  for (vector<RooRealVar*>::const_iterator siter = fail_floating_vars.begin() ;  fail_floating_vars.end() != siter; ++siter) {
    parameter_values the_range;
    RooRealVar *iter = *siter;
    the_range.start = iter->getValV();
    the_range.lo = iter->getMin();
    the_range.hi = iter->getMax();
    initial_range_for_fail.push_back(the_range);
  }

  float pavetextsize = 0.04;
  float labeltextsize = 0.05;

  RooFitResult *fitres;

  //cout << "Fitting.." << endl;

  if(myconfig->get_verbosity() > 0) fitres = sum->fitTo(*plusdatapass,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;//Minos(true) ,NumCPU(2)
  else fitres = sum->fitTo(*plusdatapass,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;

  if(simple_fits == false){
    if(convcheck(fitres,vars_for_convcheck,posis,sum,plusdatapass,myconfig)==false){
      cerr << "\033[0;31mError while fitting plusdatapass \033[0m" << endl << endl << endl;
      //return NULL;
    }
  }

  TCanvas* c1 = new TCanvas("c1","Fit pass",10,10,1280,2048) ;
  float padsplitmargin = 1e-6;
  float Left_margin  = 0.15;
  float Right_margin = 0.005;
  float Bigy_margin = 0.13;
  float padsplit = 0.25;//height of the pull histo
  c1->Divide(1,2,padsplitmargin,padsplitmargin);
  c1->cd(1);

  //Must control the margins in the individual pads and not the subpads of the canvas, because of the axis labels
  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  TPad *pad1 = new TPad("pad1","pad1",padsplitmargin,padsplit,1-padsplitmargin,1-padsplitmargin);
  pad1->Draw();
  pad1->cd();
  pad1->SetBorderSize(0);
  pad1->SetLeftMargin(Left_margin);
  pad1->SetRightMargin(Right_margin);
  pad1->SetTopMargin(Bigy_margin/(1-padsplit));
  pad1->SetBottomMargin(padsplitmargin);
  pad1->SetTickx();
  pad1->SetTicky();
  //pad1->SetLogy();
  c1->Update();

  // --- Plot ---
  RooPlot* Jpsi_Mframeplus = Jpsi_M.frame(Title(" ")) ;

  plusdatapass->plotOn(Jpsi_Mframeplus,Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1)) ;
  sum->plotOn(Jpsi_Mframeplus,Components("bkg"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines()) ;
  sum->plotOn(Jpsi_Mframeplus,Components("bkg"),LineStyle(kDashed),LineWidth(2),LineColor(4)) ;
  sum->plotOn(Jpsi_Mframeplus,LineColor(4),LineWidth(2)) ;
  RooHist* ptpull = Jpsi_Mframeplus->pullHist() ;
  int bins_gt_threshold = 0;
  if(simple_fits == false){
    double *plus_matched_pulls = ptpull->GetY();
    for(int ij = 0; ij < myconfig->get_fitbinning(); ij++){
      if(plus_matched_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in plusdatapass \033[0m" << endl << endl;
      //return NULL;
    }
  }
  plusdatapass->plotOn(Jpsi_Mframeplus,Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1)) ;

  double nsigplus = signal_yield->getVal();
  //double passnbkgplus = background_yield->getVal();

  TPaveText *ptplus = new TPaveText(0.64,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),1-pad1->GetRightMargin()-0.03,1-pad1->GetTopMargin()-0.03,"BRNDC");
  ptplus->SetName("paramBox");
  ptplus->SetBorderSize(0);ptplus->SetFillColor(kWhite);ptplus->SetTextAlign(12);ptplus->SetTextSize(pavetextsize/(1-padsplit));ptplus->SetTextFont(42);ptplus->SetFillStyle(0);
  ptplus->AddText(myconfig->get_methodlabel());
  ptplus->AddText("matched #mu^{+} probe legs");
  temp.Form( "N_{J/#psi,eff} = %.0f #pm %.0f", nsigplus,signal_yield->getError());
  ptplus->AddText(temp);
  Jpsi_Mframeplus->addObject(ptplus) ;

  TPaveText *ptparbinp = new TPaveText(pad1->GetLeftMargin()+0.03,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),0.4,1-pad1->GetTopMargin()-0.03,"BRNDC");
  ptparbinp->SetBorderSize(0);ptparbinp->SetFillColor(kWhite);ptparbinp->SetTextAlign(12);ptparbinp->SetTextSize(pavetextsize/(1-padsplit));ptparbinp->SetTextFont(42);ptparbinp->SetFillStyle(0);
  if(myconfig->get_dim1().CompareTo("ETA") == 0)temp.Form("%s #in [%g,%g]",myconfig->get_dim1_label().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi());
  else temp.Form("%s #in [%g,%g] GeV",myconfig->get_dim1_label().Data(),myconfig->get_dim1_binlo()/1000,myconfig->get_dim1_binhi()/1000);
  if(!myconfig->is_1D()){
    if(myconfig->get_dim2().CompareTo("ETA") == 0)temp.Form("%s #in [%g,%g]",myconfig->get_dim2_label().Data(),myconfig->get_dim2_binlo(),myconfig->get_dim2_binhi());
    else temp.Form("%s #in [%g,%g] GeV",myconfig->get_dim2_label().Data(),myconfig->get_dim2_binlo()/1000,myconfig->get_dim2_binhi()/1000);
  }
  ptparbinp->AddText(temp);
  //ptparbinp->AddText(myconfig->get_channellabel());
  ptparbinp->AddText(myconfig->get_samplelabel());
  Jpsi_Mframeplus->addObject(ptparbinp);

  TPaveText *pt3 = new TPaveText(pad1->GetLeftMargin(),1.002-pad1->GetTopMargin(),pad1->GetLeftMargin()+0.1,1,"BRNDC");
  pt3->SetBorderSize(0);pt3->SetFillColor(kWhite);pt3->SetFillStyle(1001);pt3->SetTextColor(kWhite);
  pt3->AddText(" ");
  Jpsi_Mframeplus->addObject(pt3);

  TPaveText *pt4 = new TPaveText(1.002-pad1->GetRightMargin(),pad1->GetBottomMargin(),1,pad1->GetBottomMargin()+0.05,"BRNDC");
  pt4->SetBorderSize(0);pt4->SetFillColor(kWhite);pt4->SetFillStyle(1001);pt4->SetTextColor(kWhite);
  pt4->AddText(" ");
  Jpsi_Mframeplus->addObject(pt4);

  Jpsi_Mframeplus->GetXaxis()->SetNdivisions(305) ;
  char yaxislabel[100];
  if(Jpsi_Mframeplus->GetMaximum() > 1e+3 && Jpsi_Mframeplus->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);
  //sprintf(yaxislabel,"weighted Events/%.0f MeV",binsize);
  Jpsi_Mframeplus->GetYaxis()->SetTitle(yaxislabel) ;
  Jpsi_Mframeplus->GetYaxis()->SetTitleSize(labeltextsize/(1-padsplit)) ;
  Jpsi_Mframeplus->GetYaxis()->SetLabelSize(labeltextsize/(1-padsplit)) ;
  Jpsi_Mframeplus->GetYaxis()->SetTitleOffset(1.2*(1-padsplit)) ;
  Jpsi_Mframeplus->GetYaxis()->SetRangeUser(0.1,Jpsi_Mframeplus->GetMaximum());
  //Jpsi_Mframeplus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeplus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeplus->GetMaximum())));
  Jpsi_Mframeplus->Draw() ;

  //reset prameter ranges
  int j = 0;
  for (vector<RooRealVar*>::const_iterator siter = fail_floating_vars.begin(); fail_floating_vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    //It is important to set the range first, because after the first fit, the parameters are constrained to tight ranges and the RooRealVar::setVal(x) checks if x is in the given range!
    iter->setRange( initial_range_for_fail.at(j).lo , initial_range_for_fail.at(j).hi );
    iter->setVal( initial_range_for_fail.at(j).start );
    if(myconfig->get_verbosity() > 0)cout << "Setting ranges of " << iter->GetName() << " back to [" << iter->getMin() << "," <<iter->getMax() << "] starting at " << iter->getVal() <<  endl;
    j++;
  }

  for (vector<RooRealVar*>::const_iterator siter = shared_vars.begin() ;  shared_vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    iter->setConstant(true);
    if(myconfig->get_verbosity() > 0) cout << "Fixing " << iter->GetName() << " at value: " << iter->getVal() <<  endl;
  }

  //signal_yield->setRange(myconfig->get_fitparams().Jpsisiglo,passnsigplus);

  RooFitResult *fitRespr;
  if(myconfig->get_verbosity() > 0) fitRespr = sum_fail->fitTo(*plusdatafail,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;
  else fitRespr = sum_fail->fitTo(*plusdatafail,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false)) ;//,EvalErrorWall(false)

  if(simple_fits == false){
    if(convcheck(fitRespr,vars_for_convcheck,reco_posis,sum_fail,plusdatafail,myconfig)==false){
      cerr << "\033[0;31mError while fitting plusdatafail\033[0m" << endl << endl << endl;
      //return NULL;
    }
    if(effz->getMax() > 1)effz->setMax(1);
    if(effz->getMax() < 0)effz->setMin(0);
  }

  /*RooAbsReal *neg_log_likelihood = sum_fail->createNLL(*plusdatafail,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false)) ;
    RooMinuit(*neg_log_likelihood).migrad() ;*/

  /*double reconsigplus = signal_yield->getVal();
    double reconbkgplus = background_yield->getVal();*/
  double epsilonplus = effz->getVal();
  /*double lo_delta_epsilonplus = fabs(effz->getErrorLo());
    double hi_delta_epsilonplus = effz->getErrorHi();*/
  double delta_epsilonplus = effz->getError();

  TCanvas* c2 = new TCanvas("c2","Fit fail",10,10,1280,2048) ;
  c2->Divide(1,2,padsplitmargin,padsplitmargin);
  c2->cd(1);

  //Must control the margins in the individual pads and not the subpads of the canvas, because of the axis labels
  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  TPad *pad5 = (TPad*)pad1->Clone("pad5");
  pad5->Draw();
  pad5->cd();
  c2->Update();

  // --- Plot ---
  RooPlot* Jpsi_Mframeplus_fail = Jpsi_M.frame(Title(" ")) ;

  plusdatafail->plotOn(Jpsi_Mframeplus_fail,Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1)) ;
  sum_fail->plotOn(Jpsi_Mframeplus_fail,Components("bkg"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines()) ;
  sum_fail->plotOn(Jpsi_Mframeplus_fail,Components("bkg"),LineStyle(kDashed),LineWidth(2),LineColor(4)) ;
  sum_fail->plotOn(Jpsi_Mframeplus_fail,LineColor(4),LineWidth(2)) ;
  RooHist* prpull = Jpsi_Mframeplus_fail->pullHist() ;
  if(simple_fits == false){
    double *plus_unmatched_pulls = prpull->GetY();
    bins_gt_threshold = 0;
    for(int ij = 0; ij < myconfig->get_fitbinning(); ij++){
      if(plus_unmatched_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in plusdatafail \033[0m" << endl << endl;
      //return NULL;
    }
  }
  plusdatafail->plotOn(Jpsi_Mframeplus_fail,Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1)) ;

  TPaveText *prplus = new TPaveText(0.60,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),1-pad1->GetRightMargin()-0.03,1-pad1->GetTopMargin()-0.03,"BRNDC");
  prplus->SetName("paramBox");
  prplus->SetBorderSize(0);prplus->SetFillColor(kWhite);prplus->SetTextAlign(12);prplus->SetTextSize(pavetextsize/(1-padsplit));prplus->SetTextFont(42);prplus->SetFillStyle(0);
  prplus->AddText(myconfig->get_methodlabel());
  prplus->AddText("unmatched #mu^{+} probe legs");
  temp.Form( "#varepsilon_{+} = %.2f #pm %.2f %%", 100*epsilonplus,100*delta_epsilonplus);
  prplus->AddText(temp);
  Jpsi_Mframeplus_fail->addObject(prplus) ;
  //add the other PaveTexts
  Jpsi_Mframeplus_fail->addObject(ptparbinp) ;
  Jpsi_Mframeplus_fail->addObject(pt3);
  Jpsi_Mframeplus_fail->addObject(pt4);

  if(Jpsi_Mframeplus_fail->GetMaximum() > 1e+3 && Jpsi_Mframeplus->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus_fail->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus_fail->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);

  Jpsi_Mframeplus_fail->GetXaxis()->SetNdivisions(305) ;
  Jpsi_Mframeplus_fail->GetYaxis()->SetTitle(yaxislabel) ;
  Jpsi_Mframeplus_fail->GetYaxis()->SetTitleSize(labeltextsize/(1-padsplit)) ;
  Jpsi_Mframeplus_fail->GetYaxis()->SetLabelSize(labeltextsize/(1-padsplit)) ;
  Jpsi_Mframeplus_fail->GetYaxis()->SetTitleOffset(1.2*(1-padsplit)) ;
  //if(meth == 1)Jpsi_Mframeplus_fail->GetYaxis()->SetRangeUser(0.1,Jpsi_Mframeplus_fail->GetMaximum());
  //else
  Jpsi_Mframeplus_fail->GetYaxis()->SetRangeUser(0.1,1.4*Jpsi_Mframeplus_fail->GetMaximum());
  //Jpsi_Mframeplus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeplus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeplus->GetMaximum())));
  Jpsi_Mframeplus_fail->Draw() ;

  /*TCanvas *c35 = new TCanvas("c35","Effz likelihood",1280,960);
    RooNLLVar *nll = new RooNLLVar("nll","nll",*sum_fail,*plusdatafail);
    RooPlot *effzframe = effz->frame(Title(" "),Range(0.98,1));
    RooAbsReal* pll = neg_log_likelihood->createProfile(effz) ;
    nll->plotOn(effzframe,ShiftToZero(),Precision(1e-4));
    pll->plotOn(effzframe,ShiftToZero(),Precision(1e-4),LineColor(2),LineStyle(kDashed));
    effzframe->Draw();
    c35->SaveAs("NLL.pdf");*/

  c2->cd(1);

  TPad *pad2 = new TPad("pad2","pad2",padsplitmargin,padsplitmargin,1-padsplitmargin,padsplit);
  pad2->SetBorderSize(0);
  pad2->SetLeftMargin(Left_margin);
  pad2->SetRightMargin(Right_margin);
  pad2->SetTopMargin(padsplitmargin);
  pad2->SetBottomMargin(padsplitmargin);
  pad2->SetTickx();
  pad2->SetTicky();

  RooPlot* ppullframe = Jpsi_M.frame(Title(" ")) ;

  Color_t pullcolor = kAzure+3;
  ptpull->SetLineColor(pullcolor);
  ptpull->SetMarkerColor(pullcolor);
  prpull->SetLineColor(pullcolor);
  prpull->SetMarkerColor(pullcolor);

  ppullframe->GetXaxis()->SetTickLength(Jpsi_Mframeplus->GetXaxis()->GetTickLength()/(padsplit/(1-padsplit)));
  ppullframe->GetXaxis()->SetNdivisions(305);
  //h18->GetYaxis()->SetTickLength(Jpsi_Mframeplus->GetYaxis()->GetTickLength()/(padsplit/(1-padsplit)));
  ppullframe->GetYaxis()->SetLabelSize(Jpsi_Mframeplus->GetYaxis()->GetLabelSize()/(padsplit/(1-padsplit)));
  ppullframe->GetYaxis()->SetTitleSize(Jpsi_Mframeplus->GetYaxis()->GetTitleSize()/(padsplit/(1-padsplit)));
  ppullframe->GetYaxis()->SetNdivisions(205);
  ppullframe->GetYaxis()->SetTitleOffset(1.2*padsplit);
  ppullframe->GetYaxis()->SetTitle("Pull ");

  TLine* zeroline = new TLine(myconfig->get_Jpsi_M_min(),0,myconfig->get_Jpsi_M_max(),0);
  zeroline->SetLineStyle(2);
  zeroline->SetLineColor(1);
  TLine* m3sigma = new TLine(myconfig->get_Jpsi_M_min(),-3,myconfig->get_Jpsi_M_max(),-3);
  m3sigma->SetLineStyle(2);
  m3sigma->SetLineColor(1);
  TLine* p3sigma = new TLine(myconfig->get_Jpsi_M_min(),3,myconfig->get_Jpsi_M_max(),3);
  p3sigma->SetLineStyle(2);
  p3sigma->SetLineColor(1);
  ppullframe->addObject(zeroline,"L") ;
  ppullframe->addObject(m3sigma,"L") ;
  ppullframe->addObject(p3sigma,"L") ;
  TPaveText *pt5 = new TPaveText(1.002-pad2->GetRightMargin(),pad2->GetBottomMargin(),1,pad2->GetBottomMargin()+0.05,"BRNDC");
  pt5->SetBorderSize(0);pt5->SetFillColor(kWhite);pt5->SetFillStyle(1001);pt5->SetTextColor(kWhite);
  pt5->AddText(" ");
  ppullframe->addObject(pt5);
  ppullframe->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);

  TPad *pad6 = (TPad*)pad2->Clone("pad6");
  pad6->Draw();
  pad6->cd();
  RooPlot* ppullframe_fail = (RooPlot*)ppullframe->Clone("ppullframe_fail");
  ppullframe_fail->addPlotable(prpull,"P") ;
  ppullframe_fail->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);
  ppullframe_fail->Draw();
  c2->Update();

  c1->cd(1);
  pad2->Draw();
  pad2->cd();
  ppullframe->addPlotable(ptpull,"P") ;
  ppullframe->Draw();
  c1->Update();

  c1->cd(2);

  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  float lowerpadsplit = padsplit+(Bigy_margin/(1-padsplit));

  TPad *pad3 = new TPad("pad3","pad3",padsplitmargin,lowerpadsplit,1-padsplitmargin,1-padsplitmargin);
  pad3->Draw();
  pad3->cd();
  pad3->SetBorderSize(0);
  pad3->SetLeftMargin(Left_margin);
  pad3->SetRightMargin(Right_margin);
  pad3->SetTopMargin(padsplitmargin);
  pad3->SetBottomMargin(padsplitmargin);
  pad3->SetTickx();
  pad3->SetTicky();
  //pad3->SetLogy();
  c1->Update();

  for (vector<RooRealVar*>::const_iterator siter = shared_vars.begin() ;  shared_vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    iter->setConstant(false);
  }

  //reset prameter ranges
  j = 0;
  for (vector<RooRealVar*>::const_iterator siter = vars.begin() ;  vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    iter->setRange( initial_range.at(j).lo , initial_range.at(j).hi );
    iter->setVal( initial_range.at(j).start );
    if(myconfig->get_verbosity() > 0)cout << "Setting ranges of " << iter->GetName() << " back to [" << iter->getMin() << "," <<iter->getMax() << "] starting at " << iter->getVal() <<  endl;
    j++;
  }

  RooFitResult *fitResmt;
  if(myconfig->get_verbosity() > 0)fitResmt = sum->fitTo(*minusdatapass,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;//Minos(true) ,NumCPU(2)
  else fitResmt = sum->fitTo(*minusdatapass,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;

  if(simple_fits == false){
    if(convcheck(fitResmt,vars_for_convcheck,posis,sum,minusdatapass,myconfig)==false){
      cerr << "\033[0;31mError while fitting minusdatapass\033[0m" << endl << endl << endl;
      //return NULL;
    }
  }

  // --- Plot ---
  RooPlot* Jpsi_Mframeminus = Jpsi_M.frame(Title(" ")) ;

  minusdatapass->plotOn(Jpsi_Mframeminus,Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1)) ;
  sum->plotOn(Jpsi_Mframeminus,Components("bkg"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines()) ;
  sum->plotOn(Jpsi_Mframeminus,Components("bkg"),LineStyle(kDashed),LineWidth(2),LineColor(4)) ;
  sum->plotOn(Jpsi_Mframeminus,LineColor(4),LineWidth(2)) ;
  RooHist* mtpull = Jpsi_Mframeminus->pullHist() ;
  if(simple_fits == false){
    bins_gt_threshold = 0;
    double *minus_matched_pulls = mtpull->GetY();
    for(int ij = 0; ij < myconfig->get_fitbinning(); ij++){
      if(minus_matched_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in minusdatapass\033[0m" << endl << endl;
      //return NULL;
    }
  }
  minusdatapass->plotOn(Jpsi_Mframeminus,Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1)) ;

  double nsigminus = signal_yield->getVal();
  //double passnbkgminus = background_yield->getVal();

  TPaveText *ptminus = new TPaveText(0.64,1-gPad->GetTopMargin()-0.03-(4*pavetextsize/(1-lowerpadsplit)),1-gPad->GetRightMargin()-0.03,1-gPad->GetTopMargin()-0.03,"BRNDC");
  ptminus->SetName("paramBox");
  ptminus->SetBorderSize(0);ptminus->SetFillColor(kWhite);ptminus->SetTextAlign(12);ptminus->SetTextSize(pavetextsize/(1-lowerpadsplit));ptminus->SetTextFont(42);ptminus->SetFillStyle(0);
  ptminus->AddText(myconfig->get_methodlabel());
  ptminus->AddText("matched #mu^{-} probe legs");
  temp.Form( "N_{J/#psi,eff} = %.0f #pm %.0f", nsigminus,signal_yield->getError());
  ptminus->AddText(temp);
  Jpsi_Mframeminus->addObject(ptminus) ;

  TPaveText *ptparbin = new TPaveText(gPad->GetLeftMargin()+0.03,1-gPad->GetTopMargin()-0.03-(4*pavetextsize/(1-lowerpadsplit)),0.4,1-gPad->GetTopMargin()-0.03,"BRNDC");
  ptparbin->SetBorderSize(0);ptparbin->SetFillColor(kWhite);ptparbin->SetTextAlign(12);ptparbin->SetTextSize(pavetextsize/(1-lowerpadsplit));ptparbin->SetTextFont(42);ptparbin->SetFillStyle(0);
  if(myconfig->get_dim1().CompareTo("ETA") == 0)temp.Form("%s #in [%g,%g]",myconfig->get_dim1_label().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi());
  else temp.Form("%s #in [%g,%g] GeV",myconfig->get_dim1_label().Data(),myconfig->get_dim1_binlo()/1000,myconfig->get_dim1_binhi()/1000);
  if(!myconfig->is_1D()){
    if(myconfig->get_dim2().CompareTo("ETA") == 0)temp.Form("%s #in [%g,%g]",myconfig->get_dim2_label().Data(),myconfig->get_dim2_binlo(),myconfig->get_dim2_binhi());
    else temp.Form("%s #in [%g,%g] GeV",myconfig->get_dim2_label().Data(),myconfig->get_dim2_binlo()/1000,myconfig->get_dim2_binhi()/1000);
  }
  ptparbin->AddText(temp);
  //ptparbin->AddText(myconfig->get_channellabel());
  ptparbin->AddText(myconfig->get_samplelabel());
  Jpsi_Mframeminus->addObject(ptparbin);

  TPaveText *pt6 = new TPaveText(1.002-pad3->GetRightMargin(),pad3->GetBottomMargin(),1,pad3->GetBottomMargin()+0.05,"BRNDC");
  pt6->SetBorderSize(0);pt6->SetFillColor(kWhite);pt6->SetFillStyle(1001);pt6->SetTextColor(kWhite);
  pt6->AddText(" ");
  Jpsi_Mframeminus->addObject(pt6);

  Jpsi_Mframeminus->GetXaxis()->SetNdivisions(305) ;
  if(Jpsi_Mframeminus->GetMaximum() > 1e+3 && Jpsi_Mframeminus->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeminus->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);
  //Jpsi_Mframeminus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeminus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeminus->GetMaximum())));
  //sprintf(yaxislabel,"weighted Events/%.0f MeV",binsize);
  Jpsi_Mframeminus->GetYaxis()->SetTitle(yaxislabel) ;
  Jpsi_Mframeminus->GetYaxis()->SetTitleSize(labeltextsize/(1-lowerpadsplit)) ;
  Jpsi_Mframeminus->GetYaxis()->SetLabelSize(labeltextsize/(1-lowerpadsplit)) ;
  Jpsi_Mframeminus->GetYaxis()->SetTitleOffset(1.2*(1-lowerpadsplit)) ;
  Jpsi_Mframeminus->GetYaxis()->SetRangeUser(0.1,Jpsi_Mframeminus->GetMaximum());
  Jpsi_Mframeminus->Draw() ;

  j = 0;
  //reset prameter ranges
  for (vector<RooRealVar*>::const_iterator siter = fail_floating_vars.begin() ;  fail_floating_vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    iter->setRange( initial_range_for_fail.at(j).lo , initial_range_for_fail.at(j).hi );
    iter->setVal( initial_range_for_fail.at(j).start );
    if(myconfig->get_verbosity() > 0)cout << "Setting ranges of " << iter->GetName() << " back to [" << iter->getMin() << "," <<iter->getMax() << "] starting at " << iter->getVal() <<  endl;
    j++;
  }

  for (vector<RooRealVar*>::const_iterator siter = shared_vars.begin() ;  shared_vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    iter->setConstant(true);
    if(myconfig->get_verbosity() > 0) cout << "Fixing " << iter->GetName() << " at value: " << iter->getVal() <<  endl;
  }

  RooFitResult *fitResmr;
  if(myconfig->get_verbosity() > 0)fitResmr = sum_fail->fitTo(*minusdatafail,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;
  else fitResmr = sum_fail->fitTo(*minusdatafail,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;//

  if(simple_fits == false){
    if(convcheck(fitResmr,vars_for_convcheck,reco_posis,sum_fail,minusdatafail,myconfig)==false){
      cerr << "\033[0;31mError while fitting minusdatafail\033[0m" << endl << endl << endl;
      //return NULL;
    }
    if(effz->getMax() > 1)effz->setMax(1);
    if(effz->getMax() < 0)effz->setMin(0);
  }

  /*RooAbsReal *neg_log_likelihood_m = sum_fail->createNLL(*minusdatafail,Strategy(2),Save(true),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false)) ;
    RooMinuit(*neg_log_likelihood_m).migrad() ;*/

  /*double reconsigminus = signal_yield->getVal();
    double reconbkgminus = background_yield->getVal();*/
  double epsilonminus = effz->getVal();
  //double lo_delta_epsilonminus = fabs(effz->getErrorLo());
  //double hi_delta_epsilonminus = effz->getErrorHi();
  double delta_epsilonminus = effz->getError();

  c2->cd(2);

  //Must control the margins in the individual pads and not the subpads of the canvas, because of the axis labels
  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  TPad *pad7 = (TPad*)pad3->Clone("pad7");
  pad7->Draw();
  pad7->cd();
  c2->Update();

  RooPlot* Jpsi_Mframeminus_fail = Jpsi_M.frame(Title(" ")) ;

  minusdatafail->plotOn(Jpsi_Mframeminus_fail,Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1)) ;
  sum_fail->plotOn(Jpsi_Mframeminus_fail,Components("bkg"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines()) ;
  sum_fail->plotOn(Jpsi_Mframeminus_fail,Components("bkg"),LineStyle(kDashed),LineWidth(2),LineColor(4)) ;
  sum_fail->plotOn(Jpsi_Mframeminus_fail,LineColor(4),LineWidth(2)) ;
  RooHist* mrpull = Jpsi_Mframeminus_fail->pullHist() ;
  if(simple_fits == false){
    double *minus_unmatched_pulls = mrpull->GetY();
    bins_gt_threshold = 0;
    for(int ij = 0; ij < myconfig->get_fitbinning(); ij++){
      if(minus_unmatched_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in minusdatafail \033[0m" << endl << endl;
      //return NULL;
    }
  }
  minusdatafail->plotOn(Jpsi_Mframeminus_fail,Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1)) ;

  TPaveText *prminus = new TPaveText(0.6,1-gPad->GetTopMargin()-0.03-(4*pavetextsize/(1-lowerpadsplit)),1-gPad->GetRightMargin()-0.03,1-gPad->GetTopMargin()-0.03,"BRNDC");
  prminus->SetName("paramBox");
  prminus->SetBorderSize(0);prminus->SetFillColor(kWhite);prminus->SetTextAlign(12);prminus->SetTextSize(pavetextsize/(1-lowerpadsplit));prminus->SetTextFont(42);prminus->SetFillStyle(0);
  prminus->AddText(myconfig->get_methodlabel());
  prminus->AddText("unmatched #mu^{-} probe legs");
  temp.Form( "#varepsilon_{-} = %.2f #pm %.2f %%", 100*epsilonminus,100*delta_epsilonminus);
  prminus->AddText(temp);
  Jpsi_Mframeminus_fail->addObject(prminus) ;
  //add the other PaveTexts
  Jpsi_Mframeminus_fail->addObject(ptparbin) ;
  Jpsi_Mframeminus_fail->addObject(pt6);

  Jpsi_Mframeminus_fail->GetXaxis()->SetNdivisions(305) ;
  if(Jpsi_Mframeminus_fail->GetMaximum() > 1e+3 && Jpsi_Mframeminus_fail->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeminus_fail->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeminus_fail->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);
  //Jpsi_Mframeminus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeminus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeminus->GetMaximum())));
  //sprintf(yaxislabel,"weighted Events/%.0f MeV",binsize);
  Jpsi_Mframeminus_fail->GetYaxis()->SetTitle(yaxislabel) ;
  Jpsi_Mframeminus_fail->GetYaxis()->SetTitleSize(labeltextsize/(1-lowerpadsplit)) ;
  Jpsi_Mframeminus_fail->GetYaxis()->SetLabelSize(labeltextsize/(1-lowerpadsplit)) ;
  Jpsi_Mframeminus_fail->GetYaxis()->SetTitleOffset(1.2*(1-lowerpadsplit)) ;
  //if(meth == 1)Jpsi_Mframeminus_fail->GetYaxis()->SetRangeUser(0.1,Jpsi_Mframeminus_fail->GetMaximum());
  //else
  Jpsi_Mframeminus_fail->GetYaxis()->SetRangeUser(0.1,1.4*Jpsi_Mframeminus_fail->GetMaximum());
  Jpsi_Mframeminus_fail->Draw() ;

  c2->Update();
  c2->cd(2);

  TPad *pad4 = new TPad("pad4","pad4",padsplitmargin,padsplitmargin,1-padsplitmargin,lowerpadsplit);
  pad4->SetBorderSize(0);
  pad4->SetLeftMargin(Left_margin);
  pad4->SetRightMargin(Right_margin);
  pad4->SetTopMargin(padsplitmargin);
  pad4->SetBottomMargin(Bigy_margin/lowerpadsplit);
  pad4->SetTickx();
  pad4->SetTicky();

  mtpull->SetLineColor(pullcolor);
  mtpull->SetMarkerColor(pullcolor);
  mrpull->SetLineColor(pullcolor);
  mrpull->SetMarkerColor(pullcolor);

  RooPlot* mpullframe = Jpsi_M.frame(Title(" ")) ;

  mpullframe->GetXaxis()->SetTitleOffset(1.05);
  mpullframe->GetXaxis()->SetTitleSize(labeltextsize/lowerpadsplit);
  mpullframe->GetXaxis()->SetLabelSize(labeltextsize/lowerpadsplit);
  mpullframe->GetXaxis()->SetTitle("M_{inv}(#mu^{+}#mu^{-}) (GeV)");
  mpullframe->GetXaxis()->SetTickLength(Jpsi_Mframeplus->GetXaxis()->GetTickLength()/(lowerpadsplit/(1-lowerpadsplit)));
  //h18->GetYaxis()->SetTickLength(Jpsi_Mframeplus->GetYaxis()->GetTickLength()/(lowerpadsplit/(1-lowerpadsplit)));
  mpullframe->GetYaxis()->SetLabelSize(labeltextsize/lowerpadsplit);
  mpullframe->GetYaxis()->SetTitleSize(labeltextsize/lowerpadsplit);
  mpullframe->GetXaxis()->SetNdivisions(305);
  mpullframe->GetYaxis()->SetNdivisions(205);
  mpullframe->GetYaxis()->SetTitleOffset(1.2*lowerpadsplit);
  mpullframe->GetYaxis()->SetTitle("Pull ");

  TLine* mzeroline = new TLine(myconfig->get_Jpsi_M_min(),0,myconfig->get_Jpsi_M_max(),0);
  mzeroline->SetLineStyle(2);
  mzeroline->SetLineColor(1);
  TLine* mm3sigma = new TLine(myconfig->get_Jpsi_M_min(),-3,myconfig->get_Jpsi_M_max(),-3);
  mm3sigma->SetLineStyle(2);
  mm3sigma->SetLineColor(1);
  TLine* mp3sigma = new TLine(myconfig->get_Jpsi_M_min(),3,myconfig->get_Jpsi_M_max(),3);
  mp3sigma->SetLineStyle(2);
  mp3sigma->SetLineColor(1);
  mpullframe->addObject(mzeroline,"L") ;
  mpullframe->addObject(mm3sigma,"L") ;
  mpullframe->addObject(mp3sigma,"L") ;
  TPaveText *pt7 = new TPaveText(1.002-pad4->GetRightMargin(),pad4->GetBottomMargin(),1,pad4->GetBottomMargin()+0.05,"BRNDC");
  pt7->SetBorderSize(0);pt7->SetFillColor(kWhite);pt7->SetFillStyle(1001);pt7->SetTextColor(kWhite);
  pt7->AddText(" ");
  mpullframe->addObject(pt7);
  mpullframe->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);

  TPad *pad8 = (TPad*)pad4->Clone("pad8");
  pad8->Draw();
  pad8->cd();
  RooPlot* mpullframe_fail = (RooPlot*)mpullframe->Clone("mpullframe_fail");
  mpullframe_fail->addPlotable(mrpull,"P") ;
  mpullframe_fail->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);
  mpullframe_fail->Draw();
  c2->Update();

  c1->cd(2);
  pad4->Draw();
  pad4->cd();
  mpullframe->addPlotable(mtpull,"P") ;
  mpullframe->Draw();
  c1->Update();

  temp = myconfig->get_dumpdir()+"Fits";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);

  if(myconfig->is_1D())temp.Form("%sFits/raw_%s_%s_%s_%g_%g_pass.pdf",myconfig->get_dumpdir().Data(),
                                 myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi());
  else temp.Form("%sFits/raw_%s_%s_%s_%g_%g_%g_%g_pass.pdf",myconfig->get_dumpdir().Data(),
                 myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi(),myconfig->get_dim2_binlo(),myconfig->get_dim2_binhi());
  if(!myconfig->is_debug()) c1->SaveAs(temp);

  if(myconfig->is_1D())temp.Form("%sFits/raw_%s_%s_%s_%g_%g_fail.pdf",myconfig->get_dumpdir().Data(),
                                 myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi());
  else temp.Form("%sFits/raw_%s_%s_%s_%g_%g_%g_%g_fail.pdf",myconfig->get_dumpdir().Data(),
                 myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi(),myconfig->get_dim2_binlo(),myconfig->get_dim2_binhi());
  if(!myconfig->is_debug()) c2->SaveAs(temp);

  delete c1;
  delete c2;

  //this is the important stuff
  double* efficiencies;
  efficiencies = new double[4];
  efficiencies[0] = epsilonplus;
  if(epsilonplus + delta_epsilonplus < 1 && epsilonplus - delta_epsilonplus > 0)efficiencies[1] = delta_epsilonplus;
  else {
    if(epsilonplus + delta_epsilonplus > 1){
      efficiencies[1] = delta_epsilonplus;
      if(myconfig->get_verbosity() > 0)cout << "efficiency + upper uncertainty exceeds 1, be careful" << endl;
    }
    if(epsilonplus - delta_epsilonplus < 0){
      cout << "efficiency - lower uncertainty smaller than 0, this shouldn't happen - return NULL!" << endl;
      return NULL;
    }
  }
  efficiencies[2] = epsilonminus;
  if(epsilonminus + delta_epsilonminus < 1 && epsilonminus - delta_epsilonminus > 0)efficiencies[3] = delta_epsilonminus;
  else {
    if(epsilonminus + delta_epsilonminus > 1){
      efficiencies[3] = delta_epsilonplus;
      if(myconfig->get_verbosity() > 0)cout << "efficiency + upper uncertainty exceeds 1, be careful" << endl;
    }
    if(epsilonminus - delta_epsilonminus < 0){
      cout << "efficiency - lower uncertainty smaller than 0, this shouldn't happen - return NULL!" << endl;
      return NULL;
    }
  }

  double *asymmetries;
  asymmetries = new double[2];
  asymmetries[0] = (epsilonplus-epsilonminus)/(epsilonplus+epsilonminus);
  asymmetries[1] = 2*sqrt((pow(epsilonplus,2)*pow(delta_epsilonminus,2)+pow(epsilonminus,2)*pow(delta_epsilonplus,2))/(pow((epsilonplus+epsilonminus),4)));
  cout << "\033[0;33m Raw Asymmetry = ( " << 100*asymmetries[0] << " +- "<< 100*asymmetries[1] << " ) %\033[0m" << endl;

  delete ptplus;
  delete ptparbinp;
  delete ptminus;
  delete ptparbin;
  delete prplus;
  delete prminus;
  delete zeroline;
  delete m3sigma;
  delete p3sigma;
  delete mzeroline;
  delete mm3sigma;
  delete mp3sigma;
  delete pt3;
  delete pt4;
  delete pt5;
  delete pt6;
  delete pt7;
  delete efficiencies;
  return asymmetries;
}

double* Fit_Jpsi_sim(RooDataHist *plusdatapass, RooDataHist *plusdatafail, RooDataHist *minusdatapass, RooDataHist *minusdatafail, configuration* myconfig){

  double binsize = (myconfig->get_Jpsi_M_max() - myconfig->get_Jpsi_M_min())/myconfig->get_fitbinning();
  double pull_threshold = 4.5;

  // --- Observable ---
  RooRealVar *Jpsi_M = new RooRealVar("J_psi_1S_M","M_{inv}(#mu^{+}#mu^{-})",myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max(),"MeV") ;

  // --- Signal J/psis ---
  // Define RooRealVars as pointers to get the correct values after convergence checks (C++: "arguments passed by reference")
  RooRealVar *Jpsi_sigmean_plus = new RooRealVar("#mu_{J/#psi,+}","#mu_{J/#psi}",myconfig->get_fitparams().MJpsistart,myconfig->get_fitparams().MJpsilo,myconfig->get_fitparams().MJpsihi) ;
  RooRealVar *Jpsi_sigmean_minus = new RooRealVar("#mu_{J/#psi,-}","#mu_{J/#psi}",myconfig->get_fitparams().MJpsistart,myconfig->get_fitparams().MJpsilo,myconfig->get_fitparams().MJpsihi) ;
  RooRealVar *Jpsi_sigma1 = new RooRealVar("#sigma_{J/#psi,1}","#sigma_{J/#psi,1}",myconfig->get_fitparams().sigma1start,myconfig->get_fitparams().sigma1lo,myconfig->get_fitparams().sigma1hi) ;//,myconfig->get_fitparams().sigma1lo,myconfig->get_fitparams().sigma1hi
  RooGaussian *Jpsi_gauss1_plus = new RooGaussian("Jpsi_gauss1_plus","Gaussian PDF 1 plus",*Jpsi_M,*Jpsi_sigmean_plus,*Jpsi_sigma1) ;
  RooGaussian *Jpsi_gauss1_minus = new RooGaussian("Jpsi_gauss1_minus","Gaussian PDF 1 minus",*Jpsi_M,*Jpsi_sigmean_minus,*Jpsi_sigma1) ;

  RooRealVar *CB_n = new RooRealVar("CB_n","CB parameter n",myconfig->get_fitparams().fixCBn) ;
  RooRealVar *CB_alpha = new RooRealVar("CB_alpha","CB parameter alpha",2.0,1.0,6.0) ;
  RooCBShape *Jpsi_CB_plus = new RooCBShape("Jpsi_CB_plus","Crystal Ball shape plus",*Jpsi_M,*Jpsi_sigmean_plus,*Jpsi_sigma1,*CB_alpha,*CB_n);
  RooCBShape *Jpsi_CB_minus = new RooCBShape("Jpsi_CB_minus","Crystal Ball shape minus",*Jpsi_M,*Jpsi_sigmean_minus,*Jpsi_sigma1,*CB_alpha,*CB_n);

  RooRealVar *widthfactor;
  if(myconfig->get_shared_params() == 1)
    widthfactor = new RooRealVar("widthfactor","widthfactor",1.6,1.3,2.2);
  else
    widthfactor = new RooRealVar("widthfactor","widthfactor",myconfig->get_fitparams().fixwf);
  RooFormulaVar *Jpsi_sigma2 = new RooFormulaVar("#sigma_{J/#psi,2}","@0*@1",RooArgList(*Jpsi_sigma1,*widthfactor));
  RooGaussian *Jpsi_gauss2_plus = new RooGaussian("Jpsi_gauss2_plus","Gaussian PDF 2",*Jpsi_M,*Jpsi_sigmean_plus,*Jpsi_sigma2) ;
  RooGaussian *Jpsi_gauss2_minus = new RooGaussian("Jpsi_gauss2_minus","Gaussian PDF 2",*Jpsi_M,*Jpsi_sigmean_minus,*Jpsi_sigma2) ;

  RooRealVar *signalfraction_plus = new RooRealVar("f_{G_{1},G_{2},+}","signalfraction plus",myconfig->get_fitparams().sfracstart,myconfig->get_fitparams().sfraclo,myconfig->get_fitparams().sfrachi);
  RooRealVar *signalfraction_minus = new RooRealVar("f_{G_{1},G_{2},-}","signalfraction minus",myconfig->get_fitparams().sfracstart,myconfig->get_fitparams().sfraclo,myconfig->get_fitparams().sfrachi);

  RooArgList *signal_pdfs_plus = new RooArgList(*Jpsi_gauss1_plus);RooArgList *signal_pdfs_minus = new RooArgList(*Jpsi_gauss1_minus);

  if(myconfig->get_sigmodel() == 1){
    signal_pdfs_plus->add(*Jpsi_gauss2_plus);
    signal_pdfs_minus->add(*Jpsi_gauss2_minus);
  }
  if(myconfig->get_sigmodel() == 2){
    signal_pdfs_plus->add(*Jpsi_CB_plus);
    signal_pdfs_minus->add(*Jpsi_CB_minus);
  }
  RooAddPdf *signal_plus_pass = new RooAddPdf("signal_plus_pass","signal plus pass",*signal_pdfs_plus,*signalfraction_plus);
  RooAddPdf *signal_minus_pass = new RooAddPdf("signal_minus_pass","signal minus pass",*signal_pdfs_minus,*signalfraction_minus);

  RooAddPdf *signal_plus_fail = new RooAddPdf("signal_plus_fail","signal plus fail",*signal_pdfs_plus,*signalfraction_plus);
  RooAddPdf *signal_minus_fail = new RooAddPdf("signal_minus_fail","signal minus fail",*signal_pdfs_minus,*signalfraction_minus);

  // --- Signal extended parameters ---

  RooRealVar *average_eff = new RooRealVar("average_eff","average_eff",myconfig->get_fitparams().effstart,myconfig->get_fitparams().efflo,myconfig->get_fitparams().effhi);
  RooRealVar *asymmetry = new RooRealVar("asymmetry","asymmetry",myconfig->get_fitparams().Astart,myconfig->get_fitparams().Alo,myconfig->get_fitparams().Ahi);

  RooRealVar *total_signal_yield_plus = new RooRealVar("N_{J/#psi,+}","total Jpsi signal yield",myconfig->get_fitparams().Jpsisig1start,myconfig->get_fitparams().Jpsisiglo,myconfig->get_fitparams().Jpsisighi);
  RooFormulaVar *efficiency_plus = new RooFormulaVar("#varepsilon_{+}","@0*(1+@1)",RooArgList(*average_eff,*asymmetry));
  RooFormulaVar *signal_yield_plus_pass = new RooFormulaVar("N^{pass}_{J/#psi,+}","@0*@1",RooArgList(*total_signal_yield_plus,*efficiency_plus));
  RooFormulaVar *signal_yield_plus_fail = new RooFormulaVar("N^{fail}_{J/#psi,+}","@0*(1-@1)",RooArgList(*total_signal_yield_plus,*efficiency_plus));

  RooRealVar *total_signal_yield_minus = new RooRealVar("N_{J/#psi,-}","total Jpsi signal yield",myconfig->get_fitparams().Jpsisig1start,myconfig->get_fitparams().Jpsisiglo,myconfig->get_fitparams().Jpsisighi);
  RooFormulaVar *efficiency_minus = new RooFormulaVar("#varepsilon_{-}","@0*(1-@1)",RooArgList(*average_eff,*asymmetry));
  RooFormulaVar *signal_yield_minus_pass = new RooFormulaVar("N^{pass}_{J/#psi,-}","@0*@1",RooArgList(*total_signal_yield_minus,*efficiency_minus));
  RooFormulaVar *signal_yield_minus_fail = new RooFormulaVar("N^{fail}_{J/#psi,-}","@0*(1-@1)",RooArgList(*total_signal_yield_minus,*efficiency_minus));

  // --- Background ---

  RooRealVar *exp_tau_plus_pass = new RooRealVar("exp_tau_plus_pass","exponential parameter tau pass",myconfig->get_fitparams().taustart,myconfig->get_fitparams().taulo,myconfig->get_fitparams().tauhi);
  RooRealVar *exp_tau_plus_fail = new RooRealVar("exp_tau_plus_fail","exponential parameter tau fail",myconfig->get_fitparams().taustart,myconfig->get_fitparams().taulo,myconfig->get_fitparams().tauhi);
  RooRealVar *exp_tau_minus_pass = new RooRealVar("exp_tau_minus_pass","exponential parameter tau pass",myconfig->get_fitparams().taustart,myconfig->get_fitparams().taulo,myconfig->get_fitparams().tauhi);
  RooRealVar *exp_tau_minus_fail = new RooRealVar("exp_tau_minus_fail","exponential parameter tau fail",myconfig->get_fitparams().taustart,myconfig->get_fitparams().taulo,myconfig->get_fitparams().tauhi);

  RooRealVar *cheba_plus_pass = new RooRealVar("cheba_plus_pass","chebBackground parameter a",myconfig->get_fitparams().chebastart,myconfig->get_fitparams().chebalo,myconfig->get_fitparams().chebahi ) ;
  RooRealVar *cheba_plus_fail = new RooRealVar("cheba_plus_fail","chebBackground parameter a",myconfig->get_fitparams().chebastart,myconfig->get_fitparams().chebalo,myconfig->get_fitparams().chebahi ) ;
  RooRealVar *cheba_minus_pass = new RooRealVar("cheba_minus_pass","chebBackground parameter a",myconfig->get_fitparams().chebastart,myconfig->get_fitparams().chebalo,myconfig->get_fitparams().chebahi );
  RooRealVar *cheba_minus_fail = new RooRealVar("cheba_minus_fail","chebBackground parameter a",myconfig->get_fitparams().chebastart,myconfig->get_fitparams().chebalo,myconfig->get_fitparams().chebahi );

  RooRealVar *chebb_plus_pass = new RooRealVar("chebb_plus_pass","chebBackground parameter b",myconfig->get_fitparams().chebbstart,myconfig->get_fitparams().chebblo,myconfig->get_fitparams().chebbhi) ;
  RooRealVar *chebb_plus_fail = new RooRealVar("chebb_plus_fail","chebBackground parameter b",myconfig->get_fitparams().chebbstart,myconfig->get_fitparams().chebblo,myconfig->get_fitparams().chebbhi) ;
  RooRealVar *chebb_minus_pass = new RooRealVar("chebb_minus_pass","chebBackground parameter b",myconfig->get_fitparams().chebbstart,myconfig->get_fitparams().chebblo,myconfig->get_fitparams().chebbhi);
  RooRealVar *chebb_minus_fail = new RooRealVar("chebb_minus_fail","chebBackground parameter b",myconfig->get_fitparams().chebbstart,myconfig->get_fitparams().chebblo,myconfig->get_fitparams().chebbhi);

  RooRealVar *chebc_plus_pass = new RooRealVar("chebc_plus_pass","chebBackground parameter c",myconfig->get_fitparams().chebcstart,myconfig->get_fitparams().chebclo,myconfig->get_fitparams().chebchi) ;
  RooRealVar *chebc_plus_fail = new RooRealVar("chebc_plus_fail","chebBackground parameter c",myconfig->get_fitparams().chebcstart,myconfig->get_fitparams().chebclo,myconfig->get_fitparams().chebchi) ;
  RooRealVar *chebc_minus_pass = new RooRealVar("chebc_minus_pass","chebBackground parameter c",myconfig->get_fitparams().chebcstart,myconfig->get_fitparams().chebclo,myconfig->get_fitparams().chebchi);
  RooRealVar *chebc_minus_fail = new RooRealVar("chebc_minus_fail","chebBackground parameter c",myconfig->get_fitparams().chebcstart,myconfig->get_fitparams().chebclo,myconfig->get_fitparams().chebchi);

  RooExponential *exponential_plus_pass = new RooExponential("exponential_plus_pass","Exponential Background PDF pass",*Jpsi_M,*exp_tau_plus_pass);
  RooExponential *exponential_plus_fail = new RooExponential("exponential_plus_fail","Exponential Background PDF fail",*Jpsi_M,*exp_tau_plus_fail);
  RooExponential *exponential_minus_pass = new RooExponential("exponential_minus_pass","Exponential Background PDF pass",*Jpsi_M,*exp_tau_minus_pass);
  RooExponential *exponential_minus_fail = new RooExponential("exponential_minus_fail","Exponential Background PDF fail",*Jpsi_M,*exp_tau_minus_fail);

  RooArgList *bkg_params_plus_pass = new RooArgList(*cheba_plus_pass);
  RooArgList *bkg_params_minus_pass = new RooArgList(*cheba_minus_pass);
  RooArgList *bkg_params_plus_fail = new RooArgList(*cheba_plus_fail);
  RooArgList *bkg_params_minus_fail = new RooArgList(*cheba_minus_fail);

  if(myconfig->get_bkgmodel() == 2){
    bkg_params_plus_pass->add(*chebb_plus_pass);
    bkg_params_plus_fail->add(*chebb_plus_fail);
    bkg_params_minus_pass->add(*chebb_minus_pass);
    bkg_params_minus_fail->add(*chebb_minus_fail);
  }
  if(myconfig->get_bkgmodel() == 3){
    bkg_params_plus_pass->add(*chebb_plus_pass);
    bkg_params_plus_pass->add(*chebc_plus_pass);

    bkg_params_plus_fail->add(*chebb_plus_fail);
    bkg_params_plus_fail->add(*chebc_plus_fail);

    bkg_params_minus_pass->add(*chebb_minus_pass);
    bkg_params_minus_pass->add(*chebc_minus_pass);

    bkg_params_minus_fail->add(*chebb_minus_fail);
    bkg_params_minus_fail->add(*chebc_minus_fail);
  }

  RooChebychev *chebbkg_plus_pass = new RooChebychev("chebbkg_plus_pass","cheb Background PDF pass plus",*Jpsi_M,*bkg_params_plus_pass) ;
  RooChebychev *chebbkg_plus_fail = new RooChebychev("chebbkg_plus_fail","cheb Background PDF fail plus",*Jpsi_M,*bkg_params_plus_fail) ;
  RooChebychev *chebbkg_minus_pass = new RooChebychev("chebbkg_minus_pass","cheb Background PDF pass minus",*Jpsi_M,*bkg_params_minus_pass) ;
  RooChebychev *chebbkg_minus_fail = new RooChebychev("chebbkg_minus_fail","cheb Background PDF fail minus",*Jpsi_M,*bkg_params_minus_fail) ;

  RooRealVar *background_yield_plus_pass = new RooRealVar("N_{bkg,pass,+}","background yield pass",myconfig->get_fitparams().bkgstart,myconfig->get_fitparams().bkglo,myconfig->get_fitparams().bkghi);
  RooRealVar *background_yield_plus_fail = new RooRealVar("N_{bkg,fail,+}","background yield fail",myconfig->get_fitparams().bkgstart,myconfig->get_fitparams().bkglo,myconfig->get_fitparams().bkghi);
  RooRealVar *background_yield_minus_pass = new RooRealVar("N_{bkg,pass,-}","background yield pass",myconfig->get_fitparams().bkgstart,myconfig->get_fitparams().bkglo,myconfig->get_fitparams().bkghi);
  RooRealVar *background_yield_minus_fail = new RooRealVar("N_{bkg,fail,-}","background yield fail",myconfig->get_fitparams().bkgstart,myconfig->get_fitparams().bkglo,myconfig->get_fitparams().bkghi);


  RooAddPdf* sum_plus_pass; RooAddPdf* sum_plus_fail; RooAddPdf* sum_minus_pass; RooAddPdf* sum_minus_fail;
  if(myconfig->get_bkgmodel() == 0)
    sum_plus_pass = new RooAddPdf("sum_plus_pass","Full Model pdf plus pass",RooArgList(*signal_plus_pass,*exponential_plus_pass),RooArgList(*signal_yield_plus_pass,*background_yield_plus_pass));
  else
    sum_plus_pass = new RooAddPdf("sum_plus_pass","Full Model pdf plus pass",RooArgList(*signal_plus_pass,*chebbkg_plus_pass),RooArgList(*signal_yield_plus_pass,*background_yield_plus_pass));
  if(myconfig->get_bkgmodel() == 0)
    sum_plus_fail = new RooAddPdf("sum_plus_fail","Full Model pdf plus fail",RooArgList(*signal_plus_fail,*exponential_plus_fail),RooArgList(*signal_yield_plus_fail,*background_yield_plus_fail));
  else if(myconfig->get_shared_params() == 2)
    sum_plus_fail = new RooAddPdf("sum_plus_fail","Full Model pdf plus fail",RooArgList(*Jpsi_gauss1_plus,*chebbkg_plus_fail),RooArgList(*signal_yield_plus_fail,*background_yield_plus_fail));
  else
    sum_plus_fail = new RooAddPdf("sum_plus_fail","Full Model pdf plus fail",RooArgList(*signal_plus_fail,*chebbkg_plus_fail),RooArgList(*signal_yield_plus_fail,*background_yield_plus_fail));
  if(myconfig->get_bkgmodel() == 0)
    sum_minus_pass = new RooAddPdf("sum_minus_pass","Full Model pdf minus pass",RooArgList(*signal_minus_pass,*exponential_minus_pass),RooArgList(*signal_yield_minus_pass,*background_yield_minus_pass));
  else
    sum_minus_pass = new RooAddPdf("sum_minus_pass","Full Model pdf minus pass",RooArgList(*signal_minus_pass,*chebbkg_minus_pass),RooArgList(*signal_yield_minus_pass,*background_yield_minus_pass));
  if(myconfig->get_bkgmodel() == 0)
    sum_minus_fail = new RooAddPdf("sum_minus_fail","Full Model pdf minus fail",RooArgList(*signal_minus_fail,*exponential_minus_fail),RooArgList(*signal_yield_minus_fail,*background_yield_minus_fail));
  else if(myconfig->get_shared_params() == 2)
    sum_minus_fail = new RooAddPdf("sum_minus_fail","Full Model pdf minus fail",RooArgList(*Jpsi_gauss1_minus,*chebbkg_minus_fail),RooArgList(*signal_yield_minus_fail,*background_yield_minus_fail));
  else
    sum_minus_fail = new RooAddPdf("sum_minus_fail","Full Model pdf minus fail",RooArgList(*signal_minus_fail,*chebbkg_minus_fail),RooArgList(*signal_yield_minus_fail,*background_yield_minus_fail));

  RooCategory *Match = new RooCategory("Match","Match") ;
  Match->defineType("plus_pass") ;
  Match->defineType("plus_fail") ;
  Match->defineType("minus_pass");
  Match->defineType("minus_fail");

  // Construct combined dataset in (Jpsi_M,Match)
  RooDataHist *combData = new RooDataHist("combData","combined data",*Jpsi_M,Index(*Match),
                                          Import("plus_pass",*plusdatapass),Import("plus_fail",*plusdatafail),Import("minus_pass",*minusdatapass),Import("minus_fail",*minusdatafail));

  // Construct a simultaneous pdf using category Charge as index
  RooSimultaneous *simPdf = new RooSimultaneous("simPdf","simultaneous pdf",*Match) ;

  // Associate the pass and fail models with the match category
  simPdf->addPdf(*sum_plus_pass,"plus_pass") ;
  simPdf->addPdf(*sum_plus_fail,"plus_fail") ;
  simPdf->addPdf(*sum_minus_pass,"minus_pass");
  simPdf->addPdf(*sum_minus_fail,"minus_fail");

  if (myconfig->get_verbosity() > 0)simPdf->fitTo(*combData,SumW2Error(erroroption),InitialHesse(true),Extended(true),Minos(RooArgSet(*asymmetry,*average_eff)));
  else simPdf->fitTo(*combData,SumW2Error(erroroption),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));

  double tracking_asymmetry = asymmetry->getVal();
  double unc_asymmetry = asymmetry->getError();

  float pavetextsize = 0.04;
  float labeltextsize = 0.05;

  TCanvas* c1 = new TCanvas("c1","Fit pass",10,10,1280,2048) ;
  float padsplitmargin = 1e-6;
  float Left_margin  = 0.15;
  float Right_margin = 0.005;
  float Bigy_margin = 0.13;
  float padsplit = 0.25;//height of the pull histo
  c1->Divide(1,2,padsplitmargin,padsplitmargin);
  c1->cd(1);

  //Must control the margins in the individual pads and not the subpads of the canvas, because of the axis labels
  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  TPad *pad1 = new TPad("pad1","pad1",padsplitmargin,padsplit,1-padsplitmargin,1-padsplitmargin);
  pad1->Draw();
  pad1->cd();
  pad1->SetBorderSize(0);
  pad1->SetLeftMargin(Left_margin);
  pad1->SetRightMargin(Right_margin);
  pad1->SetTopMargin(Bigy_margin/(1-padsplit));
  pad1->SetBottomMargin(padsplitmargin);
  pad1->SetTickx();
  pad1->SetTicky();
  //pad1->SetLogy();
  c1->Update();

  // --- Plot ---
  RooPlot* Jpsi_Mframeplus = Jpsi_M->frame(Title(" ")) ;

  combData->plotOn(Jpsi_Mframeplus,Cut("Match==Match::plus_pass"),Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1));
  simPdf->plotOn(Jpsi_Mframeplus,Slice(*Match,"plus_pass"),ProjWData(*Match,*combData,true),Components("chebbkg_plus_pass"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines());
  simPdf->plotOn(Jpsi_Mframeplus,Slice(*Match,"plus_pass"),ProjWData(*Match,*combData,true),Components("chebbkg_plus_pass"),LineStyle(kDashed),LineWidth(2),LineColor(4));
  simPdf->plotOn(Jpsi_Mframeplus,Slice(*Match,"plus_pass"),ProjWData(*Match,*combData,true),LineWidth(2),LineColor(4));
  RooHist* ptpull = Jpsi_Mframeplus->pullHist() ;
  int bins_gt_threshold = 0;
  if(simple_fits == false){
    double *plus_matched_pulls = ptpull->GetY();
    for(int ij = 0; ij < myconfig->get_fitbinning(); ij++){
      if(plus_matched_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in plusdatapass \033[0m" << endl << endl;
      //terminate();
    }
  }
  combData->plotOn(Jpsi_Mframeplus,Cut("Match==Match::plus_pass"),Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1));


  TPaveText *ptplus = new TPaveText(0.64,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),1-pad1->GetRightMargin()-0.03,1-pad1->GetTopMargin()-0.03,"BRNDC");
  ptplus->SetName("paramBox");
  ptplus->SetBorderSize(0);ptplus->SetFillColor(kWhite);ptplus->SetTextAlign(12);ptplus->SetTextSize(pavetextsize/(1-padsplit));ptplus->SetTextFont(42);ptplus->SetFillStyle(0);
  ptplus->AddText(myconfig->get_methodlabel());
  ptplus->AddText("matched #mu^{+} probe legs");
  temp.Form( "N^{total}_{J/#psi,eff} = %.0f #pm %.0f", total_signal_yield_plus->getVal(),total_signal_yield_plus->getError());
  ptplus->AddText(temp);
  Jpsi_Mframeplus->addObject(ptplus) ;

  TPaveText *ptparbinp = new TPaveText(pad1->GetLeftMargin()+0.03,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),0.4,1-pad1->GetTopMargin()-0.03,"BRNDC");
  ptparbinp->SetBorderSize(0);ptparbinp->SetFillColor(kWhite);ptparbinp->SetTextAlign(12);ptparbinp->SetTextSize(pavetextsize/(1-padsplit));ptparbinp->SetTextFont(42);ptparbinp->SetFillStyle(0);
  if(myconfig->get_dim1().CompareTo("ETA") == 0)temp.Form("%s #in [%g,%g]",myconfig->get_dim1_label().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi());
  else temp.Form("%s #in [%g,%g] GeV",myconfig->get_dim1_label().Data(),myconfig->get_dim1_binlo()/1000,myconfig->get_dim1_binhi()/1000);
  if(!myconfig->is_1D()){
    if(myconfig->get_dim2().CompareTo("ETA") == 0)temp.Form("%s #in [%g,%g]",myconfig->get_dim2_label().Data(),myconfig->get_dim2_binlo(),myconfig->get_dim2_binhi());
    else temp.Form("%s #in [%g,%g] GeV",myconfig->get_dim2_label().Data(),myconfig->get_dim2_binlo()/1000,myconfig->get_dim2_binhi()/1000);
  }
  ptparbinp->AddText(temp);
  //ptparbinp->AddText(myconfig->get_channellabel());
  ptparbinp->AddText(myconfig->get_samplelabel());
  Jpsi_Mframeplus->addObject(ptparbinp);

  TPaveText *pt3 = new TPaveText(pad1->GetLeftMargin(),1.002-pad1->GetTopMargin(),pad1->GetLeftMargin()+0.1,1,"BRNDC");
  pt3->SetBorderSize(0);pt3->SetFillColor(kWhite);pt3->SetFillStyle(1001);pt3->SetTextColor(kWhite);
  pt3->AddText(" ");
  Jpsi_Mframeplus->addObject(pt3);

  TPaveText *pt4 = new TPaveText(1.002-pad1->GetRightMargin(),pad1->GetBottomMargin(),1,pad1->GetBottomMargin()+0.05,"BRNDC");
  pt4->SetBorderSize(0);pt4->SetFillColor(kWhite);pt4->SetFillStyle(1001);pt4->SetTextColor(kWhite);
  pt4->AddText(" ");
  Jpsi_Mframeplus->addObject(pt4);

  Jpsi_Mframeplus->GetXaxis()->SetNdivisions(305) ;
  char yaxislabel[100];
  if(Jpsi_Mframeplus->GetMaximum() > 1e+3 && Jpsi_Mframeplus->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);
  //sprintf(yaxislabel,"weighted Events/%.0f MeV",binsize);
  Jpsi_Mframeplus->GetYaxis()->SetTitle(yaxislabel) ;
  Jpsi_Mframeplus->GetYaxis()->SetTitleSize(labeltextsize/(1-padsplit)) ;
  Jpsi_Mframeplus->GetYaxis()->SetLabelSize(labeltextsize/(1-padsplit)) ;
  Jpsi_Mframeplus->GetYaxis()->SetTitleOffset(1.2*(1-padsplit)) ;
  Jpsi_Mframeplus->GetYaxis()->SetRangeUser(0.1,Jpsi_Mframeplus->GetMaximum());
  //Jpsi_Mframeplus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeplus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeplus->GetMaximum())));
  Jpsi_Mframeplus->Draw() ;


  TCanvas* c2 = new TCanvas("c2","Fit fail",10,10,1280,2048) ;
  c2->Divide(1,2,padsplitmargin,padsplitmargin);
  c2->cd(1);

  //Must control the margins in the individual pads and not the subpads of the canvas, because of the axis labels
  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  TPad *pad5 = (TPad*)pad1->Clone("pad5");
  pad5->Draw();
  pad5->cd();
  c2->Update();

  // --- Plot ---
  RooPlot* Jpsi_Mframeplus_fail = Jpsi_M->frame(Title(" ")) ;

  combData->plotOn(Jpsi_Mframeplus_fail,Cut("Match==Match::plus_fail"),Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1));
  simPdf->plotOn(Jpsi_Mframeplus_fail,Slice(*Match,"plus_fail"),ProjWData(*Match,*combData,true),Components("chebbkg_plus_fail"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines());
  simPdf->plotOn(Jpsi_Mframeplus_fail,Slice(*Match,"plus_fail"),ProjWData(*Match,*combData,true),Components("chebbkg_plus_fail"),LineStyle(kDashed),LineWidth(2),LineColor(4));
  simPdf->plotOn(Jpsi_Mframeplus_fail,Slice(*Match,"plus_fail"),ProjWData(*Match,*combData,true),LineWidth(2),LineColor(4));

  RooHist* prpull = Jpsi_Mframeplus_fail->pullHist() ;
  if(simple_fits == false){
    double *plus_unmatched_pulls = prpull->GetY();
    bins_gt_threshold = 0;
    for(int ij = 0; ij < myconfig->get_fitbinning(); ij++){
      if(plus_unmatched_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in plusdatafail \033[0m" << endl << endl;
      //terminate();
    }
  }
  //plusdatafail->plotOn(Jpsi_Mframeplus_fail,Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1)) ;
  combData->plotOn(Jpsi_Mframeplus_fail,Cut("Match==Match::plus_fail"),Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1));

  TPaveText *prplus = new TPaveText(0.60,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),1-pad1->GetRightMargin()-0.03,1-pad1->GetTopMargin()-0.03,"BRNDC");
  prplus->SetName("paramBox");
  prplus->SetBorderSize(0);prplus->SetFillColor(kWhite);prplus->SetTextAlign(12);prplus->SetTextSize(pavetextsize/(1-padsplit));prplus->SetTextFont(42);prplus->SetFillStyle(0);
  prplus->AddText(myconfig->get_methodlabel());
  prplus->AddText("unmatched #mu^{+} probe legs");
  temp.Form( "<#varepsilon> = %.2f #pm %.2f %%", 100*average_eff->getVal(),100*average_eff->getError());
  prplus->AddText(temp);
  Jpsi_Mframeplus_fail->addObject(prplus) ;
  //add the other PaveTexts
  Jpsi_Mframeplus_fail->addObject(ptparbinp) ;
  Jpsi_Mframeplus_fail->addObject(pt3);
  Jpsi_Mframeplus_fail->addObject(pt4);

  if(Jpsi_Mframeplus_fail->GetMaximum() > 1e+3 && Jpsi_Mframeplus->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus_fail->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus_fail->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);

  Jpsi_Mframeplus_fail->GetXaxis()->SetNdivisions(305) ;
  Jpsi_Mframeplus_fail->GetYaxis()->SetTitle(yaxislabel) ;
  Jpsi_Mframeplus_fail->GetYaxis()->SetTitleSize(labeltextsize/(1-padsplit)) ;
  Jpsi_Mframeplus_fail->GetYaxis()->SetLabelSize(labeltextsize/(1-padsplit)) ;
  Jpsi_Mframeplus_fail->GetYaxis()->SetTitleOffset(1.2*(1-padsplit)) ;
  //if(meth == 1)Jpsi_Mframeplus_fail->GetYaxis()->SetRangeUser(0.1,Jpsi_Mframeplus_fail->GetMaximum());
  //else
  Jpsi_Mframeplus_fail->GetYaxis()->SetRangeUser(0.1,1.4*Jpsi_Mframeplus_fail->GetMaximum());
  //Jpsi_Mframeplus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeplus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeplus->GetMaximum())));
  Jpsi_Mframeplus_fail->Draw() ;

  c2->cd(1);

  TPad *pad2 = new TPad("pad2","pad2",padsplitmargin,padsplitmargin,1-padsplitmargin,padsplit);
  pad2->SetBorderSize(0);
  pad2->SetLeftMargin(Left_margin);
  pad2->SetRightMargin(Right_margin);
  pad2->SetTopMargin(padsplitmargin);
  pad2->SetBottomMargin(padsplitmargin);
  pad2->SetTickx();
  pad2->SetTicky();

  RooPlot* ppullframe = Jpsi_M->frame(Title(" ")) ;

  Color_t pullcolor = kAzure+3;
  ptpull->SetLineColor(pullcolor);
  ptpull->SetMarkerColor(pullcolor);
  prpull->SetLineColor(pullcolor);
  prpull->SetMarkerColor(pullcolor);

  ppullframe->GetXaxis()->SetTickLength(Jpsi_Mframeplus->GetXaxis()->GetTickLength()/(padsplit/(1-padsplit)));
  ppullframe->GetXaxis()->SetNdivisions(305);
  //h18->GetYaxis()->SetTickLength(Jpsi_Mframeplus->GetYaxis()->GetTickLength()/(padsplit/(1-padsplit)));
  ppullframe->GetYaxis()->SetLabelSize(Jpsi_Mframeplus->GetYaxis()->GetLabelSize()/(padsplit/(1-padsplit)));
  ppullframe->GetYaxis()->SetTitleSize(Jpsi_Mframeplus->GetYaxis()->GetTitleSize()/(padsplit/(1-padsplit)));
  ppullframe->GetYaxis()->SetNdivisions(205);
  ppullframe->GetYaxis()->SetTitleOffset(1.2*padsplit);
  ppullframe->GetYaxis()->SetTitle("Pull ");

  TLine* zeroline = new TLine(myconfig->get_Jpsi_M_min(),0,myconfig->get_Jpsi_M_max(),0);
  zeroline->SetLineStyle(2);
  zeroline->SetLineColor(1);
  TLine* m3sigma = new TLine(myconfig->get_Jpsi_M_min(),-3,myconfig->get_Jpsi_M_max(),-3);
  m3sigma->SetLineStyle(2);
  m3sigma->SetLineColor(1);
  TLine* p3sigma = new TLine(myconfig->get_Jpsi_M_min(),3,myconfig->get_Jpsi_M_max(),3);
  p3sigma->SetLineStyle(2);
  p3sigma->SetLineColor(1);
  ppullframe->addObject(zeroline,"L") ;
  ppullframe->addObject(m3sigma,"L") ;
  ppullframe->addObject(p3sigma,"L") ;
  TPaveText *pt5 = new TPaveText(1.002-pad2->GetRightMargin(),pad2->GetBottomMargin(),1,pad2->GetBottomMargin()+0.05,"BRNDC");
  pt5->SetBorderSize(0);pt5->SetFillColor(kWhite);pt5->SetFillStyle(1001);pt5->SetTextColor(kWhite);
  pt5->AddText(" ");
  ppullframe->addObject(pt5);
  ppullframe->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);

  TPad *pad6 = (TPad*)pad2->Clone("pad6");
  pad6->Draw();
  pad6->cd();
  RooPlot* ppullframe_fail = (RooPlot*)ppullframe->Clone("ppullframe_fail");
  ppullframe_fail->addPlotable(prpull,"P") ;
  ppullframe_fail->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);
  ppullframe_fail->Draw();
  c2->Update();

  c1->cd(1);
  pad2->Draw();
  pad2->cd();
  ppullframe->addPlotable(ptpull,"P") ;
  ppullframe->Draw();
  c1->Update();

  c1->cd(2);

  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  float lowerpadsplit = padsplit+(Bigy_margin/(1-padsplit));

  TPad *pad3 = new TPad("pad3","pad3",padsplitmargin,lowerpadsplit,1-padsplitmargin,1-padsplitmargin);
  pad3->Draw();
  pad3->cd();
  pad3->SetBorderSize(0);
  pad3->SetLeftMargin(Left_margin);
  pad3->SetRightMargin(Right_margin);
  pad3->SetTopMargin(padsplitmargin);
  pad3->SetBottomMargin(padsplitmargin);
  pad3->SetTickx();
  pad3->SetTicky();
  //pad3->SetLogy();
  c1->Update();

  // --- Plot ---
  RooPlot* Jpsi_Mframeminus = Jpsi_M->frame(Title(" ")) ;

  combData->plotOn(Jpsi_Mframeminus,Cut("Match==Match::minus_pass"),Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1));
  simPdf->plotOn(Jpsi_Mframeminus,Slice(*Match,"minus_pass"),ProjWData(*Match,*combData,true),Components("chebbkg_minus_pass"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines());
  simPdf->plotOn(Jpsi_Mframeminus,Slice(*Match,"minus_pass"),ProjWData(*Match,*combData,true),Components("chebbkg_minus_pass"),LineStyle(kDashed),LineWidth(2),LineColor(4));
  simPdf->plotOn(Jpsi_Mframeminus,Slice(*Match,"minus_pass"),ProjWData(*Match,*combData,true),LineWidth(2),LineColor(4));
  RooHist* mtpull = Jpsi_Mframeminus->pullHist() ;
  if(simple_fits == false){
    bins_gt_threshold = 0;
    double *minus_matched_pulls = mtpull->GetY();
    for(int ij = 0; ij < myconfig->get_fitbinning(); ij++){
      if(minus_matched_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in minusdatapass\033[0m" << endl << endl;
      //return NULL;
    }
  }
  combData->plotOn(Jpsi_Mframeminus,Cut("Match==Match::minus_pass"),Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1));

  TPaveText *ptminus = new TPaveText(0.64,1-gPad->GetTopMargin()-0.03-(4*pavetextsize/(1-lowerpadsplit)),1-gPad->GetRightMargin()-0.03,1-gPad->GetTopMargin()-0.03,"BRNDC");
  ptminus->SetName("paramBox");
  ptminus->SetBorderSize(0);ptminus->SetFillColor(kWhite);ptminus->SetTextAlign(12);ptminus->SetTextSize(pavetextsize/(1-lowerpadsplit));ptminus->SetTextFont(42);ptminus->SetFillStyle(0);
  ptminus->AddText(myconfig->get_methodlabel());
  ptminus->AddText("matched #mu^{-} probe legs");
  temp.Form("N^{total}_{J/#psi,eff} = %.0f #pm %.0f", total_signal_yield_minus->getVal(),total_signal_yield_minus->getError());
  ptminus->AddText(temp);
  Jpsi_Mframeminus->addObject(ptminus) ;

  TPaveText *ptparbin = new TPaveText(gPad->GetLeftMargin()+0.03,1-gPad->GetTopMargin()-0.03-(4*pavetextsize/(1-lowerpadsplit)),0.4,1-gPad->GetTopMargin()-0.03,"BRNDC");
  ptparbin->SetBorderSize(0);ptparbin->SetFillColor(kWhite);ptparbin->SetTextAlign(12);ptparbin->SetTextSize(pavetextsize/(1-lowerpadsplit));ptparbin->SetTextFont(42);ptparbin->SetFillStyle(0);
  if(myconfig->get_dim1().CompareTo("ETA") == 0)temp.Form("%s #in [%g,%g]",myconfig->get_dim1_label().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi());
  else temp.Form("%s #in [%g,%g] GeV",myconfig->get_dim1_label().Data(),myconfig->get_dim1_binlo()/1000,myconfig->get_dim1_binhi()/1000);
  if(!myconfig->is_1D()){
    if(myconfig->get_dim2().CompareTo("ETA") == 0)temp.Form("%s #in [%g,%g]",myconfig->get_dim2_label().Data(),myconfig->get_dim2_binlo(),myconfig->get_dim2_binhi());
    else temp.Form("%s #in [%g,%g] GeV",myconfig->get_dim2_label().Data(),myconfig->get_dim2_binlo()/1000,myconfig->get_dim2_binhi()/1000);
  }
  ptparbin->AddText(temp);
  //ptparbin->AddText(myconfig->get_channellabel());
  ptparbin->AddText(myconfig->get_samplelabel());
  Jpsi_Mframeminus->addObject(ptparbin);

  TPaveText *pt6 = new TPaveText(1.002-pad3->GetRightMargin(),pad3->GetBottomMargin(),1,pad3->GetBottomMargin()+0.05,"BRNDC");
  pt6->SetBorderSize(0);pt6->SetFillColor(kWhite);pt6->SetFillStyle(1001);pt6->SetTextColor(kWhite);
  pt6->AddText(" ");
  Jpsi_Mframeminus->addObject(pt6);

  Jpsi_Mframeminus->GetXaxis()->SetNdivisions(305) ;
  if(Jpsi_Mframeminus->GetMaximum() > 1e+3 && Jpsi_Mframeminus->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeminus->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeplus->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);
  //Jpsi_Mframeminus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeminus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeminus->GetMaximum())));
  //sprintf(yaxislabel,"weighted Events/%.0f MeV",binsize);
  Jpsi_Mframeminus->GetYaxis()->SetTitle(yaxislabel) ;
  Jpsi_Mframeminus->GetYaxis()->SetTitleSize(labeltextsize/(1-lowerpadsplit)) ;
  Jpsi_Mframeminus->GetYaxis()->SetLabelSize(labeltextsize/(1-lowerpadsplit)) ;
  Jpsi_Mframeminus->GetYaxis()->SetTitleOffset(1.2*(1-lowerpadsplit)) ;
  Jpsi_Mframeminus->GetYaxis()->SetRangeUser(0.1,Jpsi_Mframeminus->GetMaximum());
  Jpsi_Mframeminus->Draw() ;

  c2->cd(2);

  //Must control the margins in the individual pads and not the subpads of the canvas, because of the axis labels
  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  TPad *pad7 = (TPad*)pad3->Clone("pad7");
  pad7->Draw();
  pad7->cd();
  c2->Update();

  RooPlot* Jpsi_Mframeminus_fail = Jpsi_M->frame(Title(" ")) ;

  combData->plotOn(Jpsi_Mframeminus_fail,Cut("Match==Match::minus_fail"),Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1));
  simPdf->plotOn(Jpsi_Mframeminus_fail,Slice(*Match,"minus_fail"),ProjWData(*Match,*combData,true),Components("chebbkg_minus_fail"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines());
  simPdf->plotOn(Jpsi_Mframeminus_fail,Slice(*Match,"minus_fail"),ProjWData(*Match,*combData,true),Components("chebbkg_minus_fail"),LineStyle(kDashed),LineWidth(2),LineColor(4));
  simPdf->plotOn(Jpsi_Mframeminus_fail,Slice(*Match,"minus_fail"),ProjWData(*Match,*combData,true),LineWidth(2),LineColor(4));
  RooHist* mrpull = Jpsi_Mframeminus_fail->pullHist() ;
  if(simple_fits == false){
    double *minus_unmatched_pulls = mrpull->GetY();
    bins_gt_threshold = 0;
    for(int ij = 0; ij < myconfig->get_fitbinning(); ij++){
      if(minus_unmatched_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in minusdatafail \033[0m" << endl << endl;
      //terminate();
    }
  }
  combData->plotOn(Jpsi_Mframeminus_fail,Cut("Match==Match::minus_fail"),Binning(myconfig->get_fitbinning()),MarkerColor(1),LineColor(1));

  TPaveText *prminus = new TPaveText(0.6,1-gPad->GetTopMargin()-0.03-(4*pavetextsize/(1-lowerpadsplit)),1-gPad->GetRightMargin()-0.03,1-gPad->GetTopMargin()-0.03,"BRNDC");
  prminus->SetName("paramBox");
  prminus->SetBorderSize(0);prminus->SetFillColor(kWhite);prminus->SetTextAlign(12);prminus->SetTextSize(pavetextsize/(1-lowerpadsplit));prminus->SetTextFont(42);prminus->SetFillStyle(0);
  prminus->AddText(myconfig->get_methodlabel());
  prminus->AddText("unmatched #mu^{-} probe legs");
  temp.Form( "A = %.2f #pm %.2f %%", 100*tracking_asymmetry,100*unc_asymmetry);
  prminus->AddText(temp);
  Jpsi_Mframeminus_fail->addObject(prminus) ;
  //add the other PaveTexts
  Jpsi_Mframeminus_fail->addObject(ptparbin) ;
  Jpsi_Mframeminus_fail->addObject(pt6);

  Jpsi_Mframeminus_fail->GetXaxis()->SetNdivisions(305) ;
  if(Jpsi_Mframeminus_fail->GetMaximum() > 1e+3 && Jpsi_Mframeminus_fail->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeminus_fail->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(Jpsi_Mframeminus_fail->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);
  //Jpsi_Mframeminus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeminus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeminus->GetMaximum())));
  //sprintf(yaxislabel,"weighted Events/%.0f MeV",binsize);
  Jpsi_Mframeminus_fail->GetYaxis()->SetTitle(yaxislabel) ;
  Jpsi_Mframeminus_fail->GetYaxis()->SetTitleSize(labeltextsize/(1-lowerpadsplit)) ;
  Jpsi_Mframeminus_fail->GetYaxis()->SetLabelSize(labeltextsize/(1-lowerpadsplit)) ;
  Jpsi_Mframeminus_fail->GetYaxis()->SetTitleOffset(1.2*(1-lowerpadsplit)) ;
  //if(meth == 1)Jpsi_Mframeminus_fail->GetYaxis()->SetRangeUser(0.1,Jpsi_Mframeminus_fail->GetMaximum());
  //else
  Jpsi_Mframeminus_fail->GetYaxis()->SetRangeUser(0.1,1.4*Jpsi_Mframeminus_fail->GetMaximum());
  Jpsi_Mframeminus_fail->Draw() ;

  c2->Update();
  c2->cd(2);

  TPad *pad4 = new TPad("pad4","pad4",padsplitmargin,padsplitmargin,1-padsplitmargin,lowerpadsplit);
  pad4->SetBorderSize(0);
  pad4->SetLeftMargin(Left_margin);
  pad4->SetRightMargin(Right_margin);
  pad4->SetTopMargin(padsplitmargin);
  pad4->SetBottomMargin(Bigy_margin/lowerpadsplit);
  pad4->SetTickx();
  pad4->SetTicky();

  mtpull->SetLineColor(pullcolor);
  mtpull->SetMarkerColor(pullcolor);
  mrpull->SetLineColor(pullcolor);
  mrpull->SetMarkerColor(pullcolor);

  RooPlot* mpullframe = Jpsi_M->frame(Title(" ")) ;

  mpullframe->GetXaxis()->SetTitleOffset(1.05);
  mpullframe->GetXaxis()->SetTitleSize(labeltextsize/lowerpadsplit);
  mpullframe->GetXaxis()->SetLabelSize(labeltextsize/lowerpadsplit);
  mpullframe->GetXaxis()->SetTitle("M_{inv}(#mu^{+}#mu^{-}) (GeV)");
  mpullframe->GetXaxis()->SetTickLength(Jpsi_Mframeplus->GetXaxis()->GetTickLength()/(lowerpadsplit/(1-lowerpadsplit)));
  //h18->GetYaxis()->SetTickLength(Jpsi_Mframeplus->GetYaxis()->GetTickLength()/(lowerpadsplit/(1-lowerpadsplit)));
  mpullframe->GetYaxis()->SetLabelSize(labeltextsize/lowerpadsplit);
  mpullframe->GetYaxis()->SetTitleSize(labeltextsize/lowerpadsplit);
  mpullframe->GetXaxis()->SetNdivisions(305);
  mpullframe->GetYaxis()->SetNdivisions(205);
  mpullframe->GetYaxis()->SetTitleOffset(1.2*lowerpadsplit);
  mpullframe->GetYaxis()->SetTitle("Pull ");

  TLine* mzeroline = new TLine(myconfig->get_Jpsi_M_min(),0,myconfig->get_Jpsi_M_max(),0);
  mzeroline->SetLineStyle(2);
  mzeroline->SetLineColor(1);
  TLine* mm3sigma = new TLine(myconfig->get_Jpsi_M_min(),-3,myconfig->get_Jpsi_M_max(),-3);
  mm3sigma->SetLineStyle(2);
  mm3sigma->SetLineColor(1);
  TLine* mp3sigma = new TLine(myconfig->get_Jpsi_M_min(),3,myconfig->get_Jpsi_M_max(),3);
  mp3sigma->SetLineStyle(2);
  mp3sigma->SetLineColor(1);
  mpullframe->addObject(mzeroline,"L") ;
  mpullframe->addObject(mm3sigma,"L") ;
  mpullframe->addObject(mp3sigma,"L") ;
  TPaveText *pt7 = new TPaveText(1.002-pad4->GetRightMargin(),pad4->GetBottomMargin(),1,pad4->GetBottomMargin()+0.05,"BRNDC");
  pt7->SetBorderSize(0);pt7->SetFillColor(kWhite);pt7->SetFillStyle(1001);pt7->SetTextColor(kWhite);
  pt7->AddText(" ");
  mpullframe->addObject(pt7);
  mpullframe->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);

  TPad *pad8 = (TPad*)pad4->Clone("pad8");
  pad8->Draw();
  pad8->cd();
  RooPlot* mpullframe_fail = (RooPlot*)mpullframe->Clone("mpullframe_fail");
  mpullframe_fail->addPlotable(mrpull,"P") ;
  mpullframe_fail->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);
  mpullframe_fail->Draw();
  c2->Update();

  c1->cd(2);
  pad4->Draw();
  pad4->cd();
  mpullframe->addPlotable(mtpull,"P") ;
  mpullframe->Draw();
  c1->Update();

  temp = myconfig->get_dumpdir()+"Fits";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);

  if(myconfig->is_1D())temp.Form("%sFits/raw_%s_%s_%s_%g_%g_pass.pdf",myconfig->get_dumpdir().Data(),
                                 myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi());
  else temp.Form("%sFits/raw_%s_%s_%s_%g_%g_%g_%g_pass.pdf",myconfig->get_dumpdir().Data(),
                 myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi(),myconfig->get_dim2_binlo(),myconfig->get_dim2_binhi());
  if(!myconfig->is_debug()) c1->SaveAs(temp);


  if(myconfig->is_1D())temp.Form("%sFits/raw_%s_%s_%s_%g_%g_fail.pdf",myconfig->get_dumpdir().Data(),
                                 myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi());
  else temp.Form("%sFits/raw_%s_%s_%s_%g_%g_%g_%g_fail.pdf",myconfig->get_dumpdir().Data(),
                 myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi(),myconfig->get_dim2_binlo(),myconfig->get_dim2_binhi());
  if(!myconfig->is_debug()) c2->SaveAs(temp);

  delete c1;
  delete c2;

  //this is the important stuff
  double* the_asymmetries;
  the_asymmetries = new double[2];
  the_asymmetries[0] = tracking_asymmetry;
  the_asymmetries[1] = unc_asymmetry;
  cout << "\033[0;33m Raw Asymmetry = ( " << 100*the_asymmetries[0] << " +- "<< 100*the_asymmetries[1] << " ) %\033[0m" << endl;


  delete ptplus;
  delete ptparbinp;
  delete ptminus;
  delete ptparbin;
  delete prplus;
  delete prminus;
  delete zeroline;
  delete m3sigma;
  delete p3sigma;
  delete mzeroline;
  delete mm3sigma;
  delete mp3sigma;
  delete pt3;
  delete pt4;
  delete pt5;
  delete pt6;
  delete pt7;
  return the_asymmetries;
}

bool convcheck(RooFitResult *fitres, vector<RooRealVar*> vars_for_convcheck, vector<RooRealVar*> posis, RooAbsPdf *model, RooDataHist *data, configuration* myconfig){

  /*
Checks if the fitresult of model to data given the parameters vars_for_convcheck and the parameters of special interest posis converged and refits automatically in case not.
The non-convergence considered here has 3 main symptomps: either someting went wrong in the calculation of error matrices/covariances, the fit din't converge at all, or the parameters reached their limits.

The first 2 problems (case A) can be tackled by simply refitting the data with different, more careful options (minimizers) of the fit strategy.
To be more accurate: if something went wrong in the calculation of error matrices/covariances, it is likely to be fixed by a more careful treatment of the errors by the MINOS routine,
which takes rather long and sometimes has problems with fits where the start parameters are not chosen carefully. In our case, the starting parameters are the result from the previous fit, so this should be fine.
If the fit didn't converge at all, new start paramers should be chosen, which is automatically the case here (it is assumed that the non-converged fit is at least a bit better than the model with only the start paraeters).
The third problem (parameters reached limits, case B) can be tackled by defining a new range for these parameters. Currently all other parameters are also constrained to make the fit run faster.

It turned out that an alternating query of the cases gives best performance. I.e. check A first, if it failed, check then B first, if that failed, check again A first.

Currenty 3 levels of checks and refits are implemented (maybe this could be solved more general, i.e. give the number of refits as an option and be able to refit forever and a day ;) )
    */

  double maxedm = 1e-3;
  if ((fitres -> status() != 0  && fitres->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
    cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: this fit did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;
    //Refit I case A: status or uncertainties failed. This can be severe. Try to refit with different options
    if ((fitres -> status() != 0 && fitres->edm() > maxedm)|| errcheck(posis) == false){
      cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options\033[0m" << endl << endl;
      RooFitResult *fitres_refitstat;
      if(myconfig->get_verbosity() > 0)fitres_refitstat = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),InitialHesse(true),Extended(true),NumCPU(cpus_for_fit));
      else fitres_refitstat = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),InitialHesse(true),Extended(true),NumCPU(cpus_for_fit),PrintLevel(-1),Verbose(false));
      if ((fitres_refitstat -> status() != 0  && fitres_refitstat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
        cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 1st attempt to re-fit (statuscheck or bad uncertainties) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;
        //Refit II case AB: parameter reached limit. Will be solved by manipulation of the parameterranges
        if (paratlimit(vars_for_convcheck) == false){
          cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
          setnewparranges(vars_for_convcheck);
          RooFitResult *fitres_refitstatlim;
          if(myconfig->get_verbosity() > 0)fitres_refitstatlim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;
          else fitres_refitstatlim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
          if ((fitres_refitstatlim -> status() != 0 && fitres_refitstatlim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
            cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 2nd attempt to re-fit (parameter reached limits) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;

            //Refit III case ABA: status or uncertainties failed. This can be severe. Try to refit with different options
            if ((fitres_refitstatlim -> status() != 0 && fitres_refitstatlim->edm() > maxedm) || errcheck(posis) == false){
              cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options\033[0m" << endl << endl;
              RooFitResult *fitres_refitstatlimstat;
              if(myconfig->get_verbosity() > 0)fitres_refitstatlimstat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),SumW2Error(erroroption),NumCPU(cpus_for_fit));
              else fitres_refitstatlimstat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refitstatlimstat -> status() != 0 && fitres_refitstatlimstat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
            //Refit III case ABB: parameter reached limit. Will be solved by manipulation of the parameterranges
            else{
              cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
              setnewparranges(vars_for_convcheck);
              RooFitResult *fitres_refitstat2lim;
              if(myconfig->get_verbosity() > 0) fitres_refitstat2lim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;
              else fitres_refitstat2lim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
              if ((fitres_refitstat2lim -> status() != 0 && fitres_refitstat2lim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
          }
          return true;
        }
        //Refit II case AA: status or uncertainties failed. This can be severe. Try to refit with different options
        else{
          cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options. 2nd try...\033[0m" << endl << endl;
          RooFitResult *fitres_refit2stat;
          if(myconfig->get_verbosity() > 0) fitres_refit2stat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),SumW2Error(erroroption),NumCPU(cpus_for_fit));
          else fitres_refit2stat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
          if ((fitres_refit2stat -> status() != 0  && fitres_refit2stat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
            cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 2nd attempt to re-fit (statuscheck or bad uncertainties) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;
            //Refit III case AAB: parameter reached limit. Will be solved by manipulation of the parameterranges
            if (paratlimit(vars_for_convcheck) == false){
              cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
              setnewparranges(vars_for_convcheck);
              RooFitResult *fitres_refit2statlim;
              if(myconfig->get_verbosity() > 0) fitres_refit2statlim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;
              else fitres_refit2statlim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
              if ((fitres_refit2statlim -> status() != 0 && fitres_refit2statlim->edm() > maxedm)|| errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
            //Refit III case AAA: status or uncertainties failed. This can be severe. Try to refit with different options
            else{
              cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options. Last try!\033[0m" << endl << endl;
              RooFitResult *fitres_refit3stat;
              if(myconfig->get_verbosity() > 0) fitres_refit3stat = model->fitTo(*data,Save(true),Minimizer("Minuit","migradimproved"),SumW2Error(erroroption),Strategy(2),NumCPU(cpus_for_fit));
              else fitres_refit3stat = model->fitTo(*data,Save(true),Minimizer("Minuit","migradimproved"),SumW2Error(erroroption),Strategy(2),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refit3stat -> status() != 0 && fitres_refit3stat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
          }
          return true;
        }
      }
      return true;
    }
    //Refit I case B: parameter reached limit. Will be solved by manipulation of the parameterranges
    else{
      cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
      setnewparranges(vars_for_convcheck);
      RooFitResult *fitres_refitlim;
      if(myconfig->get_verbosity() > 0) fitres_refitlim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;
      else fitres_refitlim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
      if ((fitres_refitlim -> status() != 0 && fitres_refitlim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
        cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 1st attempt to re-fit (parameter reached limits) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;

        //Refit II case BA: status or uncertainties failed. This can be severe. Try to refit with different options
        if ((fitres_refitlim -> status() != 0 && fitres_refitlim->edm() > maxedm) || errcheck(posis) == false){
          cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options. 2nd try...\033[0m" << endl << endl;
          RooFitResult *fitres_refitlimstat;
          if(myconfig->get_verbosity() > 0) fitres_refitlimstat = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),InitialHesse(true),Extended(true),NumCPU(cpus_for_fit));
          else fitres_refitlimstat = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),InitialHesse(true),Extended(true),NumCPU(cpus_for_fit),PrintLevel(-1),Verbose(false));
          if ((fitres_refitlimstat -> status() != 0 && fitres_refitlimstat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
            cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 2nd attempt to re-fit (statuscheck or bad uncertainties) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;

            //Refit III case BAB: parameter reached limit. Will be solved by manipulation of the parameterranges
            if (paratlimit(vars_for_convcheck) == false){
              cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
              setnewparranges(vars_for_convcheck);
              RooFitResult *fitres_refitlimstatlim;
              if(myconfig->get_verbosity() > 0) fitres_refitlimstatlim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;
              else fitres_refitlimstatlim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
              if ((fitres_refitlimstatlim -> status() != 0 && fitres_refitlimstatlim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
            //Refit III case BAA: status or uncertainties failed. This can be severe. Try to refit with different options
            else{
              cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options. Last try!\033[0m" << endl << endl;
              RooFitResult *fitres_refitlim2stat;
              if(myconfig->get_verbosity() > 0) fitres_refitlim2stat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),SumW2Error(erroroption),NumCPU(cpus_for_fit));
              else fitres_refitlim2stat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refitlim2stat -> status() != 0 && fitres_refitlim2stat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
          }
          return true;
        }
        //Refit II case BB: parameter reached limit. Will be solved by manipulation of the parameterranges
        else{
          cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
          setnewparranges(vars_for_convcheck);
          RooFitResult *fitres_refit2lim;
          if(myconfig->get_verbosity() > 0) fitres_refit2lim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit)) ;
          else fitres_refit2lim = model->fitTo(*data,Save(true),Strategy(2),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
          if ((fitres_refit2lim -> status() != 0 && fitres_refit2lim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
            cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 2nd attempt to re-fit (parameter reached limits) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;

            //Refit III case BBA: status or uncertainties failed. This can be severe. Try to refit with different options
            if ((fitres_refit2lim -> status() != 0 && fitres_refit2lim->edm() > maxedm) || errcheck(posis) == false){
              cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options\033[0m" << endl << endl;
              RooFitResult *fitres_refit2limstat;
              if(myconfig->get_verbosity() > 0) fitres_refit2limstat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),SumW2Error(erroroption),NumCPU(cpus_for_fit));
              else fitres_refit2limstat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refit2limstat -> status() != 0 && fitres_refit2limstat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
            //Refit III case BBB: parameter reached limit. Will be solved by manipulation of the parameterranges
            else{
              cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
              setnewparranges(vars_for_convcheck);
              RooFitResult *fitres_refit3lim;
              if(myconfig->get_verbosity() > 0) fitres_refit3lim = model->fitTo(*data,Save(true),Strategy(1),SumW2Error(erroroption),NumCPU(cpus_for_fit));
              else fitres_refit3lim = model->fitTo(*data,Save(true),Strategy(1),SumW2Error(erroroption),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refit3lim -> status() != 0 && fitres_refit3lim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
          }
          return true;
        }
      }
      return true;
    }
  }
  return true;
}

bool errcheck(vector<RooRealVar*> posis){
  cout << endl << "\033[0;34mChecking size of the uncertainties for the parameters of special interest...";

  // Rule of thumb could be Number of free parameters ...
  // Dealing with weighted spectra, so the error from SumW2Errors can get large
  double insaneerrors = 32.0;

  // Figure of merit is the Poissonian standard deviation sqrt(N)
  for (vector<RooRealVar*>::const_iterator siter = posis.begin() ;  posis.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    if( (iter->getError() > insaneerrors*sqrt(iter->getVal()) ) || (iter->getError() < (1/insaneerrors)*sqrt(iter->getVal()) ) ){
      cout << endl << "\033[0;31mParameter: "<< iter->GetName() << " has a unexpectedly high or low uncertainty ( Delta N/sqrt(N) = " << iter->getError()/sqrt(iter->getVal()) << " )\033[0m" << endl;
      return false;
    }
  }
  cout << " OK\033[0m" << endl;
  return true;
}

bool paratlimit(vector<RooRealVar*> vars){

  double convergencecriterion = 0.001;

  cout << endl << "Checking for convergence..." << endl<<endl;

  for (vector<RooRealVar*>::const_iterator siter = vars.begin() ;  vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    if(iter->isConstant())continue;
    cout << "Checking parameter "<< iter->GetName() << " in range [" << iter->getMin() + convergencecriterion*(iter->getMax() - iter->getMin()) << "," << iter->getMax() - convergencecriterion*(iter->getMax() - iter->getMin()) << "]";

    // if the value with its errors reaches the lower edge or the higher edge (+epsilon = convergencecriterion*total range) of the range, the fit has not converged
    if( ( (iter->getVal() - iter->getError()) < iter->getMin() + convergencecriterion*(iter->getMax() - iter->getMin()) ) ||
        ((iter->getVal() + iter->getError()) > iter->getMax() - convergencecriterion*(iter->getMax() - iter->getMin()) ) ){
      cout << endl << "\033[0;31mConvergence check failed at: "<< iter->GetName() << " = " << iter->getVal() << " +- " << iter->getError() << "\033[0m" << endl;
      return false;
    }
    else cout << " ... OK " << endl;
  }
  cout << endl << "\033[0;32mConvergence check passed! " << "\033[0m" << endl << endl;
  return true;
}

void setnewparranges(vector<RooRealVar*> vars){
  //double errfornewrange = 0.02;
  //double errangecheck = 2; // size of error compared to the initial range
  double newrange = 6; // sigma
  double convergencecriterion = 0.001;
  cout << "Assigning new ranges for some parameter..." << endl;
  for (vector<RooRealVar*>::const_iterator siter = vars.begin() ;  vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    //This checks if a parameter reached it's lower limit
    if( (iter->getVal() - iter->getError()) < iter->getMin() + convergencecriterion*(iter->getMax() - iter->getMin())){
      //If this guy has reached it's limit, expand its limits towards lower values
      iter->setRange(iter->getMin()-newrange*iter->getError(),iter->getMax());//be generous ;)
      cout << "New range for parameter " << iter->GetName() << " = [" << iter->getMin() <<"," << iter->getMax() << "]" << endl;
      //continue;
    }
    if ((iter->getVal() + iter->getError()) > iter->getMax() - convergencecriterion*(iter->getMax() - iter->getMin()) ){
      //If this guy has reached it's limit, expand its limits towards higher values
      iter->setRange(iter->getMin(),iter->getMax()+newrange*iter->getError());//be generous ;)
      cout << "New range for parameter " << iter->GetName() << " = [" << iter->getMin() <<"," << iter->getMax() << "]" << endl;
      //continue;
    }
    //This is to constrain the other parameters
    // if the error of the parameter is > 0 and reasonably small compared to the initial parameter range, set a new range
    /*if (iter->getError() > 0 && iter->getError()/iter->getVal() > errfornewrange && iter->getError() < errangecheck*(iter->getMax() - iter->getMin())){
            iter->setRange(iter->getVal()-newrange*iter->getError(),iter->getVal()+newrange*iter->getError());
            cout << "New range for parameter " << iter->GetName() << " = [" << iter->getMin() <<"," << iter->getMax() << "]" << endl;
        }*/
  }
  return;
}

#endif
