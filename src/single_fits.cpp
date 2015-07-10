/*!  

  \date July, 1st 2015
  \author Marian Stahl
  \brief If individual fits didn't converge...

*/

#include <iostream>     // std::cout
#include <algorithm>    // std::copy
#include <vector>       // std::vector

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TStopwatch.h"
#include "Riostream.h"
#include "TH1.h"

#include "RooRealVar.h"
#include "RooDataHist.h"

#include "../include/Fit_Jpsi.h"
#include "../include/MyStyle.h"
#include "../include/configuration.h"
#include "../include/get_binning.h"

using namespace RooFit;//std namespace set with ROOT

TString temp;//global string buffer
void get_asymmetries(configuration* myconfig);
int main(int argc, char **argv){
  int dataset = atoi(argv[1]);
  const TString systematic(argv[2]);
  double xbinlo = atof(argv[3]);
  double xbinhi = atof(argv[4]);
  configuration* myconfig = new configuration(dataset,systematic);
  myconfig->set_dim1_bin(xbinlo,xbinhi);
  if(argc == 6){
    double ybinlo = atof(argv[5]);
    double ybinhi = atof(argv[6]);
    myconfig->set_dim2_bin(ybinlo,ybinhi);
  }
  get_asymmetries(myconfig);
  delete myconfig;
  return 0;
}

void get_asymmetries(configuration* myconfig){

  TStopwatch *clock = new TStopwatch();  clock->Start(1);
  MyStyle();
  cout << "computing raw tracking asymmetries with the " << myconfig->get_method() << " method using " << myconfig->get_sample() << " data." << endl;

  //Initialize probe muon plus/minus pass/fail histograms for J/psi mass fits
  TH1D::AddDirectory(false);//don't append them to a file

  TH1D Jpsi_pass_p("pass_hist_p"," ; ; ",myconfig->get_fitbinning(),myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max());
  TH1D Jpsi_fail_p("fail_hist_p"," ; ; ",myconfig->get_fitbinning(),myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max());
  TH1D Jpsi_pass_m("pass_hist_m"," ; ; ",myconfig->get_fitbinning(),myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max());
  TH1D Jpsi_fail_m("fail_hist_m"," ; ; ",myconfig->get_fitbinning(),myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max());

  //get the tag and probe sample
  temp = myconfig->get_taptupledir() + "weighted_trackeff_" + myconfig->get_method() + ".root";
  TFile *fJpsi = new TFile(temp,"read");
  TTree *tap_tree = (TTree*)gDirectory->Get(myconfig->get_sample());
  Double_t var_dim1,var_dim2 = 0,Jpsi_mass,sWeight_matched,sWeight_total;
  Bool_t matched;
  Int_t charge;

  tap_tree->SetBranchStatus("*",0); //disable all branches
  //now switch on the ones we need
  tap_tree->SetBranchStatus(myconfig->get_dim1_alias(),1);
  if(!myconfig->is_1D())tap_tree->SetBranchStatus(myconfig->get_dim2_alias(),1);
  tap_tree->SetBranchStatus("Jpsi_mass",1);
  tap_tree->SetBranchStatus("matched",1);
  tap_tree->SetBranchStatus("charge",1);
  tap_tree->SetBranchStatus("sWeight_matched",1);
  tap_tree->SetBranchStatus("sWeight_total",1);
  tap_tree->SetBranchAddress(myconfig->get_dim1_alias(),&var_dim1);
  if(!myconfig->is_1D())tap_tree->SetBranchAddress(myconfig->get_dim2_alias(),&var_dim2);
  tap_tree->SetBranchAddress("Jpsi_mass",&Jpsi_mass);
  tap_tree->SetBranchAddress("matched",&matched);
  tap_tree->SetBranchAddress("charge",&charge);
  tap_tree->SetBranchAddress("sWeight_matched",&sWeight_matched);
  tap_tree->SetBranchAddress("sWeight_total",&sWeight_total);

  UInt_t ntap = tap_tree->GetEntries();
  cout << "Entries in tag-and-probe tree: " << ntap << endl;

  double sum_tap_sweights = 0;
  for (UInt_t i = 0; i < ntap;i++) {
    tap_tree->GetEntry(i);
    sum_tap_sweights += sWeight_total;
    //Fill the mass histograms in bins of kinematics
    if(var_dim1 >= myconfig->get_dim1_binlo() && var_dim1 < myconfig->get_dim1_binhi()){
      if( myconfig->is_1D() || (var_dim2 >= myconfig->get_dim2_binlo() && var_dim2 < myconfig->get_dim2_binhi()) ){
        if(charge > 0) matched ? Jpsi_pass_p.Fill(Jpsi_mass) : Jpsi_fail_p.Fill(Jpsi_mass);
        else matched ? Jpsi_pass_m.Fill(Jpsi_mass) : Jpsi_fail_m.Fill(Jpsi_mass);        
      }
    }
  }
  cout << "Tag-and-probe sum of sWeights " << sum_tap_sweights << endl;

  tap_tree->SetDirectory(0);
  fJpsi->Close();
  delete fJpsi;

  RooRealVar Jpsi_M("J_psi_1S_M","M_{inv}(#mu^{+}#mu^{-})",myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max(),"GeV") ;
  RooDataHist *Jpsi_plusprobe_pass = new RooDataHist("Jpsi_plusprobe_pass","Jpsi_plusprobe_pass",Jpsi_M,Import(Jpsi_pass_p)) ;
  RooDataHist *Jpsi_plusprobe_fail = new RooDataHist("Jpsi_plusprobe_fail","Jpsi_plusprobe_fail",Jpsi_M,Import(Jpsi_fail_p)) ;
  RooDataHist *Jpsi_minusprobe_pass = new RooDataHist("Jpsi_minusprobe_pass","Jpsi_minusprobe_pass",Jpsi_M,Import(Jpsi_pass_m)) ;
  RooDataHist *Jpsi_minusprobe_fail = new RooDataHist("Jpsi_minusprobe_fail","Jpsi_minusprobe_fail",Jpsi_M,Import(Jpsi_fail_m)) ;

  double *single_raw_asyms = Fit_Jpsi(Jpsi_plusprobe_pass,Jpsi_plusprobe_fail,Jpsi_minusprobe_pass,Jpsi_minusprobe_fail,myconfig);

  //save asymmetries in an ASCII file
  if(!myconfig->is_debug()) {
    ofstream *ofile = new ofstream;
    if(myconfig->is_1D())temp.Form("ATrackRaw/%s_%s_%s_%g_%g",myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),
                                   myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi());
    else temp.Form("ATrackRaw/%s_%s_%s_%g_%g_%g_%g",myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),
                   myconfig->get_dim1_binlo(),myconfig->get_dim1_binhi(),myconfig->get_dim2_binlo(),myconfig->get_dim2_binhi());
    ofile->open(temp);
    *ofile << single_raw_asyms[0] << "\t" << single_raw_asyms[1] << endl;
    ofile->close();
    delete ofile;
  }
  delete Jpsi_plusprobe_pass;
  delete Jpsi_plusprobe_fail;
  delete Jpsi_minusprobe_pass;
  delete Jpsi_minusprobe_fail;

  cout << "\033[0;32mFinished! Elapsed time: \033[0m" << endl;
  clock->Stop();clock->Print();delete clock;
  return;
}
