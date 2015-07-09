/*!  

  \date July, 1st 2015
  \author Marian Stahl
  \brief Calculate raw tracking asymmetries in 1 or 2D split by year, polarity and tag-and-probe method

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
  configuration* myconfig = new configuration(dataset,systematic);
  get_asymmetries(myconfig);
  delete myconfig;
  return 0;
}

void get_asymmetries(configuration* myconfig){

  TStopwatch *clock = new TStopwatch();  clock->Start(1);
  MyStyle();
  cout << "computing raw tracking asymmetries with the " << myconfig->get_method() << " method using " << myconfig->get_sample() << " data." << endl;

  //get Binning from ASCII file and store it in the binning vectors of the base-class
  // defined in ../include/get_binning.h
  get_binning(myconfig);

  unsigned int nbins_dim1 = myconfig->get_dim1_binedges().size()-1;
  unsigned int nbins_dim2 = myconfig->get_dim2_binedges().size()-1;
  //copy binning-vectors to arrays to deal with root histograms
  double* dim1_binedges = new double[nbins_dim1+1];//copy(dim1_bev.begin(),dim1_bev.end(),dim1_binedges);<-- why does this not work?
  double* dim2_binedges = new double[nbins_dim2+1];  
  dvec dim1_bev = myconfig->get_dim1_binedges();
  dvec dim2_bev = myconfig->get_dim2_binedges();
  copy(dim1_bev.begin(),dim1_bev.end(),dim1_binedges);
  copy(dim2_bev.begin(),dim2_bev.end(),dim2_binedges);

  //Initialize probe muon plus/minus pass/fail histograms for J/psi mass fits
  TH1D::AddDirectory(false);//don't append them to a file
  vector< vector<TH1D> > Jpsi_pass_p(nbins_dim1,vector<TH1D>(nbins_dim2));//a nbins_dim1 x nbins_dim2 TH1D matrix
  vector< vector<TH1D> > Jpsi_fail_p(nbins_dim1,vector<TH1D>(nbins_dim2));
  vector< vector<TH1D> > Jpsi_pass_m(nbins_dim1,vector<TH1D>(nbins_dim2));
  vector< vector<TH1D> > Jpsi_fail_m(nbins_dim1,vector<TH1D>(nbins_dim2));
  for(unsigned int bini = 0; bini < nbins_dim1; bini++){
    for(unsigned int binj = 0; binj < nbins_dim2; binj++){
      temp.Form("pass_hist_p_%d_%d",bini,binj);
      Jpsi_pass_p[bini][binj] = TH1D(temp," ; ; ",myconfig->get_fitbinning(),myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max());
      temp.Form("fail_hist_p_%d_%d",bini,binj);      
      Jpsi_fail_p[bini][binj] = TH1D(temp," ; ; ",myconfig->get_fitbinning(),myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max());
      temp.Form("pass_hist_m_%d_%d",bini,binj);
      Jpsi_pass_m[bini][binj] = TH1D(temp," ; ; ",myconfig->get_fitbinning(),myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max());
      temp.Form("fail_hist_m_%d_%d",bini,binj);
      Jpsi_fail_m[bini][binj] = TH1D(temp," ; ; ",myconfig->get_fitbinning(),myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max());
    }    
  }

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
    for(unsigned int bini = 0; bini < nbins_dim1; bini++){
      if(myconfig->is_debug() && i < 20)cout << dim1_bev.at(bini) << " " << dim1_bev.at(bini+1) << endl;
      if(var_dim1 >= dim1_bev.at(bini) && var_dim1 < dim1_bev.at(bini+1)){

        for(unsigned int binj = 0; binj < nbins_dim2; binj++){
          if( myconfig->is_1D() || (var_dim2 >= dim2_bev.at(binj) && var_dim2 < dim2_bev.at(binj+1)) ){

            if(charge > 0) matched ? Jpsi_pass_p[bini][binj].Fill(Jpsi_mass) : Jpsi_fail_p[bini][binj].Fill(Jpsi_mass);
            else matched ? Jpsi_pass_m[bini][binj].Fill(Jpsi_mass) : Jpsi_fail_m[bini][binj].Fill(Jpsi_mass);
            if(myconfig->is_debug() && i < 20)cout << "bin found" << endl;
            break;

          }
        }
      }
    }
  }
  cout << "Tag-and-probe sum of sWeights " << sum_tap_sweights << endl;

  tap_tree->SetDirectory(0);
  fJpsi->Close();
  delete fJpsi;

  //fit every of the nbins_dim1 x nbins_dim2 probe plus/minus pass/fail J/psi histograms
  for(unsigned int bini = 0; bini < nbins_dim1; bini++){
    for(unsigned int binj = 0; binj < nbins_dim2; binj++){

    RooRealVar Jpsi_M("J_psi_1S_M","M_{inv}(#mu^{+}#mu^{-})",myconfig->get_Jpsi_M_min(),myconfig->get_Jpsi_M_max(),"GeV") ;
    RooDataHist *Jpsi_plusprobe_pass = new RooDataHist("Jpsi_plusprobe_pass","Jpsi_plusprobe_pass",Jpsi_M,Import(Jpsi_pass_p[bini][binj])) ;
    RooDataHist *Jpsi_plusprobe_fail = new RooDataHist("Jpsi_plusprobe_fail","Jpsi_plusprobe_fail",Jpsi_M,Import(Jpsi_fail_p[bini][binj])) ;
    RooDataHist *Jpsi_minusprobe_pass = new RooDataHist("Jpsi_minusprobe_pass","Jpsi_minusprobe_pass",Jpsi_M,Import(Jpsi_pass_m[bini][binj])) ;
    RooDataHist *Jpsi_minusprobe_fail = new RooDataHist("Jpsi_minusprobe_fail","Jpsi_minusprobe_fail",Jpsi_M,Import(Jpsi_fail_m[bini][binj])) ;

    myconfig->set_dim1_bin(bini);
    if(!myconfig->is_1D())myconfig->set_dim2_bin(binj);
    double *single_raw_asyms = Fit_Jpsi(Jpsi_plusprobe_pass,Jpsi_plusprobe_fail,Jpsi_minusprobe_pass,Jpsi_minusprobe_fail,myconfig);

    //save asymmetries in an ASCII file
    if(!myconfig->is_debug()) {
      ofstream *ofile = new ofstream;
      if(myconfig->is_1D())temp.Form("ATrackRaw/%s_%s_%s_%g_%g",myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),
                                     dim1_bev.at(bini),dim1_bev.at(bini+1));
      else temp.Form("ATrackRaw/%s_%s_%s_%g_%g_%g_%g",myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),
                     dim1_bev.at(bini),dim1_bev.at(bini+1),dim2_bev.at(binj),dim2_bev.at(binj+1));
      ofile->open(temp);
      *ofile << single_raw_asyms[0] << "\t" << single_raw_asyms[1] << endl;
      ofile->close();
      delete ofile;
    }
    delete Jpsi_plusprobe_pass;
    delete Jpsi_plusprobe_fail;
    delete Jpsi_minusprobe_pass;
    delete Jpsi_minusprobe_fail;
    }
  }

  cout << "\033[0;32mFinished! Elapsed time: \033[0m" << endl;
  clock->Stop();clock->Print();delete clock;
  return;
}
