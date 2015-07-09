/*!

  \date July, 1st 2015
  \author Marian Stahl
  \brief Example script to calculate mupi tracking asymmetries for asls in 1D

*/

#include <stdio.h>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TObject.h"
#include "Riostream.h"
#include "TH1.h"
#include "TStopwatch.h"

#include "../include/AsymCalc.h"
#include "../include/configuration.h"
#include "../include/MyStyle.h"
#include "../include/get_binning.h"
#include "../include/make_tables.h"
#include "../include/make_plots.h"

TString temp;
int nbinsCH = 5000;//number of bins for histo to compute mu-pi correlation in

//dvec, dmtx and dten3 typedef'd in configuration.h, short for vectors, matrices and dim. 3 tensors with doubles
//three_vectors is also defined in configuration.h
three_vectors combine_asymmetries(configuration *myconfig);
dmtx get_correlations(configuration *myconfig);
double fill_histograms_for_mupi_correlation(TH1D *mu_hist, TH1D *pi_hist, configuration *myconfig);

int main(int argc, char **argv){

  const TString jobname(argv[1]);
  TStopwatch *clock = new TStopwatch();
  clock->Start(1);

  //Construct base-class. Defined in "../include/configuration.h"
  configuration *myconfig = new configuration(jobname);
  //myconfig->set_verbosity(0);

  vector<three_vectors> results;
  dvec center_VT;dvec unc_VT;dvec center_Long; dvec unc_Long;
  for(int mode = 0; mode < myconfig->get_nchannels(); mode++){
    try {if(!myconfig->set_channel(mode))throw 3;}catch(int e){cout << "Exception number " << e << "\n Something wrong with setter of sample, channel, method! Terminating..." << endl;terminate();}
    results.push_back(combine_asymmetries(myconfig));
    center_VT.push_back(results.at(mode).combined_results.at(0));
    unc_VT.push_back(results.at(mode).combined_results.at(1));
    center_Long.push_back(results.at(mode).combined_results.at(2));
    unc_Long.push_back(results.at(mode).combined_results.at(3));
  }

  dmtx correlations = get_correlations(myconfig);

  if(myconfig->get_verbosity() > 0){
    cout << endl << endl << "Correlation matrix Dalitz-regions:" << endl << endl;
    for(unsigned int i = 0; i < correlations.at(0).size(); i++){
      for(unsigned int j = 0; j < correlations.at(0).size(); j++){
        printf("%.4f \t", correlations[i][j]);
      }
      cout << " " << endl;
    }
  }

  double comb_center = corr_comb(center_VT,unc_VT,correlations);
  double comb_unc = corr_comb_unc(unc_VT,correlations);
  double Long_center = corr_comb(center_Long,unc_Long,correlations);
  double Long_unc = corr_comb_unc(unc_Long,correlations);
  cout << endl << string(128, '*') << endl << endl << "Integrated " << myconfig->get_firstParticle() << myconfig->get_secondParticle() << " tracking asymmetries for "
       << myconfig->get_bw() << endl << endl;
  printf("Long\t\t\t( %+.3f +- %.3f ) %%\n", Long_center, Long_unc);
  printf("VELO + T Station\t\033[0;32m( %+.3f +- %.3f ) %%\033[0m\n\n", comb_center, comb_unc);

  ofstream resultfile;
  temp = myconfig->get_dumpdir() + "Results/FullDP_" + myconfig->get_bw() + "_raw_VT";
  resultfile.open(temp);
  resultfile << comb_center << "\t" << comb_unc;
  resultfile.close();

  //defined in "../include/make_tables.h"
  make_latex_table(results,myconfig);
  make_python_dictionary(results,myconfig);

  clock->Stop();
  clock->Print();

  delete myconfig;
  return 0;

}

three_vectors combine_asymmetries(configuration *myconfig){

  //get Binning from ASCII file and store it in the binning vectors of the base-class
  //defined in "../include/get_binning.h"
  get_binning(myconfig);

  //copy binning-vectors to arrays to deal with root histograms
  unsigned int nbins_dim1 = myconfig->get_dim1_binedges().size()-1;
  //copy binning-vectors to arrays to deal with root histograms
  double* dim1_binedges = new double[nbins_dim1+1];//copy(dim1_bev.begin(),dim1_bev.end(),dim1_binedges);<-- why does this not work?
  dvec dim1_bev = myconfig->get_dim1_binedges();
  copy(dim1_bev.begin(),dim1_bev.end(),dim1_binedges);
  //Tensors with all the raw asymmetries and their uncertainties (split by sample, method and dim1_bin)
  dten3 A_raw(myconfig->get_nmethods()+1,dmtx(myconfig->get_nsamples(),dvec(nbins_dim1)));
  dten3 delta_A_raw(myconfig->get_nmethods()+1,dmtx(myconfig->get_nsamples(),dvec(nbins_dim1)));

  for(int i = 0; i <= myconfig->get_nmethods(); i++){//using VELO+TStation as the 4th method
    if(i == myconfig->get_nmethods()){//combine VELO and TStation
      for(int j = 0; j < myconfig->get_nsamples(); j++){
        for(unsigned int k = 0; k < nbins_dim1; k++){
          A_raw[i][j][k] = A_raw[0][j][k]+A_raw[1][j][k];
          delta_A_raw[i][j][k] = add_asymmetry_unc(delta_A_raw[0][j][k],delta_A_raw[1][j][k]);
        }
      }
    }
    else{//Fill raw asymmetries using usual tap methods
      for(int j = 0; j < myconfig->get_nsamples(); j++){
        try {if(!myconfig->set_method(i))throw 3;}catch(int e){cout << "Exception number " << e << "\n Something wrong with setter of sample, channel, method! Terminating..." << endl;terminate();}
        try {if(!myconfig->set_sample(j))throw 3;}catch(int e){cout << "Exception number " << e << "\n Something wrong with setter of sample, channel, method! Terminating..." << endl;terminate();}
        for(unsigned int k = 0; k < nbins_dim1; k++){
          //Read the raw asymmetries from disc and fill the A_raw tensor
          ifstream infile;
          temp.Form("ATrackRaw/%s_%s_%s_%g_%g",myconfig->get_method().Data(),myconfig->get_sample().Data(),myconfig->get_bw().Data(),
                    myconfig->get_dim1_binedges().at(k),myconfig->get_dim1_binedges().at(k+1));
          infile.open(temp);
          if(!infile)cout << "\033[0;31mCould not read " << temp << " ! Check if job ran!" << endl;
          double asym,dasym;
          infile >> asym >> dasym;
          asym *= 100;
          dasym *= 100;
          A_raw[i][j][k] = asym;
          delta_A_raw[i][j][k] = dasym;
          infile.close();
        }
      }
    }
  }
  //only need to print and plot the raw asymmetry once...
  //defined in ../include/make_plots.h
  if( myconfig->get_channel().CompareTo(myconfig->get_channels().at(0)) == 0 && myconfig->is_1D() )
    make_1D_plots(A_raw,delta_A_raw,myconfig);

  //now we need loops over samples and bins, get signal muons and pions and their correlation
  vector<TH1D*> sig_mu_hists;
  vector<TH1D*> sig_pi_hists;
  dvec mu_pi_correlations;
  double sum_total_sweights = 0;

  for(int ypol = 0; ypol < myconfig->get_nsamples(); ypol++){
    try {if(!myconfig->set_sample(ypol))throw 3;}catch(int e){cout << "Exception number " << e << "\n Something wrong with setter of sample, channel, method! Terminating..." << endl;terminate();}

    //these histograms will be used later
    temp.Form("signal_%s_%d",myconfig->get_channel().Data(),ypol);
    TH1D *signal_muon = new TH1D(temp,";;",nbins_dim1,dim1_binedges);
    TH1D *signal_pion = new TH1D(*signal_muon);
    //these are just to get the correlation
    TH1D *signal_muon_fc = new TH1D(temp,";;",nbinsCH,myconfig->get_dim1_binedges().at(0),myconfig->get_dim1_binedges().at(nbins_dim1));
    TH1D *signal_pion_fc = new TH1D(*signal_muon_fc);

    double sum_sweights = fill_histograms_for_mupi_correlation(signal_muon,signal_pion,myconfig);
    sum_sweights = fill_histograms_for_mupi_correlation(signal_muon_fc,signal_pion_fc,myconfig);//yes, sum_sweights is the same as in line above

    double rho_denominator_mu = 0, rho_denominator_pi = 0, rho_numerator = 0;
    signal_muon->Scale(1/sum_sweights);
    signal_pion->Scale(1/sum_sweights);
    if(myconfig->is_debug())cout << "signal_muon.Integral() " << signal_muon->Integral() << endl;
    if(myconfig->is_debug())cout << "signal_pion.Integral() " << signal_pion->Integral() << endl;
    sig_mu_hists.push_back(move(signal_muon));
    sig_pi_hists.push_back(move(signal_pion));
    sum_total_sweights += sum_sweights;

    for(int bin = 1; bin < signal_muon_fc->GetNbinsX(); bin++){
      rho_numerator += signal_muon_fc->GetBinContent(bin)*signal_pion_fc->GetBinContent(bin);
      rho_denominator_mu += pow(signal_muon_fc->GetBinContent(bin),2);
      rho_denominator_pi += pow(signal_pion_fc->GetBinContent(bin),2);
    }
    mu_pi_correlations.push_back(rho_numerator/sqrt(rho_denominator_mu*rho_denominator_pi));

    delete signal_muon_fc;
    delete signal_pion_fc;
  }
  delete dim1_binedges;

  // BEGIN calculate A_mupi
  cout << endl << string(128, '*') << endl << endl << myconfig->get_firstParticle() << myconfig->get_secondParticle() << " tracking asymmetries for "
       << myconfig->get_bw() << " in the " << myconfig->get_channel() << " channel" << endl << endl;
  //initialize Amupi results integrated over phase space for every method and sample
  dmtx A_mupi(myconfig->get_nmethods() + 1,dvec(myconfig->get_nsamples()));
  dmtx delta_A_mupi(myconfig->get_nmethods() + 1,dvec(myconfig->get_nsamples()));

  for(int i = 0; i < myconfig->get_nmethods() + 1; i++){//VELO + T as 4th method
    for(int j = 0; j < myconfig->get_nsamples(); j++){

      dvec A_mu_b;dvec delta_A_mu_b;
      dvec A_pi_b;dvec delta_A_pi_b;
      //dvec A_mupi_b;dvec delta_A_mupi_b;
      if(myconfig->is_debug())cout << static_cast<TH1D*>(sig_mu_hists.at(j))->GetName() << endl;

      for(unsigned int k = 0; k < nbins_dim1; k++){

        double mu_binc = static_cast<TH1D*>(sig_mu_hists.at(j))->GetBinContent(k+1);
        double pi_binc = static_cast<TH1D*>(sig_pi_hists.at(j))->GetBinContent(k+1);
        if(i == 0 && myconfig->is_debug())cout << "mu_binc\t" << mu_binc << "\t\tpi_binc\t" << pi_binc << endl;

        A_mu_b.push_back(A_raw[i][j][k]*mu_binc);delta_A_mu_b.push_back(delta_A_raw[i][j][k]*mu_binc);
        A_pi_b.push_back(A_raw[i][j][k]*pi_binc);delta_A_pi_b.push_back(delta_A_raw[i][j][k]*pi_binc);
        //one can easily write out the phase space dependent mu-pi asymmetries here and use make_1D_plots(A_raw,delta_A_raw,myconfig) again later.
        //A_mupi_b.push_back(A_raw[i][j][k]*(mu_binc-pi_binc));delta_A_mupi_b.push_back(delta_A_raw[i][j][k]*(mu_binc-pi_binc));

      }
      A_mupi[i][j] = summation(A_mu_b) - summation(A_pi_b);
      dvec tmp_unc_for_comb = {summation_unc(delta_A_mu_b),summation_unc(delta_A_pi_b)};
      dmtx mupi_correlation_matrix = {{1,mu_pi_correlations.at(j)},{mu_pi_correlations.at(j),1}};
      delta_A_mupi[i][j] = corr_comb_unc(tmp_unc_for_comb,mupi_correlation_matrix);
    }
  }
  // END calculate A_mupi

  // BEGIN combine years and polarity
  dvec combined_A_mupi;dvec combined_delta_A_mupi;

  for(int i = 0; i < myconfig->get_nmethods() + 1; i++){
    //linear combination of polarities
    dvec linear_combination;dvec unc_linear_combination;
    for(int ji = 0; ji < myconfig->get_nchannels()/2; ji++){
      linear_combination.push_back((A_mupi[i][2*ji]+A_mupi[i][2*ji+1])/2);
      unc_linear_combination.push_back(add_asymmetry_unc(delta_A_mupi[i][2*ji],delta_A_mupi[i][2*ji])/2);
    }
    //A small trick here: the mu_A matrix at i is the vector containing the integrated A_mu for every sample as entries (I think it has been used above somewhere as well)
    combined_A_mupi.push_back(combination(linear_combination,unc_linear_combination));
    combined_delta_A_mupi.push_back(combined_unc(unc_linear_combination));
  }

  printf("\nVELO\t\t\t( %+.3f +- %.3f ) %%\n", combined_A_mupi.at(0), combined_delta_A_mupi.at(0));
  printf("T-station\t\t( %+.3f +- %.3f ) %%\n", combined_A_mupi.at(1), combined_delta_A_mupi.at(1));
  printf("Long\t\t\t( %+.3f +- %.3f ) %%\n", combined_A_mupi.at(2), combined_delta_A_mupi.at(2));
  printf("VELO + T Station\t\033[0;32m( %+.3f +- %.3f ) %%\033[0m\n", combined_A_mupi.at(3), combined_delta_A_mupi.at(3));

  ofstream resultfile;
  temp = myconfig->get_dumpdir() + "Results";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);
  temp += "/"+myconfig->get_channel() + "_" + myconfig->get_bw() + "_raw_VT";
  resultfile.open(temp);
  resultfile << combined_A_mupi.at(3) << "\t" << combined_delta_A_mupi.at(3);
  resultfile.close();

  three_vectors results;
  results.combined_results = {combined_A_mupi.at(3),combined_delta_A_mupi.at(3),combined_A_mupi.at(2),combined_delta_A_mupi.at(2)};
  results.cv_VT = A_mupi.at(3);
  results.unc_VT = delta_A_mupi.at(3);

  return results;
}

dmtx get_correlations(configuration *myconfig){

  dmtx correlations(myconfig->get_nchannels(),dvec(myconfig->get_nchannels()));
  dvec sum_square_weights;

  TH1D::SetDefaultSumw2(true);

  dvec dim1_bev = myconfig->get_dim1_binedges();
  double* dim1_binedges = new double[dim1_bev.size()];
  copy(dim1_bev.begin(),dim1_bev.end(),dim1_binedges);

  vector<TH1D*> mupi_diff_hists;
  for(int mode = 0; mode < myconfig->get_nchannels(); mode++){
    try {if(!myconfig->set_channel(mode))throw 3;}catch(int e){cout << "Exception number " << e << "\n Something wrong with setter of sample, channel, method! Terminating..." << endl;terminate();}

    TH1D *mu_hist = new TH1D("mu_hist",";;",dim1_bev.size()-1,dim1_binedges);
    TH1D *pi_hist = (TH1D*)mu_hist->Clone("pi_hist");
    TH1D *mupi_diff_hist = (TH1D*)mu_hist->Clone("mupi_hist");

    for(int ypol = 0; ypol < myconfig->get_nsamples(); ypol++){
      try {if(!myconfig->set_sample(ypol))throw 3;}catch(int e){cout << "Exception number " << e << "\n Something wrong with setter of sample, channel, method! Terminating..." << endl;terminate();}

      TH1D *mu_temp_hist = (TH1D*)mu_hist->Clone("mu_temp_hist");
      TH1D *pi_temp_hist = (TH1D*)mu_hist->Clone("pi_temp_hist");
      fill_histograms_for_mupi_correlation(mu_temp_hist,pi_temp_hist,myconfig);
      mu_hist->Add(mu_temp_hist,1);pi_hist->Add(pi_temp_hist,1);
      delete mu_temp_hist;
      delete pi_temp_hist;

    }
    if(!mupi_diff_hist->Add(mu_hist,pi_hist,1,-1)){cout << "\033[0;31mAdd operation of weighted histograms failed. TERMINATE\033[0m" << endl; terminate();}
    mupi_diff_hists.push_back(mupi_diff_hist);
    delete mu_hist;
    delete pi_hist;
  }
  delete dim1_binedges;
  //calculate correlations rho_ij = N_i*N_j/(sqrt(N_i^2*N_j^2))
  for(int i = 0; i < myconfig->get_nchannels(); i++){
    correlations[i][i] = 1.0;
    double sumsq = 0;
    for(unsigned int binx = 1;binx < dim1_bev.size() - 1; binx++ ){
      sumsq += pow(mupi_diff_hists.at(i)->GetBinContent(binx),2);
      if(i != myconfig->get_nchannels() - 1){
        for(int j = i+1; j < myconfig->get_nchannels(); j++){
          correlations[i][j] += mupi_diff_hists.at(i)->GetBinContent(binx)*mupi_diff_hists.at(j)->GetBinContent(binx);
        }
      }
    }
    sum_square_weights.push_back(sumsq);
  }
  //normalize the beast
  for(int i = 0; i < myconfig->get_nchannels(); i++){
    if(i != myconfig->get_nchannels() - 1){
      for(int j = i+1; j < myconfig->get_nchannels(); j++){
        correlations[i][j] /= sqrt(sum_square_weights.at(i)*sum_square_weights.at(j));
        correlations[j][i] = correlations[i][j];
      }
    }
  }
  return correlations;
}

double fill_histograms_for_mupi_correlation(TH1D *mu_hist, TH1D *pi_hist, configuration *myconfig){

  temp = myconfig->get_sigtupledir() + myconfig->get_sigtuplename();
  if(myconfig->is_debug())cout << temp << endl;
  gErrorIgnoreLevel = kError;
  TH1D::AddDirectory(false);
  TFile *fasls = new TFile(temp);
  TTree *sig_tree = (TTree*)gDirectory->Get(myconfig->get_sigtreename());

  bool sweight_in_friendtree = myconfig->get_sweighttuplename().Length() > 5;//at least 5 characters for the .root extension
  if(sweight_in_friendtree){
    temp = myconfig->get_sweighttupledir() + myconfig->get_sweighttuplename();
    sig_tree->AddFriend(myconfig->get_sweighttreename(),temp);
  }
  gErrorIgnoreLevel = kPrint;

  //using leaves here, because TLeaf::GetValue() always returns double irrespective of the declared type in the tree
  temp = myconfig->get_firstParticle() + "_" + myconfig->get_dim1();
  TLeaf *part1_dim1_leaf = sig_tree->GetLeaf(temp);
  temp = myconfig->get_secondParticle() + "_" + myconfig->get_dim1();
  TLeaf *part2_dim1_leaf = sig_tree->GetLeaf(temp);
  if(sweight_in_friendtree)
    temp = myconfig->get_sweighttreename() + "." + myconfig->get_sweightvarname();
  else
    temp = myconfig->get_sweightvarname();

  TLeaf *sWeight_leaf = sig_tree->GetLeaf(temp);

  UInt_t tree_entries = sig_tree->GetEntries();
  if(myconfig->is_debug())cout << tree_entries << endl;
  double sum_sWeights = 0.0;
  for (UInt_t i = 0; i < tree_entries;i++) {
    sig_tree->GetEntry(i);
    mu_hist->Fill(part1_dim1_leaf->GetValue(), sWeight_leaf->GetValue());
    pi_hist->Fill(part2_dim1_leaf->GetValue(), sWeight_leaf->GetValue());
    sum_sWeights += sWeight_leaf->GetValue();
  }

  fasls->Close();
  return sum_sWeights;
}
