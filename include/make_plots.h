#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "Riostream.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TStopwatch.h"
#include "TGraphErrors.h"

#include "configuration.h"
extern TString temp;

void make_a_plot(TH1D *Up_hist, TH1D *Down_hist, TH1D *Comb_hist, configuration *myconfig);
void make_a_plot(vector<TH1D*> hists, TH1D *Comb_hist, configuration *myconfig);
void make_a_plot(TH1D *Up_hist, TH1D *Down_hist, configuration *myconfig);
void make_a_plot(vector<TH1D*> hists, configuration *myconfig);
void make_a_plot(TH1D *Comb_hist, double cv, double ecv, configuration *myconfig);

void make_1D_plots(const dten3 &A_raw, dten3 &delta_A_raw, configuration *myconfig){

  TH1D *VT_Up_hist = nullptr, *VT_Down_hist = nullptr, *VT_comb_hist = nullptr;
  vector<TH1D*> VT_hists_by_year;

  //copy binning-vectors to arrays to deal with root histograms
  unsigned int nbins_dim1 = myconfig->get_dim1_binedges().size()-1;
  //copy binning-vectors to arrays to deal with root histograms
  double* dim1_binedges = new double[nbins_dim1+1];//copy(dim1_bev.begin(),dim1_bev.end(),dim1_binedges);<-- why does this not work?
  dvec dim1_bev = myconfig->get_dim1_binedges();
  copy(dim1_bev.begin(),dim1_bev.end(),dim1_binedges);

  if(myconfig->get_dim1().CompareTo("PT")==0){
    dim1_binedges[nbins_dim1]=10000;//upper limit for pT plots
  }
  else if(myconfig->get_dim1().CompareTo("ETA")==0){
    dim1_binedges[0]=2;//lower limit for eta plots
    dim1_binedges[nbins_dim1]=5;//upper limit for eta plots
  }
  else if(myconfig->get_dim1().CompareTo("P")==0){
    dim1_binedges[nbins_dim1]=100000;//upper limit for p plots
  }
  else cout << "come on... " << endl;

  temp = " ; " + myconfig->get_dim1_label() + " ;A_{track} (%)";
  TH1D::AddDirectory(false);
  VT_Up_hist = new TH1D("VT_Up_hist",temp,nbins_dim1,dim1_binedges);
  VT_Down_hist = new TH1D(*VT_Up_hist);
  VT_comb_hist = new TH1D(*VT_Up_hist);
  for(int ji = 0; ji < myconfig->get_nsamples()/2; ji++){
    VT_hists_by_year.emplace_back(new TH1D(*VT_Up_hist));
  }

  //reset binedges
  dim1_binedges[0]=myconfig->get_dim1_binedges().at(0);
  dim1_binedges[nbins_dim1]=myconfig->get_dim1_binedges().at(nbins_dim1);

  double A_raw_VELO = 0, delta_A_raw_VELO = 0, A_raw_TStation = 0, delta_A_raw_TStation = 0, A_raw_Long = 0, delta_A_raw_Long = 0, A_raw_VT = 0, delta_A_raw_VT = 0;
  dvec A_VELO_Comb;
  dvec delta_A_VELO_Comb;
  dvec A_TStation_Comb;
  dvec delta_A_TStation_Comb;
  dvec A_VT_Comb;
  dvec delta_A_VT_Comb;
  dvec A_Long_Comb;
  dvec delta_A_Long_Comb;

  for(unsigned int k = 0; k < nbins_dim1; k++){

    dvec tmpup;  dvec tmpup_unc;
    dvec tmpdown;dvec tmpdown_unc;
    for(int j = 0; j < myconfig->get_nsamples(); j++){
      if(j%2){tmpdown.push_back(A_raw[3][j][k]); tmpdown_unc.push_back(delta_A_raw[3][j][k]);} //VT Down data
      else{tmpup.push_back(A_raw[3][j][k]); tmpup_unc.push_back(delta_A_raw[3][j][k]);} //VT Up data
    }
    double A_VT_up = combination(tmpup,tmpup_unc);
    double delta_A_VT_up = combined_unc(tmpup_unc);
    double A_VT_down = combination(tmpdown,tmpdown_unc);
    double delta_A_VT_down = combined_unc(tmpdown_unc);

    double A_VT_comb = (A_VT_up+A_VT_down)/2;//linear combination of up and down
    double delta_A_VT_comb = sqrt(pow(delta_A_VT_up,2)+pow(delta_A_VT_down,2))/2;
    A_VT_Comb.push_back(A_VT_comb);
    delta_A_VT_Comb.push_back(delta_A_VT_comb);

    //Fill historgrams that are split by year
    for(int ji = 0; ji < myconfig->get_nsamples()/2; ji++){
      VT_hists_by_year.at(ji)->SetBinContent(k+1,(A_raw[3][2*ji][k]+A_raw[3][2*ji+1][k])/2);
      VT_hists_by_year.at(ji)->SetBinError(k+1,add_asymmetry_unc(delta_A_raw[3][2*ji][k],delta_A_raw[3][2*ji+1][k])/2);
    }

    VT_Up_hist->SetBinContent(k+1,A_VT_up);
    VT_Up_hist->SetBinError(k+1,delta_A_VT_up);
    VT_Down_hist->SetBinContent(k+1,A_VT_down);
    VT_Down_hist->SetBinError(k+1,delta_A_VT_down);
    VT_comb_hist->SetBinContent(k+1,A_VT_comb);
    VT_comb_hist->SetBinError(k+1,delta_A_VT_comb);

    //Combine samples of remaining methods
    for(int i = 0; i < myconfig->get_nmethods(); i++){
      for(int ji = 0; ji < myconfig->get_nsamples()/2; ji++){
        if(i == 0){
          A_VELO_Comb.push_back((A_raw[i][2*ji][k]+A_raw[i][2*ji+1][k])/2);
          delta_A_VELO_Comb.push_back(add_asymmetry_unc(delta_A_raw[i][2*ji][k],delta_A_raw[i][2*ji+1][k])/2);
        }
        else if(i == 1){
          A_TStation_Comb.push_back((A_raw[i][2*ji][k]+A_raw[i][2*ji+1][k])/2);
          delta_A_TStation_Comb.push_back(add_asymmetry_unc(delta_A_raw[i][2*ji][k],delta_A_raw[i][2*ji+1][k])/2);
        }
        else if(i == 2){
          A_Long_Comb.push_back((A_raw[i][2*ji][k]+A_raw[i][2*ji+1][k])/2);
          delta_A_Long_Comb.push_back(add_asymmetry_unc(delta_A_raw[i][2*ji][k],delta_A_raw[i][2*ji+1][k])/2);
        }
        else {cout << "something wrong with combining raw asymmetries..." << endl; terminate();}
      }
    }
  }

  A_raw_VELO = combination(A_VELO_Comb,delta_A_VELO_Comb);
  delta_A_raw_VELO = combined_unc(delta_A_VELO_Comb);
  A_raw_TStation = combination(A_TStation_Comb,delta_A_TStation_Comb);
  delta_A_raw_TStation = combined_unc(delta_A_TStation_Comb);
  A_raw_Long = combination(A_Long_Comb,delta_A_Long_Comb);
  delta_A_raw_Long = combined_unc(delta_A_Long_Comb);
  A_raw_VT = combination(A_VT_Comb,delta_A_VT_Comb);
  delta_A_raw_VT = combined_unc(delta_A_VT_Comb);

  cout << endl << string(128, '*') << endl << endl << "RAW TRACKING ASYMMETRIES FOR " << myconfig->get_bw() << endl << endl;
  cout << "VELO\t\t\t(" <<  A_raw_VELO << " $\\pm$ " << delta_A_raw_VELO << " ) \\% " << endl;
  cout << "T Station\t\t(" <<  A_raw_TStation << " $\\pm$ " << delta_A_raw_TStation << " ) \\% " << endl;
  cout << "Long\t\t\t(" <<  A_raw_Long << " $\\pm$ " << delta_A_raw_Long << " ) \\% " << endl;
  cout << "VELO + T Station\t(" <<  A_raw_VT << " $\\pm$ " << delta_A_raw_VT << " ) \\% " << endl << endl;

  MyStyle();

  //asymmetry_hist->GetYaxis()->SetTitleOffset(0.8);
  double yminmax = TMath::Max(TMath::Max(fabs(VT_Up_hist->GetBinContent(VT_Up_hist->GetMaximumBin())+VT_Up_hist->GetBinError(VT_Up_hist->GetMaximumBin())),
                                         fabs(VT_Up_hist->GetBinContent(VT_Up_hist->GetMinimumBin())-VT_Up_hist->GetBinError(VT_Up_hist->GetMinimumBin()))),
                              TMath::Max(fabs(VT_Down_hist->GetBinContent(VT_Down_hist->GetMaximumBin())+VT_Down_hist->GetBinError(VT_Down_hist->GetMaximumBin())),
                                         fabs(VT_Down_hist->GetBinContent(VT_Down_hist->GetMinimumBin())-VT_Down_hist->GetBinError(VT_Down_hist->GetMinimumBin()))));
  double ymin = -ceil(10*yminmax)/10;//assuming all histos have neg. entries too
  double ymax_temp = ceil(10*yminmax)/10;
  double ymax = ymax_temp+0.2*(ymax_temp-ymin);
  //Comb_hist->GetYaxis()->SetRangeUser(ymin,ymax_temp);
  VT_Up_hist->GetYaxis()->SetRangeUser(-ymax,ymax);
  VT_Down_hist->GetYaxis()->SetRangeUser(-ymax,ymax);
  VT_comb_hist->GetYaxis()->SetRangeUser(-ymax,ymax);
  for(int ji = 0; ji < myconfig->get_nsamples()/2; ji++){VT_hists_by_year.at(ji)->GetYaxis()->SetRangeUser(-ymax,ymax);}

  make_a_plot(VT_Up_hist,VT_Down_hist,VT_comb_hist,myconfig);
  make_a_plot(VT_Up_hist,VT_Down_hist,myconfig);
  make_a_plot(VT_comb_hist,A_raw_VT,delta_A_raw_VT,myconfig);
  make_a_plot(VT_hists_by_year,myconfig);
  make_a_plot(VT_hists_by_year,VT_comb_hist,myconfig);

  return;

}

void make_a_plot(TH1D *Up_hist, TH1D *Down_hist, TH1D *Comb_hist, configuration *myconfig){

  TCanvas* cVT = new TCanvas("cVT","Canvas",10,10,1280,960) ;
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetLeftMargin(0.15);

  Comb_hist->SetLineWidth(2);
  Comb_hist->SetLineColor(kGray+2);
  Comb_hist->SetMarkerSize(1.5);
  Comb_hist->SetMarkerColor(kGray+2);
  Comb_hist->GetXaxis()->SetLabelSize(0.05);
  Comb_hist->GetYaxis()->SetLabelSize(0.05);
  Comb_hist->GetXaxis()->SetTitleSize(0.05);
  Comb_hist->GetYaxis()->SetTitleSize(0.05);
  Comb_hist->Draw("e1p");

  Down_hist->SetLineWidth(2);
  Down_hist->SetLineColor(kBlue+2);
  Down_hist->SetMarkerSize(1.5);
  Down_hist->SetMarkerColor(kBlue+2);

  Up_hist->SetLineWidth(2);
  Up_hist->SetLineColor(kOrange+1);
  Up_hist->SetMarkerSize(1.5);
  Up_hist->SetMarkerColor(kOrange+1);

  TLine zeroline(Comb_hist->GetXaxis()->GetXmin(),0,Comb_hist->GetXaxis()->GetXmax(),0);
  zeroline.SetLineStyle(kDashed);
  zeroline.SetLineColor(1);
  zeroline.SetLineWidth(0.8);
  zeroline.Draw();

  Comb_hist->Draw("e1psame");
  Down_hist->Draw("e1psame");
  Up_hist->Draw("e1psame");

  TPaveText *blank_box = new TPaveText(1.003-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();

  TLegend *leg;
  double maxabs;
  fabs(Up_hist->GetBinContent(Up_hist->GetNbinsX()))>fabs(Down_hist->GetBinContent(Down_hist->GetNbinsX()))
      ? maxabs = Up_hist->GetBinContent(Up_hist->GetNbinsX())
      : maxabs = Down_hist->GetBinContent(Down_hist->GetNbinsX());
  if(maxabs<0)leg = new TLegend(0.6,0.77-gPad->GetTopMargin(),0.97-gPad->GetRightMargin(),0.97-gPad->GetTopMargin());
  else leg = new TLegend(0.6,0.05+gPad->GetBottomMargin(),0.97-gPad->GetRightMargin(),0.3+gPad->GetBottomMargin());
  leg->SetBorderSize(0);leg->SetFillColor(kWhite);leg->SetFillStyle(1001);leg->SetTextAlign(12);leg->SetTextSize(0.05);leg->SetTextFont(42);
  leg->SetHeader("VELO + T-station");
  leg->AddEntry(Up_hist,"Magnet Up","le1p");
  leg->AddEntry(Comb_hist,"Combined","le1p");
  leg->AddEntry(Down_hist,"Magnet Down","le1p");

  leg->Draw();

  temp = myconfig->get_dumpdir() + "Plots/UpDownCombined_asymmetries_" + myconfig->get_bw() + ".pdf";

  cVT->SaveAs(temp);
  delete cVT;
  return;
}

void make_a_plot(vector<TH1D*> hists, TH1D *Comb_hist, configuration *myconfig){

  if(hists.size() > 2){cout << "Sorry, but plotting is only implemented for 11 and 12 data" << endl; return;}

  TCanvas* cVT = new TCanvas("cVT","Canvas",10,10,1280,960) ;
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetLeftMargin(0.15);

  Comb_hist->SetLineWidth(2);
  Comb_hist->SetLineColor(kGray+2);
  Comb_hist->SetMarkerSize(1.5);
  Comb_hist->SetMarkerColor(kGray+2);
  Comb_hist->GetXaxis()->SetLabelSize(0.05);
  Comb_hist->GetYaxis()->SetLabelSize(0.05);
  Comb_hist->GetXaxis()->SetTitleSize(0.05);
  Comb_hist->GetYaxis()->SetTitleSize(0.05);
  Comb_hist->Draw("e1p");

  hists.at(0)->SetLineWidth(2);
  hists.at(0)->SetLineColor(kGreen+3);
  hists.at(0)->SetMarkerSize(1.5);
  hists.at(0)->SetMarkerColor(kGreen+3);

  hists.at(1)->SetLineWidth(2);
  hists.at(1)->SetLineColor(kRed+2);
  hists.at(1)->SetMarkerSize(1.5);
  hists.at(1)->SetMarkerColor(kRed+2);

  TLine zeroline(Comb_hist->GetXaxis()->GetXmin(),0,Comb_hist->GetXaxis()->GetXmax(),0);
  zeroline.SetLineStyle(kDashed);
  zeroline.SetLineColor(1);
  zeroline.SetLineWidth(0.8);
  zeroline.Draw();

  Comb_hist->Draw("e1psame");
  hists.at(0)->Draw("e1psame");
  hists.at(1)->Draw("e1psame");

  TPaveText *blank_box = new TPaveText(1.003-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();

  TLegend *leg;
  double maxabs;
  fabs(hists.at(1)->GetBinContent(hists.at(1)->GetNbinsX()))>fabs(hists.at(0)->GetBinContent(hists.at(0)->GetNbinsX()))
      ? maxabs = hists.at(1)->GetBinContent(hists.at(1)->GetNbinsX())
      : maxabs = hists.at(0)->GetBinContent(hists.at(0)->GetNbinsX());
  if(maxabs<0)leg = new TLegend(0.6,0.77-gPad->GetTopMargin(),0.97-gPad->GetRightMargin(),0.97-gPad->GetTopMargin());
  else leg = new TLegend(0.6,0.05+gPad->GetBottomMargin(),0.97-gPad->GetRightMargin(),0.3+gPad->GetBottomMargin());
  leg->SetBorderSize(0);leg->SetFillColor(kWhite);leg->SetFillStyle(1001);leg->SetTextAlign(12);leg->SetTextSize(0.05);leg->SetTextFont(42);
  leg->SetHeader("VELO + T-station");
  leg->AddEntry(hists.at(0),"2011 data","le1p");
  leg->AddEntry(hists.at(1),"2012 data","le1p");
  leg->AddEntry(Comb_hist,"Combined","le1p");

  leg->Draw();

  temp = myconfig->get_dumpdir() + "Plots/1112Combined_asymmetries_" + myconfig->get_bw() + ".pdf";

  cVT->SaveAs(temp);
  delete cVT;
  return;
}

void make_a_plot(TH1D *Up_hist, TH1D *Down_hist, configuration *myconfig){

  TCanvas* cVT = new TCanvas("cVT","Canvas",10,10,1280,960) ;
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetLeftMargin(0.15);

  Down_hist->SetLineWidth(2);
  Down_hist->SetLineColor(kBlue+2);
  Down_hist->SetMarkerSize(1.5);
  Down_hist->SetMarkerColor(kBlue+2);

  Up_hist->SetLineWidth(2);
  Up_hist->SetLineColor(kOrange+1);
  Up_hist->SetMarkerSize(1.5);
  Up_hist->SetMarkerColor(kOrange+1);

  Up_hist->GetXaxis()->SetLabelSize(0.05);
  Up_hist->GetYaxis()->SetLabelSize(0.05);
  Up_hist->GetXaxis()->SetTitleSize(0.05);
  Up_hist->GetYaxis()->SetTitleSize(0.05);

  Up_hist->Draw("e1p");

  TLine zeroline(Up_hist->GetXaxis()->GetXmin(),0,Up_hist->GetXaxis()->GetXmax(),0);
  zeroline.SetLineStyle(kDashed);
  zeroline.SetLineColor(1);
  zeroline.SetLineWidth(0.8);
  zeroline.Draw();

  Down_hist->Draw("e1psame");
  Up_hist->Draw("e1psame");

  TPaveText *blank_box = new TPaveText(1.003-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();

  TLegend *leg;
  double maxabs;
  fabs(Up_hist->GetBinContent(Up_hist->GetNbinsX()))>fabs(Down_hist->GetBinContent(Down_hist->GetNbinsX()))
      ? maxabs = Up_hist->GetBinContent(Up_hist->GetNbinsX())
      : maxabs = Down_hist->GetBinContent(Down_hist->GetNbinsX());
  if(maxabs<0)leg = new TLegend(0.6,0.77-gPad->GetTopMargin(),0.97-gPad->GetRightMargin(),0.97-gPad->GetTopMargin());
  else leg = new TLegend(0.6,0.05+gPad->GetBottomMargin(),0.97-gPad->GetRightMargin(),0.3+gPad->GetBottomMargin());
  leg->SetBorderSize(0);leg->SetFillColor(kWhite);leg->SetFillStyle(1001);leg->SetTextAlign(12);leg->SetTextSize(0.05);leg->SetTextFont(42);
  leg->SetHeader("VELO + T-station");
  leg->AddEntry(Up_hist,"Magnet Up","le1p");
  leg->AddEntry(Down_hist,"Magnet Down","le1p");

  leg->Draw();

  temp = myconfig->get_dumpdir() + "Plots/UpDown_asymmetries_" + myconfig->get_bw() + ".pdf";

  cVT->SaveAs(temp);
  delete cVT;
  return;
}

void make_a_plot(vector <TH1D*> hists, configuration *myconfig){

  if(hists.size() > 2){cout << "Sorry, but plotting is only implemented for 11 and 12 data" << endl; return;}

  TCanvas* cVT = new TCanvas("cVT","Canvas",10,10,1280,960) ;
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetLeftMargin(0.15);

  hists.at(0)->SetLineWidth(2);
  hists.at(0)->SetLineColor(kGreen+3);
  hists.at(0)->SetMarkerSize(1.5);
  hists.at(0)->SetMarkerColor(kGreen+3);

  hists.at(1)->SetLineWidth(2);
  hists.at(1)->SetLineColor(kRed+2);
  hists.at(1)->SetMarkerSize(1.5);
  hists.at(1)->SetMarkerColor(kRed+2);

  hists.at(0)->GetXaxis()->SetLabelSize(0.05);
  hists.at(0)->GetYaxis()->SetLabelSize(0.05);
  hists.at(0)->GetXaxis()->SetTitleSize(0.05);
  hists.at(0)->GetYaxis()->SetTitleSize(0.05);

  hists.at(0)->Draw("e1p");

  TLine zeroline(hists.at(0)->GetXaxis()->GetXmin(),0,hists.at(0)->GetXaxis()->GetXmax(),0);
  zeroline.SetLineStyle(kDashed);
  zeroline.SetLineColor(1);
  zeroline.SetLineWidth(0.8);
  zeroline.Draw();

  hists.at(1)->Draw("e1psame");
  hists.at(0)->Draw("e1psame");

  TPaveText *blank_box = new TPaveText(1.003-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();

  TLegend *leg;
  double maxabs;
  fabs(hists.at(0)->GetBinContent(hists.at(0)->GetNbinsX()))>fabs(hists.at(1)->GetBinContent(hists.at(1)->GetNbinsX()))
      ? maxabs = hists.at(0)->GetBinContent(hists.at(0)->GetNbinsX())
      : maxabs = hists.at(1)->GetBinContent(hists.at(1)->GetNbinsX());
  if(maxabs<0)leg = new TLegend(0.6,0.77-gPad->GetTopMargin(),0.97-gPad->GetRightMargin(),0.97-gPad->GetTopMargin());
  else leg = new TLegend(0.6,0.05+gPad->GetBottomMargin(),0.97-gPad->GetRightMargin(),0.3+gPad->GetBottomMargin());
  leg->SetBorderSize(0);leg->SetFillColor(kWhite);leg->SetFillStyle(1001);leg->SetTextAlign(12);leg->SetTextSize(0.05);leg->SetTextFont(42);
  leg->SetHeader("VELO + T-station");
  leg->AddEntry(hists.at(0),"2011 data","le1p");
  leg->AddEntry(hists.at(1),"2012 data","le1p");

  leg->Draw();

  temp = myconfig->get_dumpdir() + "Plots/1112_asymmetries_" + myconfig->get_bw() + ".pdf";

  cVT->SaveAs(temp);
  delete cVT;
  return;
}

void make_a_plot(TH1D *Comb_hist, double cv, double ecv, configuration *myconfig){

  TCanvas* cVT = new TCanvas("cVT","Canvas",10,10,1280,960) ;
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetLeftMargin(0.15);

  Color_t datacol = kBlack;
  Color_t cvcol = kAzure+2;
  Color_t bandcol = kAzure-9;
  Comb_hist->SetLineWidth(2);
  Comb_hist->SetLineColor(datacol);
  Comb_hist->SetMarkerSize(1.5);
  Comb_hist->SetMarkerColor(datacol);
  Comb_hist->GetXaxis()->SetLabelSize(0.05);
  Comb_hist->GetYaxis()->SetLabelSize(0.05);
  Comb_hist->GetXaxis()->SetTitleSize(0.05);
  Comb_hist->GetYaxis()->SetTitleSize(0.05);
  Comb_hist->Draw("e1p");


  double x[2]={Comb_hist->GetXaxis()->GetXmin(),Comb_hist->GetXaxis()->GetXmax()};
  double y[2]={cv,cv};
  double ex[2]={0};
  double ey[2]={ecv,ecv};
  TGraphErrors central_value(2,x,y,ex,ey);
  central_value.SetName("cval");
  central_value.SetLineColor(cvcol);
  central_value.SetFillColor(bandcol);
  central_value.SetLineWidth(1.2);
  central_value.Draw("e3same");
  central_value.Draw("c");

  /*TF1 *const_fit = new TF1("const_fit","pol0(0)",Comb_hist->GetXaxis()->GetXmin(),Comb_hist->GetXaxis()->GetXmax());
    const_fit->SetParNames("center");
    Comb_hist->Fit("const_fit","N", "");

    double fy[2]={const_fit->GetParameter(0),const_fit->GetParameter(0)};
    double fey[2]={const_fit->GetParError(0),const_fit->GetParError(0)};
    TGraphErrors fit_value(2,x,fy,ex,fey);
    fit_value.SetName("cval");
    fit_value.SetLineColor(kOrange-3);
    fit_value.SetFillColor(kOrange-4);
    fit_value.SetLineWidth(1.2);
    fit_value.Draw("e3same");
    fit_value.Draw("c");*/

  TLine zeroline(Comb_hist->GetXaxis()->GetXmin(),0,Comb_hist->GetXaxis()->GetXmax(),0);
  zeroline.SetLineStyle(kDashed);
  zeroline.SetLineColor(1);
  zeroline.SetLineWidth(0.8);
  zeroline.Draw();
  Comb_hist->Draw("e1psame");
  Comb_hist->Draw("axissame");

  TPaveText *blank_box = new TPaveText(1.003-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();

  TLegend *leg;
  if(Comb_hist->GetBinContent(Comb_hist->GetNbinsX())<0)leg = new TLegend(0.6,0.72-gPad->GetTopMargin(),0.97-gPad->GetRightMargin(),0.97-gPad->GetTopMargin());
  else leg = new TLegend(0.6,0.05+gPad->GetBottomMargin(),0.97-gPad->GetRightMargin(),0.3+gPad->GetBottomMargin());
  leg->SetBorderSize(0);leg->SetFillColor(kWhite);leg->SetFillStyle(1001);leg->SetTextAlign(12);leg->SetTextSize(0.05);leg->SetTextFont(42);
  leg->SetHeader("VELO + T Station");
  temp = "A_{track}(" + myconfig->get_dim1_label() + ")";
  leg->AddEntry(Comb_hist,temp,"le1p");
  temp = "#scale[0.8]{#int} d" + myconfig->get_dim1_label() + "A_{track}(" + myconfig->get_dim1_label() + ")";
  leg->AddEntry("cval",temp,"fl");

  leg->Draw();

  temp = myconfig->get_dumpdir() + "Plots/Combined_" + myconfig->get_bw() + ".pdf";
  cVT->SaveAs(temp);
  delete cVT;
  return;
}

