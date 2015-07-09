#include "TStyle.h"
#include "TROOT.h"
#include "TColor.h"
#include "TGaxis.h"
#include "Riostream.h"

using namespace std;

void MyStyle(){
    const int NCont = 99;
    gStyle->SetNumberContours(NCont);
    const Int_t nRGBs = 8;
    Double_t stops[nRGBs] = { 0.0, 0.18, 0.36, 0.5, 0.57, 0.66, 0.86, 1.00 };
    Double_t red[nRGBs]   = { 1.0, 0.0, 0.00, 0.0, 0.8, 1.0, 1.00, 0.25 };
    Double_t green[nRGBs] = { 1.0, 0.0, 0.81, 1.0, 1.0 ,1.0, 0.20, 0.00 };
    Double_t blue[nRGBs]  = { 1.0, 0.5, 1.00, 0.0, 0.1 ,0.0, 0.00, 0.00 };
    Int_t FI = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, NCont);
    Int_t MyPalette[NCont];
    for (int i=0;i<NCont;i++) MyPalette[i] = FI+i;
    gStyle->SetPalette(NCont, MyPalette);

    TGaxis::SetMaxDigits(3);
    gStyle->SetMarkerColor(1);
    gStyle->SetLineColor(1);
    gStyle->SetFillColor(1);
    gStyle->SetOptStat("");
    gStyle->SetMarkerStyle(20);
    gStyle->SetStatFont(42);
    gStyle->SetStatStyle(4000);
    gStyle->SetStatBorderSize(0);
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadRightMargin(0.08);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetTitleXOffset(1.2);
    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYOffset(1.3);
    gStyle->SetTitleYSize(0.05);
    gStyle->SetTitleFont(42);
    gStyle->SetTitleFont(42,"X");
    gStyle->SetTitleFont(42,"Y");
    gStyle->SetTitleFont(42,"Z");
    gStyle->SetLabelSize(0.05);
    gStyle->SetLabelSize(0.05,"Y");
    gStyle->SetLabelSize(0.05,"Z");
    gStyle->SetLabelFont(42);
    gStyle->SetLabelFont(42,"Y");
    gStyle->SetLabelFont(42,"Z");
    //gStyle->SetLabelOffset(0.007,"Y");
    gStyle->SetLegendFont(42);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetHistLineColor(1);
    gStyle->SetCanvasDefW(1280);
    gStyle->SetCanvasDefH(960);



    cout << "	" << endl;
    cout << "       Ei Gude!      "<< endl;
    cout << "  _________________  "<< endl;
    cout << "          |          "<< endl;
    cout << "          |          "<< endl;
    cout << "          |          "<< endl;
    cout << "         WWW         "<< endl;
    cout << "        <ô¿ô>        "<< endl;
    cout << "         '-'         "<< endl;
    cout << "          X          "<< endl;
    cout << "         '|'         "<< endl;
    cout << "        ' | '        "<< endl;
    cout << "       O  |  O       "<< endl;
    cout << "         ' '         "<< endl;
    cout << "        '   '        "<< endl;
    cout << "       '     '       "<< endl;
    cout << "      ()     ()      "<< endl;
}
