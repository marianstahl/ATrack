/*!

  \date July, 1st 2015
  \author Marian Stahl
  \class configuration
  \brief Change settings here if needed. Line 125ff : hardcoded stuff such as directories, file names, variable names, fitparameters etc.

*/

#ifndef CONFIGURATION
#define CONFIGURATION
#include "TString.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <vector>
#include <functional>
using namespace std;

typedef vector<double> dvec;
typedef vector< vector<double> > dmtx;
typedef vector< vector< vector<double> > > dten3;

struct three_vectors{
  dvec combined_results;
  dvec cv_VT;
  dvec unc_VT;
};

struct FitParameters{
  double MJpsistart; double MJpsilo; double MJpsihi;
  double sigma1start; double sigma1lo; double sigma1hi;
  double wfstart; double wflo; double wfhi;
  double fixCBn; double fixwf;
  double CBalphastart; double CBalphalo; double CBalphahi;
  double chebastart; double chebalo; double chebahi;
  double chebbstart; double chebblo; double chebbhi;
  double chebcstart; double chebclo; double chebchi;
  double taustart; double taulo; double tauhi;
  double sfracstart; double sfraclo; double sfrachi;
  double Jpsisig1start; double Jpsisiglo; double Jpsisighi;
  double bkgstart; double bkglo; double bkghi;
  double effstart;double efflo;double effhi;
  double Astart;double Alo;double Ahi;
};

class configuration{
private:

  bool debug;
  int verb_lvl;

  TString version;
  TString buzzword;
  bool oneD;
  TString varbinning;
  bool simultaneous;
  int signalfitmodel;
  int shared_parameters;
  int backgroundfitmodel;
  double dim1_binlo;
  double dim1_binhi;
  double dim2_binlo;
  double dim2_binhi;

  TString part1;
  TString part2;
  TString dim1;
  TString dim2;
  double dim1_lo_limit;
  double dim1_hi_limit;
  double dim2_lo_limit;
  double dim2_hi_limit;

  function<TString ()> dumpdir;//see std::function
  TString taptupledir;
  TString sigtupledir;
  TString sweighttupledir;
  function<TString ()> sigtuplename;
  function<TString ()> sigtreename;
  function<TString ()> sweighttuplename;
  TString sweighttreename;
  TString sweightvarname;

  int fitbinning;
  double Jpsi_M_min;
  double Jpsi_M_max;

  FitParameters fitparams;

  TString sample;
  TString samplelabel;
  TString sigsample;
  TString method;
  TString methodlabel;
  TString channel;
  TString channellabel;
  TString channeltexlabel;
  TString channelpylabel;

  vector<TString> samples;
  vector<TString> samplelabels;
  vector<TString> sigsamples;
  vector<TString> methods;
  vector<TString> methodlabels;
  vector<TString> channels;
  vector<TString> channellabels;
  vector<TString> channeltexlabels;
  vector<TString> channelpylabels;

  int nsamples;
  int nmethods;
  int nchannels;

  vector<double> dim1_binedge_vector;
  vector<double> dim2_binedge_vector;

public:

  configuration(){

    debug = false;
    verb_lvl = 1;

    //these things are configurable by the jobname
    buzzword = "simple_nsim_dg_wf_cheb3_PT_ETA";
    oneD = false;
    varbinning = "simple";
    simultaneous = false;
    signalfitmodel = 1;//0-gaussian, 1-double gaussian, 2-crystal-ball
    shared_parameters = 1;
    backgroundfitmodel = 3;//0-exponential, 1-chebychev 1st order, 2-2nd order, 3-3rd order
    dim1_binlo = 0.0;
    dim1_binhi = 0.0;
    dim2_binlo = 0.0;
    dim2_binhi = 0.0;

    //the following things are hardcoded
    part1 = "Mu";//Name of first particle as in signal tuple
    part2 = "Pi";//Name of second particle as in signal tuple
    dim1 = "PT";//Name of first dimension as in signal tuple
    dim2 = "ETA";//Name of second dimension as in signal tuple
    dim1_lo_limit = 0;//lower limit of dim1 variable
    dim1_hi_limit = 100000;//upper limit of dim1 variable
    dim2_lo_limit = 1.5;//lower limit of dim2 variable
    dim2_hi_limit = 5.5;//upper limit of dim2 variable

    //This way of initialization allows that the TString returned from function<TString ()> dumpdir to change whenever dumpdir() is called.
    //In this way, dumpdir only needs to be initialied once here and can be used later even though version is not set yet
    dumpdir = [&] () -> TString{return "/afs/cern.ch/work/m/mstahl/public/dumpdir/tracking_asymmetry_raw/" + version + "/";};//location where plots, results etc. will be dumped
    taptupledir = "/afs/cern.ch/work/m/mstahl/public/";//location of tag-and-probe tuples
    sigtupledir = "/afs/cern.ch/work/j/jadevrie/public/asls/ntuples/cutTuples_080515/";//location of signal tuples
    sweighttupledir = "/afs/cern.ch/work/m/mstahl/public/asls_sWeights/";//location of signal sWeight tuples
    sigtuplename = [&] () -> TString{return "tup_tree_" + sigsample + "_" + channel +"_cut.root";};//name of signal tuple
    sigtreename = [&] () -> TString{return "tree_" + sigsample + "_" + channel + "_cut";};//name of signal tree
    sweighttuplename = [&] () -> TString{return "sWeights_tree_" + sigsample + "_" + channel +"_cut.root";};//name of sWeight tuple (leave empty if nos Weights are used or if they are provided in signal tree)
    sweighttreename = "swTree";//name of sWeight tree
    sweightvarname = "sWeight";//name of sWeight variable

    samples = {"2011MagUp","2011MagDown","2012MagUp","2012MagDown"};//Sample names of the tag-and-probe samples and for plot naming. Very important to keep MagUp even and MagDown odd!(counting 0,1,2,3)
    samplelabels = {"2011 Magnet Up","2011 Magnet Down","2012 Magnet Up","2012 Magnet Down"};//just for plotting puposes
    sigsamples = {"20r1_Up","20r1_Down","20_Up","20_Down"};//names of the signal samples
    methods = {"VELO","TStation","Long"};//this should not be changed(if it's changed, the source code needs to be rewritten)
    methodlabels = {"VELO method","T-station method","Long method"};//just for plotting puposes
    channels = {"PhiPi","KStarK","KKPiNR"};//channels, bins etc... Basically a placeholder for another variable in wich the samples are split
    channellabels = {"#phi#pi^{#pm}","K*(892)K^{#pm}","(K^{+}K^{-}#pi^{#pm})_{NR}"};//just for plotting puposes
    channeltexlabels = {"$\\PhiPi$","$K^\\ast K$","$(KK\\pi)_{\\mathrm{NR}}$"};//for tex table label
    channelpylabels = {"PhiPi","KStK","NR"};//labels for output python dictionary

    fitbinning = 98;
    Jpsi_M_min = 2708.0;
    Jpsi_M_max = 3492.0;

    // (common) Mean Jpsi
    fitparams.MJpsistart = 3098; fitparams.MJpsilo = 3092; fitparams.MJpsihi = 3104;
    // sigma's
    fitparams.sigma1start = 20; fitparams.sigma1lo = 10; fitparams.sigma1hi = 40;
    fitparams.wfstart = 1.5; fitparams.wflo = 1.3; fitparams.wfhi = 2.1;
    fitparams.fixCBn = 2.0; fitparams.fixwf = 1.6;
    fitparams.CBalphastart = 2.1; fitparams.CBalphalo = 0.5; fitparams.CBalphahi = 6.0;
    //Background
    fitparams.chebastart = -0.3; fitparams.chebalo = -0.9; fitparams.chebahi = 0.9;
    fitparams.chebbstart = -0.03; fitparams.chebblo = -0.3; fitparams.chebbhi = 0.3;
    fitparams.chebcstart = 0.002; fitparams.chebclo = -0.1; fitparams.chebchi = 0.1;
    fitparams.taustart = -0.001; fitparams.taulo = -0.1; fitparams.tauhi = 0.1;
    //normalizations
    fitparams.sfracstart = 0.5; fitparams.sfraclo = 0.25; fitparams.sfrachi = 0.75;
    fitparams.Jpsisig1start = 15000; fitparams.Jpsisiglo = 500; fitparams.Jpsisighi = 3.5e+5;
    fitparams.bkgstart = 10000; fitparams.bkglo = 100; fitparams.bkghi = 2e+6;
    //efficiency
    fitparams.effstart = 0.97 ;fitparams.efflo = 0.8 ;fitparams.effhi = 1.0;
    //asymmetry(for simultaneous fits)
    fitparams.Astart = 0.0 ;fitparams.Alo = -0.05 ;fitparams.Ahi = 0.05;

    nsamples = static_cast<int>(samples.size());
    nmethods = static_cast<int>(methods.size());
    nchannels = static_cast<int>(channels.size());

    cout << "Template configuration set:" << endl
         << "\t simple variable binning" << endl
         << "\t sequential fits" << endl
         << "\t double gaussian signal" << endl
         << "\t shared parameter mode 1" << endl
         << "\t 3rd order chebychev background" << endl
         << "\t PT as first dimension" << endl
         << "\t ETA as second dimension" << endl
         << "No sample, tap method and Dalitz region defined yet!" << endl;
  }  

  configuration(const TString jobname) : configuration(){

    TString vers = (TString)((TObjString*)(jobname.Tokenize("-")->At(0)))->String();
    TString systematic = (TString)((TObjString*)(jobname.Tokenize("-")->At(1)))->String();
    this->version = vers;

    cout << endl << string(128, '*') << endl << "Running systematic check: " << systematic << endl << string(128, '*') << endl;
    this->buzzword = systematic;
    //split the string at every "_"
    TObjArray* strar = systematic.Tokenize("_");
    if(!(strar->GetEntries() == 7))
      cout << "\033[0;31mSystematic string cannot be parsed, running with default configuration\033[0m" << endl;

    if(strar->GetEntries() == 7){

      this->varbinning = ((TObjString*)strar->At(0))->String();
      if(this->varbinning.CompareTo("simple") != 0)cout << "default changed to " << this->varbinning << " binning" << endl;

      if(((TObjString*)strar->At(1))->String().CompareTo("nsim") == 0)
        this->simultaneous = false;
      else if(((TObjString*)strar->At(1))->String().CompareTo("sim") == 0){
        this->simultaneous = true;
        cout << "default changed to simultaneous fits" << endl;
      }
      else cout << "Option for simultaneous fits could not be parsed => Running simultaneous fits" << endl;

      if(((TObjString*)strar->At(2))->String().CompareTo("sg") == 0){
        this->signalfitmodel = 0;
        cout << "default changed to single gaussian signal model" << endl;
      }
      else if(((TObjString*)strar->At(2))->String().CompareTo("dg") == 0)
        this->signalfitmodel = 1;
      else if(((TObjString*)strar->At(2))->String().CompareTo("cb") == 0){
        this->signalfitmodel = 2;
        cout << "default changed to crystal ball signal model" << endl;
      }
      else cout << "Option for signal model could not be parsed => Running double gaussian fits" << endl;

      if(((TObjString*)strar->At(3))->String().CompareTo("def") == 0){
        this->shared_parameters = 0;
        cout << "default changed to shared parameter mode 0: fixed widthfactor" << endl;
      }
      else if(((TObjString*)strar->At(3))->String().CompareTo("wf") == 0)
        this->shared_parameters = 1;
      else if(((TObjString*)strar->At(3))->String().CompareTo("sgf") == 0){
        this->shared_parameters = 2;
        cout << "default changed to shared parameter mode 2: single gaussian for fail" << endl;
      }
      else cout << "Option for shared parameter mode could not be parsed => Running default mode" << endl;

      if(((TObjString*)strar->At(4))->String().CompareTo("exp") == 0){
        this->backgroundfitmodel = 0;
        cout << "default changed to exponential background model" << endl;
      }
      else if(((TObjString*)strar->At(4))->String().CompareTo("cheb1") == 0){
        this->backgroundfitmodel = 1;
        cout << "default changed to 1st order chebychev background model" << endl;
      }
      else if(((TObjString*)strar->At(4))->String().CompareTo("cheb2") == 0){
        this->backgroundfitmodel = 2;
        cout << "default changed to 2nd order chebychev background model" << endl;
      }
      else if(((TObjString*)strar->At(4))->String().CompareTo("cheb3") == 0)
        this->backgroundfitmodel = 3;
      else cout << "Option for background model could not be parsed => Running 3rd order chebychev fits" << endl;

      this->dim1 = ((TObjString*)strar->At(5))->String();
      if(this->dim1.CompareTo("PT") != 0)cout << "default first dimenstion changed to " << this->dim1 << endl;
      if(this->dim1.CompareTo("ETA") == 0){
        this->dim1_lo_limit = 1.5; this->dim1_hi_limit = 5.5;
      }
      if(this->dim1.CompareTo("P") == 0)this->dim1_hi_limit = 1e+6;

      this->dim2 = ((TObjString*)strar->At(6))->String();
      if(this->dim2.CompareTo("1D") == 0){cout << "1D binned method chosen" << endl;oneD = true;}
      else if(this->dim2.CompareTo("ETA") != 0)cout << "default second dimenstion changed to " << this->dim2 << endl;
      if(this->dim2.CompareTo("PT") == 0){
        this->dim2_lo_limit = 0; this->dim2_hi_limit = 1e+5;
      }
      if(this->dim2.CompareTo("P") == 0){
        this->dim2_lo_limit = 0; this->dim2_hi_limit = 1e+6;
      }
    }
  }

  configuration(const int dataset, const TString jobname) : configuration(jobname){

    //this sets all default values and they can be used from now on
    int meth = static_cast<int>((dataset-0.1)/static_cast<double>(nsamples));
    int ypol = static_cast<int>((dataset-1)%nsamples);

    bool good_job;
    try{
      good_job = this->set_sample(ypol);
      good_job &= this->set_method(meth);
      if(!good_job)throw 1;
    }
    catch(int e){
      cout << "Exception number " << e << "\n Please specify a valid subjob number! Terminating..." << endl; terminate();
    }

    if(verb_lvl > 0){
      cout << "Dump directory\t: " << dumpdir() << endl;
      cout << "Signal sample\t: " << sigtuplename() << endl;
      cout << "Signal tree\t: " << sigtreename() << endl;
      cout << "SWeight sample\t: " << sweighttuplename() << endl;
    }

  }

  ~configuration(){}

  TString short_systematic(){

    TString bw;

    if(this->varbinning.CompareTo("simple") != 0){
      bw += " " + varbinning;
    }

    if(this->simultaneous){
      bw += " simultaneous fits";
    }

    if(this->signalfitmodel == 0 ){
      bw += " single Gaussian";
    }
    if(this->signalfitmodel == 2 ){
      bw += " Crystal Ball";
    }

    if(this->shared_parameters == 0){
      bw += " fixed width factor";
    }
    if(this->shared_parameters == 2){
      bw += " sG in fail";
    }

    if(this->backgroundfitmodel == 0){
      bw += " exponential bg";
    }
    if(this->backgroundfitmodel == 1){
      bw += " 1st order Chebychev";
    }
    if(this->backgroundfitmodel == 2){
      bw += " 2nd order Chebychev";
    }

    if(this->dim1.CompareTo("PT") != 0 || this->dim2.CompareTo("ETA") != 0){
      bw += " " + this->dim1 + " " + this->dim2;
    }

    if(bw.Length() == 0){
      bw = "default";
    }

    return bw;
  }

  void set_verbosity(int lvl, bool dbg = 0){
    verb_lvl = lvl;
    debug = dbg;
  }

  bool is_debug(){return debug;}
  int get_verbosity(){return verb_lvl;}
  double get_dim1_lo(){return dim1_lo_limit;}
  double get_dim2_lo(){return dim2_lo_limit;}
  double get_dim1_hi(){return dim1_hi_limit;}
  double get_dim2_hi(){return dim2_hi_limit;}

  TString get_bw(){return buzzword;}

  TString get_binning_str(){return varbinning;}
  bool is_simultaneous(){return simultaneous;}
  int get_sigmodel(){return signalfitmodel;}
  int get_shared_params(){return shared_parameters;}
  int get_bkgmodel(){return backgroundfitmodel;}

  TString get_firstParticle(){return part1;}
  TString get_secondParticle(){return part2;}
  TString get_dim1(){return dim1;}
  TString get_dim2(){return dim2;}
  bool is_1D(){return oneD;}

  TString get_dumpdir(){return dumpdir();}
  TString get_sigtupledir(){return sigtupledir;}
  TString get_taptupledir(){return taptupledir;}
  TString get_sweighttupledir(){return sweighttupledir;}
  TString get_sigtuplename(){return sigtuplename();}
  TString get_sigtreename(){return sigtreename();}
  TString get_sweighttuplename(){return sweighttuplename();}
  TString get_sweighttreename(){return sweighttreename;}
  TString get_sweightvarname(){return sweightvarname;}

  template<class intlike>
  bool set_sample(intlike ypol){
    if(ypol < samples.size()){
      sample = samples.at(ypol);
      samplelabel = samplelabels.at(ypol);
      sigsample = sigsamples.at(ypol);
      return true;
    }
    else{cout << "Valid year & polartiy combinations:";
      for(unsigned int i = 0; i < samples.size(); i++){
        if(i == samples.size() - 1) cout << " " << i << " - " << samples.at(i) << "." << endl;
        else cout << " " << i << " - " << samples.at(i) << ", ";
      }
      return false;
    }
  }
  bool set_sample(TString new_sample){
    sample = new_sample;
    bool is_valid = false;
    for(auto it_sample : samples){
      is_valid |= new_sample.CompareTo(it_sample) == 0;
    }
    return is_valid;
  }

  template<class intlike>
  bool set_method(intlike meth){
    if(meth == 0){Jpsi_M_min = 2904.0; Jpsi_M_max = 3296.0;}
    if(meth < methods.size()){
      method = methods.at(meth);
      methodlabel = methodlabels.at(meth);
      return true;
    }
    else {cout << "Valid methods:";
      for(unsigned int i = 0; i < methods.size(); i++){
        if(i == methods.size() - 1) cout << " " << i << " - " << methods.at(i) << "." << endl;
        else cout << " " << i << " - " << methods.at(i) << ", ";
      }
      return false;
    }
  }
  bool set_method(TString new_method){
    method = new_method;
    bool is_valid = false;
    for(auto it_method : methods){
      is_valid |= new_method.CompareTo(it_method) == 0;
    }
    return is_valid;
  }

  template<class intlike>
  bool set_channel(intlike mode){
    if(mode < channels.size()){
      channel = channels.at(mode);
      channellabel = channellabels.at(mode);
      return true;
    }
    else{cout << "Valid channels:";
      for(unsigned int i = 0; i < channels.size(); i++){
        if(i == channels.size() - 1) cout << " " << i << " - " << channels.at(i) << "." << endl;
        else cout << " " << i << " - " << channels.at(i) << ", ";
      }
      return false;
    }
  }
  bool set_channel(TString new_channel){
    channel = new_channel;
    bool is_valid = false;
    for(auto it_channel : channels){
      is_valid |= new_channel.CompareTo(it_channel) == 0;
    }
    return is_valid;
  }

  TString get_sample(){return sample;}
  TString get_samplelabel(){return samplelabel;}
  TString get_sigsample(){return sigsample;}
  TString get_method(){return method;}
  TString get_methodlabel(){return methodlabel;}
  TString get_channel(){return channel;}
  TString get_channellabel(){return channellabel;}

  void set_fitbinning(int nbins){fitbinning = nbins; return;}
  int get_fitbinning(){return fitbinning;}
  double get_Jpsi_M_min(){return Jpsi_M_min;}
  double get_Jpsi_M_max(){return Jpsi_M_max;}

  TString get_dim1_alias(){
    if(dim1.CompareTo("PT") == 0)return "PT";
    else if(dim1.CompareTo("ETA") == 0)return "eta";
    else if(dim1.CompareTo("P") == 0)return "P";
    else return "";
  }

  TString get_dim1_label(){
    if(dim1.CompareTo("PT") == 0)return "p_{T} (GeV)";
    else if(dim1.CompareTo("ETA") == 0)return "#eta";
    else if(dim1.CompareTo("P") == 0)return "p (GeV)";
    else return "";
  }

  TString get_dim2_alias(){
    if(dim2.CompareTo("PT") == 0)return "PT";
    else if(dim2.CompareTo("ETA") == 0)return "eta";
    else if(dim2.CompareTo("P") == 0)return "P";
    else return "";
  }

  TString get_dim2_label(){
    if(dim2.CompareTo("PT") == 0)return "p_{T} (GeV)";
    else if(dim2.CompareTo("ETA") == 0)return "#eta";
    else if(dim2.CompareTo("P") == 0)return "p (GeV)";
    else return "";
  }

  void fill_dim1_binvector(double binedge){dim1_binedge_vector.push_back(binedge);}
  void fill_dim2_binvector(double binedge){dim2_binedge_vector.push_back(binedge);}
  vector<double> get_dim1_binedges(){return dim1_binedge_vector;}
  vector<double> get_dim2_binedges(){return dim2_binedge_vector;}
  void clear_binning(){dim1_binedge_vector.clear();dim2_binedge_vector.clear();return;}

  void set_dim1_bin(double binlo, double binhi){
    dim1_binlo = binlo;dim1_binhi = binhi;
  }
  template<class intlike>
  void set_dim1_bin(intlike binnumber){
    dim1_binlo = dim1_binedge_vector.at(binnumber);dim1_binhi = dim1_binedge_vector.at(binnumber + 1);
  }
  template<class intlike>
  void set_dim2_bin(intlike binnumber){
    dim2_binlo = dim2_binedge_vector.at(binnumber);dim2_binhi = dim2_binedge_vector.at(binnumber + 1);
  }

  double get_dim1_binlo(){return dim1_binlo;}
  double get_dim1_binhi(){return dim1_binhi;}
  double get_dim2_binlo(){return dim2_binlo;}
  double get_dim2_binhi(){return dim2_binhi;}

  FitParameters get_fitparams(){return fitparams;}

  int get_nsamples(){return nsamples;}
  int get_nmethods(){return nmethods;}
  int get_nchannels(){return nchannels;}
  vector<TString> get_channels(){return channels;}

  TString get_channeltexlabel(){return channeltexlabel;}
  TString get_channelpylabel(){return channelpylabel;}

};
#endif
