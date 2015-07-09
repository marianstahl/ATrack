#include "configuration.h"

extern TString temp;

void make_latex_table(const vector<three_vectors> &results, configuration *myconfig){

  ofstream table;
  temp = myconfig->get_dumpdir() + "Latex";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);
  temp += "/MuPiTrackingAsymmetries_" + myconfig->get_bw() + "_VT";
  table.open(temp);
  table << "\\begin{table}[]" << endl;
  table << "\t\\begin{center}" << endl;
  //table << "\t\t\\resizebox*{!}{\textheight}{" << endl;
  table << "\t\t\\begin{tabular}{|l| r@{ $\\pm$ }l | r@{ $\\pm$ }l | r@{ $\\pm$ }l |}" << endl;
  table << "\t\t\t\\hline" << endl;
  table << "\t\t\t\\multirow{2}{*}{Sample} & \\multicolumn{6}{c|}{$\\Amupi$ [\%]} \\\\ \\cline{2-7}" << endl;
  for(int ypol = -1; ypol < myconfig->get_nsamples(); ypol++){
    if(ypol > -1)try {if(!myconfig->set_sample(ypol))throw 3;}catch(int e){cout << "Exception number " << e << "\n Something wrong with setter of sample, channel, method! Terminating..." << endl;terminate();}
    for(int mode = 0; mode < myconfig->get_nchannels(); mode++){
      try {if(!myconfig->set_channel(mode))throw 3;}catch(int e){cout << "Exception number " << e << "\n Something wrong with setter of sample, channel, method! Terminating..." << endl;terminate();}
      if(ypol == -1){
        if(mode == 0){temp.Form("\t\t\t  & \\multicolumn{2}{c|}{%s}",myconfig->get_channeltexlabel().Data());table << temp;}
        else if(mode < myconfig->get_nchannels() - 1){temp.Form("& \\multicolumn{2}{c|}{%s}",myconfig->get_channeltexlabel().Data());table << temp;}
        else {temp.Form("& \\multicolumn{2}{c||}{%s} \\\\ \\hline",myconfig->get_channeltexlabel().Data());table << temp << endl;}
      }
      else{
        if(mode == 0){temp.Form("\t\t\t %s & %+.2f & %.2f",myconfig->get_samplelabel().Data(),results.at(mode).cv_VT.at(ypol),results.at(mode).unc_VT.at(ypol));table << temp;}
        else if(mode < myconfig->get_nchannels() - 1){temp.Form("& %+.2f & %.2f",results.at(mode).cv_VT.at(ypol),results.at(mode).unc_VT.at(ypol));table << temp;}
        else {temp.Form("& %+.2f & %.2f \\\\",results.at(mode).cv_VT.at(ypol),results.at(mode).unc_VT.at(ypol));table << temp << endl;}
      }
    }
  }
  table << "\t\t\t\\hline" << endl;
  table << "\t\t\\end{tabular}" << endl;
  table << "\t\t\\caption[Default $\\Amupi$ from $\\Jpsimupmum$]{$\\Amupi$ from $\\Jpsimupmum$ for different data taking periods, magnet polarities and Dalitz regions.}"<< endl;
  table << "\t\t\\label{tab:mupi_tr_asym_" << myconfig->get_bw() << "}" << endl;
  table << "\t\\end{center}" << endl;
  table << "\\end{table}" << endl;
  table.close();
  return;
}

void make_python_dictionary(const vector<three_vectors> &results, configuration *myconfig){

  ofstream pyfile;
  temp = myconfig->get_dumpdir() + "Python";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);
  temp += "/Marian_" + myconfig->get_bw() + ".py";
  pyfile.open(temp);
  pyfile << "Marian = {'AmupiTrack':" << endl;
  for(int ypol = 0; ypol < myconfig->get_nsamples(); ypol++){
    try {if(!myconfig->set_sample(ypol))throw 3;}catch(int e){cout << "Exception number " << e << "\n Something wrong with setter of sample, channel, method! Terminating..." << endl;terminate();}
    pyfile << "\t'"<< myconfig->get_sample() <<"' : {" << endl;
    for(int mode = 0; mode < myconfig->get_nchannels() ; mode++){
      myconfig->set_channel(mode);
      pyfile << "\t\t'"<< myconfig->get_channelpylabel() <<"' : [" << results.at(mode).cv_VT.at(ypol) << ", " << results.at(mode).unc_VT.at(ypol) << ", 0.0]";
      if(mode < myconfig->get_nchannels() -1 )pyfile << "," << endl;
      else pyfile << endl;
    }
    if(ypol < myconfig->get_nsamples() - 1 )pyfile << "\t}," << endl;
    else pyfile << "\t}\n}" << endl;
  }
  pyfile.close();
  return;
}
