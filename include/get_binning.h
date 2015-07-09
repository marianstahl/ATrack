#include "configuration.h"
extern TString temp;
void get_binning(configuration *myconfig){
  myconfig->clear_binning();
  ifstream dim1bins_in;
  temp = "Binning/" + myconfig->get_binning_str() + "_variable_binning_" + myconfig->get_dim1();
  if(myconfig->get_verbosity() > 1)cout << "Reading " << temp << endl;
  dim1bins_in.open(temp);
  ifstream dim2bins_in;
  temp = "Binning/" + myconfig->get_binning_str() + "_variable_binning_" + myconfig->get_dim2();
  if(!myconfig->is_1D())dim2bins_in.open(temp);
  if(!dim1bins_in || (!dim2bins_in && !myconfig->is_1D())){cout << "Variable binning not found!" << endl;return;}

  while (dim1bins_in.good()) {
    double binedge_dim1;
    dim1bins_in >> binedge_dim1;
    myconfig->fill_dim1_binvector(binedge_dim1);
  }
  dim1bins_in.close();
  //fill vector with whole dim1 range to be able to initialize histgrams correctly
  if(myconfig->is_1D()){myconfig->fill_dim2_binvector(myconfig->get_dim1_lo());myconfig->fill_dim2_binvector(myconfig->get_dim1_hi());}
  else{
    while (dim2bins_in.good()) {
      double binedge_dim2;
      dim2bins_in >> binedge_dim2;
      myconfig->fill_dim2_binvector(binedge_dim2);
    }
    dim2bins_in.close();
  }
}
