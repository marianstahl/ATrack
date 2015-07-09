#include <vector>
#include "TMath.h"
#include "TMatrixD.h"
#include "Riostream.h"

using namespace std;

inline double asymmetry_calc(double eff_f, double eff_fbar){
    return (eff_f-eff_fbar)/(eff_f+eff_fbar);
}

inline double asymmetry_unc_calc(double eff_f, double eff_fbar, double delta_eff_f, double delta_eff_fbar){
    return 2*sqrt((pow(eff_f,2)*pow(delta_eff_fbar,2)+pow(eff_fbar,2)*pow(delta_eff_f,2))/(pow((eff_f+eff_fbar),4)));
}
inline void pdgrounding(double &central_value, double &stat_unc, double &syst_unc){
    if((stat_unc == 0 || stat_unc != stat_unc || syst_unc != syst_unc) && syst_unc == 0 )return;
    double uncertainty = TMath::Max(stat_unc,syst_unc);
    int exponent = floor(TMath::Log10(uncertainty));
    double significant_digit = pow(10,-exponent+2)*uncertainty;
    if(significant_digit < 354){
        central_value = (round(pow(10,-exponent+1)*central_value))*pow(10,exponent-1);
        stat_unc = (round(pow(10,-exponent+1)*stat_unc))*pow(10,exponent-1);
        syst_unc = (round(pow(10,-exponent+1)*syst_unc))*pow(10,exponent-1);
    }
    else{
        central_value = (round(pow(10,-exponent)*central_value))*pow(10,exponent);
        stat_unc = (round(pow(10,-exponent)*stat_unc))*pow(10,exponent);
        syst_unc = (round(pow(10,-exponent)*syst_unc))*pow(10,exponent);
    }
}
inline void pdgrounding(double &central_value, double &stat_unc){
    if(stat_unc == 0 || stat_unc != stat_unc)return;
    double uncertainty = stat_unc;
    int exponent = floor(TMath::Log10(uncertainty));
    double significant_digit = pow(10,-exponent+2)*uncertainty;
    if(significant_digit < 354){
        central_value = (round(pow(10,-exponent+1)*central_value))*pow(10,exponent-1);
        stat_unc = (round(pow(10,-exponent+1)*stat_unc))*pow(10,exponent-1);
    }
    else{
        central_value = (round(pow(10,-exponent)*central_value))*pow(10,exponent);
        stat_unc = (round(pow(10,-exponent)*stat_unc))*pow(10,exponent);
    }
}
inline void pdgrounding(vector<double> &central_values, vector<double> &stat_uncertainties){
    if(central_values.size() != stat_uncertainties.size())return;
    for(unsigned int i = 0; i < central_values.size(); i++){
        if(stat_uncertainties.at(i) == 0 || stat_uncertainties.at(i) != stat_uncertainties.at(i))return;
        pdgrounding(central_values.at(i),stat_uncertainties.at(i));
    }
}
inline void rounding_to_ref(double &central_value, double &stat_unc, double reference){
    if((stat_unc == 0 || stat_unc != stat_unc || reference != reference) && reference == 0 )return;
    double uncertainty = reference;
    int exponent = floor(TMath::Log10(uncertainty));
    double significant_digit = pow(10,-exponent+2)*uncertainty;
    if(significant_digit < 354){
        central_value = (round(pow(10,-exponent+1)*central_value))*pow(10,exponent-1);
        stat_unc = (round(pow(10,-exponent+1)*stat_unc))*pow(10,exponent-1);
    }
    else{
        central_value = (round(pow(10,-exponent)*central_value))*pow(10,exponent);
        stat_unc = (round(pow(10,-exponent)*stat_unc))*pow(10,exponent);
    }
}

inline double summation(vector<double> x){
    double sum = 0;
    for(unsigned int i = 0; i < x.size(); i++){
        sum += x.at(i);
    }
    return sum;
}

inline double summation_unc(vector<double> delta_x){
    double unc = 0;
    for(unsigned int i = 0; i < delta_x.size(); i++){
        unc += pow(delta_x.at(i),2);
    }
    return sqrt(unc);
}

inline double combination(vector<double> x, vector<double> delta_x){
    double sum_w_i = 0, x_bar = 0;
    for(unsigned int i = 0; i < x.size(); i++){
        double w_i = 1/pow(delta_x.at(i),2);
        x_bar += x.at(i)*w_i;
        sum_w_i += w_i;
    }
    return x_bar/sum_w_i;
}

inline double combined_unc(vector<double> delta_x){
    double sum_w_i = 0;
    for(unsigned int i = 0; i < delta_x.size(); i++){
        sum_w_i += 1/pow(delta_x.at(i),2);
    }
    return 1/sqrt(sum_w_i);
}

inline double red_chi2_combination(vector<double> x, vector<double> delta_x){
    double x_bar = combination(x,delta_x);
    double chi2 = 0;
    for(unsigned int i = 0; i < x.size(); i++){
        double w_i = 1/pow(delta_x.at(i),2);
        chi2 += w_i*pow((x_bar-x.at(i)),2);
    }
    return chi2/(x.size()-1);
}

double corr_comb(vector<double> x, vector<double> delta_x, vector< vector<double> > correlations){
    int dimension = static_cast<int>(delta_x.size());
    TMatrixD covariance_matrix(dimension,dimension);
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            covariance_matrix[i][j] = correlations[i][j]*delta_x[i]*delta_x[j];
        }
    }
    TMatrixD inv_cov = covariance_matrix.Invert();
    double sum_inv_cov = 0;
    for(int j = 0; j < dimension; j++){
        for(int k = 0; k < dimension; k++){
            sum_inv_cov += fabs(inv_cov[j][k]);
        }
    }
    double x_bar = 0, sum_w_i = 0;
    for(int i = 0; i < dimension; i++){
        double sum_col_inv_cov = 0;
        for(int j = 0; j < dimension; j++){
            sum_col_inv_cov += fabs(inv_cov[i][j]);
        }
        double w_i = sum_col_inv_cov/sum_inv_cov;
        x_bar += w_i*x[i];
        sum_w_i += w_i;
    }
    if(fabs(sum_w_i-1)>0.00001)cout << "something wrong with weights in the correlated data combination" << endl;
    return x_bar;
}

vector<double> corr_weights(vector<double> delta_x, vector< vector<double> > correlations){
    int dimension = static_cast<int>(delta_x.size());
    TMatrixD covariance_matrix(dimension,dimension);
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            covariance_matrix[i][j] = correlations[i][j]*delta_x[i]*delta_x[j];
        }
    }
    TMatrixD inv_cov = covariance_matrix.Invert();
    double sum_inv_cov = 0;
    for(int j = 0; j < dimension; j++){
        for(int k = 0; k < dimension; k++){
            sum_inv_cov += fabs(inv_cov[j][k]);
        }
    }
    double sum_w_i = 0;
    vector<double> weights;
    for(int i = 0; i < dimension; i++){
        double sum_col_inv_cov = 0;
        for(int j = 0; j < dimension; j++){
            sum_col_inv_cov += fabs(inv_cov[i][j]);
        }
        double w_i = sum_col_inv_cov/sum_inv_cov;
        sum_w_i += w_i;
        weights.push_back(w_i);
    }
    if(fabs(sum_w_i-1)>0.00001)cout << "something wrong with weights in the correlated data combination" << endl;
    return weights;
}

double corr_comb_unc(vector<double> delta_x, vector< vector<double> > correlations){
    int dimension = static_cast<int>(delta_x.size());
    TMatrixD covariance_matrix(dimension,dimension);
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            covariance_matrix[i][j] = correlations[i][j]*delta_x[i]*delta_x[j];
        }
    }
    //need to make a copy of the original matrix before inverting, else it will be overwritten
    TMatrixD copy_cov(covariance_matrix);
    TMatrixD inv_cov = copy_cov.Invert();
    double sum_inv_cov = 0;
    for(int j = 0; j < dimension; j++){
        for(int k = 0; k < dimension; k++){
            sum_inv_cov += fabs(inv_cov[j][k]);
        }
    }
    vector<double> w;
    double sum_w_i = 0;
    for(int i = 0; i < dimension; i++){
        double sum_col_inv_cov = 0;
        for(int j = 0; j < dimension; j++){
            sum_col_inv_cov += fabs(inv_cov[i][j]);
        }
        double w_i = sum_col_inv_cov/sum_inv_cov;
        w.push_back(w_i);
        sum_w_i += w_i;
        //cout << w_i << endl;
    }
    if(fabs(sum_w_i-1)>0.00001)cout << "something wrong with weights in the correlated data combination" << endl;
    double delta_x_bar = 0;
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            delta_x_bar += w[i]*w[j]*covariance_matrix[i][j];
        }
    }
    return sqrt(delta_x_bar);
}

double corr_red_chi2(vector<double> x,vector<double> delta_x, vector< vector<double> > correlations){
    int dimension = static_cast<int>(delta_x.size());
    TMatrixD covariance_matrix(dimension,dimension);
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            covariance_matrix[i][j] = correlations[i][j]*delta_x[i]*delta_x[j];
        }
    }
    TMatrixD inv_cov = covariance_matrix.Invert();
    double x_bar = corr_comb(x,delta_x,correlations);
    double chi2 = 0;
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            chi2 += (x_bar-x[i])*inv_cov[i][j]*(x_bar-x[j]);
        }
    }
    return chi2/dimension-1;
}

/*inline double combine_asymmetry(double A1, double A2, double wA1, double wA2){
    return (wA1*A1+wA2*A2)/(wA1+wA2);
}*/
inline double combine_asymmetry_unc(double dA1, double dA2, double wA1, double wA2){
    return sqrt(pow(wA1/(wA1+wA2)*dA1,2)+pow(wA2/(wA1+wA2)*dA2,2));
}
inline double add_asymmetry_unc(double dA1, double dA2){
    return sqrt(pow(dA1,2)+pow(dA2,2));
}
