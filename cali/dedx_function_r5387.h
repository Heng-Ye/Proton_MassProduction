#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <string>
#include <TImage.h>
#include <iomanip>
#include <TSpline.h>
#include <TText.h>
#include <TFrame.h>
#include <TH3.h>
#include <TH2.h>
#include "TMinuit.h"


double calib_factor =1.029e-3;
double Rho = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
double betap = 0.212;//(kV/cm)(g/cm^2)/MeV
double alpha = 0.93;//parameter from ArgoNeuT experiment at 0.481kV/cm 
double Wion = 23.6e-6;//parameter from ArgoNeuT experiment at 0.481kV/cm
double Efield = 0.4867;//kV/cm protoDUNE electric filed
double normalisation_factor=0.977;//for plane 2
const int Z=18; //Atomic number of Argon
const double A=39.948; // g/mol Atomic mass of Argon
const double I=188.0e-6; // ev
const double K=0.307; // Mev.cm^2 / mol
const double Mmu=105.658; // Mev for Mu
const double Me=0.51; // Mev for electron
const double rho=1.396;//g/cm^3


////getting the variable Efield using data driven maps
TFile *ef=new TFile("$DUNE_PARDATA_DIR/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");
TH3F *zpos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Pos");
TH3F *zneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Neg");
TH3F *zpos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Pos");
TH3F *zneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Neg");

TH3F *xpos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_X_Pos");
TH3F *xneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_X_Neg");
TH3F *xpos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_X_Pos");
TH3F *xneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_X_Neg");

TH3F *ypos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Y_Pos");
TH3F *yneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Y_Neg");
TH3F *ypos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Y_Pos");
TH3F *yneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Y_Neg");



double tot_Ef(double xval,double yval,double zval){
  double E0v=0.4867;
  if(xval>=0){
    double ex=E0v+E0v*xpos->Interpolate(xval,yval,zval);
    double ey=E0v*ypos->Interpolate(xval,yval,zval);
    double ez=E0v*zpos->Interpolate(xval,yval,zval);
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
if(xval<0){
  //else if {
  double ex=E0v+E0v*xneg->Interpolate(xval,yval,zval);
    double ey=E0v*yneg->Interpolate(xval,yval,zval);
    double ez=E0v*zneg->Interpolate(xval,yval,zval);
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
 else
   return E0v;
}

float zoffsetbd(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_bd->Interpolate(xval,yval,zval);
  }
  //if(xval<0){
  else {
    return zneg_bd->Interpolate(xval,yval,zval);
  }
}

float xoffsetbd(float xval,float yval,float zval){
  if(xval>=0){
    return xpos_bd->Interpolate(xval,yval,zval);
  }
  //if(xval<0){
  else {
    return xneg_bd->Interpolate(xval,yval,zval);
  }
}

float yoffsetbd(float xval,float yval,float zval){
  if(xval>=0){
    return ypos_bd->Interpolate(xval,yval,zval);
  }
  //if(xval<0){
  else {
    return yneg_bd->Interpolate(xval,yval,zval);
  }
}

/////////////////////////
//////////////////////Importing X fractional corrections//////////////////////
TFile my_file0("/dune/app/users/apaudel/code/prod4calib_factors/r5387/YZcalo_r5387.root");
TFile my_file2("/dune/app/users/apaudel/code/prod4calib_factors/r5387/Xcalo_r5387.root");
TH1F *X_correction_hist = (TH1F*)my_file2.Get("dqdx_X_correction_hist_2");
TH2F *YZ_correction_neg_hist_2=(TH2F*)my_file0.Get("correction_dqdx_ZvsY_negativeX_hist_2");
TH2F *YZ_correction_pos_hist_2=(TH2F*)my_file0.Get("correction_dqdx_ZvsY_positiveX_hist_2");
 ////////////////////////////////////////////////////////////////////////////////// 

Double_t dedx_function_r5387(double dqdx, double x, double y, double z){
  double Ef=tot_Ef(x,y,z);
  double Cx=X_correction_hist->GetBinContent(X_correction_hist->FindBin(x));
  double Cyz=0.0;
  if(x<0){
    Cyz=YZ_correction_neg_hist_2->GetBinContent(YZ_correction_neg_hist_2->FindBin(z,y));
  }
  //if(x>=0){
  else {
    Cyz=YZ_correction_pos_hist_2->GetBinContent(YZ_correction_pos_hist_2->FindBin(z,y));
  }
 //std::cout<<"(x,y,z):("<<x<<","<<y<<","<<","<<z<<")"<<std::endl; 
 double corrected_dqdx=dqdx*Cx*Cyz*normalisation_factor/calib_factor;
 return (exp(corrected_dqdx*(betap/(Rho*Ef)*Wion))-alpha)/(betap/(Rho*Ef));
}

Double_t corrected_dqdx_function_r5387(double dqdx, double x, double y, double z){
  double Ef=tot_Ef(x,y,z);
  double Cx=X_correction_hist->GetBinContent(X_correction_hist->FindBin(x));
  double Cyz=0.0;
  if(x<0){
    Cyz=YZ_correction_neg_hist_2->GetBinContent(YZ_correction_neg_hist_2->FindBin(z,y));
  }
  //if(x>=0){
  else {
    Cyz=YZ_correction_pos_hist_2->GetBinContent(YZ_correction_pos_hist_2->FindBin(z,y));
  }
 double corrected_dqdx=dqdx*Cx*Cyz*normalisation_factor/calib_factor;
 return corrected_dqdx;
}

////////////////////////////////////////////////////////////////////////////////
