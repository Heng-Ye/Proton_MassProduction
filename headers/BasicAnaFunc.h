/*
#include <stdexcept>      // std::out_of_range
#include <vector>
#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TProfile2D.h>
//#include <iostream>
#include <fstream>
#include <string>
#include "TCanvas.h" 
#include "TGraph.h" 
#include "TGraphSmooth.h" 
#include "TVectorD.h"
#include "TParameter.h"
#include "TGraphErrors.h"
#include <TProfile3D.h>
#include <TProfile2D.h>
#include "TVector3.h"

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
#include <cassert>
*/

//using namespace std;
//using namespace ROOT::Math;


//===================================================================================================//


double p2ke(double p) {
	double ke=m_proton*(-1.+sqrt(1.+(p*p)/(m_proton*m_proton)));
	return ke;
}

double ke2p(double ke) { //input ke unit: GeV
	double p=m_proton*sqrt(-1+pow((1+ke/m_proton),2));
	return p;
}

Double_t fitg(Double_t* x,Double_t *par) {
	double m=par[0];
	double s=par[1];
	double n=par[2];

	Double_t g=n*TMath::Exp(-(x[0]-m)*(x[0]-m)/(2*s*s));
	//Double_t g=n/(s*sqrt(2*3.14159))*TMath::Exp(-(x[0]-m)*(x[0]-m)/(2*s*s));
	return (g);
}

TF1* VFit(TH1D* h, Int_t col) {
	//pre-fit parameters
	float pre_mean=h->GetBinCenter(h->GetMaximumBin());
	float pre_max=h->GetMaximum();
	float pre_rms=h->GetRMS();
	cout<<"mean: "<<pre_mean<<endl;
	cout<<"rms: "<<pre_rms<<endl;
	cout<<"max: "<<pre_max<<endl;
	cout<<""<<endl;

	//1st fitting
	TF1 *gg=new TF1("gg",fitg,pre_mean-1.*pre_rms,pre_mean+1.*pre_rms,3);
	gg->SetParameter(0,pre_mean);
	gg->SetParameter(1,pre_rms);
	gg->SetParameter(2,pre_max);
	//if (pre_rms>1.0e+06) { gg->SetParLimits(1,0,100); }

	gg->SetLineColor(col);
	gg->SetLineStyle(2);
	h->Fit("gg","remn");

	//2nd fitting
	TF1 *g=new TF1("g",fitg,gg->GetParameter(0)-3.*gg->GetParameter(1),gg->GetParameter(0)+3.*gg->GetParameter(1),3);
	//TF1 *g=new TF1("g",fitg,gg->GetParameter(0)-1,gg->GetParameter(0)+.5,3);
	g->SetParameter(0,gg->GetParameter(0));
	g->SetParameter(1,gg->GetParameter(1));
	g->SetParameter(2,gg->GetParameter(2));

	g->SetParLimits(0,gg->GetParameter(0)-3*gg->GetParameter(1), gg->GetParameter(0)+3*gg->GetParameter(1));
	double sss=gg->GetParameter(1); if (sss<0) sss=-sss;
	g->SetParLimits(1,0,5.*sss);
	g->SetParLimits(2,0,10.*sqrt(pre_max));

	g->SetLineColor(col);
	g->SetLineStyle(2);
	g->SetLineWidth(2);

	h->Fit("g","rem+");
	return g;
}

Double_t govg(Double_t* x,Double_t *par) {
	//g1
	double m1=par[0];
	double s1=par[1];
	//double n1=1.;
	double g1=-(x[0]-m1)*(x[0]-m1)/(2*s1*s1);

	//g2
	double m2=par[2];
	double s2=par[3];
	//double n2=1.;
	double g2=-(x[0]-m2)*(x[0]-m2)/(2*s2*s2);

	//g2/g1
	double g_ov_g=0; 
	g_ov_g=TMath::Exp(g2-g1);
	if (m1==m2&&s1==s2) g_ov_g=1;

	return g_ov_g;
}

//Read file of csda range versus kinetic energy---------------------------------//
TFile *fke_csda=new TFile("../proton_ke_csda_converter.root");
TGraph *csda_range_vs_ke_sm=(TGraph *)fke_csda->Get("csda_range_vs_ke_sm");
TGraph *ke_vs_csda_range_sm=(TGraph *)fke_csda->Get("ke_vs_csda_range_sm");

//Read file of csda range versus momentum -----------------------------------------------------//
TString conv_path="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/conversion/";
TFile *fmom_csda=new TFile(Form("%s/proton_mom_csda_converter.root",conv_path.Data()));
TGraph *csda_range_vs_mom_sm=(TGraph *)fmom_csda->Get("csda_range_vs_mom_sm");
TGraph *mom_vs_csda_range_sm=(TGraph *)fmom_csda->Get("mom_vs_csda_range_sm");


