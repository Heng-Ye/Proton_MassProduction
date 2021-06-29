#define ProtonMassProduction_run5387_cxx
#include "ProtonMassProduction_run5387.h"

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
#include <TGraph.h>
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

#include "./cali/dedx_function_r5387.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicAnaFunc.h"

#include <cassert>

using namespace std;
using namespace ROOT::Math;

void ProtonMassProduction_run5387::Loop() {

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	//histograms for mass production --------------------//
	//beam momentum (from beamline inst.)
	//start x/y/z position [after SCE correction]
	//end x/y/z position [after SCE correction]
	//proton track length
	//angle between beam track and TPC track (cosine_theta)

	//theta_x = acos(dcosx) etc
	//theta_y = acos(dcosx) etc
	//theta_z = acos(dcosx) etc

	int n_posz=600;
	double posz_min=-100.;
	double posz_max=200.;
	int n_posx=400;
	double posx_min=-200;
	double posx_max=200;
	int n_posy=1200;
	double posy_min=-500;
	double posy_max=700;

	TH1D *h1d_primtrk_stz_noSCECorr=new TH1D("h1d_primtrk_stz_noSCECorr","",2*n_posz,posz_min,posz_max); 
	h1d_primtrk_stz_noSCECorr->Sumw2();

	TH1D *h1d_primtrk_stx=new TH1D("h1d_primtrk_stx","",n_posx,posx_min,posx_max); h1d_primtrk_stx->SetName("h1d_primtrk_stx");
	TH1D *h1d_primtrk_sty=new TH1D("h1d_primtrk_sty","",n_posy,posy_min,posy_max); h1d_primtrk_sty->SetName("h1d_primtrk_sty");
	TH1D *h1d_primtrk_stz=new TH1D("h1d_primtrk_stz","",2*n_posz,posz_min,posz_max); h1d_primtrk_stz->SetName("h1d_primtrk_stz");
	h1d_primtrk_stx->Sumw2();
	h1d_primtrk_sty->Sumw2();
	h1d_primtrk_stz->Sumw2();

	TH1D *h1d_primtrk_endx=new TH1D("h1d_primtrk_endx","",n_posx,posx_min,posx_max); h1d_primtrk_endx->SetName("h1d_primtrk_endx");
	TH1D *h1d_primtrk_endy=new TH1D("h1d_primtrk_endy","",n_posy,posy_min,posy_max); h1d_primtrk_endy->SetName("h1d_primtrk_endy");
	TH1D *h1d_primtrk_endz=new TH1D("h1d_primtrk_endz","",n_posz,posz_min,posz_max); h1d_primtrk_endz->SetName("h1d_primtrk_endz");
	h1d_primtrk_endx->Sumw2();
	h1d_primtrk_endy->Sumw2();
	h1d_primtrk_endz->Sumw2();

	//track length
	int n_norm_primtrklen=200;
	double norm_primtrklen_min=0;
	double norm_primtrllen_max=2;
	TH1D *trklen=new TH1D("trklen","",150,0,150);
	TH1D *h1d_norm_primtrklen=new TH1D("h1d_norm_primtrklen","",n_norm_primtrklen,norm_primtrklen_min,norm_primtrllen_max);
	trklen->Sumw2();
	h1d_norm_primtrklen->Sumw2();

	//proton kinetic energies
	int n_ke=40;
	double ke_min=0.;
	double ke_max=800;
	TH1D *h1d_ke_beam=new TH1D("h1d_ke_beam","",n_ke, ke_min, ke_max);
	TH1D *h1d_ke_trklen=new TH1D("h1d_ke_trklen","",n_ke, ke_min, ke_max);
	h1d_ke_beam->GetXaxis()->SetTitle("MeV");

	//Cosine_beamdirection_primtrkdirection
	int n_cos=300;
	double cos_min=-1.5;
	double cos_max=1.5;
	TH1D *h1d_cosine_beam_primtrk=new TH1D("h1d_cosine_beam_primtrk","",n_cos,cos_min,cos_max);

	//xy-distribution
	TH2D *h2d_x_y_1sthit_before_sce_xycut=new TH2D("h2d_x_y_1sthit_before_sce_xycut","", 70,-60,10,60,390,450);

	//track direction
        TH1D *hreco_beam_angleX = new TH1D(Form("hreco_beam_angleX"), "", 180, 0, 180); //theta_x: angle betweeb track direction and x-axis, unit in degree
      	hreco_beam_angleX->Sumw2();



	//start entering event loop
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //evt loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

		size_t j=0;
		//if (Cut(jentry) < 0) continue;

		//get parameters here
		double primtrkx=-99;
		double primtrky=-99;
		double primtrkz=-99;
		double primtrkendx=-99;
		double primtrkendy=-99;
		double primtrkendz=-99;
		double primtrk_len=-99;

		double tmp_cosine_beam_primtrk=-99; 

		double mom_beam=-99; 

		//std::cout<<"ck 0"<<std::endl;
		//check if the vector is empty, if empty, skip the event
		bool IsbeamPosEmpty=beamPosz->empty();
		bool IsBeamEmpty=beamMomentum->empty();
		bool IsPrimtrk_hitEmpty=primtrk_hitz->empty();
		bool IsPrimtrk_start=primtrk_startz->empty();
		bool IsPrimtrk_dqdxEmpty=primtrk_dqdx->empty();
		bool IsPrimtrklen=primtrklen->empty();
		bool IsPrimtrk_range=primtrk_range->empty();
		bool IsPrimtrk_resrange=primtrk_resrange->empty();
		bool IsPrimtrk_pitch=primtrk_pitch->empty();
		bool IsEmpty=false;
		bool IsRecoBeamDirectionFlip=false;
		if (IsPrimtrk_start||IsPrimtrk_hitEmpty||IsBeamEmpty||IsbeamPosEmpty||IsPrimtrk_dqdxEmpty||IsPrimtrklen||IsPrimtrk_range||IsPrimtrk_resrange||IsPrimtrk_pitch) { //if any of these containers are empty
			IsEmpty=true;
			continue;
		}
		int sz=primtrk_hitz->at(j).size();
		if (sz==0) { //non-zero container check
			continue;
		} //non-zero container check


		//check if pandora track flipped and access start/end x/y/z
		if (primtrk_hitz->at(j)[0]>primtrk_hitz->at(j)[primtrk_hitz->at(j).size()-1]) { //flip the start/end point if the direction flipped
			primtrkx=primtrk_hitx->at(j)[primtrk_hitx->at(j).size()-1];
			primtrky=primtrk_hity->at(j)[primtrk_hity->at(j).size()-1];
			primtrkz=primtrk_hitz->at(j)[primtrk_hitz->at(j).size()-1];
			primtrkendx=primtrk_hitx->at(j)[0];
			primtrkendy=primtrk_hity->at(j)[0];
			primtrkendz=primtrk_hitz->at(j)[0];
		} //flip the start/end point if the direction flipped
		else {
			primtrkendx=primtrk_hitx->at(j)[primtrk_hitx->at(j).size()-1];
			primtrkendy=primtrk_hity->at(j)[primtrk_hity->at(j).size()-1];
			primtrkendz=primtrk_hitz->at(j)[primtrk_hitz->at(j).size()-1];
			primtrkx=primtrk_hitx->at(j)[0];
			primtrky=primtrk_hity->at(j)[0];
			primtrkz=primtrk_hitz->at(j)[0];
		}
		


		//std::cout<<"ck 2"<<std::endl;
		primtrk_len=primtrk_range->at(j);
		mom_beam=beamMomentum->at(j); //beam momentum [in unit of GeV/c]
		double ke_beam=p2ke(mom_beam); //GeV
		double ke_beam_MeV=1000.*ke_beam; //MeV
		double ke_trklen=1000.*ke_vs_csda_range_sm->Eval(primtrk_len); //MeV		

		double csda_val=csda_range_vs_mom_sm->Eval(mom_beam); //the expected value for the stopping protons
		double norm_trklen_csda=primtrk_len/csda_val; //normalized track length

		tmp_cosine_beam_primtrk=cosine_beam_primtrk; //cosine between beam and tpc track
		if (tmp_cosine_beam_primtrk<0) { //flip sign to fix the numerical bug 
			tmp_cosine_beam_primtrk=-1.*tmp_cosine_beam_primtrk;
		} //flip sign to fix the numerical bug


		size_t k=0; //old data format
		//calculate track direction (track end - start)
    		if (!primtrk_dqdx->at(k).empty()){    
      			TVector3 pt0(primtrkx, primtrky, primtrkz);
      			TVector3 pt1(primtrkendx, primtrkendy, primtrkendz);
      			TVector3 dir = pt1 - pt0;
      			dir = dir.Unit();
			hreco_beam_angleX->Fill(acos(dir.X())*180/TMath::Pi()); //theta_x
			//theta_y: dir.Y()

		}

		//xy distribution
		double x0_tmp=0;
		double y0_tmp=0;
		//double z0_tmp=0;
		if ((primtrk_startz->at(-1+primtrk_startz->size()))>(primtrk_startz->at(0))) { //check if Pandora flip the sign
			x0_tmp=primtrk_startx->at(0);
			y0_tmp=primtrk_starty->at(0);
			//z0_tmp=primtrk_startz->at(0);
		} //check if Pandora flip the sign
		else {
			x0_tmp=primtrk_startx->at(-1+primtrk_startx->size());
			y0_tmp=primtrk_starty->at(-1+primtrk_starty->size());
			//z0_tmp=primtrk_starty->at(-1+primtrk_startz->size());
		}


		//fill histograms -----------------------------------//
		h1d_primtrk_stx->Fill(primtrkx);
		h1d_primtrk_sty->Fill(primtrky);
		h1d_primtrk_stz->Fill(primtrkz);

		h1d_primtrk_endx->Fill(primtrkendx);
		h1d_primtrk_endy->Fill(primtrkendy);
		h1d_primtrk_endz->Fill(primtrkendz);

		h1d_primtrk_stz_noSCECorr->Fill(primtrk_startz->at(0));
		//y:primtrk_starty

		trklen->Fill(primtrk_len);
		h1d_norm_primtrklen->Fill(norm_trklen_csda);

		h1d_cosine_beam_primtrk->Fill(tmp_cosine_beam_primtrk);
		h1d_ke_beam->Fill(ke_beam_MeV);
		h1d_ke_trklen->Fill(ke_trklen);
		h2d_x_y_1sthit_before_sce_xycut->Fill(x0_tmp,y0_tmp);

	} //evt loop


	//Time stamp and run number of this data run	
   	TParameter<Int_t>* run_num=new TParameter<Int_t>("run_num",0.);
   	TParameter<Double_t>* evt_time=new TParameter<Double_t>("evt_time",0.);
   	run_num->SetVal(run); //run #
   	evt_time->SetVal(evttime); //event start time


	TFile *fout = new TFile(Form("massprod_run%d.root",run),"RECREATE");
     	run_num->Write();
     	evt_time->Write();	

	h1d_primtrk_stz_noSCECorr->Write();

	h1d_primtrk_stx->Write();
	h1d_primtrk_sty->Write();
	h1d_primtrk_stz->Write();

	h1d_primtrk_endx->Write();
	h1d_primtrk_endy->Write();
	h1d_primtrk_endz->Write();

	trklen->Write();
	h1d_norm_primtrklen->Write(); 
	h1d_cosine_beam_primtrk->Write();

	h1d_ke_beam->Write();
	h1d_ke_trklen->Write();
	
	h2d_x_y_1sthit_before_sce_xycut->Write();

	hreco_beam_angleX->Write();

	fout->Close();

}
