#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLine.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;



void go(){ 

	
    float phi[2], ach[9], theta[2], z0[2];
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true_combined.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);	

    const Long64_t entries = chain.GetEntries();
    //printf("%lld\n",entries);

	float Goodphi[36] = { 0.1, 0.26, 0.34, 0.7,    0.76, 0.97, 1.03, 1.41,    1.46, 1.66, 1.73, 2.06,    2.18, 2.35, 2.43, 2.8,    2.85, 3.05, 3.13, 3.49,    3.55, 3.76, 3.83, 4.17,    \
		4.26, 4.48, 4.55, 4.9,    4.98, 5.19, 5.26, 5.63,    5.67, 5.87, 5.95, 0.05 };	

    TCanvas *c1 = new TCanvas("c1","A",200,10,500,300); 
    c1->Divide(3,3);
	for (int j = 1; j < 10; j++) {
		c1->cd(j);
		TProfile *hprof = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);
		float ztr = 0.;
		for (int i = 0; i < entries; i++) {
			chain.GetEntry(i);
			ztr = 13. / tan(theta[0]) + z0[0];
			if ((ztr > -11.) && (ztr < 11.)) {
				if ((j == 9) && (phi[0]) < 1.) {
					hprof->Fill(phi[0]+ 2.*TMath::Pi(), ach[j-1]);
				}
				else {
					hprof->Fill(phi[0], ach[j-1]);
				}
			}
		}
		hprof->Draw();
		TLine *line1 = new TLine(Goodphi[4 * (j - 1)  ], 0, Goodphi[4 * (j - 1)  ], 18);
		TLine *line2 = new TLine(Goodphi[4 * (j - 1)+1], 0, Goodphi[4 * (j - 1)+1], 18);
		TLine *line3 = new TLine(Goodphi[4 * (j - 1)+2], 0, Goodphi[4 * (j - 1)+2], 18);
		if (j == 9) {
			TLine *line4 = new TLine(Goodphi[4 * (j - 1) + 3] + 2.*TMath::Pi(), 0, Goodphi[4 * (j - 1) + 3] + 2.*TMath::Pi(), 18);
			line4->SetLineColor(kRed);
			line4->Draw();
		}
		else {
			TLine *line4 = new TLine(Goodphi[4 * (j - 1) + 3], 0, Goodphi[4 * (j - 1) + 3], 18);
			line4->SetLineColor(kRed);
			line4->Draw();
		}
		line1->SetLineColor(kRed);
		line2->SetLineColor(kRed);
		line3->SetLineColor(kRed);
		line1->Draw();
		line2->Draw();
		line3->Draw();
		


		/*TProfile *hprof1 = new TProfile(Form("hprof%d",j),"Profile of ampl[j] versus ztr[0])",40,-20,20,0,20);
		float ztr = 0.;
		for (int i=0; i<entries; i++){
			chain.GetEntry(i);
			ztr = 13. / tan(theta[0]) + z0[0];
			if(j == 9){
				if (((phi[0] > Goodphi[4 * (j-1)]) && (phi[0] < Goodphi[4 * (j-1) + 1])) || ((phi[0] > Goodphi[4 * (j-1) + 2]) || (phi[0] < Goodphi[4 * (j-1) + 3]))) {
					hprof1->Fill(ztr, ach[j-1]);
				}
			}
			else{
				if (((phi[0] > Goodphi[4 * (j-1)]) && (phi[0] < Goodphi[4 * (j-1) + 1])) || ((phi[0] > Goodphi[4 * (j-1) + 2]) && (phi[0] < Goodphi[4 * (j-1) + 3]))){
					hprof1->Fill(ztr, ach[j-1]);
				}
			}
		}
		hprof1->Draw();
		TLine *line1 = new TLine(-11, 0,-11, 11);
		TLine *line2 = new TLine(11, 0, 11, 11);
		line1->SetLineColor(kRed);
		line2->SetLineColor(kRed);
		line1->Draw();
		line2->Draw();*/
	}
    
//_____________________________________________________________________________дальше двумерное__________________________________________________________________________________________\\
	
	
	
	
	
	/*TCanvas *c2 = new TCanvas("c2","B",200,10,500,300);
	TProfile2D *hprofile = new TProfile2D("hprof2d","Profile of pz versus px and py",100,0.0,6.3,40,-13,13,0,20);
        //TH2* h2 = new TH2F("h2", "h2 title", 100, 0.0, 6.3, 40, -13, 13);
	//TH2* h3 = new TH2F("h3", "h3 title", 100, 0.0, 6.3, 40, -13, 13);
	float ztr = 0.;
	for (int i = 0; i < entries; i++) {
		chain.GetEntry(i);
		const int j = (const int)(fabs(phi[0] - 0.05) * 9. / (2.*TMath::Pi()));
		ztr = 13. / tan(theta[0]) + z0[0];
		if (j == 8) {
			if ((ztr > -11.) && (ztr < 11.) && (((phi[0] > Goodphi[4 * (j)]) && (phi[0] < Goodphi[4 * (j) + 1])) || ((phi[0] > Goodphi[4 * (j) + 2]) || (phi[0] < Goodphi[4 * (j) + 3])))) {
				hprofile->Fill(phi[0],ztr,ach[j],1);
			}
		}
		else {
			if ((ztr > -11.) && (ztr < 11.) && (((phi[0] > Goodphi[4 * (j)]) && (phi[0] < Goodphi[4 * (j) + 1])) || ((phi[0] > Goodphi[4 * (j) + 2]) && (phi[0] < Goodphi[4 * (j) + 3])))) {
				hprofile->Fill(phi[0],ztr,ach[j],1);
			}
		}
	}
	hprofile->Draw("LEGO2");*/
    

}











 
