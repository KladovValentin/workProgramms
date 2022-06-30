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
#include <cmath>
#include "TLine.h"
#include "TChain.h"
#include "TString.h"
#include "TF1.h"
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;




float Goodphi[36] = { 0.1, 0.26, 0.34, 0.7,    0.76, 0.97, 1.03, 1.41,    1.46, 1.66, 1.73, 2.06,    2.18, 2.35, 2.43, 2.8,    2.85, 3.05, 3.13, 3.49,    3.55, 3.76, 3.83, 4.17,    4.26, 4.48, 4.55, 4.9,    4.98, 5.19, 5.26, 5.63,    5.67, 5.87, 5.95, 0.05 + 2.*TMath::Pi() };
float phi[40], ach[40], theta[40], z0[40], d0[40], energy[40], beam;
int cosm, nc, nn, charge[40], nch, schr[40];	
const char *cond0 = "(cosm == 0) && (nc == 2) && (nn==0)";
const char *cond1 = "(charge[0] == 1) && (charge[1] == 1) && (energy[0]/beam>=0.8) && (energy[1]/beam>=0.7) && (nch == 9)";
const char *cond2 = "(fabs(theta[0] + theta[1] - TMath::Pi()) <= 0.1) && (fabs(fabs(phi[0] - phi[1]) - TMath::Pi()) <= 0.07) && (sin(theta[0]) >=0.2)";
const char *cond3 = "(fabs(z0[0] - z0[1]) <= 1.0) && (fabs(d0[0] - d0[1]) <= 0.5)";




void go(){ 


	TChain chain("t1");
	chain.Add("/work/users/konctbel/calibs/R007-001/output/ntuples/col/exp2017_1000_col.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("energy", &energy);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("d0", &d0);

	chain.SetBranchAddress("cosm", &cosm);
	chain.SetBranchAddress("nc", &nc);
	chain.SetBranchAddress("nn", &nn);
	chain.SetBranchAddress("charge", &charge);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("beam", &beam);
	chain.SetBranchAddress("schr", &schr);
	
    const Long64_t entries = chain.GetEntries();
    //printf("%lld\n",entries);
	int counter = 0, scount = 0;
	TCanvas *c1 = new TCanvas("c1", "A", 200, 10, 500, 300);
	//c1->Divide(3, 3);  

	int j = 1;
	c1->cd(1);
	TProfile *hprof1 = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	/*j += 1;
	c1->cd(2);
	TProfile *hprof2 = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	j += 1;
	c1->cd(3);
	TProfile *hprof3 = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	j += 1;
	c1->cd(4);
	TProfile *hprof4 = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	j += 1;
	c1->cd(5);
	TProfile *hprof5 = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	j += 1;
	c1->cd(6);
	TProfile *hprof6 = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	j += 1;
	c1->cd(7);
	TProfile *hprof7 = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	j += 1;
	c1->cd(8);
	TProfile *hprof8 = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	j += 1;
	c1->cd(9);
	TProfile *hprof9 = new TProfile(Form("hprof%d", j), "Profile of ampl[j] versus phi[0]", 100, (float)(j - 2)*(2.*TMath::Pi()) / 9., (float)(j + 1)*(2.*TMath::Pi()) / 9., 0, 20);*/

	float ztr;
	for (int i = 1; i < entries; i++) {
		chain.GetEntry(i);

		if ((cosm == 0) && (nc == 2) && (nn == 0)) {
			//cout << beam << "	" << charge[0] << "	" << nch << "	" << energy[0] << "	" << fabs(theta[0] + theta[1] - TMath::Pi()) << endl;
			if ((charge[0] == 1) && (charge[1] == 1) && (energy[0] / beam >= 0.8) && (energy[1] / beam >= 0.7) /*&& (nch == 9)*/ && \
				(fabs(theta[0] + theta[1] - TMath::Pi()) <= 0.12) && (fabs(fabs(phi[0] - phi[1]) - TMath::Pi()) <= 0.1) && (sin(theta[0]) >= 0.2) && \
				(fabs(z0[0] - z0[1]) <= 2.0) && (fabs(d0[0] - d0[1]) <= 1.0)) {
				counter += 1;

				ztr = 13. / tan(theta[0]) + z0[0];
				if ((ztr > -11.) && (ztr < 11.) && (((phi[0] > 0.05) && (phi[0] < 0.28)) || ((phi[0] > 0.32) && (phi[0] < 0.75)))) {
					scount = 0;
					c1->cd(1);
					while (schr[1 + scount] != 0) { scount += 1; }
					hprof1->Fill(phi[0], ach[scount]); scount += 1;

					/*c1->cd(2);
					while (schr[1 + scount] != 0) { scount += 1; }
					hprof2->Fill(phi[0], ach[scount]); scount += 1;

					c1->cd(3);
					while (schr[1 + scount] != 0) { scount += 1; }
					hprof3->Fill(phi[0], ach[scount]); scount += 1;

					c1->cd(4);
					while (schr[1 + scount] != 0) { scount += 1; }
					hprof4->Fill(phi[0], ach[scount]); scount += 1;

					c1->cd(5);
					while (schr[1 + scount] != 0) { scount += 1; }
					hprof5->Fill(phi[0], ach[scount]); scount += 1;

					c1->cd(6);
					while (schr[1 + scount] != 0) { scount += 1; }
					hprof6->Fill(phi[0], ach[scount]); scount += 1;

					c1->cd(7);
					while (schr[1 + scount] != 0) { scount += 1; }
					hprof7->Fill(phi[0], ach[scount]); scount += 1;

					c1->cd(8);
					while (schr[1 + scount] != 0) { scount += 1; }
					hprof8->Fill(phi[0], ach[scount]); scount += 1;

					c1->cd(9);
					while (1 + scount < nch) {
						while (schr[1 + scount] != 0) { scount += 1; }
					}
					if ((phi[0]) < 1.) {
						hprof9->Fill(phi[0] + 2.*TMath::Pi(), ach[scount]);
					}
					else {
						hprof9->Fill(phi[0], ach[scount]);
					}*/
				}
			}
		}


	}
	//cout << counter << endl;

	c1->cd(1); hprof1->Draw();
	cout << hprof1->GetMean(2) << endl;
	/*	c1->cd(2); hprof2->Draw();
	c1->cd(3); hprof3->Draw();
	c1->cd(4); hprof4->Draw();
	c1->cd(5); hprof5->Draw();
	c1->cd(6); hprof6->Draw();
	c1->cd(7); hprof7->Draw();
	c1->cd(8); hprof8->Draw();
	c1->cd(9); hprof9->Draw();*/

				/*TLine *line1 = new TLine(Goodphi[4 * (j - 1)], 0, Goodphi[4 * (j - 1)], 18);
				TLine *line2 = new TLine(Goodphi[4 * (j - 1) + 1], 0, Goodphi[4 * (j - 1) + 1], 18);
				TLine *line3 = new TLine(Goodphi[4 * (j - 1) + 2], 0, Goodphi[4 * (j - 1) + 2], 18);
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
				line3->Draw();*/
				
				/*TProfile *hprof1 = new TProfile(Form("hprof%d",j),"Profile of ampl[j] versus ztr[0])",40,-20,20,0,20);
				float ztr = 0.;
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
				hprof1->Draw();
				TLine *line1 = new TLine(-11, 0,-11, 11);
				TLine *line2 = new TLine(11, 0, 11, 11);
				line1->SetLineColor(kRed);
				line2->SetLineColor(kRed);
				line1->Draw();
				line2->Draw();*/
				
    
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
    
	//pscp C:\disk\select_ach.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/select_ach.cpp
}











 
