#include <fstream> 
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TMath.h"
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TSpline.h"
#include "TTree.h"
#include "TH1D.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TCollection.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TLegend.h"
#include "TKey.h"
#include "TGraph.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TLine.h"
#include "TString.h"
#include "TLatex.h"
#include "TAxis.h"

#include "/work/users/kladov/snd2k/R006-003/maindir/a1a2a3a4.cpp"
#include "/work/users/kladov/snd2k/R006-003/maindir/gran.cpp"
//#include "/work/users/kladov/snd2k/R006-003/maindir/gran1.cpp"
#include "/work/users/kladov/snd2k/R006-003/model/gran1.cpp"
#include "/work/users/kladov/snd2k/R006-003/model/truegran1m.cpp"
#include "/work/users/kladov/snd2k/R006-003/model/counteri.cpp"
#include "/work/users/kladov/snd2k/R006-003/maindir/goodruns.cpp"
#include "/work/users/kladov/snd2k/R006-003/maindir/truegran.cpp"
#include "/work/users/kladov/snd2k/R006-003/maindir/truegran1.cpp"
#include "/work/users/kladov/snd2k/R006-003/maindir/yfm.cpp"
#include "/work/users/kladov/snd2k/R006-003/maindir/shiftersfm.cpp"
using namespace std;


float Goodphi[36] = { 0.08, 0.27, 0.33, 0.71,    0.78, 0.98, 1.04, 1.42,    1.48, 1.68, 1.75, 2.12,    2.21, 2.39, 2.45, 2.85,    2.92, 3.11, 3.17, 3.53,    3.61, 3.78, 3.85, 4.23,    \
	4.32, 4.5, 4.56, 4.92,    4.99, 5.17, 5.24, 5.61,    5.67, 5.87, 5.93, 6.28 };
//double a1[9][14], a2[9][14], a3[9][14], a4[9][14];
//float gran[4][9][14];
//double gran1[4][9][14];
//double pik[9][14];
//double truegran[9][4];
//double truegran1m[9][4];
//double yfm[9][14][12];
float zobl[15] = { -11.6013, -10.0367, -8.47203, -6.90737, -5.34273, -3.77807, -2.21343, -0.648766, 0.9159, 2.48053, 4.0452, 5.60983, 7.1745, 8.73913, 9.6};
float egran[5] = { 600., 780., 900., 955., 1010. };
float phi[2], phis[2], ach[40], theta[2], z0[2], schr[40], tch[40], tchr[40], energy[40], beam, eton;
int region[40];

const char* cond = "(cosm == 0) && (nc == 2) && (nn==0) && (charge[0] == 1) && (charge[1] == 1) && (energy[0]/beam>=0.8) && (energy[1]/beam>=0.7) && (nch >= 9) && \
((theta[0] + theta[1] - 3.14159) <= 0.12) && ((theta[0] + theta[1] - 3.14159) >= -0.12) && (sin(theta[0]) >= 0.41) && \
((((phi[0] - phi[1]) - 3.14159 <= 0.1) && (((phi[0] - phi[1]) - 3.14159) >= -0.1)) || (((phi[0] - phi[1]) + 3.14159 <= 0.1) && (((phi[0] - phi[1]) + 3.14159) >= -0.1))) && \
((z0[0] - z0[1]) <= 1.0) && ((z0[0] - z0[1]) >= -1.0) && ((d0[0] - d0[1]) <= 0.5) &&  ((d0[0] - d0[1]) >= -0.5)";
int mybin = 100;
int klmno[5], nch, run;
double bound[5];
double prob = 0.;
double probmax = 0.;
double eff;
float pedestali1[25][100];
float ach1[40];
float ped[14][9];
float pedz[9];
float ped1[9][14][9];
float ped2[9];
double pede[14][9][4];
double meanample[14];
double meanamplm[14];
double meanampleerr[14];
double meanamplmerr[14];
int wichcounter1[10];
/*double countc1[14];
double countc2[14];
double counta1[9][14];
double counta2[9][14];
double countz1[9];
double countz2[9];
double ecountc1[14];
double ecountc2[14];
double ecounta1[9][14];
double ecounta2[9][14];
double ecountz1[9];
double ecountz2[9];*/
float PI = 3.14159;

TString str = "a";
TString elem;

double p3g0(double x, double xr, double sr, std::vector<double> p);
double aerogel(size_t n, double x, double* par);
double shifter(double x, double* par);
double p3g(double x, double xr, double sr, std::vector<double> p, double smin, double smax);
double scphi(double* x, double* par);


void copy1() {
	vector<TString> ans;
	
	TString basedir("/work/users/kladov/snd2k/R007-001/second/");
	TString dir("/work/users/kladov/snd2k/R006-003/model/selected/");
	TString files = gSystem->GetFromPipe("ls " + basedir + "*.root");

	Ssiz_t from = 0;
	
	while (files.Tokenize(elem, from, "\n")) {
		ans.push_back(elem);
	}
	unsigned int vector_size = ans.size();
	for (int i = 0; i < vector_size; i++) {
		str = ans[i].Copy();
		TString tok;
		Ssiz_t from1 = 1;
		vector<TString> tokens;
		cout << str << endl;

		TFile* filecomb = TFile::Open(str);
		TTree* originalTree = (TTree*)filecomb->Get("h1");

		while (str.Tokenize(tok, from1, "/")) {
			tokens.push_back(tok);
		}

		cout << tokens[6] << endl;
		cout << dir + "true3_" + tokens[6] << endl;
		TFile* ouput = TFile::Open(dir + "true3_" + tokens[6], "RECREATE");
		ouput->cd();
		TTree* selectedTree = originalTree->CopyTree(cond);
		selectedTree->Print();
		selectedTree->Write();
		filecomb->Close();
		ouput->Close();
	}
}

void timespectra() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true_exp**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("run", &run);

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;


	TH1* histo = new TH1F("histo", "time spectre;time;", 1500, 0, 300.);
	

	int Count = 0, count1 = 0;
	float ztr = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.25 / tan(theta[0]) + z0[0];
		if ((ztr > zobl[0]) && (ztr < zobl[14]) && (run > 27225) && (goodrnumb[run - 25001] == 1) && (goodreff[run - 25001] == 1)) {
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9; j++) {
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					while (1 + scount < nch) { scount += 1; }
				}
				else
					while (schr[1 + scount] > 0.5) { scount += 1; }

				for (int a = scount; a >= scount; a--) {
					if ((schr[scount] > 0.5))
						histo->Fill(tch[a]);
				}
				scount += 1;
				scount1 = scount;
			}
		}

		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			Count = 0;
		}
	}
	histo->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	TLine* line11 = new TLine(110, 0, 110, 1.1 * histo->GetMaximum());
	line11->SetLineColor(kGreen);
	line11->Draw();
	c->Update();

	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/timespectre.root", "RECREATE");
	histo->Write("Tspectre");
	//MyFile->Write();
	MyFile->Close();

	
	cin.get();
}

void draweffspectre() {
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/runsi.root");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0;
	while ((key = (TKey*)next()) && counterr< 1) {
		TString name(key->GetName());
		TString type = "sensor";
		if (name.Contains("spec")) {
			counterr +=1;
			//cin.get();
			int counter;
			sscanf(name.Data(), "sens%d_eff_spec", &counter);

			TH1F* h_amp_phi = (TH1F*)f->Get(Form("sens%d_eff_spec", counter));
			h_amp_phi->SetTitle(Form("sens%d_eff_spec;eff", counter));
			//h_amp_phi->Rebin(10);
			h_amp_phi->Draw();
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->SetFillColor(0);
			c->SetFrameFillColor(0);
			c->Update();
			TLine* line11 = new TLine(0.8, 0, 0.8, 1.1 * h_amp_phi->GetMaximum());
			line11->SetLineColor(kGreen);
			line11->Draw();
			c->Update();
			cin.get();
		}
	}
}

void raspramp()	 {
	TChain chain("h1");
	chain.Add("/work/users/kladov/snd2k/R006-003/model/selected/true3_ee**.root");
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("schr", &schr);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("tchr", &tchr);
	chain.SetBranchAddress("region", &region);

	int scount = 0;
	int scount1 = 0;

	//hists for ach or eff
	TH1* ha[9][14][9];
	TH1* hn[9][14][9];
	char namea[40];
	char namen[40];
	char title[200];
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				sprintf(namea, "hprofa%d", 14 * 9 * j + 9 * k + i);
				sprintf(namen, "hprofn%d", 14 * 9 * j + 9 * k + i);
				sprintf(title, "sens%d,z%d,phi%d;eff;run", i + 1, k + 1, j + 1);
				ha[j][k][i] = new TH1F(namea, title, 550, -5, 50);
				hn[j][k][i] = new TH1F(namen, title, 100, -0.1, 2.);
			}
		}
	}



	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	float ztr1 = 0.;
	double lgran[4];
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.25 / tan(theta[0]) + z0[0];
		ztr1 = 12.25 / tan(theta[1]) + z0[1];
		schr[nch] = 0.;
		if (1) {
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9; j++) {
				//wichcounter1[j] = scount1;
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					while (1 + scount < nch) { scount += 1; }
				}
				else
					while (schr[1 + scount] > 0.5) { scount += 1; }
				lgran[0] = truegran1m[j][0];// +0.03;
				lgran[1] = truegran1m[j][1];// -0.045;
				lgran[2] = truegran1m[j][1];// +0.045;
				lgran[3] = truegran1m[j][3];// -0.03;

				//ach izm i eff izm 9X14X9
				for (int k = 0; k < 14; k++) {
					if ((ztr > zobl[k]) && (ztr < zobl[k + 1]) && (region[0] == 1)) {
						int it = 9;
						for (int i = 0; i < 3; i++) {
							if ((phi[0] >= (lgran[0] + i * (lgran[1] - lgran[0]) / 3.)) && (phi[0] < (lgran[0] + (i + 1) * (lgran[1] - lgran[0]) / 3.)))
								it = i;
							if ((phi[0] >= (lgran[2] + i * (lgran[3] - lgran[2]) / 6.)) && (phi[0] < (lgran[2] + (i + 1) * (lgran[3] - lgran[2]) / 6.)))
								it = i + 3;
							if ((phi[0] >= (lgran[2] + (i + 3) * (lgran[3] - lgran[2]) / 6.)) && (phi[0] < (lgran[2] + (i + 4) * (lgran[3] - lgran[2]) / 6.)))
								it = i + 6;
						}
						if (it < 9) {
							if (ach[j] > 0.2)
								hn[0][k][j]->Fill(1);
							else 
								hn[0][k][j]->Fill(0);
							ha[0][k][j]->Fill(ach[j]);
						}
					}
					if ((ztr1 > zobl[k]) && (ztr1 < zobl[k + 1]) && (region[1] == 1)) {
						int it = 9;
						for (int i = 0; i < 3; i++) {
							if ((phi[1] >= (lgran[0] + i * (lgran[1] - lgran[0]) / 3.)) && (phi[1] < (lgran[0] + (i + 1) * (lgran[1] - lgran[0]) / 3.)))
								it = i;
							if ((phi[1] >= (lgran[2] + i * (lgran[3] - lgran[2]) / 6.)) && (phi[1] < (lgran[2] + (i + 1) * (lgran[3] - lgran[2]) / 6.)))
								it = i + 3;
							if ((phi[1] >= (lgran[2] + (i + 3) * (lgran[3] - lgran[2]) / 6.)) && (phi[1] < (lgran[2] + (i + 4) * (lgran[3] - lgran[2]) / 6.)))
								it = i + 6;
						}
						if (it < 9) {
							if (it < 9) {
								if (ach[j] > 0.2)
									hn[0][k][j]->Fill(1);
								else
									hn[0][k][j]->Fill(0);
								ha[0][k][j]->Fill(ach[j]);
							}
						}
					}
				}
				scount += 1;
				scount1 = scount;
			}
		}

		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}

	vector<float> errx[9], effect[9], erry[9], ampl[9];
	//effotach

	char titlea[200];
	char titlen[200];
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/model/random.root", "RECREATE");
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				sprintf(titlea, "ach_sens%d,z%d,phi%d", i + 1, k + 1, j + 1);
				sprintf(titlen, "eff_sens%d,z%d,phi%d", i + 1, k + 1, j + 1);
				//cout << h[j][k][i]->GetMean(1) << endl;
				hn[j][k][i]->Write(titlen);
				ha[j][k][i]->Write(titlea);
			}
		}
	}
	MyFile->Close();

	TFile* ff = new TFile("/work/users/kladov/snd2k/R006-003/model/random.root");
	int schit = 0;
	for (int j = 0; j < 1; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				TH1F* ha = (TH1F*)ff->Get(Form("ach_sens%d,z%d,phi%d", i + 1, k + 1, j + 1));
				TH1F* hn = (TH1F*)ff->Get(Form("eff_sens%d,z%d,phi%d", i + 1, k + 1, j + 1));
				if (hn->GetMean(1) > 0.2) {
					schit = schit + 1;
					ampl[0].push_back(ha->GetMean(1));
					effect[0].push_back(hn->GetMean(1));
					errx[0].push_back(1. * ha->GetMeanError(1));
					erry[0].push_back(1. * hn->GetMeanError(1));
				}
			}
		}
	}
	cout << "col_tochek" << schit << endl;
	ff->Close();
	int numbrr = 0;
	double parametr1[9];
	double counter1[9];
	double parametr1Err[9];
	double nolliki[9];
	double chisquare1[9];
	double chisquare2[9];
	double degrof1[9];
	double degrof2[9];
	double erre[9];
	TFile* MyFile1 = new TFile("/work/users/kladov/snd2k/R006-003/model/effotachfitted.root", "RECREATE");
	for (int i = 0; i < 1; i++) {
		counter1[i] = (double)i + 1;
		nolliki[i] = 0.;
		numbrr = ampl[i].size();
		TGraphErrors* gr = new TGraphErrors(numbrr, &ampl[i][0], &effect[i][0], &errx[i][0], &erry[i][0]);
		gr->SetTitle(";ach,pe;eff");
		gr->Draw("A*");
		TF1* func = new TF1("func", "[1]-exp(-x/[0])", 2., gr->GetXaxis()->GetXmax());
		TF1* funcx = new TF1("funcx", "1-exp(-x/[0])", 2., gr->GetXaxis()->GetXmax());
		func->SetLineColor(kRed);
		funcx->SetLineColor(kGreen);
		func->SetParameter(0, 1.);
		func->SetParameter(1, 1.);
		func->SetParLimits(0, 0.5, 10.);
		func->SetParLimits(1, 0.5, 1.5);
		funcx->SetParameter(0, 1.);
		funcx->SetParLimits(0, 0.5, 10.);
		gStyle->SetOptFit(1111);

		gr->Fit(funcx, "", "", 2., gr->GetXaxis()->GetXmax());
		TF1* f1 = gr->GetFunction("funcx");
		cout << Form("1-counter%d", i + 1) << endl << "Chisquare= " << f1->GetChisquare() << "	" << "NDF= " << f1->GetNDF() << endl << endl << endl;
		chisquare2[i] = f1->GetChisquare();
		degrof2[i] = 2 * f1->GetNDF();

		gr->Fit(func, "", "", 2., gr->GetXaxis()->GetXmax());
		TF1* f = gr->GetFunction("func");
		parametr1[i] = func->GetParameter(1);
		erre[i] = func->GetParError(1);
		parametr1Err[i] = func->GetParError(1);
		cout << parametr1[i] << endl << endl;
		cout << Form("[1]-counter%d", i + 1) << endl << "Chisquare= " << f->GetChisquare() << "	" << "NDF= " << f->GetNDF() << endl;
		chisquare1[i] = f->GetChisquare();
		degrof1[i] = 2 * f->GetNDF();


		TF1* func1 = new TF1("func1", "[1]-exp(-x/[0])", 0, gr->GetXaxis()->GetXmax());
		//Double_t par[2];
		//func->GetParameters(par);
		func1->SetParameters(func->GetParameters());
		func1->SetLineColor(kRed);
		func1->Draw("same");

		TF1* funcx1 = new TF1("funcx1", "1.-exp(-x/[0])", 0, gr->GetXaxis()->GetXmax());
		//Double_t par1[1];
		//funcx->GetParameters(par1);
		funcx1->SetParameters(funcx->GetParameters());
		funcx1->SetLineColor(kGreen);
		funcx1->Draw("same");

		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->SetFillColor(0);
		c->Update();
		TLegend* legend = new TLegend(0.1, 0.83, 0.4, 0.98);
		//legend->SetHeader("The Legend Title", "C");
		legend->AddEntry("func1", "fit function [1]-exp(-#frac{x}{[0]})", "l");
		//legend->AddEntry("func0", "function 1-exp(-x)", "l");
		legend->AddEntry("funcx1", "fit function 1-exp(-#frac{x}{[0]})", "l");
		legend->Draw();
		c->Update();
		MyFile1->cd();
		c->Write(Form("effotachfit&d", i + 1));
		//cin.get();
		//c->Close();
		//canv->Close();
	}
	/*MyFile1->Close();
	TGraphErrors* gr = new TGraphErrors(9, counter1, parametr1, nolliki, erre);
	gr->SetMarkerStyle(21);
	gr->SetMarkerSize(2.0);
	gr->SetTitle(";counter;max efficiency");
	gr->GetYaxis()->SetRangeUser(0.997,1.005);
	gr->Draw("ap");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	cin.get();*/

}

void effcomp1() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_exp**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("tchr", &tchr);
	chain.SetBranchAddress("schr", &schr);
	chain.SetBranchAddress("region", &region);

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;


	char namen[20];
	char title[100];
	TH1* hn[9][14][9];
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				sprintf(namen, "hprofn%d", 14 * 9 * j + 9 * k + i);
				sprintf(title, "sens%d,z%d,phi%d;eff;run", i + 1, k + 1, j + 1);
				hn[j][k][i] = new TH1F(namen, title, 100, -0.1, 2.);
			}
		}
	}

	/*TH1* h[9][14][9];
	char name[20];
	char title[100];
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				sprintf(name, "hprof%d", 14 * 9 * j + 9 * k + i);
				sprintf(title, "sens%d,z%d,phi%d;ach", i + 1, k + 1, j + 1);
				h[j][k][i] = new TH1F(name, title, 100, -1, 1);
			}
		}
	}

	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedest.root");
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				TH1* hped = (TH1F*)f->Get(Form("ped_sens%d,z%d,phi%d", i + 1, k + 1, j + 1));
				ped1[j][k][i] = hped->GetMean(1);
			}
		}
	}*/
	double lgran[4];
	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.25 / tan(theta[0]) + z0[0];
		schr[nch] = 0.;
		if ((ztr > zobl[0]) && (ztr < zobl[14]) && (run > 27225) && (goodrnumb[run - 25001] == 1) && (goodreff[run - 25001] == 1) && (region[0] == 1)) {
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9; j++) {
				//wichcounter1[j] = scount1;
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					while (1 + scount < nch) { scount += 1; }
				}
				else
					while (schr[1 + scount] > 0.5) { scount += 1; }


				//lgran[0] = max(truegran1[j][0], truegran1m[j][0]) + 0.03;
				//lgran[1] = min(truegran1[j][1], truegran1m[j][1]) - 0.045;
				//lgran[2] = max(truegran1[j][1], truegran1m[j][1]) + 0.045;
				//lgran[3] = min(truegran1[j][3], truegran1m[j][3]) - 0.03;
				lgran[0] = truegran1[j][0] + 0.03;
				lgran[1] = truegran1[j][1] - 0.045;
				lgran[2] = truegran1[j][1] + 0.045;
				lgran[3] = truegran1[j][3] - 0.03;

				//ach izm i eff izm 9X14X9
				for (int k = 0; k < 14; k++) {
					if ((ztr > zobl[k]) && (ztr < zobl[k + 1])) {
						int it = 9;
						for (int i = 0; i < 3; i++) {
							if ((phi[0] >= (lgran[0] + i * (lgran[1] - lgran[0]) / 3.)) && (phi[0] < (lgran[0] + (i + 1) * (lgran[1] - lgran[0]) / 3.)))
								it = 0;
							if ((phi[0] >= (lgran[2] + i * (lgran[3] - lgran[2]) / 6.)) && (phi[0] < (lgran[2] + (i + 1) * (lgran[3] - lgran[2]) / 6.)))
								it = 0 + 1;
							if ((phi[0] >= (lgran[2] + (i + 3) * (lgran[3] - lgran[2]) / 6.)) && (phi[0] < (lgran[2] + (i + 4) * (lgran[3] - lgran[2]) / 6.)))
								it = 0 + 2;
						}
						if (it < 9) {
							if ((schr[1 + scount1] < 0.5) || ((schr[1 + scount1] > 0.5) && (tchr[1 + scount1] < 110)))
								hn[0][k][j]->Fill(schr[1 + scount1]);
						}
					}
				}

				/*for (int k = 0; k < 14; k++) {
					if ((ztr > zobl[k]) && (ztr < zobl[k + 1])) {
						int it = 9;
						for (int i = 0; i < 3; i++) {
							if ((phi[0] >= (truegran1[j][0]+0.03 + i * (truegran1[j][1]-0.045 - truegran1[j][0]-0.03) / 3.)) && (phi[0] < (truegran1[j][0]+0.03 + (i + 1) * (truegran1[j][1]-0.045 - truegran1[j][0]-0.03) / 3.)))
								it = i;
							if ((phi[0] >= (truegran1[j][1]+0.045 + i * (truegran1[j][3]-0.03 - truegran1[j][1]-0.045) / 6.)) && (phi[0] < (truegran1[j][1]+0.045 + (i + 1) * (truegran1[j][3]-0.03 - truegran1[j][1]-0.045) / 6.)))
								it = i + 3;
							if ((phi[0] >= (truegran1[j][1]+0.045 + (i + 3) * (truegran1[j][3]-0.03 - truegran1[j][1]-0.045) / 6.)) && (phi[0] < (truegran1[j][1]+0.045 + (i + 4) * (truegran1[j][3]-0.03 - truegran1[j][1]-0.045) / 6.)))
								it = i + 6;
						}
						if (it < 9)
							h[it][k][j]->Fill(ach[scount1]);
					}
				}*/

				scount += 1;
				scount1 = scount;
			}
		}

		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}

	vector<float> errx[9], effect[9], erry[9], ampl[9];
	//effotach


	char titlen[200];
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/random1.root", "RECREATE");
	for (int j = 0; j < 1; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				sprintf(titlen, "exp_eff_sens%d,z%d,phi%d", i + 1, k + 1, j + 1);
				//cout << h[j][k][i]->GetMean(1) << endl;
				hn[j][k][i]->Write(titlen); 
			}
		}
	}
	MyFile->Close();
}
void effcomp2() {
	TChain chain("h1");
	chain.Add("/work/users/kladov/snd2k/R006-003/model/selected/true3_ee**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("beam", &beam);
	chain.SetBranchAddress("tch", &tch);
	chain.SetBranchAddress("region", &region);

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;
	int counterr = 0;

	char namen[20];
	char title[100];
	TH1* hn[9][14][9];
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				sprintf(namen, "hprofn%d", 14 * 9 * j + 9 * k + i);
				sprintf(title, "sens%d,z%d,phi%d;eff;run", i + 1, k + 1, j + 1);
				hn[j][k][i] = new TH1F(namen, title, 100, -0.1, 2.);
			}
		}
	}

	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	double lgran[4];
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.25 / tan(theta[0]) + z0[0];
		if ((ztr > zobl[0]) && (ztr < zobl[14]) && (region[0] == 1)) {
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9; j++) {
				//wichcounter1[j] = scount1;
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					while (1 + scount < nch) { scount += 1; }
				}
				else
					while (schr[1 + scount] > 0.5) { scount += 1; }

				//lgran[0] = max(truegran1[j][0], truegran1m[j][0]) + 0.03;
				//lgran[1] = min(truegran1[j][1], truegran1m[j][1]) - 0.045;
				//lgran[2] = max(truegran1[j][1], truegran1m[j][1]) + 0.045;
				//lgran[3] = min(truegran1[j][3], truegran1m[j][3]) - 0.03;
				lgran[0] = truegran1m[j][0] + 0.03;
				lgran[1] = truegran1m[j][1] - 0.045;
				lgran[2] = truegran1m[j][1] + 0.045;
				lgran[3] = truegran1m[j][3] - 0.03;

				//ach izm i eff izm 9X14X9
				for (int k = 0; k < 14; k++) {
					if ((ztr > zobl[k]) && (ztr < zobl[k + 1])) {
						int it = 9;
						for (int i = 0; i < 3; i++) {
							if ((phi[0] >= (lgran[0] + i * (lgran[1] - lgran[0]) / 3.)) && (phi[0] < (lgran[0] + (i + 1) * (lgran[1] - lgran[0]) / 3.)))
								it = 0;
							if ((phi[0] >= (lgran[2] + i * (lgran[3] - lgran[2]) / 6.)) && (phi[0] < (lgran[2] + (i + 1) * (lgran[3] - lgran[2]) / 6.)))
								it = 0 + 1;
							if ((phi[0] >= (lgran[2] + (i + 3) * (lgran[3] - lgran[2]) / 6.)) && (phi[0] < (lgran[2] + (i + 4) * (lgran[3] - lgran[2]) / 6.)))
								it = 0 + 2;
						}
						if (it < 9) {
							if (ach[j] > 0.2)
								hn[0][k][j]->Fill(1);
							else {
								hn[0][k][j]->Fill(0);
								counterr += 1;
							}
						}
					}
				}

				scount += 1;
				scount1 = scount;
			}
		}

		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}
	cout << counterr << endl;

	vector<float> errx[9], effect[9], erry[9], ampl[9];
	//effotach


	char titlen[200];
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/random2.root", "RECREATE");
	for (int j = 0; j < 1; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				sprintf(titlen, "mod_eff_sens%d,z%d,phi%d", i + 1, k + 1, j + 1);
				//cout << h[j][k][i]->GetMean(1) << endl;
				hn[j][k][i]->Write(titlen);
			}
		}
	}
	MyFile->Close();
}

void effcalib() {
	vector<float> erre[9], effm[9], errm[9], effe[9];
	TFile* f1 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/random1.root");
	TFile* f2 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/random2.root");
	int schit = 0;
	for (int j = 0; j < 1; j++) {
		for (int k = 1; k < 13; k++) {
			for (int i = 0; i < 9; i++) {
				TH1F* he = (TH1F*)f1->Get(Form("exp_eff_sens%d,z%d,phi%d", i + 1, k + 1, j + 1));
				TH1F* hm = (TH1F*)f2->Get(Form("mod_eff_sens%d,z%d,phi%d", i + 1, k + 1, j + 1));
				//if (he->GetMean(1) > 0.6) {
					schit = schit + 1;
					effe[i].push_back(he->GetMean(1));
					effm[i].push_back(hm->GetMean(1));
					erre[i].push_back(he->GetMeanError(1));
					errm[i].push_back(hm->GetMeanError(1));
				//}
			}
		}
	}
	int numbrr = 0;
	//TFile* MyFile1 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/effoteff.root", "RECREATE");
	for (int i = 0; i < 9; i++) {
		//counter1[i] = (double)i + 1;
		//nolliki[i] = 0.;
		numbrr = effe[i].size();
		TGraphErrors* gr = new TGraphErrors(numbrr, &effe[i][0], &effm[i][0], &erre[i][0], &errm[i][0]);
		gr->SetTitle(";effe;effm");
		gr->GetXaxis()->SetRangeUser(0.9, 1.01);
		gr->SetMaximum(1.01);
		gr->SetMinimum(0.9);
		gr->GetXaxis()->SetLimits(0.9, 1.01);
		gr->GetYaxis()->SetRangeUser(0.9, 1.01);
		gr->Draw("AP");
		TF1* func = new TF1("func", "1+[0]*(x-[1])", 0.8, 1);
		func->SetLineColor(kRed);
		func->SetParameter(0, 1.); 
		func->SetParameter(1, 1.);
		func->SetParLimits(0, 0.8, 1.5);
		func->SetParLimits(1, 1, 1);
		gStyle->SetOptFit(1111);


		gr->Fit(func, "", "", 0.8, gr->GetXaxis()->GetXmax());

		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->SetFillColor(0);
		c->SetFrameFillColor(0);
		c->Update();
		cin.get();
		//c->Close();
	}

}


void zraspr() {
	TChain chain("h1");
	chain.Add("/work/users/kladov/snd2k/R006-003/model/selected/true3_ee**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("beam", &beam);
	chain.SetBranchAddress("tch", &tch);
	chain.SetBranchAddress("region", &region);

	int scount = 0;
	int scount1 = 0;

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;


	char name[20];
	char title[100];

	

	TProfile* hprof[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprof%d", i + 1);
		sprintf(title, "sensor%d;phi;ach", i + 1);
		hprof[i] = new TProfile(name, title, 300, -15, 15);
	}
	TProfile* hprof1 = new TProfile("hpprof", "all;ztr;ach", 300, -15, 15);

	TProfile* hprof0[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprof0%d", i + 1);
		sprintf(title, "sensor0%d;phi;ach", i + 1);
		hprof0[i] = new TProfile(name, title, 300, -15, 15);
	}
	TProfile* hprof01 = new TProfile("hpprof0", "all;ztr;ach", 300, -15, 15);

	TH1* hr = new TH1F("histr", "recount right;z,cm", 100, 5., 15.);
	TH1* hl = new TH1F("histl", "recount left;z,cm", 100, -16., -6.);


	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	float ztr1 = 0.;
	float zin = 0.;
	float zin1 = 0.;
	double lgran[4];
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.0 / tan(theta[0]) + z0[0];
		ztr1 = 12.0 / tan(theta[1]) + z0[1];
		zin = 10.5 / tan(theta[0]) + z0[0];
		zin1 = 10.5 / tan(theta[1]) + z0[1];
		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			//wichcounter1[j] = scount1;
			if (j == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
			}

			//lgran[0] = max(truegran1[j][0], truegran1m[j][0]) + 0.03;
			//lgran[1] = min(truegran1[j][1], truegran1m[j][1]) - 0.045;
			//lgran[2] = max(truegran1[j][1], truegran1m[j][1]) + 0.045;
			//lgran[3] = min(truegran1[j][3], truegran1m[j][3]) - 0.03;
			lgran[0] = truegran1m[j][0];// +0.03;
			lgran[1] = truegran1m[j][1];// -0.045;
			lgran[2] = truegran1m[j][1];// +0.045;
			lgran[3] = truegran1m[j][3];// -0.03;

			pedz[j] = 0.;
			if ((((phi[0] > lgran[0]) && (phi[0] < lgran[1])) || ((phi[0] > lgran[2]) && (phi[0] < lgran[3]))) && (region[0] == 1)) {
				if (1) {
					ach1[scount] = (ach[scount] - pedz[j]) * sin(theta[0]);
					hprof0[j]->Fill(ztr, ach1[scount]);
					hprof01->Fill(ztr, ach1[scount]);
					if (zin < 0) {
						if (((11.8 - fabs(zin)) * fabs(tan(theta[0])) < 3.0) && ((11.8 - fabs(zin)) * fabs(tan(theta[0])) >= 0.0)) {
							ach1[scount] = (ach[scount] - pedz[j]) * 3.0 * fabs(cos(theta[0])) / (11.8 - fabs(zin));
							hl->Fill((zin - 11.8) / 2.);
							hprof[j]->Fill((zin - 11.8) / 2., ach1[scount]);
							hprof1->Fill((zin - 11.8) / 2., ach1[scount]);
							//hprof[j]->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
							//hprof1->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
						}
						else if ((11.8 - fabs(zin)) * fabs(tan(theta[0])) >= 3.0) {
							ach1[scount] = (ach[scount] - pedz[j]) * sin(theta[0]);
							hprof[j]->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
							hprof1->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
						}
						else {
							ach1[scount] = ach[scount] - pedz[j];
							hprof[j]->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
							hprof1->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
						}
					}
					if (zin >= 0) {
						if (((10.8 - fabs(zin)) * fabs(tan(theta[0])) < 3.0) && ((10.8 - fabs(zin)) * fabs(tan(theta[0])) >= 0.0)) {
							ach1[scount] = (ach[scount] - pedz[j]) * 3.0 * fabs(cos(theta[0])) / (10.8 - fabs(zin));
							hr->Fill((zin + 10.8) / 2.);
							hprof[j]->Fill((zin + 10.8) / 2., ach1[scount]);
							hprof1->Fill((zin + 10.8) / 2., ach1[scount]);
							//hprof[j]->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
							//hprof1->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
						}
						else if ((10.8 - fabs(zin)) * fabs(tan(theta[0])) >= 3.0) {
							ach1[scount] = (ach[scount] - pedz[j]) * sin(theta[0]);
							hprof[j]->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
							hprof1->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
						}
						else {
							ach1[scount] = ach[scount] - pedz[j];
							hprof[j]->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
							hprof1->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
						}
					}
				}
			}
			if ((((phi[1] > lgran[0]) && (phi[1] < lgran[1])) || ((phi[1] > lgran[2]) && (phi[1] < lgran[3]))) && (region[1] == 1)) {
				if (1) {
					ach1[scount] = (ach[scount] - pedz[j]) * sin(theta[1]);
					hprof0[j]->Fill(ztr1, ach1[scount]);
					hprof01->Fill(ztr1, ach1[scount]);
					if (zin1 < 0) {
						if (((11.8 - fabs(zin1)) * fabs(tan(theta[1])) < 3.0) && ((11.8 - fabs(zin1)) * fabs(tan(theta[1])) >= 0.0)) {
							ach1[scount] = (ach[scount] - pedz[j]) * 3.0 * fabs(cos(theta[1])) / (11.8 - fabs(zin1));
							hl->Fill((zin1 - 11.8) / 2.);
							hprof[j]->Fill((zin1 - 11.8) / 2., ach1[scount]);
							hprof1->Fill((zin1 - 11.8) / 2., ach1[scount]);
							//hprof[j]->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
							//hprof1->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
						}
						else if ((11.8 - fabs(zin1)) * fabs(tan(theta[1])) >= 3.0) {
							ach1[scount] = (ach[scount] - pedz[j]) * sin(theta[1]);
							hprof[j]->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
							hprof1->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
						}
						else {
							ach1[scount] = ach[scount] - pedz[j];
							hprof[j]->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
							hprof1->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
						}
					}
					if (zin1 >= 0) {
						if (((10.8 - fabs(zin1)) * fabs(tan(theta[1])) < 3.0) && ((10.8 - fabs(zin1)) * fabs(tan(theta[1])) >= 0.0)) {
							ach1[scount] = (ach[scount] - pedz[j]) * 3.0 * fabs(cos(theta[1])) / (10.8 - fabs(zin1));
							hr->Fill((zin1 + 10.8) / 2.);
							hprof[j]->Fill((zin1 + 10.8) / 2., ach1[scount]);
							hprof1->Fill((zin1 + 10.8) / 2., ach1[scount]);
							//hprof[j]->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
							//hprof1->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
						}
						else if ((10.8 - fabs(zin1)) * fabs(tan(theta[1])) >= 3.0) {
							ach1[scount] = (ach[scount] - pedz[j]) * sin(theta[1]);
							hprof[j]->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
							hprof1->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
						}
						else {
							ach1[scount] = ach[scount] - pedz[j];
							hprof[j]->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
							hprof1->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
						}
					}
				}
			}
			scount += 1;
			scount1 = scount;
		}
		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}


	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/model/zprofiles.root", "RECREATE");
	for (Int_t i = 0; i < 9; i++) {
		sprintf(title, "sensor%d",i + 1);
		hprof[i]->Write(title);
		sprintf(title, "sensor0%d", i + 1);
		hprof0[i]->Write(title);
	}
	hprof1->Write("allsensors");
	hprof01->Write("allsensors0");
	hr->Write("zinr");
	hl->Write("zinl");
	MyFile->Close();

}

void zsravn() {
	//float zobl[15] = { -11.6013, -10.0367, -8.47203, -6.90737, -5.34273, -3.77807, -2.21343, -0.648766, 0.9159, 2.48053, 4.0452, 5.60983, 7.1745, 8.73913, 9.6 };
	TFile* MyFile1 = new TFile("/work/users/kladov/snd2k/R006-003/model/zprofiles.root");
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/zprofiles.root");
	TProfile* he = (TProfile*)MyFile->Get("allsensors");
	TProfile* he0 = (TProfile*)MyFile->Get("allsensors0");
	he->SetLineColor(2);
	he0->SetLineColor(1);
	he->SetTitle("blue - modeling, red - experiment;z,sm;mean amplitude, pe");
	he->Draw();
	TProfile* hm = (TProfile*)MyFile1->Get("allsensors");
	TProfile* hm0 = (TProfile*)MyFile1->Get("allsensors0");
	hm->SetLineColor(4);
	hm0->SetLineColor(3);
	hm->Draw("same");
	hm0->Draw("same");
	he0->Draw("same");
	
	/*TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	TLine* linez[15];
	for (int i = 0; i < 15; i++) {
		linez[i] = new TLine(zobl[i], 0, zobl[i], 1.05 * he->GetMaximum());
		linez[i]->SetLineColor(kGreen);
		linez[i]->Draw();
	}
	c->Update();
	MyFile1->Close();
	MyFile->Close();

	TFile* Myf = new TFile("/work/users/kladov/snd2k/R006-003/model/achspectreme.root");
	for (int i = 10; i < 11; i++) {
		TH1* he = (TH1F*)Myf->Get(Form("e_zobl%d", i + 1));
		TH1* hm = (TH1F*)Myf->Get(Form("m_zobl%d", i + 1));
		hm->SetLineColor(4);
		he->SetLineColor(2);
		he->DrawNormalized("", 1);
		hm->DrawNormalized("same", 1);
		cout << he->GetMean(1) << endl << hm->GetMean();
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
	}*/
}

void linesforz() {
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/model/zprofiles.root");

	TProfile* hprof = (TProfile*)f->Get("allsensors");
	TProfile* hprof0 = (TProfile*)f->Get("allsensors0");
	TH1* hr = (TH1F*)f->Get("zinr");
	TH1* hl = (TH1F*)f->Get("zinl");
	hr->SetLineColor(kRed);
	hl->SetLineColor(kRed);
	hprof->SetTitle("amplitude distribution by z;z,ñm;amplitude,pe");
	hprof->SetLineColor(kBlack);
	hprof0->SetLineColor(kBlue);
	hprof->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	/*TLine* linez[15];
	for (int i = 0; i < 15; i++) {
		linez[i] = new TLine(zobl[i], 0, zobl[i], 1.05 * hprof->GetMaximum());
		linez[i]->SetLineColor(kGreen);
		linez[i]->Draw();
	}*/
	c->Update();
	TLine* line11 = new TLine(-11.8, 0, -11.8, 7);
	TLine* line12 = new TLine(10.8, 0, 10.8, 7);
	line11->SetLineColor(kGreen);
	line12->SetLineColor(kGreen);
	line11->Draw();
	line12->Draw();
	hr->DrawNormalized("same", 50);
	hl->DrawNormalized("same", 50);
	hprof0->Draw("same");
	c->Update();
}

void raspr() {

	TChain chain("h1");
	chain.Add("/work/users/kladov/snd2k/R006-003/model/selected/true3_ee**.root");
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("beam", &beam);
	chain.SetBranchAddress("tch", &tch);

	int scount = 0;
	int scount1 = 0;

	//raspr po phi
	/*TProfile* hprof[14][9];
	char name[20];
	char title[100];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "hprof%d", 9*j + i + 1);
			sprintf(title, "zobl%d,sensor%d;phi;ach", j + 1,i + 1);
			hprof[j][i] = new TProfile(name,title,200, (float)(i - 1.) * (2. * PI) / 9., (float)(i + 2.) * (2. * PI) / 9.);
		}
	}*/
	TProfile* hprof[9];
	char name[20];
	char title[100];
	for (int j = 0; j < 9; j++) {
		sprintf(name, "hprof%d", j + 1);
		sprintf(title, "counter%d;phi;ach", j + 1);
		hprof[j] = new TProfile(name, title, 200, (float)(j - 1.) * (2. * PI) / 9., (float)(j + 2.) * (2. * PI) / 9.);
	}

	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	float ztr1 = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//cout << beam << endl;
		ztr = 12.25 / tan(theta[0]) + z0[0];
		ztr1 = 12.25 / tan(theta[1]) + z0[1];
		if (1) {
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9; j++) {
				//wichcounter1[j] = scount1;
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
				}

				if ((phi[0] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[0] < (float)(j + 2) * (2. * PI) / 9.) && (ztr > zobl[0]) && (ztr < zobl[14]))
					hprof[j]->Fill(phi[0], ach[j] * sin(theta[0]));
				if ((phi[1] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[1] < (float)(j + 2) * (2. * PI) / 9.) && (ztr1 > zobl[0]) && (ztr1 < zobl[14]))

				scount += 1;
				scount1 = scount;
			}
		}
		
		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}

	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/model/profiles.root", "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "zobl%d,sensor%d", m + 1, i + 1);
			hprof[m][i]->Write(title);
		}
	}
	MyFile->Close();*/

	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/model/profilesmean.root", "RECREATE");
	for (int m = 0; m < 9; m++) {
		sprintf(title, "counter%d", m + 1);
		hprof[m]->Write(title);
	}
	MyFile->Close();

}

double compare1() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("tchr", &tchr);
	chain.SetBranchAddress("schr", &schr);
	chain.SetBranchAddress("region", &region);

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;
	
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			ecountc1[i] = 0;
			ecountc2[i] = 0;
			ecounta1[j][i] = 0;
			ecounta2[j][i] = 0;
			ecountz1[j] = 0;
			ecountz2[j] = 0;
		}
	}

	char name[20];
	char title[100];
	TProfile* hprof[9];
	TH1* meancampl[14];
	TH1* meanzampl[9];
	TH1* h[9][14];
	for (int j = 0; j < 9; j++) {
		sprintf(name, "histzm%d", j + 1);
		sprintf(title, "meanz,count%d;ach,pe", j + 1);
		meanzampl[j] = new TH1F(name, title, 500, -5., 45.);
	}
	for (int i = 0; i < 14; i++) {
		sprintf(name, "histcm%d", i + 1);
		sprintf(title, "meanc,obl%d;ach,pe", i + 1);
		meancampl[i] = new TH1F(name, title, 500, -5., 45.);
	}
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			sprintf(name, "allhistsm%d", 9 * i + j + 1);
			sprintf(title, "all,obl%d,count%d;ach,pe", i + 1, j + 1);
			h[j][i] = new TH1F(name, title, 500, -5, 45);
		}
	}

	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/phintegr_pedest.root");
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			TH1* hped = (TH1F*)f->Get(Form("ped_phintegr_zobl%d,sensor%d", j + 1, i + 1));
			ped[j][i] = hped->GetMean(1);
		}
	}
	double lgran[4];
	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	float ztr1 = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.25 / tan(theta[0]) + z0[0];
		ztr1 = 12.25 / tan(theta[1]) + z0[1];
		if ((run > 27225) && (goodrnumb[run - 25001] == 1) && (goodreff[run - 25001] == 1)) {
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9; j++) {
				//wichcounter1[j] = scount1;
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					while (1 + scount < nch) { scount += 1; }
				}
				else
					while (schr[1 + scount] > 0.5) { scount += 1; }

				//lgran[0] = max(truegran1[j][0], truegran1m[j][0])+0.03;
				//lgran[1] = min(truegran1[j][1], truegran1m[j][1])-0.045;
				//lgran[2] = max(truegran1[j][1], truegran1m[j][1])+0.045;
				//lgran[3] = min(truegran1[j][3], truegran1m[j][3])-0.03;
				lgran[0] = truegran1[j][0];// +0.03;
				lgran[1] = truegran1[j][1];// -0.045;
				lgran[2] = truegran1[j][1];// +0.045;
				lgran[3] = truegran1[j][3];// -0.03;

				if ((((phi[0] > lgran[0]) && (phi[0] < lgran[1])) || ((phi[0] > lgran[2]) && (phi[0] < lgran[3]))) && (region[0] == 1)) {
					for (int i = 0; i < 14; i++) {
						if ((ztr > zobl[i]) && (ztr < zobl[i + 1]) && ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110)))) {
							ach1[scount] = ach[scount] - ped[i][j];
							meancampl[i]->Fill(ach1[scount]);
							meanzampl[j]->Fill(ach1[scount]);
							h[j][i]->Fill(ach1[scount]);
							if (ach[scount] - ped[i][j] >= 0.2) {
								if (ztr > zobl[2])
									ecountz1[j] += 1;
								ecounta1[j][i] += 1;
								ecountc1[i] += 1;
							}
							if (ach[scount] - ped[i][j] < 0.2) {
								if (ztr > zobl[2])
									ecountz2[j] += 1;
								ecounta2[j][i] += 1;
								ecountc2[i] += 1;
							}
						}
					}
				}
				if ((((phi[1] > lgran[0]) && (phi[1] < lgran[1])) || ((phi[1] > lgran[2]) && (phi[1] < lgran[3]))) && (region[1] == 1)) {
					for (int i = 0; i < 14; i++) {
						if ((ztr1 > zobl[i]) && (ztr1 < zobl[i + 1]) && ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110)))) {
							ach1[scount] = ach[scount] - ped[i][j];
							meancampl[i]->Fill(ach1[scount]);
							meanzampl[j]->Fill(ach1[scount]);
							h[j][i]->Fill(ach1[scount]);
							if (ach[scount] - ped[i][j] >= 0.2) {
								if (ztr1 > zobl[2])
									ecountz1[j] += 1;
								ecounta1[j][i] += 1;
								ecountc1[i] += 1;
							}
							if (ach[scount] - ped[i][j] < 0.2) {
								if (ztr1 > zobl[2])
									ecountz2[j] += 1;
								ecounta2[j][i] += 1;
								ecountc2[i] += 1;
							}
						}
					}
				}
				scount += 1;
				scount1 = scount;
			}
		}
		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}
	{
		TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/model/achspectreme.root", "RECREATE");
		for (int m = 0; m < 14; m++) {
			sprintf(title, "e_zobl%d", m + 1);
			meancampl[m]->Write(title);
		}
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++) {
				sprintf(title, "e_all_zobl%d_count%d", i + 1, j + 1);
				h[j][i]->Write(title);
			}
		}
		for (int j = 0; j < 9; j++) {
			sprintf(title, "e_count%d", j + 1);
			meanzampl[j]->Write(title);
		}
		MyFile->Close();
	}

	{
		ofstream fout;
		fout.open("/work/users/kladov/snd2k/R006-003/model/counteri.cpp");
		fout << "double ecounta1[9][14] = { ";
		for (int j = 0; j < 8; j++) {
			for (int i = 0; i < 14; i++) {
				fout << ecounta1[j][i] << ", ";
			}
			fout << (char)92 << endl;
		}
		for (int i = 0; i < 13; i++) {
			fout << ecounta1[8][i] << ", ";
		}
		fout << ecounta1[8][13] << " };" << endl;

		fout << "double ecounta2[9][14] = { ";
		for (int j = 0; j < 8; j++) {
			for (int i = 0; i < 14; i++) {
				fout << ecounta2[j][i] << ", ";
			}
			fout << (char)92 << endl;
		}
		for (int i = 0; i < 13; i++) {
			fout << ecounta2[8][i] << ", ";
		}
		fout << ecounta2[8][13] << " };" << endl;

		fout << "double ecountc1[14] = { ";
		for (int i = 0; i < 13; i++) {
			fout << ecountc1[i] << ", ";
		}
		fout << ecountc1[13] << " };" << endl;
		fout << "double ecountc2[14] = { ";
		for (int i = 0; i < 13; i++) {
			fout << ecountc2[i] << ", ";
		}
		fout << ecountc2[13] << " };" << endl;

		fout << "double ecountz1[9] = { ";
		for (int j = 0; j < 8; j++) {
			fout << ecountz1[j] << ", ";
		}
		fout << ecountz1[8] << " };" << endl;
		fout << "double ecountz2[9] = { ";
		for (int j = 0; j < 8; j++) {
			fout << ecountz2[j] << ", ";
		}
		fout << ecountz2[8] << " };" << endl;
	}

	/*for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			cout << ecountc1[i] << endl;
			cout << ecountc2[i] << endl;
			cout << ecounta1[j][i] << endl;
			cout << ecounta2[j][i] << endl;
			cout << ecountz1[j] << endl;
			cout << ecountz2[j] << endl;
		}
	}*/
	return 1;

	/*char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* fm = new TFile("/work/users/kladov/snd2k/R006-003/model/profiles.root");
	TFile* fe = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles.root");
	TIter next(fe->GetListOfKeys());
	TKey* key;
	int counterr = 0, ennn = 0;
	vector<double> meanache[14];
	vector<double> meanachm[14];
	int ae[14]; double be[14];
	int am[14]; double bm[14];
	double lgranm[4];
	double lgrane[4];
	while ((key = (TKey*)next())) {
		counterr += 1;
		TString name(key->GetName());
		TString type = "zobl";
		if (name.BeginsWith(type) && name.Contains("sensor")) {
			//cin.get();
			int obl, counter;
			sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
			//counter = 7;
			//obl = (counterr-1)/4+1;
			TProfile* h_amp_phie = (TProfile*)fe->Get(Form("zobl%d,sensor%d", obl, counter));
			TProfile* h_amp_phim = (TProfile*)fm->Get(Form("zobl%d,sensor%d", obl, counter));
			h_amp_phie->SetLineColor(3);
			h_amp_phim->SetLineColor(4);

			h_amp_phie->Draw();
			h_amp_phim->Draw("same");
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();

			lgranm[0] = truegran1m[counter - 1][0] + 0.04;
			lgranm[1] = truegran1m[counter - 1][1] - 0.05;
			lgranm[2] = truegran1m[counter - 1][1] + 0.05;
			lgranm[3] = truegran1m[counter - 1][3] - 0.04;
			lgrane[0] = truegran1[counter - 1][0] + 0.04;
			lgrane[1] = truegran1[counter - 1][1] - 0.05;
			lgrane[2] = truegran1[counter - 1][1] + 0.05;
			lgrane[3] = truegran1[counter - 1][3] - 0.04;

			TLine* line11 = new TLine(lgran[0], 0, lgran[0], gPad->GetUymax());
			TLine* line12 = new TLine(lgran[1], 0, lgran[1], gPad->GetUymax());
			TLine* line13 = new TLine(lgran[2], 0, lgran[2], gPad->GetUymax());
			TLine* line14 = new TLine(lgran[3], 0, lgran[3], gPad->GetUymax());
			line11->SetLineColor(kRed);
			line12->SetLineColor(kRed);
			line13->SetLineColor(kRed);
			line14->SetLineColor(kRed);
			line11->Draw();
			line12->Draw();
			line13->Draw();
			line14->Draw();
			c->Update();
			cout << "ready" << endl;
			//cin.get();

			
			double maxphi;
			for (int i = 20; i < 180; i++) {
				maxphi = (h_amp_phie->GetXaxis()->GetBinLowEdge(i) + h_amp_phie->GetXaxis()->GetBinUpEdge(i)) / 2.;
				if (((maxphi> lgranm[0]) && (maxphi < lgranm[1])) || ((maxphi > lgranm[2]) && (maxphi < lgranm[3]))) 
					meanachm[obl - 1].push_back(h_amp_phim->GetBinContent(i));
				if (((maxphi > lgrane[0]) && (maxphi < lgrane[1])) || ((maxphi > lgrane[2]) && (maxphi < lgrane[3]))) 
					meanache[obl - 1].push_back(h_amp_phie->GetBinContent(i));
			}
			
		}
	}
	

	for (int k = 0; k < 14; k++) {
		be[k] = 0;
		bm[k] = 0;
		ae[k] = meanache[k].size();
		//cout << a << endl;
		for (int i = 0; i < ae[k]; i++)
			be[k] = be[k] + meanache[k][i];
		//cout << b << endl;
		am[k] = meanachm[k].size();
		for (int i = 0; i < am[k]; i++)
			bm[k] = bm[k] + meanachm[k][i];
		cout << Form("obl	%d",k+1) << endl;
		//cout << be[k] << endl << bm[k] << endl;
		//cout << "exp -> " << (be[k] / (double)ae[k]) << endl;
		//cout << "mod -> " << (bm[k] / (double)am[k]) << endl;
		cout << (bm[k] / (double)am[k]) / (be[k] / (double)ae[k]) << endl;
	}
	for (int k = 0; k < 14; k++) {
		cout << (bm[k] / (double)am[k]) / (be[k] / (double)ae[k]) << " , ";
	}
	cout << endl;
	
	return 3;*/
}

double compare2() {
	TChain chain("h1");
	chain.Add("/work/users/kladov/snd2k/R006-003/model/selected/true3_ee**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("beam", &beam);
	chain.SetBranchAddress("tch", &tch);
	chain.SetBranchAddress("region", &region);

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			countc1[i] = 0;
			countc2[i] = 0;
			counta1[j][i] = 0;
			counta2[j][i] = 0;
			countz1[j] = 0;
			countz2[j] = 0;
		}
	}

	char name[20];
	char title[100];
	TProfile* hprof[9];
	TH1* meancampl[14];
	TH1* meanzampl[9];
	TH1* h[9][14];
	for (int j = 0; j < 9; j++) {
		sprintf(name, "histzm%d", j + 1);
		sprintf(title, "meanz,count%d;ach,pe", j + 1);
		meanzampl[j] = new TH1F(name, title, 500, -5., 45.);
	}
	for (int i = 0; i < 14; i++) {
		sprintf(name, "histcm%d", i + 1);
		sprintf(title, "meanc,obl%d;ach,pe", i + 1);
		meancampl[i] = new TH1F(name, title, 500, -5., 45.);
	}
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			sprintf(name, "allhistsm%d", 9 * i + j + 1);
			sprintf(title, "all,obl%d,count%d;ach,pe", i + 1, j + 1);
			h[j][i] = new TH1F(name, title, 500, -5., 45.);
		}
	}
	
	double lgran[4];
	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	float ztr1 = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.25 / tan(theta[0]) + z0[0];
		ztr1 = 12.25 / tan(theta[1]) + z0[1];
		if (1) {
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9; j++) {
				//wichcounter1[j] = scount1;
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					while (1 + scount < nch) { scount += 1; }
				}
				else
					while (schr[1 + scount] > 0.5) { scount += 1; }

				//lgran[0] = max(truegran1[j][0], truegran1m[j][0]) + 0.03;
				//lgran[1] = min(truegran1[j][1], truegran1m[j][1]) - 0.045;
				//lgran[2] = max(truegran1[j][1], truegran1m[j][1]) + 0.045;
				//lgran[3] = min(truegran1[j][3], truegran1m[j][3]) - 0.03;
				lgran[0] = truegran1m[j][0];// +0.03;
				lgran[1] = truegran1m[j][1];// -0.045;
				lgran[2] = truegran1m[j][1];// +0.045;
				lgran[3] = truegran1m[j][3];// -0.03;

				if ((((phi[0] > lgran[0]) && (phi[0] < lgran[1])) || ((phi[0] > lgran[2]) && (phi[0] < lgran[3]))) && (region[0] == 1)) {
					for (int i = 0; i < 14; i++) {
						if ((ztr > zobl[i]) && (ztr < zobl[i + 1])/* && ((ach[j] > 0.5) && (tch[j] < 1000))*/) {
							meancampl[i]->Fill(ach[j]);
							meanzampl[j]->Fill(ach[j]);
							h[j][i]->Fill(ach[j]);
							if (ach[j] >= 0.2) {
								if (ztr > zobl[2])
									countz1[j] += 1;
								counta1[j][i] += 1;	
								countc1[i] += 1;
							}
							if (ach[j] < 0.2) {
								if (ztr > zobl[2])
									countz2[j] += 1;
								counta2[j][i] += 1;
								countc2[i] += 1;
							}

						}
					}
				}
				if ((((phi[1] > lgran[0]) && (phi[1] < lgran[1])) || ((phi[1] > lgran[2]) && (phi[1] < lgran[3]))) && (region[1] == 1)) {
					for (int i = 0; i < 14; i++) {
						if ((ztr1 > zobl[i]) && (ztr1 < zobl[i + 1])/* && ((ach[j] > 0.5) && (tch[j] < 1000))*/) {
							meancampl[i]->Fill(ach[j]);
							meanzampl[j]->Fill(ach[j]);
							h[j][i]->Fill(ach[j]);
							if (ach[j] >= 0.2) {
								if (ztr1 > zobl[2])
									countz1[j] += 1;
								counta1[j][i] += 1;
								countc1[i] += 1;
							}
							if (ach[j] < 0.2) {
								if(ztr1>zobl[2])
									countz2[j] += 1;
								counta2[j][i] += 1;
								countc2[i] += 1;
							}

						}
					}
				}

				scount += 1;
				scount1 = scount;
			}
		}
		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}
	{
		TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/model/achspectreme.root", "UPDATE");
		for (int m = 0; m < 14; m++) {
			sprintf(title, "m_zobl%d", m + 1);
			meancampl[m]->Write(title);
		}
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++) {
				sprintf(title, "m_all_zobl%d_count%d", i + 1, j + 1);
				h[j][i]->Write(title);
			}
		}
		for (int j = 0; j < 9; j++) {
			sprintf(title, "m_count%d", j + 1);
			meanzampl[j]->Write(title);
		}
		MyFile->Close();
	}

	{
		ofstream fout;
		fout.open("/work/users/kladov/snd2k/R006-003/model/counteri.cpp", ios::out | ios::app);
		fout << "double counta1[9][14] = { ";
		for (int j = 0; j < 8; j++) {
			for (int i = 0; i < 14; i++) {
				fout << counta1[j][i] << ", ";
			}
			fout << (char)92 << endl;
		}
		for (int i = 0; i < 13; i++) {
			fout << counta1[8][i] << ", ";
		}
		fout << counta1[8][13] << " };" << endl;

		fout << "double counta2[9][14] = { ";
		for (int j = 0; j < 8; j++) {
			for (int i = 0; i < 14; i++) {
				fout << counta2[j][i] << ", ";
			}
			fout << (char)92 << endl;
		}
		for (int i = 0; i < 13; i++) {
			fout << counta2[8][i] << ", ";
		}
		fout << counta2[8][13] << " };" << endl;

		fout << "double countc1[14] = { ";
		for (int i = 0; i < 13; i++) {
			fout << countc1[i] << ", ";
		}
		fout << countc1[13] << " };" << endl;
		fout << "double countc2[14] = { ";
		for (int i = 0; i < 13; i++) {
			fout << countc2[i] << ", ";
		}
		fout << countc2[13] << " };" << endl;

		fout << "double countz1[9] = { ";
		for (int j = 0; j < 8; j++) {
			fout << countz1[j] << ", ";
		}
		fout << countz1[8] << " };" << endl;
		fout << "double countz2[9] = { ";
		for (int j = 0; j < 8; j++) {
			fout << countz2[j] << ", ";
		}
		fout << countz2[8] << " };" << endl;
	}

	/*for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			cout << countc1[i] << endl;
			cout << countc2[i] << endl;
			cout << counta1[j][i] << endl;
			cout << counta2[j][i] << endl;
			cout << countz1[j] << endl;
			cout << countz2[j] << endl;
		}
	}*/
	return 2;
}

void compare() {
	//compare1();
	//compare2();
	//cout << "popr = " << compare1() << endl;
	TFile* Myf = new TFile("/work/users/kladov/snd2k/R006-003/model/achspectreme.root");
	for (int i = 0; i < 14; i++) {
		TH1* meanamplocke = (TH1F*)Myf->Get(Form("e_zobl%d", i + 1));
		TH1* meanamplockm = (TH1F*)Myf->Get(Form("m_zobl%d", i + 1));
		meanample[i] = meanamplocke->GetMean(1);
		meanampleerr[i] = meanamplocke->GetMeanError(1);
		meanamplm[i] = meanamplockm->GetMean(1);
		meanamplmerr[i] = meanamplockm->GetMeanError(1);
	}

	cout << "popr = " << endl;
	double popr = 0.;
	cout << (meanamplm[0] / meanample[0]) << "    " << meanamplmerr[0] / meanample[0] << "	"/* << meanampleerr[0] */ << endl;
	for (int i = 1; i < 14; i++) {
		cout << (meanamplm[i] / meanample[i]) << "    " <<  meanamplmerr[i] / meanample[i] << "	"/* << meanampleerr[i] */<< endl;
		popr = popr + meanamplm[i] / meanample[i];
	}
	cout << "mean	" << popr/13. << endl;
	double oblast[14];
	double noll[14];
	double koef[14] = { 1.0707799, 1.07214, 1.08057, 1.08138, 1.07637, 1.07359, 1.07303, 1.08661, 1.0774, 1.07921, 1.08933, 1.08917, 1.07887, 1.06691 };
	double koeferr[14];
	double lokk = 0.;
	for (int i = 0; i < 14; i++) {
		oblast[i] = (zobl[i]+zobl[i+1])/2.;
		noll[i] = 0.;
		lokk = koef[i] * (meanamplm[i+0] / meanample[i+0]);
		koef[i] = 100.-100./(lokk);
		koeferr[i] = 100.*lokk * (meanamplmerr[i+0] + meanampleerr[i+0]) / meanamplm[i+0];
		cout << lokk << ", ";
	}
	TGraphErrors* gr = new TGraphErrors(14, &oblast[0], &koef[0], &noll[0], &koeferr[0]);
	TF1* fn = new TF1("fn", "[0]", -10, 10);
	fn->SetParameter(0, 7.);
	fn->SetParLimits(0,2.,11.);
	gr->SetLineColor(1);
	gr->SetMarkerStyle(21);
	gr->SetMarkerColor(1);
	gr->GetYaxis()->SetRangeUser(1., 12.);
	gr->SetTitle("secondary particles contribution;z,sm;percent");
	gr->Draw("APL");
	fn->SetLineColor(kRed);
	gr->Fit("fn", "", "", -9, 9);
	cout << fn->GetParameter(0) << endl;

	/*for (int i = 1; i < 13; i++) {
		for (int j = 0; j < 9; j++) {
			cout << Form("obl %d__count %d", i + 1, j + 1) << endl;
			cout << counta2[j][i]/ (counta1[j][i] + counta2[j][i]) << " +- " << pow(((counta1[j][i] * counta2[j][i])/ pow((counta1[j][i] + counta2[j][i]),3.)),0.5) << endl;
			cout << ecounta2[j][i]/ (ecounta1[j][i] + ecounta2[j][i]) << " +- " << pow(((ecounta1[j][i] * ecounta2[j][i])/ pow((ecounta1[j][i] + ecounta2[j][i]),3.)),0.5) << endl;
			cout << (counta2[j][i] / (counta1[j][i] + counta2[j][i]) - ecounta2[j][i] / (ecounta1[j][i] + ecounta2[j][i])) / pow(((counta1[j][i] * counta2[j][i]) / pow((counta1[j][i] + counta2[j][i]), 3.)), 0.5) << endl << endl << endl;
		}
	}*/
	/*for (int j = 0; j < 9; j++) {
		cout << Form("count %d", j + 1) << endl;
		cout << countz2[j] / (countz1[j] + countz2[j]) << " +- " << pow(((countz1[j] * countz2[j]) / pow((countz1[j] + countz2[j]), 3.)), 0.5) << endl;
		cout << ecountz2[j] / (ecountz1[j] + ecountz2[j]) << " +- " << pow(((ecountz1[j] * ecountz2[j]) / pow((ecountz1[j] + ecountz2[j]), 3.)), 0.5) << endl;
		cout << (countz2[j] / (countz1[j] + countz2[j]) - ecountz2[j] / (ecountz1[j] + ecountz2[j])) / pow(((countz1[j] * countz2[j]) / pow((countz1[j] + countz2[j]), 3.)), 0.5) << endl << endl << endl;
	}*/
	double effm[14];
	double efferrm[14];
	double effe[14];
	double efferre[14];
	for (int j = 0; j < 14; j++) {
		cout << Form("obl %d", j + 1) << endl;
		effm[j] = 100. * (1.-countc2[j] / (countc1[j] + countc2[j]));
		efferrm[j] = 120. * pow(((countc1[j] * countc2[j]) / pow((countc1[j] + countc2[j]), 3.)), 0.5);
		cout << effm[j] << " +- " << efferrm[j] << endl;
		effe[j] = 100. * (1. - ecountc2[j] / (ecountc1[j] + ecountc2[j]));
		efferre[j] = 150.*pow(((ecountc1[j] * ecountc2[j]) / pow((ecountc1[j] + ecountc2[j]), 3.)), 0.5);
		cout << effe[j] << " +- " << efferre[j] << endl;
		//cout << ((countc2[j] / (countc1[j] + countc2[j])) - (ecountc2[j] / (ecountc1[j] + ecountc2[j]))) / pow(((countc1[j] * countc2[j]) / pow((countc1[j] + countc2[j]), 3.)), 0.5) << endl << endl << endl;
		cout << effm[j]/effe[j]-1 << "	+-	" << (efferrm[j]+ efferre[j])/ effe[j] << endl << endl << endl;
	}
	/*double counterr[9];
	double nolll[9];
	for (int j = 0; j < 9; j++) {
		noll[j] = 0.;
		counterr[j] = (double)j + 1;
		cout << Form("counter %d", j + 1) << endl;
		effm[j] = 1. - countz2[j] / (countz1[j] + countz2[j]);
		efferrm[j] = pow(((countz1[j] * countz2[j]) / pow((countz1[j] + countz2[j]), 3.)), 0.5);
		cout << effm[j] << " +- " << efferrm[j] << endl;
		effe[j] = 1. - ecountz2[j] / (ecountz1[j] + ecountz2[j]);
		efferre[j] = pow(((ecountz1[j] * ecountz2[j]) / pow((ecountz1[j] + ecountz2[j]), 3.)), 0.5);
		cout << effe[j] << " +- " << efferre[j] << endl;
		//cout << ((countc2[j] / (countc1[j] + countc2[j])) - (ecountc2[j] / (ecountc1[j] + ecountc2[j]))) / pow(((countc1[j] * countc2[j]) / pow((countc1[j] + countc2[j]), 3.)), 0.5) << endl << endl << endl;
		cout << effm[j] / effe[j] - 1 << "	+-	" << (efferrm[j] + efferre[j]) / effe[j] << endl << endl << endl;
	}*/
	//cin.get();
	//TGraphErrors* greffe = new TGraphErrors(9, counterr, effe, nolll, efferre);
	//TGraphErrors* greffm = new TGraphErrors(9, counterr, effm, nolll, efferrm);
	/*TGraphErrors* greffe = new TGraphErrors(13, &oblast[1], &effe[1], &noll[1], &efferre[1]);
	TGraphErrors* greffm = new TGraphErrors(13, &oblast[1], &effm[1], &noll[1], &efferrm[1]);
	greffe->SetLineColor(2);
	greffm->SetLineColor(4);
	greffe->SetMarkerStyle(21);
	greffm->SetMarkerStyle(20);
	greffe->SetMarkerColor(1);
	greffm->SetMarkerColor(1);
	greffe->GetYaxis()->SetRangeUser(94, 99);
	greffe->SetTitle("detection efficiency;z,sm;efficiency,%");
	greffe->Draw("APL");
	greffm->Draw("PLsame");*/
}


std::vector<double> m4(std::vector<double> x, std::vector<double> y)
{
	double x1 = x[0];
	double x2 = x[1];
	double x3 = x[2];
	double x4 = x[3];
	double y1 = y[0];
	double y2 = y[1];
	double y3 = y[2];
	double y4 = y[3];
	double d = (x1 - x2) * (x1 - x3) * (x2 - x3) * (x1 - x4) * (x2 - x4) * (x3 - x4);
	std::vector<double> p(4, 0);
	p[0] = (x1 * (x1 - x3) * x3 * (x1 - x4) * (x3 - x4) * x4 * y2 +
		x2 * pow(x4, 2) * (-pow(x3, 3) * y1 + pow(x3, 2) * x4 * y1 + pow(x1, 2) * (x1 - x4) * y3) +
		pow(x1, 2) * x2 * pow(x3, 2) * (-x1 + x3) * y4 +
		pow(x2, 3) * (x4 * (-pow(x3, 2) * y1 + x3 * x4 * y1 + x1 * (x1 - x4) * y3) +
			x1 * x3 * (-x1 + x3) * y4) +
		pow(x2, 2) * (x1 * x4 * (-pow(x1, 2) + pow(x4, 2)) * y3 +
			pow(x3, 3) * (x4 * y1 - x1 * y4) + x3 * (-pow(x4, 3) * y1 + pow(x1, 3) * y4))) / d;
	p[1] = (pow(x1, 2) * (x1 - x4) * pow(x4, 2) * (y2 - y3) +
		pow(x3, 3) * (pow(x4, 2) * (y1 - y2) + pow(x1, 2) * (y2 - y4)) +
		pow(x2, 2) * (pow(x4, 3) * (y1 - y3) + pow(x1, 3) * (y3 - y4) + pow(x3, 3) * (-y1 + y4)) +
		pow(x3, 2) * (pow(x4, 3) * (-y1 + y2) + pow(x1, 3) * (-y2 + y4)) +
		pow(x2, 3) * (pow(x4, 2) * (-y1 + y3) + pow(x3, 2) * (y1 - y4) + pow(x1, 2) * (-y3 + y4))) / d;
	p[2] = (-x1 * (x1 - x4) * x4 * (x1 + x4) * (y2 - y3) +
		x3 * (pow(x4, 3) * (y1 - y2) + pow(x1, 3) * (y2 - y4)) +
		pow(x3, 3) * (-x4 * y1 - x1 * y2 + x4 * y2 + x1 * y4) +
		pow(x2, 3) * (x4 * y1 + x1 * y3 - x4 * y3 - x1 * y4 + x3 * (-y1 + y4)) +
		x2 * (pow(x4, 3) * (-y1 + y3) + pow(x3, 3) * (y1 - y4) + pow(x1, 3) * (-y3 + y4))) / d;
	p[3] = (x1 * (x1 - x4) * x4 * (y2 - y3) +
		pow(x3, 2) * (x4 * y1 + x1 * y2 - x4 * y2 - x1 * y4) +
		pow(x2, 2) * (-x4 * y1 - x1 * y3 + x4 * y3 + x3 * (y1 - y4) + x1 * y4) +
		x2 * (pow(x4, 2) * (y1 - y3) + pow(x1, 2) * (y3 - y4) + pow(x3, 2) * (-y1 + y4)) +
		x3 * (pow(x4, 2) * (-y1 + y2) + pow(x1, 2) * (-y2 + y4))) / d;
	return p;
}

double p3g0(double x, double xr, double sr, std::vector<double> p)
{
	double p0, p1, p2, p3;
	double arg;
	p0 = p[0];
	p1 = p[1];
	p2 = p[2];
	p3 = p[3];
	arg = (x - xr) / sqrt(2.0) / sr;
	return exp(-arg * arg) * sr / sqrt(2.0 * PI) *
		(p1 + p2 * (x + xr) + p3 * (x * x + x * xr + xr * xr + 2 * sr * sr)) +
		(1 + erf((arg))) / 2 *
		(p0 + x * (p1 + x * (p2 + p3 * x)) + (p2 + 3 * p3 * x) * sr * sr);
}

double aerogel(size_t n, double x, double* par)
{
	std::vector<double> xv(4);
	std::vector<double> yv(4);
	for (size_t i = 0; i < 4; i++)
		yv[i] = par[4 * n + i];
	double s1, s2;
	switch (n) {
	case 0:
		xv[0] = par[12];
		xv[3] = par[13] - 1.5 / 120.0;
		s1 = par[17];
		s2 = par[18];
		break;
	case 1:
		xv[0] = par[13] + 1.5 / 120.0;
		xv[3] = par[14];
		s1 = par[18];
		s2 = par[18];
		break;
	case 2:
		s1 = par[18];
		s2 = par[19];
		xv[0] = par[14];
		xv[3] = par[15];
		break;
	}
	xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
	xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
	double smin = xv[0];
	double smax = xv[3];
	std::vector<double> pv = m4(xv, yv);
	return p3g(x, xv[0], s1, pv, smin, smax) -
		p3g(x, xv[3], s2, pv, smin, smax);
}

double shifter(double x, double* par)
{
	double ds = 1.5 / 120.0;
	return par[16] * (erf((x - (par[13] - ds)) / sqrt(2.0) / par[18]) -
		erf((x - (par[13] + ds)) / sqrt(2.0) / par[18])) / 2;
}

double p3g(double x, double xr, double sr, std::vector<double> p, double smin, double smax)
{
	/*   double scal = PI/180;
	   double er=10.821691, der=2.40385*scal;
	   if (((fabs(smin-0)<10*scal)||
			(fabs(smin-120*scal)<10*scal)||
			(fabs(smin-240*scal)<10*scal)||
			(fabs(smin-360*scal)<10*scal))&&
		   ((xr-smin)<der))
		  return p3gl(x,xr,sr,er,smin,der,p);
	   if (((fabs(smax-0)<10*scal)||
			(fabs(smax-120*scal)<10*scal)||
			(fabs(smax-240*scal)<10*scal)||
			(fabs(smax-360*scal)<10*scal))&&
		   ((smax-xr)<der))
		  return p3gl(x,xr,sr,er,smax,-der,p);*/
	return p3g0(x, xr, sr, p);
}

double scphi(double* x, double* par)
{
	double f = 0;
	for (size_t i = 0; i < 3; i++)
		f += aerogel(i, x[0], par);
	f += shifter(x[0], par);
	return f;
}

void drawsomefit() {
	TCanvas* ci = new TCanvas("c1","c1",1500,1000);
	ci->SetFillColor(0);
	ci->SetFrameFillColor(0);
	ci->Divide(3,3);
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	double scal = PI / 180;
	int interobl[9] = {1,1,1,2,2,4,4,5,14};
	int intercount[9] = {1,3,6,1,9,3,4,7,1};
	double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles.root");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0;
	while ((key = (TKey*)next()) && counterr < 9) {
		counterr += 1;
		TString name(key->GetName());
		TString type = "zobl";
		if (name.BeginsWith(type) && name.Contains("sensor")) {
			//cin.get();
			int obl, counter;
			sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
			counter = intercount[counterr-1];
			obl = interobl[counterr-1];
			TProfile* h_amp_phi = (TProfile*)f->Get(Form("zobl%d,sensor%d", obl, counter));
			h_amp_phi->SetTitle(Form("counter%d,z_area%d;#phi,rad;amplitude,p.e.", counter, obl));
			gStyle->SetOptStat(11);
			ci->cd(counterr);
			h_amp_phi->Draw();
			std::vector<double> par(20);
			std::vector<double> parc(20);
			for (int o = 0; o < 4; o++)
				par[o] = 4.5;
			for (int o = 4; o < 8; o++)
				par[o] = 3.5;
			for (int o = 8; o < 12; o++)
				par[o] = 2.5;

			//pik________________
			double a = 0., b = 0.;
			int maxbin; double maxphi;
			for (int i = 50; i < 150; i++) {
				b = h_amp_phi->GetBinContent(i);
				if (b > a) {
					a = b;
					maxbin = i;
				}
			}
			maxphi = (h_amp_phi->GetXaxis()->GetBinLowEdge(maxbin) + h_amp_phi->GetXaxis()->GetBinUpEdge(maxbin)) / 2.;
			par[13] = maxphi;
			par[12] = par[13] - (13.2) * scal;
			par[14] = maxphi + (13.5) * scal;
			par[15] = maxphi + (23.5) * scal;

			if (counter == 7) {
				par[12] = maxphi - 13.7 * scal;
				par[14] = maxphi + 13.4 * scal;
				par[15] = maxphi + 23.5 * scal;
			}
			if (counter == 6) {
				par[15] = maxphi + 23.9 * scal;
			}
			if (counter == 8) {
				par[15] = maxphi + 23.9 * scal;
				par[12] = maxphi - 13.6 * scal;
			}
			if (counter == 3) {
				par[14] = maxphi + 13. * scal;
			}
			if (counter == 4)
				par[15] = maxphi + 25 * scal;
			if (counter == 9)
				par[12] = maxphi - 14 * scal;
			if (counter == 1 || counter == 2)
				par[15] = maxphi + 24 * scal;

			par[16] = a;
			par[17] = (0.6) * scal;
			par[18] = (0.4) * scal;
			par[19] = (0.6) * scal;

			parc[12] = par[12];
			parc[13] = par[13];
			parc[14] = par[14];
			parc[15] = par[15];

			//TCanvas* ci = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			ci->Update();
			TLine* line11 = new TLine(par[12], 0, par[12], gPad->GetUymax());
			TLine* line12 = new TLine(par[14], 0, par[14], gPad->GetUymax());
			TLine* line13 = new TLine(par[13], 0, par[13], gPad->GetUymax());
			TLine* line14 = new TLine(par[15], 0, par[15], gPad->GetUymax());
			line11->SetLineColor(kRed);
			line12->SetLineColor(kRed);
			line13->SetLineColor(kRed);
			line14->SetLineColor(kRed);
			line11->Draw();
			line12->Draw();
			line13->Draw();
			line14->Draw();
			ci->Update();
			cout << "ready" << endl;
			//cin.get();

			TF1* fn = new TF1("scphi", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
			//names
			for (size_t i = 0; i < 12; i++)
				fn->SetParName(i, Form("y%d", i));
			fn->SetParName(12, "x1");
			fn->SetParName(13, "x2");
			fn->SetParName(14, "x3");
			fn->SetParName(15, "x4");
			fn->SetParName(16, "ys");
			fn->SetParName(17, "#sigma_{1}");
			fn->SetParName(18, "#sigma_{s}");
			fn->SetParName(19, "#sigma_{2}");
			fn->SetParameters(par.data());
			fn->SetLineColor(kRed);
			fn->SetNpx(1000);

			//for(size_t i=0;i<12;i++)
		   //   fn->FixParameter(i,fn->GetParameter(i));
		   //fn->FixParameter(12, fn->GetParameter(12));
		   //fn->FixParameter(13,fn->GetParameter(13));
			//fn->FixParameter(14, fn->GetParameter(14));
			//fn->FixParameter(15, fn->GetParameter(15));
			//fn->FixParameter(16,fn->GetParameter(16));
			fn->FixParameter(17, fn->GetParameter(17));
			fn->FixParameter(18, fn->GetParameter(18));
			fn->FixParameter(19, fn->GetParameter(19));
			/*fn->SetParLimits(14, (counter-1)*41*scal+dphi+27.*scal, (counter-1)*41*scal+dphi+31*scal);
			fn->SetParLimits(12, par[14]-29*scal-0.01, par[14]-29*scal+0.05);
			fn->SetParLimits(13, par[14]-14*scal-0.02, par[14]-14*scal+0.02);
			fn->SetParLimits(15, par[14]+ 11*scal-0.02, par[14]+11*scal+0.02);*/
			//fn->SetParLimits(14, parc[14]-0.01, parc[14]+0.01);
			fn->SetParLimits(12, parc[12] - 0.02, parc[12] + 0.);
			fn->SetParLimits(13, parc[13] - 0.01, parc[13] + 0.01);
			fn->SetParLimits(14, parc[14] - 0.01, parc[14] + 0.01);
			fn->SetParLimits(15, parc[15] - 0.0, parc[15] + 0.03);
			fn->ReleaseParameter(16);
			fn->SetParLimits(16, a, a + 5. + 1. * obl);


			for (int o = 0; o < 12; o++) {
				fn->ReleaseParameter(o);
			}

			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			std::vector<double> xv(12);
			std::vector<double> yv(12);
			std::vector<int> binv(14);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;
			float yyy = 0.;
			float yyy1 = 0.;
			float yyy2 = 0.;
			float yyy4 = 0.;

			binv[12] = h_amp_phi->GetXaxis()->FindFixBin(xv[0]);
			yyy1 = h_amp_phi->GetBinContent(binv[12] + 3);
			binv[13] = h_amp_phi->GetXaxis()->FindFixBin(xv[11]);
			yyy2 = h_amp_phi->GetBinContent(binv[13] - 3);
			for (int o = 0; o < 3; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 0)
					fn->SetParLimits(o, yyy - 1., yyy1 + 0.1);
				else
					fn->SetParLimits(o, yyy - 0.2, yyy + 0.4);
			}
			fn->SetParLimits(3, yyy + 0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 3., yyy + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(3, yyy + 2., yyy + 10.);
				else
					fn->SetParLimits(3, yyy + 5., yyy + 20.);
			}
			for (int o = 5; o < 12; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 11)
					fn->SetParLimits(o, yyy - 1., yyy2 + 0.1);
				else
					fn->SetParLimits(o, yyy - 0.4, yyy + 0.4);
			}
			yyy4 = (h_amp_phi->GetBinContent(binv[5]) + h_amp_phi->GetBinContent(binv[5] - 1) + h_amp_phi->GetBinContent(binv[5] + 1)) / 3.;
			fn->SetParLimits(4, yyy4 + 0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 2., yyy4 + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 8.);
				else
					fn->SetParLimits(4, yyy4 + 5., yyy4 + 20.);
			}

			h_amp_phi->Fit("scphi");
			//fn->GetParameters(&par);
			ci->Update();
			cout << "ready" << endl;
			//cin.get();


			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;

			binv[12] = h_amp_phi->GetXaxis()->FindFixBin(xv[0]);
			yyy1 = h_amp_phi->GetBinContent(binv[12] + 5);
			binv[13] = h_amp_phi->GetXaxis()->FindFixBin(xv[11]);
			yyy2 = h_amp_phi->GetBinContent(binv[13] - 3);
			for (int o = 0; o < 3; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 0)
					fn->SetParLimits(o, yyy - 1., yyy1 + 0.1);
				else
					fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
			}
			fn->SetParLimits(3, yyy + 0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 3., yyy + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(3, yyy + 2., yyy + 8.);
				else
					fn->SetParLimits(3, yyy + 5., yyy + 20.);
			}
			for (int o = 5; o < 12; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 11)
					fn->SetParLimits(o, yyy - 1., yyy2 + 0.3);
				else
					fn->SetParLimits(o, yyy - 0.3, yyy + 0.15);
			}
			yyy4 = (h_amp_phi->GetBinContent(binv[5]) + h_amp_phi->GetBinContent(binv[5] - 1) + h_amp_phi->GetBinContent(binv[5] + 1)) / 3.;
			fn->SetParLimits(4, yyy4 + 0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 2., yyy4 + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 8.);
				else
					fn->SetParLimits(4, yyy4 + 5., yyy4 + 20.);
			}

			/*fn->ReleaseParameter(17);
			fn->ReleaseParameter(18);
			fn->ReleaseParameter(19);
			fn->SetParLimits(17, 2.0e-3, 1.4e-2);
			fn->SetParLimits(18, 1.0e-3, 1.2e-2);
			fn->SetParLimits(19, 2.0e-3, 2.0e-2);*/

			h_amp_phi->Fit("scphi");
			ci->Update();
			cout << "ready" << endl;
			//cin.get();

			fn->FixParameter(12, fn->GetParameter(12));
			fn->FixParameter(13, fn->GetParameter(13));
			fn->FixParameter(14, fn->GetParameter(14));
			fn->FixParameter(15, fn->GetParameter(15));


			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;

			binv[12] = h_amp_phi->GetXaxis()->FindFixBin(xv[0]);
			yyy1 = h_amp_phi->GetBinContent(binv[12] + 5);
			binv[13] = h_amp_phi->GetXaxis()->FindFixBin(xv[11]);
			yyy2 = h_amp_phi->GetBinContent(binv[13] - 5);
			for (int o = 0; o < 3; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 0)
					fn->SetParLimits(o, yyy - 1., yyy1 + 0.3);
				else
					fn->SetParLimits(o, yyy - 0.2, yyy + 0.15);
			}
			cout << yyy << endl;
			fn->SetParLimits(3, yyy + 0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 2., yyy + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(3, yyy + 2., yyy + 8.);
				else
					fn->SetParLimits(3, yyy + 5., yyy + 20.);
			}
			for (int o = 5; o < 12; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 11)
					fn->SetParLimits(o, yyy - 1., yyy2 + 0.3);
				else
					fn->SetParLimits(o, yyy - 0.25, yyy + 0.15);
			}
			yyy4 = (h_amp_phi->GetBinContent(binv[5]) + h_amp_phi->GetBinContent(binv[5] - 1) + h_amp_phi->GetBinContent(binv[5] + 1)) / 3.;
			fn->SetParLimits(4, yyy4 + 0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 2., yyy4 + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 8.);
				else
					fn->SetParLimits(4, yyy4 + 7., yyy4 + 20.);
			}

			fn->ReleaseParameter(17);
			fn->ReleaseParameter(18);
			fn->ReleaseParameter(19);
			fn->SetParLimits(17, 6.0e-3, 1.4e-2);
			fn->SetParLimits(18, 1.0e-3, 1.2e-2);
			fn->SetParLimits(19, 6.0e-3, 2.0e-2);



			h_amp_phi->Fit("scphi");
			ci->Update();


			/*TLine* line1 = new TLine(a1[counter - 1][obl - 1], 0, a1[counter - 1][obl - 1], gPad->GetUymax());
			TLine* line2 = new TLine(a2[counter - 1][obl - 1], 0, a2[counter - 1][obl - 1], gPad->GetUymax());
			TLine* line3 = new TLine(a3[counter - 1][obl - 1], 0, a3[counter - 1][obl - 1], gPad->GetUymax());
			TLine* line4 = new TLine(a4[counter - 1][obl - 1], 0, a4[counter - 1][obl - 1], gPad->GetUymax());*/
			TLine* line1 = new TLine(fn->GetParameter(12), 0, fn->GetParameter(12), gPad->GetUymax());
			TLine* line2 = new TLine(fn->GetParameter(13), 0, fn->GetParameter(13), gPad->GetUymax());
			TLine* line3 = new TLine(fn->GetParameter(14), 0, fn->GetParameter(14), gPad->GetUymax());
			TLine* line4 = new TLine(fn->GetParameter(15), 0, fn->GetParameter(15), gPad->GetUymax());
			line1->SetLineColor(kGreen);
			line2->SetLineColor(kGreen);
			line3->SetLineColor(kGreen);
			line4->SetLineColor(kGreen);
			line1->Draw();
			line2->Draw();
			line3->Draw();
			line4->Draw();
			ci->Update();
			cout << "ready" << endl;
			//cin.get();


			//drawing
			{
				//fn->Draw("same");
				for (size_t c = 0; c < 3; c++) {
					TF1* fna = new TF1(Form("aerogel%d", c), scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
					fna->SetParameters(fn->GetParameters());
					for (size_t i = 0; i < 12; i++)
						fna->SetParameter(i, 0);
					for (size_t i = c * 4; i < c * 4 + 4; i++)
						fna->SetParameter(i, fn->GetParameter(i));
					fna->SetParameter(16, 0);
					fna->SetLineColor(kBlue);
					fna->SetLineWidth(0);
					fna->SetNpx(1000);
					fna->Draw("same");
					//
				}
				//
				TF1* fns = new TF1("shifter", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
				fns->SetParameters(fn->GetParameters());
				for (size_t i = 0; i < 12; i++)
					fns->SetParameter(i, 0);
				fns->SetLineColor(kGreen);
				fns->SetLineWidth(0);
				fns->SetNpx(1000);
				fns->Draw("same");
				h_amp_phi->GetXaxis()->SetRangeUser(par[12] - 10 * scal, par[15] + 10 * scal);
				h_amp_phi->GetYaxis()->SetRangeUser(0, 1.1 * h_amp_phi->GetMaximum());
				//
				//c1->Write();
				//TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
				ci->Update();
				cin.get();
				//
				//size_t xxx;
				//cin>>xxx;
			}
		}
	}
}

void testforfit() {
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	std::vector<double> xv1(12);
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/model/profiles.root");
	TFile* f1 = new TFile("/work/users/kladov/snd2k/R006-003/model/fittedprofiles.root", "RECREATE");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0, ennn = 0;
	while ((key = (TKey*)next())) {
		counterr += 1; ennn += 1;
		if (ennn == 5)
			ennn = 1;
		TString name(key->GetName());
		TString type = "zobl";
		if (name.BeginsWith(type) && name.Contains("sensor")) {
			//cin.get();
			int obl, counter, enbeam;
			//sscanf(name.Data(), "zobl%d,sensor%d,en%d", &obl, &counter, &enbeam);
			sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
			TDirectory* dir = f1->mkdir(Form("zobl%d,sensor%d", obl, counter));
			//counter = 7;
			//obl = (counterr-1)/4+1;
			enbeam = ennn;
			TProfile* h_amp_phi = (TProfile*)f->Get(Form("zobl%d,sensor%d", obl, counter));

			h_amp_phi->Draw();
			std::vector<double> par(20);
			std::vector<double> parc(20);
			for (int o = 0; o < 4; o++)
				par[o] = 4.5;
			for (int o = 4; o < 8; o++)
				par[o] = 3.5;
			for (int o = 8; o < 12; o++)
				par[o] = 2.5;

			//pik________________
			double a = 0., b = 0.;
			int maxbin; double maxphi;
			for (int i = 50; i < 150; i++) {
				b = h_amp_phi->GetBinContent(i);
				if (b > a) {
					a = b;
					maxbin = i;
				}
			}
			maxphi = (h_amp_phi->GetXaxis()->GetBinLowEdge(maxbin) + h_amp_phi->GetXaxis()->GetBinUpEdge(maxbin)) / 2.;
			par[13] = maxphi;
			par[12] = par[13] - (14.0) * scal;
			par[14] = maxphi + (13.5) * scal;
			par[15] = maxphi + (24.0) * scal;

			if (counter == 7) {
				par[12] = maxphi - 13.8 * scal;
				par[14] = maxphi + 13.5 * scal;
				par[15] = maxphi + 23.8 * scal;
			}
			if (counter == 6) {
				par[15] = maxphi + 23.9 * scal;
				if (obl == 11)
					par[15] = maxphi + 24.3 * scal;
			}
			if (counter == 8) {
				par[15] = maxphi + 24.7 * scal;
				if (obl == 13 || obl == 14)
					par[15] = maxphi + 24.2 * scal;
				if (obl == 5)
					par[15] = maxphi + 25.0 * scal;
				par[12] = maxphi - 14.3 * scal;
				par[14] = maxphi + 13.1 * scal;
			}
			if (counter == 3) {
				par[14] = maxphi + 13. * scal;
			}
			if (counter == 4)
				par[15] = maxphi + 25 * scal;
			if (counter == 9) {
				par[12] = maxphi - 14.5 * scal;
				par[15] = maxphi + 23.8 * scal;
			}
			if (counter == 1 || counter == 2)
				par[15] = maxphi + 24.4 * scal;


			//________________
			for (int k = 0; k < 4; k++)
				//par[12+k] = truegran1[counter-1][k];
				par[12 + k] = gran[k][counter - 1][obl - 1];
			par[12] = par[12] + 0.003;
			par[15] = par[15] - 0.003;
			
			if (counter == 8 && obl == 3)
				par[15] = par[15] - 0.01;
			if (counter == 8 && obl == 2)
				par[12] = par[12] - 0.003;
			if (counter == 7 && obl == 3)
				par[15] = par[15] - 0.003;
			if (counter == 7 && obl == 4)
				par[15] = par[15] + 0.003;
			if (counter == 6 && (obl == 6 || obl == 7 || obl == 8 || obl == 9))
				par[13] = par[13] + 0.01;
			
			par[14] = truegran[counter - 1][2];
			par[16] = a;
			par[17] = (0.6) * scal;
			par[18] = (0.4) * scal;
			par[19] = (0.5) * scal;

			if (counter == 5) {
				par[15] = par[15] - 0.01;
				par[19] = (0.9) * scal;
			}
			if (counter == 2) {
				par[15] = par[15] - 0.005;
				if (obl == 10 || obl == 11) {
					par[15] = par[15] - 0.005+(obl-11.)/100.;
					par[19] = (0.3) * scal;
				}
			}
			if (counter == 9) {
				par[12] = par[12] - 0.001;
				if (obl == 14) {
					par[12] = par[12] + 0.002;
					par[15] = par[15] - 0.006;
				}
				if (obl == 10) {
					par[12] = par[12]-0.002;
					par[15] = par[15] - 0.002;
				}
				par[17] = (0.4) * scal;
				par[19] = (0.4) * scal;
			}
			if (counter == 7 && obl == 1) {
				par[14] = par[14] - 0.01;
				par[12] = par[12] - 0.004;
			}
			//________________

			//_________
			//for (int k = 0; k < 4; k++)
			//	par[12 + k] = truegran1[counter - 1][k];
			//par[12] = par[12] + 0.003;
			//par[15] = par[15] - 0.003;
			//__________

			parc[12] = par[12];
			parc[13] = par[13];
			parc[14] = par[14];
			parc[15] = par[15];

			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			TLine* line11 = new TLine(par[12], 0, par[12], gPad->GetUymax());
			TLine* line12 = new TLine(par[14], 0, par[14], gPad->GetUymax());
			TLine* line13 = new TLine(par[13], 0, par[13], gPad->GetUymax());
			TLine* line14 = new TLine(par[15], 0, par[15], gPad->GetUymax());
			line11->SetLineColor(kRed);
			line12->SetLineColor(kRed);
			line13->SetLineColor(kRed);
			line14->SetLineColor(kRed);
			line11->Draw();
			line12->Draw();
			line13->Draw();
			line14->Draw();
			c->Update();
			cout << "ready" << endl;
			//cin.get();

			TF1* fn = new TF1("scphi", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
			//names
			for (size_t i = 0; i < 12; i++)
				fn->SetParName(i, Form("y%d", i));
			fn->SetParName(12, "x1");
			fn->SetParName(13, "x2");
			fn->SetParName(14, "x3");
			fn->SetParName(15, "x4");
			fn->SetParName(16, "ys");
			fn->SetParName(17, "#sigma_{1}");
			fn->SetParName(18, "#sigma_{s}");
			fn->SetParName(19, "#sigma_{2}");
			fn->SetParameters(par.data());
			fn->SetLineColor(kRed);
			fn->SetNpx(1000);

			//for(size_t i=0;i<12;i++)
		   //   fn->FixParameter(i,fn->GetParameter(i));
			/*fn->FixParameter(12, fn->GetParameter(12));
			fn->FixParameter(13, fn->GetParameter(13));
			fn->FixParameter(14, fn->GetParameter(14));
			fn->FixParameter(15, fn->GetParameter(15));*/
			//fn->FixParameter(16,fn->GetParameter(16));
			fn->FixParameter(17, fn->GetParameter(17));
			fn->FixParameter(18, fn->GetParameter(18));
			fn->FixParameter(19, fn->GetParameter(19));
			fn->SetParLimits(12, parc[12] - 0.02, parc[12] + 0.015);
			fn->SetParLimits(13, parc[13] - 0.015, parc[13] + 0.015);
			fn->SetParLimits(14, parc[14] - 0.01, parc[14] + 0.01);
			fn->SetParLimits(15, parc[15] - 0.015, parc[15] + 0.02);
			fn->ReleaseParameter(16);
			if(counter == 8)
				fn->SetParLimits(16, a-4., a + 15.+1.*obl);
			if (counter == 6)
				fn->SetParLimits(16, a-1., a + 25. + 1. * obl);
			if (counter == 7) {
				fn->SetParLimits(16, a, a + 25. + 1. * obl);   
				//if (obl == 1)
				//	fn->SetParLimits(16, a-1., a + 25. + 1. * obl);
			}
			else
				fn->SetParLimits(16, a-2., a + 15.+1.*obl);

			
			for (int o = 0; o < 12; o++) {
				fn->ReleaseParameter(o);
			}

			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			std::vector<double> xv(12);
			std::vector<double> yv(12);
			std::vector<int> binv(14);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;
			float yyy = 0.;
			float yyy1 = 0.;
			float yyy2 = 0.;
			float yyy4 = 0.;

			binv[12] = h_amp_phi->GetXaxis()->FindFixBin(xv[0]);
			yyy1 = h_amp_phi->GetBinContent(binv[12] + 3);
			binv[13] = h_amp_phi->GetXaxis()->FindFixBin(xv[11]);
			yyy2 = h_amp_phi->GetBinContent(binv[13] - 3);
			for (int o = 0; o < 3; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o]-1) + h_amp_phi->GetBinContent(binv[o]+1))/3.;
				if (o == 0)
					fn->SetParLimits(o, yyy - 0.5, yyy1 + 0.1);
				else
					fn->SetParLimits(o, yyy - 0.2, yyy + 0.3);
			}
			fn->SetParLimits(3, yyy+0.5, yyy + 8.);
			if(counter == 9)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			if (obl == 14) {
				if(counter<5)
					fn->SetParLimits(3, yyy + 2., yyy + 10.);
				else if (counter == 9)
					fn->SetParLimits(3, yyy + 6., yyy + 20.);
				else if (counter == 8)
					fn->SetParLimits(3, yyy + 3., yyy + 20.);
				else
					fn->SetParLimits(3, yyy + 2., yyy + 20.);
			}
			for (int o = 5; o < 12; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o]-1) + h_amp_phi->GetBinContent(binv[o]+1))/3.;
				if (o == 11)
					fn->SetParLimits(o, yyy-0.5, yyy2 + 0.1);
				else
					fn->SetParLimits(o, yyy-0.3, yyy + 0.3);
			}
			yyy4 = (h_amp_phi->GetBinContent(binv[5]) + h_amp_phi->GetBinContent(binv[5] - 1) + h_amp_phi->GetBinContent(binv[5] + 1)) / 3.;
			fn->SetParLimits(4, yyy4+0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 8.);
				else if (counter == 9)
					fn->SetParLimits(4, yyy4 + 6., yyy4 + 20.);
				else if (counter == 8)
					fn->SetParLimits(4, yyy4 + 3., yyy4 + 20.);
				else
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 20.);
			}
			
			h_amp_phi->Fit("scphi","","", par[12] - 15 * scal, par[15] + 15 * scal);
			//fn->GetParameters(&par);
			c->Update();
			cout << "ready" << endl;
			//cin.get();


			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;

			binv[12] = h_amp_phi->GetXaxis()->FindFixBin(xv[0]);
			yyy1 = h_amp_phi->GetBinContent(binv[12]+3);
			binv[13] = h_amp_phi->GetXaxis()->FindFixBin(xv[11]);
			yyy2 = h_amp_phi->GetBinContent(binv[13]-3);
			for (int o = 0; o < 3; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o]-1) + h_amp_phi->GetBinContent(binv[o]+1))/3.;
				if(o == 0)
					fn->SetParLimits(o, yyy - 0.5, yyy1 +0.1);
				else
					fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
			}
			fn->SetParLimits(3, yyy+0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(3, yyy + 2., yyy + 8.);
				else if (counter == 9)
					fn->SetParLimits(3, yyy + 6., yyy + 20.);
				else if (counter == 8)
					fn->SetParLimits(3, yyy + 3., yyy + 20.);
				else
					fn->SetParLimits(3, yyy + 2., yyy + 20.);
			}
			for (int o = 5; o < 12; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o]-1) + h_amp_phi->GetBinContent(binv[o]+1))/3.;
				if (o == 11)
					fn->SetParLimits(o, yyy - 1., yyy2 + 0.3);
				else
					fn->SetParLimits(o, yyy - 0.3, yyy + 0.2);
			}
			yyy4 = (h_amp_phi->GetBinContent(binv[5]) + h_amp_phi->GetBinContent(binv[5] - 1) + h_amp_phi->GetBinContent(binv[5] + 1)) / 3.;
			fn->SetParLimits(4, yyy4+0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 8.);
				else if (counter == 9)
					fn->SetParLimits(4, yyy4 + 6., yyy4 + 20.);
				else if (counter == 8)
					fn->SetParLimits(4, yyy4 + 3., yyy4 + 20.);
				else
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 20.);
			}

			/*fn->ReleaseParameter(17);
			fn->ReleaseParameter(18);
			fn->ReleaseParameter(19);
			fn->SetParLimits(17, 2.0e-3, 1.4e-2);
			fn->SetParLimits(18, 1.0e-3, 1.2e-2);
			fn->SetParLimits(19, 2.0e-3, 2.0e-2);*/

			h_amp_phi->Fit("scphi", "","", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();
			cout << "ready" << endl;
			//cin.get();

			fn->FixParameter(12, fn->GetParameter(12));
			fn->FixParameter(13, fn->GetParameter(13));
			fn->FixParameter(14, fn->GetParameter(14));
			fn->FixParameter(15, fn->GetParameter(15));


			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;

			binv[12] = h_amp_phi->GetXaxis()->FindFixBin(xv[0]);
			yyy1 = h_amp_phi->GetBinContent(binv[12] + 3);
			binv[13] = h_amp_phi->GetXaxis()->FindFixBin(xv[11]);
			yyy2 = h_amp_phi->GetBinContent(binv[13] - 3);
			for (int o = 0; o < 3; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 0) {
					if (counter == 8 && obl > 4)
						fn->SetParLimits(o, yyy + 2., yyy1 + 0.3);
					else
						fn->SetParLimits(o, yyy - 0.5, yyy1 + 0.3);
				}
				else {
					if(counter == 7 && obl <3)
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.0);
					if (counter == 8)
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
					else
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
				}
			}
			cout << yyy << endl;
			fn->SetParLimits(3, yyy+0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			if (counter == 8)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			if (counter == 7 && obl == 1)
				fn->SetParLimits(3, yyy, yyy + 1.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(3, yyy + 2., yyy + 8.);
				else if (counter == 9)
					fn->SetParLimits(3, yyy + 6., yyy + 20.);
				else if (counter == 8)
					fn->SetParLimits(3, yyy + 5., yyy + 20.);
				else if (counter == 6)
					fn->SetParLimits(3, yyy + 6., yyy + 20.);
				else	
					fn->SetParLimits(3, yyy + 2., yyy + 20.);
			}
			for (int o = 5; o < 12; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 11)
					fn->SetParLimits(o, yyy - 0.5, yyy2 + 0.3);
				else {
					if (counter == 7 && obl < 3)
						fn->SetParLimits(o, yyy - 0.25, yyy + 0.1);
					if (counter == 8)
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
					else
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
				}
			}
			yyy4 = (h_amp_phi->GetBinContent(binv[5]) + h_amp_phi->GetBinContent(binv[5] - 1) + h_amp_phi->GetBinContent(binv[5] + 1)) / 3.;
			fn->SetParLimits(4, yyy4+0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);
			if (counter == 8)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 8.);
				else if (counter == 9)
					fn->SetParLimits(4, yyy4 + 6., yyy4 + 20.);
				else if (counter == 8)
					fn->SetParLimits(4, yyy4 + 6., yyy4 + 20.);
				else
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 20.);
			}

			fn->ReleaseParameter(17);
			fn->ReleaseParameter(18);
			fn->ReleaseParameter(19);
			fn->SetParLimits(17, 6.0e-3, 1.5e-2);
			fn->SetParLimits(18, 3.0e-3, 1.2e-2);
			fn->SetParLimits(19, 6.0e-3, 1.4e-2);
			if(counter == 6)
				fn->SetParLimits(18, 3.0e-3, 10.0e-3);
			if (counter == 5)
				fn->SetParLimits(19, 1.0e-2, 2.0e-2);
			if (counter == 8)
				fn->SetParLimits(18, 4.3e-3, 1.2e-2);
			if (counter == 7)
				fn->SetParLimits(18, 4.3e-3, 1.2e-2);

			

			h_amp_phi->Fit("scphi", "","", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();
			

			/*std::vector<double> xv2(4);
			std::vector<double> yv2(4);
			fn->GetParameters(&par[0]);
			xv1[0] = truegran1[counter - 1][0];
			xv1[3] = truegran1[counter - 1][1] - 1.5 / 120.0;
			xv1[1] = xv1[0] + (xv1[3] - xv1[0]) / 3;
			xv1[2] = xv1[3] - (xv1[3] - xv1[0]) / 3;
			xv1[4] = truegran1[counter - 1][1] + 1.5 / 120.0;
			xv1[7] = truegran1[counter - 1][2];
			xv1[5] = xv1[4] + (xv1[7] - xv1[4]) / 3;
			xv1[6] = xv1[7] - (xv1[7] - xv1[4]) / 3;
			xv1[8] = truegran1[counter - 1][2];
			xv1[11] = truegran1[counter - 1][3];
			xv1[9] = xv1[8] + (xv1[11] - xv1[8]) / 3;
			xv1[10] = xv1[11] - (xv1[11] - xv1[8]) / 3;

			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;
			for (int n = 0; n < 3; n++) {
				for (size_t i = 0; i < 4; i++) {
					yv2[i] = par[4 * n + i];
					xv2[i] =  xv[4 * n + i];
				}
				std::vector<double> pv1 = m4(xv2, yv2);
				//cout << pv1[0] << "	" << pv1[1] << "	" << pv1[2] << "	" << pv1[3] << endl;
				double x = xv1[1 + 4 * n];
				//cout << pv1[0] + pv1[1] * x + pv1[2] * pow(x, 2) + pv1[3] * pow(x, 3) << endl;
				for (int i = 0; i < 4; i++) {
					double x1 = xv1[i + 4 * n];
					yfm[counter - 1][obl - 1][4 * n + i] = pv1[0] + pv1[1] * x1 + pv1[2] * pow(x1, 2) + pv1[3] * pow(x1, 3);
					//cout << yfm[counter - 1][obl - 1][4 * n + i] << endl;
				}
			}
			//for (int i = 0; i < 12; i++)
			//	yfm[counter-1][obl-1][i] = scphi(&xv1[i],&par[0]);
			
			
			*/
			for (int ii = 0; ii < 4; ii++)
				gran1[ii][counter - 1][obl - 1] = fn->GetParameter(12 + ii);
			pik[counter-1][obl-1] = fn->GetParameter(16);

			/*a1[counter - 1][obl - 1] = fn->GetParameter(12) + 0.02;
			a4[counter - 1][obl - 1] = fn->GetParameter(15) - 0.02;
			a2[counter - 1][obl - 1] = fn->GetParameter(13) - 0.035;
			a3[counter - 1][obl - 1] = fn->GetParameter(13) + 0.035;
			if (obl == 14 || counter == 6) {
				a2[counter - 1][obl - 1] = fn->GetParameter(13) - 0.055;
				a3[counter - 1][obl - 1] = fn->GetParameter(13) + 0.055;
			}
			if (obl == 14 && counter < 4) {
				a2[counter - 1][obl - 1] = fn->GetParameter(13) - 0.045;
				a3[counter - 1][obl - 1] = fn->GetParameter(13) + 0.045;
			}*/

			/*TLine* line1 = new TLine(a1[counter - 1][obl - 1], 0, a1[counter - 1][obl - 1], gPad->GetUymax());
			TLine* line2 = new TLine(a2[counter - 1][obl - 1], 0, a2[counter - 1][obl - 1], gPad->GetUymax());
			TLine* line3 = new TLine(a3[counter - 1][obl - 1], 0, a3[counter - 1][obl - 1], gPad->GetUymax());
			TLine* line4 = new TLine(a4[counter - 1][obl - 1], 0, a4[counter - 1][obl - 1], gPad->GetUymax());*/
			TLine* line1 = new TLine(fn->GetParameter(12), 0, fn->GetParameter(12), gPad->GetUymax());
			TLine* line2 = new TLine(fn->GetParameter(13), 0, fn->GetParameter(13), gPad->GetUymax());
			TLine* line3 = new TLine(fn->GetParameter(14), 0, fn->GetParameter(14), gPad->GetUymax());
			TLine* line4 = new TLine(fn->GetParameter(15), 0, fn->GetParameter(15), gPad->GetUymax());
			line1->SetLineColor(kGreen);
			line2->SetLineColor(kGreen);
			line3->SetLineColor(kGreen);
			line4->SetLineColor(kGreen);
			line1->Draw();
			line2->Draw();
			line3->Draw();
			line4->Draw();
			c->Update();
			cout << "ready" << endl;
			//cin.get();
			
			
			//drawing
			{
				//fn->Draw("same");
				for (size_t c = 0; c < 3; c++) {
					TF1* fna = new TF1(Form("aerogel%d", c), scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
					fna->SetParameters(fn->GetParameters());
					for (size_t i = 0; i < 12; i++)
						fna->SetParameter(i, 0);
					for (size_t i = c * 4; i < c * 4 + 4; i++)
						fna->SetParameter(i, fn->GetParameter(i));
					fna->SetParameter(16, 0);
					fna->SetLineColor(kBlue);
					fna->SetLineWidth(0);
					fna->SetNpx(1000);
					fna->Draw("same");
					//
					dir->cd();
					sprintf(titlea, "zobl%d,sensor%d,aerogel%d", obl, counter, c + 1);
					fna->Write(titlea);
				}
				//
				TF1* fns = new TF1("shifter", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
				fns->SetParameters(fn->GetParameters());
				for (size_t i = 0; i < 12; i++)
					fns->SetParameter(i, 0);
				fns->SetLineColor(kGreen);
				fns->SetLineWidth(0);
				fns->SetNpx(1000);
				fns->Draw("same");
				h_amp_phi->GetXaxis()->SetRangeUser(par[12] - 15 * scal, par[15] + 15 * scal);
				h_amp_phi->GetYaxis()->SetRangeUser(0, 1.5 * h_amp_phi->GetMaximum());
				//
				//c1->Write();
				//TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
				c->Update();
				//
				sprintf(titles, "zobl%d,sensor%d,shifter", obl, counter);
				sprintf(titlef, "zobl%d,sensor%d,full", obl, counter);
				sprintf(titlep, "zobl%d,sensor%d,profile", obl, counter);
				sprintf(title, "zobl%d,sensor%d,canvas", obl, counter);
				fns->Write(titles);
				fn->Write(titlef);
				h_amp_phi->Write(titlep);
				c->Write(title);
				f1->cd();
				c->Write(title);
				//size_t xxx;
				//cin>>xxx;
			}
			//cin.get();
		}
	}
	f1->Write();
	f1->Close();
	/*ofstream fout;
	fout.open("/work/users/kladov/snd2k/R006-003/maindir/a1a2a3a4.cpp");
	{
		fout << "double a1[9][14] = { ";
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 14; j++) {
				fout << a1[i][j] << ", ";
			}
			fout << (char)92 << endl;
		}
		for (int j = 0; j < 13; j++) {
			fout << a1[8][j] << ", ";
		}
		fout << a1[8][13] << " };" << endl;
	}
	{
		fout << "double a2[9][14] = { ";
			for (int i = 0; i < 8; i++) {
				for (int j = 0; j < 14; j++) {
					fout << a2[i][j] << ", ";
				}
				fout << (char)92 << endl;
			}
		for (int j = 0; j < 13; j++) {
			fout << a2[8][j] << ", ";
		}
		fout << a2[8][13] << " };" << endl;
	}
	{
		fout << "double a3[9][14] = { ";
			for (int i = 0; i < 8; i++) {
				for (int j = 0; j < 14; j++) {
					fout << a3[i][j] << ", ";
				}
				fout << (char)92 << endl;
			}
		for (int j = 0; j < 13; j++) {
			fout << a3[8][j] << ", ";
		}
		fout << a3[8][13] << " };" << endl;
	}
	{
		fout << "double a4[9][14] = { ";
			for (int i = 0; i < 8; i++) {
				for (int j = 0; j < 14; j++) {
					fout << a4[i][j] << ", ";
				}
				fout << (char)92 << endl;
			}
		for (int j = 0; j < 13; j++) {
			fout << a4[8][j] << ", ";
		}
		fout << a4[8][13] << " };" << endl;
	}*/
	ofstream fout1;
	fout1.open("/work/users/kladov/snd2k/R006-003/model/gran1.cpp");
	{
		fout1 << "double gran1[4][9][14] = { ";
		for (int k = 0; k < 3; k++) {
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 14; j++) {
					fout1 << gran1[k][i][j] << ", ";
				}
				fout1 << (char)92 << endl;
			}
		}
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 14; j++) {
				fout1 << gran1[3][i][j] << ", ";
			}
			fout1 << (char)92 << endl;
		}
		for (int j = 0; j < 13; j++) {
			fout1 << gran1[3][8][j] << ", ";
		}
		fout1 << gran1[3][8][13] << " };" << endl;
		
		fout1 << "double pik[9][14] = { ";
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 14; j++) {
				fout1 << pik[i][j] << ", ";
			}
			fout1 << (char)92 << endl;
		}
		for (int j = 0; j < 13; j++) {
			fout1 << pik[8][j] << ", ";
		}
		fout1 << pik[8][13] << " };" << endl;
	}
	
	/*ofstream fout2;
	fout2.open("/work/users/kladov/snd2k/R006-003/maindir/yfm.cpp");
	{
		fout2 << "double yfm[9][14][12] = { ";
		for (int k = 0; k < 8; k++) {
			for (int i = 0; i < 14; i++) {
				for (int j = 0; j < 12; j++) {
					fout2 << yfm[k][i][j] << ", ";
				}
				fout2 << (char)92 << endl;
			}
		}
		for (int i = 0; i < 13; i++) {
			for (int j = 0; j < 12; j++) {
				fout2 << yfm[8][i][j] << ", ";
			}
			fout2 << (char)92 << endl;
		}
		for (int j = 0; j < 11; j++) {
			fout2 << yfm[8][13][j] << ", ";
		}
		fout2 << yfm[8][13][11] << " };" << endl;
	}*/


}

void testforfitmean() {
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	std::vector<double> xv1(12);
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/model/profilesmean.root");
	TFile* f1 = new TFile("/work/users/kladov/snd2k/R006-003/model/fittedprofilesmean.root", "RECREATE");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0;
	while ((key = (TKey*)next()) && counterr < 9) {
		counterr += 1;
		TString name(key->GetName());
		TString type = "counter";
		if (name.BeginsWith(type) && name.Contains("counter")) {
			//cin.get();
			int counter;
			//sscanf(name.Data(), "zobl%d,sensor%d,en%d", &obl, &counter, &enbeam);
			sscanf(name.Data(), "counter%d",&counter);
			TDirectory* dir = f1->mkdir(Form("counter%d",counter));
			//counter = 7;
			//obl = (counterr-1)/4+1;
			TProfile* h_amp_phi = (TProfile*)f->Get(Form("counter%d", counter));

			h_amp_phi->Draw();
			std::vector<double> par(20);
			std::vector<double> parc(20);
			for (int o = 0; o < 4; o++)
				par[o] = 4.5;
			for (int o = 4; o < 8; o++)
				par[o] = 3.5;
			for (int o = 8; o < 12; o++)
				par[o] = 2.5;

			//pik________________
			double a = 0., b = 0.;
			int maxbin; double maxphi;
			for (int i = 50; i < 150; i++) {
				b = h_amp_phi->GetBinContent(i);
				if (b > a) {
					a = b;
					maxbin = i;
				}
			}
			maxphi = (h_amp_phi->GetXaxis()->GetBinLowEdge(maxbin) + h_amp_phi->GetXaxis()->GetBinUpEdge(maxbin)) / 2.;
			par[13] = maxphi;


			//________________
			for (int k = 0; k < 4; k++)
				//par[12+k] = truegran1[counter-1][k];
				par[12 + k] = gran[k][counter - 1][5 - 1];
			par[12] = par[12] + 0.003;
			par[15] = par[15] - 0.003;


			par[14] = truegran[counter - 1][2];
			par[16] = a;
			par[17] = (0.6) * scal;
			par[18] = (0.4) * scal;
			par[19] = (0.5) * scal;

			if (counter == 5) {
				par[15] = par[15] - 0.03;
				par[12] = par[12] - 0.01;
				par[19] = (0.9) * scal;
			}
			if (counter == 7) {
				par[15] = par[15] - 0.01;
				par[12] = par[12] - 0.01;
				par[19] = (0.9) * scal;
			}
			if (counter == 8) {
				par[12] = par[12] - 0.01;
				par[13] = par[13] - 0.005;
				par[15] = par[15] - 0.008;
			}
			if (counter == 6) {
				par[12] = par[12] - 0.03;
				par[15] = par[15] - 0.02;
			}
			if (counter == 2) {
				par[15] = par[15] - 0.01;
			}
			if (counter == 1) {
				par[12] = par[12] - 0.01;
			}
			if (counter == 1) {
				par[15] = par[15] - 0.007;
			}
			if (counter == 4) {
				par[15] = par[15] - 0.015;
			}
			if (counter == 9) {
				par[12] = par[12] - 0.004;
				par[15] = par[15] - 0.005;
				par[17] = (0.4) * scal;
				par[19] = (0.4) * scal;
			}
			//________________

			//_________
			//for (int k = 0; k < 4; k++)
			//	par[12 + k] = truegran1[counter - 1][k];
			//par[12] = par[12] + 0.003;
			//par[15] = par[15] - 0.003;
			//__________

			parc[12] = par[12];
			parc[13] = par[13];
			parc[14] = par[14];
			parc[15] = par[15];

			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			TLine* line11 = new TLine(par[12], 0, par[12], gPad->GetUymax());
			TLine* line12 = new TLine(par[14], 0, par[14], gPad->GetUymax());
			TLine* line13 = new TLine(par[13], 0, par[13], gPad->GetUymax());
			TLine* line14 = new TLine(par[15], 0, par[15], gPad->GetUymax());
			line11->SetLineColor(kRed);
			line12->SetLineColor(kRed);
			line13->SetLineColor(kRed);
			line14->SetLineColor(kRed);
			line11->Draw();
			line12->Draw();
			line13->Draw();
			line14->Draw();
			c->Update();
			cout << "ready" << endl;
			//cin.get();

			TF1* fn = new TF1("scphi", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
			//names
			for (size_t i = 0; i < 12; i++)
				fn->SetParName(i, Form("y%d", i));
			fn->SetParName(12, "x1");
			fn->SetParName(13, "x2");
			fn->SetParName(14, "x3");
			fn->SetParName(15, "x4");
			fn->SetParName(16, "ys");
			fn->SetParName(17, "#sigma_{1}");
			fn->SetParName(18, "#sigma_{s}");
			fn->SetParName(19, "#sigma_{2}");
			fn->SetParameters(par.data());
			fn->SetLineColor(kRed);
			fn->SetNpx(1000);

			//for(size_t i=0;i<12;i++)
		   //   fn->FixParameter(i,fn->GetParameter(i));
			/*fn->FixParameter(12, fn->GetParameter(12));
			fn->FixParameter(13, fn->GetParameter(13));
			fn->FixParameter(14, fn->GetParameter(14));
			fn->FixParameter(15, fn->GetParameter(15));*/
			//fn->FixParameter(16,fn->GetParameter(16));
			fn->FixParameter(17, fn->GetParameter(17));
			fn->FixParameter(18, fn->GetParameter(18));
			fn->FixParameter(19, fn->GetParameter(19));
			fn->SetParLimits(12, parc[12] - 0.02, parc[12] + 0.015);
			fn->SetParLimits(13, parc[13] - 0.015, parc[13] + 0.015);
			if(counter == 5)
				fn->SetParLimits(13, parc[13]-0.02, parc[13] + 0.01);
			fn->SetParLimits(14, parc[14] - 0.03, parc[14] + 0.02);
			fn->SetParLimits(15, parc[15] - 0.015, parc[15] + 0.02);
			if(counter == 9)
				fn->SetParLimits(15, parc[15] - 0.005, parc[15] + 0.1);
			fn->ReleaseParameter(16);
			fn->SetParLimits(16, a-1, a + 15.);
			if(counter == 6)
				fn->SetParLimits(16, a-1, a + 15.);


			for (int o = 0; o < 12; o++) {
				fn->ReleaseParameter(o);
			}

			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			std::vector<double> xv(12);
			std::vector<double> yv(12);
			std::vector<int> binv(14);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;
			float yyy = 0.;
			float yyy1 = 0.;
			float yyy2 = 0.;
			float yyy4 = 0.;

			binv[12] = h_amp_phi->GetXaxis()->FindFixBin(xv[0]);
			yyy1 = h_amp_phi->GetBinContent(binv[12] + 3);
			binv[13] = h_amp_phi->GetXaxis()->FindFixBin(xv[11]);
			yyy2 = h_amp_phi->GetBinContent(binv[13] - 3);
			for (int o = 0; o < 3; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 0)
					fn->SetParLimits(o, yyy - 0.5, yyy1 + 0.1);
				else
					fn->SetParLimits(o, yyy - 0.2, yyy + 0.3);
			}
			fn->SetParLimits(3, yyy + 0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			for (int o = 5; o < 12; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 11)
					fn->SetParLimits(o, yyy - 0.5, yyy2 + 0.1);
				else
					fn->SetParLimits(o, yyy - 0.3, yyy + 0.3);
			}
			yyy4 = (h_amp_phi->GetBinContent(binv[5]) + h_amp_phi->GetBinContent(binv[5] - 1) + h_amp_phi->GetBinContent(binv[5] + 1)) / 3.;
			fn->SetParLimits(4, yyy4 + 0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);

			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			//fn->GetParameters(&par);
			c->Update();
			cout << "ready" << endl;
			//cin.get();


			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;

			binv[12] = h_amp_phi->GetXaxis()->FindFixBin(xv[0]);
			yyy1 = h_amp_phi->GetBinContent(binv[12] + 3);
			binv[13] = h_amp_phi->GetXaxis()->FindFixBin(xv[11]);
			yyy2 = h_amp_phi->GetBinContent(binv[13] - 3);
			for (int o = 0; o < 3; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 0)
					fn->SetParLimits(o, yyy - 0.5, yyy1 + 0.1);
				else
					fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
			}
			fn->SetParLimits(3, yyy + 0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			
			for (int o = 5; o < 12; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 11)
					fn->SetParLimits(o, yyy - 1., yyy2 + 0.3);
				else
					fn->SetParLimits(o, yyy - 0.3, yyy + 0.2);
			}
			yyy4 = (h_amp_phi->GetBinContent(binv[5]) + h_amp_phi->GetBinContent(binv[5] - 1) + h_amp_phi->GetBinContent(binv[5] + 1)) / 3.;
			fn->SetParLimits(4, yyy4 + 0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);

			/*fn->ReleaseParameter(17);
			fn->ReleaseParameter(18);
			fn->ReleaseParameter(19);
			fn->SetParLimits(17, 2.0e-3, 1.4e-2);
			fn->SetParLimits(18, 1.0e-3, 1.2e-2);
			fn->SetParLimits(19, 2.0e-3, 2.0e-2);*/

			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();
			cout << "ready" << endl;
			//cin.get();

			fn->FixParameter(12, fn->GetParameter(12));
			fn->FixParameter(13, fn->GetParameter(13));
			fn->FixParameter(14, fn->GetParameter(14));
			fn->FixParameter(15, fn->GetParameter(15));


			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;

			binv[12] = h_amp_phi->GetXaxis()->FindFixBin(xv[0]);
			yyy1 = h_amp_phi->GetBinContent(binv[12] + 3);
			binv[13] = h_amp_phi->GetXaxis()->FindFixBin(xv[11]);
			yyy2 = h_amp_phi->GetBinContent(binv[13] - 3);
			for (int o = 0; o < 3; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 0) {
					fn->SetParLimits(o, yyy - 0.5, yyy1 + 0.3);
				}
				else {
					if (counter == 8)
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
					else
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
				}
			}
			cout << yyy << endl;
			fn->SetParLimits(3, yyy + 0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			if (counter == 8)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			
			for (int o = 5; o < 12; o++) {
				binv[o] = h_amp_phi->GetXaxis()->FindFixBin(xv[o]);
				yyy = (h_amp_phi->GetBinContent(binv[o]) + h_amp_phi->GetBinContent(binv[o] - 1) + h_amp_phi->GetBinContent(binv[o] + 1)) / 3.;
				if (o == 11)
					fn->SetParLimits(o, yyy - 0.5, yyy2 + 0.3);
				else {
					if (counter == 7)
						fn->SetParLimits(o, yyy - 0.25, yyy + 0.1);
					if (counter == 8)
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
					else
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
				}
			}
			yyy4 = (h_amp_phi->GetBinContent(binv[5]) + h_amp_phi->GetBinContent(binv[5] - 1) + h_amp_phi->GetBinContent(binv[5] + 1)) / 3.;
			fn->SetParLimits(4, yyy4 + 0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);
			if (counter == 8)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);
			if (counter == 8)
				fn->SetParLimits(4, yyy4 + 0.0, yyy4 + 10.);

			fn->ReleaseParameter(17);
			fn->ReleaseParameter(18);
			fn->ReleaseParameter(19);
			fn->SetParLimits(17, 6.0e-3, 1.5e-2);
			fn->SetParLimits(18, 3.0e-3, 1.2e-2);
			fn->SetParLimits(19, 6.0e-3, 1.4e-2);
			if (counter == 6)
				fn->SetParLimits(18, 3.0e-3, 10.0e-3);
			if (counter == 5)
				fn->SetParLimits(19, 1.0e-2, 1.6e-2);
			if (counter == 8)
				fn->SetParLimits(18, 4.3e-3, 1.2e-2);
			if (counter == 7)
				fn->SetParLimits(18, 4.3e-3, 1.2e-2);



			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();


			/*std::vector<double> xv2(4);
			std::vector<double> yv2(4);
			fn->GetParameters(&par[0]);
			xv1[0] = truegran1[counter - 1][0];
			xv1[3] = truegran1[counter - 1][1] - 1.5 / 120.0;
			xv1[1] = xv1[0] + (xv1[3] - xv1[0]) / 3;
			xv1[2] = xv1[3] - (xv1[3] - xv1[0]) / 3;
			xv1[4] = truegran1[counter - 1][1] + 1.5 / 120.0;
			xv1[7] = truegran1[counter - 1][2];
			xv1[5] = xv1[4] + (xv1[7] - xv1[4]) / 3;
			xv1[6] = xv1[7] - (xv1[7] - xv1[4]) / 3;
			xv1[8] = truegran1[counter - 1][2];
			xv1[11] = truegran1[counter - 1][3];
			xv1[9] = xv1[8] + (xv1[11] - xv1[8]) / 3;
			xv1[10] = xv1[11] - (xv1[11] - xv1[8]) / 3;

			par[12] = fn->GetParameter(12);
			par[13] = fn->GetParameter(13);
			par[14] = fn->GetParameter(14);
			par[15] = fn->GetParameter(15);
			xv[0] = par[12];
			xv[3] = par[13] - 1.5 / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = par[13] + 1.5 / 120.0;
			xv[7] = par[14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = par[14];
			xv[11] = par[15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;
			for (int n = 0; n < 3; n++) {
				for (size_t i = 0; i < 4; i++) {
					yv2[i] = par[4 * n + i];
					xv2[i] =  xv[4 * n + i];
				}
				std::vector<double> pv1 = m4(xv2, yv2);
				//cout << pv1[0] << "	" << pv1[1] << "	" << pv1[2] << "	" << pv1[3] << endl;
				double x = xv1[1 + 4 * n];
				//cout << pv1[0] + pv1[1] * x + pv1[2] * pow(x, 2) + pv1[3] * pow(x, 3) << endl;
				for (int i = 0; i < 4; i++) {
					double x1 = xv1[i + 4 * n];
					yfm[counter - 1][obl - 1][4 * n + i] = pv1[0] + pv1[1] * x1 + pv1[2] * pow(x1, 2) + pv1[3] * pow(x1, 3);
					//cout << yfm[counter - 1][obl - 1][4 * n + i] << endl;
				}
			}
			//for (int i = 0; i < 12; i++)
			//	yfm[counter-1][obl-1][i] = scphi(&xv1[i],&par[0]);


			*/
			for (int ii = 0; ii < 4; ii++)
				truegran1m[counter - 1][ii] = fn->GetParameter(12 + ii);
			//pik[counter - 1][5 - 1] = fn->GetParameter(16);

			/*a1[counter - 1][obl - 1] = fn->GetParameter(12) + 0.02;
			a4[counter - 1][obl - 1] = fn->GetParameter(15) - 0.02;
			a2[counter - 1][obl - 1] = fn->GetParameter(13) - 0.035;
			a3[counter - 1][obl - 1] = fn->GetParameter(13) + 0.035;
			if (obl == 14 || counter == 6) {
				a2[counter - 1][obl - 1] = fn->GetParameter(13) - 0.055;
				a3[counter - 1][obl - 1] = fn->GetParameter(13) + 0.055;
			}
			if (obl == 14 && counter < 4) {
				a2[counter - 1][obl - 1] = fn->GetParameter(13) - 0.045;
				a3[counter - 1][obl - 1] = fn->GetParameter(13) + 0.045;
			}*/

			/*TLine* line1 = new TLine(a1[counter - 1][obl - 1], 0, a1[counter - 1][obl - 1], gPad->GetUymax());
			TLine* line2 = new TLine(a2[counter - 1][obl - 1], 0, a2[counter - 1][obl - 1], gPad->GetUymax());
			TLine* line3 = new TLine(a3[counter - 1][obl - 1], 0, a3[counter - 1][obl - 1], gPad->GetUymax());
			TLine* line4 = new TLine(a4[counter - 1][obl - 1], 0, a4[counter - 1][obl - 1], gPad->GetUymax());*/
			TLine* line1 = new TLine(fn->GetParameter(12), 0, fn->GetParameter(12), gPad->GetUymax());
			TLine* line2 = new TLine(fn->GetParameter(13), 0, fn->GetParameter(13), gPad->GetUymax());
			TLine* line3 = new TLine(fn->GetParameter(14), 0, fn->GetParameter(14), gPad->GetUymax());
			TLine* line4 = new TLine(fn->GetParameter(15), 0, fn->GetParameter(15), gPad->GetUymax());
			TLine* line21 = new TLine(truegran1[counter - 1][0], 0, truegran1[counter - 1][0], gPad->GetUymax());
			TLine* line22 = new TLine(truegran1[counter - 1][1], 0, truegran1[counter - 1][1], gPad->GetUymax());
			TLine* line23 = new TLine(truegran1[counter - 1][2], 0, truegran1[counter - 1][2], gPad->GetUymax());
			TLine* line24 = new TLine(truegran1[counter - 1][3], 0, truegran1[counter - 1][3], gPad->GetUymax());
			line1->SetLineColor(kGreen);
			line2->SetLineColor(kGreen);
			line3->SetLineColor(kGreen);
			line4->SetLineColor(kGreen);
			line21->SetLineColor(kBlue);
			line22->SetLineColor(kBlue);
			line23->SetLineColor(kBlue);
			line24->SetLineColor(kBlue);
			line1->Draw();
			line2->Draw();
			line3->Draw();
			line4->Draw();
			line21->Draw();
			line22->Draw();
			line23->Draw();
			line24->Draw();
			c->Update();
			cout << "ready" << endl;
			//cin.get();


			//drawing
			{
				//fn->Draw("same");
				for (size_t c = 0; c < 3; c++) {
					TF1* fna = new TF1(Form("aerogel%d", c), scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
					fna->SetParameters(fn->GetParameters());
					for (size_t i = 0; i < 12; i++)
						fna->SetParameter(i, 0);
					for (size_t i = c * 4; i < c * 4 + 4; i++)
						fna->SetParameter(i, fn->GetParameter(i));
					fna->SetParameter(16, 0);
					fna->SetLineColor(kBlue);
					fna->SetLineWidth(0);
					fna->SetNpx(1000);
					fna->Draw("same");
					//
					dir->cd();
					sprintf(titlea, "sensor%d,aerogel%d",counter, c + 1);
					fna->Write(titlea);
				}
				//
				TF1* fns = new TF1("shifter", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
				fns->SetParameters(fn->GetParameters());
				for (size_t i = 0; i < 12; i++)
					fns->SetParameter(i, 0);
				fns->SetLineColor(kGreen);
				fns->SetLineWidth(0);
				fns->SetNpx(1000);
				fns->Draw("same");
				h_amp_phi->GetXaxis()->SetRangeUser(par[12] - 15 * scal, par[15] + 15 * scal);
				h_amp_phi->GetYaxis()->SetRangeUser(0, 1.5 * h_amp_phi->GetMaximum());
				//
				//c1->Write();
				//TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
				c->Update();
				//
				sprintf(titles, "sensor%d,shifter", counter);
				sprintf(titlef, "sensor%d,full", counter);
				sprintf(titlep, "sensor%d,profile", counter);
				sprintf(title, "sensor%d,canvas", counter);
				fns->Write(titles);
				fn->Write(titlef);
				h_amp_phi->Write(titlep);
				c->Write(title);
				f1->cd();
				c->Write(title);
				//size_t xxx;
				//cin>>xxx;
			}
			//cin.get();
		}
	}
	f1->Write();
	f1->Close();
	ofstream fout1;
	fout1.open("/work/users/kladov/snd2k/R006-003/model/truegran1m.cpp");
	{
		fout1 << "double truegran1m[9][4] = { ";
		for (int k = 0; k < 8; k++) {
			for (int i = 0; i < 4; i++) {
				fout1 << truegran1m[k][i] << ", ";
			}
			fout1 << (char)92 << endl;
		}
		for (int i = 0; i < 3; i++) {
			fout1 << truegran1m[8][i] << ", ";
		}
		fout1 << truegran1m[8][3] << " };" << endl;
	}



}

void bordersplot() {
	char title[100];
	TFile* f1 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/bordersfromfit.root", "RECREATE");
	double scal1 = 180. / PI;
	double oblast[14];
	for (int i = 0; i < 14; i++)
		oblast[i] = (double)i + 1;
	double maximum = 0., minimum = 0.;
	cout << "rasbros" << endl << "counter	|	1gran	2gran	3gran	4gran" << endl;
	TGraph* graph[4][9];
	TGraph* grpik[9];
	TF1* fn = new TF1("fn", "[0]", 1, 14);
	fn->SetLineColor(kRed);
	for (int j = 0; j < 9; j++) {
		cout << "  " << j + 1 << "	|";
		for (int i = 0; i < 4; i++) {
			maximum = 0.; minimum = 10.;
			graph[i][j] = new TGraph(14, oblast, gran1[i][j]);
			graph[i][j]->SetLineColor(1);
			graph[i][j]->SetTitle(Form("counter%d,gran%d;zobl;phi",j+1,i+1));
			graph[i][j]->SetMarkerStyle(21);
			graph[i][j]->SetMarkerSize(1);
			graph[i][j]->SetMarkerColor(kBlack);
			//graph[i][j]->Draw("APL");
			sprintf(title, "sens%d,gran%d", j+1, i+1);
			graph[i][j]->Write(title);
			graph[i][j]->Fit(fn);
			cout << fn->GetParameter(0) << endl;
			truegran1[j][i] = fn->GetParameter(0);
			for (int k = 0; k < 14; k++) {
				if (gran1[i][j][k] > maximum) maximum = gran1[i][j][k];
				if (gran1[i][j][k] < minimum) minimum = gran1[i][j][k];
			}
			cout << "	" << Form("%.1f",(maximum - minimum)*scal1);
			//TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			//c->Update();
			//cin.get();
		}
		cout << endl;
		grpik[j] = new TGraph(14, oblast, pik[j]);
	}
	f1->Close();
	ofstream fout1;
	fout1.open("/work/users/kladov/snd2k/R006-003/maindir/truegran1.cpp");
	{
		fout1 << "double truegran1[9][4] = { ";
		for (int k = 0; k < 8; k++) {
			for (int i = 0; i < 4; i++) {
				fout1 << truegran1[k][i] << ", ";
			}
			fout1 << (char)92 << endl;
		}
		for (int i = 0; i < 3; i++) {
			fout1 << truegran1[8][i] << ", ";
		}
		fout1 << truegran1[8][3] << " };" << endl;
	}


	/*graph[0][0]->SetLineColor(1);
	graph[0][0]->GetYaxis()->SetRangeUser(0.,6.35);
	//graph[0][0]->GetHistogram()->SetMinimum(0.);
	//graph[0][0]->GetHistogram()->SetMinimum(0.);
	graph[0][0]->SetTitle("counterx;oblz;granitza");
	graph[0][0]->SetMarkerStyle(21);
	graph[0][0]->SetMarkerSize(1);
	graph[0][0]->SetMarkerColor(kBlack);
	graph[0][0]->Draw("APL");
	for (int j = 6; j < 7; j++) {
		for (int i = 0; i < 4; i++) {
			graph[i][j]->SetLineColor(i+1);
			graph[i][j]->SetMarkerStyle(21);
			graph[i][j]->SetMarkerSize(1);
			graph[i][j]->SetMarkerColor(kBlack);
			graph[i][j]->Draw("samePL");
		}
	}*/

	/*graph[2][1]->SetLineColor(1);
	graph[2][1]->GetYaxis()->SetRangeUser(0., 6.35);
	//graph[0][0]->GetHistogram()->SetMinimum(0.);
	graph[2][1]->SetTitle("counterx;oblz;granitza");
	graph[2][1]->SetMarkerStyle(21);
	graph[2][1]->SetMarkerSize(1);
	graph[2][1]->SetMarkerColor(kBlack);
	graph[2][1]->Draw("APL");*/


	/*grpik[0]->SetLineColor(1);
	grpik[0]->GetHistogram()->SetMaximum(60.);
	grpik[0]->GetHistogram()->SetMinimum(0.);
	grpik[0]->SetTitle("counterx;oblz;granitza");
	grpik[0]->SetMarkerStyle(21);
	grpik[0]->SetMarkerSize(1);
	grpik[0]->SetMarkerColor(kBlack);
	grpik[0]->Draw("APL");
	for (int j = 0; j < 9; j++) {
		grpik[j]->SetLineColor(1);
		grpik[j]->SetMarkerStyle(21);
		grpik[j]->SetMarkerSize(1);
		grpik[j]->SetMarkerColor(kBlack);
		grpik[j]->Draw("samePL");
	}*/
}

void modelfile() {
	/*double scal1 = 180. / PI;
	ofstream fout1;
	std::vector<double> xv1(12);
	int n = 14;
	fout1.open("/work/users/kladov/snd2k/R006-003/maindir/map.txt");
	for (int j = 0; j < 9; j++) {
		xv1[0] = truegran1[j][0];
		xv1[3] = truegran1[j][1] - 1.5 / 120.0;
		xv1[1] = xv1[0] + (xv1[3] - xv1[0]) / 3;
		xv1[2] = xv1[3] - (xv1[3] - xv1[0]) / 3;
		xv1[4] = truegran1[j][1] + 1.5 / 120.0;
		xv1[7] = truegran1[j][2];
		xv1[5] = xv1[4] + (xv1[7] - xv1[4]) / 3;
		xv1[6] = xv1[7] - (xv1[7] - xv1[4]) / 3;
		xv1[8] = truegran1[j][2];
		xv1[11] = truegran1[j][3];
		xv1[9] = xv1[8] + (xv1[11] - xv1[8]) / 3;
		xv1[10] = xv1[11] - (xv1[11] - xv1[8]) / 3;
		fout1 << Form("Counter %d", j+1) << endl << Form("Shifter %f %f %f",shift[j][0],shift[j][1],shift[j][2]) << endl << "AmplitudeCorrection 1" << endl << "Aerogel" << endl;
		for (int k = 0; k < 3; k++) {
			fout1 << n << " " << xv1[0+4*k]*scal1 << " " << xv1[1+4*k]*scal1 << " " << xv1[2+4*k]*scal1 << " " << xv1[3+4*k]*scal1 << endl;
			for (int i = 0; i < 14; i++) {
				fout1 << zobl[i] << " " << yfm[j][i][0+4*k] << " " << yfm[j][i][1+4*k] << " " << yfm[j][i][2+4*k] << " " << yfm[j][i][3+4*k] << endl;
			}
		}
	}*/
	
	double oblast[14];
	for (int i = 0; i < 14; i++)
		oblast[i] = (double)i + 1;
	double par[9][3];
	TGraph* grpik[9];
	for (int j = 0; j < 9; j++) {
		grpik[j] = new TGraph(14, oblast, pik[j]);
	}
	TF1* fn = new TF1("fn", "[0]*exp(-x/[2])+[1]*exp(x/[2])", 1, 14);
	//TF1* fn = new TF1("fn", "[0]+[1]*x", 1, 14);
	fn->SetParameters(1.,2.,3.);
	fn->SetParLimits(0,0.1,1000);
	fn->SetParLimits(1, 0.1, 1000);
	fn->SetParLimits(2, 0.01, 10);
	fn->SetLineColor(kRed);
	for (int j = 0; j < 9; j++) {
		grpik[j]->SetLineColor(1);
		grpik[j]->SetMarkerStyle(21);
		grpik[j]->SetMarkerSize(1);
		grpik[j]->SetMarkerColor(kBlack);
		grpik[j]->Draw("APL");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		grpik[j]->Fit("fn", "", "", 1,14);
		c->Update();
		cin.get();
		par[j][0] = fn->GetParameter(0);
		par[j][1] = fn->GetParameter(1);
		par[j][2] = fn->GetParameter(2);
	}
	/*ofstream fout;
	fout.open("/work/users/kladov/snd2k/R006-003/maindir/shiftersfm.cpp");
	fout << "double shift[9][3] = { ";
	for (int k = 0; k < 8; k++) {
		for (int i = 0; i < 3; i++) {
			fout << par[k][i] << ", ";
		}
		fout << (char)92 << endl;
	}
	for (int i = 0; i < 2; i++) {
		fout << par[8][i] << ", ";
	}
	fout << par[8][2] << " };" << endl;*/
}


void checkmod() {
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	std::vector<double> xv1(12);
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/model/profiles.root");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0, ennn = 0;
	while ((key = (TKey*)next()) && counterr < 9) {
		counterr += 1; ennn += 1;
		if (ennn == 5)
			ennn = 1;
		TString name(key->GetName());
		TString type = "zobl";
		if (name.BeginsWith(type) && name.Contains("sensor")) {
			//cin.get();
			int obl, counter, enbeam;
			//sscanf(name.Data(), "zobl%d,sensor%d,en%d", &obl, &counter, &enbeam);
			sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
			counter = counterr;
			obl = 2;
			enbeam = ennn;
			TProfile* h_amp_phi = (TProfile*)f->Get(Form("zobl%d,sensor%d", obl, counter));

			h_amp_phi->Draw();
			std::vector<double> par(20);
			for (int i = 0; i < 12; i++)
				par[i] = yfm[counter - 1][obl - 1][i];
			for (int i = 0; i < 4; i++)
				par[i+12] = truegran1[counter-1][i];
			par[16] = shift[counter - 1][0] * exp(-obl/shift[counter - 1][2]) + shift[counter - 1][1] * exp(obl/shift[counter - 1][2]);
			par[17] = (0.5) * scal;
			par[18] = (0.4) * scal;
			par[19] = (0.4) * scal;
			for (int i = 0; i < 20; i++)
				cout << Form("par[%d] = %f", i, par[i]) << endl;
			TF1* fn = new TF1("scphi", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 20);
			//names
			for (size_t i = 0; i < 12; i++)
				fn->SetParName(i, Form("y%d", i));
			fn->SetParName(12, "x1");
			fn->SetParName(13, "x2");
			fn->SetParName(14, "x3");
			fn->SetParName(15, "x4");
			fn->SetParName(16, "ys");
			fn->SetParName(17, "#sigma_{1}");
			fn->SetParName(18, "#sigma_{s}");
			fn->SetParName(19, "#sigma_{2}");
			fn->SetParameters(par.data());
			fn->SetLineColor(kRed);
			fn->SetNpx(1000);
			fn->Draw("same");
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update(); 
			TLine* line1 = new TLine(truegran1[counter - 1][0] + 0.03, 0, truegran1[counter - 1][0] + 0.03, 10);
			TLine* line2 = new TLine(truegran1[counter - 1][1] - 0.045, 0, truegran1[counter - 1][1] - 0.045, 10);
			TLine* line3 = new TLine(truegran1[counter - 1][1] + 0.045, 0, truegran1[counter - 1][1] + 0.045, 10);
			TLine* line4 = new TLine(truegran1[counter - 1][3] - 0.03, 0, truegran1[counter - 1][3] - 0.03, 10);
			line1->SetLineColor(kRed);
			line2->SetLineColor(kRed);
			line3->SetLineColor(kRed);
			line4->SetLineColor(kRed);
			line1->Draw();
			line2->Draw();
			line3->Draw();
			line4->Draw();
			c->Update();
			cin.get();
		}
	}
}

void go(){ 
	
	
//	copy1();
//	raspr();
//	testforfit();
//	testforfitmean();
//	bordersplot();
//	modelfile();
//	checkmod();
//	raspr();
//	raspramp();
//	findruns();
//	testt();	
//	calorimetr();
//	timespectra();
//	draweffspectre();
//	drawcalorimetr();
//	zraspr();
//	linesforz();
//	zsravn();
//	drawsomefit();
//	ampltime();
//	drawamptime();
//	compare();
//	checkmod();
//	effcomp1();
//	effcomp2();
//	effcalib();
//	compare1();
//	linesa();


	//pscp D:\programs\select_ach1mod.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/select_ach1mod.cpp
	//pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/select_ach1.cpp D:\programs\select_ach1.cpp
}

