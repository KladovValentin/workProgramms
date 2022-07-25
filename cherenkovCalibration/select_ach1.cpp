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
#include "/work/users/kladov/snd2k/R006-003/maindir/gran1.cpp"
#include "/work/users/kladov/snd2k/R006-003/maindir/goodruns.cpp"
#include "/work/users/kladov/snd2k/R006-003/maindir/truegran.cpp"
#include "/work/users/kladov/snd2k/R006-003/maindir/truegran1.cpp"
#include "/work/users/kladov/snd2k/R006-003/model/truegran1m.cpp"
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
//double truegran1[9][4];
//double yfm[9][14][12];
float zobl[15] = { -11.7013, -10.0367, -8.47203, -6.90737, -5.34273, -3.77807, -2.21343, -0.648766, 0.9159, 2.48053, 4.0452, 5.60983, 7.1745, 8.73913, 9.6};
float egran[5] = { 600., 780., 900., 955., 1010. };
float phi[2], phis[2], ach[40], theta[2], z0[2], schr[40], tchr[40], energy[40], beam, eton;
int region[40];

const char *cond = "(cosm == 0) && (nc == 2) && (nn==0) && (charge[0] == 1) && (charge[1] == 1) && (energy[0]/beam>=0.8) && (energy[1]/beam>=0.7) && (nch >= 9) && \
((theta[0] + theta[1] - 3.14159) <= 0.12) && ((theta[0] + theta[1] - 3.14159) >= -0.12) && (sin(theta[0]) >= 0.41) && \
((((phi[0] - phi[1]) - 3.14159 <= 0.1) && (((phi[0] - phi[1]) - 3.14159) >= -0.1)) || (((phi[0] - phi[1]) + 3.14159 <= 0.1) && (((phi[0] - phi[1]) + 3.14159) >= -0.1))) && \
((z0[0] - z0[1]) <= 1.0) && ((z0[0] - z0[1]) >= -1.0) && ((d0[0] - d0[1]) <= 0.5) &&  ((d0[0] - d0[1]) >= -0.5)";

const char* condkaons = "(cosm == 0) && (nc == 2) && (nn==0) && (charge[0] == 1) && (charge[1] == 1) && (energy[0]/beam<0.8) && (energy[1]/beam<0.7) && (nch >= 9) && \
(energy[0]/beam>0.2) && (energy[1]/beam>0.2)     &&\
((theta[0] + theta[1] - 3.14159) <= 0.12) && ((theta[0] + theta[1] - 3.14159) >= -0.12) && (sin(theta[0]) >= 0.41) && \
((((phi[0] - phi[1]) - 3.14159 <= 0.1) && (((phi[0] - phi[1]) - 3.14159) >= -0.1)) || (((phi[0] - phi[1]) + 3.14159 <= 0.1) && (((phi[0] - phi[1]) + 3.14159) >= -0.1))) && \
((z0[0] - z0[1]) <= 1.0) && ((z0[0] - z0[1]) >= -1.0) && ((d0[0] - d0[1]) <= 0.2) &&  ((d0[0] - d0[1]) >= -0.2)";
int mybin = 100;
int klmno[5], nch, run, nc, nn, cosm;
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
int wichcounter1[10];
float PI = 3.14159;

TString str = "a";
TString elem;

double p3g0(double x, double xr, double sr, std::vector<double> p);
double aerogel(size_t n, double x, double* par);
double shifter(double x, double* par);
double p3g(double x, double xr, double sr, std::vector<double> p, double smin, double smax);
double scphi(double* x, double* par);


size_t fdate(size_t run){
	TString cmd;
	cmd = Form("mysql -h sndfarm09.sndonline -u daqread -pdaqread rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"", run);
	cmd = cmd + " | " + Form("fgrep -e %d", run);
	cmd = cmd + " | " + Form("sed 's/%d//g'", run);
	cmd = Form("date -d \"$(%s)\" +\"%s%s%s\"", cmd.Data(), "%Y", "%m", "%d");
	return std::atoi(gSystem->GetFromPipe(cmd.Data()));
}
size_t ftime(size_t run){
	TString cmd;
	cmd = Form("mysql -h sndfarm09.sndonline -u daqread -pdaqread rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"", run);
	cmd = cmd + " | " + Form("fgrep -e %d", run);
	cmd = cmd + " | " + Form("sed 's/%d//g'", run);
	cmd = Form("date -d \"$(%s)\" +\"%s\"", cmd.Data(), "%s");
	return std::atoi(gSystem->GetFromPipe(cmd.Data())) - gStyle->GetTimeOffset() - 6 * 3600;
}

void testt() {
	int nnn = 0; int nnn1 = 0;
	int run = 0;
	for (int i = 2224; i < 5001; i++) {
		run = i + 25001;
		if (goodreff[i] == 1 && goodrnumb[i] == 1)
			nnn += 1;
		if ((/* ((run > 25950) && (run < 26150)) || ((run > 26175) && (run < 26275)) || ((run > 26300) && (run < 26825)) || */ ((run > 27225) && (run < 27275)) || ((run > 27325) && (run < 27725)) || ((run > 27850) && (run < 28425)) || ((run > 28450) && (run < 28575)) || ((run > 28625) && (run < 28850)) || ((run > 28875) && (run < 29475))))
			nnn1 += 1;
	}
	cout << nnn << endl;
	cout << nnn1 << endl;
}

void calorimetr() {
	TChain chain("t1");
	chain.Add("/work/users/konctbel/calibs/R007-001/output/ntuples/MHAD2017/*.root");
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	chain.SetBranchAddress("energy", &energy);
	chain.SetBranchAddress("beam", &beam);
	chain.SetBranchAddress("eton", &eton);
	chain.SetBranchAddress("nc", &nc);
	chain.SetBranchAddress("cosm", &cosm);
	chain.SetBranchAddress("nn", &nn);

	

	TH1* histo = new TH1F("histo", "energy0", 100, 0, 1.2);
	TH1* hist1 = new TH1F("hist1", "energy1", 100, 0, 1.2);
	TH1* hist2 = new TH1F("hist2", "eton", 100, 0, 1.2);
	
	
	int Count = 0;
	int count1 = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		if (beam == 1000 && nc == 2 && cosm == 0 && nn == 0) {
			histo->Fill(energy[0] / beam);
			hist2->Fill(eton);
			hist1->Fill(energy[1] / beam);
		}
		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d00k entries", count1) << endl;
			Count = 0;
		}
	}
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/calorimeter.root", "RECREATE");
	histo->Write("energy0");
	hist1->Write("energy1");
	hist2->Write("eton");
	MyFile->Close();

	/*histo->SetFillColor(kRed);
	histo->SetTitle("Energy deposition;Energy/beam");
	histo->Draw("LF2");
	hist1->SetFillColor(kBlue);
	hist1->SetTitle("Energy deposition;Energy/beam");
	hist1->Draw("LF2same");*/
	hist2->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	/*TLegend* legend = new TLegend(0.1, 0.83, 0.4, 0.98);
	//legend->SetHeader("The Legend Title", "C");
	legend->AddEntry("histo", "first particle", "f");
	legend->AddEntry("hist1", "second particle", "f");
	legend->Draw();*/
	//c->Update();
	TLine* line11 = new TLine(0.8, 0, 0.8, 1.05*histo->GetMaximum());
	TLine* line12 = new TLine(0.75, 0, 0.75, 1.05*histo->GetMaximum());
	line11->SetLineColor(kGreen);
	line12->SetLineColor(kGreen);
	//line11->Draw();
	line12->Draw();
	c->Update();
	cin.get();

}

void drawcalorimetr() {
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/calorimeter.root");

	
	TH1F* histo = (TH1F*)f->Get("energy0");
	TH1F* hist1 = (TH1F*)f->Get("energy1");
	TH1F* hist2 = (TH1F*)f->Get("eton");
	histo->SetTitle("Energy deposition;Energy/beam");
	hist2->SetTitle("Energy deposition;total deposition/2beam");
	hist2->Draw();
	/*histo->SetLineColor(kRed);
	histo->Draw();
	hist1->SetLineColor(kBlue);
	hist1->Draw("same");
	//TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);*/
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	c->Update();
	/*TLegend* legend = new TLegend(0.1, 0.83, 0.4, 0.98);
	//legend->SetHeader("The Legend Title", "C");
	legend->AddEntry("histo", "first particle", "f");
	legend->AddEntry("hist1", "second particle", "f");
	legend->Draw();
	c->Update();*/

	TLine* line11 = new TLine(0.8, 0, 0.8, 1.1 * histo->GetMaximum());
	TLine* line12 = new TLine(0.75, 0, 0.75, 1.05 * histo->GetMaximum());
	line11->SetLineColor(kGreen);
	line12->SetLineColor(kGreen);
	//line11->Draw();
	line12->Draw();
	c->Update();
}

void copy1() {
	vector<TString> ans;
	
	TString basedir("/work/users/konctbel/calibs/R007-001/output/ntuples/MHAD2017/*.root");
	TString dir("/work/users/kladov/snd2k/R006-003/selected1/");
	TString files = gSystem->GetFromPipe("ls " + basedir);

	Ssiz_t from = 0;
	
	while (files.Tokenize(elem, from, "\n")) {
		ans.push_back(elem);
	}
	unsigned int vector_size = ans.size();
	for (int i = 0; i < 23; i++) {
		str = ans[i].Copy();
		TString tok;
		Ssiz_t from1 = 1;
		vector<TString> tokens;
		TFile* filecomb = TFile::Open(str);
		TTree* originalTree = (TTree*)filecomb->Get("t1");

		while (str.Tokenize(tok, from1, "/")) {
			tokens.push_back(tok);
		}
		cout << tokens[8] << endl;
		TFile* ouput = TFile::Open(dir + "true1__" + tokens[8], "RECREATE");
		ouput->cd();
		TTree* selectedTree = originalTree->CopyTree(cond);
		selectedTree->Print();
		selectedTree->Write();
		filecomb->Close();
		ouput->Close();
	}
	for (int i = 24; i < vector_size; i++) {
		str = ans[i].Copy();
		TString tok;
		Ssiz_t from1 = 1;
		vector<TString> tokens;
		TFile* filecomb = TFile::Open(str);
		TTree* originalTree = (TTree*)filecomb->Get("t1");

		while (str.Tokenize(tok, from1, "/")) {
			tokens.push_back(tok);
		}
		cout << tokens[8] << endl;
		TFile* ouput = TFile::Open(dir + "true1__" + tokens[8], "RECREATE");
		ouput->cd();
		TTree* selectedTree = originalTree->CopyTree(cond);
		selectedTree->Print();
		selectedTree->Write();
		filecomb->Close();
		ouput->Close();
	}
}

void copy2() {
	vector<TString> ans;

	TString basedir("/work/users/konctbel/calibs/R007-001/output/ntuples/MHAD2017/*.root");
	TString dir("/work/users/kladov/snd2k/R006-003/selected1/");
	TString files = gSystem->GetFromPipe("ls " + basedir);

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
		TFile* filecomb = TFile::Open(str);
		TTree* originalTree = (TTree*)filecomb->Get("t1");

		while (str.Tokenize(tok, from1, "/")) {
			tokens.push_back(tok);
		}
		cout << tokens[8] << endl;
		TFile* ouput = TFile::Open(dir + "truekaons_" + tokens[8], "RECREATE");
		ouput->cd();
		TTree* selectedTree = originalTree->CopyTree(condkaons);
		selectedTree->Print();
		selectedTree->Write();
		filecomb->Close();
		ouput->Close();
	}
}

void dopeff() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/truekaons_**.root");
	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("schr", &schr);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("tchr", &tchr);
	chain.SetBranchAddress("energy", &energy);
	chain.SetBranchAddress("beam", &beam);
	int scount = 0;
	int scount1 = 0;

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;


	char name[20];
	char title[100];
	TH1* histdpek = new TH1I("histdpek", "do_por_eff_kaon_exp;eff", 100, -1., 2.);

	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	float ztr1 = 0.;
	int ph0 = 0;
	int ph1 = 0;
	int aka0 = 0;
	int aka1 = 0;
	/*&& \
		(energy[0] / beam < 0.9) && (energy[1] / beam < 0.9)*/
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.25 / tan(theta[0]) + z0[0];
		ztr1 = 12.25 / tan(theta[1]) + z0[1];
		if ((ztr > zobl[0]) && (ztr < zobl[14]) && (ztr1 > zobl[0]) && (ztr1 < zobl[14]) && (run > 27225) && (goodrnumb[run - 25001] == 1) && (goodreff[run - 25001] == 1) && (energy[0]/beam > 0.3) && (energy[1]/beam > 0.3) && (energy[0]/beam < 0.7) && (energy[1]/beam < 0.6)) {
			scount = 0;
			scount1 = 0;
			aka0 = 10;
			aka1 = 10;
			schr[nch] = 0.;
			for (int j = 0; j < 9; j++) {
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					if (phi[1] < 1.)
						phi[1] = phi[1] + 2. * PI;
					while (1 + scount < nch) { scount += 1; }
				}
				else
					while (schr[1 + scount] > 0.5) { scount += 1; }
				wichcounter1[j] = scount;

				if ((((phi[0] > truegran1[j][0] + 0.03) && (phi[0] < truegran1[j][1] - 0.045)) || ((phi[0] > truegran1[j][1] + 0.045) && (phi[0] < truegran1[j][3] - 0.03))) && (sin(theta[0])>0.9)) {
					ph0 = scount;
					/*if ((schr[scount] < 0.5))
						aka0 = 0;
					else
						aka0 = 1;*/
					aka0 = (int)schr[1 + scount1];
					//if (1+scount1>=nch)
					//	cout << schr[1+scount1] << endl;
				}
				if ((((phi[1] > truegran1[j][0] + 0.03) && (phi[1] < truegran1[j][1] - 0.045)) || ((phi[1] > truegran1[j][1] + 0.045) && (phi[1] < truegran1[j][3] - 0.03))) && (sin(theta[1])>0.9)) {
					ph1 = scount;
					/*if ((schr[scount] < 0.5))
						aka1 = 0;
					else
						aka1 = 1;*/
					aka1 = (int)schr[1 + scount1];
				}
				scount += 1;
				scount1 = scount;
			}
			//cout << aka0 << endl;
			if (aka0 == 0 && aka1 == 0)
				histdpek->Fill(0);
			else if ((aka0 == 0 && aka1 == 1) || (aka1 == 0 && aka0 == 1))
				histdpek->Fill(1);
		}
		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}
	histdpek->Draw();
	cout << "eff =	" << histdpek->GetMean(1) << endl;

}

void ampltime() {
	//chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_exp2017**.root");
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");
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

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;

	//hists for ach or eff
	TProfile* ha[10];
	char name[20];
	char title[100];
	for (int i = 0; i < 10; i++) {
		sprintf(name, "hprof%d", i);
		sprintf(title, "sens%d;ach;run", i);
		ha[i] = new TProfile(name,title, 50, 27000, 30000);
	}
	TProfile* he[10];
	for (int i = 0; i < 10; i++) {
		sprintf(name, "hprofeff%d", i);
		sprintf(title, "sens%d;eff;run", i);
		he[i] = new TProfile(name, title, 50, 27000, 30000);
	}


	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/phintegr_pedest.root");
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			TH1* hped = (TH1F*)f->Get(Form("ped_phintegr_zobl%d,sensor%d", j + 1, i + 1));
			ped[j][i] = hped->GetMean(1);
		}
	}
	TFile* f1 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedestz.root");
	for (int i = 0; i < 9; i++) {
		TH1* hped = (TH1F*)f1->Get(Form("ped_forz,sensor%d", i + 1));
		pedz[i] = hped->GetMean(1);
	}
	double lgran[4];
	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	float ztr1 = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12. / tan(theta[0]) + z0[0];
		ztr1 = 12. / tan(theta[1]) + z0[1];
		schr[nch] = 0;
		if ((run > 27225) && (goodrnumb[run - 25001] == 1) && (goodreff[run - 25001] == 1)) {
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9	; j++) {
				//wichcounter1[j] = scount1;
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					while (1 + scount < nch) { scount += 1; }
				}
				else
					while (schr[1 + scount] > 0.5) { scount += 1; }

				lgran[0] = truegran1[j][0];// +0.03;
				lgran[1] = truegran1[j][1];// -0.045;
				lgran[2] = truegran1[j][1];// +0.045;
				lgran[3] = truegran1[j][3];// -0.03;

				for (int k = 1; k < 13; k++) {
					if ((ztr > zobl[k]) && (ztr < zobl[k + 1]) && (region[0] == 1)) {
						if (((phi[0] > lgran[0]) && (phi[0] < lgran[1])) || ((phi[0] > lgran[2]) && (phi[0] < lgran[3]))) {
							if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
								ach1[scount] = ach[scount] - pedz[j];
								ha[0]->Fill(run, ach1[scount]);
								ha[j+1]->Fill(run,ach1[scount]);
								he[0]->Fill(run, schr[scount1+1]);
								he[j+1]->Fill(run, schr[scount1+1]);
							}
						}
					}
					if ((ztr1 > zobl[k]) && (ztr1 < zobl[k + 1]) && (region[1] == 1)) {
						if (((phi[1] > lgran[0]) && (phi[1] < lgran[1])) || ((phi[1] > lgran[2]) && (phi[1] < lgran[3]))) {
							if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
								ach1[scount] = ach[scount] - pedz[j];
								ha[0]->Fill(run, ach1[scount]);
								ha[j+1]->Fill(run, ach1[scount]);
								he[0]->Fill(run, schr[scount1 + 1]);
								he[j + 1]->Fill(run, schr[scount1 + 1]);
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

	/*int aaa = 0;
	vector<double> date, ampl, ddate, dampl;
	for (int i = 1; i < 26; i++) {
		if (ha->GetBinContent(i) > 2) {
			aaa = fdate(27000 - 60 + 120 * i) - 20170000;
			date.push_back((27000 - 60 + 120 * i));
			ampl.push_back(ha->GetBinContent(i));
			dampl.push_back(ha->GetBinError(i));
			ddate.push_back(0);
		}
	}
	int numbb = ampl.size();
	TGraphErrors* gr = new TGraphErrors(numbb, date.data(), ampl.data(), ddate.data(), dampl.data());
	//gr->GetXaxis()->SetTimeDisplay(1);
	//gr->GetXaxis()->SetTimeFormat("%y/%m/%d");
	gr->SetMarkerSize(0.5);
	gr->SetMarkerStyle(20);
	gr->Draw("AP");*/

	//ha->Draw();
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/achvsrun.root", "RECREATE");
	for (int i = 1; i < 10; i++) {
		sprintf(title, "avrsens%d", i);
		ha[i]->Write(title);
		sprintf(title, "evrsens%d", i);
		he[i]->Write(title);
	}
	ha[0]->Write("amplitudevsrun");
	he[0]->Write("effvsrun");
	MyFile->Close();

	/*TF1* fn = new TF1("scphi", "[0]+[1]*x", 27225, 29500);
	gr->Fit("scphi");
	cout << "0:" << fn->GetParameter(0) << endl;
	cout << "1:" << fn->GetParameter(1) << endl;
	double a = fn->GetParameter(1);
	double b = fn->GetParameter(0);
	cout << a*29500+b << endl;
	cout << "done!" << endl;*/
}

void drawamptime() {
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/achvsrun.root");
	TProfile* ha = (TProfile*)MyFile->Get("amplitudevsrun");
	TProfile* he = (TProfile*)MyFile->Get("effvsrun");
	int aaa = 0;
	vector<double> date, ampl, ddate, dampl, effic, deffic;
	for (int i = 1; i < 51; i++) {
		if (((i<20) && (ha->GetBinContent(i) >5.1)) || ((i >= 20) && (ha->GetBinContent(i) > 2.))) {
			ampl.push_back(1-exp(-ha->GetBinContent(i)/1.53));
			effic.push_back(he->GetBinContent(i));
			date.push_back(ftime(27000 - 30 + 60 * i));
			dampl.push_back(exp(-ha->GetBinContent(i)/1.5)*ha->GetBinError(i));
			deffic.push_back(he->GetBinError(i));
			ddate.push_back(0.);
		}
	}
	int numbb = ampl.size();

	//ha->Draw();
	TGraphErrors* gr = new TGraphErrors(numbb, &date[0], &ampl[0], &ddate[0], &dampl[0]);
	TGraphErrors* gr1 = new TGraphErrors(numbb, &date[0], &effic[0], &ddate[0], &deffic[0]);
	gr->SetTitle("amplitude vs run; run; mean amplitude-eff, p.e.");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr1->GetXaxis()->SetTimeDisplay(1);
	//gr->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d-%m}");
	gr->SetMarkerSize(0.5);
	gr->SetMarkerStyle(20);
	gr1->SetMarkerSize(0.5);
	gr1->SetMarkerStyle(21);
	gr->GetYaxis()->SetRangeUser(0.95,0.99);
	gr->SetLineColor(3);
	gr1->SetLineColor(4);
	gr->Draw("AP");
	gr1->Draw("same");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	c->SetGridy();
	c->Update();
	TF1* fn = new TF1("scphi", "[0]+[1]*x", ftime(27225), ftime(29300));
	/*gr->Fit("scphi","","", ftime(27225), ftime(29400));
	cout << "0:" << fn->GetParameter(0) << endl;
	cout << "1:" << fn->GetParameter(1) << endl;
	double a = fn->GetParameter(1);
	double b = fn->GetParameter(0);
	cout << a  << endl;
	cout << b << endl;*/
	//TF1* fn = new TF1("scphi", "[0]+[1]*x", 27225, 29500);
	/*TF1* fn = new TF1("scphi", "[0]+[1]*exp(-x/[2])", 27225, 29500);
	fn->SetLineColor(2);
	fn->SetParameters(4.,0.1,10000.);
	fn->SetParLimits(0,3.,10.);
	fn->SetParLimits(1,0.001,10.);
	fn->SetParLimits(2,5000.,100000.);
	gr->Fit("scphi");
	c->Update();
	cout << "0:" << fn->GetParameter(0) << endl;
	cout << "1:" << fn->GetParameter(1) << endl;
	cout << "2:" << fn->GetParameter(2) << endl;
	double a = fn->GetParameter(1);
	double b = fn->GetParameter(0);
	double d = fn->GetParameter(2);
	//cout << a * 27100 + b << endl;
	//cout << a * 29464 + b << endl;
	cout << b-a*exp(27100./d) << endl;
	cout << b-a*exp(29464./d) << endl;*/
	cout << "done!" << endl;
	MyFile->Close();
}

void timespectra() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");
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


	TH1* histo = new TH1F("histo", "time spectre;time,TDC chanels;", 300, 0, 300);
	TH1* histi = new TH1F("histi", "time spectre;time,TDC chanels;", 300, 0, 300);

	int Count = 0, count1 = 0;
	float ztr = 0.;
	float ztr1 = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.25 / tan(theta[0]) + z0[0];
		ztr = 12.25 / tan(theta[1]) + z0[1];
		schr[nch] = 0.;
		if ((run > 27225) && (goodrnumb[run - 25001] == 1) && (goodreff[run - 25001] == 1)) {
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

				if (scount - scount1 > 0.5 && (ztr > zobl[0]) && (ztr < zobl[14]) && (region[0] == 1)) {
					for (int a = scount; a >= scount; a--) {
						if ((((phi[0] > truegran1[j][0] + 0.03) && (phi[0] < truegran1[j][1] - 0.045)) || ((phi[0] > truegran1[j][1] + 0.045) && (phi[0] < truegran1[j][3] - 0.03))) && (schr[a] > 0.5) && j == 6)
							histo->Fill(tchr[a]);
					}
				}
				if (scount - scount1 > 0.5 && (ztr1 > zobl[0]) && (ztr1 < zobl[14]) && (region[1] == 1)) {
					for (int a = scount; a >= scount; a--) {
						if ((((phi[1] > truegran1[j][0] + 0.03) && (phi[1] < truegran1[j][1] - 0.045)) || ((phi[1] > truegran1[j][1] + 0.045) && (phi[1] < truegran1[j][3] - 0.03))) && (schr[a] > 0.5) && j == 6)
							histo->Fill(tchr[a]);
					}
				}
				if (scount - scount1 > 1.5 && (ztr > zobl[0]) && (ztr < zobl[14]) && (region[0] == 1)) {
					for (int a = scount-1; a >= scount-1; a--) {
						if ((((phi[0] > truegran1[j][0] + 0.03) && (phi[0] < truegran1[j][1] - 0.045)) || ((phi[0] > truegran1[j][1] + 0.045) && (phi[0] < truegran1[j][3] - 0.03))) && (schr[a] > 0.5) && j == 6)
							histi->Fill(tchr[a]);
					}
				}
				if (scount - scount1 > 1.5 && (ztr1 > zobl[0]) && (ztr1 < zobl[14]) && (region[1] == 1)) {
					for (int a = scount-1; a >= scount-1; a--) {
						if ((((phi[1] > truegran1[j][0] + 0.03) && (phi[1] < truegran1[j][1] - 0.045)) || ((phi[1] > truegran1[j][1] + 0.045) && (phi[1] < truegran1[j][3] - 0.03))) && (schr[a] > 0.5) && j == 6)
							histi->Fill(tchr[a]);
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
			Count = 0;
		}
	}
	histo->SetLineColor(2);
	histi->SetLineColor(13);
	histo->Draw();
	histi->Draw("same");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	TLine* line11 = new TLine(110, 0, 110, 1.1 * histo->GetMaximum());
	line11->SetLineColor(kGreen);
	line11->Draw();
	c->Update();

	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/timespectre.root", "RECREATE");
	histo->Write("TspectreMain");
	histi->Write("TspectreSec");
	//MyFile->Write();
	MyFile->Close();

	
	cin.get();
}

void drawtsp() {
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/timespectre.root");
	TH1* hma = (TH1F*)MyFile->Get("TspectreMain");
	TH1* hse = (TH1F*)MyFile->Get("TspectreSec");
	hma->SetLineColor(2);
	hse->SetLineColor(13);
	hma->Draw();
	hse->Draw("same");


}

void findruns() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");
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


	TH1* histo = new TH1I("runs","runs",5001,25000,30001);
	TProfile* hprof[9];
	char name[20];
	char title[100];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprof%d", i);
		sprintf(title, "sens%d;eff;run", i + 1);
		hprof[i] = new TProfile(name, title, 5001, 25000, 30001);
	}

	int Count = 0, count1 = 0;
	double lgran[4];
	float ztr = 0.;
	float ztr1 = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.25 / tan(theta[0]) + z0[0];
		ztr1 = 12.25 / tan(theta[1]) + z0[1];
		histo->Fill(run);
		schr[nch] = 0.;
		//if ((ztr > zobl[0]) && (ztr < zobl[14])) {
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

				lgran[0] = truegran1[j][0] + 0.03;
				lgran[1] = truegran1[j][1] - 0.045;
				lgran[2] = truegran1[j][1] + 0.045;
				lgran[3] = truegran1[j][3] - 0.03;

				if ((((phi[0] > lgran[0]) && (phi[0] < lgran[1])) || ((phi[0] > lgran[2]) && (phi[0] < lgran[3]))) && (ztr > zobl[0]) && (ztr < zobl[14]) && (region[0] == 1))
					hprof[j]->Fill(run, schr[1 + scount1]);
				if ((((phi[1] > lgran[0]) && (phi[1] < lgran[1])) || ((phi[1] > lgran[2]) && (phi[1] < lgran[3]))) && (ztr1 > zobl[0]) && (ztr1 < zobl[14]) && (region[1] == 1))
					hprof[j]->Fill(run, schr[1 + scount1]);
				scount += 1;
				scount1 = scount;
			}
		//}

		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			Count = 0;
		}
	}

	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/runsi.root", "RECREATE");
	TH1* histi[9];
	for(int j = 0; j<9;j++){
		sprintf(name, "sens%d_eff_spec", j+1);
		sprintf(title, "sens%d_eff_spec", j + 1);
		histi[j] = new TH1F(name, title, 100, 0., 1.);
		for (int i = 2; i <= 5001; i++) {
			if(histo->GetBinContent(i) > 500)
				histi[j]->Fill(hprof[j]->GetBinContent(i));
		}
		histi[j]->Write(title);
		sprintf(title, "sens%d_eff_run", j + 1);
		hprof[j]->Write(title);
	}
	histo->Write("runsnumb");
	//MyFile->Write();
	MyFile->Close();

	int ed = 1;
	int no = 0;
	int countt1 = 0; int countt2 = 0;
	ofstream fout;
	fout.open("/work/users/kladov/snd2k/R006-003/maindir/goodruns.cpp");
	{
		int k = 0;
		fout << "int goodreff[5001] = { ";
		for (int i = 2; i <= 5001; i++) {
			k = 0;
			for (int j = 0; j < 9; j++) {
				if (hprof[j]->GetBinContent(i) > 0.8)
					k += 1;
			}
			if (k == 9) {
				fout << ed << ",";
				countt1 += 1;
			}
			else if (k < 9)
				fout << no << ",";
			if ((i-1) % 100 == 0)
				fout << (char)92 << endl;
		}
		fout << no << " };" << endl;
	}
	{
		int k = 0;
		fout << "int goodrnumb[5001] = { ";
		for (int i = 2; i <= 5001; i++) {
			k = 0;
			if (histo->GetBinContent(i) > 2700)
				k = 1;
			if (k == 1) {
				fout << ed << ",";
				countt2 += 1;
			}
			else if (k == 0)
				fout << no << ",";
			if ((i - 1) % 100 == 0)
				fout << (char)92 << endl;
		}
		fout << no << " };" << endl;
	}
	cout << countt1 << endl << countt2 << endl;

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

void raspramp() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");

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

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;

	//hists for ach or eff
	/*TH1* ha[9][14][9];
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


	//info ab ped 9X14X9
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedest.root");
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				TH1* hped = (TH1F*)f->Get(Form("ped_sens%d,z%d,phi%d", i + 1, k + 1, j + 1));
				ped1[j][k][i] = hped->GetMean(1);
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
				lgran[0] = truegran1[j][0];// +0.03;
				lgran[1] = truegran1[j][1];// -0.045;
				lgran[2] = truegran1[j][1];// +0.045;
				lgran[3] = truegran1[j][3];// -0.03;
				
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
							if((schr[1+scount1] < 0.5) || ((schr[1+scount1] > 0.5) && (tchr[1+scount1] < 110)))
								hn[it][k][j]->Fill(schr[1 + scount1]); //it
							if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
								ach1[scount] = ach[scount] - ped1[it][k][j];
								ha[it][k][j]->Fill(ach1[scount]); //it
							}
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
							if ((schr[1 + scount1] < 0.5) || ((schr[1 + scount1] > 0.5) && (tchr[1 + scount1] < 110)))
								hn[it][k][j]->Fill(schr[1 + scount1]); //it
							if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
								ach1[scount] = ach[scount] - ped1[it][k][j];
								ha[it][k][j]->Fill(ach1[scount]); //it
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

	
	//effotach

	char titlea[200];
	char titlen[200];
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/random.root", "RECREATE");
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
	MyFile->Close();*/
	
	vector<float> errx[9], effect[9], erry[9], ampl[9];
	TFile* ff = new TFile("/work/users/kladov/snd2k/R006-003/maindir/random.root");
	int schit = 0;
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				TH1F* ha = (TH1F*)ff->Get(Form("ach_sens%d,z%d,phi%d", i + 1, k + 1, j + 1));
				TH1F* hn = (TH1F*)ff->Get(Form("eff_sens%d,z%d,phi%d", i + 1, k + 1, j + 1));
				if (hn->GetMean(1) > 0.2) {
					schit = schit + 1;
					ampl[i].push_back(ha->GetMean(1));
					effect[i].push_back(hn->GetMean(1));
					errx[i].push_back(1. * ha->GetMeanError(1));
					erry[i].push_back(1. * hn->GetMeanError(1));
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
	//TGraphErrors* gr[9];
	//TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/effotach.root", "RECREATE");
	//TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/effotach.root");
	TFile* MyFile1 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/effotachfitted.root", "RECREATE");
	for (int i = 0; i < 9; i++) {
		counter1[i] = (double)i + 1;
		nolliki[i] = 0.;
		numbrr = ampl[i].size();
		TGraphErrors* gr = new TGraphErrors(numbrr, &ampl[i][0], &effect[i][0], &errx[i][0], &erry[i][0]);
		gr->SetTitle(";ach,pe;eff");
		gr->GetYaxis()->SetRangeUser(0.8,1.01);
		gr->GetXaxis()->SetRangeUser(1., gr->GetXaxis()->GetXmax());
		gr->Draw("A*");
		TF1* func = new TF1("func", "[1]-[0]/(2.*[2]*x)*(exp(-x*(1.-[2])/[0])-exp(-x*(1.+[2])/[0]))", 2., gr->GetXaxis()->GetXmax());
		TF1* funcx = new TF1("funcx", "[1]-[0]/(2.*[2]*x)*(exp(-x*(1.-[2])/[0])-exp(-x*(1.+[2])/[0]))", 2., gr->GetXaxis()->GetXmax());
		func->SetLineColor(kRed);
		funcx->SetLineColor(kGreen);
		func->SetParameter(0, 1.);
		func->SetParameter(1, 1.);
		func->SetParameter(2, 0.2);
		func->SetParLimits(0, 1, 1);
		func->SetParLimits(1, 1, 1.0);
		func->SetParLimits(2, 0.01, 1.);
		funcx->SetParameter(0, 1.);
		funcx->FixParameter(1, 1.);
		funcx->SetParLimits(0, 1, 1);
		funcx->SetParameter(2, 0.2);
		funcx->SetParLimits(2, 0.01, 1.);
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


		TF1* func1 = new TF1("func1", "[1]-[0]/(2.*[2]*x)*(exp(-x*(1.-[2])/[0])-exp(-x*(1.+[2])/[0]))", 2., gr->GetXaxis()->GetXmax());
		//Double_t par[2];
		//func->GetParameters(par);
		func1->SetParameters(func->GetParameters());
		func1->SetLineColor(kRed);
		func1->Draw("same");

		TF1* funcx1 = new TF1("funcx1", "[1]-[0]/(2.*[2]*x)*(exp(-x*(1.-[2])/[0])-exp(-x*(1.+[2])/[0]))", 2., gr->GetXaxis()->GetXmax());
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
		legend->AddEntry("func1", "fit function a-exp(-#frac{x}{b})", "l");
		//legend->AddEntry("func0", "function 1-exp(-x)", "l");
		legend->AddEntry("funcx1", "fit function 1-exp(-#frac{x}{b})", "l");
		legend->Draw();
		c->Update();
		MyFile1->cd();
		c->Write(Form("effotachfit&d", i + 1));
		cin.get();
		//c->Close();
		//canv->Close();
	}
	/*MyFile1->Close();
	TGraphErrors* gr = new TGraphErrors(9, counter1, parametr1, nolliki, erre);
	gr->SetMarkerStyle(21);
	gr->SetMarkerSize(2.0);
	gr->SetTitle(";counter;max efficiency");
	gr->Draw("ap");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	cin.get();*/

	cout << "done!" << endl;

}

void drawnvsphi() {
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/fittedprofiles.root");
	//TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/nvsphi.root");
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/nvsphio19.root");
	

	for (int i = 0; i < 1; i++) {
		//TH1* h = (TH1*)MyFile->Get(Form("zobl%d,sensor%d", i + 1, 6));
		TH1* h = (TH1*)MyFile->Get(Form("sensor%d", 6));
		//TH1* h = (TH1*)MyFile->Get(Form("sensor%d", 6));
		h->GetYaxis()->SetRangeUser(0.,5.);
		h->SetLineColor(i + 1);
		h->DrawNormalized("", 2000);
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
	}

	TDirectory* d1 = (TDirectory*)f->Get(Form("zobl%d,sensor6", 3));
	TF1* tf1 = (TF1*)d1->Get(Form("zobl%d,sensor%d,full", 3, 6));
	tf1->SetLineColor(2);
	tf1->Draw("same");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	
}

void zraspr() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");

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

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;


	char name[20];
	char title[100];

	/*TH1* h[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprof%d", i + 1);
		sprintf(title, "ped,sensor%d;ach", i + 1);
		h[i] = new TH1F(name, title, 100, -1, 1);
	}*/
	
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

	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedestz.root");
	for (int i = 0; i < 9; i++) {
		TH1* hped = (TH1F*)f->Get(Form("ped_forz,sensor%d",i + 1));
		pedz[i] = hped->GetMean(1);
	}
	TH1* hr = new TH1F("histr", "recount right;z,cm", 100, 5., 15.);
	TH1* hl = new TH1F("histl", "recount left;z,cm", 100, -16., -6.);

	int kolotobr;
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
		schr[nch] = 0.;
		if ((run > 27225) && (goodrnumb[run - 25001] == 1) && (goodreff[run - 25001] == 1) ) {
			kolotobr += 1;
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

				lgran[0] = truegran1[j][0];// +0.03;
				lgran[1] = truegran1[j][1];// -0.045;
				lgran[2] = truegran1[j][1];// +0.045;
				lgran[3] = truegran1[j][3];// -0.03;

				/*if ((((phi[0] > lgran[0]) && (phi[0] < lgran[1])) || ((phi[0] > lgran[2]) && (phi[0] < lgran[3]))) && (region[0] == 1))
					h[j]->Fill(ach[scount1]);
				if ((((phi[1] > lgran[0]) && (phi[1] < lgran[1])) || ((phi[1] > lgran[2]) && (phi[1] < lgran[3]))) && (region[1] == 1))
					h[j]->Fill(ach[scount1]);*/
				

				if ((((phi[0] > lgran[0]) && (phi[0] < lgran[1])) || ((phi[0] > lgran[2]) && (phi[0] < lgran[3]))) && (region[0] == 1)) {
					if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
						ach1[scount] = (ach[scount] - pedz[j])*sin(theta[0]);
						hprof0[j]->Fill(ztr, ach1[scount]);
						hprof01->Fill(ztr, ach1[scount]);
						if (zin < 0) {
							if (((11.77 - fabs(zin)) * fabs(tan(theta[0])) < 3.0) && ((11.77 - fabs(zin)) * fabs(tan(theta[0])) >= 0.0)) {
								ach1[scount] = (ach[scount] - pedz[j]) * 3.0 * fabs(cos(theta[0])) / (11.77 - fabs(zin));
								hl->Fill(zin + 1.5 / tan(theta[0]));
								//hprof[j]->Fill((zin - 11.77) / 2., ach1[scount]);
								//hprof1->Fill((zin - 11.77) / 2., ach1[scount]);
								hprof[j]->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
								hprof1->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
							}
							else if ((11.77 - fabs(zin)) * fabs(tan(theta[0])) >= 3.0) {
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
								hr->Fill(zin + 1.5 / tan(theta[0]));
								//hprof[j]->Fill((zin + 10.8) / 2., ach1[scount]);
								//hprof1->Fill((zin + 10.8) / 2., ach1[scount]);
								hprof[j]->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
								hprof1->Fill(zin + 1.5 / tan(theta[0]), ach1[scount]);
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
					if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
						ach1[scount] = (ach[scount] - pedz[j])*sin(theta[1]);
						hprof0[j]->Fill(ztr1, ach1[scount]);
						hprof01->Fill(ztr1, ach1[scount]);
						if (zin1 < 0) {
							if (((11.77 - fabs(zin1)) * fabs(tan(theta[1])) < 3.0) && ((11.77 - fabs(zin1)) * fabs(tan(theta[1])) >= 0.0)) {
								ach1[scount] = (ach[scount] - pedz[j]) * 3.0 * fabs(cos(theta[1])) / (11.77 - fabs(zin1));
								hl->Fill(zin1 + 1.5 / tan(theta[1]));
								//hprof[j]->Fill((zin1 - 11.77) / 2., ach1[scount]);
								//hprof1->Fill((zin1 - 11.77) / 2., ach1[scount]);
								hprof[j]->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
								hprof1->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
							}
							else if ((11.77 - fabs(zin1)) * fabs(tan(theta[1])) >= 3.0) {
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
								hr->Fill(zin1 + 1.5 / tan(theta[1]));
								//hprof[j]->Fill((zin1 + 10.8) / 2., ach1[scount]);
								//hprof1->Fill((zin1 + 10.8) / 2., ach1[scount]);
								hprof[j]->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
								hprof1->Fill(zin1 + 1.5 / tan(theta[1]), ach1[scount]);
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
		}
		Count += 1;
		if (Count == 100000) {
			count1 += 1;
			cout << Form("obrabotano %d*100k entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}
	//cout << endl << kolotobr << endl << endl;
	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedestz.root", "RECREATE");
	for (int j = 0; j < 9; j++) {
		sprintf(title, "ped_forz,sensor%d", j + 1);
		h[j]->Write(title);
	}
	MyFile->Close();*/

	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/zprofiles.root", "RECREATE");
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

void linesforz() {

	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/zprofiles.root");

	TProfile* hprof = (TProfile*)f->Get("allsensorsmid");
	TProfile* hprof0 = (TProfile*)f->Get("allsensorsmid30");
	//TH1* hr = (TH1F*)f->Get("zinr");
	//TH1* hl = (TH1F*)f->Get("zinl");
	//hr->SetLineColor(kRed);
	//hl->SetLineColor(kRed);
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
	TLine* line11 = new TLine(-11.77, 0, -11.77, 7);
	TLine* line12 = new TLine(10.8, 0, 10.8, 7);
	line11->SetLineColor(kGreen);
	line12->SetLineColor(kGreen);
	line11->Draw();
	line12->Draw();
	//hr->DrawNormalized("same", 50);
	//hl->DrawNormalized("same", 50);
	hprof0->Draw("same");
	c->Update();
}

void raspr() {

	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("schr", &schr);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("tchr", &tchr);
	chain.SetBranchAddress("beam", &beam);
	chain.SetBranchAddress("region", &region);

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;

	//hist difference phi phis

	//TH1* histpps = new TH1F("sss", "differ;phi-phis", 100, -0.2, 0.2);
	//TProfile* profpps = new TProfile("sss", "differ;phi;phi-phis", 1200, 0., 6.3);




	//raspr po phi
	TProfile* hprof[14][9];
	char name[20];
	char title[100];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {	
			sprintf(name, "hprof%d", 9*j + i + 1);
			sprintf(title, "zobl%d,sensor%d;phi;ach", j + 1,i + 1);
			hprof[j][i] = new TProfile(name,title,  400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
		}
	}

	TProfile* hprofd[9][5];
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 5; j++) {
			sprintf(name, "hprofd%d", 5 * i + j + 1);
			sprintf(title, "data_range%d,sensor%d;phi;ach", j + 1, i + 1);
			hprofd[i][j] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
		}
	}
	
	TH1* hii[9][14];
	for (int a = 0; a < 14; a++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "histog%d", 9 * a + i + 1);
			sprintf(title, "#phi,sensor%d,obl%d;#phi", i + 1, a + 1);
			hii[i][a] = new TH1F(name, title, 400, truegran1[i][1] - 2., truegran1[i][1] + 2.);
		}
	}

	/*TProfile* hprof[9];
	char name[20];
	char title[100];
	for (int j = 0; j < 9; j++) {
		sprintf(name, "hprof%d", j + 1);
		sprintf(title, "counter%d;phi;ach", j + 1);
		hprof[j] = new TProfile(name, title, 200, (float)(j - 1.) * (2. * PI) / 9., (float)(j + 2.) * (2. * PI) / 9.);
	}*/

	/*TProfile* hprof[9];
	char name[20];
	char title[100];
	for (int j = 0; j < 9; j++) {
		sprintf(name, "hprof%d", j+1);
		sprintf(title, "sensor%d;phi;ach", j + 1);
		hprof[j] = new TProfile(name, title, 200, (float)(j - 1) * (2. * PI) / 9., (float)(j + 2) * (2. * PI) / 9.);
	}
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedest1.root");
	for (int j = 0; j < 9; j++) {
		TH1* hped = (TH1F*)f->Get(Form("ped_sens%d", j + 1));
		ped2[j] = hped->GetMean(1);
	}*/


	//hists for pedestal integr po phi
	/*TH1* h[14][9];
	char name[20];
	char title[100];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "hprof%d", 9*j + i + 1);
			sprintf(title, "ped,zobl%d,sensor%d;ach", j + 1, i + 1);
			h[j][i] = new TH1F(name, title, 100, -1, 1);
		}
	}*/

	/*TH1* h[14][9][4];
	char name[20];
	char title[100];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			for (int k = 0; k < 4; k++) {
				sprintf(name, "hprof%d", 4*9*j + 4*i + k + 1);
				sprintf(title, "ped,zobl%d,sensor%d,en%d;ach", j + 1, i + 1,k+1);
				h[j][i][k] = new TH1F(name, title, 100, -1, 1);
			}
		}
	}*/

	/*TH1* h[9];
	char name[20];
	char title[100];
	for (int j = 0; j < 14; j++) {
		sprintf(name, "hprof%d", j + 1);
		sprintf(title, "ped,sens%d;ach", j + 1);
		h[j] = new TH1F(name, title, 100, -1, 1);
	}*/

	//hists for pedestal or ach spectra 9X14X9
	/*TH1* h0[9][14][9];
	//char name[20];
	//char title[100];
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				sprintf(name, "hprof0%d", 14*9*j + 9*k + i);
				sprintf(title, "sens%d,z%d,phi%d;ach", i+1, k+1, j+1);
				h0[j][k][i] = new TH1F(name, title, 500, -5, 5);
			}
		}
	}*/

	//hists for ach or eff
	/*TH1* ha[9][14][9];
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
				ha[j][k][i] = new TH1F(namea, title, 500, -1, 49);
				hn[j][k][i] = new TH1F(namen, title, 100, -0.1, 2.);
			}
		}
	}*/

	//TH1* h1 = new TH1I("h1", "col otobrannih", 200, 25000, 30000);
	
	//info ab ped integrated phi
	TFile* f0 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/phintegr_pedest.root");
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			TH1* hped = (TH1F*)f0->Get(Form("ped_phintegr_zobl%d,sensor%d", j + 1, i + 1));
			ped[j][i] = hped->GetMean(1);
		}
	}

	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedestz.root");
	for (int i = 0; i < 9; i++) {
		TH1* hped1 = (TH1F*)f->Get(Form("ped_forz,sensor%d", i + 1));
		pedz[i] = hped1->GetMean(1);
	}

	//info ab ped 25X9
	/*TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedestvse.root");
	TProfile2D* hped = (TProfile2D*)f->Get("pedestvse");
	for (int k = 0; k < 25; k++) {
		for (int i = 0; i < 100; i++)
			pedestali1[k][i] = hped->GetBinContent(i+1,k+1);
	}*/

	//info ab ped 9X14X9
	/*TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedest.root");
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				TH1* hped = (TH1F*)f->Get(Form("ped_sens%d,z%d,phi%d", i + 1, k + 1, j + 1));
				ped1[j][k][i] = hped->GetMean(1);
			}
		} 
	}*/

	int Count = 0, count1 = 0;
	float ant;
	float ztr = 0.;
	float ztr1 = 0.;
	float zin = 0.;
	float zin1 = 0.;
	double lgran[4];
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//cout << beam << endl;
		ztr = 12.0 / tan(theta[0]) + z0[0];
		ztr1 = 12.0 / tan(theta[1]) + z0[1];
		zin = 10.5 / tan(theta[0]) + z0[0];
		zin1 = 10.5 / tan(theta[1]) + z0[1];
		schr[nch] = 0.;
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
				lgran[0] = truegran1[j][0] + 0.03;
				lgran[1] = truegran1[j][1] - 0.045;
				lgran[2] = truegran1[j][1] + 0.045;
				lgran[3] = truegran1[j][3] - 0.03;
				//if (((ztr > zobl[13]) && (ztr < zobl[14])))
				//profpps->Fill(phi[0], phi[0] - phis[0]);
				//histpps->Fill(phi[0]-phis[0]);

				//eff
				/*for (int k = 0; k < 14; k++) {
					for (int i = 0; i < 9; i++) {
						if ((ztr > zobl[k]) && (ztr < zobl[k+1]) && (phi[0] >= Goodphi[4 * j] + i * (Goodphi[4 * j + 3] - Goodphi[4 * j]) / 9.) && (phi[0] < Goodphi[4 * j] + (i + 1) * (Goodphi[4 * j + 3] - Goodphi[4 * j]) / 9.) && ((phi[0] < Goodphi[4 * j + 1]) || (phi[0] >= Goodphi[4 * j + 2]))) {
							if ((schr[1 + scount1] < 0.5) || ((schr[1 + scount1] > 0.5) && (tchr[1 + scount1] < 110)))
								hn[i][k][j]->Fill(run, schr[1 + scount1]);
							if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
								TH1* hped = (TH1F*)f->Get(Form("ped_sens%d,z%d,phi%d", j + 1, k + 1, i + 1));
								float ped = hped->GetMean(1);
								ach[scount] = ach[scount] - ped;
								hprofa[i][k][j]->Fill(run, ach[scount]);
							}
						}
					}
				}*/

				//ach*sin ot phi
				if ((phi[0] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[0] < (float)(j + 2) * (2. * PI) / 9.) && ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110)))) {
					if (zin < 0) {
						if (((11.77 - fabs(zin)) * fabs(tan(theta[0])) < 3.3) && ((11.77 - fabs(zin)) * fabs(tan(theta[0])) > 0.0)) {
							for (int i = 0; i < 14; i++) {
								//if (((zin - 11.77) / 2. >= zobl[i]) && ((zin - 11.77) / 2. < zobl[i + 1])) {
								if ((ztr >= zobl[i]) && (ztr < zobl[i + 1])) {
									ach1[scount] = (ach[scount] - ped[i][j]) * 3.3 * fabs(cos(theta[0])) / (11.77 - fabs(zin));
									hprof[i][j]->Fill(phi[0], ach1[scount]);
									hii[j][i]->Fill(phi[0]);
								}
							}
						}
						else if ((11.77 - fabs(zin)) * fabs(tan(theta[0])) >= 3.3) {
							for (int i = 0; i < 14; i++) {
								if ((ztr >= zobl[i]) && (ztr < zobl[i + 1])) {
									ach1[scount] = (ach[scount] - ped[i][j]) * sin(theta[0]);
									hprof[i][j]->Fill(phi[0], ach1[scount]);
									hii[j][i]->Fill(phi[0]);
								}
							}
						}
					}
					if (zin >= 0) {
						if (((10.8 - fabs(zin)) * fabs(tan(theta[0])) < 3.3) && ((10.8 - fabs(zin)) * fabs(tan(theta[0])) > 0.0)) {
							for (int i = 0; i < 14; i++) {
								//if (((zin + 10.8) / 2. > zobl[i]) && ((zin + 10.8) / 2. < zobl[i + 1])) {
								if ((ztr >= zobl[i]) && (ztr < zobl[i + 1])) {
									ach1[scount] = (ach[scount] - ped[i][j]) * 3.3 * fabs(cos(theta[0])) / (10.8 - fabs(zin));
									hprof[i][j]->Fill(phi[0], ach1[scount]);
									hii[j][i]->Fill(phi[0]);
								}
							}
						}
						else if ((10.8 - fabs(zin)) * fabs(tan(theta[0])) >= 3.3) {
							for (int i = 0; i < 14; i++) {
								if ((ztr >= zobl[i]) && (ztr < zobl[i + 1])) {
									ach1[scount] = (ach[scount] - ped[i][j]) * sin(theta[0]);
									hprof[i][j]->Fill(phi[0], ach1[scount]);
									hii[j][i]->Fill(phi[0]);
								}
							}
						}
					}

				}
				if ((phi[1] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[1] < (float)(j + 2) * (2. * PI) / 9.) && ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110)))) {
					if (zin1 < 0) {
						if (((11.77 - fabs(zin1)) * fabs(tan(theta[1])) < 3.3) && ((11.77 - fabs(zin1)) * fabs(tan(theta[1])) >= 0.0)) {
							for (int i = 0; i < 14; i++) {
								//if (((zin1 - 11.77) / 2. > zobl[i]) && ((zin1 - 11.77) / 2. < zobl[i + 1])) {
								if ((ztr1 >= zobl[i]) && (ztr1 < zobl[i + 1])) {
									ach1[scount] = (ach[scount] - ped[i][j]) * 3.3 * fabs(cos(theta[1])) / (11.77 - fabs(zin1));
									hprof[i][j]->Fill(phi[1], ach1[scount]);
									hii[j][i]->Fill(phi[1]);
								}
							}
						}
						else if ((11.77 - fabs(zin1)) * fabs(tan(theta[1])) >= 3.3) {
							for (int i = 0; i < 14; i++) {
								if ((ztr1 >= zobl[i]) && (ztr1 < zobl[i + 1])) {
									ach1[scount] = (ach[scount] - ped[i][j]) * sin(theta[1]);
									hprof[i][j]->Fill(phi[1], ach1[scount]);
									hii[j][i]->Fill(phi[1]);
								}
							}
						}
					}
					if (zin1 >= 0) {
						if (((10.8 - fabs(zin1)) * fabs(tan(theta[1])) < 3.3) && ((10.8 - fabs(zin1)) * fabs(tan(theta[1])) >= 0.0)) {
							for (int i = 0; i < 14; i++) {
								//if (((zin1 + 10.8) / 2. > zobl[i]) && ((zin1 + 10.8) / 2. < zobl[i + 1])) {
								if ((ztr1 >= zobl[i]) && (ztr1 < zobl[i + 1])) {
									ach1[scount] = (ach[scount] - ped[i][j]) * 3.3 * fabs(cos(theta[1])) / (10.8 - fabs(zin1));
									hprof[i][j]->Fill(phi[1], ach1[scount]);
									hii[j][i]->Fill(phi[1]);
								}
							}
						}
						else if ((10.8 - fabs(zin1)) * fabs(tan(theta[1])) >= 3.3) {
							for (int i = 0; i < 14; i++) {
								if ((ztr1 >= zobl[i]) && (ztr1 < zobl[i + 1])) {
									ach1[scount] = (ach[scount] - ped[i][j]) * sin(theta[1]);
									hprof[i][j]->Fill(phi[1], ach1[scount]);
									hii[j][i]->Fill(phi[1]);
								}
							}
						}
					}
				}
				/*if ((phi[0] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[0] < (float)(j + 2) * (2. * PI) / 9.)) {
					for (int i = 0; i < 14; i++) {
						if ((ztr > zobl[i]) && (ztr < zobl[i + 1])) {
							if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
								ach1[scount] = ach[scount] * sin(theta[0]) - ped[i][j];
								hprof[j]->Fill(phi[0], ach1[scount]);
							}
						}
					}
				}*/

				/*if ((phi[0] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[0] < (float)(j + 2) * (2. * PI) / 9.) && ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110)))) {
					if ((ztr > zobl[1]) && (ztr < zobl[14])) {
						ach1[scount] = (ach[scount] - pedz[j]) * sin(theta[0]);
						int i = (run - 27225) / 450;
						hprofd[j][i]->Fill(phi[0], ach1[scount]);
					}
				}
				if ((phi[1] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[1] < (float)(j + 2) * (2. * PI) / 9.) && ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110)))) {
					if ((ztr1 > zobl[1]) && (ztr1 < zobl[14])) {
						ach1[scount] = (ach[scount] - pedz[j]) * sin(theta[1]);
						int i = (run - 27225) / 450;
						hprofd[j][i]->Fill(phi[1], ach1[scount]);
					}
				}*/

				/*if ((phi[0] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[0] < (float)(j + 2) * (2. * PI) / 9.)) {
					for (int i = 0; i < 14; i++) {
						for (int k = 0; k < 4; k++) {
							if ((ztr > zobl[i]) && (ztr < zobl[i + 1]) && (beam > egran[k]) && (beam < egran[k+1])) {
								if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
									//if (((phi[0] >= Goodphi[4 * j]) && phi[0] < Goodphi[4 * j + 1]) || ((phi[0] >= Goodphi[4 * j + 2]) && phi[0] < Goodphi[4 * j + 3])) {
									ach1[scount] = ach[scount] - pede[i][j][k];
									hprof[i][j][k]->Fill(phi[0], ach1[scount] * sin(theta[0]));
									//}
								}
							}
						}
					}
				}*/

				//ped int po phi
				/*if ((phi[0] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[0] < (float)(j + 2) * (2. * PI) / 9.)) {
					if (zin < 0) {
						if (((11.77 - fabs(zin)) * fabs(tan(theta[0])) < 3.0) && ((11.77 - fabs(zin)) * fabs(tan(theta[0])) >= 0.0)) {
							for (int i = 0; i < 14; i++) {
								//if (((zin - 11.77) / 2. > zobl[i]) && ((zin - 11.77) / 2. < zobl[i + 1])) {
								if ((ztr > zobl[i]) && (ztr < zobl[i + 1])) {
									h[i][j]->Fill(ach[scount1]);
								}
							}
						}
						else if ((11.77 - fabs(zin)) * fabs(tan(theta[0])) >= 3.0) {
							for (int i = 0; i < 14; i++) {
								if ((ztr > zobl[i]) && (ztr < zobl[i + 1])) {
									h[i][j]->Fill(ach[scount1]);
								}
							}
						}
					}
					if (zin >= 0) {
						if (((10.8 - fabs(zin)) * fabs(tan(theta[0])) < 3.0) && ((10.8 - fabs(zin)) * fabs(tan(theta[0])) >= 0.0)) {
							for (int i = 0; i < 14; i++) {
								//if (((zin + 10.8) / 2. > zobl[i]) && ((zin + 10.8) / 2. < zobl[i + 1])) {
								if ((ztr > zobl[i]) && (ztr < zobl[i + 1])) {
									h[i][j]->Fill(ach[scount1]);
								}
							}
						}
						else if ((10.8 - fabs(zin)) * fabs(tan(theta[0])) >= 3.0) {
							for (int i = 0; i < 14; i++) {
								if ((ztr > zobl[i]) && (ztr < zobl[i + 1])) {
									h[i][j]->Fill(ach[scount1]);
								}
							}
						}
					}
				}
				if ((phi[1] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[1] < (float)(j + 2) * (2. * PI) / 9.)) {
					if (zin1 < 0) {
						if (((11.77 - fabs(zin1)) * fabs(tan(theta[1])) < 3.0) && ((11.77 - fabs(zin1)) * fabs(tan(theta[1])) >= 0.0)) {
							for (int i = 0; i < 14; i++) {
								//if (((zin1 - 11.77) / 2. > zobl[i]) && ((zin1 - 11.77) / 2. < zobl[i + 1])) {
								if ((ztr1 > zobl[i]) && (ztr1 < zobl[i + 1])) {
									h[i][j]->Fill(ach[scount1]);
								}
							}
						}
						else if ((11.77 - fabs(zin1)) * fabs(tan(theta[1])) >= 3.0) {
							for (int i = 0; i < 14; i++) {
								if ((ztr1 > zobl[i]) && (ztr1 < zobl[i + 1])) {
									h[i][j]->Fill(ach[scount1]);
								}
							}
						}
					}
					if (zin1 >= 0) {
						if (((10.8 - fabs(zin1)) * fabs(tan(theta[1])) < 3.0) && ((10.8 - fabs(zin1)) * fabs(tan(theta[1])) >= 0.0)) {
							for (int i = 0; i < 14; i++) {
								//if (((zin1 + 10.8) / 2. > zobl[i]) && ((zin1 + 10.8) / 2. < zobl[i + 1])) {
								if ((ztr1 > zobl[i]) && (ztr1 < zobl[i + 1])) {
									h[i][j]->Fill(ach[scount1]);
								}
							}
						}
						else if ((10.8 - fabs(zin1)) * fabs(tan(theta[1])) >= 3.0) {
							for (int i = 0; i < 14; i++) {
								if ((ztr1 > zobl[i]) && (ztr1 < zobl[i + 1])) {
									h[i][j]->Fill(ach[scount1]);
								}
							}
						}
					}
				}*/
				/*if ((phi[0] >= (float)(j - 1) * (2. * PI) / 9. + 0.3) && (phi[0] < (float)(j + 2) * (2. * PI) / 9. - 0.3)) {
					for (int i = 0; i < 14; i++) {
						for (int k = 0; k < 4; k++) {
							if ((ztr > zobl[i]) && (ztr < zobl[i + 1]) && (beam > egran[k]) && (beam < egran[k + 1])) {
								h[i][j][k]->Fill(ach[scount1]);
							}
						}
					}
				}*/

				//pedest 9X14X9 hist
				/*for (int k = 0; k < 14; k++) {
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
						if (it < 9)
							h0[it][k][j]->Fill(ach[scount1]);
					}
					if ((ztr1 > zobl[k]) && (ztr1 < zobl[k + 1]) && (region[1] == 1)) {
						int it = 9;
						for (int i = 0; i < 3; i++) {
							if ((phi[0] >= (lgran[0] + i * (lgran[1] - lgran[0]) / 3.)) && (phi[0] < (lgran[0] + (i + 1) * (lgran[1] - lgran[0]) / 3.)))
								it = i;
							if ((phi[0] >= (lgran[2] + i * (lgran[3] - lgran[2]) / 6.)) && (phi[0] < (lgran[2] + (i + 1) * (lgran[3] - lgran[2]) / 6.)))
								it = i + 3;
							if ((phi[0] >= (lgran[2] + (i + 3) * (lgran[3] - lgran[2]) / 6.)) && (phi[0] < (lgran[2] + (i + 4) * (lgran[3] - lgran[2]) / 6.)))
								it = i + 6;
						}
						if (it < 9)
							h0[it][k][j]->Fill(ach[scount1]);
					}
				}*/


				/*if ((phi[0] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[0] < (float)(j + 2) * (2. * PI) / 9.)) {
					h[j]->Fill(ach[scount1]);
				}*/

				/*if ((phi[0] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[0] < (float)(j + 2) * (2. * PI) / 9.)) {
					if ((schr[scount] < 0.5) || ((schr[scount] > 0.5) && (tchr[scount] < 110))) {
						//if (((phi[0] >= Goodphi[4 * j]) && phi[0] < Goodphi[4 * j + 1]) || ((phi[0] >= Goodphi[4 * j + 2]) && phi[0] < Goodphi[4 * j + 3])) {
						ach1[scount] = ach[scount] - ped2[j];
						hprof[j]->Fill(phi[0], ach1[scount] * sin(theta[0]));
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



	
	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles.root", "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "zobl%d,sensor%d", m + 1, i + 1);
			hprof[m][i]->Write(title);
		}
	}
	MyFile->Close();*/
	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profilesdate.root", "RECREATE");
	for (int m = 0; m < 5; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "date%d,sensor%d", m + 1, i + 1);
			hprofd[i][m]->Write(title);
		}
	}
	MyFile->Close();*/
	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles.root", "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "zobl%d,sensor%d", m + 1, i + 1);
			hprof[m][i]->Write(title);
		}
	}
	MyFile->Close();*/
	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profilesmean.root", "RECREATE");
	for (Int_t i = 0; i < 9; i++) {
		sprintf(title, "sensor%d", i + 1);
		hprof[i]->Write(title);
	}
	MyFile->Close();*/

	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles1.root", "RECREATE");
	for (Int_t i = 0; i < 9; i++) {
		sprintf(title, "sensor%d",i + 1);
		hprof[i]->Write(title);
	}
	MyFile->Close();*/
	
	/*TFile* MyFile1 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedest.root", "RECREATE");
	for (int j = 0; j < 9; j++) {
		for (int k = 0; k < 14; k++) {
			for (int i = 0; i < 9; i++) {
				sprintf(title, "ped_sens%d,z%d,phi%d", i + 1, k + 1, j + 1);
				//cout << h[j][k][i]->GetMean(1) << endl;
				h0[j][k][i]->Write(title);
			}
		}
	}
	MyFile1->Close();*/

	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/pedest1.root", "RECREATE");
	for (int j = 0; j < 9; j++) {
		sprintf(title, "ped_sens%d", j + 1);
		//cout << h[j][k][i]->GetMean(1) << endl;
		h[j]->Write(title);
	}
	MyFile->Close();*/

	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/phintegr_pedest.root", "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "ped_phintegr_zobl%d,sensor%d", m + 1, i + 1);
			//cout << h[m][i]->GetMean(1) << endl;
			h[m][i]->Write(title);
		}
	}
	MyFile->Close();*/
	/*TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/phintegr_pedest_en.root", "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			for (int k = 0; k < 4; k++) {
				sprintf(title, "ped_phintegr_zobl%d,sensor%d,en%d", m + 1, i + 1,k+1);
				//cout << h[m][i]->GetMean(1) << endl;
				h[m][i][k]->Write(title);
			}
		}
	}
	MyFile->Close();*/

	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/nvsphi.root", "RECREATE");
	for (int a = 0; a < 9; a++) {
		for (int i = 0; i < 14; i++) {
			sprintf(title, "zobl%d,sensor%d", i + 1, a + 1);
			hii[a][i]->Write(title);
		}
	}
	MyFile->Close();
	cout << "done!" << endl;

}

void nvsphiorig() {
	TChain chain("t1");
	//chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected2019/true1__**exp00**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("run", &run);
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	char name[20];
	char title[100];
	TH1* hii[9][1];

	for (int l = 0; l < 1; l++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "histog%d",9*l + i + 1);
			sprintf(title, "#phi,sensor%d;#phi", i + 1);
			hii[i][l] = new TH1F(name, title, 400, truegran1[i][1] - 2., truegran1[i][1] + 2.);
		}
	}
	int tl = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		for(int i = 0; i < 9; i++){
			hii[i][0]->Fill(phi[0]);
			hii[i][0]->Fill(phi[1]);
		}
	}
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/nvsphio19.root", "RECREATE");
	for (int l = 0; l < 1; l++) {
		for (int a = 0; a < 9; a++) {
			sprintf(title, "sensor%d", a + 1);
			hii[a][l]->Write(title);
		}
	}
	MyFile->Close();
	cout << "done!" << endl;
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

Double_t gg(Double_t* x, Double_t* par)
{
	double a1 = par[0];
	double m1 = par[1];
	double s1 = par[2];
	double dx1 = (x[0] - m1) / s1;
	double f1 = a1 * exp(-dx1 * dx1 / 2);///std::sqrt(2*TMath::Pi())/s1;
	double a2 = par[3];
	double m2 = par[4];
	double s2 = par[5];
	double dx2 = (x[0] - m2) / s2;
	double f2 = a2 * exp(-dx2 * dx2 / 2);///std::sqrt(2*TMath::Pi())/s2;
	return f1 + f2;
}

Double_t gp2x(Double_t* x, Double_t* par)
{
	double xr = par[0];
	double sr = par[1];
	double a = par[4] / 2;
	double b = par[3] - 2 * a * xr;
	double c = par[2] - a * xr * xr - b * xr;
	double dx = (x[0] - xr) / sr;
	double f = (c + b * x[0] + a * (x[0] * x[0] + sr * sr)) * (1 + erf(dx / sqrt(2.0))) / 2 +
		(b + a * (x[0] + xr)) * sr / sqrt(2 * 3.1415927) * exp(-dx * dx / 2);
	return f;
}

double scphi(double* x, double* par)
{
	double f = 0;
	for (size_t i = 0; i < 3; i++)
		f += aerogel(i, x[0], par);
	f += shifter(x[0], par);
	return f;
}

void linesa() {
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles1.root");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0;
	while ((key = (TKey*)next())) {
		TString name(key->GetName());
		TString type = "sensor";
		if (name.BeginsWith(type) && name.Contains("sensor")) {
			//cin.get();
			int counter;
			sscanf(name.Data(), "sensor%d", &counter);

			TProfile* h_amp_phi = (TProfile*)f->Get(Form("sensor%d", counter));

			h_amp_phi->Draw();

			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			int j = counter - 1;
			TLine* line11 = new TLine(Goodphi[4*j], 0, Goodphi[4*j], gPad->GetUymax());
			TLine* line12 = new TLine(Goodphi[4*j+1], 0, Goodphi[4*j+1], gPad->GetUymax());
			TLine* line13 = new TLine(Goodphi[4*j+2], 0, Goodphi[4*j+2], gPad->GetUymax());
			TLine* line14 = new TLine(Goodphi[4*j+3], 0, Goodphi[4*j+3], gPad->GetUymax());
			line11->SetLineColor(kRed);
			line12->SetLineColor(kRed);
			line13->SetLineColor(kRed);
			line14->SetLineColor(kRed);
			line11->Draw();
			line12->Draw();
			line13->Draw();
			line14->Draw();
			c->Update();
			cin.get();
		}
	}
}


void drawsomefit() {
	TCanvas* ci = new TCanvas("c1","c1",1500,1000);
	ci->SetFillColor(0);
	ci->SetFrameFillColor(0);
	//ci->Divide(3,3);
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	double scal = PI / 180;
	int interobl[9] = {14,1,1,2,2,4,4,5,14};
	int intercount[9] = {8,3,6,1,9,3,4,7,1};
	double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles.root");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0;
	while ((key = (TKey*)next()) && counterr < 1) {
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
			h_amp_phi->SetTitle(Form("counter%d, z #in (%.1f, %.1f);#phi,rad;amplitude,p.e.", counter, zobl[obl-1], zobl[obl]));
			gStyle->SetOptStat(11);
			//ci->cd(counterr);
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
				par[12] = par[12] - 0.003;
			if (counter == 7 && obl == 4)
				par[15] = par[15] + 0.005;
			if (counter == 7 && obl == 6)
				par[15] = par[15] + 0.003;
			if (counter == 6 && (obl == 6 || obl == 7 || obl == 8 || obl == 9))
				par[13] = par[13] + 0.01;
			if (counter == 8 && obl == 1)
				par[12] = par[12] - 0.003;
			if (counter == 8 && obl == 12)
				par[12] = par[12] + 0.003;


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
					par[15] = par[15] - 0.005 + (obl - 11.) / 100.;
					par[19] = (0.3) * scal;
				}
				if (obl == 4)
					par[15] = par[15] - 0.003;
			}
			if (counter == 1 && obl == 12)
				par[12] = par[12] + 0.005;
			if (counter == 9) {
				par[12] = par[12] - 0.001;
				if (obl == 14) {
					par[12] = par[12] + 0.002;
					par[15] = par[15] - 0.006;
				}
				if (obl == 10) {
					par[12] = par[12] - 0.002;
					par[15] = par[15] - 0.002;
				}
				par[17] = (0.4) * scal;
				par[19] = (0.4) * scal;
			}
			if (counter == 7 && obl == 1) {
				par[14] = par[14] - 0.015;
				par[12] = par[12] - 0.002;
			}
			if (counter == 7 && obl == 2) {
				par[15] = par[15] + 0.005;
				par[12] = par[12] + 0.003;
			}
			if (counter == 7 && obl == 11) {
				par[15] = par[15] - 0.005;
				par[12] = par[12] - 0.003;
			}
			if (counter == 7 && obl == 4)
				par[19] = (0.4) * scal;
			if (counter == 7 && obl == 5)
				par[15] = par[15] + 0.005;
			if (counter == 7 && obl == 7)
				par[15] = par[15] + 0.005;
			if (counter == 7 && obl == 10)
				par[12] = par[12] - 0.005;
			if (counter == 8 && obl == 3)
				par[15] = par[15] - 0.005;
			if (counter == 8 && obl == 5)
				par[12] = par[12] + 0.005;
			if (counter == 8 && obl == 9)
				par[12] = par[12] - 0.003;
			if (counter == 9 && obl == 13)
				par[15] = par[15] + 0.005;
			if (counter == 9 && obl == 7)
				par[15] = par[15] + 0.003;
			if (counter == 9 && obl == 1)
				par[15] = par[15] + 0.003;
			if (counter == 6 && (obl == 6 || obl == 7))
				par[13] = par[13] + 0.01;
			if (counter == 8 && (obl != 12 && obl != 1 && obl != 2 && obl != 3 && obl != 4 && obl != 14 && obl != 5 && obl != 6)) {
				par[12] = par[12] + 0.001;
				if (obl == 13)
					par[12] = par[12] + 0.006;
				par[17] = (0.6) * scal;
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
			if (counter == 8)
				fn->SetParLimits(16, a - 2., a + 15. + 1. * obl);
			if (counter == 6)
				fn->SetParLimits(16, a - 2., a + 25. + 1. * obl);
			if (counter == 7) {
				fn->SetParLimits(16, a, a + 25. + 1. * obl);
				//if (obl == 1)
				//	fn->SetParLimits(16, a-1., a + 25. + 1. * obl);
			}
			else
				fn->SetParLimits(16, a - 2., a + 15. + 1. * obl);


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
				if (o == 0) {
					fn->SetParLimits(o, yyy - 0.5, yyy1 + 0.1);
					if (counter == 8)
						fn->SetParLimits(o, yyy + 0.5, yyy1 + 0.5);
				}
				else
					fn->SetParLimits(o, yyy - 0.2, yyy + 0.3);
			}
			fn->SetParLimits(3, yyy + 0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			if (obl == 14) {
				if (counter < 5)
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
					if (counter == 8 && obl > 4)
						fn->SetParLimits(o, yyy + 2., yyy1 + 0.3);
					if (counter == 8 && obl == 1)
						fn->SetParLimits(o, yyy + 1., yyy1 + 2.3);
					else
						fn->SetParLimits(o, yyy - 0.5, yyy1 + 0.3);
				}
				else {
					if (counter == 7 && obl < 3)
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.0);
					if (counter == 8)
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
					if (counter == 3 && obl == 8)
						fn->SetParLimits(o, yyy - 1.2, yyy + 0.2);
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
			if (counter == 7 && obl == 1)
				fn->SetParLimits(3, yyy, yyy + 1.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(3, yyy + 2., yyy + 8.);
				else if (counter == 9)
					fn->SetParLimits(3, yyy + 1., yyy + 20.);
				else if (counter == 8)
					fn->SetParLimits(3, yyy + 2., yyy + 20.);
				else if (counter == 7)
					fn->SetParLimits(3, yyy + 1., yyy + 20.);
				else if (counter == 6)
					fn->SetParLimits(3, yyy + 5., yyy + 20.);
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
			fn->SetParLimits(4, yyy4 + 0.5, yyy4 + 8.);
			if (counter == 9)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);
			if (counter == 8)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 8.);
				else if (counter == 9)
					fn->SetParLimits(4, yyy4 + 1., yyy4 + 20.);
				else if (counter == 8)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 20.);
				else if (counter == 7)
					fn->SetParLimits(4, yyy4 + 1., yyy4 + 20.);
				else
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 20.);
			}

			fn->ReleaseParameter(17);
			fn->ReleaseParameter(18);
			fn->ReleaseParameter(19);
			fn->SetParLimits(17, 6.0e-3, 1.5e-2);
			fn->SetParLimits(18, 3.0e-3, 1.2e-2);
			fn->SetParLimits(19, 6.0e-3, 1.4e-2);
			if (counter == 6)
				fn->SetParLimits(18, 3.0e-3, 10.0e-3);
			if (counter == 5)
				fn->SetParLimits(19, 1.0e-2, 2.0e-2);
			if (counter == 8)
				fn->SetParLimits(18, 4.3e-3, 1.2e-2);
			if (counter == 8 && obl == 1)
				fn->SetParLimits(17, 4.3e-3, 1.1e-2);
			if (counter == 7)
				fn->SetParLimits(18, 4.3e-3, 1.2e-2);



			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();


			std::vector<double> xv2(4);
			std::vector<double> yv2(4);
			fn->GetParameters(&par[0]);
			

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
				//size_t xxx;
				//cin>>xxx;
			}
			//cin.get();
		}
	}
}

void testfit() {
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	std::vector<double> xv1(12);
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles.root");
	TFile* f1 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/fittedprofiles.root", "RECREATE");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0, ennn = 0;
	while ((key = (TKey*)next())) {
		counterr += 1; ennn += 1;
		if (ennn == 5)
			ennn = 1;
		TString name(key->GetName());
		TString type = "zobl";
		if (/*name.BeginsWith(type) && */name.Contains("sensor")) {
			//cin.get();
			int obl, counter, enbeam, date;
			//sscanf(name.Data(), "zobl%d,sensor%d,en%d", &obl, &counter, &enbeam);
			sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
			//sscanf(name.Data(), "date%d,sensor%d", &date, &counter);
			TDirectory* dir = f1->mkdir(Form("zobl%d,sensor%d", obl, counter));
			//counter = 7;
			//obl = counterr;
			//date = counterr;
			//TDirectory* dir = f1->mkdir(Form("date%d,sensor%d", date, counter));
			TProfile* h_amp_phi = (TProfile*)f->Get(Form("zobl%d,sensor%d", obl, counter));
			//TProfile* h_amp_phi = (TProfile*)f->Get(Form("date%d,sensor%d", date, counter));

			h_amp_phi->Draw();
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			std::vector<double> par(20);
			std::vector<double> parc(20);


			double xmin = (counter * 40. - 35.) * scal + dphi;
			double xmax = (counter * 40. - 15.) * scal + dphi;
			TF1* fxc = new TF1("fxc", gg, xmin, xmax, 6);
			double xm = h_amp_phi->GetBinCenter(h_amp_phi->GetMaximumBin());
			fxc->SetParameter(0, h_amp_phi->GetMaximum());
			fxc->SetParameter(1, xm);
			fxc->SetParameter(2, 0.5 * scal);
			fxc->SetParameter(3, 5);
			fxc->SetParameter(4, xm);
			fxc->SetParameter(5, 10 * scal);
			fxc->FixParameter(0, h_amp_phi->GetMaximum());
			fxc->FixParameter(1, xm);
			fxc->FixParameter(2, 0.5 * scal);
			fxc->SetParLimits(4, xm - dphi, xm + dphi);
			fxc->SetParLimits(5, 5 * scal, 20 * scal);
			h_amp_phi->Fit(fxc, "", "", xmin, xmax);
			c->Update();
			fxc->ReleaseParameter(0);
			fxc->ReleaseParameter(1);
			fxc->ReleaseParameter(2);
			fxc->SetParLimits(0, 0, 100);
			fxc->SetParLimits(1, (counter * 40. - 30.) * scal + dphi, (counter * 40. - 20.) * scal + dphi);
			fxc->SetParLimits(2, 0.1 * scal, 0.6 * scal);
			h_amp_phi->Fit(fxc, "", "", xmin, xmax);
			c->Update();
			//
			xmin = ((counter - 1) * 40. - 10.) * scal + dphi;
			xmax = fxc->GetParameter(1) - 3 * fxc->GetParameter(2);//((counter-1)*40.+10.)*scal+dphi;
			TF1* fxl = new TF1("fxl", gp2x, xmin, xmax, 5);
			fxl->SetParameter(0, ((counter - 1) * 40 + 2.5) * scal);
			fxl->SetParameter(1, 0.5 * scal);
			fxl->SetParameter(2, 5);
			fxl->SetParameter(3, 0);
			fxl->SetParameter(4, 0);
			h_amp_phi->Fit(fxl, "", "", xmin, xmax);
			c->Update();
			//
			xmin = fxc->GetParameter(1) + 3 * fxc->GetParameter(2);//(counter*40.-10.)*scal+dphi;
			xmax = (counter * 40. + 10.) * scal + dphi;
			TF1* fxr = new TF1("fxr", gp2x, xmin, xmax, 5);
			fxr->SetParameter(0, (counter * 40 + 2.5) * scal);
			fxr->SetParameter(1, -0.5 * scal);
			fxr->SetParameter(2, 5);
			fxr->SetParameter(3, 0);
			fxr->SetParameter(4, 0);
			h_amp_phi->Fit(fxr, "", "", xmin, xmax);
			c->Update();
			//
			par[12] = fxl->GetParameter(0);
			par[13] = fxc->GetParameter(1);
			par[14] = par[13] + 12.5 * scal;
			par[15] = fxr->GetParameter(0);
			//
			xmin = par[12];
			xmax = par[13];
			par[0] = 2 * fxl->Eval(xmin);
			par[1] = fxl->Eval(xmin + 0.33 * (xmax - xmin));
			par[2] = fxl->Eval(xmax - 0.33 * (xmax - xmin));
			par[3] = fxl->Eval(xmax);
			//
			xmin = par[13];
			xmax = par[14];
			par[4] = fxr->Eval(xmin);
			par[5] = fxr->Eval(xmin + 0.33 * (xmax - xmin));
			par[6] = fxr->Eval(xmax - 0.33 * (xmax - xmin));
			par[7] = fxr->Eval(xmax);
			//
			xmin = par[14];
			xmax = par[15];
			par[8] = fxr->Eval(xmin);
			par[9] = fxr->Eval(xmin + 0.33 * (xmax - xmin));
			par[10] = fxr->Eval(xmax - 0.33 * (xmax - xmin));
			par[11] = 2 * fxr->Eval(xmax);
			//
			par[16] = fxc->GetParameter(0);
			par[17] = fxl->GetParameter(1);
			par[18] = fxc->GetParameter(2);
			par[19] = fabs(fxr->GetParameter(1));


			parc[12] = par[12];
			parc[13] = par[13];
			parc[14] = par[14];
			parc[15] = par[15];

			
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
			std::vector<double> xv(12);

			//______fit___
			for (size_t i = 0; i < 12; i++)
				fn->FixParameter(i, fn->GetParameter(i));
			fn->FixParameter(12, fn->GetParameter(12));
			fn->FixParameter(13, fn->GetParameter(13));
			fn->FixParameter(14, fn->GetParameter(14));
			fn->FixParameter(15, fn->GetParameter(15));
			fn->SetParLimits(16, fxc->GetParameter(0) - 0.3, fxc->GetParameter(0) + 9.);
			fn->FixParameter(17, fn->GetParameter(17));
			fn->FixParameter(18, fn->GetParameter(18));
			fn->FixParameter(19, fn->GetParameter(19));

			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();

			for (size_t i = 0; i < 3; i++)
				fn->ReleaseParameter(i);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();
			for (size_t i = 3; i < 8; i++)
				fn->ReleaseParameter(i);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();
			for (size_t i = 8; i < 12; i++)
				fn->ReleaseParameter(i);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();

			fn->ReleaseParameter(12);
			fn->ReleaseParameter(13);
			fn->ReleaseParameter(15);
			fn->SetParLimits(12, fn->GetParameter(12) - dphi, fn->GetParameter(12) + dphi);
			fn->SetParLimits(13, fn->GetParameter(13) - dphi, fn->GetParameter(13) + dphi);
			fn->SetParLimits(15, fn->GetParameter(15) - dphi, fn->GetParameter(15) + dphi);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();

			fn->ReleaseParameter(17);
			fn->ReleaseParameter(18);
			fn->ReleaseParameter(19);
			fn->SetParLimits(17, 2.0e-3, 1.4e-2);
			fn->SetParLimits(18, 1.0e-3, 1.2e-2);
			fn->SetParLimits(19, 2.0e-3, 2.0e-2);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();

			for (size_t i = 0; i < 2; i++) {
				xmin = fn->GetParameter(12) - 3 * fn->GetParameter(17);
				xmax = fn->GetParameter(15) + 3 * fn->GetParameter(19);
				h_amp_phi->Fit("scphi", "", "", xmin, xmax);
			}


			//___writing__

			std::vector<double> xv2(4);
			std::vector<double> yv2(4);
			fn->GetParameters(&par[0]);
			xv1[0] = truegran1[counter - 1][0];
			xv1[3] = truegran1[counter - 1][1];
			xv1[1] = xv1[0] + (xv1[3] - xv1[0]) / 3;
			xv1[2] = xv1[3] - (xv1[3] - xv1[0]) / 3;
			xv1[4] = truegran1[counter - 1][1];
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
					xv2[i] = xv[4 * n + i];
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



			for (int ii = 0; ii < 4; ii++)
				gran1[ii][counter - 1][obl - 1] = fn->GetParameter(12 + ii);
			pik[counter - 1][obl - 1] = fn->GetParameter(16);

			a1[counter - 1][obl - 1] = fn->GetParameter(12) + 0.02;
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
			}

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
					//sprintf(titlea, "date%d,sensor%d,aerogel%d", date, counter, c + 1);
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
				/*sprintf(titles, "date%d,sensor%d,shifter", date, counter);
				sprintf(titlef, "date%d,sensor%d,full", date, counter);
				sprintf(titlep, "date%d,sensor%d,profile", date, counter);
				sprintf(title, "date%d,sensor%d,canvas", date, counter);*/
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
	}
	ofstream fout1;
	fout1.open("/work/users/kladov/snd2k/R006-003/maindir/gran1.cpp");
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

	ofstream fout2;
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

void testforfit() {
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	std::vector<double> xv1(12);
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profilesdate.root");
	TFile* f1 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/fittedprofilesdate.root", "RECREATE");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0, ennn = 0;
	double grandata[4][9][5];
	while ((key = (TKey*)next()) && counterr < 45) {
		counterr += 1; ennn += 1;
		if (ennn == 5)
			ennn = 1;
		TString name(key->GetName());
		TString type = "zobl";
		if (/*name.BeginsWith(type) && */name.Contains("sensor")) {
			//cin.get();
			int obl, counter, enbeam, date;
			//sscanf(name.Data(), "zobl%d,sensor%d,en%d", &obl, &counter, &enbeam);
			//sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
			sscanf(name.Data(), "date%d,sensor%d", &date, &counter);
			//TDirectory* dir = f1->mkdir(Form("zobl%d,sensor%d", obl, counter));
			//counter = 6;
			obl = 5;
			//date = counterr;
			TDirectory* dir = f1->mkdir(Form("date%d,sensor%d", date, counter));
			//TProfile* h_amp_phi = (TProfile*)f->Get(Form("zobl%d,sensor%d", obl, counter));
			TProfile* h_amp_phi = (TProfile*)f->Get(Form("date%d,sensor%d", date, counter));

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
			par[13] = maxphi;
			if (counter == 8 && obl == 3)
				par[15] = par[15] - 0.01;
			if (counter == 8 && obl == 2)
				par[12] = par[12] - 0.003;
			if (counter == 7 && obl == 3)
				par[12] = par[12] - 0.003;
			if (counter == 7 && obl == 4)
				par[15] = par[15] + 0.005;
			if (counter == 7 && obl == 6)
				par[15] = par[15] + 0.003;
			if (counter == 6 && (obl == 6 || obl == 7 || obl == 8 || obl == 9))
				par[13] = par[13] + 0.01;
			if (counter == 8 && obl == 1)
				par[12] = par[12] - 0.003;
			if (counter == 8 && obl == 12)
				par[12] = par[12] + 0.003;
			
			
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
				if (obl == 4)
					par[15] = par[15] - 0.003;
			}
			if (counter == 1 && obl == 12)
				par[12] = par[12] + 0.005;
			if (counter == 8 && obl == 1)
				par[15] = par[15] - 0.002;
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
				par[14] = par[14] - 0.015;
				par[12] = par[12] - 0.002;
				par[15] = par[15] + 0.01;
				par[13] = par[13] + 0.002;
			}
			if (counter == 7 && obl == 2) {
				par[15] = par[15] + 0.003;
				par[12] = par[12] + 0.005;
				par[14] = par[14] - 0.002;
			}
			if (counter == 7 && obl == 11) {
				par[15] = par[15] - 0.005;
				par[12] = par[12] - 0.003;
			}
			if (counter == 7 && obl == 4) 
				par[19] = (0.4) * scal;
			if (counter == 7 && obl == 5)
				par[15] = par[15] + 0.005;
			if (counter == 7 && obl == 7)
				par[15] = par[15] + 0.005;
			if (counter == 7 && obl == 10)
				par[12] = par[12] - 0.005;
			if (counter == 7 && obl == 14)
				par[15] = par[15] - 0.01;
			if (counter == 8 && obl == 3)
				par[15] = par[15] - 0.005;
			if (counter == 8 && obl == 5)
				par[12] = par[12] + 0.005;
			if (counter == 8 && obl == 9)
				par[12] = par[12] - 0.003;
			if (counter == 9 && obl == 13)
				par[15] = par[15] + 0.005;
			if (counter == 9 && obl == 7)
				par[15] = par[15] + 0.003;
			if (counter == 9 && obl == 14)
				par[13] = par[13] - 0.003;
			if (counter == 9 && obl == 1)
				par[15] = par[15] - 0.006;
			if (counter == 9 && obl == 4)
				par[13] = par[13] - 0.002;
			if (counter == 9 && obl == 9)
				par[13] = par[13] - 0.002;
			if (counter == 9 && obl == 10)
				par[13] = par[13] - 0.002;
			if (counter == 9 && obl == 11)
				par[13] = par[13] - 0.002;
			if (counter == 6 && (obl == 6 || obl == 7))
				par[13] = par[13] + 0.01;
			if (counter == 8 && (obl != 12 && obl != 1 && obl != 2 && obl != 3 && obl != 4 && obl != 14 && obl != 5 && obl != 6)) {
				par[12] = par[12] + 0.001;
				if (obl == 13)
					par[12] = par[12] + 0.006;
				par[17] = (0.6) * scal;
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
				fn->SetParLimits(16, a-2., a + 15.+1.*obl);
			if (counter == 6)
				fn->SetParLimits(16, a-2., a + 25. + 1. * obl);
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
				if (o == 0) {
					fn->SetParLimits(o, yyy - 0.5, yyy1 + 0.1);
					if (counter == 8)
						fn->SetParLimits(o, yyy + 0.5, yyy1 + 0.5);
				}
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
					if (counter == 8 && obl == 1)
						fn->SetParLimits(o, yyy + 1., yyy1 + 2.3);
					else
						fn->SetParLimits(o, yyy - 0.5, yyy1 + 0.3);
				}
				else {
					if(counter == 7 && obl <3)
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.0);
					if (counter == 8)
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
					if (counter == 3 && obl == 8)
						fn->SetParLimits(o, yyy - 1.2, yyy + 0.2);
					else
						fn->SetParLimits(o, yyy - 0.2, yyy + 0.2);
				}
			}
			cout << yyy << endl;
			fn->SetParLimits(3, yyy+0.5, yyy + 8.);
			if (counter == 9)
				fn->SetParLimits(3, yyy + 0.5, yyy + 10.);
			if (counter == 8)
				fn->SetParLimits(3, yyy + 1.5, yyy + 10.);
			if (counter == 7 && obl == 1)
				fn->SetParLimits(3, yyy, yyy + 1.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(3, yyy + 2., yyy + 8.);
				else if (counter == 9)
					fn->SetParLimits(3, yyy + 0.5, yyy + 6.);
				else if (counter == 8)
					fn->SetParLimits(3, yyy + 2., yyy + 20.);
				else if (counter == 7)
					fn->SetParLimits(3, yyy + 1., yyy + 20.);
				else if (counter == 6)
					fn->SetParLimits(3, yyy + 5., yyy + 20.);
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
				fn->SetParLimits(4, yyy4 + 0.5, yyy4 + 10.);
			if (counter == 8)
				fn->SetParLimits(4, yyy4 + 1.5, yyy4 + 10.);
			if (obl == 14) {
				if (counter < 5)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 8.);
				else if (counter == 9)
					fn->SetParLimits(4, yyy4 + 0.5, yyy4 + 20.);
				else if (counter == 8)
					fn->SetParLimits(4, yyy4 + 2., yyy4 + 20.);
				else if (counter == 7)
					fn->SetParLimits(4, yyy4 + 1., yyy4 + 20.);
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
			if (counter == 8 && obl == 1)
				fn->SetParLimits(17, 4.3e-3, 1.1e-2);
			if (counter == 7)
				fn->SetParLimits(18, 4.3e-3, 1.2e-2);

			

			h_amp_phi->Fit("scphi", "","", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();
			

			std::vector<double> xv2(4);
			std::vector<double> yv2(4);
			fn->GetParameters(&par[0]);
			xv1[0] = truegran1[counter - 1][0];
			xv1[3] = truegran1[counter - 1][1];
			xv1[1] = xv1[0] + (xv1[3] - xv1[0]) / 3;
			xv1[2] = xv1[3] - (xv1[3] - xv1[0]) / 3;
			xv1[4] = truegran1[counter - 1][1];
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
			
			grandata[0][counter - 1][date - 1] = par[12];
			grandata[1][counter - 1][date - 1] = par[13];
			grandata[2][counter - 1][date - 1] = par[14];
			grandata[3][counter - 1][date - 1] = par[15];

			for (int ii = 0; ii < 4; ii++)
				gran1[ii][counter - 1][obl - 1] = fn->GetParameter(12 + ii);
			pik[counter-1][obl-1] = fn->GetParameter(16);

			a1[counter - 1][obl - 1] = fn->GetParameter(12) + 0.02;
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
			}

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
					//sprintf(titlea, "zobl%d,sensor%d,aerogel%d", obl, counter, c + 1);
					sprintf(titlea, "date%d,sensor%d,aerogel%d", date, counter, c + 1);
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
				/*sprintf(titles, "zobl%d,sensor%d,shifter", obl, counter);
				sprintf(titlef, "zobl%d,sensor%d,full", obl, counter);
				sprintf(titlep, "zobl%d,sensor%d,profile", obl, counter);
				sprintf(title, "zobl%d,sensor%d,canvas", obl, counter);*/
				sprintf(titles, "date%d,sensor%d,shifter", date, counter);
				sprintf(titlef, "date%d,sensor%d,full", date, counter);
				sprintf(titlep, "date%d,sensor%d,profile", date, counter);
				sprintf(title, "date%d,sensor%d,canvas", date, counter);
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
	}
	ofstream fout1;
	fout1.open("/work/users/kladov/snd2k/R006-003/maindir/gran1.cpp");
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
	
	ofstream fout2;
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
	
	double data[5];
	for (int i = 0; i < 5; i++)
		data[i] = (double)(i+1);
	for (int i = 0; i < 9; i++) {
		TGraph* graph = new TGraph(5, data, grandata[0][i]);
		graph->SetMarkerSize(1.5);
		graph->Draw("APL");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
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
	/*ofstream fout1;
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
	}*/


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
	double scal1 = 180. / PI;
	ofstream fout1;
	std::vector<double> xv1(12);
	int n = 14; 
	double koef[14] = { 1.0707799, 1.07214, 1.08057, 1.08138, 1.07637, 1.07359, 1.07303, 1.08661, 1.0774, 1.07921, 1.08933, 1.08917, 1.07887, 1.07091 };
	fout1.open("/work/users/kladov/snd2k/R006-003/maindir/map.txt");
	for (int j = 0; j < 9; j++) {
		xv1[0] = truegran1[j][0];
		xv1[3] = truegran1[j][1];
		xv1[1] = xv1[0] + (xv1[3] - xv1[0]) / 3;
		xv1[2] = xv1[3] - (xv1[3] - xv1[0]) / 3;
		xv1[4] = truegran1[j][1];
		xv1[7] = truegran1[j][2];
		xv1[5] = xv1[4] + (xv1[7] - xv1[4]) / 3;
		xv1[6] = xv1[7] - (xv1[7] - xv1[4]) / 3;
		xv1[8] = truegran1[j][2];
		xv1[11] = truegran1[j][3];
		xv1[9] = xv1[8] + (xv1[11] - xv1[8]) / 3;
		xv1[10] = xv1[11] - (xv1[11] - xv1[8]) / 3;
		fout1 << Form("Counter %d", j+1) << endl << Form("Shifter %f %f %f",shift[j][0]/koef[0],shift[j][1]/koef[0],shift[j][2]) << endl << "AmplitudeCorrection 1" << endl << "Aerogel" << endl;
		for (int k = 0; k < 3; k++) {
			fout1 << n << " " << xv1[0+4*k]*scal1 << " " << xv1[1+4*k]*scal1 << " " << xv1[2+4*k]*scal1 << " " << xv1[3+4*k]*scal1 << endl;
			for (int i = 0; i < 14; i++) {
				fout1 << zobl[i]-0.0 << " " << yfm[j][i][0+4*k]/koef[i] << " " << yfm[j][i][1+4*k]/koef[i] << " " << yfm[j][i][2+4*k]/koef[i] << " " << yfm[j][i][3+4*k]/koef[i] << endl;
			}
		}
	}

	/*double oblast[14];
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
		//cin.get();
		par[j][0] = fn->GetParameter(0);
		par[j][1] = fn->GetParameter(1);
		par[j][2] = fn->GetParameter(2);
	}
	ofstream fout;
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
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles.root");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0, ennn = 0;
	while ((key = (TKey*)next()) && counterr < 14) {
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
			counter = 7;
			obl = counterr;
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
			cin.get();
		}
	}
}

void checkbord() {
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	std::vector<double> xv1(12);
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profilesmean.root");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0;
	while ((key = (TKey*)next())) {
		counterr += 1;
		TString name(key->GetName());
		TString type = "sensor";
		if (name.BeginsWith(type) && name.Contains("sensor")) {
			//cin.get();
			int obl, counter;
			//sscanf(name.Data(), "zobl%d,sensor%d,en%d", &obl, &counter, &enbeam);
			sscanf(name.Data(), "sensor%d",&counter);
			//counter = 7;
			//obl = counterr;
			TProfile* h_amp_phi = (TProfile*)f->Get(Form("sensor%d",counter));

			h_amp_phi->Draw();
			TLine* line1 = new TLine(truegran1m[counter - 1][0], 0, truegran1m[counter - 1][0], gPad->GetUymax());
			TLine* line2 = new TLine(truegran1m[counter - 1][1], 0, truegran1m[counter - 1][1], gPad->GetUymax());
			TLine* line3 = new TLine(truegran1m[counter - 1][2], 0, truegran1m[counter - 1][2], gPad->GetUymax());
			TLine* line4 = new TLine(truegran1m[counter - 1][3], 0, truegran1m[counter - 1][3], gPad->GetUymax());
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
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			cin.get();
		}
	}
}

void go(){ 

//	copy1();
//	copy2();
//	findruns();
//	raspr();
//	zraspr(); 
//	linesforz();
//	raspramp();
//	timespectra();
//	linesa();
//	testforfit();
//	testfit();
//	bordersplot();
//	modelfile();
//	checkmod();
//	checkbord();
//	testt();	
//	calorimetr();
//	draweffspectre();
//	drawcalorimetr();
//	drawsomefit();
//	ampltime();
//	drawamptime();
//	dopeff();
//	nvsphiorig();
	drawnvsphi();

	//pscp D:\programs\select_ach1.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/select_ach1.cpp
	//pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/select_ach1.cpp D:\programs\select_ach1.cpp
}

