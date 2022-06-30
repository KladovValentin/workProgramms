#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TH1.h"
#include "TMath.h"
#include <cmath>
#include "TLine.h"
#include "TChain.h"
#include "TSystem.h"
#include "TString.h"
#include "TF1.h"
#include "TSystemDirectory.h"
#include "TCollection.h"
#include "TList.h"
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector> 
using namespace std;

void m4(double x[], double y[]) {
	cout << (pow(x[0],2))*x[1] << "    " << x[1] << endl;
	cout << TMath::Erf(y[0]) << "    " << y[1] << endl;
}


void raspr() {

	TCanvas* c1 = new TCanvas("c1", "A", 200, 10, 500, 300);

	TProfile* hprof[18];
	char name[20];
	char title[100];
	for (Int_t i = 0; i < 18; i++) {
		sprintf(name, "hprof%d", i);
		sprintf(title, "title of hprof%d", i);
		hprof[i] = new TProfile(name, title, 100, 25, 30);
	}
	hprof[1]->Fill(27, 30);
	hprof[1]->Draw();



	/*TProfile *hprof = new TProfile("hprof", "Profile of ampl[1] versus phi[0]", 100, (float)(1 - 2)*(2.*TMath::Pi()) / 9., (float)(1 + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	TObjArray *face = new TObjArray(1000);
	for (int j = 0; j < 8; j++) {
		face->Add(hprof);
	}
	for (int i = 0; i < 8; i++) {
		TProfile *e;
		e = (TProfile*)face->At(i);
		e->Fill(0.1*i, 1.*i);
	}
	TProfile *e;
	e = (TProfile*)face->At(5);
	e->Draw();*/

	/*TCanvas *c1 = new TCanvas("c1", "A", 200, 10, 500, 300);
	c1->Divide(3, 3);
	TProfile  *hprof[8];
	for (int j = 0; j < 8; j++) {
		TProfile *hprof[j] = new TProfile(Form("hprof%d", j), "Profile of ampl[1] versus phi[0];phi;ach", 100, (float)(1 - 2)*(2.*TMath::Pi()) / 9., (float)(1 + 1)*(2.*TMath::Pi()) / 9., 0, 20);
	}
	c1->cd(1);
	hprof[1]->Fill(0.5, 12);
	//hprof[1]->SetTitle("; phi; ach");
	hprof[1]->Draw();*/
}



void go(){ 
	double x1[4] = { 0.1, 0.27, 0.33, 0.71 };
	double y1[4] = { 0.2, 0.17, 0.43, 0.61 };
	
//	copy1();
//	copy2();
	raspr();
//	m4(x1,y1);
  
	//pscp D:\programs\test.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/test.cpp
	//pscp pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/select_ach1.cpp C:\disk\select_ach1.cpp
}

