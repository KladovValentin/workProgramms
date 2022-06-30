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








void go(){ 
	float energy[40];
	int nc, nn, cosm;


	TChain chain("t1");
	chain.Add("/work/users/konctbel/calibs/R007-001/output/ntuples/col/exp2017_860_col.root");

	chain.SetBranchAddress("energy", &energy);
	chain.SetBranchAddress("nc", &nc);
	chain.SetBranchAddress("nn", &nn);
	chain.SetBranchAddress("cosm", &cosm);


	for (int i = 1; i < 10; i++) {
		chain.GetEntry(i);

		cout << cosm << " " << nc << " " << " " << nn << endl;
		if ((cosm == 0) && (nc == 2) && (nn == 0)) {
			cout << energy[0] << endl;
		}
	}
   	//pscp C:\disk\kaketo.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/kaketo.cpp
}











 
