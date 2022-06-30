
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TChain.h"
#include "TString.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
using namespace std;



int copy() {

	const char *cond = "(cosm == 0) && (nc == 2) && (nn==0) && (charge[0] == 1) && (charge[1] == 1) && (energy[0]>=900) && (energy[1] >=900) && (nch == 9) && (theta[0] >=0.1) && (theta[0] <=3.0) && ((13. / tan(theta[0]) + z0[0]) > -18.) && ((13. / tan(theta[0]) + z0[0]) < 18.) && ((13. / tan(theta[1]) + z0[1]) > -18.) && ((13. / tan(theta[1]) + z0[1]) < 18.)";
	
	/*TChain ch("t1");
	ch.Add("/work/users/kladov/snd2k/R006-003/output/ntuples/ee/*.root");	
	ch.Merge("/work/users/kladov/snd2k/R006-003/output/ntuples/ee1/combined.root");*/
	TFile* filecomb = TFile::Open("/work/users/kladov/snd2k/R006-003/output/ntuples/ee1/combined.root");
	TTree* originalTree = (TTree*)filecomb->Get("t1");
	TFile* ouput = TFile::Open("/work/users/kladov/snd2k/R006-003/selected1/true_combined.root", "RECREATE");
	ouput->cd();
	TTree* selectedTree = originalTree->CopyTree(cond);
	selectedTree->Print();
	selectedTree->Write();
	filecomb->Close();
	ouput->Close();


	
	
	
	
	/*char s0[256];
	TString s;
	ifstream input("/work/users/kladov/snd2k/R006-003/fwk/filelist.txt");
	ofstream fout;
	fout.open("/work/users/kladov/snd2k/R006-003/fwk/filelistout.dat");   
	//for (int i = 1; i <= 120; i++) {
		input.getline(s0, 256, '\r\n');
	    s = s0;
	    s = s + "_wmix.root";	  	 
	    TString directory1 = "/work/users/kladov/snd2k/R006-003/output/ntuples/ee/";
	    TString directory2 = "/work/users/kladov/snd2k/R006-003/selected1/";
	    TString fullfile1 = directory1 + s;
	    TString fullfile2 = directory2 + "true_" + s;
	    cout << fullfile1 << endl;
	    const char *oldfile = (const char*)fullfile1;
	    const char *newfile = (const char*)fullfile2;
		fout << fullfile2 << endl;
		
		const char *oldfile = "/work/users/kladov/snd2k/R006-003/output/ntuples/ee/ee_bhwide_1000_10114-1_wmix.root";
		const char *newfile = "/work/users/kladov/snd2k/R006-003/selected1/true_ee_bhwide_1000_10114-1_wmix.root";
		TFile* file = TFile::Open(oldfile);
		TTree* originalTree = (TTree*)file->Get("t1");
		TFile* ouput = TFile::Open(newfile, "RECREATE");
		ouput->cd();
		TTree* selectedTree = originalTree->CopyTree(cond);
		selectedTree->Print();
		selectedTree->Write();
		ouput->Close();
		
		input.get();
	//}
		input.close();
		fout.close();
		*/











	    //strcpy(oldfile, fullfile1.c_str());
	    //strcpy(newfile, fullfile2.c_str());
	    //cout << newfile << endl;

	    /*char directory1[] = "/work/users/kladov/snd2k/R006-003/output/ntuples/ee/";
	    //char directory2[] = "/work/users/kladov/snd2k/R006-003/selected/";
	    //char filename[] = "ee_bhwide_990_13008-2_wmix.root";
	    char oldfile[128] = "a";
	    strcpy(oldfile,directory1);
	    strcat(oldfile,filename);

	    char newfile[128] = "b";
	    strcpy(newfile,directory2);
	    strcat(newfile,"true_");
	    strcat(newfile,filename);
	    */
		/*char myfile[] = "a";
		char directory1[] = "/work/users/kladov/snd2k/R006-003/output/ntuples/ee/";
		char directory2[] = "/work/users/kladov/snd2k/R006-003/selected1/";*/
		//char oldfile[256] = "a";
		//char newfile[256] = "b";

		/*strcpy(myfile, s.c_str());
		strcpy(oldfile, directory1);
		strcat(oldfile, myfile);
		strcat(oldfile, "_wmix.root");
		strcpy(newfile, directory2);
		strcat(newfile, "true_");
		strcat(newfile, myfile);
		strcat(newfile, "_wmix.root");*/
   
    return 0;
}
