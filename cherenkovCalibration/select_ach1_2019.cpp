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

using namespace std;


float Goodphi[36] = { 0.08, 0.27, 0.33, 0.71,    0.78, 0.98, 1.04, 1.42,    1.48, 1.68, 1.75, 2.12,    2.21, 2.39, 2.45, 2.85,    2.92, 3.11, 3.17, 3.53,    3.61, 3.78, 3.85, 4.23,    \
	4.32, 4.5, 4.56, 4.92,    4.99, 5.17, 5.24, 5.61,    5.67, 5.87, 5.93, 6.28 };
float zobl[15] = { -11.7013, -10.0367, -8.47203, -6.90737, -5.34273, -3.77807, -2.21343, -0.648766, 0.9159, 2.48053, 4.0452, 5.60983, 7.1745, 8.73913, 9.6};
float egran[5] = { 600., 780., 900., 955., 1010. };
float phi[2], phis[2], mcphi[2], ach[40], achr[40], theta[2], z0[2], schr[40], tch[40], tchr[40], energy[40], beam, eton;
float d2phi[2], dphirho[2], d2rho[2], d2z0[2], d2cosTh[2], dz0cosTh[2], Dtheta[2], Dphi[2], energyerr[2];
int klmno[5], nch, run, nc, nn, cosm, eventtime, act, region[40];
int eventn;
double shwidth;
int calibsimple[224] = { 37361	, 37405	,37425	,37444	,38183	,38211	,38212	,38245	,38288	,38359	,
38422	,38463	,38500	,38518	,38518	,38554	,38601	,38619	,38621	,38639	,38665	,38686	,38729	,38769	,38797	,38800	,
38804	,38804	,38872	,38896	,38928	,38966	,39027	,39067	,39105	,39224	,39224	,39224	,39225	,39269	,39385	,39471	,
39527	,39588	,39599	,39635	,39734	,39767	,39816	,39866	,39912	,39934	,39980	,39980	,40023	,40028	,40053	,40091	,
40112	,40352	,40370	,40384	,40405	,40462	,40480	,40480	,40480	,40502	,40502	,40502	,40503	,40534	,40567	,40567	,
40567	,40567	,40590	,40590	,40590	,40622	,40655	,40679	,40704	,40716	,40742	,40765	,40765	,40767	,40783	,40808	,
40827	,40842	,40856	,40885	,40909	,40925	,40951	,40966	,40980	,41009	,41024	,41046	,41068	,41091	,41103	,41125	,
41150	,41180	,41198	,41220	,41236	,41257	,41279	,41303	,41320	,41346	,41374	,41403	,41425	,41446	,41477	,41499	,
41508	,41514	,41521	,41541	,41561	,41580	,41595	,41609	,41621	,41658	,41682	,41711	,41728	,41754	,41776	,41810	,
41836	,41866	,41886	,41911	,41936	,41952	,41976	,42002	,42014	,42034	,42055	,42078	,42143	,42178	,42209	,42238	,
42245	,42273	,42288	,42311	,42338	,42371	,42396	,42400	,42410	,42438	,42459	,42467	,42496	,42520	,42548	,42579	,
42602	,42628	,42654	,42669	,42700	,42714	,42752	,42765	,42772	,42789	,42807	,42916	,42932	,42963	,42989	,43013	,
43028	,43047	,43067	,43080	,43096	,43122	,43160	,43185	,43205	,43220	,43222	,43236	,43254	,43275	,43285	,43297	,
43314	,43334	,43351	,43383	,43408	,43430	,43452	,43493	,43526	,43563	,43587	,43657	,43692	,43747	,43787	,43814	,
43836	,43878	,43896	,43921	,43981	,44008 };

const char *cond = "(cosm == 0) && (nc == 2) && (nn==0) && (charge[0] == 1) && (charge[1] == 1) && (energy[0]/beam>=0.8) && (energy[1]/beam>=0.7) && (nch >= 9) && (act == 0) &&\
((theta[0] + theta[1] - 3.14159) <= 0.12) && ((theta[0] + theta[1] - 3.14159) >= -0.12) && (sin(theta[0]) >= 0.41) && \
((((phi[0] - phi[1]) - 3.14159 <= 0.1) && (((phi[0] - phi[1]) - 3.14159) >= -0.1)) || (((phi[0] - phi[1]) + 3.14159 <= 0.1) && (((phi[0] - phi[1]) + 3.14159) >= -0.1))) && \
((z0[0] - z0[1]) <= 1.0) && ((z0[0] - z0[1]) >= -1.0) && ((d0[0] - d0[1]) <= 0.5) &&  ((d0[0] - d0[1]) >= -0.5)";

const char* condkaons = "(cosm == 0) && (nc == 2) && (nn==0) && (charge[0] == 1) && (charge[1] == 1) && (energy[0]/beam<0.8) && (energy[1]/beam<0.7) && (nch >= 9) && \
(energy[0]/beam>0.2) && (energy[1]/beam>0.2)     &&\
((theta[0] + theta[1] - 3.14159) <= 0.12) && ((theta[0] + theta[1] - 3.14159) >= -0.12) && (sin(theta[0]) >= 0.41) && \
((((phi[0] - phi[1]) - 3.14159 <= 0.1) && (((phi[0] - phi[1]) - 3.14159) >= -0.1)) || (((phi[0] - phi[1]) + 3.14159 <= 0.1) && (((phi[0] - phi[1]) + 3.14159) >= -0.1))) && \
((z0[0] - z0[1]) <= 1.0) && ((z0[0] - z0[1]) >= -1.0) && ((d0[0] - d0[1]) <= 0.2) &&  ((d0[0] - d0[1]) >= -0.2)";
int mybin = 100;
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

struct EntriesInRun {
	int run;
	float energy;
	int entries;
};

double p3g0(double x, double xr, double sr, std::vector<double> p);
double aerogel(size_t n, double x, double* par);
double shifter(double x, double* par);
double p3g(double x, double xr, double sr, std::vector<double> p, double smin, double smax);
double scphi(double* x, double* par);


//for amptime function to draw vs time and not vs run
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
	return (size_t)(std::atoi(gSystem->GetFromPipe(cmd.Data())) - gStyle->GetTimeOffset() - 6 * 3600);
}

//for fit of amptime with an array of connected lines (probably, not needed, and certainly they shouldn't be connected as changes in amplitude can be caused by abrupt changes in detector)
double lineotrezok(int n, double x, double* par)
{
	double a = par[3 * n + 1];
	double b = par[3 * n + 2];
	return a * x + b;
}
double linesall(double* x, double* par)
{
	double xout = 1;
	int nmax = 7;
	//cout << nmax << endl;
	for (int n = 0; n < nmax; n++) {
		if ( (x[0] > par[3*n]) && (x[0] <= par[3*(n+1)]) )
			xout = lineotrezok(n, x[0], par);
	}
	return xout;
}

//drawing EMC deposition of particles in original col stream to determine selection criteria
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
	//TH1F* hist1 = (TH1F*)f->Get("energy1");
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

class CalibrExp {
public:
	TString basedirExp;
	TString basedirMod;
	string workingDir;
	string dirSelectedExp;
	string prefixSelectedExp;
	string dirSelectedMod;
	string prefixSelectedMod;
	float timeCut;
	double aerBord[2];
	double aerRLowBord;
	double aerWidth;
	double achCut;
	int entriesInGoodRun;
	double effInGoodRun;

	int MinimRun, MaximRun;
	vector<bool> goodrun;
	double truegran1[9][4];
	double truegran1Err[9][4];
	double gran1[4][9][14];
	double gran1Err[4][9][14];
	double pik[9][14];
	double meanshir[9][3];
	double truegran1m[9][4];
	double truegran1mErr[9][4];
	double meanshirm[9][3];
	double shift[9][3];
	double fitPar[9][14][21];
	double yfm[9][14][12];

	//expModComparison
	double ecounta1[9][14];
	double ecounta2[9][14];
	double ecountc1[14];
	double ecountc2[14];
	double ecountz1[9];
	double ecountz2[9];
	double counta1[9][14];
	double counta2[9][14];
	double countc1[14];
	double countc2[14];
	double countz1[9];
	double countz2[9];

	double coef[14][9];

	vector<double> bins;

	CalibrExp() {
		//specify basedirExp
		//mkdir basedirMod, workingDir, dirSelectedExp, dirSelectedMod
		//copy files /work/users/kladov/snd2k/R006-003/maindir/2017/*.dat and /work/users/kladov/snd2k/R006-003/maindir/2017/*.root to workingDir
		//change TimeCut (or TineIsGood function)

		//generate map with first block in go()
		//---setup the required soft for working with the database (if it is the first time) (can copy from konctbel guide for me) 
		//write map to database: 1) in file map.kumac find readmapx; 2) change there map path to workingDir/map.txt; 3) paw -> exec map#readmapx -> y -> exit; 4) check calibration clbixlist accmap MCTEST

		//find modeling files with conserved hits .mod.gz for runs in required range (from file workingDir/minmaxrun.dat)
		//execute listoffiles.cpp script getlist() (pass there basedirMod, directory with .mod.gz files and min/max run cuts) and copy the output lines (in file) in runrecobase
		//make sure that files fwk/simreco-neu.. and fwk/histcomon-col.... are correct
		//./runrecobase to start modeling processing

		//do h2root for them, command can be found in a console output of getlist()
		//compare the modeling with the experiment using second block in go()

		//check the time spectre and get the amp-time depentence and store it in workingDir/runAmp.dat via third block of functions in go()

		//finally, profiles can be found in workingDir/profiles.root; the corresponding map with area-regarding coefficients in a certain text form can be found in workingDir/map.txt
		//maybe, I will write a script to multiply the amplitudes in this map on the amp-time dependence

		//basedirExp = "/work/users/konctbel/calibs/R007-001/output/ntuples/MHAD2020/col/*col.root";	//directory with exp col stream files to which apply calibration
		//basedirExp = "/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2019/*col.root";			//directory with exp col stream files to which apply calibration
		basedirExp = "/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2019/eecalib/*col.root";	//directory with exp col stream files to which apply calibration
		basedirMod = "/work/users/kladov/snd2k/R007-001/2019/";										//directory with mod files after map creation

		workingDir = "/work/users/kladov/snd2k/R006-003/maindir/2019/";								//specify where to store temporary files
		dirSelectedExp = "/work/users/kladov/snd2k/R006-003/selected2019/";							//where to store selected Baba events
		dirSelectedMod = workingDir + "model/selected/";											//where to store selected Baba events for modeling
		prefixSelectedExp = "true1__";
		prefixSelectedMod = "true_";
		timeCut = 1100;
		//timeCut = 88;
		aerBord[0] = 11.7; //11.7
		aerBord[1] = 10.8; //10.8
		aerRLowBord = 10.5;
		aerWidth = 3.4; //3.4
		achCut = 70.;
		entriesInGoodRun = 2000;
		effInGoodRun = 0.81;
	}

	bool TimeIsGood(int signnumb);

	void readMinMaxRun() {
		ifstream fin1;
		fin1.open((workingDir + "minmaxrun.dat").c_str());
		fin1 >> MinimRun >> MaximRun;
		fin1.close();
	}
	void writeMinMaxRun() {
		ofstream fout1;
		fout1.open((workingDir + "minmaxrun.dat").c_str());
		fout1 << MinimRun << "	" << MaximRun << endl;
		fout1.close();
	}

	void readGoodRuns() {
		readMinMaxRun();
		ifstream fin1;
		fin1.open((workingDir + "goodRuns.dat").c_str());
		for (size_t i = 0; i < size_t(MaximRun - MinimRun); i++) {
			int grt = 0;
			fin1 >> grt;
			if (goodrun.size() <= i)
				goodrun.push_back((bool)grt);
			else
				goodrun[i] = (bool)grt;
			fin1.get();
		}
		fin1.close();
		cout << goodrun.size() << "	" << MaximRun - MinimRun << endl;
	}
	void writeGoodRuns() {
		ofstream fout1;
		fout1.open((workingDir + "goodRuns.dat").c_str());
		for (size_t i = 0; i < goodrun.size(); i++) {
			cout << (int)goodrun[i] << " ";
			fout1 << (int)goodrun[i] << endl;
		}
		fout1.close();
	}

	void readTrueGran1() {
		ifstream fin1;
		fin1.open((workingDir + "truegran1.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++)
				fin1 >> truegran1[i][j];
			fin1.get();
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++)
				fin1 >> truegran1Err[i][j];
			fin1.get();
		}
		fin1.close();
	}
	void writeTrueGran1() {
		ofstream fout1;
		fout1.open((workingDir + "truegran1.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++)
				fout1 << truegran1[i][j] << "	";
			fout1 << endl;
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++)
				fout1 << truegran1Err[i][j] << "	";
			fout1 << endl;
		}
		fout1.close();
	}

	void readFitPar() {
		ifstream fin1;
		fin1.open((workingDir + "FitPar.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++) {
				for (int k = 0; k < 21; k++)
					fin1 >> fitPar[i][j][k];
				fin1.get();
			}
			fin1.get();
		}
		fin1.close();
	}
	void writeFitPar() {
		ofstream fout1;
		fout1.open((workingDir + "FitPar.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++) {
				for (int k = 0; k < 21; k++)
					fout1 << fitPar[i][j][k] << "	";
				fout1 << endl;
			}
			fout1 << endl;
		}
		fout1.close();
	}

	void readGran1() {
		ifstream fin1;
		fin1.open((workingDir + "gran1.dat").c_str());
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 9; j++) {
				for (int k = 0; k < 14; k++){
					fin1 >> gran1[i][j][k];
				}
				fin1.get();
			}
			fin1.get();
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 9; j++) {
				for (int k = 0; k < 14; k++)
					fin1 >> gran1Err[i][j][k];
				fin1.get();
			}
			fin1.get();
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fin1 >> pik[i][j];
			fin1.get();
		}
		fin1.close();
	}
	void writeGran1() {
		ofstream fout1;
		fout1.open((workingDir + "gran1.dat").c_str());
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 9; j++) {
				for (int k = 0; k < 14; k++)
					fout1 << gran1[i][j][k] << "	";
				fout1 << endl;
			}
			fout1 << endl;
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 9; j++) {
				for (int k = 0; k < 14; k++)
					fout1 << gran1Err[i][j][k] << "	";
				fout1 << endl;
			}
			fout1 << endl;
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fout1 << pik[i][j] << "	";
			fout1 << endl;
		}
		fout1.close();
	}

	void readMeanshir() {
		ifstream fin1;
		fin1.open((workingDir + "meanshir.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 3; j++)
				fin1 >> meanshir[i][j];
			fin1.get();
		}
		fin1.close();
	}
	void writeMeanshir() {
		ofstream fout1;
		fout1.open((workingDir + "meanshir.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 3; j++)
				fout1 << fabs(meanshir[i][j]) << "	";
			fout1 << endl;
		}
		fout1.close();
	}

	void readTruegran1m() {
		ifstream fin1;
		fin1.open((workingDir + "truegran1m.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++)
				fin1 >> truegran1m[i][j];
			fin1.get();
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++)
				fin1 >> truegran1mErr[i][j];
			fin1.get();
		}
		fin1.close();
	}
	void writeTruegran1m() {
		ofstream fout1;
		fout1.open((workingDir + "truegran1m.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++)
				fout1 << truegran1m[i][j] << "	";
			fout1 << endl;
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++)
				fout1 << truegran1mErr[i][j] << "	";
			fout1 << endl;
		}
		fout1.close();
	}

	void readMeanshirm() {
		ifstream fin1;
		fin1.open((workingDir + "meanshirm.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 3; j++)
				fin1 >> meanshirm[i][j];
			fin1.get();
		}
		fin1.close();
	}
	void writeMeanshirm() {
		ofstream fout1;
		fout1.open((workingDir + "meanshirm.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 3; j++)
				fout1 << fabs(meanshirm[i][j]) << "	";
			fout1 << endl;
		}
		fout1.close();
	}

	void readShift() {
		ifstream fin1;
		fin1.open((workingDir + "shift.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 3; j++)
				fin1 >> shift[i][j];
			fin1.get();
		}
		fin1.close();
	}
	void writeShift() {
		ofstream fout1;
		fout1.open((workingDir + "shift.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 3; j++)
				fout1 << shift[i][j] << "	";
			fout1 << endl;
		}
		fout1.close();
	}

	void readComparisone() {
		ifstream fin1;
		fin1.open((workingDir + "ecounteri.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fin1 >> ecounta1[i][j];
			fin1.get();
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fin1 >> ecounta2[i][j];
			fin1.get();
		}
		for (int i = 0; i < 14; i++)
			fin1 >> ecountc1[i];
		fin1.get();
		for (int i = 0; i < 14; i++)
			fin1 >> ecountc2[i];
		fin1.get();
		for (int i = 0; i < 9; i++)
			fin1 >> ecountz1[i];
		fin1.get();
		for (int i = 0; i < 9; i++)
			fin1 >> ecountz2[i];
		fin1.get();
		fin1.close();
	}
	void writeComparisone() {
		ofstream fout1;
		fout1.open((workingDir + "ecounteri.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fout1 << ecounta1[i][j] << "	";
			fout1 << endl;
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fout1 << ecounta2[i][j] << "	";
			fout1 << endl;
		}
		for (int i = 0; i < 14; i++)
			fout1 << ecountc1[i] << "	";
		fout1 << endl;
		for (int i = 0; i < 14; i++)
			fout1 << ecountc2[i] << "	";
		fout1 << endl;
		for (int i = 0; i < 9; i++)
			fout1 << ecountz1[i] << "	";
		fout1 << endl;
		for (int i = 0; i < 9; i++)
			fout1 << ecountz2[i] << "	";
		fout1 << endl;
		fout1.close();
	}

	void readComparisonm() {
		ifstream fin1;
		fin1.open((workingDir + "mcounteri.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fin1 >> counta1[i][j];
			fin1.get();
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fin1 >> counta2[i][j];
			fin1.get();
		}
		for (int i = 0; i < 14; i++)
			fin1 >> countc1[i];
		fin1.get();
		for (int i = 0; i < 14; i++)
			fin1 >> countc2[i];
		fin1.get();
		for (int i = 0; i < 9; i++)
			fin1 >> countz1[i];
		fin1.get();
		for (int i = 0; i < 9; i++)
			fin1 >> countz2[i];
		fin1.get();
		fin1.close();
	}
	void writeComparisonm() {
		ofstream fout1;
		fout1.open((workingDir + "mcounteri.dat").c_str());
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fout1 << counta1[i][j] << "	";
			fout1 << endl;
		}
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 14; j++)
				fout1 << counta2[i][j] << "	";
			fout1 << endl;
		}
		for (int i = 0; i < 14; i++)
			fout1 << countc1[i] << "	";
		fout1 << endl;
		for (int i = 0; i < 14; i++)
			fout1 << countc2[i] << "	";
		fout1 << endl;
		for (int i = 0; i < 9; i++)
			fout1 << countz1[i] << "	";
		fout1 << endl;
		for (int i = 0; i < 9; i++)
			fout1 << countz2[i] << "	";
		fout1 << endl;
		fout1.close();
	}

	void readMcoef() {	//coefficients from current modeling database
		ifstream fin1;
		fin1.open((workingDir + "coef.dat").c_str());
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++)
				fin1 >> coef[i][j];
			fin1.get();
		}
		fin1.close();
	}
	void writeMcoef() {
		ofstream fout1;
		fout1.open((workingDir + "coef.dat").c_str());
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++)
				fout1 << coef[i][j] << "	";
			fout1 << endl;
		}
		fout1.close();
	}

	void readNcoef() {	//new coefficients from comparison with modeling to pass to a new database
		ifstream fin1;
		fin1.open((workingDir + "coefNew.dat").c_str());
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++)
				fin1 >> coef[i][j];
			fin1.get();
		}
		fin1.close();
	}
	void writeNcoef() {
		ofstream fout1;
		fout1.open((workingDir + "coefNew.dat").c_str());
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++)
				fout1 << coef[i][j] << "	";
			fout1 << endl;
		}
		fout1.close();
	}

	double recalculateZ(double zin, double theta, double ach, double pedz);
	void copy1();
	void copymod();
	void findrunborders();
	bool RunIsGood(int run);
	void findruns(int amount, double effThreshold);
	void zraspr();
	void linesforz(std::string name, const char* title1, const char* title2);
	void raspr();
	void rasprmod(bool transformPYN);
	void testforfit(string basefile, string pName);
	void ampltime();
	void bordersplot();
	void modelfile(bool transformPYN);
	vector<int> getCalibRuns(string calName, bool read);
	void drawamptime();
	void enpoint();
	void achspec();
	void zrasprmod();
	void timespectra();
	void drawtsp();
	void nvsphiorig();
	void phiminusphis();
	void checkmod();
	void checkbord();
	void dopeff();
	double compare1();
	double compare2(bool transformPYN);
	void compare();
	void compareAmpSpectr();
	void compareAmpSpectrG();
	void compareEffSpectr();
	void compareZAmpSpectr();
	void phiShift();
	void comparePhiShifts();
};


bool CalibrExp::TimeIsGood(int signnumb) {
	//return (((schr[signnumb] < 0.5) || ((schr[signnumb] > 0.5) && ((tchr[signnumb] > 87) && (tchr[signnumb] < 104)))) && (eventtime == 0));
	return (((schr[signnumb] < 0.5 && tchr[signnumb] > 400 ) || ((schr[signnumb] > 0.5) && (tchr[signnumb] < 10*timeCut) && (tchr[signnumb] > 87))));
	//return (((schr[signnumb] < 0.5) || ((schr[signnumb] > 0.5) && (tchr[signnumb] > timeCut))));
	//return (tchr[signnumb] < timeCut);
}

bool CalibrExp::RunIsGood(int run) {
	bool isrungood = (run - (MinimRun + 1)) < 0 ? false : goodrun[run - MinimRun - 1];
	TString dirSelectedExpT = (TString)dirSelectedExp;
	isrungood = isrungood && (run > 27250 || !dirSelectedExpT.EndsWith("2017/"));
	return isrungood;
}

double CalibrExp::recalculateZ(double zin, double theta, double ach, double pedz) {
	double amplitude = -1.0;
	int aerBordInd = -1;
	if (zin >= 0)
		aerBordInd = 1;
	else
		aerBordInd = 0;
	if (((aerBord[aerBordInd] - fabs(zin)) * fabs(tan(theta)) < aerWidth) && (fabs(aerBord[aerBordInd] - fabs(zin)) >= 0.01)) {

		amplitude = (ach - pedz) * aerWidth * fabs(cos(theta)) / fabs(aerBord[aerBordInd] - fabs(zin));
	}
	else if ((aerBord[aerBordInd] - fabs(zin)) * fabs(tan(theta)) >= aerWidth) {

		amplitude = (ach - pedz) * sin(theta);
	}
	return amplitude;
}

//copy only e+e- entries (cond above) from col stream, located in "basedir", to a new directory "dir", 
//corresponding to selected year, in order to not analyse all col entries every time in the following calculations
void CalibrExp::copy1() {
	vector<TString> ans;
	///base directory for files
	TString files = gSystem->GetFromPipe("ls " + basedirExp);

	Ssiz_t from = 0;
	
	while (files.Tokenize(elem, from, "\n")) {
		ans.push_back(elem);
	}
	size_t vector_size = ans.size();
	for (size_t i = 0; i < vector_size; i++) {
		str = ans[i].Copy();
		TString tok;
		Ssiz_t from1 = 1;
		vector<TString> tokens;
		TFile* filecomb = TFile::Open(str);
		TTree* originalTree = (TTree*)filecomb->Get("t1");

		while (str.Tokenize(tok, from1, "/")) {
			tokens.push_back(tok);
		}
		TString fname = tokens.back();
		cout << fname << endl;
		TFile* ouput = TFile::Open(dirSelectedExp + prefixSelectedExp + fname, "RECREATE");
		ouput->cd();
		TTree* selectedTree = originalTree->CopyTree(cond);
		selectedTree->Print();
		selectedTree->Write();
		ouput->Close();
		filecomb->Close();
	}
	/*for (int i = 98; i < vector_size; i++) {
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
	}*/
}

void CalibrExp::copymod() {
	vector<TString> ans;
	TString files = gSystem->GetFromPipe("ls " + basedirMod + "*0.root");

	Ssiz_t from = 0;

	while (files.Tokenize(elem, from, "\n")) {
		ans.push_back(elem);
	}
	size_t vector_size = ans.size();
	for (size_t i = 0; i < vector_size; i++) {
		str = ans[i].Copy();
		TString tok;
		Ssiz_t from1 = 1;
		vector<TString> tokens;

		TFile* filecomb = TFile::Open(str);
		TTree* originalTree = (TTree*)filecomb->Get("h1");

		while (str.Tokenize(tok, from1, "/")) {
			tokens.push_back(tok);
		}
		TString fname = tokens.back();
		cout << dirSelectedMod + "true_" + fname << endl;
		TFile* ouput = TFile::Open(dirSelectedMod + prefixSelectedMod + fname, "RECREATE");
		ouput->cd();
		TTree* selectedTree = originalTree->CopyTree(cond);
		selectedTree->Print();
		selectedTree->Write();
		ouput->Close();
		filecomb->Close();
	}
}

//find min and max run in the data to make distribution in this borders
void CalibrExp::findrunborders() {
	TChain chain("t1");
	std::string chadds = dirSelectedExp + prefixSelectedExp + "**exp00**.root";
	const char* chadd = chadds.c_str();
	chain.Add(chadd);
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	cout << "|";
	for (int i = 0; i < entries / 1000000; i++)
		cout << " ";
	cout << "|" << endl;
	cout << " ";
	chain.SetBranchAddress("run", &run);
	MinimRun = 100000;
	MaximRun = 0;

	int Count = 0, count1 = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		if (run < MinimRun)
			MinimRun = run;
		if (run > MaximRun)
			MaximRun = run;

		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			//cout << Form("obrabotano %d M entries", count1) << endl;
			cout << "a";
			Count = 0;
		}
	}
	cout << endl;
	MinimRun = MinimRun - 1;
	MaximRun = MaximRun + 1;
	cout << "run range:  " << MinimRun << " - " << MaximRun << endl;
	writeMinMaxRun();
}

//find runs that contain (events > "amount") and (mean efficiency > "effThreshold"), conditions on phi(+shir) and z
void CalibrExp::findruns(int amount, double effThreshold) {
	readMinMaxRun();
	readTrueGran1();
	readMeanshir();
	TChain chain("t1");
	std::string chadds = dirSelectedExp + prefixSelectedExp + "**exp00**.root";
	chain.Add(chadds.c_str());
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
	chain.SetBranchAddress("beam", &beam);
	
	TH1* histo = new TH1I("runs","runs", MaximRun - MinimRun, MinimRun, MaximRun);		//profile for amount of events in run

	TProfile* hprof[9];
	char name[20];
	char title[100];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprof%d", i);
		sprintf(title, "sens%d;eff;run", i + 1);
		hprof[i] = new TProfile(name, title, MaximRun - MinimRun, MinimRun, MaximRun);	//profile for mean efficiency in run 
	}

	int Count = 0, count1 = 0;
	double lgran[4];
	float ztr = 0.;
	float ztr1 = 0.;
	int scount = 0;
	int scount1 = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12. / tan(theta[0]) + z0[0];
		ztr1 = 12. / tan(theta[1]) + z0[1];
		histo->Fill(run);
		schr[nch] = 0.;	// to know that nch-1 is the last for 9's counter
		if ((ztr > zobl[0]) && (ztr < zobl[14])) {
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9; j++) {
				//wichcounter1[j] = scount1;
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					if (phi[1] < 1.)
						phi[1] = phi[1] + 2. * PI;
					while (1 + scount < nch) { scount += 1; }
				}
				else
					while (schr[1 + scount] > 0.5) { scount += 1; }

				lgran[0] = truegran1[j][0] + 3. * meanshir[j][0];
				lgran[1] = truegran1[j][1] - 3. * meanshir[j][1];// -0.045;
				lgran[2] = truegran1[j][1] + 3. * meanshir[j][1];// +0.045;
				lgran[3] = truegran1[j][3] - 3. * meanshir[j][2];// -0.03;

				if ((((phi[0] > lgran[0]) && (phi[0] < lgran[1])) || ((phi[0] > lgran[2]) && (phi[0] < lgran[3]))) && (ztr > zobl[0]) && (ztr < zobl[14])) {
					cout << schr[1 + scount1] << endl;
					hprof[j]->Fill(run, schr[1 + scount1]);
				}
				if ((((phi[1] > lgran[0]) && (phi[1] < lgran[1])) || ((phi[1] > lgran[2]) && (phi[1] < lgran[3]))) && (ztr1 > zobl[0]) && (ztr1 < zobl[14]))
					hprof[j]->Fill(run, schr[1 + scount1]);
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

	TFile* MyFile = new TFile((workingDir + "runsi.root").c_str(), "RECREATE");
	for(int j = 0; j<9;j++){
		sprintf(title, "sens%d_eff_run", j + 1);
		hprof[j]->Write(title);
	}
	histo->Write("runsnumb");
	MyFile->Close();

	int countt1 = 0;
	int k = 0;
	for (int i = 2; i <= MaximRun - MinimRun; i++) {
		k = 0;
		for (int j = 0; j < 9; j++) {
			//cout << "eff	" << hprof[j]->GetBinContent(i) << endl;
			if (hprof[j]->GetBinContent(i) > effThreshold)
				k += 1;
		}
		//cout << "number of entries	" << histo->GetBinContent(i) << endl;
		if (k == 9 && (histo->GetBinContent(i) > amount)) {
			goodrun.push_back(1);
			countt1 += 1;
		}
		else
			goodrun.push_back(0);
	}
	writeGoodRuns();
	cout << "amount of good runs	" << countt1 << endl;
}

//make profile with mean ach vs z on the inner cillinder and medium cillinder with recalculation on borders according to intersection line segment lenght and without
void CalibrExp::zraspr() {
	readGoodRuns();
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
	readTruegran1m();
	TChain chain("t1");
	std::string chadds = dirSelectedExp + prefixSelectedExp + "**exp00**.root";
	chain.Add(chadds.c_str());
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
	chain.SetBranchAddress("eventtime", &eventtime);

	int scount = 0;
	int scount1 = 0;
	char name[20];
	char title[100];

	TH1* h[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprofped%d", i + 1);
		sprintf(title, "ped,sensor%d;ach", i + 1);
		h[i] = new TH1F(name, title, 100, -1, 1);				//hists for pedestals, 9 counters 
	}
	
	TProfile* hprof[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprof%d", i + 1);
		sprintf(title, "sensor%d;z;ach", i + 1);
		hprof[i] = new TProfile(name, title, 300, -15, 15);		//profiles for ach vs ztr(12.0), 9 counters, with recalculation at borders
	}
	TProfile* hprof1[2];
	hprof1[0] = new TProfile("hpprofl", "all;ztr;ach", 300, -15, 15);
	hprof1[1] = new TProfile("hpprofr", "all;ztr;ach", 300, -15, 15);

	TProfile* hprof0[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprof0%d", i + 1);
		sprintf(title, "sensor0%d;z;ach", i + 1);
		hprof0[i] = new TProfile(name, title, 300, -15, 15);	//profiles for ach vs ztr(12.0), 9 counters, without recalculation
	}
	TProfile* hprof01[2];
	hprof01[0] = new TProfile("hpprof0l", "all;ztr;ach", 300, -15, 15);
	hprof01[1] = new TProfile("hpprof0r", "all;ztr;ach", 300, -15, 15);


	TFile* f = new TFile((workingDir + "pedestz.root").c_str());
	for (int i = 0; i < 9; i++) {
		TH1* hped = (TH1F*)f->Get(Form("ped_forz,sensor%d",i + 1));		//info about pedestals
		pedz[i] = hped->GetMean(1);
	}
	TH1* hr = new TH1F("histr", "recount right;z,cm", 100, 5., 15.);	//z, where ach is recalculated
	TH1* hl = new TH1F("histl", "recount left;z,cm", 100, -16., -6.);

	int Count = 0, count1 = 0;
	float ztr[2];
	ztr[0] = 0.; ztr[1] = 0.;
	float zin[2];
	zin[0] = 0.; zin[1] = 0.;
	double lgran[4];
	int maxampla = 0;
	double maxach = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr[0] = 12.0 / tan(theta[0]) + z0[0];
		ztr[1] = 12.0 / tan(theta[1]) + z0[1];
		zin[0] = aerRLowBord / tan(theta[0]) + z0[0];
		zin[1] = aerRLowBord / tan(theta[1]) + z0[1];
		schr[nch] = 0.;
		TString dirSelectedExpT = (TString)dirSelectedExp;
		if ((!RunIsGood(run)) || (eventtime != 0)/* && (run > 39420 && run < 39580) run < 39350*/) {
			continue;
		}
		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			if (j == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
				if (phi[1] < 1.)
					phi[1] = phi[1] + 2. * PI;
				//while (1 + scount < nch) { scount += 1; }
			}
			//else
				while (schr[1 + scount] > 0.5) { scount += 1; }

			lgran[0] = max(truegran1[j][0], truegran1m[j][0]) + 3. * max(meanshir[j][0], meanshirm[j][0]);
			lgran[1] = min(truegran1[j][1], truegran1m[j][1]) - 3. * max(meanshir[j][1], meanshirm[j][1]);// -0.045;
			lgran[2] = max(truegran1[j][1], truegran1m[j][1]) + 3. * max(meanshir[j][1], meanshirm[j][1]);// +0.045;
			lgran[3] = min(truegran1[j][3], truegran1m[j][3]) - 3. * max(meanshir[j][2], meanshirm[j][2]);// -0.03;


			maxach = -1000.;
			for (int f = scount1; f <= scount; f++) {
				if (ach[f] > maxach) {
					maxach = ach[f];
					maxampla = f;
				}
			}
			int aerBordInd = 0;

				
			for (int prtc = 0; prtc < 2; prtc++) {

				if (!(((phi[prtc] > lgran[0]) && (phi[prtc] < lgran[1])) || ((phi[prtc] > lgran[2]) && (phi[prtc] < lgran[3]))) || !TimeIsGood(maxampla))
					continue;
					
				//for pedestals
				h[j]->Fill(ach[scount1]);

				//for amplitude vs z
				//fill 0 hists with standart sin
				hprof01[prtc]->Fill(zin[prtc], (ach[maxampla] - pedz[j]) * sin(theta[prtc]));
				hprof0[j]->Fill(zin[prtc], (ach[maxampla] - pedz[j]) * sin(theta[prtc]));

				//fill with recalculated
				double fillAmplitude = recalculateZ(zin[prtc], theta[prtc], ach[maxampla], pedz[j]);
				if (fillAmplitude == -1.0 || ach[maxampla] > achCut)
					continue;

				hprof[j]->Fill(zin[prtc], fillAmplitude);
				hprof1[prtc]->Fill(zin[prtc], fillAmplitude);

			}
			scount += 1;
			scount1 = scount;
		
		}
		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("obrabotano %d M entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}
	TF1* fn = new TF1("scphi", "([0]+[1]*x)", -12, -7);
	fn->SetParameters(3, 0.2);
	hprof01[1]->Draw();
	hprof01[1]->Fit("scphi","","",-11,-9);
	cout << -fn->GetParameter(0) / fn->GetParameter(1) << endl;
	TF1* fn1 = new TF1("scphi1", "([0]+[1]*x)+[2]*exp(-((x-[3])/[4])**2)", 8, 12);
	fn1->SetParameters(5, -0.3, 1, 10.5, 1);
	hprof01[1]->Fit("scphi1", "", "", 8.5, 10.8);
	cout << -fn1->GetParameter(0) / fn1->GetParameter(1) << endl;


	TFile* MyFile = new TFile((workingDir + "pedestz.root").c_str(), "RECREATE");
	for (int j = 0; j < 9; j++) {
		sprintf(title, "ped_forz,sensor%d", j + 1);
		h[j]->Write(title);
	}
	MyFile->Close();

	TFile* MyFile1 = new TFile((workingDir + "zprofiles.root").c_str(), "RECREATE");
	for (Int_t i = 0; i < 9; i++) {
		sprintf(title, "sensor%d",i + 1);
		hprof[i]->Write(title);
		sprintf(title, "sensor0%d", i + 1);
		hprof0[i]->Write(title);
	}
	hprof1[0]->Write("allsensorsr");
	hprof1[1]->Write("allsensorsl");
	hprof01[0]->Write("allsensors0l");
	hprof01[1]->Write("allsensors0r");
	hr->Write("zinr");
	hl->Write("zinl");
	MyFile1->Close();	
}

//same as zraspr but without conditions on run and nch is always 9, need to rewrite or, better, unite it with zraspr
void CalibrExp::zrasprmod() {
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
	readTruegran1m();
	TChain chain("h1");
	std::string chadds = dirSelectedMod + prefixSelectedMod + "ee**.root";
	chain.Add(chadds.c_str());

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
	float ztr[2] = { 0.,0. };
	float zin[2] = { 0.,0. };
	double lgran[4];
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr[0] = 12.0 / tan(theta[0]) + z0[0];
		ztr[1] = 12.0 / tan(theta[1]) + z0[1];
		zin[0] = aerRLowBord / tan(theta[0]) + z0[0];
		zin[1] = aerRLowBord / tan(theta[1]) + z0[1];
		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			if (j == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
				if (phi[1] < 1.)
					phi[1] = phi[1] + 2. * PI;
			}

			lgran[0] = max(truegran1[j][0], truegran1m[j][0]) + 3. * max(meanshir[j][0], meanshirm[j][0]);
			lgran[1] = min(truegran1[j][1], truegran1m[j][1]) - 3. * max(meanshir[j][1], meanshirm[j][1]);// -0.045;
			lgran[2] = max(truegran1[j][1], truegran1m[j][1]) + 3. * max(meanshir[j][1], meanshirm[j][1]);// +0.045;
			lgran[3] = min(truegran1[j][3], truegran1m[j][3]) - 3. * max(meanshir[j][2], meanshirm[j][2]);// -0.03;

			pedz[j] = 0.;
			int aerBordInd = 0;

			for (int prtc = 0; prtc < 2; prtc++) {

				if (!(((phi[prtc] > lgran[0]) && (phi[prtc] < lgran[1])) || ((phi[prtc] > lgran[2]) && (phi[prtc] < lgran[3]))))
					continue;

				//fill 0 hists with standart sin
				hprof01->Fill(zin[prtc], (ach[j] - 0.0) * sin(theta[prtc]));
				hprof0[j]->Fill(zin[prtc], (ach[j] - 0.0) * sin(theta[prtc]));

				//fill with recalculated
				double fillAmplitude = recalculateZ(zin[prtc], theta[prtc], ach[j], 0);
				if (fillAmplitude == -1.0 || ach[j] > achCut)
					continue;
				hprof[j]->Fill(zin[prtc], fillAmplitude);
				hprof1->Fill(zin[prtc], fillAmplitude);

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


	TFile* MyFile1 = new TFile((workingDir + "zprofilesmod.root").c_str(), "RECREATE");
	for (Int_t i = 0; i < 9; i++) {
		sprintf(title, "sensor%d", i + 1);
		hprof[i]->Write(title);
		sprintf(title, "sensor0%d", i + 1);
		hprof0[i]->Write(title);
	}
	hprof1->Write("allsensors");
	hprof01->Write("allsensors0");
	hr->Write("zinr");
	hl->Write("zinl");
	MyFile1->Close();
}

//draw z profiles from one file "name" with titles 1 and 2
void CalibrExp::linesforz(std::string name, const char* title1, const char* title2) {
	TFile* f = new TFile((workingDir + name).c_str());

	TProfile* hprof = (TProfile*)f->Get(title1);
	TProfile* hprof0 = (TProfile*)f->Get(title2);
	//TH1* hr = (TH1F*)f->Get("zinr");
	//TH1* hl = (TH1F*)f->Get("zinl");
	//hr->SetLineColor(kRed);
	//hl->SetLineColor(kRed);
	hprof->SetTitle("amplitude distribution by z;z,cm;amplitude,pe");
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
	TLine* line11 = new TLine(-aerBord[0], 0, -aerBord[0], 7);
	TLine* line12 = new TLine(aerBord[1], 0, aerBord[1], 7);
	line11->SetLineColor(kGreen);
	line12->SetLineColor(kGreen);
	line11->Draw();
	line12->Draw();
	//hr->DrawNormalized("same", 50);
	//hl->DrawNormalized("same", 50);
	hprof0->Draw("same");
	c->Update();
}

//cleaned, maybe get parameters for recalculation on borders from zraspr fit 
//dont take entries with ztr > borders, because it ruin statistics - make big corner influence on first zone
void CalibrExp::raspr() {
	/// <summary>
	/// ~400 for fit, normal sim statistic - ~100 bins
	/// </summary>
	
	readGoodRuns();
	const char* treename = "t1";
	TChain chain(treename);

	std::string chadds = dirSelectedExp + prefixSelectedExp + "**exp00**.root";
	chain.Add(chadds.c_str());

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
	chain.SetBranchAddress("eventtime", &eventtime);

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;

	//raspr po phi
	TProfile* hprof[14][9];
	char name[20];
	char title[100];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {	
			sprintf(name, "hprofrp%d", 9*j + i + 1);
			sprintf(title, "zobl%d,sensor%d;phi;ach", j + 1,i + 1);
			hprof[j][i] = new TProfile(name,title, 200, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
			//hprof[j][i] = new TProfile(name,title, 100, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
			//hprof[j][i] = new TProfile(name,title, 4800, 0, 7.0);
		}
	}
	TProfile* hprofMean[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprofrpZm%d", i + 1);
		sprintf(title, "sensor%d;phi;ach",i + 1);
		hprofMean[i] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
	}
	//eff vs phi
	TProfile* hprofeff[14][9];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "hprofeff%d", 9 * j + i + 1);
			sprintf(title, "zobl%d,sensor%d;phi;eff", j + 1, i + 1);
			hprofeff[j][i] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
		}
	}
	TProfile* hprofEffMean[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprofeffZm%d", i + 1);
		sprintf(title, "sensor%d;phi;ach", i + 1);
		hprofEffMean[i] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
	}

	//hists for zasel po phi
	TH1* hza[14][9];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "hzas%d", 9*j + i + 1);
			sprintf(title, "zas,zobl%d,sensor%d;ach", j + 1, i + 1);
			hza[j][i] = new TH1F(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
		}
	}

	TProfile* hprof0 = new TProfile("tempamp", "tempamp counter 8 ", 800, (float)(7 - 1) * (2. * PI) / 9., (float)(7 + 2) * (2. * PI) / 9.);
	TProfile* hprofeff0 = new TProfile("tempeff", "tempeff count 8", 100, (float)(7 - 1) * (2. * PI) / 9., (float)(7 + 2) * (2. * PI) / 9.);

	//hists for pedestal integr po phi
	TH1* h[14][9];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "hprof%d", 9*j + i + 1);
			sprintf(title, "ped,zobl%d,sensor%d;ach", j + 1, i + 1);
			h[j][i] = new TH1F(name, title, 100, -1, 1);
		}
	}

	//info ab ped integrated phi
	TFile* f00 = new TFile((workingDir + "phintegr_pedest.root").c_str());
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			TH1* hped = (TH1F*)f00->Get(Form("ped_phintegr_zobl%d,sensor%d", j + 1, i + 1));
			ped[j][i] = hped->GetMean(1);
		}
	}
	TFile* f0 = new TFile((workingDir + "pedestz.root").c_str());
	for (int i = 0; i < 9; i++) {
		TH1* hped1 = (TH1F*)f0->Get(Form("ped_forz,sensor%d", i + 1));
		pedz[i] = hped1->GetMean(1);
	}

	int Count = 0, count1 = 0;
	float ztr[2] = { 0.,0. };
	float zin[2] = { 0.,0. };
	int maxampla = 0;
	double maxach = 0.;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//cout << beam << endl;
		ztr[0] = 12.0 / tan(theta[0]) + z0[0];
		ztr[1] = 12.0 / tan(theta[1]) + z0[1];
		zin[0] = aerRLowBord / tan(theta[0]) + z0[0];
		zin[1] = aerRLowBord / tan(theta[1]) + z0[1];
		schr[nch] = 0.;
		if ((!RunIsGood(run))/* || eventtime != 0*/) {
			continue;
		}
		scount = 0;     //2017: 4001, 25500, 29501    ; 2018: 6001, 29800, 35801       ; 2019: 38000-43001       ; 2020: 45040-46541
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			//wichcounter1[j] = scount1;
				//counters are shifted in oreder to determine first one and direction of count, I think. Because of this 9th counter have the tail of signal around 0, I transfer it on >2pi
				//nch is >= 9 in 2017+, so for one counter I choose from scount1 to scount 2 the highest amplitude 
			if (j == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
				if (phi[1] < 1.)
					phi[1] = phi[1] + 2. * PI;
				while (1 + scount < nch) { scount += 1; }
			}
			else
				while (schr[1 + scount] > 0.5) { scount += 1; }
			
			maxach = -100.;
			maxampla = scount1;
			for (int f = scount1; f <= scount; f++) {
				if (ach[f] > maxach) {
					maxach = ach[f];
					maxampla = f;
				}
			}
			
			//ach*sin ot phi 
				//parameters: good t in timeisgood, aerBord - assumed left and right borders of aerogel, aerWidth - width of aerogel to find intersection with its end, achCut - treshold to prevent random high signals
			int aerBordInd = 0;
			int zoblInd = -1;
			
			for (size_t pind = 0; pind < 2; pind++) {
				zoblInd = -1;
				for (int i = 0; i < 14; i++) {
					if ((ztr[pind] >= zobl[i]) && (ztr[pind] < zobl[i + 1]))
						zoblInd = i;
				}
				//if(tchr[scount1]<200){
					//	maxampla = scount1;
					//}
					//ach*sin ot phi 
				if ((phi[pind] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[pind] < (float)(j + 2) * (2. * PI) / 9.) && (TimeIsGood(maxampla)) && (TimeIsGood(scount1)) && zoblInd != -1) {

					//hprof[zoblInd][j]->Fill(phi[pind], (ach[maxampla] - ped[zoblInd][j]) * sin(theta[pind]));

					//double fillAmplitude = recalculateZ(zin[pind], theta[pind], ach[maxampla], ach[scount1]);
					double fillAmplitude = recalculateZ(zin[pind], theta[pind], ach[maxampla], ped[zoblInd][j]);
					if (fillAmplitude != -1.0 && ach[maxampla] < achCut && (TimeIsGood(maxampla))) {
						hprof[zoblInd][j]->Fill(phi[pind], fillAmplitude);
						if (zoblInd > 2 && zoblInd < 12)
							hprofMean[j]->Fill(phi[pind], fillAmplitude);
					}

					//efficiency
					double effWhatFill = 0.;
					if (ach[maxampla] >= 0.2)
						//if (schr[scount1+1] >= 1)
						effWhatFill = 1;
					//else if (ach[maxampla] - ach[scount1] < 0.2)
					else
						//if (schr[scount1+1] < 1)
						effWhatFill = 0;
					hprofeff[zoblInd][j]->Fill(phi[pind], effWhatFill);
					if (zoblInd > 2 && zoblInd < 12)
						hprofEffMean[j]->Fill(phi[pind], effWhatFill);

					//pedestals
					h[zoblInd][j]->Fill(ach[scount1]);
				}
				
				//ped int po phi
				//if ((phi[pind] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[pind] < (float)(j + 2) * (2. * PI) / 9.) && zoblInd != -1 && (TimeIsGood(maxampla))) {
				if (((phi[pind] >= (float)(j + 1.5) * (2. * PI) / 9. && phi[pind] < (float)(j + 2) * (2. * PI) / 9.) || (phi[pind] >= (float)(j - 1.0) * (2. * PI) / 9. && phi[pind] < (float)(j - 0.5) * (2. * PI) / 9.)) && zoblInd != -1 && (TimeIsGood(maxampla))) {
					if (((aerBord[aerBordInd] - fabs(zin[pind])) * fabs(tan(theta[pind])) < aerWidth) && ((aerBord[aerBordInd] - fabs(zin[pind])) >= 0.01)) {
						//h[zoblInd][j]->Fill(ach[maxampla]); //scount1
					}
					else if ((aerBord[aerBordInd] - fabs(zin[pind])) * fabs(tan(theta[pind])) >= aerWidth) {
						//h[zoblInd][j]->Fill(ach[maxampla]); //scount1
					}
				}

				//zasel po phi
				if ((phi[pind] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[pind] < (float)(j + 2) * (2. * PI) / 9.) && zoblInd != -1) {
					hza[zoblInd][j]->Fill(phi[pind]);
				}
			}

			scount += 1;
			scount1 = scount;
		}
		
		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("obrabotano %d M entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}

	TFile* MyFile = new TFile((workingDir + "profiles.root").c_str(), "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "zobl%d,sensor%d", m + 1, i + 1);
			hprof[m][i]->Write(title);
		}
	}
	for (Int_t i = 0; i < 9; i++) {
		sprintf(title, "mean,sensor%d", i + 1);
		hprofMean[i]->Write(title);
	}
	MyFile->Close();

	TFile* MyFileEff = new TFile((workingDir + "profilesEff.root").c_str(), "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "zobl%d,sensor%d", m + 1, i + 1);
			hprofeff[m][i]->Write(title);
		}
	}
	for (Int_t i = 0; i < 9; i++) {
		sprintf(title, "mean,sensor%d", i + 1);
		hprofEffMean[i]->Write(title);
	}
	MyFileEff->Close();

	TFile* MyFileZasel = new TFile((workingDir + "Zaselphi.root").c_str(), "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "zobl%d,sensor%d", m + 1, i + 1);
			hza[m][i]->Write(title);
		}
	}
	MyFileZasel->Close();

	TFile* MyFile1 = new TFile((workingDir + "phintegr_pedest.root").c_str(), "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "ped_phintegr_zobl%d,sensor%d", m + 1, i + 1);
			h[m][i]->Write(title);
		}
	}
	MyFile1->Close();

	cout << "done!" << endl;
}

//in modeling time and runs are always good, but phi distribution can be shifted, 
//because of little statistics I only make 9 distributions for counters with mean by z amplitude in z[1]<z<z[13] to determine borders in fit
void CalibrExp::rasprmod(bool transformPYN) {
	readTrueGran1();
	readTruegran1m();

	TChain chain("h1");
	std::string chadds = dirSelectedMod + prefixSelectedMod + "ee**.root";
	chain.Add(chadds.c_str());
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
	chain.SetBranchAddress("mcphi", &mcphi);

	int scount = 0;
	int scount1 = 0;

	//raspr po phi
	TProfile* hprof[14][9];
	char name[20];
	char title[100];
	for (int m = 0; m < 14; m++) {
		for (int j = 0; j < 9; j++) {
			sprintf(name, "hprofmod%d", m*9 + j + 1);
			sprintf(title, "zobl%d,counter%d;phi;ach", m + 1, j + 1);
			hprof[m][j] = new TProfile(name, title, 200, (float)(j - 1.) * (2. * PI) / 9., (float)(j + 2.) * (2. * PI) / 9.);
		}
	}
	TProfile* hprofMean[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprofrpZm%d", i + 1);
		sprintf(title, "sensor%d;phi;ach", i + 1);
		hprofMean[i] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
	}
	//eff vs phi
	TProfile* hprofeff[14][9];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "hprofeffmod%d", 9 * j + i + 1);
			sprintf(title, "zobl%d,sensor%d;phi;eff", j + 1, i + 1);
			hprofeff[j][i] = new TProfile(name, title, 100, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
		}
	}
	TProfile* hprofEffMean[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprofeffZm%d", i + 1);
		sprintf(title, "sensor%d;phi;ach", i + 1);
		hprofEffMean[i] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
	}

	int Count = 0, count1 = 0;
	float ztr[2] = { 0.,0. };
	float zin[2] = { 0.,0. };
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//cout << beam << endl;
		ztr[0] = 12.0 / tan(theta[0]) + z0[0];
		ztr[1] = 12.0 / tan(theta[1]) + z0[1];
		zin[0] = aerRLowBord / tan(theta[0]) + z0[0];
		zin[1] = aerRLowBord / tan(theta[1]) + z0[1];

		/*if(fabs(mcphi[0]-phi[0])<0.2 && fabs(mcphi[1]-phi[1])<0.2){
			phi[0] = mcphi[0];
			phi[1] = mcphi[1];
		}
		else if(fabs(mcphi[0]-phi[1])<0.2 && fabs(mcphi[1]-phi[0])<0.2){
			phi[0] = mcphi[1];
			phi[1] = mcphi[0];
		}*/

		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			//wichcounter1[j] = scount1;
			if (j == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
				if (phi[1] < 1.)
					phi[1] = phi[1] + 2. * PI;
			}
			int aerBordInd = 0;
			int zoblInd = -1;
			for (size_t pind = 0; pind < 2; pind++) {
				if (zin[pind] >= 0)
					aerBordInd = 1;
				else
					aerBordInd = 0;
				zoblInd = -1;
				for (int i = 0; i < 14; i++) {
					if ((ztr[pind] >= zobl[i]) && (ztr[pind] < zobl[i + 1]))
						zoblInd = i;
				}
				
				if ((phi[pind] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[pind] < (float)(j + 2) * (2. * PI) / 9.) && zoblInd != -1) {
					//hprof[zoblInd][j]->Fill(phi[pind], ach[j] * sin(theta[pind]));
					// 
					//transform angle
					double fillPhiWhat = phi[pind];
					double newPhi = 0;
					if (phi[pind] < truegran1m[j][1])
						newPhi = truegran1[j][0] + (phi[pind] - truegran1m[j][0]) * (truegran1[j][1] - truegran1[j][0]) / (truegran1m[j][1] - truegran1m[j][0]);
					else if (phi[pind] >= truegran1m[j][1])
						newPhi = truegran1[j][1] + (phi[pind] - truegran1m[j][1]) * (truegran1[j][3] - truegran1[j][1]) / (truegran1m[j][3] - truegran1m[j][1]);
					if (transformPYN)
						fillPhiWhat = newPhi;

					double fillAmplitude = recalculateZ(zin[pind], theta[pind], ach[j], 0);
					if (fillAmplitude != -1.0 && ach[j] < achCut) {
						hprof[zoblInd][j]->Fill(fillPhiWhat, fillAmplitude);
						if (zoblInd > 2 && zoblInd < 12)
							hprofMean[j]->Fill(fillPhiWhat, fillAmplitude);
					}

					//eff
					double effWhatFill = 0.;
					if (ach[j] >= 0.2)
					//if (tch[j] < 900)
						effWhatFill = 1.;
					if (ach[j] < 0.2)
					//if (tch[j] > 900)
						effWhatFill = 0.;
					hprofeff[zoblInd][j]->Fill(fillPhiWhat, effWhatFill);
					if (zoblInd > 2 && zoblInd < 12)
						hprofEffMean[j]->Fill(fillPhiWhat, effWhatFill);
				}
			}
			scount += 1;
			scount1 = scount;
		}
		

		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("obrabotano %d M entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}


	TFile* MyFile = new TFile((workingDir + "profilesmmod.root").c_str(), "RECREATE");
	for (int i = 0; i < 14; i++) {
		for (int m = 0; m < 9; m++) {
			sprintf(title, "zobl%d,sensor%d", i + 1, m + 1);
			hprof[i][m]->Write(title);
		}
	}
	for (Int_t i = 0; i < 9; i++) {
		sprintf(title, "mean,sensor%d", i + 1);
		hprofMean[i]->Write(title);
	}
	MyFile->Close();

	TFile* MyFileEff = new TFile((workingDir + "profilesEffMod.root").c_str(), "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "zobl%d,sensor%d", m + 1, i + 1);
			hprofeff[m][i]->Write(title);
		}
	}
	for (Int_t i = 0; i < 9; i++) {
		sprintf(title, "mean,sensor%d", i + 1);
		hprofEffMean[i]->Write(title);
	}
	MyFileEff->Close();
}

//attempts to test registration efficiency of kaons
void CalibrExp::dopeff() {
	readGoodRuns();
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


	TH1* histdpek = new TH1I("histdpek", "do_por_eff_kaon_exp;eff", 100, -1., 2.);

	int Count = 0, count1 = 0;
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
		if ((ztr > zobl[0]) && (ztr < zobl[14]) && (ztr1 > zobl[0]) && (ztr1 < zobl[14]) && (run > 27225) && RunIsGood(run) && (energy[0]/beam > 0.3) && (energy[1]/beam > 0.3) && (energy[0]/beam < 0.7) && (energy[1]/beam < 0.6)) {
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

void CalibrExp::ampltime() {
	readGoodRuns();
	readTrueGran1();
	readMeanshir();
	TChain chain("t1");

	std::string chadds = dirSelectedExp + prefixSelectedExp + "**exp00**.root";
	chain.Add(chadds.c_str());

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("achr", &achr);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("schr", &schr);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("tchr", &tchr);
	chain.SetBranchAddress("region", &region);
	chain.SetBranchAddress("beam", &beam);
	chain.SetBranchAddress("act", &act);
	chain.SetBranchAddress("eventtime", &eventtime);

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	int scount = 0;
	int scount1 = 0;

	//profiles for ach or eff
	TProfile* ha[10];
	char name[20];
	char title[100];
	for (int i = 0; i < 10; i++) {
		sprintf(name, "hprof%d", i);
		sprintf(title, "sens%d;run;ach,phe", i);
		if (bins.size() > 10)
			ha[i] = new TProfile(name,title, bins.size()-2, &bins[0]); 
		else
			ha[i] = new TProfile(name,title, (MaximRun - MinimRun) / 10, MinimRun, MaximRun);  //2017: 4001, 25500, 29501    ; 2018: 6001, 29800, 35801       ; 2019: 38000-43001       ; 2020: 45040-46541
	}
	TProfile* he[10];
	for (int i = 0; i < 10; i++) {
		sprintf(name, "hprofeff%d", i);
		sprintf(title, "sens%d;run;eff", i);
		if (bins.size() > 10)
			he[i] = new TProfile(name, title, bins.size() - 2, &bins[0]);
		else
			he[i] = new TProfile(name, title, (MaximRun - MinimRun)/10, MinimRun, MaximRun);
	}
	//ach spectre to find out mean value
	TH1* ampCounter[9];
	for (int j = 0; j < 9; j++) {
		sprintf(name, "adistr%d", j + 1);
		sprintf(title, "amplitude spectrum, sensor%d;ach", j + 1);
		ampCounter[j] = new TH1F(name, title, 5000, -1, 499);
	}

	
	TH1* phidistrgood = new TH1F("phidistrgood","phidistrgood", 100, (float)(7 - 1) * (2. * PI) / 9., (float)(7 + 2) * (2. * PI) / 9.);
	TH1* phidistrzero = new TH1F("phidistrzero","phidistrzero", 100, (float)(7 - 1) * (2. * PI) / 9., (float)(7 + 2) * (2. * PI) / 9.);
	TH1* thetadistrgood = new TH1F("thetadistrgood","thetadistrgood", 300, -15, 15);	
	TH1* phidistrgood1 = new TH1F("phidistrgood1", "phidistrgood", 100, (float)(7 - 1) * (2. * PI) / 9., (float)(7 + 2) * (2. * PI) / 9.);
	TH1* phidistrzero1 = new TH1F("phidistrzero1", "phidistrzero", 100, (float)(7 - 1) * (2. * PI) / 9., (float)(7 + 2) * (2. * PI) / 9.);
	TH1* thetadistrgood1 = new TH1F("thetadistrgood1", "thetadistrgood", 300, -15, 15);

	TH1* achsp = new TH1F("achsp","achsp",5000,-1,499);
	TH1* achspGood = new TH1F("achspGood","achspGood",5000,-1,499);
	TH1* achsp1 = new TH1F("achsp1", "achsp", 5000, -1, 499);
	TH1* achspGood1 = new TH1F("achspGood1", "achspGood", 5000, -1, 499);
	TH1* achsp2 = new TH1F("achsp2", "achsp", 5000, -1, 499);
	TH1* achspGood2 = new TH1F("achspGood2", "achspGood", 5000, -1, 499);

	TH1* tchrdistr = new TH1F("tchrdistr","tchrdistr", 3000, -1000, 2000);
	TH1* tchrdistr1 = new TH1F("tchrdistr1","tchrdistr", 3000, -1000, 2000);
	TH1* tchrdistr2 = new TH1F("tchrdistr2","tchrdistr", 3000, -1000, 2000);

	TH1* thetadistrzero = new TH1F("thetadistrzero","thetadistrzero", 300, -15, 15);
	TH1* thetadistrzero1 = new TH1F("thetadistrzero1", "thetadistrzero", 300, -15, 15);
	TH1* thetadistrzero2 = new TH1F("thetadistrzero2", "thetadistrzero", 300, -15, 15);

	TH1* achspZero = new TH1F("achspZero","achspGood",20000,-100,100);
	TH1* achspZero1 = new TH1F("achspZero1","achspGood",20000,-100,100);
	TH1* achspZero2 = new TH1F("achspZero2","achspGood",20000,-100,100);

	TH1* ampspectre[9][14];
	for (int j = 0; j < 9; j++) {
		for (int i = 0; i < 14; i++) {
			sprintf(name, "hist%d", 14 * j + i + 1);
			sprintf(title, "amplitude spectrum, sensor%d, run range %d;ach", j + 1, i + 1);
			ampspectre[j][i] = new TH1F(name, title, 5000, -1, 499);
		}
	}

	//info ab ped integrated phi
	TFile* MFileIP = new TFile((workingDir + "phintegr_pedest.root").c_str());
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			TH1* hped = (TH1F*)MFileIP->Get(Form("ped_phintegr_zobl%d,sensor%d", j + 1, i + 1));
			ped[j][i] = hped->GetMean(1);
		}
	}
	TFile* MFileP = new TFile((workingDir + "pedestz.root").c_str());
	for (int i = 0; i < 9; i++) {
		TH1* hped1 = (TH1F*)MFileP->Get(Form("ped_forz,sensor%d", i + 1));
		pedz[i] = hped1->GetMean(1);
	}
	vector<float> rench;
	vector<float> eench;
	float beamlast = 0.;
	double lgran[4];
	int Count = 0, count1 = 0;
	float ztr[2] = { 0., 0.};
	int maxampla = 0;
	double maxach = 0.;
	int runrange = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr[0] = 12. / tan(theta[0]) + z0[0];
		ztr[1] = 12. / tan(theta[1]) + z0[1];
		schr[nch] = 0;
		if (fabs(beam - beamlast) > 1.7) {  //when energy change
			rench.push_back((float)run);
			eench.push_back(beam);
		}
		beamlast = beam;
		if (!RunIsGood(run)/* || (eventtime != 0)*/) {
			continue;
		}
			
		scount = 0;     //2017: 4001, 25500, 29501    ; 2018: 6001, 29800, 35801       ; 2019: 38000-43001       ; 2020: 45040-46541
		scount1 = 0;

		for (int j = 0; j < 9; j++) {
			//wichcounter1[j] = scount1;
			if (j == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
				if (phi[1] < 1.)
					phi[1] = phi[1] + 2. * PI;
			}
			//else
			while (schr[1 + scount] > 0.5) { scount += 1; }

			lgran[0] = truegran1[j][0] + 3. * meanshir[j][0];// +0.03;
			lgran[1] = truegran1[j][1] - 3. * meanshir[j][1];// -0.045;
			lgran[2] = truegran1[j][1] + 3. * meanshir[j][1];// +0.045;
			lgran[3] = truegran1[j][3] - 3. * meanshir[j][2];// -0.03;


			maxach = 0.;
			//find maxach and its index in amp signal array for this counter
			for (int f = scount1; f <= scount; f++) {
				if (ach[f] > maxach) {
					maxach = ach[f];
					maxampla = f;
				}
			}
			//split into run ranges to check amplitude distributions
			for (int i = 0; i < 14; i++) {
				if ((run > MinimRun + i*(MaximRun-MinimRun)/14) && (run < MinimRun + (i+1)*(MaximRun - MinimRun) / 14))
					runrange = i;
			}
			for (size_t pind = 0; pind < 2; pind++) {
				int zoblInd = -1;
				for (int k = 1; k < 13; k++) {
					if ((ztr[pind] > zobl[k]) && (ztr[pind] < zobl[k + 1])) {
						zoblInd = k;
					}
				}
				if ((((phi[pind] > lgran[0]) && (phi[pind] < lgran[1])) || ((phi[pind] > lgran[2]) && (phi[pind] < lgran[3]))) && (TimeIsGood(maxampla)) && (TimeIsGood(scount1)) && zoblInd != -1 && ach[maxampla] < achCut) {
					ach1[maxampla] = ach[maxampla]*sin(theta[pind]) - ach[scount1];
					bool eff = (bool)schr[scount1 + 1] && ach1[maxampla]>0.2;
					he[0]->Fill(run, (int)eff);
					he[j + 1]->Fill(run, (int)eff);
					ampspectre[j][runrange]->Fill(ach1[scount]);
					if((schr[maxampla] > 0.5)) {
						ha[0]->Fill(run, ach1[maxampla]);
						ha[j + 1]->Fill(run, ach1[maxampla]);

						/*if (run > 39420 && run < 39580 && j == 7) {
								achsp->Fill(ach1[maxampla]);
								tchrdistr->Fill(tchr[maxampla]);//scount1 for efficiency
								if(schr[scount1+1]==0){
									thetadistrzero->Fill(ztr[pind]);
									achspZero->Fill(ach[scount1]);
								}
								if (ach1[maxampla] > 0.2)
									achspGood->Fill(ach1[maxampla]);
							}
						if (run > 39000 && run < 39350 && j == 7){
								achsp1->Fill(ach1[maxampla]);
								tchrdistr1->Fill(tchr[maxampla]);
								if(schr[scount1+1]==0){
									thetadistrzero1->Fill(ztr[pind]);
									achspZero1->Fill(ach[scount1]);
								}
								if (ach1[maxampla] > 0.2)
									achspGood1->Fill(ach1[maxampla]);
							}
						if (run > 39650 && run < 39850 && j == 7) {
								achsp2->Fill(ach1[maxampla]);
								tchrdistr2->Fill(tchr[maxampla]);
								if(schr[scount1+1]==0){
									thetadistrzero2->Fill(ztr[pind]);
									achspZero2->Fill(ach[scount1]);
								}
								if (ach1[maxampla] > 0.2)
									achspGood2->Fill(ach1[maxampla]);
							}*/
						ampCounter[j]->Fill(ach1[maxampla] * sin(theta[pind]));
					}
				}
			}



			scount += 1;
			scount1 = scount;
		}

		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("obrabotano %d M entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}


	TFile* MyFile = new TFile((workingDir + "achvsrun.root").c_str(), "RECREATE");
	for (int i = 1; i < 10; i++) {
		sprintf(title, "avrsens%d", i);
		ha[i]->Write(title);
		sprintf(title, "evrsens%d", i);
		he[i]->Write(title);
	}
	ha[0]->Write("amplitudevsrun");
	he[0]->Write("effvsrun");
	MyFile->Close();

	TFile* MyFile1 = new TFile((workingDir + "achspectrum.root").c_str(), "RECREATE");
	for (int j = 0; j < 9; j++) {
		for (int i = 0; i < 14; i++) {
			sprintf(title, "aspecs%dr%d", j + 1, i + 1);
			ampspectre[j][i]->Write(title);
		}
	}
	MyFile1->Close();

	/*TFile* MyFileT = new TFile("/work/users/kladov/snd2k/R007-002/tempSpectra.root", "RECREATE");
	achsp->Write("achsp");
	achspGood->Write("achspGood");
	phidistrgood->Write("phidistrgood");
	phidistrzero->Write("phidistrzero");
	thetadistrgood->Write("thetadistrgood");
	thetadistrzero->Write("thetadistrzero");
	achsp1->Write("achsp1");
	achspGood1->Write("achspGood1");
	phidistrgood1->Write("phidistrgood1");
	phidistrzero1->Write("phidistrzero1");
	thetadistrgood1->Write("thetadistrgood1");
	thetadistrzero1->Write("thetadistrzero1");
	MyFileT->Close();*/
	/*TFile* MyFileT = new TFile("/work/users/kladov/snd2k/R007-002/tempSpectra1.root", "RECREATE");
	achsp->Write("achsp");
	achspGood->Write("achspGood");
	achsp1->Write("achsp1");
	achspGood1->Write("achspGood1");
	achsp2->Write("achsp2");
	achspGood2->Write("achspGood2");*/
	TFile* MyFileT = new TFile("/work/users/kladov/snd2k/R007-002/tempSpectra2.root", "RECREATE");
	/*tchrdistr->Write("tchrdistr");
	tchrdistr1->Write("tchrdistr1");
	tchrdistr2->Write("tchrdistr2");
	thetadistrzero->Write("thetadistrzero");
	thetadistrzero1->Write("thetadistrzero1");
	thetadistrzero2->Write("thetadistrzero2");
	achspZero->Write("achspZero");
	achspZero1->Write("achspZero1");
	achspZero2->Write("achspZero2");*/
	MyFileT->Close();

	TFile* MyFile2 = new TFile((workingDir + "achSpectCount.root").c_str(), "RECREATE");
	for (int j = 0; j < 9; j++) {
		sprintf(title, "aspecs%d", j + 1);
		ampCounter[j]->Write(title);
	}
	MyFile2->Close();
	/*TF1* fn = new TF1("scphi", "[0]+[1]*x", 27225, 29500);
	gr->Fit("scphi");
	cout << "0:" << fn->GetParameter(0) << endl;
	cout << "1:" << fn->GetParameter(1) << endl;
	double a = fn->GetParameter(1);
	double b = fn->GetParameter(0);
	cout << a*29500+b << endl;
	cout << "done!" << endl;*/

	ofstream fout;
	fout.open((workingDir + "ench.cpp").c_str());
	{
		for (size_t i = 0; i < rench.size(); i++) {
			fout << rench[i] << "	";
			fout << eench[i] << endl;
		}
	}
	fout.close();
}

void achspectresravn(){
	TFile* Myf = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2019/achspectrum.root");
	TH1* he0 = (TH1F*)Myf->Get(Form("aspecs%dr%d", 5, 0 + 1));
	he0->SetLineColor(1);
	he0->DrawNormalized("", 1);
	cout << he0->GetMean(1) << endl;
	for (int i = 1; i < 6; i++) {
		TH1* he = (TH1F*)Myf->Get(Form("aspecs%dr%d", 5, i + 1));
		he->SetLineColor(i+1);
		he->DrawNormalized("same", 1);
		cout << he->GetMean(1) << endl;
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
	}
}

vector<int> CalibrExp::getCalibRuns(string calName, bool read) {
	readMinMaxRun();
	TString cmd;
	cmd = "clbixlist " + calName + " CURRENT > " + calName + ".list";
	if (read) gSystem->Exec(cmd.Data());
	//cmd = "cat /work/users/kladov/snd2k/R007-002/" + calName + ".list | sed '/^-/ d' | sed 's#^.*\\([0-9]\\{4\\}-[0-9]\\{2\\}-[0-9]\\{2\\} [0-9]\\{2\\}:[0-9]\\{2\\}:[0-9]\\{2\\}\\).*$#\\1#'";
	cmd = "cat /work/users/kladov/snd2k/R007-002/" + calName + ".list | sed '/^-/ d' ";
	TString s = gSystem->GetFromPipe(cmd.Data());

	vector<TString> ans;
	TString elem;
	Ssiz_t from0 = 0;
	while (s.Tokenize(elem, from0, "\n")) {
		ans.push_back(elem);
	}
	TString str = "a";
	TString str1 = "b";
	vector<int> output;
	unsigned int vector_size = ans.size();
	if (vector_size != 0) {
		for (int i = 1; i < vector_size - 1; i++) {
			//cout << ans[i] << endl;
			str = ans[i].Copy();
			str1 = ans[i].Copy();
			TString tok;
			TString tok1;
			Ssiz_t from = 0;
			Ssiz_t from1 = 0;
			vector<TString> tokens;
			vector<TString> tokens1;
			//while (str.Tokenize(tok, from, "[_-]")) {
			while (str.Tokenize(tok, from, "[ ]")) {
				//cout << tok << endl;
				tokens.push_back(tok);
			}
			while (str1.Tokenize(tok1, from1, "[	]")) {
				//cout << tok1 << endl;
				tokens1.push_back(tok1);
			}
			const char* a0 = (const char*)tokens[0];
			int a = atoi(a0);
			const char* b0 = (const char*)tokens[1];
			int b = atoi(b0);
			//cout << a << endl;
			if (a < MaximRun && a > MinimRun)
				output.push_back(a);
		}
	}
	else
		output.clear();
	return output;
}

void CalibrExp::drawamptime() {
	readMinMaxRun();
	double re, ee;
	vector<float> rench;
	vector<float> eench;
	ifstream fin;
	fin.open((workingDir + "ench.cpp").c_str());
	{
		int j = 0;
		while (fin >> re >> ee){
			rench.push_back(re);
			cout << re << "	";
			eench.push_back(ee/100);
			cout << ee / 100 << endl;
			j++;
		}
	}
	fin.close();

	TFile* MyFile = new TFile((workingDir + "achvsrun.root").c_str());
	TFile* MyFile2 = new TFile((workingDir + "achSpectCount.root").c_str());

	//new bins finding with usage of average (by counters) ampl-time graph + filling file runAmp with amplitudes vs run
	TProfile* ha0 = (TProfile*)MyFile->Get("amplitudevsrun");
	TProfile* he0 = (TProfile*)MyFile->Get("effvsrun");
	he0->GetYaxis()->SetRangeUser(0.9, 1.);
	he0->Draw();
	int Nbins = ha0->GetNbinsX();
	int Nbinse = he0->GetNbinsX();
	ha0->Draw();
	bins.clear();
	ofstream foutra0;
	foutra0.open((workingDir + "runAmp0.dat").c_str());
	int binsI = 0;
	int bigE = 0;
	for (int i = 1; i < Nbins + 1; i++) {
		foutra0 << ha0->GetBinCenter(i) << "	" << ha0->GetBinWidth(i) << "	" << ha0->GetBinContent(i) << "	" << ha0->GetBinError(i) << endl;
		if (ha0->GetBinError(i) > 0.05)
			bigE += 1;
		else {
			if (bigE % 2 == 0)
				for (int x = 0; x < bigE + 1; x = x + 2)
					bins.push_back(ha0->GetBinLowEdge(i - bigE + x));
			else
				for (int x = 0; x < bigE + 1; x = x + 2)
					bins.push_back(ha0->GetBinCenter(i - bigE + x));
			bigE = 0;
		}
	}
	foutra0.close();

	//calculate consistency with drops of amplitude
	//naimvg, 
	/*vector<string> calibnames;
	ifstream in1;
	in1.open("caliblist.dat");
	string calname;
	while (in1 >> calname)
		calibnames.push_back(calname);
	in1.close();
	vector< vector<int> > calibsBorders;
	vector<int> allcalibruns;
	vector< vector<double> > allcalibrunslist;
	vector<int> dropCoincidense;
	vector<int> dropCoincidenseE;
	//calibsBorders.push_back(getCalibRuns("naimvg",false));
	//calibsBorders.push_back(getCalibRuns("dcpagen",false));
	for (size_t k = 0; k < calibnames.size(); k++) {
		TString temptstr = "ls " + calibnames[k] + ".list";
		TString s = gSystem->GetFromPipe(temptstr.Data());
		cout << s << endl;
		cout << calibnames[k] + ".list" << endl;
		if (s == calibnames[k] + ".list") {
			calibsBorders.push_back(getCalibRuns(calibnames[k], false));
			cout << "a" << endl;
		}
		else
			calibsBorders.push_back(getCalibRuns(calibnames[k], true));
		allcalibruns.push_back(calibsBorders[k].size());
		dropCoincidense.push_back(0);
		dropCoincidenseE.push_back(0);
		vector<double> tempvd;
		tempvd.push_back(2);
		allcalibrunslist.push_back(tempvd);
		allcalibrunslist[k].clear();
		for (size_t j = 0; j < calibsBorders[k].size(); j++) {
			for (int i = 1; i < Nbins + 1; i++)
				if (ha0->GetBinCenter(i - 1) < calibsBorders[k][j] && ha0->GetBinCenter(i) >= calibsBorders[k][j])
					if ((1.5 * (ha0->GetBinError(i - 1) + ha0->GetBinError(i)) < fabs(ha0->GetBinContent(i) - ha0->GetBinContent(i - 1))) && (ha0->GetBinContent(i - 1) > 2 && ha0->GetBinContent(i) > 2))
						dropCoincidense[k] += 1;
			for (int i = 1; i < Nbinse + 1; i++)
				if (he0->GetBinCenter(i - 1) < calibsBorders[k][j] && he0->GetBinCenter(i) >= calibsBorders[k][j])
					//if ((2 * (he0->GetBinError(i - 1) + he0->GetBinError(i)) < fabs(he0->GetBinContent(i) - he0->GetBinContent(i - 1))) && (he0->GetBinContent(i - 1) > 0.5 && he0->GetBinContent(i) > 0.5)) {
					if ((0.003 < fabs(he0->GetBinContent(i) - he0->GetBinContent(i - 1))) && (he0->GetBinContent(i - 1) > 0.9 && he0->GetBinContent(i) > 0.9)) {
						dropCoincidenseE[k] += 1;
						allcalibrunslist[k].push_back(calibsBorders[k][j]);
					}
		}
	}
	int alldrops = 0;
	int alldropsE = 0;
	for (int i = 1; i < Nbins + 1; i++)
		if ((1.5*(ha0->GetBinError(i - 1) + ha0->GetBinError(i)) < fabs(ha0->GetBinContent(i) - ha0->GetBinContent(i - 1))) && (ha0->GetBinContent(i - 1) > 2 && ha0->GetBinContent(i) > 2))
			alldrops += 1;
	for (int i = 1; i < Nbinse + 1; i++)
		//if ((2 * (he0->GetBinError(i - 1) + he0->GetBinError(i)) < fabs(he0->GetBinContent(i) - he0->GetBinContent(i - 1))) && (he0->GetBinContent(i - 1) > 0.5 && he0->GetBinContent(i) > 0.5))
		if ((0.003 < fabs(he0->GetBinContent(i) - he0->GetBinContent(i - 1))) && (he0->GetBinContent(i - 1) > 0.9 && he0->GetBinContent(i) > 0.9))
			alldropsE += 1;
	cout << "all drops	" << alldrops << endl;
	for (size_t k = 0; k < calibnames.size(); k++)
		cout << calibnames[k] << "		" << allcalibruns[k] << "	" << dropCoincidense[k] << endl;
	cout << endl << endl;
	cout << "all dropsE	" << alldropsE << endl;
	for (size_t k = 0; k < calibnames.size(); k++)
		cout << calibnames[k] << "		" << allcalibruns[k] << "	" << dropCoincidenseE[k] << endl;

	he0->Draw();
	he0->GetYaxis()->SetRangeUser(0.9,1.0);
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	vector<TLine*> line1;
	for (size_t k = 0; k < calibnames.size(); k++) {
		for (size_t i = 0; i < allcalibrunslist[k].size(); i++) {
			line1.push_back(new TLine(allcalibrunslist[k][i], 0, allcalibrunslist[k][i], 1.1 * he0->GetMaximum()));
			line1.back()->SetLineColor((2+k)%5);
			line1.back()->Draw();
			c->Update();
		}
	}*/


	//get profiles (1-9) and hists from files
	TProfile* ha[9];
	TProfile* he[9];
	double meanAmp[9];
	for (int k = 0; k < 9; k++) {
		meanAmp[k] = ((TProfile*)MyFile2->Get(Form("aspecs%d", k + 1)))->GetMean(1);
		ha[k] = (TProfile*)MyFile->Get(Form("avrsens%d", k + 1));
		he[k] = (TProfile*)MyFile->Get(Form("evrsens%d", k + 1));
	}

	//writing coefficients to the file
	ofstream foutra;
	foutra.open((workingDir + "runAmp.dat").c_str());
	foutra << "run	" << "bin width";
	for (int k = 0; k < 9; k++) {
		foutra << Form("	count%d", k);
	}
	foutra << endl;
	for (size_t j = 1; j < Nbins + 1; j++) {
		//cout << ha[0]->GetBinCenter(j) << "	" << ha[0]->GetBinWidth(j);
		foutra << ha[0]->GetBinCenter(j) << "	" << ha[0]->GetBinWidth(j);
		for (int k = 0; k < 9; k++) {
			//cout << "	" << ha[k]->GetBinContent(j) / meanAmp[k];
			foutra << "	" << ha[k]->GetBinContent(j) / meanAmp[k];
		}
		foutra << endl;
	}
	foutra.close();



	//drawing
	for (int k = 0; k < 9; k++) {
		vector<double> date, ampl, ddate, dampl, effic, deffic;
		binsI = 0;
		bigE = 0;
		for (size_t i = 1; i < Nbins + 1; i++) {
			ampl.push_back(ha[k]->GetBinContent(i));
			effic.push_back(100*he[k]->GetBinContent(i));
			date.push_back(ha[k]->GetBinCenter(i));
			dampl.push_back(ha[k]->GetBinError(i));
			deffic.push_back(100*he[k]->GetBinError(i));
			ddate.push_back(ha[k]->GetBinWidth(i)/2.);
		}
		int numbb = ampl.size();
		//ha[k]->Draw();
		TGraphErrors* gr = new TGraphErrors(numbb, &date[0], &ampl[0], &ddate[0], &dampl[0]);
		TGraphErrors* gr1 = new TGraphErrors(numbb, &date[0], &effic[0], &ddate[0], &deffic[0]);
		TGraph* gr2 = new TGraph(rench.size(), &rench[0], &eench[0]);
		gr->SetTitle("amplitude vs run; run; mean amplitude-eff, p.e.");
		//gr->GetXaxis()->SetTimeDisplay(1);
		//gr1->GetXaxis()->SetTimeDisplay(1);
		//gr->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d-%m}");
		gr->SetMarkerSize(0.5);
		gr->SetMarkerStyle(20);
		gr2->SetMarkerSize(0.5);
		gr1->SetMarkerSize(0.5);
		gr2->SetMarkerStyle(21);
		gr1->SetMarkerStyle(21);
		gr1->GetYaxis()->SetRangeUser(85,100);
		gr->GetYaxis()->SetRangeUser(0.0,7);
		gr->SetLineColor(3);
		gr2->SetLineColor(4);
		gr1->SetLineColor(2);
		gr->Draw("AP");
		//gr2->Draw("same");
		//gr1->Draw("same");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->SetFillColor(0);
		c->SetFrameFillColor(0);
		c->SetGridy();
		c->Update();
		cout << "done!" << endl;
		cin.get();
	}
	MyFile->Close();

}

void CalibrExp::timespectra() {
	readGoodRuns();
	readMeanshir();
	readTrueGran1();
	TChain chain("t1");
	std::string chadds = dirSelectedExp + prefixSelectedExp + "**exp00**.root";
	chain.Add(chadds.c_str());
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


	TH1* histo = new TH1F("histo", "time spectre;time,TDC chanels;", 500, -150, 350);
	TH1* histi = new TH1F("histi", "time spectre;time,TDC chanels;", 500, -150, 350);
	TH1* histp = new TH1F("histp", "time spectre;time,TDC chanels;", 1000, -150, 850);

	int Count = 0, count1 = 0;
	float ztr[2] = { 0.,0. };
	int maxampla = 0;
	int predmaxampla = 0;
	double maxach = 0.;
	double predmaxach = 0.;
	double lgran[4];
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr[0] = 12.0 / tan(theta[0]) + z0[0];
		ztr[1] = 12.0 / tan(theta[1]) + z0[1];
		schr[nch] = 0.;
		if (!RunIsGood(run)) {
			continue;
		}
		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			if (j == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
				if (phi[1] < 1.)
					phi[1] = phi[1] + 2. * PI;
			}
			while (schr[1 + scount] > 0.5) { scount += 1; }

			lgran[0] = truegran1[j][0] + 3. * meanshir[j][0];
			lgran[1] = truegran1[j][1] - 3. * meanshir[j][1];// -0.045;
			lgran[2] = truegran1[j][1] + 3. * meanshir[j][1];// +0.045;
			lgran[3] = truegran1[j][3] - 3. * meanshir[j][2];// -0.03;

			predmaxach = -1.;
			maxach = -100.;
			predmaxampla = scount1;
			for (int f = scount1; f <= scount; f++) {
				if (ach[f] > maxach) {
					maxach = ach[f];
					maxampla = f;
				}
				if (ach[f] < maxach-0.0001 && ach[f] > predmaxach) {
					predmaxach = ach[f];
					predmaxampla = f;
				}
			}
			for (size_t partI = 0; partI < 2; partI++) {
				//scount 1 is always 0, so if we have some signal other than 0 
				if (scount - scount1 > 0.5 && (ztr[partI] > zobl[1]) && (ztr[partI] < zobl[13])) {
					for (int a = maxampla; a >= maxampla; a--) {
						if ((((phi[partI] > lgran[0]) && (phi[partI] < lgran[1])) || ((phi[partI] > lgran[2]) && (phi[partI] < lgran[3]))) && (schr[a] > 0.5))
							histo->Fill(tchr[a]);
							//histo->Fill(tchr[scount1]);
					}
				}
				// if it for example 0 1 2 0 - second by amplitude signal exists
				if (scount - scount1 > 1.5 && (ztr[partI] > zobl[1]) && (ztr[partI] < zobl[13])) {
					for (int a = predmaxampla; a >= predmaxampla; a--) {
						if ((((phi[partI] > lgran[0]) && (phi[partI] < lgran[1])) || ((phi[partI] > lgran[2]) && (phi[partI] < lgran[3]))) && (schr[a] > 0.5))
							histi->Fill(tchr[a]);
					}
				}
				if ((ztr[partI] > zobl[1]) && (ztr[partI] < zobl[13])) {
					if ((((phi[partI] > lgran[0]) && (phi[partI] < lgran[1])) || ((phi[partI] > lgran[2]) && (phi[partI] < lgran[3]))) && (schr[scount1] < 0.5))
						histp->Fill(tchr[scount1]);
				}
			}
			scount += 1;
			scount1 = scount;
		}
		

		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("processed %d*M entries", count1) << endl;
			Count = 0;
		}
	}
	histo->SetLineColor(2);
	histi->SetLineColor(13);
	histo->Draw();
	//histi->Draw("same");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	TLine* line11 = new TLine(110, 0, 110, 1.1 * histo->GetMaximum());
	line11->SetLineColor(kGreen);
	line11->Draw();
	c->Update();

	TFile* MyFile = new TFile((workingDir + "timespectre.root").c_str(), "RECREATE");
	histo->Write("TspectreMain");
	histi->Write("TspectreSec");
	histp->Write("TspectrePed");
	//MyFile->Write();
	MyFile->Close();
}

void CalibrExp::drawtsp() {

	// get and draw
	TFile* MyFile = new TFile((workingDir + "timespectre.root").c_str());
	TH1* hma = (TH1F*)MyFile->Get("TspectreMain");
	TH1* hse = (TH1F*)MyFile->Get("TspectreSec");
	hma->SetLineColor(2);
	hse->SetLineColor(13);
	hma->Draw();
	hse->Draw("same");

	// fit around maximum bin (width is about 10-20 ADC) with Gauss
	double center = hma->GetBinCenter(hma->GetMaximumBin());
	TF1* f1 = new TF1("f1", "[0]*exp(- (x-[1])*(x-[1]) / ([2]*[2]))", center-5, center+5);
	f1->SetParameters(hma->GetMaximum(),100,3);
	hma->Fit("f1", "", "", center-5,center+5);

	// find maximum values on the sides of main pike (mean - width - 0::15) 
	double maxRight = 0;
	size_t maxRightInd = 0;
	double maxLeft = 0;
	size_t maxLeftInd = 0;
	for (size_t i = 0; i < 15; i++) {
		if (hma->GetBinContent(hma->GetMaximumBin() - (int)f1->GetParameter(2) - i) > maxLeft) {
			maxLeft = hma->GetBinContent(hma->GetMaximumBin() - 2.5*(int)f1->GetParameter(2) - i);
			maxLeftInd = i;
		}
		if (hma->GetBinContent(hma->GetMaximumBin() + (int)f1->GetParameter(2) + i) > maxRight) {
			maxRight = hma->GetBinContent(hma->GetMaximumBin() + 2.5*(int)f1->GetParameter(2) + i);
			maxRightInd = i;
		}
	}
	// if right maximum is higher than left one - then we should take right side and cut left (it is considered to occur before it can kinematically), so tcut is mean-2*widthw
	double tCut = f1->GetParameter(1);
	int lowerOrHigher = 0;
	if (maxRight > maxLeft) {
		lowerOrHigher = 1;
		tCut -= 2 * (f1->GetParameter(2));
	}
	else {
		lowerOrHigher = -1;
		tCut += 2 * (f1->GetParameter(2));
	}
	cout << lowerOrHigher << "	" << tCut << endl;

	TLine* l1 = new TLine(tCut, 0, tCut, f1->GetParameter(0));
	l1->Draw("same");

	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	cin.get();

	TH1* hsp = (TH1F*)MyFile->Get("TspectrePed");
	hsp->Draw();
	c->Update();

}

void CalibrExp::enpoint() {
	readGoodRuns();
	TChain chain("t1");

	std::string chadds = dirSelectedExp + "/true1__**exp00**.root";
	chain.Add(chadds.c_str());

	//chain.Add("/work/users/kladov/snd2k/R006-003/selected2018/true1__**exp00**.root");
	const Long64_t entries = chain.GetEntries();
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("beam", &beam);
	TH1* histenpo = new TH1F("enpo", "enpo", 800, 300, 1100);

	EntriesInRun mass1[MaximRun - MinimRun];
	for (int i = 0; i < MaximRun - MinimRun; i++) {
		mass1[i].run = i+ MinimRun;
		mass1[i].energy = 0;
		mass1[i].entries = 0;
	}
	int entrall = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//histo->Fill(run);
		if ((run > 26814)/* && RunIsGood(run)*/) {
			histenpo->Fill(beam);
			mass1[run - MinimRun].run = run;
			if (mass1[run - MinimRun].energy > 1 && ((mass1[run - MinimRun].energy < beam-0.5) || (mass1[run - MinimRun].energy > beam+0.5)))
				cout << "energy is not the same in run " << run << endl;
			mass1[run - MinimRun].energy = beam;
			mass1[run - MinimRun].entries += 1;
			entrall += 1;
		}
	}
	vector<double> enpoint[2];
	for (int i = 0; i < 800; i++) {
		if (histenpo->GetBinContent(i) > 0) {
			enpoint[0].push_back(i + 299);
			enpoint[1].push_back(histenpo->GetBinContent(i));
		}
	}
	for (size_t i = 0; i < enpoint[0].size(); i++)
		cout << enpoint[0][i] << "	" << (((int)(enpoint[1][i] * (6000000.0 / entrall)))/10000)*10000 << endl;
	ofstream fout;
	fout.open("/work/users/kladov/snd2k/R006-003/maindir/2018/entriesinrune.txt");
	int SumEntrForMod = 0;
	for (int i = 0; i < MaximRun - MinimRun; i++) {
		if ((mass1[i].entries * (2000000.0 / entrall)) >= 1) {
			fout << mass1[i].run << "	";
			fout << mass1[i].energy << "	";
			fout << (int)(mass1[i].entries * (2000000.0 / entrall)) << endl;
			SumEntrForMod += (int)(mass1[i].entries * (2000000.0 / entrall));
		}
	}
	cout << SumEntrForMod << endl;
	fout.close();
}

void CalibrExp::nvsphiorig() {
	TChain chain("t1");
	//chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected2018/true1__**exp00**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("schr", &schr);
	chain.SetBranchAddress("tchr", &tchr);
	chain.SetBranchAddress("event", &eventn);
	chain.SetBranchAddress("region", &region);
	chain.SetBranchAddress("eventtime", &eventtime);
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	char name[20];
	char title[100];
	TH1* hii[9][4];

	for (int l = 0; l < 4; l++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "histog%d", 9 * l + i + 1);
			sprintf(title, "#phi,sensor%d,gran%d;#phi", i + 1,l+1);
			hii[i][l] = new TH1F(name, title, 1200, truegran1[i][l] - 1.3, truegran1[i][l] + 1.3);
		}
	}
	TH1* hii1 = new TH1F("profile1", "zasel;#phi", 6300, 0. - 0.01, 2. * PI + 0.01);
	ofstream fout3;
	fout3.open((workingDir + "bad_entries_87_104.txt").c_str());
	fout3 << "run" << "	" << "event" << endl;
	float ztr = 0.;
	float ztr1 = 0.;
	double lgran[4];
	int maxampla = 0;
	double maxach = 0.;
	int scount = 0;
	int scount1 = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.0 / tan(theta[0]) + z0[0];
		ztr1 = 12.0 / tan(theta[1]) + z0[1];
		schr[nch] = 0.;
		scount = 0;
		scount1 = 0;
		for (int i = 0; i < 9; i++) {
			if (i == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
				if (phi[1] < 1.)
					phi[1] = phi[1] + 2. * PI;
				while (1 + scount < nch) { scount += 1; }
			}
			else
				while (schr[1 + scount] > 0.5) { scount += 1; }
			maxach = 0.;
			for (int f = scount1; f <= scount; f++) {
				if (ach[f] > maxach) {
					maxach = ach[f];
					maxampla = f;
				}
			}
			lgran[0] = truegran1[i][0]-0.05;// +0.03;
			lgran[1] = truegran1[i][1];// -0.045;
			lgran[2] = truegran1[i][1];// +0.045;
			lgran[3] = truegran1[i][3]+0.05;// -0.03;
			if (((phi[0] > lgran[0]) && (phi[0] < lgran[3]))){
				if ((ztr > zobl[0] && ztr < zobl[14]) && ((schr[maxampla] < 0.5) /*|| ((schr[maxampla] > 0.5)*/ && /*(tchr[maxampla] > 87) && *//*(tchr[maxampla] < 110) &&*/ (eventtime == 0) && (run>33000) && (run < 34000))) {
					for (int l = 0; l < 4; l++) {
						hii[i][l]->Fill(phi[0]);
					}
					hii1->Fill(phi[0]);
					//fout3 << run << "	" << eventn << endl;
				}
			}
			if (((phi[1] > lgran[0]) && (phi[1] < lgran[3]))) {
				if ((ztr1 > zobl[0] && ztr1 < zobl[14]) && ((schr[maxampla] < 0.5)/* || ((schr[maxampla] > 0.5)*/ &&/*(tchr[maxampla] > 87) && *//*(tchr[maxampla] < 110) && */(eventtime == 0) && (run > 33000) && (run < 34000))) {
					for (int l = 0; l < 4; l++) {
						hii[i][l]->Fill(phi[1]);
					}
					hii1->Fill(phi[1]);
					//fout3 << run << "	" << eventn << endl;
				}
			}
			scount += 1;
			scount1 = scount;
		}
	}
	fout3.close();
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2018/nvsphio.root", "RECREATE");
	for (int l = 0; l < 4; l++) {
		for (int a = 0; a < 9; a++) {
			sprintf(title, "sensor%d,gran%d", a + 1, l+1);
			hii[a][l]->Write(title);
		}
	}
	hii1->Write("zasel");
	MyFile->Close();
	cout << "done!" << endl;
}

void drawnvsphi() {
	//TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2018/fittedprofiles.root");
	//TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/nvsphi.root");
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2018/nvsphio.root");


	for (int i = 0; i < 1; i++) {
		//TH1* h = (TH1*)MyFile->Get(Form("zobl%d,sensor%d", i + 1, 6));
		TH1* h = (TH1*)MyFile->Get(Form("zasel"));
		//TH1* h = (TH1*)MyFile->Get(Form("sensor%d", 6));
		//h->GetYaxis()->SetRangeUser(0., 5.);
		h->SetLineColor(i + 1);
		h->Draw();
		//h->DrawNormalized("", 2500);
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
	}

	/*TDirectory* d1 = (TDirectory*)f->Get(Form("zobl%d,sensor8", 2));
	TF1* tf1 = (TF1*)d1->Get(Form("zobl%d,sensor%d,full", 2, 8));
	//TProfile *tf1 = (TProfile*)f->Get(Form("zobl%d,sensor%d", 2, 8));
	tf1->SetLineColor(2);
	tf1->Draw("same");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();*/

}

void CalibrExp::phiminusphis() {
	TChain chain("t1");
	//chain.Add("/work/users/kladov/snd2k/R006-003/selected1/true1_**exp00**.root");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected2020/true1__**exp00**.root");

	chain.SetBranchAddress("phi", &phi);
	chain.SetBranchAddress("phis", &phis);
	chain.SetBranchAddress("run", &run);
	chain.SetBranchAddress("tchr", &tchr);
	chain.SetBranchAddress("ach", &ach);
	chain.SetBranchAddress("schr", &schr);
	chain.SetBranchAddress("theta", &theta);
	chain.SetBranchAddress("z0", &z0);
	chain.SetBranchAddress("nch", &nch);
	chain.SetBranchAddress("eventtime", &eventtime);

	chain.SetBranchAddress("d2phi", &d2phi);
	chain.SetBranchAddress("dphirho", &dphirho);
	chain.SetBranchAddress("d2rho", &d2rho);
	chain.SetBranchAddress("d2z0", &d2z0);
	chain.SetBranchAddress("d2cosTh", &d2cosTh);
	chain.SetBranchAddress("dz0cosTh", &dz0cosTh);

	chain.SetBranchAddress("Dtheta", &Dtheta);
	chain.SetBranchAddress("Dphi", &Dphi);
	chain.SetBranchAddress("energyerr", &energyerr);

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	char name[20];
	char title[100];

	sprintf(name, "profile");
	sprintf(title, "#phi-phis;#phi");
	TProfile* hii = new TProfile(name, title, 6300, 0.-0.01, 2.*PI+0.01);
	TProfile* hii1 = new TProfile("profile1", "#phi-phis1;#phi", 6300, 0.-0.01, 2.*PI+0.01);

	//info ab phi-phis
	/*TFile* f11 = new TFile((workingDir + "phiminphis.root").c_str());
	TProfile* pphi = (TProfile*)f11->Get("phi-phis_vs_phis");*/

	/*TFile* f12 = new TFile((workingDir + "phiminphisphi0.root").c_str());
	TProfile* pphi0 = (TProfile*)f12->Get("phi-phis_vs_phi");*/
	int scount = 0;
	int scount1 = 0;
	int smplc = 1;
	float ztr = 0.;
	float ztr1 = 0.;
	int maxampla = 0;
	double maxach = 0.;
	double lgran[4];
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr = 12.0 / tan(theta[0]) + z0[0];
		ztr1 = 12.0 / tan(theta[1]) + z0[1];
		scount = 0;
		scount1 = 0;
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
			maxach = 0.;
			for (int f = scount1; f <= scount; f++) {
				if (ach[f] > maxach) {
					maxach = ach[f];
					maxampla = f;
				}
			}
			lgran[0] = truegran1[j][0] - 0.05;// +0.03;
			lgran[1] = truegran1[j][1];// -0.045;
			lgran[2] = truegran1[j][1];// +0.045;
			lgran[3] = truegran1[j][3] + 0.05;// -0.03;
			if (((phi[0] > lgran[0]) && (phi[0] < lgran[3]))) {
				//if (fabs((phi[0] - phis[0]) - (pphi->GetBinContent(pphi->FindBin(phis[0])))) < 0.01)
				//if ((d2phi[0] < 0.00003) && (dphirho[0] < 0.00012) && (d2rho[0] < 0.0007) && (d2z0[0] < 0.1) && (d2cosTh[0] < 0.015) && (dz0cosTh[0] > -0.01) && (energyerr[0] > 31.3 && energyerr[0] < 32.9) && (Dphi[0] < 0.02)) {
				if (ztr > zobl[0] && ztr < zobl[14]) {
					if (/*(schr[maxampla] < 0.5) || ((schr[maxampla] > 0.5) && */(/*(tchr[maxampla] > 87) && (tchr[maxampla] < 104) && */(eventtime == 0))) {
						hii->Fill(phi[0], phi[0] - phis[0]);
						smplc += 1;
					}
				}
				//}
			}
			hii1->Fill(phi[0], phi[0] - phis[0]);
			if (((phi[1] > lgran[0]) && (phi[1] < lgran[3]))) {
				//if (fabs((phi[1] - phis[1]) - (pphi->GetBinContent(pphi->FindBin(phis[1])))) < 0.01)
				//if ((d2phi[1] < 0.00003) && (dphirho[1] < 0.00012) && (d2rho[1] < 0.0007) && (d2z0[1] < 0.1) && (d2cosTh[1] < 0.015) && (dz0cosTh[1] > -0.01) && (energyerr[1] > 31.3 && energyerr[1] < 32.9) && (Dphi[1] < 0.02)){
				if ((ztr1 > zobl[0] && ztr1 < zobl[14]) && (/*(schr[maxampla] < 0.5) || ((schr[maxampla] > 0.5) && */(/*(tchr[maxampla] > 87) && (tchr[maxampla] < 104) && */(eventtime == 0)))) {
					hii->Fill(phi[1], phi[1] - phis[1]);
				}
				//}
			}
			hii1->Fill(phi[1], phi[1] - phis[1]);
			//if (smplc / 100000 == 0)
				//cout << smplc << endl;
			scount += 1;
			scount1 = scount;
		}
	}

	/*double binnu[6300];
	double phimphis[6300];
	for (int i = 0; i < 6300; i++) {
		binnu[i] = pphi0->GetBinCenter(i+1);
		phimphis[i] = pphi0->GetBinContent(i+1);
	}
	for (int j = 0; j < 100; j++) {
		for (int i = 1; i < 6299; i++) {
			if (fabs(phimphis[i - 1] - phimphis[i + 1]) > 0.002)
				phimphis[i] = (phimphis[i + 1] + phimphis[i - 1]) / 2.;
		}
	}
	phimphis[0] = 0.;
	phimphis[6299] = 0.;
	TGraph* gr1 = new TGraph(6300, binnu, phimphis);
	gr1->SetLineColor(2);*/

	hii->Draw();
	hii1->SetLineColor(2);
	hii1->Draw("same");
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2020/phiminphisphi0.root", "RECREATE");
	sprintf(title, "phi-phis_vs_phi");
	hii->Write(title);
	//pphi0->Draw();
	//gr1->Draw("same");
	MyFile->Close();
	cout << "done!" << endl;
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
	/*TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R006-003/selected2019/true1_**exp00**.root");

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


	//info ab ped 
	TFile* f0 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2019/phintegr_pedest.root");
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			TH1* hped = (TH1F*)f0->Get(Form("ped_phintegr_zobl%d,sensor%d", j + 1, i + 1));
			ped[j][i] = hped->GetMean(1);
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
		if ((run < 40100) && (goodrnumb[run - 38001] == 1) && (goodreff[run - 38001] == 1)) {
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
								ach1[scount] = ach[scount] - ped[k][j];
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
								ach1[scount] = ach[scount] - ped[k][j];
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
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2019/random.root", "RECREATE");
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
	TFile* ff = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2019/random.root");
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

void CalibrExp::achspec() {
	readGoodRuns();
	TChain chain("t1");

	std::string chadds = dirSelectedExp + "/true1__**exp00**.root";
	chain.Add(chadds.c_str());

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
	char name[20];
	char title[100];
	TH1* h[10];	//ach spectra, 9 counters + 10 for average
	for (int i = 0; i < 10; i++) {
		sprintf(name, "hprofped%d", i + 1);
		sprintf(title, "ped,sensor%d;ach", i + 1);
		h[i] = new TH1F(name, title, 5000, -1, 499);
	}
	TH1* hg[10]; //amount of all events to map and ach spectra
	for (int i = 0; i < 10; i++) {
		sprintf(name, "hprofgped%d", i + 1);
		sprintf(title, "ped,sensor%d;ach", i + 1);
		hg[i] = new TH1F(name, title, 5000, -1, 499);
	}
	TH1* hig[10]; //time spectra, 9 counters + 10 for average
	for (int i = 0; i < 10; i++) {
		sprintf(name, "hprofgt%d", i + 1);
		sprintf(title, "tchr,sensor%d;tds chanals", i + 1);
		hig[i] = new TH1F(name, title, 500, -1, 499);
	}
	TH1* hi[10][30];	//time spectra divided into 30 run ranges, 9 counters + 10 for average
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 30; j++) {
			sprintf(name, "hproft%d", 30*(i)+j+1);
			sprintf(title, "tchr,sensor%d,runr%d;tchr", i + 1, j+1);
			hi[i][j] = new TH1F(name, title, 500, -1, 499);
		}
	}
	TH1* hiz[30];	//where zero counts are sitting vs ztr(12)
	for (int i = 0; i < 30; i++) {
		sprintf(name, "hprofz%d", i + 1);
		sprintf(title, "zasel0,runr%d;z", i + 1);
		hiz[i] = new TH1F(name, title, 300, -15, 15);
	}
	//TH1* hii1 = new TH1F("profile1", "zasel;#phi", 6300, 0. - 0.01, 2. * PI + 0.01);
	TH1* hii2 = new TH1F("profile2", "zasel;#phi", 6300, 0. - 0.01, 2. * PI + 0.01);

	float ztr[2];
	ztr[0] = 0.; ztr[1] = 0.;
	double lgran[4];
	int maxampla = 0;
	double maxach = 0.;
	int Count = 0;
	int count1 = 0;
	int runrange = 0;
	int countentr[30];
	for (int i = 0; i < 30; i++)
		countentr[i] = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr[0] = 12.0 / tan(theta[0]) + z0[0];
		ztr[1] = 12.0 / tan(theta[1]) + z0[1];
		schr[nch] = 0.;
		if ((run > 27250) && RunIsGood(run)) { 
			scount = 0;
			scount1 = 0;
			for (int j = 0; j < 9; j++) {
				if (j == 8) {
					if (phi[0] < 1.)
						phi[0] = phi[0] + 2. * PI;
					if (phi[1] < 1.)
						phi[1] = phi[1] + 2. * PI;
					//while (1 + scount < nch) { scount += 1; }
				}
				//else
					while (schr[1 + scount] > 0.5) { scount += 1; }


				lgran[0] = truegran1[j][0] + 3. * meanshir[j][0];// +0.03;
				lgran[1] = truegran1[j][1] - 3. * meanshir[j][1];// -0.045;
				lgran[2] = truegran1[j][1] + 3. * meanshir[j][1];// +0.045;
				lgran[3] = truegran1[j][3] - 3. * meanshir[j][2];// -0.03;

				maxach = -1000.;
				for (int f = scount1; f <= scount; f++) {
					if (ach[f] > maxach) {
						maxach = ach[f];
						maxampla = f;
					}
				}
				for (int i = 0; i < 30; i++) {
					if ((run > MinimRun + i * (MaximRun - MinimRun) / 30.) && (run < MinimRun + (i + 1) * (MaximRun - MinimRun) / 30.))
						runrange = i;
				}

				for (int prtc = 0; prtc < 2; prtc++) {
					if ((((phi[prtc] > lgran[0]) && (phi[prtc] < lgran[1])) || ((phi[prtc] > lgran[2]) && (phi[prtc] < lgran[3]))) && (ztr[prtc] > zobl[1] && ztr[prtc] < zobl[13])) {
						//bool rng = (rand() % 100) < 4;
						//if (rng == 1)
							//schr[maxampla] = 0;
						if (((schr[maxampla] > 0.5))) {
							hi[j][runrange]->Fill(tchr[maxampla]);
							hi[9][runrange]->Fill(tchr[maxampla]);
							hig[j]->Fill(tchr[maxampla]);
							hig[9]->Fill(tchr[maxampla]);
						}
						if ((schr[maxampla] < 0.5) || ((schr[maxampla] > 0.5) && (tchr[maxampla] < 110))) {
							hg[j]->Fill(ach[maxampla]);
							hg[9]->Fill(ach[maxampla]);
							countentr[runrange] += 1;
						}
						if (schr[maxampla] < 0.5) {
							hii2->Fill(phi[prtc]);
							hiz[runrange]->Fill(ztr[prtc]);
						}
					}
				}

				scount += 1;
				scount1 = scount;
			}
		}
		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("obrabotano %d*M entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}
	/*cout << (double)(hiz1->GetEntries()) << endl;
	cout << (double)(h[9]->GetEntries()) << endl;
	cout << (double)(hiz2->GetEntries()) << endl;
	cout << (double)(hg[9]->GetEntries()) << endl;

	cout << "countentr 2,3,4,5,6" << endl;
	cout << countentr2 << endl;
	cout << countentr3 << endl;
	cout << countentr4 << endl;
	cout << countentr5 << endl;
	cout << countentr6 << endl;

	cout << "norms 2,3,4,5,6" << endl;
	cout << ((double)(hiz2->GetEntries())) / countentr2 << endl;
	cout << ((double)(hiz3->GetEntries())) / countentr3 << endl;
	cout << ((double)(hiz4->GetEntries())) / countentr4 << endl;
	cout << ((double)(hiz5->GetEntries())) / countentr5 << endl;
	cout << ((double)(hiz6->GetEntries())) / countentr6 << endl;

	double norm1 = ((double)(hiz1->GetEntries()))/ ((double)(h[9]->GetEntries()));
	double norm2 = ((double)(hiz2->GetEntries()))/ ((double)(hg[9]->GetEntries()));
	hiz1->SetLineColor(2);
	hiz1->DrawNormalized("",norm1);
	hiz2->DrawNormalized("same",norm2);
	hiz2->DrawNormalized("", ((double)(hiz2->GetEntries())) / countentr2);
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	hiz3->DrawNormalized("", ((double)(hiz3->GetEntries())) / countentr3);
	c->Update();
	hiz4->DrawNormalized("", ((double)(hiz4->GetEntries())) / countentr4);
	c->Update();
	hiz5->DrawNormalized("", ((double)(hiz5->GetEntries())) / countentr5);
	c->Update();
	hiz6->DrawNormalized("", ((double)(hiz6->GetEntries())) / countentr6);
	c->Update();*/

	cout << countentr[0] << endl;
	hiz[0]->DrawNormalized("", ((double)(hiz[0]->GetEntries())) / countentr[0]);
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	for (int i = 1; i < 30; i++) {
		cout << countentr[i] << endl;
		hiz[i]->DrawNormalized("", ((double)(hiz[i]->GetEntries())) / countentr[i]);
		c->Update();
		cin.get();
	}

	cout << "done" << endl;
	TFile* MyFile = new TFile((workingDir + "tchranom.root").c_str(), "RECREATE");

	for (int j = 0; j < 10; j++) {
		//sprintf(title, "ach_anom,sensor%d", j + 1);
		//h[j]->Write(title);
		sprintf(title, "ach_norm,sensor%d", j + 1);
		hg[j]->Write(title);
		for (int i = 0; i < 30; i++) {
			sprintf(title, "t,sensor%d,runr%d", j + 1, i + 1);
			hi[j][i]->Write(title);
		}
		sprintf(title, "tch_norm,sensor%d", j + 1);
		hig[j]->Write(title);
	}
	sprintf(title, "zasp_norm,sensor%d", 9 + 1);
	hii2->Write(title);
	for (int i = 0; i < 30; i++) {
		sprintf(title, "zasel0,runr%d", i + 1);
		hiz[i]->Write(title);
	}
	MyFile->Close();
	//h[9]->Draw();
	//hg[4]->SetLineColor(2);
	//hg[4]->Draw("same");
}

//for fit, starts here
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
		(p1 + p2 * (x + xr) + p3 * (x * x + x * xr + xr * xr + 2. * sr * sr)) +
		(1. + erf((arg))) / 2. *
		(p0 + x * (p1 + x * (p2 + p3 * x)) + (p2 + 3. * p3 * x) * sr * sr);
}

double aerogel(size_t n, double x, double* par)
{
	std::vector<double> xv(4);
	std::vector<double> yv(4);
	for (size_t i = 0; i < 4; i++)
		yv[i] = par[4 * n + i];
	double s1 = 0, s2 = 0;
	switch (n) {
	case 0:
		xv[0] = par[12];
		xv[3] = par[13] - par[20] / 120.0;
		s1 = par[17];
		s2 = par[18];
		break;
	case 1:
		xv[0] = par[13] + par[20] / 120.0;
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
	if(x > smin - 3 * s1 && x < smax + 3 * s2){
		return p3g(x, xv[0], s1, pv, smin, smax) -
			p3g(x, xv[3], s2, pv, smin, smax);
	}
	else{
		return 0.;
	}
}

double shifter(double x, double* par)
{
	double ds = par[20] / 120.0;
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
	smin = smin + 1;
	smax = smax + 1;
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
	//par[17]=par[18];
	//par[19]=par[18];
	for (size_t i = 0; i < 3; i++)
		f += aerogel(i, x[0], par);
	f += shifter(x[0], par);
	return f;
}
double scphi2(double* x, double* par)
{
	double f = 0;
	for (size_t i = 0; i < 3; i++)
		f += aerogel(i, x[0], par);
	return f;
}
//ends here

void linesa() {
	//char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	//double scal = PI / 180;
	//double dphi = 2.4 * scal;
	TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/profiles1.root");
	TIter next(f->GetListOfKeys());
	TKey* key;
	while ((key = (TKey*)next())) {
		TString name(key->GetName());
		TString type = "sensor";
		if (name.BeginsWith(type) && name.Contains("sensor8")) {
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

void CalibrExp::testforfit(string basefile, string pName ) {
	readTrueGran1();
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	double shir[9][14][3];
	std::vector<double> xv1(12);
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* f = new TFile((workingDir + basefile).c_str());
	TFile* f1 = new TFile((workingDir + "fitted" + basefile).c_str(), "RECREATE");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0;
	string pNameFullsPred = "a";
	while ((key = (TKey*)next()) && counterr < 500) {
		counterr += 1;
		TString name(key->GetName());
		if (/*name.BeginsWith(type) && */name.Contains(pName)) {
			int obl, counter;

			// ____find the name of a file 
			string pNameFulls;
			const char* pNameFull = "a";
			if (pName == "zobl") {
				sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
				pNameFulls = Form("zobl%d,sensor%d",  obl,  counter);
				pNameFull = pNameFulls.c_str();
			}
			else if (pName == "mean") {
				sscanf(name.Data(), "mean,sensor%d", &counter);
				pNameFulls = Form("mean,sensor%d",  counter);
				pNameFull = pNameFulls.c_str();
				obl = 1;
			}

			// ____check if this name is the same as pred
			if (pNameFulls == pNameFullsPred)
				continue;
			pNameFullsPred = pNameFulls;

			TDirectory* dir = f1->mkdir(pNameFull);
			TProfile* h_amp_phi = (TProfile*)f->Get(pNameFull);

			h_amp_phi->Draw();
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			std::vector<double> par(21);
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
			xmax = fxc->GetParameter(1) - 5 * fxc->GetParameter(2);//((counter-1)*40.+10.)*scal+dphi;
			TF1* fxl = new TF1("fxl", gp2x, xmin, xmax, 5);
			fxl->SetParameter(0, ((counter - 1) * 40 + 4.5) * scal);
			fxl->SetParameter(1, 0.5 * scal);
			fxl->SetParLimits(1, 0.1 * scal, 2.2 * scal);
			fxl->SetParameter(2, 1.0);
			fxl->SetParameter(3, 1.0);
			fxl->SetParameter(4, 2.0);
			h_amp_phi->Fit(fxl, "", "", xmin, xmax);
			c->Update();
			//cin.get();
			//cin.get();
			//
			xmin = fxc->GetParameter(1) + 6 * fxc->GetParameter(2);//(counter*40.-10.)*scal+dphi;
			xmax = (counter * 40. + 10.) * scal + dphi;
			TF1* fxr = new TF1("fxr", gp2x, xmin, xmax, 5);
			fxr->SetParameter(0, (counter * 40 + 2.5) * scal);
			fxr->SetParLimits(0, ((counter-0.5) * 40 + 2.5) * scal, ((counter+0.5) * 40 + 2.5) * scal);
			fxr->SetParameter(1, -0.5 * scal);
			fxr->SetParLimits(1, -1.2 * scal, -0.3 * scal);
			fxr->SetParameter(2, 5);
			fxr->SetParameter(3, 0);
			fxr->SetParameter(4, 0);
			h_amp_phi->Fit(fxr, "", "", xmin, xmax);
			c->Update();
			//cin.get();
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
			par[16] = (fxc->GetParameter(0) + fxc->GetParameter(3));
			par[17] = fxl->GetParameter(1)/1.2;
			par[18] = fxc->GetParameter(2)/1.5;
			par[19] = fabs(fxr->GetParameter(1))/1.2;
			//
			par[20] = 1.5;

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

			TF1* fn = new TF1("scphi", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 21);
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
			fn->SetParName(20, "shifterWidth");
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
			fn->FixParameter(16, fn->GetParameter(16));
			fn->SetParLimits(16, fxc->GetParameter(0) - 0.3, fxc->GetParameter(0) + 19.);
			fn->FixParameter(17, fn->GetParameter(17));
			fn->FixParameter(18, fn->GetParameter(18));
			fn->FixParameter(19, fn->GetParameter(19));
			fn->FixParameter(20, fn->GetParameter(20));


			fn->ReleaseParameter(20);
			fn->SetParLimits(20, 1.4, 4.0);

			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();
			//fn->FixParameter(16, fn->GetParameter(16));
			//fn->ReleaseParameter(16);
			//fn->SetParLimits(16, fn->GetParameter(16)-1.0, fn->GetParameter(16)*5.0);

			for (size_t i = 0; i < 4; i++)
				fn->ReleaseParameter(i);
			fn->SetParLimits(0,0.5,5.);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();
			for (size_t i = 4; i < 8; i++)
				fn->ReleaseParameter(i);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();
			for (size_t i = 8; i < 12; i++)
				fn->ReleaseParameter(i);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();

			//fn->ReleaseParameter(16);
			fn->ReleaseParameter(14);

			
			fn->ReleaseParameter(12);
			fn->ReleaseParameter(13);
			fn->ReleaseParameter(15);
			fn->SetParLimits(12, fn->GetParameter(12) - 2*dphi, fn->GetParameter(12) + 2*dphi);
			fn->SetParLimits(13, fn->GetParameter(13) - dphi, fn->GetParameter(13) + dphi);
			fn->SetParLimits(14, fn->GetParameter(14) - 0.5*dphi, fn->GetParameter(14) + 0.5*dphi);
			fn->SetParLimits(15, fn->GetParameter(15) - 2*dphi, fn->GetParameter(15) + 2*dphi);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();

			fn->ReleaseParameter(17);
			fn->ReleaseParameter(18);
			fn->ReleaseParameter(19);
			fn->SetParLimits(17, 2.0e-3, 1.4e-2);
			fn->SetParLimits(18, 1.0e-3, 1.2e-2);
			fn->SetParLimits(19, 1.0e-3, 2.0e-2);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();

			fn->ReleaseParameter(16);
			fn->SetParLimits(16, fn->GetParameter(16)-10.0, fn->GetParameter(16)*5.0);

			for (size_t i = 0; i < 2; i++) {
				xmin = fn->GetParameter(12) - 3 * fn->GetParameter(17);
				xmax = fn->GetParameter(15) + 3 * fn->GetParameter(19);
				h_amp_phi->Fit("scphi", "", "", xmin, xmax);
			}

			//___writing__

			fn->GetParameters(&par[0]);
			for (size_t i = 0; i < 21; i++)
				fitPar[counter - 1][obl - 1][i] = par[i];
			
			shir[counter-1][obl-1][0] = fn->GetParameter(17);
			shir[counter-1][obl-1][1] = fn->GetParameter(18);
			shir[counter-1][obl-1][2] = fn->GetParameter(19);

			double granerrorst[4];
			granerrorst[0] = shir[counter-1][obl-1][0]/2.;
			granerrorst[1] = (shir[counter-1][obl-1][0] + shir[counter-1][obl-1][1])/4.;
			granerrorst[2] = (shir[counter-1][obl-1][1] + shir[counter-1][obl-1][2])/4.;
			granerrorst[3] = shir[counter-1][obl-1][2]/2.;

			for (int ii = 0; ii < 4; ii++) {
				gran1[ii][counter - 1][obl - 1] = fn->GetParameter(12 + ii);
				//gran1Err[ii][counter - 1][obl - 1] = fn->GetParError(12 + ii);
				gran1Err[ii][counter - 1][obl - 1] = granerrorst[ii];
			}
			pik[counter-1][obl-1] = fn->GetParameter(16);

			cout << "ready" << endl;
			//cin.get();
			
			
			//drawing
			{
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

				//fn->Draw("same");
				for (size_t cc = 0; cc < 3; cc++) {
					TF1* fna = new TF1(Form("aerogel%d", cc), scphi, par[12] - 15 * scal, par[15] + 15 * scal, 21);
					fna->SetParameters(fn->GetParameters());
					for (size_t i = 0; i < 12; i++)
						fna->SetParameter(i, 0);
					for (size_t i = cc * 4; i < cc * 4 + 4; i++)
						fna->SetParameter(i, fn->GetParameter(i));
					fna->SetParameter(16, 0);
					fna->SetLineColor(kBlue);
					fna->SetLineWidth(0);
					fna->SetNpx(1000);
					fna->Draw("same");
					//
					dir->cd();
					cout << "a" << endl;
					sprintf(titlea, "zobl%d,sensor%d,aerogel%d", obl, counter, cc + 1);
					//sprintf(titlea, "date%d,sensor%d,aerogel%d", date, counter, c + 1);
					cout << "b" << endl;
					fna->Write(titlea);
					cout << "c" << endl;
					c->Update();
					cout << "d" << endl;
				}
				c->Update();
				//
				TF1* fns = new TF1("shifter", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 21);
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
			}
			//cin.get();
		}
	}
	cout << "I am here" << endl;
	f1->Write();
	f1->Close();
	if (basefile == "profiles.root") {
		writeGran1();
		for (size_t a = 0; a < 9; a++) {
			meanshir[a][0] = 0.;
			meanshir[a][1] = 0.;
			meanshir[a][2] = 0.;
			for (size_t i = 0; i < 14; i++) {
				meanshir[a][0] += shir[a][i][0]/14.;
				meanshir[a][1] += shir[a][i][1]/14.;
				meanshir[a][2] += shir[a][i][2]/14.;
			}
		}
		writeMeanshir();
		writeFitPar();
	}
	if (basefile == "profilesmmod.root") {
		double oblast[14];
		double oblastErr[14];
		for (int i = 0; i < 14; i++) {
			oblast[i] = (double)i + 1;
			oblastErr[i] = 0.;
		}
		for (size_t i = 0; i < 9; i++) {
			for (size_t j = 0; j < 4; j++) {
				truegran1m[i][j] = gran1[j][i][0];
				/*TGraphErrors* graph = new TGraphErrors(14, oblast, gran1[j][i], oblastErr, gran1Err[j][i]);
				TF1* fn = new TF1("fn", "[0]", 1, 14);
				graph->Fit(fn);
				truegran1m[i][j] = fn->GetParameter(0);
				truegran1mErr[i][j] = fn->GetParError(0);*/
				//for(size_t m = 2; m < 12; m++)
				//	truegran1m[i][j] += gran1[j][i][m]/10.;
			}
		}
		writeTruegran1m();

		for (size_t a = 0; a < 9; a++) {
			meanshirm[a][0] = 0.;
			meanshirm[a][1] = 0.;
			meanshirm[a][2] = 0.;
			for (size_t i = 2; i < 12; i++) {
				meanshirm[a][0] += shir[a][0][0] / 10.;
				meanshirm[a][1] += shir[a][0][1] / 10.;
				meanshirm[a][2] += shir[a][0][2] / 10.;
			}
		}
		writeMeanshirm();

		/*//info about difference in grans 
		double scal1 = 180. / PI;
		double maximum[9][4];
		vector<int> badfit1[9][4];
		for (int j = 0; j < 9; j++) {
			for (int i = 0; i < 4; i++) {
				maximum[j][i] = 0.;
				for (int k = 0; k < 14; k++) {
					if (fabs(gran1[i][j][k]-truegran1m[j][i]) > maximum[j][i])
						maximum[j][i] = fabs(gran1[i][j][k]-truegran1m[j][i]);
					if (fabs(gran1[i][j][k]-truegran1m[j][i])*scal1>0.9){
						badfit1[j][i].push_back(k+1);
					}
				}
			}
		}
		//cout what above
		cout << "rasbros" << endl << "counter	|	1gran	2gran	4gran" << endl;
		for (int j = 0; j < 9; j++) {
			cout << "  " << j + 1 << "	|";
			for (int i = 0; i < 4; i++) {
				if(i!=2){
					cout << "	" << Form("%.1f",(maximum[j][i])*scal1) << Form("/%.1f- ",(gran1Err[i][j][0])*scal1);
					for(size_t k = 0; k < badfit1[j][i].size(); k++){
						cout << badfit1[j][i][k] << " ";
					}
				}
			}
			cout << endl;
		}
		cout << endl;*/
	}

}

//find mean by z borders of aerogel parts truegran1 from obtained in fit gran1.
void CalibrExp::bordersplot() {
	readGran1();
	readFitPar();
	char title[100];
	TFile* f1 = new TFile((workingDir + "bordersfromfit.root").c_str(), "RECREATE");
	double scal1 = 180. / PI;
	double oblast[14];
	double oblastErr[14];
	for (int i = 0; i < 14; i++) {
		oblast[i] = (double)i + 1;
		oblastErr[i] = 0.;
	}
	double maximum[9][4];
	double minimum[9][4];
	vector<int> badfit1[9][4];
	TGraphErrors* graph[4][9];
	TF1* fn = new TF1("fn", "[0]", 1, 14);
	fn->SetLineColor(kRed);
	for (int j = 0; j < 9; j++) {
		cout << "  " << j + 1 << "	|";
		for (int i = 0; i < 4; i++) {
			maximum[j][i] = 0.; minimum[j][i] = 10.;
			graph[i][j] = new TGraphErrors(14, oblast, gran1[i][j],oblastErr,gran1Err[i][j]);
			graph[i][j]->SetLineColor(1);
			graph[i][j]->SetTitle(Form("counter%d,gran%d;zobl;phi",j+1,i+1));
			graph[i][j]->SetMarkerStyle(21);
			graph[i][j]->SetMarkerSize(1);
			graph[i][j]->SetMarkerColor(kBlack);
			graph[i][j]->Draw("APL");
			sprintf(title, "sens%d,gran%d", j+1, i+1);
			graph[i][j]->Write(title);
			graph[i][j]->Fit(fn);
			cout << fn->GetParameter(0) << endl;
			truegran1[j][i] = fn->GetParameter(0);
			truegran1Err[j][i] = fn->GetParError(0);
			for (int k = 0; k < 14; k++) {
				if (fabs(gran1[i][j][k]-truegran1[j][i]) > maximum[j][i]) maximum[j][i] = fabs(gran1[i][j][k]-truegran1[j][i]);
				if (fabs(gran1[i][j][k]-truegran1[j][i])*scal1>0.9){
					badfit1[j][i].push_back(k+1);
				}
				if (gran1[i][j][k] < minimum[j][i]) minimum[j][i] = gran1[i][j][k];
			}
			cout << "	" << Form("%.1f",(maximum[j][i])*scal1);
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			//cin.get();
		}
		
	}
	f1->Close();
	writeTrueGran1();

	//info about changes in 
	cout << "rasbros" << endl << "counter	|	1gran	2gran	4gran" << endl;
	for (int j = 0; j < 9; j++) {
		cout << "  " << j + 1 << "	|";
		for (int i = 0; i < 4; i++) {
			if(i!=2){
				cout << "	" << Form("%.1f",(maximum[j][i])*scal1) << Form("/%.1f- ",(gran1Err[i][j][0])*scal1);
				for(size_t k = 0; k < badfit1[j][i].size(); k++){
					cout << badfit1[j][i][k] << " ";
				}
			}
		}
		cout << endl;
	}
	cout << endl;

	
	double dphi = 2.4 * scal1;
	TFile* f = new TFile((workingDir + "profiles.root").c_str());
	TIter next(f->GetListOfKeys());
	TKey* key;
	double maxdiffArrL[9][14];
	double maxdiffArrR[9][14];
	double maxdiffArrM[9][14];
	int counterr = 0;
	while ((key = (TKey*)next()) && counterr < 500) {
		counterr += 1;
		TString name(key->GetName());
		TString type = "zobl";
		if (/*name.BeginsWith(type) && */name.Contains("sensor")) {
			int obl, counter;
			sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
			TDirectory* dir = f1->mkdir(Form("zobl%d,sensor%d", obl, counter));
			TProfile* h_amp_phi = (TProfile*)f->Get(Form("zobl%d,sensor%d", obl, counter));
			std::vector<double> par;
			for(size_t k =0;k<20;k++){
				par.push_back(fitPar[counter-1][obl-1][k]);
			}
			TF1* fn = new TF1("scphi", scphi, par[12] - 15 * scal1, par[15] + 15 * scal1, 20);
			double maxdiffL = 0;
			double maxdiffR = 0; 
			double maxdiffM = 0;
			fn->SetParameters(par.data());
			for(size_t i = 1; i < h_amp_phi->GetNbinsX()-1; i++){
				double x0 = h_amp_phi->GetXaxis()->GetBinCenter(i);
				double x1 = h_amp_phi->GetXaxis()->GetBinCenter(i-1);
				double x2 = h_amp_phi->GetXaxis()->GetBinCenter(i+1);
				double erri = h_amp_phi->GetBinError(i);
				double hy0 = h_amp_phi->GetBinContent(i);
				double hy1 = h_amp_phi->GetBinContent(i-1);
				double hy2 = h_amp_phi->GetBinContent(i+1);
				double fy0 = scphi(&x0,&par[0]);
				double fy1 = scphi(&x1,&par[0]);
				double fy2 = scphi(&x2,&par[0]);
				double bord1 = gran1[0][counter-1][obl-1] + 3.*gran1Err[0][counter-1][obl-1];
				double bord2 = gran1[1][counter-1][obl-1] - 3.*gran1Err[1][counter-1][obl-1];
				double bord3 = gran1[1][counter-1][obl-1] + 3.*gran1Err[1][counter-1][obl-1];
				double bord4 = gran1[3][counter-1][obl-1] - 3.*gran1Err[3][counter-1][obl-1];

				if(fabs(x0-gran1[0][counter-1][obl-1]) < 3.*gran1Err[0][counter-1][obl-1] && maxdiffL < fabs(hy0 - fy0)/erri ) {
					maxdiffL = fabs(hy0 - fy0)/erri;
				}
				if(fabs(x0-gran1[3][counter-1][obl-1]) < 3.*gran1Err[3][counter-1][obl-1] && maxdiffR < fabs(hy0 - fy0)/erri ) {
					maxdiffR = fabs(hy0 - fy0)/erri;
				}
				if(((x0 > bord1 && x0 < bord2) || (x0 > bord3 && x0 < bord4))  ) {
					double locmeandiff = (hy0 - fy0 + hy1 - fy1 + hy2 - fy2)/(3.*erri);
					if (maxdiffM < fabs(locmeandiff))
						maxdiffM = fabs(locmeandiff);
				}
			}
			maxdiffArrL[counter-1][obl-1] = maxdiffL;
			maxdiffArrR[counter-1][obl-1] = maxdiffR;
			maxdiffArrM[counter-1][obl-1] = maxdiffM;
		}
	}
	cout << "rasbrosL" << endl << "counter	|	z1	z2	z3	z4	z5	z6	z7	z8	z9	z10	z11	z12	z13	z14" << endl;
	for (int i = 0; i < 9; i++) {
		cout << "  " << i + 1 << "	|";
		for (int j = 0; j < 14; j++) {
			cout << "	" << Form("%.1f", maxdiffArrL[i][j]);
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 9; i++) {
		cout << "  " << i + 1 << "	|";
		for (int j = 0; j < 14; j++) {
			cout << "	" << Form("%.1f", maxdiffArrR[i][j]);
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 9; i++) {
		cout << "  " << i + 1 << "	|";
		for (int j = 0; j < 14; j++) {
			cout << "	" << Form("%.1f", maxdiffArrM[i][j]);
		}
		cout << endl;
	}
	f->Close();
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

//find shifter array from pik, stored in gran1 file and also use it, Yfm and truegran to create file with parameters for modeling database
void CalibrExp::modelfile(bool transformPYN) {
	readTrueGran1();
	readTruegran1m();
	double scal = PI / 180;
	//shifter fit
	readGran1();
	readNcoef();
	double oblast[14];
	for (int i = 0; i < 14; i++)
		oblast[i] = (double)i + 1;
	TGraph* grpik[9];
	for (int j = 0; j < 9; j++) {
		grpik[j] = new TGraph(14, oblast, pik[j]);
	}
	TF1* fn = new TF1("fn", "[0]*exp(-x/[2])+[1]*exp(x/[2])", 1, 14);
	//TF1* fn = new TF1("fn", "[0]+[1]*x", 1, 14);
	fn->SetParameters(1., 2., 3.);
	fn->SetParLimits(0, 0.1, 1000);
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
		grpik[j]->Fit("fn", "", "", 2, 14);
		c->Update();
		//cin.get();
		shift[j][0] = fn->GetParameter(0);
		shift[j][1] = fn->GetParameter(1);
		shift[j][2] = fn->GetParameter(2);
	}
	writeShift();
	readShift();
	readTrueGran1();
	readFitPar();
	double scal1 = 180. / PI;
	ofstream fout1;
	std::vector<double> xv(12);
	std::vector<double> xv1(12);
	std::vector<double> xvt(12);
	std::vector<double> xv2(4);
	std::vector<double> yv2(4);
	int n = 14; 
	fout1.open((workingDir + "map.txt").c_str());
	for (int j = 0; j < 9; j++) {
		xvt[0] = truegran1[j][0];
		xvt[3] = truegran1[j][1];
		xvt[1] = xvt[0] + (xvt[3] - xvt[0]) / 3;
		xvt[2] = xvt[3] - (xvt[3] - xvt[0]) / 3;
		xvt[4] = truegran1[j][1];
		xvt[7] = truegran1[j][2];
		xvt[5] = xvt[4] + (xvt[7] - xvt[4]) / 3;
		xvt[6] = xvt[7] - (xvt[7] - xvt[4]) / 3;
		xvt[8] = truegran1[j][2];
		xvt[11] = truegran1[j][3];
		xvt[9] = xvt[8] + (xvt[11] - xvt[8]) / 3;
		xvt[10] = xvt[11] - (xvt[11] - xvt[8]) / 3;
		//Yfm calculation
		xv1[0] = truegran1m[j][0];
		xv1[3] = truegran1m[j][1];
		xv1[1] = xv1[0] + (xv1[3] - xv1[0]) / 3;
		xv1[2] = xv1[3] - (xv1[3] - xv1[0]) / 3;
		xv1[4] = truegran1m[j][1];
		xv1[7] = truegran1m[j][2];
		xv1[5] = xv1[4] + (xv1[7] - xv1[4]) / 3;
		xv1[6] = xv1[7] - (xv1[7] - xv1[4]) / 3;
		xv1[8] = truegran1m[j][2];
		xv1[11] = truegran1m[j][3];
		xv1[9] = xv1[8] + (xv1[11] - xv1[8]) / 3;
		xv1[10] = xv1[11] - (xv1[11] - xv1[8]) / 3;
		for (size_t k = 0; k < 14; k++) {
			xv[0] = fitPar[j][k][12];
			xv[3] = fitPar[j][k][13] - fitPar[j][k][20] / 120.0;
			xv[1] = xv[0] + (xv[3] - xv[0]) / 3;
			xv[2] = xv[3] - (xv[3] - xv[0]) / 3;
			xv[4] = fitPar[j][k][13] + fitPar[j][k][20] / 120.0;
			xv[7] = fitPar[j][k][14];
			xv[5] = xv[4] + (xv[7] - xv[4]) / 3;
			xv[6] = xv[7] - (xv[7] - xv[4]) / 3;
			xv[8] = fitPar[j][k][14];
			xv[11] = fitPar[j][k][15];
			xv[9] = xv[8] + (xv[11] - xv[8]) / 3;
			xv[10] = xv[11] - (xv[11] - xv[8]) / 3;
			for (int n = 0; n < 3; n++) {
				for (size_t i = 0; i < 4; i++) {
					yv2[i] = fitPar[j][k][4 * n + i];
					xv2[i] = xv[4 * n + i];
				}
				vector<double> pv1 = m4(xv2, yv2);
				//cout << pv1[0] << "	" << pv1[1] << "	" << pv1[2] << "	" << pv1[3] << endl;
				//cout << pv1[0] + pv1[1] * x + pv1[2] * pow(x, 2) + pv1[3] * pow(x, 3) << endl;
				for (int i = 0; i < 4; i++) {
					double x1 = xvt[i + 4 * n];
					if (transformPYN)
						x1 = xv1[i + 4 * n];
					//just polinoms
					yfm[j][k][4 * n + i] = pv1[0] + pv1[1] * x1 + pv1[2] * pow(x1, 2) + pv1[3] * pow(x1, 3);
					//yfm[j][k][4 * n + i] = 2-i%2;
					//cout << yfm[j][k][4 * n + i] / fitPar[j][k][4 * n + i] << endl;
					//aerogel with exponents (fit f without shifters)
					//yfm[j][k][4 * n + i] = aerogel(n, x1, fitPar[j][k]);
					//cout << yfm[counter - 1][obl - 1][4 * n + i] << endl;
				}
			}
			TF1* fn = new TF1("scphi", scphi, fitPar[j][k][12] - 15 * scal, fitPar[j][k][15] + 15 * scal, 21);
			for (size_t i = 0; i < 21; i++)
				fn->SetParameter(i, fitPar[j][k][i]);
			cout << fitPar[j][k][16] << endl;;
			fn->Draw();
			TF1* fn1 = new TF1("scphi1", scphi, fitPar[j][k][12] - 15 * scal, fitPar[j][k][15] + 15 * scal, 21);
			for (size_t i = 0; i < 12; i++)
				fn1->SetParameter(i, yfm[j][k][i]);
			for (size_t i = 12; i < 16; i++)
				fn1->SetParameter(i, truegran1[j][i-12]);	
			for (size_t i = 16; i < 20; i++)
				fn1->SetParameter(i, fitPar[j][k][i]);
			fn1->SetParameter(20, 0);
			fn1->SetLineColor(2);
			fn1->Draw("same");
			TLine* line1 = new TLine(fitPar[j][k][12], 0, fitPar[j][k][12], 5);
			TLine* line2 = new TLine(fitPar[j][k][13], 0, fitPar[j][k][13], 5);
			TLine* line3 = new TLine(truegran1[j][0], 0, truegran1[j][0], 5);
			TLine* line4 = new TLine(truegran1[j][1], 0, truegran1[j][1], 5);
			line3->SetLineColor(2);
			line4->SetLineColor(2);
			line1->Draw("same");
			line2->Draw("same");
			line3->Draw("same");
			line4->Draw("same");
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			//cin.get();

		}
		//write in file
		fout1 << Form("Counter %d", j+1) << endl << Form("Shifter %f %f %f",shift[j][0], shift[j][1], shift[j][2]) << endl << "AmplitudeCorrection 1" << endl << "Aerogel" << endl;
		for (int k = 0; k < 3; k++) {
			fout1 << n << " " << xvt[0+4*k]*scal1 << " " << xvt[1+4*k]*scal1 << " " << xvt[2+4*k]*scal1 << " " << xvt[3+4*k]*scal1 << endl;
			for (int i = 0; i < 14; i++) {
				fout1 << zobl[i]-0.0 << " " << yfm[j][i][0+4*k]/coef[i][j] << " " << yfm[j][i][1+4*k]/ coef[i][j] << " " << yfm[j][i][2+4*k]/ coef[i][j] << " " << yfm[j][i][3+4*k]/ coef[i][j] << endl;
			}
		}
	}
	fout1.close();
	writeMcoef();
}

void CalibrExp::checkmod() {
	//char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	std::vector<double> xv1(12);
	double scal = PI / 180;
	//double dphi = 2.4 * scal;
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

void CalibrExp::checkbord() {
	readTrueGran1();
	//char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	std::vector<double> xv1(12);
	//double scal = PI / 180;
	//double dphi = 2.4 * scal;
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
			int counter;
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

double CalibrExp::compare1() {
	readGoodRuns();
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
	readTruegran1m();
	TChain chain("t1");
	std::string chadds = dirSelectedExp + prefixSelectedExp + "**exp00**.root";
	chain.Add(chadds.c_str());

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

	TFile* f00 = new TFile((workingDir + "phintegr_pedest.root").c_str());
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			TH1* hped = (TH1F*)f00->Get(Form("ped_phintegr_zobl%d,sensor%d", j + 1, i + 1));
			ped[j][i] = hped->GetMean(1);
		}
	}
	int Count = 0, count1 = 0;
	float ztr[2] = { 0.,0. };
	float zin[2] = { 0.,0. };
	int maxampla = 0;
	double maxach = 0.;
	double lgran[4];
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr[0] = 12.0 / tan(theta[0]) + z0[0];
		ztr[1] = 12.0 / tan(theta[1]) + z0[1];
		zin[0] = aerRLowBord / tan(theta[0]) + z0[0];
		zin[1] = aerRLowBord / tan(theta[1]) + z0[1];
		schr[nch] = 0.;
		if ((!RunIsGood(run))/* || eventtime != 0*/) {
			continue;
		}
		scount = 0;     //2017: 4001, 25500, 29501    ; 2018: 6001, 29800, 35801       ; 2019: 38000-43001       ; 2020: 45040-46541
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			//wichcounter1[j] = scount1;
			if (j == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
				if (phi[1] < 1.)
					phi[1] = phi[1] + 2. * PI;
				while (1 + scount < nch) { scount += 1; }
			}
			else
				while (schr[1 + scount] > 0.5) { scount += 1; }

			maxach = -100.;
			maxampla = scount1;
			for (int f = scount1; f <= scount; f++) {
				if (ach[f] > maxach) {
					maxach = ach[f];
					maxampla = f;
				}
			}

			lgran[0] = max(truegran1[j][0], truegran1m[j][0]) + 3. * max(meanshir[j][0], meanshirm[j][0]);
			lgran[1] = min(truegran1[j][1], truegran1m[j][1]) - 5. * max(meanshir[j][1], meanshirm[j][1]);// -0.045;
			lgran[2] = max(truegran1[j][1], truegran1m[j][1]) + 5. * max(meanshir[j][1], meanshirm[j][1]);// +0.045;
			lgran[3] = min(truegran1[j][3], truegran1m[j][3]) - 3. * max(meanshir[j][2], meanshirm[j][2]);// -0.03;

			
			for (size_t pind = 0; pind < 2; pind++) {
				int zoblInd = -1;
				for (int i = 0; i < 14; i++) {
					if ((ztr[pind] >= zobl[i]) && (ztr[pind] < zobl[i + 1]))
						zoblInd = i;
				}
				if ((((phi[pind] > lgran[0]) && (phi[pind] < lgran[1])) || ((phi[pind] > lgran[2]) && (phi[pind] < lgran[3]))) && (TimeIsGood(maxampla)) && (TimeIsGood(scount1)) && zoblInd != -1 && ach[maxampla] < achCut) {
					ach1[maxampla] = ach[maxampla] - ped[zoblInd][j];
					//ach1[maxampla] = ach[maxampla] - ach[scount1];
					meancampl[zoblInd]->Fill(ach1[maxampla]);
					if (zoblInd > 2 && zoblInd < 12)
						meanzampl[j]->Fill(ach1[maxampla]);

					double fillAmplitude = recalculateZ(zin[pind], theta[pind], ach[maxampla], ped[zoblInd][j]);
					if (fillAmplitude != -1.0 && ach[maxampla] < achCut )
						h[j][zoblInd]->Fill(fillAmplitude);
					if (ach[maxampla] >= 0.2) {
						if (zoblInd > 2 && zoblInd < 12)
							ecountz1[j] += 1;
						ecounta1[j][zoblInd] += 1;
						ecountc1[zoblInd] += 1;
					}
					if (ach[maxampla] < 0.2) {
						if (zoblInd > 2 && zoblInd < 12)
							ecountz2[j] += 1;
						ecounta2[j][zoblInd] += 1;
						ecountc2[zoblInd] += 1;
					}
				}
			}
			scount += 1;
			scount1 = scount;
		}
		
		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("obrabotano %d M entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}
	{
		TFile* MyFile = new TFile((workingDir + "achspectreme.root").c_str(), "RECREATE");
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
	writeComparisone();

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

double CalibrExp::compare2(bool transformPYN) {
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
	readTruegran1m();
	TChain chain("h1");
	std::string chadds = dirSelectedMod + prefixSelectedMod + "ee**.root";
	chain.Add(chadds.c_str());
	
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
	float ztr[2] = { 0.,0. };
	float zin[2] = { 0.,0. };
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		ztr[0] = 12.0 / tan(theta[0]) + z0[0];
		ztr[1] = 12.0 / tan(theta[1]) + z0[1];
		zin[0] = aerRLowBord / tan(theta[0]) + z0[0];
		zin[1] = aerRLowBord / tan(theta[1]) + z0[1];
		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			//wichcounter1[j] = scount1;
			if (j == 8) {
				if (phi[0] < 1.)
					phi[0] = phi[0] + 2. * PI;
				if (phi[1] < 1.)
					phi[1] = phi[1] + 2. * PI;
			}

			lgran[0] = max(truegran1[j][0], truegran1m[j][0]) + 3. * max(meanshir[j][0], meanshirm[j][0]);
			lgran[1] = min(truegran1[j][1], truegran1m[j][1]) - 5. * max(meanshir[j][1], meanshirm[j][1]);// -0.045;
			lgran[2] = max(truegran1[j][1], truegran1m[j][1]) + 5. * max(meanshir[j][1], meanshirm[j][1]);// +0.045;
			lgran[3] = min(truegran1[j][3], truegran1m[j][3]) - 3. * max(meanshir[j][2], meanshirm[j][2]);// -0.03;
			
			for (size_t pind = 0; pind < 2; pind++) {
				int zoblInd = -1;
				for (int i = 0; i < 14; i++) {
					if ((ztr[pind] >= zobl[i]) && (ztr[pind] < zobl[i + 1]))
						zoblInd = i;
				}
				if ((((phi[pind] > lgran[0]) && (phi[pind] < lgran[1])) || ((phi[pind] > lgran[2]) && (phi[pind] < lgran[3]))) && zoblInd != -1 && ach[j]< achCut) {

					meancampl[zoblInd]->Fill(ach[j]);
					if (zoblInd > 2 && zoblInd < 12)
						meanzampl[j]->Fill(ach[j]);

					double fillAmplitude = recalculateZ(zin[pind], theta[pind], ach[j], 0);
					if (fillAmplitude != -1.0 && ach[j] < achCut)
						h[j][zoblInd]->Fill(fillAmplitude);
					if (ach[j] >= 0.2) {
						if (zoblInd > 2 && zoblInd < 12)
							countz1[j] += 1;
						counta1[j][zoblInd] += 1;
						countc1[zoblInd] += 1;
					}
					if (ach[j] < 0.2) {
						if (zoblInd > 2 && zoblInd < 12)
							countz2[j] += 1;
						counta2[j][zoblInd] += 1;
						countc2[zoblInd] += 1;
					}
				}
			}

			scount += 1;
			scount1 = scount;
		}
		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("obrabotano %d M entries", count1) << endl;
			Count = 0;
		}
	}

	{
		TFile* MyFile = new TFile((workingDir + "achspectreme.root").c_str(), "UPDATE");
		for (int i = 0; i < 14; i++) {
			sprintf(title, "m_zobl%d", i + 1);
			meancampl[i]->Write(title);
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
	writeComparisonm();
	return 2;
}

void CalibrExp::compare() {
	//read counts for efficiency
	//compare1();
	//compare2();
	readComparisonm();
	readComparisone();
	readMcoef();
	//cout << "popr = " << compare1() << endl;
	double meanample[14][9];
	double meanamplm[14][9];
	double meanampleerr[14][9];
	double meanamplmerr[14][9];

	TFile* Myf = new TFile((workingDir + "achspectreme.root").c_str());
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			TH1* meanamplocke = (TH1F*)Myf->Get(Form("e_all_zobl%d_count%d", i + 1, j + 1));
			TH1* meanamplockm = (TH1F*)Myf->Get(Form("m_all_zobl%d_count%d", i + 1, j + 1));
			meanample[i][j] = meanamplocke->GetMean(1);
			meanampleerr[i][j] = meanamplocke->GetMeanError(1);
			meanamplm[i][j] = meanamplockm->GetMean(1);
			meanamplmerr[i][j] = meanamplockm->GetMeanError(1);
		}
	}

	//Pure comparison between modeling and experiment, mean by counters, 14 zobl
	cout << "popr = " << endl;
	double popr = 0.;
	for (int i = 0; i < 14; i++) {
		cout << (meanamplm[i][0] / meanample[i][0]) << "    " << sqrt(pow(meanamplmerr[i][0] / meanample[i][0],2) + pow(meanampleerr[i][0]/ meanample[i][0] * meanamplm[i][0] / meanample[i][0],2)) << "	"/* << meanampleerr[i] */ << endl;
		popr = popr + meanamplm[i][0] / meanample[i][0];
	}
	cout << "mean	" << popr / 14. << endl;

	double oblast[14];
	double noll[14];
	//coeff which were transfered in map.txt before modeling

	double koefnew[9][14];
	double koefnew100[9][14];
	double koeferr[9][14];
	double koeferr100[9][14];
	double lokk = 0.;
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			oblast[i] = (zobl[i] + zobl[i + 1]) / 2.;
			noll[i] = 0.;
			if (coef[i][j] != coef[0][j])
				cout << "ABCDEF" << endl;
			koefnew[j][i] = coef[i][j] * (meanamplm[i][j] / meanample[i][j]);
			//koefnew[i][j] = 100. - 100. / (lokk);
			koefnew100[j][i] = (koefnew[j][i] -1.)*100.;
			koeferr[j][i] = koefnew[j][i] * (meanamplmerr[i][j]/ meanamplm[i][j] + meanampleerr[i][j] / meanample[i][j]) ;
			koeferr100[j][i] = koeferr[j][i] * 100.;
			cout << koefnew[j][i] << ", ";
		}
		cout << endl;
	}
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			coef[i][j] = koefnew[j][i];
		}
	}
	for (int j = 0; j < 9; j++) {
		TGraphErrors* gr = new TGraphErrors(14, &oblast[0], &koefnew100[j][0], &noll[0], &koeferr100[j][0]);
		TF1* fn = new TF1("fn", "[0]", -12, 11);
		fn->SetParameter(0, 7.);
		fn->SetParLimits(0, 2., 11.);
		gr->SetLineColor(1);
		gr->SetMarkerStyle(21);
		gr->SetMarkerColor(1);
		//gr->GetYaxis()->SetRangeUser(1., 12.);
		gr->SetTitle("secondary particles contribution;z,sm;percent");
		gr->Draw("AP");
		fn->SetLineColor(kRed);
		gr->Fit("fn", "", "", -9, 9);
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
		cout << fn->GetParameter(0) << endl;
	}

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
	//efficiency comparison with threshold ach = 0.2
	double effm[14];
	double efferrm[14];
	double effe[14];
	double efferre[14];
	for (int j = 0; j < 14; j++) {
		cout << Form("obl %d", j + 1) << endl;
		effm[j] = 100. * (1. - countc2[j] / (countc1[j] + countc2[j]));
		efferrm[j] = 120. * pow(((countc1[j] * countc2[j]) / pow((countc1[j] + countc2[j]), 3.)), 0.5);
		cout << effm[j] << " +- " << efferrm[j] << endl;
		effe[j] = 100. * (1. - ecountc2[j] / (ecountc1[j] + ecountc2[j]));
		efferre[j] = 150. * pow(((ecountc1[j] * ecountc2[j]) / pow((ecountc1[j] + ecountc2[j]), 3.)), 0.5);
		cout << effe[j] << " +- " << efferre[j] << endl;
		//cout << ((countc2[j] / (countc1[j] + countc2[j])) - (ecountc2[j] / (ecountc1[j] + ecountc2[j]))) / pow(((countc1[j] * countc2[j]) / pow((countc1[j] + countc2[j]), 3.)), 0.5) << endl << endl << endl;
		cout << effm[j] / effe[j] - 1 << "	+-	" << (efferrm[j] + efferre[j]) / effe[j] << endl << endl << endl;
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

	writeNcoef();
}

void CalibrExp::compareAmpSpectr() {
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
	readTruegran1m();
	readGran1();
	readMcoef();
	readFitPar();
	double scal = PI / 180;
	TFile* f1 = new TFile((workingDir + "profilesmmod.root").c_str());
	//TFile* f1 = new TFile((workingDir + "profilesEffMod.root").c_str());
	TFile* f2 = new TFile((workingDir + "profiles.root").c_str());
	TFile* f3 = new TFile((workingDir + "Zaselphi.root").c_str());
	cin.get();
	double koefnewp[9][14];
	double koefnewpErr[9][14];
	double koefnewpM[9];
	double koefnewpMErr[9];
	double koefnewh[9][14];
	double modelGran[9][4] = { 3.96954, 17.4603, 30.8096, 42.9946, 43.6686, 58.0747, 70.5603, 82.8407, 83.7598, 97.8676, 110.007, 122.273, 124.125, 138.028, 150.658, 163.491, 164.016, 178.457, 192.73, 203.046, 203.326, 217.694, 231.488, 242.669, 244.352, 258.508, 271.534, 282.903, 283.639, 298.182, 310.422, 322.809, 323.23, 337.838, 351.363, 362 };
	for (size_t obl = 0; obl < 14; obl++) {
		//TProfile* hprof10 = (TProfile*)f1->Get(Form("zobl%d,sensor%d", zobl + 1, 1));
		//TProfile* hprof20 = (TProfile*)f2->Get(Form("zobl%d,sensor%d", zobl + 1, 1));
		//hprof10->SetLineColor(2);
		//hprof10->Draw();
		//hprof20->Draw("same");
		for (size_t counter = 0; counter < 9; counter++) {
			TF1* fn = new TF1("scphi", scphi, fitPar[counter][obl][12] - 15 * scal, fitPar[counter][obl][12] + 15 * scal, 20);
			for (size_t i = 0; i < 21; i++)
				fn->SetParameter(i, fitPar[counter][obl][i]);


			double lgran[4];
			lgran[0] = max(truegran1[counter][0], truegran1m[counter][0]) + 3. * max(meanshir[counter][0], meanshirm[counter][0]);
			lgran[1] = min(truegran1[counter][1], truegran1m[counter][1]) - 5. * max(meanshir[counter][1], meanshirm[counter][1]);// -0.045;
			lgran[2] = max(truegran1[counter][1], truegran1m[counter][1]) + 5. * max(meanshir[counter][1], meanshirm[counter][1]);// +0.045;
			lgran[3] = min(truegran1[counter][3], truegran1m[counter][3]) - 3. * max(meanshir[counter][2], meanshirm[counter][2]);// -0.03;
			//TDirectory* d1 = (TDirectory*)f2->Get(Form("zobl%d,sensor%d", obl+1, counter+1));
			cout << "asdasdasd" << endl;
			//TF1* tf1 = (TF1*)d1->Get(Form("zobl%d,sensor%d,full", 2, 8));
			//TProfile* hprof1 = (TProfile*)f1->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));
			//TProfile* hprof1 = (TProfile*)f1->Get(Form("mean,sensor%d", counter + 1));
			TProfile* hprof1 = (TProfile*)f1->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));
			TProfile* hprof2 = (TProfile*)f2->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));
			//TProfile* hprof2 = (TProfile*)f2->Get(Form("mean,sensor%d", counter + 1));
			TProfile* hz = (TProfile*)f3->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));
			/*hprof1->SetLineColor(2);
			hprof1->GetXaxis()->SetRangeUser(3./4.*hprof1->GetBinCenter(1)+1./4.*hprof1->GetBinCenter(hprof1->GetNbinsX()-1), 1./4.*hprof1->GetBinCenter(1)+3./4.*hprof1->GetBinCenter(hprof1->GetNbinsX()-1));
			hprof1->Draw();
			hprof2->Draw("same");
			//hz->DrawNormalized("same",800);
			TLine* lineet1 = new TLine(truegran1[counter][0],0,truegran1[counter][0],10);
			TLine* lineel1 = new TLine(gran1[0][counter][obl],0,gran1[0][counter][obl],10);
			//TLine* linemt1 = new TLine(modelGran[counter][0]*PI/180.,0,modelGran[counter][0]*PI/180.,10);
			TLine* linemt1 = new TLine(truegran1m[counter][0],0,truegran1m[counter][0],10);
			TLine* lineet2 = new TLine(truegran1[counter][3],0,truegran1[counter][3],10);
			TLine* lineel2 = new TLine(gran1[3][counter][obl],0,gran1[3][counter][obl],10);
			//TLine* linemt2 = new TLine(modelGran[counter][3]*PI/180.,0,modelGran[counter][3]*PI/180.,10);
			TLine* linemt2 = new TLine(truegran1m[counter][3],0,truegran1m[counter][3],10);
			linemt1->SetLineColor(2);
			linemt2->SetLineColor(2);
			lineet1->SetLineColor(4);
			lineet2->SetLineColor(4);
			lineet1->Draw("same");
			lineet2->Draw("same");
			linemt1->Draw("same");
			linemt2->Draw("same");
			lineel1->Draw("same");
			lineel2->Draw("same");
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			cin.get();*/

			
			TLine* lineg1 = new TLine(lgran[0],0.7,lgran[0],1.3);
			TLine* lineg2 = new TLine(lgran[1],0.7,lgran[1],1.3);
			TLine* lineg3 = new TLine(lgran[2],0.7,lgran[2],1.3);
			TLine* lineg4 = new TLine(lgran[3],0.7,lgran[3],1.3);

			vector<double> x;
			vector<double> dx;
			vector<double> y;
			vector<double> dy;
			vector<double> z;
			vector<double> dz;
			
			for(size_t i = 2; i < hprof2->GetNbinsX()-1; i++){
				if(hprof2->GetBinEntries(i)>0 && ((hprof2->GetBinCenter(i) > lgran[0] && hprof2->GetBinCenter(i) < lgran[1]) || (hprof2->GetBinCenter(i) > lgran[2] && hprof2->GetBinCenter(i) < lgran[3]))){
					x.push_back(hprof2->GetBinCenter(i));
					dx.push_back(hprof2->GetBinWidth(i));

					size_t modBin = hprof1->FindBin(x.back());
					if (modBin != i)
						cout << "bin is not the same for some reason" << endl;
					double modA = 0;
					if (x.back() < hprof1->GetBinCenter(modBin))
						modA = hprof1->GetBinContent(modBin - 1) + (hprof1->GetBinContent(modBin) - hprof1->GetBinContent(modBin - 1)) / (hprof1->GetBinCenter(modBin) - hprof1->GetBinCenter(modBin - 1)) * (x.back()- hprof1->GetBinCenter(modBin - 1));
					else if (x.back() >= hprof1->GetBinCenter(modBin))
						modA = hprof1->GetBinContent(modBin) + (hprof1->GetBinContent(modBin+1) - hprof1->GetBinContent(modBin)) / (hprof1->GetBinCenter(modBin+1) - hprof1->GetBinCenter(modBin)) * (x.back() - hprof1->GetBinCenter(modBin));

					//y.push_back(modA /hprof2->GetBinContent(i));
					//dy.push_back(sqrt(pow(hprof1->GetBinError(i)/hprof2->GetBinContent(i),2) + pow(y.back()*hprof2->GetBinError(i)/hprof2->GetBinContent(i),2)));
					y.push_back(hprof1->GetBinContent(i) / hprof2->GetBinContent(modBin));
					dy.push_back(sqrt(pow(hprof1->GetBinError(i) / hprof2->GetBinContent(modBin), 2) + pow(y.back() * hprof2->GetBinError(modBin) / hprof2->GetBinContent(modBin), 2)));

					z.push_back(modA / scphi(&x.back(),&fitPar[counter][obl][0]));
					dz.push_back(sqrt(pow(hprof1->GetBinError(modBin) / scphi(&x.back(), &fitPar[counter][obl][0]), 2) + pow(z.back() * 0. / scphi(&x.back(), &fitPar[counter][obl][0]), 2)));
					//z.push_back(hprof1->GetBinContent(i) / scphi(&x.back(),&fitPar[counter][obl][0]));
					//dz.push_back(sqrt(pow(hprof1->GetBinError(i) / scphi(&x.back(), &fitPar[counter][obl][0]), 2) + pow(z.back() * 0. / scphi(&x.back(), &fitPar[counter][obl][0]), 2)));
				}
			}

			TGraph* graph = new TGraphErrors(x.size(), &x[0], &y[0], &dx[0], &dy[0]);
			TGraph* graph1 = new TGraphErrors(x.size(), &x[0], &z[0], &dx[0], &dz[0]);
			graph->SetLineColor(1);
			graph1->SetLineColor(1);
			graph->SetTitle(Form("counter%d, zobl%d;x;coeff",counter+1,obl+1));
			graph->SetMarkerStyle(21);
			graph->SetMarkerSize(1);
			graph->SetMarkerColor(kBlack);
			graph->GetYaxis()->SetRangeUser(0.75, 1.25);
			graph1->SetMarkerStyle(21);
			graph1->SetMarkerSize(1);
			graph1->SetMarkerColor(kRed);
			graph1->GetYaxis()->SetRangeUser(0.75, 1.25);
			graph->Draw("AP");
			//graph1->Draw("sameP");
			lineg1->Draw("same");
			lineg2->Draw("same");
			lineg3->Draw("same");
			lineg4->Draw("same");
			TF1* f1 = new TF1("f1", "[0]", lgran[0], lgran[3]);
			f1->SetParameter(0,1.0);
			graph->Fit("f1", "", "", lgran[0], lgran[3]);
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			//cin.get();
			koefnewp[counter][obl] = coef[obl][counter] * f1->GetParameter(0);
			koefnewpErr[counter][obl] = f1->GetParError(0);
			koefnewpM[counter] = coef[obl][counter] * f1->GetParameter(0);
			koefnewpMErr[counter] = coef[obl][counter] * f1->GetParError(0);
		}
	}
	readNcoef();
	double noll[14];
	double oblast[14];
	double koefnewhErr[9][14];
	for (size_t obl = 0; obl < 14; obl++) {
		for (size_t counter = 0; counter < 9; counter++) {
			noll[obl] = 0;
			oblast[obl] = (zobl[obl] + zobl[obl + 1]) / 2.;
			koefnewh[counter][obl] = coef[obl][counter];
			koefnewhErr[counter][obl] = 0.01;
		}
	}
	double oblast9[9] = {1,2,3,4,5,6,7,8,9};
	double nol9[9] = {0,0,0,0,0,0,0,0,0};
	TGraphErrors* gro = new TGraphErrors(9, &oblast9[0], &koefnewpM[0], &nol9[0], &koefnewpMErr[0]);
	gro->SetLineColor(1);
	gro->SetMarkerStyle(21);
	gro->SetMarkerColor(1);
	gro->Draw("AP");

	for (int j = 0; j < 9; j++) {
		TGraphErrors* gr = new TGraphErrors(14, &oblast[0], &koefnewh[j][0], &noll[0], &koefnewhErr[j][0]);
		TGraphErrors* gr1 = new TGraphErrors(14, &oblast[0], &koefnewp[j][0], &noll[0], &koefnewpErr[j][0]);
		TF1* fn = new TF1("fn", "[0]", -12, 11);
		fn->SetParameter(0, 7.);
		fn->SetParLimits(0, 2., 11.);
		gr->SetLineColor(1);
		gr->SetMarkerStyle(21);
		gr->SetMarkerColor(1);
		
		gr1->SetLineColor(2);
		gr1->SetMarkerStyle(21);
		gr1->SetMarkerColor(2);
		//gr->GetYaxis()->SetRangeUser(1., 12.);
		gr->SetTitle("secondary particles contribution;z,sm;percent");
		
		gr->Draw("AP");
		gr1->Draw("Psame");

		fn->SetLineColor(kRed);
		//gr->Fit("fn", "", "", -9, 9);
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
		//cout << fn->GetParameter(0) << endl;
	}
	for (size_t obl = 0; obl < 14; obl++) {
		for (size_t counter = 0; counter < 9; counter++) {
			coef[obl][counter] = (koefnewp[counter][obl] + koefnewh[counter][obl])/2.;
		}
	}
	writeNcoef();
}

void CalibrExp::compareAmpSpectrG() {
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
	readTruegran1m();
	readFitPar();
	std::vector<double> xvt(12);
	xvt[0] = truegran1[4][0];
	xvt[3] = truegran1[4][1];
	xvt[1] = xvt[0] + (xvt[3] - xvt[0]) / 3;
	xvt[2] = xvt[3] - (xvt[3] - xvt[0]) / 3;
	xvt[4] = truegran1[4][1];
	xvt[7] = truegran1[4][2];
	xvt[5] = xvt[4] + (xvt[7] - xvt[4]) / 3;
	xvt[6] = xvt[7] - (xvt[7] - xvt[4]) / 3;
	xvt[8] = truegran1[4][2];
	xvt[11] = truegran1[4][3];
	xvt[9] = xvt[8] + (xvt[11] - xvt[8]) / 3;
	xvt[10] = xvt[11] - (xvt[11] - xvt[8]) / 3;
	TFile* f1 = new TFile((workingDir + "profilesmmod.root").c_str());
	TFile* f2 = new TFile((workingDir + "profiles.root").c_str());
	TProfile* hm0 = (TProfile*)f1->Get(Form("zobl%d,sensor%d", 2, 5));
	TProfile* he0 = (TProfile*)f2->Get(Form("zobl%d,sensor%d", 2, 5));
	hm0->SetLineColor(2);
	//he0->Draw();
	hm0->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	for (size_t obl = 5; obl < 6; obl++) {
		double scal = PI / 180;
		TF1* fn = new TF1(Form("scphi%d",obl), scphi, fitPar[4][obl][12] - 15 * scal, fitPar[4][obl][15] + 15 * scal, 20);
		for (size_t i = 0; i < 12; i++)
			fn->SetParameter(i, 2 - i % 2);
		for (size_t i = 12; i < 16; i++)
			fn->SetParameter(i, truegran1[4][i - 12]);
		for (size_t i = 16; i < 20; i++)
			fn->SetParameter(i, fitPar[4][obl][i]);
		fn->SetParameter(20, 0);
		//fn->Draw("same");
		TProfile* hm = (TProfile*)f1->Get(Form("zobl%d,sensor%d", obl + 1, 5));
		TProfile* he = (TProfile*)f2->Get(Form("zobl%d,sensor%d", obl + 1, 5));
		hm->SetLineColor(2);
		he->Draw("same");
		hm->Draw("same");
		c->Update();
		cin.get();
	}
	TLine* lineg1 = new TLine(truegran1[4][0], 0, truegran1[4][0], 5.0);
	TLine* lineg2 = new TLine(truegran1m[4][0], 0, truegran1m[4][0], 5.0);
	TLine* lineg3 = new TLine(xvt[0], 0, xvt[0], 5.0);
	TLine* lineg4 = new TLine(xvt[1], 0, xvt[1], 5.0);
	TLine* lineg5 = new TLine(xvt[2], 0, xvt[2], 5.0);
	TLine* lineg6 = new TLine(xvt[3], 0, xvt[3], 5.0);
	lineg2->SetLineColor(2);
	lineg1->Draw("same");
	lineg2->Draw("same");
	lineg3->Draw("same");
	lineg4->Draw("same");
	lineg5->Draw("same");
	lineg6->Draw("same");
	c->Update();
}

void CalibrExp::comparePhiShifts() {
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
	readTruegran1m();
	double phiBE[9];
	double phiBEE[9];
	double shiftBM[9];
	double shiftBME[9];
	double shiftBE[9];
	double shiftBEE[9];

	double phiSE[27];
	double phiSEE[27];
	double shift[27];
	double shiftE[27];
	for (size_t i = 0; i < 9; i++) {
		phiBE[i] = truegran1[i][2] * 180. / 3.14159;
		phiBEE[i] = truegran1mErr[i][2] * 180. / 3.14159 / 2.;
		if (i != 8) {
			shiftBM[i] = (truegran1m[i + 1][0] - truegran1m[i][3]) * 180. / 3.14159;
			shiftBE[i] = (truegran1[i + 1][0] - truegran1[i][3]) * 180. / 3.14159;
			shiftBME[i] = max(truegran1mErr[i + 1][0], truegran1mErr[i][3]) * 180. / 3.14159;
			shiftBEE[i] = max(truegran1Err[i + 1][0], truegran1Err[i][3]) * 180. / 3.14159;
		}
		else if (i == 8) {
			shiftBM[i] = (truegran1m[0][0] - truegran1m[i][3] ) * 180. / 3.14159 + 360;
			shiftBE[i] = (truegran1[0][0] - truegran1[i][3]) * 180. / 3.14159 + 360;
			shiftBME[i] = max(truegran1mErr[0][0], truegran1mErr[i][3]) * 180. / 3.14159;
			shiftBEE[i] = max(truegran1Err[0][0], truegran1Err[i][3]) * 180. / 3.14159;
		}

		phiSE[3*i] = truegran1[i][0] * 180. / 3.14159;
		shift[3*i] = (truegran1m[i][0] - truegran1[i][0]) * 180. / 3.14159;
		phiSEE[3*i] = truegran1Err[i][0] * 180. / 3.14159 / 2.;
		shiftE[3*i] = max(truegran1mErr[i][0], truegran1Err[i][0]) * 180. / 3.14159 / 2.;
		phiSE[3 * i+1] = truegran1[i][1] * 180. / 3.14159;
		shift[3 * i+1] = (truegran1m[i][1] - truegran1[i][1]) * 180. / 3.14159;
		phiSEE[3 * i+1] = truegran1Err[i][1] * 180. / 3.14159 / 2.;
		shiftE[3 * i+1] = max(truegran1mErr[i][1], truegran1Err[i][1]) * 180. / 3.14159 / 2.;
		phiSE[3 * i+2] = truegran1[i][3] * 180. / 3.14159;
		shift[3 * i+2] = (truegran1m[i][3] - truegran1[i][3]) * 180. / 3.14159;
		phiSEE[3 * i+2] = truegran1Err[i][3] * 180. / 3.14159 / 2.;
		shiftE[3 * i+2] = max(truegran1mErr[i][3], truegran1Err[i][3]) * 180. / 3.14159 / 2.;
	}
	TGraphErrors* gro = new TGraphErrors(27, &phiSE[0], &shift[0], &phiSEE[0], &shiftE[0]);
	TGraphErrors* grbm = new TGraphErrors(9, &phiBE[0], &shiftBM[0], &phiBEE[0], &shiftBME[0]);
	TGraphErrors* grbe = new TGraphErrors(9, &phiBE[0], &shiftBE[0], &phiBEE[0], &shiftBEE[0]);
	grbm->SetLineColor(2);
	grbe->SetLineColor(3);
	//gro->Draw("AP");
	grbm->Draw("AP");
	grbe->Draw("Psame");
}

void CalibrExp::compareEffSpectr() {
	TFile* f1 = new TFile((workingDir + "profilesEffMod.root").c_str());
	TFile* f2 = new TFile((workingDir + "profilesEff.root").c_str());
	for (size_t zobl = 3; zobl < 4; zobl++) {
		for (size_t counter = 0; counter < 9; counter++) {
			TProfile* hprof1 = (TProfile*)f1->Get(Form("mean,sensor%d", counter + 1));
			TProfile* hprof2 = (TProfile*)f2->Get(Form("mean,sensor%d", counter + 1));
			hprof1->SetLineColor(2);
			hprof1->SetLineWidth(3);
			hprof2->SetLineWidth(2);
			//hprof1->Rebin(3);
			//hprof2->Rebin(2);
			//hprof1->Divide(hprof2);
			hprof1->Draw();
			hprof1->GetYaxis()->SetRangeUser(0.7, 1.1);
			hprof2->Draw("same");
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			cin.get();
		}
	}
}

void CalibrExp::compareZAmpSpectr() {
	TFile* f1 = new TFile((workingDir + "zprofilesmod.root").c_str());
	TFile* f2 = new TFile((workingDir + "zprofiles.root").c_str());
	for (size_t counter = 0; counter < 9; counter++) {
		TProfile* hprof1 = (TProfile*)f1->Get(Form("sensor0%d", counter + 1));
		TProfile* hprof2 = (TProfile*)f2->Get(Form("sensor0%d", counter + 1));
		hprof1->SetLineColor(2);
		hprof1->Draw();
		hprof2->Draw("same");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
	}
	
}

void CalibrExp::phiShift() {
	ifstream ifile1("/work/users/kladov/snd2k/R007-002/mapdb.dat");
	double xfromdb[9][4][3];
	for(size_t j = 0; j < 9; j++){
		for(size_t i = 0; i < 4; i++){
			for(size_t k = 0; k < 3; k++){
				ifile1 >> xfromdb[j][i][k];
			}
		}
	}
	ifile1.close();

	readTrueGran1();
	readTruegran1m();
	readMeanshir();
	readMeanshirm();
	double dbord[4][9];
	double dbord1[4][9];
	double dbordErr[4][9];
	double dbord1Err[4][9];
	double modelGran[9][4] = { 3.96954, 17.4603, 30.8096, 42.9946, 43.6686, 58.0747, 70.5603, 82.8407, 83.7598, 97.8676, 110.007, 122.273, 124.125, 138.028, 150.658, 163.491, 164.016, 178.457, 192.73, 203.046, 203.326, 217.694, 231.488, 242.669, 244.352, 258.508, 271.534, 282.903, 283.639, 298.182, 310.422, 322.809, 323.23, 337.838, 351.363, 362 };
	for (size_t j = 0; j < 4; j++) {
		for (size_t i = 0; i < 9; i++) {
			//cout << truegran1[i][j] - truegran1m[i][j] << endl;
			dbord[j][i] = (truegran1[i][j]- truegran1m[i][j]) * 180. / PI;
			dbord1[j][i] = modelGran[i][j]- truegran1m[i][j] * 180. / PI;
			//dbordErr[j][i] = sqrt(pow(truegran1Err[i][j],2) + pow(truegran1mErr[i][j], 2));
			dbordErr[j][i] = 180. / PI * sqrt(pow(meanshir[i][(int)(j/1.5)]/2.,2) + pow(meanshirm[i][(int)(j/1.5)]/2., 2));
			dbord1Err[j][i] = 180. / PI * meanshirm[i][(int)(j/1.5)]/2.;
			cout << dbord[j][i] << "	" << dbordErr[j][i] << endl;
		}
		cout << endl;
	}
	
	for (size_t i = 0; i < 9; i++) {
		for (size_t j = 0; j < 4; j++) {
			cout << truegran1m[i][j]* 180. / PI << endl;
		}
		cout << endl;
	}
	double counters[9] = { 1.,2.,3.,4.,5.,6.,7.,8.,9. };
	double countersErr[9] = { 0.,0.,0.,0.,0.,0.,0.,0.,0. };
	for (size_t j = 0; j < 4; j++) {
		TGraphErrors* gr = new TGraphErrors(9, counters, dbord[j], countersErr, dbordErr[j]);
		gr->SetMarkerSize(0.9);
		gr->SetMarkerStyle(20);
		gr->SetLineColor(2);
		//gr->Draw("AP");
		TGraphErrors* gr1 = new TGraphErrors(9, counters, dbord1[j], countersErr, dbord1Err[j]);
		gr1->SetMarkerSize(0.9);
		gr1->SetMarkerStyle(20);
		gr1->SetLineColor(3);
		gr1->Draw("AP");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
	}
}

void phisravn() {
	CalibrExp a9;
	CalibrExp a7;
	a7.basedirExp = "/work/users/konctbel/calibs/R007-001/output/ntuples/MHAD2017/col/*col.root";	//directory with exp col stream files to which apply calibration
	a7.basedirMod = "/work/users/kladov/snd2k/R007-001/2017/";										//directory with mod files after map creation
	a7.workingDir = "/work/users/kladov/snd2k/R006-003/maindir/2017/";								//specify where to store temporary files
	a7.dirSelectedExp = "/work/users/kladov/snd2k/R006-003/selected2017/";							//where to store selected Baba events
	a7.dirSelectedMod = a7.workingDir + "model/selected/";											//where to store selected Baba events for modeling
	a9.readTrueGran1();
	a9.readTruegran1m();
	a9.readMeanshir();
	a9.readMeanshirm();
	a7.readTrueGran1();
	a7.readTruegran1m();
	a7.readMeanshir();
	a7.readMeanshirm();

	double dbord[4][9];
	double dbord1[4][9];
	double dbordErr[4][9];
	double dbord1Err[4][9];
	double modelGran[9][4] = { 3.96954, 17.4603, 30.8096, 42.9946, 43.6686, 58.0747, 70.5603, 82.8407, 83.7598, 97.8676, 110.007, 122.273, 124.125, 138.028, 150.658, 163.491, 164.016, 178.457, 192.73, 203.046, 203.326, 217.694, 231.488, 242.669, 244.352, 258.508, 271.534, 282.903, 283.639, 298.182, 310.422, 322.809, 323.23, 337.838, 351.363, 362 };
	for (size_t j = 0; j < 4; j++) {
		for (size_t i = 0; i < 9; i++) {
			//cout << truegran1[i][j] - truegran1m[i][j] << endl;
			dbord[j][i] = (a9.truegran1[i][j] - a7.truegran1[i][j]) * 180. / PI;
			//dbord1[j][i] = modelGran[i][j] - truegran1m[i][j] * 180. / PI;
			dbord1[j][i] = (a9.truegran1m[i][j] - a7.truegran1m[i][j]) * 180. / PI;
			//dbordErr[j][i] = sqrt(pow(truegran1Err[i][j],2) + pow(truegran1mErr[i][j], 2));
			dbordErr[j][i] = 180. / PI * sqrt(pow(a9.meanshir[i][(int)(j / 1.5)] / 2., 2) + pow(a7.meanshir[i][(int)(j / 1.5)] / 2., 2));
			dbord1Err[j][i] = 180. / PI * sqrt(pow(a9.meanshirm[i][(int)(j / 1.5)] / 2., 2) + pow(a7.meanshirm[i][(int)(j / 1.5)] / 2., 2));
			//dbord1Err[j][i] = 180. / PI * meanshirm[i][(int)(j / 1.5)] / 2.;
			cout << dbord[j][i] << "	" << dbordErr[j][i] << endl;
		}
		cout << endl;
	}

	double counters[9] = { 1.,2.,3.,4.,5.,6.,7.,8.,9. };
	double countersErr[9] = { 0.,0.,0.,0.,0.,0.,0.,0.,0. };
	for (size_t j = 0; j < 4; j++) {
		TGraphErrors* gr = new TGraphErrors(9, counters, dbord[j], countersErr, dbordErr[j]);
		gr->SetMarkerSize(0.9);
		gr->SetMarkerStyle(20);
		gr->SetLineColor(2);
		gr->Draw("AP");
		TGraphErrors* gr1 = new TGraphErrors(9, counters, dbord1[j], countersErr, dbord1Err[j]);
		gr1->SetMarkerSize(0.9);
		gr1->SetMarkerStyle(20);
		gr1->SetLineColor(3);
		gr1->Draw("same");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
	}

}

void go(){ 
	CalibrExp ab;
	int cicleV = 1;
	int chooseV = 0;
	shwidth = 1.5;
	while (cicleV == 1) {
		cout << "11 - copy; 12 - findborders; 13 - findruns; 14 - zProf(2); 15- drawZ; 16 - phiProf(2); 17 - fit; 18 - meanBorders; 19 - generateMap (191 - without shifts)" << endl;
		cout << "1 - to launch full block one (as above) - maybe I will ask someone and make it as job" << endl;
		cout << "21 - copyMod; 22 - phiProfMod (221 - with shifts); 23 - fitMod; 24 - zProfMod; 25 - drawZ; 26 - compareExp; 27 - compareMod (271 - with shifts); 28 - compare" << endl;
		cout << "2 to launch full modeling analisys block" << endl;
		cout << "31 - getTspectre; 32- draw it; 33 - one iteration of amptime; 34 (48) - compareAmpSp; 35 - EffSp; 36 - zAmpSp; 37 - phiShift, 48 - borders and profiles comparison" << endl;

		cin >> chooseV;
		cout << endl;
		switch (chooseV)
		{
		case 0:
			cicleV = 0;
			break;
		case 1: {
			shwidth = 1.5;
			ab.copy1();
			ab.findrunborders();
			ab.findruns(ab.entriesInGoodRun, ab.effInGoodRun);
			ab.zraspr();
			ab.zraspr();
			ab.linesforz("zprofiles.root", "allsensorsl", "allsensors0r");
			ab.raspr();
			ab.raspr();
			cin.get();
			ab.testforfit("profiles.root", "zobl");
			ab.bordersplot();
			ab.modelfile(true);
			break;
		}
		case 11:
			ab.copy1();
			break;
		case 12:
			ab.findrunborders();
			break;
		case 13:
			ab.findruns(ab.entriesInGoodRun, ab.effInGoodRun);
			break;
		case 14:
			ab.zraspr();
			ab.zraspr();
			break;
		case 15:
			ab.linesforz("zprofiles.root", "allsensorsl", "allsensors0r");
			break;
		case 16:
			ab.raspr();
			ab.raspr();
			break;
		case 17:
			shwidth=1.5;
			ab.testforfit("profiles.root","zobl");
			break;
		case 18:
			ab.bordersplot();
			break;
		case 19:
			ab.modelfile(true);
			break;
		case 191:
			ab.modelfile(false);
			break;
		case 2: {
			shwidth = 1.5;
			ab.copymod();
			ab.rasprmod(false);
			ab.testforfit("profilesmmod.root","mean");
			ab.zrasprmod(); 
			ab.linesforz("zprofilesmod.root","allsensors","allsensors0");
			ab.compare1();
			ab.compare2(false);
			ab.compare();
		}
		case 21:
			ab.copymod();
			break;
		case 22:
			ab.rasprmod(false);
			break;
		case 221:
			ab.rasprmod(true);
			break;
		case 23:
			shwidth = 1.5;
			ab.testforfit("profilesmmod.root","mean");
			break;
		case 24:
			ab.zrasprmod();
			break;
		case 25:
			ab.linesforz("zprofilesmod.root", "allsensors", "allsensors0");
			break;
		case 26:
			ab.compare1();
			break;
		case 27:
			ab.compare2(false);
			break;
		case 271:
			ab.compare2(true);
			break;
		case 28:
			ab.compare();
			break;

		case 31:
			ab.timespectra();
			break;
		case 32:
			ab.drawtsp();
			break;
		case 33:
			ab.ampltime();
			ab.drawamptime();
			break;
		case 34:
			ab.compareAmpSpectr();
			break;
		case 35:
			ab.compareEffSpectr();
			break;
		case 36:
			ab.compareZAmpSpectr();
			break;
		case 37:
			ab.phiShift();
			break;
		case 48:
			ab.compareAmpSpectrG();
			break;
		case 50:
			ab.comparePhiShifts();
			break;
		default:
			break;
		}
	}
//	raspramp();
//	linesa();
//	checkmod(); 
//	checkbord();
//	testt();	
//	calorimetr();
//	draweffspectre();
//	drawcalorimetr();
//	drawsomefit();
//	dopeff();
//	achspectresravn();
//	enpoint();
//	nvsphiorig();
//	drawnvsphi();
//	achspec();	
//	phiminusphis();
	
// variables that needt to be set: basedir for col stream files after "otpis`"; 

	// can be done: fit z distribution to better recalculate borders
	// think about pik amp vs zarea fitting - distribution of events inside is not plain, so maybe point should be shifted from the middle according to event amount

	//pscp D:\workProgramms\cherenkovCalibration\select_ach1_2019.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/select_ach1_2019.cpp
	//pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/select_ach1_2019.cpp D:\workProgramms\cherenkovCalibration\select_ach1_2019.cpp
}

