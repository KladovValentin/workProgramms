#include <fstream> 
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>

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

float zobl[15] = { -11.7013, -10.0367, -8.47203, -6.90737, -5.34273, -3.77807, -2.21343, -0.648766, 0.9159, 2.48053, 4.0452, 5.60983, 7.1745, 8.73913, 9.6};

float phi[2], phis[2], mcphi[2], ach[40], achr[40], theta[2], z0[2], schr[40], tch[40], tchr[40], energy[40], beam, eton, eventtime;	/// normal float event parameters
int nch, run, nc, nn, cosm, act, region[40], eventn;																					/// normal int   event parameters
float d2phi[2], dphirho[2], d2rho[2], d2z0[2], d2cosTh[2], dz0cosTh[2], Dtheta[2], Dphi[2], energyerr[2];								/// errors for   event parameters

const char *cond = "(cosm == 0) && (nc == 2) && (nn==0) && (charge[0] == 1) && (charge[1] == 1) && (energy[0]/beam>=0.8) && (energy[1]/beam>=0.7) && (nch >= 9) && (act == 0) &&\
((theta[0] + theta[1] - 3.14159) <= 0.12) && ((theta[0] + theta[1] - 3.14159) >= -0.12) && (sin(theta[0]) >= 0.41) && \
((((phi[0] - phi[1]) - 3.14159 <= 0.1) && (((phi[0] - phi[1]) - 3.14159) >= -0.1)) || (((phi[0] - phi[1]) + 3.14159 <= 0.1) && (((phi[0] - phi[1]) + 3.14159) >= -0.1))) && \
((z0[0] - z0[1]) <= 1.0) && ((z0[0] - z0[1]) >= -1.0) && ((d0[0] - d0[1]) <= 0.5) &&  ((d0[0] - d0[1]) >= -0.5)";


float pedestali1[25][100], ped[14][9], pedz[9], ped1[9][14][9], ped2[9], pede[14][9][4];
float PI = 3.14159;

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


class CalibrExp {
public:
	TString basedirExp;
	TString basedirMod;
	string workingDir;
	string dirSelectedExp;
	string prefixSelectedExp;
	string dirSelectedMod;
	string prefixSelectedMod;
	pair<float,float> timeCut;
	pair<double,double> aerBord;
	pair<double,double> aerBordM;
	double aerRLowBord;
	double aerWidth;
	double aerWidthM;
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
	double zBordMod[9][2];
	double zBordExp[9][2];

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

	double shwidth;

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

		//basedirExp = "/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2019/eecalibT/*col.root";	//directory with exp col stream files to which apply calibration
		basedirExp = "/online/users2/kladov/R007-002/output/ntuples/MHAD2019/eecalibT/*col.root";		//directory with exp col stream files to which apply calibration
		//basedirExp = "/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2019/eecalibT/*col.root";	//directory with exp col stream files to which apply calibration
		basedirMod = "/online/users2/kladov/R007-002/2019/";											//directory with mod files after map creation
		//basedirMod = "/work/users/kladov/snd2k/R007-002/2019/";										//directory with mod files after map creation

		workingDir = "/online/users2/kladov/R007-002/2019/workDir/";								//specify where to store temporary files
		//dirSelectedExp = "/work/users/kladov/snd2k/R006-003/selected2019T/";						//where to store selected Baba events
		dirSelectedExp = "/online/users2/kladov/R007-002/2019/selectedExp/";						//where to store selected Baba events
		//dirSelectedMod = "/online/users2/kladov/snd2kOld/R006-003/maindir/2017/model/selected/";	//where to store selected Baba events for modeling
		dirSelectedMod = "/online/users2/kladov/R007-002/2019/selectedMod/";						//where to store selected Baba events for modeling
		prefixSelectedExp = "true1__";
		prefixSelectedMod = "true_";
		timeCut  = make_pair(84, 994);			//(50, 113)		(2017);		(80, 1100)		(2019);		(80, 1100)		(2020)
		aerBord  = make_pair(11.7, 10.8);		//(11.7, 10.8)	(2017);		(12.34, 10.3)	(2019);		(12.05, 10.57)	(2020)
		aerBordM = make_pair(11.95, 10.45);		//(11.75,10.6)	(2017);		(11.75,10.6)	(2019);
		
		aerRLowBord = 10.5;
		aerWidth = 3.4; //3.4
		aerWidthM = 3.4; //3.4
		//
		achCut = 70.;
		entriesInGoodRun = 2000;
		effInGoodRun = 0.81;
		shwidth = 2.5;
	}

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
		fout1 << std::fixed;
		fout1 << setprecision(4);
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
		fout1 << std::fixed;
		fout1 << setprecision(4);
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
		fout1 << std::fixed;
		fout1 << setprecision(4);
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
		fout1 << std::fixed;
		fout1 << setprecision(4);
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
		fout1 << std::fixed;
		fout1 << setprecision(4);
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

	void readZBord(string fileName) {
		ifstream fin1;
		fin1.open((workingDir + fileName).c_str());
		for (int i = 0; i < 9; i++) {
			if (fileName == "zBordMod.dat")
				fin1 >> zBordMod[i][0] >> zBordMod[i][1];
			else if (fileName == "zBordExp.dat")
				fin1 >> zBordExp[i][0] >> zBordExp[i][1];
			fin1.get();
		}
		fin1.close();
	}
	void writeZBord(string fileName) {
		ofstream fout1;
		fout1.open((workingDir + fileName).c_str());
		for (int i = 0; i < 9; i++) {
			if (fileName == "zBordMod.dat") 
				fout1 << zBordMod[i][0] << "	" << zBordMod[i][1];
			else if (fileName == "zBordExp.dat")
				fout1 << zBordExp[i][0] << "	" << zBordExp[i][1];
			fout1 << endl;
		}
		fout1.close();
	}

	bool TimeIsGood(int signnumb);
	double recalculateZ(double zin, double theta, double ach, double pedz, int state);
	void copy(TString in, TString out, const char* tname);
	void copy1();
	void copymod();
	bool RunIsGood(int run);
	void findruns(int amount, double effThreshold);
	void zraspr();
	void linesforz(std::string name, const char* title1, const char* title2);
	void raspr();
	void rasprmod(bool transformPYN, bool shiftZYN);
	void testforfit(string basefile, string pName);
	void ampltime();
	void bordersplot();
	void modelfile(bool transformPYN);
	void drawamptime();
	void zrasprmod(bool shiftZYN);
	void timespectra();
	void drawtsp();
	double compare1();
	double compare2(bool transformPYN);
	void compare();
	void compareAmpSpectr();
	void compareAmpSpectrG();
	void compareAmpDistributions();
	void compareEffSpectr();
	void compareZAmpSpectr();
	void comparePhiShifts();
	vector<double> defineLGrans(int j);
	void shiftPhi2Pi();
	void getlist(int basedirInd, string hbookdir);
	double shiftZ(double z);
};


int findZOblInd(double z) {
	for (int i = 0; i < 14; i++) {
		if ((z >= zobl[i]) && (z < zobl[i + 1]))
			return i;
	}
	return -1;
}

int findMaxAmpInd(int s1, int s, float* a) {
	double maxach = -100.;
	int maxampla = s1;
	if (s1 == s)
		return 0;
	for (int f = s1+1; f <= s; f++) {
		if (a[f] > maxach) {
			maxach = a[f];
			maxampla = f;
		}
	}
	return maxampla;
}

double divisionError(double divident, double divider, double dividentErr, double dividerErr) {
	//just summing up dispersions 
	return sqrt(pow(dividentErr / divider, 2) + pow(dividerErr / divider * divident / divider, 2));
}

bool CalibrExp::TimeIsGood(int signnumb) {
	//return (((schr[signnumb] < 0.5) || ((schr[signnumb] > 0.5) && ((tchr[signnumb] > 87) && (tchr[signnumb] < 104)))) && (eventtime == 0));
	return (signnumb == 0 || (schr[signnumb] < 0.5) || ((tchr[signnumb] < timeCut.second) && (tchr[signnumb] > timeCut.first))); 
	//return (((schr[signnumb] < 0.5) || ((schr[signnumb] > 0.5) && (tchr[signnumb] > timeCut))));
	//return (tchr[signnumb] < timeCut);
}

bool CalibrExp::RunIsGood(int run) {
	bool isrungood = (run - (MinimRun + 1)) < 0 ? false : goodrun[run - MinimRun - 1];
	TString dirSelectedExpT = (TString)dirSelectedExp;
	isrungood = isrungood && (run > 27250 || !dirSelectedExpT.EndsWith("2017/"));
	return isrungood;
}

double CalibrExp::shiftZ(double z) {
	return z + (-aerBord.first - (-aerBordM.first)) + (z - (- aerBordM.first)) * (aerBord.second - aerBordM.second - (-aerBord.first) + (-aerBordM.first)) / (aerBordM.second - (-aerBordM.first));
}

double CalibrExp::recalculateZ(double zin, double theta, double ach, double ped, int state) {
	//dont take entries with ztr > borders, because it ruin statistics - make big corner influence on first zone
	double amplitude = -1.0;
	double border = aerBord.second;
	if (zin < 0)
		border = aerBord.first;

	double width = aerWidth;
	if (state == 2) {
		border = aerBordM.second;
		if (zin < 0)
			border = aerBordM.first;
	}

	if (((border - fabs(zin)) * fabs(tan(theta)) < aerWidth) && (fabs(border - fabs(zin)) >= 0.01)) {

		amplitude = (ach - ped) * aerWidth * fabs(cos(theta)) / fabs(border - fabs(zin));
	}
	else if ((border - fabs(zin)) * fabs(tan(theta)) >= aerWidth) {

		amplitude = (ach - ped) * sin(theta);
	}
	return amplitude;
}

vector<double> CalibrExp::defineLGrans(int j) {
	vector<double> lgransl(4);
	lgransl[0] = max(truegran1[j][0], truegran1m[j][0]) + 2. * max(meanshir[j][0], meanshirm[j][0]);
	lgransl[1] = min(truegran1[j][1], truegran1m[j][1]) - 4. * max(meanshir[j][1], meanshirm[j][1]);// -0.045;
	lgransl[2] = max(truegran1[j][1], truegran1m[j][1]) + 4. * max(meanshir[j][1], meanshirm[j][1]);// +0.045;
	lgransl[3] = min(truegran1[j][3], truegran1m[j][3]) - 2. * max(meanshir[j][2], meanshirm[j][2]);// -0.03;
	return lgransl;
}

bool insideLGrans(double angle, double* lgran) {
	return ((angle > lgran[0]) && (angle < lgran[1])) || ((angle > lgran[2]) && (angle < lgran[3]));
	//return ((angle > (lgran[0]+lgran[1])/2) && (angle < lgran[1])) || ((angle > lgran[2]) && (angle < (lgran[3]+lgran[2]*3)/4.));
	//return ((angle > lgran[2]) && (angle < lgran[3]));
}

void CalibrExp::shiftPhi2Pi() {
	if (phi[0] < 1.)
		phi[0] = phi[0] + 2. * PI;
	if (phi[1] < 1.)
		phi[1] = phi[1] + 2. * PI;
}

void CalibrExp::copy(TString in, TString out, const char *tname) {
	/// copy only e+e- entries (cond above) from col stream, located in "basedir", to a new directory "dir"
	vector<TString> ans;
	TString files = gSystem->GetFromPipe("ls " + in); ///base directory for files

	Ssiz_t from = 0;
	TString elem;

	while (files.Tokenize(elem, from, "\n")) {
		ans.push_back(elem);
	}
	size_t vector_size = ans.size();
	for (size_t i = 0; i < vector_size; i++) {
		TString str = ans[i].Copy();
		TString tok;
		Ssiz_t from1 = 1;
		vector<TString> tokens;
		TFile* filecomb = TFile::Open(str);
		TTree* originalTree = (TTree*)filecomb->Get(tname);
		while (str.Tokenize(tok, from1, "/")) {
			tokens.push_back(tok);
		}
		TString fname = tokens.back();
		cout << out + fname << endl;
		TFile* ouput = TFile::Open(out + fname, "RECREATE");
		ouput->cd();
		TTree* selectedTree = originalTree->CopyTree(cond);
		selectedTree->Print();
		selectedTree->Write();
		ouput->Close();
		filecomb->Close();
	}
}

void CalibrExp::copy1() {
	copy(basedirExp, dirSelectedExp + prefixSelectedExp, "t1");
}

void CalibrExp::copymod() {
	copy(basedirMod + "*0.root", dirSelectedMod + prefixSelectedMod, "h1");
}



void CalibrExp::findruns(int amount, double effThreshold) {
	//find runs that contain (events > "amount") and (mean efficiency > "effThreshold"), conditions on phi(+shir) and z
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
	MinimRun = 100000;
	MaximRun = 0;
	int scount = 0;
	int scount1 = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		histo->Fill(run);
		if (run < MinimRun)
			MinimRun = run;
		if (run > MaximRun)
			MaximRun = run;
		schr[nch] = 0.;	// to know that nch-1 is the last for 9's counter
		
		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			//wichcounter1[j] = scount1;
			if (j == 8) { shiftPhi2Pi(); }
			while (schr[1 + scount] > 0.5) { scount += 1; }

			vector<double> lgran = defineLGrans(j);

			for (int prtc = 0; prtc < 2; prtc++) {
				double ztr = 12. / tan(theta[prtc]) + z0[prtc];
				if (insideLGrans(phi[prtc], &lgran[0]) && (ztr > zobl[0]) && (ztr < zobl[14]))
					hprof[j]->Fill(run, schr[1 + scount1]);
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

	TFile* MyFile = new TFile((workingDir + "runsi.root").c_str(), "RECREATE");
	for(int j = 0; j<9;j++){
		sprintf(title, "sens%d_eff_run", j + 1);
		hprof[j]->Write(title);
	}
	histo->Write("runsnumb");
	MyFile->Close();

	MinimRun = MinimRun - 1;
	MaximRun = MaximRun + 1;
	cout << "run range:  " << MinimRun << " - " << MaximRun << endl;
	writeMinMaxRun();

	readMinMaxRun();
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

void CalibrExp::zraspr() {
	//make profile with mean ach vs z on the inner cillinder and medium cillinder with recalculation on borders according to intersection line segment lenght and without
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
	float ztr[2] = {0,0};
	float zin[2] = {0,0};
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		schr[nch] = 0.;
		ach[0] = 0;
		if ((!RunIsGood(run)) || (eventtime != 0)) {
			continue;
		}
		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			if (j == 8) { shiftPhi2Pi(); }
			while (schr[1 + scount] > 0.5) { scount += 1; }

			int maxampla = findMaxAmpInd(scount1, scount, &ach[0]);
			vector<double> lgran = defineLGrans(j);
				
			for (int prtc = 0; prtc < 2; prtc++) {
				ztr[prtc] = 12.0 / tan(theta[prtc]) + z0[prtc];
				zin[prtc] = aerRLowBord / tan(theta[prtc]) + z0[prtc];

				if (!insideLGrans(phi[prtc], &lgran[0]) || !TimeIsGood(maxampla))
					continue;
					
				//for pedestals
				h[j]->Fill(ach[scount1]);

				//for amplitude vs z
				//fill 0 hists with standart sin
				hprof01[prtc]->Fill(zin[prtc], (ach[maxampla]) * sin(theta[prtc]));
				hprof0[j]->Fill(zin[prtc], (ach[maxampla]) * sin(theta[prtc]));

				//fill with recalculated
				double fillAmplitude = recalculateZ(zin[prtc], theta[prtc], ach[maxampla], 0, 1);
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
	
	TF1* fn = new TF1("scphi", "[0]*(x-[1])", -12, -7);
	fn->SetParameters(hprof01[1]->GetMaximum() / aerWidth, -aerBord.first);
	hprof01[1]->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	hprof01[1]->Fit("scphi", "", "", -aerBord.first + 0.2, -aerBord.first + aerWidth);
	double newBord = fn->GetParameter(1);
	hprof01[1]->Fit("scphi", "", "", newBord + aerWidth/4., newBord + aerWidth);
	cout << fn->GetParameter(1) << endl;
	c->Update();
	cin.get();
	cin.get();


	TF1* fn1 = new TF1("scphi1", "[0]*([1]-x)+[2]*exp(-((x-[3])/[4])**2)", 8, 12);
	fn1->SetParameters(10, aerBord.second, 1, 10.5, 1);
	fn1->SetParLimits(2,0,0.5* hprof01[1]->GetMaximum());
	fn1->SetParLimits(3, aerBord.second - aerWidth, aerBord.second);
	fn1->SetParLimits(4, 0, 2*aerWidth);
	hprof01[1]->Fit("scphi1", "", "", 8.5, 10.8);
	hprof01[1]->Fit("scphi1", "", "", fn1->GetParameter(1) - aerWidth / 1.5, fn1->GetParameter(1));
	cout << fn1->GetParameter(1) << endl;


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

void CalibrExp::zrasprmod(bool shiftZYN) {
	//same as zraspr but without conditions on run and nch is always 9, need to rewrite or, better, unite it with zraspr
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
		sprintf(title, "sensor%d;z;ach", i + 1);
		hprof[i] = new TProfile(name, title, 300, -15, 15);
	}
	TProfile* hprof1 = new TProfile("hpprof", "all;ztr;ach", 300, -15, 15);

	TProfile* hprof0[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprof0%d", i + 1);
		sprintf(title, "sensor0%d;z;ach", i + 1);
		hprof0[i] = new TProfile(name, title, 300, -15, 15);
	}
	TProfile* hprof01 = new TProfile("hpprof0", "all;ztr;ach", 300, -15, 15);

	TH1* hr = new TH1F("histr", "recount right;z,cm", 100, 5., 15.);
	TH1* hl = new TH1F("histl", "recount left;z,cm", 100, -16., -6.);

	int Count = 0, count1 = 0;
	float ztr[2] = { 0.,0. };
	float zin[2] = { 0.,0. };
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			if (j == 8) { shiftPhi2Pi(); }
			vector<double> lgran = defineLGrans(j);

			for (int prtc = 0; prtc < 2; prtc++) {
				ztr[prtc] = 12.0 / tan(theta[prtc]) + z0[prtc];
				//double radia = aerRLowBord * (1 + 4.14159 / 180. * 0.370066 * cos(phi[prtc] - 6.19112 * 3.14159 / 180.) + 3.14159 / 180. * 0.270066 * sin(phi[prtc] - 6.19112 * 3.14159 / 180.));
				double radia = aerRLowBord * (1 + 3.14159 / 180. * 0.370066 * cos(phi[prtc] - 6.19112 * 3.14159 / 180.));
				zin[prtc] = aerRLowBord / tan(theta[prtc]) + z0[prtc];
				if (shiftZYN) {
					ztr[prtc] = shiftZ(ztr[prtc]);
					zin[prtc] = shiftZ(zin[prtc]);
				}

				double fillPhiWhat = phi[prtc];
				double newPhi = 0;
				if (phi[prtc] < truegran1m[j][1])
					newPhi = truegran1[j][0] + (phi[prtc] - truegran1m[j][0]) * (truegran1[j][1] - truegran1[j][0]) / (truegran1m[j][1] - truegran1m[j][0]);
				else if (phi[prtc] >= truegran1m[j][1])
					newPhi = truegran1[j][1] + (phi[prtc] - truegran1m[j][1]) * (truegran1[j][3] - truegran1[j][1]) / (truegran1m[j][3] - truegran1m[j][1]);
				//fillPhiWhat = newPhi;
				//fillPhiWhat = fillPhiWhat - 3.14159 / 180. * (0.370066 * cos(fillPhiWhat - 6.19112 * 3.14159 / 180.) - 0.403304);

				if (!insideLGrans(fillPhiWhat, &lgran[0]))
					continue;

				//fill 0 hists with standart sin
				hprof01->Fill(zin[prtc], (ach[j]) * sin(theta[prtc]));
				hprof0[j]->Fill(zin[prtc], (ach[j]) * sin(theta[prtc]));

				//fill with recalculated
				double fillAmplitude = recalculateZ(zin[prtc], theta[prtc], ach[j], 0, 2-(int)shiftZYN);
				if (fillAmplitude == -1.0 || ach[j] > achCut)
					continue;
				hprof[j]->Fill(zin[prtc], fillAmplitude);
				hprof1->Fill(zin[prtc], fillAmplitude);

			}
			scount += 1;
			scount1 = scount;
		}
		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("obrabotano %d M entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
	}

	//TF1* fn = new TF1("scphi", "[0]*(1. / (exp( - (x-[1])/[2] + 3) + 1 ) )", -12, -7);
	TF1* fn = new TF1("scphi", "[0]*(x-[1])", -12, -7);
	for (size_t i = 0; i < 9; i++) {
		fn->SetParameters(hprof01->GetMaximum() / aerWidth, -aerBord.first);
		hprof0[i]->Draw();
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		hprof0[i]->Fit("scphi", "", "", -aerBord.first + 0.2, -aerBord.first + aerWidth);
		double newBord = fn->GetParameter(1);
		hprof0[i]->Fit("scphi", "", "", newBord + aerWidth / 4., newBord + aerWidth);
		cout << fn->GetParameter(1) << endl;
		c->Update();
		//cin.get();
	}
	//cin.get();
	
	//TF1* fn1 = new TF1("scphi1", "([0]+[1]*x)+[2]*exp(-((x-[3])/[4])**2)", 8, 12);
	TF1* fn1 = new TF1("scphi1", "[0]*([1]-x)", 8, 12);
	fn1->SetParameters(hprof01->GetMaximum() / aerWidth, aerBord.second);
	hprof01->Fit("scphi1", "", "", aerBord.second - aerWidth / 2., aerBord.second);
	hprof01->Fit("scphi1", "", "", fn1->GetParameter(1) - aerWidth / 2., fn1->GetParameter(1));
	cout << fn1->GetParameter(1) << endl;
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	//cin.get();

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

void CalibrExp::linesforz(std::string name, const char* title1, const char* title2) {
	//draw z profiles from one file "name" with titles 1 and 2
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
	for (int i = 0; i < 15; i++) {
		TLine* line = new TLine(zobl[i], 0, zobl[i], 1.05 * hprof->GetMaximum());
		line->Draw();
	}
	c->Update();
	double lb = name == "zprofilesmod.root" ? -aerBordM.first : -aerBord.first;
	double rb = name == "zprofilesmod.root" ? aerBordM.second : aerBord.second;
	TLine* line11 = new TLine(lb, 0, lb, 7);
	TLine* line12 = new TLine(rb, 0, rb, 7);
	line11->SetLineColor(kGreen);
	line12->SetLineColor(kGreen);
	line11->Draw();
	line12->Draw();
	//hr->DrawNormalized("same", 50);
	//hl->DrawNormalized("same", 50);
	hprof0->Draw("same");
	c->Update();
}

void CalibrExp::raspr() {
	/// cleaned, maybe get parameters for recalculation on borders from zraspr fit 
	/// ~400 for fit, normal sim statistic - ~100 bins
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
	readTruegran1m();
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

	//amp vs phi, raspr po phi and amp distributions
	TProfile* hprof[14][9];
	TProfile* hprofeff[14][9];
	TH1* hdist[14][9];
	TH1* hza[14][9];
	char name[20];
	char title[100];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {	
			sprintf(name, "hprofrp%d", 9*j + i + 1);
			sprintf(title, "zobl%d,sensor%d;#phi;ach", j + 1,i + 1);
			hprof[j][i] = new TProfile(name,title, 800, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);

			sprintf(name, "hprofeff%d", 9 * j + i + 1);
			sprintf(title, "zobl%d,sensor%d;#phi;eff", j + 1, i + 1);
			hprofeff[j][i] = new TProfile(name, title, 200, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);

			sprintf(name, "hdist%d", 9 * j + i + 1);
			sprintf(title, "zobl%d,sensor%d;ach", j + 1, i + 1);
			hdist[j][i] = new TH1F(name, title, 250, -5, 45);

			sprintf(name, "hzas%d", 9 * j + i + 1);
			sprintf(title, "zas,zobl%d,sensor%d;#phi", j + 1, i + 1);
			hza[j][i] = new TH1F(name, title, 800, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
		}
	}

	TProfile* hprofMean[9];
	TProfile* hprofEffMean[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprofrpZm%d", i + 1);
		sprintf(title, "sensor%d;#phi;ach",i + 1);
		hprofMean[i] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);

		sprintf(name, "hprofeffZm%d", i + 1);
		sprintf(title, "sensor%d;#phi;ach", i + 1);
		hprofEffMean[i] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
	}

	TH1* htheta0 = new TH1F("theta0","theta;#theta",100,0,180);
	TH1* htheta = new TH1F("theta","theta;#theta",100,0,180);

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
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		Count += 1;
		if (Count == 1000000) {
			count1 += 1;
			cout << Form("obrabotano %d M entries", count1) << endl;
			//cout << Form("obrabotano %d*10k entries	or %d%", count1, (count1*1000000)/entries) << endl;
			Count = 0;
		}
		//cout << beam << endl;
		schr[nch] = 0.;
		ach[0] = 0;
		if ((!RunIsGood(run) || eventtime != 0)) {
			continue;
		}
		scount = 0;   
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			//wichcounter1[j] = scount1;
			//counters are shifted in oreder to determine first one and direction of count, I think. Because of this 9th counter have the tail of signal around 0, I transfer it on >2pi
			//nch is >= 9 in 2017+, so for one counter I choose from scount1 to scount 2 the highest amplitude 
			if (j == 8) { shiftPhi2Pi(); }
			while (schr[1 + scount] > 0.5) { scount += 1; }
			
			int maxampla = findMaxAmpInd(scount1, scount, &ach[0]);
			vector<double> lgran = defineLGrans(j);
			
			for (size_t pind = 0; pind < 2; pind++) {
				ztr[pind] = 12.0 / tan(theta[pind]) + z0[pind];
				zin[pind] = aerRLowBord / tan(theta[pind]) + z0[pind];
				int zoblInd = findZOblInd(ztr[pind]);

				if ((phi[pind] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[pind] < (float)(j + 2) * (2. * PI) / 9.) && (TimeIsGood(maxampla)) &&/* (TimeIsGood(scount1)) && */ zoblInd != -1) {
				//if (insideLGrans(phi[pind], &lgran[0]) && (TimeIsGood(maxampla))/* && (TimeIsGood(scount1)) */ && zoblInd != -1) {
					//zasel po phi
					hza[zoblInd][j]->Fill(phi[pind]);

					htheta->Fill(theta[0]*180./3.14159);
					htheta->Fill(theta[1]*180./3.14159);
					if (ach[maxampla] < 0.2) {
						htheta0->Fill(theta[0]*180./3.14159);
						htheta0->Fill(theta[1]*180./3.14159);
					}

					double fillAmplitude = recalculateZ(zin[pind], theta[pind], ach[maxampla], 0, 1);
					if (fillAmplitude != -1.0 && ach[maxampla] < achCut) {
						if (insideLGrans(phi[pind], &lgran[0]))
							hdist[zoblInd][j]->Fill(fillAmplitude);
						hprof[zoblInd][j]->Fill(phi[pind], fillAmplitude);
						if (zoblInd > 2 && zoblInd < 12)
							hprofMean[j]->Fill(phi[pind], fillAmplitude);
					}

					//efficiency
					double effWhatFill = 0.;
					if (ach[maxampla] >= 0.2)
					//if (schr[scount1+1] >= 1)
						effWhatFill = 1;

					hprofeff[zoblInd][j]->Fill(phi[pind], effWhatFill);
					if (zoblInd >= 2 && zoblInd <= 12)
						hprofEffMean[j]->Fill(phi[pind], effWhatFill);

					//pedestals
					//if (insideLGrans(phi[pind], &lgran[0]))
					h[zoblInd][j]->Fill(ach[scount1]);
				}

			}

			scount += 1;
			scount1 = scount;
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

	TFile* MyFileD = new TFile((workingDir + "distributions.root").c_str(), "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "zobl%d,sensor%d", m + 1, i + 1);
			hdist[m][i]->Write(title);
		}
	}
	MyFileD->Close();

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
		htheta0->Write("theta0");
		htheta->Write("theta");
	}
	MyFileEff->Close();

	TFile* MyFileZasel = new TFile((workingDir + "zaselphi.root").c_str(), "RECREATE");
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

double hevisideFunction(double x) {
	return x > 0 ? 1 : 0;
}
vector< vector<double> > transformProfileByPhiDistr(TProfile* profile, TH1* phiDistr) {
	vector<double> xBin, entriesInBin, amplInBin, amplInBinCopy, xErr, amplError, nolErr;
	vector< pair<double, double> > minmaxPairs;
	vector< vector<double> > extremums;
	for (size_t i = 0; i < (size_t)phiDistr->GetNbinsX(); i++) {
		if (phiDistr->GetBinContent(i) < 3)
			continue;
		xBin.push_back((double)phiDistr->GetBinCenter(i));
		xErr.push_back((double)phiDistr->GetBinWidth(i) / 2.);
		entriesInBin.push_back((double)phiDistr->GetBinContent(i));
		amplInBin.push_back((double)profile->GetBinContent(i));
		amplInBinCopy.push_back((double)profile->GetBinContent(i));
		amplError.push_back((double)profile->GetBinError(i));

		nolErr.push_back(0.0);
	}
	double meanEntriesAmountBin = (double)phiDistr->GetEntries() / (double)xBin.size();
	double binWidth = phiDistr->GetBinWidth(250);


	bool makingPair = false;
	bool makingPairLow = false;
	double locMaximum = meanEntriesAmountBin;
	double locMinimum = meanEntriesAmountBin;
	int locMaxBin = 0;
	int locMinBin = 0;
	int startingBinPair = 0;
	/*for (size_t i = 0; i < xBin.size(); i++) {
		if (entriesInBin[i] > meanEntriesAmountBin * 1.15 && !makingPair) {
			makingPair = true;
			startingBinPair = i;
		}
		if (entriesInBin[i] < meanEntriesAmountBin / 1.2 && makingPair)
			makingPairLow = true;
		if ((entriesInBin[i] > meanEntriesAmountBin / 1.2 && makingPairLow))
			minmaxPairs.push_back(make_pair(xBin[locMaxBin], xBin[locMinBin]));
		if ((entriesInBin[i] > meanEntriesAmountBin / 1.2 && makingPairLow) || ((int)i-startingBinPair)>20 && makingPair) {
			cout << i << endl;
			makingPair = false;
			makingPairLow = false;
			locMaximum = meanEntriesAmountBin;
			locMinimum = meanEntriesAmountBin;
		}

		if (!makingPair)
			continue;

		if (entriesInBin[i] > locMaximum) { locMaximum = entriesInBin[i]; locMaxBin = int(i); }
		if (entriesInBin[i] < locMinimum) { locMinimum = entriesInBin[i]; locMinBin = int(i); }
	}*/

	//TCanvas* c = new TCanvas("c", "c", 1200, 800);
	//c->SetFillColor(0);
	//c->SetBorderMode(0);
	//TPad* pad1 = new TPad("pad1", "", 0, 0, 1, 1, 0, 4, 0);
	//phiDistr->Draw();
	//cin.get();
	//cin.get();
	for (size_t i = 0; i < xBin.size(); i++) {
		if (entriesInBin[i] < meanEntriesAmountBin / 1.3) {
			TF1* fxc = new TF1("fxc", "[9] + [0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[8])", xBin[i] - 0.04, xBin[i] + 0.06);
			fxc->SetParameters(meanEntriesAmountBin*0.3, xBin[i] - 0.02,0.02, -meanEntriesAmountBin * 0.3, xBin[i]+0.005, 0.02, meanEntriesAmountBin * 0.3, xBin[i] + 0.03, 0.02, meanEntriesAmountBin);
			fxc->SetParLimits(0, meanEntriesAmountBin*0.05, meanEntriesAmountBin * 0.5);
			fxc->SetParLimits(3, -meanEntriesAmountBin*0.5, -meanEntriesAmountBin * 0.05);
			fxc->SetParLimits(6, meanEntriesAmountBin*0.01, meanEntriesAmountBin * 0.5);
			fxc->SetParLimits(1, xBin[i] - 0.04, xBin[i]);
			fxc->SetParLimits(4, xBin[i] - 0.04, xBin[i] + 0.06);
			fxc->SetParLimits(7, xBin[i], xBin[i] + 0.06);

			if (fabs((xBin[i] * 180. / 3.14159) - (int)(xBin[i] * 180. / 3.14159)) > 0.15 && fabs((xBin[i] * 180. / 3.14159) - (int)(xBin[i] * 180. / 3.14159)) < 0.35)
				continue;

			//if (fabs((xBin[i] * 180. / 3.14159) - (int)(xBin[i] * 180. / 3.14159)) > 0.35)
			//	fxc->FixParameter(6, 0);
			phiDistr->Fit(fxc, "", "", xBin[i] - 0.04, xBin[i] + 0.06);
			//c->Update();
			//cin.get();
			vector<double> arr;
			arr.push_back(fxc->GetParameter(1));
			arr.push_back(fxc->GetParameter(4));
			arr.push_back(fxc->GetParameter(7));
			//if (fabs((arr[1] * 180. / 3.14159) - (int)(arr[1] * 180. / 3.14159)) > 0.35)
			//	arr[2] = arr[1];
			if (fxc->GetParameter(6)< meanEntriesAmountBin * 0.05)
				arr[2] = 0;
			extremums.push_back(arr);
			i = i + (size_t)(0.06 / binWidth);
		}
	}


	double maxNdiff = meanEntriesAmountBin;
	double prevMax = 0.;
	for (size_t k = 0; k < 1; k++) {
		for (size_t i = 0; i < xBin.size(); i++) {
			for (size_t j = 0; j < xBin.size(); j++) {
				if (fabs(entriesInBin[j] - entriesInBin[i]) > prevMax) { prevMax = fabs(entriesInBin[j] - entriesInBin[i]); }
			}
		}
		maxNdiff = prevMax;
		prevMax = 0.;
		for (size_t i = 0; i < xBin.size(); i++) {
			if ((amplInBin[i] < 0.5))
				continue;
			for (int sign = -1; sign < 0; sign = sign + 2) {
				vector<double> deltaEntr;
				vector<double> deltaAmp;
				for (size_t j = 0; j < xBin.size(); j++) {
					deltaEntr.push_back(0.);
					if (fabs(xBin[j] - xBin[i]) < -0.05 && sign*(i-j)>0)
						deltaEntr[j] = (entriesInBin[i] - entriesInBin[j]) * fabs(pow((entriesInBin[i] - entriesInBin[j]) / maxNdiff, 5)) * TMath::Gaus(xBin[i]-sign*0.02, xBin[j], 0.015) * hevisideFunction(entriesInBin[i] - meanEntriesAmountBin) * hevisideFunction(meanEntriesAmountBin - entriesInBin[j]) / 2.;

					deltaAmp.push_back((amplInBin[i] > 0.5) ? deltaEntr[j] * (amplInBin[i] - amplInBin[j]) / (entriesInBin[i] - deltaEntr[j]) : 0.);
				}
				for (size_t j = 0; j < xBin.size(); j++) {
					entriesInBin[j] += deltaEntr[j];
					entriesInBin[i] -= deltaEntr[j];
					amplInBin[i] += deltaAmp[j];
				}
			}
		}
		/*TGraph* grDistr = new TGraph(xBin.size(), &xBin[0], &entriesInBin[0]);
		phiDistr->Draw();
		grDistr->SetMarkerSize(0.5);
		grDistr->SetMarkerStyle(20);
		grDistr->Draw("Psame");
		TCanvas* c1 = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c1->SetFillColor(0);
		c1->SetFrameFillColor(0);
		c1->Update();
		cin.get();*/
		if (k%1==0)
			cout << "next round " << k << endl;
	}
	TGraph* grDistr = new TGraph(xBin.size(), &xBin[0], &entriesInBin[0]);
	//for (size_t i = 0; i < xBin.size(); i++)
	//	cout << entriesInBin[i] << endl;;
	/*grDistr->SetMarkerSize(0.5);
	grDistr->SetMarkerStyle(20);
	//grDistr->Draw("Psame");
	//c->SetFillColor(0);
	//c->SetFrameFillColor(0);
	for (size_t i = 0; i < extremums.size(); i++) {
		vector<double> arr = extremums[i];
		cout << arr[0] << "	" << arr[1] << "	" << arr[2] << endl;
		TLine* lh = new TLine(arr[0], 0, arr[0], meanEntriesAmountBin*1.3);
		TLine* ll = new TLine(arr[1], 0, arr[1], meanEntriesAmountBin*1.3);
		TLine* ls = new TLine(arr[2], 0, arr[2], meanEntriesAmountBin*1.3);
		lh->SetLineColor(kBlack);
		ll->SetLineColor(kBlack);
		ls->SetLineColor(kBlack);
		lh->Draw("same");
		ll->Draw("same");
		ls->Draw("same");
		//c->Update();
	}
	c->Update();
	cin.get();
	cin.get();*/
	TGraphErrors* grAmp = new TGraphErrors(xBin.size(), &xBin[0], &amplInBin[0], &xErr[0], &amplError[0]);
	//profile->Draw();
	//grAmp->Draw("Psame");
	//c->Update();
	//return grAmp;
	//c->Delete();
	return extremums;
}
void checkTrans() {
	TFile* fA = new TFile("/online/users2/kladov/R007-002/2019/workDir/profiles.root");
	TFile* fD = new TFile("/online/users2/kladov/R007-002/2019/workDir/zaselphi.root");
	TProfile* h_amp_phi = (TProfile*)fA->Get(Form("zobl%d,sensor%d", 2, 9));
	TH1* h_phi_distr = (TH1F*)fD->Get(Form("zobl%d,sensor%d", 2, 9));
	//TGraphErrors* gr_amp_phi = transformProfileByPhiDistr(h_amp_phi, h_phi_distr);
	vector< vector<double> > gr_amp_phi = transformProfileByPhiDistr(h_amp_phi, h_phi_distr);
}

void CalibrExp::rasprmod(bool transformPYN, bool shiftZYN) {
	//in modeling time and runs are always good, but phi distribution can be shifted, 
	//because of little statistics I only make 9 distributions for counters with mean by z amplitude in z[1]<z<z[13] to determine borders in fit
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
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
	TH1* hdist[14][9];
	char name[20];
	char title[100];
	for (int m = 0; m < 14; m++) {
		for (int j = 0; j < 9; j++) {
			sprintf(name, "hprofmod%d", m*9 + j + 1);
			sprintf(title, "zobl%d,counter%d;phi;ach", m + 1, j + 1);
			hprof[m][j] = new TProfile(name, title, 200, (float)(j - 1.) * (2. * PI) / 9., (float)(j + 2.) * (2. * PI) / 9.);
			sprintf(name, "hdist%d", 9 * m + j + 1);
			sprintf(title, "zobl%d,sensor%d;ach", m + 1, j + 1);
			hdist[m][j] = new TH1F(name, title, 250, -5, 45);
		}
	}
	TProfile* hprofMean[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprofrpZm%d", i + 1);
		sprintf(title, "sensor%d;#phi;ach", i + 1);
		hprofMean[i] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
	}
	//eff vs phi
	TProfile* hprofeff[14][9];
	for (int j = 0; j < 14; j++) {
		for (int i = 0; i < 9; i++) {
			sprintf(name, "hprofeffmod%d", 9 * j + i + 1);
			sprintf(title, "zobl%d,sensor%d;phi;eff", j + 1, i + 1);
			hprofeff[j][i] = new TProfile(name, title, 200, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
		}
	}
	TProfile* hprofEffMean[9];
	for (int i = 0; i < 9; i++) {
		sprintf(name, "hprofeffZm%d", i + 1);
		sprintf(title, "sensor%d;phi;ach", i + 1);
		hprofEffMean[i] = new TProfile(name, title, 400, (float)(i - 1) * (2. * PI) / 9., (float)(i + 2) * (2. * PI) / 9.);
	}

	TH1* htheta0 = new TH1F("theta0", "theta;#theta", 100, 0, 180);
	TH1* htheta = new TH1F("theta", "theta;#theta", 100, 0, 180);

	int Count = 0, count1 = 0;
	float ztr[2] = { 0.,0. };
	float zin[2] = { 0.,0. };
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//cout << beam << endl;
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
			if (j == 8) { shiftPhi2Pi(); }

			vector<double> lgran = defineLGrans(j);

			for (size_t pind = 0; pind < 2; pind++) {
				ztr[pind] = 12.0 / tan(theta[pind]) + z0[pind];
				zin[pind] = aerRLowBord / tan(theta[pind]) + z0[pind];
				if (shiftZYN) {
					ztr[pind] = shiftZ(ztr[pind]);
					zin[pind] = shiftZ(zin[pind]);
				}
				int zoblInd = findZOblInd(ztr[pind]);
				
				if ((phi[pind] >= (float)(j - 1) * (2. * PI) / 9.) && (phi[pind] < (float)(j + 2) * (2. * PI) / 9.) && zoblInd != -1) {
				//if (insideLGrans(phi[pind], &lgran[0]) && zoblInd != -1) {
					//hprof[zoblInd][j]->Fill(phi[pind], ach[j] * sin(theta[pind]));
					// 
					//transform angle
					htheta->Fill(theta[0]*180./3.14159);
					htheta->Fill(theta[1]*180./3.14159);
					if (ach[j] < 0.2) {
						htheta0->Fill(theta[0]*180./3.14159);
						htheta0->Fill(theta[1]*180./3.14159);
					}

					double fillPhiWhat = phi[pind];
					double newPhi = 0;
					if (phi[pind] < truegran1m[j][1])
						newPhi = truegran1[j][0] + (phi[pind] - truegran1m[j][0]) * (truegran1[j][1] - truegran1[j][0]) / (truegran1m[j][1] - truegran1m[j][0]);
					else if (phi[pind] >= truegran1m[j][1])
						newPhi = truegran1[j][1] + (phi[pind] - truegran1m[j][1]) * (truegran1[j][3] - truegran1[j][1]) / (truegran1m[j][3] - truegran1m[j][1]);
					if (transformPYN)
						//fillPhiWhat = newPhi;
						fillPhiWhat = fillPhiWhat - 3.14159/180. * (0.370066*cos( fillPhiWhat - 6.19112*3.14159/180. ) - 0.403304);
					//fillPhiWhat = phi[pind];

					double fillAmplitude = recalculateZ(zin[pind], theta[pind], ach[j], 0, 2-(int)shiftZYN);
					if (fillAmplitude != -1.0 && ach[j] < achCut) {
						if (insideLGrans(phi[pind], &lgran[0]))
							hdist[zoblInd][j]->Fill(fillAmplitude);
						hprof[zoblInd][j]->Fill(fillPhiWhat, fillAmplitude);
						if (zoblInd > 2 && zoblInd < 12) {
							hprofMean[j]->Fill(fillPhiWhat, fillAmplitude);
						}
					}

					//eff
					double effWhatFill = 0.;
					if (ach[j] >= 0.2)
					//if (tch[j] < 900)
						effWhatFill = 1.;
					hprofeff[zoblInd][j]->Fill(fillPhiWhat, effWhatFill);
					if (zoblInd >= 2 && zoblInd <= 12)
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

	TFile* MyFileD = new TFile((workingDir + "distributionsmod.root").c_str(), "RECREATE");
	for (int m = 0; m < 14; m++) {
		for (Int_t i = 0; i < 9; i++) {
			sprintf(title, "zobl%d,sensor%d", m + 1, i + 1);
			hdist[m][i]->Write(title);
		}
	}
	MyFileD->Close();

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
		htheta0->Write("theta0");
		htheta->Write("theta");
	}
	MyFileEff->Close();
}

void CalibrExp::ampltime() {
	/// get profiles of eff and ach vs run, conditions as in raspr except "insideGrans" - region condition
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
			ha[i] = new TProfile(name,title, (MaximRun - MinimRun) / 10, MinimRun, MaximRun); 
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
	int Count = 0, count1 = 0;
	float ztr[2] = { 0., 0.};
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);

		if (fabs(beam - beamlast) > 1.7) {  //when energy change
			rench.push_back((float)run);
			eench.push_back(beam);
		}
		beamlast = beam;

		schr[nch] = 0.;
		ach[0] = 0;
		if ((!RunIsGood(run))) {
			continue;
		}
			
		scount = 0;
		scount1 = 0;

		for (int j = 0; j < 9; j++) {
			if (j == 8) { shiftPhi2Pi(); }
			while (schr[1 + scount] > 0.5) { scount += 1; }

			int maxampla = findMaxAmpInd(scount1, scount, &ach[0]);
			vector<double> lgran = defineLGrans(j);

			for (size_t pind = 0; pind < 2; pind++) {
				ztr[pind] = 12.0 / tan(theta[pind]) + z0[pind];
				int zoblInd = findZOblInd(ztr[pind]);
				if (insideLGrans(phi[pind], &lgran[0])  && (TimeIsGood(maxampla)) && zoblInd >= 2 && zoblInd <= 11 && ach[maxampla] < achCut ) {
					
					double ach1 = ach[maxampla]*sin(theta[pind]);
					bool eff = (bool)schr[scount1 + 1] && ach1>0.2;
					he[0]->Fill(run, (int)eff);
					he[j + 1]->Fill(run, (int)eff);

					ha[0]->Fill(run, ach1);
					ha[j + 1]->Fill(run, ach1);
					ampCounter[j]->Fill(ach1);
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

	TFile* MyFile2 = new TFile((workingDir + "achSpectCount.root").c_str(), "RECREATE");
	for (int j = 0; j < 9; j++) {
		sprintf(title, "aspecs%d", j + 1);
		ampCounter[j]->Write(title);
	}
	MyFile2->Close();

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

void CalibrExp::drawamptime() {
	/// derive new bin sequence from mean ach vs run, draw, write coeffs in .dat and profiles in .root
	readMinMaxRun();
	double re, ee;
	vector<float> rench;
	vector<float> eench;
	ifstream fin;
	fin.open((workingDir + "ench.cpp").c_str());
	{
		int j = 0;
		while (fin >> re >> ee) {
			rench.push_back(re);
			cout << re << "	";
			eench.push_back(ee / 100);
			cout << ee / 100 << endl;
			j++;
		}
	}
	fin.close();

	TFile* MyFile = new TFile((workingDir + "achvsrun.root").c_str());       // profiles vs run
	TFile* MyFile2 = new TFile((workingDir + "achSpectCount.root").c_str()); // hists for mean

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
	for (size_t j = 2; j < (size_t)Nbins + 1; j++) {
		//cout << ha[0]->GetBinCenter(j) << "	" << ha[0]->GetBinWidth(j);
		foutra << ha[0]->GetBinCenter(j) << "	" << ha[0]->GetBinWidth(j);
		for (int k = 0; k < 9; k++) {
			foutra << "	" << ha[k]->GetBinContent(j) / meanAmp[k];
		}
		foutra << endl;
	}
	foutra.close();

	TFile* MyFile1 = new TFile((workingDir + "amptime.root").c_str(), "RECREATE");

	//drawing
	for (int k = 0; k < 9; k++) {
		vector<double> date, ampl, ddate, dampl, effic, deffic;
		binsI = 0;
		bigE = 0;
		for (size_t i = 1; i < (size_t)Nbins + 1; i++) {
			ampl.push_back(ha[k]->GetBinContent(i));
			effic.push_back(100 * he[k]->GetBinContent(i));
			date.push_back(ha[k]->GetBinCenter(i));
			dampl.push_back(ha[k]->GetBinError(i));
			deffic.push_back(100 * he[k]->GetBinError(i));
			ddate.push_back(ha[k]->GetBinWidth(i) / 2.);
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
		gr1->GetYaxis()->SetRangeUser(85, 100);
		gr->GetYaxis()->SetRangeUser(0.0, 7);
		gr->SetLineColor(3);
		gr2->SetLineColor(4);
		gr1->SetLineColor(2);
		gr->Draw("AP");

		gr->Write(Form("amp profile sensor%d", k + 1));
		
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
	MyFile1->Close();

	MyFile->Close();
}

void CalibrExp::timespectra() {
	/// get time distributions for maximum ach, for pre-max and 0-s, with conditions as in raspr
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
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		schr[nch] = 0.;
		ach[0] = 0;
		if (!RunIsGood(run)) {
			continue;
		}
		scount = 0;
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			if (j == 8) { shiftPhi2Pi(); }
			while (schr[1 + scount] > 0.5) { scount += 1; }

			int maxampla = findMaxAmpInd(scount1, scount, &ach[0]);
			vector<double> lgran = defineLGrans(j);

			predmaxach = -1.;
			predmaxampla = scount1;
			for (int f = scount1; f <= scount; f++) {
				if (ach[f] < ach[maxampla] - 0.0001 && ach[f] > predmaxach) {
					predmaxach = ach[f];
					predmaxampla = f;
				}
			}

			for (size_t partI = 0; partI < 2; partI++) {
				ztr[partI] = 12.0 / tan(theta[partI]) + z0[partI];
				int zoblInd = findZOblInd(ztr[partI]);
				//scount 1 is always 0, so if we have some signal other than 0 
				if (scount - scount1 > 0.5 && (zoblInd > 1) && (zoblInd < 12)) {
					if (insideLGrans(phi[partI], &lgran[0]) && (schr[maxampla] > 0.5))
						histo->Fill(tchr[maxampla]);
				}
				// if it for example 0 1 2 0 - second by amplitude signal exists
				if (scount - scount1 > 1.5 && (ztr[partI] > zobl[1]) && (ztr[partI] < zobl[13])) {
					if (insideLGrans(phi[partI], &lgran[0]) && (schr[predmaxampla] > 0.5))
						histi->Fill(tchr[predmaxampla]);
				}
				if ((ztr[partI] > zobl[1]) && (ztr[partI] < zobl[13])) {
					if (insideLGrans(phi[partI], &lgran[0]) && (ach[maxampla] < 2000))
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
			maxLeft = hma->GetBinContent((Int_t)(hma->GetMaximumBin() - 2.5*f1->GetParameter(2) - i));
			maxLeftInd = i;
		}
		if (hma->GetBinContent(hma->GetMaximumBin() + (int)f1->GetParameter(2) + i) > maxRight) {
			maxRight = hma->GetBinContent((Int_t)(hma->GetMaximumBin() + 2.5*f1->GetParameter(2) + i));
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
	f += par[21]*TMath::Gaus(x[0],par[22],par[23]);
	f += par[24]*TMath::Gaus(x[0],par[25],par[26]);
	return f;
}
double scphi2(double* x, double* par)
{
	double f = 0;
	for (size_t i = 0; i < 3; i++)
		f += aerogel(i, x[0], par);
	f += shifter(x[0], par);
	return f;
}
//ends here


void CalibrExp::testforfit(string basefile, string pName ) {
	readTrueGran1();
	char title[100], titlea[100], titles[100], titlep[100], titlef[100];
	double shir[9][14][3];
	std::vector<double> xv1(12);
	double scal = PI / 180;
	double dphi = 2.4 * scal;
	TFile* f = new TFile((workingDir + basefile).c_str());
	TFile* fD = new TFile((workingDir + "zaselphi.root").c_str());
	TFile* f1 = new TFile((workingDir + "fitted" + basefile).c_str(), "RECREATE");
	TIter next(f->GetListOfKeys());
	TKey* key;
	int counterr = 0;
	string pNameFullsPred = "a";
	while ((key = (TKey*)next()) && counterr < 500) {
		counterr += 1;
		TString name(key->GetName());
		if (!name.Contains(pName))
			continue;
		
		// Fit procedure for one profile, maybe I should make it as a separate function
		{
			int obl, counter;

			// ____find the name of a file 
			string pNameFulls;
			const char* pNameFull = "a";
			if (pName == "zobl") {
				sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
				//obl = 4;
				//counter = 3;
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
			TH1* h_phi_distr = (TH1F*)fD->Get(pNameFull);
			//TGraphErrors* h_amp_phi = transformProfileByPhiDistr(h_amp_phiP, h_phi_distr);
			
			//TProfile* h_amp_phi = (TProfile*)h_amp_phiP->Clone(Form("zobl%d,sensor%d", obl, counter));
			double shiftPhi = 0;
			double shiftPhiR = 0;
			double shiftSource = 0;
			if (basefile.find("mod") == std::string::npos) {
				vector< vector<double> > shifts = transformProfileByPhiDistr(h_amp_phi, h_phi_distr);
				for (size_t i = 0; i < shifts.size(); i++) {
					if (fabs(h_amp_phi->GetBinCenter(h_amp_phi->GetMaximumBin()) - shifts[i][1]) > 0.025)
						continue;
					shiftPhi = shifts[i][0];
					shiftPhiR = shifts[i][2];
					shiftSource = shifts[i][1];
				}
			}

			h_amp_phi->Draw();
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			std::vector<double> par(27);
			
			// Fit Central shifter pick
			double xmin = (counter * 40. - 35.) * scal + dphi;
			double xmax = (counter * 40. - 15.) * scal + dphi;
			TF1* fxc = new TF1("fxc", gg, xmin, xmax, 6);
			double xm = h_amp_phi->GetBinCenter(h_amp_phi->GetMaximumBin());
			fxc->SetParameters(h_amp_phi->GetMaximum(), xm, 0.5*scal, 5, xm, 10 * scal);
			fxc->FixParameter(0, fxc->GetParameter(0));
			fxc->FixParameter(1, fxc->GetParameter(1));
			fxc->FixParameter(2, fxc->GetParameter(2));
			fxc->SetParLimits(4, xm - dphi, xm + dphi);
			fxc->SetParLimits(5, 5 * scal, 20 * scal);
			h_amp_phi->Fit(fxc, "", "", xmin, xmax);
			c->Update();
			for (size_t i = 0; i < 3; i++)
				fxc->ReleaseParameter(i);
			fxc->SetParLimits(0, 0, 100);
			fxc->SetParLimits(1, (counter * 40. - 30.) * scal + dphi, (counter * 40. - 20.) * scal + dphi);
			fxc->SetParLimits(2, 0.1 * scal, 0.6 * scal);
			h_amp_phi->Fit(fxc, "", "", xmin, xmax);
			c->Update();

			// Fit Left border
			xmin = fxc->GetParameter(1) - 0.32; //((counter - 1) * 40. - 10.) * scal + dphi;
			xmax = fxc->GetParameter(1) - 9 * fxc->GetParameter(2);//((counter-1)*40.+10.)*scal+dphi;
			TF1* fxl = new TF1("fxl", gp2x, xmin, xmax, 5);
			size_t binTresholdUp = 1;
			while (h_amp_phi->GetBinContent(binTresholdUp) < 0.5)
				binTresholdUp++;
			fxl->SetParameters(h_amp_phi->GetBinCenter(binTresholdUp), 0.5 * scal, 1.0, 1.0, 2.0);
			fxl->SetParLimits(1, 0.1 * scal, 2.2 * scal);
			h_amp_phi->Fit(fxl, "", "", xmin, xmax);
			c->Update();

			// Fit Right border
			xmin = fxc->GetParameter(1) + 9 * fxc->GetParameter(2);//(counter*40.-10.)*scal+dphi;
			xmax = fxc->GetParameter(1) + 0.5; // (counter * 40. + 12.)* scal + dphi;
			TF1* fxr = new TF1("fxr", gp2x, xmin, xmax, 5);
			fxr->SetParameters(fxc->GetParameter(1) + 0.43, -0.5 * scal, 5, 0, 0);
			fxr->SetParLimits(0, fxc->GetParameter(1) + 0.38, fxc->GetParameter(1) + 0.48);
			fxr->SetParLimits(1, -1.2 * scal, -0.3 * scal);
			h_amp_phi->Fit(fxr, "", "", xmin, xmax);
			c->Update();

			// Set Parameters
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
			par[20] = shwidth;

			c->Update();
			for (size_t i = 0; i < 4; i++) {
				TLine* linet = new TLine(par[12+i], 0, par[12+i], gPad->GetUymax());
				linet->SetLineColor(kRed);
				//linet->Draw();
				//c->Update();
			}
			//cin.get();

			TF1* fn = new TF1("scphi", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 27);
			// Set Names
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
			// Set Values
			fn->SetParameters(par.data());
			fn->SetLineColor(kRed);
			fn->SetNpx(1000);
			std::vector<double> xv(12);

			//______fit___
			
			//adding gaus around shifter
			fn->FixParameter(21, 0);
			fn->SetParameter(22, shiftPhi);
			fn->SetParLimits(22, shiftPhi - 0.0025, shiftPhi + 0.0025);
			fn->SetParLimits(23, 0.0005, max((-shiftPhi + shiftSource) / 2., 0.0005));
			fn->FixParameter(24, 0);
			fn->SetParameter(25, shiftPhiR);
			fn->SetParLimits(25, shiftPhiR - 0.0025, shiftPhiR + 0.0025);
			fn->SetParLimits(26, 0.0005, max((shiftPhiR - shiftSource) / 2., 0.0005));
			if (shiftSource < par[13] && shiftSource > 0)
				fn->SetParLimits(22, shiftPhi - 0.0025, shiftPhi + 0.007);
			if (shiftSource > par[13] && shiftSource > 0)
				fn->SetParLimits(25, shiftPhiR - 0.007, shiftPhiR + 0.0025);

			for (size_t i = 0; i < 20; i++)
				fn->FixParameter(i, fn->GetParameter(i));
			fn->SetParLimits(16, fxc->GetParameter(0) - 0.3, fxc->GetParameter(0) + 19.);

			fn->ReleaseParameter(20);
			fn->SetParLimits(20, 1.4, 4.0);

			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();

			if (shiftPhi > 0) {
				fn->ReleaseParameter(21);
				fn->SetParLimits(21, 0.1, 10.);
			}

			for (size_t i = 0; i < 4; i++)
				fn->ReleaseParameter(i);
			fn->SetParLimits(0,0.5,5.);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			//c->Update();
			if (shiftPhiR > 0) {
				fn->ReleaseParameter(24);
				fn->SetParLimits(24, 0.01,10.);
			}
			for (size_t i = 4; i < 8; i++)
				fn->ReleaseParameter(i);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			//c->Update();
			for (size_t i = 8; i < 12; i++)
				fn->ReleaseParameter(i);
			h_amp_phi->Fit("scphi", "", "", par[12] - 15 * scal, par[15] + 15 * scal);
			c->Update();

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
			fn->SetParLimits(17, 1.0e-3, 1.4e-2);
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
			c->Update();

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
					TF1* fna = new TF1(Form("aerogel%d", cc), scphi2, par[12] - 15 * scal, par[15] + 15 * scal, 21);
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
					sprintf(titlea, "zobl%d,sensor%d,aerogel%d", obl, counter, cc + 1);
					//sprintf(titlea, "date%d,sensor%d,aerogel%d", date, counter, c + 1);
					fna->Write(titlea);
					c->Update();
				}
				c->Update();
				//
				TF1* fns = new TF1("shifter", scphi, par[12] - 15 * scal, par[15] + 15 * scal, 27);
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
			}
			//cin.get();
		}
	}
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
				truegran1mErr[i][j] = gran1Err[j][i][0];//(shir[i][0][0]+shir[i][0][1]+shir[i][0][2])/3.;
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


void CalibrExp::bordersplot() {
	//find mean by z borders of aerogel parts truegran1 from obtained in fit gran1.
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
	double truegran1ErrVer2[9][4];
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

			TH1* hGr1 = new TH1F("hGr1","border distribution",20, truegran1[j][i]-0.05, truegran1[j][i] + 0.05);
			for (int k = 0; k < 14; k++) {
				hGr1->Fill(gran1[i][j][k]);
				if (fabs(gran1[i][j][k]-truegran1[j][i]) > maximum[j][i]) maximum[j][i] = fabs(gran1[i][j][k]-truegran1[j][i]);
				if (fabs(gran1[i][j][k]-truegran1[j][i])*scal1>0.9){
					badfit1[j][i].push_back(k+1);
				}
				if (gran1[i][j][k] < minimum[j][i]) minimum[j][i] = gran1[i][j][k];
			}
			cout << "	" << Form("%.1f",(maximum[j][i])*scal1);
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			//hGr1->Draw();
			truegran1ErrVer2[j][i] = hGr1->GetRMS();
			//c->Update();
			//cin.get();
			hGr1->Delete();
		}
		
	}
	f1->Close();
	writeTrueGran1();

	//info about changes in borders
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

	//info about difference in fit-profile
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
		if (!name.Contains("zobl"))
			continue;

		int obl, counter;
		sscanf(name.Data(), "zobl%d,sensor%d", &obl, &counter);
		TProfile* h_amp_phi = (TProfile*)f->Get(Form("zobl%d,sensor%d", obl, counter));
		std::vector<double> par;
		for(size_t k =0;k<21;k++)
			par.push_back(fitPar[counter-1][obl-1][k]);
		TF1* fn = new TF1("scphi", scphi2, par[12] - 15/scal1, par[15] + 15/scal1, 21);
		double maxdiffL = 0;
		double maxdiffR = 0; 
		double maxdiffM = 0;
		fn->SetParameters(par.data());
		for(size_t i = 1; i < (size_t)h_amp_phi->GetNbinsX()-1; i++){
			double x0 = h_amp_phi->GetXaxis()->GetBinCenter(i);
			double x1 = h_amp_phi->GetXaxis()->GetBinCenter(i-1);
			double x2 = h_amp_phi->GetXaxis()->GetBinCenter(i+1);
			double erri = h_amp_phi->GetBinError(i);
			double hy0 = h_amp_phi->GetBinContent(i);
			double hy1 = h_amp_phi->GetBinContent(i-1);
			double hy2 = h_amp_phi->GetBinContent(i+1);
			double fy0 = scphi2(&x0,&par[0]);
			double fy1 = scphi2(&x1,&par[0]);
			double fy2 = scphi2(&x2,&par[0]);
			double bord1 = gran1[0][counter-1][obl-1] + 3.*gran1Err[0][counter-1][obl-1];
			double bord2 = gran1[1][counter-1][obl-1] - 3.*gran1Err[1][counter-1][obl-1];
			double bord3 = gran1[1][counter-1][obl-1] + 3.*gran1Err[1][counter-1][obl-1];
			double bord4 = gran1[3][counter-1][obl-1] - 3.*gran1Err[3][counter-1][obl-1];

			//cout << fy0 << endl;
			if (x0>0.4 && fabs(x0 - gran1[0][counter - 1][obl - 1]) < 2. * gran1Err[0][counter - 1][obl - 1] && maxdiffL < fabs(hy0 - fy0) / erri) {
				maxdiffL = fabs(hy0 - fy0) / erri;
			}
			if(fabs(x0-gran1[3][counter-1][obl-1]) < 2.*gran1Err[3][counter-1][obl-1] && maxdiffR < fabs(hy0 - fy0)/erri )
				maxdiffR = fabs(hy0 - fy0)/erri;
			if(((x0 > bord1 && x0 < bord2) || (x0 > bord3 && x0 < bord4))  ) {
				double locmeandiff = ((hy0 - fy0) + (hy1 - fy1) + (hy2 - fy2))/(3.*erri);
				if (maxdiffM < fabs(locmeandiff))
					maxdiffM = fabs(locmeandiff);
			}
		}
		maxdiffArrL[counter-1][obl-1] = maxdiffL;
		maxdiffArrR[counter-1][obl-1] = maxdiffR;
		maxdiffArrM[counter-1][obl-1] = maxdiffM;
		
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

	//graph of fluctuations in borders
	double tg1t[9];
	double tg1tErr[9];
	double counterInd[9];
	double counterIndErr[9];
	for (int i = 0; i < 9; i++) {
		counterInd[i] = (double)i + 1;
		counterIndErr[i] = 0.;
		tg1t[i] = truegran1[i][1];
		tg1tErr[i] = truegran1ErrVer2[i][1];
	}
	TGraphErrors* grSh = new TGraphErrors(9, counterInd, tg1t, counterIndErr, tg1tErr);
	grSh->SetMarkerStyle(21);
	grSh->SetMarkerSize(1);
	grSh->Draw("AP");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	cin.get();
	cin.get();
	TF1* fngrSh = new TF1("fngrSh", "[0]+[1]*x", 0, 10);
	grSh->Fit("fngrSh");
	for (int i = 0; i < 9; i++)
		tg1t[i] = tg1t[i]-(fngrSh->GetParameter(0)+counterInd[i]* fngrSh->GetParameter(1));
	TGraphErrors* grShCorr = new TGraphErrors(9, counterInd, tg1t, counterIndErr, tg1tErr);
	grShCorr->SetMarkerStyle(21);
	grShCorr->SetMarkerSize(1);
	grShCorr->Draw("AP");
}


void CalibrExp::modelfile(bool transformPYN) {
	//find shifter array from pik, stored in gran1 file and also use it, Yfm and truegran to create file with parameters for modeling database
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
		if (j == 4)
			xv1[4] = truegran1m[j][1]-0.007;
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
						//x1 = xv1[i + 4 * n];
						x1 = x1 + xv1[4] - xvt[4];
					//just polinoms
					yfm[j][k][4 * n + i] = pv1[0] + pv1[1] * x1 + pv1[2] * pow(x1, 2) + pv1[3] * pow(x1, 3);
					//yfm[j][k][4 * n + i] = 2-i%2;
					//cout << yfm[j][k][4 * n + i] / fitPar[j][k][4 * n + i] << endl;
					//aerogel with exponents (fit f without shifters)
					//yfm[j][k][4 * n + i] = aerogel(n, x1, fitPar[j][k]);
					//cout << yfm[counter - 1][obl - 1][4 * n + i] << endl;
				}
			}
			TF1* fn = new TF1("scphi", scphi2, fitPar[j][k][12] - 15 * scal, fitPar[j][k][15] + 15 * scal, 21);
			for (size_t i = 0; i < 21; i++)
				fn->SetParameter(i, fitPar[j][k][i]);
			fn->Draw();
			TF1* fn1 = new TF1("scphi1", scphi2, fitPar[j][k][12] - 15 * scal, fitPar[j][k][15] + 15 * scal, 21);
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
				fout1 << zobl[i] << " " << yfm[j][i][0 + 4 * k] / coef[i][j] << " " << yfm[j][i][1 + 4 * k] / coef[i][j] << " " << yfm[j][i][2 + 4 * k] / coef[i][j] << " " << yfm[j][i][3 + 4 * k] / coef[i][j] << endl;
			}
		}
	}
	fout1.close();
	writeMcoef();
}


void CalibrExp::getlist(int basedirInd, string hbookdir) {
	readMinMaxRun();
	ofstream fout;
	fout.open((workingDir + "simflist2.fwi").c_str());

	vector<TString> ans;
	TString basedir;
	if (basedirInd == 1)
		basedir = basedirMod + "*.hbook";
	else
		basedir = "/online/simulation/MC/R006-004/ee/output/*.mod.gz";
	TString files = gSystem->GetFromPipe("ls " + basedir);
	TString elem;
	Ssiz_t from = 0;
	while (files.Tokenize(elem, from, "\n")) {
		ans.push_back(elem);
	}
	int vector_size = ans.size();
	TString str = "a";
	TString str1 = "b";
	fout << "source /work/snd2000/root/setup2k.sh i386-SL5-opt-debug" << endl << "source/etc/sysconfig/gridengine_nsu" << endl << endl;
	for (int i = 0; i < vector_size; i++) {
		str = ans[i].Copy();
		str1 = ans[i].Copy();
		TString tok;
		TString tok1;
		Ssiz_t from = 0;
		Ssiz_t from1 = 0;
		vector<TString> tokens;
		vector<TString> tokens1;
		while (str.Tokenize(tok, from, "[/.-]")) {
			//cout << tok << endl;
			tokens.push_back(tok);
		}
		while (str1.Tokenize(tok1, from1, "[/.]")) {
			//cout << tok1 << endl;
			tokens1.push_back(tok1);
		}

		//cout << tokens1[6] << endl;

		const char* a0 = (const char*)tokens[tokens.size()-4];
		int a = atoi(a0);
		//cout << a << endl;
		if (a > 27250 && a > MinimRun && a < MaximRun) {
			//cout << a << endl;
			fout << ".mainrelease/Offline/submit.sh -q clusters,1440 FAKERUN=" << a << " MODFILENAME=" << tokens1[tokens1.size()-2] << "     MODFILEDIR=/online/simulation/MC/R006-004/ee/output/ HBOOKDIR=./" << hbookdir << "/  SimRecApp fwk/simreco_col_point.fwk" << endl;
			cout << "h2root " << basedirMod << tokens1[tokens1.size()-2] << ".hbook " << basedirMod << tokens1[tokens1.size() - 2] << ".root" << (char)59 << " ";
		}

	}
	cout << endl;
	fout.close();
}


//odbsolent
double CalibrExp::compare1() {
	/// Functions to yield coefficients for the map by comparing mean values of amplitude for 1 counter and 1 z area
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
	TH1* h1[9][14];
	double h[9][14];
	double hcount[9][14];
	double hsquare[9][14];
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
			h1[j][i] = new TH1F(name, title, 5000, -50, 2*achCut);
			h[j][i] = 0.;
			hcount[j][i] = -0.5;
			hsquare[j][i] = 0.;
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
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		schr[nch] = 0.;
		ach[0] = 0;
		if ((!RunIsGood(run))/* || eventtime != 0*/) {
			continue;
		}
		scount = 0;     //2017: 4001, 25500, 29501    ; 2018: 6001, 29800, 35801       ; 2019: 38000-43001       ; 2020: 45040-46541
		scount1 = 0;
		for (int j = 0; j < 9; j++) {
			//wichcounter1[j] = scount1;
			if (j == 8) { shiftPhi2Pi(); }
			while (schr[1 + scount] > 0.5) { scount += 1; }

			int maxampla = findMaxAmpInd(scount1, scount, &ach[0]);
			vector<double> lgran = defineLGrans(j);

			for (size_t pind = 0; pind < 2; pind++) {
				ztr[pind] = 12.0 / tan(theta[pind]) + z0[pind];
				zin[pind] = aerRLowBord / tan(theta[pind]) + z0[pind];
				int zoblInd = findZOblInd(ztr[pind]);

				if ( insideLGrans(phi[pind], &lgran[0]) && (TimeIsGood(maxampla)) && zoblInd != -1) {

					double ach1 = ach[maxampla];
					meancampl[zoblInd]->Fill(ach[maxampla]);
					if (zoblInd > 2 && zoblInd < 12)
						meanzampl[j]->Fill(ach[maxampla]);

					double fillAmplitude = recalculateZ(zin[pind], theta[pind], ach[maxampla], 0, 1);
					int n = count1 * 1000000 + Count;
					if (fillAmplitude != -1.0 && ach[maxampla] < achCut) {
						h1[j][zoblInd]->Fill(ach[maxampla]);
						//cout << fillAmplitude << endl;
						//h[j][zoblInd] = h[j][zoblInd] * hcount[j][zoblInd] / (hcount[j][zoblInd] + 1) + fillAmplitude / (hcount[j][zoblInd] + 1);
						//cout << h[j][zoblInd] << endl;
						h[j][zoblInd] += fillAmplitude;
						hcount[j][zoblInd] += 1;
						hsquare[j][zoblInd] += pow(fillAmplitude,2);
					}

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
				h1[j][i]->Write(title);
			}
		}
		for (int j = 0; j < 9; j++) {
			sprintf(title, "e_count%d", j + 1);
			meanzampl[j]->Write(title);
		}
		MyFile->Close();

		ofstream fout1;
		fout1.open((workingDir + "meanAchExp.dat").c_str());
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++)
				fout1 << h[j][i]/hcount[j][i] << "	";
			fout1 << endl;
		}
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++)
				fout1 << sqrt( hsquare[j][i] / hcount[j][i] - pow(h[j][i] / hcount[j][i],2) ) / sqrt(hcount[j][i]) << "	";
			fout1 << endl;
		}
		fout1.close();
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
	/// same as previous, just extracting mean values for modeling
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
	TH1* h1[9][14];
	double h[9][14];
	double hcount[9][14];
	double hsquare[9][14];
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
			h1[j][i] = new TH1F(name, title, 5000, -50., 2*achCut);
			h[j][i] = 0.;
			hcount[j][i] = -0.5;
			hsquare[j][i] = 0.;
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
			if (j == 8) { shiftPhi2Pi(); }

			vector<double> lgran = defineLGrans(j);

			for (size_t pind = 0; pind < 2; pind++) {
				int zoblInd = findZOblInd(ztr[pind]);

				//if ((((phi[pind] > lgran[0]) && (phi[pind] < lgran[1])) || ((phi[pind] > lgran[2]) && (phi[pind] < lgran[3]))) && zoblInd != -1) {
				if (insideLGrans(phi[pind], &lgran[0]) && zoblInd != -1) {

					meancampl[zoblInd]->Fill(ach[j]);
					if (zoblInd > 2 && zoblInd < 12)
						meanzampl[j]->Fill(ach[j]);

					double fillAmplitude = recalculateZ(zin[pind], theta[pind], ach[j], 0, 2);
					int n = count1 * 1000000 + Count;
					if (fillAmplitude != -1.0 && ach[j] < achCut) {
						h1[j][zoblInd]->Fill(ach[j]);
						//h[j][zoblInd] = h[j][zoblInd] * hcount[j][zoblInd] / (hcount[j][zoblInd]+1) + fillAmplitude / (hcount[j][zoblInd] + 1);
						h[j][zoblInd] += fillAmplitude;
						hcount[j][zoblInd] += 1;
						hsquare[j][zoblInd] += pow(fillAmplitude, 2);
					}

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
				h1[j][i]->Write(title);
			}
		}
		for (int j = 0; j < 9; j++) {
			sprintf(title, "m_count%d", j + 1);
			meanzampl[j]->Write(title);
		}
		MyFile->Close();

		ofstream fout1;
		fout1.open((workingDir + "meanAchMod.dat").c_str());
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++)
				fout1 << h[j][i]/hcount[j][i] << "	";
			fout1 << endl;
		}
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 9; j++)
				fout1 << sqrt(hsquare[j][i] / hcount[j][i] - pow(h[j][i] / hcount[j][i], 2)) / sqrt(hcount[j][i]) << "	";
			fout1 << endl;
		}
		fout1.close();
	}
	writeComparisonm();
	return 2;
}
void CalibrExp::compare() {
	/// comparing these mean values
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

	double he[9][14];
	double hce[9][14];
	double hm[9][14];
	double hcm[9][14];
	ifstream fin1;
	fin1.open((workingDir + "meanAchExp.dat").c_str());
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++)
			fin1 >> he[j][i];
		fin1.get();
	}
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++)
			fin1 >> hce[j][i];
		fin1.get();
	}
	fin1.close();
	ifstream fin2;
	fin2.open((workingDir + "meanAchMod.dat").c_str());
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++)
			fin2 >> hm[j][i];
		fin2.get();
	}
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++)
			fin2 >> hcm[j][i];
		fin2.get();
	}
	fin2.close();

	TFile* Myf = new TFile((workingDir + "achspectreme.root").c_str());
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 9; j++) {
			/*TH1* meanamplocke = (TH1F*)Myf->Get(Form("e_all_zobl%d_count%d", i + 1, j + 1));
			TH1* meanamplockm = (TH1F*)Myf->Get(Form("m_all_zobl%d_count%d", i + 1, j + 1));
			meanample[i][j] = meanamplocke->GetMean(1);
			meanampleerr[i][j] = meanamplocke->GetMeanError(1);
			meanamplm[i][j] = meanamplockm->GetMean(1);
			meanamplmerr[i][j] = meanamplockm->GetMeanError(1);*/
			meanample[i][j] = he[j][i];
			meanampleerr[i][j] = hce[j][i];
			meanamplm[i][j] = hm[j][i];
			meanamplmerr[i][j] = hcm[j][i];
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
		//gr->Fit("fn", "", "", -9, 9);
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


TDirectory* createDir(const char* name, TFile* file) {
	/// function to delete and create dir
	if (file->FindKey(name))
		file->rmdir(name);
	return file->mkdir(name);
}

void CalibrExp::compareAmpSpectr() {
	/// a function to extract map coefficients by dividing phi profiles and fitting resulting dependence with a line, lying around 1
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

	TFile* MyFile = new TFile((workingDir + "resultingComparisons.root").c_str(), "UPDATE");
	TDirectory* dir = createDir("amplitude phi profiles ratio fit", MyFile);
	dir->cd();

	TCanvas* c = new TCanvas("c", "c", 1200, 800);
	c->SetFillColor(0);
	c->SetBorderMode(0);
	TPad* pad1 = new TPad("pad1", "", 0, 0, 1, 1, 0, 4, 0);

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


			vector<double> lgran = defineLGrans(counter);
			//TF1* tf1 = (TF1*)d1->Get(Form("zobl%d,sensor%d,full", 2, 8));
			//TProfile* hprof1 = (TProfile*)f1->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));
			//TProfile* hprof1 = (TProfile*)f1->Get(Form("mean,sensor%d", counter + 1));
			TProfile* hprof1 = (TProfile*)f1->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));
			TProfile* hprof2 = (TProfile*)f2->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));
			hprof2->Rebin(4);
			//TProfile* hprof2 = (TProfile*)f2->Get(Form("mean,sensor%d", counter + 1));
			TH1* hz = (TH1F*)f3->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));

			
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
			TH1* hRelation = new TH1F("hRelation", "ratio", 100, 0.75, 1.25);
			for(size_t i = 2; i < (size_t)hprof2->GetNbinsX()-1; i++){
				if (hprof2->GetBinEntries(i) <= 0 || !insideLGrans(hprof2->GetBinLowEdge(i), &lgran[0]) || !insideLGrans(hprof2->GetBinLowEdge(i) + hprof2->GetBinWidth(i), &lgran[0]))
					continue;
				{
				x.push_back(hprof2->GetBinCenter(i));
				//cout << "asdasd" << endl;
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
				hRelation->Fill(y.back());

				z.push_back(modA / scphi(&x.back(),&fitPar[counter][obl][0]));
				dz.push_back(sqrt(pow(hprof1->GetBinError(modBin) / scphi(&x.back(), &fitPar[counter][obl][0]), 2) + pow(z.back() * 0. / scphi(&x.back(), &fitPar[counter][obl][0]), 2)));
				//z.push_back(hprof1->GetBinContent(i) / scphi(&x.back(),&fitPar[counter][obl][0]));
				//dz.push_back(sqrt(pow(hprof1->GetBinError(i) / scphi(&x.back(), &fitPar[counter][obl][0]), 2) + pow(z.back() * 0. / scphi(&x.back(), &fitPar[counter][obl][0]), 2)));
				}
			}
			cout << (hRelation->GetMean()-1)*100 << "		" << hRelation->GetRMS() * 100 << endl;
			hRelation->Delete();
			TGraph* graph = new TGraphErrors(x.size(), &x[0], &y[0], &dx[0], &dy[0]);
			TGraph* graph1 = new TGraphErrors(x.size(), &x[0], &z[0], &dx[0], &dz[0]);
			graph->SetLineColor(1);
			graph1->SetLineColor(1);
			graph->SetTitle(Form("counter%d, zobl%d m/e;#phi;ratio",counter+1,obl+1));
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
			//hz->DrawNormalized("same",400);
			lineg1->Draw("same");
			lineg2->Draw("same");
			lineg3->Draw("same");
			lineg4->Draw("same");
			TF1* f1 = new TF1("f1", "[0]", lgran[0], lgran[3]);
			f1->SetParameter(0,1.0);
			graph->Fit("f1", "", "", lgran[0], lgran[3]);
			c->Update();
			//cin.get();

			graph->Write(Form("ratio zarea%d sensor%d", obl + 1, counter + 1));

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


	TDirectory* dir1 = createDir("map coefficients", MyFile);
	dir1->cd();
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
		gr1->SetTitle("map coefficients (m/e);z,cm;coeff");
		
		//gr->Draw("AP");
		gr1->Draw("AP");

		gr1->Write(Form("coefficients sensor%d", j + 1));

		fn->SetLineColor(kRed);
		//gr->Fit("fn", "", "", -9, 9);
		c->Update();
		//cin.get();
		//cout << fn->GetParameter(0) << endl;
	}
	MyFile->Close();
	for (size_t obl = 0; obl < 14; obl++) {
		for (size_t counter = 0; counter < 9; counter++) {
			coef[obl][counter] = (koefnewp[counter][obl] + koefnewp[counter][obl])/2.;
			//coef[obl][counter] = 1.1;
		}
	}
	writeNcoef();
	c->Delete();
}
///
void CalibrExp::compareEffSpectr() {
	readMeanshir();
	readTrueGran1();
	readMeanshirm();
	readTruegran1m();

	TFile* fmEff = new TFile((workingDir + "profilesEffMod.root").c_str());
	TFile* feEff = new TFile((workingDir + "profilesEff.root").c_str());

	TFile* fmAmp = new TFile((workingDir + "profilesmmod.root").c_str());
	TFile* feAmp = new TFile((workingDir + "profiles.root").c_str());

	TFile* MyFile = new TFile((workingDir + "resultingComparisons.root").c_str(), "UPDATE");
	TDirectory* dir = createDir("efficiency and amp phi profiles ratios", MyFile);
	TDirectory* dir1 = createDir("efficiency profiles on top of one another", MyFile);
	dir->cd();

	TCanvas* c = new TCanvas("c", "c", 1200, 800);
	c->SetFillColor(0);
	c->SetBorderMode(0);
	TPad* pad1 = new TPad("pad1", "", 0, 0, 1, 1, 0, 4, 0);

	cout << fixed;
	cout.precision(3);
	cout << "max-min(%)	error(%)	mean m/e	stdValue(%)" << endl;
	for (size_t zobl = 0; zobl < 1; zobl++) {
		for (size_t counter = 0; counter < 9; counter++) {
			TProfile* hpmEff = (TProfile*)fmEff->Get(Form("mean,sensor%d", counter + 1));
			TProfile* hpeEff = (TProfile*)feEff->Get(Form("mean,sensor%d", counter + 1));
			TProfile* hpmAmp = (TProfile*)fmAmp->Get(Form("mean,sensor%d", counter + 1));
			TProfile* hpeAmp = (TProfile*)feAmp->Get(Form("mean,sensor%d", counter + 1));
			hpmEff->Rebin(2);
			hpeEff->Rebin(2);
			hpmAmp->Rebin(2);
			hpeAmp->Rebin(2);
			hpmEff->SetLineColor(2);
			hpmEff->SetLineWidth(3);
			hpeEff->SetLineWidth(2);

			vector<double> lgran = defineLGrans(counter);
			vector<double> xa, xae, a, ae, xe, xee, e, ee, ue, uee;
			double maxValue = 0; double minValue = 2; 
			double maxErr = 0;
			double minErr = 0;
			TH1* hRelation = new TH1F("hRelation","ratio",100,0.75,1.25);
			for (size_t i = 2; i < (size_t)hpmAmp->GetNbinsX() - 1; i++) {
				if (hpeAmp->GetBinContent(i) < 0.1 || !insideLGrans(hpmAmp->GetBinCenter(i), &lgran[0]))
					continue;

				xa.push_back(hpmAmp->GetBinCenter(i));
				xae.push_back(hpmAmp->GetBinWidth(i));

				a.push_back(hpmAmp->GetBinContent(i) / hpeAmp->GetBinContent(i));
				ae.push_back(divisionError(hpmAmp->GetBinContent(i), hpeAmp->GetBinContent(i), hpmAmp->GetBinError(i), hpeAmp->GetBinError(i)));
			}
			for (size_t i = 2; i < (size_t)hpmEff->GetNbinsX() - 1; i++) {
				if (hpeEff->GetBinContent(i) < 0.5 || !insideLGrans(hpmEff->GetBinCenter(i), &lgran[0]))
					continue;

				xe.push_back(hpmEff->GetBinCenter(i));
				xee.push_back(hpmEff->GetBinWidth(i));

				ue.push_back((1- hpmEff->GetBinContent(i)) / (1- hpeEff->GetBinContent(i)));
				uee.push_back(divisionError(1- hpmEff->GetBinContent(i), 1- hpeEff->GetBinContent(i), hpmEff->GetBinError(i), hpeEff->GetBinError(i)));

				e.push_back((hpmEff->GetBinContent(i)) / (hpeEff->GetBinContent(i)));
				ee.push_back(divisionError(hpmEff->GetBinContent(i), hpeEff->GetBinContent(i), hpmEff->GetBinError(i), hpeEff->GetBinError(i)));

				if (e.back() > maxValue) {
					maxValue = e.back();
					maxErr = ee.back();
				}
				if (e.back() < minValue) {
					minValue = e.back();
					minErr = ee.back();
				}
				hRelation->Fill(e.back());
			}
			TGraphErrors* gra = new TGraphErrors(xa.size(), &xa[0], &a[0], &xae[0], &ae[0]);
			TGraphErrors* grue = new TGraphErrors(xe.size(), &xe[0], &ue[0], &xee[0], &uee[0]);
			TGraphErrors* gre = new TGraphErrors(xe.size(), &xe[0], &e[0], &xee[0], &ee[0]);

			gra->SetMarkerStyle(22);
			gra->SetMarkerColor(2);
			gre->SetMarkerStyle(22);
			grue->SetMarkerStyle(22);

			gra->SetName("gra");
			gre->SetName("gre");

			gre->Draw("AP");
			cout << (maxValue-minValue)/2.*100 << "		" << sqrt(pow(minErr,2)+pow(maxErr, 2))/2.*100 << "		" << hRelation->GetMean() << "		" << hRelation->GetRMS()*100 << endl;
			gra->Draw("Psame");

			double maxY = std::max(TMath::MaxElement(xa.size(), gra->GetY()), TMath::MaxElement(xe.size(), gre->GetY()));
			double minY = std::min(TMath::MinElement(xa.size(), gra->GetY()), TMath::MinElement(xe.size(), gre->GetY()));

			gre->GetYaxis()->SetRangeUser(minY, maxY);

			TLegend* legend = new TLegend(0, 0, 0.2, 0.2);
			legend->AddEntry("gra", "Amp relat m/e", "lep");
			legend->AddEntry("gre", "Eff relat m/e", "lep");
			legend->Draw();

			//TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			//cin.get();

			dir->cd();
			c->Write(Form("ratio mod/exp sensor%d", counter + 1));
			grue->SetTitle("undefficiency relation m/e;#phi");
			grue->Draw("AP");
			c->Update();
			c->Write(Form("uneff ratio mod/exp sensor%d", counter + 1));

			hpmEff->Draw();
			hpmEff->GetYaxis()->SetRangeUser(0.8, 1.1);
			hpeEff->Draw("same");
			c->Update();
			dir1->cd();
			c->Write(Form("profiles sensor%d", counter + 1));
			hRelation->Write(Form("effRelation, sensor%d", counter + 1));
			hRelation->Delete();

		}
	}
	c->Delete();
	MyFile->Close();
}
/// 
void CalibrExp::compareZAmpSpectr() {
	TFile* f1 = new TFile((workingDir + "zprofilesmod.root").c_str());
	TFile* f2 = new TFile((workingDir + "zprofiles.root").c_str());

	TFile* MyFile = new TFile((workingDir + "resultingComparisons.root").c_str(), "UPDATE");
	TDirectory* dir = createDir("amplitude z profiles ratio", MyFile);
	dir->cd();

	for (size_t counter = 0; counter < 9; counter++) {
		TProfile* hprof1 = (TProfile*)f1->Get(Form("sensor%d", counter + 1));
		TProfile* hprof2 = (TProfile*)f2->Get(Form("sensor%d", counter + 1));
		hprof1->SetLineColor(2);
		hprof1->Draw();
		hprof2->Draw("same");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		for (size_t i = 0; i < 15; i++) {
			TLine* lineg = new TLine(zobl[i], 0, zobl[i], 5.0);
			lineg->SetLineColor(4);
			lineg->Draw("same");
		}
		c->Update();
		//cin.get();

		c->Write(Form("ratio sensor%d", counter + 1));
	}
	MyFile->Close();
}
/// 
void CalibrExp::compareAmpDistributions() {
	TFile* f1 = new TFile((workingDir + "distributionsmod.root").c_str());
	TFile* f2 = new TFile((workingDir + "distributions.root").c_str());

	TFile* MyFile = new TFile((workingDir + "resultingComparisons.root").c_str(), "UPDATE");
	TDirectory* dir = createDir("amplitude distribution comparison", MyFile);
	dir->cd();

	for (size_t obl = 0; obl < 14; obl++) {
		for (size_t counter = 0; counter < 9; counter++) {
			TH1* hm = (TProfile*)f1->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));
			TH1* he = (TProfile*)f2->Get(Form("zobl%d,sensor%d", obl + 1, counter + 1));
			hm->SetName("modelHist");	hm->SetStats(0);	hm->SetLineColor(2);
			he->SetName("experHist");	he->SetStats(0);	he->SetLineWidth(2);

			he->DrawNormalized("E",hm->GetEntries());
			hm->DrawNormalized("same", hm->GetEntries());

			float meanMValue = hm->GetMean();
			float meanEValue = he->GetMean();
			TLine* lm = new TLine(meanMValue,0,meanMValue, gPad->GetUymax());
			TLine* le = new TLine(meanEValue,0,meanEValue, gPad->GetUymax());
			lm->SetLineColor(2);
			lm->Draw();
			le->Draw();

			TLegend* legend = new TLegend(0.8, 0.8, 1.0, 1.0);
			legend->AddEntry("modelHist", "Modeling", "lep");
			legend->AddEntry("experHist", "Experiment", "lep");
			legend->Draw();

			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			//gPad->SetLogy();
			c->Update();
			//cin.get();

			c->Write(Form("distrs zarea%d sensor%d", obl + 1, counter + 1));
		}
	}
	MyFile->Close();
}

void CalibrExp::compareAmpSpectrG() {
	/// draw images of amplitude phi profiles on top of one another to give a clear understanding of how close modeling is to exp
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

	TFile* MyFile = new TFile((workingDir + "resultingComparisons.root").c_str(), "UPDATE");
	TDirectory* dir = createDir("amplitude profiles on top of one another", MyFile);
	dir->cd();

	for (size_t obl = 0; obl < 14; obl++) {
		for (size_t cc = 0; cc < 9; cc++) {
			double scal = PI / 180;
			TF1* fn = new TF1(Form("scphi%d", obl), scphi2, fitPar[4][obl][12] - 15 * scal, fitPar[4][obl][15] + 15 * scal, 21);
			for (size_t i = 0; i < 12; i++)
				fn->SetParameter(i, 2 - i % 2);
			for (size_t i = 12; i < 16; i++)
				fn->SetParameter(i, truegran1[4][i - 12]);
			for (size_t i = 16; i < 20; i++)
				fn->SetParameter(i, fitPar[4][obl][i]);
			fn->SetParameter(20, 0);
			//fn->Draw("same");
			TProfile* hm = (TProfile*)f1->Get(Form("zobl%d,sensor%d", obl+1,cc+1));
			TProfile* he = (TProfile*)f2->Get(Form("zobl%d,sensor%d", obl+1,cc+1));
			hm->SetLineColor(2);
			he->Draw();
			hm->Draw("same");
			c->Update();
			
			TLine* lineg1 = new TLine(truegran1[cc][0], 0, truegran1[cc][0], 5.0);
			TLine* lineg2 = new TLine(truegran1m[cc][0], 0, truegran1m[cc][0], 5.0);
			TLine* lineg3 = new TLine(truegran1[cc][1], 0, truegran1[cc][1], 5.0);
			TLine* lineg4 = new TLine(truegran1m[cc][1], 0, truegran1m[cc][1], 5.0);
			TLine* lineg5 = new TLine(truegran1[cc][2], 0, truegran1[cc][2], 5.0);
			TLine* lineg6 = new TLine(truegran1m[cc][2], 0, truegran1m[cc][2], 5.0);
			TLine* lineg7 = new TLine(truegran1[cc][3], 0, truegran1[cc][3], 5.0);
			TLine* lineg8 = new TLine(truegran1m[cc][3], 0, truegran1m[cc][3], 5.0);
			lineg2->SetLineColor(2);
			lineg4->SetLineColor(2);
			lineg6->SetLineColor(2);
			lineg8->SetLineColor(2);
			lineg1->Draw("same");
			lineg2->Draw("same");
			lineg3->Draw("same");
			lineg4->Draw("same");
			lineg5->Draw("same");
			lineg6->Draw("same");
			lineg7->Draw("same");
			lineg8->Draw("same");
			c->Update();

			c->Write(Form("canvas zArea%d sensor%d",obl+1, cc + 1));

			//cin.get();
		}
	}
	MyFile->Close();
	/*TLine* lineg1 = new TLine(truegran1[4][0], 0, truegran1[4][0], 5.0);
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
	lineg6->Draw("same");*/
	c->Update();
}

void CalibrExp::comparePhiShifts() {
	/// draws distance between borders for exp and modeling, so we can check how severe or not systematics of angle are 
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

	TFile* MyFile = new TFile((workingDir + "resultingComparisons.root").c_str(), "UPDATE");
	TDirectory* dir = createDir("phi shifts", MyFile);
	dir->cd();

	for (size_t i = 0; i < 9; i++) {
		phiBE[i] = truegran1[i][2] * 180. / 3.14159;
		phiBEE[i] = truegran1mErr[i][2] * 180. / 3.14159 / 2.;
		if (i != 8) {
			shiftBM[i] = (truegran1m[i][3] - truegran1m[i][0]) * 180. / 3.14159;
			shiftBE[i] = (truegran1[i][3] - truegran1[i][0]) * 180. / 3.14159;
			shiftBME[i] = max(truegran1mErr[i][0], truegran1mErr[i][3]) * 180. / 3.14159;
			shiftBEE[i] = max(truegran1Err[i][0], truegran1Err[i][3]) * 180. / 3.14159;
		}
		else if (i == 8) {
			shiftBM[i] = (truegran1m[i][3] - truegran1m[i][0]) * 180. / 3.14159;
			shiftBE[i] = (truegran1[i][3] - truegran1[i][0]) * 180. / 3.14159;
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

		phiSEE[3 * i] = 0.4;
		phiSEE[3 * i+1] = 0.4;
		phiSEE[3 * i+2] = 0.4;
	}
	TGraphErrors* gro = new TGraphErrors(27, &phiSE[0], &shift[0], &phiSEE[0], &shiftE[0]);
	TGraphErrors* grbm = new TGraphErrors(9, &phiBE[0], &shiftBM[0], &phiBEE[0], &shiftBME[0]);
	TGraphErrors* grbe = new TGraphErrors(9, &phiBE[0], &shiftBE[0], &phiBEE[0], &shiftBEE[0]);
	grbm->SetLineColor(2);
	grbe->SetLineColor(3);
	grbm->SetMarkerStyle(22);
	grbm->SetMarkerSize(1);
	grbm->SetMarkerColor(kBlack);
	grbe->SetMarkerStyle(21);
	grbe->SetMarkerSize(1);
	grbe->SetMarkerColor(kBlack);
	//gro->Draw("AP");
	grbm->Draw("AP");
	//grbe->Draw("Psame");

	grbm->Write("between borders modeling");
	grbe->Write("between borders experiment");
	
	gro->SetTitle("difference between borders mod-exp;#phi#circ;#phi#circ");
	gro->Draw("AP");
	TF1* f1 = new TF1("f1", "[0]*cos((x-[1])*3.14159/180.)+[2]", 0, 360);
	f1->SetParameters(0.5,0,-0.5);
	gro->Fit("f1");
	cout << fixed;
	cout.precision(2);
	string sign1 = -f1->GetParameter(1) < 0 ? " - " : " + ";
	string sign2 = f1->GetParameter(2) < 0 ? " - " : " + ";
	cout << "fit function is: " << f1->GetParameter(0) << "*cos(phi" << sign1 << fabs(f1->GetParameter(1)) << ")" << sign2 << fabs(f1->GetParameter(2)) << endl;

	gro->Write("between mod and exp");

	MyFile->Close();
}


void compareShitRandTemp() {
	/// compare theta distributions for exp and modeling
	TFile* f1 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2019/profilesEffMod.root");
	TFile* f2 = new TFile("/work/users/kladov/snd2k/R006-003/maindir/2019/profilesEff.root");
	TH1* hm0 = (TH1F*)f1->Get("theta0");
	TH1* he0 = (TH1F*)f2->Get("theta0");
	TH1* hm = (TH1F*)f1->Get("theta");
	TH1* he = (TH1F*)f2->Get("theta");

	vector<double> xa, xae, a, ae, xe, xee, e, ee;
	for (size_t i = 2; i < (size_t)hm0->GetNbinsX() - 1; i++) {
		if (he0->GetBinContent(i) == 0)
			continue;

		xa.push_back(hm0->GetBinCenter(i));
		xae.push_back(hm0->GetBinWidth(i));
		
		a.push_back(hm0->GetBinContent(i) / he0->GetBinContent(i));
		ae.push_back(divisionError(hm0->GetBinContent(i), he0->GetBinContent(i), hm0->GetBinError(i), he0->GetBinError(i)));
	}
	for (size_t i = 2; i < (size_t)hm->GetNbinsX() - 1; i++) {
		if (he->GetBinContent(i) == 0)
			continue;

		xe.push_back(hm->GetBinCenter(i));
		xee.push_back(hm->GetBinWidth(i));

		e.push_back(hm->GetBinContent(i) / he->GetBinContent(i));
		ee.push_back(divisionError(hm->GetBinContent(i), he->GetBinContent(i), hm->GetBinError(i), he->GetBinError(i)));
	}
	TGraphErrors* gra = new TGraphErrors(xa.size(), &xa[0], &a[0], &xae[0], &ae[0]);
	TGraphErrors* gre = new TGraphErrors(xe.size(), &xe[0], &e[0], &xee[0], &ee[0]);
	gra->SetMarkerStyle(22);
	gra->SetMarkerColor(2);
	gre->SetMarkerStyle(22);

	gre->Draw("AP");
	//gre->GetYaxis()->SetRangeUser(0.2, 1.5);
	gra->Draw("Psame");
}

void go(){ 
	CalibrExp ab;
	int cicleV = 1;
	int chooseV = 0;
	while (cicleV == 1) {
		cout << "11 - select BhaBha from raw .col stream ; 13 - find good runs for this dataset (eff and event number);" << endl;
		cout << "14 - get z profiles; 15- draw z profiles; 16 - get phi prof; 17 - fit phi profiles; 18 - average borders values; " << endl;
		cout << "19 - generate map (191 - without shifter shifts); 20 - get list of modeling files (jobs) / h2root command" << endl;
		cout << "1  - launch full block one (as above) - maybe I will ask someone and make it to be a job" << endl << endl;
		cout << "21 - apply the same BhaBha conditons for modeling; 22 - get phi profiles for mod; 23 - fit mod; 24 - get z profiles mod; 25 - draw z profiles mod;" << endl;
		cout << "2  - launch full modeling analisys block" << endl;
		cout << "31 - getTspectre; 32- draw it; 33 - one iteration of amptime; 34 - get profile coefficients; 35 - compare EffSp; 36 - compare zAmpSp; 37 - borders and profiles comparison; 38 - phi shifts; 39 - amp distr" << endl;

		cin >> chooseV;
		cout << endl;
		switch (chooseV)
		{
		case 0:
			cicleV = 0;
			break;
		case 1: {
			ab.shwidth = 1.5;
			ab.copy1();
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
		case 13:
			ab.findruns(ab.entriesInGoodRun, ab.effInGoodRun);
			break;
		case 14:
			ab.zraspr();
			break;
		case 15:
			ab.linesforz("zprofiles.root", "allsensorsl", "allsensors0r");
			break;
		case 16:
			ab.raspr();
			break;
		case 17: {
			ab.shwidth = 1.5;
			ab.testforfit("profiles.root", "zobl");
			break;
		}
		case 18: {
			ab.shwidth = 1.5;
			ab.bordersplot();
			break;
		}
		case 19:
			ab.modelfile(true);
			break;
		case 191:
			ab.modelfile(false);
			break;
		case 20: {
			int basedirInd = 0;
			string year = "2019";
			cout << "input basedir ind, 1 - from basedirmod, other - from /online/simulation:	";
			cin >> basedirInd;
			cout << "input where to store (basedirmod):	";
			cin >> year;
			ab.getlist(basedirInd, year);
			break;
		}
		case 2: {
			ab.shwidth = 1.5;
			ab.copymod();
			ab.rasprmod(false,false);
			ab.testforfit("profilesmmod.root","mean");
			ab.zrasprmod(false); 
			ab.linesforz("zprofilesmod.root","allsensors","allsensors0");
			break;
		}
		case 21:
			ab.copymod();
			break;
		case 22:
			ab.rasprmod(false,false);
			break;
		case 221:
			ab.rasprmod(true,false);
			break;
		case 23: {
			ab.shwidth = 1.5;
			ab.testforfit("profilesmmod.root", "mean");
			break;
		}
		case 24: {
			ab.zrasprmod(false);
			ab.linesforz("zprofilesmod.root", "allsensors", "allsensors0");
			break;
		}
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
		case 33: {
			ab.ampltime();
			ab.drawamptime();
			break;
		}
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
			ab.compareAmpSpectrG();
			break;
		case 38:
			ab.comparePhiShifts();
			break;
		case 39:
			ab.compareAmpDistributions();
			break;
		case 499:
			checkTrans();
			break;
		default:
			break;
		}
	}

/*
pscp D:\workProgramms\cherenkovCalibration\select_ach1_2019.cpp kladov@sndxt1.inp.nsk.su:/sweet/home/kladov/select_ach1_2019.cpp
12081998Vkl

*/
	// cp select_ach1_2019.cpp /online/users2/kladov/maindir/select_ach1_2019.cpp
	//pscp kladov@sndhw3.inp.nsk.su:/online/users2/kladov/maindir/select_ach1_2019.cpp D:\workProgramms\cherenkovCalibration\select_ach1_2019.cpp
}

