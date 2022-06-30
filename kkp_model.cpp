#include <fstream> 
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <map>
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


double temparrayE2011[41] = {510.5, 525, 537.5, 550, 562.5, 575, 587.5, 600, 612.5, 625, 637.5, 650, 662.5, 675, 687.5, 700, 712.5, 725, 737.5, 750, 762.5, 775, 787.5, 800, 812.5, 825, 837.5, 850, 862.5, 875, 887.5, 900, 912.5, 925, 935, 945, 950, 962.5, 975, 987.5, 1000};
//double temparrayE2011[25] = {712.5, 725, 737.5, 750, 762.5, 775, 787.5, 800, 812.5, 825, 837.5, 850, 862.5, 875, 887.5, 900, 912.5, 925, 935, 945, 950, 962.5, 975, 987.5, 1000};
double temparraydE2011[41] = { 454.382, 478.319, 541.344, 524.745, 529.335, 522.558, 524.538, 534.044, 531.015, 532.518, 534.64, 505.007, 527.997, 505.063, 526.018, 496.186, 552.671, 519.782, 534.653, 1074.09, 527.142, 540.323, 541.752, 525.711, 538.466, 512.586, 537.793, 526.935, 547.567, 522.779, 542.942, 492.381, 558.21, 526.207, 543.09, 534.995, 474.281, 494.361, 481.722, 500.974, 480.053 };
//double temparraydE2011[41] = {667.491, 676.499, 673.06, 667.055, 672.703, 669.067, 668.95, 671.122, 669.824, 668.622, 674.237, 658.593, 672.157, 659.678, 667.867, 653.578, 681.1, 667.45, 674.208, 662.667, 668.973, 671.764, 676.031, 673.469, 672.878, 664.36, 675.005, 666.756, 676.974, 670.569, 679.697, 671.171, 681.507, 676.257, 683.036, 676.602, 669.235, 669.884, 668.304, 676.59, 668.943};
//double temparraydE2011[25] = {681.1, 667.45, 674.208, 662.667, 668.973, 671.764, 676.031, 673.469, 672.878, 664.36, 675.005, 666.756, 676.974, 670.569, 679.697, 671.171, 681.507, 676.257, 683.036, 676.602, 669.235, 669.884, 668.304, 676.59, 668.943};

//double temparrayE2011[25] = {712.5, 725, 737.5, 750, 762.5, 775, 787.5, 800, 812.5, 825, 837.5, 850, 862.5, 875, 887.5, 900, 912.5, 925, 935, 945, 950, 962.5, 975, 987.5, 1000};
//double temparraydE2011[25] = {628.104, 612.092, 621.343, 614.365, 616.937, 619.699, 623.654, 618.042, 620.558, 606.827, 621.242, 612.074, 622.695, 612.71, 624.236, 628.149, 625.437, 615.705, 624.96, 617.568, 597.489, 606.168, 599.665, 603.833, 598.85};
//double temparraydE2011[25] = {658.211, 643, 652.637, 644.916, 647.679, 650.025, 653.886, 649.248, 651.457, 637.787, 651.793, 641.716, 652.657, 642.615, 654.763, 664.697, 654.299, 646.084, 656.114, 648.224, 633.268, 639.662, 634.477, 637.33, 633.334};

double temparrayE2012[25] = {505.3, 508, 509, 509.16, 510, 511, 512, 512.299, 513, 640, 680, 720, 760, 800, 840, 860, 880, 900, 920, 936, 950, 960, 970, 980, 990};
double temparraydE2012[25] = { 271.106, 256.058, 246.169, 223.335, 209.316, 278.774, 199.951, 553.503, 194.777, 499.562, 546.499, 534.244, 537.816, 537.936, 539.296, 532.78, 517.295, 531.966, 506.613, 526.644, 523.712, 521.012, 509.934, 491.007, 457.576 };
//double temparraydE2012[25] = {663.131, 661.548, 661.114, 661.348, 660.537, 657.828, 657.256, 657.551, 656.726, 659.901, 674.224, 675.999, 679.889, 680.174, 679.888, 679.717, 680.06, 681.697, 682.396, 682.911, 682.312, 687.014, 682.002, 678.519, 678.87};
		
double temparrayE2017[25] = {641, 650, 675, 700, 725, 750, 775, 800, 825, 840, 937.5, 938.299, 938.899, 939.6, 940.2, 940.799, 942, 950, 960, 970, 980, 990, 1000, 1003, 1003.5};
double temparraydE2017[25] = { 479.5, 469.341, 567.696, 645.29, 520.767, 565.233, 567.633, 606.596, 597.607, 454.387, 517.799, 429.772, 442.747, 428.076, 447.762, 585.478, 451.692, 504.218, 405.881, 388.087, 476.957, 483.759, 451.058, 489.836, 504.781 };
//double temparraydE2017[25] = {678.093, 678.551, 685.889, 690.533, 685.975, 687.116, 688.196, 689.187, 690.504, 683.958, 689.022, 683.551, 685.126, 683.597, 687.259, 693.601, 687.821, 689.2, 684.379, 682.569, 687.17, 689.363, 686.795, 626.044, 669.004};

//double energiesMod2011[24] = {712.5, 725, 737.5, 750, 762.5, 775, 787.5, 800, 812.5, 825, 837.5, 850, 862.5, 875, 887.5, 900, 912.5, 925, 935, 945, 950, 962.5, 975, 987.5};
double energiesMod2011[34] = { 600, 612.5, 625, 637.5, 650, 662.5, 675, 687.5, 700, 712.5, 725, 737.5, 750, 762.5, 775, 787.5, 800, 812.5, 825, 837.5, 850, 862.5, 875, 887.5, 900, 912.5, 925, 935, 945, 950, 962.5, 975, 987.5 , 1000 };
double efficienciesMod2011[34] = {0.0606788, 0.105458, 0.129497, 0.157937, 0.168477, 0.172677, 0.177356, 0.168957, 0.148457, 0.135577, 0.136217, 0.139097, 0.138489, 0.127017, 0.118818, 0.116038, 0.106218, 0.105838, 0.102558, 0.101218, 0.0952781, 0.0938381, 0.0887582, 0.0814184, 0.0754585, 0.0739385, 0.0680186, 0.0670987, 0.0641987, 0.0593588, 0.0591588, 0.0555389, 0.0572189 , 0.0540989 };
//double efficienciesMod2011[24] = {0.124595, 0.116131, 0.112509, 0.102764, 0.0982206, 0.0951253, 0.0899502, 0.0836739, 0.0862153, 0.0787446, 0.078604, 0.0773094, 0.0783167, 0.0760487, 0.0705346, 0.0717419, 0.0723049, 0.0695546, 0.0701623, 0.0700201, 0.0677261, 0.0651169, 0.0628661, 0.0631473};

double energiesMod2012[16] = {640, 680, 720, 760, 800, 840, 860, 880, 900, 920, 936, 950, 960, 970, 980, 990};
double efficienciesMod2012[16] = {0.130117, 0.159937, 0.150717, 0.111978, 0.0890582, 0.0874383, 0.0833983, 0.0711786, 0.0732185, 0.0685786, 0.0625187, 0.0628987, 0.0564589, 0.0588788, 0.0589788, 0.0579988};

double energiesMod2017[24] = {641, 650, 675, 700, 725, 750, 775, 800, 825, 840, 936, 937.5, 938.899, 938.299, 939.6, 940.2, 940.799, 942, 950, 960, 970, 980, 990, 1000};
double efficienciesMod2017[24] = {0.0591988, 0.0709986, 0.100398, 0.0965981, 0.136377, 0.131317, 0.121538, 0.115758, 0.113598, 0.112478, 0.0817184, 0.0800392, 0.0807384, 0.0785984, 0.0810784, 0.0793584, 0.0818784, 0.0800784, 0.0806184, 0.0773185, 0.0774985, 0.0759185, 0.0752185, 0.0752185, };

struct Entryvar {  //parameters of entry
	float phi[40], phis[40], ach[40], achr[40], theta[40], z0[40], schr[40], tch[40], tchr[40], energy[40], beam, eton, x0[40], y0[40], d0[40], dExs[40], dExnC[40], dExn[40];  //kinematic, counters parameters float
	float d2phi[40], dphirho[40], d2rho[40], d2z0[40], d2cosTh[40], dz0cosTh[40], Dtheta[40], Dphi[40], energyerr[40];  //errors afrer reconstruction
	int nch, run, nc, nn, cosm, eventtime, act, region[40];     //parameters int
	float amplitude[40], amplitudeE[40]; //cherenkov general amplitude
	float x2ikf1, ppkf1[11], thetakf1[11], phikf1[11];  //kinfit
	int npkf1[11], ipkf1[11];   //kinfit 
	float x2ikf2, ppkf2[11], thetakf2[11], phikf2[11];  //kinfit
	int npkf2[11], ipkf2[11];   //kinfit 
	float mctheta[40], mcphi[40], mcpx[40], mcpy[40], mcpz[40], mce[40];	//montecarlo modeling
	int mcpdg[40], nmc;	//montecarlo modeling
};
//Entryvar sEntr;

//"(x2ikf1 < 50 && ((((ppkf1[0] < 350 && dExs[ipkf1[0]-1] > 600. * (ppkf1[0]*ppkf1[0] + 320.*320.)/(ppkf1[0]*ppkf1[0]) * (log((ppkf1[0]*ppkf1[0])/(320.*320.)) - (ppkf1[0]*ppkf1[0])/(ppkf1[0]*ppkf1[0] + 320.*320.) + 6.5)) || ((ppkf1[0] >= 350) && (amplitude[ipkf1[0]-1] > 500 || amplitude[ipkf1[0]-1] < 0.22) && region[ipkf1[0]-1] == 1)) && ((ppkf1[1] < 350 && dExs[ipkf1[1]-1] < 600. * (ppkf1[1]*ppkf1[1] + 320.*320.)/(ppkf1[1]*ppkf1[1]) * (log((ppkf1[1]*ppkf1[1])/(320.*320.)) - (ppkf1[1]*ppkf1[1])/(ppkf1[1]*ppkf1[1] + 320.*320.) + 6.5)) || ((ppkf1[1] >= 350) && (amplitude[ipkf1[1]-1] < 500 && amplitude[ipkf1[1]-1] > 0.22) && region[ipkf1[1]-1] == 1))) ||  (((ppkf1[0] < 350 && dExs[ipkf1[0]-1] < 600. * (ppkf1[0]*ppkf1[0] + 320.*320.)/(ppkf1[0]*ppkf1[0]) * (log((ppkf1[0]*ppkf1[0])/(320.*320.)) - (ppkf1[0]*ppkf1[0])/(ppkf1[0]*ppkf1[0] + 320.*320.) + 6.5)) || ((ppkf1[0] >= 350) && (amplitude[ipkf1[0]-1] < 500 && amplitude[ipkf1[0]-1] > 0.22) && region[ipkf1[0]-1] == 1)) && ((ppkf1[1] < 350 && dExs[ipkf1[1]-1] > 600. * (ppkf1[1]*ppkf1[1] + 320.*320.)/(ppkf1[1]*ppkf1[1]) * (log((ppkf1[1]*ppkf1[1])/(320.*320.)) - (ppkf1[1]*ppkf1[1])/(ppkf1[1]*ppkf1[1] + 320.*320.) + 6.5)) || ((ppkf1[1] >= 350) && (amplitude[ipkf1[1]-1] > 500 || amplitude[ipkf1[1]-1] < 0.22) && region[ipkf1[1]-1] == 1)))) && eton>0.3 && eton < 0.8 && (d0[ipkf1[0] - 1] > -0.5) && (d0[ipkf1[0] - 1] < 0.5) && (z0[ipkf1[0] - 1] > -10) && (z0[ipkf1[0] - 1] < 10) && (d0[ipkf1[1] - 1] > -0.5) && (d0[ipkf1[1] - 1] < 0.5) && (z0[ipkf1[1] - 1] > -10) && (z0[ipkf1[1] - 1] < 10) && (z0[ipkf1[0] - 1] - z0[ipkf1[1] - 1]) < 1.5 && (z0[ipkf1[0] - 1] - z0[ipkf1[1] - 1]) > -1.5 && (theta[ipkf1[0] - 1] > 30 * 3.14159 / 180. && theta[ipkf1[0] - 1] < 150 * 3.14159 / 180.) && (theta[ipkf1[1] - 1] > 30 * 3.14159 / 180. && theta[ipkf1[1] - 1] < 150 * 3.14159 / 180.))"

using namespace std;
double PI = 3.14159;
double scalrad = 180. / PI;
double normcoeff = 671.395;

//_________section of ee analisys for luminosity



class luminCalculation {
public:
	Entryvar sEntr;
	float phi[40], phis[40], ach[40], achr[40], theta[40], z0[40], schr[40], tch[40], tchr[40], energy[40], beam, eton, x0[40], y0[40], d0[40], dExs[40], dExnC[40], dExn[40];  //kinematic, counters parameters float
	float d2phi[40], dphirho[40], d2rho[40], d2z0[40], d2cosTh[40], dz0cosTh[40], Dtheta[40], Dphi[40], energyerr[40];  //errors afrer reconstruction
	int nch, run, nc, nn, cosm, eventtime, act, cact , region[40], col;     //parameters int
	float amplitude[40], amplitudeE[40]; //cherenkov general amplitude
	float x2ikf1, ppkf1[11], thetakf1[11], phikf1[11];  //kinfit
	int npkf1[11], ipkf1[11];   //kinfit 
	float x2ikf2, ppkf2[11], thetakf2[11], phikf2[11];  //kinfit
	int npkf2[11], ipkf2[11];   //kinfit 
	float mctheta[40], mcphi[40], mcpx[40], mcpy[40], mcpz[40], mce[40];	//montecarlo modeling
	int mcpdg[40], nmc;	//montecarlo modeling
	float dEx1C[10], dEx2C[10], dEx3C[10], dEx4C[10], dEx5C[10], dEx6C[10], dEx7C[10], dEx8C[10], dEx9C[10];

	double thetacut;
	double energycutl;
	double energycutr;
	double dphicut;
	double dthetacut;
	double d0cut;
	double z0cut;
	TProfile* dedxm;
	TProfile* dedxe;

	vector<double> energyPoints;
	vector<int> countModEventsEE;

	vector<double> meandExn;
	vector<double> meandExnErr;
	vector<double> meandExnC;
	vector<double> meandExnCErr;
	vector<double> dExnpopr;
	vector<double> dExnpoprErr;

	vector<double> effem;
	vector<double> efftm;
	vector<double> effstm;
	vector<double> effpm;
	vector<double> effdm;
	vector<double> effzm;
	vector<double> effcm;
	vector<double> deffem;
	vector<double> defftm;
	vector<double> deffstm;
	vector<double> deffpm;
	vector<double> deffdm;
	vector<double> deffzm;
	vector<double> deffcm;

	vector<double> effee;
	vector<double> effte;
	vector<double> effste;
	vector<double> effpe;
	vector<double> effde;
	vector<double> effze;
	vector<double> effce;
	vector<double> deffee;
	vector<double> deffte;
	vector<double> deffste;
	vector<double> deffpe;
	vector<double> deffde;
	vector<double> deffze;
	vector<double> deffce;

	vector<double> popre;
	vector<double> poprt;
	vector<double> poprst;
	vector<double> poprp;
	vector<double> poprd;
	vector<double> poprz;
	vector<double> poprc;
	vector<double> dpopre;
	vector<double> dpoprt;
	vector<double> dpoprst;
	vector<double> dpoprp;
	vector<double> dpoprd;
	vector<double> dpoprz;
	vector<double> dpoprc;


	vector< vector<int> > countSelectedEventsEE[8];
	vector<double> energyMod;
	vector<double> energyExp;
	vector<TH1*> dEdxMod1;
	vector<TH1*> dEdxExp1;
	vector< vector<int> > entriesExp[8];
	vector< vector<double> > efficiency[8];
	vector< vector<double> > efficiencyErr[8];
	vector<double> energyErr;
	vector<double> luminosity;
	vector<double> luminosityErr;

	vector< vector<double> > lumforsystematic[8];
	vector< vector<double> > lumforsystematicErr[8];

	vector<double> baseparam;
	double varsizes[8];
	vector<double> varwidths;
	vector< vector<double> > parameters[8];

	luminCalculation(){
		col = 0;
		thetacut = 45.;
		energycutl = 0.25; //!!! 0.25
		energycutr = 0.25; //0.2 - strong cuts for systematic analysis
		dphicut = 0.16;
		dthetacut = 0.25;
		d0cut = 1.0;   //wide cuts to ignore discrepancies in mod with exp
		z0cut = 3.0;
		//dedxm = new TProfile("dedxm","mod dExn vs energy", 100,400,1400);
		//dedxe = new TProfile("dedxe","exp dExnC vs energy", 100,400,1400);
		vector<double> vtd;
		vector<int> vti;
		for(size_t i = 0; i < 7; i++)
			varsizes[i] = 21;
		varsizes[7] = 2;
		varsizes[3] = 21;
		baseparam.push_back(energycutl);
		baseparam.push_back(energycutr);
		baseparam.push_back(dthetacut);
		baseparam.push_back(thetacut);
		baseparam.push_back(dphicut);
		baseparam.push_back(d0cut); 
		baseparam.push_back(z0cut);
		baseparam.push_back(0);

		varwidths.push_back(0.75);
		varwidths.push_back(0.75);
		varwidths.push_back(0.4);
		varwidths.push_back(5.);
		varwidths.push_back(0.4);
		varwidths.push_back(2.0); 
		varwidths.push_back(3.0);
		varwidths.push_back(1);

		for(size_t i = 0; i < 8; i++){
			for(size_t j = 0; j<varsizes[i]; j++){
				countSelectedEventsEE[i].push_back(vti);
				entriesExp[i].push_back(vti);
				efficiency[i].push_back(vtd);
				efficiencyErr[i].push_back(vtd);
				lumforsystematic[i].push_back(vtd);
				lumforsystematicErr[i].push_back(vtd);
				
				parameters[i].push_back(baseparam);
			}
		}
		for(size_t i = 0; i < 8; i++){
			for(size_t j = 1; j<varsizes[i]; j++){
				if(i==(size_t)3){
					if (j<5)
						parameters[i][j][i] += -2.5*((int)j-5);
					else
						parameters[i][j][i] += -2.5*((int)j-4);
				}
				else
					parameters[i][j][i]+=varwidths[i]*(((double)j)/(varsizes[i]-1.));
				//if(i!=0)
				//	parameters[i][j][0]-=0.15;
				//if(i!=1)
				//	parameters[i][j][1]-=0.15;
			}
		}
	}
	//selection criteria
	bool dtheta(double cut) {
		return (fabs(theta[0] + theta[1] - 3.14159) <= cut);
	}
	bool dphi(double cut) {
		return ((fabs(fabs(phi[0] - phi[1]) - 3.14159) <= cut));
	}
	bool thetag(double cut) {
		return (theta[0] >= (90.-cut)/scalrad && theta[0] <= (90.+cut)/scalrad) && (theta[1] >= (90.-cut)/scalrad && theta[1] <= (90.+cut)/scalrad);
		//return (theta[0] >= 0.7 && theta[0] <= 2.44) && (theta[1] >= 0.7 && theta[1] <= 2.44);
	}
	bool energy0cut(double cutl, double cutr,int indexp) {
		return (energy[indexp] / beam >= 1-cutl && energy[indexp] / beam <= 1+cutl);
	}
	bool energy1cut(double cutl, double cutr,int indexp) {
		return (energy[indexp] / beam >= 1-cutr && energy[indexp] / beam <= 1+cutr);
	}
	bool dd0cut(double cut) {
		return (fabs(d0[0] - d0[1]) <= cut);
	}
	bool dz0cut(double cut) {
		return (fabs(z0[0] - z0[1]) <= cut);
	}
	bool partnumb(){
		return (nc >= 2) && (nn <= 3) && (nc <= 3);
	}
	bool actb(int cutyn){
		//return ((bool)cutyn) || (cact == 0);
		return true;
	}

	bool isItEE(double* cutvector) {
		int indET = rand() % 2;
		int indET2 = abs(indET-1);
		return (col == 1) && partnumb() && energy0cut(cutvector[0],cutvector[1],indET) && energy1cut(cutvector[0],cutvector[1],indET2) && actb((int)cutvector[7]) && dtheta(cutvector[2]) && thetag(cutvector[3]) && dphi(cutvector[4]) && dd0cut(cutvector[5]) && dz0cut(cutvector[6]);
	}
	bool isItEEEn0() {
		int indET = 0;//rand() % 2;
		int indET2 = abs(indET-1);
		return (col == 1) && partnumb() && energy1cut(baseparam[0],baseparam[1],indET2) && actb((int)baseparam[7]) && dtheta(baseparam[2]) && thetag(baseparam[3]) && dphi(baseparam[4]) && dd0cut(baseparam[5]) && dz0cut(baseparam[6]);
	}
	bool isItEEEn1() {
		int indET = 1;//rand() % 2;
		int indET2 = abs(indET-1);
		return (col == 1) && partnumb() && energy1cut(baseparam[0],baseparam[1],indET2) && actb((int)baseparam[7]) && dtheta(baseparam[2]) && thetag(baseparam[3]) && dphi(baseparam[4]) && dd0cut(baseparam[5]) && dz0cut(baseparam[6]);
	}
	bool isItEETheta() {
		int indET = rand() % 2;
		int indET2 = abs(indET-1);
		return (col == 1) && partnumb() &&  energy0cut(baseparam[0],baseparam[1],indET) && energy1cut(baseparam[0],baseparam[1],indET2)  && actb((int)baseparam[7]) && dtheta(baseparam[2]) && dphi(baseparam[4]) && dd0cut(baseparam[5]) && dz0cut(baseparam[6]);
	}
	bool isItEESTheta() {
		int indET = rand() % 2;
		int indET2 = abs(indET-1);
		return (col == 1) && partnumb() &&  energy0cut(baseparam[0],baseparam[1],indET) && energy1cut(baseparam[0],baseparam[1],indET2)  && actb((int)baseparam[7]) && thetag(baseparam[3]) && dphi(baseparam[4]) && dd0cut(baseparam[5]) && dz0cut(baseparam[6]);
	}
	bool isItEEPhi() {
		int indET = rand() % 2;
		int indET2 = abs(indET-1);
		return (col == 1) && partnumb() &&  energy0cut(baseparam[0],baseparam[1],indET) && energy1cut(baseparam[0],baseparam[1],indET2)  && actb((int)baseparam[7]) && dtheta(baseparam[2]) && thetag(baseparam[3]) && dd0cut(baseparam[5]) && dz0cut(baseparam[6]);
	}
	bool isItEEDd0() {
		int indET = rand() % 2;
		int indET2 = abs(indET-1);
		return (col == 1) && partnumb() &&  energy0cut(baseparam[0],baseparam[1],indET) && energy1cut(baseparam[0],baseparam[1],indET2)  && actb((int)baseparam[7]) && dtheta(baseparam[2]) && thetag(baseparam[3]) && dphi(baseparam[4]) && dz0cut(baseparam[6]);
	}
	bool isItEEDz0() {
		int indET = rand() % 2;
		int indET2 = abs(indET-1);
		return (col == 1) && partnumb() &&  energy0cut(baseparam[0],baseparam[1],indET) && energy1cut(baseparam[0],baseparam[1],indET2)  && actb((int)baseparam[7]) && dtheta(baseparam[2]) && thetag(baseparam[3]) && dphi(baseparam[4]) && dd0cut(baseparam[5]);
	}
	bool isItEEcact() {
		int indET = rand() % 2;
		int indET2 = abs(indET-1);
		return (col == 1) && partnumb() &&  energy0cut(baseparam[0],baseparam[1],indET) && energy1cut(baseparam[0],baseparam[1],indET2)  && dtheta(baseparam[2]) && thetag(baseparam[3]) && dphi(baseparam[4]) && dd0cut(baseparam[5]) && dz0cut(baseparam[6]);
	}

	void effEeEnergy(char* name, char* add, char* file);
	void pointComparisonAndSort();
	void luminosityCalc();
	void differenceForLum();
	void drawdEdxspectrum();
};

//deriving all parameters from mod and exp (change files and t1<->h1)
void luminCalculation::effEeEnergy(char *name, char *add, char *file) {
	TChain chain(name);
	chain.Add(add);
	//chain.Add("/work/users/konctbel/calibs/R007-001/output/ntuples/MHAD2017/col/*col.root");
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	TString b = add;
	{
		chain.SetBranchAddress("beam", &beam);
		chain.SetBranchAddress("energy", &energy);
		chain.SetBranchAddress("nc", &nc);
		chain.SetBranchAddress("nn", &nn);
		chain.SetBranchAddress("z0", &z0);
		chain.SetBranchAddress("d0", &d0);
		chain.SetBranchAddress("phi", &phi);
		chain.SetBranchAddress("theta", &theta);
		chain.SetBranchAddress("col", &col);
		chain.SetBranchAddress("cosm", &cosm);
		chain.SetBranchAddress("act", &act);
		chain.SetBranchAddress("cact", &cact);
		chain.SetBranchAddress("dExn", &dExn);
		chain.SetBranchAddress("dExnC", &dExnC);
		if (b.Contains("MHAD")){
			chain.SetBranchAddress("dEx1C", &dEx1C);
			chain.SetBranchAddress("dEx2C", &dEx2C);
			chain.SetBranchAddress("dEx3C", &dEx3C);
			chain.SetBranchAddress("dEx4C", &dEx4C);
			chain.SetBranchAddress("dEx5C", &dEx5C);
			chain.SetBranchAddress("dEx6C", &dEx6C);
			chain.SetBranchAddress("dEx7C", &dEx7C);
			chain.SetBranchAddress("dEx8C", &dEx8C);
			chain.SetBranchAddress("dEx9C", &dEx9C);
		}

	}
	map<double, TH1*> hE0m;
	map<double, TH1*> hE1m;
	map<double, TH1*> hThetam;
	map<double, TH1*> hSThetam;
	map<double, TH1*> hPhim;
	map<double, TH1*> hDd0m;
	map<double, TH1*> hDz0m;
	TProfile* dedxT = new TProfile("dedxm","mod dExn vs energy", 100,400,1400);

	vector<double> alle;
	vector<double> allt;
	vector<double> allst;
	vector<double> allp;
	vector<double> alld;
	vector<double> allz;
	vector<double> allc;
	vector<double> passe;
	vector<double> passt;
	vector<double> passst;
	vector<double> passp;
	vector<double> passd;
	vector<double> passz;
	vector<double> passc;
	
	vector<int> ndExncount;

	vector<TH1*> hE0;// = new TH1F("hE0", "energy0/beam", 100, 0, 2);
	vector<TH1*> hE1;// = new TH1F("hE1", "energy1/beam", 100, 0, 2);
	vector<TH1*> hTheta;// = new TH1F("hTheta", "theta", 100, 0, 3.1416);
	vector<TH1*> hSTheta;// = new TH1F("hSTheta", "theta1 + theta2 - Pi", 100, -0.5, 0.5);
	vector<TH1*> hPhi;// = new TH1F("hPhi", "phi1 - phi2 - pi", 100, -0.5, 0.5);
	vector<TH1*> hDd0;// = new TH1F("hDd0", "d01 - d02", 100, -2.5, 2.5);
	vector<TH1*> hDz0;// = new TH1F("hDz0", "z01 - z02", 100, -15, 15);
	vector<TH1*> hdExn;// = new TH1F("hDz0", "z01 - z02", 100, -15, 15);


	TFile* MyFile = new TFile(file, "RECREATE");
	
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		bool newEnPoint = true;
		int pointIndex = energyPoints.size();
		for (size_t i = 0; i < energyPoints.size(); i++) {
			if (energyPoints[i] == (double)beam) {
				newEnPoint = false;
				pointIndex = i;
			}
		}
		if (newEnPoint) {
			energyPoints.push_back((double)beam);
			countModEventsEE.push_back(1);
			for(size_t i = 0; i < 8; i++){
				for(size_t j = 0; j<varsizes[i]; j++){
					countSelectedEventsEE[i][j].push_back(0);
				}
			}

			hE0.push_back(new TH1F(Form("hE0%f", energyPoints.back()), "energy0/beam", 100, 0, 2));
			hE1.push_back(new TH1F(Form("hE1%f", energyPoints.back()), "energy1/beam", 100, 0, 2));
			hTheta.push_back(new TH1F(Form("hTheta%f", energyPoints.back()), "theta", 100, 0, 3.1416));
			hSTheta.push_back(new TH1F(Form("hSTheta%f", energyPoints.back()), "theta1 + theta2 - Pi", 100, -0.5, 0.5));
			hPhi.push_back(new TH1F(Form("hPhi%f", energyPoints.back()), "phi1 - phi2 - pi", 100, -0.5, 0.5));
			hDd0.push_back(new TH1F(Form("hDd0%f", energyPoints.back()), "d01 - d02", 100, -1.5, 1.5));
			hDz0.push_back(new TH1F(Form("hDz0%f", energyPoints.back()), "z01 - z02", 100, -10, 10));
			hdExn.push_back(new TH1F(Form("hdExn%f", energyPoints.back()), ";dE/dx", 100, 0, 10));

			alle.push_back(0.);
			allt.push_back(0.);
			allst.push_back(0.);
			allp.push_back(0.);
			alld.push_back(0.);
			allz.push_back(0);
			allc.push_back(0);
			passe.push_back(0.);
			passt.push_back(0);
			passst.push_back(0);
			passp.push_back(0);
			passd.push_back(0);
			passz.push_back(0);
			passc.push_back(0);
			ndExncount.push_back(0);
			if (b.Contains("eemod")) {
				meandExn.push_back(0);
				cout << Form("dEdxMod1%d", (int)(energyPoints.back()) + pointIndex) << endl;
				dEdxMod1.push_back(new TH1F(Form("dEdxMod1%d", (int)energyPoints.back()+ pointIndex), ";dE/dx", 1000, 0, 5000));
			}
			if (b.Contains("MHAD")) {
				meandExnC.push_back(0);
				dEdxExp1.push_back(new TH1F(Form("dEdxExp1%d", (int)energyPoints.back() + pointIndex), ";dE/dx", 1000, 0, 10));
			}
		}
		else
			countModEventsEE[pointIndex] += 1;

		int indET = rand() % 2;
		if (indET == 0 && isItEEEn0() || indET == 1 && isItEEEn1()){
			hE0[pointIndex]->Fill(energy[indET] / beam);
			hE1[pointIndex]->Fill(energy[indET] / beam);
			alle[pointIndex]+=1;
			//if (isItEE(&parameters[0][0][0]))
			if (energy0cut(baseparam[0],baseparam[1],indET))
				passe[pointIndex]+=1;
		}

		if (isItEETheta()) {
			hTheta[pointIndex]->Fill(theta[0]);
			hTheta[pointIndex]->Fill(theta[1]);
			allt[pointIndex]+=1;
			if (isItEE(&parameters[0][0][0])){
				passt[pointIndex]+=1;
			}
		}
		if (isItEESTheta()){
			hSTheta[pointIndex]->Fill(theta[0] + theta[1] - 3.14159);
			allst[pointIndex]+=1;
			if (isItEE(&parameters[0][0][0]))
				passst[pointIndex]+=1;
		}
		if (isItEEPhi()){
			hPhi[pointIndex]->Fill(fabs(phi[0] - phi[1]) - 3.14159);
			allp[pointIndex]+=1;
			if (isItEE(&parameters[0][0][0]))
				passp[pointIndex]+=1;
		}
		if (isItEEDd0()){
			hDd0[pointIndex]->Fill(d0[0] - d0[1]);
			alld[pointIndex]+=1;
			if (isItEE(&parameters[0][0][0]))
				passd[pointIndex]+=1;
		}
		if (isItEEDz0()){
			hDz0[pointIndex]->Fill(z0[0] - z0[1]);
			allz[pointIndex]+=1;
			if (isItEE(&parameters[0][0][0]))
				passz[pointIndex]+=1;
		}
		if (isItEEcact()){
			allc[pointIndex]+=1;
			if (isItEE(&parameters[0][0][0]))
				passc[pointIndex]+=1;
		}

		
		if (isItEE(&parameters[0][0][0])) {
			//cout << b << endl;
			if (b.Contains("eemod")){
				int closestPind = -1;
				double mindist = 1000;
				for (size_t i = 0; i < 41; i++) {
			
					if (fabs(beam - temparrayE2011[i]) < mindist) {
						mindist = fabs(beam - temparrayE2011[i]);
						closestPind = i;
					}
				}
				hdExn[pointIndex]->Fill(dExn[0]/temparraydE2011[closestPind]);
				hdExn[pointIndex]->Fill(dExn[1]/temparraydE2011[closestPind]);
				dEdxMod1[pointIndex]->Fill(dExn[0]);
				dEdxMod1[pointIndex]->Fill(dExn[1]);
				//if(dExn[0]/680. < 3.5 && dExn[1]/680. < 3.5){
					//cout << energy[0] << "	" << dExn[0] << endl;
					dedxT->Fill(beam, dExn[0]);
					dedxT->Fill(beam, dExn[1]);
					ndExncount[pointIndex]+=1;
					meandExn[pointIndex] = meandExn[pointIndex] + ((dExn[0]+dExn[1])/2. - meandExn[pointIndex])/(ndExncount[pointIndex]);
				//}
			}
			if (b.Contains("MHAD")){
				double partialded[2][9] = {dEx1C[0], dEx2C[0], dEx3C[0], dEx4C[0], dEx5C[0], dEx6C[0], dEx7C[0], dEx8C[0], dEx9C[0], \
										   dEx1C[1], dEx2C[1], dEx3C[1], dEx4C[1], dEx5C[1], dEx6C[1], dEx7C[1], dEx8C[1], dEx9C[1]};
				double dExnTrue[2] = { 0., 0. };
				double maxValue[2] = { -1, -1 };
				size_t maxIndex[2] = { 0, 0 };
				for (size_t j = 0; j < 2; j++) {
					for (size_t i = 0; i < 9; i++) {
						dExnTrue[j] += partialded[j][i];
						if (partialded[j][i] > maxValue[j]) { maxValue[j] = partialded[j][i];  maxIndex[j] = i; } //find maximum
					}
					dExnTrue[j] = (dExnTrue[j] - partialded[j][maxIndex[j]])/8.; //subtract maximum value
				}

				hdExn[pointIndex]->Fill(dExnTrue[0]);
				hdExn[pointIndex]->Fill(dExnTrue[1]);
				dEdxExp1[pointIndex]->Fill(dExnTrue[0]);
				dEdxExp1[pointIndex]->Fill(dExnTrue[1]);
				//cout << energy[0] << "	" << dExnC[0] << endl;
				//dedxT->Fill(beam, dExnC[0]);
				//dedxT->Fill(beam, dExnC[1]);
				//if(dExnTrue[0] > 0.01 && dExnTrue[1] > 0.01){
					dedxT->Fill(beam, dExnTrue[0]);
					dedxT->Fill(beam, dExnTrue[1]);
					ndExncount[pointIndex]+=1;
					meandExnC[pointIndex] = meandExnC[pointIndex] + ((dExnTrue[0] + dExnTrue[1]) / 2. - meandExnC[pointIndex]) / (ndExncount[pointIndex]);
				//}
			}
		}
		
		for(size_t i = 0; i < 8; i++){
			//cout << parameters[i][0][i] << endl;
			for(size_t j = 0; j<varsizes[i]; j++){
				//cout << col << endl;
				//cout << isItEE(&parameters[0][0][0]) << endl;
				if (isItEE(&parameters[i][j][0])) {
					countSelectedEventsEE[i][j][pointIndex] += 1;
					//if(i==7)
					//	cout << "<<ASDAD" << endl;
				}
			}
		}
		if (e % 1000000 == 0)
			cout << e/1000000 << " M" << endl;
	}

	for (size_t i = 0; i < energyPoints.size(); i++) {
		hE0[i]->Write(Form("hE0_%.3f", energyPoints[i]));
		hE1[i]->Write(Form("hE1_%.3f", energyPoints[i]));
		hTheta[i]->Write(Form("hTheta_%.3f", energyPoints[i]));
		hSTheta[i]->Write(Form("hSTheta_%.3f", energyPoints[i]));
		hPhi[i]->Write(Form("hPhi_%.3f", energyPoints[i]));
		hDd0[i]->Write(Form("hDd0_%.3f", energyPoints[i]));
		hDz0[i]->Write(Form("hDz0_%.3f", energyPoints[i]));
		hdExn[i]->Write(Form("hdExn_%.3f", energyPoints[i]));
	}
	MyFile->Close();

	//erase all elements from vectors in bad points with zero selected events to evade nan and 0s
	vector<int> eraseInd;
	for (size_t i = 0; i < energyPoints.size(); i++) {
		if(countSelectedEventsEE[0][0][i] == 0){
			eraseInd.push_back((int)i);
		}
	}
	cout << eraseInd.size() << endl;
	if(eraseInd.size()>=1){
		for (int i = eraseInd.size()-1; i >-1; i--) {
			cout << i << endl;
			cout << eraseInd[i] << endl;
			energyPoints.erase(energyPoints.begin()+eraseInd[i]);
			alle.erase(alle.begin()+eraseInd[i]);
			allt.erase(allt.begin()+eraseInd[i]);
			allst.erase(allst.begin()+eraseInd[i]);
			allp.erase(allp.begin()+eraseInd[i]);
			alld.erase(alld.begin()+eraseInd[i]);
			allz.erase(allz.begin()+eraseInd[i]);
			allc.erase(allc.begin()+eraseInd[i]);
			passe.erase(passe.begin()+eraseInd[i]);
			passt.erase(passt.begin()+eraseInd[i]);
			passst.erase(passst.begin()+eraseInd[i]);
			passp.erase(passp.begin()+eraseInd[i]);
			passd.erase(passd.begin()+eraseInd[i]);
			passz.erase(passz.begin()+eraseInd[i]);
			passc.erase(passc.begin()+eraseInd[i]);

			countModEventsEE.erase(countModEventsEE.begin()+eraseInd[i]);
			for(size_t k = 0; k < 8; k++){
				for(size_t j = 0; j<varsizes[k]; j++){
					countSelectedEventsEE[k][j].erase(countSelectedEventsEE[k][j].begin()+eraseInd[i]);
				}
			}
		}
	}

	for (size_t i = 0; i < energyPoints.size(); i++) {
		if (b.Contains("eemod")) {
			meandExnErr.push_back(200. / sqrt(ndExncount[i]));   //taken from histogramm for modeling
			cout << energyPoints[i] << endl;
			cout << dEdxMod1[i]->GetEntries() << "	" << i << endl;
			dEdxMod1[i]->Draw();
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
			//cin.get();
		}
		if (b.Contains("MHAD")) {
			meandExnCErr.push_back(0.4 / sqrt(ndExncount[i]));
			dEdxExp1[i]->Draw();
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->Update();
		}
	}

	/*for (size_t i = 0; i < energyPoints.size(); i++) {
                energyMod.push_back(energyPoints[i]);
                efficiency.push_back((double)countSelectedEventsEE[i] / (double)countModEventsEE[i]);
                energyErr.push_back(0.);
                efficiencyErr.push_back(sqrt((double)countSelectedEventsEE[i]) / (double)countModEventsEE[i]);
		cout << energyPoints[i] << ", ";
	}
	cout << endl;
	for (size_t i = 0; i < energyPoints.size(); i++)
		cout << efficiency[i] << ", ";
	for (size_t i = 0; i < energyPoints.size(); i++)
		cout << countSelectedEventsEE[i] << ", ";
	cout << endl;
	TGraphErrors* gr = new TGraphErrors(energyPoints.size(), &energyPoints[0], &efficiency[0], &energyErr[0], &efficiencyErr[0]);
	gr->SetTitle("efficiencyEE/energy");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->Draw("AP");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	cin.get();*/
	cout << "b" << endl;
	if (b.Contains("eemod")){
		dedxm = dedxT;

		for (size_t i = 0; i < energyPoints.size(); i++) {
			if(alle[i]<1)
				cout << "m " << energyPoints[i] << endl;
			if(passe[i]<1)
				cout << "m " << energyPoints[i] << endl;

			effem.push_back((double)passe[i] / (double)alle[i]);
			deffem.push_back(sqrt(effem.back()*(1-effem.back())/(double)alle[i]));
			if(deffem.back()>1){
				cout << passe[i] << endl;
				cout << alle[i] << endl;
			}

			efftm.push_back((double)passt[i] / (double)allt[i]);
			defftm.push_back(sqrt(efftm.back()*(1-efftm.back())/(double)allt[i]));

			effstm.push_back((double)passst[i] / (double)allst[i]);
			deffstm.push_back(sqrt(effstm.back()*(1-effstm.back())/(double)allst[i]));

			effpm.push_back((double)passp[i] / (double)allp[i]);
			deffpm.push_back(sqrt(effpm.back()*(1-effpm.back())/(double)allp[i]));

			effdm.push_back((double)passd[i] / (double)alld[i]);
			deffdm.push_back(sqrt(effdm.back()*(1-effdm.back())/(double)alld[i]));

			effzm.push_back((double)passz[i] / (double)allz[i]);
			deffzm.push_back(sqrt(effzm.back()*(1-effzm.back())/(double)allz[i]));

			effcm.push_back((double)passc[i] / (double)allc[i]);
			deffcm.push_back(sqrt(effcm.back()*(1-effcm.back())/(double)allc[i]));
		}
	}
	if (b.Contains("MHAD")){
		dedxe = dedxT;

		for (size_t i = 0; i < energyPoints.size(); i++) {

			if(alle[i]<1)
				cout << "e " << energyPoints[i] << endl;
			if(passe[i]<1)
				cout << "e " << energyPoints[i] << endl;	

			effee.push_back((double)passe[i] / (double)alle[i]);
			deffee.push_back(sqrt(effee.back()*(1-effee.back())/(double)alle[i]));

			effte.push_back((double)passt[i] / (double)allt[i]);
			deffte.push_back(sqrt(effte.back()*(1-effte.back())/(double)allt[i]));

			effste.push_back((double)passst[i] / (double)allst[i]);
			deffste.push_back(sqrt(effste.back()*(1-effste.back())/(double)allst[i]));

			effpe.push_back((double)passp[i] / (double)allp[i]);
			deffpe.push_back(sqrt(effpe.back()*(1-effpe.back())/(double)allp[i]));

			effde.push_back((double)passd[i] / (double)alld[i]);
			deffde.push_back(sqrt(effde.back()*(1-effde.back())/(double)alld[i]));

			effze.push_back((double)passz[i] / (double)allz[i]);
			deffze.push_back(sqrt(effze.back()*(1-effze.back())/(double)allz[i]));

			effce.push_back((double)passc[i] / (double)allc[i]);
			deffce.push_back(sqrt(effce.back()*(1-effce.back())/(double)allc[i]));
		}
	}
}

//compare images for exp and mod
void compareSelectionCuts() {
	TFile* fm = new TFile("/work/users/kladov/snd2k/R007-001/2019/lumdistrMod.root");
	vector<double> energym;
	vector<double> energye;
	vector<TH1*> hE0m;
	vector<TH1*> hE1m;
	vector<TH1*> hThetam;
	vector<TH1*> hSThetam;
	vector<TH1*> hPhim;
	vector<TH1*> hDd0m;
	vector<TH1*> hDz0m;
	vector<TH1*> hE0e;
	vector<TH1*> hE1e;
	vector<TH1*> hThetae;
	vector<TH1*> hSThetae;
	vector<TH1*> hPhie;
	vector<TH1*> hDd0e;
	vector<TH1*> hDz0e;
	TKey* key;

	TIter next(fm->GetListOfKeys());
	while ((key = (TKey*)next())) {
		TString name(key->GetName());
		if (name.Contains("hE0")) {
			double tempen = 0.;
			sscanf(name.Data(), "hE0%f", &tempen);
			energym.push_back(tempen);
			hE0m.push_back((TH1F*)fm->Get(Form("hE0%f", tempen)));
			hE1m.push_back((TH1F*)fm->Get(Form("hE1%f", tempen)));
			hThetam.push_back((TH1F*)fm->Get(Form("hTheta%f", tempen)));
			hSThetam.push_back((TH1F*)fm->Get(Form("hSTheta%f", tempen)));
			hPhim.push_back((TH1F*)fm->Get(Form("hPhi%f", tempen)));
			hDd0m.push_back((TH1F*)fm->Get(Form("hDd0%f", tempen)));
			hDz0m.push_back((TH1F*)fm->Get(Form("hDz0%f", tempen)));
		}
	}
	fm->Close();
	TFile* fe = new TFile("/work/users/kladov/snd2k/R007-001/2017/lumdistrExp.root");
	TIter nexte(fe->GetListOfKeys());
	while ((key = (TKey*)nexte())) {
		TString name(key->GetName());
		if (name.Contains("hE0")) {
			double tempen = 0.;
			sscanf(name.Data(), "hE0%f", &tempen);
			energye.push_back(tempen);
			hE0e.push_back((TH1F*)fe->Get(Form("hE0%f", tempen)));
			hE1e.push_back((TH1F*)fe->Get(Form("hE1%f", tempen)));
			hThetae.push_back((TH1F*)fe->Get(Form("hTheta%f", tempen)));
			hSThetae.push_back((TH1F*)fe->Get(Form("hSTheta%f", tempen)));
			hPhie.push_back((TH1F*)fe->Get(Form("hPhi%f", tempen)));
			hDd0e.push_back((TH1F*)fe->Get(Form("hDd0%f", tempen)));
			hDz0e.push_back((TH1F*)fe->Get(Form("hDz0%f", tempen)));
		}
	}
	fe->Close();

	for (size_t i = 0; i < energye.size(); i++) {
		size_t closestPInd;
		double mindist = 1000;
		for (size_t j = 0; j < energym.size(); j++) {
			if (fabs(energym[j] - energye[i]) < mindist) {
				mindist = fabs(energym[j] - energye[i]);
				closestPInd = j;
			}
		}
		hE0e[i]->DrawNormalized("",1);
		hE0m[closestPInd]->SetLineColor(2);
		hE0m[closestPInd]->DrawNormalized("same", 1);

		//hE0e[i]->GetYaxis()->SetLogY();
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->cd(1);
		gPad->SetLogy();
		c->SetFillColor(0);
		c->SetFrameFillColor(0);
		c->SetBorderMode(0);
		c->SetFrameBorderMode(0);
		c->SetGrid();
		c->Update();
		cin.get();

	}

}

//calculated cross section in mathematica
double crossSectNb(double en) {
	double e3 = en / 1000.;
	//return 470.251 / e3 / e3 * (1 + 0.0122513 * e3) * 1.017 / 1.145 / 1.0097; //it is from fit of modeling files info //t35
	
	//return 470.251 / e3 / e3 * (1 + 0.0122513 * e3);	//for t35 from fit
	//return 687.358 / e3 / e3 * (1 + 0.0148081 * e3);	//for t30 from fit
	
	//return 470.251 / e3 / e3 * (1 + 0.0122513 * e3) * 0.998;	//2011
	return 470.251 / e3 / e3 * (1 + 0.0122513 * e3) * 0.994;	//2012
	//return 687.358 / e3 / e3 * (1 + 0.0148081 * e3) * 1.0273;	//2017

	//return 689.0899 / e3 / e3 * (1 + 0.0122513 * e3); //it is for t30, but not excactly (need to fit)
	//return 465.251 / e3 / e3 * (1 + 0.0132513 * e3);
	//return 610.251 / e3 / e3 * (1 + 0.0522513 * e3);
}


//calculate luminosity using copied arrays from effEeEnergy output
void luminCalculation::pointComparisonAndSort() {
	/*double energyExp1[40] = { 850.2, 851, 850.222, 850, 860, 870, 880, 882, 890, 891.173, 901, 910, 921.695, 930, 936, 942, 950, 960, 970, 980, 990, \
		1003, 1000, 1003.5, 938.899, 939.6, 940.2, 938.299, 937.5, 940.799, 840, 825, 800, 775, 750, 725, 700, 675, 650, 641 };
	int entriesExp1[40] = { 89761, 18906, 200481, 360017, 554728, 686451, 463121, 227999, 122528, 778226, 662033, 689867, 737398, 702633, 556285, 625710, 463656, \
		561717, 702032, 448714, 444683, 246, 1275290, 672, 815347, 905017, 1008181, 1304764, 1006829, 829104, 451308, 652348, 483015, 435068, 667688, 914909, 801415, 1197150, 1005018, 3370375 };
	double energyMod1[33] = { 1000, 641, 650, 675, 700, 725, 750, 775, 800, 825, 840, 850, 860, 870, 880, 890, 900, 910, \
		920, 930, 936, 937.5, 938.899, 938.299, 939.6, 940.2, 940.799, 942, 950, 960, 970, 980, 990 };
	double efficiency1[33] = { 0.519377, 0.510032, 0.50329, 0.522935, 0.503535, 0.505585, 0.520161, 0.520081, 0.518721, 0.52015, 0.51956, 0.53901, 0.53827, \
		0.534149, 0.531942, 0.538274, 0.540115, 0.543004, 0.538042, 0.539975, 0.533841, 0.519006, 0.524675, 0.509678, 0.51885, 0.518145, 0.519235, 0.51939, 0.52474, 0.528581, 0.52471, 0.523305, 0.525741 };
	vector<double> energyExp;
	vector<int> entriesExp;
	vector<double> energyMod;
	vector<double> efficiency;
	for (int i = 0; i < 40; i++) {
		energyExp.push_back(energyExp1[i]);
		entriesExp.push_back(entriesExp1[i]);
	}
	for (int i = 0; i < 33; i++) {
		energyMod.push_back(energyMod1[i]);
		efficiency.push_back(efficiency1[i]);
	}*/
	ofstream ofile("/work/users/kladov/snd2k/R007-001/2011/luminosity.dat");

	//calculating of lum and popr (syst)
	for (size_t i = 0; i < energyExp.size(); i++) {
		size_t closestPInd;
		double mindist = 1000;
		for (size_t j = 0; j < energyMod.size(); j++) {
			if (fabs(energyMod[j] - energyExp[i]) < mindist) {
				mindist = fabs(energyMod[j] - energyExp[i]);
				closestPInd = j;
			}
		}
		luminosity.push_back((double)entriesExp[0][0][i]/efficiency[0][0][closestPInd]/crossSectNb(energyExp[i]));
		luminosityErr.push_back(/*luminosity[i]/sqrt((double)countModEventsEE[i]) + */luminosity[i] * (efficiencyErr[0][0][closestPInd]/ efficiency[0][0][closestPInd]/* + 0.02*/));
		//luminosityErr.push_back(luminosity[i]/sqrt((double)entriesExp[0][0][i]) + luminosity[i] * (efficiencyErr[0][0][closestPInd]/ efficiency[0][0][closestPInd]/* + 0.02*/));
		for(size_t k = 0; k < 8; k++){
			for(size_t j = 0; j<varsizes[k]; j++){
				lumforsystematic[k][j].push_back((double)entriesExp[k][j][i]/efficiency[k][j][closestPInd]/crossSectNb(energyExp[i]));
				lumforsystematicErr[k][j].push_back(/*lumforsystematic[k][j][i]/sqrt((double)countModEventsEE[i]) + */lumforsystematic[k][j][i] * (efficiencyErr[k][j][closestPInd]/ efficiency[k][j][closestPInd]/* + 0.02*/));
				//lumforsystematicErr[k][j].push_back(lumforsystematic[k][j][i]/sqrt((double)entriesExp[k][j][i]) + lumforsystematic[k][j][i] * (efficiencyErr[k][j][closestPInd]/ efficiency[k][j][closestPInd]/* + 0.02*/));
			}
		}

		popre.push_back(effee[i]/effem[closestPInd]);
		dpopre.push_back(sqrt(pow(deffee[i]/effem[closestPInd],2) + pow(effee[i]/effem[closestPInd] * deffem[closestPInd]/effem[closestPInd],2)));
		//dpopre.push_back(0);

		poprt.push_back(effte[i]/efftm[closestPInd]);
		dpoprt.push_back(sqrt(pow(deffte[i]/efftm[closestPInd],2) + pow(effte[i]/efftm[closestPInd] * defftm[closestPInd]/efftm[closestPInd],2)));
		//dpoprt.push_back(0);

		poprst.push_back(effste[i]/effstm[closestPInd]);
		dpoprst.push_back(sqrt(pow(deffste[i]/effstm[closestPInd],2) + pow(effste[i]/effstm[closestPInd] * deffstm[closestPInd]/effstm[closestPInd],2)));
		//dpoprst.push_back(0);

		poprp.push_back(effpe[i]/effpm[closestPInd]);
		dpoprp.push_back(sqrt(pow(deffpe[i]/effpm[closestPInd],2) + pow(effpe[i]/effpm[closestPInd] * deffpm[closestPInd]/effpm[closestPInd],2)));
		//dpoprp.push_back(0);

		poprd.push_back(effde[i]/effdm[closestPInd]);
		dpoprd.push_back(sqrt(pow(deffde[i]/effdm[closestPInd],2) + pow(effde[i]/effdm[closestPInd] * deffdm[closestPInd]/effdm[closestPInd],2)));
		//dpoprd.push_back(0);

		poprz.push_back(effze[i]/effzm[closestPInd]);
		dpoprz.push_back(sqrt(pow(deffze[i]/effzm[closestPInd],2) + pow(effze[i]/effzm[closestPInd] * deffzm[closestPInd]/effzm[closestPInd],2)));
		//dpoprz.push_back(0);

		poprc.push_back(effce[i]/effcm[closestPInd]);
		dpoprc.push_back(sqrt(pow(deffce[i]/effcm[closestPInd],2) + pow(effce[i]/effcm[closestPInd] * deffcm[closestPInd]/effcm[closestPInd],2)));
		//dpoprc.push_back(0);
	}
	double tempEnergy = 0;
	double tempLuminosity = 0;
	double tempLuminosityErr = 0;
	double tp = 0;
	double dtp = 0;

	//sort by energy
	for (size_t i = 0; i < energyExp.size(); i++) {
		for (size_t j = 0; j < energyExp.size() - i - 1; j++) {
			if (energyExp[j] > energyExp[j + 1]) {
				tempEnergy = energyExp[j];
				energyExp[j] = energyExp[j + 1];
				energyExp[j + 1] = tempEnergy;

				tempLuminosity = luminosity[j];
				luminosity[j] = luminosity[j + 1];
				luminosity[j + 1] = tempLuminosity;

				tempLuminosityErr = luminosityErr[j];
				luminosityErr[j] = luminosityErr[j + 1];
				luminosityErr[j + 1] = tempLuminosityErr;

				for(size_t k = 0; k < 8; k++){
					for(size_t l = 0; l<varsizes[k]; l++){
						tempLuminosity = lumforsystematic[k][l][j];
						lumforsystematic[k][l][j] = lumforsystematic[k][l][j + 1];
						lumforsystematic[k][l][j + 1] = tempLuminosity;

						tempLuminosityErr = lumforsystematicErr[k][l][j];
						lumforsystematicErr[k][l][j] = lumforsystematicErr[k][l][j + 1];
						lumforsystematicErr[k][l][j + 1] = tempLuminosityErr;
					}
				}

				tp = popre[j];
				popre[j] = popre[j + 1];
				popre[j + 1] = tp;
				dtp = dpopre[j];
				dpopre[j] = dpopre[j + 1];
				dpopre[j + 1] = dtp;

				tp = poprt[j];
				poprt[j] = poprt[j + 1];
				poprt[j + 1] = tp;
				dtp = dpoprt[j];
				dpoprt[j] = dpoprt[j + 1];
				dpoprt[j + 1] = dtp;
				
				tp = poprst[j];
				poprst[j] = poprst[j + 1];
				poprst[j + 1] = tp;
				dtp = dpoprst[j];
				dpoprst[j] = dpoprst[j + 1];
				dpoprst[j + 1] = dtp;
				
				tp = poprp[j];
				poprp[j] = poprp[j + 1];
				poprp[j + 1] = tp;
				dtp = dpoprp[j];
				dpoprp[j] = dpoprp[j + 1];
				dpoprp[j + 1] = dtp;
				
				tp = poprd[j];
				poprd[j] = poprd[j + 1];
				poprd[j + 1] = tp;
				dtp = dpoprd[j];
				dpoprd[j] = dpoprd[j + 1];
				dpoprd[j + 1] = dtp;
				
				tp = poprz[j];
				poprz[j] = poprz[j + 1];
				poprz[j + 1] = tp;
				dtp = dpoprz[j];
				dpoprz[j] = dpoprz[j + 1];
				dpoprz[j + 1] = dtp;

				tp = poprc[j];
				poprc[j] = poprc[j + 1];
				poprc[j + 1] = tp;
				dtp = dpoprc[j];
				dpoprc[j] = dpoprc[j + 1];
				dpoprc[j + 1] = dtp;


				//dExn
				tp = dExnpopr[j];
				dExnpopr[j] = dExnpopr[j + 1];
				dExnpopr[j + 1] = tp;
				dtp = dExnpoprErr[j];
				dExnpoprErr[j] = dExnpoprErr[j + 1];
				dExnpoprErr[j + 1] = dtp;

			}
	
		}
	}
	ofile << endl;
	vector<double> enerrexp;
	for (size_t i = 0; i < energyExp.size(); i++) {
		enerrexp.push_back(0.);
		cout << energyExp[i] << "	" << luminosity[i] << "	" << luminosityErr[i] << endl;
		ofile << energyExp[i] << "	" << luminosity[i] << "	" << luminosityErr[i] << endl;
	}
	ofile.close();

	TGraphErrors* gre = new TGraphErrors(energyExp.size()-2, &energyExp[0], &popre[0], &enerrexp[0], &dpopre[0]);
	TGraphErrors* grt = new TGraphErrors(energyExp.size()-2, &energyExp[0], &poprt[0], &enerrexp[0], &dpoprt[0]);
	TGraphErrors* grst = new TGraphErrors(energyExp.size()-2, &energyExp[0], &poprst[0], &enerrexp[0], &dpoprst[0]);
	TGraphErrors* grp = new TGraphErrors(energyExp.size()-2, &energyExp[0], &poprp[0], &enerrexp[0], &dpoprp[0]);
	TGraphErrors* grd = new TGraphErrors(energyExp.size()-2, &energyExp[0], &poprd[0], &enerrexp[0], &dpoprd[0]);
	TGraphErrors* grz = new TGraphErrors(energyExp.size()-2, &energyExp[0], &poprz[0], &enerrexp[0], &dpoprz[0]);
	TGraphErrors* grc = new TGraphErrors(energyExp.size()-2, &energyExp[0], &poprc[0], &enerrexp[0], &dpoprc[0]);
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R007-001/2011/lumsgraphs.root", "RECREATE");

	gre->SetMarkerColor(4);
	gre->SetTitle("efficiency exp / efficiency mod;beam energy");
	gre->SetMarkerStyle(21);
	gre->Draw("AP");
	TF1* f1 = new TF1("f1", "[0]", 600, 1100);
	f1->SetParameter(0,1.0);
	gre->Fit("f1", "", "", 600, 1100);
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	c->SetBorderMode(0);
	c->SetFrameBorderMode(0);
	c->SetGrid();
	c->Update();
	cin.get();
	gre->Write("e");

	grt->SetMarkerColor(4);
	grt->SetTitle("efficiency exp / efficiency mod;beam energy");
	grt->SetMarkerStyle(21);
	grt->Draw("AP");
	f1->SetParameter(0,1.0);
	grt->Fit("f1", "", "", 600, 1100);
	c->Update();
	cin.get();
	grt->Write("t");
	
	grst->SetMarkerColor(4);
	grst->SetTitle("efficiency exp / efficiency mod;beam energy");
	grst->SetMarkerStyle(21);
	grst->Draw("AP");
	f1->SetParameter(0,1.0);
	grst->Fit("f1", "", "", 600, 1100);
	c->Update();
	cin.get();
	grst->Write("st");

	grp->SetMarkerColor(4);
	grp->SetTitle("efficiency exp / efficiency mod;beam energy");
	grp->SetMarkerStyle(21);
	grp->Draw("AP");
	f1->SetParameter(0,1.0);
	grp->Fit("f1", "", "", 600, 1100);
	c->Update();
	cin.get();
	grp->Write("p");
	
	grd->SetMarkerColor(4);
	grd->SetTitle("efficiency exp / efficiency mod;beam energy");
	grd->SetMarkerStyle(21);
	grd->Draw("AP");
	f1->SetParameter(0,1.0);
	grd->Fit("f1", "", "", 600, 1100);
	c->Update();
	cin.get();
	grd->Write("d");
	
	grz->SetMarkerColor(4);
	grz->SetTitle("efficiency exp / efficiency mod;beam energy");
	grz->SetMarkerStyle(21);
	grz->Draw("AP");
	f1->SetParameter(0,1.0);
	grz->Fit("f1", "", "", 600, 1100);
	c->Update();
	cin.get();
	grz->Write("z");

	grc->SetMarkerColor(4);
	grc->SetTitle("efficiency exp / efficiency mod;beam energy");
	grc->SetMarkerStyle(21);
	grc->Draw("AP");
	f1->SetParameter(0,1.0);
	grc->Fit("f1", "", "", 600, 1100);
	c->Update();
	//cin.get();
	grc->Write("c");

	MyFile->Close();
}

//all process summed up
void luminCalculation::luminosityCalc() {
	//effEeEnergy("h1", "/work/users/kladov/snd2k/R007-001/2019/*0.root", "/work/users/kladov/snd2k/R007-001/2019/lumdistrMod.root");
	
	//effEeEnergy("t1", "/work/users/kladov/snd2k/R007-002/output/ntuples/eemod/2017n/*odch*x.root", "/work/users/kladov/snd2k/R007-001/2017/lumdistrMod.root");
	effEeEnergy("t1", "/work/users/kladov/snd2k/R007-002/output/ntuples/eemod/2012n/dEdx/*x.root", "/work/users/kladov/snd2k/R007-001/2012/lumdistrMod.root");
	//effEeEnergy("t1", "/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2017/*col.root", "/work/users/kladov/snd2k/R007-001/2017/lumdistrMod.root");
	
	//effEeEnergy("t1", "/work/users/kladov/snd2k/R007-002/output/ntuples/eemod/2017/*x.root", "/work/users/kladov/snd2k/R007-001/2017/lumdistrMod.root");
	//effEeEnergy("t1", "/work/users/kladov/snd2k/R007-001/2011/eeCorr/*x.root", "/work/users/kladov/snd2k/R007-001/2011/lumdistr1Mod.root");
	for (size_t i = 0; i < energyPoints.size(); i++) {
		energyMod.push_back(energyPoints[i]);
		energyErr.push_back(0.);
		for(size_t k = 0; k < 8; k++){
			for(size_t l = 0; l<varsizes[k]; l++){
				efficiency[k][l].push_back((double)countSelectedEventsEE[k][l][i] / (double)countModEventsEE[i]);
				efficiencyErr[k][l].push_back(sqrt(efficiency[k][l].back() * (1.-efficiency[k][l].back()) / (double)countModEventsEE[i]));// +efficiency.back() * 0.02);
			}
		}
	}

	TGraphErrors* gr = new TGraphErrors(energyMod.size(),&energyMod[0], &efficiency[3][0][0], &energyErr[0], &efficiencyErr[3][0][0]);
	gr->SetMarkerColor(4);
	gr->SetTitle("ee scatterig registration efficiency vs energy from modeling files");
	gr->SetMarkerStyle(21);
	gr->Draw("AP");
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R007-001/2012/efficiency_en.root", "RECREATE");
	gr->Write("effen");
	MyFile->Close();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	//cin.get();

	energyPoints.clear();
	countModEventsEE.clear();
	for(size_t k = 0; k < 8; k++){
		for(size_t l = 0; l<varsizes[k]; l++){
			countSelectedEventsEE[k][l].clear();
		}
	}
	//effEeEnergy("t1", "/work/users/konctbel/calibs/R007-001/output/ntuples/MHAD2019/col/*col.root","/work/users/kladov/snd2k/R007-001/2019/lumdistrExp.root");
	effEeEnergy("t1", "/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2012/*col.root","/work/users/kladov/snd2k/R007-001/2012/lumdistrExp.root");
	//effEeEnergy("t1", "/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2017/*col.root","/work/users/kladov/snd2k/R007-001/2017/lumdistrExp.root");
	//effEeEnergy("t1", "/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2011/*col.root","/work/users/kladov/snd2k/R007-001/2011/lumdistr1Exp.root");
	for (size_t i = 0; i < energyPoints.size(); i++) {
		energyExp.push_back(energyPoints[i]);
		for(size_t k = 0; k < 8; k++){
			for(size_t l = 0; l<varsizes[k]; l++){
				entriesExp[k][l].push_back(countSelectedEventsEE[k][l][i]);
			}
		}
	}
	cout << "expProcessed" << endl;

	vector<double> enparticle;
	vector<double> enparticleErr;
	vector<double> coeff;
	vector<double> coeffErr;
	dedxm->Draw();
	c->Update();
	//cin.get();
	dedxe->Draw();
	c->Update();

	for (size_t i = 0; i < energyExp.size(); i++) {
		cout << "a" << endl;
		size_t closestPInd;
		double mindist = 1000;
		for (size_t j = 0; j < energyMod.size(); j++) {
			if (fabs(energyMod[j] - energyExp[i]) < mindist) {
				mindist = fabs(energyMod[j] - energyExp[i]);
				closestPInd = j;
			}
		}
		cout << "b" << endl;
		cout << dEdxMod1.size() << endl;
		cout << closestPInd << endl;
		dEdxMod1[closestPInd]->SetLineColor(2);
		dEdxMod1[closestPInd]->Draw();
		cout << "b0" << endl;
		double maximumCenter = dEdxMod1[closestPInd]->GetBinCenter(dEdxMod1[closestPInd]->GetMaximumBin());
		double rMS = dEdxMod1[closestPInd]->GetRMS();
		cout << "ba" << endl;
		TF1* f1 = new TF1("f1", "[0]*Landau(-x,[1],[2])", maximumCenter-rMS, maximumCenter + rMS);
		f1->SetParameters(dEdxMod1[closestPInd]->GetMaximum(), maximumCenter, rMS);
		dEdxMod1[closestPInd]->Fit("f1", "", "", maximumCenter - rMS, maximumCenter + rMS);
		
		//TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->Update();
		cin.get();
		cout << "c" << endl;
		dEdxExp1[i]->SetLineColor(2);
		dEdxExp1[i]->Draw();
		maximumCenter = dEdxExp1[i]->GetBinCenter(dEdxExp1[i]->GetMaximumBin());
		rMS = dEdxExp1[i]->GetRMS();
		TF1* f2 = new TF1("f1", "[0]*Landau(-x,[1],[2])", maximumCenter - rMS, maximumCenter + rMS);
		f2->SetParameters(dEdxExp1[i]->GetMaximum(), maximumCenter, rMS);
		dEdxExp1[i]->Fit("f2", "", "", maximumCenter - rMS, maximumCenter + rMS);
		c->Update();
		cin.get();
	}

	//dEx for beam dependence 
	vector<double> nulvector;
	for (size_t i = 0; i < energyExp.size(); i++) {
		size_t closestPInd;
		double mindist = 1000;
		for (size_t j = 0; j < energyMod.size(); j++) {
			if (fabs(energyMod[j] - energyExp[i]) < mindist) {
				mindist = fabs(energyMod[j] - energyExp[i]);
				closestPInd = j;
			}
		}
		dExnpopr.push_back(meandExn[closestPInd]/meandExnC[i]);
		dExnpoprErr.push_back(meandExnErr[closestPInd]/meandExnC[i] + dExnpopr.back()*meandExnCErr[i]/meandExnC[i]);
		nulvector.push_back(0.);
	}
	//pointComparisonAndSort();
	//dEx for energy dependence
	for(size_t i = 2; i < dedxm->GetNbinsX()-1; i++){
		if(dedxe->GetBinEntries(i)>0){
			enparticle.push_back(dedxm->GetBinCenter(i));
			enparticleErr.push_back(1000./100.);
			coeff.push_back(dedxm->GetBinContent(i)/dedxe->GetBinContent(i));
			coeffErr.push_back(sqrt(pow(dedxm->GetBinError(i)/dedxe->GetBinContent(i),2) + pow(dedxm->GetBinContent(i)/dedxe->GetBinContent(i)*dedxe->GetBinError(i)/dedxe->GetBinContent(i),2)));
			cout << enparticle.back() << "	" << coeff.back() << "	+-	" << coeffErr.back() << endl;
		}
	}


	//resulting correction output
	vector<double> tempEnSorted;
	for(size_t i = 0; i < energyExp.size(); i++)
		tempEnSorted.push_back(energyExp[i]);
	//sort by energy
	double tp = 0;
	for (size_t i = 0; i < tempEnSorted.size(); i++) {
		for (size_t j = 0; j < tempEnSorted.size() - i - 1; j++) {
			if (tempEnSorted[j] > tempEnSorted[j + 1]) {
				tp = tempEnSorted[j];
				tempEnSorted[j] = tempEnSorted[j + 1];
				tempEnSorted[j + 1] = tp;

				tp = dExnpopr[j];
				dExnpopr[j] = dExnpopr[j + 1];
				dExnpopr[j + 1] = tp;

				tp = dExnpoprErr[j];
				dExnpoprErr[j] = dExnpoprErr[j + 1];
				dExnpoprErr[j + 1] = tp;
			}
		}
	}
	//TGraphErrors* gr1 = new TGraphErrors(enparticle.size(), &enparticle[0], &coeff[0], &enparticleErr[0], &coeffErr[0]);
	TGraphErrors* gr1 = new TGraphErrors(tempEnSorted.size()-1, &tempEnSorted[1], &dExnpopr[1], &nulvector[1], &dExnpoprErr[1]);
	gr1->SetMarkerColor(4);
	gr1->SetTitle("dedx correction");
	gr1->SetMarkerStyle(21);
	gr1->Draw("AP");
	TF1* f1 = new TF1("f1", "[0]", 500, 1100);
	f1->SetParameter(0,950);
	gr1->Fit("f1", "", "", 500, 1100);
	c->Update();
	cout << f1->GetParameter(0) << "	" << f1->GetParError(0) << endl;
	//cout
	for(size_t i = 0; i < tempEnSorted.size(); i++){
		cout << tempEnSorted[i] << ", ";
	}
	cout << endl;
	for(size_t i = 0; i < tempEnSorted.size(); i++){
		cout << dExnpopr[i] << ", ";
	}
	cout << endl;

	//cin.get();
	//pointComparisonAndSort();
}

double pep2lumt(double m, double L, double t, double dm){
    double x, s, em, ep, k;
    double c, pi, scal, alpha;
    alpha=1.0/137;
    pi=3.1415927;
    scal=pi/180;
    c=cos(t);
    em=9;
    ep=3.1;
    s=4*em*ep;
    x=1.-m*m/s;
    k=1e+06*dm;
    double pep2lumt=k*alpha/pi/x*m/s*L*((2.-2.*x+x*x)*log((1+c)/(1-c))-x*x*c);
    return pep2lumt;
}

double babarLum(double m, double L, double dm){
    double pep2lum=pep2lumt(m,L,0.35,dm)+pep2lumt(m,L,3.1415927-2.4,dm);
    return pep2lum;
}

void drawLumBabar(){
	TF1* f1 = new TF1("f1","babarLum(x,535.,0.025)",0,2);
	f1->Draw();
}

void drawLum() {
	vector<double> luminosityEn;
	vector<double> luminosity;
	vector<double> luminosityt;
	vector<double> luminosityErr;
	vector<double> enErr;
	ifstream ifile("/work/users/kladov/snd2k/R007-001/2011/luminosity.dat");
	//ifstream ifile("/work/users/kladov/snd2k/R007-001/2022/luminosity1.dat");
	vector<double> luminosityEn1;
	vector<double> luminosity1;
	vector<double> luminosityt1;
	vector<double> luminosityErr1;
	vector<double> enErr1;
	ifstream ifile1("/work/users/kladov/snd2k/R007-001/2012/luminosity.dat");
	vector<double> luminosityEn2;
	vector<double> luminosity2;
	vector<double> luminosityt2;
	vector<double> luminosityErr2;
	vector<double> enErr2;
	ifstream ifile2("/work/users/kladov/snd2k/R007-001/2017/luminosity.dat");

	TH1* lumh = new TH1F("lumh", "luminosity density", 80, 0+1./80, 2+1./80);

	float sumlum = 0;
	while (ifile1.get() != EOF) {
		double a, b, c;
		ifile1 >> a >> b >> c;
		//if (luminosityEn1.size() == 0 || a != luminosityEn1.back()) {
			luminosityEn1.push_back(2*a);
			luminosityt1.push_back(b/1000.);
			double suml =0;
			for (size_t j = 0; j < luminosityEn1.size(); j++)
				suml+=luminosityt1[j];
			//luminosity1.push_back(suml);
			luminosity1.push_back(b);
			luminosityErr1.push_back(c);
			enErr1.push_back(0.);
			cout << a * 2 << "	" << b << "	" << c << endl;
			lumh->Fill(2./1000.*a,b);
			sumlum+=b/1000.;
		//}
	}
	ifile1.close();
	while (ifile2.get() != EOF) {
		double a, b, c;
		ifile2 >> a >> b >> c;
		//if (luminosityEn2.size() == 0 || a != luminosityEn2.back()) {
			luminosityEn2.push_back(2*a);
			luminosityt2.push_back(b/1000.);
			double suml =0;
			for (size_t j = 0; j < luminosityEn2.size(); j++)
				suml+=luminosityt2[j];
			//luminosity2.push_back(suml);
			luminosity2.push_back(b);
			luminosityErr2.push_back(c);
			enErr2.push_back(0.);
			cout << a * 2 << "	" << b << "	" << c << endl;
			lumh->Fill(2./1000.*a,b);
			sumlum+=b/1000.;
		//}

	}
	ifile2.close();
	while (ifile.get() != EOF) {
		double a, b, c;	
		ifile >> a >> b >> c;
		//for (size_t j = 0; j < luminosityEn.size(); j++) {
			//if (luminosityEn1[j] == a && (luminosityEn.size() == 0 || a != luminosityEn.back())) {
				luminosityEn.push_back(2*a);
				luminosityt.push_back(b/1000.);
				double suml =0;
				for (size_t j = 0; j < luminosityEn.size(); j++)
					suml+=luminosityt[j];
				//luminosity.push_back(suml);
				luminosity.push_back(b);
				luminosityErr.push_back(c);
				//luminosityErr.push_back(0.);
				enErr.push_back(0.);
				cout << a * 2 << "	" << b << "	" << "0" << endl;
				lumh->Fill(2./1000.*a,b);
				sumlum+=b/1000.;
			//}
		//}
	}
	cout << "sumlum	=	" << sumlum << endl; 
	ifile.close();
	TGraphErrors* gr1 = new TGraphErrors(luminosityEn1.size(), &luminosityEn1[0], &luminosity1[0], &enErr1[0], &luminosityErr1[0]);
	TGraphErrors* gr2 = new TGraphErrors(luminosityEn2.size(), &luminosityEn2[0], &luminosity2[0], &enErr2[0], &luminosityErr2[0]);
	TGraphErrors* gr = new TGraphErrors(luminosityEn.size(), &luminosityEn[0], &luminosity[0], &enErr[0], &luminosityErr[0]);
	gr1->SetLineColor(1);
	gr1->SetMarkerColor(2);
	gr1->SetMarkerStyle(21);
	gr1->SetTitle("Integrated luminosity; beam energy, MeV; IL, pb^{-1}");
	
	gr->SetLineColor(1);
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->SetTitle("Integrated luminosity; beam energy, MeV; IL, pb^{-1}");
	
	gr2->SetLineColor(1);
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(21);
	gr2->SetTitle("Integrated luminosity; beam energy, MeV; IL, pb^{-1}");

	/*gr2->Draw("AP");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	c->SetBorderMode(0);
	c->SetFrameBorderMode(0);
	c->SetGrid();
	c->Update();
	gr->Draw("Psame");
	gr1->Draw("Psame");
	c->Update();*/
	lumh->Draw();
	TF1* f1 = new TF1("f1","babarLum(x,535.,0.025)",0,2);
	f1->Draw("same");
	TF1* f2 = new TF1("f2","babarLum(x,232.,0.025)",0,2);
	f2->Draw("same");
}



void sumpoints(){
	ifstream ifile("/work/users/kladov/snd2k/R007-001/2019/luminosity.dat");
	ifstream ifile1("/work/users/kladov/snd2k/R007-001/2020/luminosity.dat");
	ifstream ifile2("/work/users/kladov/snd2k/R007-001/2021/luminosity.dat");
	ifstream ifile3("/work/users/kladov/snd2k/R007-001/2022/luminosity.dat");
	ofstream ofile("/work/users/kladov/snd2k/R007-001/2022/luminosity1.dat");
	vector<double> luminosityEn;
	vector<double> luminosity;
	vector<double> enErr;
	while (ifile1.get() != EOF) {
		double a, b;
		ifile1 >> a >> b;
		bool nP = true;
		for(size_t i = 0; i < luminosityEn.size();i++){
			if (luminosityEn[i] == a){
				nP = false;
				luminosity[i] += (b);
			}
		}
		if (nP){
			luminosityEn.push_back(a);
			luminosity.push_back(b);
		}
	}
	while (ifile.get() != EOF) {
		double a, b;
		ifile >> a >> b;
		bool nP = true;
		for(size_t i = 0; i < luminosityEn.size();i++){
			if (luminosityEn[i] == a){
				nP = false;
				luminosity[i] += (b);
			}
		}
		if (nP){
			luminosityEn.push_back(a);
			luminosity.push_back(b);
		}
	}
	while (ifile2.get() != EOF) {
		double a, b;
		ifile2 >> a >> b;
		bool nP = true;
		for(size_t i = 0; i < luminosityEn.size();i++){
			if (luminosityEn[i] == a){
				nP = false;
				luminosity[i] += (b);
			}
		}
		if (nP){
			luminosityEn.push_back(a);
			luminosity.push_back(b);
		}
	}
	while (ifile3.get() != EOF) {
		double a, b;
		ifile3 >> a >> b;
		bool nP = true;
		for(size_t i = 0; i < luminosityEn.size();i++){
			if (luminosityEn[i] == a){
				nP = false;
				luminosity[i] += (b);
			}
		}
		if (nP){
			luminosityEn.push_back(a);
			luminosity.push_back(b);
		}
	}
	ifile.close();
	ifile1.close();
	ifile2.close();
	ifile3.close();

	double tempEnergy;
	for (size_t i = 0; i < luminosityEn.size(); i++) {
		for (size_t j = 0; j < luminosityEn.size() - i - 1; j++) {
			if (luminosityEn[j] > luminosityEn[j + 1]) {
				tempEnergy = luminosityEn[j];
				luminosityEn[j] = luminosityEn[j + 1];
				luminosityEn[j + 1] = tempEnergy;

				tempEnergy = luminosity[j];
				luminosity[j] = luminosity[j + 1];
				luminosity[j + 1] = tempEnergy;
			}
		}
	}
	for(size_t i = 0; i < luminosityEn.size();i++){
		ofile << luminosityEn[i] << "	" << luminosity[i] << endl;
	}
	ofile.close();
}

void luminCalculation::differenceForLum() {
	/*vector<double> luminosityEn;
	vector<double> luminosity;
	vector<double> luminosityErr;
	vector<double> enErr;
	vector<double> luminosityEn1;
	vector<double> luminosity1;
	vector<double> luminosityErr1;
	vector<double> enErr1;*/
	//ifstream ifile("/work/users/kladov/snd2k/R007-001/2011/luminosity1.dat");
	//ifstream ifile1("/work/users/kladov/snd2k/R007-001/2011/luminosity2.dat");

	/*
	while (ifile1.get() != EOF) {
		double a, b, c;
		ifile1 >> a >> b >> c;
		if (luminosityEn1.size() == 0 || a != luminosityEn1.back()) {
			luminosityEn1.push_back(a);
			luminosity1.push_back(b);
			luminosityErr1.push_back(c);
			enErr1.push_back(0.);
			cout << a * 2 << "	" << b << "	" << c << endl;
		}
	}
	ifile1.close();
	while (ifile.get() != EOF) {
		double a, b, c;
		ifile >> a >> b >> c;
		for (size_t j = 0; j < luminosityEn1.size(); j++) {
			if (luminosityEn1[j] == a && (luminosityEn.size() == 0 || a != luminosityEn.back())) {
				luminosityEn.push_back(a);
				luminosity.push_back(b);
				luminosityErr.push_back(c);
				cout << a * 2 << "	" << b << "	" << c << endl;
			}
		}
	}
	*/
	vector<double> meandiff[8];
	vector<double> maxdiff[8];
	cout << "energyExp" << "	" << "meandiff" << "	" << "maxdiff" << endl;
	
	cout << std::fixed;
	cout << setprecision(3);
	TFile* MyFile = new TFile("/work/users/kladov/snd2k/R007-001/2011/lumsgraphsTheta.root", "RECREATE");
	for (size_t i = 0; i < energyExp.size(); i++) {
		for(size_t a = 0; a < 8; a++){
			meandiff[a].push_back(0);
			maxdiff[a].push_back(0);
		}
		double diff[21];
		double differr[21];
		double tcut[21];
		double tcuterr[21];
		for(size_t tt = 0; tt<21;tt++){
			//tcut[tt]=0.75-0.075*tt;
			if (tt == 0)
				tcut[tt] = 45;
			else if (tt>0 && tt<5)
				tcut[tt] = 45 + 2.5*((int)tt-5);
			else
				tcut[tt] = 45 + 2.5*((int)tt-4);
			tcuterr[tt]=0.0;
		}
		//cout << fabs(luminosity[i] - luminosity1[i])/ luminosity[i] << "	" << luminosityErr[i]/ luminosity[i] << endl;
		for(size_t k = 0; k < 8; k++){
			for(size_t l = 0; l<varsizes[k]; l++){
				meandiff[k][i] += pow(lumforsystematic[k][l][i] - luminosity[i],2)/ pow(luminosity[i],2) / (double)varsizes[k];
				if (fabs(lumforsystematic[k][l][i] - luminosity[i])/ luminosity[i] > maxdiff[k][i])
					maxdiff[k][i] = fabs(lumforsystematic[k][l][i] - luminosity[i])/ luminosity[i];
				if (k==3){
					diff[l]=(lumforsystematic[k][l][i] - luminosity[i])/luminosity[i]*100;
					differr[l]=sqrt((pow(lumforsystematicErr[k][l][i],2) + pow(luminosityErr[i],2))/pow(luminosity[i],2))*100;
				}
			}
		}
		
		TGraphErrors* gr = new TGraphErrors(21, tcut, diff, tcuterr, differr);
		gr->SetTitle("TGraphErrors Example");
		gr->SetMarkerColor(4);
		gr->SetMarkerStyle(21);
		gr->Draw("AP");
		gr->Write(Form("t%d",(int)i+1));
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->SetFillColor(0);
		c->SetFrameFillColor(0);
		c->SetBorderMode(0);
		c->SetFrameBorderMode(0);
		c->SetGrid();
		c->Update();
		cin.get();

		cout << energyExp[i] << "	";
		for(size_t a = 0; a < 8; a++){
			cout << maxdiff[a][i]*100. << "	";
		}
		cout << endl;
	}
	
	MyFile->Close();

	ofstream ofile("/work/users/kladov/snd2k/R007-001/2011/luminositySyst.dat");
	for (size_t i = 0; i < energyExp.size(); i++) {
		for(size_t k = 0; k < 8; k++){
			for(size_t l = 0; l<varsizes[k]; l++){
				ofile << energyExp[i] << "	" << lumforsystematic[k][l][i] << "	";
			}
			ofile << endl;
		}
		ofile << endl;
	}
	ofile.close();
}

void luminCalculation::drawdEdxspectrum(){
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2017/*col.root");
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	{
		chain.SetBranchAddress("beam", &beam);
		chain.SetBranchAddress("energy", &energy);
		chain.SetBranchAddress("nc", &nc);
		chain.SetBranchAddress("nn", &nn);
		chain.SetBranchAddress("z0", &z0);
		chain.SetBranchAddress("d0", &d0);
		chain.SetBranchAddress("phi", &phi);
		chain.SetBranchAddress("theta", &theta);
		chain.SetBranchAddress("col", &col);
		chain.SetBranchAddress("cosm", &cosm);
		chain.SetBranchAddress("act", &act);
		chain.SetBranchAddress("cact", &cact);
		chain.SetBranchAddress("dExn", &dExn);
		chain.SetBranchAddress("dExnC", &dExnC);
		

		chain.SetBranchAddress("dEx1C", &dEx1C);
		chain.SetBranchAddress("dEx2C", &dEx2C);
		chain.SetBranchAddress("dEx3C", &dEx3C);
		chain.SetBranchAddress("dEx4C", &dEx4C);
		chain.SetBranchAddress("dEx5C", &dEx5C);
		chain.SetBranchAddress("dEx6C", &dEx6C);
		chain.SetBranchAddress("dEx7C", &dEx7C);
		chain.SetBranchAddress("dEx8C", &dEx8C);
		chain.SetBranchAddress("dEx9C", &dEx9C);
	}
	TH1* h = new TH1F("h", "dEdxC distribution", 100, 0, 15);
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		if (isItEE(&parameters[0][0][0])) {
			if (1){
				double partialded0[9] = {dEx1C[0], dEx2C[0], dEx3C[0], dEx4C[0], dEx5C[0], dEx6C[0], dEx7C[0], dEx8C[0], dEx9C[0]};
				double partialded1[9] = {dEx1C[1], dEx2C[1], dEx3C[1], dEx4C[1], dEx5C[1], dEx6C[1], dEx7C[1], dEx8C[1], dEx9C[1]};
				double dExnTrue0 = 0.;
				double dExnTrue1 = 0.;
				int normcount0 = 0;
				int normcount1 = 0;
				for(size_t i = 0; i<9; i++){
					if(partialded0[i] < 3.5 && partialded0[i] > 0.01){
						normcount0+=1;
						dExnTrue0+=partialded0[i];
					}
					if(partialded1[i] < 3.5 && partialded1[i] > 0.01){
						normcount1+=1;
						dExnTrue1+=partialded1[i];
					}
				}
				if (normcount0 > 0)
					dExnTrue0 = dExnTrue0/normcount0;
				else
					dExnTrue0 = 0;
				if (normcount1 > 0)
					dExnTrue1 = dExnTrue1/normcount1;
				else
					dExnTrue1 = 0;

				//cout << energy[0] << "	" << dExnC[0] << endl;
				//dedxT->Fill(beam, dExnC[0]);
				//dedxT->Fill(beam, dExnC[1]);
				if(dExnTrue0>=0.01 && dExnTrue1>=0.01){
					h->Fill(dExnTrue0);
					h->Fill(dExnTrue1);
				}
				//h->Fill(dExnC[0]);
				//h->Fill(dExnC[1]);
				//}
			}
		}
	}
	h->Draw();
}



class kkpDistr {
public:
	Entryvar sEntr;
	float phi[40], phis[40], ach[40], achr[40], theta[40], z0[40], schr[40], tch[40], tchr[40], energy[40], beam, eton, x0[40], y0[40], d0[40], dExs[40], dExnC[40], dExn[40];  //kinematic, counters parameters float
	float d2phi[40], dphirho[40], d2rho[40], d2z0[40], d2cosTh[40], dz0cosTh[40], Dtheta[40], Dphi[40], energyerr[40];  //errors afrer reconstruction
	int nch, run, nc, nn, cosm, eventtime, act, region[40], charge[40];     //parameters int
	float amplitude[40], amplitudeE[40]; //cherenkov general amplitude
	float x2ikf1, x2kf1, ppkf1[11], thetakf1[11], phikf1[11];  //kinfit
	int npkf1[11], ipkf1[11];   //kinfit 
	float x2ikf2, x2kf2, ppkf2[11], thetakf2[11], phikf2[11];  //kinfit
	int npkf2[11], ipkf2[11];   //kinfit 
	float x2ikf3, x2kf3, ppkf3[11], thetakf3[11], phikf3[11];  //kinfit
	int npkf3[11], ipkf3[11];   //kinfit 
	float mctheta[40], mcphi[40], mcpx[40], mcpy[40], mcpz[40], mce[40];	//montecarlo modeling
	int mcpdg[40], nmc;	//montecarlo modeling
	float dEx1C[10], dEx2C[10], dEx3C[10], dEx4C[10], dEx5C[10], dEx6C[10], dEx7C[10], dEx8C[10], dEx9C[10];
	vector<int> countModEvents;
	vector<int> countExpEvents;
	bool kf10, kf11, kf20, kf21;

	int hetonInd;
	int beamInd;
	int closestPind;

	double tempdExn[4];
	double tempamp[4];
	double tempampE[4];

	vector<double> energiesEXP;
	vector<double> dEdxpopr;
	int year;
	vector<double> energiesMod;
	vector<double> efficienciesModBeam;
	vector<double> countModEventsBeam;


	vector<double> pareffmce[5];
	vector<double> pareffmceErr[5];

	vector<bool> conditions;

	kkpDistr() {
		year = 2011;
	}
	void definededpopr(){
		switch (year){
			case 2011:
				for(size_t i =0; i<41; i++){
					energiesEXP.push_back(temparrayE2011[i]);
					dEdxpopr.push_back(temparraydE2011[i]);
				}
				break;
			case 2012:
				for(size_t i =0; i<25; i++){
					energiesEXP.push_back(temparrayE2012[i]);
					dEdxpopr.push_back(temparraydE2012[i]);
				}
				break;
			case 2017:
				for(size_t i =0; i<25; i++){
					energiesEXP.push_back(temparrayE2017[i]);
					dEdxpopr.push_back(temparraydE2017[i]);
				}
				break;
		}
	}

	void findhetonInd(double entemp, vector<pair<double, double> > eranges){
		hetonInd = -1;
		for (size_t i = 0; i < eranges.size(); i++)
			if (2. * entemp >= eranges[i].first && 2. * entemp < eranges[i].second)
				hetonInd = i;
	}
	void findbeamInd(double entemp){
		beamInd = -1;
		for (size_t i = 0; i < energiesEXP.size(); i++)
			if(fabs(entemp - energiesEXP[i])<0.4)
				beamInd = i;
	}
	void findclosestPind(double entemp){
		
		closestPind = -1;
		double mindist = 1000;
		for (size_t i = 0; i < energiesEXP.size(); i++) {
			
			if (fabs(entemp - energiesEXP[i]) < mindist) {
				mindist = fabs(entemp - energiesEXP[i]);
				closestPind = i;
			}
		}
	}

	bool indexCheck(string regime){
		return ((mce[1] < 1000 && regime == "mod" && closestPind != -1) || (regime == "exp" && beamInd !=-1)) && hetonInd!=-1;
	}
	
	void normalizedExn(string regime){
		normcoeff = dEdxpopr[beamInd];
		if(regime == "mod")
			normcoeff = dEdxpopr[closestPind];

		tempdExn[0] = dExn[ipkf1[0] - 1] / normcoeff;
		tempdExn[1] = dExn[ipkf1[1] - 1] / normcoeff;
		tempdExn[2] = dExn[ipkf2[0] - 1] / normcoeff;
		tempdExn[3] = dExn[ipkf2[1] - 1] / normcoeff;

		if(regime=="exp"){
			double partialded[4][9] = {dEx1C[ipkf1[0]-1], dEx2C[ipkf1[0]-1], dEx3C[ipkf1[0]-1], dEx4C[ipkf1[0]-1], dEx5C[ipkf1[0]-1], dEx6C[ipkf1[0]-1], dEx7C[ipkf1[0]-1], dEx8C[ipkf1[0]-1], dEx9C[ipkf1[0]-1] , \
									   dEx1C[ipkf1[1]-1], dEx2C[ipkf1[1]-1], dEx3C[ipkf1[1]-1], dEx4C[ipkf1[1]-1], dEx5C[ipkf1[1]-1], dEx6C[ipkf1[1]-1], dEx7C[ipkf1[1]-1], dEx8C[ipkf1[1]-1], dEx9C[ipkf1[1]-1] , \
									   dEx1C[ipkf2[0]-1], dEx2C[ipkf2[0]-1], dEx3C[ipkf2[0]-1], dEx4C[ipkf2[0]-1], dEx5C[ipkf2[0]-1], dEx6C[ipkf2[0]-1], dEx7C[ipkf2[0]-1], dEx8C[ipkf2[0]-1], dEx9C[ipkf2[0]-1] , \
									   dEx1C[ipkf2[1]-1], dEx2C[ipkf2[1]-1], dEx3C[ipkf2[1]-1], dEx4C[ipkf2[1]-1], dEx5C[ipkf2[1]-1], dEx6C[ipkf2[1]-1], dEx7C[ipkf2[1]-1], dEx8C[ipkf2[1]-1], dEx9C[ipkf2[1]-1]};
			double dExnTrue[4] = { 0, 0, 0, 0 };
			int normcount10 = 0;
			int normcount11 = 0;
			int normcount20 = 0;
			int normcount21 = 0;
			/*for (size_t i = 0; i<9; i++) {
				if(partialded10[i] < 15. && partialded10[i] > 0.001){
					normcount10 += 1;
					dExnTrue10 += partialded10[i];
				}
				if(partialded11[i] < 15. && partialded11[i] > 0.001){
					normcount11+=1;
					dExnTrue11+=partialded11[i];
				}
				if(partialded20[i] < 15. && partialded20[i] > 0.001){
					normcount20+=1;
					dExnTrue20+=partialded20[i];
				}
				if(partialded21[i] < 15. && partialded21[i] > 0.001){
					normcount21+=1;
					dExnTrue21+=partialded21[i];
				}
			}*/
			double maxValue[4] = { 0, 0, 0, 0 };
			size_t maxIndex[4] = { 0, 0, 0, 0 };
			for (size_t j = 0; j < 4; j++) {
				for (size_t i = 0; i < 9; i++) {
					dExnTrue[j] += partialded[j][i];
					if (partialded[j][i] > maxValue[j]) { maxValue[j] = partialded[j][i];  maxIndex[j] = i; } //find maximum
				}
				dExnTrue[j] -= partialded[j][maxIndex[j]]; //subtract maximum value
				tempdExn[j] = dExnTrue[j];
			}

			/*if (normcount10 > 0)
				dExnTrue10 = dExnTrue10/(double)normcount10;
			else
				dExnTrue10 = 0;
			if(normcount11 > 0)
				dExnTrue11 = dExnTrue11/(double)normcount11;
			else
				dExnTrue11 = 0;
			if(normcount20 > 0)
				dExnTrue20 = dExnTrue20/(double)normcount20;
			else
				dExnTrue20 = 0;
			if(normcount21 > 0)
				dExnTrue21 = dExnTrue21/(double)normcount21;
			else
				dExnTrue21 = 0;
			*/

			//tempdExn[ipkf1[0] - 1] = dExnC[ipkf1[0] - 1];
			//tempdExn[ipkf1[1] - 1] = dExnC[ipkf1[1] - 1];
			//tempdExn[ipkf2[0] - 1] = dExnC[ipkf2[0] - 1];
			//tempdExn[ipkf2[1] - 1] = dExnC[ipkf2[1] - 1];
		}
	}
	void normalizeAch(){
		tempamp[0] = amplitude[ipkf1[0] - 1]/amplitudeE[ipkf1[0] - 1];
		tempamp[1] = amplitude[ipkf1[1] - 1]/amplitudeE[ipkf1[1] - 1];
		tempamp[2] = amplitude[ipkf2[0] - 1]/amplitudeE[ipkf2[0] - 1]; 
		tempamp[3] = amplitude[ipkf2[1] - 1]/amplitudeE[ipkf2[1] - 1];

		tempampE[0] = fabs(amplitude[ipkf1[0] - 1] - amplitudeE[ipkf1[0] - 1]) / sqrt(amplitudeE[ipkf1[0] - 1]);
		tempampE[1] = fabs(amplitude[ipkf1[1] - 1] - amplitudeE[ipkf1[1] - 1]) / sqrt(amplitudeE[ipkf1[1] - 1]);
		tempampE[2] = fabs(amplitude[ipkf2[0] - 1] - amplitudeE[ipkf2[0] - 1]) / sqrt(amplitudeE[ipkf2[0] - 1]);
		tempampE[3] = fabs(amplitude[ipkf2[1] - 1] - amplitudeE[ipkf2[1] - 1]) / sqrt(amplitudeE[ipkf2[1] - 1]);
		/*tempampE[0] = amplitudeE[ipkf1[0] - 1];
		tempampE[1] = amplitudeE[ipkf1[1] - 1];
		tempampE[2] = amplitudeE[ipkf2[0] - 1]; 
		tempampE[3] = amplitudeE[ipkf2[1] - 1];*/
		if (ppkf1[0] < 300){
			tempamp[0] = 0.0;
			tempampE[0] = 1.0;
		}
		if (ppkf1[1] < 300){
			tempamp[1] = 0.0;
			tempampE[1] = 1.0;
		}
		if (ppkf2[0] < 300){
			tempamp[2] = 0.0;
			tempampE[2] = 1.0;
		}
		if (ppkf2[1] < 300){
			tempamp[3] = 0.0;
			tempampE[3] = 1.0;
		}
		if(x2kf3>99)
			x2kf3=100;
	}


	//
	double kpiInvMassMod(){
		double pxp = mcpx[0];
		double pyp = mcpy[0];
		double pzp = mcpz[0];
		double pp = sqrt(pow(pxp,2)+pow(pyp,2)+pow(pzp,2));
		double ep = mce[0];
		
		double pxk = mcpx[3];
		double pyk = mcpy[3];
		double pzk = mcpz[3];
		double kp = sqrt(pow(pxk,2)+pow(pyk,2)+pow(pzk,2));
		double ek = mce[3];

		double kstpx = pxk + pxp;
		double kstpy = pyk + pyp;
		double kstpz = pzk + pzp;

		double kstp = sqrt(pow(kstpx, 2) + pow(kstpy, 2) + pow(kstpz, 2));
		double kpiInvMass = sqrt(pow(ek + ep, 2) - pow(kstp, 2));

		return kpiInvMass;
	}
	//same
	double kspiInvMassMod(){
		double pxg[4], pyg[4], pzg[4], eg[4];
		double pG[4], thetaG[4], phiG[4];
		for (size_t i = 0; i < 4; i++) {
			pG[i] = mce[6+i];
			thetaG[i] = mctheta[6+i];
			phiG[i] = mcphi[6+i];

			pxg[i] = mcpx[6+i];
			pyg[i] = mcpy[6+i];
			pzg[i] = mcpz[6+i];
			eg[i] = pG[i];
		}
		double pxp = mcpx[0];
		double pyp = mcpy[0];
		double pzp = mcpz[0];
		double pp = sqrt(pow(pxp,2)+pow(pyp,2)+pow(pzp,2));
		double ep = mce[0];

		double p0px1 = pxg[0] + pxg[1];
		double p0py1 = pyg[0] + pyg[1];
		double p0pz1 = pzg[0] + pzg[1];

		double p0px2 = pxg[2] + pxg[3];
		double p0py2 = pyg[2] + pyg[3];
		double p0pz2 = pzg[2] + pzg[3];

		double kstpx = p0px1 + p0px2 + pxp;
		double kstpy = p0py1 + p0py2 + pyp;
		double kstpz = p0pz1 + p0pz2 + pzp;

		double kstp = sqrt(pow(kstpx, 2) + pow(kstpy, 2) + pow(kstpz, 2));
		double kspiInvMass = sqrt(pow(eg[0] + eg[1] + eg[2] + eg[3] + ep, 2) - pow(kstp, 2));

		return kspiInvMass;
	}

	pair<double,double> bwKstar(double s){
		s = s/1000./1000.;
		double msPi = 0.13957;
		double msK  = 0.493677;
		double msKst= 0.89166;
		double gKst = 0.0508;

		double w = pow( (1.-pow(msK+msPi,2)/s) * (1.-pow(msK-msPi,2)/s) / (1.-pow(msK+msPi,2)/pow(msKst,2)) / (1.-pow(msK-msPi,2)/pow(msKst,2)) ,3);
		if (w>=0)
			w = sqrt(w);
		else
			w = 0;

		pair<double,double> bwKstarD = make_pair(1. - s/pow(msKst,2), -gKst * s / pow(msKst,3) * w);
		pair<double,double> bwKstar;
		bwKstar.first = bwKstarD.first / (pow(bwKstarD.first,2) + pow(bwKstarD.second,2));
		bwKstar.second = -bwKstarD.second / (pow(bwKstarD.first,2) + pow(bwKstarD.second,2));
		return bwKstar;
	}
	pair<double,double> bwPhi(double s,double mphi, double gphi){
		s = s/1000./1000.;
		double msPhi = mphi;
		double msK  = 0.495;
		double msKst= 0.89166;
		double gPhi = gphi;

		double w = pow( (1.-pow(msKst+msK,2)/s) * (1.-pow(msKst-msK,2)/s) / (1.-pow(msKst+msK,2)/pow(msPhi,2)) / (1.-pow(msKst-msK,2)/pow(msPhi,2)) ,3);
		if (w>=0)
			w = sqrt(w);
		else
			w = 0;

		pair<double,double> bwPhiD = make_pair(1. - s/pow(msPhi,2), -gPhi * s / pow(msPhi,3) * w);
		pair<double,double> bwPhi;
		bwPhi.first = bwPhiD.first / (pow(bwPhiD.first,2) + pow(bwPhiD.second,2));
		bwPhi.second = -bwPhiD.second / (pow(bwPhiD.first,2) + pow(bwPhiD.second,2));
		
		return bwPhi;
	}
	pair<double,double> bwRho(double s, double mrho, double grho){
		s = s/1000./1000.;
		double msRho = mrho;
		double msK  = 0.495;
		double msKst= 0.89166;
		double gRho = grho;

		double w = pow( (1.-pow(msKst+msK,2)/s) * (1.-pow(msKst-msK,2)/s) / (1.-pow(msKst+msK,2)/pow(msRho,2)) / (1.-pow(msKst-msK,2)/pow(msRho,2)) ,3);
		if (w>=0)
			w = sqrt(w);
		else
			w = 0;

		pair<double,double> bwRhoD = make_pair(1. - s/pow(msRho,2), -gRho * s / pow(msRho,3) * w);
		pair<double,double> bwRho;
		bwRho.first = bwRhoD.first / (pow(bwRhoD.first,2) + pow(bwRhoD.second,2));
		bwRho.second = -bwRhoD.second / (pow(bwRhoD.first,2) + pow(bwRhoD.second,2));
		return bwRho;
	}


	double recalculateWeight(){
		double previousWweight = 0;
		double currentWeight = 0;

		double couplPhi0 = 17.0923;
		double couplRho0 = couplPhi0*0.67;

		pair<double,double> c1T;
		pair<double,double> c2T;
		pair<double,double> c0Phi;
		pair<double,double> c0Rho;
		pair<double,double> cfacPhi;
		pair<double,double> cfacRho;


		pair<double,double> bwKstarL = bwKstar(pow(kpiInvMassMod(),2));
		c1T.first = bwKstarL.first;
		c1T.second = bwKstarL.second;

		bwKstarL = bwKstar(pow(kspiInvMassMod(),2));
		c2T.first = bwKstarL.first;
		c2T.second = bwKstarL.second;

		pair<double,double> bwPhiL = bwPhi(4*beam*beam, 1.68, 0.15);
		c0Phi.first = bwPhiL.first * couplPhi0;
		c0Phi.second = bwPhiL.second * couplPhi0;
		cfacPhi.first = c0Phi.first*(c1T.first+c2T.first) - c0Phi.second*(c1T.second+c2T.second);
		cfacPhi.second = c0Phi.first*(c1T.second+c2T.second) + c0Phi.second*(c1T.first+c2T.first);

		pair<double,double> bwRhoL = bwRho(4*beam*beam, 1.465, 0.4);
		c0Rho.first = bwRhoL.first * couplRho0;
		c0Rho.second = bwRhoL.second * couplRho0;
		cfacRho.first = c0Rho.first*(c1T.first-c2T.first) - c0Rho.second*(c1T.second-c2T.second);
		cfacRho.second = c0Rho.first*(c1T.second-c2T.second) + c0Rho.second*(c1T.first-c2T.first);

		previousWweight = pow(cfacPhi.first + cfacRho.first,2) + pow(cfacPhi.second + cfacRho.second,2);

		couplPhi0 = 17.0923;
		couplRho0 = couplPhi0 * 0.67;

		bwKstarL = bwKstar(pow(kpiInvMassMod(),2));
		c1T.first = bwKstarL.first;
		c1T.second = bwKstarL.second;

		bwKstarL = bwKstar(pow(kspiInvMassMod(),2));
		c2T.first = bwKstarL.first;
		c2T.second = bwKstarL.second;

		bwPhiL = bwPhi(4*beam*beam, 1.68, 0.15);
		c0Phi.first = bwPhiL.first * couplPhi0;
		c0Phi.first = bwPhiL.first * couplPhi0;
		c0Phi.second = bwPhiL.second * couplPhi0;
		cfacPhi.first = c0Phi.first*(c1T.first+c2T.first) - c0Phi.second*(c1T.second+c2T.second);
		cfacPhi.second = c0Phi.first*(c1T.second+c2T.second) + c0Phi.second*(c1T.first+c2T.first);

		bwRhoL = bwRho(4*beam*beam, 1.465, 0.4);
		c0Rho.first = bwRhoL.first * couplRho0;
		c0Rho.second = bwRhoL.second * couplRho0;
		cfacRho.first = c0Rho.first*(c1T.first-c2T.first) - c0Rho.second*(c1T.second-c2T.second);
		cfacRho.second = c0Rho.first*(c1T.second-c2T.second) + c0Rho.second*(c1T.first-c2T.first);

		currentWeight = pow(cfacPhi.first + cfacRho.first,2) + pow(cfacPhi.second + cfacRho.second,2);


		if(previousWweight!=0){
			//if (currentWeight/previousWweight > 10)
			//	cout << currentWeight/previousWweight << endl;
			return currentWeight/previousWweight;
		}
		else{
			return 0;
		}
	}


	//ks0 inv mass calculation based on kinfit parameters
	double ks0InvMass(int kfvers) {
		double pxg[4], pyg[4], pzg[4], eg[4];
		double pG[4], thetaG[4], phiG[4];
		for (size_t i = 0; i < 4; i++) {
			pG[i] = kfvers < 1.5 ? ppkf1[i + 2] : ppkf2[i + 2];
			thetaG[i] = kfvers < 1.5 ? thetakf1[i + 2] : thetakf2[i + 2];
			phiG[i] = kfvers < 1.5 ? phikf1[i + 2] : phikf2[i + 2];

			pxg[i] = pG[i] * sin(thetaG[i]) * cos(phiG[i]);
			pyg[i] = pG[i] * sin(thetaG[i]) * sin(phiG[i]);
			pzg[i] = pG[i] * cos(thetaG[i]);
			eg[i] = pG[i];
		}

		double p0px1 = pxg[0] + pxg[1];
		double p0py1 = pyg[0] + pyg[1];
		double p0pz1 = pzg[0] + pzg[1];
		//double p0p1 = sqrt(pow(p0px1, 2) + pow(p0py1, 2) + pow(p0pz1, 2));

		double p0px2 = pxg[2] + pxg[3];
		double p0py2 = pyg[2] + pyg[3];
		double p0pz2 = pzg[2] + pzg[3];
		//double p0p2 = sqrt(pow(p0px2, 2) + pow(p0py2, 2) + pow(p0pz2, 2));

		double k0px = p0px1 + p0px2;
		double k0py = p0py1 + p0py2;
		double k0pz = p0pz1 + p0pz2;
		double k0p = sqrt(pow(k0px, 2) + pow(k0py, 2) + pow(k0pz, 2));
		double ks0invmass = sqrt(pow(eg[0] + eg[1] + eg[2] + eg[3], 2) - pow(k0p, 2));

		return ks0invmass;
	}
	//energy deposition of all detected particles except 6 defined by kinFit, used to exclude multi pion final states
	double eGammas(int kfvers) {
		double egammas = 0;
		double ip[4];
		for (size_t j = 0; j < 4; j++){
			ip[j] = kfvers < 1.5 ? ipkf1[2+j] - 1 : ipkf2[2+j] - 1;
		}
		for (size_t i = 0; i < (nn + nc); i++){
			if (charge[i]==0 && i != ip[0] && i != ip[1] && i != ip[2] && i != ip[3] )
				egammas += energy[i];
		}
		egammas = egammas / (2. * beam);
		return egammas;
	}
	//energy deposition of charged kaon and pion (indicies from kinfit), used to exclude processes like eta gamma
	double eKPi(int kfvers) {
		double energyK  = kfvers < 1.5 ? energy[ipkf1[0] - 1] : energy[ipkf2[0] - 1];
		double energyPi = kfvers < 1.5 ? energy[ipkf1[1] - 1] : energy[ipkf2[1] - 1];
		return (energyK + energyPi) / (2. * beam);
	}

	double anglekspi(int kfvers, int p0ind) {
		double pxg[4], pyg[4], pzg[4], eg[4];
		double pG[4], thetaG[4], phiG[4];
		for (size_t i = 0; i < 4; i++) {
			pG[i] = kfvers < 1.5 ? ppkf1[i + 2] : ppkf2[i + 2];
			thetaG[i] = kfvers < 1.5 ? thetakf1[i + 2] : thetakf2[i + 2];
			phiG[i] = kfvers < 1.5 ? phikf1[i + 2] : phikf2[i + 2];

			pxg[i] = pG[i] * sin(thetaG[i]) * cos(phiG[i]);
			pyg[i] = pG[i] * sin(thetaG[i]) * sin(phiG[i]);
			pzg[i] = pG[i] * cos(thetaG[i]);
			eg[i] = pG[i];
		}
		/*if(nmc>10)
			return 1000;
		for (size_t i = 0; i < 4; i++) {
			pG[i] = mce[6+i];
			thetaG[i] = mctheta[6+i];
			phiG[i] = mcphi[6+i];

			pxg[i] = pG[i] * sin(thetaG[i]) * cos(phiG[i]);
			pyg[i] = pG[i] * sin(thetaG[i]) * sin(phiG[i]);
			pzg[i] = pG[i] * cos(thetaG[i]);
			eg[i] = pG[i];
		}*/

		double p0px1 = pxg[0] + pxg[1];
		double p0py1 = pyg[0] + pyg[1];
		double p0pz1 = pzg[0] + pzg[1];
		double p0p1 = sqrt(pow(p0px1, 2) + pow(p0py1, 2) + pow(p0pz1, 2));

		double p0px2 = pxg[2] + pxg[3];
		double p0py2 = pyg[2] + pyg[3];
		double p0pz2 = pzg[2] + pzg[3];
		double p0p2 = sqrt(pow(p0px2, 2) + pow(p0py2, 2) + pow(p0pz2, 2));

		/*double p0px1 = mcpx[5];
		double p0py1 = mcpy[5];
		double p0pz1 = mcpz[5];
		double p0p1 = sqrt(pow(p0px1, 2) + pow(p0py1, 2) + pow(p0pz1, 2));

		double p0px2 = mcpx[4];
		double p0py2 = mcpy[4];
		double p0pz2 = mcpz[4];
		double p0p2 = sqrt(pow(p0px2, 2) + pow(p0py2, 2) + pow(p0pz2, 2));*/

		double k0px = p0px1 + p0px2;
		double k0py = p0py1 + p0py2;
		double k0pz = p0pz1 + p0pz2;
		double k0p = sqrt(pow(k0px, 2) + pow(k0py, 2) + pow(k0pz, 2));
	
		double pcol1 = (p0px1*k0px + p0py1*k0py + p0pz1*k0pz)/k0p;
		double ptrans1 = pow(p0p1,2) - pow(pcol1,2);
		double pcol2 = (p0px2*k0px + p0py2*k0py + p0pz2*k0pz)/k0p;
		double ptrans2 = pow(p0p2,2) - pow(pcol2,2);
		double mp0 = 135.;
		double enp01 = eg[0] + eg[1];
		double enp02 = eg[2] + eg[3];
		double enk01 = eg[0] + eg[1] + eg[2] + eg[3];

		double beta = k0p/enk01;
		double gamma = 1./sqrt(1.-pow(beta,2));

		double pcol11 = gamma * (pcol1 - beta * enp01);
		double pcol21 = gamma * (pcol2 - beta * enp02);
		double ptrans11 = ptrans1;
		double ptrans21 = ptrans2;

		double angle1 = 90;
		if(pcol11 != 0)
			angle1 = atan(ptrans11/pcol11)*180./3.14159;
		if(angle1<0)
			angle1 = angle1 + 180.;
		angle1 = pcol11/sqrt(pow(pcol11,2)+ptrans11);
		double angle2 = 90;
		if(pcol21 != 0)
			angle2 = atan(ptrans21/pcol21)*180./3.14159;
		if(angle2<0)
			angle2 = angle2 + 180.;
		angle2 = pcol21/sqrt(pow(pcol21,2)+ptrans21);
		
		//return acos(pcol1/k0p)*180./3.14159;
		if(p0ind == 1)
			return angle1;
		if(p0ind == 2)
			return angle2;
	}

	double gammagammaInvMasskf(int kfvers) {
		double pxg[2], pyg[2], pzg[2], eg[2];
		double pG[2], thetaG[2], phiG[2];
		for (size_t i = 0; i < 2; i++) {
			pG[i] = kfvers < 1.5 ? ppkf1[i + 2] : ppkf2[i + 2];
			thetaG[i] = kfvers < 1.5 ? thetakf1[i + 2] : thetakf2[i + 2];
			phiG[i] = kfvers < 1.5 ? phikf1[i + 2] : phikf2[i + 2];

			pxg[i] = pG[i] * sin(thetaG[i]) * cos(phiG[i]);
			pyg[i] = pG[i] * sin(thetaG[i]) * sin(phiG[i]);
			pzg[i] = pG[i] * cos(thetaG[i]);
			eg[i] = pG[i];
		}
		
		double p0px = pxg[0] + pxg[1];
		double p0py = pyg[0] + pyg[1];
		double p0pz = pzg[0] + pzg[1];
		double ep0 = eg[0] + eg[1];
		
		return sqrt(pow(ep0, 2) - pow(p0px, 2) - pow(p0py, 2) - pow(p0pz, 2));
	}
	
	double gammagammaInvMass(pair<int, int> ind) {
		double pxg[2], pyg[2], pzg[2], eg[2];
		double pG[2], thetaG[2], phiG[2];

		pG[0] = ppkf1[ind.first];
		eg[0] = energy[ind.first];
		pxg[0] = pG[0] * sin(theta[ind.first]) * cos(phi[ind.first]);
		pyg[0] = pG[0] * sin(theta[ind.first]) * sin(phi[ind.first]);
		pzg[0] = pG[0] * cos(theta[ind.first]);
		pG[1] = energy[ind.second];
		eg[1] = energy[ind.second];
		pxg[1] = pG[1] * sin(theta[ind.second]) * cos(phi[ind.second]);
		pyg[1] = pG[1] * sin(theta[ind.second]) * sin(phi[ind.second]);
		pzg[1] = pG[1] * cos(theta[ind.second]);

		double p0px = pxg[0] + pxg[1];
		double p0py = pyg[0] + pyg[1];
		double p0pz = pzg[0] + pzg[1];
		double ep0 = eg[0] + eg[1];

		return sqrt(pow(ep0, 2) - pow(p0px, 2) - pow(p0py, 2) - pow(p0pz, 2));
	}

	double extraPartInvMass(){
		vector< pair<int, int> > gammasInd;
		vector<int> excessiveGInd;
		for(size_t i = 0; i < nc+nn; i++){
			if(charge[i]==0 && x2kf1<20 && ipkf1[2]-1 != i && ipkf1[3]-1 != i && ipkf1[4]-1 != i && ipkf1[5]-1 != i){
				excessiveGInd.push_back((int)i);
			}
			else if(charge[i]==0 && x2kf2<20 && ipkf2[2]-1 != i && ipkf2[3]-1 != i && ipkf2[4]-1 != i && ipkf2[5]-1 != i){
				excessiveGInd.push_back((int)i);
			}
		}
		size_t imax = max(min(10,(int)excessiveGInd.size()) - 1,0);
		for(size_t i = 0; i < imax; i++){
			for(size_t j = i+1; j < min(10,(int)excessiveGInd.size()); j++){
				gammasInd.push_back(make_pair(excessiveGInd[i],excessiveGInd[j]));
				//cout << excessiveGInd.size() <<endl;
			}
		}
		//for(size_t i = 0; i < min(gammasInd.size(),1); i++)
		if (gammasInd.size() > 0)
			return gammagammaInvMass(gammasInd[0]);
		else if (gammasInd.size()==0){
			return 0;
		}
		return 0;
	}

	int mixFix();
	bool isITKaonOrPion(int ipkf, int kfx, int dyn);
	double itIsKaon(int ipkf, int kfx);
	double itIsPion(int ipkf, int kfx);

	//some conditions on event parameters, standing "before" kinfit
	bool central(int kfid, int kfvers) {
		double d = kfvers < 1.5 ? d0[ipkf1[kfid] - 1] : d0[ipkf2[kfid] - 1];
		double z = kfvers < 1.5 ? z0[ipkf1[kfid] - 1] : z0[ipkf2[kfid] - 1];
		return (fabs(d) < 0.5 && fabs(z) < 8);
	}
	bool unionVertex(int ind1, int ind2, int kfvers) {
		double z01 = kfvers < 1.5 ? z0[ipkf1[ind1] - 1] : z0[ipkf2[ind1] - 1];
		double z02 = kfvers < 1.5 ? z0[ipkf1[ind2] - 1] : z0[ipkf2[ind2] - 1];
		return fabs(z01 - z02) < 2;
	}
	bool thetaxBetween(int ind, double theta1, double theta2, int kfvers) {
		double thetaI = kfvers < 1.5 ? theta[ipkf1[ind] - 1] : theta[ipkf2[ind] - 1];
		return (thetaI > theta1 * 3.14159 / 180. && thetaI < theta2 * 3.14159 / 180.);
	}
	bool egmtht(int kfvers){
		double eg[4];
		for (size_t i = 0; i < 4; i++)
			eg[i] = kfvers < 1.5 ? ppkf1[i + 2] : ppkf2[i + 2];
		return eg[0]>40 && eg[1]>40 && eg[2]>40 && eg[3]>40;
	}

	void changeOn4pi(){
		for(size_t i = 0; i < 6; i++){
			ipkf1[i] = ipkf3[i];
			ipkf2[i] = ipkf3[i];
			ppkf1[i] = ppkf3[i];
			ppkf2[i] = ppkf3[i];
			phikf1[i] = phikf3[i];
			phikf2[i] = phikf3[i];
			thetakf1[i] = phikf3[i];
			thetakf2[i] = phikf3[i];
		}
		x2kf1 = x2kf3;
		x2kf2 = x2kf3;
		x2ikf1 = x2kf3;
		x2ikf2 = x2kf3;
	}

	vector<bool> formConditions(int model, vector<size_t> ommit){
		vector<bool> temp;
		switch(model){
			case 1:
				temp.push_back(x2kf1 < 20 && x2kf1 < x2kf3+3);
				temp.push_back(central(0, 1));
				temp.push_back(central(1, 1));
				temp.push_back(unionVertex(0, 1, 1) );
				//temp.push_back(thetaxBetween(0, 45, 135, 1) && thetaxBetween(1, 45, 135, 1));
				temp.push_back(eGammas(1) < 0.6 * pow(500./beam,3));
				//temp.push_back(eKPi(1) > 0.1);
				temp.push_back(ppkf1[0] > 0 && ppkf1[1] > 0);
				temp.push_back(nn <=7 && nc <=3);
				//temp.push_back(tempdExn[0]>=0.001 && tempdExn[1]>=0.001 && tempdExn[0]<15.0 && tempdExn[1]<15.0);
				//temp.push_back(extraPartInvMass() < 110 || extraPartInvMass() > 160);
				//temp.push_back(tempdExn[0]<=3.5 && tempdExn[1]<=3.5);
			break;
			case 2:
				temp.push_back(x2kf2 < 20 && x2kf2 < x2kf3 + 3);
				temp.push_back(central(0, 2));
				temp.push_back(central(1, 2));
				temp.push_back(unionVertex(0, 1, 2) );
				//temp.push_back(thetaxBetween(0, 45, 135, 2) && thetaxBetween(1, 45, 135, 2));
				temp.push_back(eGammas(2) < 0.6 * pow(500./beam,3));
				//temp.push_back(eKPi(2) > 0.1);
				temp.push_back(ppkf2[0] > 0 && ppkf2[1] > 0);
				temp.push_back(nn <=7 && nc <=3);
				//temp.push_back(tempdExn[2]>=0.001 && tempdExn[3]>=0.001 && tempdExn[2]<15.0 && tempdExn[3]<15.0);
				//temp.push_back(extraPartInvMass() < 110 || extraPartInvMass() > 160);
				//temp.push_back(tempdExn[2]<=3.5 && tempdExn[3]<=3.5);
			break;
		}
		for(size_t i = 0; i<ommit.size();i++){
			temp[ommit[i]] = true;
		}
		return temp;
	}

	vector<bool> formConditions(int model){
		vector<size_t> abcd;
		abcd.clear();
		return formConditions(model, abcd);
	}

	bool goodEventReturnResult(){
		bool result = true;
		for(size_t i = 0; i < conditions.size(); i++)
			result = result && conditions[i];
		return result;
	}
	bool goodKkpEntry1Cut(int ind){
		vector<size_t> temp;
		temp.push_back((int)ind);
		conditions = formConditions(1,temp);
		return goodEventReturnResult();
	}
	bool goodKkpEntry2Cut(int ind){
		vector<size_t> temp;
		temp.push_back((int)ind);
		conditions = formConditions(2,temp);
		return goodEventReturnResult();
	}

	bool goodKkpEntry1(){
		conditions = formConditions(1);
		//return goodKkpEntry1Cut(0);
		return goodEventReturnResult();
	}
	bool goodKkpEntry2(){
		conditions = formConditions(2);
		//goodKkpEntry2Cut(0);
		return goodEventReturnResult();
	}
	
	bool goodEventCheck(string regime){
		//normalizing aach and dedx
		kf10 = ppkf1[0] < 300 || region[ipkf1[0] - 1] == 1 && tempamp[0] < 30 && amplitudeE[ipkf1[0] - 1] < 1000 && amplitudeE[ipkf1[0] - 1] > 0.01;//&& amplitude[ipkf1[0] - 1] > 0.01;
		kf11 = ppkf1[1] < 300 || region[ipkf1[1] - 1] == 1 && tempamp[1] < 30 && amplitudeE[ipkf1[1] - 1] < 1000 && amplitudeE[ipkf1[1] - 1] > 0.01;//&& amplitude[ipkf1[1] - 1] > 0.01;
		kf20 = ppkf2[0] < 300 || region[ipkf2[0] - 1] == 1 && tempamp[2] < 30 && amplitudeE[ipkf2[0] - 1] < 1000 && amplitudeE[ipkf2[0] - 1] > 0.01;//&& amplitude[ipkf2[0] - 1] > 0.01;
		kf21 = ppkf2[1] < 300 || region[ipkf2[1] - 1] == 1 && tempamp[3] < 30 && amplitudeE[ipkf2[1] - 1] < 1000 && amplitudeE[ipkf2[1] - 1] > 0.01;//&& amplitude[ipkf2[1] - 1] > 0.01;
		return ((goodKkpEntry1() && kf10 && kf11) || (goodKkpEntry2() && kf20 && kf21)) /*&& eton > 0.4 && eton < 0.85*/;
	}


	void someAnalisis();
	void effEnergy();
	void birdfix();
	void drawdEvsp(vector<pair<double, double> > eranges);
	void landaufit();
	void genTableML(char *add, char *outputFile, int appYN);
	void genTableMLExp(char *add, char *outdistr, char *txtfile, int ind, vector<pair<double, double> > eranges, string regime);
	vector<TH1*> imvmass(char* add, char* name, char* title, int nbins, double lb, double hb, vector<pair<double, double> > eranges, string regime);
	void calculateEffAndNoE();
	void invmass2part();
};

pair<size_t,size_t> findclosest(double en1, vector<double> en2){
	double mindistl = 100000;
	double mindistr = 100000;
	size_t mindistIndl = 500;
	size_t mindistIndr = 500;
	for(size_t i = 0; i < en2.size(); i++){
		if(fabs(en1-en2[i])<mindistl && en2[i] <= en1){
			mindistl = fabs(en1-en2[i]);
			mindistIndl = i;
		}
		if(fabs(en1-en2[i])<mindistr && en2[i] >= en1){
			mindistr = fabs(en1-en2[i]);
			mindistIndr = i;
		}
	}
	if (mindistIndl == 500)
		mindistIndl = mindistIndr;
	if (mindistIndr == 500)
		mindistIndr = mindistIndl;
	return make_pair(mindistIndl,mindistIndr);
}

double energygaps[13] = {1200, 1395, 1499, 1523, 1590, 1655, 1690, 1755, 1795, 1805, 1875, 1945, 2005};//1499 1523
//double energygaps[11] = {1200, 1355, 1495, 1605, 1705, 1745, 1795, 1820, 1865, 1882, 2007};

//usefull function for pion-kaon separation in track chamber
double betaBloch(double p, double m) {
	return 0.1* (p * p + m * m) / (p * p) * (log((p * p) / (m * m)) - (p * p) / (p * p + m * m) + 4.8); //5.3 f.e. to cut pi, 600
}
double pionLine(double p) {
	return 600. * (p*p + 139.6*139.6)/(p*p) * (log((p*p)/(139.6*139.6)) - (p*p)/(p*p + 139.6*139.6) + 6.5);
}
double sigmapion(double p) {
	return p; // have to find sigma dependence on momentum by fitting in some ranges with mc data usage
}
double kaonLine(double p) {
	return 600. * (p*p + 493.7*493.7)/(p*p) * (log((p*p)/(493.7*493.7)) - (p*p)/(p*p + 493.7*493.7) + 6.5);
}
double sigmakaon(double p) {
	return p; // have to find sigma dependence on momentum by fitting in some ranges with mc data usage
}
double meanLine(double p) {
	return betaBloch(p,355); //375 to cut pi
}

double sigma(float p, float mean, float width) {
	return 1./(1.+exp(-(p-mean)/width));
}

double eff(int pcode, float p, float ampE) {  //pcode - 0 if kaon hypithesis (poisson distribution), 1 if pion
	if (pcode == 0)
		return 1. * (0.2+0.1*ampE);
	if (pcode == 1)
		return sigma(p, 350., 100.) * (0.2 + 0.1 * ampE);
}


//First try to fix mixing charged particles of kinfit, is not used now, because I have 2 models and separation function,which is described right below
int kkpDistr::mixFix() {
	if (x2ikf1 < 100)
	{
		double pp2[2] = { sqrt(ppkf1[0] * ppkf1[0] - 139.6 * 139.6 + 493.7 * 493.7) , sqrt(ppkf1[1] * ppkf1[1] + 139.6 * 139.6 - 493.7 * 493.7) };
		if (fabs(dExs[ipkf1[1] - 1] - pionLine(ppkf1[1])) * fabs(dExs[ipkf1[1] - 1] - pionLine(ppkf1[1])) +
			fabs(dExs[ipkf1[0] - 1] - kaonLine(ppkf1[0])) * fabs(dExs[ipkf1[0] - 1] - kaonLine(ppkf1[0])) >
			fabs(dExs[ipkf1[1] - 1] - kaonLine(pp2[1])) * fabs(dExs[ipkf1[1] - 1] - kaonLine(pp2[1])) +
			fabs(dExs[ipkf1[0] - 1] - pionLine(pp2[0])) * fabs(dExs[ipkf1[0] - 1] - pionLine(pp2[0])))
		{//почему я сравниваю не поменяный дедх, а потом пишу поменяный? бред какой то
			ppkf1[0] = pp2[1]; ppkf1[1] = pp2[0];
			float dexs = dExs[ipkf1[1] - 1];
			dExs[ipkf1[1] - 1] = dExs[ipkf1[0] - 1];
			dExs[ipkf1[0] - 1] = dexs;
			return 1;
		}
	}
	return 0;
}

//differenciate charged particles, ipkf - particle index (kinfit), kfx - kinfit model
//если каон (выше линии в дедх или плохой сигнал счетчика), то True, если пион то False
bool kkpDistr::isITKaonOrPion(int ipkf, int kfx, int dyn) {
	int ip = ipkf1[ipkf] - 1;
	float p= ppkf1[ipkf];
	if (kfx == 2) {
		ip = ipkf2[ipkf] - 1;
		p  = ppkf2[ipkf];
	}
	float amp; amp = amplitude[ip];
	float dE; dE = dExnC[ip];
	if (dyn == 1){
		dE = dExn[ip];
		dE = dE/605.;
	}
	int reg; reg = region[ip];
	return (p < 350 && dE > meanLine(p)) ||
		((p >= 350) && (amp > 500 || amp < 0.25) && reg==1);
}

//functions for probability estimations of particle being kaon or pion (first try), have to read some papers to do it correctly.
/*int itIsKaon(int ipkf, int kfx) {
	int ip = ipkf1[ipkf] - 1;
	float amplitudeE;
	float p = ppkf1[ipkf];
	if (kfx == 2) {
		ip = ipkf2[ipkf] - 1;
		p = ppkf2[ipkf];
	}
	float amp; amp = amplitude[ip];
	float dE; dE = dExs[ip];
	int reg; reg = region[ip];

	int probab = (kaonLine(p) - meanLine(p) - fabs(dE - kaonLine(p)))/(kaonLine(p) - meanLine(p));
	if (p >= 350) {
		probab = pow(amplitudeE, amp) * exp(-amplitudeE);
		float amptemp = amp;
		while (amptemp >= 1)
			probab = probab / amptemp;
		probab 
	}
	if (probab >= 0)
		return probab;
	return 0;
}*/
/*int itIsPion(int ipkf, int kfx) {
	int ip = ipkf1[ipkf] - 1;
	float p = ppkf1[ipkf];
	if (kfx == 2) {
		ip = ipkf2[ipkf] - 1;
		p = ppkf2[ipkf];
	}
	float amp; amp = amplitude[ip];
	float dE; dE = dExs[ip];
	int reg; reg = region[ip];

	int probab = (meanLine(p) - pionLine(p) - fabs(dE - pionLine(p))) / (meanLine(p) - pionLine(p));
	if (probab >= 0)
		return probab;
	return 0;
}*/

//functions for probability estimations of particle being kaon or pion (second try), have to read some papers to do it correctly.
double kkpDistr::itIsKaon(int ipkf, int kfx) {   // index from kinfit and version of kinfit, returns ratio of probabilities in two hypotheses
	int ip = ipkf1[ipkf] - 1;
	float p = ppkf1[ipkf];
	if (kfx == 2) {
		ip = ipkf2[ipkf] - 1;
		p = ppkf2[ipkf];
	}
	float amp; amp = amplitude[ip];
	float ampE; ampE = amplitudeE[ip];
	float dE; dE = dExs[ip];
	int reg; reg = (int)region[ip];
	double probab = 0;
	if (p < 350) {
		probab = TMath::Landau(dE, kaonLine(p), sigmakaon(p)) / TMath::Landau(dE, pionLine(p), sigmapion(p));
	}
	if (p >= 350) {
		if (amp > 0.2){
			probab = eff(0,p,ampE)/eff(1,p,ampE);
		}
		if (amp <= 0.2) {
			probab =(1-eff(0,p,ampE))/(1-eff(1,p,ampE));
		}
	}
	if (probab >= 0)
		return probab;
	return 0;
}
double kkpDistr::itIsPion(int ipkf, int kfx) {
	int ip = ipkf1[ipkf] - 1;
	float p = ppkf1[ipkf];
	if (kfx == 2) {
		ip = ipkf2[ipkf] - 1;
		p = ppkf2[ipkf];
	}
	float amp; amp = amplitude[ip];
	float ampE; ampE = amplitudeE[ip];
	float dE; dE = dExs[ip];
	int reg; reg = region[ip];
	double probab = 0;
	if (p < 350) {
		probab = TMath::Landau(dE, pionLine(p), sigmapion(p)) / TMath::Landau(dE, kaonLine(p), sigmakaon(p));
	}
	if (p >= 350) {
		if (amp > 0.2) {
			probab = eff(1, p, ampE) / eff(0, p, ampE);
		}
		if (amp <= 0.2) {
			probab = (1 - eff(1, p, ampE)) / (1 - eff(0, p, ampE));
		}
	}
	if (probab >= 0)
		return probab;
	return 0;
}

//trying to create a class, that holds vector of histograms and can fill them with some cuts from some data, 
//but it returns segmentation violation depending on certain amount of branches used
/*class hwf {
	//const int n;
	//object = array of hists/profiles (or vector better?) пока 1 гистограмма
public:
	TH1F* hist;
	//object 2 = fuctions to fill (пока не буду делать этого)

	//some constructor
	hwf(const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xup) {
		*hist = TH1F(name, title, nbins, xlow, xup);
	}

	//filling function
	void fillhists(char* data) {

		TChain* chain = new TChain("t1");
		chain->Add(data);
		const Long64_t entries = chain->GetEntries();
		cout << entries << endl;

		chain->SetBranchAddress("beam", &beam);
		chain->SetBranchAddress("eton", &eton);
		chain->SetBranchAddress("nc", &nc);
		chain->SetBranchAddress("z0", &z0);
		chain->SetBranchAddress("d0", &d0);
		chain->SetBranchAddress("phi", &phi);
		chain->SetBranchAddress("phis", &phis);
		chain->SetBranchAddress("theta", &theta);
		chain->SetBranchAddress("region", &region);

		chain->SetBranchAddress("x2ikf1", &x2ikf1);
		chain->SetBranchAddress("ipkf1", &ipkf1);
		chain->SetBranchAddress("ppkf1", &ppkf1);
		chain->SetBranchAddress("phikf1", &phikf1);
		chain->SetBranchAddress("thetakf1", &thetakf1);

		chain->SetBranchAddress("amplitude", &amplitude);
		chain->SetBranchAddress("dExs", &dExs);

		chain->SetBranchAddress("x2ikf2", &x2ikf2);
		chain->SetBranchAddress("ipkf2", &ipkf2);
		chain->SetBranchAddress("ppkf2", &ppkf2);
		chain->SetBranchAddress("phikf2", &phikf2);
		chain->SetBranchAddress("thetakf2", &thetakf2);

	}

};*/

//function to select kkp events, maybe will be usefull for experiment of 2017-2020, but probably not, and it should be updatad seriously
void copy() {
	TString str = "a";
	TString elem;
	vector<TString> ans;
	TString basedir("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/*.root");
	TString dir("/work/users/kladov/snd2k/R007-001/kkpi/selectedkkpiCharged/");
	TString files = gSystem->GetFromPipe("ls " + basedir);

	Ssiz_t from = 0;

	while (files.Tokenize(elem, from, "\n")) {
		ans.push_back(elem);
	}
	int vector_size = ans.size();
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
		TFile* ouput = TFile::Open(dir + "true__" + tokens[8], "RECREATE");
		ouput->cd();
		TTree* selectedTree = originalTree->CopyTree("(cosm == 0) && (nc == 4) && (nn==0)");
		selectedTree->Print();
		selectedTree->Write();
		ouput->Close();
		filecomb->Close();
	}
}

//scan all processes (listed in directory), use selection criteria and store some histogramm in a file (for example, inv mass distributions)
void kkpDistr::someAnalisis() {
	vector<TString> ans;
	//what directory
	TString files = gSystem->GetFromPipe("ls /work/users/kladov/snd2k/R007-001/output/ntuples/all_KsNF/*");
	TString elem;
	Ssiz_t from0 = 0;
	while (files.Tokenize(elem, from0, "\n")) {
		ans.push_back(elem);
	}

	//where to store
	TFile* MyFile1 = new TFile("/work/users/kladov/snd2k/R007-001/kkpi/invmassfonovWC.root", "RECREATE"); 
	int choice = 0;
	int mycount = 0;
	unsigned int vector_size = ans.size();
	for (unsigned int i = 0; i < vector_size; i++) {
		TChain chain("t1");
		TString str = ans[i].Copy();
		chain.Add(str);
		//TFile* filecomb = TFile::Open(str);
		const Long64_t entries = chain.GetEntries();
		cout << entries << endl;

		mycount = 0;
		{
			chain.SetBranchAddress("beam", &beam);
			chain.SetBranchAddress("eton", &eton);
			chain.SetBranchAddress("energy", &energy);
			chain.SetBranchAddress("nc", &nc);
			chain.SetBranchAddress("nn", &nn);
			chain.SetBranchAddress("z0", &z0);
			chain.SetBranchAddress("d0", &d0);
			chain.SetBranchAddress("x0", &x0);
			chain.SetBranchAddress("y0", &y0);
			chain.SetBranchAddress("phi", &phi);
			chain.SetBranchAddress("phis", &phis);
			chain.SetBranchAddress("theta", &theta);
			chain.SetBranchAddress("region", &region);
			chain.SetBranchAddress("amplitude", &amplitude);
			chain.SetBranchAddress("dExs", &dExs);

			chain.SetBranchAddress("x2ikf1", &x2ikf1);
			chain.SetBranchAddress("ipkf1", &ipkf1);
			chain.SetBranchAddress("ppkf1", &ppkf1);
			chain.SetBranchAddress("phikf1", &phikf1);
			chain.SetBranchAddress("thetakf1", &thetakf1);

			chain.SetBranchAddress("x2ikf2", &x2ikf2);
			chain.SetBranchAddress("ipkf2", &ipkf2);
			chain.SetBranchAddress("ppkf2", &ppkf2);
			chain.SetBranchAddress("phikf2", &phikf2);
			chain.SetBranchAddress("thetakf2", &thetakf2);
		}
		//TH1* h = new TH1F("histx", "x2afterallcond", 25, 0, 100);
		TH1* h = new TH1F("invks01", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265);

		for (int e = 0; e < entries; e++) {
			chain.GetEntry(e);
			choice = 0;
			//_! and < for process, !_ and > for background
			if (x2ikf1 < 50 && isITKaonOrPion(0, 1, 1) && !isITKaonOrPion(1, 1, 1))
				choice = 1;
			if (x2ikf2 < 50 && isITKaonOrPion(0, 2, 1) && !isITKaonOrPion(1, 2, 1) && (choice == 0 || x2ikf2 < x2ikf1))
				choice = 2;
			if (choice == 2) {
				x2ikf1 = x2ikf2;
				for (size_t li = 0; li < 6; li++) {
					ipkf1[li] = ipkf2[li];
					ppkf1[li] = ppkf2[li];
					phikf1[li] = phikf2[li];
					thetakf1[li] = thetakf2[li];
				}
			}
			double ks0invmass = ks0InvMass(1);
			double egammas = eGammas(1);
			double ekpi = eKPi(1);

			bool goodKkpEntry = (choice != 0 && central(0,1) && central(1,1) && unionVertex(0, 1, 1) && eton > 0.2 && eton < 1.0 && thetaxBetween(0, 30, 150, 1) && thetaxBetween(1, 30, 150, 1) && nc >= 2);
			if (goodKkpEntry && ppkf1[0] > 10 && ppkf1[1] > 10 /*&& (ks0invmass< 495. + 35 && ks0invmass > 495. - 35)*/ && egammas < 0.15 && ekpi>0.15) {
				//h->Fill(x2ikf1);
				h->Fill(ks0invmass);
				mycount++;
			}
		}
		MyFile1->cd();

		TString tok;
		Ssiz_t from1 = 0;
		vector<TString> tokens;

		while (str.Tokenize(tok, from1, "[/.]")) {
			tokens.push_back(tok);
		}
		cout << mycount << endl;
		if (mycount >= 20) {
			h->Write(tokens[8]);
		}
	}
	MyFile1->Close();
}

//diferent manipulations with efficiency, in main loop it has only selection criteria and filling efficiency profile with different x values
void kkpDistr::effEnergy() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp_KsNF/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/2etg/*.root");
	//chain1.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/all/4pi_wrc*");
	{
		chain.SetBranchAddress("beam", &beam);
		chain.SetBranchAddress("eton", &eton);
		chain.SetBranchAddress("energy", &energy);
		chain.SetBranchAddress("nc", &nc);
		chain.SetBranchAddress("nn", &nn);
		chain.SetBranchAddress("z0", &z0);
		chain.SetBranchAddress("d0", &d0);
		chain.SetBranchAddress("x0", &x0);
		chain.SetBranchAddress("y0", &y0);
		chain.SetBranchAddress("phi", &phi);
		chain.SetBranchAddress("phis", &phis);
		chain.SetBranchAddress("theta", &theta);
		chain.SetBranchAddress("region", &region);
		chain.SetBranchAddress("amplitude", &amplitude);
		chain.SetBranchAddress("dExs", &dExs);

		chain.SetBranchAddress("x2ikf1", &x2ikf1);
		chain.SetBranchAddress("ipkf1", &ipkf1);
		chain.SetBranchAddress("ppkf1", &ppkf1);
		chain.SetBranchAddress("phikf1", &phikf1);
		chain.SetBranchAddress("thetakf1", &thetakf1);

		chain.SetBranchAddress("x2ikf2", &x2ikf2);
		chain.SetBranchAddress("ipkf2", &ipkf2);
		chain.SetBranchAddress("ppkf2", &ppkf2);
		chain.SetBranchAddress("phikf2", &phikf2);
		chain.SetBranchAddress("thetakf2", &thetakf2);

		chain.SetBranchAddress("mctheta", &mctheta);
		chain.SetBranchAddress("mcphi", &mcphi);
		chain.SetBranchAddress("mce", &mce);
		chain.SetBranchAddress("mcpx", &mcpx);
		chain.SetBranchAddress("mcpy", &mcpy);
		chain.SetBranchAddress("mcpz", &mcpz);
		chain.SetBranchAddress("mcpdg", &mcpdg);
		chain.SetBranchAddress("nmc", &nmc);
	}
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	double bins[31];
	bins[0] = 0;
	for (size_t i = 1; i < 31; i++)
		bins[i] = 2 + i + bins[i-1];
	vector<TProfile*> prof;
	//TProfile* prof = new TProfile("effenergy", "eff vs energy",30, bins);
	for(size_t i = 0; i < 12; i++)
		prof.push_back(new TProfile(Form("effenergy%d",i), "eff vs energy", 30, bins));
	//TH1* h = new TH1F("x2kf1", "x2i kskp", 50, 0, 200);
	//TH1* h1 = new TH1F("x2kf1pi", "x2i 4pi", 50, 0, 200);
	int choice = 0;
	double counteff = 0;
	double countall = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//if ((2 * beam >= energygaps[4] && 2 * beam < energygaps[5])) {
		size_t erangel = 0;
		for (size_t i = 0; i < sizeof(energygaps) / sizeof(energygaps[0]); i++) {
			if ((2 * beam >= energygaps[i] && 2 * beam < energygaps[i + 1])) {
				erangel = i;
			}
		}
			countall += 1;
			choice = 0;
			if (x2ikf1 < 50 && isITKaonOrPion(0, 1, 1) && !isITKaonOrPion(1, 1, 1))
				choice = 1;
			if (x2ikf2 < 50 && isITKaonOrPion(0, 2, 1) && !isITKaonOrPion(1, 2, 1) && (choice == 0 || x2ikf2 < x2ikf1))
				choice = 2;

			if (choice == 2) {
				x2ikf1 = x2ikf2;
				for (size_t li = 0; li < 6; li++) {
					ipkf1[li] = ipkf2[li];
					ppkf1[li] = ppkf2[li];
					phikf1[li] = phikf2[li];
					thetakf1[li] = thetakf2[li];
				}
			}
			//double ks0invmass = ks0InvMass(1);
			double egammas = eGammas(1);
			double ekpi = eKPi(1);

			bool goodKkpEntry = (choice != 0 && central(0,1) && central(1,1) && unionVertex(0, 1, 1) && eton > 0.2 && eton < 1.0 && thetaxBetween(0, 30, 150, 1) && thetaxBetween(1, 30, 150, 1) && nc >= 2);
			if (goodKkpEntry && ppkf1[0] > 10 && ppkf1[1] > 10 && egammas < 0.15 && ekpi>0.15) {
				prof[erangel]->Fill(mce[1], 1); //for dependence on ISR photon
				//prof->Fill(2*beam, 1);
				counteff += 1;
			}
			else
				prof[erangel]->Fill(mce[1], 0);
		//}
		/*if (nmc >= 10)
			prof->Fill(2. * beam, 1);
		if (nmc < 9)
			prof->Fill(2. * beam, 0);*/
		/*if (x2ikf1 <= 50)
			prof->Fill(2. * beam, 1);
		if (x2ikf1 > 50)
			prof->Fill(2. * beam, 0);*/
	}

	//double en[12] = { 1300, 1438, 1509, 1550, 1625, 1678, 1725, 1765, 1800, 1854, 1910, 1975 };
	//double enerr[12] = { 80, 35, 10, 25, 25, 2.5, 15, 8, 3, 20, 25, 20 };
	double mean[12];
	double meanerr[12];

	//prof->Scale(countall/counteff);
	prof[0]->Scale(1./prof[0]->GetBinContent(1));
	prof[0]->SetTitle("efficiency multiplier vs ISR photon energy; Energy, MeV; multiplier");
	prof[0]->Draw();
	prof[0]->GetYaxis()->SetRangeUser(0, 1.2);
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	c->SetBorderMode(0);
	c->SetFrameBorderMode(0);
	c->SetGrid();
	c->Update();

	//fit with 2 lines of efficiency dependence on E ISR
	//TF1* f1 = new TF1("f1", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	//f1->SetParameters(50, 1, -0.003, -0.01);
	TF1* f1 = new TF1("f1", "[2]*(TMath::Erf(1./[0]*([1]-x))+[3])", 0, 600);
	//TF1* f1 = new TF1("f1", "1./(1+TMath::Exp((x-[1])/[0]))", 0, 300);
	f1->SetParameters(50, 150,1.1,1);
	prof[0]->Fit("f1", "", "", 0, 300);
	c->Update();
	vector<double> par(4);
	f1->GetParameters(&par[0]);
	mean[0] = par[3];
	meanerr[0] = f1->GetParError(3);
	cout << "parameters:\n" << par[0] << ", " << par[1] << ", " << par[2] << ", " << par[3] << endl;
	for (size_t i = 1; i < 12; i++) {
		prof[i]->Scale(1. / prof[i]->GetBinContent(1));
		f1->SetParameters(50, 150, 1.1,1);
		prof[i]->Fit("f1", "", "", 0, 600);
		prof[i]->SetLineColor(i);
		prof[i]->Draw("same");

		c->Update();
		cin.get();
		vector<double> par(4);
		f1->GetParameters(&par[0]);
		mean[i] = par[3];
		meanerr[i] = f1->GetParError(3);
		cout << "parameters:\n" << par[0] << ", " << par[1] << ", " << par[2] << ", " << par[3] << endl;
	}

	/*TGraphErrors* gr = new TGraphErrors(12, en, mean, enerr, meanerr);
	gr->SetTitle("TGraphErrors Example");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->Draw("AP");
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	c->SetBorderMode(0);
	c->SetFrameBorderMode(0);
	c->SetGrid();
	c->Update();

	TF1* f2 = new TF1("f2", "[0]+x*[1]", 1520, 2000);
	f2->SetParameters(-10,0.1);
	gr->Fit("f2", "", "", 1520, 2000);
	cout << f2->GetParameter(0) << "	" << f2->GetParameter(1) << endl;*/

	//all graphs in one for eff/Eisr by energy points
	/*TF1* f2 = new TF1("f2", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	TF1* f3 = new TF1("f3", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	TF1* f4 = new TF1("f4", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	TF1* f5 = new TF1("f5", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	TF1* f6 = new TF1("f6", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	TF1* f7 = new TF1("f7", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	TF1* f8 = new TF1("f8", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	TF1* f9 = new TF1("f9", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	TF1* f10 = new TF1("f10", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	TF1* f11 = new TF1("f11", "x<[0] ? [1]+[2]*(x-[0]) : [1]+[3]*(x-[0])", 0, 200);
	f2->SetParameters(35.2359, 0.966004, -0.00115771, -0.0113549);
	f3->SetParameters(28.75, 0.970494, -0.00138306, -0.00620383);
	f4->SetParameters(40.7501, 0.997323, -0.000304604, -0.00841768);
	f5->SetParameters(54.75, 0.994797, -0.000339034, -0.00767365);
	f6->SetParameters(54.75, 0.993215, -0.000577483, -0.00850181);
	f7->SetParameters(54.1692, 1.0039, -0.000771498, -0.00789003);
	f8->SetParameters(34.5002, 1.16736, 0.00329028, -0.00674025);
	f9->SetParameters(54.7501, 1.04613, -0.00124571, -0.00749069);
	f10->SetParameters(63.115, 1.11227, -0.00117442, -0.00708491);
	f11->SetParameters(44.1119, 1.33755, 0.00229626, -0.00761448);
	f2->SetLineColor(2);
	f3->SetLineColor(3);
	f4->SetLineColor(4);
	f5->SetLineColor(5);
	f6->SetLineColor(6);
	f7->SetLineColor(7);
	f8->SetLineColor(8);
	f9->SetLineColor(9);
	f10->SetLineColor(10);
	f11->SetLineColor(11);
	f2->Draw("same");
	f3->Draw("same");
	f4->Draw("same");
	f5->Draw("same");
	f6->Draw("same");
	f7->Draw("same");
	f8->Draw("same");
	f9->Draw("same");
	f10->Draw("same");
	f11->Draw("same");*/

	//cout << prof->Integral(0,20) << "	" << prof->Integral(21,100);

	//real efficiency
	/*double en[12] = { 1300, 1438, 1509, 1550, 1625, 1678, 1725, 1765, 1800, 1854, 1910, 1975 };
	double enerr[12] = { 80, 35, 10, 25, 25, 2.5, 15, 8, 3, 20, 25, 20 };
	double eff[12] = { 3.84, 3.51, 3.35, 3.27422, 3.12232, 2.9758, 2.9437, 2.60861, 2.74384, 2.79617, 2.58806, 2.43196 };
	double efferr[12] = { 0.064, 0.058, 0.056, 0.08, 0.0395681, 0.0546186, 0.0384114, 0.0511319, 0.0524544, 0.0334872, 0.029405, 0.0312301 };
	TGraphErrors* gr = new TGraphErrors(12, en, eff, enerr, efferr);
	gr->SetTitle("TGraphErrors Example");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->Draw("Psame");
	TF1* f1 = new TF1("f1", "[0]+[1]*x+[2]*x**2", 1200, 2010);
	gr->Fit("f1", "", "", 1200, 2010);*/
	//h1->SetLineColor(2);
	//h->DrawNormalized("",h->Integral()/entries);
	//h1->DrawNormalized("",h1->Integral()/entries1);
	//h1->Draw();
}

//now it only draws the theta distribution of dEdxs for some root files
void kkpDistr::birdfix() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/selected2011/*.root");
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	TProfile* p0 = new TProfile("devstheta","devstheta",50,0,3.1416);
	chain.SetBranchAddress("dExs", &dExs);
	chain.SetBranchAddress("theta", &theta);
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		p0->Fill(theta[0], dExs[0]);
		p0->Fill(theta[1], dExs[1]);
	}
	p0->Draw();
}

double invmass2p(double e1, double e2, double p1x, double p2x, double p1y, double p2y, double p1z, double p2z){
	return sqrt(pow(e1+e2,2) - pow(p1x+p2x,2) - pow(p1y+p2y,2) - pow(p1z+p2z,2) ); 
}

void kkpDistr::invmass2part(){
	TChain chain("t1");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2017/kkp_KsNF/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2017/*.root");
	chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/4pi/2011/*x.root");
	{
		chain.SetBranchAddress("beam", &beam);
		chain.SetBranchAddress("eton", &eton);
		chain.SetBranchAddress("energy", &energy);
		chain.SetBranchAddress("nc", &nc);
		chain.SetBranchAddress("nn", &nn);
		chain.SetBranchAddress("z0", &z0);
		chain.SetBranchAddress("d0", &d0);
		chain.SetBranchAddress("x0", &x0);
		chain.SetBranchAddress("y0", &y0);
		chain.SetBranchAddress("phi", &phi);
		chain.SetBranchAddress("phis", &phis);
		chain.SetBranchAddress("theta", &theta);
		chain.SetBranchAddress("region", &region);
		chain.SetBranchAddress("amplitude", &amplitude);
		chain.SetBranchAddress("amplitudeE", &amplitudeE);
		chain.SetBranchAddress("dExs", &dExs);
		chain.SetBranchAddress("dExnC", &dExnC);
		chain.SetBranchAddress("dExn", &dExn);

		chain.SetBranchAddress("x2ikf1", &x2ikf1);
		chain.SetBranchAddress("ipkf1", &ipkf1);
		chain.SetBranchAddress("ppkf1", &ppkf1);
		chain.SetBranchAddress("phikf1", &phikf1);
		chain.SetBranchAddress("thetakf1", &thetakf1);

		chain.SetBranchAddress("x2ikf2", &x2ikf2);
		chain.SetBranchAddress("ipkf2", &ipkf2);
		chain.SetBranchAddress("ppkf2", &ppkf2);
		chain.SetBranchAddress("phikf2", &phikf2);
		chain.SetBranchAddress("thetakf2", &thetakf2);

		chain.SetBranchAddress("mctheta", &mctheta);
		chain.SetBranchAddress("mcphi", &mcphi);
		chain.SetBranchAddress("mce", &mce);
		chain.SetBranchAddress("mcpx", &mcpx);
		chain.SetBranchAddress("mcpy", &mcpy);
		chain.SetBranchAddress("mcpz", &mcpz);
		chain.SetBranchAddress("mcpdg", &mcpdg);
	}
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	int choice = 0;
	TH1* h0 = new TH1F("kspi","ks + pi inv mass distr",100,0,2000);
	TH1* h01 = new TH1F("kpi","k + pi inv mass distr",100,0,2000);
	//TH2* h1 = new TH2D("devspp","dEdx vs p,",200,0,600,400,-50,350);
	TH1* h2 = new TH1F("ach","amp ch, aa",100,-0.5,4.5);
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		choice = 0;
		int de = 1;
		if (x2ikf1 < 50 && isITKaonOrPion(0, 1, de) && !isITKaonOrPion(1, 1, de) && ppkf1[0]>400 && ppkf1[1] > 40)
			choice = 1;
		if (x2ikf2 < 50 && isITKaonOrPion(0, 2, de) && !isITKaonOrPion(1, 2, de) && ppkf2[0]>400 && ppkf2[1] > 40 && (choice == 0 || x2ikf2 < x2ikf1))
			choice = 2;
		bool kf10 = region[ipkf1[0] - 1] == 1 && amplitude[ipkf1[0] - 1] < 1000 && amplitudeE[ipkf1[0] - 1] < 1000 && amplitudeE[ipkf1[0] - 1] > 0.01;//&& amplitude[ipkf1[0] - 1] > 0.01;
		bool kf11 = region[ipkf1[1] - 1] == 1 && amplitude[ipkf1[1] - 1] < 1000 && amplitudeE[ipkf1[1] - 1] < 1000 && amplitudeE[ipkf1[1] - 1] > 0.01;//&& amplitude[ipkf1[1] - 1] > 0.01;
		bool kf20 = region[ipkf2[0] - 1] == 1 && amplitude[ipkf2[0] - 1] < 1000 && amplitudeE[ipkf2[0] - 1] < 1000 && amplitudeE[ipkf2[0] - 1] > 0.01;//&& amplitude[ipkf2[0] - 1] > 0.01;
		bool kf21 = region[ipkf2[1] - 1] == 1 && amplitude[ipkf2[1] - 1] < 1000 && amplitudeE[ipkf2[1] - 1] < 1000 && amplitudeE[ipkf2[1] - 1] > 0.01;//&& amplitude[ipkf2[1] - 1] > 0.01;
		h01->Fill(invmass2p( mce[0], mce[3] + mce[4] , mcpx[0] , mcpx[3]+mcpx[4] , mcpy[0], mcpy[3] + mcpy[4] , mcpz[0] , mcpz[3] + mcpz[4] ) );
		if ((goodKkpEntry1() && ((kf10 || ppkf1[0] < 300) && (kf11 || ppkf1[1] < 300))) && eton > 0.4 && eton < 0.85 && choice==1 && ks0InvMass(1)>440 && ks0InvMass(1)<530) {
			//if (fabs(mctheta[0] - thetakf1[1]) < 0.05 && fabs(mctheta[3] - thetakf1[0]) < 0.05) { //0,1
				//double p1x = ppkf1[1]*sin(thetakf1[1])*sin(phikf1[1]) + ppkf1[0]*sin(thetakf1[0])*sin(phikf1[0]);
				//double p1y = ppkf1[1]*sin(thetakf1[1])*cos(phikf1[1]) + ppkf1[0]*sin(thetakf1[0])*cos(phikf1[0]);
				//double p1z = ppkf1[1]*cos(thetakf1[1]) + ppkf1[0]*cos(thetakf1[0]);
				//double e1 = sqrt( pow(ppkf1[1],2) + pow(139.57039,2)) + sqrt( pow(ppkf1[0],2) + pow(493.677,2));

				double p2x = ppkf1[1]*sin(thetakf1[1])*sin(phikf1[1]) + ppkf1[2]*sin(thetakf1[2])*sin(phikf1[2]) + ppkf1[3]*sin(thetakf1[3])*sin(phikf1[3]) + ppkf1[4]*sin(thetakf1[4])*sin(phikf1[4]) + ppkf1[5]*sin(thetakf1[5])*sin(phikf1[5]);
				double p2y = ppkf1[1]*sin(thetakf1[1])*cos(phikf1[1]) + ppkf1[2]*sin(thetakf1[2])*cos(phikf1[2]) + ppkf1[3]*sin(thetakf1[3])*cos(phikf1[3]) + ppkf1[4]*sin(thetakf1[4])*cos(phikf1[4]) + ppkf1[5]*sin(thetakf1[5])*cos(phikf1[5]);
				double p2z = ppkf1[1]*cos(thetakf1[1]) + ppkf1[2]*cos(thetakf1[2]) + ppkf1[3]*cos(thetakf1[3]) + ppkf1[4]*cos(thetakf1[4]) + ppkf1[5]*cos(thetakf1[5]);
				double e2 = sqrt( pow(ppkf1[1],2) + pow(139.57039,2)) + ppkf1[2] + ppkf1[3] + ppkf1[4] + ppkf1[5];

				h2->Fill(amplitude[ipkf1[0]]);
				//h01->Fill(sqrt(pow(e1,2) - pow(p1x,2) - pow(p1y,2) - pow(p1z,2) ) );
				//h01->Fill(invmass2p( sqrt( pow(ppkf1[1],2) + pow(139.57039,2)), sqrt( pow(ppkf1[0],2) + pow(493.677,2)) , ppkf1[1] , ppkf1[0] , thetakf1[1], thetakf1[0] , phikf1[1] , phikf1[0] ) );
				//h0->Fill( sqrt(pow(e2,2) - pow(p2x,2) - pow(p2y,2) - pow(p2z,2) ) );
			//}
		}
		if ((goodKkpEntry2() && ((kf20 || ppkf2[0] < 300) && (kf21 || ppkf2[1] < 300))) && eton > 0.4 && eton < 0.85 && choice==2 && ks0InvMass(2)>440 && ks0InvMass(2)<530) {
			//if (fabs(mctheta[0] - thetakf2[1]) < 0.05 && fabs(mctheta[3] - thetakf2[0]) < 0.05) { //0,1
				double p1x = ppkf2[1]*sin(thetakf2[1])*sin(phikf2[1]) + ppkf2[0]*sin(thetakf2[0])*sin(phikf2[0]);
				double p1y = ppkf2[1]*sin(thetakf2[1])*cos(phikf2[1]) + ppkf2[0]*sin(thetakf2[0])*cos(phikf2[0]);
				double p1z = ppkf2[1]*cos(thetakf2[1]) + ppkf2[0]*cos(thetakf2[0]);
				double e1 = sqrt( pow(ppkf2[1],2) + pow(139.57039,2)) + sqrt( pow(ppkf2[0],2) + pow(493.677,2));

				double p2x = ppkf2[1]*sin(thetakf2[1])*sin(phikf2[1]) + ppkf2[2]*sin(thetakf2[2])*sin(phikf2[2]) + ppkf2[3]*sin(thetakf2[3])*sin(phikf2[3]) + ppkf2[4]*sin(thetakf2[4])*sin(phikf2[4]) + ppkf2[5]*sin(thetakf2[5])*sin(phikf2[5]);
				double p2y = ppkf2[1]*sin(thetakf2[1])*cos(phikf2[1]) + ppkf2[2]*sin(thetakf2[2])*cos(phikf2[2]) + ppkf2[3]*sin(thetakf2[3])*cos(phikf2[3]) + ppkf2[4]*sin(thetakf2[4])*cos(phikf2[4]) + ppkf2[5]*sin(thetakf2[5])*cos(phikf2[5]);
				double p2z = ppkf2[1]*cos(thetakf2[1]) + ppkf2[2]*cos(thetakf2[2]) + ppkf2[3]*cos(thetakf2[3]) + ppkf2[4]*cos(thetakf2[4]) + ppkf2[5]*cos(thetakf2[5]);
				double e2 = sqrt( pow(ppkf2[1],2) + pow(139.57039,2)) + ppkf2[2] + ppkf2[3] + ppkf2[4] + ppkf2[5];

				h2->Fill(amplitude[ipkf2[0]]);
				//h01->Fill(sqrt(pow(e1,2) - pow(p1x,2) - pow(p1y,2) - pow(p1z,2) ) );
				//h0->Fill( sqrt(pow(e2,2) - pow(p2x,2) - pow(p2y,2) - pow(p2z,2) ) );
			//}
		}
	}
	h0->SetLineColor(2);
	h01->Draw();
	//h0->Draw("same");
	//h2->Draw();
}

//drawing different tipes of dE vs P distributions (with compare of mc and kfi, for example).
void kkpDistr::drawdEvsp(vector<pair<double, double> > eranges) {
	cout << eranges.size() << endl;
	TChain chain("t1");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/kskp*");
	//chain.Add("/work/users/konctbel/calibs/R007-001/ksknpp_ks22p0_wrc_nemcc_odch-1000-10043-126-50000_wmix.root");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/exp2011/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/2etg/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp_ks0notfixed/*.root");
	
	chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2011IA/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2017In/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/4pi/2017/*x.root");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2011/kkp_KsNFn/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/exp2011_ks0notfixed/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/exp2011/*.root");
	string regime = "mod";
	{
		chain.SetBranchAddress("beam", &beam);
		chain.SetBranchAddress("eton", &eton);
		chain.SetBranchAddress("energy", &energy);
		chain.SetBranchAddress("nc", &nc);
		chain.SetBranchAddress("nn", &nn);
		chain.SetBranchAddress("z0", &z0);
		chain.SetBranchAddress("d0", &d0);
		chain.SetBranchAddress("x0", &x0);
		chain.SetBranchAddress("y0", &y0);
		chain.SetBranchAddress("phi", &phi);
		chain.SetBranchAddress("phis", &phis);
		chain.SetBranchAddress("theta", &theta);
		chain.SetBranchAddress("region", &region);
		chain.SetBranchAddress("amplitude", &amplitude);
		chain.SetBranchAddress("amplitudeE", &amplitudeE);
		chain.SetBranchAddress("dExs", &dExs);
		chain.SetBranchAddress("dExnC", &dExnC);
		chain.SetBranchAddress("dExn", &dExn);
		chain.SetBranchAddress("charge", &charge);

		chain.SetBranchAddress("x2ikf1", &x2ikf1);
		chain.SetBranchAddress("x2kf1", &x2kf1);
		chain.SetBranchAddress("ipkf1", &ipkf1);
		chain.SetBranchAddress("ppkf1", &ppkf1);
		chain.SetBranchAddress("phikf1", &phikf1);
		chain.SetBranchAddress("thetakf1", &thetakf1);

		chain.SetBranchAddress("x2ikf2", &x2ikf2);
		chain.SetBranchAddress("x2kf2", &x2kf2);
		chain.SetBranchAddress("ipkf2", &ipkf2);
		chain.SetBranchAddress("ppkf2", &ppkf2);
		chain.SetBranchAddress("phikf2", &phikf2);
		chain.SetBranchAddress("thetakf2", &thetakf2);

		chain.SetBranchAddress("x2ikf3", &x2ikf3);
		chain.SetBranchAddress("x2kf3", &x2kf3);

		if (regime == "exp") {
			chain.SetBranchAddress("dEx1C", &dEx1C);
			chain.SetBranchAddress("dEx2C", &dEx2C);
			chain.SetBranchAddress("dEx3C", &dEx3C);
			chain.SetBranchAddress("dEx4C", &dEx4C);
			chain.SetBranchAddress("dEx5C", &dEx5C);
			chain.SetBranchAddress("dEx6C", &dEx6C);
			chain.SetBranchAddress("dEx7C", &dEx7C);
			chain.SetBranchAddress("dEx8C", &dEx8C);
			chain.SetBranchAddress("dEx9C", &dEx9C);
		}

		if (regime == "mod") {
			chain.SetBranchAddress("mctheta", &mctheta);
			chain.SetBranchAddress("mcphi", &mcphi);
			chain.SetBranchAddress("mce", &mce);
			chain.SetBranchAddress("mcpx", &mcpx);
			chain.SetBranchAddress("mcpy", &mcpy);
			chain.SetBranchAddress("mcpz", &mcpz);
			chain.SetBranchAddress("mcpdg", &mcpdg);
		}
	}
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	int choice = 0;
	TH2* h0 = new TH2F("devspk","dEdx vs p,",800,0,800,300,0,10);
	TH2* h01 = new TH2F("devspk1","dEdx vs p,",800,0,800,300,0,10);
	TH2* h1 = new TH2F("ach1vspp","ach vs p,",250,300,800,20000,0,10);
	TProfile* p1 = new TProfile("pach1vspp","ach vs p,",250,300,800);
	TH2* h2 = new TH2F("ach0vspp","ach vs p,",250,300,800,20000,0,10);
	TProfile* p2 = new TProfile("pach0vspp","ach vs p,",250,300,800);
	TProfile* prof = new TProfile("effenergy", "eff vs energy",200, 0, 50);
	//
	ifstream ifilem("/work/users/kladov/snd2k/R007-001/2011/testMarksMod.txt");
	int linecount = 0;
	int nx2 = 0;
	int nkf = 0;
	int nge = 0;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//find gap of energy where we are
		findhetonInd(beam,eranges);
		//find excact beam position for exp
		findbeamInd(beam);
		//find the closest beam position for modeling
		findclosestPind(beam);
		if(indexCheck(regime)){
			normalizedExn(regime);
			normalizeAch();
			if(goodEventCheck(regime)){
				linecount+=1;
				int goodEv = 0;
				int whatKF = 0;
				double x2ikfML = 0;	
				ifilem >> goodEv >> whatKF >> x2ikfML;
				if (goodEv == 1) {
					int nntrue =0;
					for(size_t i = 0; i<nn+nc;i++){
						if(charge[i]==0 && energy[i] > 20)
							nntrue+=1;
					}
					if (whatKF == 1 && (goodKkpEntry1() && kf10 && kf11)/* && (ppkf1[0]>200 || tempdExn[0] > meanLine(ppkf1[0]))*/ /*&& (ppkf1[0]<300 || tempamp[0]< min(0.4,0.3/(double)amplitudeE[ipkf1[0] - 1])) && (ppkf1[1]<350 || tempamp[1]> 0.1)*/ && nntrue<7) {
					//if ((goodKkpEntry1() && (kf10 || ppkf1[0] < 300) && (kf11 || ppkf1[1] < 300)) && (ppkf1[0]>200 || dExn[ipkf1[0] - 1]/605.>meanLine(ppkf1[0])) && ppkf1[0]>50) {
						h0->Fill(ppkf1[0], tempdExn[0]);
						h01->Fill(ppkf1[1], tempdExn[1]);
						h2->Fill(ppkf1[0],tempamp[0]);
						p2->Fill(ppkf1[0],tempamp[0]);
						h1->Fill(ppkf1[1],tempamp[1]);
						p1->Fill(ppkf1[1],tempamp[1]);
						if (fabs(x2kf1 - x2ikfML) > 0.01)
							cout << x2ikfML << "	" << linecount << endl;
					}
					if (whatKF == 2 && (goodKkpEntry2() && kf20 && kf21) /*&& (ppkf2[0]>200 || tempdExn[2] > meanLine(ppkf2[0])) *//*&& (ppkf2[0]<300 || tempamp[2]< min(0.4,0.3/(double)amplitudeE[ipkf2[0] - 1])) && (ppkf2[1]<350 || tempamp[3]> 0.1)*/ && nntrue<7) {
					//if ((goodKkpEntry2() && (kf20 || ppkf2[0] < 300) && (kf21 || ppkf2[1] < 300)) && (ppkf2[0]>200 || dExn[ipkf2[0] - 1]/605.>meanLine(ppkf2[0])) && ppkf2[0]>50) {
						h0->Fill(ppkf2[0], tempdExn[2]);
						h01->Fill(ppkf2[1], tempdExn[3]);
						h2->Fill(ppkf2[0],tempamp[2]);
						p2->Fill(ppkf2[0],tempamp[2]);
						h1->Fill(ppkf2[1],tempamp[3]);
						p1->Fill(ppkf2[1],tempamp[3]);
						if (fabs(x2kf2 - x2ikfML) > 0.01)
							cout << x2ikfML << endl;
					}
				}
			}
			/*if ((goodKkpEntry1() && ((kf10 || ppkf1[0] < 300) && (kf11 || ppkf1[1] < 300))) && eton > 0.4 && eton < 0.85 && ks0InvMass(1)>440 && ks0InvMass(1)<530) {
			//if ((goodKkpEntry1() && ((kf10 && kf11))) && eton > 0.4 && eton < 0.85 && ks0InvMass(1)>440 && ks0InvMass(1)<530) {
				//if (fabs(mctheta[0] - thetakf1[1]) < 0.05 && fabs(mctheta[3] - thetakf1[0]) < 0.05) { //0,1
					//h1->Fill(x2ikf1,x2ikf3);
					prof->Fill(x2ikf1,x2ikf3);
					h01->Fill(ppkf1[1], dExn[ipkf1[1] - 1]/605.);
					h0->Fill(ppkf1[0], dExn[ipkf1[0] - 1]/605.);
					if(ppkf1[1]<300)
						amplitude[ipkf1[1]-1]=0;
					if(ppkf1[0]<300)
						amplitude[ipkf1[0]-1]=0;
					if(amplitudeE[ipkf1[0]-1]>0.001)
						h2->Fill(ppkf1[0],amplitude[ipkf1[0]-1]/amplitudeE[ipkf1[0]-1]);
					if(amplitudeE[ipkf1[1]-1]>0.001)
						h1->Fill(ppkf1[1],amplitude[ipkf1[1]-1]/amplitudeE[ipkf1[1]-1]);
				//}
			}
			if ((goodKkpEntry2() && ((kf20 || ppkf2[0] < 300) && (kf21 || ppkf2[1] < 300))) && eton > 0.4 && eton < 0.85 && ks0InvMass(2)>440 && ks0InvMass(2)<530) {
			//if ((goodKkpEntry2() && ((kf20 && kf21))) && eton > 0.4 && eton < 0.85 && ks0InvMass(2)>440 && ks0InvMass(2)<530) {
				//if (fabs(mctheta[0] - thetakf2[1]) < 0.05 && fabs(mctheta[3] - thetakf2[0]) < 0.05) { //0,1
					//h1->Fill(x2ikf2,x2ikf3);
					prof->Fill(x2ikf2,x2ikf3);
					h01->Fill(ppkf2[1], dExn[ipkf2[1] - 1])/605.;
					h0->Fill(ppkf2[0], dExn[ipkf2[0] - 1]/605.);
					if(ppkf2[1]<300)
						amplitude[ipkf2[1]-1]=0;
					if(ppkf2[0]<300)
						amplitude[ipkf2[0]-1]=0;
					if(amplitudeE[ipkf2[0]-1]>0.001)
						h2->Fill(ppkf2[0],amplitude[ipkf2[0]-1]/amplitudeE[ipkf2[0]-1]);
					if(amplitudeE[ipkf2[1]-1]>0.001)
						h1->Fill(ppkf2[1],amplitude[ipkf2[1]-1]/amplitudeE[ipkf2[1]-1]);
				//}
			}*/
		}
	}
	cout << "nx2 " << nx2 << ";	nkf " << nkf << ";	nge " << nkf << endl << endl;
	ifilem.close();
	//dedx
	/*h01->SetMarkerStyle(8);
	h01->SetMarkerSize(0.4);
	h01->SetMarkerColor(4);

	cout << linecount << endl;
	h0->SetMarkerStyle(8);
	h0->SetMarkerSize(0.4);
	h0->SetMarkerColor(2);
	h0->SetTitle("dE/dx vs p; P (MeV/c)");

	h0->Draw();
	h01->Draw("same");*/

	//ach
	h1->SetMarkerStyle(8);
	h1->SetMarkerSize(0.4);
	h1->SetMarkerColor(4);
	h2->SetMarkerStyle(8);
	h2->SetMarkerSize(0.4);
	h2->SetMarkerColor(2);
	h2->Draw();
	h1->Draw("same");

	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	c->SetBorderMode(0);
	c->SetFrameBorderMode(0);
	c->Update();

	/*TF1* f1 = new TF1("f1","meanLine(x)",100,200);
	f1->Draw("same");
	c->Update();*/

	/*TFile* MyFile1 = new TFile("/work/users/kladov/snd2k/R007-001/kkpi/dEvspT.root", "RECREATE");
	h0->Write("kaons");
	h1->Write("pions");
	MyFile1->Close();*/
}

//fit for dEdx histogramm distribution with p in certain range, writen long ago, have to change it later when decide to do this properly (after radcorrections)
void kkpDistr::landaufit() {
	TChain chain("t1");
	chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/kskp*");
	{
		chain.SetBranchAddress("eton", &eton);
		chain.SetBranchAddress("phi", &phi);
		chain.SetBranchAddress("phis", &phis);
		chain.SetBranchAddress("theta", &theta);
		chain.SetBranchAddress("region", &region);
		chain.SetBranchAddress("amplitude", &amplitude);
		chain.SetBranchAddress("amplitudeE", &amplitudeE);
		chain.SetBranchAddress("dExs", &dExs);
		chain.SetBranchAddress("x2ikf1", &x2ikf1);
		chain.SetBranchAddress("ipkf1", &ipkf1);
		chain.SetBranchAddress("ppkf1", &ppkf1);
		chain.SetBranchAddress("phikf1", &phikf1);
		chain.SetBranchAddress("thetakf1", &thetakf1);
		chain.SetBranchAddress("mctheta", &mctheta);
	}
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;

	TH1* h5 = new TH1F("de30", "dEdx, p 300-350", 50, 0, 30000);
	TH1* h4 = new TH1F("de25", "dEdx, p 250-300", 50, 0, 30000);
	TH1* h3 = new TH1F("de20", "dEdx, p 200-250", 50, 0, 30000);
	TH1* h2 = new TH1F("de15", "dEdx, p 150-200", 50, 0, 30000);
	TH1* h1 = new TH1F("de10", "dEdx, p 100-150", 50, 0, 30000);
	TH2* h = new TH2D("devspk", "dEdx vs p", 100, 0, 600, 100, 0, 30000);
	TProfile* pr = new TProfile("pr","pi_eff_vs_p",100,0,1000);
	TProfile* pr1 = new TProfile("pr1","pi_eff_vs_aE",100,0,10);
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//int fwt = 0; //fill what?
		if (mctheta[2] - thetakf1[0] < 0.07 && mctheta[2] - thetakf1[0] > -0.07 && x2ikf1 < 50) {
			h->Fill(ppkf1[0], dExs[ipkf1[0] - 1]);
			if (ppkf1[0] > 300 && ppkf1[0] < 350)
				h5->Fill(dExs[ipkf1[0] - 1]);
			if (ppkf1[0] > 250 && ppkf1[0] < 300)
				h4->Fill(dExs[ipkf1[0] - 1]);
			if (ppkf1[0] > 200 && ppkf1[0] < 250)
				h3->Fill(dExs[ipkf1[0] - 1]);
			if (ppkf1[0] > 150 && ppkf1[0] < 200)
				h2->Fill(dExs[ipkf1[0] - 1]);
			if (ppkf1[0] > 100 && ppkf1[0] < 150)
				h1->Fill(dExs[ipkf1[0] - 1]);
		}
		if (mctheta[0] - thetakf1[1] < 0.07 && mctheta[0] - thetakf1[1] > -0.07 && x2ikf1 < 50 && region[ipkf1[1]-1] == 1) {
			if (amplitude[ipkf1[1] - 1] > 0.2)
				pr->Fill(ppkf1[1], 1);
			else
				pr->Fill(ppkf1[1], 0);
			if (ppkf1[1] > 380) {
				if (amplitude[ipkf1[1] - 1] > 0.2)
					pr1->Fill(amplitudeE[ipkf1[1] - 1], 1);
				else
					pr1->Fill(amplitudeE[ipkf1[1] - 1], 0);
			}
		}
	}
	/*TF1* f1 = new TF1("f1", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0, 30000);
	h1->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	f1->SetParameter(1, kaonLine(125));
	f1->SetParameter(2, 0.5*kaonLine(125));
	f1->SetParLimits(2, 0.1*kaonLine(125), 1.0*kaonLine(125));
	h1->Fit(f1,"","",0,30000);
	c->Update();
	cin.get();
	h2->Draw();
	f1->SetParameter(1, kaonLine(175));
	f1->SetParameter(2, 0.5 * kaonLine(175));
	f1->SetParLimits(2, 0.1 * kaonLine(125), 1.0 * kaonLine(125));
	h2->Fit(f1, "", "", 0, 30000);
	c->Update();
	cin.get();
	h3->Draw();
	f1->SetParameter(1, kaonLine(225));
	f1->SetParameter(2, 0.5 * kaonLine(225));
	f1->SetParLimits(2, 0.1 * kaonLine(125), 1.0 * kaonLine(125));
	h3->Fit(f1, "", "", 0, 30000);
	c->Update();
	cin.get();
	h4->Draw();
	f1->SetParameter(1, kaonLine(275));
	f1->SetParameter(2, 0.5 * kaonLine(275));
	f1->SetParLimits(2, 0.05 * kaonLine(125), 1.5 * kaonLine(125));
	h4->Fit(f1, "", "", 0, 30000);
	c->Update();
	cin.get();
	h5->Draw();
	f1->SetParameter(1, kaonLine(325));
	f1->SetParameter(2, 0.5 * kaonLine(325));
	f1->SetParLimits(2, 0.05 * kaonLine(125), 1.5 * kaonLine(125));
	h5->Fit(f1, "", "", 0, 30000);
	c->Update();*/
	/*TF1* f1 = new TF1("f1", "[0] * (x*x + [1]*[1])/(x*x) * (log((x*x)/([1]*[1])) - (x*x)/(x*x + [1]*[1]) + [2])", 0, 30000);
	h->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	f1->SetParameter(0, 600);
	f1->SetParameter(1, 493);
	f1->SetParLimits(1, 480, 520);
	f1->SetParameter(2, 6.5);
	f1->SetParLimits(2, 6, 7);
	h->Fit(f1, "", "", 125, 600);*/
	/*pr->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->Update();
	TF1* f1 = new TF1("f1", "[0]+[1]/(1.+exp(-(x-[2])/[3]))", 0, 30000);
	f1->SetParameter(0, 0.03);
	f1->SetParameter(1, 0.9);
	f1->SetParLimits(1, 0.5, 1.);
	f1->SetParameter(2, 285);
	f1->SetParLimits(2, 250, 310);
	f1->SetParameter(3, 30);
	pr->Fit(f1, "", "", 140, 500);*/
	pr1->Draw();
}

double fitEffMcE(double *x, double *par) {
	return par[2]*(TMath::Erf(1./fabs(par[0]+(x[0]-par[1])*par[4])*(par[1]-x[0]))+par[3]);
}

//finction that can return different histograms from different sources (see parameters), now it returns invariant mass distribution in certain range of beam energies
//also it normalise histogramm in the end
//need to change energy cut to use predefined earray
vector<TH1*> kkpDistr::imvmass(char *add, char* name, char* title,int nbins,double lb, double hb, vector<pair<double,double> > eranges, string regime) {
	cout << eranges.size() << endl;
	TChain chain("t1");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/kskp*");
	chain.Add(add);
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/exp/*");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/exp2011/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-001/output/ntuples/2etg/*.root");
	{
		chain.SetBranchAddress("beam", &beam);
		chain.SetBranchAddress("eton", &eton);
		chain.SetBranchAddress("energy", &energy);
		chain.SetBranchAddress("nc", &nc);
		chain.SetBranchAddress("nn", &nn);
		chain.SetBranchAddress("z0", &z0);
		chain.SetBranchAddress("d0", &d0);
		chain.SetBranchAddress("x0", &x0);
		chain.SetBranchAddress("y0", &y0);
		chain.SetBranchAddress("phi", &phi);
		chain.SetBranchAddress("phis", &phis);
		chain.SetBranchAddress("theta", &theta);
		chain.SetBranchAddress("region", &region);
		chain.SetBranchAddress("amplitude", &amplitude);
		chain.SetBranchAddress("amplitudeE", &amplitudeE);
		chain.SetBranchAddress("dExs", &dExs);
		chain.SetBranchAddress("dExnC", &dExnC);
		chain.SetBranchAddress("dExn", &dExn);
		chain.SetBranchAddress("charge", &charge);

		chain.SetBranchAddress("x2ikf1", &x2ikf1);
		chain.SetBranchAddress("x2kf1", &x2kf1);
		chain.SetBranchAddress("ipkf1", &ipkf1);
		chain.SetBranchAddress("ppkf1", &ppkf1);
		chain.SetBranchAddress("phikf1", &phikf1);
		chain.SetBranchAddress("thetakf1", &thetakf1);

		chain.SetBranchAddress("x2ikf2", &x2ikf2);
		chain.SetBranchAddress("x2kf2", &x2kf2);
		chain.SetBranchAddress("ipkf2", &ipkf2);
		chain.SetBranchAddress("ppkf2", &ppkf2);
		chain.SetBranchAddress("phikf2", &phikf2);
		chain.SetBranchAddress("thetakf2", &thetakf2);

		chain.SetBranchAddress("x2ikf3", &x2ikf3);
		chain.SetBranchAddress("x2kf3", &x2kf3);
		chain.SetBranchAddress("ipkf3", &ipkf3);
		chain.SetBranchAddress("ppkf3", &ppkf3);
		chain.SetBranchAddress("phikf3", &phikf3);
		chain.SetBranchAddress("thetakf3", &thetakf3);

		if (regime == "exp") {
			chain.SetBranchAddress("dEx1C", &dEx1C);
			chain.SetBranchAddress("dEx2C", &dEx2C);
			chain.SetBranchAddress("dEx3C", &dEx3C);
			chain.SetBranchAddress("dEx4C", &dEx4C);
			chain.SetBranchAddress("dEx5C", &dEx5C);
			chain.SetBranchAddress("dEx6C", &dEx6C);
			chain.SetBranchAddress("dEx7C", &dEx7C);
			chain.SetBranchAddress("dEx8C", &dEx8C);
			chain.SetBranchAddress("dEx9C", &dEx9C);
		}

		if (regime == "mod") {
			chain.SetBranchAddress("mctheta", &mctheta);
			chain.SetBranchAddress("mcphi", &mcphi);
			chain.SetBranchAddress("mce", &mce);
			chain.SetBranchAddress("mcpx", &mcpx);
			chain.SetBranchAddress("mcpy", &mcpy);
			chain.SetBranchAddress("mcpz", &mcpz);
			chain.SetBranchAddress("mcpdg", &mcpdg);
		}
	}
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	TH1* h = new TH1F("invmass", "k+pi inv mass distr", 100, 500, 1500);
	TH1* h1 = new TH1F("cs", "k+pi inv mass distr", 30, 1400, 2000);
	TH1* h2 = new TH1F("invks0", "ks0 inv mass distr", 20, 497.611-65, 497.611+65);
	TH2F* hDeDx = new TH2F("hDeDx","dEdxn / 605",500,0,500,100,0,25);
	vector<TH1*> heton;
	if(regime == "mod")
		countModEvents.clear();
	if(regime == "exp")
		countExpEvents.clear();
	vector<TProfile*> effvsmce;
	double xEdges[26];
	for(size_t j = 0; j <26; j++)
		xEdges[j] = 0 + ((double)(j*j));

	//cout << (int)eranges.size() << endl;
	for (size_t i = 0; i < eranges.size(); i++) {
		heton.push_back(new TH1F(Form("invmassks%d",i), title, nbins, lb, hb));
		countModEvents.push_back(0);
		countExpEvents.push_back(0);
		
		//effvsmce.push_back(new TProfile(Form("effvsmce%d",i),"efficiency vs ISR photon energy",25,xEdges));
	}

	//TH1* heton = new TH1F("energyyy", "ks0 inv mass", 100, 0.0, 1.0);
	TH1* hx2kf = new TH1F("hx2kf", "x2 distribution", 200, 0, 200);
	TH1* hdedxloc = new TH1F("hdedxloc", "distribution;dE/dx", 200, 0, 10);
	TH1* hdedxlock = new TH1F("hdedxlock", "distribution;dE/dx", 200, 0, 10);
	TH1* hdedxlocp = new TH1F("hdedxlocp", "distribution;dE/dx", 200, 0, 10);
	TH1* hTrueEton = new TH1F("hTrueEton", "distribution;dE/dx", 100, 0, 2);

	TH1* hthetakp = new TH1F("theta", "theta between K and p0 in K0 system distr", 180, 0, 180);
	TH1* hthetakpmc = new TH1F("thetamc", "theta between K and p0 in K0 system distr mc", 180, 0, 180);

	TProfile* effvsp = new TProfile("effvsp","efficiency vs p",100,0,1000);

	TH1* hpsumr = new TH1F("hpsumr", "total pr distr", 100, 0, 200);
	TH1* hpsumz = new TH1F("hpsumz", "total pz distr", 200, -200, 200);
	TH1* hpsum = new TH1F("hpsum", "total p distr", 100, 0, 200);

	TH1* chi2 = new TH1F("chi2", "chi2", 50, 0, 100);

	TH1* hcru = new TH1F("hcru", "cru distribution", 100, 0, 10);
	int choice = 0;
	int linecount = 0;
	ifstream ifile("/work/users/kladov/snd2k/R007-001/2011/testMarksExp.txt");
	ifstream ifilem("/work/users/kladov/snd2k/R007-001/2011/testMarksMod.txt");
	vector<double> energyPoints;
	vector<double> countEventsBeam;
	vector<double> countSelectedEventsBeam;
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		double cru = 1;
		if (regime == "mod")
			cru = recalculateWeight();
		
		//changeOn4pi();

		//find where energy is
		findhetonInd(beam,eranges);
		findbeamInd(beam);
		if(regime == "exp" && beamInd == -1){
			cout << beam << endl;
		}
		findclosestPind(beam);
		if(indexCheck(regime)){
			
			//__________________________________//
			//for efficiency in mod by beam 
			bool newEnPoint = true;
			int pointIndex = energyPoints.size();
			for (size_t i = 0; i < energyPoints.size(); i++) {
				if (energyPoints[i] == (double)beam) {
					newEnPoint = false;
					pointIndex = i;
					//countEventsBeam[pointIndex]+=1;
				}
			}
			if (newEnPoint) {
				energyPoints.push_back((double)beam);
				countEventsBeam.push_back(1);
				countSelectedEventsBeam.push_back(0);
				effvsmce.push_back(new TProfile(Form("effvsmce%d",effvsmce.size()+1),"efficiency vs ISR photon energy",25,xEdges));
			}
			countModEvents[hetonInd] += 1;
			countExpEvents[hetonInd] += 1;
			//efficiency vs p
			int fillyn = 0;
			//__________________________________//
			countEventsBeam[pointIndex]+=cru;
			hcru->Fill(cru);
			normalizedExn(regime);
			normalizeAch();

			/*if(x2kf1<20){
				//countEventsBeam[pointIndex]+=1;
				//if(eGammas(1) < 0.6 * pow(500./beam,3)){
				//if( goodEventCheck(regime) && goodKkpEntry1() ){
					countSelectedEventsBeam[pointIndex]+=1;
				//}
			}
			else if(x2kf2<20){
				//countEventsBeam[pointIndex]+=1;
				//if(eGammas(2) < 0.6 * pow(500./beam,3)){
				//if( goodEventCheck(regime) && goodKkpEntry2() ){
					countSelectedEventsBeam[pointIndex]+=1;
				//}
			}*/

			if(goodEventCheck(regime)){
				linecount++;
				
				//if (cru > 0.01 && cru < 100)
				//countEventsBeam[pointIndex]+=1;
				
				int goodEv = 0;
				int whatKF = 0;
				double x2ikfML = 0;
				if (regime == "exp")
					ifile >> goodEv >> whatKF >> x2ikfML;
				if (regime == "mod")
					ifilem >> goodEv >> whatKF >> x2ikfML;

				//if (whatKF == 3)
				//	countEventsBeam[pointIndex]-=1;
				if (goodEv == 1) {
					
					//if (whatKF == 1 && (goodKkpEntry1() && kf10 && kf11)/* && (ppkf1[0]>200 || tempdExn[0] > meanLine(ppkf1[0]))*/ /*&& (ppkf1[0]<300 || tempamp[0]< min(0.4,0.3/(double)amplitudeE[ipkf1[0] - 1])) && (ppkf1[1]<350 || tempamp[1]> 0.1)*/)
					//	countSelectedEventsBeam[pointIndex]+=1;
					//if (whatKF == 2 && (goodKkpEntry2() && kf20 && kf21)/* && (ppkf2[0]>200 || tempdExn[2] > meanLine(ppkf2[0]))*/ /*&& (ppkf2[0]<300 || tempamp[2]< min(0.4,0.3/(double)amplitudeE[ipkf2[0]-1])) && (ppkf2[1]<350 || tempamp[1]> 0.1)*/)
					//	countSelectedEventsBeam[pointIndex]+=1;
					int nntrue =0;
					for(size_t i = 0; i<nn+nc;i++){
						if(charge[i]==0 && energy[i] > 20)
							nntrue+=1;
					}
					if (whatKF == 1 && (goodKkpEntry1() && kf10 && kf11) && (ppkf1[0]<300 || tempamp[0]< min(0.4,0.3/(double)amplitudeE[ipkf1[0] - 1])) /*&& (ppkf1[1]<350 || tempamp[1]> 0.1)*/) {
						hx2kf->Fill(x2kf1);
						if(ks0InvMass(1)>440 && ks0InvMass(1)<560)
							hTrueEton->Fill(eton);
						if (ppkf1[0]>250 && ppkf1[0]<300){
							hdedxloc->Fill(tempdExn[0]);
							hdedxlock->Fill(tempdExn[0]);
						}
						if (ppkf1[1]>250 && ppkf1[1]<300){
							hdedxloc->Fill(tempdExn[1]);
							hdedxlocp->Fill(tempdExn[1]);
						}

						if (regime == "mod")
							effvsmce[pointIndex]->Fill(mce[1],1);
						effvsp->Fill(ppkf1[0],1);
						fillyn = 1;
						//if (cru > 0.01 && cru < 100){
						countSelectedEventsBeam[pointIndex]+=cru;
						heton[hetonInd]->Fill(ks0InvMass(1),cru);
						//}
						if (fabs(x2kf1 - x2ikfML) > 0.01)
							cout << "something went wrong" << endl;
					}
					if (whatKF == 2 && (goodKkpEntry2() && kf20 && kf21) && (ppkf2[0]<300 || tempamp[2]< min(0.4,0.3/(double)amplitudeE[ipkf2[0]-1])) /*&& (ppkf2[1]<350 || tempamp[1]> 0.1)*/) {
						hx2kf->Fill(x2kf2);
						if(ks0InvMass(2)>440 && ks0InvMass(2)<560)
							hTrueEton->Fill(eton);
						if (ppkf2[0]>250 && ppkf2[0]<300){
							hdedxloc->Fill(tempdExn[2]);
							hdedxlock->Fill(tempdExn[2]);
						}
						if (ppkf2[1]>250 && ppkf2[1]<300){
							hdedxloc->Fill(tempdExn[3]);
							hdedxlocp->Fill(tempdExn[3]);
						}

						if (regime == "mod")
							effvsmce[pointIndex]->Fill(mce[1],1);
						effvsp->Fill(ppkf2[0],1);
						fillyn = 1;
						//if (cru > 0.01 && cru < 100){
						countSelectedEventsBeam[pointIndex]+=cru;
						heton[hetonInd]->Fill(ks0InvMass(2),cru);
						//}
						if (fabs(x2kf2 - x2ikfML) > 0.01)
							cout << "something went wrong" << endl;
					}
				}


				if (fillyn == 0 && (goodKkpEntry1() && kf10 && kf11)){
					effvsp->Fill(ppkf1[0],0);
				}
				else if (fillyn == 0 && (goodKkpEntry2() && kf20 && kf21)){
					effvsp->Fill(ppkf2[0],0);
				}
				
			}
			
			if (fillyn == 0 ){
				if (regime == "mod")
					effvsmce[pointIndex]->Fill(mce[1],0);
			}
		}


		if(indexCheck(regime)){
			normalizedExn(regime);
			normalizeAch();
			if(goodEventCheck(regime)){
				//cout << 2.*beam << "	" << eranges[0].second << endl;
				choice = 0;
				//_! and < for process, !_ and > for background
				int modechoice = 0;
				if (regime=="mod")
					modechoice = 1;
				if (x2ikf1 < 40 && isITKaonOrPion(0, 1, modechoice) && !isITKaonOrPion(1, 1, modechoice) /*&& ( ((x2ikf3 > x2ikf1-15 && x2ikf1<30) || x2ikf3>15))*/) 
					choice = 1;
				if (x2ikf2 < 40 && isITKaonOrPion(0, 2, modechoice) && !isITKaonOrPion(1, 2, modechoice) && (choice == 0 || x2ikf2 < x2ikf1) /*&& ( ((x2ikf3 > x2ikf2-15 && x2ikf2<30) || x2ikf3>15))*/) 
					choice = 2;
				/*if (choice != 0){
					if(choice == 1 && (goodKkpEntry1() && (kf10 || ppkf1[0] < 300) && (kf11 || ppkf1[1] < 300))){
						//cout << hetonInd << endl;
						heton[hetonInd]->Fill(ks0InvMass(1));
						//cout << hetonInd << endl;
						//if(hetonInd==0){
						//	cout << ks0InvMass(1) << endl;
						//}
						//cout << heton[hetonInd]->GetEntries() << endl;
					}
					else if(choice == 2 && goodKkpEntry2() && (kf20 || ppkf2[0] < 300) && (kf21 || ppkf2[1] < 300)){
						heton[hetonInd]->Fill(ks0InvMass(2));
					}
				}*/

				/*double pxk = ppkf1[0] * sin(thetakf1[0]) * cos(phikf1[0]);
					double pxp = ppkf1[1] * sin(thetakf1[1]) * cos(phikf1[1]);
					double pyk = ppkf1[0] * sin(thetakf1[0]) * sin(phikf1[0]);
					double pyp = ppkf1[1] * sin(thetakf1[1]) * sin(phikf1[1]);
					double pzk = ppkf1[0] * cos(thetakf1[0]);
					double pzp = ppkf1[1] * cos(thetakf1[1]);
					double ek = sqrt(ppkf1[0] * ppkf1[0] + 493.7 * 493.7);
					double ep = sqrt(ppkf1[1] * ppkf1[1] + 139.6 * 139.6);



					double invMasskp = sqrt(pow(ek + ep, 2) - pow(pxk + pxp, 2) - pow(pyk + pyp, 2) - pow(pzk + pzp, 2));
					//double invMassksp = sqrt(pow(eg1 + eg2 + eg3 + eg4 + ep, 2) - pow(pxg1 + pxg2 + pxg3 + pxg4 + pxp, 2) - pow(pyg1 + pyg2 + pyg3 + pyg4 + pyp, 2) - pow(pzg1 + pzg2 + pzg3 + pzg4 + pzp, 2));
					h1->Fill(2. * beam);
					h->Fill(invMasskp);*/

				//psum
				/*{
						double psumr = 0;
						double psumx = 0;
						double psumy = 0;
						double psumz = 0;
						double psum = 0;
						for (int i = 0; i < 6; i++) {
							psumx += ppkf1[i] * sin(thetakf1[i]) * cos(phikf1[i]);
							psumy += ppkf1[i] * sin(thetakf1[i]) * sin(phikf1[i]);
							psumz += ppkf1[i] * cos(thetakf1[i]);
						}
						psumr = sqrt(pow(psumx, 2) + pow(psumy, 2));
						hpsumr->Fill(psumr);
						hpsumz->Fill(psumz);
						hpsum->Fill(sqrt(pow(psumz, 2) + pow(psumr, 2)));
						if (sqrt(pow(psumz, 2) + pow(psumr, 2)) > 5)
							chi2->Fill(x2ikf1);
					}*/
					//angle kspi
				//{
					/*double cospk = (p0px1 * k0px + p0py1 * k0py + p0pz1 * k0pz) / (p0p1 * k0p);
						double sinpk = sqrt(1. - pow(cospk, 2));

						double cospk2 = (p0px2 * k0px + p0py2 * k0py + p0pz2 * k0pz) / (p0p2 * k0p);
						double sinpk2 = sqrt(1. - pow(cospk2, 2));

						double p0pcoll = p0p1 * cospk;
						double p0ppe = p0p1 * sinpk;

						double p0pcoll2 = p0p2 * cospk2;
						double p0ppe2 = p0p2 * sinpk2;

						double k0beta = 1. / sqrt(1. + 497.611 * 497.611 / (k0p * k0p));
						double k0gamma = k0p / 497.611 / k0beta;

						double p0pcolldash = k0gamma * (p0pcoll - k0beta * (eg1 + eg2));
						double p0pcolldash2 = k0gamma * (p0pcoll2 - k0beta * (eg3 + eg4));

						double thetadash = p0pcolldash / sqrt(pow(p0pcolldash, 2) + pow(p0ppe, 2));
						double thetadash2 = p0pcolldash2 / sqrt(pow(p0pcolldash2, 2) + pow(p0ppe2, 2));*/

					//heton->Fill(thetadash);
						//heton->Fill(thetadash2);
						// ks0 inv mass

					//if(thetadash<5 || thetadash2 <5)
						//	chi2->Fill(x2ikf1);
					
					// for modeling angle (kkpi only)
					/*double p0p1mc = sqrt(pow(mcpx[3], 2) + pow(mcpy[3], 2) + pow(mcpz[3], 2));
						double p0p1mc = sqrt(pow(mce[3], 2) - pow(139., 2));
						double p0p2mc = sqrt(pow(mcpx[4], 2) + pow(mcpy[4], 2) + pow(mcpz[4], 2));
						double p0p2mc = sqrt(pow(mce[4], 2) - pow(139., 2));
						double k0pmc = sqrt(pow(mcpx[1], 2) + pow(mcpy[1], 2) + pow(mcpz[1], 2));
						double k0pmc = sqrt(pow(mce[1], 2) - pow(497.6, 2));
						double cospkmc = (mcpx[1] * mcpx[3] + mcpy[1] * mcpy[3] + mcpz[1] * mcpz[3]) / (p0p1mc * k0pmc);
						double sinpkmc = sqrt(1. - pow(cospkmc, 2));
						double cospk2mc = (mcpx[1] * mcpx[4] + mcpy[1] * mcpy[4] + mcpz[1] * mcpz[4]) / (p0p2mc * k0pmc);
						double sinpk2mc = sqrt(1. - pow(cospk2mc, 2));
						double p0pcollmc = p0p1mc * cospkmc;
						double p0ppemc = p0p1mc * sinpkmc;
						double p0pcoll2mc = p0p2mc * cospk2mc;
						double p0ppe2mc = p0p2mc * sinpk2mc;
						double k0betamc = 1. / sqrt(1. + 497.611 * 497.611 / (k0pmc * k0pmc));
						double k0gammamc = k0pmc / 497.611 / k0betamc;
						double p0pcolldashmc = k0gammamc * (p0pcollmc - k0betamc * (mce[3]));
						double p0pcolldash2mc = k0gammamc * (p0pcoll2mc - k0betamc * (mce[4]));
						double thetadashmc = acos(p0pcolldashmc / sqrt(pow(p0pcolldashmc, 2) + pow(p0ppemc, 2))) * 180. / 3.14159;
						double thetadash2mc = acos(p0pcolldash2mc / sqrt(pow(p0pcolldash2mc, 2) + pow(p0ppe2mc, 2))) * 180. / 3.14159;
						hthetakpmc->Fill(thetadashmc);
						hthetakpmc->Fill(thetadash2mc);*/
					

				//}
			}
		}
	}


	cout << heton[0]->GetEntries() << endl;
	cout << "linecount	" << linecount << endl;
	ifile.close();
	ifilem.close();
	char* infoOname = "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelIDistrKsKPiE.root";
	if (regime == "mod") {
		infoOname = "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelIDistrKsKPiM.root";
	}
	TFile* MyFile = new TFile(infoOname, "RECREATE");
	//hDeDx->Write("dedxn0");
	hx2kf->Write("hx2kf");
	hdedxloc->Write("hdedxloc");
	hdedxlocp->Write("hdedxlocp");
	hdedxlock->Write("hdedxlock");
	hTrueEton->Write("hTrueEton");
	MyFile->Close();
	//hpsumz->Draw();
	/*hpsumr->SetLineColor(4);
	hpsumr->Draw("same");
	hpsum->SetLineColor(2);
	hpsum->Draw("same");*/
	//h2->Draw();
	//cout << heton->Integral(1,50) << endl;
	//chi2->DrawAt
	//hthetakp->Draw();
	//hthetakpmc->SetLineColor(2);
	//hthetakpmc->Draw("same");
	//if(erange != 3)
	//	heton->Scale(1590.8/countModEvents);
	for (size_t i = 0; i < eranges.size(); i++) {
		cout << "modeled events	" << eranges[i].first << "	" << countModEvents[i] << endl;
	}
	cout << heton[0]->GetEntries() << endl;
	if (regime == "mod"){
		for(size_t i = 0; i<energyPoints.size();i++){
			if(energyPoints[i]>0){
				energiesMod.push_back(energyPoints[i]);
				if(countEventsBeam[i]>0)
					efficienciesModBeam.push_back(countSelectedEventsBeam[i]/countEventsBeam[i]);
				else
					efficienciesModBeam.push_back(0);
				countModEventsBeam.push_back(countEventsBeam[i]);
			}
		}
	}
	cout << "here" << endl;
	//heton[0]->Draw();
	//TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	//c->Update();
	//cin.get();
	
	/*hcru->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	c->SetBorderMode(0);
	c->SetFrameBorderMode(0);
	c->SetGrid();
	c->Update();*/
	//cin.get();

	//effvsp->Draw();
	/*
	if (regime == "mod"){
		TFile* fout = new TFile("/work/users/kladov/snd2k/R007-001/2011/mceEff.root", "RECREATE");
		//for (size_t i = 0; i < eranges.size(); i++) {
		for (size_t i = 0; i < energiesMod.size(); i++) {
			//effvsmce[i]->Scale(1./effvsmce[i]->GetBinContent(1));
			effvsmce[i]->Draw();
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			c->SetFillColor(0);
			c->SetFrameFillColor(0);
			c->SetBorderMode(0);
			c->SetFrameBorderMode(0);
			c->SetGrid();
			c->Update();

			//TF1* f1 = new TF1("f1", "[2]*(TMath::Erf(1./[0]*([1]-x))+[3])", 0, 600);
			//TF1* f1 = new TF1("f1", "[2]*(TMath::Erf(1./fabs([0]+(x-[1])*[4])*([1]-x))+[3])", 0, 600);
			TF1* f1 = new TF1("f1", fitEffMcE, 0, 600,5);
			//TF1* f1 = new TF1("f1", "exp(-(x-[1])/[0])", 0, 600);
			f1->SetParameters(50, 150,0.5,1.1,0.5);
			f1->SetParLimits(0,10,100);
			f1->SetParLimits(1,50,200);
			f1->SetParLimits(2,0.01,1.5);
			f1->SetParLimits(3,0.5,1.5);
			f1->SetParLimits(4,0.05,1.0);
			//f1->FixParameter(2,0.5);
			f1->FixParameter(3,1.);
			effvsmce[i]->Fit("f1", "", "", 0, 500);
			c->Update();
			effvsmce[i]->Scale(0.5/f1->GetParameter(2));
			f1->FixParameter(2,0.5);
			effvsmce[i]->Fit("f1", "", "", 0, 500);
			vector<double> par(5);
			vector<double> parErr;
			f1->GetParameters(&par[0]);
			par[2] = 0.5;
			for(size_t j = 0; j<5; j++){
				parErr.push_back(f1->GetParError(j));
				pareffmce[j].push_back(par[j]);
				if(j == 2)
					pareffmceErr[j].push_back(parErr[j]/f1->GetParameter(j));
				else
					pareffmceErr[j].push_back(parErr[j]);
			}
			c->Update();
			//cin.get();
			effvsmce[i]->Write(Form("effmce%d",i));
		}
		fout->Close();
	}
	*/

	return heton;
}

void kkpDistr::genTableML(char *add, char *outputFile, int appYN) {
	TChain chain("t1");
	chain.Add(add);
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/4pi/2017/*x.root");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2019/*.root");
	{
		chain.SetBranchAddress("beam", &beam);
		chain.SetBranchAddress("eton", &eton);
		chain.SetBranchAddress("energy", &energy);
		chain.SetBranchAddress("nc", &nc);
		chain.SetBranchAddress("nn", &nn);
		chain.SetBranchAddress("z0", &z0);
		chain.SetBranchAddress("d0", &d0);
		chain.SetBranchAddress("x0", &x0);
		chain.SetBranchAddress("y0", &y0);
		chain.SetBranchAddress("phi", &phi);
		chain.SetBranchAddress("phis", &phis);
		chain.SetBranchAddress("theta", &theta);
		chain.SetBranchAddress("region", &region);
		chain.SetBranchAddress("amplitude", &amplitude);
		chain.SetBranchAddress("amplitudeE", &amplitudeE);
		chain.SetBranchAddress("dExs", &dExs);
		chain.SetBranchAddress("dExnC", &dExnC);
		chain.SetBranchAddress("dExn", &dExn);
		chain.SetBranchAddress("charge", &charge);

		chain.SetBranchAddress("x2ikf1", &x2ikf1);
		chain.SetBranchAddress("x2kf1", &x2kf1);
		chain.SetBranchAddress("ipkf1", &ipkf1);
		chain.SetBranchAddress("ppkf1", &ppkf1);
		chain.SetBranchAddress("phikf1", &phikf1);
		chain.SetBranchAddress("thetakf1", &thetakf1);

		chain.SetBranchAddress("x2ikf2", &x2ikf2);
		chain.SetBranchAddress("x2kf2", &x2kf2);
		chain.SetBranchAddress("ipkf2", &ipkf2);
		chain.SetBranchAddress("ppkf2", &ppkf2);
		chain.SetBranchAddress("phikf2", &phikf2);
		chain.SetBranchAddress("thetakf2", &thetakf2);

		chain.SetBranchAddress("x2ikf3", &x2ikf3);
		chain.SetBranchAddress("x2kf3", &x2kf3);
		chain.SetBranchAddress("ipkf3", &ipkf3);
		chain.SetBranchAddress("ppkf3", &ppkf3);
		chain.SetBranchAddress("phikf3", &phikf3);
		chain.SetBranchAddress("thetakf3", &thetakf3);

		chain.SetBranchAddress("mctheta", &mctheta);
		chain.SetBranchAddress("mcphi", &mcphi);
		chain.SetBranchAddress("mce", &mce);
		chain.SetBranchAddress("mcpx", &mcpx);
		chain.SetBranchAddress("mcpy", &mcpy);
		chain.SetBranchAddress("mcpz", &mcpz);
		chain.SetBranchAddress("mcpdg", &mcpdg);
	}
	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	int good = 0;
	int all = 0;
	ofstream ofile;
	if(appYN == 0)
		ofile.open(outputFile);
	else
		ofile.open(outputFile,std::ios_base::app);
	ofile << std::fixed;
	ofile << setprecision(4);
	TH2* h2 = new TH2F("ach0vspp","ach vs p,",250,300,800,2000,0,10);
	TH2* h1 = new TH2F("ach1vspp","ach vs p,",250,300,800,2000,0,10);
	
	TH1* hp0invmass = new TH1F("hp0invmass",";#pi^{0} invariant mass, MeV; N entries",100,50,250);
	TH1* hdedxloc = new TH1F("hdedxloc", "distribution;dE/dx", 200, 0, 10);
	TH1* hdedxlocp = new TH1F("hdedxlocp", "distribution;dE/dx", 200, 0, 10);
	TH1* hdedxlock = new TH1F("hdedxlock", "distribution;dE/dx", 200, 0, 10);

	TH1* hachloc = new TH1F("hachloc", "distribution;Amplitude, p.e.", 100, 0, 10);
	TH1* hachlocp = new TH1F("hachlocp", "distribution;Amplitude, p.e.", 100, 0, 10);
	TH1* hachlock = new TH1F("hachlock", "distribution;Amplitude, p.e.", 100, 0, 10);
	TProfile* pEffChi2 = new TProfile("pChi2","chi squared;P_{K}, MeV/c;Eff_{#chi^{2}<20}",70,100,800);

	TProfile* pAchPP = new TProfile("pAchPP",";P, MeV/c;amp / amp(ee)",100,0,1000);
	TProfile* pAchPK = new TProfile("pAchPK",";P, MeV/c;amp / amp(ee)",75,250,1000);

	TProfile* pdEdxP = new TProfile("pdEdxP",";P, MeV/c;dE/dx",100,0,1000);
	TProfile* pdEdxK = new TProfile("pdEdxK",";P, MeV/c;dE/dx",88,120,1000);
	srand (time(NULL));
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//changeOn4pi();
		bool p1 = !isITKaonOrPion(1, 1, 1);
		bool k1 = isITKaonOrPion(0, 1, 1);
		bool p2 = !isITKaonOrPion(1, 2, 1);
		bool k2 = isITKaonOrPion(0, 2, 1);

		findclosestPind(beam);
		if (mce[1] < 1000 && closestPind != -1){
			normalizedExn("mod");
			normalizeAch();
			if(goodEventCheck("mod")){
				int pionTag = 1; 
				if(appYN == 2)
					pionTag = 0;
				int kaonTag = appYN; 
				if(appYN == 2)
					kaonTag = 0;
				int tag = 0;
				int tag2 = 1;
				if (appYN == 1){
					tag = 2;
					tag2 = 2;
				}
				if (appYN == 2){
					tag = 3;
					tag2 = 3;
				}

				if (goodKkpEntry1() && kf10 && kf11) {
					if ((appYN == 1) || (appYN == 2) || (appYN == 0 && fabs(mctheta[0] - thetakf1[1]) < 0.03 && fabs(mctheta[3] - thetakf1[0]) < 0.03 && fabs(mcphi[0] - phikf1[1]) < 0.01 && fabs(mcphi[3] - phikf1[0]) < 0.01 && fabs(sqrt(pow(mcpx[0],2)+pow(mcpy[0],2)+pow(mcpz[0],2))-ppkf1[1]) < 20)) {
						//ofile << ppkf1[1] / 1000. << "	" << tempdExn[1] << "	" << tempamp[1] << "	" << pionTag << endl;
						//ofile << ppkf1[0] / 1000. << "	" << tempdExn[0] << "	" << tempamp[0] << "	" << kaonTag << endl;
						ofile << ppkf1[0] / 1000. << "	" << tempdExn[0] << "	" << tempamp[0] <<  "	" << tempampE[0] << "	" << ppkf1[1] / 1000. << "	" << tempdExn[1] << "	" << tempamp[1] <<  "	" << tempampE[1] <<  "	" << x2kf1 << "	" << x2kf3 << "	" << beam/750. << "	" << tag << endl;
						//if (x2kf2<1000)
						//	ofile << ppkf2[0] / 1000. << "	" << tempdExn[2] << "	" << tempamp[2] << "	" << ppkf2[1] / 1000. << "	" << tempdExn[3] << "	" << tempamp[3] <<  "	" << x2kf2 << "	" << tag2 << endl;
						
						//if (tempamp[0]<0.5)
						pAchPK->Fill(ppkf1[0],tempamp[0]);
						pAchPP->Fill(ppkf1[1],tempamp[1]);
						hp0invmass->Fill(gammagammaInvMasskf(1));

						pdEdxK->Fill(ppkf1[0],tempdExn[0]);
						pdEdxP->Fill(ppkf1[1],tempdExn[1]);


						if (ppkf1[0]>300 && ppkf1[0]<350){
							hdedxloc->Fill(tempdExn[0]);
							hdedxlock->Fill(tempdExn[0]);
							hachloc->Fill(tempamp[0]);
							hachlock->Fill(tempamp[0]);
						}
						if (ppkf1[1]>300 && ppkf1[1]<350){
							hdedxloc->Fill(tempdExn[1]);
							hdedxlocp->Fill(tempdExn[1]);
							hachloc->Fill(tempamp[1]);
							hachlocp->Fill(tempamp[1]);
						}
						
						int iSecret = rand() % 100;
						if(iSecret >= 80){
							h2->Fill(ppkf1[0],tempamp[0]);
							h1->Fill(ppkf1[1],tempamp[1]);
						}
						all += 1;
						if (!k1)
							good += 1;
						//if this particle which is good pion by mc angles comparison is ok by my separation criteria - +1 to efficiency
						//but I am more interested in possibility that pion will be marked as kaon
						//all += 1;
						//if (!p1)
						//	good += 1;
						if(x2kf1<x2kf3)
							pEffChi2->Fill(ppkf1[0],1);
						else
							pEffChi2->Fill(ppkf1[0],0);
					}
					if ((appYN == 0 && fabs(mctheta[0] - thetakf1[0]) < 0.05 && fabs(mctheta[3] - thetakf1[1]) < 0.05 && fabs(mcphi[0] - phikf1[0]) < 0.01 && fabs(mcphi[3] - phikf1[1]) < 0.01)) {
						ofile << ppkf1[0] / 1000. << "	" << tempdExn[0] << "	" << tempamp[0] <<  "	" << tempampE[0] << "	" << ppkf1[1] / 1000. << "	" << tempdExn[1] << "	" << tempamp[1] <<  "	" << tempampE[1] << "	" << x2kf1 << "	" << x2kf3 << "	" << beam/750. << "	" << tag2 << endl;
					}
				}
				if (goodKkpEntry2() && kf20 && kf21) {
					if ((appYN == 1) || (appYN == 2) || (appYN == 0 && fabs(mctheta[0] - thetakf2[1]) < 0.03 && fabs(mctheta[3] - thetakf2[0]) < 0.03 && fabs(mcphi[0] - phikf2[1]) < 0.01 && fabs(mcphi[3] - phikf2[0]) < 0.01 && fabs(sqrt(pow(mcpx[0],2)+pow(mcpy[0],2)+pow(mcpz[0],2))-ppkf2[1]) < 20)) {
						//ofile << ppkf2[1] / 1000. << "	" << tempdExn[3] << "	" << tempamp[3] << "	" << pionTag << endl;
						//ofile << ppkf2[0] / 1000. << "	" << tempdExn[2] << "	" << tempamp[2] << "	" << kaonTag << endl;
						ofile << ppkf2[0] / 1000. << "	" << tempdExn[2] << "	" << tempamp[2] <<  "	" << tempampE[2] << "	" << ppkf2[1] / 1000. << "	" << tempdExn[3] << "	" << tempamp[3]<<  "	" << tempampE[3] << "	" << x2kf2 << "	" << x2kf3 << "	" << beam/750. << "	" << tag << endl;
						//if (x2kf1<1000)
						//	ofile << ppkf1[0] / 1000. << "	" << tempdExn[0] << "	" << tempamp[0] << "	" << ppkf1[1] / 1000. << "	" << tempdExn[1] << "	" << tempamp[1] << "	" << x2kf1 << "	" << tag2 << endl;
						
						//if (tempamp[2]<0.5)
						pAchPK->Fill(ppkf2[0],tempamp[2]);
						pAchPP->Fill(ppkf2[1],tempamp[3]);
						hp0invmass->Fill(gammagammaInvMasskf(2));

						pdEdxK->Fill(ppkf2[0],tempdExn[2]);
						pdEdxP->Fill(ppkf2[1],tempdExn[3]);

						
						if (ppkf2[0]>300 && ppkf2[0]<350){
							hdedxloc->Fill(tempdExn[2]);
							hdedxlock->Fill(tempdExn[2]);
							hachloc->Fill(tempamp[2]);
							hachlock->Fill(tempamp[2]);
						}
						if (ppkf2[1]>300 && ppkf2[1]<350){
							hdedxloc->Fill(tempdExn[3]);
							hdedxlocp->Fill(tempdExn[3]);
							hachloc->Fill(tempamp[3]);
							hachlocp->Fill(tempamp[3]);
						}

						int iSecret = rand() % 100;
						if(iSecret >= 80){
							h2->Fill(ppkf2[0],tempamp[2]);
							h1->Fill(ppkf2[1],tempamp[3]);
						}
						//if this particle which is good kaon by mc angles comparison is ok by my separation criteria - +1 to efficiency
						all += 1;
						if (!k2)
							good += 1;
						//if this particle which is good pion by mc angles comparison is ok by my separation criteria - +1 to efficiency
						//all += 1;
						//if (!p2)
						//	good += 1;
						if(x2kf2<x2kf3)
							pEffChi2->Fill(ppkf2[0],1);
						else
							pEffChi2->Fill(ppkf2[0],0);
					}
					if ((appYN == 0 && fabs(mctheta[0] - thetakf2[0]) < 0.05 && fabs(mctheta[3] - thetakf2[1]) < 0.05 && fabs(mcphi[0] - phikf2[0]) < 0.01 && fabs(mcphi[3] - phikf2[1]) < 0.01)) {
						ofile << ppkf2[0] / 1000. << "	" << tempdExn[2] << "	" << tempamp[2] <<  "	" << tempampE[2] << "	" << ppkf2[1] / 1000. << "	" << tempdExn[3] << "	" << tempamp[3] <<  "	" << tempampE[3] <<  "	" << x2kf2 << "	"  << x2kf3 << "	" << beam/750. << "	" << tag2 << endl;
					}
				}
			}
		}
	}
	ofile.close();
	cout << (double)good / (double)all << endl;

	if (appYN == 0){
		TFile* MyFile = new TFile("/work/users/kladov/snd2k/R007-001/kkpi/2011/trainDistrKsKPi.root", "RECREATE");
		hdedxloc->Write("hdedxloc");
		hdedxlocp->Write("hdedxlocp");
		hdedxlock->Write("hdedxlock");

		hachloc->Write("hachloc");
		hachlocp->Write("hachlocp");
		hachlock->Write("hachlock");

		pEffChi2->Write("pEffChi2");

		pAchPP->Write("pAchPP");
		pAchPK->Write("pAchPK");

		pdEdxP->Write("pdEdxP");
		pdEdxK->Write("pdEdxK");

		hp0invmass->Write("hp0invmass");

		MyFile->Close();
	}

	/*h1->SetMarkerStyle(8);
	h1->SetMarkerSize(0.4);
	h1->SetMarkerColor(4);
	h1->SetLineColor(4);
	h2->SetMarkerStyle(8);
	h2->SetMarkerSize(0.4);
	h2->SetMarkerColor(2);
	h2->SetLineColor(2);
	h2->Draw();
	h1->Draw("same");

	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	c->SetFillColor(0);
	c->SetFrameFillColor(0);
	c->SetBorderMode(0);
	c->SetFrameBorderMode(0);
	c->Update();*/
}

void kkpDistr::genTableMLExp(char *add, char *outdistr, char *txtfile, int ind, vector<pair<double, double> > eranges, string regime) {
	TChain chain("t1");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2017/kkp_KsNFT/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2017/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2019/kkp_KsNF/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2019/*.root");
	//chain.Add("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2017/*.root");
	chain.Add(add);
	{
		chain.SetBranchAddress("beam", &beam);
		chain.SetBranchAddress("eton", &eton);
		chain.SetBranchAddress("energy", &energy);
		chain.SetBranchAddress("nc", &nc);
		chain.SetBranchAddress("nn", &nn);
		chain.SetBranchAddress("z0", &z0);
		chain.SetBranchAddress("d0", &d0);
		chain.SetBranchAddress("x0", &x0);
		chain.SetBranchAddress("y0", &y0);
		chain.SetBranchAddress("phi", &phi);
		chain.SetBranchAddress("phis", &phis);
		chain.SetBranchAddress("theta", &theta);
		chain.SetBranchAddress("region", &region);
		chain.SetBranchAddress("amplitude", &amplitude);
		chain.SetBranchAddress("amplitudeE", &amplitudeE);
		chain.SetBranchAddress("dExs", &dExs);
		chain.SetBranchAddress("dExnC", &dExnC);
		chain.SetBranchAddress("dExn", &dExn);
		chain.SetBranchAddress("charge", &charge);

		chain.SetBranchAddress("x2ikf1", &x2ikf1);
		chain.SetBranchAddress("x2kf1", &x2kf1);
		chain.SetBranchAddress("ipkf1", &ipkf1);
		chain.SetBranchAddress("ppkf1", &ppkf1);
		chain.SetBranchAddress("phikf1", &phikf1);
		chain.SetBranchAddress("thetakf1", &thetakf1);

		chain.SetBranchAddress("x2ikf2", &x2ikf2);
		chain.SetBranchAddress("x2kf2", &x2kf2);
		chain.SetBranchAddress("ipkf2", &ipkf2);
		chain.SetBranchAddress("ppkf2", &ppkf2);
		chain.SetBranchAddress("phikf2", &phikf2);
		chain.SetBranchAddress("thetakf2", &thetakf2);

		chain.SetBranchAddress("x2ikf3", &x2ikf3);
		chain.SetBranchAddress("x2kf3", &x2kf3);
		chain.SetBranchAddress("ipkf3", &ipkf3);
		chain.SetBranchAddress("ppkf3", &ppkf3);
		chain.SetBranchAddress("phikf3", &phikf3);
		chain.SetBranchAddress("thetakf3", &thetakf3);

		if (regime == "exp") {
			chain.SetBranchAddress("dEx1C", &dEx1C);
			chain.SetBranchAddress("dEx2C", &dEx2C);
			chain.SetBranchAddress("dEx3C", &dEx3C);
			chain.SetBranchAddress("dEx4C", &dEx4C);
			chain.SetBranchAddress("dEx5C", &dEx5C);
			chain.SetBranchAddress("dEx6C", &dEx6C);
			chain.SetBranchAddress("dEx7C", &dEx7C);
			chain.SetBranchAddress("dEx8C", &dEx8C);
			chain.SetBranchAddress("dEx9C", &dEx9C);
		}

		if (regime == "mod") {
			chain.SetBranchAddress("mctheta", &mctheta);
			chain.SetBranchAddress("mcphi", &mcphi);
			chain.SetBranchAddress("mce", &mce);
			chain.SetBranchAddress("mcpx", &mcpx);
			chain.SetBranchAddress("mcpy", &mcpy);
			chain.SetBranchAddress("mcpz", &mcpz);
			chain.SetBranchAddress("mcpdg", &mcpdg);
		}
	}

	TH1* hInvMass = new TH1F("invks01", "ks0 inv mass distr", 50, 497.611 - 265, 497.611 + 265);
	TH1* hChi2 = new TH1F("Chi2", "chi squared", 200, -100, 100);
	TProfile* pEffChi2 = new TProfile("pEffChi2","chi squared;P, MeV/c;#chi^{2}",70,100,800);
	TH1* hEgammas = new TH1F("hegammas", "energy dep of not used in kinfit particles", 100, 0, 1.5);
	TH1* hEKPi = new TH1F("heKPi", "energy dep of k and pi charged", 100, 0, 1.5);
	TH1* hEton = new TH1F("heton", "eton", 100, 0, 1.5);
	TH1* hetonT = new TH1F("hetonT", "eton", 100, 0, 1.5);
	TH1* hd0 = new TH1F("hd0", "d0", 250, -2.5, 2.5);
	TH1* hz0 = new TH1F("hz0", "z0", 250, -25, 25);
	TH1* hdz0 = new TH1F("hdz0", "z01 - z00", 250, -25, 25);
	TH1* htheta = new TH1F("htheta", "theta", 250, 0, 180);
	TH1* hthetaKsP0 = new TH1F("hthetaKsP0", "thetaKsP0", 50, -1, 1);
	TH1* hp0invmass = new TH1F("hp0invmass", "p0invmass", 50, 0, 500);
	TH2F* hDeDx = new TH2F("hDeDx","dEdxn / 605",500,0,500,100,0,25);
	TH1* hchoice = new TH1D("hchoice","choice",3,0,3);
	TProfile* pampphi = new TProfile("pampphi","amplitude profile",20,2,3);
	TProfile* pampphi2 = new TProfile("pampphi2","amplitude profile",30,2,3);
	TProfile* pampphi3 = new TProfile("pampphi3","amplitude profile",30,2,3);
	TProfile* pampphi4 = new TProfile("pampphi4","amplitude profile",30,2,3);
	TProfile* pampEphi = new TProfile("pampEphi","amplitudeE profile",30,2,3);
	TProfile* effvspa = new TProfile("effvspa","efficiency vs p",100,0,1000);
	TFile* MyFile = new TFile(outdistr, "RECREATE");

	

	vector<int> countSelEventsBeam;
	vector<int> countAllEventsBeam;
	for(size_t i = 0; i<energiesEXP.size(); i++){
		countAllEventsBeam.push_back(0);
		countSelEventsBeam.push_back(0);
	}

	vector<TH1*> hegmax;
	vector<TH1*> hegmean;
	for (size_t i = 0; i < eranges.size(); i++){
		hegmax.push_back(new TH1F(Form("hegmax%d",i), "energy dep of maximum energy photon except 4 from kinfit", 100, 0, 1.5));
		hegmean.push_back(new TH1F(Form("hegmean%d",i), "energy dep each photon from kinfit", 50, 0, 400));
	}

	const Long64_t entries = chain.GetEntries();
	cout << entries << endl;
	int good = 0;
	int all = 0;
	int all1 = 0;
	int selected = 0;
	int good1 = 0;
	int good2 = 0;
	int bad = 0;
	int choice = 0;
	ofstream ofile(txtfile);
	//ofstream ofile("/work/users/kladov/snd2k/R007-001/2017/testSetMod.txt",std::ios_base::app);
	ofile << std::fixed;
	ofile << setprecision(4);
	int linecount = 0;
	int nEntries[11] = {0,0,0,0,0,0,0,0,0,0,0};
	double nol = 0;
	
	for (int e = 0; e < entries; e++) {
		chain.GetEntry(e);
		//changeOn4pi();
		//hthetaKsP0->Fill((anglekspi(1,1)));
		//hthetaKsP0->Fill((anglekspi(1,2)));
		if (x2kf1<20){
			hthetaKsP0->Fill(anglekspi(1,1));
			hthetaKsP0->Fill(anglekspi(1,2));
		}
		if (x2kf2<20){
			hthetaKsP0->Fill(anglekspi(2,1));
			hthetaKsP0->Fill(anglekspi(2,2));
		}
		//find where energy is
		findhetonInd(beam,eranges);
		findbeamInd(beam);
		findclosestPind(beam);
		if(beamInd <0 && regime == "exp")
			cout << beam << endl;
		if(indexCheck(regime)){
			all1 += 1;
			normalizedExn(regime);
			normalizeAch();
			if(goodEventCheck(regime)){
				all += 1;
				linecount++;
				selected += 1;
				if (goodKkpEntry1() && kf10 && kf11) {
					ofile << ppkf1[0] / 1000. << "	" << tempdExn[0] << "	" << tempamp[0] <<  "	" << tempampE[0] << "	";
					ofile << ppkf1[1] / 1000. << "	" << tempdExn[1] << "	" << tempamp[1] <<  "	" << tempampE[1] << "	" << x2kf1 << "	";
				}
				else
					ofile << "0.2000	3.5000	0.0000	1.0000	0.5000	0.5000	0.2000	1.0000	100.0000	";
				ofile << x2kf3 << "	" << beam/750. << "	" << ind << endl;
				if (goodKkpEntry2() && kf20 && kf21) {
					ofile << ppkf2[0] / 1000. << "	" << tempdExn[2] << "	" << tempamp[2] <<  "	" << tempampE[2] << "	";
					ofile << ppkf2[1] / 1000. << "	" << tempdExn[3] << "	" << tempamp[3] <<  "	" << tempampE[3] << "	" << x2kf2 << "	";
				}
				else
					ofile << "0.2000	3.5000	0.0000	1.0000	0.5000	0.5000	0.2000	1.0000	100.0000	";
				ofile << x2kf3 << "	" << beam/750. << "	" << ind << endl;
			

				//Fill  hp0inv mass with all possible pairs of excessive photons (excluding kinfit)
				hp0invmass->Fill(extraPartInvMass());
				
			}
		}
		
		if(indexCheck(regime)){
			//normalizedExn(regime);
			//normalizeAch();
			if(goodEventCheck(regime)){
				if (goodKkpEntry1() && kf10 && kf11) {
					if(x2kf1<x2kf3)
						pEffChi2->Fill(ppkf1[0],1);
					else
						pEffChi2->Fill(ppkf1[0],0);
					//hthetaKsP0->Fill(anglekspi(1));
					//random shit
					if (choice == 1 && ks0InvMass(1) > 400 && ks0InvMass(1) < 600){
						good1 += 1;
						int indET = 2 + rand() % 4; //2-3-4-5
						//cout << ppkf1[indET] << endl;
						hegmean[0]->Fill(ppkf1[indET]);
						//find first particle by energy that is not used in kinfit and its charge == 0, then write hegmax with it
						if(x2ikf1<50 && hetonInd < hegmax.size() && nn>=4 && nc>=2){
							int egmeind = 0;
							int it = 0;
							bool excessive = false;
							while(it < nc+nn && egmeind == 0){
								int equality = 0;
								for(size_t j = 0; j < 6; j++)
									if(ipkf1[j]-1 == it)
										equality = 1;
								if (equality == 0 && charge[it] == 0){
									egmeind = it;
									excessive = true;
								}
								it+=1;
							}
							if (excessive)
								hegmax[hetonInd]->Fill(energy[egmeind]/(2.*beam));
							else
								hegmax[hetonInd]->Fill(0.);
						}
					}
				}
				if (goodKkpEntry2() && kf20 && kf21) {
					if(x2kf2<x2kf3)
						pEffChi2->Fill(ppkf2[0],1);
					else
						pEffChi2->Fill(ppkf2[0],0);
					//hthetaKsP0->Fill(anglekspi(2));
					//random shit
					if (choice == 2 && ks0InvMass(2) > 400 && ks0InvMass(2) < 600){
						good2 += 1;
						int indET = 2 + rand() % 4; //2-3-4-5
						hegmean[0]->Fill(ppkf2[indET]);
					}
				}
			}
			
			

			//random shit

			//filling vector of temp efficiencies for test
			//if(regime == "mod")
			//	countAllEventsBeam[closestPind]+=1;
			int fillyn = -1;
			if(regime == "exp")
				countAllEventsBeam[beamInd]+=1;
			if(regime == "mod"){
				//if (beam>900)
				//	cout << closestPind << "	" << beam << endl;
				//countSelEventsBeam[closestPind]+=1;
				countAllEventsBeam[closestPind]+=1;
				if(goodEventCheck(regime)){
					countSelEventsBeam[closestPind]+=1;
					fillyn = 0;
					if(kf10 && kf11){
						effvspa->Fill(beam,1);
						fillyn = 1;
					}
				}
			}
			if(fillyn == 0)
				effvspa->Fill(beam,0);
			if(regime == "exp")
				countSelEventsBeam[beamInd]+=1;

			choice = 0;
			//_! and < for process, !_ and > for background
			if (x2ikf1 < 50 && isITKaonOrPion(0, 1, 1) && !isITKaonOrPion(1, 1, 1))
				choice = 1;
			if (x2ikf2 < 50 && isITKaonOrPion(0, 2, 1) && !isITKaonOrPion(1, 2, 1) && (choice == 0 || x2ikf2 < x2ikf1))
				choice = 2;


			if(x2ikf1<50)
				hchoice->Fill(choice);
			if (choice == 1 && x2ikf1<50){
				//if (x2ikf1<50 && ppkf1[1]>400 && region[ipkf1[1]-1]==1 && amplitude[ipkf1[1]-1]>0.4 && amplitude[ipkf1[1]-1]<100){
				//	pampphi->Fill(phikf1[1],amplitude[ipkf1[1]-1]);
				//	pampEphi->Fill(phikf1[1],amplitudeE[ipkf1[1]-1]);
				//}
				hEton->Fill(eton);
				//hEgammas->Fill(eGammas(1));
				hEKPi->Fill(eKPi(1));
				htheta->Fill(theta[ipkf1[0] - 1] * 180. / 3.14159);
				htheta->Fill(theta[ipkf1[1] - 1] * 180. / 3.14159);

				hInvMass->Fill(ks0InvMass(1));
			}

			if (x2ikf1 < 50 && (goodKkpEntry1() && ((kf10 || ppkf1[0] < 300) && (kf11 || ppkf1[1] < 300))) && ks0InvMass(1) > 400 && ks0InvMass(1) < 600 && eton > 0.4 && eton < 0.85 && choice==1) {
				if (x2ikf1<50 && ppkf1[1]>400 && region[ipkf1[1]-1]==1 && amplitude[ipkf1[1]-1]>0.4 && amplitude[ipkf1[1]-1]<100){
					pampphi->Fill(phikf1[1],amplitude[ipkf1[1]-1]);
					pampEphi->Fill(phikf1[1],amplitudeE[ipkf1[1]-1]);
				}
			}
			if (choice==1) {
				if (x2ikf1<50 && ppkf1[1]>400 && region[ipkf1[1]-1]==1 && amplitude[ipkf1[1]-1]>0.4 && amplitude[ipkf1[1]-1]<100){
					pampphi2->Fill(phikf1[1],amplitude[ipkf1[1]-1]);
				}
			}
			if (x2ikf1 < 50 && (goodKkpEntry1() && ((kf10 || ppkf1[0] < 300) && (kf11 || ppkf1[1] < 300))) && ks0InvMass(1) > 400 && ks0InvMass(1) < 600 && eton > 0.4 && eton < 0.85) {
				if (x2ikf1<50 && ppkf1[1]>400 && region[ipkf1[1]-1]==1 && amplitude[ipkf1[1]-1]>0.4 && amplitude[ipkf1[1]-1]<100){
					pampphi3->Fill(phikf1[1],amplitude[ipkf1[1]-1]);
				}
			}
			if (x2ikf1<50 && ppkf1[1]>400 && region[ipkf1[1]-1]==1 && amplitude[ipkf1[1]-1]>0.4 && amplitude[ipkf1[1]-1]<100){
				pampphi4->Fill(phikf1[1],amplitude[ipkf1[1]-1]);
			}


			//x2ikf1 - have to cut seriously(1/3 mod, 1/16 exp), choice (neseccary) - also (-10% mod,), eton (strong cut), 
			//thetabetween (?-seems like it cuts nothing, so why I would cut it? if I will completely remove this - nothing will change I suppose, maybe systematic will rise (but how)),
			//centrals, union vertex, egammas (it is generally coincide with isr photon I guess, but can cut 2etagamma (5 photons)), eKpi (try without it),
			//nc<, nn<, ppkf1>, kf(region) (can try without it)
			if(x2kf1 < 20 || x2kf2 < 20)
				nEntries[1]+=1;
			if(x2kf1 < 20 && ppkf1[0] > 10 && ppkf1[1] > 5 && nn <=7 && nc <=3 || x2kf2<20 && ppkf2[0] > 10 && ppkf2[1] > 5 && nn <=7 && nc <=3)
				nEntries[2]+=1;
			if(x2kf1 < 20 && ppkf1[0] > 10 && ppkf1[1] > 5 && nn <=7 && nc <=3 && eGammas(1) < 0.2 && eKPi(1) > 0.1 || x2kf2<20 && ppkf2[0] > 10 && ppkf2[1] > 5 && nn <=7 && nc <=3 && eGammas(2) < 0.2 && eKPi(2) > 0.1)
				nEntries[3]+=1;
			if(goodKkpEntry1() || goodKkpEntry2())
				nEntries[4]+=1;
			if (((goodKkpEntry1() && kf10 && kf11) || (goodKkpEntry2() && kf20 && kf21)))
				nEntries[5]+=1;
			if(((goodKkpEntry1() && kf10 && kf11) || (goodKkpEntry2() && kf20 && kf21)) && eton > 0.4 && eton < 0.85)
				nEntries[6]+=1;
			/*if (x2kf1 < 20) {
				nEntries[1] += 1;
				if(choice == 1){
					nEntries[2]+=1;
					if(eton<0.85 && eton>0.4){
						nEntries[3]+=1;
						if (ppkf1[0] > 10 && ppkf1[1] > 10) {
							nEntries[4] += 1;
							if(nc<=3 && nn<=7){
								nEntries[5]+=1;
								if (central(0, 1) && central(1, 1) && unionVertex(0, 1, 1)) {
									nEntries[6] += 1;
									if (eGammas(1) < 0.2 ) {
										nEntries[7] += 1;
										if(eKPi(1)>0.1){
											nEntries[8] += 1;
											if(thetaxBetween(0, 30, 150, 1) && thetaxBetween(1, 30, 150, 1)){
												nEntries[9] += 1;
												if (((kf11 || ppkf1[1] < 300) && (kf10 || ppkf1[0] < 300))) {
													nEntries[10] += 1;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}*/

			if (x2kf1<999)
				hChi2->Fill(x2kf1);
			/*if (x2kf1<999 && x2kf1<x2kf2){
				pChi2->Fill(ppkf1[0],x2kf1);
			}
			if (x2kf2<999 && x2kf2<x2kf1){
				pChi2->Fill(ppkf2[0],x2kf2);
			}*/
			hetonT->Fill(eton);
			if(x2kf1<15 && beam<700 && nc<=3 && nn<=7){
				hEgammas->Fill(eGammas(1));
				hd0->Fill(d0[ipkf1[1]-1]);
				hd0->Fill(d0[ipkf1[0]-1]);
				hz0->Fill(z0[ipkf1[1]-1]);
				hz0->Fill(z0[ipkf1[0]-1]);
				hdz0->Fill(z0[ipkf1[0]-1]- z0[ipkf1[1] - 1]);
			}
			if(x2kf2<=15 && beam<700 && nc<=3 && nn<=7)
				hEgammas->Fill(eGammas(2));
			if (x2ikf1 < 50 && (goodKkpEntry1() && ((kf10 || ppkf1[0] < 300) && (kf11 || ppkf1[1] < 300))) && ks0InvMass(1) > 400 && ks0InvMass(1) < 600 && eton > 0.4 && eton < 0.85 && nc >= 2) {
				hDeDx->Fill(ppkf1[0],dExn[ipkf1[0]-1]/605.);
			}
			if (x2ikf2 < 50 && (goodKkpEntry2() && ((kf20 || ppkf2[0] < 300) && (kf21 || ppkf2[1] < 300))) && ks0InvMass(2) > 400 && ks0InvMass(2) < 600 && eton > 0.4 && eton < 0.85 && nc >= 2) {
				hDeDx->Fill(ppkf2[0],dExn[ipkf2[0]-1]/605.);
				//if(region[ipkf2[1]-1] == 1 && amplitude[ipkf2[1]-1]>0.4 && amplitude[ipkf2[1]-1] < 500 && ppkf2[1]>400){
				//	pampphi->Fill(phikf2[1],amplitude[ipkf2[1]-1]);
				//	pampEphi->Fill(phikf2[1],amplitudeE[ipkf2[1]-1]);
				//}
				if (choice == 2){
					hEgammas->Fill(eGammas(2));
					hEKPi->Fill(eKPi(2));
					hd0->Fill(d0[ipkf2[1] - 1]);
					hd0->Fill(d0[ipkf2[0] - 1]);
					hz0->Fill(z0[ipkf2[1] - 1]);
					hz0->Fill(z0[ipkf2[0] - 1]);
					hdz0->Fill(z0[ipkf2[0] - 1] - z0[ipkf2[1] - 1]);
					htheta->Fill(theta[ipkf2[0] - 1] * 180. / 3.14159);
					htheta->Fill(theta[ipkf2[1] - 1] * 180. / 3.14159);

					hInvMass->Fill(ks0InvMass(2));
				}
			}
		}
	}
	cout << "linecount" << linecount << endl;
	ofile.close();
	cout << all << "	" << good1 << "	" << good2 << "	" << all - good1 - good2 << endl;
	cout << selected << "	" << all1 << "	" << (double)selected / (double)all1 << endl << endl;

	for (size_t i = 0; i < 7; i++)
		cout << (double)nEntries[i]/(double)nEntries[0]*100. << endl;

	hInvMass->Write("hInvMass");
	hChi2->Write("hChi2");
	pEffChi2->Write("pEffChi2");
	hetonT->Write("hetonT");
	hEgammas->Write("hEgammas");
	hEKPi->Write("hEKPi");
	hEton->Write("hEton");
	hd0->Write("hd0");
	hz0->Write("hz0");
	hdz0->Write("hdz0");
	htheta->Write("htheta");
	hthetaKsP0->Write("hthetaKsP0");
	hp0invmass->Write("p0invmass");
	hDeDx->Write("hdEdx");
	hchoice->Write("hchoice");
	pampphi->Write("pampphi");
	pampphi2->Write("pampphi2");
	pampphi3->Write("pampphi3");
	pampphi4->Write("pampphi4");
	pampEphi->Write("pampEphi");
	for (size_t i = 0; i < eranges.size(); i++){
		hegmax[i]->Write(Form("hegmax%d",i));
		hegmean[i]->Write(Form("hegmean%d",i));
	}
	MyFile->Close();

	vector<double> tempRandEff;
	for(size_t i = 0; i<energiesEXP.size(); i++){
		cout << energiesEXP[i] << "	";
		cout << countAllEventsBeam[i] << "	";
		cout << countSelEventsBeam[i] << endl;
		if(countAllEventsBeam[i]!=0)
			tempRandEff.push_back((double)countSelEventsBeam[i]/(double)countAllEventsBeam[i] * 100);
		else
			tempRandEff.push_back(0.);
	}
	//cin.get();
	effvspa->Draw();
	/*TGraph* gr = new TGraph(energiesEXP.size(),&energiesEXP[0],&tempRandEff[0]);
	gr->SetTitle("detection efficiency;beam energy, MeV");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->Draw("AP");*/
}

double fitfunction(double *x, double *par) {
	return par[0] + x[0] * par[1] + pow(x[0],2.) * par[2];
}

//just draw final efficiency/cross section graphs using manual array filling
void draweffen() {
	double en[12] = { 1300, 1438, 1509, 1550, 1625, 1678, 1725, 1765, 1800, 1854, 1910, 1975};
	double enerr[12] = { 80, 35, 10, 25, 25, 2.5, 15, 8, 3, 20, 25, 20};
	double eff[12] = { 0.0209077, 0.0424553, 0.0411024, 0.0375454, 0.0345898, 0.0352175, 0.035, 0.0353802, 0.033944, 0.0325318, 0.0294251, 0.0273688};
	double efferr[12] = { 0.000209908, 0.000397473, 0.000562783, 0.000544652, 0.000376667, 0.00054754, 0.001, 0.000412039, 0.000591606, 0.00027107, 0.000176018, 0.00031863 };
	double luminosity[12] = { 10381.2, 3955.99, 1655.25, 1139.1, 3285.64, 1377.89, 5839.5, 5193.13, 2247.29, 9534.01, 32089.3, 9284.82 };
	double nEntr[12] = { 0., 53.6182, 21.5545, 65.2062, 159.964, 65.2508, 132.773, 119.826, 43.4956, 71.6009, 242.103, 26.9307 };
	double nEntrErr[12] = { 3.36557, 12.6933, 10.121, 12.1374, 19.7295, 14.8155, 37.0766, 32.2492, 20.9728, 38.5198, 39.1053, 16.3883 };
	double crossSect[12], crossSectErr[12];
	for (size_t i = 0; i < 12; i++) {
		crossSect[i] = nEntr[i] / eff[i] / luminosity[i];
		crossSectErr[i] = crossSect[i] * (nEntrErr[i] / nEntr[i] + efferr[i] / eff[i]);
	}
	double en2011[14] = {1436.82, 1500, 1525, 1550, 1575, 1600, 1625, 1650, 1675, 1700, 1725, 1762.1, 1840.2, 1947.87};
	double enErr2011[14] = {28.7801, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12.4935, 26.3546, 42.8463};
	double crossSect2011[14] = {0.281965, 0.895183, 1.32283, 1.83505, 2.36817, 2.80429, 3.02308, 2.98551, 2.75274, 2.42461, 2.08032, 1.63036, 0.986036, 0.532670};
	double crossSectErr2011[14] = {0.281965, 0.895183, 1.32283, 1.83505, 2.36817, 2.80429, 3.02308, 2.98551, 2.75274, 2.42461, 2.08032, 1.63036, 0.986036, 0.532670};
	//TGraphErrors* gr = new TGraphErrors(12,en,eff,enerr,efferr);
	TGraphErrors* gr = new TGraphErrors(12,en,crossSect,enerr,crossSectErr);
	gr->SetTitle("TGraphErrors Example");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->Draw("AP");
}

//fitting cross section usefull stuff
double PhSp2(double en, double m1, double m2) {
	if (en <= (m1 + m2))
		return 0.;
	else {
		double tm2 = m1 / en;
		double tm3 = m2 / en;
		return 0.5 * en * sqrt((1.0 - pow(tm2 - tm3,2)) * (1.0 - pow(tm2 + tm3,2)));
	}
}

double funcBWK(double en) {
	int n = 100;
	double msKc = 495.;
	double mpi = 135.;
	double h = (pow(en - msKc,2) - pow(msKc+mpi,2)) / n;
	double mKst = 892.;
	double gKst = 50.;
	double pi = 3.14159;
	double sum = 0;
	for (size_t i = 0; i <= (size_t)n; i++) {
		double m = sqrt( pow(msKc+mpi,2) + h * double(i));
		double grec = gKst * pow((PhSp2(m,mpi,msKc)/PhSp2(mKst,mpi,msKc)),3) * pow(mKst/m,2);
		sum = sum + mKst * grec / pi / (pow( pow(m, 2) - pow(mKst, 2), 2) + pow(grec * m, 2)) * pow(PhSp2(en, m, msKc), 3) / pow(PhSp2(en, mKst, msKc), 3);
	}
	double m = sqrt( pow(msKc+mpi,2) + h * double(0));
	double grec = gKst * pow((PhSp2(m,mpi,msKc)/PhSp2(mKst,mpi,msKc)),3) * pow(mKst/m,2);
	sum = sum - mKst * grec / pi / (pow( pow(m, 2) - pow(mKst, 2), 2) + pow(grec * m, 2)) * pow(PhSp2(en, m, msKc), 3) / pow(PhSp2(en, mKst, msKc), 3) / 2.0;

	m = sqrt( pow(msKc+mpi,2) + h * double(n));
	grec = gKst * pow((PhSp2(m,mpi,msKc)/PhSp2(mKst,mpi,msKc)),3) * pow(mKst/m,2);
	sum = sum - mKst * grec / pi / (pow( pow(m, 2) - pow(mKst, 2), 2) + pow(grec * m, 2)) * pow(PhSp2(en, m, msKc), 3) / pow(PhSp2(en, mKst, msKc), 3) / 2.0;

	return sum * h;
}

double crSAmplitude(double en){
	double phphi1680 = 0.;
	double msPhi1680 = 1680.;
	double gPhi1680 = 250.;
	double sPhi1680 = 6.;

	double phrho1450 = 0.;
	double msRho1450 = 1550.;
	double gRho1450 = 300.;
	double sRho1450 = 2.0;

	double pi = 3.14159;
	double deg2rad = pi/180.;
	double mKst = 892.;
	double msKc = 495.;
	double meVnb = 0.389379; //Calculated cross sections are often given in terms of electronvolts, via the conversion h2c2/GeV2 = 0.3894 mb
	double wdFKstKc1 = 0.;
	pair<double,double> ampPhi1680 = make_pair(0.,0.);
	double wdRKstKc1 = 0.;
	pair<double,double> ampRho1450 = make_pair(0.,0.);

	pair<double,double> xf1 = make_pair(cos(deg2rad*phphi1680),sin(deg2rad*phphi1680));
	pair<double,double> xr1 = make_pair(cos(deg2rad*phrho1450),sin(deg2rad*phrho1450));

	double wsp0F1=PhSp2(msPhi1680,mKst,msKc);
	if (wsp0F1 != 0){
		wdFKstKc1 = gPhi1680*pow((PhSp2(en,mKst,msKc)/wsp0F1),3) * pow(msPhi1680/en,2);
	}
	if (wdFKstKc1 != 0){
		double c = -pow(msPhi1680,2) + pow(en,2);
		double d = +(en * wdFKstKc1);
		double deter = pow(c,2) + pow(d,2);
		double a = xf1.first * c + xf1.second * d;
		double b = -xf1.first * d + xf1.second * c;
		double amp = sqrt(sPhi1680 * pow(msPhi1680,5) * gPhi1680 * wdFKstKc1 / (12.*pi) / meVnb);

		ampPhi1680 = make_pair(amp/deter * a, amp/deter * b);
	}

	double wsp0R1=PhSp2(msRho1450,mKst,msKc);
	if (wsp0R1 != 0)
		wdRKstKc1 = gRho1450*pow((PhSp2(en,mKst,msKc)/wsp0R1),3) * pow(msRho1450/en,2);
	if (wdRKstKc1 != 0){
		double c = pow(msRho1450,2) - pow(en,2);
		double d = -(en * wdRKstKc1);
		double deter = pow(c,2) + pow(d,2);
		double a = xr1.first * c + xr1.second * d;
		double b = -xr1.first * d + xr1.second * c;
		double amp = sqrt(sRho1450 * pow(msRho1450,5) * gRho1450 * wdRKstKc1 / (12.*pi) / meVnb);

		ampRho1450 = make_pair(amp/deter * a, amp/deter * b);
	}
		//ampRho1450 = sqrt(sRho1450 * pow(msRho1450,5) * gRho1450 * wdRKstKc1 / (12.*pi) / meVnb) * xr1 / dcmplx(msRho1450**2-e22,-(en*wdRKstKc1));


	pair<double,double> amp = make_pair(ampPhi1680.first + ampRho1450.first, ampPhi1680.second + ampRho1450.second);
	double ppcs;
	ppcs=12.*pi*meVnb* (pow(amp.first,2) + pow(amp.second,2))/pow(en,3) * funcBWK(en);
	if(ppcs <= 0)
		ppcs = 0.;
	return ppcs;

}



void drawBWK() {
	//TF1* fa1 = new TF1("fa1", "funcBWK(x)", 1400, 2000);
	TF1* fa1 = new TF1("fa1", "crSAmplitude(x)", 1400, 2000);
	fa1->Draw();
}

//function for drawing fitted with snd RADCOR library cross section, _thr -1 , _exp - 2; can change what to draw with chosing push_back argument
//definition of values can be found here https://www.inp.nsk.su/images/preprint/1999_103.pdf
void drawfittcross(char* infile1, char* infile2) {
	double enc1[15] = {1289.76, 1436.8, 1500, 1525, 1550, 1575, 1600, 1625, 1650, 1675, 1700, 1725, 1762.1, 1840.2, 1947.9};
	double encerr1[15] = {58.99, 28.78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12.49, 26.35, 42.85};
	double crossc1[15] = {0, 0.339469, 0.475308, 1.59304, 2.33625, 1.9604, 2.22457, 2.87796, 3.09978, 2.76731, 2.50083, 2.05777, 2.20684, 1.05526, 0.39898};
	double crosscerr1[15] = {0, 0.162064, 0.0652043, 0.209812, 0.243483, 0.265761, 0.314915, 0.307244, 0.276515, 0.302423, 0.33667, 0.251882, 0.283739, 0.208137, 0.126557};
	
	double enc2[14] = {1436.8, 1500, 1525, 1550, 1575, 1600, 1625, 1650, 1675, 1700, 1725, 1762.1, 1840.2, 1947.9 };
	double encerr2[14] = {28.78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12.49, 26.35, 42.85 };
	double crossc2[14] = {0.180663, 0.98034, 1.60153, 2.64386, 3.51816, 3.50217, 3.25667, 2.99725, 3.11927, 2.41277, 2.46896, 2.22681, 1.21504, 0.382014 };
	double crosscerr2[14] = {0.0766402, 0.112167, 0.237348, 0.303451, 0.311275, 0.335448, 0.293908, 0.340387, 0.458649, 0.372432, 0.326841, 0.291852, 0.334811, 0.0801421 };
	double babarEn[] = {};
	ifstream myfile(infile1);
	ifstream myfile1(infile2);
	string line;
	int np;
	double e, cre, cr, c, cl;
	getline(myfile, line);
	cout << line << endl;
	vector<double> en, cross, enerr, crosserr;
	while(myfile >> np >> e >> cre >> cr >> c >> cl){
		if (np == 1 && c != 0) {
			cout << np << "	" << e << "	" << cre << endl;
			en.push_back(e);
			cross.push_back(cr/c);
			enerr.push_back(0);
			crosserr.push_back(0);
		}
	}
	double erre, cexp, ecup, ecdn, rcr, rc;
	vector<double> en1, cross1, enerr1, crosserr1;
	vector<double> en12, cross12, enerr12, crosserr12;
	vector<double> en17, cross17, enerr17, crosserr17;
	getline(myfile1, line);
	while (myfile1 >> np >> e >> erre >> cexp >> ecup >> ecdn >> rcr >> rc) {
		if (np == 1) {
			cout << np << "	" << e << "	" << cexp << endl;
			en1.push_back(e);
			cross1.push_back(cexp*rc);
			enerr1.push_back(erre);
			crosserr1.push_back((ecup+ecdn)/2.*rc);
		}
		if (np == 2) {
			cout << np << "	" << e << "	" << cexp << endl;
			en12.push_back(e);
			cross12.push_back(cexp*rc);
			enerr12.push_back(erre);
			crosserr12.push_back((ecup+ecdn)/2.*rc);
		}
		if (np == 3) {
			cout << np << "	" << e << "	" << cexp << endl;
			en17.push_back(e);
			cross17.push_back(cexp*rc);
			enerr17.push_back(erre);
			crosserr17.push_back((ecup+ecdn)/2.*rc);
		}
	}
	myfile.close();
	myfile1.close();
	/*
	for(size_t i = 0; i < 14; i++){
		cross1[i] = crossc2[i]/crossc1[i];
		crosserr1[i] = sqrt(pow(crosscerr2[i]/crossc1[i],2) + pow(crosscerr1[i]/crossc1[i],2));
	}*/

	/*TGraphErrors* grrad = new TGraphErrors(en.size(), &en[0], &cross[0], &enerr[0], &crosserr[0]);
	grrad->SetTitle("radiative correction");
	grrad->SetName("gr1");
	grrad->SetMarkerColor(2);
	grrad->SetMarkerStyle(21);
	grrad->Draw("APL");*/

	ifstream babarf("/work/users/kladov/snd2k/R007-001/table.dat");
	vector<double> enb, crossb, enerrb, crosserrb;
	getline(babarf, line);
	while (babarf >> e >> erre >> erre >> cexp >> ecdn >> ecdn) {
		cout << e << "	" << cexp << endl;
		enb.push_back(e*1000);
		crossb.push_back(cexp);
		enerrb.push_back((erre-e)*1000);
		crosserrb.push_back(fabs(ecdn));
	}
	babarf.close();

	TGraphErrors* gr1 = new TGraphErrors(en1.size(), &en1[0], &cross1[0], &enerr1[0], &crosserr1[0]);
	gr1->SetTitle("2011");
	gr1->SetName("gr1");
	gr1->SetMarkerColor(2);
	gr1->SetMarkerStyle(21);
	TGraphErrors* gr12 = new TGraphErrors(en12.size(), &en12[0], &cross12[0], &enerr12[0], &crosserr12[0]);
	gr12->SetTitle("2012");
	gr12->SetName("gr12");
	gr12->SetMarkerColor(3);
	gr12->SetMarkerStyle(22);
	TGraphErrors* gr17 = new TGraphErrors(en17.size(), &en17[0], &cross17[0], &enerr17[0], &crosserr17[0]);
	gr17->SetTitle("2017");
	gr17->SetName("gr17");
	gr17->SetMarkerColor(4);
	gr17->SetMarkerStyle(23);
	TGraphErrors* grb = new TGraphErrors(enb.size(), &enb[0], &crossb[0], &enerrb[0], &crosserrb[0]);
	grb->SetTitle("babar");
	grb->SetName("grb");
	grb->SetMarkerColor(2);
	grb->SetMarkerStyle(21);
	
	TGraphErrors* gr = new TGraphErrors(15, enc1, crossc1, encerr1, crosscerr1);
	gr->SetTitle("2011");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	grb->Draw("AP");
	gr12->Draw("Psame");
	gr17->Draw("Psame");
	gr1->Draw("Psame");
	//gr->Draw("Psame");
	TCanvas* can = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
	legend->AddEntry("gr1","2011","lep");
	legend->AddEntry("gr12","2012","lep");
	legend->AddEntry("gr17","2017","lep");
	legend->AddEntry("grb","BaBar","lep");
	legend->Draw();
	can->Update();

	cout << "double enc2[" << en1.size() << "] = {";
	for(size_t i = 0; i < en1.size(); i++){
		cout << en1[i] << ", ";
	}
	cout << "};" << endl << "double encerr2[" << en1.size() << "] = {";
	for(size_t i = 0; i < en1.size(); i++){
		cout << enerr1[i] << ", ";
	}
	cout << "};" << endl << "double crossc2[" << en1.size() << "] = {";
	for(size_t i = 0; i < en1.size(); i++){
		cout << cross1[i] << ", ";
	}
	cout << "};" << endl << "double crosscerr2[" << en1.size() << "] = {";
	for(size_t i = 0; i < en1.size(); i++){
		cout << crosserr1[i] << ", ";
	}
	cout << "};" << endl;
}



double crsBB(double x) {
	if (x < 1400)
		return 0.1;
	if (x >= 1400 && x < 1640)
		return 4. / 240. * (x - 1400);
	if (x >= 1640 && x < 1850)
		return 3.5 / 210. * (1850 - x) + 0.5;
	else
		return 0.5;
}

double fouPiCrosssect(double en){
	if (en<1470.)
		return 6.7 + (33.-6.7)/470. * (en-1000.);
	if (en>1470.)
		return 33. - 23./(1900.-1470.) * (en-1470.);
}
double twok2piCrosssect(double en){
	if(en<1600.)
		return 0;
	else
		return 0.7/400.*(en-1600);
}
double ksKPiPi0Crosssect(double en){
	if(en<1550.)
		return 0;
	else
		return 2./350.*(en-1550);
}


//make ks0 inv mass distr for modeling and experiment in energy ranges "energygaps", fit them - > find efficiency, number of enries, integrated luminosity in range, store exp hists with fit to the file "spectres", display output arrays for fit
void kkpDistr::calculateEffAndNoE() {
	vector<double> efficiencyArr;
	vector<double> efficiencyErrArr;
	vector<double> entriesExpArr;
	vector<double> entriesExpErrArr;
	vector<double> luminosityArr;
	vector<double> luminosityErrArr;
	vector<double> luminosityTArr;
	vector<double> luminosityErrTArr;
	vector<double> meanEnArr;
	vector<double> meanEnTArr;
	vector<double> meanEnErrArr;
	vector<double> meanEnErrTArr;

	vector<double> luminosityEn;
	vector<double> luminosity;
	vector<double> luminosityErr;
	
	vector<double> meanEff;
	vector<double> meanEffErr;

	vector<double> nEntrieset[3][5][12];
	vector<double> nEntriesetErr[3][5][12];
	vector<double> chi2NDF[3][5][12];



	double crosssectvar[14][5] = { 0.64336, 3.00583, 2.96723, 2.16193, 0.670964 ,
									0.700812, 2.81528, 2.78754, 2.10175, 0.628002 ,
									0.630943, 3.00127, 2.85307, 2.10808, 0.599928 ,
									0.68644, 2.80075, 2.65093, 2.02698, 0.604753 ,
									0.613889, 3.07613, 2.94808, 2.13294, 0.640811 ,
									0.674373, 2.86997, 2.74815, 2.08579, 0.632016 ,
									0.568987, 2.73546, 2.55269, 1.96181, 0.561838 ,
									0.618841, 2.64016, 2.49995, 2.0366, 0.589529 ,
									0.671849, 2.88517, 2.77172, 2.23393, 0.543576 ,
									0.719923, 2.70113, 2.82625, 2.09024, 0.552712 ,
									0.574229, 2.82879, 2.63593, 2.03702, 0.569287 ,
									0.612098, 2.72402, 2.57228, 2.08331, 0.606202 ,
									0.50644, 2.64921, 2.59956, 1.94529, 0.428749 ,
									0.548634, 2.52885, 2.54577, 1.83083, 0.446011 };



	//double baseCrSect[9] = {0.131411, 0.825262, 1.78216, 2.82785, 2.71616, 2.6189, 2.45985, 2.63786, 1.07576};
	//double baseCrSect[9] = {0.131718, 0.786979, 1.82338, 2.82674, 2.74565, 2.67606, 2.30665, 2.60589, 1.10024};
	//double baseCrSect[9] = {0.0396138, 0.86488, 2.38423, 3.08366, 2.96012, 2.88019, 1.54994, 2.25272, 0.452545};
	double baseCrSect[12] = {0.227621, 1.09922, 2.03196, 1.23158, 3.14614, 2.58523, 0.755855, 0.704579, 0.721263, 0.885269, 0.29374, 0.366902};
	//double baseCrSect[9] = {0.164876, 0.802952, 2.08926, 2.67804, 2.65615, 2.7048, 2.30513, 2.52432, 1.10234};
	//double baseCrSectErr[9] = {0.054685, 0.176039, 0.302482, 0.37705, 0.496256, 0.414917, 0.427455, 0.373468, 0.160493};
	double baseCrSectErr[12] = {0.056702, 0.308974, 0.442277, 0.37505, 0.440877, 0.491549, 0.144883, 0.171468, 0.161611, 0.188838, 0.17491, 0.129518};

	//double entries4pi[9] = {97.0802, 78.461, 39.599, 36.0705, 19.1618, 39.6009, 37.5789, 43.4471, 115.638};
	//double entries4pi[9] = {170.632, 77.4467, 47.6104, 47.2324, 23.8192, 41.7883, 44.9313, 52.5379, 136.799};
	//double entries4pi[9] = {677.148, 197.395, 114.85, 95.2913, 46.5738, 70.5837, 60.3931, 71.3482, 163.433};
	double entries4pi[9] = {677.148, 197.395, 114.85, 95.2913, 46.5738, 70.5837, 60.3931, 71.3482, 163.433};
	//double entrieskskpp0[9] = {0, 0, 0, 6.59587, 7.57025, 20.6056, 31.4171, 58.7507, 352.948};
	//double entrieskskpp0[15] = {0, 0, 0, 0, 0, 0.400469, 0.711906, 1.30864, 1.66758, 2.03551, 2.42558, 3.14834, 7.17717, 19.2752, 46.5599};
	//double entrieskskpp0[15] = {0, 0, 0, 0, 0, 0.321543, 0.574611, 1.05575, 1.34471, 1.64731, 1.958, 2.53098, 5.73339, 15.1881, 35.7772};
	//double entrieskskpp0[15] = {0, 0, 0, 0, 0, 0.70363, 1.22432, 2.13675, 2.58815, 2.98438, 3.40088, 4.22852, 8.993, 20.4276, 39.4013};
	double entrieskskpp0[8] = { 0, 0, 0, 1.24716, 1.35717, 3.52858, 5.25401, 9.9807 };
	//double entrieskskpp0[8] = {0, 0, 2.52579, 5.09155, 3.42369, 7.46554, 16.5081, 32.8375};
	//double entrieskskpp0[12] = {0, 0, 0, 1.88623, 4.83941, 4.22317, 48.4439, 34.5174, 37.4401, 34.2308, 25.0402, 38.9704};
	//double entries2kc2pi[9] = {0, 0, 0, 0, 0.115448, 0.589727, 1.0897, 2.5459, 23.5336};
	double entries2kc2pi[9] = {0, 0, 0, 0, 0.017435, 0.0952805, 0.212486, 0.60844, 6.8576};

	double n4pidivbyeff[15] = {67520.3, 52298.4, 19495.1, 12707, 14531.7, 12986.4, 11541.6, 11913.5, 10413.1, 9710.55, 9129.25, 9514.04, 16660.4, 22297, 21495.8};
	double n4pidivbyeff80[15] = {66440, 51711.4, 19030.4, 12537.3, 14283.7, 12712.7, 11305.1, 11361.1, 9800.86, 9436.95, 8774.02, 9445.55, 15777.2, 21115.9, 19761.7};
	double n4pidivbyeff25[15] = {63051.1, 49750.4, 18168.4, 11948.8, 13796.6, 12529.4, 10491.8, 11108.7, 8884.81, 8981.18, 8359.99, 8508.77, 14665.8, 18941.3, 18021.9};
	double n4pidivbyeff25loose[15] = {67477.7, 52806, 19162.6, 12481, 15059.5, 13204.7, 11174.9, 11885.1, 9639.05, 9170.37, 9149.74, 9187.46, 16556.4, 20411.3, 19680.3};
	double n4pidivbyeff80loose[15] = {70724.8, 54809.8, 20082.2, 13086.7, 15565.3, 13356.9, 12055.9, 12220.3, 10747.7, 9831.4, 9526.4, 10198.6, 17538.3, 22771.3, 21532.7};
	
	
	double ne4pi[8] = { 182746, 47122.1, 35309, 30476.1, 15053.7, 26248.4, 27176.7, 34556.6 };



	double effbase[15] = {0.0253735, 0.0327809, 0.0347007, 0.0326093, 0.0317334, 0.0309954, 0.0288474, 0.0287814, 0.0275634, 0.0273235, 0.0256675, 0.0251995, 0.0232206, 0.0199766, 0.0158896};
	double effbase2012[14] = {0,0.176216, 0.140637, 0.130117, 0.099278, 0.0973381, 0.0899982, 0.0777184, 0.0774785, 0.0722186, 0.0651387, 0.0655787, 0.0591988, 0.0591588};
	for(size_t i = 0; i < 15; i++){
		cout << n4pidivbyeff80loose[i]/n4pidivbyeff25loose[i] << "	" << (n4pidivbyeff80loose[i]/n4pidivbyeff25loose[i])/(n4pidivbyeff80[i]/n4pidivbyeff25[i]) << endl;
	}

	ifstream ifile("/work/users/kladov/snd2k/R007-001/2011/luminosity.dat");
	double thnentr = 0;
	double lbord = 1200;
	double hbord = 0;
	vector<pair<double, double> > eranges;
	vector<double> enmod;
	vector<double> effmod;
	for(size_t i = 0; i < sizeof(energiesMod2011)/sizeof(energiesMod2011[0]); i++){
		enmod.push_back(energiesMod2011[i]);
		effmod.push_back(efficienciesMod2011[i]);
	}
	while (ifile.get() != EOF) {
		double a, b, c;
		ifile >> a >> b >> c;
		luminosityEn.push_back(a);
		luminosity.push_back(b);
		luminosityErr.push_back(c/3.);
		cout << a*2 << "	" << b << "	" << c << endl;
		pair<size_t,size_t> modind = findclosest(a,enmod);
		double efftemp = effmod[modind.first];
		if (modind.first != modind.second)
			efftemp += (effmod[modind.second] - effmod[modind.first]) * (a - enmod[modind.first]) / (enmod[modind.second] - enmod[modind.first]);
		thnentr += b * crsBB(a * 2) * efftemp/10.*3.;
		if (thnentr >= 60) {
			thnentr = 0;
			hbord = 2*a + 0.1;
			eranges.push_back(make_pair(lbord, hbord));
			lbord = hbord;
		}
	}
	if (eranges.back().second < 1980) {
		eranges.push_back(make_pair(eranges.back().second, 2010));
	}
	//eranges.clear();
	//eranges.push_back(make_pair(1530, 1760));
	ifile.close();
	for (size_t k = 0; k < eranges.size(); k++) {
		cout << eranges[k].first << "-----" << eranges[k].second << endl;
	}
	//for (size_t k = 0; k < luminosityEn.size(); k++)
	//for (size_t k = 0; k < 12; k++)
		//eranges.push_back(make_pair(2. * luminosityEn[k] - 30, 2. * luminosityEn[k] + 30));
		//eranges.push_back(make_pair(energygaps[k], energygaps[k+1]));
	
	vector<TH1*> hkskpe = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2011/kkp_KsNFn/*.root", "invks04", "ks0 inv mass distr", 60, 497.611 - 265, 497.611 + 265, eranges, "exp");
	//vector<TH1*> hkskpe = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2012/kkp_KsNF/*.root", "invks04", "ks0 inv mass distr", 60, 497.611 - 265, 497.611 + 265, eranges, "exp");
	//vector<TH1*> hkskpe = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2017/kkp_KsNFnn/*.root", "invks04", "ks0 inv mass distr", 60, 497.611 - 265, 497.611 + 265, eranges, "exp");
	vector<TH1*> hksinvmass[5];
	for(size_t k = 0; k < eranges.size(); k++) {
		for(size_t i = 0; i < 5; i++) {
			cout << k << "	" << i << endl;
			hksinvmass[i].push_back((TH1*)hkskpe[k]->Clone(Form("hksinvmass%d",i)));
		}
	}
	cout << "b" << endl;
	//vector<TH1*> hkskpe = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/exp2011_ks0notfixed/*.root", "invks04", "ks0 inv mass distr", 100, 497.611 - 265, 497.611 + 265, eranges, "exp");
	vector<TH1*> hkskp = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2011IA/*.root", "invks03", "ks0 inv mass distr", 80, 497.611 - 265, 497.611 + 265, eranges, "mod");
	//vector<TH1*> hkskp = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2012IA/*.root", "invks03", "ks0 inv mass distr", 80, 497.611 - 265, 497.611 + 265, eranges, "mod");
	//vector<TH1*> hkskp = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2017IA/*.root", "invks03", "ks0 inv mass distr", 80, 497.611 - 265, 497.611 + 265, eranges, "mod");
	//vector<TH1*> hkskp = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/4pi/2011A/*x.root", "invks03", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, eranges, "mod");
	//vector<TH1*> hkskp = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/kskpp0/*x.root", "invks03", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, eranges, "mod");
	//vector<TH1*> hkskp = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/2kc2p0/*x.root", "invks03", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, eranges, "mod");
	//vector<TH1*> hkskp = imvmass("/work/users/kladov/snd2k/R007-002/output/ntuples/kskleta/*x.root", "invks03", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, eranges, "mod");
	//vector<TH1*> hkskp = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/all_KsNF/et2p0g_wrc**x.root", "invks03", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, eranges, "mod");
	
	vector<double> noleff;
	TFile* fout = new TFile("/work/users/kladov/snd2k/R007-001/2011/spectres.root", "RECREATE");
	//TFile* fout = new TFile("/work/users/kladov/snd2k/R007-001/2011/spectres4pi.root", "RECREATE");
	//TFile* fout = new TFile("/work/users/kladov/snd2k/R007-001/2011/spectresKsKpipi0.root", "RECREATE");
	//TFile* fout = new TFile("/work/users/kladov/snd2k/R007-001/2011/spectres2kc2p0.root", "RECREATE");
	//TFile* fout = new TFile("/work/users/kladov/snd2k/R007-001/2011/spectreskskleta.root", "RECREATE");

	//______couting and drawing efficiencies in mod beam points
	{
		vector<double> effModErr;
		cout << "double energiesMod2011[" << energiesMod.size() << "] = {";
		for(size_t i = 0; i < energiesMod.size(); i++){
			cout << energiesMod[i];
			if(i<energiesMod.size()-1)
				cout << ", ";
			else
				cout << "};" << endl;
			effModErr.push_back(sqrt(efficienciesModBeam[i]*(1.-efficienciesModBeam[i])/countModEventsBeam[i]));
			noleff.push_back(0);
			//if(effbase2012[i]!=0){
			//	efficienciesModBeam[i]=(efficienciesModBeam[i]/effbase2012[i]-1)*100.;
			//	effModErr[i] = effModErr[i]/effbase2012[i]*100.;
			//}
		}
		cout << "double efficienciesMod2011[" << energiesMod.size() << "] = {";
		for(size_t i = 0; i < energiesMod.size(); i++){
			cout << efficienciesModBeam[i];
			if(i<energiesMod.size()-1)
				cout << ", ";
			else
				cout << "};" << endl;
		}
		TGraphErrors* gre = new TGraphErrors(energiesMod.size(), &energiesMod[0], &efficienciesModBeam[0], &noleff[0], &effModErr[0]);
		gre->SetTitle("detection efficiency;beam energy, MeV");
		gre->SetMarkerColor(4);
		gre->SetMarkerStyle(21);
		gre->Draw("AP");
		gre->Write("efficiency");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->SetFillColor(0);
		c->SetFrameFillColor(0);
		c->SetBorderMode(0);
		c->SetFrameBorderMode(0);
		c->SetGrid();
		c->Update();
		cin.get();
	}
	//________recalculating of efficiency with lum weights (for eranges)
	
	
	vector<double> nEntries4pi;
	vector<double> nEntries4piErr;
	vector<double> popr4pi;
	vector<double> popr4piErr;
	for (size_t k = 0; k < eranges.size(); k++) {
		double lum1 = 0;
		double modEventsSum = 0;
		int countPointsInErange = 0;
		double entriesNumber = 0;
		pair<size_t,size_t>  closestInd;
		meanEff.push_back(0);
		meanEffErr.push_back(sqrt((double)hkskp[k]->GetEntries() * (1. - (double)hkskp[k]->GetEntries()/(double)countModEvents[k])) / (double)countModEvents[k]);
		if((double)hkskp[k]->GetEntries()==0){
			meanEffErr[k]=meanEffErr[k-1];
		}
		nEntries4pi.push_back((double)hkskpe[k]->GetEntries());
		nEntries4piErr.push_back(sqrt((double)hkskpe[k]->GetEntries() * (1 - (double)hkskpe[k]->GetEntries()/(double)countExpEvents[k])) / (double)countExpEvents[k]);
		if(nEntries4piErr[k]==0){
			nEntries4piErr[k]=nEntries4piErr[k-1];
		}
		//cout << meanEffErr[k] << "	" << nEntries4piErr[k] << endl;
		double sumWeightPoints = 0.;
		for (size_t j = 0; j < luminosityEn.size(); j++) {
			if (2 * luminosityEn[j] >= eranges[k].first && 2 * luminosityEn[j] < eranges[k].second) {
				countPointsInErange+=1;
				lum1 += luminosity[j];
				closestInd = findclosest(luminosityEn[j], energiesMod);
				//cout << closestInd.first << "	ASDSAFASFAS	" << closestInd.second << endl;
				double tempEffPoint = efficienciesModBeam[closestInd.first];
				if(closestInd.first != closestInd.second)
					tempEffPoint += (efficienciesModBeam[closestInd.second] - efficienciesModBeam[closestInd.first]) * (luminosityEn[j] - energiesMod[closestInd.first]) / (energiesMod[closestInd.second] - energiesMod[closestInd.first]);
				//meanEff.back() += luminosity[j] * efficienciesModBeam[closestInd] / countModEventsBeam[closestInd];

				double tempWeightPoint = luminosity[j] * tempEffPoint * crsBB(2 * luminosityEn[j]);
				sumWeightPoints += tempWeightPoint;

				meanEff.back() += tempWeightPoint * tempEffPoint;
				//entriesNumber += luminosity[j] * tempEffPoint * twok2piCrosssect(2.*luminosityEn[j]);
				//entriesNumber += luminosity[j] * tempEffPoint * fouPiCrosssect(2.*luminosityEn[j]);
				//entriesNumber += luminosity[j] * tempEffPoint * 0.5;
				entriesNumber += luminosity[j] * tempEffPoint * ksKPiPi0Crosssect(2.*luminosityEn[j])*3./10.;
			}
		}
		//meanEff.back() = meanEff.back() / (double)countPointsInErange * modEventsSum / lum1;
		meanEff.back() = meanEff.back() / sumWeightPoints;
		//cout << "meaneff " << meanEff.back() << "	" << lum1 << "	" << k << " / " << eranges.size() << endl;
		//cout << entriesNumber << ", ";
		popr4pi.push_back(nEntries4pi[k]/meanEff[k] / ne4pi[k]);
		popr4piErr.push_back(sqrt(2.)*sqrt(pow(nEntries4piErr[k]/meanEff[k] / ne4pi[k],2) + pow(popr4pi[k]*meanEffErr[k]/meanEff[k],2) ) );
		cout << nEntries4pi[k]/meanEff[k]/ ne4pi[k] << ", ";
		//cout << nEntries4pi[k]/meanEff[k] << ", ";
	}
	cout << endl;
	//________luminosity erange calculation, energy as arithmetic mean, and then error as standart dispersion
	for (size_t k = 0; k < eranges.size(); k++) {
		double backnorm = 0;
		double enNorm = 0;
		pair<size_t, size_t>  closestInd;
		luminosityTArr.push_back(0.);
		luminosityErrTArr.push_back(0.);
		meanEnTArr.push_back(0.);
		meanEnErrTArr.push_back(0.);
		double sumWeightPoints = 0.;
		for (size_t j = 0; j < luminosityEn.size(); j++) {
			if (2 * luminosityEn[j] >= eranges[k].first && 2 * luminosityEn[j] < eranges[k].second) {
				closestInd = findclosest(luminosityEn[j], energiesMod);
				double tempEffPoint = efficienciesModBeam[closestInd.first];
				if (closestInd.first != closestInd.second)
					tempEffPoint += (efficienciesModBeam[closestInd.second] - efficienciesModBeam[closestInd.first]) * (luminosityEn[j] - energiesMod[closestInd.first]) / (energiesMod[closestInd.second] - energiesMod[closestInd.first]);
				double tempWeightPoint = luminosity[j] * tempEffPoint * crsBB(2 * luminosityEn[j]);
				sumWeightPoints += tempWeightPoint;

				luminosityTArr.back() += luminosity[j];
				luminosityErrTArr.back() += pow(luminosityErr[j],2);
				meanEnTArr.back() += 2 * luminosityEn[j] * tempWeightPoint;
				enNorm+=luminosity[j];
				backnorm+=meanEff[k]*luminosity[j]*fouPiCrosssect(luminosityEn[j]);

			}
		}
		meanEnTArr.back() = meanEnTArr.back() / sumWeightPoints;
		luminosityErrTArr.back() = sqrt(luminosityErrTArr.back());
		for (size_t j = 0; j < luminosityEn.size(); j++) {
			if (2 * luminosityEn[j] >= eranges[k].first && 2 * luminosityEn[j] < eranges[k].second) {
				closestInd = findclosest(luminosityEn[j], energiesMod);
				double tempEffPoint = efficienciesModBeam[closestInd.first];
				if (closestInd.first != closestInd.second)
					tempEffPoint += (efficienciesModBeam[closestInd.second] - efficienciesModBeam[closestInd.first]) * (luminosityEn[j] - energiesMod[closestInd.first]) / (energiesMod[closestInd.second] - energiesMod[closestInd.first]);
				double tempWeightPoint = luminosity[j] * tempEffPoint * crsBB(2 * luminosityEn[j]);
				meanEnErrTArr.back() += pow(2 * luminosityEn[j] - meanEnTArr.back(), 2) * tempWeightPoint;
			}
		}
		meanEnErrTArr.back() = sqrt(meanEnErrTArr.back() / sumWeightPoints);
	}

	

	TGraphErrors* gre = new TGraphErrors(eranges.size(), &meanEnTArr[0], &popr4pi[0], &meanEnErrTArr[0], &popr4piErr[0]);
	gre->SetTitle("detection efficiency vs ISR photon energy parameter %d dependence on 2 beam;  2 * beam energy, MeV");
	gre->SetMarkerColor(4);
	gre->SetMarkerStyle(21);
	gre->Draw("AP");
	//c->Update();

	//________drawing and fitting efficiency dependence on mce[1]
	/*for(size_t i = 0; i < 5; i++){
		//TGraphErrors* gre = new TGraphErrors(meanEnArr.size(), &meanEnArr[0], &pareffmce[i][0], &meanEnErrArr[0], &pareffmceErr[i][0]);
		TGraphErrors* gre = new TGraphErrors(energiesMod.size(), &energiesMod[0], &pareffmce[i][0], &noleff[0], &pareffmceErr[i][0]);
		gre->SetTitle(Form("detection efficiency vs ISR photon energy parameter %d dependence on 2 beam;  beam energy, MeV",i+1));
		gre->SetMarkerColor(4);
		gre->SetMarkerStyle(21);
		gre->Draw("AP");
		gre->Write("efficiency");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		c->SetFillColor(0);
		c->SetFrameFillColor(0);
		c->SetBorderMode(0);
		c->SetFrameBorderMode(0);
		c->SetGrid();
		c->Update();
		

		TF1* f1 = new TF1("f1", "[0]+[1]*x", 650, 1000);
		f1->SetParameters(100.,0.1);
		gre->Fit("f1", "", "", 700, 1000);
		c->Update();
		cout << f1->GetParameter(0) << " + " << f1->GetParameter(1) << " * x" << endl; 
		cin.get();
	}*/


	
	//files for background getting
	TFile* f2kc2p0 = new TFile("/work/users/kladov/snd2k/R007-001/2011/spectres2kc2p0.root");
	TFile* f4pi = new TFile("/work/users/kladov/snd2k/R007-001/2011/spectres4pin.root");
	//TFile* fkskpipi0 = new TFile("/work/users/kladov/snd2k/R007-001/2011/spectres4pin.root");
	TFile* fkskpipi0 = new TFile("/work/users/kladov/snd2k/R007-001/2011/spectresKsKpipi0all.root");
	
	//for (size_t k = 0; k < luminosityEn.size(); k++) {//12
	for (size_t k = 0; k < eranges.size(); k++) {//12
	//for (size_t k = 7; k < 8; k++) {//12
		cout << "a" << endl;
		//TH1* hkskp = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp_KsNF/2017/*.root", "invks03", "ks0 inv mass distr", 200, 497.611 - 265, 497.611 + 265, 2.*luminosityEn[k]-2.5, 2.*luminosityEn[k] + 2.5, "mod");
		cout << hkskp[k]->GetEntries() << endl;
		cout << hkskpe[k]->GetEntries() << endl;
		int h4pisource = k;
		if(h4pisource == 4)
			h4pisource = 3;
		if(h4pisource == 7 || h4pisource>8)
			h4pisource = 8;
		TH1* h2kc2p0 = (TH1*)f2kc2p0->Get(Form("modeling%d",h4pisource));
		TH1* h4pi = (TH1*)f4pi->Get(Form("modeling%d",h4pisource));
		TH1* hkskpipi0 = (TH1*)fkskpipi0->Get(Form("modeling%d",0));
		//if (hkskp[k]->GetEntries() != 0 && hkskpe[k]->GetEntries() >= 0 && k != 4 && k != 7) { //20
		
		if (hkskpe[k]->GetEntries() > 0) { //20
			double entriesInKsKPihist = (double)hkskpipi0->GetEntries();
			h2kc2p0->Scale(entries2kc2pi[h4pisource]/(double)h2kc2p0->GetEntries());
			h4pi->Scale(entries4pi[h4pisource]/(double)h4pi->GetEntries());
			hkskpipi0->Scale(entrieskskpp0[k]/entriesInKsKPihist);

			//if(h2kc2p0->GetEntries()>0)
			//	hkskpipi0->Add(h2kc2p0,1);
			if(hkskpipi0->GetEntries()>0)
				h4pi->Add(h4pi,hkskpipi0,0,1);
			//h4pi->Rebin(2);
			h4pi->Draw();
			TF1* f2 = new TF1("f2", "[0]+[1]*x+[2]*x*x", 300, 700);
			//h4pi->Fit("f2", "", "", 300, 750);
			vector<double> parBack(3);
			//f2->GetParameters(&parBack[0]);

			fout->cd();
			TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
			gStyle->SetOptTitle(kFALSE);
			gStyle->SetOptStat(0);
			c->Update();
			//cin.get();

			//TH1* hkskp = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp_KsNF/2011/*.root", "invks03", "ks0 inv mass distr", 200, 497.611 - 265, 497.611 + 265, k, "mod");
			//Fit for kskp modeling pik with 3 gauss, fixed parameters
			//TF1* f1 = new TF1("f1", "[0]*(exp(-0.5*((x-[1]-[9])/[2])**2)+[3]*exp(-0.5*((x-[4]-[9])/[5])**2)+[6]*exp(-0.5*((x-[7]-[9])/[8])**2))*(265./100.)/53.3188", 300, 700);
			TF1* f1 = new TF1("f1", "([0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)+[6]*exp(-0.5*((x-[7])/[8])**2))*(265./40.)", 200, 700);
			f1->SetParameters(50.0, 497.611, 20., 25., 500, 50., 25., 495, 50.);
			f1->SetParLimits(1, 400., 550.);
			f1->SetParLimits(4, 450., 550.);
			f1->SetParLimits(7, 450., 550.);
			f1->SetParLimits(0, 0., 1000.);
			f1->SetParLimits(3, 0., 1000.);
			f1->SetParLimits(6, 0., 1000.);
			f1->SetParLimits(2, 10., 50.);
			f1->SetParLimits(5, 10., 50.);
			f1->SetParLimits(8, 10., 50.);
			//f1->FixParameter(1, 490.5);
			//f1->FixParameter(2, 14.57);
			//f1->FixParameter(3, 0.113);
			//f1->FixParameter(4, 461.3);
			//f1->FixParameter(5, 8.503);
			//f1->FixParameter(6, 0.182);
			//f1->FixParameter(7, 493.8);
			//f1->FixParameter(8, 31.54);
			//f1->SetParameter(9, 0);
			//f1->SetParLimits(9, -10, 10);
			vector<double> par(9);
			hkskp[k]->Draw();

			//set dedault paameters if no modeling in this point
			if (hkskp[k]->GetEntries() == 0) {
				f1->FixParameter(0, 1.);
				f1->FixParameter(1, 490.5);
				f1->FixParameter(2, 14.57);
				f1->FixParameter(3, 0.113);
				f1->FixParameter(4, 461.3);
				f1->FixParameter(5, 8.503);
				f1->FixParameter(6, 0.182);
				f1->FixParameter(7, 493.8);
				f1->FixParameter(8, 31.54);
			}
			hkskp[k]->Fit("f1", "", "", 350, 650);
			f1->GetParameters(&par[0]);
			f1->SetParameter(0, max(par[0], max(par[3], par[6])));
			f1->SetParameter(3, min(par[0], min(par[3], par[6])));
			f1->SetParameter(6, min(par[0], min(par[3], par[6])));
			hkskp[k]->Fit("f1", "", "", 350, 650);
			hkskp[k]->Write(Form("modeling%d",k));
			//h4pi->Draw("same");
			
			f1->GetParameters(&par[0]);
			//change the order, with maximum coefficient outside
			if (par[0] < par[3] || par[0] < par[6]) {
				cout << "switchig" << endl;
				if (par[3] < par[6]) {
					double p0t = par[0];
					par[0] = par[6];
					par[3] = par[3] / par[0];
					par[6] = p0t / par[0];

					double p1t = par[7];
					par[7] = par[1];
					par[1] = p1t;

					double p2t = par[8];
					par[8] = par[2];
					par[2] = p2t;
				}
				else {
					double p0t = par[0];
					par[0] = par[3];
					par[6] = par[6] / par[0];
					par[3] = p0t / par[0];

					double p1t = par[4];
					par[4] = par[1];
					par[1] = p1t;

					double p2t = par[5];
					par[5] = par[2];
					par[2] = p2t;
				}
			}
			else {
				par[3] = par[3] / par[0];
				par[6] = par[6] / par[0];
			}

			double nEntries = par[0];
			double multiplier = (fabs(par[2]) + par[3] * fabs(par[5]) + par[6] * fabs(par[8])) * (sqrt(3.14159 * 2.));
			nEntries = nEntries * multiplier;
			nEntries = hkskp[k]->GetEntries();
			//double	effic = nEntries / (double)countModEvents[k];
			double	effic = meanEff[k];
			double	effErr = meanEffErr[k];
			double multiplierErr = f1->GetParError(2) + f1->GetParError(5) * fabs(par[3]) + f1->GetParError(3) * fabs(par[5]) + f1->GetParError(8) * fabs(par[6]) + f1->GetParError(6) * fabs(par[8]);
			//double effErr = (f1->GetParError(0) * multiplier + multiplierErr * par[0]) / (double)countModEvents;
			//double multiplier = hkskp->GetEntries()/nEntries;
			c->Update();
			//cin.get();
			//TH1* hkskpe = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/MHAD2017/kkp_KsNF/*.root", Form("invks%d",k), "ks0 inv mass distr", 20, 497.611 - 265, 497.611 + 265, 2.*luminosityEn[k] - 2.5, 2.*luminosityEn[k] + 2.5, "exp");

			//TH1* hkskpe = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/exp2011_ks0notfixed/*.root", "invks03", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, k, "exp");
			//fit for experiment
			//TF1* f1 = new TF1("f1", "[0]*(exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)+[6]*exp(-0.5*((x-[7])/[8])**2)) + [9]*([10]*exp(-0.5*((x-[11])/[12])**2)+[13]*exp(-0.5*((x-[14])/[15])**2)+[16]+[17]*x) + [18]*([19]*exp(-0.5*((x-[20])/[21])**2)+[22]*exp(-0.5*((x-[23])/[24])**2)) + [25]*(exp(-0.5*((x-[26])/[27])**2)+[28]*exp(-0.5*((x-[29])/[30])**2))", 300, 700);
			TF1* f1e = new TF1("f1e", "[0]*(exp(-0.5*((x-[1]-[12])/([2]/[13]))**2)+[3]*exp(-0.5*((x-[4]-[12])/([5]/[13]))**2)+[6]*exp(-0.5*((x-[7]-[12])/([8]/[13]))**2))*(265./30.)/1. + [14]*([9]+[10]*x+[11]*x*x) + [15]*[16]*(exp(-0.5*((x-512.473)/25.6123)**2)+1.668832*exp(-0.5*((x-554.803)/105.413)**2)+0.1819033*exp(-0.5*((x-344.423)/27.6483)**2))", 200, 800); //300 700
			TF1* f2e = new TF1("f2e", "[0]*(exp(-0.5*((x-[1]-[12])/[2])**2)+[3]*exp(-0.5*((x-[4]-[12])/[5])**2)+[6]*exp(-0.5*((x-[7]-[12])/[8])**2))*(265./30.)/1. + [9]+[10]*x+[11]*x*x + [13]*x*x*x", 200, 800); //300 700
			TF1* f3e = new TF1("f3e", "[0]*(exp(-0.5*((x-[1]-[12])/[2])**2)+[3]*exp(-0.5*((x-[4]-[12])/[5])**2)+[6]*exp(-0.5*((x-[7]-[12])/[8])**2))*(265./30.)/1. + [9]+[10]*x+[11]*x*x + [13]*exp(-((x-[14])/[15])**2)", 200, 800); //300 700
			//TF1* f1 = new TF1("f1", "[0]*(exp(-0.5*((x-[1]-[25])/[2])**2)+[3]*exp(-0.5*((x-[4]-[25])/[5])**2)+[6]*exp(-0.5*((x-[7]-[25])/[8])**2))*(265./20.)/53.3188 + [9]*([10]*exp(-0.5*((x-[11])/[12])**2)+[13]*exp(-0.5*((x-[14])/[15])**2)+[16]+[17]*x) + [18]*([19]*exp(-0.5*((x-[20])/[21])**2)+[22]*exp(-0.5*((x-[23])/[24])**2))", 300, 700);
			//f1->FixParameter(1, 490.5);
			//f1->FixParameter(2, 14.57);
			//f1->FixParameter(3, 0.113);
			//f1->FixParameter(4, 461.3);
			//f1->FixParameter(5, 8.503);
			//f1->FixParameter(6, 0.182);
			//f1->FixParameter(7, 493.8);
			//f1->FixParameter(8, 31.54);
			//f1->SetParameter(14, 0);
			//f1->SetParLimits(14, -10, 10);

			for (size_t i = 1; i < 9; i++)
				f1e->FixParameter(i, par[i]);
			//for (size_t i = 9; i < 12; i++)
			//	f1e->FixParameter(i, parBack[i-9]);
			f1e->FixParameter(14,1.);
			f1e->SetParLimits(11,-0.007,0.000001);
			f1e->SetParameter(12, 0);
			f1e->SetParLimits(12, -2, 2);
			f1e->SetParLimits(13, 0.97, 1.03);
			f1e->SetParLimits(0, 0, 500);
			double par16 =  265./40. * 87.904 * entrieskskpp0[k]/entriesInKsKPihist * 3./2.;
			f1e->FixParameter(16,par16);
			f1e->SetParameter(15,1);
			//f1e->SetParLimits(15,-1e-10,par15*1.1);
			f1e->SetParLimits(15,0.05-(1e-10),1.1);
			//f1e->FixParameter(15,0);
			
			{
				for (size_t i = 1; i < 9; i++)
					f2e->FixParameter(i, par[i]);
				f2e->SetParLimits(11,-0.00005,0.00005);
				f2e->SetParameter(12, 0);
				f2e->SetParLimits(12, -10, 10);
				f2e->SetParLimits(0, 0, 500);

				for (size_t i = 1; i < 9; i++)
					f3e->FixParameter(i, par[i]);
				f3e->SetParLimits(11,-0.00005,0.00005);
				f3e->SetParameter(12, 0);
				f3e->SetParLimits(12, -10, 10);
				f3e->SetParLimits(0, 0, 500);
				f3e->SetParameter(13,4);
				f3e->SetParLimits(13,0,10);
				f3e->SetParameter(14,350);
				f3e->SetParLimits(14,200,400);
				f3e->SetParameter(15,50);
			}

			//f1->FixParameter(10, 1.13285e+03);
			//f1->FixParameter(11, 4.95439e+02);
			//f1->FixParameter(12, 1.47172e+02);
			//f1->FixParameter(13, 1.37301e+03);
			//f1->FixParameter(14, 2.69848e+02);
			//f1->FixParameter(15, 1.11795e+02);
			//f1->FixParameter(16, -2.53885e+03);
			//f1->FixParameter(17, 3.18818e+00);

			//f1->FixParameter(19, 6.59567e+02);
			//f1->FixParameter(20, 3.64866e+02);
			//f1->FixParameter(21, 7.08758e+01);
			//f1->FixParameter(22, 5.58503e+02);
			//f1->FixParameter(23, 4.61473e+02);
			//f1->FixParameter(24, 1.03831e+02);

			//f1->SetParLimits(25,0,100);
			//f1->FixParameter(26, 503.203);
			//f1->FixParameter(27, 107.813);
			//f1->FixParameter(28, 0.728);
			//f1->FixParameter(29, 516.542);
			//f1->FixParameter(30, 27.6604);

			int rebin[5];
			//rebin[0] = (int)((600. / (double)(hkskpe[k]->GetEntries()))*((pow((eranges[k].first-1550.)/300.,2)+1))) + 1;
			//rebin[1] = (int)((600. / (double)(hkskpe[k]->GetEntries()))*((pow((eranges[k].first-1600.)/250.,2)+1))) + 1;
			//rebin[2] = (int)((400. / (double)(hkskpe[k]->GetEntries()))*((pow((eranges[k].first-1550.)/300.,2)+1))) + 1;
			//rebin[3] = (int)((500. / (double)(hkskpe[k]->GetEntries()))*((pow((eranges[k].first-1550.)/300.,2)+1))) + 1;
			//rebin[4] = (int)((600. / (double)(hkskpe[k]->GetEntries()))*((pow((eranges[k].first-1550.)/200.,2)+1))) + 1;
			rebin[0] = 3;
			rebin[1] = 2;
			rebin[2] = 5;
			rebin[3] = 4;
			rebin[4] = (int)((600. / (double)(hkskpe[k]->GetEntries()))*((pow((eranges[k].first-1550.)/300.,2)+1))) + 1;
			//int rebin = 3;
			for(size_t i = 0; i < 5; i++){
				hksinvmass[i][k]->Rebin(rebin[i]);
			}
			hkskpe[k]->Rebin(rebin[0]);
			//hkskpe[k]->Add(h4pi,-0.9);
			//hksinvmass[0]->Rebin(rebin[0]);
			//h4pi->Rebin(rebin[0]/8);


			vector<double> pare(17);
			vector<double> pare2(14);
			vector<double> pare3(16);
			hkskpe[k]->SetMarkerStyle(8);
			hkskpe[k]->SetMarkerSize(1);
			hkskpe[k]->Draw("PE1");
			//h4pi->Scale(4./3.);
			//hkskpe[k]->Fit("f1e", "", "", 300 + (meanEnArr[k]-meanEnArr[0])/6., 700);
			hkskpe[k]->Fit("f1e", "", "", 300, 700);
			hkskpe[k]->Draw("PE1");
			//h4pi->Draw("same");
			f1e->GetParameters(&pare[0]);
			double multipliere = (fabs(pare[2]) + pare[3] * fabs(pare[5]) + pare[6] * fabs(pare[8])) * (sqrt(3.14159 * 2.) / pare[13]);
			double nEntriese = pare[0] / (double)rebin[0];
			double parecout = pare[0];
			double rebincout = rebin[0];
			nEntriese = nEntriese * multipliere;
			entriesExpArr.push_back(nEntriese);
			entriesExpErrArr.push_back(f1e->GetParError(0) * multiplier / (double)rebin[0]);
			
			

			/*
			for(size_t j = 0; j < 5; j++){
				for(size_t i = 0; i < 4; i++){
					hksinvmass[j][k]->Fit("f1e", "", "", 200+i*50, 750);
					f1e->GetParameters(&pare[0]);
					nEntrieset[0][j][i+0].push_back(pare[0] / (double)rebin[j] * multiplier);
					nEntriesetErr[0][j][i+0].push_back(f1e->GetParError(0) * multiplier / (double)rebin[j]);
					chi2NDF[0][j][i+0].push_back(f1e->GetChisquare()/(double)f1e->GetNDF());

					hksinvmass[j][k]->Fit("f1e", "", "", 350, 800-i*50);
					f1e->GetParameters(&pare[0]);
					nEntrieset[0][j][i+4].push_back(pare[0] / (double)rebin[j] * multiplier);
					nEntriesetErr[0][j][i+4].push_back(f1e->GetParError(0) * multiplier / (double)rebin[j]);
					chi2NDF[0][j][i+4].push_back(f1e->GetChisquare()/(double)f1e->GetNDF());

					hksinvmass[j][k]->Fit("f1e", "", "", 200+i*50, 800-i*50);
					f1e->GetParameters(&pare[0]);
					nEntrieset[0][j][i+8].push_back(pare[0] / (double)rebin[j] * multiplier);
					nEntriesetErr[0][j][i+8].push_back(f1e->GetParError(0) * multiplier / (double)rebin[j]);
					chi2NDF[0][j][i+8].push_back(f1e->GetChisquare()/(double)f1e->GetNDF());

					hksinvmass[j][k]->Fit("f2e", "", "", 200+i*50, 750);
					f2e->GetParameters(&pare2[0]);
					nEntrieset[1][j][i+0].push_back(pare2[0] / (double)rebin[j] * multiplier);
					nEntriesetErr[1][j][i+0].push_back(f2e->GetParError(0) * multiplier / (double)rebin[j]);
					chi2NDF[1][j][i+0].push_back(f2e->GetChisquare()/(double)f2e->GetNDF());

					hksinvmass[j][k]->Fit("f2e", "", "", 350, 800-i*50);
					f2e->GetParameters(&pare2[0]);
					nEntrieset[1][j][i+4].push_back(pare2[0] / (double)rebin[j] * multiplier);
					nEntriesetErr[1][j][i+4].push_back(f2e->GetParError(0) * multiplier / (double)rebin[j]);
					chi2NDF[1][j][i+4].push_back(f2e->GetChisquare()/(double)f2e->GetNDF());

					hksinvmass[j][k]->Fit("f2e", "", "", 200+i*50, 800-i*50);
					f2e->GetParameters(&pare2[0]);
					nEntrieset[1][j][i+8].push_back(pare2[0] / (double)rebin[j] * multiplier);
					nEntriesetErr[1][j][i+8].push_back(f2e->GetParError(0) * multiplier / (double)rebin[j]);
					chi2NDF[1][j][i+8].push_back(f2e->GetChisquare()/(double)f2e->GetNDF());

					hksinvmass[j][k]->Fit("f3e", "", "", 200+i*50, 750);
					f3e->GetParameters(&pare3[0]);
					nEntrieset[2][j][i+0].push_back(pare3[0] / (double)rebin[j] * multiplier);
					nEntriesetErr[2][j][i+0].push_back(f3e->GetParError(0) * multiplier / (double)rebin[j]);
					chi2NDF[2][j][i+0].push_back(f3e->GetChisquare()/(double)f3e->GetNDF());

					hksinvmass[j][k]->Fit("f3e", "", "", 350, 800-i*50);
					f3e->GetParameters(&pare3[0]);
					nEntrieset[2][j][i+4].push_back(pare3[0] / (double)rebin[j] * multiplier);
					nEntriesetErr[2][j][i+4].push_back(f3e->GetParError(0) * multiplier / (double)rebin[j]);
					chi2NDF[2][j][i+4].push_back(f3e->GetChisquare()/(double)f3e->GetNDF());

					hksinvmass[j][k]->Fit("f3e", "", "", 200+i*50, 800-i*50);
					f3e->GetParameters(&pare3[0]);
					nEntrieset[2][j][i+8].push_back(pare3[0] / (double)rebin[j] * multiplier);
					nEntriesetErr[2][j][i+8].push_back(f3e->GetParError(0) * multiplier / (double)rebin[j]);
					chi2NDF[2][j][i+8].push_back(f3e->GetChisquare()/(double)f3e->GetNDF());
				}
			}*/
			
			

			c->Update();
			//TF1* f2 = new TF1("f2", "[0]+[1]*x+[2]*x*x", 300, 700);
			for (size_t i = 0; i < 3; i++)
				f2->FixParameter(i, pare[9 + i]*pare[14]);
			f2->SetLineColor(2);
			TF1* f3 = new TF1("f3", "[0]+[1]*x+[2]*x*x + [3]*(exp(-0.5*((x-512.473)/25.6123)**2)+1.668832*exp(-0.5*((x-554.803)/105.413)**2)+0.1819033*exp(-0.5*((x-344.423)/27.6483)**2))", 300, 700);
			for (size_t i = 0; i < 3; i++)
				f3->FixParameter(i, pare[9 + i]*pare[14]);
			f3->FixParameter(3, pare[15]*pare[16]);
			f3->SetLineColor(4);
			
			f2->Draw("Lsame");
			f3->Draw("Lsame");
			//h4pi->Draw("same");
			//hkskp[k]->DrawNormalized("same",backnorm);
			c->Update();
			//cin.get();
			

			hkskpe[k]->Write(Form("experiment%d",k));
			c->Write(Form("sp%d", k));
		
			//double nEntriese = pare[0];

			cout << "fitted" << endl;
		
			//to exclued bad areas with no modeling or no experiment entries
			//if (hkskp[k]->GetEntries() != 0 && hkskpe[k]->GetEntries() >= 20) {
			
			cout << "rebin " << rebincout << endl;
			cout << "pare[0] " << parecout << endl;
			//cout << "entries1 " << nEntrieset[0][3][3+0][k] << endl;

			cout << "multiplier check: " << multiplier << endl;
			//cout << "number of entries	" << nEntries << endl;
			cout << "Chi2 / NDF	" << f1->GetChisquare() << " / " << f1->GetNDF() << endl;

			cout << "multiplier check: " << multipliere << endl;
			cout << "Chi2 / NDF	" << f1e->GetChisquare() << " / " << f1e->GetNDF() << endl;

			cout << "efficiency	" << effic * 3. / 10. * 100. << endl;
			cout << "efficiencyErr	" << effErr * 3. / 10. * 100. << endl;

			cout << "number of entries	" << nEntriese << endl;
			cout << "number of entries error	" << f1e->GetParError(0) * multipliere << endl;

			//lum = 1476.2;
			//lumErr = 0.;
			cout << "luminosity	" << luminosityTArr[k] << endl;
			cout << "luminosity error	" << luminosityErrTArr[k] << endl;

			//cout << endl << energygaps[k] << " - " << energygaps[k + 1] << endl;
			cout << "cross section	" << ((double)nEntriese / effic / luminosityTArr[k]) * 10. / 3. << endl;

			cout << "energy Point	" << meanEnTArr[k] << endl;
			cout << "energy Point error	" << meanEnErrTArr[k] << endl;

			efficiencyArr.push_back(effic * 3. / 10.);
			efficiencyErrArr.push_back(effErr * 3. / 10.);
			meanEnArr.push_back(meanEnTArr[k]);
			meanEnErrArr.push_back(meanEnErrTArr[k]);
			luminosityArr.push_back(luminosityTArr[k]);
			luminosityErrArr.push_back(luminosityErrTArr[k]);
			//luminosityArr.push_back(lum);
			//luminosityErrArr.push_back(lumErr);
			cin.get();
			cout << "b" << endl;
		}
		else{
			cout << "asdasd" << endl;
		}
	}
	
	fout->Write();
	fout->Close();
	
	cout << meanEnArr.size() << endl;
	for (size_t i = 0; i < meanEnArr.size(); i++)
		//cout << fabs(efficiencyArr[i]/effbase[i]-1)*100 << ", ";
		cout << efficiencyArr[i] << ", ";
	cout << endl;
	for (size_t i = 0; i < meanEnArr.size(); i++)
		cout << efficiencyErrArr[i] << ", ";
	cout << endl;
	for (size_t i = 0; i < meanEnArr.size(); i++)
		cout << entriesExpArr[i] << ", ";
	cout << endl;
	for (size_t i = 0; i < meanEnArr.size(); i++)
		cout << entriesExpErrArr[i]/2. << ", ";
	cout << endl;
	for (size_t i = 0; i < meanEnArr.size(); i++)
		cout << luminosityArr[i] << ", ";
	cout << endl;
	for (size_t i = 0; i < meanEnArr.size(); i++)
		cout << luminosityErrArr[i] << ", ";
	cout << endl;
	for (size_t i = 0; i < meanEnArr.size(); i++)
		cout << meanEnArr[i] << ", ";
	cout << endl;
	for (size_t i = 0; i < meanEnArr.size(); i++)
		cout << meanEnErrArr[i] << ", ";
	cout << endl;
	vector<double> expCr;
	vector<double> expCr1;
	vector<double> expCrErr;
	vector<double> expCrErr1;
	for (size_t i = 0; i < meanEnArr.size(); i++) {
		expCr.push_back(entriesExpArr[i] / luminosityArr[i] / efficiencyArr[i]);
		expCr1.push_back((entriesExpArr[i] / luminosityArr[i] / efficiencyArr[i])/baseCrSect[i]);
		expCrErr.push_back(sqrt(pow(entriesExpErrArr[i]/entriesExpArr[i],2) + pow(luminosityErrArr[i]/luminosityArr[i],2) + pow(efficiencyErrArr[i]/ efficiencyArr[i],2))* expCr[i]);
		expCrErr1.push_back(sqrt(pow(entriesExpErrArr[i]/entriesExpArr[i],2) + pow(luminosityErrArr[i]/luminosityArr[i],2) + pow(efficiencyErrArr[i]/ efficiencyArr[i],2) + pow(baseCrSectErr[i]/baseCrSect[i],2)) * expCr[i]);
		cout << expCr[i] << ", ";
	}
	cout << endl;
	for (size_t i = 0; i < meanEnArr.size(); i++) {
		cout << expCrErr[i] << ", ";
	}
	cout << endl;
	
	/*for (size_t i = 0; i < meanEnArr.size(); i++) {
		expCrErr[i] = expCrErr[i] / crosssectvar[6][i];
	}
	for (size_t j = 0; j < 14; j++) {
		for (size_t i = 0; i < meanEnArr.size(); i++) {
			expCr[i] = crosssectvar[j][i] / crosssectvar[6][i];
		}
		TGraphErrors* gr = new TGraphErrors(meanEnArr.size(), &meanEnArr[0], &expCr[0], &meanEnErrArr[0], &expCrErr[0]);
		gr->SetTitle("cross section without radcor and some other corrections");
		gr->SetMarkerColor(4);
		gr->SetMarkerStyle(21);
		gr->Draw("AP");
		TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
		gStyle->SetOptTitle(kFALSE);
		gStyle->SetOptStat(0);
		c->Update();
		cin.get();
	}*/
	
	TGraphErrors* gr = new TGraphErrors(meanEnArr.size(), &meanEnArr[0], &expCr[0], &meanEnErrArr[0], &expCrErr[0]);
	gr->SetTitle("cross section without radcor and some other corrections");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->Draw("AP");
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	gStyle->SetOptTitle(kFALSE);
	gStyle->SetOptStat(0);
	c->Update();
	cin.get();
	
	/*TGraphErrors* gr1 = new TGraphErrors(meanEnArr.size(), &meanEnArr[0], &expCr1[0], &meanEnErrArr[0], &expCrErr1[0]);
	gr1->SetTitle("cross section without radcor and some other corrections");
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(21);
	gr1->Draw("AP");
	c->Update();*/
	/*vector<TH1*> hCrT;
	vector<double> meanCrSect;
	vector<double> meanCrSectErr;
	vector<double> nmean;
	vector<double> sumErr;
	vector<double> sumErr2;
	for (size_t k = 0; k < meanEnArr.size(); k++) {
		hCrT.push_back(new TH1F(Form("hCrT%d",k), "crSpread", 10, expCr[k]-max(0.5,expCr[k]/2.), expCr[k] + max(0.5,expCr[k]/2.)));
		meanCrSect.push_back(expCr[k]);
		meanCrSectErr.push_back(pow(expCrErr[k]/expCrErr[k],2));
		nmean.push_back(1.);
		sumErr.push_back(1./expCrErr[k]);
		sumErr2.push_back(pow(1./expCrErr[k],2));
	}
	for(size_t j = 0; j < 5; j++){
		for(size_t i = 0; i < 12; i++){
			for(size_t l = 0; l < 3; l++){
				vector<double> expCr1;
				vector<double> expCrErr1;
				for (size_t k = 0; k < meanEnArr.size(); k++) {
					expCr1.push_back(nEntrieset[l][j][i][k] / luminosityArr[k] / efficiencyArr[k]);
					expCrErr1.push_back((nEntriesetErr[l][j][i][k]/nEntrieset[l][j][i][k] + luminosityErrArr[k]/luminosityArr[k] + efficiencyErrArr[k]/ efficiencyArr[k])* expCr1[k]);
					if(chi2NDF[l][j][i][k] < 1.5){
						hCrT[k]->Fill(expCr1[k]);
						meanCrSect[k]+=expCr1[k]/expCrErr1[k];
						meanCrSectErr[k]+=pow(expCrErr1[k]/expCrErr1[k],2);
						nmean[k]+=1;
						sumErr[k]+=1./expCrErr1[k];
						sumErr2[k]+=pow(1./expCrErr1[k],2);
					}
				}
				TGraphErrors* gr1 = new TGraphErrors(meanEnArr.size(), &meanEnArr[0], &expCr1[0], &meanEnErrArr[0], &expCrErr1[0]);
				gr1->SetTitle("cross section without radcor and some other corrections");
				gr1->SetMarkerColor(4);
				gr1->SetMarkerStyle(21);
				if(chi2NDF[l][j][i][3] < 1.5)
					gr1->Draw("Psame");
			}
		}
	}
	c->Update();
	cin.get();
	for (size_t k = 0; k < meanEnArr.size(); k++) {
		hCrT[k]->Draw();
		TLine *baseLine = new TLine(expCr[k],0,expCr[k],hCrT[k]->GetMaximum());
		baseLine->Draw("same");
		c->Update();
		cin.get();
	}
	for (size_t k = 0; k < meanEnArr.size(); k++) {
		//meanCrSect[k] = meanCrSect[k]/nmean[k];
		//meanCrSect[k] = meanCrSect[k]/sumErr[k];
		meanCrSect[k] = hCrT[k]->GetMean()/baseCrSect[k];
		//meanCrSect[k] = hCrT[k]->GetMean();
		//meanCrSectErr[k] = sqrt(meanCrSectErr[k]/pow(sumErr[k],2));
		//meanCrSectErr[k] = hCrT[k]->GetRMS();
		meanCrSectErr[k] = expCrErr[k]/baseCrSect[k];
		//cout << meanCrSectErr[k]/expCr[k] << endl;
		cout << meanCrSect[k] << ", ";
	}
	cout << endl;
	//TGraphErrors* grr = new TGraphErrors(meanEnArr.size(), &meanEnArr[0], &expCr[0], &meanEnErrArr[0], &meanCrSectErr[0]);
	TGraphErrors* grr = new TGraphErrors(meanEnArr.size(), &meanEnArr[0], &meanCrSect[0], &meanEnErrArr[0], &meanCrSectErr[0]);
	grr->SetTitle("cross section without radcor and some other corrections");
	grr->SetMarkerColor(4);
	grr->SetMarkerStyle(21);
	grr->Draw("AP");

	c->Update();
	cin.get();*/
	
	/*
	f2kc2p0->Close();
	f4pi->Close();
	fkskpipi0->Close();*/
}




void drawCrossSection(){
	TF1* f1 = new TF1("f1","crsBB(x)",1000,2000);
	f1->Draw();
}


void go() {
	kkpDistr processing;
	processing.year = 2011;
	processing.definededpopr();
	//normcoeff = 671.395;
	//effEeEnergy();
	//luminosity();
	//copy();
	//someAnalisis();
	//effEnergy();
	
	//landaufit();
	//draweffen();
	//drawfittcross("pp_cs_thr.his", "pp_cs_exp.his");
	//drawBWK();
	//processing.calculateEffAndNoE();
	//processing.invmass2part();
	
	
	/*
	processing.genTableML("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2011IA/*.root","/work/users/kladov/snd2k/R007-001/2011/trainSet.txt",0);
	processing.genTableML("/work/users/kladov/snd2k/R007-002/output/ntuples/4pi/2011A/*x.root","/work/users/kladov/snd2k/R007-001/2011/trainSet.txt",1);
	processing.genTableML("/work/users/kladov/snd2k/R007-002/output/ntuples/2kc2p0/2011A/*x.root","/work/users/kladov/snd2k/R007-001/2011/trainSet.txt",2);
	ifstream ifile("/work/users/kladov/snd2k/R007-001/2011/luminosity.dat");
	double thnentr = 0;
	double lbord = 1200;
	double hbord = 0;
	vector<pair<double, double> > eranges;
	int f = 0;
	double en0 = 0.;
	vector<double> enmod;
	vector<double> effmod;
	for(size_t i = 0; i < sizeof(energiesMod2011)/sizeof(energiesMod2011[0]); i++){
		enmod.push_back(energiesMod2011[i]);
		effmod.push_back(efficienciesMod2011[i]);
	}
	while (ifile.get() != EOF) {
		double a, b, c;
		ifile >> a >> b >> c;
		if (f==0)
			en0 = a;
		cout << a * 2 << "	" << b << "	" << c << endl;
		pair<size_t,size_t> modind = findclosest(a,enmod);
		double efftemp = effmod[modind.first];
		if (modind.first != modind.second)
			efftemp += (effmod[modind.second] - effmod[modind.first]) * (a - enmod[modind.first]) / (enmod[modind.second] - enmod[modind.first]);
		thnentr += b * crsBB(a * 2) * efftemp/10.*3.;
		if (thnentr >= 60) {
			thnentr = 0;
			hbord = 2 * a + 0.1;
			eranges.push_back(make_pair(lbord, hbord));
			lbord = hbord;
		}
		f++;
	}
	if (eranges.back().second < 1980) {
		eranges.push_back(make_pair(eranges.back().second, 2010));
	}
	//eranges.clear();
	//eranges.push_back(make_pair(1530, 1760));
	ifile.close();
	//processing.drawdEvsp(eranges);
	processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2011/kkp_KsNFn/*.root", "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelDistrKsKPiE.root", "/work/users/kladov/snd2k/R007-001/2011/testSetExp.txt", 2 ,eranges, "exp");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2012/kkp_KsNF/*.root", "/work/users/kladov/snd2k/R007-001/kkpi/2012/SelDistrKsKPiE.root", "/work/users/kladov/snd2k/R007-001/2011/testSetExp.txt", 2 ,eranges, "exp");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/MHAD2017/kkp_KsNFnn/*.root", "/work/users/kladov/snd2k/R007-001/kkpi/2017/SelDistrKsKPiE.root", "/work/users/kladov/snd2k/R007-001/2011/testSetExp.txt", 2 ,eranges, "exp");
	processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2011IA/*.root", "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelDistrKsKPiM.root", "/work/users/kladov/snd2k/R007-001/2011/testSetMod.txt", 1 ,eranges, "mod");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2012IA/*.root", "/work/users/kladov/snd2k/R007-001/kkpi/2012/SelDistrKsKPiM.root", "/work/users/kladov/snd2k/R007-001/2011/testSetMod.txt", 1 ,eranges, "mod");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/kkp_KsNF/2017IA/*.root", "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelDistrKsKPiM.root", "/work/users/kladov/snd2k/R007-001/2011/testSetMod.txt", 1 ,eranges, "mod");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/4pi/2011A/*x.root", "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelDistrKsKPiM4pi1.root", "/work/users/kladov/snd2k/R007-001/2011/testSetMod.txt", 1 ,eranges, "mod");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/kskpp0/*x.root", "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelDistrKsKPiMkskpp0.root", "/work/users/kladov/snd2k/R007-001/2011/testSetMod.txt", 1 ,eranges, "mod");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/2kc2p0/*x.root", "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelDistrKsKPiM4pi.root", "/work/users/kladov/snd2k/R007-001/2011/testSetMod.txt", 1 ,eranges, "mod");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/kskleta/*x.root", "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelDistrKsKPiMkskleta.root", "/work/users/kladov/snd2k/R007-001/2011/testSetMod.txt", 1 ,eranges, "mod");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-001/output/ntuples/all_KsNF/eta3pi_wrc**x.root", "/work/users/kladov/snd2k/R007-001/kkpi/2011/SelDistreta3pi.root", "/work/users/kladov/snd2k/R007-001/2011/testSetMod.txt", 1 ,eranges, "mod");
	//processing.genTableMLExp("/work/users/kladov/snd2k/R007-002/output/ntuples/4pi/2011/*x.root", "/work/users/kladov/snd2k/R007-001/kkpi/SelDistrKsKPiB.root", "/work/users/kladov/snd2k/R007-001/2017/testSetBack.txt", 0 ,eranges, "mod");
	*/


	//compareSelectionCuts();
	luminCalculation lumcalc;
	//lumcalc.effEeEnergy("t1", "/work/users/kladov/snd2k/R007-001/2011/*x.root");
	lumcalc.luminosityCalc();
	//lumcalc.differenceForLum();
	//lumcalc.drawdEdxspectrum();
	//drawLum();
	//sumpoints();
	
	//hwf hk("etonkkp", "etonkkp", 50, 0.3, 0.8);
	//hwf he("etoneta", "etoneta", 50, 0.3, 0.8);
	
	//hk.fillhists("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp_ks0notfixed/*.root");
	//he.fillhists("/work/users/kladov/snd2k/R007-001/output/ntuples/2etg_ks0notfixed/*.root");
	/*
	hk.hist->Draw();
	he.hist->SetLineColor(2);
	he.hist->Draw("same");*/

	//TH1* hkp = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/*.root");
	//TH1* hkk = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/2kc2pi/*.root");
	//TH1* hpp = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/4pi_KsNF/*.root", "invks01", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, 1);
	//TH1* hpp2= imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/2kc2pi_KsNF/*.root", "invks02", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, 1);
	//TH1* hkskp= imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp_ks0notfixed/*.root", "invks03", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, 2);
	//TH1* hexp= imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp/exp2011_ks0notfixed/*.root", "invks04", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265, 3);
	
	/*TH1* hs = new TH1F("invks05", "ks0 inv mass distr", 30, 497.611 - 265, 497.611 + 265);

	hexp->Draw();

	hs->Add(hpp, hpp2, 30, 0.5);
	hs->Add(hkskp, 3);
	hs->SetLineColor(2);
	hs->Draw("same");

	hpp->Scale(30);
	hpp2->Scale(0.5);
	hkskp->Scale(3);
	hkskp->Draw("same");
	hpp2->SetLineColor(2);
	hpp2->Draw("same");
	hpp->SetLineColor(3);
	hpp->Draw("same");*/
	
	//TH1* hpp = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/all_KsNF/ksknppp0_wrc_nemcc_825-9426-27_wmix.root", "invks01", "ks0 inv mass distr", 25, 497.611 - 265, 497.611 + 265, 1);
	//hpp->Draw();


	//TH1* hkskp = imvmass("/work/users/kladov/snd2k/R007-001/output/ntuples/kkp_KsNF/2017/*.root", "invks03", "ks0 inv mass distr", 200, 497.611 - 265, 497.611 + 265, 4);
	
	//for kskpipi0
	/*TF1* f1 = new TF1("f1", "gaus(0)+gaus(3)", 230, 800);
	f1->SetParameters(50.0, 497.611, 20.0, 20, 500, 50);
	hkskp->Draw();
	hkskp->Fit("f1", "", "", 230, 800);
	cout << "Chi2 / NDF	" << f1->GetChisquare() << " / " << f1->GetNDF() << endl;*/

	//Fit for 4pi background or 2kc2pi
	/*TF1* f1 = new TF1("f1", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", 300, 700);
	f1->SetParameters(50,350,20,20,400,100);
	hkskp->Draw();
	hkskp->Fit("f1", "", "", 300, 700);
	cout << "Chi2 / NDF	" << f1->GetChisquare() << " / " << f1->GetNDF() << endl;*/
	

	/*hexp1->SetLineColor(1);
	hexp2->SetLineColor(2);
	hexp3->SetLineColor(3);
	hexp4->SetLineColor(4);
	hexp5->SetLineColor(5);
	hexp6->SetLineColor(6);
	TF1* f1 = new TF1("f1", "gaus(0)+[3]+[4]*(x-497.611)", 340, 740);
	f1->SetParameters(50.0, 497.611, 20.0, 50.0, -0.25);
	hexp1->Draw();
	TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
	//hexp1->Fit("f1", "", "", 340, 740);
	c->Update();
	cin.get();
	hexp2->Draw();
	//hexp2->Fit("f1", "", "", 340, 740);
	c->Update();
	cin.get();
	hexp3->Draw();
	//hexp3->Fit("f1", "", "", 340, 740);
	c->Update();
	cin.get();
	hexp4->Draw();
	//hexp4->Fit("f1", "", "", 340, 740);
	c->Update();
	cin.get();
	hexp5->Draw();
	//hexp5->Fit("f1", "", "", 340, 740);
	c->Update();
	cin.get();
	hexp6->Draw();
	//hexp6->Fit("f1", "", "", 340, 740);
	c->Update();*/
	
	//hexp->Fit("f1","","",340,740);
	//hexp->DrawNormalized("", 1);
	//cout << hkp->Integral() << "	" << hkk->Integral() << "	" << hpp->Integral() << "	" << hexp->Integral() << endl;*/
	
	//pscp D:\programs\kkp_model.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/kkpi/kkp_model.cpp
	//pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/kkpi/kkp_model.cpp D:\programs\kkp_model.cpp

	//pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/2017/testSetExp.txt D:\programs\testSetExp.txt 
	//pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/2017/testSetMod.txt D:\programs\testSetMod.txt
	//pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/2017/testSetBack.txt D:\programs\testSetBack.txt
	//pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/2017/trainSet.txt D:\programs\trainSet.txt
/*
pscp D:\workProgramms\kkp_model.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/kkpi/kkp_model.cpp
12081998Vkl

pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/2011/testSetExp.txt D:\workProgramms\testSetExp.txt
12081998Vkl
pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/2011/testSetMod.txt D:\workProgramms\testSetMod.txt
12081998Vkl

pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/2011/trainSet.txt D:\workProgramms\trainSet.txt
12081998Vkl

pscp D:\workProgramms\testMarksExp.txt kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/2011/testMarksExp.txt
12081998Vkl
pscp D:\workProgramms\testMarksMod.txt kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/2011/testMarksMod.txt
12081998Vkl

*/
}
