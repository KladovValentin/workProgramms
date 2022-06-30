#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TH1.h"
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
#include "TMath.h"
#include<algorithm> 
using namespace std;

double p[4];
double smin = 0., smax = 1.;
double fitpad[29] = { 0.1, 0.27, 0.33, 0.71,    0.78, 0.98, 1.04, 1.42,    1.48, 1.68, 1.75, 2.12,    2.21, 2.39, 2.45, 2.85,    2.92, 3.11, 3.17, 3.55,    3.61, 3.78, 3.85, 4.23,    \
	4.31, 4.5, 4.56, 4.93, 1.6};

double p3g0(double x, double xr, double sr);
double TMath::Erf(double a);
double TMath::Erf(double a);
void m4(double* x, double* y);
double p3g(double x, double xr, double sr);
double p3gl(double x, double xr, double sr, double er, double xer, double der);
double scphi(double x);




double agt0(double a) {
	//double w, a0;
	double agt00;
	/*vector<double> pix; //?
	vector<double> wp;
	pix.push_back(1.);
	pix.push_back(2.);
	pix.push_back(3.);
	if (wp[0] == 1) {
		w = pix[0];
		a0 = pix[1];
		agt00 = (w * a0 + a) / (1.0 - (1.0 - w) * exp(-a));
	}
	else {*/
		agt00 = a;
	//}
	return agt00;
}


void m4(double x[], double y[]){ //?
    double d;
    double x1, x2, x3, x4;
    double y1, y2, y3, y4;
    x1=x[0];
    x2=x[1];
    x3=x[2];
    x4=x[3];
    y1=agt0(fabs(y[0]));
    y2=agt0(fabs(y[1]));
    y3=agt0(fabs(y[2]));
    y4=agt0(fabs(y[3]));
    d=(x1-x2)*(x1-x3)*(x2-x3)*(x1-x4)*(x2-x4)*(x3-x4);
    p[0]=( x1*(x1-x3)*x3*(x1-x4)*(x3-x4)*x4*y2+ \
		x2*(pow(x4,2))*((-pow(x3,3))*y1+(pow(x3,2))*x4*y1+(pow(x1,2))*(x1-x4)*y3)+ \
		(pow(x1,2))*x2*(pow(x3,2))*(-x1+x3)*y4+ \
		(pow(x2,3))*(x4*(-(pow(x3,2))*y1+x3*x4*y1+x1*(x1-x4)*y3) + x1*x3*(-x1+x3)*y4)+ \
		(pow(x2,2))*(x1*x4*(-(pow(x1,2))+(pow(x4,2)))*y3 + (pow(x3,3))*(x4*y1-x1*y4)+x3*(-(pow(x4,3))*y1+(pow(x1,3))*y4)))/d;

    p[1]=( (pow(x1,2))*(x1-x4)*(pow(x4,2))*(y2-y3)+ \
		(pow(x3,3))*((pow(x4,2))*(y1-y2)+ (pow(x1,2))*(y2-y4))+ \
		(pow(x2,2))*((pow(x4,3))*(y1-y3)+ (pow(x1,3))*(y3-y4)+(pow(x3,3))*(-y1+y4))+ \
		(pow(x3,2))*((pow(x4,3))*(-y1+y2)+(pow(x1,3))*(-y2+y4))+ \
		(pow(x2,3))*((pow(x4,2))*(-y1+y3)+(pow(x3,2))*(y1-y4)+(pow(x1,2))*(-y3+y4)))/d;

    p[2]=(-x1*(x1-x4)*x4*(x1+x4)*(y2-y3)+ \
		x3*((pow(x4,3))*(y1-y2)+(pow(x1,3))*(y2-y4))+ \
		(pow(x3,3))*(-x4*y1-x1*y2+x4*y2+x1*y4)+ \
		(pow(x2,3))*(x4*y1+x1*y3-x4*y3-x1*y4+x3*(-y1+y4))+ \
		x2*((pow(x4,3))*(-y1+y3)+(pow(x3,3))*(y1-y4)+(pow(x1,3))*(-y3+y4)))/d;

    p[3]=(x1*(x1-x4)*x4*(y2-y3)+ \
		(pow(x3,2))*(x4*y1+x1*y2-x4*y2-x1*y4)+ \
		(pow(x2,2))*(-x4*y1-x1*y3+x4*y3+x3*(y1-y4)+x1*y4)+ \
		x2*((pow(x4,2))*(y1-y3)+(pow(x1,2))*(y3-y4)+(pow(x3,2))*(-y1+y4))+ \
		x3*((pow(x4,2))*(-y1+y2)+(pow(x1,2))*(-y2+y4)))/d;
}


double p3g(double x, double xr, double sr) {
	double p3g;
	double er, der;
	//common / rims / smin, smax;  //?
	//double smin, smax;
	er = 0.821691;
	der = 2.40385;
	if (((fabs(smin - 0) < 10) || (fabs(smin - 120) < 10) || (fabs(smin - 240) < 10) || (fabs(smin - 360) < 10)) && ((xr - smin) < der)) {
		p3g = p3gl(x, xr, sr, er, smin, der);
		return p3g;
	}
	if (((fabs(smax - 0) < 10) || (fabs(smax - 120) < 10) || (fabs(smax - 240) < 10) || (fabs(smax - 360) < 10)) && ((smax - xr) < der)) {
		p3g = p3gl(x, xr, sr, er, smax, -der);
		return p3g;
	}
	p3g = p3g0(x, xr, sr);
	return p3g;
}


double p3g0(double x, double xr, double sr){
    double p3g0;
    double p0,p1,p2,p3;
    double srr, arg, PI;
    PI=3.1415926535;
    p0=p[0];
    p1=p[1];
    p2=p[2];
    p3=p[3];
    srr=0.5+(sr-0.5)*(TMath::Erf((x-(xr-5.*sr))/sqrt(2.0)/sr)-TMath::Erf((x-(xr+5.*sr))/sqrt(2.0)/sr))/2.0;
    arg=(x-xr)/sqrt(2.0)/srr;
    p3g0=exp(-(pow(arg,2)))*srr/sqrt(2.0*PI)*(p1+p2*(x+xr)+p3*(x*x+x*xr+xr*xr+2*srr*srr))+(1.+TMath::Erf((arg)))/2.*(p0+x*(p1+x*(p2+p3*x))+(p2+3.*p3*x)*srr*srr);
    return p3g0;
}


double p3gl(double x, double xr, double sr, double er, double xer, double der){
    double p3gl;
    double p0,p1,p2,p3;
    double earg, arg, PI, srr;
    double part1, part2, part3, part4, part5;
    PI=3.1415926535;
    p0=p[0];
    p1=p[1];
    p2=p[2];
    p3=p[3];
    srr=0.5+(sr-0.5)*(TMath::Erf((x-(xr-5.*sr))/sqrt(2.0)/sr)-TMath::Erf((x-(xr+5.*sr))/sqrt(2.0)/sr))/2.0;
    srr=max(sr,0.001);
    arg=(x - xr)/(sqrt(2.0)*srr);
    earg=exp(-(pow(arg,2)));
    p3gl=(2.*(der*er*(sqrt(2.*PI)*(p0 + p2*(pow(srr,2) + pow(x,2)) + x*(p1 + p3*(3.*pow(srr,2) + pow(x,2)))) + srr*(p1 + p2*(x + xr) + p3*(2.*pow(srr,2) + pow(x,2) + x*xr + pow(xr,2)))*earg) - \
			(-1. + er)*(sqrt(2.*PI)* (p1*pow(srr,2) + 3.*p3*pow(srr,4) + p0*x + 3.*p2*pow(srr,2)*x + p1*pow(x,2) + 6.*p3*pow(srr,2)*pow(x,2) + p2*pow(x,3) + p3*pow(x,4) - \
			(p0 + p2*(pow(srr,2) + pow(x,2)) + x*(p1 + p3*(3.*pow(srr,2) + pow(x,2))))*xer) + \
			srr*(p0 + p1*(x - xer + xr) + p2*(2.*pow(srr,2) + pow(x,2) - x*xer + x*xr - xer*xr + pow(xr,2)) + \
			p3*(pow(x,2)*(x - xer) + x*(x - xer)*xr + (x - xer)*pow(xr,2) + pow(xr,3) + pow(srr,2)*(5.*x - 2.*xer + 3.*xr)))*earg)) + \
				sqrt(2.*PI)*(-(der*er*(p0 + p2*(pow(srr,2) + pow(x,2)) + x*(p1 + p3*(3.*pow(srr,2) + pow(x,2))))) + \
				(-1. + er)*(p1*pow(srr,2) + 3.*p3*pow(srr,4) + p0*x + 3.*p2*pow(srr,2)*x + p1*pow(x,2) + 6.*p3*pow(srr,2)*pow(x,2) + p2*pow(x,3) + p3*pow(x,4) - \
				(p0 + p2*(pow(srr,2) + pow(x,2)) + x*(p1 + p3*(3.*pow(srr,2) + pow(x,2))))*xer))*(1. - TMath::Erf(arg))) \
					/(2.*der*sqrt(2.*PI));
    return p3gl;
}




double scphi(double x) {
	double scphi;
	double fitfun;
	//common/hcfitd/fitpad,fitpad;  //?
	double y1, y2, y3, y4 ,y5, y6, y7, y8, y9, y10,y11, y12, w1, w2, w3;
	double s1, ss, sm, s2, ys, sig, ds, xd, sigs, sigg;
	double x1[4], x2[4], x3[4];
	double y[4];
	double sxi[12], dsi[3], ssi[3], sxim, dsxim;
	int ni, i, ind, ind1, ind2, j;
	//common/rims/ smin, smax;  //?
	//double smin, smax; //
	double ssj[5], dsj[5], ssjt, dsjt, der;
	vector<double> xi;
	vector<double> si;
	xi.push_back(1.); xi.push_back(2.); xi.push_back(3.); xi.push_back(4.); xi.push_back(5.); xi.push_back(6.); xi.push_back(7.); xi.push_back(8.); xi.push_back(9.);
	si.push_back(1.); si.push_back(2.); si.push_back(3.); si.push_back(4.);
	der=2.40385;
	ni=0;
	for (int i=0; i<3; i++){
		if (si[i] != -1000){
			ni=ni+1;
		}
	}
	ds=1.5/120.0*180./TMath::Pi();
	xd=x;
	y1=fabs(fitpad[0]);
	y2=fabs(fitpad[1]);
	y3=fabs(fitpad[2]);
	y4=fabs(fitpad[3]);
	y5=fabs(fitpad[4]);
	y6=fabs(fitpad[5]);
	y7=fabs(fitpad[6]);
	y8=fabs(fitpad[7]);
	y9=fabs(fitpad[8]);
	y10=fabs(fitpad[9]);
	y11=fabs(fitpad[10]);
	y12=fabs(fitpad[8]);
	s1=fabs(fitpad[11]);
	ss=fabs(fitpad[12]);
	sm=fabs(fitpad[13]);
	s2=fabs(fitpad[14]);
	ys=fabs(fitpad[15]);
	sig=fabs(fitpad[16]);
	sigs=fabs(fitpad[17]);
    sigg=fabs(fitpad[18]);
    sigg=sigs;
    dsi[0]=fabs(fitpad[19]/2.0);
    dsi[1]=fabs(fitpad[20]/2.0);
    dsi[2]=fabs(fitpad[21]/2.0);
    ssi[0]=si[0]+fitpad[22];
    ssi[1]=si[1]+fitpad[23];
    ssi[2]=si[2]+fitpad[24];
    w1=fabs(fitpad[25]);
    w2=fabs(fitpad[26]);
    w3=fabs(fitpad[27]);
    smin=s1;
    smax=s2;
    ind=0;
    if (((fabs(smin-0.) < 10.) || (fabs(smin-120.) < 10.) ||(fabs(smin-240.) < 10.) || (fabs(smin-360.) < 10.))){
		ind=ind+1;
		dsj[ind]=0.001;
		ssj[ind]=smin+der;
    }
    if (((fabs(smax-0.) < 10.) || (fabs(smax-120.) < 10.) || (fabs(smax-240.) < 10.) || (fabs(smax-360.) < 10.))){
		ind=ind+1;
		dsj[ind]=0.001;
		ssj[ind]=smax-der;
    }
    ni=ni+ind;
    for(int i=0; i<3; i++){
		ind=ind+1;
		dsj[ind]=dsi[i];
		ssj[ind]=ssi[i];
    }
    for(int i=0; i < ni-1; i++){
		for(int j=i+1; j < ni; j++){
			if (ssj[i] > ssj[j]){
			ssjt=ssj[i];
			dsjt=dsj[i];
			ssj[i]=ssj[j];
			dsj[i]=dsj[j];
			ssj[j]=ssjt;
			dsj[j]=dsjt;
			}
		}
    }
    fitfun=0.;
    ds=fabs(fitpad[20])/120.0*180./TMath::Pi();
	if (xd < sm) {
		x1[0] = s1;
		x1[1] = xi[1];
		x1[2] = xi[2];
		x1[3] = ss;
		y[0] = y1;
		y[1] = y2;
		y[2] = y3;
		y[3] = y4;
		m4(x1, y);
		fitfun = p3g(xd, s1, sig) - p3g(xd, ss - ds, sigs);
		
		ind = 0;
		ind1 = ind;
		ind2 = -1;
		sxi[ind1] = s1;
		for (int i = 0; i < ni; i++) {
			if ((ssj[i] - dsj[i] > s1) && (ssj[i] - dsj[i] < ss - ds)) {
				ind = ind + 1;
				sxi[ind] = ssj[i] - dsj[i];
			}
			if ((ssj[i] + dsj[i] > s1) && (ssj[i] + dsj[i] < ss - ds)) {
				ind = ind + 1;
				sxi[ind] = ssj[i] + dsj[i];
			}
			if ((s1 > ssj[i] - dsj[i]) && (s1 < ssj[i] + dsj[i])) {
				ind1 = ind;
			}
			if ((ss - ds > ssj[i] - dsj[i]) && (ss - ds < ssj[i] + dsj[i])) {
				ind2 = ind;
			}
		}
		if (ind2 == -1) {
			ind = ind + 1;
			ind2 = ind;
			sxi[ind2] = ss - ds - w1;
		}
		for (int i = ind1; i <= ind2; i = i+2) { //?
			if ((i == ind1) && (i + 1 == ind2)) {
				fitfun = fitfun + p3g(xd, sxi[i], sig) - p3g(xd, sxi[i + 1], sigs);
			}
			else if ((i == ind1) && (i + 1 != ind2)) {
				fitfun = fitfun + p3g(xd, sxi[i], sig) - p3g(xd, sxi[i + 1], sigg);
				sxim = (sxi[i + 2] + sxi[i + 1]) / 2.;
				dsxim = (sxi[i + 2] - sxi[i + 1]) / 2. / 3.;
				fitfun = fitfun + (p3g(xd, sxim, sig) - p3g(xd, sxim, sigs)) / 2. / 2. * (TMath::Erf(((xd - (sxim - dsxim)) / sqrt(2.0) / sigg)) - TMath::Erf(((xd - (sxim + dsxim)) / sqrt(2.0) / sigg)));
			}
			else if ((i != ind1) && (i + 1 == ind2)) {
				fitfun = fitfun + p3g(xd, sxi[i], sigg) - p3g(xd, sxi[i + 1], sigs);
				sxim = (sxi[i] + sxi[i - 1]) / 2.;
				dsxim = (sxi[i] - sxi[i - 1]) / 2. / 3.;
				fitfun = fitfun + (p3g(xd, sxim, sig) - p3g(xd, sxim, sigs)) / 2. / 2. * (TMath::Erf(((xd - (sxim - dsxim)) / sqrt(2.0) / sigg)) - TMath::Erf(((xd - (sxim + dsxim)) / sqrt(2.0) / sigg)));
			}
			else {
				fitfun = fitfun + p3g(xd, sxi[i], sigg) - p3g(xd, sxi[i + 1], sigg);
				sxim = (sxi[i] + sxi[i - 1]) / 2.;
				dsxim = (sxi[i] - sxi[i - 1]) / 2. / 3.;
				fitfun = fitfun + (p3g(xd, sxim, sig) - p3g(xd, sxim, sigs)) / 2. / 2. * (TMath::Erf(((xd - (sxim - dsxim)) / sqrt(2.0) / sigg)) - TMath::Erf(((xd - (sxim + dsxim)) / sqrt(2.0) / sigg)));
				sxim = (sxi[i + 2] + sxi[i + 1]) / 2.;
				dsxim = (sxi[i + 2] - sxi[i + 1]) / 2. / 3.;
				fitfun = fitfun + (p3g(xd, sxim, sig) - p3g(xd, sxim, sigs)) / 2. / 2. * (TMath::Erf(((xd - (sxim - dsxim)) / sqrt(2.0) / sigg)) - TMath::Erf(((xd - (sxim + dsxim)) / sqrt(2.0) / sigg)));
			}
		}
		x2[0] = ss;
		x2[1] = xi[4];
		x2[2] = xi[5];
		x2[3] = sm;
		y[0] = y5;
		y[1] = y6;
		y[2] = y7;
		y[3] = y8;
		m4(x2, y);
		fitfun = fitfun + p3g(xd, ss + ds, sigs) - p3g(xd, sm, sigs);
		
		ind = 0;
		ind1 = ind;
		ind2 = -1;
		sxi[ind1] = ss + ds + w2;
		for (int i = 0; i < ni; i++) {
			if ((ssj[i] - dsj[i] > ss + ds) && (ssj[i] - dsj[i] < sm)) {
				ind = ind + 1;
				sxi[ind] = ssj[i] - dsj[i];
			}
			if ((ssj[i] + dsj[i] > ss + ds) && (ssj[i] + dsj[i] < sm)) {
				ind = ind + 1;
				sxi[ind] = ssj[i] + dsj[i];
			}
			if ((ss + ds > ssj[i] - dsj[i]) && (ss + ds < ssj[i] + dsj[i])) {
				ind1 = ind;
			}
			if ((sm > ssj[i] - dsj[i]) && (sm < ssj[i] + dsj[i])) {
				ind2 = ind;
			}
		}
		if (ind2 == -1) {
			ind = ind + 1;
			ind2 = ind;
			sxi[ind2] = sm - w3 / 2.0;
		}
		for (int i = ind1; i <= ind2; i = i+2) {   //
			if (i == ind1 && i + 1 == ind2) {
				fitfun = fitfun + p3g(xd, sxi[i], sigs) - p3g(xd, sxi[i + 1], sigs);
			}
			else if (i == ind1 && i + 1 != ind2) {
				fitfun = fitfun + p3g(xd, sxi[i], sigs) - p3g(xd, sxi[i + 1], sigg);
			}
			else if (i != ind1 && i + 1 == ind2) {
				fitfun = fitfun + p3g(xd, sxi[i], sigg) - p3g(xd, sxi[i + 1], sigs);
			}
			else {
				fitfun = fitfun + p3g(xd, sxi[i], sigg) - p3g(xd, sxi[i + 1], sigg);
			}
		}

	}
	else{
		x3[0]=sm;
		x3[1]=xi[7];
		x3[2]=xi[8];
		x3[3]=s2;
		y[0]=y12;
		y[1]=y9;
		y[2]=y10;
		y[3]=y11;
		m4(x3,y);
		fitfun=fitfun+p3g(xd,sm,sigs)-p3g(xd,s2,sig); 
		fitfun=p3g(xd,s2,-sig);
		
		ind=0;
		ind1=ind;
		ind2=-1;
		sxi[ind1]=sm+w3/2.0;
		for (int i = 0; i < ni; i++) {
			if ((ssj[i] - dsj[i] > sm) && (ssj[i] - dsj[i] < s2)) {
				ind = ind + 1;
				sxi[ind] = ssj[i] - dsj[i];
			}
			if ((ssj[i] + dsj[i] > sm) && (ssj[i] + dsj[i] < s2)) {
				ind = ind + 1;
				sxi[ind] = ssj[i] + dsj[i];
			}
			if ((sm > ssj[i] - dsj[i]) && (sm < ssj[i] + dsj[i])) {
				ind1 = ind;
			}
			if ((s2 > ssj[i] - dsj[i]) && (s2 < ssj[i] + dsj[i])) {
				ind2 = ind;
			}
		}
		if (ind2 == -1){
			ind=ind+1;
			ind2=ind;
			sxi[ind2]=s2;
		}
		for(int i=ind1; i <= ind2; i=i+2){  //
			if ((i == ind1) && (i + 1 == ind2)) {
				fitfun = fitfun + p3g(xd, sxi[i], sigs) - p3g(xd, sxi[i + 1], sig);
			}
			else if (i == ind1 && i + 1 != ind2) {
				fitfun = fitfun + p3g(xd, sxi[i], sigs) - p3g(xd, sxi[i + 1], sigg);
			}
			else if (i != ind1 && i + 1 == ind2) {
				fitfun = fitfun + p3g(xd, sxi[i], sigg) - p3g(xd, sxi[i + 1], sig);
			}
			else {
				fitfun=fitfun+p3g(xd,sxi[i],sigg)-p3g(xd,sxi[i+1],sigg);
			}
		}
	}
	fitfun=fitfun+agt0(ys)/2.0*(TMath::Erf(((xd-(ss-ds))/sqrt(2.0)/sigs))-TMath::Erf(((xd-(ss+ds))/sqrt(2.0)/sigs)));
	fitfun=fitfun+agt0(0.0)-agt0(0.0)/2.0*(TMath::Erf(((xd-s1)/sqrt(2.0)/sig))-TMath::Erf(((xd-s2)/sqrt(2.0)/sig)));
	if (fitfun > 1.0e-5) {
		fitfun = fitfun / (1.0 - exp(-fitfun));
	}
	else {
		fitfun = 2.0/(2.0-fitfun);
	}
	if (fitfun > 1000){
		fitfun = 1.0;
	}
	scphi=fitfun;
	return scphi;
}


void go() {
	TCanvas* c1 = new TCanvas("c1", "A", 200, 10, 500, 300);
	TF1* f1 = new TF1("f1", "scphi(x)", 0, 10);
	f1->Draw();
}

//pscp D:\programs\funcfit.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/funcfit.cpp