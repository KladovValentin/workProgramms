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

#include "TColor.h"
#include "THistPainter.h"
#include "TFile.h"
#include "TCollection.h"
#include "TKey.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TMinuit.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TArrow.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TVector3.h"
#include "TTree.h"
#include "TLegend.h"
#include "TChain.h"
#include "TGraph.h"

#define pi TMath::Pi()
//#define scal (pi/180)

//-----------------------------------------------------------------------------
// /work/users/konctbel/snd2k/R006-003/.mainrelease/TemDatabase/addon-snd-cps.sql
// mysql -h sndfs1.sndonline -P 3307 -u daqread -pdaqread temdbase -e "select * from v_dir" | grep 192
// mysql -h sndfs1.sndonline -P 3307 -u daqread -pdaqread temdbase -e "select r.c_run, d.* from t_run r, t_datadouble d where r.c_run=39808 and d.c_stamp between r.c_stamp0 and r.c_stamp1 and d.c_channel=385"
// mysql -h sndfs1.sndonline -P 3307 -u daqread -pdaqread temdbase -e "select r.c_run, sum(d.c_value) from t_run r, t_datadouble d where r.c_run=39808 and d.c_stamp between r.c_stamp0 and r.c_stamp1 and d.c_channel=385"
// mysql -h sndfs1.sndonline -P 3307 -u daqread -pdaqread temdbase -e "select r.c_run, avg(d.c_value) from t_run r, t_datadouble d where r.c_run=39808 and d.c_stamp between r.c_stamp0 and r.c_stamp1 and d.c_channel=385"

TCanvas * stop()
{
   TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
   c->Update();
   std::cout << "Push any key!" << std::endl;
   char key;
   std::cin>>key;
   return c;
}

Double_t fDC1(size_t run)
{
   TString cmd;
   cmd = Form("mysql -h sndfarm10.sndonline -P3307 -udaqread -pdaqread temdbase -e \"select i.c_value / d.c_value from t_run r inner join t_datainteger0 i on i.c_channel = 47 and i.c_stamp = r.c_stamp1 inner join t_datadouble0 d on d.c_channel =4 and d.c_stamp = r.c_stamp1 where r.c_run= %d\"",run);
   cmd = cmd + " | " + Form("sed '2q;d'");
   return std::atof(gSystem->GetFromPipe(cmd.Data()));
}

Double_t fbeam(size_t run)
{
   TString cmd;
   //cmd = Form("mysql -h sndfarm11.sndonline -u daqread -pdaqread rundbase -e \"select c_run,c_energy from t_readout where c_run=%d\"",run);
   cmd = Form("cd /work/users/konctbel/snd2k/R999-999clb; source /work/snd2000/root/setup2k.sh; sndsql.py rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"",run);
   cmd = cmd + " | " + Form("fgrep -e %d",run);
   cmd = cmd + " | " + Form("sed 's/%d//g'",run);
   return std::atof(gSystem->GetFromPipe(cmd.Data()));
}

size_t fdate(size_t run)
{
   TString cmd;
   //cmd = Form("mysql -h sndfarm11.sndonline -u daqread -pdaqread rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"",run);
   cmd = Form("cd /work/users/konctbel/snd2k/R999-999clb; source /work/snd2000/root/setup2k.sh; sndsql.py rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"",run);
   cmd = cmd + " | " + Form("fgrep -e %d",run);
   cmd = cmd + " | " + Form("sed 's/%d//g'",run);
   cmd = Form("date -d \"$(%s)\" +\"%s%s%s\"",cmd.Data(),"%Y","%m","%d");
   return std::atoi(gSystem->GetFromPipe(cmd.Data()));
}

TString fdatetime(size_t run)
{
   TString cmd;
   //cmd = Form("mysql -h sndfarm11.sndonline -u daqread -pdaqread rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"",run);
   cmd = Form("cd /work/users/konctbel/snd2k/R999-999clb; source /work/snd2000/root/setup2k.sh; sndsql.py rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"",run);
   cmd = cmd + " | " + Form("fgrep -e %d",run);
   cmd = cmd + " | " + Form("sed 's/%d//g'",run);
   cmd = Form("date -d \"$(%s)\" +\"%s-%s-%s %s:%s:%s\"",cmd.Data(),"%Y","%m","%d","%H","%M","%S");
   return gSystem->GetFromPipe(cmd.Data());
}

size_t fdatetotime(size_t date)
{
   TString cmd;
   cmd = Form("date -d \"%d\" +\"%s\"",date,"%s");
   return std::atoi(gSystem->GetFromPipe(cmd.Data())) - gStyle->GetTimeOffset() - 6*3600;
}

size_t ftime(size_t run)
{
   TString cmd;
   //cmd = Form("mysql -h sndfarm11.sndonline -u daqread -pdaqread rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"",run);
   cmd = Form("cd /work/users/konctbel/snd2k/R999-999clb; source /work/snd2000/root/setup2k.sh; sndsql.py rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"",run);
   cmd = cmd + " | " + Form("fgrep -e %d",run);
   cmd = cmd + " | " + Form("sed 's/%d//g'",run);
   cmd = Form("date -d \"$(%s)\" +\"%s\"",cmd.Data(),"%s");
   return std::atoi(gSystem->GetFromPipe(cmd.Data())) - gStyle->GetTimeOffset() - 6*3600;
}

Double_t fpar(size_t run, size_t par)
{
   TString cmd;
   cmd = "mysql -h sndfs1.sndonline -P 3307 -u daqread -pdaqread temdbase -e \"select * from v_dir\"";
   cmd = cmd + " | " + Form("grep %d",par);
   size_t ind = std::atoi(gSystem->GetFromPipe(cmd.Data()));
   cmd = Form("mysql -h sndfs1.sndonline -P 3307 -u daqread -pdaqread temdbase -e \"select r.c_run, avg(d.c_value) from t_run r, t_datadouble d where r.c_run=%d and d.c_stamp between r.c_stamp0 and r.c_stamp1 and d.c_channel=%d\"",run,ind);
   cmd = cmd + " | " + Form("fgrep -e %d",run);
   cmd = cmd + " | " + Form("sed 's/%d//g'",run);
   return std::atof(gSystem->GetFromPipe(cmd.Data()));
}

bool myOrder(std::pair<Int_t,TString> a, std::pair<Int_t,TString> b) 
{ 
   return a.first < b.first; 
} 

std::vector< std::pair<Int_t,TString> > getCalTime(TString calName = "dcpagen", bool read = false)
{
   TString cmd;
   cmd = "cd /work/users/konctbel/calibs/R007-001/; source /work/snd2000/root/setup2k.sh; clbixlist " + calName + " CURRENT > " + calName + ".list";
   if (read) gSystem->Exec(cmd.Data());
   cmd = "cat /work/users/konctbel/calibs/R007-001/" + calName + ".list | sed '/^-/ d' | sed 's#^.*\\([0-9]\\{4\\}-[0-9]\\{2\\}-[0-9]\\{2\\} [0-9]\\{2\\}:[0-9]\\{2\\}:[0-9]\\{2\\}\\).*$#\\1#'";
   TString s = gSystem->GetFromPipe(cmd.Data());
   TObjArray *ts = s.Tokenize("\n");
   std::vector< std::pair<Int_t,TString> > res;
   for (Int_t i = 2; i < ts->GetEntries(); i++) {
      TString str = ((TObjString *)(ts->At(i)))->String();
      cmd = Form("date -d \"%s\" +\"%s\"",str.Data(),"%s");
      Int_t tm = (Int_t)(std::atoi(gSystem->GetFromPipe(cmd.Data())) - gStyle->GetTimeOffset() - 6*3600);
      res.push_back(std::make_pair(tm,str));
   }
   std::sort(res.begin(),res.end(),myOrder);
   return res;
}

bool myOrderInt(std::pair<Int_t,Int_t> a, std::pair<Int_t,Int_t> b) 
{ 
   return a.first < b.first; 
} 

std::vector< std::pair<Int_t,Int_t> > getCalRun(TString calName = "dcpagen", bool tmrn = true, bool read = false, TString usage = "CURRENT")
{
   TString fName = calName + "_" + usage + ".list";
   TString cmd;
   // reading cals list from db
   cmd = "cd /work/users/konctbel/calibs/R007-001/; source /work/snd2000/root/setup2k.sh; clbixlist " + calName + " " + usage + " > " + fName;
   if (read) gSystem->Exec(cmd.Data());
   // start run
   cmd = "cat /work/users/konctbel/calibs/R007-001/" + fName + " | sed '/^-/ d' | awk '{print $1}'";
   TString r1 = gSystem->GetFromPipe(cmd.Data());
   TObjArray *tr1 = r1.Tokenize("\n");
   // stop run
   cmd = "cat /work/users/konctbel/calibs/R007-001/" + fName + " | sed '/^-/ d' | awk '{print $2}'";
   TString r2 = gSystem->GetFromPipe(cmd.Data());
   TObjArray *tr2 = r2.Tokenize("\n");
   // save to db time
   //cmd = "cat /work/users/konctbel/calibs/R007-001/" + fName + " | sed '/^-/ d' | awk '{print $8\" \"$9}'";
   //TString s = gSystem->GetFromPipe(cmd.Data());
   //TObjArray *ts = s.Tokenize("\n");
   std::vector< std::pair<Int_t,Int_t> > res;
   for (Int_t i = 2; i < tr1->GetEntries(); i++) {
      //TString str = ((TObjString *)(ts->At(i)))->String();
      //cmd = Form("date -d \"%s\" +\"%s\"",str.Data(),"%s");
      //Int_t tm = (Int_t)(std::atoi(gSystem->GetFromPipe(cmd.Data())) - gStyle->GetTimeOffset() - 6*3600);
      Int_t run1 = (((TObjString *)(tr1->At(i)))->String()).Atoi();
      Int_t run2 = (((TObjString *)(tr2->At(i)))->String()).Atoi();
      if (tmrn) res.push_back(std::make_pair(ftime(run1),ftime(run2)));
      else res.push_back(std::make_pair(run1,run2));
   }
   std::sort(res.begin(),res.end(),myOrderInt);
   return res;
}

class DcPaGen
{
 public:
   DcPaGen(){};
   size_t run;
   size_t t;
   std::vector<Double_t> al;
   std::vector<Double_t> ar;
   std::vector<Double_t> bl;
   std::vector<Double_t> br;
   std::vector<Double_t> gl;
   std::vector<Double_t> gr;
   std::vector<Double_t> chi2l;
   std::vector<Double_t> chi2r;
};

class PcGen
{
 public:
   PcGen(){};
   size_t run;
   size_t t;
   std::vector<Double_t> a;
   std::vector<Double_t> b;
   std::vector<Double_t> chi2;
};

class DcZabs
{
 public:
   DcZabs(){};
   size_t run;
   size_t t;
   std::vector<Double_t> a0;
   std::vector<Double_t> a0err;
   std::vector<Double_t> a1;
   std::vector<Double_t> a1err;
   std::vector<Double_t> c01;
   std::vector<Double_t> chi2;
   std::vector<Double_t> ndim;
   std::vector<Double_t> z1;
   std::vector<Double_t> z2;
};

std::vector<Double_t> getCalArray(int run, TString cal, TString array, TString usage = "CURRENT")
{
   //TString cal = "dcpagen";
   TString dir = "/work/users/konctbel/calibs/R007-001/";
   TString cmd;
   TString fname = Form("%s%s",dir.Data(),cal.Data());
   if (gSystem->AccessPathName(fname)) {
      cmd = "mkdir " + fname;
      gSystem->Exec(cmd.Data());
   }
   fname = Form("%s%s/%s",dir.Data(),cal.Data(),usage.Data());
   if (gSystem->AccessPathName(fname)) {
      cmd = "mkdir " + fname;
      gSystem->Exec(cmd.Data());
   }
   fname = Form("%s%s/%s/%s%d.txt",dir.Data(),cal.Data(),usage.Data(),array.Data(),run);
   if (gSystem->AccessPathName(fname)) {
      cmd = "cd "+dir+"; source /work/snd2000/root/setup2k.sh; clbarray " + cal + " " + array + " " + usage + " " + Form("%d",run) +" > " + fname;
      gSystem->Exec(cmd.Data());
      std::cout << "Calibration array " << run << ":" << cal << ":" << array << "was stored on disk" << std::endl;
   }
   cmd = "cat " + fname + " | grep -e \"data:\" | awk '{print $2}'";
   size_t nvar = std::atoi(gSystem->GetFromPipe(cmd.Data()));
   cmd = "cat " + fname + " | grep '[0-9][0-9][0-9][0-9]:' | sed 's/\\(^.*:\\)/ /g' | sed 's/ /\\n/g' | sed '/^$/d'";
   TString s = gSystem->GetFromPipe(cmd.Data());
   TObjArray *ts = s.Tokenize("\n");
   std::vector<Double_t> data(nvar,0);
   for (Int_t i = 0; i < ts->GetEntries(); i++)
      data[i] = std::atof(((TObjString *)(ts->At(i)))->String());
   return data;
}

std::vector<DcPaGen> getDcPaGen(Int_t tStart, Int_t tStop, bool tm = true)
{
   std::vector< std::pair<Int_t,Int_t> > res = getCalRun("dcpagen",tm);
   std::vector<DcPaGen> cal;
   for(size_t i=0;i<res.size();i++) {
      if (res[i].first<tStart||res[i].first>tStop) continue;
      cal.push_back(DcPaGen());
      cal.back().run = res[i].second;
      cal.back().t = ftime(res[i].second);
      cal.back().al = getCalArray(res[i].second,"dcpagen","al");
      cal.back().ar = getCalArray(res[i].second,"dcpagen","ar");
      cal.back().bl = getCalArray(res[i].second,"dcpagen","bl");
      cal.back().br = getCalArray(res[i].second,"dcpagen","br");
      cal.back().gl = getCalArray(res[i].second,"dcpagen","gl");
      cal.back().gr = getCalArray(res[i].second,"dcpagen","gr");
      cal.back().chi2l = getCalArray(res[i].second,"dcpagen","chi2l");
      cal.back().chi2r = getCalArray(res[i].second,"dcpagen","chi2r");
   }
   return cal;
}

std::vector<PcGen> getPcGen(Int_t tStart, Int_t tStop, bool tm = true)
{
   std::vector< std::pair<Int_t,Int_t> > res = getCalRun("pcgen",tm);
   std::vector<PcGen> cal;
   for(size_t i=0;i<res.size();i++) {
      if (res[i].first<tStart||res[i].first>tStop) continue;
      cal.push_back(PcGen());
      cal.back().run = res[i].second;
      cal.back().t = ftime(res[i].second);
      cal.back().a = getCalArray(res[i].second,"pcgen","a");
      cal.back().b = getCalArray(res[i].second,"pcgen","b");
      cal.back().chi2 = getCalArray(res[i].second,"pcgen","chi2");
   }
   return cal;
}

std::vector<DcZabs> getDcZabs(Int_t tStart, Int_t tStop, bool tm = true, TString usage = "RECONSTR")
{
   TString name = "dczabs";
   std::vector< std::pair<Int_t,Int_t> > res = getCalRun(name,tm,true,usage);
   std::vector<DcZabs> cal;
   for(size_t i=0;i<res.size();i++) {
      if (res[i].first<tStart||res[i].first>tStop) continue;
      cal.push_back(DcZabs());
      cal.back().run = res[i].second;
      cal.back().t = ftime(res[i].second);
      cal.back().a0 = getCalArray(res[i].second,name,"a0",usage);
      cal.back().a0err = getCalArray(res[i].second,name,"a0err",usage);
      cal.back().a1 = getCalArray(res[i].second,name,"a1",usage);
      cal.back().a1err = getCalArray(res[i].second,name,"a1err",usage);
      cal.back().c01 = getCalArray(res[i].second,name,"c01",usage);
      cal.back().chi2 = getCalArray(res[i].second,name,"chi2",usage);
      cal.back().ndim = getCalArray(res[i].second,name,"ndim",usage);
      cal.back().z1 = getCalArray(res[i].second,name,"z1",usage);
      cal.back().z2 = getCalArray(res[i].second,name,"z2",usage);
   }
   return cal;
}

void xxx()
{
   TString calName = "dcpagen";
   std::vector< std::pair<Int_t,TString> > res = getCalTime(calName,true);
   std::vector< std::pair<Int_t,Int_t> > res1 = getCalRun(calName);
   for (size_t i=0;i<res.size();i++)
      std::cout << i << " " << res[i].first << " " << res[i].second << " " << res1[i].second << std::endl;   
   //TGraph *gr = new TGraph(res.size(),res.data(),res.data());
   //gr->GetXaxis()->SetTimeDisplay(1);
   //gr->SetMarkerSize(0.5);
   //gr->SetMarkerStyle(20);
   //gr->Draw("AP");
}

void yyy()
{
   std::vector<Int_t> res;
   for (Int_t i = 0; i < 24; i++) {
      TString cmd = Form("date -d \"2020-02-02 %02d:30:00\" +\"%s\"",i,"%s");
      std::cout << cmd << std::endl;
      res.push_back(std::atoi(gSystem->GetFromPipe(cmd.Data())) - gStyle->GetTimeOffset() - 6*3600);
   }
   TGraph *gr = new TGraph(res.size(),res.data(),res.data());
   gr->GetXaxis()->SetTimeDisplay(1);
   gr->SetMarkerSize(0.5);
   gr->SetMarkerStyle(20);
   gr->Draw("AP");
}

//-----------------------------------------------------------------------------

size_t npars = 0;
double fcn_min = 0;
size_t nmax;
size_t nps = 0;
std::vector<size_t> ls; 
std::vector<size_t> im;
std::vector<double> t;
std::vector<double> dt;
std::vector<double> theta;
std::vector< std::vector<double> > ne(10,std::vector<double>());
std::vector< std::vector<double> > sz(10,std::vector<double>());
std::vector< std::vector<double> > dsz(10,std::vector<double>());
std::vector< std::vector<double> > dedx0(10,std::vector<double>());
std::vector< std::vector<double> > ddedx0(10,std::vector<double>());
std::vector< std::vector<double> > dedxm(10,std::vector<double>());
std::vector< std::vector<double> > ddedxm(10,std::vector<double>());
std::vector< std::vector<double> > slope(10,std::vector<double>());
std::vector< std::vector<double> > dslope(10,std::vector<double>());
std::vector< std::vector<double> > scale(10,std::vector<double>());
std::vector< std::vector<double> > dscale(10,std::vector<double>());
std::vector< std::vector< std::vector<double> > > dedx(10,std::vector< std::vector<double> >());
std::vector< std::vector< std::vector<double> > > ddedx(10,std::vector< std::vector<double> >());

//______________________________________________________________________________

double dei(double z, double r, double theta, double sz)
{
   double cth = fabs(cos(theta));
   double sth = sin(theta);
   double zm = r*cth/sth/2;
   if (10*zm>sz) 
      return (erf((z+zm)/sqrt(2.)/sz)-erf((z-zm)/sqrt(2.)/sz))/2/cth;
   else
      return exp(-(z/sz)*(z/sz)/2)/sqrt(2*3.1415927)/sz*r/sth;
}

double satfsx( double dedx, double theta, 
               double sz, double demax,
               double &g0, double &g1, double &g2 )
{
   double r, zm, cth, sth;
   
   //sz = std::fabs(sz);
   cth = fabs(cos(theta));
   sth = sin(theta);
   
   r = 0.4*2;
   zm = r*cth/sth/2;
   
   double zmm = zm-3*sz;
   double zmp = zmm>0 ? 6*sz : zm+3*sz;
   size_t n1 = 30;
   double d1 = zmp/n1;
   size_t n2 = zmm>0 ? (size_t)std::min(20.,n1*zmm/zmp+1) : 0;
   double d2 = zmm>0 ? zmm/n2 : 0;
   double zi = -zm-3*sz;
   double ndei = dei(zi,r,theta,sz)/demax; 
   double f0l = demax*(1-exp(-dedx*ndei));
   double f1l = demax*(exp(-dedx*ndei)*ndei);
   double f2l = demax*(-exp(-dedx*ndei)*ndei*ndei/2);
   g0 = 0;
   g1 = 0;
   g2 = 0;
   for(size_t i=0;i<n1+n2;i++) {
      double d = i<n1 ? d1 : d2;
      zi += d;
      ndei = dei(zi,r,theta,sz)/demax;
      double f0r = demax*(1-exp(-dedx*ndei));
      double f1r = demax*(exp(-dedx*ndei)*ndei);
      double f2r = demax*(-exp(-dedx*ndei)*ndei*ndei/2);
      g0 += (f0r+f0l)*d/2;
      g1 += (f1r+f1l)*d/2;
      g2 += (f2r+f2l)*d/2;
      f0l = f0r;
      f1l = f1r;
      f2l = f2r;
   }
   g0 *= 2*sth/r;
   g1 *= 2*sth/r;
   g2 *= 2*sth/r;

   return g0;
}

//______________________________________________________________________________

class fParN
{
 public:
   fParN(size_t nm_, size_t xmin_, size_t xmax_, TString name_ = ""){nm = nm_; xmin = xmin_; xmax = xmax_; name = name_;};
   double Eval(size_t t) {return Eval((double)t,p);};
   double Eval(double x, std::vector<double> y);
   double Eval(size_t t, Double_t *par);
   bool IsInTime(size_t t) {if (t>=xmin&&t<=xmax) return true; else return false;};
   void init(size_t ind_, std::vector<double> t, std::vector<double> q, std::vector<double> dq);
   void plot(Double_t low=0, Double_t up=0);
   void function(double *par);
   
   TString name;
   TF1* pn;
   size_t nm;
   size_t ind;
   double xmin;
   double xmax;
   std::vector<Double_t> p;
   std::vector<Double_t> dp;
   std::vector<Double_t> xp;
   std::vector<Double_t> dxp;
   std::vector<Double_t> yp;
   std::vector<Double_t> dyp;
   std::vector<Double_t> xt;
   std::vector<Double_t> dxt;
   std::vector<Double_t> yt;
   std::vector<Double_t> dyt;
};

double fParN::Eval(size_t t, Double_t *par)
{
   if (!nm) return par[ind];
   std::vector<double> y(nm+1);
   for(size_t i=0;i<=nm;i++) y[i] = par[ind+i];
   return Eval((double)t,y);
}

double fParN::Eval(double x, std::vector<double> y)
{
   size_t n = y.size();
   if (!n) return 0;
   if (n==1) return y[0];
   //
   std::vector<double> xi(n,0);
   for(size_t i=0;i<n;i++) xi[i] = (double)i/(n-1);
   //
   double db[n];
   double da[n*n];
   for(size_t i=0;i<n;i++) {
      db[i] = y[i];
      for(size_t j=0;j<n;j++)
         da[n*i+j] = std::pow(xi[i],(double)j);
   }
   //
   TVectorD B(n);
   B.SetElements(db);
   TMatrixD A(n,n);
   A.SetMatrixArray(da);
   //
   TDecompLU lu(A);
   Bool_t ok;
   TVectorD px = lu.Solve(B,ok);
   double f = px(0);
   double xn = (x-xmin)/(xmax-xmin);
   for(size_t i=1;i<n;i++)
      f += px(i)*std::pow(xn,(double)i);
   return f;
}

void fParN::init(size_t ind_, std::vector<double> t, std::vector<double> q, std::vector<double> dq)
{
   ind = ind_;
   xp.resize(0);
   yp.resize(0);
   dxp.resize(0);
   dyp.resize(0);
   std::vector<Double_t> xc;
   std::vector<double> dqs(dq.size());
   std::copy(dq.begin(),dq.end(),dqs.begin());
   std::sort(dqs.begin(),dqs.end());
   double dqmin = dqs[dq.size()*0.05];
   double dqmax = dqs[dq.size()*0.95];
   for(size_t i=0;i<t.size();i++) {
      if (TMath::IsNaN(q[i])||TMath::IsNaN(dq[i])) continue;
      if (dq[i]<dqmin||dq[i]>dqmax) continue;
      xt.push_back(t[i]);
      dxt.push_back(0);
      yt.push_back(q[i]);
      dyt.push_back(dq[i]);
      if (t[i]<xmin||t[i]>xmax) continue;
      xc.push_back((t[i]-xmin)/(xmax-xmin));
      xp.push_back(t[i]);
      dxp.push_back(0);
      yp.push_back(q[i]);
      dyp.push_back(dq[i]);
   }
   TGraphErrors *gr = new TGraphErrors(xc.size(),xc.data(),yp.data(),dxp.data(),dyp.data());
   double ns = 100;
   gr->Fit("pol0");
   pn = (TF1*)gr->GetFunction("pol0");
   for(size_t s=0;s<3;s++) {
      std::vector<double> xs;
      std::vector<double> ys;
      std::vector<double> dxs;
      std::vector<double> dys;
      for(size_t i=0;i<xc.size();i++) {
         if (std::fabs(yp[i]-pn->Eval(xc[i]))>ns*dyp[i]) continue;
         xs.push_back(xc[i]);
         dxs.push_back(dxp[i]);
         ys.push_back(yp[i]);
         dys.push_back(dyp[i]);
      }
      TGraphErrors *grf = new TGraphErrors(xs.size(),xs.data(),ys.data(),dxs.data(),dys.data());
      grf->Fit(Form("pol%d",nm));
      pn = (TF1*)grf->GetFunction(Form("pol%d",nm));
      ns = 3;
   }
   p.resize(0);
   if (!nm) {p.push_back(pn->Eval(xmin)); return;}
   for(double xi=0;xi<=1;xi+=1.0/nm)
      p.push_back(pn->Eval(xi));
}

void fParN::plot(Double_t low, Double_t up)
{
   TGraphErrors *gr0 = new TGraphErrors(xt.size(),xt.data(),yt.data(),dxt.data(),dyt.data());
   gr0->GetXaxis()->SetTimeDisplay(1);
   gr0->SetMarkerSize(0.5);
   gr0->SetMarkerStyle(20);
   gr0->Draw("AP");
   gr0->SetName(name);
   if (low!=up) gr0->GetYaxis()->SetRangeUser(std::max(low,gr0->GetYaxis()->GetXmin()),
                                              std::min(up,gr0->GetYaxis()->GetXmax()));
   //
   TGraphErrors *gr = new TGraphErrors(xp.size(),xp.data(),yp.data(),dxp.data(),dyp.data());
   gr->GetXaxis()->SetTimeDisplay(1);
   gr->SetMarkerSize(0.5);
   gr->SetMarkerStyle(20);
   gr->Draw("P");
   //
   std::vector<double> f;
   for(size_t i=0;i<xp.size();i++)
      f.push_back(Eval(xp[i],p));
   TGraphErrors *grf = new TGraphErrors(xp.size(),xp.data(),f.data(),dxp.data(),dxp.data());
   grf->GetXaxis()->SetTimeDisplay(1);
   grf->SetMarkerSize(0.5);
   grf->SetMarkerStyle(20);
   grf->SetLineColor(kRed);
   grf->SetMarkerColor(kRed);
   grf->Draw("PL");
   TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
   c->SetTitle(name);
   c->SetName(name);
   c->Update();
}

void fParN::function(double *par)
{
   std::vector<double> y(nm+1);
   for(size_t i=0;i<=nm;i++) y[i] = par[ind+i];
   std::vector<double> xf;
   std::vector<double> yf;
   for(double x=xmin;x<=xmax;x+=(xmax-xmin)/100) {
      xf.push_back(x);
      //yf.push_back(Eval(x,y));
      yf.push_back(Eval((size_t)x,par));
   }
   TGraph *gr = new TGraph(xf.size(),xf.data(),yf.data());
   gr->SetMarkerSize(0.5);
   gr->SetMarkerStyle(20);
   gr->SetLineColor(kRed);
   gr->SetMarkerColor(kRed);
   gr->Draw("PL");
}

//______________________________________________________________________________

std::vector<fParN> f_sz;
std::vector<fParN> f_alpha;
std::vector<fParN> f_slope;

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   //calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   npars = 0;
   for (size_t l=0;l<ls.size();l++ ) {
      size_t layer = ls[l];
      //for (size_t r=0;r<dedx[layer].size();r++) {
      for (size_t r=0;r<nmax;r++) {
         size_t ind = 100;
         for(size_t q=0;q<f_sz.size();q++) if (f_sz[q].IsInTime(t[r])) ind = q;
         if (ind==100) std::cout << "ind = " << ind << std::endl;
         double szi = f_sz[ind].Eval(t[r],par);
         //std::cout << ind << ": " << szi << " " << f_sz[ind].xmin << " " << t[r] << " " << f_sz[ind].xmax << std::endl;;
         ind = 100;
         for(size_t q=0;q<f_alpha.size();q++) if (f_alpha[q].IsInTime(t[r])) ind = q;
         if (ind==100) std::cout << "ind = " << ind << std::endl;
         double alpha = f_alpha[ind].Eval(t[r],par);
         double demaxi = alpha/szi;
         //std::cout << "\t" << ind << ": " << alpha << " " << demaxi << " ";
         ind = 100;
         for(size_t q=0;q<f_slope.size();q++) if (f_slope[q].IsInTime(t[r])) ind = q;
         if (ind==100) std::cout << "ind = " << ind << std::endl;
         double slope = f_slope[ind].Eval(t[r],par);
         //std::cout << ind << ": " << slope << std::endl;;
         double g0, g1, g2;
         size_t i1 = 12-(im[layer]-1);
         size_t i2 = 12+(im[layer]-1);
         for (size_t i=i1;i<=i2; i++) {
            if (!ddedx[layer][r][i]) continue;
            if (ddedx[layer][r][i]<dedx[layer][r][i]/100) continue;
            double dedxi = par[nps+r] + slope*(theta[i]-TMath::Pi()/2);
            double dedxs = satfsx(dedxi, theta[i], szi, demaxi, g0, g1, g2 );
            delta  = (dedx[layer][r][i]-dedxs)/ddedx[layer][r][i];
            chisq += delta*delta;
            npars++;
         }
      }
   }
   f = chisq;
   fcn_min = f;
}

void fitall(size_t flayer = 0)
{
   TFile *f = new TFile("mhad2017y.root");
   TIter next_run(f->GetListOfKeys());
   TKey *key_run;
   TString key_run_name("TDirectoryFile");
   TString key_name("TFitResult");
   TString key_name_fun("f_fits_");
   size_t ts = fdatetotime(20170320);
   while ((key_run=(TKey*)next_run())) {
      if (!key_run_name.BeginsWith(key_run->GetClassName())) continue;
      TDirectory *d = (TDirectory*)f->Get(key_run->GetName());
      TIter next(d->GetListOfKeys());
      TKey *key;
      while ((key=(TKey*)next())) {
         TString name(key->GetName());
         //if (name.BeginsWith("TFitResult-dcpagen")&&
         if (name.BeginsWith("TFitResult-pcgen")&&
             name.Contains(key_name_fun)) {
            int run, layer;
            //sscanf(key->GetName(),"TFitResult-dcpagen%d_layer%d",&run,&layer);
            sscanf(key->GetName(),"TFitResult-pcgen%d_layer%d",&run,&layer);
            TFitResult *r = (TFitResult*)d->Get(key->GetName());
            std::cout << run << " " << layer << " " << r->Chi2()/r->Ndf() << 
               " " << r->Parameters().at(0) << " " << r->Errors().at(0) << 
               " " << r->Parameters().at(1) << " " << r->Errors().at(1) <<
               " " << r->Parameters().at(2) << " " << r->Errors().at(2) << std::endl;
            if (//run>ts||
                TMath::IsNaN(r->Parameters().at(0))||
                TMath::IsNaN(r->Parameters().at(1))||
                TMath::IsNaN(r->Parameters().at(2))||
                //r->Chi2()/r->Ndf()>100||
                //r->Parameters().at(0)<0.021||
                //r->Parameters().at(0)>0.199||
                r->Errors().at(0)<0.001) continue;
            //if (layer==flayer) t.push_back(ftime(run));
            if (layer==flayer) t.push_back(run);
            sz[layer].push_back(r->Parameters().at(0));
            dsz[layer].push_back(r->Errors().at(0));
            dedx0[layer].push_back(r->Parameters().at(1));
            ddedx0[layer].push_back(r->Errors().at(1));
            dedxm[layer].push_back(r->Parameters().at(2));
            ddedxm[layer].push_back(r->Errors().at(2));
            slope[layer].push_back(r->Parameters().at(3));
            dslope[layer].push_back(r->Errors().at(3));
            //
            //TProfile *h = (TProfile*)d->Get(Form("dcpagen%08d_layer%d",run,layer));
            TProfile *h = (TProfile*)d->Get(Form("pcgen%08d_layer%d",run,layer));
            ne[layer].push_back(h->GetEntries());
            //
            std::vector<double> dedxi;
            std::vector<double> ddedxi;
            for(size_t i=0;i<h->GetNbinsX();i++) {
               dedxi.push_back(h->GetBinContent(i));
               ddedxi.push_back(h->GetBinError(i));
            }
            dedx[layer].push_back(dedxi);
            ddedx[layer].push_back(ddedxi);
            //
            if (!theta.size())
               for(size_t i=0;i<h->GetNbinsX();i++)
                  theta.push_back(h->GetBinCenter(i));
         }
      }
   }
   //---------------------------------------
   size_t tmin = fdatetotime(20180101);
   size_t tmax = fdatetotime(20170101);
   for(size_t i=0;i<t.size();i++) {
      if (t[i]<tmin) tmin=t[i];
      if (t[i]>tmax) tmax=t[i];
   }
   //---------------------------------------
   
   nmax = dedx0[flayer].size();
   size_t nrun = dedx0[flayer].size();
   size_t npar = 0;
   //
   size_t xxxx;
   f_sz.push_back(fParN(3,tmin,fdatetotime(20170320)));
   f_sz.push_back(fParN(5,fdatetotime(20170320),tmax));
   for(size_t i=0;i<f_sz.size();i++) {
      f_sz[i].init(npar,t,sz[flayer],dsz[flayer]);
      npar += f_sz[i].nm+1;
      f_sz[i].plot();
      cin>>xxxx;
   }
   //
   f_alpha.push_back(fParN(6,tmin,tmax));
   for(size_t i=0;i<f_alpha.size();i++) {
      f_alpha[i].init(npar,t,dedxm[flayer],ddedxm[flayer]);
      npar += f_alpha[i].nm+1;
      f_alpha[i].plot();
      cin>>xxxx;
   }
   //
   f_slope.push_back(fParN(6,tmin,tmax));
   for(size_t i=0;i<f_slope.size();i++) {
      f_slope[i].init(npar,t,slope[flayer],dslope[flayer]);
      npar += f_slope[i].nm+1;
      f_slope[i].plot();
      cin>>xxxx;
   }
   npar += nrun;
   //
   TMinuit *gMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of 5 params
   gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

   // Set starting values and step sizes for parameters
   nps = 0;
   for(size_t i=0;i<f_sz.size();i++)
      for(size_t k=0;k<f_sz[i].p.size();k++)
         gMinuit->mnparm(nps++, Form("sz%d_p%d",i,k), f_sz[i].p[k], 0.1*f_sz[i].p[k], 
                         0.5*f_sz[i].p[k], 1.5*f_sz[i].p[k], ierflg);
   for(size_t i=0;i<f_alpha.size();i++)
      for(size_t k=0;k<f_alpha[i].p.size();k++)
         gMinuit->mnparm(nps++, Form("alpha%d_p%d",i,k), f_alpha[i].p[k], 0.1*f_alpha[i].p[k], 
                         0.5*f_alpha[i].p[k], 1.5*f_alpha[i].p[k], ierflg);
   for(size_t i=0;i<f_slope.size();i++)
      for(size_t k=0;k<f_slope[i].p.size();k++)
         gMinuit->mnparm(nps++, Form("slope%d_p%d",i,k), f_slope[i].p[k], 0.1*f_slope[i].p[k], 
                         0.5*f_slope[i].p[k], 1.5*f_slope[i].p[k], ierflg);
   for(size_t i=0;i<nrun;i++)
      gMinuit->mnparm(nps+i, Form("dedx0_%d",i), dedx0[flayer][i], ddedx0[flayer][i], 500,1000,ierflg);
   
   ls.resize(0);
   ls.push_back(flayer);
   im.resize(0);
   //im.push_back(7);
   //im.push_back(9);
   //im.push_back(9);
   //im.push_back(9);
   //im.push_back(8);
   //im.push_back(7);
   //im.push_back(7);
   //im.push_back(6);
   //im.push_back(6);
   //im.push_back(5);
   
   im.push_back(8);
   im.push_back(10);
   im.push_back(10);
   im.push_back(10);
   im.push_back(9);
   im.push_back(9);
   im.push_back(8);
   im.push_back(7);
   im.push_back(6);
   im.push_back(6);
   
   // Now ready for minimization step
   arglist[0] = 50000000;
   arglist[1] = 1.;
   for(size_t i=0;i<nps;i++)
      gMinuit->FixParameter(i);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   //for(size_t i=0;i<nps;i++)
   //   gMinuit->Release(i);
   //for(size_t i=0;i<nrun;i++)
   //   gMinuit->FixParameter(nps+i);
   //gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   //for(size_t i=0;i<nrun;i++)
   //   gMinuit->Release(nps+i);
   //gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   
   // Print results
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   //gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   //gMinuit->mnprin(3,amin);
   TFile *fout = TFile::Open("fitall_2017.root","RECREATE");      
   
   TCanvas *c = new TCanvas("c","c");
   std::vector<double> dt(t.size(),0);
   TGraphErrors *gr = new TGraphErrors(t.size(),t.data(),sz[flayer].data(),dt.data(),dsz[flayer].data());
   gr->GetXaxis()->SetTimeDisplay(1);
   gr->SetMarkerSize(0.5);
   gr->SetMarkerStyle(20);
   gr->Draw("AP");
   
   std::vector<double> par; 
   std::vector<double> dpar;
   for(size_t i = 0;i<npar;i++) {
      double p, dp;
      gMinuit->GetParameter(i,p,dp);
      std::cout << i << " " << p << " " << dp << std::endl;
      par.push_back(p);
      dpar.push_back(dp);
   }
   for(size_t i=0;i<f_sz.size();i++)
      f_sz[i].function(par.data());
   
   c->Update();
   c->Write();
   int xxx;
   std::cin >> xxx;
   
   TGraphErrors *gr1 = new TGraphErrors(t.size(),t.data(),dedxm[flayer].data(),dt.data(),ddedxm[flayer].data());
   gr1->GetXaxis()->SetTimeDisplay(1);
   gr1->SetMarkerSize(0.5);
   gr1->SetMarkerStyle(20);
   gr1->Draw("AP");
   
   for(size_t i=0;i<f_alpha.size();i++)
      f_alpha[i].function(par.data());
   
   c->Update();
   c->Write();   
   std::cin >> xxx;
   
   TGraphErrors *gr0 = new TGraphErrors(t.size(),t.data(),dedx0[flayer].data(),dt.data(),ddedx0[flayer].data());
   gr0->GetXaxis()->SetTimeDisplay(1);
   gr0->SetMarkerSize(0.5);
   gr0->SetMarkerStyle(20);
   gr0->Draw("AP");
   
   std::vector<double> dedx0f;
   std::vector<double> ddedx0f;
   for(size_t i=0;i<nrun;i++) {
      dedx0f.push_back(par[nps+i]);
      ddedx0f.push_back(dpar[nps+i]);
   }
   TGraphErrors *gr0f = new TGraphErrors(t.size(),t.data(),dedx0f.data(),dt.data(),ddedx0f.data());
   gr0f->GetXaxis()->SetTimeDisplay(1);
   gr0f->SetMarkerSize(0.5);
   gr0f->SetMarkerStyle(20);
   gr0f->SetLineColor(kRed);
   gr0f->SetMarkerColor(kRed);
   gr0f->Draw("P");
   
   c->Update();
   c->Write();
   std::cin >> xxx;
   
   TGraphErrors *gr2 = new TGraphErrors(t.size(),t.data(),slope[flayer].data(),dt.data(),dslope[flayer].data());
   gr2->GetXaxis()->SetTimeDisplay(1);
   gr2->SetMarkerSize(0.5);
   gr2->SetMarkerStyle(20);
   gr2->Draw("AP");
   
   for(size_t i=0;i<f_slope.size();i++)
      f_slope[i].function(par.data());
   
   c->Update();
   c->Write();
   fout->Close();
   
   std::cout << "fcn = " << fcn_min << " " << npars << " " << npar << std::endl;
}

//------------------------------------------------------------------------------------

double satfsx(double *x, double *par)
{
   double theta = x[0];
   double sz = par[0];
   double dphi = theta-TMath::Pi()/2;
   double dedx = par[1] + par[3]*dphi;// + par[4]*std::pow(std::fabs(dphi),par[5]);
   double demax = par[2]/par[0];
   double g0, g1, g2;
   double sth = sin(theta);
   double tth = sth/cos(theta);
   double b = par[4]/tth;
   //double b2 = par[5]/tth/tth;
   //double b = b1 + b2;
   double k = sth*sqrt(1 + b*b);
   return k*satfsx(dedx, theta, sz, demax, g0, g1, g2 );
}

class dEdXclb
{
 public:
   dEdXclb(){fitted = false;};
   TFitResultPtr fits(size_t mask = 2);
  TF1* plot(TFitResultPtr r);
  bool IsEmpty();

   int run;
   int layer;
   size_t t;
   size_t dt;
   size_t tmin;
   size_t tmax;
   size_t ind;
   bool fitted;
   double sz;
   double dedx0;
   double alpha;
   double slope;
   double scale;
   TProfile *h;
   TFitResult *r;
   TFitResultPtr rf;
   TH1D *h_zCorr;
};

bool dEdXclb::IsEmpty()
{
   double xmin = TMath::Pi()/25*(13-ind);
   double xmax = TMath::Pi() - xmin;
   bool empty = true;
   for(size_t i=0;i<h->GetNbinsX();i++) {
     double x = h->GetBinCenter(i);
     double y = h->GetBinContent(i);
     if (x>=xmin&&x<=xmax&&y>0) empty = false;
   }
   return empty;
}

TFitResultPtr dEdXclb::fits(size_t mask)
{
   fitted = true;
   TF1 *f = new TF1("satfsx",satfsx,0,TMath::Pi(),5);
   f->SetLineWidth(0.1);
   f->SetParName(0,"sz");
   f->SetParName(1,"dedx0");
   f->SetParName(2,"alpha");
   f->SetParName(3,"slope");
   f->SetParName(4,"scale");
   f->SetParameter(0,sz);
   f->SetParameter(1,dedx0);
   f->SetParameter(2,alpha);
   f->SetParameter(3,slope);
   f->SetParameter(4,scale);
   f->SetParLimits(1,0,20000);
   f->SetParLimits(3,0,10000);
   f->SetParLimits(4,0,10000);
   if (!(mask&1)) f->FixParameter(0,sz);
   if (!(mask&2)) f->FixParameter(1,dedx0);
   if (!(mask&4)) f->FixParameter(2,alpha);
   if (!(mask&8)) f->FixParameter(3,slope);
   if (!(mask&16)) f->FixParameter(4,scale);
   double xmin = TMath::Pi()/25*(13-ind);
   double xmax = TMath::Pi() - xmin;
   h->Fit(f,"SQ","",xmin,xmax);
   //h->Fit(f,"S","",xmin,xmax);
   //f->ReleaseParameter(4);
   //h->Fit(f,"S","",xmin,xmax);
   rf = h->Fit(f,"S","",xmin,xmax);
   return rf;
}

TF1* dEdXclb::plot(TFitResultPtr r)
{
   fitted = true;
   TF1 *f = new TF1("satfsx",satfsx,0,TMath::Pi(),5);
   f->SetLineWidth(0.1);
   f->SetParName(0,"sz");
   f->SetParName(1,"dedx0");
   f->SetParName(2,"alpha");
   f->SetParName(3,"slope");
   f->SetParName(4,"scale");
   f->SetParameter(0,r->Parameters().at(0));
   f->SetParameter(1,r->Parameters().at(1));
   f->SetParameter(2,r->Parameters().at(2));
   f->SetParameter(3,r->Parameters().at(3));
   f->SetParameter(4,r->Parameters().at(4));
   f->Draw("same");
   return f;
}

class cPars
{
 public:
   cPars(){};
   double getSz(size_t t);
   double getAlpha(size_t t);
   double getSlope(size_t t);
   double getScale(size_t t);
   
   std::vector<fParN*> sz;
   std::vector<fParN*> alpha;
   std::vector<fParN*> slope;
   std::vector<fParN*> scale;
};

double cPars::getSz(size_t t)
{
   size_t ind = 100;
   for(size_t q=0;q<sz.size();q++) if (sz[q]->IsInTime(t)) ind = q;
   if (ind==100) std::cout << "getSz: " << t << " " << sz.front()->xmin << " " << sz.back()->xmax << std::endl;
   return sz[ind]->Eval(t);
}

double cPars::getAlpha(size_t t)
{
   size_t ind = 100;
   for(size_t q=0;q<alpha.size();q++) if (alpha[q]->IsInTime(t)) ind = q;
   if (ind==100) std::cout << "getAlpha: " << t << std::endl;
   return alpha[ind]->Eval(t);
}

double cPars::getSlope(size_t t)
{
   size_t ind = 100;
   for(size_t q=0;q<slope.size();q++) if (slope[q]->IsInTime(t)) ind = q;
   if (ind==100) std::cout << "getSlope: " << t << std::endl;
   return slope[ind]->Eval(t);
}

double cPars::getScale(size_t t)
{
   size_t ind = 100;
   for(size_t q=0;q<scale.size();q++) if (scale[q]->IsInTime(t)) ind = q;
   if (ind==100) std::cout << "getScale: " << t << std::endl;
   return scale[ind]->Eval(t);
}

std::vector<dEdXclb> crun;
std::vector<dEdXclb> pcgen;

void plotRun(size_t run, std::vector<dEdXclb> &pcgen)
{
   std::vector<size_t> ih;
   for(size_t i=0;i<pcgen.size();i++)
      if (pcgen[i].run==run && pcgen[i].fitted) 
         ih.push_back(i);
   if (!ih.size()) return;
   TCanvas *c = new TCanvas("c","c",900,600);
   c->Divide(3,3);
   for(size_t i=0;i<ih.size();i++) {
      c->cd(i+1);
      pcgen[ih[i]].h->Draw();
      double par4 = ((TF1*)pcgen[ih[i]].h->FindObject("satfsx"))->GetParameter(4);
      double dpar4 = ((TF1*)pcgen[ih[i]].h->FindObject("satfsx"))->GetParError(4);
      double par5 = ((TF1*)pcgen[ih[i]].h->FindObject("satfsx"))->GetParameter(5);
      double dpar5 = ((TF1*)pcgen[ih[i]].h->FindObject("satfsx"))->GetParError(5);
      std::cout << (1-par4)/dpar4 << " " << par4 << " " << dpar4 << " | "
         << par5/dpar5 << " " << par5 << " " << dpar5 << std::endl;
   }
   c->Update();
   size_t xxx;
   cin >> xxx;
}

void reFitAll(size_t flayer = 0)
{
   im.resize(0);
   //im.push_back(7);
   //im.push_back(9);
   //im.push_back(9);
   //im.push_back(9);
   //im.push_back(8);
   //im.push_back(7);
   //im.push_back(7);
   //im.push_back(6);
   //im.push_back(6);
   //im.push_back(5);
   //
   im.push_back(8);
   im.push_back(10);
   im.push_back(10);
   im.push_back(9);
   im.push_back(9);
   im.push_back(9);
   im.push_back(9);
   im.push_back(9);
   im.push_back(8);
   im.push_back(8);
   //
   //TFile *f = new TFile("mhad2017-4.root");
   //TFile *f = new TFile("mhad2017new.root");
   //TFile *f = new TFile("mhad2017new_col-4.root");
   //TFile *f = new TFile("mhad2017new_mhad2017-4.root");
   //TFile *f = new TFile("rho-2018.root");
   //TFile *f = new TFile("mhad2019.root");
   std::vector<TFile*> f;
   // 103 изменил условия отбора. усилил
   //f.push_back( new TFile("mhad2017-4_0.000000_103.000000.root"));
   // 104 добывил условие dEdx>0
   //f.push_back( new TFile("mhad2017-4_0.000000_104.000000.root"));
   // 105 добывил fix: scale = 0
   //f.push_back( new TFile("mhad2017-4_0.000000_105.000000.root"));
   //f.push_back(new TFile("rho_2018-3_0.000000_25.000000.root"));
   //f.push_back(new TFile("rho_2018-3_25.000000_50.000000.root"));
   //f.push_back(new TFile("rho_2018-3_50.000000_75.000000.root"));
   //f.push_back(new TFile("rho_2018-3_75.000000_105.000000.root"));
   //
   //f.push_back(new TFile("mhad2010-4_0.000000_103.000000.root"));
   //f.push_back(new TFile("mhad2011-6_0.000000_103.000000.root"));
   //f.push_back(new TFile("mhad2012-6_0.000000_103.000000.root"));
   //f.push_back(new TFile("rho_2013-4_0.000000_103.000000.root"));
   //
   //f.push_back(new TFile("mhad2017-5_0.000000_103.000000.root"));
   //f.push_back(new TFile("mhad2019-2_0.000000_103.000000.root"));
   //f.push_back(new TFile("mhad2020-1_0.000000_103.000000.root"));
   //f.push_back(new TFile("rho_2018-3_0.000000_26.000000.root"));
   //f.push_back(new TFile("rho_2018-3_25.000000_51.000000.root"));
   //f.push_back(new TFile("rho_2018-3_50.000000_76.000000.root"));
   //f.push_back(new TFile("rho_2018-3_75.000000_101.000000.root"));
   //f.push_back(new TFile("rho_2013-4_0.000000_26.000000.root"));
   //f.push_back(new TFile("rho_2013-4_25.000000_51.000000.root"));
   //f.push_back(new TFile("rho_2013-4_50.000000_76.000000.root"));
   //f.push_back(new TFile("rho_2013-4_75.000000_101.000000.root"));
   //
   f.push_back(new TFile("dEdxCalibData/v0/mhad2017-5/mhad2017-5_v0_0.000000_101.000000.root"));
   //f.push_back(new TFile("dEdxCalibData/v1/mhad2017-5/mhad2017-5_v1_0.000000_101.000000.root"));
   TString pref;
   TString dir;
   double tmin = 1e20;
   double tmax = 0;
   for(size_t ifp=0;ifp<f.size();ifp++) {
      TString name(f[ifp]->GetName());
      //pref = name(0,name.Index(".root"));
      pref = name(name.Last('/')+1,name.Index(".root")-name.Last('/')-1);
      dir = name(0,name.Last('/'));
      std::cout << pref << " " << dir << std::endl;
   TIter next_run(f[ifp]->GetListOfKeys());
   TKey *key_run;
   TString key_run_name("TDirectoryFile");
   TString key_name("TFitResult");
   TString key_name_fun("f_fits_");
   while ((key_run=(TKey*)next_run())) {
      if (!key_run_name.BeginsWith(key_run->GetClassName())) continue;
      TDirectory *d = (TDirectory*)f[ifp]->Get(key_run->GetName());
      TIter next(d->GetListOfKeys());
      TKey *key;
      while ((key=(TKey*)next())) {
         TString name(key->GetName());
         if (name.BeginsWith("TFitResult-pcgen")&&
             name.Contains(key_name_fun)) {
            pcgen.push_back(dEdXclb());
            pcgen.back().r = (TFitResult*)d->Get(key->GetName());
            sscanf(key->GetName(),"TFitResult-pcgen%d_layer%d",&pcgen.back().run,&pcgen.back().layer);
            pcgen.back().h = (TProfile*)d->Get(Form("pcgen%08d_layer%d",pcgen.back().run,pcgen.back().layer));
            pcgen.back().ind = im[pcgen.back().layer];
            pcgen.back().h_zCorr = (TH1D*)d->Get(Form("pcgen%08d_zCorr",pcgen.back().run));
            std::vector<size_t> *t = (std::vector<size_t> *)d->Get(Form("time_pcgen%08d",pcgen.back().run));
            std::cout << pcgen.back().run << " " << (*t)[0] << " " << (*t)[1] << std::endl;
            pcgen.back().t = t->size()<1 ? 0 : (*t)[0];
            pcgen.back().dt = t->size()<2 ? 0 : (*t)[1];
            pcgen.back().tmin = t->size()<3 ? 0 : (*t)[2];
            pcgen.back().tmax = t->size()<4 ? 0 : (*t)[3];
            //pcgen.back().t = ftime(pcgen.back().run);
            //pcgen.back().dt = 0;
         }
         if (name.BeginsWith("TFitResult-run")&&
             name.Contains(key_name_fun)) {
            crun.push_back(dEdXclb());
            crun.back().r = (TFitResult*)d->Get(key->GetName());
            sscanf(key->GetName(),"TFitResult-run%d_layer%d",&crun.back().run,&crun.back().layer);
            crun.back().h = (TProfile*)d->Get(Form("run%08d_layer%d",crun.back().run,crun.back().layer));
            crun.back().ind = im[crun.back().layer];
            crun.back().h_zCorr = (TH1D*)d->Get(Form("run%08d_zCorr",crun.back().run));
            std::vector<size_t> *t = (std::vector<size_t> *)d->Get(Form("time_run%08d",crun.back().run));
            crun.back().t = t->size()<1 ? 0 : (*t)[0];
            crun.back().dt = t->size()<2 ? 0 : (*t)[1];
            crun.back().tmin = t->size()<3 ? 0 : (*t)[2];
            crun.back().tmax = t->size()<4 ? 0 : (*t)[3];
            //crun.back().t = ftime(crun.back().run);
            //crun.back().dt = 0;
            if (tmin>crun.back().t) tmin = crun.back().t;
            if (tmax<crun.back().t) tmax = crun.back().t;
         }
      }
   }
   }
   //
   std::vector<double> zCorr;
   std::vector<double> dzCorr;
   for(size_t i=0;i<pcgen.size();i++) {
      TFitResult *r(pcgen[i].r);
      // run/time
      size_t run = pcgen[i].run;
      size_t layer = pcgen[i].layer;
      if (layer==flayer) {
         t.push_back((pcgen[i].tmax+pcgen[i].tmin)/2);
         dt.push_back((pcgen[i].tmax-pcgen[i].tmin)/2);
         std::cout << pcgen[i].t << " " << pcgen[i].dt << std::endl;
      }
      Double_t scal = std::sqrt(r->Chi2()/r->Ndf());
      
      // parameters
      sz[layer].push_back(r->Parameters().at(0));
      dsz[layer].push_back(r->Errors().at(0)*scal+1e-6);
      dedx0[layer].push_back(r->Parameters().at(1));
      ddedx0[layer].push_back(r->Errors().at(1)*scal+1e-6);
      dedxm[layer].push_back(r->Parameters().at(2));
      ddedxm[layer].push_back(r->Errors().at(2)*scal+1e-6);
      slope[layer].push_back(r->Parameters().at(3));
      dslope[layer].push_back(r->Errors().at(3)*scal+1e-6);
      scale[layer].push_back(r->Parameters().at(4));
      dscale[layer].push_back(r->Errors().at(4)*scal+1e-6);
      if (layer==flayer) {
         //zCorr.push_back(4*(pcgen[i].h_zCorr->GetMean()-1)+1);
         zCorr.push_back(pcgen[i].h_zCorr->GetMean());
         dzCorr.push_back(pcgen[i].h_zCorr->GetMeanError());
      }
   }
   //std::vector<double> dt(t.size(),0);
   TGraphErrors *grz = new TGraphErrors(t.size(),t.data(),zCorr.data(),dt.data(),dzCorr.data());
   grz->GetXaxis()->SetTimeDisplay(1);
   grz->SetMarkerSize(0.5);
   grz->SetMarkerStyle(20);
   grz->SetLineColor(kBlue);
   grz->SetMarkerColor(kBlue);
   grz->Draw("AP");
//   stop();
   TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
   c->Update();
   c->SaveAs(Form("%s/%s_dTh.jpg",dir.Data(),pref.Data()));
   size_t xxxxx;
   //cin>>xxxxx;
   //
   //std::vector<DcPaGen> cal = getDcPaGen(t.front(),t.back());
   //std::vector<PcGen> cal = getPcGen(t.front(),t.back());
   //std::vector<DcZabs> cal = getDcZabs(t.front(),t.back(),false,"CURRENT");
   //std::vector<DcZabs> cal = getDcZabs(21803,29446,false,"CURRENT");
   //std::vector<DcZabs> cal = getDcZabs(37920,42826,false,"CURRENT");
   //std::vector<DcZabs> cal = getDcZabs(21803,26815,false,"CURRENT");
   //std::vector<DcZabs> cal = getDcZabs(21803,29446,false,"KB");
/*   for(size_t ip=0;ip<cal[0].z1.size();ip++) {
      grz->Draw("AP");
      grz->GetYaxis()->SetRangeUser(0.9,1.1);
      std::vector<double> x;
      std::vector<double> y;
      for(size_t i=0;i<cal.size();i++) {
         if (cal[i].z1.size()) {
            x.push_back(cal[i].t);
            y.push_back(std::fabs(cal[i].z1[ip]));
         }
      }
      x.push_back(t.back());
      y.push_back(y.back());
      //std::vector<double>::iterator ymax = std::max_element(y.begin(), y.end());
      //double ymax = -1e10;
      //for(size_t i=0;i<y.size();i++) if (ymax<y[i]) ymax = y[i];
      //std::cout << (*ymax) << " " << x.front() << " " << x.back() << " " << t.front() << " " << t.back() << std::endl;
      double norm = y[0];
      for(size_t i=0;i<y.size();i++)
         y[i] /= norm;
      TGraph *grc = new TGraph(x.size(),x.data(),y.data());
      grc->GetXaxis()->SetTimeDisplay(1);
      grc->SetMarkerSize(0.5);
      grc->SetMarkerStyle(20);
      grc->SetLineColor(kRed);
      grc->SetMarkerColor(kRed);
      grc->Draw("PL");
      //grc->GetXaxis()->SetRangeUser(fdatetotime(20170113),fdatetotime(20170715));
      TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
      c->Update();
      //c->SaveAs(Form("pict/dczabs/z1%d.jpg",ip));
      size_t xxxxx;
      cin>>xxxxx;
   }*/
   //---------------------------------------
   
   nmax = dedx0[flayer].size();
   size_t nrun = dedx0[flayer].size();
   size_t npar = 0;
   size_t npar0 = 0;
   //
   std::vector<cPars> pars(10,cPars());
   std::vector<cPars> par0(10,cPars());
   size_t xxxx;
   std::vector<size_t> flyr;
   flyr.push_back(flayer);
   //flyr.push_back(0);
   //flyr.push_back(1);
   //flyr.push_back(2);
   //flyr.push_back(3);
   //flyr.push_back(4);
   //flyr.push_back(5);
   //flyr.push_back(6);
   //flyr.push_back(7);
   //flyr.push_back(8);
   //flyr.push_back(9);
   double t1 = fdatetotime(20120313) + 9*3600;
   TFile *fp = new TFile(Form("%s/%s_pars_x_layer%d.root",dir.Data(),pref.Data(),flyr[0]),"RECREATE");
   for(size_t l=0;l<flyr.size();l++) {
      //
      par0[flyr[l]].alpha.push_back(new fParN(6,tmin,tmax,Form("dEdx: layer = %d",flyr[l])));
      for(size_t i=0;i<par0[flyr[l]].alpha.size();i++) {
         par0[flyr[l]].alpha[i]->init(npar0,t,dedx0[flyr[l]],ddedx0[flyr[l]]);
         npar0 += par0[flyr[l]].alpha[i]->nm+1;
         par0[flyr[l]].alpha[i]->plot(0,3000);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         c->Write();
         c->SaveAs(Form("%s/%s_dEdx_layer%d_%d.jpg",dir.Data(),pref.Data(),flyr[l],i));
         //cin>>xxxx;
      }
      if (((TString)f[0]->GetName()).Contains("mhad2012")) 
      {
         double tmd1 = fdatetotime(20120125) + 10*3600;
         double tmd2 = fdatetotime(20120215) + 10*3600;
         pars[flyr[l]].sz.push_back(new fParN(0,tmin,tmd1,Form("sz: layer = %d",flyr[l])));
         pars[flyr[l]].sz.push_back(new fParN(1,tmd1,tmd2,Form("sz: layer = %d",flyr[l])));
         pars[flyr[l]].sz.push_back(new fParN(6,tmd2,tmax,Form("sz: layer = %d",flyr[l])));
      } 
      else if (((TString)f[0]->GetName()).Contains("rho_2013")) 
      {
         double tmd1 = fdatetotime(20130130);
         double tmd2 = fdatetotime(20130412) + 10*3600;
         pars[flyr[l]].sz.push_back(new fParN(2,tmin,tmd1,Form("sz: layer = %d",flyr[l])));
         pars[flyr[l]].sz.push_back(new fParN(4,tmd1,tmd2,Form("sz: layer = %d",flyr[l])));
         pars[flyr[l]].sz.push_back(new fParN(4,tmd2,tmax,Form("sz: layer = %d",flyr[l])));
      } 
      else if (((TString)f[0]->GetName()).Contains("rho_2018")) 
      {
         double tmd1 = fdatetotime(20180107);
         double tmd2 = fdatetotime(20180305);
         pars[flyr[l]].sz.push_back(new fParN(2,tmin,tmd1,Form("sz: layer = %d",flyr[l])));
         pars[flyr[l]].sz.push_back(new fParN(4,tmd1,tmd2,Form("sz: layer = %d",flyr[l])));
         pars[flyr[l]].sz.push_back(new fParN(4,tmd2,tmax,Form("sz: layer = %d",flyr[l])));
      } 
      else
         pars[flyr[l]].sz.push_back(new fParN(6,tmin,tmax,Form("sz: layer = %d",flyr[l])));
      for(size_t i=0;i<pars[flyr[l]].sz.size();i++) {
         pars[flyr[l]].sz[i]->init(npar,t,sz[flyr[l]],dsz[flyr[l]]);
         npar += pars[flyr[l]].sz[i]->nm+1;
         pars[flyr[l]].sz[i]->plot(0,0.1);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         c->Write();
         c->SaveAs(Form("%s/%s_sz_layer%d_%d.jpg",dir.Data(),pref.Data(),flyr[l],i));
         //cin>>xxxx;
      }
      //
      if (((TString)f[0]->GetName()).Contains("mhad2011")) 
      {
         pars[flyr[l]].alpha.push_back(new fParN(6,tmin,fdatetotime(20110305),Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(6,fdatetotime(20110305),tmax,Form("alpha: layer = %d",flyr[l])));
      } 
      else if (((TString)f[0]->GetName()).Contains("mhad2012")) 
      {
         double tmd1 = fdatetotime(20120119) + 10*3600;
         double tmd2 = fdatetotime(20120215) + 10*3600;
         double tmd3 = fdatetotime(20120313) + 13*3600;
         pars[flyr[l]].alpha.push_back(new fParN(1,tmin,tmd1,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(2,tmd1,tmd2,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(2,tmd2,tmd3,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(4,tmd3,tmax,Form("alpha: layer = %d",flyr[l])));
      } 
      else if (((TString)f[0]->GetName()).Contains("rho_2013")) 
      {
         double tmd1 = fdatetotime(20130130);
         double tmd2 = fdatetotime(20130412) + 10*3600;
         double tmd3 = fdatetotime(20130625) + 10*3600;
         pars[flyr[l]].alpha.push_back(new fParN(4,tmin,tmd1,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(4,tmd1,tmd2,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(4,tmd2,tmd3,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(2,tmd3,tmax,Form("alpha: layer = %d",flyr[l])));
      } 
      else if (((TString)f[0]->GetName()).Contains("mhad2017")) 
      {
         double tmid = fdatetotime(20170319);
         pars[flyr[l]].alpha.push_back(new fParN(4,tmin,tmid,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(6,tmid,tmax,Form("alpha: layer = %d",flyr[l])));
      } 
      else if (((TString)f[0]->GetName()).Contains("rho_2018")) 
      {
         double tmd1 = fdatetotime(20180107);
         double tmd2 = fdatetotime(20180129);
         double tmd3 = fdatetotime(20180319) + 10*3600;
         double tmd4 = fdatetotime(20180325);
         double tmd5 = fdatetotime(20180426) + 10*3600;
         pars[flyr[l]].alpha.push_back(new fParN(6,tmin,tmd1,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(2,tmd1,tmd2,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(8,tmd2,tmd3,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(0,tmd3,tmd4,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(4,tmd4,tmd5,Form("alpha: layer = %d",flyr[l])));
         pars[flyr[l]].alpha.push_back(new fParN(4,tmd5,tmax,Form("alpha: layer = %d",flyr[l])));
      } 
      else
         pars[flyr[l]].alpha.push_back(new fParN(6,tmin,tmax,Form("alpha: layer = %d",flyr[l])));
      for(size_t i=0;i<pars[flyr[l]].alpha.size();i++) {
         pars[flyr[l]].alpha[i]->init(npar,t,dedxm[flyr[l]],ddedxm[flyr[l]]);
         npar += pars[flyr[l]].alpha[i]->nm+1;
         pars[flyr[l]].alpha[i]->plot(0,300);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         c->Write();
         c->SaveAs(Form("%s/%s_alpha_layer%d_%d.jpg",dir.Data(),pref.Data(),flyr[l],i));
         //cin>>xxxx;
      }
      //
      if (((TString)f[0]->GetName()).Contains("mhad2011")) 
      {
         double tmid = fdatetotime(20110424);
         pars[flyr[l]].slope.push_back(new fParN(6,tmin,tmid,Form("slope: layer = %d",flyr[l])));
         pars[flyr[l]].slope.push_back(new fParN(4,tmid,tmax,Form("slope: layer = %d",flyr[l])));
      }
      else if (((TString)f[0]->GetName()).Contains("mhad2012")) 
      {
         double tmd1 = fdatetotime(20120125) + 10*3600;
         double tmd2 = fdatetotime(20120215) + 10*3600;
         pars[flyr[l]].slope.push_back(new fParN(0,tmin,tmd1,Form("slope: layer = %d",flyr[l])));
         pars[flyr[l]].slope.push_back(new fParN(2,tmd1,tmd2,Form("slope: layer = %d",flyr[l])));
         pars[flyr[l]].slope.push_back(new fParN(6,tmd2,tmax,Form("slope: layer = %d",flyr[l])));
      } 
      else if (((TString)f[0]->GetName()).Contains("mhad2017")) 
      {
         double tmid = fdatetotime(20170319);
         pars[flyr[l]].slope.push_back(new fParN(4,tmin,tmid,Form("slope: layer = %d",flyr[l])));
         pars[flyr[l]].slope.push_back(new fParN(6,tmid,tmax,Form("slope: layer = %d",flyr[l])));
      } 
      else if (((TString)f[0]->GetName()).Contains("rho_2018")) 
      {
         double tmd1 = fdatetotime(20180107);
         double tmd2 = fdatetotime(20180305);
         pars[flyr[l]].slope.push_back(new fParN(2,tmin,tmd1,Form("slope: layer = %d",flyr[l])));
         pars[flyr[l]].slope.push_back(new fParN(4,tmd1,tmd2,Form("slope: layer = %d",flyr[l])));
         pars[flyr[l]].slope.push_back(new fParN(4,tmd2,tmax,Form("slope: layer = %d",flyr[l])));
      } 
      else
         pars[flyr[l]].slope.push_back(new fParN(8,tmin,tmax,Form("slope: layer = %d",flyr[l])));
      for(size_t i=0;i<pars[flyr[l]].slope.size();i++) {
         pars[flyr[l]].slope[i]->init(npar,t,slope[flyr[l]],dslope[flyr[l]]);
         npar += pars[flyr[l]].slope[i]->nm+1;
         pars[flyr[l]].slope[i]->plot(-100,100);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         c->Write();
         c->SaveAs(Form("%s/%s_slope_layer%d_%d.jpg",dir.Data(),pref.Data(),flyr[l],i));
         //cin>>xxxx;
      }
      //
      pars[flyr[l]].scale.push_back(new fParN(6,tmin,tmax,Form("scale: layer = %d",flyr[l])));
      for(size_t i=0;i<pars[flyr[l]].scale.size();i++) {
         pars[flyr[l]].scale[i]->init(npar,t,scale[flyr[l]],dscale[flyr[l]]);
         npar += pars[flyr[l]].scale[i]->nm+1;
         pars[flyr[l]].scale[i]->plot(0,2);
         //grz->Draw("P");
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         c->Write();
         c->SaveAs(Form("%s/%s_scale_layer%d_%d.jpg",dir.Data(),pref.Data(),flyr[l],i));
         //cin>>xxxx;
      }
   }
   //fp->Close();
   npar += nrun;
   //
   std::vector< std::vector<double> > dedx0r(10,std::vector<double>());
   std::vector< std::vector<double> > ddedx0r(10,std::vector<double>());
   std::vector< std::vector<double> > tr(10,std::vector<double>());
   std::vector< std::vector<double> > dtr(10,std::vector<double>());
   std::vector< std::vector<double> > dedx0f(10,std::vector<double>());
   std::vector< std::vector<double> > ddedx0f(10,std::vector<double>());
   std::vector< std::vector<double> > tf(10,std::vector<double>());
   std::vector< std::vector<double> > dtf(10,std::vector<double>());
   std::vector< std::vector<double> > rf(10,std::vector<double>());
   std::vector< std::vector<double> > drf(10,std::vector<double>());
   std::vector< std::vector<double> > rs(10,std::vector<double>());
   TH1D *hs = new TH1D("r","r",100,0,10);
   size_t runToPlot = 0;
   for(size_t i=0;i<crun.size();i++) {
     //if (crun[i].IsEmpty()) continue;
      size_t run = crun[i].run;
      size_t layer = crun[i].layer;
      TFitResult *r0(crun[i].r);
      dedx0r[layer].push_back(r0->Parameters().at(1));
      ddedx0r[layer].push_back(r0->Errors().at(1));
      //if (run!=runToPlot) plotRun(runToPlot,crun);
      runToPlot = run;
      //if (crun[i].layer!=flayer) continue;
      if (std::find(flyr.begin(), flyr.end(), layer) == flyr.end()) continue;
      crun[i].dedx0 = par0[layer].getAlpha(ftime(crun[i].run));
      crun[i].sz = pars[layer].getSz(ftime(crun[i].run));
      crun[i].alpha = pars[layer].getAlpha(ftime(crun[i].run));
      crun[i].slope = pars[layer].getSlope(ftime(crun[i].run));
      crun[i].scale = 1;//pars[layer].getScale(ftime(crun[i].run));
      TFitResultPtr r = crun[i].fits();
      if (TMath::IsNaN(r->Parameters().at(1))) {
         std::cout << crun[i].sz << " " << crun[i].alpha << " " << crun[i].slope << " " << crun[i].scale << std::endl;
      }
      std::cout << r->GetName() << std::endl;
      r->Write();
      dedx0f[layer].push_back(r->Parameters().at(1));
      ddedx0f[layer].push_back(r->Errors().at(1));
      tf[layer].push_back(crun[i].t);
      dtf[layer].push_back(crun[i].dt);
      rf[layer].push_back(crun[i].run);
      drf[layer].push_back(0);
      rs[layer].push_back(crun[i].r->Errors().at(1)/r->Errors().at(1));
      hs->Fill(rs[layer].back());
      //crun[i].h->Draw();
      //std::cout << crun[i].r->Parameters().at(1) << " " << crun[i].r->Errors().at(1) << std::endl;
      //TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
      //c->Update();
      //for(size_t k=0;k<crun[i].h->GetNbinsX();k++)
        // std::cout << k << " " << crun[i].h->GetBinContent(k) << " " << crun[i].h->GetBinError(k) << std::endl;
      //cin>>xxxx;
   }
   {
      //TGraphErrors *gr = new TGraphErrors(t.size(),t.data(),dedx0[flayer].data(),dt.data(),ddedx0[flayer].data());
      TGraphErrors *gr = new TGraphErrors(tf[flayer].size(),tf[flayer].data(),dedx0r[flayer].data(),dtf[flayer].data(),ddedx0r[flayer].data());
      gr->GetXaxis()->SetTimeDisplay(1);
      gr->SetMarkerSize(0.5);
      gr->SetMarkerStyle(20);
      gr->Draw("AP");
      for(size_t flayer=0;flayer<10;flayer++)
         std::cout << flayer << " " << tf[flayer].size() << " " << dtf[flayer].size() << " " << dedx0f[flayer].size() << " " << ddedx0f[flayer].size() << std::endl;
      TGraphErrors *grf = new TGraphErrors(tf[flayer].size(),tf[flayer].data(),dedx0f[flayer].data(),dtf[flayer].data(),ddedx0f[flayer].data());
      grf->GetXaxis()->SetTimeDisplay(1);
      grf->SetMarkerSize(0.5);
      grf->SetMarkerStyle(20);
      grf->SetMarkerColor(kRed);
      grf->SetLineColor(kRed);
      grf->Draw("P");
      //TFile *fp = new TFile(Form("%s/%s_dEdx0_layer%d.root",dir.Data(),pref.Data(),flyr[0]),"RECREATE");
      {
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetTitle("dEdx_fit");
         c->SetName("dEdx_fit");
         c->Update();
         c->Write();
      }
   }
   {
      TGraphErrors *gr = new TGraphErrors(rf[flayer].size(),rf[flayer].data(),dedx0r[flayer].data(),drf[flayer].data(),ddedx0r[flayer].data());
      //gr->GetXaxis()->SetTimeDisplay(1);
      gr->SetMarkerSize(0.5);
      gr->SetMarkerStyle(20);
      gr->Draw("AP");
      TGraphErrors *grf = new TGraphErrors(rf[flayer].size(),rf[flayer].data(),dedx0f[flayer].data(),drf[flayer].data(),ddedx0f[flayer].data());
      //grf->GetXaxis()->SetTimeDisplay(1);
      grf->SetMarkerSize(0.5);
      grf->SetMarkerStyle(20);
      grf->SetMarkerColor(kRed);
      grf->SetLineColor(kRed);
      grf->Draw("P");
      {
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetTitle("dEdx_fit_run");
         c->SetName("dEdx_fit_run");
         c->Update();
         c->Write();
      }
   }
   fp->Close();
   //
   //hs->Draw();
   TGraphErrors *grr = new TGraphErrors(tf[flayer].size(),tf[flayer].data(),rs[flayer].data(),dtf[flayer].data(),dtf[flayer].data());
   grr->GetXaxis()->SetTimeDisplay(1);
   grr->SetMarkerSize(0.5);
   grr->SetMarkerStyle(20);
   grr->SetMarkerColor(kRed);
   grr->SetLineColor(kRed);
   //grr->Draw("AP");
   //
   TFile *fs = new TFile(Form("%s/%s_dEdx_layer%d.root",dir.Data(),pref.Data(),flayer),"RECREATE");
   for(size_t i=0;i<crun.size();i++) {
      size_t run = crun[i].run;
      size_t layer = crun[i].layer;
      if (std::find(flyr.begin(), flyr.end(), layer) == flyr.end()) continue;
      std::cout << i << " " << run << " " << layer << std::endl;
      TString dname = Form("run%08d",run);
      TDirectory *dir = (TDirectory*)fs->Get(dname);
      if (!dir) dir = fs->mkdir(dname);
      dir->cd();
      crun[i].rf->Write();
      dir->cd("../");
   }
   fs->Close();
}

void caldEdxPrep(TString scan = "mhad2019-2_0.000000_103.000000")
{
   std::vector<TFile*> f;
   f.push_back(new TFile(Form("%s_dEdx_layer1.root",scan.Data())));
   f.push_back(new TFile(Form("%s_dEdx_layer2.root",scan.Data())));
   f.push_back(new TFile(Form("%s_dEdx_layer3.root",scan.Data())));
   f.push_back(new TFile(Form("%s_dEdx_layer4.root",scan.Data())));
   f.push_back(new TFile(Form("%s_dEdx_layer5.root",scan.Data())));
   f.push_back(new TFile(Form("%s_dEdx_layer6.root",scan.Data())));
   f.push_back(new TFile(Form("%s_dEdx_layer7.root",scan.Data())));
   f.push_back(new TFile(Form("%s_dEdx_layer8.root",scan.Data())));
   f.push_back(new TFile(Form("%s_dEdx_layer9.root",scan.Data())));
   std::map<Int_t,std::vector<Double_t> > cal;
   TString pref;
   for(size_t ifp=0;ifp<f.size();ifp++) {
      TString name(f[ifp]->GetName());
      //pref = name(0,name.Index(".root")-7);
      //pref = name(0,10);
      pref = name(0,name.Last('/')+11);
      TIter next_run(f[ifp]->GetListOfKeys());
      TKey *key_run;
      TString key_run_name("TDirectoryFile");
      TString key_name("TFitResult");
      TString key_name_fun("f_fits_");
      while ((key_run=(TKey*)next_run())) {
         if (!key_run_name.BeginsWith(key_run->GetClassName())) continue;
         TDirectory *d = (TDirectory*)f[ifp]->Get(key_run->GetName());
         TIter next(d->GetListOfKeys());
         TKey *key;
         while ((key=(TKey*)next())) {
            TString name(key->GetName());
            if (name.BeginsWith("TFitResult")) {
               TFitResult *r = (TFitResult*)d->Get(key->GetName());
               Int_t run;
               Int_t layer;
               sscanf(key->GetName(),"TFitResult-run%d_layer%d",&run,&layer);
               std::cout << name << " " << run << " " << layer << std::endl;
               if (!cal[run].size()) cal[run] = std::vector<Double_t>(45,0);
               cal[run][layer-1] = r->Parameters().at(1); // dedxe,dedx0
               cal[run][layer+8] = r->Parameters().at(2)/r->Parameters().at(0); // dedxm = alpha/sigz
               cal[run][layer+17] = r->Parameters().at(4); // scale
               cal[run][layer+26] = r->Parameters().at(0); // sigz,sz
               cal[run][layer+35] = r->Parameters().at(3); // slope
            }
         }
      }
   }
/*   ofstream out(Form("%s.txt",pref.Data()));
   std::map<Int_t,std::vector<Double_t> >::iterator it = cal.begin();
   for(;it!=cal.end();it++) {
      out << (*it).first << " ";
      for(size_t i=0;i<(*it).second.size();i++)
         out << std::scientific << (*it).second[i] << " ";
      out << std::endl;
   }
   out.close();*/
   std::vector<TString> suff;
   suff.push_back("dedxe");
   suff.push_back("dedxm");
   suff.push_back("scale");
   suff.push_back("sigz");
   suff.push_back("slope");
   for(size_t p=0;p<suff.size();p++)
   {
      ofstream out(Form("%s_%s.txt",pref.Data(),suff[p].Data()));
      std::map<Int_t,std::vector<Double_t> >::iterator it = cal.begin();
      for(;it!=cal.end();it++) {
         out << (*it).first << " ";
         for(size_t i=9*p;i<9*p+9;i++)
            out << std::scientific << (*it).second[i] << " ";
         out << std::endl;
      }
      out.close();
   }
}

void caldEdxTest(TString scan = "rho_2013-4_75.000000_101.000000_pars_x_layer1.root")
{
   TFile *f = new TFile(Form("pict/%s",scan.Data()));
   //TCanvas *c = (TCanvas*)f->Get("dEdx_fit_run");
   TCanvas *c = (TCanvas*)f->Get("alpha: layer = 1");
   c->Draw();
   //
   TString pref = scan(0,scan.Index(".root")-14);
   std::cout << Form("%s_dEdx.txt",pref.Data()) << std::endl;
   ifstream in(Form("%s_dEdx.txt",pref.Data()));
   std::vector<Double_t> run;
   std::vector<std::vector<Double_t> > dedx0(9,std::vector<Double_t>());
   std::vector<std::vector<Double_t> > alpha(9,std::vector<Double_t>());
   std::vector<std::vector<Double_t> > sz(9,std::vector<Double_t>());
   std::vector<std::vector<Double_t> > slope(9,std::vector<Double_t>());
   std::vector<std::vector<Double_t> > scale(9,std::vector<Double_t>());
   while (1) {
      Double_t par[46];
      for(size_t i=0;i<46;i++) in >> par[i];
      if (!in.good()) break;
      run.push_back(ftime(par[0]));
      for(size_t i=0;i<9;i++) {
         dedx0[i].push_back(par[i+1]);
         alpha[i].push_back(par[i+10]*par[i+28]);
         sz[i].push_back(par[i+28]);
         scale[i].push_back(par[i+19]);
         slope[i].push_back(par[i+37]);
      }
      std::cout << run.back()-run[0] << " " << par[0] << " " << alpha[0].back() << std::endl;
   }
   std::cout << run.size() << std::endl;
   Int_t layer = 1;
   TGraph *gr = new TGraph(run.size(),run.data(),alpha[layer-1].data());
   gr->GetXaxis()->SetTimeDisplay(1);
   gr->SetMarkerSize(0.5);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed);
   gr->SetLineColor(kRed);
   gr->SetLineColor(kBlue);
   gr->Draw("same");
}

void dczabs_pl()
{
  std::vector<DcZabs> cal = getDcZabs(21803,29446,false,"CURRENT");
  std::vector<DcZabs> calr = getDcZabs(21803,29446,false,"RECONSTR");
  //std::vector<DcZabs> cal = getDcZabs(21803,40000,false);
   std::cout << cal.size() << std::endl;
   double t1 = fdatetotime(20170310) + 12*3600;
   TH2D *hr = new TH2D("hr","hr",24,0,24,9,0,9);
   std::vector<TH1D*> hk;
   for(size_t l=0;l<9;l++)
     hk.push_back(new TH1D(Form("hk%d",l),"hk",100,0.8,1.2));
   for(size_t ip=0;ip<cal[0].z1.size();ip++) {
     std::vector<double> x;
     std::vector<double> y;
     double ymin = 1e5;
     double ymax = 0;
     double yb = 0;
     double ya = 0;
     size_t nb = 0;
     size_t na = 0;
     for(size_t i=0;i<cal.size();i++) {
       if (cal[i].z1.size()) {
         x.push_back(cal[i].t);
         y.push_back(std::fabs(cal[i].z1[ip]));
         if (y.back()>ymax) ymax = y.back();
         if (y.back()<ymin) ymin = y.back();
         if (x.back()<t1) {nb++;yb+=y.back();}
         else {na++;ya+=y.back();}
       }
     }
     yb/=nb;
     ya/=na;
     double k = ya/yb;
     hk[ip%9]->Fill(k);
     hr->SetBinContent(ip/9+1,ip%9+1,k);
     std::cout << ip%9 << " " << yb << " " << ya << " " << k << std::endl;
     TGraph *gr = new TGraph(x.size(),x.data(),y.data());
     gr->GetXaxis()->SetTimeDisplay(1);
     gr->SetMarkerSize(2);
     gr->SetMarkerStyle(20);
     gr->Draw("AP");
     TLine *l = new TLine();
     l->DrawLine(t1,0.9*ymin,t1,1.1*ymax);
      //
     std::vector<double> xr;
     std::vector<double> yr;
     for(size_t i=0;i<calr.size();i++) {
       if (calr[i].z1.size()) {
         xr.push_back(calr[i].t);
         yr.push_back(std::fabs(calr[i].z1[ip]));
       }
     }
     TGraph *grr = new TGraph(xr.size(),xr.data(),yr.data());
     grr->GetXaxis()->SetTimeDisplay(1);
     grr->SetMarkerSize(1);
     grr->SetMarkerStyle(20);
     grr->SetMarkerColor(kRed);
     grr->Draw("P");
      //
     TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
     c->Update();
     size_t xxx;
     cin>>xxx;
   }
   std::vector<double> x;
   std::vector<double> dx;
   std::vector<double> y;
   std::vector<double> dy;
   for(size_t l=0;l<hk.size();l++) {
     x.push_back(l);
     dx.push_back(0);
     y.push_back(hk[l]->GetMean());
     dy.push_back(hk[l]->GetMeanError());
   }
   TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
   gr->SetMarkerSize(0.5);
   gr->SetMarkerStyle(20);
   gr->Draw("AP");
   gStyle->SetPalette(57);
   hr->Draw("COLZ");
}

//---------------------------------------------------------------------

bool sortby(const std::pair<size_t,double> &a, 
            const std::pair<size_t,double> &b) 
{ 
   return (a.first < b.first); 
} 

Double_t p1p2(Double_t *x, Double_t *par)
{
   Double_t f = par[0]/x[0];
   return f;
}

Double_t p3p(Double_t *x, Double_t *par)
{
   Double_t f = par[0]/(x[0] - par[1]);
   return f;
}

class dEdxData
{
 public:
   dEdxData(){};
   
   std::vector<Double_t> x;
   std::vector<Double_t> y;
   std::vector<Double_t> dy;
   std::vector< std::vector<Double_t> > z;
   std::vector<Double_t> k;
   std::vector<Double_t> b;
   std::vector<Double_t> dk;
   std::vector<Double_t> db;
};

dEdxData dedx_data;

double fcr(double *x, double *par)
{
   double theta = x[0];
   double sz = par[0];
   double dedx = par[1];
   double demax = par[2]/par[0];
   double g0, g1, g2;
   return satfsx(dedx, theta, sz, demax, g0, g1, g2 );
}

double fcorr(const double *par)
{
   //calculate chisquare
   Double_t chi2 = 0;
   Double_t dchi2;
   for (size_t i=0;i<dedx_data.y.size(); i++) {
      double yx = 0;
      for (size_t j=0;j<dedx_data.z.size();j++)
         yx += dedx_data.z[j][i]*par[2*j] + par[2*j+1];
      dchi2 = (dedx_data.y[i]-yx)/dedx_data.dy[i];
      chi2 += dchi2*dchi2;
   }
   return chi2;
}

void dEdx0m_fit()
{
   const char *minName = "Minuit2";
   const char *algoName = "";

   double chi2;
   size_t ndf;
   const double *p;
   const double *dp;
   {
      ROOT::Math::Minimizer* min =
         ROOT::Math::Factory::CreateMinimizer(minName, algoName);
   
      // set tolerance , etc...
      min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      min->SetMaxIterations(10000);  // for GSL
      min->SetTolerance(0.00001);
      min->SetPrintLevel(1);
       
      // create funciton wrapper for minmizer
      // a IMultiGenFunction type
      size_t np = dedx_data.z.size();
      ROOT::Math::Functor f(&fcorr,2*np);   
      min->SetFunction(f);
       
      // Set the free variables to be minimized!
      for(size_t i=0;i<np;i++) {
         min->SetLimitedVariable(2*i,Form("k_%d",i),-5000/np, 10, -10000,10000);
         min->SetLimitedVariable(2*i+1,Form("b_%d",i),9500/np, 10, -10000, 10000);
      }
      
      // do the minimization
      for(size_t i=0;i<np;i++) 
         min->SetFixedVariable(2*i,Form("k_%d",i),1);
      min->Minimize();
      p = min->X();
      dp = min->Errors();
   }
   {
      ROOT::Math::Minimizer* min =
         ROOT::Math::Factory::CreateMinimizer(minName, algoName);
   
      // set tolerance , etc...
      min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      min->SetMaxIterations(10000);  // for GSL
      min->SetTolerance(0.00001);
      min->SetPrintLevel(1);
       
      // create funciton wrapper for minmizer
      // a IMultiGenFunction type
      size_t np = dedx_data.z.size();
      ROOT::Math::Functor f(&fcorr,2*np);   
      min->SetFunction(f);
       
      // Set the free variables to be minimized!
      for(size_t i=0;i<np;i++) {
         min->SetLimitedVariable(2*i,Form("k_%d",i),-5000/np, 10, -10000,10000);
         min->SetLimitedVariable(2*i+1,Form("b_%d",i),p[2*i+1],dp[2*i+1], -10000, 10000);
      }
      
      // do the minimization
      for(size_t i=0;i<np;i++) 
         min->SetFixedVariable(2*i+1,Form("b_%d",i),p[2*i+1]);
      min->Minimize();
      p = min->X();
      dp = min->Errors();
   }
   {
      ROOT::Math::Minimizer* min =
         ROOT::Math::Factory::CreateMinimizer(minName, algoName);
   
      // set tolerance , etc...
      min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      min->SetMaxIterations(10000);  // for GSL
      min->SetTolerance(0.00001);
      min->SetPrintLevel(1);
       
      // create funciton wrapper for minmizer
      // a IMultiGenFunction type
      size_t np = dedx_data.z.size();
      ROOT::Math::Functor f(&fcorr,2*np);   
      min->SetFunction(f);
       
      // Set the free variables to be minimized!
      for(size_t i=0;i<np;i++) {
         min->SetLimitedVariable(2*i,Form("k_%d",i),p[2*i],dp[2*i], -20000,20000);
         min->SetLimitedVariable(2*i+1,Form("b_%d",i),p[2*i+1],dp[2*i+1], -20000, 20000);
      }
      
      // do the minimization
      min->Minimize();
      p = min->X();
      dp = min->Errors();
      dedx_data.k.resize(0);
      dedx_data.dk.resize(0);
      dedx_data.b.resize(0);
      dedx_data.db.resize(0);
      for(size_t i=0;i<np;i++) {
         dedx_data.k.push_back(p[2*i]);
         dedx_data.dk.push_back(dp[2*i]);
         dedx_data.b.push_back(p[2*i+1]);
         dedx_data.db.push_back(dp[2*i+1]);
      }
   }
}

std::vector<double> fX;
std::vector<double> fY;
std::vector<double> fVxx;
std::vector<double> fVxy;
std::vector<double> fVyy;

double f_sz_dEdxm(const double *par)
{
   //calculate chisquare
   Double_t chi2 = 0;
   for (size_t i=0;i<fX.size(); i++) {
      Double_t dx = par[0] - fX[i];
      Double_t dy = par[1] - fY[i];
      chi2 += dx*dx*fVxx[i] + 2*dx*dy*fVxy[i] + dy*dy*fVyy[i];
   }
   return chi2;
}

void sz_dEdxm_fit(const std::vector<double> &x,
                  const std::vector<double> &y,
                  const std::vector<double> &Vxx,
                  const std::vector<double> &Vxy,
                  const std::vector<double> &Vyy)
{
   const char *minName = "Minuit2";
   const char *algoName = "";

   double chi2;
   size_t ndf;
   const double *p;
   const double *dp;
   {
      ROOT::Math::Minimizer* min =
         ROOT::Math::Factory::CreateMinimizer(minName, algoName);
   
      // set tolerance , etc...
      min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      min->SetMaxIterations(10000);  // for GSL
      min->SetTolerance(0.00001);
      min->SetPrintLevel(1);
       
      size_t n = x.size();
      fX.resize(n,0);
      fY.resize(n,0);
      fVxx.resize(n,0);
      fVxy.resize(n,0);
      fVyy.resize(n,0);
      for(size_t i=0;i<n;i++) {
         fX[i] = x[i];
         fY[i] = y[i];
         fVxx[i] = Vxx[i];
         fVxy[i] = Vxy[i];
         fVyy[i] = Vyy[i];
      }
      
      std::cout << "size = " << n << std::endl;
      // create funciton wrapper for minmizer
      // a IMultiGenFunction type
      ROOT::Math::Functor f(&f_sz_dEdxm,2);   
      min->SetFunction(f);
       
      // Set the free variables to be minimized!
      min->SetLimitedVariable(0,"sz",0.05, 0.01, 0, 0.2);
      min->SetLimitedVariable(1,"dedxm",200, 10, 10, 1000);
      
      // do the minimization
      min->Minimize();
      p = min->X();
      dp = min->Errors();
   }
}


void plot(size_t np = 0, size_t ip = 0)
{
   //TFile *f = new TFile("xxx.root");
   //TFile *f = new TFile("mhad2012.root");
   TFile *f = new TFile("mhad2017y.root");
   TIter next_run(f->GetListOfKeys());
   TKey *key_run;
   TString key_run_name("TDirectoryFile");
   TString key_name("TFitResult");
   TString key_name_fun("f_fits_");
   std::vector<double> bm;
   std::vector<double> rd;
   std::vector<double> td;
   std::vector<double> tdm;
   std::vector<double> rn;
   std::vector<double> tm;
   std::vector<size_t> pind;
   for(size_t ind=192;ind<208;ind++)
      pind.push_back(ind);
   for(size_t ind=69;ind<78;ind+=2)
      pind.push_back(ind);
   std::vector< std::vector<double> > dbp(pind.size(),std::vector<double>());
   std::vector< std::vector<double> > szm(10,std::vector<double>());
   std::vector< std::vector<double> > dszm(10,std::vector<double>());
   std::vector< std::vector<double> > dedx0m(10,std::vector<double>());
   std::vector< std::vector<double> > ddedx0m(10,std::vector<double>());
   std::vector< std::vector<double> > dedxmm(10,std::vector<double>());
   std::vector< std::vector<double> > ddedxmm(10,std::vector<double>());
   std::vector< std::vector<double> > szd(10,std::vector<double>());
   std::vector< std::vector<double> > dszd(10,std::vector<double>());
   std::vector< std::vector<double> > dedx0d(10,std::vector<double>());
   std::vector< std::vector<double> > ddedx0d(10,std::vector<double>());
   std::vector< std::vector<double> > dedxmd(10,std::vector<double>());
   std::vector< std::vector<double> > ddedxmd(10,std::vector<double>());
   std::vector< std::vector<double> > sloped(10,std::vector<double>());
   std::vector< std::vector<double> > dsloped(10,std::vector<double>());
   std::vector< std::vector<double> > szb(10,std::vector<double>());
   std::vector< std::vector<double> > dszb(10,std::vector<double>());
   std::vector< std::vector<double> > dedx0b(10,std::vector<double>());
   std::vector< std::vector<double> > ddedx0b(10,std::vector<double>());
   std::vector< std::vector<double> > dedxmb(10,std::vector<double>());
   std::vector< std::vector<double> > ddedxmb(10,std::vector<double>());
   std::vector< std::vector<double> > ne(10,std::vector<double>());
   std::vector< std::vector<TMatrixDSym> > Dmatrix(10,std::vector<TMatrixDSym>());
   while ((key_run=(TKey*)next_run())) {
      if (!key_run_name.BeginsWith(key_run->GetClassName())) continue;
      TDirectory *d = (TDirectory*)f->Get(key_run->GetName());
      TIter next(d->GetListOfKeys());
      TKey *key;
      while ((key=(TKey*)next())) {
         TString name(key->GetName());
         if (name.BeginsWith("par_run")) {
            int run, layer;
            sscanf(key->GetName(),"par_run%08d_layer%d",&run,&layer);
            std::vector<Double_t> *v;
            d->GetObject(key->GetName(),v);
            szm[layer].push_back(v->at(0));
            dedx0m[layer].push_back(v->at(2));
            dedxmm[layer].push_back(v->at(1));
         }
         if (name.BeginsWith("dpar_run")) {
            int run, layer;
            sscanf(key->GetName(),"dpar_run%08d_layer%d",&run,&layer);
            std::vector<Double_t> *v;
            d->GetObject(key->GetName(),v);
            dszm[layer].push_back(v->at(0));
            ddedx0m[layer].push_back(v->at(2));
            ddedxmm[layer].push_back(v->at(1));
         }
         if (name.BeginsWith("TFitResult-beam")&&
             name.Contains(key_name_fun)) {
            int num, layer;
            sscanf(key->GetName(),"TFitResult-beam%d_layer%d",&num,&layer);
            TFitResult *r = (TFitResult*)d->Get(key->GetName());
            if (layer==0) bm.push_back(num/10.);
            szb[layer].push_back(r->Parameters().at(0));
            dszb[layer].push_back(r->Errors().at(0));
            dedx0b[layer].push_back(r->Parameters().at(1));
            ddedx0b[layer].push_back(r->Errors().at(1));
            dedxmb[layer].push_back(r->Parameters().at(2));
            ddedxmb[layer].push_back(r->Errors().at(2));
         }
         //if (name.BeginsWith("TFitResult-day")&&
         //if (name.BeginsWith("TFitResult-dcpagen")&&
         if (name.BeginsWith("TFitResult-pcgen")&&
             name.Contains(key_name_fun)) {
            int date, layer, dummy;
            //sscanf(key->GetName(),"TFitResult-day%d_layer%d",&date,&layer);
            //sscanf(key->GetName(),"TFitResult-dcpagen%d_layer%d",&date,&layer);
            sscanf(key->GetName(),"TFitResult-pcgen%d_layer%d",&date,&layer);
            TProfile *h = (TProfile*)d->Get(Form("pcgen%08d_layer%d",date,layer));
            TFitResult *r = (TFitResult*)d->Get(key->GetName());
            if (layer==0) rd.push_back(date);
            //if (layer==0) { td.push_back(fdatetotime(date)); tdm.push_back(td.back()+100);}
            if (layer==0) { td.push_back(date); tdm.push_back(td.back()+100);}
            szd[layer].push_back(r->Parameters().at(0));
            dszd[layer].push_back(r->Errors().at(0));
            dedx0d[layer].push_back(r->Parameters().at(1));
            ddedx0d[layer].push_back(r->Errors().at(1));
            dedxmd[layer].push_back(r->Parameters().at(2));
            ddedxmd[layer].push_back(r->Errors().at(2));
            sloped[layer].push_back(r->Parameters().at(3));
            dsloped[layer].push_back(r->Errors().at(3));
         }
         if (name.BeginsWith("TFitResult-run")&&
             name.Contains(key_name_fun)) {
            int run, layer;
            sscanf(key->GetName(),"TFitResult-run%d_layer%d",&run,&layer);
            //if (run<12125) continue;
            TFitResult *r = (TFitResult*)d->Get(key->GetName());
            //if (layer==0) {
            //   for(size_t i=0;i<pind.size();i++)
            //      dbp[i].push_back(fpar(run,pind[i]));
            //}
            if (layer==0) rn.push_back(run);
            if (layer==0) { t.push_back(ftime(run)); tm.push_back(t.back()+100);}
            //if (layer==0) t.push_back(run);
            sz[layer].push_back(r->Parameters().at(0));
            dsz[layer].push_back(r->Errors().at(0));
            dedx0[layer].push_back(r->Parameters().at(1));
            ddedx0[layer].push_back(r->Errors().at(1));
            dedxm[layer].push_back(r->Parameters().at(2));
            ddedxm[layer].push_back(r->Errors().at(2));
            Dmatrix[layer].push_back(r->GetCovarianceMatrix());
            //
            TProfile *h = (TProfile*)d->Get(Form("run%08d_layer%d",run,layer));
            ne[layer].push_back(h->GetEntries());
         }
      }
   }
   //{
   //   TFile *fout = TFile::Open("dbpars_2017.root","RECREATE");      
   //   for(size_t i=0;i<dbp.size();i++)
   //      fout->WriteObjectAny(&dbp[i],"std::vector<Double_t>",Form("dbp%d",i));
   //   fout->Close();
   //}
   {
      //TFile *fin = TFile::Open("dbpars.root","READ");
      TFile *fin = TFile::Open("dbpars_2017.root","READ");
      std::vector<Double_t> *v;
      for(size_t i=0;i<dbp.size();i++) {
         fin->GetObject(Form("dbp%d",i),v);
         for(std::vector<Double_t>::iterator it = v->begin(); it != v->end(); ++it) 
            dbp[i].push_back(*it);
      }
   }
   size_t layer = 0;
   std::cout << sz[layer].size() << std::endl;
   std::vector<double> dt(t.size(),0);
   std::vector<double> dtd(td.size(),0);
   std::vector<double> dtb(bm.size(),0);
   TGraphErrors *gr;
   TGraphErrors *gr1;
   if (np==-1000) {
      TFile *fout = TFile::Open("DBpars_2017_x.root","RECREATE");      
      {
         TGraphErrors *gr = new TGraphErrors(t.size(),t.data(),rn.data(),dt.data(),dt.data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name("run%time");
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }  
      for(size_t ip=0;ip<dbp.size();ip++) {
         TGraphErrors *gr = new TGraphErrors(t.size(),t.data(),dbp[ip].data(),dt.data(),dt.data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name(ip<16?
                      Form("PS#%02d sensor #%02d",ip,ip):
                      Form("HV%d-1 Drift Chamber VVI%d",ip-15,ip-15));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }      
      for(size_t i=0;i<sz.size();i++) {
         TGraphErrors *gr = new TGraphErrors(t.size(),t.data(),sz[i].data(),dt.data(),dsz[i].data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name(Form("#sigma_{z}: layer_#%d",i));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }      
      for(size_t i=0;i<sz.size();i++) {
         TString name(Form("#sigma_{z}: layer_#%d spect",i));
         TH1D *h = new TH1D(name,name,100,0,0.1);
         for(size_t j=0;j<sz[i].size();j++)
            h->Fill(sz[i][j]);
         h->Draw();
         h->GetXaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }      
      for(size_t i=0;i<dedx0.size();i++) {
         TGraphErrors *gr = new TGraphErrors(t.size(),t.data(),dedx0[i].data(),dt.data(),ddedx0[i].data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name(Form("dEdx_{0}: layer_#%d",i));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }
      for(size_t i=0;i<dedxm.size();i++) {
         TGraphErrors *gr = new TGraphErrors(t.size(),t.data(),dedxm[i].data(),dt.data(),ddedxm[i].data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name(Form("dEdx_{m}: layer_#%d",i));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }      
      for(size_t i=0;i<dedxm.size();i++) {
         TString name(Form("dEdx_{m}: layer_#%d spect",i));
         TH1D *h = new TH1D(name,name,100,0,300);
         for(size_t j=0;j<dedxm[i].size();j++)
            h->Fill(dedxm[i][j]);
         h->Draw();
         h->GetXaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }      
      for(size_t i=0;i<szb.size();i++) {
         TGraphErrors *gr = new TGraphErrors(bm.size(),bm.data(),szb[i].data(),dtb.data(),dszb[i].data());
         //gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name(Form("beam: #sigma_{z}: layer_#%d",i));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }      
      for(size_t i=0;i<dedx0b.size();i++) {
         TGraphErrors *gr = new TGraphErrors(bm.size(),bm.data(),dedx0b[i].data(),dtb.data(),ddedx0b[i].data());
         //gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name(Form("beam: dEdx_{0}: layer_#%d",i));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }
      for(size_t i=0;i<dedxmb.size();i++) {
         TGraphErrors *gr = new TGraphErrors(bm.size(),bm.data(),dedxmb[i].data(),dtb.data(),ddedxmb[i].data());
         //gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name(Form("beam: dEdx_{m}: layer_#%d",i));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }  
      std::vector<size_t> ts;
      ts.push_back(fdatetotime(20000101));
      ts.push_back(fdatetotime(20170320));
      //ts.push_back(fdatetotime(20170525));
      ts.push_back(fdatetotime(21000101));
      for(size_t it = 1; it<ts.size(); it++) {
         size_t ts1 = ts[it-1];
         size_t ts2 = ts[it];
         std::vector<double> szl;
         std::vector<double> dszl;
         std::vector<double> ls;
         std::vector<double> dls;
         std::vector<double> chi2;
         for(size_t i=0;i<szd.size();i++) {
            std::vector<double> x;
            std::vector<double> y;
            std::vector<double> dx;
            std::vector<double> dy;
            double mean = 0.04;
            double rms = 0;
            double ns = 10;
            TF1 *f0 = new TF1("p0","[0]");
            for (size_t s=0;s<3;s++) {
               x.resize(0);
               y.resize(0);
               dx.resize(0);
               dy.resize(0);
               for(size_t j=0;j<td.size();j++) {
                  if (td[j]<ts1 || td[j]>ts2) continue;
                  if (szd[i][j]<0.001) continue;
                  if (szd[i][j]>1.0) continue;
                  if (std::fabs(szd[i][j]-mean)>ns*dszd[i][j]) continue;
                  x.push_back(td[j]);
                  dx.push_back(0);
                  y.push_back(szd[i][j]);
                  dy.push_back(dszd[i][j]);
               }
               TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
               gr->Fit("p0");
               mean = f0->GetParameter(0);
               rms = f0->GetParError(0);
               std::cout << mean << " " << rms << " " << f0->GetNDF() << " " << f0->GetChisquare() << std::endl;
               ns = 3;
            }
            TGraphErrors *gr = new TGraphErrors(td.size(),td.data(),szd[i].data(),dtd.data(),dszd[i].data());
            gr->GetXaxis()->SetTimeDisplay(1);
            gr->SetMarkerSize(0.5);
            gr->SetMarkerStyle(20);
            gr->Draw("AP");
            TGraphErrors *gr0 = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
            gr0->GetXaxis()->SetTimeDisplay(1);
            gr0->SetMarkerColor(kBlue);
            gr0->SetLineColor(kBlue);
            gr0->SetMarkerSize(0.5);
            gr0->SetMarkerStyle(20);
            gr0->Draw("P");
            gr0->Fit("p0");
            szl.push_back(f0->GetParameter(0));
            dszl.push_back(f0->GetParError(0));
            ls.push_back((double)i);
            dls.push_back(0);
            chi2.push_back(f0->GetChisquare()/f0->GetNDF());
         ns = 100;
         gr->Fit("pol0");
         TF1 *pn = (TF1*)gr->GetFunction("pol0");
         for(size_t s=0;s<3;s++) {
            std::vector<double> x;
            std::vector<double> y;
            std::vector<double> dx;
            std::vector<double> dy;
            std::cout << ">>>" << std::endl;
            for(size_t j=0;j<td.size();j++) {
               if (td[j]<ts1 || td[j]>ts2) continue;
               double tx = (td[j]-fdatetotime(20170101))/86400;
               if (std::fabs(szd[i][j]-pn->Eval(tx))>ns*dszd[i][j]) continue;
               x.push_back(tx);
               dx.push_back(dtd[j]);
               y.push_back(szd[i][j]);
               dy.push_back(dszd[i][j]);
            }
            std::cout << x.size() << std::endl;
            TGraphErrors *grf = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
            //grf->SetMarkerSize(0.5);
            //grf->SetMarkerStyle(20);
            //grf->Draw("P");
            grf->Fit("pol1");
            pn = (TF1*)grf->GetFunction("pol1");
            ns = 3;
         }
         x.resize(0);;
         y.resize(0);
         dx.resize(0);;
         dy.resize(0);
         for(size_t j=0;j<td.size();j++) {
            if (td[j]<ts1 || td[j]>ts2) continue;
            x.push_back(td[j]);
            dx.push_back(dtd[j]);
            y.push_back(pn->Eval((td[j]-fdatetotime(20170101))/86400));
            dy.push_back(0);
         }
         TGraphErrors *grfx = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
         grfx->SetMarkerSize(0.5);
         grfx->SetMarkerStyle(20);
         grfx->SetMarkerColor(kRed);
         grfx->SetLineColor(kRed);
         grfx->Draw("PL");
            TString name(Form("day: #sigma_{z}: interval_%d layer_#%d",it,i));
            gr->GetYaxis()->SetTitle(name);
            TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
            c->SetName(name);
            c->SetTitle(name);
            c->Write();
         }      
         {
            TGraphErrors *gr = new TGraphErrors(ls.size(),ls.data(),szl.data(),dls.data(),dszl.data());
            gr->SetMarkerSize(0.5);
            gr->SetMarkerStyle(20);
            gr->Draw("AP");
            TString name(Form("#sigma_{z} vs layer: interval_%d",it));
            gr->GetYaxis()->SetTitle(name);
            TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
            c->SetName(name);
            c->SetTitle(name);
            c->Write();
         }
         {
            TGraphErrors *gr = new TGraphErrors(ls.size(),ls.data(),chi2.data(),dls.data(),dls.data());
            gr->SetMarkerSize(0.5);
            gr->SetMarkerStyle(20);
            gr->Draw("AP");
            TString name(Form("#sigma_{z}: #Chi^{2} vs layer: interval_%d",it));
            gr->GetYaxis()->SetTitle(name);
            TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
            c->SetName(name);
            c->SetTitle(name);
            c->Write();
         }
      }
      for(size_t i=0;i<dedx0d.size();i++) {
         TGraphErrors *gr = new TGraphErrors(td.size(),td.data(),dedx0d[i].data(),dtd.data(),ddedx0d[i].data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name(Form("day: dEdx_{0}: layer_#%d",i));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }
      for(size_t i=0;i<dedxmd.size();i++) {
         TGraphErrors *gr = new TGraphErrors(td.size(),td.data(),dedxmd[i].data(),dtd.data(),ddedxmd[i].data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         double ns = 100;
         gr->Fit("pol0");
         TF1 *pn = (TF1*)gr->GetFunction("pol0");
         for(size_t s=0;s<3;s++) {
            std::vector<double> x;
            std::vector<double> y;
            std::vector<double> dx;
            std::vector<double> dy;
            for(size_t j=0;j<td.size();j++) {
               double tx = (td[j]-fdatetotime(20170101))/86400;
               if (std::fabs(dedxmd[i][j]-pn->Eval(tx))>ns*ddedxmd[i][j]) continue;
               x.push_back(tx);
               dx.push_back(dtd[j]);
               y.push_back(dedxmd[i][j]);
               dy.push_back(ddedxmd[i][j]);
            }
            TGraphErrors *grf = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
            //grf->SetMarkerSize(0.5);
            //grf->SetMarkerStyle(20);
            //grf->Draw("P");
            grf->Fit("pol6");
            pn = (TF1*)grf->GetFunction("pol6");
            ns = 3;
         }
         std::vector<double> x;
         std::vector<double> y;
         std::vector<double> dx;
         std::vector<double> dy;
         for(size_t j=0;j<td.size();j++) {
            x.push_back(td[j]);
            dx.push_back(dtd[j]);
            y.push_back(pn->Eval((td[j]-fdatetotime(20170101))/86400));
            dy.push_back(0);
         }
         TGraphErrors *grfx = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
         grfx->SetMarkerSize(0.5);
         grfx->SetMarkerStyle(20);
         grfx->SetMarkerColor(kRed);
         grfx->SetLineColor(kRed);
         grfx->Draw("PL");
         TString name(Form("day: dEdx_{m}: layer_#%d",i));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }  
      for(size_t i=0;i<sloped.size();i++) {
         TGraphErrors *gr = new TGraphErrors(td.size(),td.data(),sloped[i].data(),dtd.data(),dsloped[i].data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         double ns = 100;
         gr->Fit("pol0");
         TF1 *pn = (TF1*)gr->GetFunction("pol0");
         for(size_t s=0;s<3;s++) {
            std::vector<double> x;
            std::vector<double> y;
            std::vector<double> dx;
            std::vector<double> dy;
            for(size_t j=0;j<td.size();j++) {
               double tx = (td[j]-fdatetotime(20170101))/86400;
               if (std::fabs(sloped[i][j]-pn->Eval(tx))>ns*dsloped[i][j]) continue;
               x.push_back(tx);
               dx.push_back(dtd[j]);
               y.push_back(sloped[i][j]);
               dy.push_back(dsloped[i][j]);
            }
            TGraphErrors *grf = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
            //grf->SetMarkerSize(0.5);
            //grf->SetMarkerStyle(20);
            //grf->Draw("P");
            grf->Fit("pol6");
            pn = (TF1*)grf->GetFunction("pol6");
            ns = 3-s;
         }
         std::vector<double> x;
         std::vector<double> y;
         std::vector<double> dx;
         std::vector<double> dy;
         for(size_t j=0;j<td.size();j++) {
            x.push_back(td[j]);
            dx.push_back(dtd[j]);
            y.push_back(pn->Eval((td[j]-fdatetotime(20170101))/86400));
            dy.push_back(0);
         }
         TGraphErrors *grfx = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
         grfx->SetMarkerSize(0.5);
         grfx->SetMarkerStyle(20);
         grfx->SetMarkerColor(kRed);
         grfx->SetLineColor(kRed);
         grfx->Draw("PL");
         TString name(Form("day: slope: layer_#%d",i));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }  
      std::vector<double> l;
      std::vector<double> dl;
      std::vector<double> k;
      std::vector<double> dk;
      for(size_t layer=2;layer<10;layer++) {
         std::vector<double> x;
         std::vector<double> dx;
         std::vector<double> y;
         std::vector<double> dy;
         for(size_t i=1;i<dedxmd[1].size();i++) {
            x.push_back(dedxmd[1][i]);
            dx.push_back(ddedxmd[1][i]);
            y.push_back(dedxmd[layer][i]);
            dy.push_back(ddedxmd[layer][i]);
         }
         TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TF1 *f1 = new TF1("f1","[0]*x",0,1000);
         f1->SetParameter(0,1);
         gr->Fit("f1");
         l.push_back(layer);
         dl.push_back(0);
         k.push_back(f1->GetParameter(0));
         dk.push_back(f1->GetParError(0));
         TString name(Form("day: dEdx_{m}[layer_1] vs dEdx_{m}[layer_#%d]",layer));
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }
      {
         TGraphErrors *gr = new TGraphErrors(l.size(),l.data(),k.data(),dl.data(),dk.data());
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TString name("k vs layer");
         gr->GetYaxis()->SetTitle(name);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();
      }
      std::vector<double> dr;
      dr.push_back(fdatetotime(20170319));
      dr.push_back(fdatetotime(20170423));
      dr.push_back(fdatetotime(20170512));
      dr.push_back(fdatetotime(20170525));
      dr.push_back(fdatetotime(20170618));
      size_t npart = dr.size()+1;
      for(size_t layer=0;layer<10;layer++) {
         std::vector< std::vector<double> > szc(npart,std::vector<double>());
         std::vector< std::vector<double> > dedxmc(npart,std::vector<double>());
         std::vector< std::vector<double> > V11(npart,std::vector<double>());
         std::vector< std::vector<double> > V12(npart,std::vector<double>());
         std::vector< std::vector<double> > V22(npart,std::vector<double>());
         for(size_t i=0;i<sz[layer].size();i++) {
            size_t jm = 0;
            for(size_t j=0;j<dr.size();j++)
               if (t[i]>dr[j]) jm++;
            if (sz[layer][i]>0.19) continue;
            Double_t D = Dmatrix[layer][i](0,0)*Dmatrix[layer][i](1,1)-Dmatrix[layer][i](1,0)*Dmatrix[layer][i](1,0);
            if (D<1e-3) continue;
            szc[jm].push_back(sz[layer][i]);
            dedxmc[jm].push_back(dedxm[layer][i]);
            V11[jm].push_back(Dmatrix[layer][i](1,1)/D);
            V12[jm].push_back(-Dmatrix[layer][i](1,0)/D);
            V22[jm].push_back(Dmatrix[layer][i](0,0)/D);
         }
         //
         for(size_t i=0;i<npart;i++) {
            sz_dEdxm_fit(szc[i],dedxmc[i],V11[i],V12[i],V22[i]);
         }
         //
      }
      for(size_t layer=0;layer<10;layer++) {
         std::vector< std::vector<double> > szc(npart,std::vector<double>());
         std::vector< std::vector<double> > dszc(npart,std::vector<double>());
         std::vector< std::vector<double> > dedxmc(npart,std::vector<double>());
         std::vector< std::vector<double> > ddedxmc(npart,std::vector<double>());
         for(size_t i=0;i<szd[layer].size();i++) {
            size_t jm = 0;
            for(size_t j=0;j<dr.size();j++)
               if (td[i]>dr[j]) jm++;
            szc[jm].push_back(szd[layer][i]);
            dszc[jm].push_back(0*dszd[layer][i]);
            dedxmc[jm].push_back(dedxmd[layer][i]);
            ddedxmc[jm].push_back(0*ddedxmd[layer][i]);
         }
         //
         std::vector<TGraphErrors*> gr;
         gr.push_back(new TGraphErrors(szd[layer].size(),szd[layer].data(),dedxmd[layer].data(),dtd.data(),dtd.data()));
         //gr.back()->GetXaxis()->SetTimeDisplay(1);
         gr.back()->SetMarkerSize(1);
         gr.back()->SetMarkerStyle(20);
         gr.back()->Draw("AP");
         for(size_t i=0;i<szc.size();i++) {
            gr.push_back(new TGraphErrors(szc[i].size(),szc[i].data(),dedxmc[i].data(),dszc[i].data(),ddedxmc[i].data()));
            //gr.back()->GetXaxis()->SetTimeDisplay(1);
            gr.back()->SetMarkerSize(1);
            gr.back()->SetMarkerStyle(20);
            gr.back()->SetMarkerColor(i+1);
            gr.back()->SetLineColor(i+1);
            gr.back()->Draw("P");
         }
         TString name(Form("day: dEdx_{m} vs sz: layer_#%d",layer));
         gr.back()->GetYaxis()->SetTitle("dEdx_{m}");
         gr.back()->GetYaxis()->SetTitle("sz");
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->SetName(name);
         c->SetTitle(name);
         c->Write();         
      }
      fout->Close();
   }
   if (np==2000) {
      std::vector<size_t> p;
      //for(size_t i=1;i<4;i++) p.push_back(i);
      p.push_back(1);
      //p.push_back(3);
      for(size_t i=0;i<p.size();i++)
         dedx_data.z.push_back(std::vector<double>());
      for(size_t i=0;i<dedx0m[layer].size();i++) {
         if (dedx0m[layer][i]>500&&ddedx0m[layer][i]>1&&
             rn[i]>12125&&rn[i]<12225) {
            dedx_data.x.push_back(t[i]);
            dedx_data.y.push_back(dedx0m[layer][i]);
            dedx_data.dy.push_back(ddedx0m[layer][i]);
            for(size_t j=0;j<p.size();j++)
               dedx_data.z[j].push_back(dbp[p[j]][i]);
         }
      }
      dEdx0m_fit();
      std::vector<double> dx(dedx_data.x.size(),0);
      TGraphErrors *gr = new TGraphErrors(dedx_data.x.size(),
                                          dedx_data.x.data(),
                                          dedx_data.y.data(),
                                          dx.data(),
                                          dedx_data.dy.data());
      gr->SetMarkerSize(0.5);
      gr->SetMarkerStyle(20);
      gr->Draw("AP");
      for(size_t j=0;j<dedx_data.z.size();j++) {
         std::vector<double> z;
         std::cout << j << " " << dedx_data.k[j] << " " << dedx_data.b[j] << std::endl;
         for(size_t i=0;i<dedx_data.z[j].size();i++)
            z.push_back(dedx_data.z[j][i]*dedx_data.k[j]+dedx_data.b[j]);
         TGraphErrors *gr = new TGraphErrors(dedx_data.x.size(),
                                             dedx_data.x.data(),
                                             z.data(),
                                             dx.data(),
                                             dx.data());
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->SetMarkerColor(kRed);
         gr->Draw("P");
      }
/*      {
         std::vector<double> z;
         std::vector<double> dz;
         size_t j = 0;
         for(size_t i=0;i<dedx_data.z[j].size();i++) {
            double zi = dedx_data.z[j][i]*dedx_data.k[j]+dedx_data.b[j];
            z.push_back(dedx_data.y[i]/zi);
            dz.push_back(dedx_data.dy[i]/zi);
         }
         TGraphErrors *gr = new TGraphErrors(dedx_data.x.size(),
                                             dedx_data.x.data(),
                                             z.data(),
                                             dx.data(),
                                             dz.data());
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->SetMarkerColor(kRed);
         gr->Draw("AP");
      }*/
   }
   if (np==1000) {
      for(size_t i=0;i<dedx0m[layer].size();i++) {
         if (TMath::IsNaN(dedx0m[layer][i]))
            std::cout << rn[i] << " " << layer << " " << dedx0m[layer][i] << " " << szm[layer][i] << " " << ddedxmm[layer][i] 
            << " " << dedx0[layer][i] << " " << sz[layer][i] << " " << ddedxm[layer][i] << std::endl;
      }
   }
   if (np==200) {
      std::vector<double> q;
      std::vector<double> p;
      size_t ip = 1;
      double pm = 0;
      size_t im = 0;
      double xm = 1.77;
      for(double x0=0;x0<1000;x0+=10) {
         TH1D* h = new TH1D("h","h",400,0,4000);
         for(size_t i=0;i<dedx0m[layer].size();i++) {
            double x = dbp[ip][i];
            double y = dedx0m[layer][i];
            if (x&&!TMath::IsNaN(y))
                h->Fill((y-x0)*x);
                //h->Fill(y*(x0-xm)/(x0-x));
         }
         std::vector<size_t> hi;
         for(size_t i=1;i<h->GetNbinsX()-1;i++)
            if (h->GetBinContent(i-1)<h->GetBinContent(i)&&
                h->GetBinContent(i)>h->GetBinContent(i+1))
               hi.push_back(h->GetBinContent(i));
         std::sort(hi.begin(), hi.end());
         size_t max = 0;
         for(int i=hi.size()-5;i<(int)hi.size();i++)
            if (i>=0) max += hi[i];
         q.push_back(x0);
         p.push_back(max);
         if (p.back()>pm) {pm=p.back(); im=p.size()-1;}
         h->Draw();
         {
            TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
            c->Update();
            size_t xxx;
            std::cin >> xxx;
         }
      }
      TGraph *gr = new TGraph(q.size(),q.data(),p.data());
      gr->Draw("APL");
      std::cout << q[im] << " " << p[im] << std::endl;
      double x0 = q[im];
      {
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         size_t xxx;
         std::cin >> xxx;
      }
      //
      TH1D* h = new TH1D("h","h",1000,0,4000);
      for(size_t i=0;i<dedx0m[layer].size();i++) {
         double x = dbp[ip][i];
         double y = dedx0m[layer][i];
         if (x&&!TMath::IsNaN(y))
            h->Fill(y*(x0-xm)/(x0-x));
      }
      h->Draw();
         {
            TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
            c->Update();
            size_t xxx;
            std::cin >> xxx;
         }
      std::vector< std::pair<size_t,double> > xi;
      for(size_t i=1;i<h->GetNbinsX()-1;i++)
         if (h->GetBinContent(i-1)<h->GetBinContent(i)&&
             h->GetBinContent(i)>h->GetBinContent(i+1))
            xi.push_back(std::make_pair(h->GetBinContent(i),
                                        h->GetBinCenter(i)));
      sort(xi.begin(), xi.end(), sortby);
      {
         std::vector<double> x;
         std::vector<double> dx;
         std::vector<double> y;
         std::vector<double> dy;
         for(size_t i=0;i<dedx0m[layer].size();i++) {
            double xi = dbp[ip][i];
            double yi = dedx0m[layer][i];
            if (xi&&!TMath::IsNaN(yi)) {
               double k = (x0-xm)/(x0-xi);
               x.push_back(xi);
               dx.push_back(0);
               y.push_back(yi);
               dy.push_back(ddedx0m[layer][i]);
            }
         }
         TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         for(int i=xi.size()-5;i<(int)xi.size();i++)
            if (i>=0) {
               std::cout << xi[i].second << std::endl;
               double x1 = 1.73;
               double y1 = xi[i].second*(x1-x0)/(xm-x0);
               double x2 = 1.81;
               double y2 = xi[i].second*(x2-x0)/(xm-x0);
               TLine *l = new TLine();
               l->DrawLine(x1,y1,x2,y2);
            }
         {
            TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
            c->Update();
            size_t xxx;
            std::cin >> xxx;
         }
      }
      gr1 = new TGraphErrors(tm.size(),tm.data(),dedx0m[layer].data(),dt.data(),ddedx0m[layer].data());
      gr1->GetXaxis()->SetTimeDisplay(1);
      gr1->SetMarkerSize(0.5);
      gr1->SetMarkerStyle(20);
      gr1->Draw("AP");
      {
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         size_t xxx;
         std::cin >> xxx;
      }
      //
      std::vector<double> s;
      std::vector<double> ds;
      std::vector<double> y;
      std::vector<double> dy;
      for(size_t i=0;i<dedx0m[layer].size();i++) {
         double xi = dbp[1][i];
         double yi = dedx0m[layer][i];
         if (xi&&!TMath::IsNaN(yi)) {
            double k = (x0-xm)/(x0-xi);
            s.push_back(t[i]);
            ds.push_back(dt[i]);
            y.push_back(yi*k);
            dy.push_back(ddedx0m[layer][i]*k);
         }
      }
      gr1 = new TGraphErrors(tm.size(),tm.data(),dedx0m[layer].data(),dt.data(),ddedx0m[layer].data());
      gr1->GetXaxis()->SetTimeDisplay(1);
      gr1->SetMarkerSize(0.5);
      gr1->SetMarkerStyle(20);
      gr1->Draw("AP");
      gr = new TGraphErrors(s.size(),s.data(),y.data(),ds.data(),dy.data());
      gr->GetXaxis()->SetTimeDisplay(1);
      gr->SetMarkerSize(0.5);
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(kRed);
      gr->SetLineColor(kRed);
      gr->Draw("P");
   }
   if (np==100) {
      size_t j=1;
      //for(size_t j=0;j<dbp.size();j++) {
         std::vector<double> x;
         std::vector<double> dx;
         std::vector<double> y;
         std::vector<double> dy;
         for(size_t i=0;i<dedx0m[layer].size();i++) {
         //for(size_t i=0;i<100;i++) {
            if (dbp[j][i]) {
               x.push_back(dbp[j][i]);
               dx.push_back(0);
               //y.push_back(szm[layer][i]);
               //dy.push_back(dszm[layer][i]);
               //y.push_back(dedxmm[layer][i]);
               //dy.push_back(ddedxmm[layer][i]);
               y.push_back(dedx0m[layer][i]);
               dy.push_back(ddedx0m[layer][i]);
            }
         }
         gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
         //gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         //size_t xxx;
         //std::cin >> xxx;
      //}
   }
   if (np==101) {
      double tmin = t.front();
      double tmax = t.back();
      size_t nt = 10;
      double dt = (tmax-tmin)/(nt-1);
      for(size_t j=0;j<nt;j++) {
         std::vector<double> x;
         std::vector<double> dx;
         std::vector<double> y;
         std::vector<double> dy;
         for(size_t i=0;i<dedx0m[layer].size();i++) {
            if (dbp[1][i]&&
                (t[i]>tmin+j*dt)&&
                (t[i]<tmin+(j+1)*dt)) {
               x.push_back(dbp[1][i]);
               dx.push_back(0);
               y.push_back(dedx0m[layer][i]);
               dy.push_back(ddedx0m[layer][i]);
            }
         }
         gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
         //gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("APL");
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         size_t xxx;
         std::cin >> xxx;
      }
   }
   if (np==-1) {
      std::vector<size_t> ii;
      ii.push_back(1);
      ii.push_back(3);
      ii.push_back(4);
      ii.push_back(9);
      ii.push_back(11);
      ii.push_back(14);
      ii.push_back(15);
      for(size_t i=0;i<ii.size();i++) {
         ip = ii[i];
         TGraphErrors *gr = new TGraphErrors(t.size(),t.data(),dbp[ip].data(),dt.data(),dt.data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->SetMarkerColor(i);
         gr->SetLineColor(i);
         if (i) gr->Draw("P");
         else gr->Draw("AP");
      }
      double td = 86400*(int)t.front()/86400;
      TLine *l = new TLine();
      while(td<t.back()) {
         l->DrawLine(td,0,td,2);
         td += 86400;
      }
   }
   if (np==0) {
      std::vector<double> l;
      std::vector<double> dl;
      std::vector<double> szl;
      std::vector<double> dszl;
      std::vector<TH1D*> h;
      for(size_t layer=0;layer<10;layer++) {
         gr = new TGraphErrors(t.size(),t.data(),sz[layer].data(),dt.data(),dsz[layer].data());
         gr->GetXaxis()->SetTimeDisplay(1);
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->Draw("AP");
         gr1 = new TGraphErrors(tm.size(),tm.data(),szm[layer].data(),dt.data(),dszm[layer].data());
         gr1->GetXaxis()->SetTimeDisplay(1);
         gr1->SetMarkerSize(0.5);
         gr1->SetMarkerStyle(20);
         gr1->SetMarkerColor(kRed);
         gr1->SetLineColor(kRed);
         gr1->Draw("P");
         //gr1->Fit("pol0");
         //
         //TF1 *f = gr->GetFunction("pol0");
         h.push_back(new TH1D(Form("h%d",layer),Form("h%d",layer),100,0,0.1));
         for(size_t i=0;i<szm[layer].size();i++)
            h.back()->Fill(szm[layer][i]);
         h.back()->Draw();
         h.back()->Fit("gaus");
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         size_t xxx;
         std::cin >> xxx;
         l.push_back(layer);
         dl.push_back(0);
         TF1 *f = h.back()->GetFunction("gaus");
         szl.push_back(f->GetParameter(1));
         dszl.push_back(f->GetParError(1));
      }
      gr = new TGraphErrors(l.size(),l.data(),szl.data(),dl.data(),dszl.data());
      gr->SetMarkerSize(0.5);
      gr->SetMarkerStyle(20);
      gr->Draw("AP");
      //h[0]->Draw();
      //for(size_t layer=1;layer<10;layer++) {
      //   h[layer]->Draw("same");
      //   h[layer]->SetLineColor(layer);
      //}
   }
   if (np==1) {
      gr = new TGraphErrors(t.size(),t.data(),dedxm[layer].data(),dt.data(),ddedxm[layer].data());
      gr->GetXaxis()->SetTimeDisplay(1);
      gr->SetMarkerSize(0.5);
      gr->SetMarkerStyle(20);
      gr->Draw("AP");
   }
   if (np==2) {
      gr = new TGraphErrors(t.size(),t.data(),dedx0[layer].data(),dt.data(),ddedx0[layer].data());
      gr->GetXaxis()->SetTimeDisplay(1);
      gr->SetMarkerSize(0.5);
      gr->SetMarkerStyle(20);
      gr->Draw("AP");
      gr1 = new TGraphErrors(tm.size(),tm.data(),dedx0m[layer].data(),dt.data(),ddedx0m[layer].data());
      gr1->GetXaxis()->SetTimeDisplay(1);
      gr1->SetMarkerSize(0.5);
      gr1->SetMarkerStyle(20);
      gr1->SetMarkerColor(kRed);
      gr1->SetLineColor(kRed);
      gr1->Draw("P");
      for(size_t ip=0;ip<0;ip++) {
         std::vector<double> y(dbp[ip].size());
         for(size_t i=0;i<dbp[ip].size();i++)
            y[i] = 5000*(dbp[ip][i]-1.7)+700;
         TGraphErrors *gr2 = new TGraphErrors(t.size(),t.data(),y.data(),dt.data(),dt.data());
         gr2->GetXaxis()->SetTimeDisplay(1);
         gr2->SetMarkerSize(0.5);
         gr2->SetMarkerStyle(20);
         gr2->SetMarkerColor(ip%8);
         gr2->SetLineColor(ip%8);
         gr2->Draw("P");
      }
   }
   if (np==10) { 
      std::vector<double> VI11;
      std::vector<double> VI12;
      std::vector<double> VI22;
      std::vector<double> M1;
      std::vector<double> M2;
      for(size_t flayer=0;flayer<10;flayer++) {
         for(size_t c=0;c<2;c++) {
            double szm = 0.05;
            double dszm = 0.1;
            double alpham = 200;
            double dalpham = 100;
            std::vector<double> x;
            std::vector<double> dx;
            std::vector<double> y;
            std::vector<double> dy;
            const TMatrixD* m;
            const TVectorD* mv;
            double V11;
            double V12;
            double V22;
            double m1;
            double m2;
            for(size_t s=0;s<2;s++) {
               TPrincipal p(2,"ND");
               Double_t data[2];
               x.resize(0);
               dx.resize(0);
               y.resize(0);
               dy.resize(0);
               for(size_t i=0;i<sz[flayer].size();i++) {
                  if (//(c ? rn[i]>=12125 : rn[i]<12125)&&
                      sz[flayer][i]>0.01&&sz[flayer][i]<0.099&&
                      dedxm[flayer][i]>10&&dedxm[flayer][i]<990&&
                      dsz[flayer][i]>0.001&&ddedxm[flayer][i]>1&&
                      std::fabs(sz[flayer][i]-szm)<3*dszm&&
                      std::fabs(dedxm[flayer][i]-alpham)<3*dalpham) {
                     x.push_back(sz[flayer][i]);
                     dx.push_back(dsz[flayer][i]);
                     //y.push_back(dedxm[flayer][i]*sz[flayer][i]);
                     y.push_back(dedxm[flayer][i]);
                     dy.push_back(ddedxm[flayer][i]);
                     data[0] = sz[flayer][i];
                     data[1] = dedxm[flayer][i];
                     p.AddRow(data);
                  }
               }
               m = p.GetCovarianceMatrix();
               mv = p.GetMeanValues();
               V11 = (*m)(0,0);
               V12 = (*m)(1,0);
               V22 = (*m)(1,1);
               m1 = (*mv)(0);
               m2 = (*mv)(1);
               dszm = std::sqrt((*m)(0,0));
               dalpham = std::sqrt((*m)(1,1));
               szm = (*mv)(0);
               alpham = (*mv)(1);
               std::cout << szm << " " << dszm << " " << alpham << " " << dalpham << std::endl;
            }
            double D = V11*V22-V12*V12;
            VI11.push_back(V22/D);
            VI12.push_back(-V12/D);
            VI22.push_back(V11/D);
            M1.push_back(m1);
            M2.push_back(m2);
            //
            std::vector<double> d(x.size(),0);
            gr = new TGraphErrors(x.size(),x.data(),y.data(),d.data(),d.data());
            gr->SetMarkerSize(0.5);
            gr->SetMarkerStyle(20);
            gr->Draw("AP");
            gr1 = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
            TF1 *f1 = new TF1("f1",p1p2,0.01,0.1,1);
            f1->SetParName(0,"amp");
            f1->SetParameter(0,200);
            f1->Draw("same");
            gr1->Fit(f1,"S","",0.01,0.95);
            gr1->Fit("pol1");
            TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
            c->Update();
            //const TMatrixD* m = p.GetCovarianceMatrix();
            //m->Print();
            //const TVectorD* mv = p.GetMeanValues();
            //mv->Print();
            size_t xxx;
            std::cin >> xxx;
         }
      }
      TFile *fout = TFile::Open("Error_matix_2017.root","RECREATE");      
      fout->WriteObjectAny(&VI11,"std::vector<Double_t>","VI11");
      fout->WriteObjectAny(&VI12,"std::vector<Double_t>","VI12");
      fout->WriteObjectAny(&VI22,"std::vector<Double_t>","VI22");
      fout->WriteObjectAny(&M1,"std::vector<Double_t>","M1");
      fout->WriteObjectAny(&M2,"std::vector<Double_t>","M2");
      fout->Close();
   }
   //
   if (np==20) {
      TH1D* h_dsz = new TH1D(Form("h_sz%d",layer),Form("h_sz%d",layer),100,0,0.1);
      for(size_t j=0;j<dsz[layer].size();j++)
         h_dsz->Fill(sz[layer][j]);
      h_dsz->Draw();
      h_dsz->Fit("gaus","","",0.,0.095);
      double l = 0;
      double r = 0.1;
      for(size_t i=0;i<3;i++) {
         h_dsz->Fit("gaus","L","",l,r);
         TF1 *f1 = h_dsz->GetFunction("gaus");
         double mean = f1->GetParameter(1);
         double rms = f1->GetParameter(2);
         l = mean - 2*rms;
         r = mean + 2*rms;
      }
   }
   //
   if (np==21) {
      TH1D* h = new TH1D("h","h",100,0,500);
      for(size_t j=0;j<dsz[layer].size();j++)
         h->Fill(dedxm[layer][j]*sz[layer][j]);
      h->Draw();
      double l = 0;
      double r = 500;
      for(size_t i=0;i<3;i++) {
         h->Fit("gaus","L","",l,r);
         TF1 *f1 = h->GetFunction("gaus");
         double mean = f1->GetParameter(1);
         double rms = f1->GetParameter(2);
         l = mean - 2*rms;
         r = mean + 2*rms;
      }
   }
   //
   if (np==22) {
      TProfile *h = new TProfile("h","h",20,0,0.1,0,1e6);
      for(size_t j=0;j<dsz[layer].size();j++)
         h->Fill(sz[layer][j],dedxm[layer][j]);
      h->Draw();
      TH1D* h1 = new TH1D("h1","h1",20,0,0.1);
      for(size_t j=0;j<h->GetNbinsX();j++)
         h1->SetBinContent(j,h->GetBinError(j)*std::sqrt(h->GetBinEntries(j))*h->GetBinCenter(j));
      h1->Draw();
      //h->Fit("gaus","","",-2000,2000);
   }
   //
   if (np==30) {
      TGraphErrors *gr1 = new TGraphErrors(dt.size(),
                                           ne[layer].data(),
                                           dsz[layer].data(),
                                           dt.data(),
                                           dt.data());
      gr1->Draw("AP");
   }
   //
   if (np==31) {
      dedx_data.z.push_back(std::vector<double>());
      size_t ip = 1;
      std::vector<double> q;
      std::vector<double> k;
      std::vector<double> b;
      std::vector<double> dq;
      std::vector<double> dk;      
      std::vector<double> db;      
      size_t i=1;
      while(i<t.size()) {
         std::vector<double> x;
         std::vector<double> y;
         std::vector<double> z;
         std::vector<double> dx;
         std::vector<double> dy;
         std::vector<double> dz;
         if (i==0) {
            z.push_back(t[0]);
            dz.push_back(0);
            x.push_back(dbp[ip][0]);
            dx.push_back(0);
            y.push_back(dedx0m[layer][0]);
            dy.push_back(ddedx0m[layer][0]);
         }
         while(i<t.size()&&(t[i]-t[i-1])<7200) {
            if (dbp[1][i]&&!TMath::IsNaN(dedx0m[layer][i])) {
               z.push_back(t[i]);
               dz.push_back(0);
               x.push_back(dbp[ip][i]);
               dx.push_back(0);
               y.push_back(dedx0m[layer][i]);
               dy.push_back(ddedx0m[layer][i]);
            }
            i++;
         }
         i++;
         if (x.size()>1) {
            dedx_data.x.resize(0);
            dedx_data.y.resize(0);
            dedx_data.dy.resize(0);
            dedx_data.z.back().resize(0);
            for(size_t j=0;j<z.size();j++) {
               dedx_data.x.push_back(x[j]);
               dedx_data.y.push_back(y[j]);
               dedx_data.dy.push_back(dy[j]);
               dedx_data.z.back().push_back(x[j]);
            }
            dEdx0m_fit();
            //
            TCanvas c("c", "c", 1200, 300);
            c.Divide(3,1);
            c.cd(1);
            TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
            gr->SetMarkerSize(0.5);
            gr->SetMarkerStyle(20);
            gr->Draw("AP");
            //gr->Fit("pol1");
            //
            TF1 *f1 = new TF1("f1",p3p,1,2,2);
            f1->SetParName(0,"amp");
            f1->SetParName(1,"p0");
            f1->SetParameter(0,1000);
            f1->SetParameter(1,0);
            f1->SetParLimits(1,0,1.5);
            gr->Fit("f1");
      for(size_t j=0;j<dedx_data.z.size();j++) {
         std::vector<double> zz;
         std::cout << j << " " << dedx_data.k[j] << " " << dedx_data.b[j] << std::endl;
         for(size_t ii=0;ii<dedx_data.z[j].size();ii++) {
            zz.push_back(dedx_data.z[j][ii]*dedx_data.k[j]+dedx_data.b[j]);
            //std::cout << ii << " " << zz.back() << " " << dedx_data.x[ii] << std::endl;
         }
         TGraphErrors *gr = new TGraphErrors(dedx_data.x.size(),
                                             dedx_data.x.data(),
                                             zz.data(),
                                             dz.data(),
                                             dz.data());
         gr->SetMarkerSize(0.5);
         gr->SetMarkerStyle(20);
         gr->SetMarkerColor(kRed);
         gr->Draw("P");
      }
            //f1->Draw("same");
            //
            //for(size_t i=1;i<x.size();i++) {
            //   TArrow *ar = new TArrow(x[i-1],y[i-1],x[i],y[i]);
            //   ar->Draw();
            //}
            c.cd(2);
            TGraphErrors *gr1 = new TGraphErrors(z.size(),z.data(),y.data(),dz.data(),dy.data());
            gr1->GetXaxis()->SetTimeDisplay(1);
            gr1->SetMarkerSize(0.5);
            gr1->SetMarkerStyle(20);
            gr1->Draw("AP");            
            c.cd(3);
            TGraphErrors *gr2 = new TGraphErrors(z.size(),z.data(),x.data(),dz.data(),dx.data());
            gr2->GetXaxis()->SetTimeDisplay(1);
            gr2->SetMarkerSize(0.5);
            gr2->SetMarkerStyle(20);
            gr2->Draw("AP");            
            //
            /*TF1 *f = gr->GetFunction("pol1");
            q.push_back((s.back()+s.front())/2);
            dq.push_back((s.back()-s.front())/2);
            b.push_back(f->GetParameter(0));
            db.push_back(f->GetParError(0));
            k.push_back(f->GetParameter(1));
            dk.push_back(f->GetParError(1));*/
            //std::cout << b.back() << " " << k.back() << std::endl;
            //TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
            c.Update();
            size_t xxx;
            std::cin >> xxx;
         }
      }
      TGraphErrors *gr = new TGraphErrors(q.size(),q.data(),k.data(),dq.data(),dk.data());
      gr->GetXaxis()->SetTimeDisplay(1);
      gr->SetMarkerSize(0.5);
      gr->SetMarkerStyle(20);
      gr->Draw("AP");
      //TH1D* h = new TH1D("h","h",1000,0,100000);
      //for(size_t i=1;i<t.size();i++)
      //   h->Fill(t[i]-t[i-1]);
      //h->Draw();
   }
}

//-------------------------------------------------------------------

class TFitResultMy
{
 public:
   TFitResultMy(){};
   
   double chi2;
   size_t ndf;
   std::vector<Double_t> par;
   std::vector<Double_t> dpar;
};

double satf(double *x, double *par)
{
   double theta = x[0];
   double sz = par[0];
   double dedx = par[1];
   double demax = par[2];
   double g0, g1, g2;
   return satfsx(dedx, theta, sz, demax, g0, g1, g2 );
}

double satfs(double *x, double *par)
{
   double theta = x[0];
   double sz = par[0];
   double dedx = par[1];
   double demax = par[2]/par[0];
   double g0, g1, g2;
   return satfsx(dedx, theta, sz, demax, g0, g1, g2 );
}

std::vector<double> VI11;
std::vector<double> VI12;
std::vector<double> VI22;
std::vector<double> M1;
std::vector<double> M2;
void readEM()
{
   TFile *fin = TFile::Open("Error_matix.root","READ");
   std::vector<Double_t>  *v;
   fin->GetObject("VI11",v);
   for(std::vector<Double_t>::iterator it = v->begin(); it != v->end(); ++it) VI11.push_back(*it);
   fin->GetObject("VI12",v);
   for(std::vector<Double_t>::iterator it = v->begin(); it != v->end(); ++it) VI12.push_back(*it);
   fin->GetObject("VI22",v);
   for(std::vector<Double_t>::iterator it = v->begin(); it != v->end(); ++it) VI22.push_back(*it);
   fin->GetObject("M1",v);
   for(std::vector<Double_t>::iterator it = v->begin(); it != v->end(); ++it) M1.push_back(*it);
   fin->GetObject("M2",v);
   for(std::vector<Double_t>::iterator it = v->begin(); it != v->end(); ++it) M2.push_back(*it);
}

size_t imx[10] = {7, 9, 9, 9, 8, 7, 7, 6, 6, 5}; // form center like in kk.kumac

size_t flayer = 0;
size_t fInd = 0;
std::vector<double> h_dedx;
std::vector<double> h_ddedx;
std::vector<double> h_theta;
double m_sz, m_dsz;
double m_alpha, m_dalpha;

double satfm(const double *par)
{
   //calculate chisquare
   Double_t chi2 = 0;
   Double_t dchi2;
   double sz = par[0];
   double dedx = par[2];
   double demax = par[1]/par[0];
   double g0, g1, g2;
   size_t i1 = 12-(imx[flayer]-1);
   size_t i2 = 12+(imx[flayer]-1);
   for (size_t i=i1;i<=i2; i++) {
      if (!h_ddedx[i]) continue;
      double dedxs = satfsx(dedx, h_theta[i], sz, demax, g0, g1, g2 );
      dchi2  = (h_dedx[i]-dedxs)/h_ddedx[i];
      chi2 += dchi2*dchi2;
   }
   //dchi2 = (sz - m_sz)/m_dsz;
   //chi2 += dchi2*dchi2;
   //dchi2 = (par[1] - m_alpha)/m_dalpha;
   //chi2 += dchi2*dchi2;
   chi2 += (sz - M1[fInd])*(sz - M1[fInd])*VI11[fInd];
   chi2 += 2*(sz - M1[fInd])*(par[1] - M2[fInd])*VI12[fInd];
   chi2 += (par[1] - M2[fInd])*(par[1] - M2[fInd])*VI22[fInd];
   return chi2;
}

void fcnx(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   f = satfm(par);
}

TFitResultPtr dEdx_fit(TProfile *h, size_t layer)
{
   TF1 *f = new TF1("satf",satf,0,TMath::Pi(),3);
   f->SetParName(0,"sz");
   f->SetParName(1,"dedx0");
   f->SetParName(2,"demax");
   f->SetParameter(0,0.045);
   f->SetParameter(1,700);
   f->SetParameter(2,3500);
   f->FixParameter(0,0.045);
   f->SetParLimits(1,0,20000);
   f->SetParLimits(2,0,20000);
   double xmin = TMath::Pi()/25*(13-imx[layer]);
   double xmax = TMath::Pi() - xmin;
   h->Fit(f,"Q","",xmin,xmax);
   f->ReleaseParameter(0);
   f->SetParLimits(0,0.01,0.2);
   f->FixParameter(1,f->GetParameter(1));
   h->Fit(f,"Q","",xmin,xmax);
   f->ReleaseParameter(1);
   f->SetParLimits(1,0,20000);
   return h->Fit(f,"S","",xmin,xmax);
}

TFitResultPtr dEdx_fits(TProfile *h, size_t layer)
{
   TF1 *f = new TF1("satfs",satfs,0,TMath::Pi(),3);
   f->SetParName(0,"sz");
   f->SetParName(1,"dedx0");
   f->SetParName(2,"alpha");
   f->SetParameter(0,0.045);
   f->SetParameter(1,700);
   f->SetParameter(2,200);
   f->FixParameter(0,0.045);
   f->SetParLimits(1,0,20000);
   f->SetParLimits(2,0,1000);
   double xmin = TMath::Pi()/25*(13-imx[layer]);
   double xmax = TMath::Pi() - xmin;
   h->Fit(f,"Q","",xmin,xmax);
   f->ReleaseParameter(0);
   f->SetParLimits(0,0.01,0.2);
   f->FixParameter(1,f->GetParameter(1));
   h->Fit(f,"Q","",xmin,xmax);
   f->ReleaseParameter(1);
   f->SetParLimits(1,0,10000);
   return h->Fit(f,"S","",xmin,xmax);
}

TFitResultMy* dEdx_fitmo(TProfile *h, size_t run, size_t layer)
{
   flayer = layer;
   h_dedx.resize(0);
   h_ddedx.resize(0);
   h_theta.resize(0);
   for(size_t j=0;j<h->GetNbinsX();j++) {
      h_dedx.push_back(h->GetBinContent(j));
      h_ddedx.push_back(h->GetBinError(j));
      h_theta.push_back(h->GetBinCenter(j));
   }
   if (run<12125) {
      m_sz = 5.09418e-02;
      m_dsz = 1.28295e-02;
      m_alpha = 305.9;
      m_dalpha = 26.1;
   }
   else {
      m_sz = 5.09418e-02;
      m_dsz = 1.28295e-02;
      m_alpha = 188.0;
      m_dalpha = 20.4;
   }
   
   TMinuit *gMinuit = new TMinuit(3);  //initialize TMinuit with a maximum of 5 params
   gMinuit->SetFCN(fcnx);

   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

   // Set starting values and step sizes for parameters
   //gMinuit->mnparm(0, "sz", m_sz, m_dsz, 0,0,ierflg);
   //gMinuit->mnparm(1, "alpha", m_alpha, m_dalpha, 0,0,ierflg);
   //gMinuit->mnparm(2, "dedx0", 2000, 100, 0,0,ierflg);
   gMinuit->DefineParameter(0, "sz", m_sz, m_dsz, m_sz-5*m_dsz, m_sz+5*m_dsz);
   gMinuit->DefineParameter(1, "alpha", m_alpha, m_dalpha, m_alpha-5*m_dalpha, m_alpha+5*m_dalpha);
   gMinuit->DefineParameter(2, "dedx0", 2000, 100, 0, 50000);
   
   // Now ready for minimization step
   arglist[0] = 5000;
   arglist[1] = 1.;
   gMinuit->FixParameter(0);
   gMinuit->FixParameter(1);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->Release(0);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->Release(1);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   
   gMinuit->DeleteArrays();
   return NULL;
}

TFitResultMy*  dEdx_fitm(TProfile *h, size_t run, size_t layer)
{
   flayer = layer;
   h_dedx.resize(0);
   h_ddedx.resize(0);
   h_theta.resize(0);
   for(size_t j=0;j<h->GetNbinsX();j++) {
      h_dedx.push_back(h->GetBinContent(j));
      h_ddedx.push_back(h->GetBinError(j));
      h_theta.push_back(h->GetBinCenter(j));
   }
   if (run<12125) {
      fInd = 2*layer;
      m_sz = M1[fInd];
      m_dsz = 1.28295e-02;
      m_alpha = M2[fInd];
      m_dalpha = 26.1;
   }
   else {
      fInd = 2*layer+1;
      m_sz = M1[fInd];
      m_dsz = 1.28295e-02;
      m_alpha = M2[fInd];
      m_dalpha = 20.4;
   }
   
   const char * minName = "Minuit2";
   const char *algoName = "";

   double chi2;
   size_t ndf;
   const double *p;
   const double *dp;
   {
      ROOT::Math::Minimizer* min =
         ROOT::Math::Factory::CreateMinimizer(minName, algoName);
   
      // set tolerance , etc...
      min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      min->SetMaxIterations(10000);  // for GSL
      min->SetTolerance(0.00001);
      min->SetPrintLevel(1);
       
      // create funciton wrapper for minmizer
      // a IMultiGenFunction type
      ROOT::Math::Functor f(&satfm,3);   
      min->SetFunction(f);
       
      // Set the free variables to be minimized!
      min->SetLimitedVariable(0,"sz",m_sz, m_dsz, m_sz-5*m_dsz,m_sz+5*m_dsz);
      min->SetLimitedVariable(1,"alpha",m_alpha,m_dalpha, m_alpha-5*m_dalpha,m_alpha+5*m_dalpha);
      min->SetLimitedVariable(2,"dedx0",2000, 50,1000,10000);
      
      // do the minimization
      min->SetFixedVariable(0,"sz",m_sz);
      min->SetFixedVariable(1,"alpha",m_alpha);
      min->Minimize();
      p = min->X();
      dp = min->Errors();
   }
   
   {
      ROOT::Math::Minimizer* min =
         ROOT::Math::Factory::CreateMinimizer(minName, algoName);
   
      // set tolerance , etc...
      min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      min->SetMaxIterations(10000);  // for GSL
      min->SetTolerance(0.00001);
      min->SetPrintLevel(1);
       
      // create funciton wrapper for minmizer
      // a IMultiGenFunction type
      ROOT::Math::Functor f(&satfm,3);   
      min->SetFunction(f);
       
      // Set the free variables to be minimized!
      min->SetLimitedVariable(0,"sz",    p[0], dp[0], p[0] - 5*dp[0], p[0] + 5*dp[0] );
      min->SetLimitedVariable(1,"alpha", p[1], dp[1], p[1] - 5*dp[1], p[1] + 5*dp[1] );
      min->SetLimitedVariable(2,"dedx0", p[2], dp[2], p[2] - 5*dp[2], p[2] + 5*dp[2] );
      
      // do the minimization
      min->Minimize();
      p = min->X();
      dp = min->Errors();
      chi2 = min->MinValue();
      ndf = min->NFree();
   }
   
   TFitResultMy *r = new TFitResultMy();;
   r->chi2 = chi2;
   r->ndf = ndf;
   for(size_t i=0;i<ndf;i++) r->par.push_back(p[i]);
   for(size_t i=0;i<ndf;i++) r->dpar.push_back(dp[i]);
   
   return r;
}

void refit(size_t run = 0, size_t layer = 0)
{
   readEM();
   TFile *f = new TFile("xxx.root");
   TIter next_run(f->GetListOfKeys());
   TKey *key_run;
   TString key_run_name("TDirectoryFile");
   TString key_name("TProfile");
   while ((key_run=(TKey*)next_run())) {
      if (!key_run_name.BeginsWith(key_run->GetClassName())) continue;
      TDirectory *d = (TDirectory*)f->Get(key_run->GetName());
      TIter next(d->GetListOfKeys());
      TKey *key;
      while ((key=(TKey*)next())) {
         TString name(key->GetName());
         if (name.BeginsWith(Form("run%08d_layer%d",run,layer))) {
            std::cout << name.Data() << " " << run << " " << layer << std::endl;
            TProfile *h = (TProfile*)d->Get(key->GetName());
            h->Draw();
            dEdx_fit(h,layer);
            dEdx_fits(h,layer);
            dEdx_fitm(h,run,layer);
            //TFitResult *r = (TFitResult*)d->Get(key->GetName());
            //std::cout << r->Parameters().at(0) << std::endl;
         }
      }
   }
}

Double_t fcn2d(Double_t *x, Double_t *p)
{
   Double_t C = 0.7;
   Double_t Dxx = 1, Dyy = 7, Dxy = sqrt(Dxx*Dyy)*C;
   Double_t D = Dxx*Dyy-Dxy*Dxy;
   Double_t V11 = Dyy/D;
   Double_t V12 = -Dxy/D;
   Double_t V22 = Dxx/D;
   Double_t xp = x[0];
   Double_t yp = x[1];
   return exp(-0.5*(xp*xp*V11+2*xp*yp*V12+yp*yp*V22));
}

void td()
{
   TF2 *f = new TF2("gg",fcn2d,-50,50,-50,50);
   f->SetNpx(1000);
   f->SetNpy(1000);
   f->Draw();
   TPrincipal p(2,"ND");
   Double_t data[2];
   std::vector<Double_t> x;
   std::vector<Double_t> y;
   for(size_t i=0;i<1000000;i++) {
      f->GetRandom2(data[0],data[1]);
      p.AddRow(data);
      x.push_back(data[0]);
      y.push_back(data[1]);
   }
   const TMatrixD* m = p.GetCovarianceMatrix();
   const TVectorD* mv  = p.GetMeanValues();
   const TVectorD* sv  = p.GetSigmas();
   std::cout << (*m)(0,0) + (*m)(1,1) << " " << (*m)(1,0)/sqrt((*m)(0,0)*(*m)(1,1)) << std::endl;
   std::cout << (*m)(0,0) << " " << (*m)(1,0) << " " << (*m)(1,1) << std::endl;
   std::cout << (*mv)(0) << " " << (*mv)(1) << std::endl;
   std::cout << (*sv)(0) << " " << (*sv)(1) << std::endl;
   TGraph *gr = new TGraph(x.size(),x.data(),y.data());
   //gr->SetMarkerSize(0.5);
   //gr->SetMarkerStyle(20);
   gr->Draw("AP");
}


double matrconv(double x, double xmin, double xmax, std::vector<double> y)
{
   size_t n = y.size();
   if (!n) return 0;
   if (n==1) return y[0];
   //
   std::vector<double> xi(n,0);
   for(size_t i=0;i<n;i++) xi[i] = (double)i/(n-1);
   //
   double db[n];
   double da[n*n];
   for(size_t i=0;i<n;i++) {
      db[i] = y[i];
      for(size_t j=0;j<n;j++)
         da[n*i+j] = std::pow(xi[i],(double)j);
   }
   //
   TVectorD B(n);
   B.SetElements(db);
   TMatrixD A(n,n);
   A.SetMatrixArray(da);
   //
   TDecompLU lu(A);
   Bool_t ok;
   TVectorD p = lu.Solve(B,ok);
   double f = p(0);
   double xn = (x-xmin)/(xmax-xmin);
   for(size_t i=1;i<n;i++)
      f += p(i)*std::pow(xn,(double)i);
   return f;
}

void matrconv0()
{
   size_t n = 8;
   double xmin = 0;
   double xmax = 3*TMath::Pi();
   TF1 *f0 = new TF1("f0","sin(x)",xmin,xmax);
   f0->Draw();
   std::vector<double> y;
   for(size_t i=0;i<n;i++)
      y.push_back(f0->Eval(xmin+i*(xmax-xmin)/(n-1)));
   std::vector<double> xp;
   std::vector<double> yp;
   for(double x=xmin;x<=xmax;x+=(xmax-xmin)/100) {
      xp.push_back(x);
      yp.push_back(matrconv(x,xmin,xmax,y));
   }
   TGraph *gr = new TGraph(xp.size(),xp.data(),yp.data());
   gr->SetMarkerSize(0.5);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed);
   gr->Draw("P");
}

void tt()
{
   std::vector<double> date, ampl, ddate, dampl;
   for(int i=1; i<26;i++){
      std::cout << i << " " << 27000 - 60 + 120*i << " " << fdate(27000 - 60 + 120*i) << std::endl;
      date.push_back(ftime(27000 - 60 + 120*i));
      ampl.push_back(1);
      dampl.push_back(0);
      ddate.push_back(0);
   }
   TGraphErrors * gr = new TGraphErrors(date.size(),date.data(),ampl.data(),ddate.data(),dampl.data());
   gr->GetXaxis()->SetTimeDisplay(1);
   gr->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}");
   gr->SetMarkerSize(0.5);
   gr->SetMarkerStyle(20);
   gr->Draw("AP");
}

//------------------------------------------

void corr(size_t l1 = 0, size_t l2 = 0)
{
   std::vector<size_t> lyr;
   lyr.push_back(l1);
   lyr.push_back(l2);
   size_t run_min = 1e10;
   size_t run_max = 0;
   std::vector< std::map<size_t,Double_t> > dedx(lyr.size());
   std::vector<TString> files;
   files.push_back(Form("pict/mhad2017-4_0.000000_104.000000_pars_layer%d.root",lyr[0]));
   files.push_back(Form("pict/mhad2017-4_0.000000_105.000000_pars_layer%d.root",lyr[1]));
   TString pref = files.back()(0,files.back().Index(".root"));
   for(size_t i=0;i<lyr.size();i++) {
      TFile *f = new TFile(files[i]);
      TIter next(f->GetListOfKeys());
      TKey *key;
      TString key_name("TFitResult");
      while ((key=(TKey*)next())) {
         TString name(key->GetName());
         if (name.BeginsWith("TFitResult")) {
            int run, layer;
            sscanf(name.Data(),"TFitResult-run%d_layer%d",&run,&layer);
            TFitResult *r = (TFitResult*)f->Get(key->GetName());
            dedx[i][run] = r->Parameters().at(1);
            //std::cout << run << " " << layer << " " << dedx[i][run] << " " << name << std::endl;
            if (run_min>run) run_min = run;
            if (run_max<run) run_max = run;
         }
      }
   }
   std::vector<Double_t> x;
   std::vector<Double_t> y;
   Double_t mean;
   Double_t sig ;
   TF1 *f = NULL;
   for(size_t s=0;s<3;s++) {
      x.resize(0);
      y.resize(0);
      TProfile *pr = new TProfile("pr","pr",50,0,2000);
      TPrincipal p(2,"ND");
      Double_t data[2];
      for(std::map<size_t,Double_t>::iterator it = dedx[0].begin(); it != dedx[0].end(); it++) {
         size_t run = it->first;
         Double_t xi = dedx[0][run];
         Double_t yi = dedx[1][run];
         Bool_t ds = f != NULL ? std::fabs(yi-f->Eval(xi)-mean)<3*sig : true;
         if (ds && dedx[1][run]) {
            x.push_back(xi);
            y.push_back(yi);
            pr->Fill(x.back(),y.back());
            data[0] = xi;
            data[1] = yi;
            p.AddRow(data);
         }
      }
      const TMatrixD* m = p.GetCovarianceMatrix();
      const TVectorD* mv = p.GetMeanValues();
      Double_t l = (*mv)(0)-3*std::sqrt((*m)(0,0));
      Double_t r = (*mv)(0)+3*std::sqrt((*m)(0,0));
      pr->Fit("pol1","","",l,r);
      f = (TF1*)pr->GetFunction("pol1");
      //stop();
      TH1D *h = new TH1D("h","h",200,-100,100);
      for(size_t i=0;i<x.size();i++) h->Fill(y[i]-f->Eval(x[i]));
      l=h->GetMean()-3*h->GetRMS();
      r=h->GetMean()+3*h->GetRMS();
      for(size_t c=0;c<3;c++) {
         h->Fit("gaus","","",l,r);
         TF1 *g = (TF1*)h->GetFunction("gaus");
         g->SetLineColor(kRed);
         mean = g->GetParameter(1);
         sig = g->GetParameter(2);
         l=mean-3*sig;
         r=mean+3*sig;
      }
      //h->Draw();
      //stop();
   }
   x.resize(0);
   y.resize(0);
   std::vector<Double_t> tn;
   std::vector<size_t> rn;
   std::vector<Double_t> rr;
   std::vector<Double_t> x0;
   std::vector<Double_t> y0;
   std::vector<Double_t> xr;
   std::vector<Double_t> yr;
   std::vector<Double_t> rn0;
   std::vector<Double_t> e0;
   std::vector<Double_t> rn1;
   std::vector<Double_t> e1;
   TProfile *pr = new TProfile("pr","pr",100,run_min,run_max);
   for(std::map<size_t,Double_t>::iterator it = dedx[0].begin(); it != dedx[0].end(); it++) {
      size_t run = it->first;
      Double_t xi = dedx[0][run];
      Double_t yi = dedx[1][run];
      if (dedx[1][run]) {
         e0.push_back(yi);
         rn0.push_back(run);
         x0.push_back(xi);
         y0.push_back(yi);
         xr.push_back(run);
         tn.push_back(run);//ftime(run));
         Double_t yf = f->Eval(xi);
         yr.push_back(yi/xi);//yi/yf);
         Bool_t r1 = yi - yf < -3*sig;
         Bool_t r2 = yi/yf < 0.8;
         if (r1&&r2) {
            e1.push_back(yi);
            rn1.push_back(run);
            x.push_back(xi);
            y.push_back(yi);
            rn.push_back(run);
            rr.push_back(yi/yf);
         } 
         r1 = std::fabs(yi - yf) < 3*sig;
         r2 = std::fabs(yi/yf-1) < 0.2;         
         if (r1&&r2) {
            pr->Fill(xr.back(),yr.back());            
         }
      }
   }
   TGraph *gr0 = new TGraph(x0.size(),x0.data(),y0.data());
   gr0->SetMarkerSize(0.5);
   gr0->SetMarkerStyle(20);
   gr0->Draw("AP");
   TGraph *gr = new TGraph(x.size(),x.data(),y.data());
   gr->SetMarkerSize(0.5);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed);
   gr->Draw("P");
   f->Draw("same");
   f->SetLineColor(kRed);
   TLine *ln = new TLine();
   ln->DrawLine(0,0,2000,2000);
   ln->SetLineColor(kBlue);
   stop();
   pr->Draw();
   std::vector<Double_t> xf;
   std::vector<Double_t> dxf;
   std::vector<Double_t> yf;
   std::vector<Double_t> dyf;
   for(size_t i=0;i<pr->GetNbinsX();i++)
      if (pr->GetBinContent(i)&&pr->GetBinError(i)) {
         xf.push_back(pr->GetBinCenter(i));
         dxf.push_back(0);
         yf.push_back(pr->GetBinContent(i));
         dyf.push_back(pr->GetBinError(i));
      }
   TGraphErrors *grf = new TGraphErrors(xf.size(),xf.data(),yf.data(),dxf.data(),dyf.data());
   grf->SetMarkerSize(0.5);
   grf->SetMarkerStyle(20);
   grf->Draw("AP");   
   grf->Fit("pol4");
   TF1 *f1 = (TF1*)grf->GetFunction("pol4");
   f1->SetLineColor(kRed);
   stop();
   TGraph *grr = new TGraph(tn.size(),tn.data(),yr.data());
   //grr->GetXaxis()->SetTimeDisplay(1);
   grr->SetMarkerSize(0.5);
   grr->SetMarkerStyle(20);
   grr->Draw("AP");
   f1->Draw("same");
   stop();
   ofstream out(Form("%s_low_amp.txt",pref.Data()));
   for(size_t i=0;i<rn.size();i++)
      out << rn[i] << " " << rr[i] << std::endl;
   out.close();
   TGraph *ge0 = new TGraph(rn0.size(),rn0.data(),e0.data());
   ge0->SetMarkerSize(0.5);
   ge0->SetMarkerStyle(20);
   ge0->Draw("AP");
   TGraph *ge = new TGraph(rn1.size(),rn1.data(),e1.data());
   ge->SetMarkerSize(0.5);
   ge->SetMarkerStyle(20);
   ge->SetMarkerColor(kRed);
   ge->Draw("P");
   stop();
}

void runtodate(size_t rmin = 25000, size_t rmax = 25100)
{
   ofstream out(Form("run_to_time_%d_%d.txt",rmin,rmax));
   for(size_t run = rmin; run<rmax ; run++)
      out << run << " " << fdatetime(run) << std::endl;
   out.close();
}

//------------------------------------------------------------------------------------
void simcal()
{
   im.resize(0);
   im.push_back(8);
   im.push_back(10);
   im.push_back(10);
   im.push_back(9);
   im.push_back(9);
   im.push_back(9);
   im.push_back(9);
   im.push_back(9);
   im.push_back(7);
   im.push_back(7);
  TFile *f = new TFile("/work/users/konctbel/calibs/R007-001/output/ntuples/neu/ee_bhwide_t18_nemcc-1000-27769-177-250000_wmix.root");
  TTree *t = (TTree*)f->Get("t1");
  TString type = "sim";
  int num = 27769;
  int nx = 25;
  std::map<Int_t,std::vector<Double_t> > cal;
  Int_t run = 1;
  for(size_t i=0;i<10;i++){
    crun.push_back(dEdXclb());
    crun.back().layer = i;
    crun.back().ind = im[crun.back().layer];
    TString hname = Form("%s%08d_layer%d",type.Data(),num,i);
    crun.back().h = new TProfile(hname,hname,nx,0,TMath::Pi(),0,1e6);
    TString pn0 = i ? Form("dEx%d[0]",i) : "dExn[0]";
    TString pn1 = i ? Form("dEx%d[1]",i) : "dExn[1]";
    t->Draw(Form("%s:theta[0]>>%s",pn0.Data(),hname.Data()),"nc>1&&ndex>1");
    crun.back().dedx0 = 1200;
    crun.back().sz = 0.05;
    crun.back().alpha = 200;
    crun.back().slope = 0;
    crun.back().scale = 1;
    //
    crun.back().ind = 7;
    TFitResultPtr r1 = crun.back().fits(7);
    //
    crun.back().ind = im[i];
    TFitResultPtr r2 = crun.back().fits(23);
    //
    crun.back().h->Draw();
    crun.back().h->SetMarkerSize(0.5);
    crun.back().h->SetMarkerStyle(20);
    crun.back().h->GetXaxis()->SetTitle("#theta");
    crun.back().h->GetYaxis()->SetTitle("dE/dx");
    TF1* f1 = crun.back().plot(r1);
    f1->SetLineColor(kBlue);
    TF1* f2 = crun.back().plot(r2);
    f2->SetLineColor(kRed);
   TLegend* legend = new TLegend(0.35,0.2,0.65,0.4);
   legend->AddEntry(crun.back().h,"dEdx(e^{-}) (MC)");
   legend->AddEntry(f1,"scale fixed (=1)");
   legend->AddEntry(f2,"scale free");
   legend->Draw();
    stop();
    //
    if (!i) continue;
    Int_t layer = i;
    if (!cal[run].size()) cal[run] = std::vector<Double_t>(45,0);
    cal[run][layer-1] = r2->Parameters().at(1); // dedxe,dedx0
    cal[run][layer+8] = r2->Parameters().at(2)/r2->Parameters().at(0); // dedxm = alpha/sigz
    cal[run][layer+17] = r2->Parameters().at(4); // scale
    cal[run][layer+26] = r2->Parameters().at(0); // sigz,sz
    cal[run][layer+35] = r2->Parameters().at(3); // slope
  }
  //
  std::vector<TString> suff;
  suff.push_back("dedxe");
  suff.push_back("dedxm");
  suff.push_back("scale");
  suff.push_back("sigz");
  suff.push_back("slope");
  TString pref = "sim";
  for(size_t p=0;p<suff.size();p++)
    {
      ofstream out(Form("%s_%s.txt",pref.Data(),suff[p].Data()));
      std::map<Int_t,std::vector<Double_t> >::iterator it = cal.begin();
      for(;it!=cal.end();it++) {
        out << (*it).first << " ";
        for(size_t i=9*p;i<9*p+9;i++)
          out << std::scientific << (*it).second[i] << " ";
        out << std::endl;
      }
      out.close();
    }
}

//-----------------------------------------------------------------------------------------


double satfsxc( double dedx, double theta, 
                double sz, double demax, double scale, double slope,
                double &g0, double &g1, double &g2 )
{
   double dedxc = dedx + (theta-TMath::Pi()/2)*slope;
   double sth = sin(theta);
   double tth = sth/cos(theta);
   double b = scale/tth;
   double k = sth*sqrt(1 + b*b);
   return k*satfsx(dedxc,theta,sz,demax,g0,g1,g2);
}

double dedxcorr( double dedxr, double thetar, 
                 double sz, double de0, double demax, double scale, double slope )
{
   if (sz==0) return 0;
   double g0 = 0, g1, g2;
   size_t i = 0;
   double dedx = dedxr;
   double theta = thetar;
   double ym = dedxr/2;
   double tm = 0;
   while (fabs(dedxr-g0)>1e-3 && i<20 && dedx<30*de0) {
      satfsxc(dedx,theta,sz,demax,scale,slope,g0,g1,g2);
      if (g1==0) return dedx/de0;
      dedx += (dedxr-g0)/g1;
      i++;
   }
/*   while (ym<dedxr && i<10) {
      satfsxc(dedx,theta,sz,demax,scale,slope,g0,g1,g2);
      tm = -g1/2/g2;
      ym = g2*tm*tm + g1*tm + g0;
      dedx += tm;
      i++;
   }
   dedx -= tm;
   
   i = 0;
   while (fabs(dedxr-g0)>1e-3 && i<10)
   {
      double t1, t2;
      ax2_2bx_c(g2,g1/2,g0-dedxr,t1,t2);
      dedx += t2;
      satfsxc(dedx,theta,sz,demax,scale,slope,g0,g1,g2);
      i++;
   }
   
   if (dedx<0 && dedxr==0) dedx = 0;
   if (dedx<0 && dedxr>0) dedx = 0*de0;
  */ 
   return dedx/de0;
}

void dedxTest()
{
  std::vector<Double_t> dedxe;
  ifstream in1("sim_dedxe.txt");
  while (1) {
    Double_t par[10];
    for(size_t i=0;i<10;i++) in1 >> par[i];
    if (!in1.good()) break;
    for(size_t i=1;i<10;i++) dedxe.push_back(par[i]);
   }
   //
  std::vector<Double_t> dedxm;
  ifstream in2("sim_dedxm.txt");
  while (1) {
    Double_t par[10];
    for(size_t i=0;i<10;i++) in2 >> par[i];
    if (!in2.good()) break;
    for(size_t i=1;i<10;i++) dedxm.push_back(par[i]);
   }
   //
  std::vector<Double_t> sz;
  ifstream in3("sim_sigz.txt");
  while (1) {
    Double_t par[10];
    for(size_t i=0;i<10;i++) in3 >> par[i];
    if (!in3.good()) break;
    for(size_t i=1;i<10;i++) sz.push_back(par[i]);
   }
   //
  std::vector<Double_t> slope;
  ifstream in4("sim_slope.txt");
  while (1) {
    Double_t par[10];
    for(size_t i=0;i<10;i++) in4 >> par[i];
    if (!in4.good()) break;
    for(size_t i=1;i<10;i++) slope.push_back(par[i]);
   }
   //
  std::vector<Double_t> scale;
  ifstream in5("sim_scale.txt");
  while (1) {
    Double_t par[10];
    for(size_t i=0;i<10;i++) in5 >> par[i];
    if (!in5.good()) break;
    for(size_t i=1;i<10;i++) scale.push_back(par[i]);
   }
  for(size_t i=0;i<sz.size();i++)
      std::cout << sz[i] << " " << dedxe[i] << " " << dedxm[i] << " " << scale[i] << " " << slope[i] << std::endl;
   //
   TChain ch("t1");
   ch.Add("/work/users/konctbel/calibs/R007-001/output/ntuples/neu/ee_bhwide_t18_nemcc-1000-27769-177-250000_wmix.root");
   //                                                                                             
   Int_t run, event;
   Float_t beam, eton;
   ch.SetBranchAddress("run",&run);
   ch.SetBranchAddress("event",&event);
   ch.SetBranchAddress("beam",&beam);
   ch.SetBranchAddress("eton",&eton);
   Int_t trin;
   ch.SetBranchAddress("trin",&trin);
   Int_t act;
   ch.SetBranchAddress("act",&act);
   Bool_t cosm;
   ch.SetBranchAddress("cosm",&cosm);
   Int_t np, nc, nn, charge[32];
   ch.SetBranchAddress("np",&np);
   ch.SetBranchAddress("nc",&nc);
   ch.SetBranchAddress("nn",&nn);
   ch.SetBranchAddress("charge",charge);
   Float_t energy[32], theta[32], phi[32], z0[32], d0[32];
   ch.SetBranchAddress("energy",energy);
   ch.SetBranchAddress("theta",theta);
   ch.SetBranchAddress("phi",phi);
   ch.SetBranchAddress("z0",z0);
   ch.SetBranchAddress("d0",d0);
   Int_t ndex;
   ch.SetBranchAddress("ndex",&ndex);
   Float_t dExi[9][11];
   ch.SetBranchAddress("dEx1",dExi[0]);
   ch.SetBranchAddress("dEx2",dExi[1]);
   ch.SetBranchAddress("dEx3",dExi[2]);
   ch.SetBranchAddress("dEx4",dExi[3]);
   ch.SetBranchAddress("dEx5",dExi[4]);
   ch.SetBranchAddress("dEx6",dExi[5]);
   ch.SetBranchAddress("dEx7",dExi[6]);
   ch.SetBranchAddress("dEx8",dExi[7]);
   ch.SetBranchAddress("dEx9",dExi[8]);
   Int_t ndexc;
   ch.SetBranchAddress("ndexc",&ndexc);
   Float_t dExiC[9][11];
   ch.SetBranchAddress("dEx1C",dExiC[0]);
   ch.SetBranchAddress("dEx2C",dExiC[1]);
   ch.SetBranchAddress("dEx3C",dExiC[2]);
   ch.SetBranchAddress("dEx4C",dExiC[3]);
   ch.SetBranchAddress("dEx5C",dExiC[4]);
   ch.SetBranchAddress("dEx6C",dExiC[5]);
   ch.SetBranchAddress("dEx7C",dExiC[6]);
   ch.SetBranchAddress("dEx8C",dExiC[7]);
   ch.SetBranchAddress("dEx9C",dExiC[8]);
   Float_t ei[6];
   ch.SetBranchAddress("e1",&ei[0]);
   ch.SetBranchAddress("e2",&ei[1]);
   ch.SetBranchAddress("e3",&ei[2]);
   ch.SetBranchAddress("e4",&ei[3]);
   ch.SetBranchAddress("e5",&ei[4]);
   ch.SetBranchAddress("e6",&ei[5]);
   Float_t x2ikf1, x2ikf2, x2ikf3;
   ch.SetBranchAddress("x2ikf1",&x2ikf1);
   ch.SetBranchAddress("x2ikf2",&x2ikf2);
   ch.SetBranchAddress("x2ikf3",&x2ikf3);
   Int_t nimskf1;
   ch.SetBranchAddress("nimskf1",&nimskf1);
   Float_t ims1kf1[11];
   ch.SetBranchAddress("ims1kf1",ims1kf1);
   Bool_t col;
   ch.SetBranchAddress("col",&col);
   
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);

   TProfile *h = new TProfile("h","h",100,0,TMath::Pi(),0,1e10);
   TProfile *hc = new TProfile("hc","hc",100,0,TMath::Pi(),0,1e10);
   TH2F* h2 = new TH2F("h2","h2",100,0,TMath::Pi(),100,0,10);
   TH2F* h2c = new TH2F("h2c","h2c",100,0,TMath::Pi(),100,0,10);
   
   Int_t nentries = (Int_t)ch.GetEntries();
   for(int i=0;i<nentries;i++) {
      ch.GetEntry(i);
      // 
      Double_t dph = std::fabs(phi[0]-phi[1])-pi;
      Double_t dth = theta[0]+theta[1]-pi;
      // 
      bool cuts = true;
      cuts &= trin>0 && col;
      cuts &= np>1 && ndex>1;
      cuts &= charge[0]>0 && charge[1]>0;
      cuts &= std::fabs(z0[0])<10 && std::fabs(z0[1])<10;
      cuts &= std::fabs(z0[0]-z0[1])<10;
      cuts &= std::fabs(d0[0])<1 && std::fabs(d0[1])<1;
      if (!cuts) continue;
      h->Fill(theta[0],dExi[0][0]);
      h->Fill(theta[1],dExi[0][1]);
      Double_t dedxc = dedxcorr(dExi[0][0],theta[0],sz[0],dedxe[0],dedxm[0],scale[0],slope[0]);
      if (dedxc<10) hc->Fill(theta[0],dedxc);
      dedxc = dedxcorr(dExi[0][1],theta[1],sz[0],dedxe[0],dedxm[0],scale[0],slope[0]);
      if (dedxc<10) hc->Fill(theta[1],dedxc);
      //
      //h2->Fill(theta[0],dExi[0][0]);
      //h2->Fill(theta[1],dExi[0][1]);
      //h2c->Fill(theta[0],dedxcorr(dExi[0][0],theta[0],sz[0],dedxe[0],dedxm[0],scale[0],slope[0]));
      //h2c->Fill(theta[1],dedxcorr(dExi[0][1],theta[0],sz[0],dedxe[0],dedxm[0],scale[0],slope[0]));
   }
   hc->Draw();
   //h->Draw("same");
}


double drw[10] = {0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4};
double rw[10] = {2.4, 2.4, 3.1, 3.9, 4.7, 5.5, 6.3, 7.1, 7.9, 8.7};
//double zw(9) = {15.7, 15.7, 15.7, 15.7, 15.45, 14.95, 14.45, 13.95, 13.7};
double zw[10] = {14.0, 14.0, 14.0, 14.0, 14.0, 13.8, 13.25, 12.7, 12.5, 11.6};

void dedxNew()
{
   TChain ch("t1");
   ch.Add("/work/users/konctbel/calibs/R007-001/output/ntuples/neu/ee_bhwide_t18_nemcc-1000-27769-177-250000_wmix.root");
   //                                                                                             
   Int_t run, event;
   Float_t beam, eton;
   ch.SetBranchAddress("run",&run);
   ch.SetBranchAddress("event",&event);
   ch.SetBranchAddress("beam",&beam);
   ch.SetBranchAddress("eton",&eton);
   Int_t trin;
   ch.SetBranchAddress("trin",&trin);
   Int_t act;
   ch.SetBranchAddress("act",&act);
   Bool_t cosm;
   ch.SetBranchAddress("cosm",&cosm);
   Int_t np, nc, nn, charge[32];
   ch.SetBranchAddress("np",&np);
   ch.SetBranchAddress("nc",&nc);
   ch.SetBranchAddress("nn",&nn);
   ch.SetBranchAddress("charge",charge);
   Float_t energy[32], theta[32], phi[32], z0[32], d0[32];
   ch.SetBranchAddress("energy",energy);
   ch.SetBranchAddress("theta",theta);
   ch.SetBranchAddress("phi",phi);
   ch.SetBranchAddress("z0",z0);
   ch.SetBranchAddress("d0",d0);
   Int_t ndex;
   ch.SetBranchAddress("ndex",&ndex);
   Float_t dExi[9][11];
   ch.SetBranchAddress("dEx1",dExi[0]);
   ch.SetBranchAddress("dEx2",dExi[1]);
   ch.SetBranchAddress("dEx3",dExi[2]);
   ch.SetBranchAddress("dEx4",dExi[3]);
   ch.SetBranchAddress("dEx5",dExi[4]);
   ch.SetBranchAddress("dEx6",dExi[5]);
   ch.SetBranchAddress("dEx7",dExi[6]);
   ch.SetBranchAddress("dEx8",dExi[7]);
   ch.SetBranchAddress("dEx9",dExi[8]);
   Int_t ndexc;
   ch.SetBranchAddress("ndexc",&ndexc);
   Float_t dExiC[9][11];
   ch.SetBranchAddress("dEx1C",dExiC[0]);
   ch.SetBranchAddress("dEx2C",dExiC[1]);
   ch.SetBranchAddress("dEx3C",dExiC[2]);
   ch.SetBranchAddress("dEx4C",dExiC[3]);
   ch.SetBranchAddress("dEx5C",dExiC[4]);
   ch.SetBranchAddress("dEx6C",dExiC[5]);
   ch.SetBranchAddress("dEx7C",dExiC[6]);
   ch.SetBranchAddress("dEx8C",dExiC[7]);
   ch.SetBranchAddress("dEx9C",dExiC[8]);
   Float_t ei[6];
   ch.SetBranchAddress("e1",&ei[0]);
   ch.SetBranchAddress("e2",&ei[1]);
   ch.SetBranchAddress("e3",&ei[2]);
   ch.SetBranchAddress("e4",&ei[3]);
   ch.SetBranchAddress("e5",&ei[4]);
   ch.SetBranchAddress("e6",&ei[5]);
   Float_t x2ikf1, x2ikf2, x2ikf3;
   ch.SetBranchAddress("x2ikf1",&x2ikf1);
   ch.SetBranchAddress("x2ikf2",&x2ikf2);
   ch.SetBranchAddress("x2ikf3",&x2ikf3);
   Int_t nimskf1;
   ch.SetBranchAddress("nimskf1",&nimskf1);
   Float_t ims1kf1[11];
   ch.SetBranchAddress("ims1kf1",ims1kf1);
   Bool_t col;
   ch.SetBranchAddress("col",&col);
   
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);

   std::vector<TH1F*> h;
   size_t nh = 25;
   for(size_t i=0;i<nh;i++)
     h.push_back(new TH1F(Form("h%d",i),Form("h%d",i),1000,0,10000));
   Double_t dTh = TMath::Pi()/nh;
   
   Int_t nentries = (Int_t)ch.GetEntries();
   for(int i=0;i<nentries;i++) {
      ch.GetEntry(i);
      // 
      bool selected = true;
      selected &= trin > 0;
      selected &= nc >= 2;
      if (!selected) continue;
      selected &= std::fabs(z0[0]) < 7;
      selected &= std::fabs(z0[1]) < 7;
      selected &= std::fabs(z0[0]-z0[1]) < 1;
      selected &= std::fabs(d0[0]) < 1;
      selected &= std::fabs(d0[1]) < 1;
      selected &= std::fabs(std::fabs(phi[0]-phi[1])-TMath::Pi()) < 1.5*TMath::Pi()/180;
      selected &= std::fabs(theta[0]+theta[1]-TMath::Pi()) < 5*TMath::Pi()/180;
      selected &= energy[0]/beam > 0.8;
      selected &= energy[1]/beam > 0.7;
      if (!selected) continue;
      for(size_t p=0;p<2;p++) {
        size_t ii = i+1;
        double tth = std::sin(theta[p])/std::cos(theta[p]);
        double zr = z0[p] + (rw[ii]+drw[ii])/tth;
        bool s = std::fabs(zr)<zw[ii]-2;
        if (s&&dExi[0][p]) h[theta[p]/dTh]->Fill(dExi[0][p]);
      }
   }
   //
   Int_t nq = 10;
   Double_t q[10];
   for(size_t i=0;i<nq;i++) q[i] = (i+1.)/nq;
   Double_t t[nh][10];
   for(size_t i=3;i<nh-3;i++)
     h[i]->GetQuantiles(nq,t[i],q);
   for(int i=nq-2;i>=0;i--) {
     std::vector<Double_t> x;
     std::vector<Double_t> y;
     for(size_t j=3;j<nh-3;j++) {
       x.push_back(j);
       y.push_back(t[j][i]);
     }
     TGraph *gr = new TGraph(x.size(),x.data(),y.data());
     if (i==nq-2) gr->Draw("APL");
     else gr->Draw("PL");
     gr->SetMarkerSize(0.5);
     gr->SetMarkerStyle(20);
     stop();
   }
   //
   /*   h[nh/2]->Scale(1.0/h[nh/2]->Integral());
   h[nh/2]->GetCumulative()->Draw();
   for(size_t i=3;i<nh-3;i++) {
     h[i]->Scale(1.0/h[i]->Integral());
     h[i]->GetCumulative()->Draw("same");
     }*/
}

//---------------------------------------------------------------

std::pair<Double_t,Double_t> mr(std::vector<Double_t> x)
{
   std::sort(x.begin(),x.end());
   Double_t m = 0;
   Double_t n = 0;
   Double_t m2 = 0;
   for(size_t i=3;i<x.size()-3;i++) {
      n++;
      m += x[i];
      m2 += x[i]*x[i];
   }
   m /= n;
   m2 /= n;
   return std::make_pair(m,sqrt((m2-m*m)/n));
}

void thetaCorrTestplot(Int_t ip = 2)
{
   TString r = gSystem->GetFromPipe("ls exp*_corr.txt");
   TObjArray *tr = r.Tokenize("\n");
   std::vector<TGraph*> gr;
   std::map<Double_t, std::vector<Double_t> > pk;
   for (Int_t i = 0; i < tr->GetEntries(); i++) {
      TString fname = ((TObjString *)(tr->At(i)))->String();
      ifstream in(fname);
      std::vector<Double_t> vk;
      std::vector<Double_t> vf;
      while (1) {
         Double_t n, x0, y0, k, b, c;
         in >> n >> x0 >> y0 >> k >> b >> c;
         if (n==15) continue;
         if (n>18) break;
         if (!in.good()) break;
         vf.push_back((n+0.5)*0.05);
         Double_t p;
         if (ip==0) p = x0;
         if (ip==1) p = y0;
         if (ip==2) p = k;
         if (ip==3) p = b;
         if (ip==4) p = c;
         if (ip==5) p = y0-k*x0;
         vk.push_back(p);
         pk[n].push_back(p);
      }
      gr.push_back(new TGraph(vf.size(),vf.data(),vk.data()));
      gr.back()->SetMarkerSize(0.5);
      gr.back()->SetMarkerStyle(20);
   }
   std::vector<Double_t> s;
   for(size_t i=0;i<gr.size();i++) {
      Double_t xi = 0;
      for(size_t j=0;j<gr.size();j++)
         for(size_t k=0;k<gr[i]->GetN();k++) {
            Double_t x1, y1;
            gr[i]->GetPoint(k,x1,y1);
            Double_t x2, y2;
            gr[j]->GetPoint(k,x2,y2);
            Double_t dy = y1-y2;
            xi += dy*dy;
         }
      s.push_back(xi);
      std::cout << i << " " << xi << std::endl;
   }
   for(size_t i=0;i<gr.size();i++) {
      //if (s[i]>410) continue;
      if (i) gr[i]->Draw("P");
      else   gr[i]->Draw("AP");
   }
   std::vector<Double_t> x;
   std::vector<Double_t> dx;
   std::vector<Double_t> y;
   std::vector<Double_t> dy;
   std::map<Double_t, std::vector<Double_t> >::iterator it = pk.begin();
   for(;it!=pk.end();it++) {
      std::pair<Double_t,Double_t> p = mr((*it).second);
      x.push_back(((*it).first+0.5)*0.05);
      dx.push_back(0);
      y.push_back(p.first);
      dy.push_back(p.second);
   }
   TGraphErrors *grs = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
   grs->SetLineColor(kRed);
   grs->SetLineWidth(3);
   grs->Draw("same");
   grs->Fit("pol5");
   TF1* f = grs->GetFunction("pol5");
   f->SetLineColor(kBlue);
}

