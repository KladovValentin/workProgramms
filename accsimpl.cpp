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

bool myOrderInt(std::pair<Int_t,Int_t> a, std::pair<Int_t,Int_t> b) 
{ 
   return a.first < b.first; 
} 

size_t ftime(size_t run)
{
   TString cmd;
   cmd = Form("mysql -h sndfarm09.sndonline -u daqread -pdaqread rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"",run);
   cmd = cmd + " | " + Form("fgrep -e %d",run);
   cmd = cmd + " | " + Form("sed 's/%d//g'",run);
   cmd = Form("date -d \"$(%s)\" +\"%s\"",cmd.Data(),"%s");
   return std::atoi(gSystem->GetFromPipe(cmd.Data())) - (size_t)gStyle->GetTimeOffset() - 6*3600;
}

size_t fdate(size_t run)
{
   TString cmd;
   cmd = Form("mysql -h sndfarm09.sndonline -u daqread -pdaqread rundbase -e \"select c_run,c_starttime from t_readout where c_run=%d\"",run);
   cmd = cmd + " | " + Form("fgrep -e %d",run);
   cmd = cmd + " | " + Form("sed 's/%d//g'",run);
   cmd = Form("date -d \"$(%s)\" +\"%s%s%s\"",cmd.Data(),"%Y","%m","%d");
   return std::atoi(gSystem->GetFromPipe(cmd.Data()));
}

size_t fdatetotime(size_t date)
{
   TString cmd;
   cmd = Form("date -d \"%d\" +\"%s\"",date,"%s");
   return std::atoi(gSystem->GetFromPipe(cmd.Data())) - (size_t)gStyle->GetTimeOffset() - 6*3600;
}

std::vector< std::pair<Int_t,Int_t> > getCalRun(TString calName = "dcpagen", bool tmrn = true, bool read = false)
{
   TString cmd;
   // reading cals list from db
   cmd = "cd /work/users/konctbel/calibs/R007-001/; setup2k.sh; clbixlist " + calName + " CURRENT > " + calName + ".list";
   if (read) gSystem->Exec(cmd.Data());
   // start run
   cmd = "cat /work/users/konctbel/calibs/R007-001/" + calName + ".list | sed '/^-/ d' | awk '{print $1}'";
   TString r1 = gSystem->GetFromPipe(cmd.Data());
   TObjArray *tr1 = r1.Tokenize("\n");
   // stop run
   cmd = "cat /work/users/konctbel/calibs/R007-001/" + calName + ".list | sed '/^-/ d' | awk '{print $2}'";
   TString r2 = gSystem->GetFromPipe(cmd.Data());
   TObjArray *tr2 = r2.Tokenize("\n");
   // save to db time
   //cmd = "cat /work/users/konctbel/calibs/R007-001/" + calName + ".list | sed '/^-/ d' | awk '{print $8\" \"$9}'";
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

class AccSimpl
{
 public:
   AccSimpl(){};
   size_t run;
   size_t t;
   std::vector<Double_t> amp1pe;
};

class AccProxy
{
 public:
   AccProxy(){};
   size_t run;
   size_t t;
   std::vector<Double_t> a;
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
      cmd = "cd "+dir+"; setup2k.sh; clbarray " + cal + " " + array + " " + usage + " " + Form("%d",run) +" > " + fname;
      gSystem->Exec(cmd.Data());
      std::cout << "Calinration array " << run << ":" << cal << ":" << array << "was stored on disk" << std::endl;
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
   std::vector< std::pair<Int_t,Int_t> > res = getCalRun(name,tm,false);
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

std::vector<AccSimpl> getAccSimpl(Int_t tStart, Int_t tStop, bool tm = true, TString usage = "RECONSTR")
{
   TString name = "accsimpl";
   std::vector< std::pair<Int_t,Int_t> > res = getCalRun(name,tm,false);
   std::vector<AccSimpl> cal;
   for(size_t i=0;i<res.size();i++) {
      if (res[i].first<tStart||res[i].first>tStop) continue;
      cal.push_back(AccSimpl());
      cal.back().run = res[i].second;
      cal.back().t = ftime(res[i].second);
      cal.back().amp1pe = getCalArray(res[i].second,name,"amp1pe",usage);
   }
   return cal;
}

std::vector<AccProxy> getAccProxy(Int_t tStart, Int_t tStop, bool tm = true, TString usage = "RECONSTR")
{
   TString name = "accproxy";
   std::vector< std::pair<Int_t,Int_t> > res = getCalRun(name,tm,false);
   std::vector<AccProxy> cal;
   for(size_t i=0;i<res.size();i++) {
      if (res[i].first<tStart||res[i].first>tStop) continue;
      cal.push_back(AccProxy());
      cal.back().run = res[i].second;
      cal.back().t = ftime(res[i].second);
      cal.back().a = getCalArray(res[i].second,name,"a",usage);
   }
   return cal;
}

void accsimple(size_t ch = 0)
{
   std::vector<AccSimpl> cal = getAccSimpl(0,44710,false,"CURRENT");
   std::vector<double> x;
   std::vector<double> y;
   for(size_t i=0;i<cal.size();i++) {
      x.push_back(cal[i].run);
      y.push_back(cal[i].amp1pe[ch]);
   }
   TGraph *gr = new TGraph(x.size(),x.data(),y.data());
   //gr->GetXaxis()->SetTimeDisplay(1);
   gr->SetMarkerSize(1);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed);
   gr->Draw("AP");
}

void accproxy(size_t ch = 0)
{
   std::vector<AccProxy> cal = getAccProxy(19999,40272,false,"CURRENT");
   std::vector<double> x;
   std::vector<double> y;
   std::cout << cal.size() << std::endl;
   for(size_t i=0;i<cal.size();i++) {
      x.push_back(cal[i].t);
      y.push_back(cal[i].a[ch]);
   }
   TGraph *gr = new TGraph(x.size(),x.data(),y.data());
   gr->GetXaxis()->SetTimeDisplay(1);
   gr->SetMarkerSize(1);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed);
   gr->Draw("AP");
}

//-------------------------------------------------------------

class Event
{
 public:
   Event(){};
   void setAccProxy(std::vector<AccProxy> &cal);
   int run, trin;
   float beam;
   int nc;
   float energy[31], theta[31], phi[31], z0[31], d0[31];
   float thetas[31], phis[31];
   int counter[31], numparts[31], region[31];
   float amplitude[31], time[31], length[31];
   int ndex;
   float dexi[10][31];
   int nch;
   float cch[31], ach[31], tch[31];
   int nchr;
   float cchr[31], schr[31], achr[31], tchr[31], cachr[31], ctchr[31];
   size_t rtime;
   std::vector<Double_t> accproxy;
};

void Event::setAccProxy(std::vector<AccProxy> &cal)
{
   size_t i=0;
   while(i<cal.size()&&(size_t)run>cal[i].run)i++;
   if (i) i--;
   accproxy = cal[i].a;
}

//----------------------------------------------------------

class ChHit
{
 public:
   ChHit(size_t num_ = 0, TString type_ = "run"){num = num_; type = type_;};
   void init(size_t nx = 25);
   void write(TString fname = "xxx.root");
   void open(TString fname = "xxx.root");
   void close();
   void fill(Event &evt);
   //void fit();
   void plot();
   void timeCalc();

   TFile *file;
   TString type;
   std::map<size_t,TH1F*> time_spect;
   std::map<size_t,TH1F*> amp_spect_pds;
   std::map<size_t,TH1F*> amp_spect_gx;
   std::map<size_t,TH1F*> amp_spect_g0;
   std::map<size_t,TProfile*> amp_phi_gx;
   std::map<size_t,TProfile*> amp_phi_g0;
   std::map<size_t,TProfile*> amp_z0_gx;
   std::map<size_t,TProfile*> amp_z0_g0;
   TProfile* amp_gx;
   TProfile* amp_g0;
   TProfile* eff_gx;
   TProfile* eff_g0;
   std::vector<size_t> tm;
   std::vector<size_t> t;
   size_t num;
   size_t nevt;
   size_t nx;
   std::vector<TCanvas*> pict;
};

void ChHit::init(size_t nx_)
{
   nx = nx_;
   for(int i=1;i<10;i++) {
      time_spect[i] = new TH1F(Form("%s%08d_time_spect_counter%d",type.Data(),num,i),"",200,0,200);
      amp_spect_pds[i] = new TH1F(Form("%s%08d_amp_spect_pds_counter%d",type.Data(),num,i),"",200,-2,2);
      amp_spect_gx[i] = new TH1F(Form("%s%08d_amp_spect_gx_counter%d",type.Data(),num,i),"",310,-1,30);
      amp_spect_g0[i] = new TH1F(Form("%s%08d_amp_spect_g0_counter%d",type.Data(),num,i),"",310,-1,30);
      amp_phi_gx[i] = new TProfile(Form("%s%08d_amp_phi_gx_counter%d",type.Data(),num,i),"",120,(2*i-4)*TMath::Pi()/9,(2*i+2)*TMath::Pi()/9);
      amp_phi_g0[i] = new TProfile(Form("%s%08d_amp_phi_g0_counter%d",type.Data(),num,i),"",120,(2*i-4)*TMath::Pi()/9,(2*i+2)*TMath::Pi()/9);
      amp_z0_gx[i] = new TProfile(Form("%s%08d_amp_z0_gx_counter%d",type.Data(),num,i),"",200,-25,25);
      amp_z0_g0[i] = new TProfile(Form("%s%08d_amp_z0_g0_counter%d",type.Data(),num,i),"",200,-25,25);
   }
   amp_gx = new TProfile(Form("%s%08d_amp_gx",type.Data(),num),"",9,1,10);
   amp_g0 = new TProfile(Form("%s%08d_amp_g0",type.Data(),num),"",9,1,10);
   eff_gx = new TProfile(Form("%s%08d_eff_gx",type.Data(),num),"",9,1,10);
   eff_g0 = new TProfile(Form("%s%08d_eff_g0",type.Data(),num),"",9,1,10);
   nevt = 0;
   file = 0;
   t.resize(4,0);
   pict.resize(10,NULL);
}

void ChHit::timeCalc()
{
   if (!tm.size()) return;
   double tmm = 0;
   double tmin = 1e20;
   double tmax = 0;
   for(size_t i=0;i<tm.size();i++) {
      tmm += (double)tm[i]-(double)tm[0];
      if (tm[i] < tmin ) tmin = tm[i];
      if (tm[i] > tmax ) tmax = tm[i];
   }
   tmm /= tm.size();
   tmm += tm[0];
   double dtmm = 0;
   for(size_t i=0;i<tm.size();i++) {
      double dt = tmm-(double)tm[i];
      dtmm += dt*dt;
   }
   dtmm = sqrt(dtmm);
   t[0] = (size_t)tmm;
   t[1] = (size_t)dtmm;
   t[2] = (size_t)tmin;
   t[3] = (size_t)tmax;
}

void ChHit::open(TString fname)
{
   file = TFile::Open(fname,"RECREATE");
}

void ChHit::write(TString fname)
{
   timeCalc();
   file = (TFile*)gROOT->FindObject(fname);
   if (!file) open(fname);
   TDirectory *dir = file->mkdir(Form("%s%08d",type.Data(),num));
   dir->cd();
   //for(size_t i=0;i<pict.size();i++) pict[i]->Write();
   for(size_t i=1;i<10;i++) {
      time_spect[i]->Write();
      amp_spect_pds[i]->Write();
      amp_spect_gx[i]->Write();
      amp_spect_g0[i]->Write();
      amp_phi_gx[i]->Write();
      amp_phi_g0[i]->Write();
      amp_z0_gx[i]->Write();
      amp_z0_g0[i]->Write();
   }
   amp_gx->Write();
   amp_g0->Write();
   eff_gx->Write();
   eff_g0->Write();
   dir->WriteObjectAny(&t,"std::vector<size_t>",Form("time_%s%08d",type.Data(),num));
   dir->cd("../");
}

void ChHit::close()
{
   file->Close();
}

void test()
{
   TFile *fout = TFile::Open("TEST.root","RECREATE");
   
   std::vector<Double_t>  xx; // I create and fill the vector to be saved.
   for( int k = 0; k <3; k++ ) xx.push_back(k);
   fout->WriteObjectAny(&xx,"std::vector<Double_t>","xx"); // I store the vector in the TFile
   fout->ls();
   fout->Close();
   
   TFile *fin = TFile::Open("TEST.root","READ");
   std::vector<Double_t>  *yy;
   fin->GetObject("xx",yy); // I try to retrieve the vector
   for(std::vector<Double_t>::iterator it = yy->begin(); it != yy->end(); ++it) {
      std::cout << *it << '\n';
   }
}


void ChHit::fill(Event &evt)
{
   std::map<size_t,std::vector<double> > ach;
   std::map<size_t,std::vector<double> > tch;
   for(size_t i=0;i<(size_t)evt.nch;i++) {
      size_t nch = (size_t)evt.cch[i];
      ach[nch].push_back(evt.ach[i]);
      tch[nch].push_back(evt.tch[i]);
   }
   //
   std::map<size_t,std::vector<double> > achr;
   std::map<size_t,std::vector<double> > tchr;
   for(size_t i=0;i<(size_t)evt.nchr;i++) {
      size_t nchr = (size_t)evt.cchr[i];
      achr[nchr].push_back(evt.achr[i]);
      tchr[nchr].push_back(evt.tchr[i]);
   }
   //
   bool selected = true;
   selected &= evt.trin > 0;
   selected &= evt.nc >= 2;
   if (!selected) return;
   // beam spot
   selected &= std::fabs(evt.z0[0]) < 7;
   selected &= std::fabs(evt.z0[1]) < 7;
   selected &= std::fabs(evt.z0[0]-evt.z0[1]) < 1;
   selected &= std::fabs(evt.d0[0]) < 1;
   selected &= std::fabs(evt.d0[1]) < 1;
   // collinearity
   selected &= std::fabs(std::fabs(evt.phi[0]-evt.phi[1])-TMath::Pi()) < 1.5*TMath::Pi()/180;
   selected &= std::fabs(evt.theta[0]+evt.theta[1]-TMath::Pi()) < 5*TMath::Pi()/180;
   // energy
   selected &= evt.energy[0]/evt.beam > 0.7;
   if (!selected) return;
   tm.push_back(evt.rtime);
   //
   for(size_t i=1;i<10;i++) {
      for(size_t p=0;p<2;p++) {
         double dphic = evt.phi[p] - (2*i-1)*TMath::Pi()/9.0 - 2.5*TMath::Pi()/180;
         double dphis = evt.phi[p] - (2*i-1)*TMath::Pi()/9.0 + 2.5*TMath::Pi()/180;
         bool phiSelected = 
            std::fabs(dphic)<15*TMath::Pi()/180 && 
            std::fabs(dphis)>3*TMath::Pi()/180;
         double zc = evt.z0[p] + 14/std::tan(evt.theta[p]);
         bool thetaSelected = std::fabs(zc)<10;
         int dt = (int)(tchr[i].back() - evt.accproxy[3*i-2]);
         bool thereAreHits = ach[i].size()>1;
         bool mainGate = (dt>-10&&dt<8)||(!thereAreHits);
         if (phiSelected&&thetaSelected) {
            if (thereAreHits) time_spect[i]->Fill(tchr[i].back());
            amp_spect_pds[i]->Fill(ach[i].front());
            amp_spect_gx[i]->Fill(ach[i].back());
            amp_gx->Fill(i,ach[i].back());
            eff_gx->Fill(i,thereAreHits);
            if (mainGate) {
               amp_spect_g0[i]->Fill(ach[i].back());
               amp_g0->Fill(i,ach[i].back());
               eff_g0->Fill(i,thereAreHits);
            }
         }
         if (phiSelected) {
            amp_z0_gx[i]->Fill(zc,ach[i].back());
            if (mainGate) amp_z0_g0[i]->Fill(zc,ach[i].back());
         }
         if (thetaSelected) {
            double phi = evt.phi[p];
            if (i==1&&phi>TMath::Pi()) phi -= 2*TMath::Pi();
            if (i==9&&phi<TMath::Pi()) phi += 2*TMath::Pi();
            amp_phi_gx[i]->Fill(phi,ach[i].back());
            if (mainGate) amp_phi_g0[i]->Fill(phi,ach[i].back());
         }
      }
   }
   //
/*   for(size_t i=0;i<dedx_th.size();i++) {
      for(size_t p=0;p<2;p++) {
         double tth = std::sin(evt.theta[p])/std::cos(evt.theta[p]);
         double zr = evt.z0[p] + (rw[i]+drw[i])/tth;
         bool s = std::fabs(zr)<zw[i]-1;
         if (s&&evt.dexi[i][!p]>0) dedx_th.at(i)->Fill(evt.theta[p],evt.dexi[i][p]);
      }
   }
   if (evt.thetas[0]&&evt.theta[0]!=evt.thetas[0])
      zCorr->Fill(std::tan(evt.theta[0])/std::tan(evt.thetas[0]));
   if (evt.thetas[1]&&evt.theta[1]!=evt.thetas[1])
      zCorr->Fill(std::tan(evt.theta[1])/std::tan(evt.thetas[1]));*/
}

void ChHit::plot()
{
/*   for(size_t i=0;i<dedx_th.size();i++) {
      TString name(Form("layer %d",i));
      pict[i] = (TCanvas*)(gROOT->GetListOfCanvases()->At(0)->DrawClone());
      pict[i]->SetName(name);
      pict[i]->SetTitle(name);
      TLegend *legend = new TLegend(0.05,0.05,0.4,0.3);
      legend->SetHeader(name);
      dedx_th[i]->Draw();
      dedx_th[i]->SetMarkerSize(0.5);
      dedx_th[i]->SetMarkerStyle(20);
      dedx_th[i]->GetXaxis()->SetTitle("#theta");
      dedx_th[i]->GetYaxis()->SetTitle("dE/dx");
      legend->AddEntry(dedx_th[i],"ChHit profile");
      if (f_fit[i]) {
         f_fit[i]->SetLineColor(kBlack);
         f_fit[i]->SetLineWidth(7);
         f_fit[i]->Draw("same");
         legend->AddEntry(Form("f_fit_%d",i),"ChHit saturation function No 1","l");
      }
      if (f_fits[i]) {
         f_fits[i]->SetLineColor(kBlue);
         f_fits[i]->SetLineWidth(3);
         f_fits[i]->Draw("same");
         legend->AddEntry(Form("f_fits_%d",i),"ChHit saturation function No 2","l");
      }
      if (f_fitm[i]) {
         f_fitm[i]->SetLineColor(kRed);
         f_fitm[i]->SetLineWidth(3);
         f_fitm[i]->Draw("same");
         legend->AddEntry(Form("f_fitm_%d",i),"ChHit saturation function No 3","l");
      }
      if (f_fitmo[i]) {
         f_fitmo[i]->SetLineColor(kGreen);
         f_fitmo[i]->SetLineWidth(3);
         f_fitmo[i]->Draw("same");
         legend->AddEntry(Form("f_fitmo_%d",i),"ChHit saturation function No 4","l");
      }
      legend->Draw();
      dedx_th[i]->Draw("same");
      pict[i]->Update();
   }*/
}

void ampspect(TString scan = "mhad2017") {
   gROOT->Reset();
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   //
   std::vector<AccProxy> accproxy = getAccProxy(19999,40272,false,"CURRENT");
   //
   Event evt;
   TChain ch("t1");//"h1");
   //ch.Add(Form("/work/users/konctbel/KKc/exp/%s/*.root",scan.Data()));
   ch.Add("/work/users/konctbel/calibs/R007-001/output/ntuples/col/*.root");
   ch.SetBranchAddress("run",&evt.run);
   ch.SetBranchAddress("trin",&evt.trin);
   ch.SetBranchAddress("beam",&evt.beam);
   ch.SetBranchAddress("nc",&evt.nc);
   ch.SetBranchAddress("energy",evt.energy);
   ch.SetBranchAddress("theta",evt.theta);
   ch.SetBranchAddress("phi",evt.phi);
   ch.SetBranchAddress("z0",evt.z0);
   ch.SetBranchAddress("d0",evt.d0);
   ch.SetBranchAddress("thetas",evt.thetas);
   ch.SetBranchAddress("phis",evt.phis);
   ch.SetBranchAddress("counter",evt.counter);
   ch.SetBranchAddress("numparts",evt.numparts);
   ch.SetBranchAddress("amplitude",evt.amplitude);
   ch.SetBranchAddress("time",evt.time);
   ch.SetBranchAddress("length",evt.length);
   ch.SetBranchAddress("region",evt.region);
   ch.SetBranchAddress("nch",&evt.nch);
   ch.SetBranchAddress("cch",&evt.cch);
   ch.SetBranchAddress("ach",&evt.ach);
   ch.SetBranchAddress("tch",&evt.tch);
   ch.SetBranchAddress("nchr",&evt.nchr);
   ch.SetBranchAddress("cchr",&evt.cchr);
   ch.SetBranchAddress("schr",&evt.schr);
   ch.SetBranchAddress("achr",&evt.achr);
   ch.SetBranchAddress("tchr",&evt.tchr);
   ch.SetBranchAddress("cachr",&evt.cachr);
   ch.SetBranchAddress("ctchr",&evt.ctchr);
   ch.SetBranchAddress("ndex",&evt.ndex);
   ch.SetBranchAddress("dExn",&evt.dexi[0]);
   ch.SetBranchAddress("dEx1",&evt.dexi[1]);
   ch.SetBranchAddress("dEx2",&evt.dexi[2]);
   ch.SetBranchAddress("dEx3",&evt.dexi[3]);
   ch.SetBranchAddress("dEx4",&evt.dexi[4]);
   ch.SetBranchAddress("dEx5",&evt.dexi[5]);
   ch.SetBranchAddress("dEx6",&evt.dexi[6]);
   ch.SetBranchAddress("dEx7",&evt.dexi[7]);
   ch.SetBranchAddress("dEx8",&evt.dexi[8]);
   ch.SetBranchAddress("dEx9",&evt.dexi[9]);
   Int_t nevent = ch.GetEntries();
   std::map<int,int> runcal;
   std::map<int,int> runtime;
   std::map<int,int> runday;
   std::map<int,int> nevtr;
   std::map<int,ChHit> chit;
   std::map<int,int> nevtrD;
   std::map<int,ChHit> chitD; 
   std::map<int,int> nevtrE;
   std::map<int,ChHit> chitE; 
   std::map<int,int> nevtrC;
   std::map<int,ChHit> chitC; 
   //nevent = 1000000;
   int dnevent = (int)((double)nevent/100+0.5);
   for (Int_t jentry=0;jentry<nevent;jentry++) {   
      Long64_t ientry = ch.LoadTree(jentry);
      if (ientry < 0)
         break;
      ch.GetEntry(jentry);
      evt.setAccProxy(accproxy);
      if (!(jentry%dnevent)) std::cout << jentry/dnevent << "% of events processed " << std::endl;
      if (!runtime[evt.run]) runtime[evt.run] = ftime(evt.run);
      evt.rtime = runtime[evt.run];
      // split events by the run
      if(!nevtr[evt.run]) {
         nevtr[evt.run]++;
         chit[evt.run] = ChHit(evt.run,"run");
         chit[evt.run].init();
      }
      chit[evt.run].fill(evt);
      chit[evt.run].nevt++;
      // split events y the day
      if (!runday[evt.run]) runday[evt.run] = fdate(evt.run);
      size_t date = runday[evt.run];
      if(!nevtrD[date]) {
         nevtrD[date]++;
         chitD[date] = ChHit(date,"day");
         chitD[date].init();
      }
      chitD[date].fill(evt);
      chitD[date].nevt++;
      // split events by the beam energy
      size_t beam10 = (size_t)(10*evt.beam);
      if(!nevtrE[beam10]) {
         nevtrE[beam10]++;
         chitE[beam10] = ChHit(beam10,"beam");
         chitE[beam10].init();
      }
      chitE[beam10].fill(evt);
      chitE[beam10].nevt++;
   }
   TString fname = Form("acc_%s.root",scan.Data());
   std::map<int,ChHit>::iterator it;
   double cnt = 0;
   for(it = chit.begin(); it != chit.end(); it++) {
      //it->second.correction();
      //it->second.fit();
      //it->second.plot();
      it->second.write(fname);
      cnt++;
      std::cout << 100*cnt/chit.size() << "%   " << it->first << " run was fitted and saved" << std::endl;
   }
   cnt = 0;
   for(it = chitD.begin(); it != chitD.end(); it++) {
      //it->second.correction();
      //it->second.fit();
      //it->second.plot();
      it->second.write(fname);
      cnt++;
      std::cout << 100*cnt/chitD.size() << "%   " << it->first << " day was fitted and saved" << std::endl;
   }
   cnt = 0;
   for(it = chitE.begin(); it != chitE.end(); it++) {
      //it->second.correction();
      //it->second.fit();
      //it->second.plot();
      it->second.write(fname);
      cnt++;
      std::cout << 100*cnt/chitE.size() << "%   " << it->first << " beam was fitted and saved" << std::endl;
   }
   TFile *file = (TFile*)gROOT->FindObject(fname);
   if (file) file->Close();
   for(it = chit.begin(); it != chit.end(); it++)
      std::cout << it->second.type << " " << it->second.num << " " << it->second.nevt << std::endl;
   for(it = chitD.begin(); it != chitD.end(); it++)
      std::cout << it->second.type << " " << it->second.num << " " << it->second.nevt << std::endl;
   for(it = chitE.begin(); it != chitE.end(); it++)
      std::cout << it->second.type << " " << it->second.num << " " << it->second.nevt << std::endl;
}

class AccMap
{
 public:
   AccMap(){init();};
   void init();
   size_t ind;
   size_t run;
   size_t t;
   std::vector<Double_t> time;
   std::vector<Double_t> dtime;
   std::vector<Double_t> amp_g0;
   std::vector<Double_t> amp_gx;
   std::vector<Double_t> damp_g0;
   std::vector<Double_t> damp_gx;
   std::vector<Double_t> eff_g0;
   std::vector<Double_t> eff_gx;
   std::vector<Double_t> deff_g0;
   std::vector<Double_t> deff_gx;
};

void AccMap::init()
{
}


class AccProxyCal
{
 public:
   AccProxyCal(size_t r1=19999, size_t r2=40000, bool rt=false, TString usage="RECONSTR");
   void plot(size_t ch = 0);
   std::vector<AccProxy> cal;
   std::vector<TGraph*> gr;
};

AccProxyCal::AccProxyCal(size_t r1, size_t r2, bool rt, TString usage)
{
   cal = getAccProxy(r1,r2,rt,usage);
   for(size_t ch=0;ch<9;ch++) {
      std::vector<Double_t> x;
      std::vector<Double_t> y;
      for(size_t i=0;i<cal.size();i++) {
         x.push_back(cal[i].t);
         y.push_back(cal[i].a[3*ch+1]+8);
      }
      for(int i=cal.size()-1;i>=0;i--) {
         x.push_back(cal[i].t);
         y.push_back(cal[i].a[3*ch+1]-10);
      }
      gr.push_back(new TGraph(x.size(),x.data(),y.data()));
      gr.back()->GetXaxis()->SetTimeDisplay(1);
      gr.back()->SetMarkerSize(0.5);
      gr.back()->SetMarkerStyle(20);
      gr.back()->SetMarkerColor(kBlue);
      gr.back()->SetLineColor(kBlue);
      gr.back()->SetFillStyle(3000);
   }
}

void AccProxyCal::plot(size_t ch)
{
   gr[ch]->Draw("PL");   
}

void plotpar()
{
   AccProxyCal accproxy(19999,40272,false,"CURRENT");
   //
   double tmin = 1e20;
   double tmax = 0;
   std::vector<AccMap> acc;
   //
   TFile *f = new TFile("acc_mhad2017.root");
   TIter next_dir(f->GetListOfKeys());
   TKey *key_dir;
   while ((key_dir=(TKey*)next_dir())) {
      if (!TString("TDirectoryFile").BeginsWith(key_dir->GetClassName())) continue;
      TDirectory *d = (TDirectory*)f->Get(key_dir->GetName());
      TIter next(d->GetListOfKeys());
      TKey *key;
      while ((key=(TKey*)next())) {
         TString name(key->GetName());
         TString type = "day";
         if (name.BeginsWith(type)&&
             name.Contains("_amp_g0")) {
            acc.push_back(AccMap());
            sscanf(name.Data(),Form("%s%s_amp_g0",type.Data(),"%d"),&acc.back().ind);
            acc.back().t = fdatetotime(acc.back().ind);
            //
            TProfile* h_amp_g0 = (TProfile*)d->Get(Form("%s%08d_amp_g0",type.Data(),acc.back().ind));
            for(int i=1;i<=h_amp_g0->GetNbinsX();i++) {
               acc.back().amp_g0.push_back(h_amp_g0->GetBinContent(i));
               acc.back().damp_g0.push_back(h_amp_g0->GetBinError(i));
            }
            //
            TProfile* h_amp_gx = (TProfile*)d->Get(Form("%s%08d_amp_gx",type.Data(),acc.back().ind));
            for(int i=1;i<=h_amp_gx->GetNbinsX();i++) {
               acc.back().amp_gx.push_back(h_amp_gx->GetBinContent(i));
               acc.back().damp_gx.push_back(h_amp_gx->GetBinError(i));
            }
            //
            TProfile* h_eff_g0 = (TProfile*)d->Get(Form("%s%08d_eff_g0",type.Data(),acc.back().ind));
            for(int i=1;i<=h_eff_g0->GetNbinsX();i++) {
               acc.back().eff_g0.push_back(h_eff_g0->GetBinContent(i));
               acc.back().deff_g0.push_back(h_eff_g0->GetBinError(i));
            }
            //
            TProfile* h_eff_gx = (TProfile*)d->Get(Form("%s%08d_eff_gx",type.Data(),acc.back().ind));
            for(int i=1;i<=h_eff_gx->GetNbinsX();i++) {
               acc.back().eff_gx.push_back(h_eff_gx->GetBinContent(i));
               acc.back().deff_gx.push_back(h_eff_gx->GetBinError(i));
            }
            //
            for(int i=1;i<10;i++) {
               TH1F *h = (TH1F*)d->Get(Form("%s%08d_time_spect_counter%d",type.Data(),acc.back().ind,i));
               acc.back().time.push_back(h->GetBinCenter(h->GetMaximumBin()));
               acc.back().dtime.push_back(0);
            }
            //std::vector<size_t> *t = (std::vector<size_t> *)d->Get(Form("time_day%08d",acc.back().ind));
            //std::cout << acc.back().run << " " << (*t)[0] << " " << (*t)[1] << std::endl;
            //acc.back().t = t->size()<1 ? 0 : (*t)[0];
            //acc.back().dt = t->size()<2 ? 0 : (*t)[1];
            //acc.back().tmin = t->size()<3 ? 0 : (*t)[2];
            //acc.back().tmax = t->size()<4 ? 0 : (*t)[3];
         }
      }
   }
   for(size_t p=3;p<5;p++) {
      for(size_t c=0;c<9;c++) {
         std::vector<double> x;
         std::vector<double> dx;
         std::vector<double> y0;
         std::vector<double> dy0;
         std::vector<double> y1;
         std::vector<double> dy1;
         for(size_t i=0;i<acc.size();i++) {
            if (p==0) {
               y0.push_back(acc[i].amp_g0[c]);
               dy0.push_back(acc[i].damp_g0[c]);
               y1.push_back(acc[i].amp_gx[c]);
               dy1.push_back(acc[i].damp_gx[c]);
            }
            if (p==1) {
               y0.push_back(acc[i].eff_g0[c]);
               dy0.push_back(acc[i].deff_g0[c]);
               y1.push_back(acc[i].eff_gx[c]);
               dy1.push_back(acc[i].deff_gx[c]);
            }
            if (p==2) {
               y0.push_back(acc[i].eff_g0[c]);
               dy0.push_back(acc[i].deff_g0[c]);
               y1.push_back(1-std::exp(-acc[i].amp_g0[c]));
               dy1.push_back(std::exp(-acc[i].amp_g0[c])*acc[i].damp_g0[c]);
            }
            if (p==3) {
               y0.push_back(acc[i].eff_gx[c]);
               dy0.push_back(acc[i].deff_gx[c]);
               y1.push_back(1-std::exp(-acc[i].amp_gx[c]));
               dy1.push_back(std::exp(-acc[i].amp_gx[c])*acc[i].damp_gx[c]);
            }
            if (p==4) {
               y0.push_back(acc[i].time[c]);
               dy0.push_back(acc[i].dtime[c]);
            }
            x.push_back(acc[i].t);
            dx.push_back(0);
         }
         TGraphErrors *gr0 = new TGraphErrors(x.size(),x.data(),y0.data(),dx.data(),dy0.data());
         gr0->GetXaxis()->SetTimeDisplay(1);
         gr0->SetMarkerSize(0.5);
         gr0->SetMarkerStyle(20);
         gr0->SetMarkerColor(kRed);
         gr0->SetLineColor(kRed);
         gr0->Draw("AP");
         TGraphErrors *gr1 = new TGraphErrors(x.size(),x.data(),y1.data(),dx.data(),dy1.data());
         gr1->GetXaxis()->SetTimeDisplay(1);
         gr1->SetMarkerSize(0.5);
         gr1->SetMarkerStyle(20);
         gr1->SetMarkerColor(kBlue);
         gr1->SetLineColor(kBlue);
         gr1->Draw("P");
         accproxy.plot(c);
         if (!p) gr0->GetYaxis()->SetRangeUser(0.0,10);
         else gr0->GetYaxis()->SetRangeUser(0.8,1.1);
         if (p==4) gr0->GetYaxis()->SetRangeUser(0.0,200);
         TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
         c->Update();
         int xxx;
         cin>>xxx;
      }
   }
}

//------------------------------------------------------------------------

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
   double d = (x1-x2)*(x1-x3)*(x2-x3)*(x1-x4)*(x2-x4)*(x3-x4);
   std::vector<double> p(4,0);
   p[0] = (x1*(x1-x3)*x3*(x1-x4)*(x3-x4)*x4*y2+
           x2*pow(x4,2)*(-pow(x3,3)*y1+pow(x3,2)*x4*y1+pow(x1,2)*(x1-x4)*y3)+
           pow(x1,2)*x2*pow(x3,2)*(-x1+x3)*y4+
           pow(x2,3)*(x4*(-pow(x3,2)*y1+x3*x4*y1+x1*(x1-x4)*y3)+
                  x1*x3*(-x1+x3)*y4)+
           pow(x2,2)*(x1*x4*(-pow(x1,2)+pow(x4,2))*y3+
                  pow(x3,3)*(x4*y1-x1*y4)+x3*(-pow(x4,3)*y1+pow(x1,3)*y4)))/d;
   p[1] = (pow(x1,2)*(x1-x4)*pow(x4,2)*(y2-y3)+
           pow(x3,3)*(pow(x4,2)*(y1-y2)+pow(x1,2)*(y2-y4))+
           pow(x2,2)*(pow(x4,3)*(y1-y3)+pow(x1,3)*(y3-y4)+pow(x3,3)*(-y1+y4))+
           pow(x3,2)*(pow(x4,3)*(-y1+y2)+pow(x1,3)*(-y2+y4))+
           pow(x2,3)*(pow(x4,2)*(-y1+y3)+pow(x3,2)*(y1-y4)+pow(x1,2)*(-y3+y4)))/d;
   p[2] = (-x1*(x1-x4)*x4*(x1+x4)*(y2-y3)+
           x3*(pow(x4,3)*(y1-y2)+pow(x1,3)*(y2-y4))+
           pow(x3,3)*(-x4*y1-x1*y2+x4*y2+x1*y4)+
           pow(x2,3)*(x4*y1+x1*y3-x4*y3-x1*y4+x3*(-y1+y4))+
           x2*(pow(x4,3)*(-y1+y3)+pow(x3,3)*(y1-y4)+pow(x1,3)*(-y3+y4)))/d;
   p[3] = (x1*(x1-x4)*x4*(y2-y3)+
           pow(x3,2)*(x4*y1+x1*y2-x4*y2-x1*y4)+
           pow(x2,2)*(-x4*y1-x1*y3+x4*y3+x3*(y1-y4)+x1*y4)+
           x2*(pow(x4,2)*(y1-y3)+pow(x1,2)*(y3-y4)+pow(x3,2)*(-y1+y4))+
           x3*(pow(x4,2)*(-y1+y2)+pow(x1,2)*(-y2+y4)))/d;
   return p;
}

double p3g0(double x, double xr, double sr, std::vector<double> p)
{
   double p0,p1,p2,p3;
   double arg;
   p0=p[0];
   p1=p[1];
   p2=p[2];
   p3=p[3];
   arg=(x-xr)/sqrt(2.0)/sr;
   return exp(-arg*arg)*sr/sqrt(2.0*TMath::Pi())*
      (p1+p2*(x+xr)+p3*(x*x+x*xr+xr*xr+2*sr*sr))+
      (1+erf((arg)))/2*
      (p0+x*(p1+x*(p2+p3*x))+(p2+3*p3*x)*sr*sr);
}

double p3gl(double x, double xr, double sr, double er, double xer, double der, std::vector<double> p)
{
   double p0,p1,p2,p3;
   double earg, arg;
   p0=p[0];
   p1=p[1];
   p2=p[2];
   p3=p[3];
   arg=(x - xr)/(sqrt(2.0)*sr);
   earg=exp(-arg*arg);
   return (2*(der*er*(sqrt(2*TMath::Pi())*(p0 + p2*(sr*sr + x*x) +
                                  x*(p1 + p3*(3*sr*sr + x*x))) +
                      sr*(p1 + p2*(x + xr) +
                           p3*(2*sr*sr + x*x + x*xr + xr*xr))*earg) -
              (-1 + er)*(sqrt(2*TMath::Pi())*
                         (p1*sr*sr + 3*p3*pow(sr,4) + p0*x + 3*p2*sr*sr*x + p1*x*x +
                          6*p3*sr*sr*x*x + p2*pow(x,3) + p3*pow(x,4) -
                          (p0 + p2*(sr*sr + x*x) +
                           x*(p1 + p3*(3*sr*sr + x*x)))*xer) +
                         sr*(p0 + p1*(x - xer + xr) +
                              p2*(2*sr*sr + x*x - x*xer + x*xr - xer*xr + xr*xr) +
                              p3*(x*x*(x - xer) + x*(x - xer)*xr + (x - xer)*xr*xr +
                                  pow(xr,3) + sr*sr*(5*x - 2*xer + 3*xr)))*earg)) +
           sqrt(2*TMath::Pi())*(-(der*er*(p0 + p2*(sr*sr + x*x) +
                                 x*(p1 + p3*(3*sr*sr + x*x)))) +
                       (-1 + er)*(p1*sr*sr + 3*p3*pow(sr,4) + p0*x +
                                  3*p2*sr*sr*x + p1*x*x +
                                  6*p3*sr*sr*x*x + p2*pow(x,3) + p3*pow(x,4) -
                                  (p0 + p2*(sr*sr + x*x) +
                                   x*(p1 + p3*(3*sr*sr + x*x)))*xer))*
           (1 - erf(arg)))/(2.*der*sqrt(2*TMath::Pi()));
}
   
double p3g(double x, double xr, double sr, std::vector<double> p, double smin, double smax)
{
/*   double scal = TMath::Pi()/180;
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
   return p3g0(x,xr,sr,p);
}

double aerogel(size_t n, double x, double *par)
{
   std::vector<double> xv(4);
   std::vector<double> yv(4);
   for(size_t i=0;i<4;i++)
      yv[i] = par[4*n+i];
   double s1, s2;
   switch (n) {
    case 0:
      xv[0] = par[12];
      xv[3] = par[13]-1.5/120.0;
      s1 = par[17];
      s2 = par[18];
      break;
    case 1:
      xv[0] = par[13]+1.5/120.0;
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
   xv[1] = xv[0]+(xv[3]-xv[0])/3;
   xv[2] = xv[3]-(xv[3]-xv[0])/3;
   double smin = xv[0];
   double smax = xv[3];
   std::vector<double> pv = m4(xv,yv);
   return p3g(x,xv[0],s1,pv,smin,smax)-
      p3g(x,xv[3],s2,pv,smin,smax);
}

double shifter(double x, double *par)
{
   double ds = 1.5/120.0;
   return par[16]*(erf((x-(par[13]-ds))/sqrt(2.0)/par[18])-
                   erf((x-(par[13]+ds))/sqrt(2.0)/par[18]))/2;
}

double scphi(double *x, double *par)
{
   double f = 0;
   for(size_t i=0;i<3;i++)
      f += aerogel(i,x[0],par);
   f += shifter(x[0],par);
   return f;
}

void scphi_plot()
{
   double c[4] = {1,0.2,-0.03,0.005};
   double x[4] = {1,2,3,4};
   double y[4];
   TF1 *f = new TF1("p4",Form("(%f)+(%f)*x+(%f)*x*x+(%f)*x*x*x",c[0],c[1],c[2],c[3]),x[0]-1,x[3]+1);
   f->Draw();
   for(size_t i=0;i<4;i++)
      y[i] = f->Eval(x[i]);
   TGraph *gr = new TGraph(4,x,y);
   gr->SetMarkerSize(2);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed);
   gr->SetLineColor(kRed);
   gr->Draw("P");
   //
   size_t nc = 0;
   double scal = TMath::Pi()/180;
   std::vector<double> par(30);
   par[0] = 3;
   par[1] = 3.2;
   par[2] = 3.5;
   par[3] = 4;
   par[4] = 4;
   par[5] = 3.5;
   par[6] = 3.2;
   par[7] = 3;
   par[8] = 3;
   par[9] = 2.9;
   par[10] = 2.7;
   par[11] = 2.6;
   //
   par[12] = nc*40*scal;
   par[13] = par[12]+15*scal;
   par[14] = par[13]+12.5*scal;
   par[15] = par[12]+40*scal;
   //
   par[16] = 10;
   par[17] = 0.5*scal;
   par[18] = 0.5*scal;
   par[19] = 0.5*scal;
   //par[] = ;
   //
   double xmin = par[12]-0.2;
   double xmax = par[15]+0.2;
   std::vector<double> ys;
   std::vector<double> xv;
   std::vector<double> yv;
   for(double xi=xmin;xi<xmax;xi+=0.001) {
      xv.push_back(xi);
      yv.push_back(shifter(xi,par.data()));
      ys.push_back(yv.back());
   }
   TGraph *gr1 = new TGraph(xv.size(),xv.data(),yv.data());
   gr1->SetMarkerSize(0.1);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerColor(kGreen);
   gr1->SetLineColor(kGreen);
   gr1->Draw("APL");   
   for(size_t i=0;i<3;i++) {
      std::vector<double> yv;
      for(double xi=xmin;xi<xmax;xi+=0.001) {
         yv.push_back(aerogel(i,xi,par.data()));
         ys[yv.size()-1] += yv.back();
      }
      TGraph *gr = new TGraph(xv.size(),xv.data(),yv.data());
      gr->SetMarkerSize(0.1);
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(kBlue);
      gr->SetLineColor(kBlue);
      gr->Draw("PL");   
   }
   TGraph *grs = new TGraph(xv.size(),xv.data(),ys.data());
   grs->SetMarkerSize(0.1);
   grs->SetMarkerStyle(20);
   grs->SetMarkerColor(kBlack);
   grs->SetLineColor(kBlack);
   grs->Draw("PL");
   //
   TF1* fn = new TF1("scphi",scphi,xmin,xmax,20);
   fn->SetParameters(par.data());
   fn->SetLineColor(kRed);
   fn->Draw("same");
}


TFile* f = new TFile("/work/users/kladov/snd2k/R006-003/maindir/fittedprofilesdate.root");
TDirectory * d1 = (TDirectory*)f->Get("date1,sensor6");
//TF1* tf1 = (TF1*)d1->Get("date1,sensor6,full");
TProfile* tf1 = (TProfile*)d1->Get("date1,sensor6,profile");
tf1->SetLineColor(1);
tf1->Draw();
for (int i = 2; i < 6; i++) {
    TDirectory* d = (TDirectory*)f->Get(Form("date%d,sensor6", i));
    //TF1* tf = (TF1*)d->Get(Form("date%d,sensor6,full", i));
    TProfile* tf = (TProfile*)d->Get(Form("date%d,sensor6,profile",i));
    tf->SetLineColor(i);
    tf->Draw("same");
}


void scphi_fit()
{
   double scal = TMath::Pi()/180;
   TFile *f = new TFile("acc_mhad2017.root");
   TIter next_dir(f->GetListOfKeys());
   TKey *key_dir;
   while ((key_dir=(TKey*)next_dir())) {
      if (!TString("TDirectoryFile").BeginsWith(key_dir->GetClassName())) continue;
      TDirectory *d = (TDirectory*)f->Get(key_dir->GetName());
      TIter next(d->GetListOfKeys());
      TKey *key;
      while ((key=(TKey*)next())) {
         TString name(key->GetName());
         TString type = "day";
         if (name.BeginsWith(type)&&
             name.Contains("_amp_phi_gx_counter")) {
            size_t ind, counter;
            sscanf(name.Data(),Form("%s%s_amp_phi_gx_counter%s",type.Data(),"%d","%d"),&ind,&counter);
            //
            if (ind<20170401) continue; 
            TProfile* h_amp_phi = (TProfile*)d->Get(Form("%s%08d_amp_phi_gx_counter%d",type.Data(),ind,counter));
            h_amp_phi->Draw();
            std::vector<double> par(20);
            par[0] = 5;
            par[1] = 5;
            par[2] = 5;
            par[3] = 5;
            par[4] = 5;
            par[5] = 5;
            par[6] = 5;
            par[7] = 5;
            par[8] = 5;
            par[9] = 5;
            par[10] = 5;
            par[11] = 5;
            //
            double dphi = 2.4*scal;
            par[12] = (counter-1)*40*scal+dphi;
            par[13] = par[12]+15*scal;
            par[14] = par[13]+15*scal;
            par[15] = par[12]+40*scal;
            //
            par[16] = 30;
            par[17] = 0.3*scal;
            par[18] = 0.3*scal;
            par[19] = 0.3*scal;
            //
            TF1* fn = new TF1("scphi",scphi,par[12]-15*scal,par[15]+15*scal,20);
            for(size_t i=0;i<12;i++)
               fn->SetParName(i,Form("y%d",i));
            fn->SetParName(12,"x1");
            fn->SetParName(13,"x2");
            fn->SetParName(14,"x3");
            fn->SetParName(15,"x4");
            fn->SetParName(16,"ys");
            fn->SetParName(17,"#sigma_{1}");
            fn->SetParName(18,"#sigma_{s}");
            fn->SetParName(19,"#sigma_{2}");
            fn->SetParameters(par.data());
            fn->SetLineColor(kRed);
            fn->SetNpx(1000);
            //
            for(size_t i=0;i<12;i++)
               fn->FixParameter(i,fn->GetParameter(i));
            //fn->FixParameter(13,fn->GetParameter(13));
            fn->FixParameter(14,fn->GetParameter(14));
            fn->FixParameter(16,fn->GetParameter(16));
            fn->FixParameter(17,fn->GetParameter(17));
            fn->FixParameter(18,fn->GetParameter(18));
            fn->FixParameter(19,fn->GetParameter(19));
            h_amp_phi->Fit("scphi");
            for(size_t i=0;i<12;i++)
               fn->ReleaseParameter(i);
            h_amp_phi->Fit("scphi");
            fn->ReleaseParameter(16);
            h_amp_phi->Fit("scphi");
            fn->ReleaseParameter(17);
            fn->ReleaseParameter(18);
            fn->ReleaseParameter(19);
            h_amp_phi->Fit("scphi");
            //fn->Draw("same");
            for(size_t c=0;c<3;c++) {
               TF1* fna = new TF1(Form("aerogel%d",c),scphi,par[12]-15*scal,par[15]+15*scal,20);
               fna->SetParameters(fn->GetParameters());
               for(size_t i=0;i<12;i++)
                  fna->SetParameter(i,0);
               for(size_t i=c*4;i<c*4+4;i++)
                  fna->SetParameter(i,fn->GetParameter(i));
               fna->SetParameter(16,0);
               fna->SetLineColor(kBlue);
               fna->SetLineWidth(0);
               fna->SetNpx(1000);
               fna->Draw("same");
            }
            //
            TF1* fns = new TF1("shifter",scphi,par[12]-15*scal,par[15]+15*scal,20);
            fns->SetParameters(fn->GetParameters());
            for(size_t i=0;i<12;i++)
               fns->SetParameter(i,0);
            fns->SetLineColor(kGreen);
            fns->SetLineWidth(0);
            fns->SetNpx(1000);
            fns->Draw("same");
            h_amp_phi->GetXaxis()->SetRangeUser(par[12]-15*scal,par[15]+15*scal);
            h_amp_phi->GetYaxis()->SetRangeUser(0,1.5*h_amp_phi->GetMaximum());
            //
            TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
            c->Update();
            size_t xxx;
            cin>>xxx;
         }
      }
   }
}


void mucorr()
{
   double er = 0.994;
   double mr = 3.58;
   double D = 0.2;
   double d = 0;
   double mu = mr;
   for(size_t i=0;i<10;i++) {
      std::cout << i << " " << mu << " " << d << std::endl;
      mu = mr/(1-D*d/2*exp(-mu));
      d = (1-(1-er)*exp(mu))/mu;
   }
}
