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
#include "TProfile2D.h"
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
#include "TText.h"
#include "TLine.h"

void stop()
{
   TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
   c->Update();
   std::cout << "Push any key!" << std::endl;
   char key;
   std::cin>>key;
}

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

TString ftimetodate(size_t t)
{
   TString cmd;
   cmd = Form("date -d @%d",t+(size_t)gStyle->GetTimeOffset()+6*3600);
   return gSystem->GetFromPipe(cmd.Data());
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

void accsimple_pl(size_t ch = 0)
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

void accproxy_pl(size_t ch = 0)
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
   void fill(Event &evt, Int_t stage);
   //void fit();
   void plot();
   void timeCalc();
   void setCorr(Double_t corr_ = 1){corr = corr_;};

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
   std::map<size_t,TProfile2D*> amp_2d;
   std::map<size_t,TProfile2D*> len_2d;
   std::map<size_t,TProfile2D*> amp_2d_gx;
   std::map<size_t,TProfile2D*> amp_2d_g0;
   std::map<size_t,TProfile2D*> eff_2d_gx;
   std::map<size_t,TProfile2D*> eff_2d_g0;
   std::map<size_t,TProfile2D*> effc_2d_gx;
   std::map<size_t,TProfile2D*> effc_2d_g0;
   TProfile* amp_gx;
   TProfile* amp_g0;
   TProfile* eff_gx;
   TProfile* eff_g0;
   TProfile* effc_gx;
   TProfile* effc_g0;
   std::vector<size_t> tm;
   std::vector<size_t> t;
   size_t num;
   size_t nevt;
   size_t nx;
   std::vector<TCanvas*> pict;
   Double_t corr;
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
      double phi1 = (2*i-1)*TMath::Pi()/9.0 + 2.5*TMath::Pi()/180 - 20*TMath::Pi()/180;
      double phi2 = (2*i-1)*TMath::Pi()/9.0 + 2.5*TMath::Pi()/180 + 19*TMath::Pi()/180;
      amp_2d[i] = new TProfile2D(Form("%s%08d_amp_2d_counter%d",type.Data(),num,i),"",50,-25,25,13,phi1,phi2);
      len_2d[i] = new TProfile2D(Form("%s%08d_len_2d_counter%d",type.Data(),num,i),"",50,-25,25,13,phi1,phi2);
      amp_2d_gx[i] = new TProfile2D(Form("%s%08d_amp_2d_gx_counter%d",type.Data(),num,i),"",50,-25,25,13,phi1,phi2);
      amp_2d_g0[i] = new TProfile2D(Form("%s%08d_amp_2d_g0_counter%d",type.Data(),num,i),"",50,-25,25,13,phi1,phi2);
      eff_2d_gx[i] = new TProfile2D(Form("%s%08d_eff_2d_gx_counter%d",type.Data(),num,i),"",50,-25,25,13,phi1,phi2);
      eff_2d_g0[i] = new TProfile2D(Form("%s%08d_eff_2d_g0_counter%d",type.Data(),num,i),"",50,-25,25,13,phi1,phi2);
      effc_2d_gx[i] = new TProfile2D(Form("%s%08d_effc_2d_gx_counter%d",type.Data(),num,i),"",50,-25,25,13,phi1,phi2);
      effc_2d_g0[i] = new TProfile2D(Form("%s%08d_effc_2d_g0_counter%d",type.Data(),num,i),"",50,-25,25,13,phi1,phi2);
   }
   amp_gx = new TProfile(Form("%s%08d_amp_gx",type.Data(),num),"",9,1,10);
   amp_g0 = new TProfile(Form("%s%08d_amp_g0",type.Data(),num),"",9,1,10);
   eff_gx = new TProfile(Form("%s%08d_eff_gx",type.Data(),num),"",9,1,10);
   eff_g0 = new TProfile(Form("%s%08d_eff_g0",type.Data(),num),"",9,1,10);
   effc_gx = new TProfile(Form("%s%08d_effc_gx",type.Data(),num),"",9,1,10);
   effc_g0 = new TProfile(Form("%s%08d_effc_g0",type.Data(),num),"",9,1,10);
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
      amp_2d[i]->Write();
      len_2d[i]->Write();
      amp_2d_gx[i]->Write();
      amp_2d_g0[i]->Write();
      eff_2d_gx[i]->Write();
      eff_2d_g0[i]->Write();
      effc_2d_gx[i]->Write();
      effc_2d_g0[i]->Write();
   }
   amp_gx->Write();
   amp_g0->Write();
   eff_gx->Write();
   eff_g0->Write();
   effc_gx->Write();
   effc_g0->Write();
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


void ChHit::fill(Event &evt, Int_t stage)
{
   std::map<size_t,std::vector<double> > ach;
   std::map<size_t,std::vector<double> > tch;
   for(size_t i=0;i<(size_t)evt.nch;i++) {
      size_t nch = (size_t)evt.cch[i];
      ach[nch].push_back(corr*evt.ach[i]);
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
         double phi = evt.phi[p];
         if (i==1&&phi>TMath::Pi()) phi -= 2*TMath::Pi();
         if (i==9&&phi<TMath::Pi()) phi += 2*TMath::Pi();
         if (stage==0) {
            amp_2d[i]->Fill(zc,phi,ach[i].back());
            len_2d[i]->Fill(zc,phi,1.0/std::sin(evt.theta[p]));
         } else {
            amp_2d_gx[i]->Fill(zc,phi,ach[i].back());
            eff_2d_gx[i]->Fill(zc,phi,thereAreHits);
            Int_t ix = amp_2d[i]->GetXaxis()->FindBin(zc);
            Int_t iy = amp_2d[i]->GetYaxis()->FindBin(phi);
            double alpha = 
               amp_2d[i]->GetBinContent(ix,iy)/
               len_2d[i]->GetBinContent(ix,iy);
            //double effc = 1-exp(-ach[i].back());
            double effc = 1-exp(-alpha/std::sin(evt.theta[p]));
            effc_2d_gx[i]->Fill(zc,phi,effc);
            if (mainGate) {
               amp_2d_g0[i]->Fill(zc,phi,ach[i].back());
               eff_2d_g0[i]->Fill(zc,phi,thereAreHits);
               effc_2d_g0[i]->Fill(zc,phi,effc);
            }
            if (phiSelected&&thetaSelected) {
               if (thereAreHits) time_spect[i]->Fill(tchr[i].back());
               amp_spect_pds[i]->Fill(ach[i].front());
               amp_spect_gx[i]->Fill(ach[i].back());
               amp_gx->Fill(i,ach[i].back());
               eff_gx->Fill(i,thereAreHits);
               effc_gx->Fill(i,effc);
               if (mainGate) {
                  amp_spect_g0[i]->Fill(ach[i].back());
                  amp_g0->Fill(i,ach[i].back());
                  eff_g0->Fill(i,thereAreHits);
                  effc_g0->Fill(i,effc);
               }
            }
            if (phiSelected) {
               amp_z0_gx[i]->Fill(zc,ach[i].back());
               if (mainGate) amp_z0_g0[i]->Fill(zc,ach[i].back());
            }
            if (thetaSelected) {
               amp_phi_gx[i]->Fill(phi,ach[i].back());
               if (mainGate) amp_phi_g0[i]->Fill(phi,ach[i].back());
            }
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

void ampspect(TString scan = "mhad2017", Double_t corr = 1) {
   gROOT->Reset();
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   //
   std::vector<AccProxy> accproxy = getAccProxy(19999,40272,false,"CURRENT");
   //
   Event evt;
   //TChain ch("h1");
   //ch.Add(Form("/work/users/konctbel/KKc/exp/%s/*.root",scan.Data()));
   TChain ch("t1");
   ch.Add("/work/users/konctbel/calibs/R007-001/output/ntuples/col/*.root");
   //ch.Add("/work/users/konctbel/calibs/R007-001/output/ntuples/col/exp2012_[6..9]60_col.root");
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
   for (Int_t stage=0;stage<2;stage++) {
      for (Int_t jentry=0;jentry<nevent;jentry++) {   
         Long64_t ientry = ch.LoadTree(jentry);
         if (ientry < 0)
            break;
         ch.GetEntry(jentry);
         evt.setAccProxy(accproxy);
         if (!(jentry%dnevent)) std::cout << jentry/dnevent << "% of events processed " << std::endl;
         /*if (!runtime[evt.run]) runtime[evt.run] = ftime(evt.run);
         evt.rtime = runtime[evt.run];
         // split events by the run
         if(!nevtr[evt.run]) {
            nevtr[evt.run]++;
            chit[evt.run] = ChHit(evt.run,"run");
            chit[evt.run].init();
            chit[evt.run].setCorr(corr);
         }
         chit[evt.run].fill(evt,stage);
         chit[evt.run].nevt++;
         // split events y the day
         if (!runday[evt.run]) runday[evt.run] = fdate(evt.run);
         size_t date = runday[evt.run];
         if(!nevtrD[date]) {
            nevtrD[date]++;
            chitD[date] = ChHit(date,"day");
            chitD[date].init();
            chitD[date].setCorr(corr);
         }
         chitD[date].fill(evt,stage);
         chitD[date].nevt++;*/
         // split events by the beam energy
         size_t beam10 = (size_t)(10*evt.beam);
         if(!nevtrE[beam10]) {
            nevtrE[beam10]++;
            chitE[beam10] = ChHit(beam10,"beam");
            chitE[beam10].init();
            chitE[beam10].setCorr(corr);
         }
         chitE[beam10].fill(evt,stage);
         chitE[beam10].nevt++;
      }
   }
   TString fname = Form("acc_%s_corr%lf_new.root",scan.Data(),corr);
   std::map<int,ChHit>::iterator it;
   double cnt = 0;
   /*for(it = chit.begin(); it != chit.end(); it++) {
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
   }*/
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
         stop();
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

Double_t gp2x(Double_t *x, Double_t *par)
{
   double xr = par[0];
   double sr = par[1];
   double a = par[4]/2;
   double b = par[3]-2*a*xr;
   double c = par[2]-a*xr*xr-b*xr;
   double dx = (x[0]-xr)/sr;
   double f = (c+b*x[0]+a*(x[0]*x[0]+sr*sr))*(1+erf(dx/sqrt(2.0)))/2+
      (b+a*(x[0]+xr))*sr/sqrt(2*3.1415927)*exp(-dx*dx/2);
   return f;
}

Double_t gg(Double_t *x, Double_t *par)
{
   double a1 = par[0];
   double m1 = par[1];
   double s1 = par[2];
   double dx1 = (x[0]-m1)/s1;
   double f1 = a1*exp(-dx1*dx1/2);///std::sqrt(2*TMath::Pi())/s1;
   double a2 = par[3];
   double m2 = par[4];
   double s2 = par[5];
   double dx2 = (x[0]-m2)/s2;
   double f2 = a2*exp(-dx2*dx2/2);///std::sqrt(2*TMath::Pi())/s2;
   return f1+f2;
}

class PhiProfile
{
 public:
   PhiProfile(){h=NULL;};
   void phi_profile_fit();
   void phi_profile_fits();
   void phi_profile_plot();
   void add(PhiProfile h);
   void readStart();
   TString type;
   size_t t;
   size_t dt;
   size_t tmin;
   size_t tmax;
   size_t ind;
   size_t counter;
   TProfile* h;
   std::vector<Double_t> p;
   std::vector<Double_t> dp;
   TString fstart;
   std::vector<Double_t> *ps;
   std::vector<Double_t> *dps;
};

void PhiProfile::readStart()
{
   TFile *f = new TFile(fstart);
   ps = (std::vector<Double_t> *)f->Get(Form("par_%s_counter%d",type.Data(),counter));
   dps = (std::vector<Double_t> *)f->Get(Form("dpar_%s_counter%d",type.Data(),counter));
   f->Close();
}

void PhiProfile::add(PhiProfile p)
{
   type = p.type;
   ind = p.ind;
   counter = p.counter;
   if (h) h->Add(p.h);
   else h = p.h;
}

void PhiProfile::phi_profile_fits()
{
   readStart();
   double scal = TMath::Pi()/180;
   double dphi = 2.4*scal;
   //
   TF1* fn = new TF1("scphi",scphi,ps->at(12)-15*scal,ps->at(15)+15*scal,20);
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
   fn->SetParameters(ps->data());
   fn->SetLineColor(kRed);
   fn->SetNpx(1000);
   //
   for(size_t i=0;i<12;i++)
      fn->FixParameter(i,fn->GetParameter(i));
   //fn->FixParameter(12,fn->GetParameter(12));
   //fn->FixParameter(13,fn->GetParameter(13));
   fn->FixParameter(14,fn->GetParameter(14));
   //fn->FixParameter(15,fn->GetParameter(15));
   fn->FixParameter(17,fn->GetParameter(17));
   fn->FixParameter(18,fn->GetParameter(18));
   fn->FixParameter(19,fn->GetParameter(19));
   fn->SetParLimits(12,fn->GetParameter(12)-dphi,fn->GetParameter(12)+dphi);
   fn->SetParLimits(13,fn->GetParameter(13)-dphi,fn->GetParameter(13)+dphi);
   fn->SetParLimits(15,fn->GetParameter(15)-dphi,fn->GetParameter(15)+dphi);
   h->Fit("scphi");
   fn->FixParameter(12,fn->GetParameter(12));
   fn->FixParameter(13,fn->GetParameter(13));
   fn->FixParameter(15,fn->GetParameter(15));
   for(size_t i=0;i<12;i++)
      fn->ReleaseParameter(i);
   h->Fit("scphi");
   fn->ReleaseParameter(12);
   fn->ReleaseParameter(13);
   fn->ReleaseParameter(15);
   fn->SetParLimits(12,fn->GetParameter(12)-dphi,fn->GetParameter(12)+dphi);
   fn->SetParLimits(13,fn->GetParameter(13)-dphi,fn->GetParameter(13)+dphi);
   fn->SetParLimits(15,fn->GetParameter(15)-dphi,fn->GetParameter(15)+dphi);
   h->Fit("scphi");
   fn->ReleaseParameter(17);
   fn->ReleaseParameter(18);
   fn->ReleaseParameter(19);
   fn->SetParLimits(17,0.1*scal,0.9*scal);
   fn->SetParLimits(18,0.1*scal,0.9*scal);
   fn->SetParLimits(19,0.1*scal,0.9*scal);
   h->Fit("scphi");
   for(size_t i=0;i<2;i++) {
      double xmin = fn->GetParameter(12)-3*fn->GetParameter(17);
      double xmax = fn->GetParameter(15)+3*fn->GetParameter(19);
      h->Fit("scphi","","",xmin,xmax);
   }
   p.resize(0);
   dp.resize(0);
   for(size_t i=0;i<fn->GetNpar();i++) {
      p.push_back(fn->GetParameter(i));
      dp.push_back(fn->GetParError(i));
   }
}

void PhiProfile::phi_profile_fit()
{
   double scal = TMath::Pi()/180;
   //
   double dphi = 2.4*scal;
   double xmin = (counter*40.-35.)*scal+dphi;
   double xmax = (counter*40.-15.)*scal+dphi;
   TF1 *fxc = new TF1("fxc",gg,xmin,xmax,6);
   double xm = h->GetBinCenter(h->GetMaximumBin());
   fxc->SetParameter(0,h->GetMaximum());
   fxc->SetParameter(1,xm);
   fxc->SetParameter(2,0.5*scal);
   fxc->SetParameter(3,5);
   fxc->SetParameter(4,xm);
   fxc->SetParameter(5,10*scal);
   fxc->FixParameter(0,h->GetMaximum());
   fxc->FixParameter(1,xm);
   fxc->FixParameter(2,0.5*scal);
   fxc->SetParLimits(4,xm-dphi,xm+dphi);
   fxc->SetParLimits(5,5*scal,20*scal);
   h->Fit(fxc,"","",xmin,xmax);
   fxc->ReleaseParameter(0);
   fxc->ReleaseParameter(1);
   fxc->ReleaseParameter(2);
   fxc->SetParLimits(0,0,100);
   fxc->SetParLimits(1,(counter*40.-30.)*scal+dphi,(counter*40.-20.)*scal+dphi);
   fxc->SetParLimits(2,0.1*scal,0.6*scal);
   h->Fit(fxc,"","",xmin,xmax);
   //stop();
   //
   xmin = ((counter-1)*40.-10.)*scal+dphi;
   xmax = fxc->GetParameter(1)-3*fxc->GetParameter(2);//((counter-1)*40.+10.)*scal+dphi;
   TF1 *fxl = new TF1("fxl",gp2x,xmin,xmax,5);
   fxl->SetParameter(0,((counter-1)*40+2.5)*scal);
   fxl->SetParameter(1,0.5*scal);
   fxl->SetParameter(2,5);
   fxl->SetParameter(3,0);
   fxl->SetParameter(4,0);
   h->Fit(fxl,"","",xmin,xmax);
   //stop();
   xmin = fxc->GetParameter(1)+3*fxc->GetParameter(2);//(counter*40.-10.)*scal+dphi;
   xmax = (counter*40.+10.)*scal+dphi;
   TF1 *fxr = new TF1("fxr",gp2x,xmin,xmax,5);
   fxr->SetParameter(0,(counter*40+2.5)*scal);
   fxr->SetParameter(1,-0.5*scal);
   fxr->SetParameter(2,5);
   fxr->SetParameter(3,0);
   fxr->SetParameter(4,0);
   h->Fit(fxr,"","",xmin,xmax);
   //stop();
   //
   std::vector<double> par(20);
   //
   par[12] = fxl->GetParameter(0);
   par[13] = fxc->GetParameter(1);
   par[14] = par[13]+12.5*scal;
   par[15] = fxr->GetParameter(0);
   //
   xmin = par[12];
   xmax = par[13];
   par[0] = 2*fxl->Eval(xmin);
   par[1] = fxl->Eval(xmin+0.33*(xmax-xmin));
   par[2] = fxl->Eval(xmax-0.33*(xmax-xmin));
   par[3] = fxl->Eval(xmax);
   //
   xmin = par[13];
   xmax = par[14];
   par[4] = fxr->Eval(xmin);
   par[5] = fxr->Eval(xmin+0.33*(xmax-xmin));
   par[6] = fxr->Eval(xmax-0.33*(xmax-xmin));
   par[7] = fxr->Eval(xmax);
   //
   xmin = par[14];
   xmax = par[15];
   par[8] = fxr->Eval(xmin);
   par[9] = fxr->Eval(xmin+0.33*(xmax-xmin));
   par[10] = fxr->Eval(xmax-0.33*(xmax-xmin));
   par[11] = 2*fxr->Eval(xmax);
   //
   par[16] = fxc->GetParameter(0);
   par[17] = fxl->GetParameter(1);
   par[18] = fxc->GetParameter(2);
   par[19] = fabs(fxr->GetParameter(1));
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
   fn->FixParameter(12,fn->GetParameter(12));
   fn->FixParameter(13,fn->GetParameter(13));
   fn->FixParameter(14,fn->GetParameter(14));
   fn->FixParameter(15,fn->GetParameter(15));
   fn->FixParameter(17,fn->GetParameter(17));
   fn->FixParameter(18,fn->GetParameter(18));
   fn->FixParameter(19,fn->GetParameter(19));
   h->Fit("scphi");
   //stop();
   for(size_t i=0;i<3;i++)
      fn->ReleaseParameter(i);
   h->Fit("scphi");
   //stop();
   for(size_t i=3;i<8;i++)
      fn->ReleaseParameter(i);
   h->Fit("scphi");
   //stop();
   for(size_t i=8;i<12;i++)
      fn->ReleaseParameter(i);
   h->Fit("scphi");
   //stop();
   fn->ReleaseParameter(12);
   fn->ReleaseParameter(13);
   fn->ReleaseParameter(15);
   fn->SetParLimits(12,fn->GetParameter(12)-dphi,fn->GetParameter(12)+dphi);
   fn->SetParLimits(13,fn->GetParameter(13)-dphi,fn->GetParameter(13)+dphi);
   fn->SetParLimits(15,fn->GetParameter(15)-dphi,fn->GetParameter(15)+dphi);
   h->Fit("scphi");
   //stop();
   fn->ReleaseParameter(17);
   fn->ReleaseParameter(18);
   fn->ReleaseParameter(19);
   fn->SetParLimits(17,0.1*scal,0.9*scal);
   fn->SetParLimits(18,0.1*scal,0.9*scal);
   fn->SetParLimits(19,0.1*scal,0.9*scal);
   h->Fit("scphi");
   //stop();
   for(size_t i=0;i<2;i++) {
      xmin = fn->GetParameter(12)-3*fn->GetParameter(17);
      xmax = fn->GetParameter(15)+3*fn->GetParameter(19);
      h->Fit("scphi","","",xmin,xmax);
   }
   p.resize(0);
   dp.resize(0);
   for(size_t i=0;i<fn->GetNpar();i++) {
      p.push_back(fn->GetParameter(i));
      dp.push_back(fn->GetParError(i));
   }
}

void PhiProfile::phi_profile_plot()
{
   h->Draw();
   double scal = TMath::Pi()/180;
   for(size_t c=0;c<3;c++) {
      TF1* fna = new TF1(Form("aerogel%d",c),scphi,p[12]-15*scal,p[15]+15*scal,20);
      fna->SetParameters(p.data());
      for(size_t i=0;i<12;i++)
         fna->SetParameter(i,0);
      for(size_t i=c*4;i<c*4+4;i++)
         fna->SetParameter(i,p[i]);
      fna->SetParameter(16,0);
      fna->SetLineColor(kBlue);
      fna->SetLineWidth(0);
      fna->SetNpx(1000);
      fna->Draw("same");
   }
   //
   TF1* fns = new TF1("shifter",scphi,p[12]-15*scal,p[15]+15*scal,20);
   fns->SetParameters(p.data());
   for(size_t i=0;i<12;i++)
      fns->SetParameter(i,0);
   fns->SetLineColor(kGreen);
   fns->SetLineWidth(0);
   fns->SetNpx(1000);
   fns->Draw("same");
   TF1* f = new TF1("f",scphi,p[12]-15*scal,p[15]+15*scal,20);
   f->SetParameters(p.data());
   h->GetXaxis()->SetRangeUser(p[12]-5*scal,p[15]+5*scal);
   h->GetYaxis()->SetRangeUser(0,1.1*std::max(h->GetMaximum(),
                                              f->GetMaximum()));
   TText *t = new TText(p[12],h->GetMaximum(),
                        Form("Beam %6.1f  Counter %d",
                             ((double)ind)/10,counter));
   t->Draw();
   //
   //stop();
}

void scphi_fit()
{
   std::map< size_t,std::vector< std::vector<Double_t> > > pars;
   std::map< size_t,std::vector< std::vector<Double_t> > > dpars;
   std::vector<PhiProfile> hp;
   std::map<size_t,PhiProfile> hs;
   TFile *f = new TFile("acc_mhad2017.root");
   TIter next_dir(f->GetListOfKeys());
   TKey *key_dir;
   TFile *fo = new TFile("mhad2017_phi_fit.root","RECREATE");
   while ((key_dir=(TKey*)next_dir())) {
      if (!TString("TDirectoryFile").BeginsWith(key_dir->GetClassName())) continue;
      TDirectory *d = (TDirectory*)f->Get(key_dir->GetName());
      TIter next(d->GetListOfKeys());
      TKey *key;
      while ((key=(TKey*)next())) {
         TString name(key->GetName());
         TString type = "beam";
         if (name.BeginsWith(type)&&
             name.Contains("_amp_phi_gx_counter")) {
            size_t ind, counter;
            sscanf(name.Data(),Form("%s%s_amp_phi_gx_counter%s",type.Data(),"%d","%d"),&ind,&counter);
            //
            std::vector<size_t> *t = (std::vector<size_t> *)d->Get(Form("time_%s%08d",type.Data(),ind));
            std::cout << ind << " " << (*t)[0] << " " << (*t)[1] << std::endl;
            if ((*t)[0]<fdatetotime(20170309)) continue;
            //if (ind!=7750) continue;
            //if (counter<8) continue;
            //
            hp.push_back(PhiProfile());
            hp.back().t = t->size()<1 ? 0 : (*t)[0];
            hp.back().dt = t->size()<2 ? 0 : (*t)[1];
            hp.back().tmin = t->size()<3 ? 0 : (*t)[2];
            hp.back().tmax = t->size()<4 ? 0 : (*t)[3];
            hp.back().type = type;
            hp.back().ind = ind;
            hp.back().counter = counter;
            hp.back().fstart = "mhad2017_phi_profiles.root";
            hp.back().h = (TProfile*)d->Get(Form("%s%08d_amp_phi_gx_counter%d",type.Data(),ind,counter));
            //
            hs[counter].add(hp.back());
            //
            hp.back().phi_profile_fit();
            hp.back().phi_profile_plot();
            //
            hp.back().h->Write();
            fo->WriteObjectAny(&hp.back().p,"std::vector<Double_t>",Form("par_%s%08d_counter%d",type.Data(),ind,counter));
            fo->WriteObjectAny(&hp.back().dp,"std::vector<Double_t>",Form("dpar_%s%08d_counter%d",type.Data(),ind,counter));
            //
            pars[counter].push_back(hp.back().p);
            dpars[counter].push_back(hp.back().dp);
         }
      }
   }
   fo->Close();
   TFile *fs = new TFile("mhad2017_phi_profiles.root","RECREATE");
   for(size_t c=1;c<10;c++) {
      hs[c].phi_profile_fit();
      hs[c].phi_profile_plot();
      hs[c].h->Write();
      fs->WriteObjectAny(&hs[c].p,"std::vector<Double_t>",Form("par_%s_counter%d",hs[c].type.Data(),hs[c].counter));
      fs->WriteObjectAny(&hs[c].dp,"std::vector<Double_t>",Form("dpar_%s_counter%d",hs[c].type.Data(),hs[c].counter));
      stop();
   }
   fs->Close();
   for(size_t c=1;c<10;c++) {
      std::vector<double> x;
      std::vector<double> dx;
      std::vector<double> y;
      std::vector<double> dy;
      size_t n = pars[c].size();
      for(size_t i=0;i<n;i++) {
         x.push_back(i+1);
         dx.push_back(0);
         y.push_back(pars[c][i][12]);
         dy.push_back(dpars[c][i][17]);
      }
      TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
      gr->SetMarkerSize(1);
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(kBlack);
      gr->Draw("AP");
      stop();      
   }
   for(size_t c=1;c<10;c++) {
      std::vector<double> x;
      std::vector<double> dx;
      std::vector<double> y;
      std::vector<double> dy;
      size_t n = pars[c].size();
      for(size_t i=0;i<n;i++) {
         x.push_back(i+1);
         dx.push_back(0);
         y.push_back(pars[c][i][13]);
         dy.push_back(dpars[c][i][18]);
      }
      TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
      gr->SetMarkerSize(1);
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(kBlack);
      gr->Draw("AP");
      stop();      
   }
   for(size_t c=1;c<10;c++) {
      std::vector<double> x;
      std::vector<double> dx;
      std::vector<double> y;
      std::vector<double> dy;
      size_t n = pars[c].size();
      for(size_t i=0;i<n;i++) {
         x.push_back(i+1);
         dx.push_back(0);
         y.push_back(pars[c][i][15]);
         dy.push_back(dpars[c][i][19]);
      }
      TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
      gr->SetMarkerSize(1);
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(kBlack);
      gr->Draw("AP");
      stop();      
   }
}



double eff_amp(double *x, double *par)
{
   double mr = x[0];
   double D = 0.2;
   double d = par[0];
   double mu = mr;
   for(size_t i=0;i<5;i++)
      mu = mr/(1-D*d/2*exp(-mu));
   return 1-(1-mu*d)*exp(-mu);
}

double eff_amp_2(double *x, double *par)
{
   double mr = x[0]*par[0];
   double D = 0.2;
   double d = par[1];
   double mu = mr;
   for(size_t i=0;i<5;i++)
      mu = mr/(1-D*d/2*exp(-mu));
   //return 1-((1-mu*d)*exp(-mu)+par[2]);
   return 1-par[2]-exp(-x[0]*par[0]);
}


class AccMapCal
{
 public:
   AccMapCal(){h_amp = NULL; h_eff = NULL;};
   void add(TProfile2D* h_amp_, TProfile2D* h_eff_);
   void add(AccMapCal &map);
   void setTime(std::vector<size_t> *t);
   void get();
   void fit();
   void plot();
   TString type;
   size_t t;
   size_t dt;
   size_t tmin;
   size_t tmax;
   size_t ind;
   size_t counter;
   TProfile2D* h_amp;
   TProfile2D* h_eff;
   std::vector<Double_t> th;
   std::vector<Double_t> amp;
   std::vector<Double_t> eff;
   std::vector<Double_t> dth;
   std::vector<Double_t> damp;
   std::vector<Double_t> deff;
   std::vector<Double_t> thf;
   std::vector<Double_t> ampf;
   std::vector<Double_t> efff;
   std::vector<Double_t> dthf;
   std::vector<Double_t> dampf;
   std::vector<Double_t> defff;
   TF1 *f;
   TGraphErrors *gr;
   TGraphErrors *grf;
   bool goodData;
   Double_t mean;
   Double_t dmean;
};

void AccMapCal::add(TProfile2D* h_amp_, TProfile2D* h_eff_)
{
   if (!h_eff_) return;
   if (!h_amp_) return;
   if (h_amp) h_amp->Add(h_amp_);
   else h_amp = h_amp_;
   if (h_eff) h_eff->Add(h_eff_);
   else h_eff = h_eff_;
}

void AccMapCal::add(AccMapCal &map)
{
   std::cout << "AccMapCal::add: " << h_amp << " " << h_eff << " " << map.h_amp << " " << map.h_eff << std::endl;
   if (!map.h_amp) return;
   if (!map.h_eff) return;
   if (h_amp) h_amp->Add(map.h_amp);
   else h_amp = map.h_amp;
   if (h_eff) h_eff->Add(map.h_eff);
   else h_eff = map.h_eff;
}

void AccMapCal::setTime(std::vector<size_t> *v)
{
   t = v->size()<1 ? 0 : (*v)[0];
   dt = v->size()<2 ? 0 : (*v)[1];
   tmin = v->size()<3 ? 0 : (*v)[2];
   tmax = v->size()<4 ? 0 : (*v)[3];
}

void AccMapCal::get()
{
   std::vector<int> ie;
   ie.push_back(1);
   ie.push_back(5);
   ie.push_back(6);
   ie.push_back(13);
   double mean = 0;
   TH1D *hx = h_eff->ProjectionX();
   for(int i=15;i<=35;i++) {//h_amp->GetNbinsX();i++) {
      for(int j=1;j<=h_amp->GetNbinsY();j++) {
         if (std::find(ie.begin(),ie.end(),j)!=ie.end()) continue;
         if (!h_eff->GetBinError(i,j)) continue;
         th.push_back(hx->GetBinCenter(i));
         dth.push_back(0);
         amp.push_back(h_amp->GetBinContent(i,j));
         damp.push_back(h_amp->GetBinError(i,j));
         eff.push_back(h_eff->GetBinContent(i,j));
         deff.push_back(h_eff->GetBinError(i,j));
         mean += amp.back();
      }
   }
   mean /= amp.size();
   gr = new TGraphErrors(amp.size(),amp.data(),eff.data(),damp.data(),deff.data());
   goodData = amp.size()>10 && mean>0.2 ? true : false; 
}

void AccMapCal::fit()
{
 //  if (!goodData) return;
   std::cout << h_eff << " " << h_amp << std::endl;
   if (!h_eff) return;
   if (!h_amp) return;
   get();
   TString fname(Form("f_%d_%d",ind,counter));
   f = new TF1(fname,eff_amp_2,0,20,3);
   f->SetParName(0,"k");
   f->SetParName(1,"d");
   f->SetParName(2,"s");
   f->SetParameter(0,1);
   f->SetParameter(1,0);
   f->SetParameter(2,0);
   f->SetParLimits(0,0,2);
   f->SetParLimits(1,0,1);
   f->SetParLimits(2,0,1);
   f->FixParameter(1,0);
   //f->FixParameter(2,0);
   //
   double ns = 1000;
   std::vector<Double_t> x;
   std::vector<Double_t> y;
   std::vector<Double_t> z;
   std::vector<Double_t> dx;
   std::vector<Double_t> dy;
   std::vector<Double_t> dz;
   for(size_t s=0;s<3;s++) {
      x.resize(0);
      y.resize(0);
      z.resize(0);
      dx.resize(0);
      dy.resize(0);
      dz.resize(0);
      TH1D *hx = new TH1D("hx","hx",1000,0,100);
      for(size_t i=0;i<amp.size();i++)
         if (std::fabs(eff[i] - f->Eval(amp[i])) < ns*std::max(0.01,deff[i])) {
            x.push_back(amp[i]);
            dx.push_back(damp[i]);
            y.push_back(eff[i]);
            dy.push_back(deff[i]);
            z.push_back(th[i]);
            dz.push_back(dth[i]);
            hx->Fill(amp[i]);
         }
      grf = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
      grf->Fit(fname);
      ns = 3;
      mean = hx->GetMean();
      dmean = hx->GetMeanError();
   }
}

void AccMapCal::plot()
{
   std::cout << gr << " " << grf << std::endl;
   if (type=="beam") gr->SetTitle(Form("Beam:%6.1f  counter:%d  %s",(double)ind/10,counter,ftimetodate(t).Data()));
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kGray);
   gr->SetLineColor(kGray);
   gr->Draw("AP");
   //
   grf->SetMarkerStyle(20);
   grf->SetMarkerColor(kBlack);
   grf->SetLineColor(kBlack);
   grf->Draw("P");
   //
   TF1 *f0 = new TF1("f0","1-exp(-x)",0,20);
   f0->SetLineColor(kBlue);
   f0->Draw("same");
   //
   f->SetLineColor(kRed);
   f->Draw("same");
   //
   gr->GetXaxis()->SetTitle("amplitude (pe)");
   gr->GetYaxis()->SetTitle("efficiency");
   //
}

void amp_eff()
{
   std::vector<Double_t> c0(9);
   for(size_t i=0;i<9;i++) c0[i] = i+1;
   std::vector<Double_t> dc0(9,0);
   std::vector<Double_t> k0(9);
   std::vector<Double_t> dk0(9);
   std::vector<Double_t> e0(9);
   std::vector<Double_t> de0(9);
   double scal = TMath::Pi()/180;
   //
   std::vector<AccMapCal> map;
   //
   TFile *f = new TFile("acc_mhad2017.root");
   //TFile *f = new TFile("acc_mhad2017_corr0.800000.root");
   //TFile *f = new TFile("acc_mhad2017_corr1.000000_new.root");
   //TFile *f = new TFile("acc_mhad2017_corr1.100000.root");
   TIter next_dir(f->GetListOfKeys());
   TKey *key_dir;
   while ((key_dir=(TKey*)next_dir())) {
      if (!TString("TDirectoryFile").BeginsWith(key_dir->GetClassName())) continue;
      TDirectory *d = (TDirectory*)f->Get(key_dir->GetName());
      TIter next(d->GetListOfKeys());
      TKey *key;
      while ((key=(TKey*)next())) {
         TString name(key->GetName());
         TString type = "beam";
         if (name.BeginsWith(type)&&
             name.Contains("_amp_2d_g0_counter")) {
            size_t ind, counter;
            sscanf(name.Data(),Form("%s%s_amp_2d_g0_counter%s",type.Data(),"%d","%d"),&ind,&counter);
            //
            map.push_back(AccMapCal());
            map.back().type = type;
            map.back().ind = ind;
            map.back().counter = counter;
            map.back().setTime((std::vector<size_t> *)d->Get(Form("time_%s%08d",type.Data(),ind)));
            //
            map.back().add((TProfile2D*)d->Get(Form("%s%08d_amp_2d_gx_counter%d",type.Data(),ind,counter)),
                           (TProfile2D*)d->Get(Form("%s%08d_eff_2d_gx_counter%d",type.Data(),ind,counter)));
            //
            map.back().fit();
            //map.back().plot();
            //
            //stop();
         }
      }
   }
   std::vector<size_t> ie;
   ie.push_back(8510);
   ie.push_back(8502);
   ie.push_back(8500);
   ie.push_back(8600);
   ie.push_back(8700);
   ie.push_back(9500);
   ie.push_back(9800);
   ie.push_back(9360);
   ie.push_back(9388);
/*   for(size_t c=1;c<10;c++) {
      std::vector<Double_t> x;
      std::vector<Double_t> y;
      std::vector<Double_t> dx;
      std::vector<Double_t> dy;
      for(size_t i=0;i<map.size();i++) {
         if (std::find(ie.begin(),ie.end(),map[i].ind)!=ie.end()) continue;
         if (map[i].counter!=c) continue;
         if (!map[i].goodData) continue;
         x.push_back(map[i].t);
         dx.push_back(0*map[i].dt);
         y.push_back(map[i].f->GetParameter(0));
         dy.push_back(map[i].f->GetParError(0));
         map[i].plot();
         stop();
      }
      TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
      gr->GetXaxis()->SetTimeDisplay(1);
      gr->SetMarkerSize(1);
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(kBlack);
      gr->Draw("AP");
      stop();
   }*/
   //
   std::vector< std::vector<AccMapCal> > mapc(2,std::vector<AccMapCal>(9));
   for(size_t i=0;i<map.size();i++) {
      if (std::find(ie.begin(),ie.end(),map[i].ind)!=ie.end()) continue;
      if (!map[i].goodData) continue;
      size_t im = map[i].t<fdatetotime(20170309) ? 0 : 1;
      size_t ic = map[i].counter-1;
      if (mapc[im][ic].h_amp) {
         mapc[im][ic].type = "xxx";
         mapc[im][ic].ind = 0;
         mapc[im][ic].counter = map[i].counter;
      }
      mapc[im][ic].add(map[i]);
      std::cout << ">>>" << im << " " << ic << mapc[im][ic].type << " " << mapc[im][ic].h_amp << " " << mapc[im][ic].h_eff << std::endl;
   }
   size_t im = 0;
   size_t ic = 0;
   std::cout << ">>>" << mapc[im][ic].type << " " << mapc[im][ic].h_amp << " " << mapc[im][ic].h_eff << std::endl;
   std::vector<Double_t> nc(9);
   std::vector<Double_t> dnc(9);
   std::vector<Double_t> meanct;
   std::vector<Double_t> dmeanct;
   std::vector<Double_t> ampct;
   std::vector<Double_t> dampct;
   std::vector< std::vector<Double_t> > meanc(2,std::vector<Double_t>(9));
   std::vector< std::vector<Double_t> > dmeanc(2,std::vector<Double_t>(9));
   std::vector< std::vector<Double_t> > ampc(2,std::vector<Double_t>(9));
   std::vector< std::vector<Double_t> > dampc(2,std::vector<Double_t>(9));
   std::vector< std::vector<Double_t> > effc(2,std::vector<Double_t>(9));
   std::vector< std::vector<Double_t> > deffc(2,std::vector<Double_t>(9));
   for(size_t im=0;im<mapc.size();im++) {
      for(size_t ic=0;ic<9;ic++) {
         std::cout << ">>>" << mapc[im][ic].type << " " << mapc[im][ic].h_amp << " " << mapc[im][ic].h_eff << std::endl;
         mapc[im][ic].fit();
         mapc[im][ic].plot();
         nc[ic] = ic+1;
         dnc[ic] = 0;
         ampc[im][ic] = mapc[im][ic].f->GetParameter(0);
         dampc[im][ic] = mapc[im][ic].f->GetParError(0);
         effc[im][ic] = mapc[im][ic].f->GetParameter(2);
         deffc[im][ic] = mapc[im][ic].f->GetParError(2);
         meanc[im][ic] = mapc[im][ic].mean;
         dmeanc[im][ic] = mapc[im][ic].dmean;
         ampct.push_back(ampc[im][ic]);
         dampct.push_back(dampc[im][ic]);
         meanct.push_back(meanc[im][ic]);
         dmeanct.push_back(dmeanc[im][ic]);
         stop();
      }
      TGraphErrors *gra = new TGraphErrors(nc.size(),nc.data(),ampc[im].data(),dnc.data(),dampc[im].data());
      gra->SetMarkerSize(1);
      gra->SetMarkerStyle(20);
      gra->SetMarkerColor(kBlack);
      gra->Draw("AP");
      stop();
      TGraphErrors *gre = new TGraphErrors(nc.size(),nc.data(),effc[im].data(),dnc.data(),deffc[im].data());
      gre->SetMarkerSize(1);
      gre->SetMarkerStyle(20);
      gre->SetMarkerColor(kBlack);
      gre->Draw("AP");
      stop();
      TGraphErrors *grm = new TGraphErrors(meanc[im].size(),meanc[im].data(),ampc[im].data(),dmeanc[im].data(),dampc[im].data());
      grm->SetMarkerSize(1);
      grm->SetMarkerStyle(20);
      grm->SetMarkerColor(kBlack);
      grm->Draw("AP");
      stop();
   }
   TGraphErrors *grm = new TGraphErrors(meanct.size(),meanct.data(),ampct.data(),dmeanct.data(),dampct.data());
   grm->SetMarkerSize(1);
   grm->SetMarkerStyle(20);
   grm->SetMarkerColor(kBlack);
   grm->Draw("AP");
   stop();
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

//-----------------------------------------------------

Double_t g2(Double_t *x, Double_t *par)
{
   double arg = (x[0] - par[1])/par[2];
   double f = par[0]*exp(-arg*arg/2)/std::sqrt(2*TMath::Pi())/par[2];
   arg = (x[0] - par[1])/par[4];
   f += par[3]*exp(-arg*arg/2)/std::sqrt(2*TMath::Pi())/par[4];
   return f;
}

std::vector<TFitResultPtr> PdsSpectFitting(std::vector<TH1D*> &pdsSpect)
{
   std::vector<TFitResultPtr> fitRes;
   Double_t xmin = 20000;
   Double_t xmax = 0;
   for(size_t i=0;i<pdsSpect.size();i++) {
      double mean = pdsSpect.at(i)->GetBinCenter(pdsSpect.at(i)->GetMaximumBin());
      double rms = pdsSpect.at(i)->GetRMS();
      xmin = mean-20*rms<xmin ? mean-20*rms : xmin;
      xmax = mean+20*rms>xmax ? mean+20*rms : xmax;
   }
   TF1 *f = new TF1("f",g2,(size_t)xmin,(size_t)xmax,5);
   f->SetParName(0,"a1");
   f->SetParName(1,"m1");
   f->SetParName(2,"s1");
   f->SetParName(3,"a2");
   f->SetParName(4,"s2");
   f->SetLineColor(kRed);
   f->SetNpx((size_t)(xmax-xmin));
   for(size_t i=0;i<pdsSpect.size();i++) {
      double peak = pdsSpect.at(i)->GetBinCenter(pdsSpect.at(i)->GetMaximumBin());
      double rms = pdsSpect.at(i)->GetRMS();
      double xmin = peak - 3*rms;
      double xmax = peak + 3*rms;
      f->SetParameter(0,pdsSpect.at(i)->GetMaximum());
      f->SetParameter(1,peak);
      f->SetParameter(2,rms/2);
      f->SetParameter(3,0);
      f->SetParameter(4,rms);
      TFitResultPtr r;
      for(size_t s=0;s<3;s++) {
         f->ReleaseParameter(0);
         f->ReleaseParameter(2);
         f->SetParLimits(0,0,1000*pdsSpect.at(i)->GetMaximum());
         f->SetParLimits(1,peak-3*rms,peak+3*rms);
         f->SetParLimits(2,rms/10,3*rms);
         f->FixParameter(3,f->GetParameter(3));
         f->FixParameter(4,f->GetParameter(4));
         r = pdsSpect.at(i)->Fit(f,"SQL","",xmin,xmax);
         //f->FixParameter(0,f->GetParameter(0));
         f->FixParameter(2,f->GetParameter(2));
         f->ReleaseParameter(3);
         f->ReleaseParameter(4);
         f->SetParLimits(3,0,f->GetParameter(0));
         f->SetParLimits(4,f->GetParameter(2),10*f->GetParameter(2));
         r = pdsSpect.at(i)->Fit(f,"SQL","",xmin,xmax);
         //f->ReleaseParameter(0);
         f->ReleaseParameter(2);
         //f->SetParLimits(0,0,10*pdsSpect.at(i)->GetMaximum());
         f->SetParLimits(2,rms/10,3*rms);
         r = pdsSpect.at(i)->Fit(f,"SQL","",xmin,xmax);
         double mean = r->Parameters().at(1);
         double sig = std::max(r->Parameters().at(2),r->Parameters().at(4));
         xmin = mean - 5*sig;
         xmax = mean + 5*sig;
      }
      fitRes.push_back(r);
   }
   return fitRes;
}

//-----------------------------------------------------------------------------

Double_t gp2(Double_t *x, Double_t *par)
{
   double xr = par[0];
   double sr = par[1];
   double a = par[4]/2;
   double b = par[3]-2*a*xr;
   double c = par[2]-a*xr*xr-b*xr;
   double amp = std::fabs(par[5]);
   double dx = (x[0]-xr)/sr;
   double f = amp*exp(-dx*dx/2)/std::sqrt(2*TMath::Pi())/sr;
   f += (c+b*x[0]+a*(x[0]*x[0]+sr*sr))*(1+erf(dx/sqrt(2.0)))/2+
      (b+a*(x[0]+xr))*sr/sqrt(2*3.1415927)*exp(-dx*dx/2);
   return f;
}

Double_t ggpar[5];

Double_t ggp2(Double_t *x, Double_t *par0)
{
   Double_t par[6];
   Double_t k = ggpar[0]/(ggpar[0]+ggpar[3]);
   par[0] = par0[0];
   par[1] = sqrt(ggpar[2]*ggpar[2]+par0[1]*par0[1]);
   par[2] = par0[2]*k;
   par[3] = par0[3]*k;
   par[4] = par0[4]*k;
   par[5] = ggpar[0]*par0[5];
   Double_t f = gp2(x,par);
   par[1] = sqrt(ggpar[4]*ggpar[4]+par0[1]*par0[1]);
   par[2] = par0[2]*(1-k);
   par[3] = par0[3]*(1-k);
   par[4] = par0[4]*(1-k);
   par[5] = ggpar[3]*par0[5];
   f += gp2(x,par);
   return f;
}

//-----------------------------------------------------------------------------

Double_t RooNovosibirskX(Double_t x, Double_t peak, Double_t width, Double_t tail)
{
   if (TMath::Abs(tail) < 1.e-7) {
      return TMath::Exp( -0.5 * TMath::Power( ( (x - peak) / width ), 2 ));
   }
   
   Double_t arg = 1.0 - ( x - peak ) * tail / width;
   
   if (arg < 1.e-7) {
      //Argument of logarithm negative. Real continuation -> function equals zero
      return 0.0;
   }
   
   Double_t log = TMath::Log(arg);
   static const Double_t xi = 2.3548200450309494; // 2 Sqrt( Ln(4) )
   
   Double_t width_zero = ( 2.0 / xi ) * TMath::ASinH( tail * xi * 0.5 );
   Double_t width_zero2 = width_zero * width_zero;
   Double_t exponent = ( -0.5 / (width_zero2) * log * log ) - ( width_zero2 * 0.5 );
   
   return TMath::Exp(exponent) ;
}

Double_t NF(Double_t *x, Double_t *par)
{
      Double_t f = par[0]*RooNovosibirskX(x[0],par[1],par[2],par[3]);
      return f;
}

//-----------------------------------------------------------------------------

std::vector<TFitResultPtr> LedSpectFitting(std::vector<TH1D*> &ledSpect,
                                           const std::vector<TFitResultPtr> &pdsFitRes)
{
   std::vector<TFitResultPtr> fitRes;
   Double_t Xmin = 20000;
   Double_t Xmax = 0;
   for(size_t i=0;i<ledSpect.size();i++) {
      double mean = ledSpect.at(i)->GetBinCenter(ledSpect.at(i)->GetMaximumBin());
      Xmin = mean-100<Xmin ? mean-100 : Xmin;
      Xmax = mean+1000>Xmax ? mean+1000 : Xmax;
   }
   TF1 *f = new TF1("f",ggp2,(size_t)Xmin,(size_t)Xmax,6);
   f->SetParName(0,"xr");
   f->SetParName(1,"sr");
   f->SetParName(2,"a");
   f->SetParName(3,"b");
   f->SetParName(4,"c");
   f->SetParName(5,"amp");
   f->SetParLimits(1,0,10);
   f->SetParLimits(2,0,100000);
   f->SetParLimits(5,0,1);
   f->SetParameters(0,0,10,0,0,1);
   f->SetLineColor(kRed);
   f->SetNpx((size_t)(Xmax-Xmin));
   for(size_t i=0;i<ledSpect.size();i++) {
      for(size_t j=0;j<pdsFitRes.at(i)->NPar();j++) {
         ggpar[j] = pdsFitRes.at(i)->Parameters().at(j);
         std::cout << j << ":" << ggpar[j] << std::endl;
      }
      double sig = std::min(10.,ggpar[0]>ggpar[3] ? ggpar[2] : ggpar[4]);
      double peak = ggpar[1];//ledSpect.at(i)->GetBinCenter(ledSpect.at(i)->GetMaximumBin());
      double rms = ledSpect.at(i)->GetRMS();
      double xmin = peak - 5*sig;
      double xmax = peak + 15*sig;
      double mean = ledSpect.at(i)->GetMean();
      f->SetParameter(0,peak);
      TFitResultPtr r;
      f->SetParLimits(0,peak-3*sig,peak+3*sig);
      f->SetParLimits(1,0,30);
      f->SetParLimits(5,0,1);
      f->FixParameter(2,0);
      f->FixParameter(3,0);
      f->FixParameter(4,0);
      r = ledSpect.at(i)->Fit(f,"SQ","",peak-3*sig,peak+3*sig);
      double mu = -log(r->Parameters().at(5));
      double a1 = (mean-peak)/mu;
      double w = std::max(a1*(1+mu/2),10*sig);
      f->FixParameter(0,r->Parameters().at(0));
      f->FixParameter(1,r->Parameters().at(1));
      f->FixParameter(5,r->Parameters().at(5));
      f->ReleaseParameter(2);
      f->ReleaseParameter(3);
      f->ReleaseParameter(4);
      r = ledSpect.at(i)->Fit(f,"SQ","",peak+5*sig,peak+w);
      if (r) {
         f->FixParameter(2,0);
         f->FixParameter(3,0);
         f->FixParameter(4,0);
      } else {
         f->FixParameter(2,r->Parameters().at(2));
         f->FixParameter(3,r->Parameters().at(3));
         f->FixParameter(4,r->Parameters().at(4));
      }
      f->ReleaseParameter(0);
      f->ReleaseParameter(1);
      f->ReleaseParameter(5);
      f->SetParLimits(0,peak-3*sig,peak+3*sig);
      f->SetParLimits(1,0,30);
      f->SetParLimits(5,0,1);
      r = ledSpect.at(i)->Fit(f,"SQL","",peak-5*sig,peak+w);
      f->ReleaseParameter(2);
      f->ReleaseParameter(3);
      f->ReleaseParameter(4);
      r = ledSpect.at(i)->Fit(f,"SQL","",peak-5*sig,peak+w);
      mu = -log(r->Parameters().at(5));
      a1 = (mean-peak)/mu;
      w = std::max(a1*(1+mu/2),10*sig);
      r = ledSpect.at(i)->Fit(f,"SQL","",peak-5*sig,peak+w);
      if (r->Parameters().at(1)<sig/100||r->Parameters().at(1)<r->Errors().at(1)) {
         f->FixParameter(1,0);
         r = ledSpect.at(i)->Fit(f,"SQL","",peak-5*sig,peak+w);
      }
      //ledSpect[i]->Draw();
      //ledSpect[i]->GetXaxis()->SetRangeUser(xmin-100,xmax+300);
      //stop();
      fitRes.push_back(r);
   }
   return fitRes;
}

struct accproxy_cal
{
   std::vector<Double_t> t0_;
   std::vector<Double_t> dt_;
   std::vector<Double_t> thr_;
   std::vector<Double_t> dt0_;
   std::vector<Double_t> a_;
   std::vector<Double_t> da_;
};

struct accsimpl_cal
{
   std::vector<Double_t> t0;
   std::vector<Double_t> dt0;
   std::vector<Double_t> pds;
   std::vector<Double_t> dpds;
   std::vector<Double_t> a1;
   std::vector<Double_t> da1;
   std::vector<Double_t> p0;
   std::vector<Double_t> dp0;
   std::vector<Double_t> mu;
   std::vector<Double_t> dmu;
   std::vector<Double_t> asym;
   std::vector<Double_t> dasym;
};

accsimpl_cal AccSimplCalculation(const std::vector<TFitResultPtr> &pdsFitRes,
                                 const std::vector<TFitResultPtr> &ledFitRes,
                                 const std::vector<TH1D*> &pdsSpect,
                                 const std::vector<TH1D*> &ledSpect,
                                 const accproxy_cal &accProxyCal)
{
   accsimpl_cal cal;
   for(size_t i=0;i<pdsFitRes.size();i++) {
      //
      cal.t0.push_back(accProxyCal.t0_.at(i));
      cal.dt0.push_back(accProxyCal.dt0_.at(i));
      //
      cal.asym.push_back(accProxyCal.a_.at(i));
      cal.dasym.push_back(accProxyCal.da_.at(i));
      //
      /*double ds = ledFitRes.at(i)->Parameters().at(1);
      double s1 = pdsFitRes.at(i)->Parameters().at(2);
      double s2 = pdsFitRes.at(i)->Parameters().at(4);
      double kA = ledFitRes.at(i)->Parameters().at(5);
      double A1 = pdsFitRes.at(i)->Parameters().at(0);
      double A2 = pdsFitRes.at(i)->Parameters().at(3);
      double n0 = A1*s1+A2*s2;
      double S1 = std::sqrt(s1*s1+ds*ds);
      double S2 = std::sqrt(s2*s2+ds*ds);
      double n1 = kA*(A1*S1+A2*S2);*/
      //
      //double p0 = ledFitRes.at(i)->Parameters().at(5)*ledFitRes.at(i)->Parameters().at(1)/
      //   pdsFitRes.at(i)->Parameters().at(0)/pdsFitRes.at(i)->Parameters().at(2);
      double p0 = ledFitRes.at(i)->Parameters().at(5);//n1/n0;
      cal.p0.push_back(p0);
      //
      double nevt = ledSpect.at(i)->GetEntries();
      double dp0 = ledFitRes.at(i)->Errors().at(5);//std::sqrt(p0*(1-p0)/nevt);
      cal.dp0.push_back(dp0);
      //
      double mu = -log(p0);
      cal.mu.push_back(mu);
      //
      double dmu = dp0/p0;
      cal.dmu.push_back(dmu);
      //
      double a0 = pdsSpect.at(i)->GetMean() - pdsFitRes.at(i)->Parameters().at(1);
      double am = ledSpect.at(i)->GetMean() - ledFitRes.at(i)->Parameters().at(0);
      //std::cout << a0 << " " << am << " " << mu << std::endl;
      double a1 = (am-a0)/mu;
      cal.a1.push_back(a1);
      //
      double da0 = pdsSpect.at(i)->GetMeanError();
      double dam = ledSpect.at(i)->GetMeanError();
      double d1 = a1*dmu/mu;
      double d2 = dam/mu;
      double d3 = da0/mu;
      double da1 = std::sqrt(d1*d1+d2*d2+d3*d3);
      cal.da1.push_back(da1);
      //
      cal.pds.push_back(a0);
      cal.dpds.push_back(da0);
   }
   return cal;
}

void PrintAccSimpl(const accsimpl_cal &cal)
{
   std::cout << "t0(dt0)              pds(pds)              a1(da1)           mu(dmu))           aysm(daysm)" << std::endl;
   for(size_t i=0; i<cal.a1.size(); i++)
      std::cout << cal.t0.at(i) << "(" << cal.dt0.at(i) << ") "
      << cal.pds.at(i) << "(" << cal.dpds.at(i) << ") "
      << cal.a1.at(i) << "(" << cal.da1.at(i) << ") "
      << cal.mu.at(i) << "(" << cal.dmu.at(i) << ") "
      << cal.asym.at(i) << "(" << cal.dasym.at(i) << ") "
      << std::endl;
}

void oldcal(TString fname = "onel_270317_v1.res")
{
   std::vector<TH1D*> pdsSpect;
   std::vector<TH1D*> ledSpect;
   for(size_t i=1;i<10;i++) {
      TString hname = Form("pds%d",i);
      if (gROOT->FindObject(hname.Data())) gROOT->Delete(hname.Data());
      pdsSpect.push_back(new TH1D(hname.Data(),hname.Data(),4000,0,20000));
      hname = Form("led%d",i);
      if (gROOT->FindObject(hname.Data())) gROOT->Delete(hname.Data());
      ledSpect.push_back(new TH1D(hname.Data(),hname.Data(),4000,0,20000));
   }
   //ifstream file("/work/users/konctbel/TESTS/flash/onel_070317_v0.res");
   //ifstream file(Form("/work/users/konctbel/TESTS/flash/%s",fname.Data()));
   fname = Form("/work/users/konctbel/TESTS/fullcal/%s",fname.Data());
   //fname = "/work/users/konctbel/TESTS/fullcal/Simple_Calibration/Cal_170214_1/onel_170214_1.res";
   //fname = "/work/users/konctbel/TESTS/fullcal/Simple_Calibration/Cal_170524_1/onel_170524_1.res";
   //fname = "/work/users/konctbel/TESTS/flash/onel_090617.res";
   //fname = "/work/users/konctbel/TESTS/flash/onel_270317_v5.res";
   TString fnamer = fname(fname.Last('/')+1,fname.Last('.')-fname.Last('/')-1) + ".root";
   ifstream file(fname.Data());
/*   size_t k = 0;
   double x;
   while (1) {
      file >> x;
      if (file.eof()) break;
      if (k%2) ledSpect[(k%18)/2]->Fill(x);
      else pdsSpect[(k%18)/2]->Fill(x);
      k++;
   }*/
   char str[10000];
   while (!file.getline(str, sizeof(str)).eof()) {
      size_t count = 0;
      std::vector<size_t> ind;
      for(size_t i = 0;i<strlen(str);i++)
         if (str[i]=='.') {count++; ind.push_back(i-6);}
      if (count!=18&&count!=0) continue;
      if (!count) {
         ind.push_back(0);
         for(size_t i = 0;i<strlen(str);i++)
            if (str[i]==' ') {count++; ind.push_back(i);}
         if (count!=18) continue;
      }
      for(size_t i=0;i<18;i++) {
         double x;
         sscanf(&(str[ind[i]]),"%lf",&x);
         if (i%2) ledSpect[i/2]->Fill(x);
         else pdsSpect[i/2]->Fill(x);
      }
   }
   //
   std::vector<TFitResultPtr> pdsFitRes = PdsSpectFitting(pdsSpect);
   std::vector<TFitResultPtr> ledFitRes = LedSpectFitting(ledSpect,pdsFitRes);
   //
   // accproxy
   accproxy_cal cal;
   cal.t0_.resize(9,0);
   cal.dt_.resize(9,16);
   cal.thr_.resize(9,3*19);
   cal.dt0_.resize(9,0);
   cal.a_.resize(9,0);
   cal.da_.resize(9,0);
   //
   accsimpl_cal accSimplCal = AccSimplCalculation(pdsFitRes,ledFitRes,pdsSpect,ledSpect,cal);
   PrintAccSimpl(accSimplCal);
   //
   TFile *f = new TFile(fnamer.Data(),"RECREATE");
   for(size_t i=0;i<pdsSpect.size();i++)
      pdsSpect.at(i)->Write();
   for(size_t i=0;i<ledSpect.size();i++)
      ledSpect.at(i)->Write();
   for(size_t i=0;i<pdsFitRes.size();i++)
      pdsFitRes.at(i)->Write();
   for(size_t i=0;i<ledFitRes.size();i++)
      ledFitRes.at(i)->Write();
   f->Close();   
}

//-----------------------------------------------------------

class AccSimpleCal
{
 public:
   AccSimpleCal(TString rfile_){rfile = rfile_;init();};
   void readData();
   void readU();
   void init();
   void prepData();
   void printU();
   TString rfile;
   std::vector<Double_t> U;
   std::vector<TH1D*> pdsSpect;
   std::vector<TH1D*> ledSpect;
   std::vector<TFitResult*> pdsFitRes;
   std::vector<TFitResult*> ledFitRes;
   std::vector<Double_t> p0;
   std::vector<Double_t> dp0;
   std::vector<Double_t> mu;
   std::vector<Double_t> dmu;
   std::vector<Double_t> a1;
   std::vector<Double_t> da1;
};

void AccSimpleCal::prepData()
{
   p0.resize(0);
   dp0.resize(0);
   mu.resize(0);
   dmu.resize(0);
   a1.resize(0);
   da1.resize(0);
   for(size_t i=0;i<ledSpect.size();i++) {
      Double_t p0_ = ledFitRes.at(i)->Parameters().at(5);
      p0.push_back(p0_);
      //
      Double_t nevt = ledSpect.at(i)->GetEntries();
      Double_t dp0_ = ledFitRes.at(i)->Errors().at(5);
      dp0.push_back(dp0_);
      //
      Double_t mu_ = -log(p0_);
      mu.push_back(mu_);
      //
      Double_t dmu_ = dp0_/p0_;
      dmu.push_back(dmu_);
      //
      Double_t a0 = pdsSpect.at(i)->GetMean() - pdsFitRes.at(i)->Parameters().at(1);
      Double_t am = ledSpect.at(i)->GetMean() - ledFitRes.at(i)->Parameters().at(0);
      Double_t a1_ = (am-a0)/mu_;
      a1.push_back(a1_);
      //
      Double_t da0 = pdsSpect.at(i)->GetMeanError();
      Double_t dam = ledSpect.at(i)->GetMeanError();
      Double_t d1 = a1_*dmu_/mu_;
      Double_t d2 = dam/mu_;
      Double_t d3 = da0/mu_;
      Double_t da1_ = std::sqrt(d1*d1+d2*d2+d3*d3);
      da1.push_back(da1_);
   }
}
   
void AccSimpleCal::readData()
{
   TFile *f = new TFile(rfile);
   pdsSpect.resize(0);
   for(Int_t i=1;i<10;i++)
      pdsSpect.push_back((TH1D*)f->Get(Form("pds%d",i)));
   ledSpect.resize(0);
   for(Int_t i=1;i<10;i++)
      ledSpect.push_back((TH1D*)f->Get(Form("led%d",i)));
   pdsFitRes.resize(0);
   for(Int_t i=1;i<10;i++)
      pdsFitRes.push_back((TFitResult*)f->Get(Form("TFitResult-pds%d-f",i)));
   ledFitRes.resize(0);
   for(Int_t i=1;i<10;i++)
      ledFitRes.push_back((TFitResult*)f->Get(Form("TFitResult-led%d-f",i)));
   TText *t1 = (TText*)f->Get("Voltages");
   TText *t2 = (TText*)f->Get("Voltages+");
   TString s;
   if (strlen(t1->GetTitle()))  s = t1->GetTitle();
   if (strlen(t2->GetTitle()))  s = t2->GetTitle();
   TObjArray *x = s.Tokenize(" ");
   for (Int_t j = 0; j < x->GetEntries(); j++) 
      U.push_back(atof(((TObjString *)(x->At(j)))->String()));
//   f->Close();
}

void AccSimpleCal::readU()
{
   TString dfile = rfile(0,rfile.Last('.'));
   TString cmd(Form("grep -A 14 -e '%s' /work/users/konctbel/TESTS/flash/hvchange_history.txt |grep -e AG0 |awk '{print $9}'",dfile.Data()));
   TString str(gSystem->GetFromPipe(cmd.Data()));
   TObjArray *x = str.Tokenize("\n");
   U.resize(0);
   for (Int_t i = 0; i < x->GetEntries(); i++) U.push_back(atof(((TObjString *)(x->At(i)))->String()));
}

void AccSimpleCal::init()
{
//   readU();
   readData();
   prepData();
}

void AccSimpleCal::printU()
{
   std::cout << rfile.Data() << std::endl;
   for (Int_t i = 0; i < U.size(); i++)
      std::cout << "U[" << i << "] = " << U[i] << " V" << std::endl;
}

class AccSimpleCalU
{
 public:
   AccSimpleCalU(){};
   void add(AccSimpleCal *C);
   void add(AccSimpleCalU C);
   void plotP0(size_t ch);
   void plotMu(size_t ch);
   TGraphErrors* plotA1(size_t ch);
   std::vector<AccSimpleCal*> c;
};

void AccSimpleCalU::add(AccSimpleCal *C)
{
   c.push_back(C);
}

void AccSimpleCalU::add(AccSimpleCalU C)
{
   for(size_t i=0;i<C.c.size();i++)
      c.push_back(C.c[i]);
}

void AccSimpleCalU::plotP0(size_t ch)
{
   std::vector<Double_t> x;
   std::vector<Double_t> y;
   std::vector<Double_t> dx;
   std::vector<Double_t> dy;
   for(size_t i=0;i<c.size();i++) {
      x.push_back(std::fabs(c[i]->U[ch]));
      dx.push_back(0);
      y.push_back(c[i]->p0[ch]);
      dy.push_back(c[i]->dp0[ch]);
   }
   TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
   gr->SetMarkerSize(1);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kBlack);
   gr->Draw("AP");
}

void AccSimpleCalU::plotMu(size_t ch)
{
   std::vector<Double_t> x;
   std::vector<Double_t> y;
   std::vector<Double_t> dx;
   std::vector<Double_t> dy;
   for(size_t i=0;i<c.size();i++) {
      x.push_back(std::fabs(c[i]->U[ch]));
      dx.push_back(0);
      y.push_back(c[i]->mu[ch]);
      dy.push_back(c[i]->dmu[ch]);
   }
   TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
   gr->SetMarkerSize(1);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kBlack);
   gr->Draw("AP");
}
   
TGraphErrors* AccSimpleCalU::plotA1(size_t ch)
{
   std::vector<Double_t> x;
   std::vector<Double_t> y;
   std::vector<Double_t> dx;
   std::vector<Double_t> dy;
   for(size_t i=0;i<c.size();i++) {
      if (c[i]->a1[ch]<0) continue;
      if (c[i]->a1[ch]>1000) continue;
      if (c[i]->da1[ch]>100) continue;
      x.push_back(std::fabs(c[i]->U[ch]));
      dx.push_back(0);
      y.push_back(c[i]->a1[ch]);
      dy.push_back(c[i]->da1[ch]);
   }
   TGraphErrors *gr = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
   gr->SetMarkerSize(1);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kBlack);
   gr->Draw("AP");
   TF1* f = new TF1("f","([1]-[2])*pow(x/3000,[0])+[2]",2000,4000);
   f->SetParName(0,"a");
   f->SetParName(1,"b");
   f->SetParName(2,"c");
   f->SetParameters(20,200,0);
   f->SetParLimits(0,0,100);
   f->SetLineColor(kRed);
   f->SetLineWidth(0);
   gr->Fit(f);
   return gr;
}

void calAccSimpleU()
{
   std::vector<AccSimpleCalU> cal;
   cal.push_back(AccSimpleCalU());
   for(size_t i=0;i<8;i++)
      cal.back().add(new AccSimpleCal(Form("onel_060417_v%d.root",i)));
   cal.push_back(AccSimpleCalU());
   for(size_t i=0;i<11;i++)
      cal.back().add(new AccSimpleCal(Form("onel_200417_v%d.root",i)));
   //cal.push_back(AccSimpleCalU());
   //for(size_t i=1;i<14;i++)
   //   cal.back().add(new AccSimpleCal(Form("hvtest%d.root",i)));
   //cal.push_back(AccSimpleCalU());
   //for(size_t i=1;i<17;i++)
   //   cal.back().add(new AccSimpleCal(Form("onel_280417_%d.root",i)));
   //cal.push_back(AccSimpleCalU());
   //for(size_t i=1;i<17;i++)
   //   cal.back().add(new AccSimpleCal(Form("onel_280417_%d_v1.root",i)));
   cal.push_back(AccSimpleCalU());
   for(size_t i=1;i<17;i++)
      cal.back().add(new AccSimpleCal(Form("onel_040517_%d.root",i)));
   cal.push_back(AccSimpleCalU());
   for(size_t i=1;i<17;i++)
      cal.back().add(new AccSimpleCal(Form("onel_170517_%d.root",i)));
   cal.push_back(AccSimpleCalU());
   for(size_t i=1;i<17;i++)
      cal.back().add(new AccSimpleCal(Form("onel_070617_%d.root",i)));
   cal.push_back(AccSimpleCalU());
   for(size_t i=1;i<17;i++)
      cal.back().add(new AccSimpleCal(Form("onel_090617_%d.root",i)));
   cal.push_back(AccSimpleCalU());
   for(size_t i=1;i<17;i++)
      cal.back().add(new AccSimpleCal(Form("onel_140617_%d.root",i)));
   //
   AccSimpleCalU calS;
   for(size_t j=0;j<cal.size();j++)
      calS.add(cal[j]);
   std::vector<TGraphErrors*> gr;
   for(size_t i=0;i<9;i++) {
      gr.push_back(calS.plotA1(i));
      stop();
   }
   //13/03/2017 10:40 26954 1.13
   //Double_t u[9] = {2967, 2935, 2867, 3317, 2991, 2945, 2812, 2878, 3200};
   //27/03/2017 11:00 27139 1.13
   //Double_t u[9] = {2965, 2851, 2867, 3151, 3022, 2945, 2901, 2901, 3102};
   //06/04/2071 11:40 27366 1.13
   Double_t u[9] = {2992, 2856, 2898, 3178, 3046, 2945, 2906, 2926, 3111};
   //
   std::vector<AccSimpleCal> cx;
   //cx.push_back(AccSimpleCal("onel_170214_1.root"));
   cx.push_back(AccSimpleCal("onel_170524_1.root"));
   cx.push_back(AccSimpleCal("onel_170602_1.root"));
   cx.push_back(AccSimpleCal("onel_170602_2.root"));
   cx.push_back(AccSimpleCal("onel_170605_1.root"));
   //cx.push_back(AccSimpleCal("onel_170605_2.root"));
   cx.push_back(AccSimpleCal("onel_170607_1.root"));
   //cx.push_back(AccSimpleCal("onel_170609_1.root"));
   cx.push_back(AccSimpleCal("onel_170609_2.root"));
   cx.push_back(AccSimpleCal("onel_170609_3.root"));
   cx.push_back(AccSimpleCal("onel_170609_4.root"));
   cx.push_back(AccSimpleCal("onel_170609_5.root"));
   cx.push_back(AccSimpleCal("onel_170609_6.root"));
   cx.push_back(AccSimpleCal("onel_170614_1.root"));
   cx.push_back(AccSimpleCal("onel_170614_2.root"));
   //cx.push_back(AccSimpleCal("onel_170627_1.root"));
   //cx.push_back(AccSimpleCal("onel_170927_1.root"));
   //
   for(size_t ch=0;ch<9;ch++) {
      std::cout << gr[ch]->GetFunction("f")->Eval(u[ch]) << std::endl;
      Double_t a1 = gr[ch]->GetFunction("f")->Eval(u[ch]);
      size_t nx = cx.size();
      std::vector<Double_t> x(nx);
      std::vector<Double_t> y(nx,a1);
      std::vector<Double_t> dx(nx,0);
      std::vector<Double_t> dy(nx,0);
      for(size_t i=0;i<cx.size();i++) {
         x[i] = i+1;
         dx[i] = 0;
         std::cout << i << " " << cx[i].a1[ch] << " " << cx[i].da1[ch] << std::endl;
         if (cx[i].a1[ch]<0) continue;
         if (cx[i].a1[ch]>1000) continue;
         if (cx[i].da1[ch]>100) continue;
         y[i] = cx[i].a1[ch];
         dy[i] = cx[i].da1[ch];
      }
      TGraphErrors *grx = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
      grx->SetMarkerSize(1);
      grx->SetMarkerStyle(20);
      grx->SetMarkerColor(kBlack);
      grx->Draw("AP");
      TLine *l = new TLine();
      l->DrawLine(x.front(),a1,x.back(),a1);
      stop();
   }
   //
   for(size_t i=0;i<9;i++) {
      //cal.plotP0(i);
      //cal.plotMu(i);
      std::vector<TGraphErrors*> gr;
      for(size_t j=0;j<cal.size();j++) {
         gr.push_back(cal[j].plotA1(i));
         gr.back()->SetMarkerColor(j);
         //if (j) gr.back()->Draw("P");
         //else gr.back()->Draw("AP");
      }
      gr[0]->Draw("AP");
      gr[1]->Draw("P");
      gr[2]->Draw("P");
      gr[3]->Draw("P");
      gr[4]->Draw("P");
      gr[5]->Draw("P");
      gr[6]->Draw("P");
      stop();
   }
   
}

void WriteText(TString name, TString s)
{
   TText t(0,0,s.Data());
   t.SetName(name);
   t.Write();
}

void hvfind()
{
   //TString s = gSystem->GetFromPipe("ls -l /work/users/konctbel/TESTS/fullcal/Cal_2017*/*.res");
   TString s = gSystem->GetFromPipe("ls -l /work/users/konctbel/TESTS/fullcal/Simple_Calibration/Cal_171*/*.res");
   TObjArray *t = s.Tokenize("\n");
   for(int i=0;i<t->GetEntries();i++) {
      TString s = ((TObjString *)(t->At(i)))->String();
      TString month = gSystem->GetFromPipe(Form("echo %s | awk '{print $6}'",s.Data()));
      TString day = gSystem->GetFromPipe(Form("echo %s | awk '{print $7}'",s.Data()));
      TString year = gSystem->GetFromPipe(Form("echo %s | awk '{print $8}'",s.Data()));
      TString fnamed = gSystem->GetFromPipe(Form("echo %s | awk '{print $9}'",s.Data()));
      TString fname = fnamed(fnamed.Last('/')+1,fnamed.Sizeof()-fnamed.Last('/')-1);
      TString dir = fnamed(0,fnamed.Last('/'));
      TString ncal = dir(dir.Last('_')+1,dir.Sizeof()-dir.Last('_')-1);
      TString version;
      if (fname.Contains("hvtest"))
         version = fname(fname.Last('t')+1,fname.Last('.')-fname.Last('t')-1);
      else if (fname.Contains("_v"))
         version = fname(fname.Last('v')+1,fname.Last('.')-fname.Last('v')-1);
      else
         version = fname(fname.Last('_')+1,fname.Last('.')-fname.Last('_')-1);
      TString nopt = "unknown";
      if (fname.Contains("1.05")) nopt = "1.05";
      if (fname.Contains("1.13")) nopt = "1.13";
      std::cout << month << " " << day << " " << year << " " << fname << " " << version << " " << nopt << " " << dir << " " << ncal << std::endl;
      TString hvag = dir+"/hvag_new.dat";
      if (gSystem->AccessPathName(hvag)) hvag = dir+"/hvag.dat";
      std::cout << hvag << std::endl;
      TString u = gSystem->GetFromPipe(Form("sed '%sq;d' %s",version.Data(),hvag.Data()));
      std::cout << u << std::endl;
      //
      TString dfile = fname(0,fname.Last('.'))+".dat";
      TString cmd(Form("grep -A 14 -e '%s' /work/users/konctbel/TESTS/flash/hvchange_history.txt |grep -e AG0 |awk '{print $8}'",dfile.Data()));
      s = gSystem->GetFromPipe(cmd.Data());
      TObjArray *x = s.Tokenize("\n");
      std::vector<Double_t> U;
      for (Int_t j = 0; j < x->GetEntries(); j++) U.push_back(atof(((TObjString *)(x->At(j)))->String()));
      TString uf = "";
      for(size_t j=0; j<U.size(); j++) uf += Form("%d ",(int)U[j]);
      //
      std::vector<TH1D*> pdsSpect;
      std::vector<TH1D*> ledSpect;
      for(size_t i=1;i<10;i++) {
         TString hname = Form("pds%d",i);
         if (gROOT->FindObject(hname.Data())) gROOT->Delete(hname.Data());
         pdsSpect.push_back(new TH1D(hname.Data(),hname.Data(),4000,0,20000));
         hname = Form("led%d",i);
         if (gROOT->FindObject(hname.Data())) gROOT->Delete(hname.Data());
         ledSpect.push_back(new TH1D(hname.Data(),hname.Data(),4000,0,20000));
      }
      TString fnamer = fname(0,fname.Last('.')) + ".root";
      std::cout << fnamer << std::endl;
      ifstream file(fnamed.Data());
      char str[10000];
      while (!file.getline(str, sizeof(str)).eof()) {
         size_t count = 0;
         std::vector<size_t> ind;
         for(size_t i = 0;i<strlen(str);i++)
            if (str[i]=='.') {count++; ind.push_back(i-6);}
         if (count!=18&&count!=0) continue;
         if (!count) {
            ind.push_back(0);
            for(size_t i = 0;i<strlen(str);i++)
               if (str[i]==' ') {count++; ind.push_back(i);}
            if (count!=18) continue;
         }
         for(size_t i=0;i<18;i++) {
            double x;
            sscanf(&(str[ind[i]]),"%lf",&x);
            if (i%2) ledSpect[i/2]->Fill(x);
            else pdsSpect[i/2]->Fill(x);
         }
      }
      //
      std::vector<TFitResultPtr> pdsFitRes = PdsSpectFitting(pdsSpect);
      std::vector<TFitResultPtr> ledFitRes = LedSpectFitting(ledSpect,pdsFitRes);
      //
      // accproxy
      accproxy_cal cal;
      cal.t0_.resize(9,0);
      cal.dt_.resize(9,16);
      cal.thr_.resize(9,3*19);
      cal.dt0_.resize(9,0);
      cal.a_.resize(9,0);
      cal.da_.resize(9,0);
      //
      accsimpl_cal accSimplCal = AccSimplCalculation(pdsFitRes,ledFitRes,pdsSpect,ledSpect,cal);
      PrintAccSimpl(accSimplCal);
      //
      TFile *f = new TFile(fnamer.Data(),"RECREATE");
      for(size_t i=0;i<pdsSpect.size();i++)
         pdsSpect.at(i)->Write();
      for(size_t i=0;i<ledSpect.size();i++)
         ledSpect.at(i)->Write();
      for(size_t i=0;i<pdsFitRes.size();i++)
         pdsFitRes.at(i)->Write();
      for(size_t i=0;i<ledFitRes.size();i++)
         ledFitRes.at(i)->Write();
      WriteText("month",month);
      WriteText("day",day);
      WriteText("year",year);
      WriteText("fname",fname);
      WriteText("dir",dir);
      WriteText("ncal",ncal);
      WriteText("version",version);
      WriteText("nopt",nopt);
      WriteText("hvag",hvag);
      WriteText("Voltages",u);
      WriteText("Voltages+",uf);
      f->Close();
      //stop();
   }
}

//----------------------------------
void p0test(size_t ch = 1)
{
   //AccSimpleCal cal("onel_170517_16.root");
   AccSimpleCal cal("onel_170614_2.root");
   cal.prepData();
   size_t ic = ch - 1;
   std::cout << cal.p0[ic] << " " << cal.dp0[ic] << std::endl;
   cal.ledSpect[ic]->Draw();
   double l = cal.ledFitRes[ic]->Parameters().at(0) - 10*cal.ledFitRes[ic]->Parameters().at(1);
   double r = cal.ledFitRes[ic]->Parameters().at(0) + 3*cal.ledFitRes[ic]->Parameters().at(1);
   double q = cal.ledSpect[ic]->Integral(cal.ledSpect[ic]->GetXaxis()->FindBin(l),
                                         cal.ledSpect[ic]->GetXaxis()->FindBin(r))/
      cal.ledSpect[ic]->Integral();
   std::cout << q << " " << l << " " << r << std::endl;
   // accproxy
/*   accproxy_cal cal0;
   cal0.t0_.resize(9,0);
   cal0.dt_.resize(9,16);
   cal0.thr_.resize(9,3*19);
   cal0.dt0_.resize(9,0);
   cal0.a_.resize(9,0);
   cal0.da_.resize(9,0);
   //
   accsimpl_cal accSimplCal = AccSimplCalculation(cal.pdsFitRes,
                                                  cal.ledFitRes,
                                                  cal.pdsSpect,
                                                  cal.ledSpect,
                                                  cal0);
   PrintAccSimpl(accSimplCal);*/
}
