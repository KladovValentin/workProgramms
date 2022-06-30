#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TChain.h"
#include "TString.h"
#include "TList.h"
#include "TCollection.h"
#include "TSystemDirectory.h"
#include "TSystem.h"
#include <string>
#include <iostream>
#include <stdio.h>
#include <fstream> 
#include <vector> 
//#include "/work/users/kladov/snd2k/R006-003/maindir/2017/minmaxrun.cpp"
using namespace std;


TString str2 = "c";
int apred = 0;
char simb = (char)34;  //"//
int counter = 1;

int getlist() {
	ofstream fout;
	fout.open("/work/users/kladov/snd2k/R006-003/fwk/simflist2.fwi");
	//fout.open("/work/users/kladov/snd2k/R007-001/runrecot.sh");

	/*vector<TString> ans;
	TString basedir("/online/gridback/MC/R006-001/ee/output");  // /online/gridback/MC/R006-001/ee/output // /online/simulation/MC/R006-004/ee/output
	TString input("/sweet/home/kladov/");
	TSystemDirectory dir(input, basedir);
	TList *files = dir.GetListOfFiles();
	cout << " basedir " << basedir << " input " << input << endl;
	if (files) {
		TIter next(files);
		TSystemFile *file;
		TString fname;

		while ((file = (TSystemFile*)next())) {

			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(".psy.gz")) {
				ans.push_back(fname);
			}
		}

	}*/
	
	//2c4n_ps_nrc, 2etg_wrc, 2kc2p0_wrc, 2kcp0_knkst, 2kcp0_kpkst, 2kcp0_kstk, 4pi_wrc, et2p0g_wrc, etapg_nrc, klknpp_wrc, 
	//klknppp0_wrc, klkppn_wrc, klkppnp0_wrc, kskl2p0_wrc, kskleta_wrc, ksknpp_ks22p0_wrc, ksknppp0_wrc, kskppn_ks22p0_wrc, kskppnp0_wrc, ompi4_wrc.
	//important: histcommon_col: условие на кинфиты -, поменять в simreco-neu_wmc.fwk имя папки и создать ее output/ntuples/2etg, fwk/simflist.fwi сюда вставить результат скрипта, выбрать строки в runrecokonct
	TString str = "a";
	TString str1 = "b";
	
	vector<TString> ans;
	//TString basedir("/online/gridback/MC/R006-001/ee/output");
	//TString basedir("/online/simulation/MC/R006-004/ee/output/");
	//TString basedir("/online/simulation/MC/R006-004/kketa/output/");
	//TString basedir("/online/simulation/MC/R006-004/kkpi/output");
	TString basedir("/online/simulation/MC/R006-004/2k2pi/output");
	TString simbs = simb;
	//TString files = gSystem->GetFromPipe("find " + basedir + " -name " + simbs + "*t30_nemcc_odch*" + simbs + " | sort");
	//TString files = gSystem->GetFromPipe("find " + basedir + " -name " + simbs + "*kskleta_wrc_nemcc*" + simbs + " | sort");
	//TString files = gSystem->GetFromPipe("find " + basedir + " -name " + simbs + "*kskppn_ks22p0_wrc*" + simbs + " | sort");
	//TString files = gSystem->GetFromPipe("find " + basedir + " -name " + simbs + "*kskcpc_ks22p0_wrc_nemcc_odch-*-50000.psy.gz" + simbs + " | sort");
	//TString files = gSystem->GetFromPipe("find " + basedir + " -name " + simbs + "*kskcpc_ks22p0_wrc_nemcc_odch*" + simbs + " | sort");
	//TString files = gSystem->GetFromPipe("find " + basedir + " -name " + simbs + "*kskppn_ks22p0_wrc_nemcc_odch*" + simbs + " | sort");
	TString files = gSystem->GetFromPipe("find " + basedir + " -name " + simbs + "2kc2p0_wrc_nemcc_odch*" + simbs + " | sort");
	//TString files = gSystem->GetFromPipe("find /online/simulation/MC/ -name " + simbs + "*4pi_wrc*" + simbs + " ! -name " + simbs + "*sim*" + simbs + " | sort");
	//TString files = gSystem->GetFromPipe("find /online/simulation/MC/ /online/gridback/MC/ -name " + simbs + "*9426*" + simbs + " ! -name " + simbs + "*MHAD*" + simbs + " -print"); //-825-
	//TString files = gSystem->GetFromPipe("ls /online/simulation/MC/R007-001/kkpi/output/kskp*");
	TString elem;
	Ssiz_t from0 = 0;
	while (files.Tokenize(elem, from0, "\n")) {
		if (elem.EndsWith(".psy.gz"))
			ans.push_back(elem);
	}


	fout << "flist = [" << (char)92 << endl;	
	TString process_name = "a";
	unsigned int vector_size = ans.size();
	for (int i = 0; i < vector_size; i++) {
		str = ans[i].Copy();
		str1 = ans[i].Copy();
		str2 = ans[i].Copy();
		TString tok;
		TString tok1;
		TString tok2;
		Ssiz_t from = 0;
		Ssiz_t from1 = 0;
		Ssiz_t from2 = 0;
		vector<TString> tokens;
		vector<TString> tokens1;
		vector<TString> tokens2;
		//while (str.Tokenize(tok, from, "[_-]")) {
		while (str.Tokenize(tok, from, "[/.-]")) {
			//cout << tok << endl;
			tokens.push_back(tok);
		}
		while (str1.Tokenize(tok1, from1, "[/.]")) {
			//cout << tok1 << endl;
			tokens1.push_back(tok1);
		}
		const char *a0 = (const char*)tokens[9];
		int a = atoi(a0);
		const char* b0 = (const char*)tokens[8];
		int b = atoi(b0);

		if (a == apred) {
			counter = counter + 1;
		}
		else {
			counter = 1;
		}

		char count[11];
		sprintf(count, "%d", counter);
		TString c = count;
		TString indyr = "/" + tokens[0] + "/" + tokens[1] + "/" + tokens[2] + "/" + tokens[3] + "-" + tokens[4] + "/" + tokens[5] + "/" + tokens[6] + simb;
		TString infname = simb + tokens1[6] + "." + tokens[12] + "." + tokens[13] + simb;
		TString outfname = simb + tokens[7]+"_"+ tokens[8]+"-"+tokens[9]+"-"+c + simb;
		TString nen = tokens[11];
		//cout << b << endl;

		//if ((a >= 7842) && (a <= 10988) && counter <2) {
		if ((a >= 11285) && (a <= 13219) && counter <2) {
		//if ((a >= 7842) && (a <= 13219) && counter <2) {
		//if ((a >= 37920) && (a <= 42826) && counter < 2 ) {
		//if ((a >= 25525) && (a <= 29461) && counter < 2 ) {
		//if ((a >= 27251) && (a <= 29461) && counter < 2 ) {
		//if ((a >= 9426-5) && (a <= 9426+5)) {
			//fout << "( " << simb << indyr << " , " << infname << " , " << a << " , " << outfname << " , 1 , " << nen << " , 9380 , 9450 )" << "," << endl;
			//if(b>750 && b<800 )
				fout << "( " << simb << indyr << " , " << infname << " , " << a << " , " << outfname << " , 1 , " << nen << " , " << a-25 << " , " << a+25 << " )" << "," << endl;
			cout << b << "-" << a << endl;
			process_name = tokens[7];
		}
		//cout << endl;
			
		apred = a;
	}
	fout << "]" << endl;
	fout.close();
	
	/*vector<TString> ans;
	TString basedir = "/work/users/kladov/snd2k/R007-001/2019/*.hbook";
	//TString basedir = "/online/simulation/MC/R006-004/ee/output/*.mod.gz";
	TString files = gSystem->GetFromPipe("ls " + basedir);
	TString elem;
	Ssiz_t from = 0;
	while (files.Tokenize(elem, from, "\n")) {
		ans.push_back(elem);
	}
	int vector_size = ans.size();
	TString str = "a";
	TString str1 = "b";
	fout << "source /work/snd2000/root/setup2k.sh i386-SL5-debug" << endl << "source/etc/sysconfig/gridengine_nsu" << endl << endl;
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
			//cout << tok << endl;
			tokens1.push_back(tok1);
		}

		//cout << tokens1[6] << endl;

		const char* a0 = (const char*)tokens[9];
		int a = atoi(a0);
		if (a > 27250 && a > 37919 && a < 42827) {
			fout << ".mainrelease/Offline/submit.sh -q clusters,1440 FAKERUN=" << a << " MODFILENAME=" << tokens1[6] << "     MODFILEDIR=/online/simulation/MC/R006-004/ee/output/ HBOOKDIR=./2019/  SimRecApp fwk/simreco_col_point.fwk" << endl;
			cout << "h2root " << "2019/" << tokens1[6] << ".hbook " << "2019/" << tokens1[6] << ".root" << (char)59 << " ";
		}
		
	}
	cout << endl;
	fout.close();*/

	return 1;
}

void splitlines(){
	ifstream fin;
	ofstream fout;
	fin.open("fwk/recflistT.fwi");
	fout.open("fwk/recflist_my.fwi",std::ios_base::app);
	string name;
	while(getline(fin,name)){
		TString a = name;
		vector<TString> ans;
		TString elem;
		Ssiz_t from0 = 0;
		while (a.Tokenize(elem, from0, ",")) {
			ans.push_back(elem);
		}

		TString tok;
		Ssiz_t from = 0;
		vector<TString> tokens;
		if(ans.size()>=5){
			const char *a1 = (const char*)ans[2];
			int b1 = atoi(a1);
			const char *a2 = (const char*)ans[3];
			int b2 = atoi(a2);
			//cout << b1 << "	" << b2 << "	";

			while (ans[4].Tokenize(tok, from, "[_]")) {
				tokens.push_back(tok);
				//cout << tok << "	";
			}
			if(tokens.size() >=3 ){
				int c1 = b1+(b2-b1)/3;
				int c2 = b1+2*(b2-b1)/3;
				fout << ans[0] << "," << ans[1] << ", " << b1 << " , " << c1 << " , " << simb << "000" << b1 << "_000" << c1 << "_" << tokens[2] << "," << endl;
				fout << ans[0] << "," << ans[1] << ", " << c1+1 << " , " << c2 << " , " << simb << "000" << c1+1 << "_000" << c2 << "_" << tokens[2] << "," << endl;
				fout << ans[0] << "," << ans[1] << ", " << c2+1 << " , " << b2 << " , " << simb << "000" << c2+1 << "_000" << b2 << "_" << tokens[2] << "," << endl;
			}
		}
		/*for(size_t i = 0; i < ans.size(); i++){
			cout << ans[i] << "	; ";
		}
		cout << endl;*/
	}
	fin.close();
	fout.close();
}

void getCalibRuns(string calName, bool read) {
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
	unsigned int vector_size = ans.size();
	for (int i = 1; i < vector_size-1; i++) {
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
		cout << a << endl;
	}
}

//pscp D:\programs\listoffiles.cpp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R006-003/maindir/listoffiles.cpp
//pscp kladov@sndxt1.inp.nsk.su:/work/users/kladov/snd2k/R007-001/fwk/simflist.fwi D:\simflist.txt


