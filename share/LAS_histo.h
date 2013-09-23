#ifndef _LAS_HISTO_H_
#define _LAS_HISTO_H_

#include <string>
#include <vector>


// histo  for LAS_stabil
// 
// 04.08  VZ



using namespace std;

class TH1F;
class TH2F;
class TGraph;
class TGraphErrors;
class LAS_stabil;

class LAS_histo
{
 public:
LAS_histo();
LAS_histo(LAS_stabil* las);
~LAS_histo();
void setlas(LAS_stabil* alas) {thelas=alas;};
void book();
void fill();
int write(string fnama);
void reset();
//
int fillprofiles(string fname);

string resultsav(string fname, double av);
int fillresults(LAS_stabil* alas);
int fillstab(string fname);
void setdebug(int d) {debug=d;};

 int bh;
 int bstab;
 int bstabt;
 int bstabtt;
 int debug;
 int nb;
 int btime;
 int brun;
 int nav;
 LAS_stabil* thelas;
// control histos
// initial profiles
 TH2F* h2_tobprof[8];
 TH2F* h2_tibprof[8];
 TH2F* h2_tecpprof[8];
 TH2F* h2_tecmprof[8];
 // means, rms and residuals
 TH2F* h2_tobrms[8];
 TH2F* h2_tibrms[8];
 TH2F* h2_tecprms[8];
 TH2F* h2_tecmrms[8];
 TH2F* h2_tobmean[8];
 TH2F* h2_tibmean[8];
 TH2F* h2_tecpmean[8];
 TH2F* h2_tecmmean[8];
 TH2F* h2_tobres[8];
 TH2F* h2_tibres[8];
 TH2F* h2_tecpres[8];
 TH2F* h2_tecmres[8];
 TH2F* h2_tobnorm[8];
 TH2F* h2_tibnorm[8];
 TH2F* h2_tecpnorm[8];
 TH2F* h2_tecmnorm[8];
 TH1F* h1_tobstat[8];
 TH1F* h1_tibstat[8];
 TH1F* h1_tecpstat[8];
 TH1F* h1_tecmstat[8];

 TH2F* h2_tecpmeanint[2][8];
 TH2F* h2_tecmmeanint[2][8];



 TH2F* h2_tobstat;
 TH2F* h2_tibstat;
 TH2F* h2_tecpstat;
 TH2F* h2_tecmstat;

 // tecinternal
 

 // AT params
TH1F* h1_ata[8];
TH1F* h1_atb[8];
TH1F* h1_atch2[8];
TH1F* h1_atn[8];
// results
TH1F* h1_tobpar[5];
TH1F* h1_tobepar[5];
TH1F* h1_tobch2;
TH1F* h1_tibpar[5];
TH1F* h1_tibepar[5];
TH1F* h1_tibch2;
TH1F* h1_tecppar[3];
TH1F* h1_tecpepar[3];
TH1F* h1_tecpch2;
TH1F* h1_tecmpar[3];
TH1F* h1_tecmepar[3];
TH1F* h1_tecmch2;
TH1F* h1_ch2;
TH1F* h1_sign;

TH1F* h1_tecpintpar[3];
TH1F* h1_tecpeintpar[3];
TH1F* h1_tecmintpar[3];
TH1F* h1_tecmeintpar[3];
 
 TGraphErrors* g_tbstab[5];
 TGraphErrors* g_tecpstab[3];
 TGraphErrors* g_tecmstab[3];
 TGraphErrors* g_ch2;
 TGraphErrors* g_sign;
 TGraphErrors* g_ch2tb;
 TGraphErrors* g_ch2tecp;
 TGraphErrors* g_ch2tecm;
 
 TGraphErrors* g_tecpintstab[3];
 TGraphErrors* g_tecmintstab[3];

 
 TGraphErrors* gt_tbstab[5];
 TGraphErrors* gt_tecpstab[3];
 TGraphErrors* gt_tecmstab[3];
 TGraphErrors* gt_ch2;
 TGraphErrors* gt_sign;
 TGraphErrors* gt_ch2tb;
 TGraphErrors* gt_ch2tecp;
 TGraphErrors* gt_ch2tecm;

 TGraphErrors* gt_tecpintstab[3];
 TGraphErrors* gt_tecmintstab[3];
 
 TGraphErrors* gtt_tbstab[5];
 TGraphErrors* gtt_tecpstab[3];
 TGraphErrors* gtt_tecmstab[3];
 TGraphErrors* gtt_ch2;
 TGraphErrors* gtt_sign;

 TGraphErrors* gtt_tecpintstab[3];
 TGraphErrors* gtt_tecmintstab[3];


};


#endif
