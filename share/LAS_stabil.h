#ifndef _LAS_STABIL_H_
#define _LAS_STABIL_H_

#include <string>
#include <vector>
#include "TMatrixT.h"

// calculate residuals for the main algn parameters 
// TIB, TECP, TECM in respect to alignment tubes attached to TOB
// 25.07  VZ



using namespace std;

class TH1F;
class TH2F;
class TGraph;
class  TRandom;
class LAS_histo;

class LAS_stabil 
{
 public:
LAS_stabil();
~LAS_stabil();
void init_geom();
int readflist(const string& sflist, const string& spref);
int alignAll();
int alignTIB();
int alignTOB();
int alignTEC(int sgn);
int alignATtob();
int alignATtib();

int alignblockATtob();

int alignTECint(int sgn);
int alignTECintRDC();

void setdebug(int d) {debug=d;};
void setrefrunn(int rnn) {nrefrun=rnn;};
void setrefrun(int rn) {refrunn=rn;};
int  writeres(const string& fout);
int  sumch2();

//tools
void printh(TMatrixT<double>& h);
void clearh(TMatrixT<double>* h);
void printa(double* a, int n);
void cleara(double* a, int n);
void cleara(int* a, int n);
void clearv(vector<double*>* vv);
void clearv(vector<int*>* vv);
int  printres();

 void enablecorat(int d) {bcor_at=d;}
 
int getnruns() {return nruns;}
void setnruns(int n) {if(n<nruns) nruns=n;}
vector<double>*  gettime() {return &vtime;   };
vector<double>*  getrun() {return &vruns;   };


vector<double*>*  getpartb() {return &vpar_tb;   };
vector<double*>*  getpartecp() {return &vpar_tecp; };
vector<double*>*  getpartecm() {return &vpar_tecm; };
vector<double*>*  getparertb() {return &vparer_tb;   };
vector<double*>*  getparertecp() {return &vparer_tecp; };
vector<double*>*  getparertecm() {return &vparer_tecm; };

vector<double*>*  getpartecpint() {return &vpar_tecpint; };
vector<double*>*  getpartecmint() {return &vpar_tecmint; };
vector<double*>*  getparertecpint() {return &vparer_tecpint; };
vector<double*>*  getparertecmint() {return &vparer_tecmint; };


vector<double>*  getsign() {return &vsign; };
vector<double>*  getallch2() {return &vch2; };
vector<double*>*  getatch2() {return &vch2_at; };
vector<double>*  gettbch2() {return &vch2_tb; };
vector<double>*  gettecpch2() {return &vch2_tecp; };
vector<double>*  gettecmch2() {return &vch2_tecm; };
vector<int*>*  getatndf() {return &vndf_at; };
vector<int>*  gettbndf() {return &vndf_tb; };
vector<int>*  gettecpndf() {return &vndf_tecp; };
vector<int>*  gettecmndf() {return &vndf_tecm; };

// MC 
void mcenablemean(int op) {bmc_mean=op; bmc+=op; };
void mcenableat(int op)    {bmc_at=op;  bmc+=op;};
void mcsmear(double* d, double s);
void mcoffset(double* d1, double d2);
void updateresmc(int s);
void updaterestb(int s);
void updaterestecp(int s);
void updaterestecpint(int s);
void updaterestecmint(int s);
void updaterestecm(int s);
void calculatermsmc();
void calculatesign();
// histos
void histenable(int i);
void sethisto(LAS_histo* ah) {thehisto=ah;};
LAS_histo* gethisto() {return thehisto;}

void setbat(int d) {battob=d; }
void settecin(int d) {btecint=d; }

//TEC internal
void tecintsummary();

// offset all par to means
 void offsetpar();

 private:
int debug;
vector<string>  sfiles;
int nruns;
int nrefrun;
int refrunn;
int battob;

vector<double> vtime;
vector<double>  vruns;
vector<double*> varo;
vector<double*> vbro;
vector<double*> vaaro;
vector<int*> vnp;



vector<double*>  vpar_tb;
vector<double*>  vpar_tecp;
vector<double*>  vpar_tecm;

vector<double*>  vpar_tecpint;
vector<double*>  vpar_tecmint;
vector<double> vch2_tecpint;
vector<double> vch2_tecmint;
vector<int> vndf_tecpint;
vector<int> vndf_tecmint;



 vector<double*> vch2_at;
 vector<double> vch2_tb;
 vector<double> vch2_tecp;
 vector<double> vch2_tecm;
 vector<int*>  vndf_at;
 vector<int> vndf_tb;
 vector<int> vndf_tecp;
 vector<int> vndf_tecm;
 vector<double>  vch2;
 vector<double>  vsign;

// running  mean and rms
 int utb;
 int utecp;
 int utecm;
 int utecpint;
 int utecmint;



vector<double*>  vparer_tb;
vector<double*>  vparer_tecp;
vector<double*>  vparer_tecm;


vector<double*>  vparm_tb;
vector<double*>  vparm_tecp;
vector<double*>  vparm_tecm;

vector<double*>  vparm_tecpint;
vector<double*>  vparm_tecmint;
vector<double*>  vparer_tecpint;
vector<double*>  vparer_tecmint;

// TEC int
vector<double*>  vpar_tecpinta;
vector<double*>  vpar_tecminta;
vector<double*>  vparer_tecpinta;
vector<double*>  vparer_tecminta;



// status
 int stat_at;
 int stat_tib;
 int stat_tecp;
 int stat_tecm;
vector<int*> stat;
// geom
double izp_tob[6],izp_tib[6],izp_tecp[9],izp_tecm[9];
double phib[8];
double pch_tob,pch_tib,rad_tob,rad_tib, pch_tec,rad_tec;
int im[8],ii[8][6];
 double pch_tecint[2], rad_tecint[2];
// smearing for MC
 int bmc;
 int bmc_mean;
 int bmc_at;
 int bcor_at;
 int btecint;

 TRandom * thernd;
 //common sigma profiles (in pitches)
double maxrms;
double mcsig_mean;
double mcsigtob_mean;
double mcsigtib_mean;
double mcsigtec_mean;

// individual
 double** mcsigmtob;// for each semsor
 double** mcsigmtib;
 double** mcsigmtecp;
 double**  mcsigmtecm;
// sigma AT params
double mcsig_ata;
double mcsig_atb;

double mcsig_atec;
double mcsig_btec;
// histo
int bhist;
LAS_histo* thehisto;
// selection 
 double Cdiffmax;
 double Cmin;
 double Cmax;
 double Crms;
 int Cnpm;
 double Cch2at;
 int   Cblock;
};


#endif
