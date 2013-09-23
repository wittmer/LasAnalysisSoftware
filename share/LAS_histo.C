#include <iostream>
#include <ctime>
#include <string>
#include <functional>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMatrixT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"

#define __AVECROOT__
#include "Avec.h"

#include "LAS_histo.h"
#include "LAS_stabil.h"
#include "LASGlobalData.h"
#include "LASGlobalDataLoop.h"
#include "LAS_vectorfloat_tools.h"
#include "LAS_basic_tools.h"
#include "LAS_globaldata_tools.h"
#include "LAS_RDC_tools.h"

using namespace std;

//
LAS_histo::LAS_histo() {
  debug=1;
  thelas=0;
  btime=1;
  brun=1;
  nav=1;
  bh=0;
  bstab=0;
  bstabt=0;
  bstabtt=0;
}
//
LAS_histo::LAS_histo(LAS_stabil* alas) {
  debug=0;
  thelas=alas;
 btime=1;
 brun=1;
  bstab=0;
  bstabt=0;
  bstabtt=0;
  book();
}

//
LAS_histo::~LAS_histo() {
 
}

//
void LAS_histo::book() {
  if(debug) cout<<" LAS_histo::book>>>>"<<endl;
  nb=8;
  int ndmax1=6;
  int ndmax2=4;
  int ndmax3=9;
  double xm1=10;
  double xm2=500;
  double xrms=100;
  double xr1=-50;
  double xr2=50;

 double xpr1=-100;
  double xpr2=100;
  int nb=100;
  int nbm=8;
  // AT related
 char si[2];
  for(int i=0;i<nbm;i++) {
sprintf(si,"%d",i);
string ss1(si);
string ss="b"+ss1;
// cout<<" book "<<ss<<endl;
 { const char* name=("h2_"+ss+"_tobrms").c_str();    h2_tobrms[i]= new TH2F(name,name,ndmax1,0,ndmax1,nb,0,xrms); } 
 { const char* name=("h2_"+ss+"_tibrms").c_str();    h2_tibrms[i]= new TH2F(name,name,ndmax1,0,ndmax1,nb,0,xrms); } 
 { const char* name=("h2_"+ss+"_tecprms").c_str();   h2_tecprms[i]= new TH2F(name,name,ndmax2,0,ndmax2,nb,0,xrms); } 
 { const char* name=("h2_"+ss+"_tecmrms").c_str();   h2_tecmrms[i]= new TH2F(name,name,ndmax2,0,ndmax2,nb,0,xrms); } 

 { const char* name=("h2_"+ss+"_tobmean").c_str();  
   //  cout<<" book "<<name <<endl;
 h2_tobmean[i]= new TH2F(name,name,ndmax1,0,ndmax1,nb,xm1,xm2); } 

 { const char* name=("h2_"+ss+"_tibmean").c_str();   h2_tibmean[i]= new TH2F(name,name,ndmax1,0,ndmax1,nb,xm1,xm2); } 
 { const char* name=("h2_"+ss+"_tecpmean").c_str();  h2_tecpmean[i]= new TH2F(name,name,ndmax2,0,ndmax2,nb,xm1,xm2); } 
 { const char* name=("h2_"+ss+"_tecmmean").c_str();  h2_tecmmean[i]= new TH2F(name,name,ndmax2,0,ndmax2,nb,xm1,xm2); } 

 { const char* name=("h2_"+ss+"_tobres").c_str();    h2_tobres[i]= new TH2F(name,name,ndmax1,0,ndmax1,nb,xr1,xr2); } 
 { const char* name=("h2_"+ss+"_tibres").c_str();    h2_tibres[i]= new TH2F(name,name,ndmax1,0,ndmax1,nb,xr1,xr2); } 
 { const char* name=("h2_"+ss+"_tecpres").c_str();   h2_tecpres[i]= new TH2F(name,name,ndmax2,0,ndmax2,nb,xr1,xr2); } 
 { const char* name=("h2_"+ss+"_tecmres").c_str();   h2_tecmres[i]= new TH2F(name,name,ndmax2,0,ndmax2,nb,xr1,xr2); } 

 { const char* name= ("h2_"+ss+"_tobnorm").c_str();   h2_tobnorm[i]= new TH2F(name,name,ndmax1,0,ndmax1,nb,0,1000); } 
 { const char* name= ("h2_"+ss+"_tibnorm").c_str();   h2_tibnorm[i]= new TH2F(name,name,ndmax1,0,ndmax1,nb,0,1000); } 
 { const char* name= ("h2_"+ss+"_tecpnorm").c_str();  h2_tecpnorm[i]= new TH2F(name,name,ndmax2,0,ndmax2,nb,0,1000); } 
 { const char* name= ("h2_"+ss+"_tecmnorm").c_str();  h2_tecmnorm[i]= new TH2F(name,name,ndmax2,0,ndmax2,nb,0,1000); } 

 { const char* name=("h1_"+ss+"_tobstat").c_str();    h1_tobstat[i]= new TH1F(name,name,ndmax1,0,ndmax1); } 
 { const char* name=("h1_"+ss+"_tibstat").c_str();    h1_tibstat[i]= new TH1F(name,name,ndmax1,0,ndmax1); } 
 { const char* name=("h1_"+ss+"_tecpstat").c_str();   h1_tecpstat[i]= new TH1F(name,name,ndmax2,0,ndmax2); } 
 { const char* name=("h1_"+ss+"_tecmstat").c_str();   h1_tecmstat[i]= new TH1F(name,name,ndmax2,0,ndmax2); } 


 { const char* name=("h1_"+ss+"_ata").c_str();       h1_ata[i]= new TH1F(name,name,nb,-0.0001,0.0001 ); } 
 { const char* name=("h1_"+ss+"_atb").c_str();       h1_atb[i]= new TH1F(name,name, nb,-0.1,0.1); } 
 { const char* name= ("h1_"+ss+"_atch2").c_str();    h1_atch2[i]= new TH1F(name,name, nb,0,10); } 
 { const char* name= ("h1_"+ss+"_atn").c_str();    h1_atn[i]= new TH1F(name,name, 10,0,10); } 
   
 { const char* name=("h2_"+ss+"_tecpmeanintr4").c_str();   h2_tecpmeanint[0][i]= new TH2F(name,name,ndmax2,0,ndmax3,nb,xm1,xm2); } 
 { const char* name=("h2_"+ss+"_tecmmeanintr4").c_str();   h2_tecmmeanint[0][i]= new TH2F(name,name,ndmax2,0,ndmax3,nb,xm1,xm2); } 

 { const char* name=("h2_"+ss+"_tecpmeanintr6").c_str();   h2_tecpmeanint[1][i]= new TH2F(name,name,ndmax2,0,ndmax3,nb,xm1,xm2); } 
 { const char* name=("h2_"+ss+"_tecmmeanintr6").c_str();   h2_tecmmeanint[1][i]= new TH2F(name,name,ndmax2,0,ndmax3,nb,xm1,xm2); } 

}//i

  // tib param related
  int npartib=5;
 for(int i=0;i<npartib;i++) {
sprintf(si,"%d",i);
string ss1(si);
string ss="_p"+ss1;
 { const char* name=("h1_"+ss+"_tib").c_str();      h1_tibpar[i]=  new TH1F(name,name,nb,xpr1, xpr2); } 
 { const char* name=("h1_"+ss+"_tibe").c_str();     h1_tibepar[i]=  new TH1F(name,name,nb,0, xr2); } 
{ const char* name=("h1_"+ss+"_tob").c_str();      h1_tobpar[i]=  new TH1F(name,name,nb,xpr1, xpr2); } 
 { const char* name=("h1_"+ss+"_tobe").c_str();     h1_tobepar[i]=  new TH1F(name,name,nb,0, xr2); } 
 }//i

  // tecp param related
  int npartec=3;
 for(int i=0;i<npartec;i++) {
sprintf(si,"%d",i);
string ss1(si);
string ss="_p"+ss1;
 { const char* name=("h1_"+ss+"_tecp").c_str();     h1_tecppar[i]=  new TH1F(name,name,nb,xpr1, xpr2); } 
 { const char* name=("h1_"+ss+"_tecpe").c_str();    h1_tecpepar[i]=  new TH1F(name,name,nb,0, xr2); } 
 { const char* name=("h1_"+ss+"_tecm").c_str();     h1_tecmpar[i]=  new TH1F(name,name,nb,xpr1, xpr2); } 
 { const char* name=("h1_"+ss+"_tecme").c_str();    h1_tecmepar[i]=  new TH1F(name,name,nb,0, xr2); } 

 { const char* name=("h1_"+ss+"_tecpint").c_str();     h1_tecpintpar[i]=  new TH1F(name,name,nb,xpr1, xpr2); } 
 { const char* name=("h1_"+ss+"_tecpeint").c_str();    h1_tecpeintpar[i]=  new TH1F(name,name,nb,0, xr2); } 
 { const char* name=("h1_"+ss+"_tecmint").c_str();     h1_tecmintpar[i]=  new TH1F(name,name,nb,xpr1, xpr2); } 
 { const char* name=("h1_"+ss+"_tecmeint").c_str();    h1_tecmeintpar[i]=  new TH1F(name,name,nb,0, xr2); } 
 
 }//i

{ const char* name="h1_tobch2";      h1_tobch2=  new TH1F(name,name,nb,-1, xpr2); } 
{ const char* name="h1_tibch2";      h1_tibch2=  new TH1F(name,name,nb,-1, xpr2); } 
{ const char* name="h1_tecpch2";     h1_tecpch2=  new TH1F(name,name,nb,-1, xpr2); } 
{ const char* name="h1_tecmch2";     h1_tecmch2=  new TH1F(name,name,nb,-1, xpr2); } 
{ const char* name="h1_ch2";         h1_ch2=  new TH1F(name,name,nb,0, xpr2); } 
{ const char* name="h1_sign";        h1_sign=  new TH1F(name,name,nb,0, 10); } 

{ const char* name= "h2_tobstat";   h2_tobstat= new TH2F(name,name,ndmax1,0,ndmax1,nbm,0,nbm); } 
{ const char* name= "h2_tibstat";   h2_tibstat= new TH2F(name,name,ndmax1,0,ndmax1,nbm,0,nbm); }
{ const char* name= "h2_tecpstat";  h2_tecpstat= new TH2F(name,name,ndmax2,0,ndmax2,nbm,0,nbm); }
{ const char* name= "h2_tecmstat";  h2_tecmstat= new TH2F(name,name,ndmax2,0,ndmax2,nbm,0,nbm); }

 
if(debug) cout<<" >>>>LAS_histo::book done"<<endl;
 bh=1;
 return;
}
//
void LAS_histo::fill() {

}
//
int LAS_histo::fillprofiles(string fname) {

  return 0;
}
// read out txt file from LAS_stabil and analyze
// format with MC errors
string LAS_histo::resultsav(string fname, double avv) {

 if(debug) cout<<" LAS histo resultsav>>> "<<fname<<endl;
ifstream infile (fname.c_str());
 if(!infile)  {
    cout<<" *********cannot open  file "<<fname<<endl;;
    return "1";
  }

 nav=int(avv);

vector<double> vruns;
vector<double> vtime;
vector<double> vch2;
vector<double> vsign;
vector<double*> vpar_tib;
vector<double*> vparer_tib;
vector<double*> vpar_tecp;
vector<double*> vparer_tecp;
vector<double*> vpar_tecm;
vector<double*> vparer_tecm;
vruns.clear();
vtime.clear();
vch2.clear();
vsign.clear();
vpar_tib.clear(); // assure there are no arrays inside!
vparer_tib.clear(); 
vpar_tecp.clear();
vparer_tecp.clear();
vpar_tecm.clear();
vparer_tecm.clear();

vector<double*> vpar_tecpint;
vector<double*> vparer_tecpint;
vector<double*> vpar_tecmint;
vector<double*> vparer_tecmint;
vpar_tecpint.clear();
vparer_tecpint.clear();
vpar_tecmint.clear();
vparer_tecmint.clear();

// list of bad and good runs
 vector<double> vbadruns;
 vector<double> vgoodruns;
 vbadruns.clear();
 vgoodruns.clear();
//
 int nptib=5;
 int nptec=3;
 double atime;
 double arun;
 // create for first time and recreate for each new averaged run
 double* atib=new double[nptib];
 double* atibe=new double[nptib];
 double* atecp=new double[nptec];
 double* atecpe=new double[nptec];
 double* atecm=new double[nptec];
 double* atecme=new double[nptec];
 double ach2;
 double asign;

 double* atecpint=new double[nptec];
 double* atecpeint=new double[nptec];
 double* atecmint=new double[nptec];
 double* atecmeint=new double[nptec];

 // this can stay
 double* atib0=new double[nptib];
 double* atibe0=new double[nptib];
 double* atecp0=new double[nptec];
 double* atecpe0=new double[nptec];
 double* atecm0=new double[nptec];
 double* atecme0=new double[nptec];

 double* atecpint0=new double[nptec];
 double* atecpeint0=new double[nptec];
 double* atecmint0=new double[nptec];
 double* atecmeint0=new double[nptec];

 double ach20=0;
 double asign0=0;
 double atime0=0;
 double arun0=0;

  for(int i=0;i<nptib;i++) {
      atib[i]=0.;
      atibe[i]=0.;
    }
  for(int i=0;i<nptec;i++) {
     atecp[i]=0.;
      atecpe[i]=0.;
      atecm[i]=0.;
      atecme[i]=0.;

      atecpint[i]=0.;
      atecpeint[i]=0.;
      atecmint[i]=0.;
      atecmeint[i]=0.;

    }

int nentry=0;
int nr0=nav; // number of runs to average, can be 1
int nr=nr0;
if(debug>10) cout<<" start averaging "<<nav<<endl;
int ne=0;
 
if(nav==0) return "2";

while(infile
>>atime0
>>arun0
>>atib0[0]
>>atibe0[0]
>>atib0[1]
>>atibe0[1]
>>atib0[2]
>>atibe0[2]
>>atib0[3]
>>atibe0[3]
>>atib0[4]
>>atibe0[4]
>>atecp0[0]
>>atecpe0[0]
>>atecp0[1]
>>atecpe0[1]
>>atecp0[2]
>>atecpe0[2]
>>atecm0[0]
>>atecme0[0]
>>atecm0[1]
>>atecme0[1]
>>atecm0[2]
>>atecme0[2]
>>ach20
>>asign0
>>atecpint0[0]
>>atecpeint0[0]
>>atecpint0[1]
>>atecpeint0[1]
>>atecpint0[2]
>>atecpeint0[2]
>>atecmint0[0]
>>atecmeint0[0]
>>atecmint0[1]
>>atecmeint0[1]
>>atecmint0[2]
>>atecmeint0[2]


) {



  //select 
  double sum=0;
  double sume=0;
 for(int i=0;i<nptib;i++) {
   sum+=atib0[i];
  sume+=atibe0[i];
 }

 for(int i=0;i<nptec;i++) {
   sum+=atecp0[i];
  sume+=atecpe0[i];
 sum+=atecm0[i];
  sume+=atecme0[i];
 }
 // selections
 bool bbadr=(fabs(sum)>500.||fabs(sume)>500.)&&(ach20>10||ach20<0.01)||(asign0>10);

 bool bgoodr=!bbadr&&ach20<5;

 if(bbadr) vbadruns.push_back(arun0);
 if(bgoodr) vgoodruns.push_back(arun0);

 if((fabs(sum)>500.)) continue;
 if((fabs(sume)>500.)) continue;
if(ach20>100.||ach20<0.01) continue;

  if(debug>10) { 
    cout<<nentry<<" run: "<<arun0<<" nr= "<<nr
	<<" "<<atime0
	<<" "<<arun0
	<<" "<<atib0[0]
	<<" "<<atibe0[0]
	<<" "<<atib0[1]
	<<" "<<atibe0[1]
	<<" "<<atib0[2]
	<<" "<<atibe0[2]
	<<" "<<atib0[3]
	<<" "<<atibe0[3]
	<<" "<<atib0[4]
	<<" "<<atibe0[4]
	<<" "<<atecp0[0]
	<<" "<<atecpe0[0]
	<<" "<<atecp0[1]
	<<" "<<atecpe0[1]
	<<" "<<atecp0[2]
	<<" "<<atecpe0[2]
	<<" "<<atecm0[0]
	<<" "<<atecme0[0]
	<<" "<<atecm0[1]
	<<" "<<atecme0[1]
	<<" "<<atecm0[2]
	<<" "<<atecme0[2]
	<<" "<<ach20
	<<" "<<asign0
<<" "<<atecpint0[0]
<<" "<<atecpeint0[0]
<<" "<<atecpint0[1]
<<" "<<atecpeint0[1]
<<" "<<atecpint0[2]
<<" "<<atecpeint0[2]
<<" "<<atecmint0[0]
<<" "<<atecmeint0[0]
<<" "<<atecmint0[1]
<<" "<<atecmeint0[1]
<<" "<<atecmint0[2]
<<" "<<atecmeint0[2]

 <<endl;
  }

if(nr==0) {
  if(debug>10) cout<<nentry<<" new meas "<<ne<<"  nr="<<nr<<endl;
    // new slice
  atime=atime/nr0;
  arun=arun/nr0;
  ach2=ach2/nr0;
  asign=asign/nr0;
    for(int i=0;i<nptib;i++) {
      // average
      atib[i]=atib[i]/nr0;
      atibe[i]=sqrt(atibe[i]/nr0);
      if(debug>10) cout<<ne<<" average  tib "<<i<<" "<<atib[i]<<" "<<atibe[i]<<endl;
    }
  for(int i=0;i<nptec;i++) {
      atecp[i]=atecp[i]/nr0;
      atecpe[i]=sqrt(atecpe[i])/nr0;
      atecm[i]=atecm[i]/nr0;
      atecme[i]=sqrt(atecme[i])/nr0;

      atecpint[i]=atecpint[i]/nr0;
      atecpeint[i]=sqrt(atecpeint[i])/nr0;
      atecmint[i]=atecmint[i]/nr0;
      atecmeint[i]=sqrt(atecmeint[i])/nr0;

    }
 
  // fill vectors
  vruns.push_back(arun);
  vtime.push_back(atime);
  vch2.push_back(ach2);
  vsign.push_back(asign);
  // here are the arrays of pointers
  vpar_tib.push_back(atib);
  vparer_tib.push_back(atibe);
  vpar_tecp.push_back(atecp);
  vparer_tecp.push_back(atecpe);
  vpar_tecm.push_back(atecm);
  vparer_tecm.push_back(atecme);

  vpar_tecpint.push_back(atecpint);
  vparer_tecpint.push_back(atecpeint);
  vpar_tecmint.push_back(atecmint);
  vparer_tecmint.push_back(atecmeint);

  if(debug>10) cout<<ne<<" fill vectors run="<< vruns.at(ne)<<" time="<< vtime.at(ne)<<" vpar_tib0="<<vpar_tib.at(ne)[0]<<endl;

  atime=0;
  arun=0;
  ach2=0;
  asign=0;
  nr=nr0;// restore nr
  // create and clear for the next run
atib=new double[nptib];
atibe=new double[nptib];
atecp=new double[nptec];
atecpe=new double[nptec];
atecm=new double[nptec];
atecme=new double[nptec];

atecpint=new double[nptec];
atecpeint=new double[nptec];
atecmint=new double[nptec];
atecmeint=new double[nptec];


for(int i=0;i<nptec;i++) {
      atecp[i]=0.;
      atecpe[i]=0.;
      atecm[i]=0.;
      atecme[i]=0.;
   
      atecpint[i]=0.;
      atecpeint[i]=0.;
      atecmint[i]=0.;
      atecmeint[i]=0.;


    }
  for(int i=0;i<nptib;i++) {
      atib[i]=0.;
      atibe[i]=0.;
    }

if(debug>10) cout<<ne<<" ++fill vectors run="<< vruns.at(ne)<<" time="<< vtime.at(ne)<<" vpar_tib0="<<vpar_tib.at(ne)[0]<<endl;

  ne++;
 }//  nr==0

  if(nr>0) {
    atime+=atime0;
    arun+=arun0;
    ach2+=ach20;
    asign+=asign0;
 for(int i=0;i<nptib;i++) {
    atib[i]+=atib0[i];
    atibe[i]+=atibe0[i]*atibe0[i];
 }
 for(int i=0;i<nptec;i++) {
    atecp[i]+=atecp0[i];
    atecpe[i]+=atecpe0[i]*atecpe0[i];
    atecm[i]+=atecm0[i];
    atecme[i]+=atecme0[i]*atecme0[i];
    //
   atecpint[i]+=atecpint0[i];
    atecpeint[i]+=atecpeint0[i]*atecpeint0[i];
    atecmint[i]+=atecmint0[i];
    atecmeint[i]+=atecmeint0[i]*atecmeint0[i];


 }
 if(debug>10) cout<<" ne="<<ne<<" nr="<<nr<<"  atib[0]="<<atib[0]<<" atime="<<atime<<"  arun="<<arun<<endl;
   nr--;
  }// nr>0
 
  nentry++;
 
 }// while runs

if(debug>10) cout<<" 0=fill vectors run="<< vruns.at(0)<<" time="<< vtime.at(0)<<" vpar_tib0="<<vpar_tib.at(0)[0]<<endl;

// list of bad runs
 string foutbad="listbad"+fname;
fstream fileb(foutbad.c_str(), std::ios::out);
 for(unsigned int i=0; i<vbadruns.size(); i++) {
   fileb<< vbadruns[i]<<endl;
 }
 fileb.close();
// list of good runs
 string foutgood="listgood"+fname;
fstream fileg(foutgood.c_str(), std::ios::out);
 for(unsigned int i=0; i<vbadruns.size(); i++) {
   fileg<< vgoodruns[i]<<endl;
 }
 fileg.close();


if(debug>10) cout<<" --finish averaging "<<vruns.size()<<" "<<nentry<<endl;
 int nruns= vruns.size();


 // fill into file the averaged

 string fout="av"+fname;
fstream file(fout.c_str(), std::ios::out);
 if(!file) {  cout<<"****cant create fout "<<endl;return "3";}
  file << setprecision(6);

 for(int irun=0; irun<nruns; irun++) {
   // store only nonzero results
double sum=vpar_tib.at(irun)[0]
    +vpar_tib.at(irun)[1]
    +vpar_tib.at(irun)[2]
    +vpar_tib.at(irun)[3]
    +vpar_tib.at(irun)[4]
    +vpar_tecp.at(irun)[0]
    +vpar_tecp.at(irun)[1]
    +vpar_tecp.at(irun)[2]
    +vpar_tecm.at(irun)[0]
    +vpar_tecm.at(irun)[1]
  +vpar_tecm.at(irun)[2];

  if(sum==0) continue;

    file<< vtime.at(irun)
    <<" "<<vruns.at(irun)
    <<" "<<vpar_tib.at(irun)[0]
    <<" "<<vparer_tib.at(irun)[0]
    <<" "<<vpar_tib.at(irun)[1]
    <<" "<<vparer_tib.at(irun)[1]
    <<" "<<vpar_tib.at(irun)[2]
    <<" "<<vparer_tib.at(irun)[2]
    <<" "<<vpar_tib.at(irun)[3]
    <<" "<<vparer_tib.at(irun)[3]
    <<" "<<vpar_tib.at(irun)[4]
    <<" "<<vparer_tib.at(irun)[4]
  
    <<" "<<vpar_tecp.at(irun)[0]
    <<" "<<vparer_tecp.at(irun)[0]
    <<" "<<vpar_tecp.at(irun)[1]
    <<" "<<vparer_tecp.at(irun)[1]
    <<" "<<vpar_tecp.at(irun)[2]
    <<" "<<vparer_tecp.at(irun)[2]
   
    <<" "<<vpar_tecm.at(irun)[0]
    <<" "<<vparer_tecm.at(irun)[0]
    <<" "<<vpar_tecm.at(irun)[1]
    <<" "<<vparer_tecm.at(irun)[1]
    <<" "<<vpar_tecm.at(irun)[2]
    <<" "<<vparer_tecm.at(irun)[2]

    <<" "<<vch2.at(irun)
    <<" "<<vsign.at(irun)


 <<" "<<vpar_tecpint.at(irun)[0]
    <<" "<<vparer_tecpint.at(irun)[0]
    <<" "<<vpar_tecpint.at(irun)[1]
    <<" "<<vparer_tecpint.at(irun)[1]
    <<" "<<vpar_tecpint.at(irun)[2]
    <<" "<<vparer_tecpint.at(irun)[2]
   
    <<" "<<vpar_tecmint.at(irun)[0]
    <<" "<<vparer_tecmint.at(irun)[0]
    <<" "<<vpar_tecmint.at(irun)[1]
    <<" "<<vparer_tecmint.at(irun)[1]
    <<" "<<vpar_tecmint.at(irun)[2]
    <<" "<<vparer_tecmint.at(irun)[2]
<<endl;


 }// iruns

 file.close();

 if(debug) cout<<" >>>>LAS_histo::resutlsav "<<fout<<endl;
  return fout;
}


// read saved txt file and fill stability graphs only  
int LAS_histo::fillstab(string fname) {

 if(debug) cout<<" LAS histo fillstab>>> "<<fname<<endl;
ifstream infile (fname.c_str());
 if(!infile)  {
    cout<<" *********cannot open  file "<<fname<<endl;;
    return 1;
  }


vector<double> vruns;
vector<double> vtime;
vector<double> vch2;
vector<double> vsign;
vector<double*> vpar_tib;
vector<double*> vparer_tib;
vector<double*> vpar_tecp;
vector<double*> vparer_tecp;
vector<double*> vpar_tecm;
vector<double*> vparer_tecm;

vector<double*> vpar_tecpint;
vector<double*> vparer_tecpint;
vector<double*> vpar_tecmint;
vector<double*> vparer_tecmint;

vruns.clear();
vtime.clear();
vch2.clear();
vsign.clear();
vpar_tib.clear(); // assure there are no arrays inside!
vparer_tib.clear(); 
vpar_tecp.clear();
vparer_tecp.clear();
vpar_tecm.clear();
vparer_tecm.clear();

vpar_tecpint.clear();
vparer_tecpint.clear();
vpar_tecmint.clear();
vparer_tecmint.clear();
//
 int nptib=5;
 int nptec=3;
 double atime;
 double arun;
 // create for first time and recreate for each new averaged run
 double* atib=new double[nptib];
 double* atibe=new double[nptib];
 double* atecp=new double[nptec];
 double* atecpe=new double[nptec];
 double* atecm=new double[nptec];
 double* atecme=new double[nptec];

 double* atecpint=new double[nptec];
 double* atecpeint=new double[nptec];
 double* atecmint=new double[nptec];
 double* atecmeint=new double[nptec];

 double ach2;
 double asign;

  for(int i=0;i<nptib;i++) {
      atib[i]=0.;
      atibe[i]=0.;
    }
  for(int i=0;i<nptec;i++) {
     atecp[i]=0.;
      atecpe[i]=0.;
      atecm[i]=0.;
      atecme[i]=0.;
    }

int nentry=0;

if(nav==0) return 2;

while(infile
>>atime
>>arun
>>atib[0]
>>atibe[0]
>>atib[1]
>>atibe[1]
>>atib[2]
>>atibe[2]
>>atib[3]
>>atibe[3]
>>atib[4]
>>atibe[4]
>>atecp[0]
>>atecpe[0]
>>atecp[1]
>>atecpe[1]
>>atecp[2]
>>atecpe[2]
>>atecm[0]
>>atecme[0]
>>atecm[1]
>>atecme[1]
>>atecm[2]
>>atecme[2]
>>ach2
>>asign
>>atecpint[0]
>>atecpeint[0]
>>atecpint[1]
>>atecpeint[1]
>>atecpint[2]
>>atecpeint[2]
>>atecmint[0]
>>atecmeint[0]
>>atecmint[1]
>>atecmeint[1]
>>atecmint[2]
>>atecmeint[2]
) {


  if(debug>1) { 
    cout<<nentry<<" run: "<<arun
	<<" "<<atime
	<<" "<<arun
	<<" "<<atib[0]
	<<" "<<atib[0]
	<<" "<<atib[1]
	<<" "<<atibe[1]
	<<" "<<atib[2]
	<<" "<<atibe[2]
	<<" "<<atib[3]
	<<" "<<atibe[3]
	<<" "<<atib[4]
	<<" "<<atibe[4]
	<<" "<<atecp[0]
	<<" "<<atecpe[0]
	<<" "<<atecp[1]
	<<" "<<atecpe[1]
	<<" "<<atecp[2]
	<<" "<<atecpe[2]
	<<" "<<atecm[0]
	<<" "<<atecme[0]
	<<" "<<atecm[1]
	<<" "<<atecme[1]
	<<" "<<atecm[2]
	<<" "<<atecme[2]
	<<" "<<ach2
 <<endl;
    }

  // make some selection
  // remove zeros


  // remove large errors





  // fill vectors
  vruns.push_back(arun);
  vtime.push_back(atime);
  vch2.push_back(ach2);
  vsign.push_back(asign);
  // here are the arrays of pointers
  vpar_tib.push_back(atib);
  vparer_tib.push_back(atibe);
  vpar_tecp.push_back(atecp);
  vparer_tecp.push_back(atecpe);
  vpar_tecm.push_back(atecm);
  vparer_tecm.push_back(atecme);

  vpar_tecpint.push_back(atecpint);
  vparer_tecpint.push_back(atecpeint);
  vpar_tecmint.push_back(atecmint);
  vparer_tecmint.push_back(atecmeint);


  // create and clear for the next run
atib=new double[nptib];
atibe=new double[nptib];
atecp=new double[nptec];
atecpe=new double[nptec];
atecm=new double[nptec];
atecme=new double[nptec];


atecpint=new double[nptec];
atecpeint=new double[nptec];
atecmint=new double[nptec];
atecmeint=new double[nptec];


for(int i=0;i<nptec;i++) {
      atecp[i]=0.;
      atecpe[i]=0.;
      atecm[i]=0.;
      atecme[i]=0.;

      atecpint[i]=0.;
      atecpeint[i]=0.;
      atecmint[i]=0.;
      atecmeint[i]=0.;
    }
  for(int i=0;i<nptib;i++) {
      atib[i]=0.;
      atibe[i]=0.;
    }

  if(debug>10) cout<<" nentry"<<nentry<<" atime="<<atime<<"  arun="<<arun<<"  atib[0]="<<atib[0]<<endl; 
  nentry++;
 
 }// while runs



if(debug>2) cout<<" --finish reading "<<vruns.size()<<" "<<nentry<<endl;
 int nruns= vruns.size();


  // for graphs
double dtime[nruns];
 double drun[nruns];
  
 double sign[nruns];
 double ch2[nruns];
 double ch2tib[nruns];
 double ch2tecp[nruns];
 double ch2tecm[nruns];

 double dx_tib[nruns]; 
 double dy_tib[nruns]; 
 double rx_tib[nruns]; 
 double ry_tib[nruns]; 
 double rz_tib[nruns];  
 double dxer_tib[nruns]; 
 double dyer_tib[nruns]; 
 double rxer_tib[nruns]; 
 double ryer_tib[nruns]; 
 double rzer_tib[nruns];  

 double dx_tecp[nruns]; 
 double dy_tecp[nruns]; 
 double rz_tecp[nruns]; 
 double dxer_tecp[nruns]; 
 double dyer_tecp[nruns]; 
 double rzer_tecp[nruns]; 

 double dx_tecm[nruns]; 
 double dy_tecm[nruns]; 
 double rz_tecm[nruns]; 
 double dxer_tecm[nruns]; 
 double dyer_tecm[nruns]; 
 double rzer_tecm[nruns]; 


 double dx_tecpint[nruns]; 
 double dy_tecpint[nruns]; 
 double rz_tecpint[nruns]; 
 double dxer_tecpint[nruns]; 
 double dyer_tecpint[nruns]; 
 double rzer_tecpint[nruns]; 

 double dx_tecmint[nruns]; 
 double dy_tecmint[nruns]; 
 double rz_tecmint[nruns]; 
 double dxer_tecmint[nruns]; 
 double dyer_tecmint[nruns]; 
 double rzer_tecmint[nruns]; 



 double dum[nruns];

// fill grphs
 for(int irun=0; irun<nruns; irun++) {
 dum[irun]=0.;

dtime[irun]=vtime[irun];
drun[irun]=vruns[irun];
//tib
 
dx_tib[irun]=vpar_tib.at(irun)[0];       if(fabs(dx_tib[irun])>100) dx_tib[irun]=100; 
if(debug>10) cout<<irun<<" fill  graph dxtib="<<dx_tib[irun]<<" dtime="<<dtime[irun]<<endl;
dy_tib[irun]=vpar_tib.at(irun)[1];        if(fabs(dy_tib[irun])>100) dy_tib[irun]=100;
rx_tib[irun]=vpar_tib.at(irun)[2];        if(fabs(rx_tib[irun])>100) rx_tib[irun]=100;
ry_tib[irun]=vpar_tib.at(irun)[3];        if(fabs(ry_tib[irun])>100) ry_tib[irun]=100;
rz_tib[irun]=vpar_tib.at(irun)[4];       if(fabs(rz_tib[irun])>100) rz_tib[irun]=100;
 if(vparer_tib.size()>0) {
 dxer_tib[irun]=vparer_tib.at(irun)[0];  if(dxer_tib[irun]>100) dxer_tib[irun]=100;
dyer_tib[irun]=vparer_tib.at(irun)[1];   if(dyer_tib[irun]>100) dyer_tib[irun]=100;
rxer_tib[irun]=vparer_tib.at(irun)[2];   if(rxer_tib[irun]>100) rxer_tib[irun]=100;
ryer_tib[irun]=vparer_tib.at(irun)[3];   if(ryer_tib[irun]>100) ryer_tib[irun]=100;
rzer_tib[irun]=vparer_tib.at(irun)[4];   if(rzer_tib[irun]>100) rzer_tib[irun]=100;
 }
//tecp 
 dx_tecp[irun]=vpar_tecp.at(irun)[0];     if(fabs(dx_tecp[irun])>100) dx_tecp[irun]=100;
 dy_tecp[irun]=vpar_tecp.at(irun)[1];      if(fabs(dy_tecp[irun])>100) dy_tecp[irun]=100;
 rz_tecp[irun]=vpar_tecp.at(irun)[2];      if(fabs(rz_tecp[irun])>100) rz_tecp[irun]=100;
 if(vparer_tecp.size()>0) {
dxer_tecp[irun]=vparer_tecp.at(irun)[0];  if(dxer_tecp[irun]>100) dxer_tecp[irun]=100;
dyer_tecp[irun]=vparer_tecp.at(irun)[1];  if(dyer_tecp[irun]>100) dyer_tecp[irun]=100;
rzer_tecp[irun]=vparer_tecp.at(irun)[2];  if(rzer_tecp[irun]>100) rzer_tecp[irun]=100;
 }
//tecm
dx_tecm[irun]=vpar_tecm.at(irun)[0];      if(fabs(dx_tecm[irun])>100) dx_tecm[irun]=100;
dy_tecm[irun]=vpar_tecm.at(irun)[1];      if(fabs(dy_tecm[irun])>100) dy_tecm[irun]=100;
rz_tecm[irun]=vpar_tecm.at(irun)[2];     if(fabs(rz_tecm[irun])>100) rz_tecm[irun]=100;
 if(vparer_tecm.size()>0) {
dxer_tecm[irun]=vparer_tecm.at(irun)[0]; if(dxer_tecm[irun]>100) dxer_tecm[irun]=100;
dyer_tecm[irun]=vparer_tecm.at(irun)[1]; if(dyer_tecm[irun]>100) dyer_tecm[irun]=100;
rzer_tecm[irun]=vparer_tecm.at(irun)[2]; if(rzer_tecm[irun]>100) rzer_tecm[irun]=100;
 }
 // tecpint
 dx_tecpint[irun]=vpar_tecpint.at(irun)[0];     if(fabs(dx_tecpint[irun])>100) dx_tecpint[irun]=100;
 dy_tecpint[irun]=vpar_tecpint.at(irun)[1];      if(fabs(dy_tecpint[irun])>100) dy_tecpint[irun]=100;
 rz_tecpint[irun]=vpar_tecpint.at(irun)[2];      if(fabs(rz_tecpint[irun])>100) rz_tecpint[irun]=100;
 if(vparer_tecpint.size()>0) {
dxer_tecpint[irun]=vparer_tecpint.at(irun)[0];  if(dxer_tecpint[irun]>100) dxer_tecpint[irun]=100;
dyer_tecpint[irun]=vparer_tecpint.at(irun)[1];  if(dyer_tecpint[irun]>100) dyer_tecpint[irun]=100;
rzer_tecpint[irun]=vparer_tecpint.at(irun)[2];  if(rzer_tecpint[irun]>100) rzer_tecpint[irun]=100;
 }
//tecmint
dx_tecmint[irun]=vpar_tecmint.at(irun)[0];      if(fabs(dx_tecmint[irun])>100) dx_tecmint[irun]=100;
dy_tecmint[irun]=vpar_tecmint.at(irun)[1];      if(fabs(dy_tecmint[irun])>100) dy_tecmint[irun]=100;
rz_tecmint[irun]=vpar_tecmint.at(irun)[2];     if(fabs(rz_tecmint[irun])>100) rz_tecmint[irun]=100;
 if(vparer_tecmint.size()>0) {
dxer_tecmint[irun]=vparer_tecmint.at(irun)[0]; if(dxer_tecmint[irun]>100) dxer_tecmint[irun]=100;
dyer_tecmint[irun]=vparer_tecmint.at(irun)[1]; if(dyer_tecmint[irun]>100) dyer_tecmint[irun]=100;
rzer_tecmint[irun]=vparer_tecmint.at(irun)[2]; if(rzer_tecmint[irun]>100) rzer_tecmint[irun]=100;
 }




ch2[irun]=vch2[irun];       if(ch2[irun]>1000) ch2[irun]=1000;
sign[irun]=vsign[irun];       if(sign[irun]>100) sign[irun]=100;

 }// iruns



 //make offset in time for easier presenetations
 double dtime1[nruns];
 double dtime2[nruns];
double dtime3[nruns];
double dtime4[nruns];
double dtime5[nruns];
 double off=1000.;
 for(int i=0; i<nruns; i++) {
   dtime1[i]=dtime[i]+off;
   dtime2[i]=dtime[i]+2.*off;
   dtime3[i]=dtime[i]+3.*off;
   dtime4[i]=dtime[i]+4.*off;
   dtime5[i]=dtime[i]+5.*off;
 }



 if(btime) {
   bstabtt=1;
gtt_tbstab[0]=new TGraphErrors(nruns,dtime,dx_tib,dum,dxer_tib);
gtt_tbstab[1]=new TGraphErrors(nruns,dtime1,dy_tib,dum,dyer_tib);
gtt_tbstab[2]=new TGraphErrors(nruns,dtime2,rx_tib,dum,rxer_tib);
gtt_tbstab[3]=new TGraphErrors(nruns,dtime3,ry_tib,dum,ryer_tib);
gtt_tbstab[4]=new TGraphErrors(nruns,dtime4,rz_tib,dum,rzer_tib);


gtt_tecpstab[0]=new TGraphErrors(nruns,dtime,dx_tecp,dum,dxer_tecp);
gtt_tecpstab[1]=new TGraphErrors(nruns,dtime1,dy_tecp,dum,dyer_tecp);
gtt_tecpstab[2]=new TGraphErrors(nruns,dtime2,rz_tecp,dum,rzer_tecp);

gtt_tecmstab[0]=new TGraphErrors(nruns,dtime,dx_tecm,dum,dxer_tecm);
gtt_tecmstab[1]=new TGraphErrors(nruns,dtime1,dy_tecm,dum,dyer_tecm);
gtt_tecmstab[2]=new TGraphErrors(nruns,dtime2,rz_tecm,dum,rzer_tecm);


gtt_tecpintstab[0]=new TGraphErrors(nruns,dtime,dx_tecpint,dum,dxer_tecpint);
gtt_tecpintstab[1]=new TGraphErrors(nruns,dtime1,dy_tecpint,dum,dyer_tecpint);
gtt_tecpintstab[2]=new TGraphErrors(nruns,dtime2,rz_tecpint,dum,rzer_tecpint);

gtt_tecmintstab[0]=new TGraphErrors(nruns,dtime,dx_tecmint,dum,dxer_tecmint);
gtt_tecmintstab[1]=new TGraphErrors(nruns,dtime1,dy_tecmint,dum,dyer_tecmint);
gtt_tecmintstab[2]=new TGraphErrors(nruns,dtime2,rz_tecmint,dum,rzer_tecmint);



gtt_sign=new TGraphErrors(nruns,dtime,sign,dum,dum);

gtt_ch2=new TGraphErrors(nruns,dtime,ch2,dum,dum);
 } 


 delete []  atib;
 delete []  atibe;
 delete []  atecp;
 delete []  atecpe;
 delete [] atecm;
 delete [] atecme;


 if(debug) cout<<" LAS_histo:avresults >>>"<< endl;

  return 0;
}






// use LAS_stabil directly to access files
int LAS_histo::fillresults(LAS_stabil* alas) {
 if(debug) cout<<"  LAS_histo::fill resutls>>>> "<<endl;

  if(alas==0) return 1;

vector<double*>  vpar_tib=   *(alas->getpartb());
vector<double*>  vpar_tecp=   *(alas->getpartecp());
vector<double*>  vpar_tecm=   *(alas->getpartecm());
vector<double*>  vparer_tib=   *(alas->getparertb());
vector<double*>  vparer_tecp=  *(alas->getparertecp()); 
vector<double*>  vparer_tecm =  *(alas->getparertecm()); 
vector<double>  vch2  =    *(alas->getallch2());
vector<double>  vsign  =    *(alas->getsign());
vector<double>  vch2tib  =    *(alas->gettbch2());
vector<double>  vch2tecp  =    *(alas->gettecpch2());
vector<double>  vch2tecm  =    *(alas->gettecmch2());
vector<double>  vtime =    *(alas->gettime());
vector<double>  vrun =    *(alas->getrun());


vector<double*>  vpar_tecpint=   *(alas->getpartecpint());
vector<double*>  vpar_tecmint=   *(alas->getpartecmint());
vector<double*>  vparer_tecpint=   *(alas->getparertecpint());
vector<double*>  vparer_tecmint=   *(alas->getparertecmint());

int nptib=5;
int nptec=3;

 int nruns=vrun.size();

 double dtime[nruns];
 double drun[nruns];
  
 double ch2[nruns];
 double sign[nruns];
 double ch2tib[nruns];
 double ch2tecp[nruns];
 double ch2tecm[nruns];

 double dx_tib[nruns]; 
 double dy_tib[nruns]; 
 double rx_tib[nruns]; 
 double ry_tib[nruns]; 
 double rz_tib[nruns];  
 double dxer_tib[nruns]; 
 double dyer_tib[nruns]; 
 double rxer_tib[nruns]; 
 double ryer_tib[nruns]; 
 double rzer_tib[nruns];  

 double dx_tecp[nruns]; 
 double dy_tecp[nruns]; 
 double rz_tecp[nruns]; 
 double dxer_tecp[nruns]; 
 double dyer_tecp[nruns]; 
 double rzer_tecp[nruns]; 

 double dx_tecm[nruns]; 
 double dy_tecm[nruns]; 
 double rz_tecm[nruns]; 
 double dxer_tecm[nruns]; 
 double dyer_tecm[nruns]; 
 double rzer_tecm[nruns]; 

double dx_tecpint[nruns]; 
 double dy_tecpint[nruns]; 
 double rz_tecpint[nruns]; 
 double dxer_tecpint[nruns]; 
 double dyer_tecpint[nruns]; 
 double rzer_tecpint[nruns]; 

 double dx_tecmint[nruns]; 
 double dy_tecmint[nruns]; 
 double rz_tecmint[nruns]; 
 double dxer_tecmint[nruns]; 
 double dyer_tecmint[nruns]; 
 double rzer_tecmint[nruns]; 


 double dum[nruns];

for(int irun=0;irun<nruns; irun++) {

 
for(int k=0;k<nptib; k++) {
  h1_tibpar[k]->Fill(vpar_tib.at(irun)[k]);
  if(vparer_tib.size()>0)
  h1_tibepar[k]->Fill(vparer_tib.at(irun)[k]);
}

for(int k=0;k<nptec; k++) {
  h1_tecppar[k]->Fill(vpar_tecp.at(irun)[k]);
 if(vparer_tecp.size()>0)
  h1_tecpepar[k]->Fill(vparer_tecp.at(irun)[k]);
  h1_tecmpar[k]->Fill(vpar_tecm.at(irun)[k]);
 if(vparer_tecm.size()>0)
  h1_tecmepar[k]->Fill(vparer_tecm.at(irun)[k]);


  h1_tecpintpar[k]->Fill(vpar_tecpint.at(irun)[k]);
 if(vparer_tecpint.size()>0)
  h1_tecpeintpar[k]->Fill(vparer_tecpint.at(irun)[k]);
  h1_tecmintpar[k]->Fill(vpar_tecmint.at(irun)[k]);
 if(vparer_tecmint.size()>0)
  h1_tecmeintpar[k]->Fill(vparer_tecmint.at(irun)[k]);


}

h1_ch2->Fill(vch2.at(irun));
h1_sign->Fill(vsign.at(irun));
// grphs
dum[irun]=0.;

dtime[irun]=vtime[irun];
drun[irun]=vrun[irun];
//tib
dx_tib[irun]=vpar_tib.at(irun)[0];        if(fabs(dx_tib[irun])>100) dx_tib[irun]=100; 
dy_tib[irun]=vpar_tib.at(irun)[1];        if(fabs(dy_tib[irun])>100) dy_tib[irun]=100;
rx_tib[irun]=vpar_tib.at(irun)[2];        if(fabs(rx_tib[irun])>100) rx_tib[irun]=100;
ry_tib[irun]=vpar_tib.at(irun)[3];        if(fabs(ry_tib[irun])>100) ry_tib[irun]=100;
rz_tib[irun]=vpar_tib.at(irun)[4];       if(fabs(rz_tib[irun])>100) rz_tib[irun]=100;
 if(vparer_tib.size()>0) {
 dxer_tib[irun]=vparer_tib.at(irun)[0];  if(dxer_tib[irun]>100) dxer_tib[irun]=100;
dyer_tib[irun]=vparer_tib.at(irun)[1];   if(dyer_tib[irun]>100) dyer_tib[irun]=100;
rxer_tib[irun]=vparer_tib.at(irun)[2];   if(rxer_tib[irun]>100) rxer_tib[irun]=100;
ryer_tib[irun]=vparer_tib.at(irun)[3];   if(ryer_tib[irun]>100) ryer_tib[irun]=100;
rzer_tib[irun]=vparer_tib.at(irun)[4];   if(rzer_tib[irun]>100) rzer_tib[irun]=100;
 }
//tecp 
 dx_tecp[irun]=vpar_tecp.at(irun)[0];     if(fabs(dx_tecp[irun])>100) dx_tecp[irun]=100;
 dy_tecp[irun]=vpar_tecp.at(irun)[1];      if(fabs(dy_tecp[irun])>100) dy_tecp[irun]=100;
 rz_tecp[irun]=vpar_tecp.at(irun)[2];      if(fabs(rz_tecp[irun])>100) rz_tecp[irun]=100;
 if(vparer_tecp.size()>0) {
dxer_tecp[irun]=vparer_tecp.at(irun)[0];  if(dxer_tecp[irun]>100) dxer_tecp[irun]=100;
dyer_tecp[irun]=vparer_tecp.at(irun)[1];  if(dyer_tecp[irun]>100) dyer_tecp[irun]=100;
rzer_tecp[irun]=vparer_tecp.at(irun)[2];  if(rzer_tecp[irun]>100) rzer_tecp[irun]=100;
 }
//tecm
dx_tecm[irun]=vpar_tecm.at(irun)[0];      if(fabs(dx_tecm[irun])>100) dx_tecm[irun]=100;
 dy_tecm[irun]=vpar_tecm.at(irun)[1];      if(fabs(dy_tecm[irun])>100) dy_tecm[irun]=100;
 rz_tecm[irun]=vpar_tecm.at(irun)[2];     if(fabs(rz_tecm[irun])>100) rz_tecm[irun]=100;
 if(vparer_tecm.size()>0) {
dxer_tecm[irun]=vparer_tecm.at(irun)[0]; if(dxer_tecm[irun]>100) dxer_tecm[irun]=100;
dyer_tecm[irun]=vparer_tecm.at(irun)[1]; if(dyer_tecm[irun]>100) dyer_tecm[irun]=100;
rzer_tecm[irun]=vparer_tecm.at(irun)[2]; if(rzer_tecm[irun]>100) rzer_tecm[irun]=100;
 }


 // tecpint
 dx_tecpint[irun]=vpar_tecpint.at(irun)[0];     if(fabs(dx_tecpint[irun])>100) dx_tecpint[irun]=100;
 dy_tecpint[irun]=vpar_tecpint.at(irun)[1];      if(fabs(dy_tecpint[irun])>100) dy_tecpint[irun]=100;
 rz_tecpint[irun]=vpar_tecpint.at(irun)[2];      if(fabs(rz_tecpint[irun])>100) rz_tecpint[irun]=100;
 if(vparer_tecpint.size()>0) {
dxer_tecpint[irun]=vparer_tecpint.at(irun)[0];  if(dxer_tecpint[irun]>100) dxer_tecpint[irun]=100;
dyer_tecpint[irun]=vparer_tecpint.at(irun)[1];  if(dyer_tecpint[irun]>100) dyer_tecpint[irun]=100;
rzer_tecpint[irun]=vparer_tecpint.at(irun)[2];  if(rzer_tecpint[irun]>100) rzer_tecpint[irun]=100;
 }
//tecmint
dx_tecmint[irun]=vpar_tecmint.at(irun)[0];      if(fabs(dx_tecmint[irun])>100) dx_tecmint[irun]=100;
dy_tecmint[irun]=vpar_tecmint.at(irun)[1];      if(fabs(dy_tecmint[irun])>100) dy_tecmint[irun]=100;
rz_tecmint[irun]=vpar_tecmint.at(irun)[2];     if(fabs(rz_tecmint[irun])>100) rz_tecmint[irun]=100;
 if(vparer_tecmint.size()>0) {
dxer_tecmint[irun]=vparer_tecmint.at(irun)[0]; if(dxer_tecmint[irun]>100) dxer_tecmint[irun]=100;
dyer_tecmint[irun]=vparer_tecmint.at(irun)[1]; if(dyer_tecmint[irun]>100) dyer_tecmint[irun]=100;
rzer_tecmint[irun]=vparer_tecmint.at(irun)[2]; if(rzer_tecmint[irun]>100) rzer_tecmint[irun]=100;
 }


sign[irun]=vsign[irun];       if(sign[irun]>100) sign[irun]=100;
ch2[irun]=vch2[irun];       if(ch2[irun]>1000) ch2[irun]=1000;
ch2tib[irun]=vch2tib[irun];  if(ch2tib[irun]>1000) ch2tib[irun]=1000;
ch2tecp[irun]=vch2tecp[irun]; if(ch2tecp[irun]>1000) ch2tecp[irun]=1000;
ch2tecm[irun]=vch2tecm[irun]; if(ch2tecm[irun]>1000) ch2tecm[irun]=1000;
 }//irun
// fill graphs
 if(btime) {
   bstabt=1;
gt_tbstab[0]=new TGraphErrors(nruns,dtime,dx_tib,dum,dxer_tib);
gt_tbstab[1]=new TGraphErrors(nruns,dtime,dy_tib,dum,dyer_tib);
gt_tbstab[2]=new TGraphErrors(nruns,dtime,rx_tib,dum,rxer_tib);
gt_tbstab[3]=new TGraphErrors(nruns,dtime,ry_tib,dum,ryer_tib);
gt_tbstab[4]=new TGraphErrors(nruns,dtime,rz_tib,dum,rzer_tib);


gt_tecpstab[0]=new TGraphErrors(nruns,dtime,dx_tecp,dum,dxer_tecp);
gt_tecpstab[1]=new TGraphErrors(nruns,dtime,dy_tecp,dum,dyer_tecp);
gt_tecpstab[2]=new TGraphErrors(nruns,dtime,rz_tecp,dum,rzer_tecp);

gt_tecmstab[0]=new TGraphErrors(nruns,dtime,dx_tecm,dum,dxer_tecm);
gt_tecmstab[1]=new TGraphErrors(nruns,dtime,dy_tecm,dum,dyer_tecm);
gt_tecmstab[2]=new TGraphErrors(nruns,dtime,rz_tecm,dum,rzer_tecm);

gt_tecpintstab[0]=new TGraphErrors(nruns,dtime,dx_tecpint,dum,dxer_tecpint);
gt_tecpintstab[1]=new TGraphErrors(nruns,dtime,dy_tecpint,dum,dyer_tecpint);
gt_tecpintstab[2]=new TGraphErrors(nruns,dtime,rz_tecpint,dum,rzer_tecpint);

gt_tecmintstab[0]=new TGraphErrors(nruns,dtime,dx_tecmint,dum,dxer_tecmint);
gt_tecmintstab[1]=new TGraphErrors(nruns,dtime,dy_tecmint,dum,dyer_tecmint);
gt_tecmintstab[2]=new TGraphErrors(nruns,dtime,rz_tecmint,dum,rzer_tecmint);


gt_sign=new TGraphErrors(nruns,dtime,sign,dum,dum);
gt_ch2=new TGraphErrors(nruns,dtime,ch2,dum,dum);
gt_ch2tb=new TGraphErrors(nruns,dtime,ch2tib,dum,dum);
gt_ch2tecp=new TGraphErrors(nruns,dtime,ch2tecp,dum,dum);
gt_ch2tecm=new TGraphErrors(nruns,dtime,ch2tecm,dum,dum);
 } 

 if(brun){
  bstab=1;
g_tbstab[0]=new TGraphErrors(nruns,drun,dx_tib,dum,dxer_tib);
g_tbstab[1]=new TGraphErrors(nruns,drun,dy_tib,dum,dyer_tib);
g_tbstab[2]=new TGraphErrors(nruns,drun,rx_tib,dum,rxer_tib);
g_tbstab[3]=new TGraphErrors(nruns,drun,ry_tib,dum,ryer_tib);
g_tbstab[4]=new TGraphErrors(nruns,drun,rz_tib,dum,rzer_tib);


g_tecpstab[0]=new TGraphErrors(nruns,drun,dx_tecp,dum,dxer_tecp);
g_tecpstab[1]=new TGraphErrors(nruns,drun,dy_tecp,dum,dyer_tecp);
g_tecpstab[2]=new TGraphErrors(nruns,drun,rz_tecp,dum,rzer_tecp);

g_tecmstab[0]=new TGraphErrors(nruns,drun,dx_tecm,dum,dxer_tecm);
g_tecmstab[1]=new TGraphErrors(nruns,drun,dy_tecm,dum,dyer_tecm);
g_tecmstab[2]=new TGraphErrors(nruns,drun,rz_tecm,dum,rzer_tecm);


g_tecpintstab[0]=new TGraphErrors(nruns,drun,dx_tecpint,dum,dxer_tecpint);
g_tecpintstab[1]=new TGraphErrors(nruns,drun,dy_tecpint,dum,dyer_tecpint);
g_tecpintstab[2]=new TGraphErrors(nruns,drun,rz_tecpint,dum,rzer_tecpint);

g_tecmintstab[0]=new TGraphErrors(nruns,drun,dx_tecmint,dum,dxer_tecmint);
g_tecmintstab[1]=new TGraphErrors(nruns,drun,dy_tecmint,dum,dyer_tecmint);
g_tecmintstab[2]=new TGraphErrors(nruns,drun,rz_tecmint,dum,rzer_tecmint);


g_sign=new TGraphErrors(nruns,drun,sign,dum,dum);

g_ch2=new TGraphErrors(nruns,drun,ch2,dum,dum);
g_ch2tb=new TGraphErrors(nruns,drun,ch2tib,dum,dum);
g_ch2tecp=new TGraphErrors(nruns,drun,ch2tecp,dum,dum);
g_ch2tecm=new TGraphErrors(nruns,drun,ch2tecm,dum,dum);


 } // btime


 if(debug) cout<<"  >>>>>LAS_histo::filleresults "<<endl;
  return 0;
}


// clear all histo
void LAS_histo::reset() {

}
//
int LAS_histo::write(string fname) {
  if(debug) cout<<"  LAS_histo::write histo>>> "<<fname<<endl;
TFile*  af   = new TFile(fname.c_str(), "RECREATE" ); 
af->cd();
// store means and rmss
TDirectory* dir=af->mkdir("Positions");
dir->cd();

TDirectory* dir1=dir->mkdir("Rms");
dir->cd();
TDirectory* dir2=dir->mkdir("Mean");
dir->cd();
TDirectory* dir3=dir->mkdir("Residuals");
dir->cd();
TDirectory* dir4=dir->mkdir("Statisics");
dir->cd();
TDirectory* dir5=dir->mkdir("ATparameters");

af->cd();
TDirectory* dir10=af->mkdir("Results");

 int nbm=8;
 int npartib=5;
int npartec=3;
 if(bh) {
for(int i=0;i<nbm;i++) {

dir1->cd();
 h2_tobrms[i]->Write();
 h2_tibrms[i]->Write();
 h2_tecprms[i]->Write();
 h2_tecmrms[i]->Write();
dir->cd();
 
dir2->cd();
 h2_tobmean[i]->Write();
 // cout<<" wite "<<h2_tobmean[i]->GetName();
 h2_tibmean[i]->Write();
 h2_tecpmean[i]->Write();
 h2_tecmmean[i]->Write();

 h2_tecpmeanint[0][i]->Write();
 h2_tecmmeanint[0][i]->Write();
 h2_tecpmeanint[1][i]->Write();
 h2_tecmmeanint[1][i]->Write();


dir->cd();
dir3->cd();
 h2_tobres[i]->Write();
 h2_tibres[i]->Write();
 h2_tecpres[i]->Write();
 h2_tecmres[i]->Write();
dir->cd();
dir4->cd();
h2_tobnorm[i]->Write();
h2_tibnorm[i]->Write();
h2_tecpnorm[i]->Write();
h2_tecmnorm[i]->Write();

h1_tobstat[i]->Write();
h1_tibstat[i]->Write();
h1_tecpstat[i]->Write();
h1_tecmstat[i]->Write();
dir->cd();
dir5->cd();
h1_ata[i]->Write();
h1_atb[i]->Write();
h1_atch2[i]->Write();
h1_atn[i]->Write();
 }
// results

dir10->cd();

 for(int i=0;i<npartib;i++) {
   h1_tibpar[i]->Write();
   h1_tibepar[i]->Write();
   h1_tobpar[i]->Write();
   h1_tobepar[i]->Write();
 }


 for(int i=0;i<npartec;i++) {
h1_tecppar[i]->Write();
h1_tecpepar[i]->Write();
h1_tecmpar[i]->Write();
h1_tecmepar[i]->Write();

h1_tecpintpar[i]->Write();
h1_tecpeintpar[i]->Write();
h1_tecmintpar[i]->Write();
h1_tecmeintpar[i]->Write();


  }

 h1_tobch2->Write();
 h1_tibch2->Write();
 h1_tecpch2->Write();
 h1_tecmch2->Write();
 h1_ch2->Write();
 h1_sign->Write();

h2_tobstat->Write();
h2_tibstat->Write();
h2_tecpstat->Write();
h2_tecmstat->Write();

 }// if booked

 if(debug>5) cout<<" write graph "<<endl;
 // graphs

dir10->cd();

char si[1];
string ss;
 string name;


 if(bstab) {
for(int i=0;i<npartib;i++) {
sprintf(si,"%d",i);
string ss(si);
name="g_tb_"+ss+"p";
 if(g_tbstab[i])
g_tbstab[i]->Write(name.c_str()); 
 }

for(int i=0;i<npartec;i++) {
sprintf(si,"%d",i);
string ss(si);
name="g_tecp_"+ss+"p";
 if(g_tecpstab[i])
g_tecpstab[i]->Write(name.c_str());
name="g_tecm_"+ss+"p";
 if(g_tecmstab[i])
g_tecmstab[i]->Write(name.c_str()); 
 
name="g_tecpint_"+ss+"p";
 if(g_tecpintstab[i])
g_tecpintstab[i]->Write(name.c_str());
name="g_tecmint_"+ss+"p";
 if(g_tecmintstab[i])
g_tecmintstab[i]->Write(name.c_str()); 



}
 if(g_sign) g_sign->Write("g_sign");
 if(g_ch2) g_ch2->Write("g_ch2");
 if(g_ch2tb) g_ch2tb->Write("g_ch2tb");
 if(g_ch2tecp) g_ch2tecp->Write("g_ch2tecp");
 if(g_ch2tecm) g_ch2tecm->Write("g_ch2tecm");

 }// bstab

if(bstabt) {
for(int i=0;i<npartib;i++) {
sprintf(si,"%d",i);
string ss(si);
name="gt_tb_"+ss+"p";
 if(gt_tbstab[i])
gt_tbstab[i]->Write(name.c_str()); 
 }
for(int i=0;i<npartec;i++) {
sprintf(si,"%d",i);
string ss(si);
name="gt_tecp_"+ss+"p";
 if(gt_tecpstab[i])
gt_tecpstab[i]->Write(name.c_str());
name="gt_tecm_"+ss+"p";
 if(gt_tecmstab[i])
gt_tecmstab[i]->Write(name.c_str()); 
 
name="gt_tecpint_"+ss+"p";
 if(gt_tecpintstab[i])
gt_tecpintstab[i]->Write(name.c_str());
name="gt_tecmint_"+ss+"p";
 if(gt_tecmintstab[i])
gt_tecmintstab[i]->Write(name.c_str()); 

}
 if(gt_sign) gt_sign->Write("gt_sign");
 if(gt_ch2) gt_ch2->Write("gt_ch2");
 if(gt_ch2tb) gt_ch2tb->Write("gt_ch2tb");
 if(gt_ch2tecp) gt_ch2tecp->Write("gt_ch2tecp");
 if(gt_ch2tecm) gt_ch2tecm->Write("gt_ch2tecm");

 }// bstabt


if(bstabtt) {
for(int i=0;i<npartib;i++) {
sprintf(si,"%d",i);
string ss(si);
name="gtt_tb_"+ss+"p";
 if(gtt_tbstab[i])
gtt_tbstab[i]->Write(name.c_str()); 
 }
for(int i=0;i<npartec;i++) {
sprintf(si,"%d",i);
string ss(si);
name="gtt_tecp_"+ss+"p";
 if(gtt_tecpstab[i])
gtt_tecpstab[i]->Write(name.c_str());
name="gtt_tecm_"+ss+"p";
 if(gtt_tecmstab[i])
gtt_tecmstab[i]->Write(name.c_str()); 
 

name="gtt_tecpint_"+ss+"p";
 if(gtt_tecpintstab[i])
gtt_tecpintstab[i]->Write(name.c_str());
name="gtt_tecmint_"+ss+"p";
 if(gtt_tecmintstab[i])
gtt_tecmintstab[i]->Write(name.c_str()); 


}

 if(gtt_sign) gtt_sign->Write("gtt_sign");
 if(gtt_ch2) gtt_ch2->Write("gtt_ch2");
 }// bstabtt



af->cd();
af->Close();
 if(debug) cout<<"  LAS>>>_histo::write done  "<<endl;
 return 0;
}

