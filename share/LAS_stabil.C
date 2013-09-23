#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <functional>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMatrixT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"

#define __AVECROOT__
#include "Avec.h"

#include "LAS_stabil.h"
#include "LAS_histo.h"
#include "LASGlobalData.h"
#include "LASGlobalDataLoop.h"
#include "LAS_vectorfloat_tools.h"
#include "LAS_basic_tools.h"
#include "LAS_globaldata_tools.h"
#include "LAS_RDC_tools.h"

using namespace std;

//
LAS_stabil::LAS_stabil() {
init_geom();
  nruns=-1;
  debug=0;
  nrefrun=0;
  refrunn=1;
  bhist=0; 
  bmc_mean=0;
  bmc_at=0;
  bcor_at=1;
  stat_at=0;
  stat_tib=0;
  stat_tecp=0;
  stat_tecm=0;
 
  battob=0;
  btecint=0; // 1 use the local TECint, 2 use RDC
 mcsig_mean=0.5; // just eg
 mcsigtob_mean=0.5; // just eg
 mcsigtib_mean=0.5; // just eg
 mcsigtec_mean=0.5; // just eg

 // mcsig_ata=3.0e-5;//from 1e-5 3e-5
 // mcsig_atb=0.01; //from 0.004 to 0.01

  mcsig_ata=1.0e-5;//from 1e-5 3e-5
  mcsig_atb=0.005; //from 0.004 to 0.01
  // for internal TEC
  mcsig_atec=0.5e-5;
  mcsig_btec=0.001;

 maxrms=2.;
  //
 thernd=new TRandom();
 //
thehisto=0;


//default selection
//Cdiffmax =0.05;// 40um 
Cdiffmax =0.02;// 40um 
Cmin=100.;
Cmax=400.;
Crms=1.0;
Cnpm=3;
Cch2at=10.0;
Cblock=5;
}

//
LAS_stabil::~LAS_stabil() {
  if(thehisto) delete thehisto;
}

//geometry constants
void LAS_stabil::init_geom() {
 
  pch_tib = 0.120;
  rad_tib = 514.;
  pch_tob = 0.183;
  rad_tob = 600.;

  pch_tec = 0.126;
  rad_tec = 564.;


  pch_tecint[0] = 0.126;
  rad_tecint[0] = 564.;
  pch_tecint[1] = 0.189;
  rad_tecint[1] = 840.;


  izp_tob[0] = 1040.;
  izp_tob[1] =  580.;
  izp_tob[2] =  220.;
  izp_tob[3] = -140.;
  izp_tob[4] = -500.;
  izp_tob[5] = -860.;

  izp_tib[0] =  620.;
  izp_tib[1] =  380.;
  izp_tib[2] =  180.;
  izp_tib[3] = -100.;
  izp_tib[4] = -340.;
  izp_tib[5] = -540.;
 
  izp_tecp[0] = 1322.5;
  izp_tecp[1] = 1462.5;
  izp_tecp[2] = 1602.5;
  izp_tecp[3] = 1742.5;
  izp_tecp[4] = 1882.5;
  izp_tecp[5] = 2057.5;
  izp_tecp[6] = 2247.5;
  izp_tecp[7] = 2452.5;
  izp_tecp[8] = 2667.5;

  izp_tecm[0] = -1322.5;
  izp_tecm[1] = -1462.5;
  izp_tecm[2] = -1602.5;
  izp_tecm[3] = -1742.5;
  izp_tecm[4] = -1882.5;
  izp_tecm[5] = 2057.5;
  izp_tecm[6] = 2247.5;
  izp_tecm[7] = 2452.5;
  izp_tecm[8] = 2667.5;

  ii[0][0] = -1;
  ii[0][1] = 1;
  ii[0][2] = 1;
  ii[0][3] = -1;
  ii[0][4] = -1;
  ii[0][5] = 1;
  ii[1][0] = -1;
  ii[1][1] = 1;
  ii[1][2] = 1;
  ii[1][3] = -1;
  ii[1][4] = -1;
  ii[1][5] = 1;
  ii[2][0] = -1;
  ii[2][1] = 1;
  ii[2][2] = 1;
  ii[2][3] = -1;
  ii[2][4] = -1;
  ii[2][5] = 1;
  ii[3][0] = -1;
  ii[3][1] = 1;
  ii[3][2] = 1;
  ii[3][3] = -1;
  ii[3][4] = -1;
  ii[3][5] = 1;
  ii[4][0] = -1;
  ii[4][1] = 1;
  ii[4][2] = 1;
  ii[4][3] = -1;
  ii[4][4] = -1;
  ii[4][5] = 1;
  ii[5][0] = -1;
  ii[5][1] = 1;
  ii[5][2] = 1;
  ii[5][3] = -1;
  ii[5][4] = -1;
  ii[5][5] = 1;
  ii[6][0] = -1;
  ii[6][1] = 1;
  ii[6][2] = 1;
  ii[6][3] = -1;
  ii[6][4] = -1;
  ii[6][5] = 1;
  ii[7][0] = -1;
  ii[7][1] = 1;
  ii[7][2] = 1;
  ii[7][3] = -1;
  ii[7][4] = -1;
  ii[7][5] = 1;

  phib[0] = 0.393;
  phib[1] = 1.290;
  phib[2] = 1.889;
  phib[3] = 2.786;
  phib[4] = 3.684;
  phib[5] = 4.282;
  phib[6] = 5.180;
  phib[7] = 5.778;
  return;
}

// prepare file list
int LAS_stabil::readflist(const string& filelist, const string& spref) {
  if(debug) cout<<" read filelist>>> "<<filelist<<endl;
ifstream infile (filelist.c_str());
 if(!infile)  {
    cout<<" *********cannot open  file "<<filelist<<endl;;
    return 1;
  }

sfiles.clear();
string sfile;
while(infile>>sfile) {
   string sf=spref+sfile;
   // check file
 int flag=0;
  TFile f_res(sf.c_str());
    if(!f_res.IsOpen()){
      cout << "Could not open file " << sf  << endl;
      return 2;
    }
  Avec unixTime=avec_get("unixTime", f_res);
    if(unixTime.empty()){
      if(debug>2) cout << "Could not find unix Time " <<sf<< endl;
      flag=3;
    }
 
    if(flag==0) {
    LASGlobalData<double> *mean_ptr=0;    
    f_res.GetObject("positions_strips", mean_ptr);   
    LASGlobalData<double> *rms_ptr=0;   
    f_res.GetObject("rms_av_0", rms_ptr);   
    if(mean_ptr == 0){
      if(debug>2) cout << "Could not find mean " <<sf<< endl;
      flag=4;
    }
    if(rms_ptr == 0){
      if(debug>2) cout << "Could not find rms " << sf<<endl;
	flag=5;
     }

 Avec block_nr=avec_get("block_nr", f_res);
 if(block_nr.empty()) flag=6;
 else{
 if(debug>2) cout<<" blocknrsize="<<block_nr.size()<<endl; 
if(block_nr.size()<Cblock) flag=7;
 }
    } // flag==0
    if(flag) continue;

    else {
sfiles.push_back(sf);
if(debug>2) cout<<" added= "<<sf.c_str()<<endl;
vtime.push_back(unixTime[1]);
// not yet stored

//Avec runnum=avec_get("runnumber", f_res);
// if(runnum.empty()) {
   if(debug) cout<<"extract runnumber from fixed name ";
   // if not yet implimented, decode from the name
string sfo=sfile;
//int n1= sfo.find("RUN",0);
int n1= sfo.find("run",0);
sfo.erase(0,n1+3);
int n2= sfo.find("_",0);
 sfo.erase(n2,sfo.size());
 int runn=atoi(sfo.c_str());
 if(debug>2) cout<<sfo<<" run= "<<runn<<endl;
vruns.push_back(runn);
// } else {

//vruns.push_back(runnum[1]);
// }
    }// else
 }// while
nruns=sfiles.size();
 if(nruns==0) return 10;
// if refernce is in
 if(refrunn>0) {
 bool bfref=false; 
for(int irun = 0; irun < nruns;  ++irun){
  if(debug>10) cout<<" search ref run "<<vruns[irun]<<" "<<refrunn<<endl;
  if(vruns[irun]==refrunn) bfref=true;
 }

 if(!bfref) {cout<<"****** reference run "<<refrunn<<" not found after selecion in "<<filelist<<endl; return 100;}
 }// if refrun is defined (negative means averaging)

  if(debug>1) {
    cout<<" list to be analyzed "<<endl;
for(int irun = 0; irun < nruns;  ++irun){
  cout<<" file ="<<sfiles[irun]<<" "<<vruns[irun]<<endl;
     }
  }

 if(debug) cout<<" >>>read filelist done"<<filelist<<endl;
  return 0;
}


//
int LAS_stabil::alignATtob() {
  if(debug) cout<<" alignATtob>>>"<<endl;

 if(debug>5) {
    cout<<" list to be analyzed alignAT "<<endl;
for(int irun = 0; irun < nruns;  ++irun){
  cout<<" file ="<<sfiles[irun]<<endl;
 }
 }
  //
 clearv(&varo);
 clearv(&vbro);
 clearv(&vaaro);

clearv(&vch2_at);
clearv(&vndf_at);
  // main loop 
 // tob
 int nbmax=8;
 int ndmax=6;

int  nav[nbmax][ndmax];
double mean0[nbmax][ndmax];
double rms0[nbmax][ndmax];
double mean[nbmax][ndmax];
double rms[nbmax][ndmax];

 if(refrunn<0) {
 for(int i=0;i<nbmax;i++)  
   for(int j=0;j<ndmax;j++) {
    mean0[i][j]=0; 
    nav[i][j]=0; 
       }
 }

LASGlobalDataLoop* loop=0;


// first loop , is the reference, independent on time

for(int irun = 0; irun < nruns;  ++irun){
 if(debug>2)  cout<<" 1read file ="<<sfiles[irun]<<endl;
 TFile f_res(sfiles[irun].c_str());
LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);   
// read LAS
loop=new LASGlobalDataLoop(LASGlobalDataLoop::TOB);
    do{
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
  if(bn>nbmax) continue;
  if(dn>ndmax) continue;  
  if(refrunn<0) {
    double m= loop->GetEntry<double>( *mean_ptr );
    double r=loop->GetEntry<double>( *rms_ptr );
    	if(m > Cmin && m < Cmax   && r <Crms){  
mean0[bn][dn] =mean0[bn][dn]+m;
nav[bn][dn]++;
	}
  }
 // if(irun==nrefrun) {  
 if(vruns[irun]==refrunn&&refrunn>0) {
   if(debug>20) cout<< "found ref run "<<sfiles[irun]<<endl;
    mean0[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms0[bn][dn] = loop->GetEntry<double>( *rms_ptr );
    if(debug>10) cout<< " mean0["<<bn<<"]["<<dn<<"]="<<mean0[bn][dn]<<endl;;
 }//ref
 }while(loop->next()); 
    delete loop;
 }//irun

 if(refrunn<0) {
 for(int bn=0; bn<nbmax;bn++) 
    for(int dn=0; dn<ndmax;dn++) 
      { if(nav[bn][dn]!=0) mean0[bn][dn]=mean0[bn][dn]/nav[bn][dn];
      } 
}

    // process all(including reference)
    // selections
// Cdiffmax =0.04;// 40um 
// Cmin=100.;
// Cmax=400.;
// Crms=0.5;
int nminp=Cnpm;

double sx2,sx,sxz,sz2,sz,sn;
 sx2=sx=sxz=sz2=sz=sn=0;
 int np=0;
 double* a;
 double* aa;
 double* b;
 int* p;
 double* ch2_at;
 int* ndf_at;

 // analysis loop
for(int irun = 0; irun < nruns;  ++irun){
  // if(debug>2)  cout<<" anal  file ="<<sfiles[irun]<<endl;
 cout<<" KKKanal  file ="<<sfiles[irun]<<endl;
TFile f_res(sfiles[irun].c_str());

// extract prameters per run
// possible per block
LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);  
LASGlobalData<int> *norm=0;   
f_res.GetObject("norm_0", norm);   

 
// read LAS
// container for params
a=new double[nbmax];
aa=new double[nbmax];
b=new double[nbmax];
p=new int[nbmax];
ch2_at=new double[nbmax];
ndf_at=new int[nbmax];

 cleara(a,nbmax);
 cleara(aa,nbmax);
 cleara(b,nbmax);
 cleara(p,nbmax);
 cleara(ch2_at,nbmax);
 cleara(ndf_at,nbmax);
 stat_at=0;

 int bn1=-1;
 int dn1=-1;
loop=new LASGlobalDataLoop(LASGlobalDataLoop::TOB);
    do{
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
     // boundaries
     if(bn>nbmax) continue;
     if(dn>ndmax) continue;
     if(bn1<0) bn1=bn;   
     if(dn1<0) dn1=dn;   
  
    mean[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms[bn][dn] = loop->GetEntry<double>( *rms_ptr );
    int  anorm = loop->GetEntry<int>( *norm );
  
 
  //hist
    if(bhist) {
      if(debug>11) cout<<" fill histoAT "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[bn][dn]<<" rms="<<rms[bn][dn]<<endl;
      //cout<<" KKKfill histoATtob "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[bn][dn]<<" rms="<<rms[bn][dn]<<endl;
	   thehisto->h2_tobmean[bn]->Fill(dn,mean[bn][dn]);
	   thehisto->h2_tobrms[bn]->Fill(dn,rms[bn][dn]);
	   thehisto->h2_tobnorm[bn]->Fill(dn,anorm);
    }


 if( rms[bn][dn]==0)
      rms[bn][dn]=mcsig_mean;

    if(anorm==0) anorm=1;
       if(anorm<10) 
           rms[bn][dn]=maxrms;

     rms[bn][dn]/=sqrt(anorm);


  // if MC
    if(bmc_mean) {
 mcsmear(&mean[bn][dn],rms[bn][dn]); 

    }
    if(debug>10) cout<<" read "<<bn<<dn<<" bn1="<<bn1<<" mean="<<mean[bn][dn]<<" "<<rms[bn][dn]<<endl;
   
   // new beam(assuming that the disks are stored inside the beams
    // calculate params from the previous one
    if(bn!=bn1||(bn==nbmax-1&&dn==ndmax-1)) {
	if(debug>8) cout<<" get params and change beam "<<bn<<" bn1="<<bn1<<" np="<<np<<" sxz="<<sxz<<" "<<"sx="<<sx<<" sz2="<<sz2<<" sn="<<sn<<" sz="<<sz<<endl;
	// cout<<" KKKKKKget params and change beam "<<bn<<" bn1="<<bn1<<" np="<<np<<" sxz="<<sxz<<" "<<"sx="<<sx<<" sz2="<<sz2<<" sn="<<sn<<" sz="<<sz<<endl;



	if(bn==nbmax-1&&dn==ndmax-1) { 
         bn1=bn; 
	  if(debug>8) cout<<" last beam "<<endl;
	}//last beam

      a[bn1] = 0.;
      b[bn1] = 0.;
      aa[bn1] = 0.;
      p[bn1]=np;
      if(np > nminp){
	a[bn1]= (sxz*sn-sx*sz)/double(sz2*sn-sz*sz);
	b[bn1] = (sx-a[bn]*sz)/sn;
	aa[bn1] = sxz/sz2;
    
	if(debug>5) cout<<" AT params  beam="<<bn1<<" a="<<a[bn1]<<" b="<<b[bn1]<<" aa="<<aa[bn1]<<" p= "<< p[bn1]<<endl;   
	// calculate ch2 for each beam
        double   chisq = 0.;
        int ndf=0;
        for(int idn = 0; idn < ndmax; ++idn){
	  double ddiff =  (mean[bn1][idn]-mean0[bn1][idn])*ii[bn1][idn]*pch_tob; // bug 2208 was ii[dn]
	 double rres = ddiff - a[bn1]*izp_tob[idn] - b[bn1];
         double err = rms[bn1][idn]*rms[bn1][idn];
        if(mean[bn1][idn] > Cmin && mean[bn1][idn] < Cmax   && fabs(ddiff) < Cdiffmax   && rms[bn1][idn]<Crms ){  
          chisq+=rres*rres/err;
	  ndf++;
	 }
	}//idn
	if(debug>5) cout<<" calulate AT ch2 for beam "<<bn1<<" ch2 "<<chisq<<" "<<ndf<<endl;
	ch2_at[bn1]=chisq;
        ndf-=2; // a,b
        ndf_at[bn1]=ndf;
	if(ndf==0){ if(debug) cout<<" ndf<3 continue"<<endl; // should never happen
	  ndf_at[bn1]=-1;
	}
        double ch2n=chisq/(ndf);
	if(ch2n!=ch2n) ch2n=1000;
	if(bhist) {
	     thehisto->h1_ata[bn1]->Fill(a[bn1]);
	     thehisto->h1_atb[bn1]->Fill(b[bn1]);
	     thehisto->h1_atch2[bn1]->Fill(ch2n);
	     thehisto->h1_atn[bn1]->Fill(p[bn1]);
	     //     if(p[bn1]>5) cout<<"KKKKKKKKK shit "<<p[bn1]<<endl;
	     // cout<<" KKKfill histoAT params "<<bn1<<"  a= "<<a[bn1]<<"  b="<<b[bn1]<<" np="<<p[bn1]<<" ch2= "<<ch2n<<endl;

	}
	if(ch2n>Cch2at||ch2n<0) {
	  // set parameters to zero if ch2 is bad
        a[bn1] = 0.;
        b[bn1] = 0.;
        aa[bn1] = 0.;
	}

  // clear
      sx2 = 0.;
      sx = 0.;
      sxz = 0.;
      sz = 0.;
      sz2 = 0.;
      sn = 0.;
      np = 0;
      bn1=bn; // start new beam
      } // npmin
 else {stat_at++;
      sx2 = 0.;
      sx = 0.;
      sxz = 0.;
      sz = 0.;
      sz2 = 0.;
      sn = 0.;
      np = 0;
      bn1=bn; // start new beam
}// not enough points
      }//ifbn

//residuals loop over dn
 
	double er=rms[bn][dn]*rms[bn][dn];
	double diff =  (mean[bn][dn]-mean0[bn][dn])*ii[bn][dn]*pch_tob;
	if(mean[bn][dn] > Cmin && mean[bn][dn] < Cmax   && fabs(diff) < Cdiffmax && rms[bn][dn]<Crms  ){  

 if(bhist) {
    thehisto->h1_tobstat[bn]->Fill(dn);
    thehisto->h2_tobstat->Fill(dn,bn);
    thehisto->h2_tobres[bn]->Fill(dn,diff);
  }
	  np = np + 1;  
	  sn = sn + 1/er;
	  sz2 = sz2 + izp_tob[dn]*izp_tob[dn]/er;
	  sz = sz + izp_tob[dn]/er;
	  sxz = sxz + izp_tob[dn]*diff/er;
	  sx = sx + diff/er;
	  if(debug>8) cout<<" select prof  bn: "<<bn<<"  dn: "<<dn<<"  sn: "<<sn<<"  sz2: "<<sz2<<"  izp[dn]: "<<izp_tob[dn]<<"  sz: "<<sz<<"  err: "<<er<<" np "<<np<<endl;   
	}// if
      if(debug>20) cout<<"irun: "<<irun<<"  np: "<<np<<"  bn: "<<bn<<"  sn: "<<sn<<"  sz2: "<<sz2<<"  sz: "<<sz<<endl;

}while(loop->next()); 
    delete loop;

    if(bmc_at) {
      if(debug>5) cout<<"smear AT params "<<endl;
      for(int i=0;i<nbmax;i++) {
	mcsmear(&b[i],mcsig_ata);
	mcsmear(&b[i],mcsig_atb);
      }
    }//bmc

    // fill ch2
 vch2_at.push_back(ch2_at);
 vndf_at.push_back(ndf_at);
 // fill res

 varo.push_back(a);
 vaaro.push_back(aa);
 vbro.push_back(b); 
 vnp.push_back(p); 

 

 }// irun

 if(debug>1) {
    cout<<" ===TOB AT params "<<endl;
    for(unsigned int  i=0;i<nruns;i++) {
      cout<<sfiles[i]<<" "<<vruns[i]<<" time="<<vtime[i]<<endl;
      for(int  j=0;j<nbmax;j++)
	cout<<"beam="<<j<<" a= "<<varo.at(i)[j]<<" b= "<<vbro.at(i)[j]<<" np="<<vnp.at(i)[j]<<endl;
      }

  }// debug
  if(debug) cout<<" >>>alignAT done"<<endl;
  return 0;
}


int LAS_stabil::alignblockATtob() {
  if(debug) cout<<" alignblockATtob>>>"<<endl;

 if(debug>5) {
    cout<<" list to be analyzed alignAT "<<endl;
for(int irun = 0; irun < nruns;  ++irun){
  cout<<" file ="<<sfiles[irun]<<endl;
 }
 }
  //
 clearv(&varo);
 clearv(&vbro);
 clearv(&vaaro);

clearv(&vch2_at);
clearv(&vndf_at);
  // main loop 
 // tob
 int nbmax=8;
 int ndmax=6;

int  nav[nbmax][ndmax];
double mean0[nbmax][ndmax];
double rms0[nbmax][ndmax];
double mean[nbmax][ndmax];
double rms[nbmax][ndmax];

 if(refrunn<0) {
 for(int i=0;i<nbmax;i++)  
   for(int j=0;j<ndmax;j++) {
    mean0[i][j]=0; 
    nav[i][j]=0; 
       }
 }

LASGlobalDataLoop* loop=0;


// first loop , is the reference, independent on time

for(int irun = 0; irun < nruns;  ++irun){


 if(debug>2)  cout<<" 1read file ="<<sfiles[irun]<<endl;
 TFile f_res(sfiles[irun].c_str());
LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);   
// read LAS
loop=new LASGlobalDataLoop(LASGlobalDataLoop::TOB);
    do{
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
  if(bn>nbmax) continue;
  if(dn>ndmax) continue;  
  if(refrunn<0) {
    double m= loop->GetEntry<double>( *mean_ptr );
    double r=loop->GetEntry<double>( *rms_ptr );
    	if(m > Cmin && m < Cmax   && r <Crms){  
mean0[bn][dn] =mean0[bn][dn]+m;
nav[bn][dn]++;
	}
  }
 // if(irun==nrefrun) {  
 if(vruns[irun]==refrunn&&refrunn>0) {
   if(debug>20) cout<< "found ref run "<<sfiles[irun]<<endl;
    mean0[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms0[bn][dn] = loop->GetEntry<double>( *rms_ptr );
    if(debug>10) cout<< " mean0["<<bn<<"]["<<dn<<"]="<<mean0[bn][dn]<<endl;;
 }//ref
 }while(loop->next()); 
    delete loop;
 

}//irun

 if(refrunn<0) {
 for(int bn=0; bn<nbmax;bn++) 
    for(int dn=0; dn<ndmax;dn++) 
      { if(nav[bn][dn]!=0) mean0[bn][dn]=mean0[bn][dn]/nav[bn][dn];
      } 
}

    // process all(including reference)
    // selections
// Cdiffmax =0.04;// 40um 
// Cmin=100.;
// Cmax=400.;
// Crms=0.5;
int nminp=Cnpm;

double sx2,sx,sxz,sz2,sz,sn;
 sx2=sx=sxz=sz2=sz=sn=0;
 int np=0;
 double* a;
 double* aa;
 double* b;
 int* p;
 double* ch2_at;
 int* ndf_at;

 // analysis loop
for(int irun = 0; irun < nruns;  ++irun){
  // if(debug>2)  cout<<" anal  file ="<<sfiles[irun]<<endl;
 cout<<" KKKanal  file ="<<sfiles[irun]<<endl;
TFile f_res(sfiles[irun].c_str());

// extract prameters per run
// possible per block
LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);  
LASGlobalData<int> *norm=0;   
f_res.GetObject("norm_0", norm);   

 
// read LAS
// container for params
a=new double[nbmax];
aa=new double[nbmax];
b=new double[nbmax];
p=new int[nbmax];
ch2_at=new double[nbmax];
ndf_at=new int[nbmax];

 cleara(a,nbmax);
 cleara(aa,nbmax);
 cleara(b,nbmax);
 cleara(p,nbmax);
 cleara(ch2_at,nbmax);
 cleara(ndf_at,nbmax);
 stat_at=0;

 int bn1=-1;
 int dn1=-1;
loop=new LASGlobalDataLoop(LASGlobalDataLoop::TOB);
    do{
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
     // boundaries
     if(bn>nbmax) continue;
     if(dn>ndmax) continue;
     if(bn1<0) bn1=bn;   
     if(dn1<0) dn1=dn;   
  
    mean[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms[bn][dn] = loop->GetEntry<double>( *rms_ptr );
    int  anorm = loop->GetEntry<int>( *norm );
  
 
  //hist
    if(bhist) {
      if(debug>11) cout<<" fill histoAT "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[bn][dn]<<" rms="<<rms[bn][dn]<<endl;
      //cout<<" KKKfill histoATtob "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[bn][dn]<<" rms="<<rms[bn][dn]<<endl;
	   thehisto->h2_tobmean[bn]->Fill(dn,mean[bn][dn]);
	   thehisto->h2_tobrms[bn]->Fill(dn,rms[bn][dn]);
	   thehisto->h2_tobnorm[bn]->Fill(dn,anorm);
    }


 if( rms[bn][dn]==0)
      rms[bn][dn]=mcsig_mean;

    if(anorm==0) anorm=1;
       if(anorm<10) 
           rms[bn][dn]=maxrms;

     rms[bn][dn]/=sqrt(anorm);


  // if MC
    if(bmc_mean) {
 mcsmear(&mean[bn][dn],rms[bn][dn]); 

    }
    if(debug>10) cout<<" read "<<bn<<dn<<" bn1="<<bn1<<" mean="<<mean[bn][dn]<<" "<<rms[bn][dn]<<endl;
   
   // new beam(assuming that the disks are stored inside the beams
    // calculate params from the previous one
    if(bn!=bn1||(bn==nbmax-1&&dn==ndmax-1)) {
	if(debug>8) cout<<" get params and change beam "<<bn<<" bn1="<<bn1<<" np="<<np<<" sxz="<<sxz<<" "<<"sx="<<sx<<" sz2="<<sz2<<" sn="<<sn<<" sz="<<sz<<endl;
	// cout<<" KKKKKKget params and change beam "<<bn<<" bn1="<<bn1<<" np="<<np<<" sxz="<<sxz<<" "<<"sx="<<sx<<" sz2="<<sz2<<" sn="<<sn<<" sz="<<sz<<endl;



	if(bn==nbmax-1&&dn==ndmax-1) { 
         bn1=bn; 
	  if(debug>8) cout<<" last beam "<<endl;
	}//last beam

      a[bn1] = 0.;
      b[bn1] = 0.;
      aa[bn1] = 0.;
      p[bn1]=np;
      if(np > nminp){
	a[bn1]= (sxz*sn-sx*sz)/double(sz2*sn-sz*sz);
	b[bn1] = (sx-a[bn]*sz)/sn;
	aa[bn1] = sxz/sz2;
    
	if(debug>5) cout<<" AT params  beam="<<bn1<<" a="<<a[bn1]<<" b="<<b[bn1]<<" aa="<<aa[bn1]<<" p= "<< p[bn1]<<endl;   
	// calculate ch2 for each beam
        double   chisq = 0.;
        int ndf=0;
        for(int idn = 0; idn < ndmax; ++idn){
	  double ddiff =  (mean[bn1][idn]-mean0[bn1][idn])*ii[bn1][idn]*pch_tob; // bug 2208 was ii[dn]
	 double rres = ddiff - a[bn1]*izp_tob[idn] - b[bn1];
         double err = rms[bn1][idn]*rms[bn1][idn];
        if(mean[bn1][idn] > Cmin && mean[bn1][idn] < Cmax   && fabs(ddiff) < Cdiffmax   && rms[bn1][idn]<Crms ){  
          chisq+=rres*rres/err;
	  ndf++;
	 }
	}//idn
	if(debug>5) cout<<" calulate AT ch2 for beam "<<bn1<<" ch2 "<<chisq<<" "<<ndf<<endl;
	ch2_at[bn1]=chisq;
        ndf-=2; // a,b
        ndf_at[bn1]=ndf;
	if(ndf==0){ if(debug) cout<<" ndf<3 continue"<<endl; // should never happen
	  ndf_at[bn1]=-1;
	}
        double ch2n=chisq/(ndf);
	if(ch2n!=ch2n) ch2n=1000;
	if(bhist) {
	     thehisto->h1_ata[bn1]->Fill(a[bn1]);
	     thehisto->h1_atb[bn1]->Fill(b[bn1]);
	     thehisto->h1_atch2[bn1]->Fill(ch2n);
	     thehisto->h1_atn[bn1]->Fill(p[bn1]);
	     //     if(p[bn1]>5) cout<<"KKKKKKKKK shit "<<p[bn1]<<endl;
	     // cout<<" KKKfill histoAT params "<<bn1<<"  a= "<<a[bn1]<<"  b="<<b[bn1]<<" np="<<p[bn1]<<" ch2= "<<ch2n<<endl;

	}
	if(ch2n>Cch2at||ch2n<0) {
	  // set parameters to zero if ch2 is bad
        a[bn1] = 0.;
        b[bn1] = 0.;
        aa[bn1] = 0.;
	}

  // clear
      sx2 = 0.;
      sx = 0.;
      sxz = 0.;
      sz = 0.;
      sz2 = 0.;
      sn = 0.;
      np = 0;
      bn1=bn; // start new beam
      } // npmin
 else {stat_at++;
      sx2 = 0.;
      sx = 0.;
      sxz = 0.;
      sz = 0.;
      sz2 = 0.;
      sn = 0.;
      np = 0;
      bn1=bn; // start new beam
}// not enough points
      }//ifbn

//residuals loop over dn
 
	double er=rms[bn][dn]*rms[bn][dn];
	double diff =  (mean[bn][dn]-mean0[bn][dn])*ii[bn][dn]*pch_tob;
	if(mean[bn][dn] > Cmin && mean[bn][dn] < Cmax   && fabs(diff) < Cdiffmax && rms[bn][dn]<Crms  ){  

 if(bhist) {
    thehisto->h1_tobstat[bn]->Fill(dn);
    thehisto->h2_tobstat->Fill(dn,bn);
    thehisto->h2_tobres[bn]->Fill(dn,diff);
  }
	  np = np + 1;  
	  sn = sn + 1/er;
	  sz2 = sz2 + izp_tob[dn]*izp_tob[dn]/er;
	  sz = sz + izp_tob[dn]/er;
	  sxz = sxz + izp_tob[dn]*diff/er;
	  sx = sx + diff/er;
	  if(debug>8) cout<<" select prof  bn: "<<bn<<"  dn: "<<dn<<"  sn: "<<sn<<"  sz2: "<<sz2<<"  izp[dn]: "<<izp_tob[dn]<<"  sz: "<<sz<<"  err: "<<er<<" np "<<np<<endl;   
	}// if
      if(debug>20) cout<<"irun: "<<irun<<"  np: "<<np<<"  bn: "<<bn<<"  sn: "<<sn<<"  sz2: "<<sz2<<"  sz: "<<sz<<endl;

}while(loop->next()); 
    delete loop;

    if(bmc_at) {
      if(debug>5) cout<<"smear AT params "<<endl;
      for(int i=0;i<nbmax;i++) {
	mcsmear(&b[i],mcsig_ata);
	mcsmear(&b[i],mcsig_atb);
      }
    }//bmc

    // fill ch2
 vch2_at.push_back(ch2_at);
 vndf_at.push_back(ndf_at);
 // fill res

 varo.push_back(a);
 vaaro.push_back(aa);
 vbro.push_back(b); 
 vnp.push_back(p); 

 

 }// irun

 if(debug>1) {
    cout<<" ===TOB AT params "<<endl;
    for(unsigned int  i=0;i<nruns;i++) {
      cout<<sfiles[i]<<" "<<vruns[i]<<" time="<<vtime[i]<<endl;
      for(int  j=0;j<nbmax;j++)
	cout<<"beam="<<j<<" a= "<<varo.at(i)[j]<<" b= "<<vbro.at(i)[j]<<" np="<<vnp.at(i)[j]<<endl;
      }

  }// debug
  if(debug) cout<<" >>>alignAT done"<<endl;
  return 0;
}










//
int LAS_stabil::alignATtib() {
  if(debug) cout<<" alignATtib>>>"<<endl;

 if(debug>5) {
    cout<<" list to be analyzed alignAT "<<endl;
for(int irun = 0; irun < nruns;  ++irun){
  cout<<" file ="<<sfiles[irun]<<endl;
 }
 }
  //
 clearv(&varo);
 clearv(&vbro);
 clearv(&vaaro);

clearv(&vch2_at);
clearv(&vndf_at);
  // main loop 
 // tob
 int nbmax=8;
 int ndmax=6;

int  nav[nbmax][ndmax];
double mean0[nbmax][ndmax];
double rms0[nbmax][ndmax];
double mean[nbmax][ndmax];
double rms[nbmax][ndmax];

 if(refrunn<0) {
 for(int i=0;i<nbmax;i++)  
   for(int j=0;j<ndmax;j++) {
    mean0[i][j]=0; 
    nav[i][j]=0; 
       }
 }

LASGlobalDataLoop* loop=0;


// first loop , is the reference, independent on time

for(int irun = 0; irun < nruns;  ++irun){
 if(debug>2)  cout<<" 1read file ="<<sfiles[irun]<<endl;
 TFile f_res(sfiles[irun].c_str());
LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);   
// read LAS
loop=new LASGlobalDataLoop(LASGlobalDataLoop::TIB);
    do{
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
  if(bn>nbmax) continue;
  if(dn>ndmax) continue;  
  if(refrunn<0) {
    double m= loop->GetEntry<double>( *mean_ptr );
    double r=loop->GetEntry<double>( *rms_ptr );
    	if(m > Cmin && m < Cmax   && r <Crms){  
mean0[bn][dn] =mean0[bn][dn]+m;
nav[bn][dn]++;
	}
  }
 // if(irun==nrefrun) {  
 if(vruns[irun]==refrunn&&refrunn>0) {
   if(debug>20) cout<< "found ref run "<<sfiles[irun]<<endl;
    mean0[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms0[bn][dn] = loop->GetEntry<double>( *rms_ptr );
    if(debug>10) cout<< " mean0["<<bn<<"]["<<dn<<"]="<<mean0[bn][dn]<<endl;;
 }//ref
 }while(loop->next()); 
    delete loop;
 }//irun

 if(refrunn<0) {
 for(int bn=0; bn<nbmax;bn++) 
    for(int dn=0; dn<ndmax;dn++) 
      { if(nav[bn][dn]!=0) mean0[bn][dn]=mean0[bn][dn]/nav[bn][dn];
      } 
}

    // process all(including reference)
    // selections
// Cdiffmax =0.04;// 40um 
// Cmin=100.;
// Cmax=400.;
// Crms=0.5;
int nminp=Cnpm;

double sx2,sx,sxz,sz2,sz,sn;
 sx2=sx=sxz=sz2=sz=sn=0;
 int np=0;
 double* a;
 double* aa;
 double* b;
 int* p;
 double* ch2_at;
 int* ndf_at;

 // analysis loop
for(int irun = 0; irun < nruns;  ++irun){
 if(debug>2)  cout<<" anal  file ="<<sfiles[irun]<<endl;

TFile f_res(sfiles[irun].c_str());

// extract prameters per run
// possible per block
LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);  
LASGlobalData<int> *norm=0;   
f_res.GetObject("norm_0", norm);   

 
// read LAS
// container for params
a=new double[nbmax];
aa=new double[nbmax];
b=new double[nbmax];
p=new int[nbmax];
ch2_at=new double[nbmax];
ndf_at=new int[nbmax];

 cleara(a,nbmax);
 cleara(aa,nbmax);
 cleara(b,nbmax);
 cleara(p,nbmax);
 cleara(ch2_at,nbmax);
 cleara(ndf_at,nbmax);
 stat_at=0;

 int bn1=-1;
 int dn1=-1;
loop=new LASGlobalDataLoop(LASGlobalDataLoop::TIB);
    do{
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
     // boundaries
     if(bn>nbmax) continue;
     if(dn>ndmax) continue;
     if(bn1<0) bn1=bn;   
     if(dn1<0) dn1=dn;   
  
    mean[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms[bn][dn] = loop->GetEntry<double>( *rms_ptr );
    int  anorm = loop->GetEntry<int>( *norm );
  
 
  //hist
    if(bhist) {
      if(debug>11) cout<<" fill histoAT "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[bn][dn]<<" rms="<<rms[bn][dn]<<endl;
	   thehisto->h2_tibmean[bn]->Fill(dn,mean[bn][dn]);
	   thehisto->h2_tibrms[bn]->Fill(dn,rms[bn][dn]);
	   thehisto->h2_tibnorm[bn]->Fill(dn,anorm);
    }


 if( rms[bn][dn]==0)
      rms[bn][dn]=mcsig_mean;

    if(anorm==0) anorm=1;
       if(anorm<10) 
           rms[bn][dn]=maxrms;

     rms[bn][dn]/=sqrt(anorm);


  // if MC
    if(bmc_mean) {
 mcsmear(&mean[bn][dn],rms[bn][dn]); 

    }
    if(debug>10) cout<<" read "<<bn<<dn<<" bn1="<<bn1<<" mean="<<mean[bn][dn]<<" "<<rms[bn][dn]<<endl;
   
   // new beam(assuming that the disks are stored inside the beams
    // calculate params from the previous one
    if(bn!=bn1||(bn==nbmax-1&&dn==ndmax-1)) {
	if(debug>8) cout<<" get params and change beam "<<bn<<" bn1="<<bn1<<" np="<<np<<" sxz="<<sxz<<" "<<"sx="<<sx<<" sz2="<<sz2<<" sn="<<sn<<" sz="<<sz<<endl;

	if(bn==nbmax-1&&dn==ndmax-1) { 
         bn1=bn; 
	  if(debug>8) cout<<" last beam "<<endl;
	}//last beam

      a[bn1] = 0.;
      b[bn1] = 0.;
      aa[bn1] = 0.;
      p[bn1]=np;
      if(np > nminp){
	a[bn1]= (sxz*sn-sx*sz)/double(sz2*sn-sz*sz);
	b[bn1] = (sx-a[bn]*sz)/sn;
	aa[bn1] = sxz/sz2;
    
	if(debug>5) cout<<" AT params  beam="<<bn1<<" a="<<a[bn1]<<" b="<<b[bn1]<<" aa="<<aa[bn1]<<" p= "<< p[bn1]<<endl;   
	// calculate ch2 for each beam
        double   chisq = 0.;
        int ndf=0;
        for(int idn = 0; idn < ndmax; ++idn){
	  double ddiff =  (mean[bn1][idn]-mean0[bn1][idn])*ii[bn1][idn]*pch_tib; //bug 2208, was ii[dn]
	 double rres = ddiff - a[bn1]*izp_tib[idn] - b[bn1];
         double err = rms[bn1][idn]*rms[bn1][idn];
        if(mean[bn1][idn] > Cmin && mean[bn1][idn] < Cmax   && fabs(ddiff) < Cdiffmax   && rms[bn1][idn]<Crms ){  
          chisq+=rres*rres/err;
	  ndf++;
	 }
	}//idn
	if(debug>5) cout<<" calulate AT ch2 for beam "<<bn1<<" ch2 "<<chisq<<" "<<ndf<<endl;
	ch2_at[bn1]=chisq;
        ndf-=2; // a,b
        ndf_at[bn1]=ndf;
	if(ndf==0){ if(debug) cout<<" ndf<3 continue"<<endl; // should never happen
	  ndf_at[bn1]=-1;
	}
        double ch2n=chisq/(ndf);
	if(bhist) {
	     thehisto->h1_ata[bn1]->Fill(a[bn1]);
	     thehisto->h1_atb[bn1]->Fill(b[bn1]);
	     thehisto->h1_atch2[bn1]->Fill(ch2n);
	     thehisto->h1_atn[bn1]->Fill(p[bn1]);
	}
	if(ch2n>Cch2at) {
	  // set parameters to zero if ch2 is bad
        a[bn1] = 0.;
        b[bn1] = 0.;
        aa[bn1] = 0.;
	}

  // clear
      sx2 = 0.;
      sx = 0.;
      sxz = 0.;
      sz = 0.;
      sz2 = 0.;
      sn = 0.;
      np = 0;
      bn1=bn; // start new beam
      } // npmin
 else {stat_at++;
      sx2 = 0.;
      sx = 0.;
      sxz = 0.;
      sz = 0.;
      sz2 = 0.;
      sn = 0.;
      np = 0;
      bn1=bn; //

}// not enough points}
      }//ifbn

//residuals loop over dn
 
	double er=rms[bn][dn]*rms[bn][dn];
	double diff =  (mean[bn][dn]-mean0[bn][dn])*ii[bn][dn]*pch_tib;
	if(mean[bn][dn] > Cmin && mean[bn][dn] < Cmax   && fabs(diff) < Cdiffmax && rms[bn][dn]<Crms  ){  

 if(bhist) {
    thehisto->h1_tibstat[bn]->Fill(dn);
    thehisto->h2_tibstat->Fill(dn,bn);
    thehisto->h2_tibres[bn]->Fill(dn,diff);
  }
	  np = np + 1;  
	  sn = sn + 1/er;
	  sz2 = sz2 + izp_tib[dn]*izp_tib[dn]/er;
	  sz = sz + izp_tib[dn]/er;
	  sxz = sxz + izp_tib[dn]*diff/er;
	  sx = sx + diff/er;
	  if(debug>8) cout<<" select prof  bn: "<<bn<<"  dn: "<<dn<<"  sn: "<<sn<<"  sz2: "<<sz2<<"  izp[dn]: "<<izp_tib[dn]<<"  sz: "<<sz<<"  err: "<<er<<" np "<<np<<endl;   
	}// if
      if(debug>20) cout<<"irun: "<<irun<<"  np: "<<np<<"  bn: "<<bn<<"  sn: "<<sn<<"  sz2: "<<sz2<<"  sz: "<<sz<<endl;

}while(loop->next()); 
    delete loop;

    if(bmc_at) {
      if(debug>5) cout<<"smear AT params "<<endl;
      for(int i=0;i<nbmax;i++) {
	mcsmear(&a[i],mcsig_ata);
	mcsmear(&b[i],mcsig_atb);
      }
    }//bmc

    // fill ch2
 vch2_at.push_back(ch2_at);
 vndf_at.push_back(ndf_at);
 // fill res
 varo.push_back(a);
 vaaro.push_back(aa);
 vbro.push_back(b); 
 vnp.push_back(p); 

 

 }// irun

 if(debug>1) {
    cout<<" ===TOB AT params "<<endl;
    for(unsigned int  i=0;i<nruns;i++) {
      cout<<sfiles[i]<<" "<<vruns[i]<<" time="<<vtime[i]<<endl;
      for(int  j=0;j<nbmax;j++)
	cout<<"beam="<<j<<" a= "<<varo.at(i)[j]<<" b= "<<vbro.at(i)[j]<<" np="<<vnp.at(i)[j]<<endl;
      }

  }// debug
  if(debug) cout<<" >>>alignATtib done"<<endl;
  return 0;
}


int LAS_stabil::alignTIB() {
  if(debug) cout<<"    alignTIB>>>"<<endl;
  //tib config
int nbmax=8;
int ndmax=6;
int npar=5;//dx,dy,rx,ry,rz

int  nav[nbmax][ndmax];
double mean0[nbmax][ndmax];
double rms0[nbmax][ndmax];
double mean[nbmax][ndmax];
double rms[nbmax][ndmax];


if(refrunn<0) {
 for(int i=0;i<nbmax;i++)  
   for(int j=0;j<ndmax;j++) {
    mean0[i][j]=0; 
    nav[i][j]=0; 
       }
 }

 LASGlobalDataLoop* loop=0;
// get reference 
for(int irun = 0; irun < nruns;  ++irun){
   
    TFile f_res(sfiles[irun].c_str());
if(debug>2)  cout<<" read file ="<<sfiles[irun]<<endl;
 
     LASGlobalData<double> *mean_ptr=0;   
    f_res.GetObject("positions_strips", mean_ptr);   
    LASGlobalData<double> *rms_ptr=0;   
    f_res.GetObject("rms_av_0", rms_ptr);   
 
loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TIB);
    do{
      int bn = loop->get_beam();
      int dn = loop->get_zpos(); 
  if(refrunn<0) {
  double m= loop->GetEntry<double>( *mean_ptr );
  double r=loop->GetEntry<double>( *rms_ptr );
    	if(m > Cmin && m < Cmax   && r <Crms){  
mean0[bn][dn] =mean0[bn][dn]+m;
nav[bn][dn]++;
	}
  }// refrun<0
 if(vruns[irun]==refrunn&&refrunn>0) {
   if(debug>10) cout<<" TIB found ref "<<vruns[irun]<<endl;
       // reference 
    mean0[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms0[bn][dn] = loop->GetEntry<double>( *rms_ptr );
     }
    }while(loop->next());   
  } // irun

 if(refrunn<0) {
 for(int bn=0; bn<nbmax;bn++) 
    for(int dn=0; dn<ndmax;dn++) 
       if(nav[bn][dn]!=0) mean0[bn][dn]=mean0[bn][dn]/nav[bn][dn];     
}


double bs[npar];
double res[nbmax][ndmax];
TMatrixT<double> h(npar,npar);
// clean vpartib
clearv(&vpar_tb);


vch2_tb.clear();
vndf_tb.clear();


// runs loop  
  for(int irun = 0; irun < nruns;  ++irun){ 

TFile f_res(sfiles[irun].c_str());
if(debug>2)  cout<<" analyze file TIB ="<<sfiles[irun]<<endl;
// clear each run
 cleara(bs,npar);
 clearh(&h);


LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);   
LASGlobalData<int> *norm=0;   
f_res.GetObject("norm_0", norm); 
 


loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TIB);
    do{
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
     // boundaries
     if(bn>nbmax||bn<0) continue;
     if(dn>ndmax||dn<0) continue;

    mean[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms[bn][dn] = loop->GetEntry<double>( *rms_ptr );
    int  anorm = loop->GetEntry<int>( *norm );
 

    //hist
    if(bhist) {
 if(debug>11) cout<<" fill histoTIB "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[bn][dn]<<" rms="<<rms[bn][dn]<<endl;
       thehisto->h2_tibmean[bn]->Fill(dn,mean[bn][dn]);
       thehisto->h2_tibrms[bn]->Fill(dn,rms[bn][dn]);
       thehisto->h2_tibnorm[bn]->Fill(dn,anorm);
    }

 if( rms[bn][dn]==0)
      rms[bn][dn]=mcsig_mean;

    if(anorm==0) anorm=1;
       if(anorm<10) 
           rms[bn][dn]=maxrms;

     rms[bn][dn]/=sqrt(anorm);


    if(debug>11) cout<<"  after  norm "<<bn<<dn<<" mean="<<mean[bn][dn]<<" "<<rms[bn][dn]<<endl;
      // new beam(assuming that the disks are stored inside the beams
  
//residuals loop over dn
	
	if(rms[bn][dn]==0.) {
	  if(debug>10) cout<<" ***rms are 0  put them to 1 "<<rms[bn][dn]<<endl;
          rms[bn][dn]=mcsig_mean;
	}
        if(anorm<10) 
        rms[bn][dn]=maxrms;

	// smear if MC
	if(bmc_mean) mcsmear(&mean[bn][dn],rms[bn][dn]);
	double er=  rms[bn][dn]* rms[bn][dn];

	double diff = (mean[bn][dn]-mean0[bn][dn])*pch_tib;
	if(bcor_at)
  	res[bn][dn]= diff - varo.at(irun)[bn]*izp_tib[dn] - vbro.at(irun)[bn];
        else 
	  res[bn][dn]= diff;


	if(debug>15) cout<<"bn="<<bn<<"  dn="<<dn<<" diff="<<diff<<" res="<<res[bn][dn]<<endl;
  
	if( mean[bn][dn] > Cmin && mean[bn][dn] < Cmax   && fabs(res[bn][dn]) <  Cdiffmax ){ 

  if(bhist) {
    thehisto->h1_tibstat[bn]->Fill(dn);
    thehisto->h2_tibstat->Fill(dn,bn);
    thehisto->h2_tibres[bn]->Fill(dn,diff);
  }
	  double dphi = res[bn][dn]/rad_tib;
	 	  
	  bs[0] = bs[0] + dphi*sin(phib[bn])*rad_tib/er;
	  bs[1] = bs[1] + dphi*cos(phib[bn])*rad_tib/er;
	  bs[2] = bs[2] + dphi*cos(phib[bn])*rad_tib*izp_tib[dn]/er;
	  bs[3] = bs[3] + dphi*sin(phib[bn])*rad_tib*izp_tib[dn]/er;
	  bs[4] = bs[4] + dphi/er;
// fill matrix
	  if(debug>10) cout<<bn<<dn<<" fill h "<<er<<endl;
	  if(debug>20) printh(h);

	  h(0,0) = h(0,0) - sin(phib[bn])*sin(phib[bn])/er;
	  h(0,1) = h(0,1) + sin(phib[bn])*cos(phib[bn])/er;
	  h(0,2) = h(0,2) + sin(phib[bn])*cos(phib[bn])*izp_tib[dn]/er;
	  h(0,3) = h(0,3) + sin(phib[bn])*sin(phib[bn])*izp_tib[dn]/er;
	  h(0,4) = h(0,4) + sin(phib[bn])*rad_tib/er;
	  
	  h(1,0) = h(1,0) - sin(phib[bn])*cos(phib[bn])/er;
	  h(1,1) = h(1,1) + cos(phib[bn])*cos(phib[bn])/er;
	  h(1,2) = h(1,2) + cos(phib[bn])*cos(phib[bn])*izp_tib[dn]/er;
	  h(1,3) = h(1,3) + sin(phib[bn])*cos(phib[bn])*izp_tib[dn]/er;
	  h(1,4) = h(1,4) + cos(phib[bn])*rad_tib/er;
	  
	  h(2,0) = h(2,0) - sin(phib[bn])*cos(phib[bn])*izp_tib[dn]/er;
	  h(2,1) = h(2,1) + cos(phib[bn])*cos(phib[bn])*izp_tib[dn]/er;
	  h(2,2) = h(2,2) + cos(phib[bn])*cos(phib[bn])*izp_tib[dn]*izp_tib[dn]/er;
	  h(2,3) = h(2,3) + sin(phib[bn])*cos(phib[bn])*izp_tib[dn]*izp_tib[dn]/er;
	  h(2,4) = h(2,4) + cos(phib[bn])*izp_tib[dn]*rad_tib/er;
	  
	  h(3,0) = h(3,0) - sin(phib[bn])*sin(phib[bn])*izp_tib[dn]/er;
	  h(3,1) = h(3,1) + sin(phib[bn])*cos(phib[bn])*izp_tib[dn]/er;
	  h(3,2) = h(3,2) + sin(phib[bn])*cos(phib[bn])*izp_tib[dn]*izp_tib[dn]/er;
	  h(3,3) = h(3,3) + sin(phib[bn])*sin(phib[bn])*izp_tib[dn]*izp_tib[dn]/er;
	  h(3,4) = h(3,4) + sin(phib[bn])*izp_tib[dn]*rad_tib/er;
	  
	  h(4,0) = h(4,0) - sin(phib[bn])/rad_tib/er;
	  h(4,1) = h(4,1) + cos(phib[bn])/rad_tib/er;
	  h(4,2) = h(4,2) + cos(phib[bn])*izp_tib[dn]/rad_tib/er;
	  h(4,3) = h(4,3) + sin(phib[bn])*izp_tib[dn]/rad_tib/er;
	  h(4,4) = h(4,4) + 1./er;
	}// ifselect
}while(loop->next()); 
    delete loop;

    if(debug>5) {
      cout<<irun<<" analTIB bs= ";   
      printa(bs,npar);
    }      

    if(debug>30) {
     cout << " TIB h before invertion----  "  << endl; 
     printh(h); 
    }//debug
    
// find solution
    h.Invert();

  if(debug>30) {
 cout << " TIB h after invertion----  "  << endl; 
   printh(h); 
  }//debug

 
  // fill results
 double* ps=new double[npar];
 
    for(int ip = 0; ip < npar; ++ip){
      ps[ip] = 0.;
      for(int ix = 0; ix < npar; ++ix){
	ps[ip] = ps[ip] + h(ip,ix)*bs[ix];
      }//ix
    }//ip


// calculate chi2
 int  ndf = 0;
  double  chisq = 0.;
    for(int bn = 0; bn < nbmax;  ++bn){
      for(int dn = 0; dn < ndmax;  ++dn){
	if(rms[bn][dn]==0.) {// should not occur
	  if(debug>10) cout<<" ***rms are 0  put them to 1 "<<rms[bn][dn]<<endl;
	  rms[bn][dn]=mcsig_mean;
	}
	double er=rms[bn][dn]*rms[bn][dn]/rad_tib/rad_tib; // take error of mean

	if( mean[bn][dn] > Cmin && mean[bn][dn] < Cmax   && fabs(res[bn][dn]) <  Cdiffmax ){ 
	  double rr=res[bn][dn]/rad_tib - (cos(phib[bn])*ps[1]/rad_tib - sin(phib[bn])*ps[0]/rad_tib + ps[4]) + ps[2]*izp_tib[dn]/rad_tib*cos(phib[bn]) + ps[3]*izp_tib[dn]/rad_tib*sin(phib[bn]);	  
	  chisq += rr*rr/er;	  
	  ndf++;
	}//if
      }//dn
    }//bn
    ndf-=npar;
    if(fabs(chisq)>1000.) chisq=10000.;
    vch2_tb.push_back(chisq);
    vndf_tb.push_back((ndf));

	if(bhist) {           
	  if(ndf==0) ndf=-1;
	  double ch2n= chisq/ndf;
	  thehisto->h1_tibch2->Fill(ch2n);
	}
    // scale to microns ans microrads 
  for(int ip = 0; ip < npar; ++ip){
       ps[ip]*=1000.;
    }
    ps[2]*=1000.;
    ps[3]*=1000.;
    ps[4]*=1000.;

  if(debug>1) {
      cout<<irun <<" "<<vruns[irun]<<" analTIB  final  ps =  "; 
      printa(ps,npar);
    }  
 
 // fill resutls
 vpar_tb.push_back(ps);
 
    if(debug>1) {
      cout<<irun <<" "<<vruns[irun]<<" analTIB  final  ps =  "; 
      printa(ps,npar);
    }  

  }//irun

  if(debug>1) {
    cout<<" ====all TIB params "<<vpar_tb.size()<<" "<<sfiles.size()<<endl;

    for(unsigned int irun=0;irun<nruns;irun++) {
      cout<<sfiles[irun]<<" "<<vruns[irun]<<" time="<<vtime[irun]<<endl;
	cout<<" dx= "<<vpar_tb.at(irun)[0]<<" dy= "<<vpar_tb.at(irun)[1]<<" Rx= "<<vpar_tb.at(irun)[2]<<" Ry= "<<vpar_tb.at(irun)[3]<<" Rz= "<<vpar_tb.at(irun)[4]<<endl;
      }

  }// debug

  return 0;
}


int LAS_stabil::alignTOB() 
{
  if(debug) cout<<"    alignTOB>>>"<<endl;
  //tob config
int nbmax=8;
int ndmax=6;
int npar=5;//dx,dy,rx,ry,rz

int  nav[nbmax][ndmax];
double mean0[nbmax][ndmax];
double rms0[nbmax][ndmax];
double mean[nbmax][ndmax];
double rms[nbmax][ndmax];


if(refrunn<0) {
 for(int i=0;i<nbmax;i++)  
   for(int j=0;j<ndmax;j++) {
    mean0[i][j]=0; 
    nav[i][j]=0; 
       }
 }

 LASGlobalDataLoop* loop=0;
// get reference 
for(int irun = 0; irun < nruns;  ++irun){
   
    TFile f_res(sfiles[irun].c_str());
if(debug>2)  cout<<" read file ="<<sfiles[irun]<<endl;
 
     LASGlobalData<double> *mean_ptr=0;   
    f_res.GetObject("positions_strips", mean_ptr);   
    LASGlobalData<double> *rms_ptr=0;   
    f_res.GetObject("rms_av_0", rms_ptr);   
 
loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TOB);
    do{
      int bn = loop->get_beam();
      int dn = loop->get_zpos(); 
  if(refrunn<0) {
  double m= loop->GetEntry<double>( *mean_ptr );
  double r=loop->GetEntry<double>( *rms_ptr );
    	if(m > Cmin && m < Cmax   && r <Crms){  
mean0[bn][dn] =mean0[bn][dn]+m;
nav[bn][dn]++;
	}
  }// refrun<0
 if(vruns[irun]==refrunn&&refrunn>0) {
   if(debug>10) cout<<" TIB found ref "<<vruns[irun]<<endl;
       // reference 
    mean0[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms0[bn][dn] = loop->GetEntry<double>( *rms_ptr );
     }
    }while(loop->next());   
  } // irun

 if(refrunn<0) {
 for(int bn=0; bn<nbmax;bn++) 
    for(int dn=0; dn<ndmax;dn++) 
       if(nav[bn][dn]!=0) mean0[bn][dn]=mean0[bn][dn]/nav[bn][dn];     
}


double bs[npar];
double res[nbmax][ndmax];
TMatrixT<double> h(npar,npar);
// clean vpartib
clearv(&vpar_tb);


vch2_tb.clear();
vndf_tb.clear();


// runs loop  
  for(int irun = 0; irun < nruns;  ++irun){ 

TFile f_res(sfiles[irun].c_str());
if(debug>2)  cout<<" analyze file TIB ="<<sfiles[irun]<<endl;
// clear each run
 cleara(bs,npar);
 clearh(&h);


LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);   
LASGlobalData<int> *norm=0;   
f_res.GetObject("norm_0", norm); 
 


loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TOB);
    do{
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
     // boundaries
     if(bn>nbmax||bn<0) continue;
     if(dn>ndmax||dn<0) continue;

    mean[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms[bn][dn] = loop->GetEntry<double>( *rms_ptr );
    int  anorm = loop->GetEntry<int>( *norm );
 

    //hist
    if(bhist) {
 if(debug>11) cout<<" fill histoTIB "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[bn][dn]<<" rms="<<rms[bn][dn]<<endl;
       thehisto->h2_tobmean[bn]->Fill(dn,mean[bn][dn]);
       thehisto->h2_tobrms[bn]->Fill(dn,rms[bn][dn]);
       thehisto->h2_tobnorm[bn]->Fill(dn,anorm);
    }

 if( rms[bn][dn]==0)
      rms[bn][dn]=mcsig_mean;

    if(anorm==0) anorm=1;
       if(anorm<10) 
           rms[bn][dn]=maxrms;

     rms[bn][dn]/=sqrt(anorm);


    if(debug>11) cout<<"  after  norm "<<bn<<dn<<" mean="<<mean[bn][dn]<<" "<<rms[bn][dn]<<endl;
      // new beam(assuming that the disks are stored inside the beams
  
//residuals loop over dn
	
	if(rms[bn][dn]==0.) {
	  if(debug>10) cout<<" ***rms are 0  put them to 1 "<<rms[bn][dn]<<endl;
          rms[bn][dn]=mcsig_mean;
	}
        if(anorm<10) 
        rms[bn][dn]=maxrms;

	// smear if MC
	if(bmc_mean) mcsmear(&mean[bn][dn],rms[bn][dn]);
	double er=  rms[bn][dn]* rms[bn][dn];

	double diff = (mean[bn][dn]-mean0[bn][dn])*pch_tob;
	if(bcor_at)
  	res[bn][dn]= diff - varo.at(irun)[bn]*izp_tob[dn] - vbro.at(irun)[bn];
        else 
	  res[bn][dn]= diff;


	if(debug>15) cout<<"bn="<<bn<<"  dn="<<dn<<" diff="<<diff<<" res="<<res[bn][dn]<<endl;
  
	if( mean[bn][dn] > Cmin && mean[bn][dn] < Cmax   && fabs(res[bn][dn]) <  Cdiffmax ){ 

  if(bhist) {
    thehisto->h1_tobstat[bn]->Fill(dn);
    thehisto->h2_tobstat->Fill(dn,bn);
    thehisto->h2_tobres[bn]->Fill(dn,diff);
  }
	  double dphi = res[bn][dn]/rad_tib;
	 	  
	  bs[0] = bs[0] + dphi*sin(phib[bn])*rad_tob/er;
	  bs[1] = bs[1] + dphi*cos(phib[bn])*rad_tob/er;
	  bs[2] = bs[2] + dphi*cos(phib[bn])*rad_tob*izp_tob[dn]/er;
	  bs[3] = bs[3] + dphi*sin(phib[bn])*rad_tob*izp_tob[dn]/er;
	  bs[4] = bs[4] + dphi/er;
// fill matrix
	  if(debug>10) cout<<bn<<dn<<" fill h "<<er<<endl;
	  if(debug>20) printh(h);

	  h(0,0) = h(0,0) - sin(phib[bn])*sin(phib[bn])/er;
	  h(0,1) = h(0,1) + sin(phib[bn])*cos(phib[bn])/er;
	  h(0,2) = h(0,2) + sin(phib[bn])*cos(phib[bn])*izp_tob[dn]/er;
	  h(0,3) = h(0,3) + sin(phib[bn])*sin(phib[bn])*izp_tob[dn]/er;
	  h(0,4) = h(0,4) + sin(phib[bn])*rad_tob/er;
	  
	  h(1,0) = h(1,0) - sin(phib[bn])*cos(phib[bn])/er;
	  h(1,1) = h(1,1) + cos(phib[bn])*cos(phib[bn])/er;
	  h(1,2) = h(1,2) + cos(phib[bn])*cos(phib[bn])*izp_tob[dn]/er;
	  h(1,3) = h(1,3) + sin(phib[bn])*cos(phib[bn])*izp_tob[dn]/er;
	  h(1,4) = h(1,4) + cos(phib[bn])*rad_tob/er;
	  
	  h(2,0) = h(2,0) - sin(phib[bn])*cos(phib[bn])*izp_tob[dn]/er;
	  h(2,1) = h(2,1) + cos(phib[bn])*cos(phib[bn])*izp_tob[dn]/er;
	  h(2,2) = h(2,2) + cos(phib[bn])*cos(phib[bn])*izp_tob[dn]*izp_tob[dn]/er;
	  h(2,3) = h(2,3) + sin(phib[bn])*cos(phib[bn])*izp_tob[dn]*izp_tob[dn]/er;
	  h(2,4) = h(2,4) + cos(phib[bn])*izp_tob[dn]*rad_tob/er;
	  
	  h(3,0) = h(3,0) - sin(phib[bn])*sin(phib[bn])*izp_tob[dn]/er;
	  h(3,1) = h(3,1) + sin(phib[bn])*cos(phib[bn])*izp_tob[dn]/er;
	  h(3,2) = h(3,2) + sin(phib[bn])*cos(phib[bn])*izp_tob[dn]*izp_tob[dn]/er;
	  h(3,3) = h(3,3) + sin(phib[bn])*sin(phib[bn])*izp_tob[dn]*izp_tob[dn]/er;
	  h(3,4) = h(3,4) + sin(phib[bn])*izp_tob[dn]*rad_tob/er;
	  
	  h(4,0) = h(4,0) - sin(phib[bn])/rad_tob/er;
	  h(4,1) = h(4,1) + cos(phib[bn])/rad_tob/er;
	  h(4,2) = h(4,2) + cos(phib[bn])*izp_tob[dn]/rad_tob/er;
	  h(4,3) = h(4,3) + sin(phib[bn])*izp_tob[dn]/rad_tob/er;
	  h(4,4) = h(4,4) + 1./er;
	}// ifselect
}while(loop->next()); 
    delete loop;

    if(debug>5) {
      cout<<irun<<" analTOB bs= ";   
      printa(bs,npar);
    }      

    if(debug>30) {
     cout << " TIB h before invertion----  "  << endl; 
     printh(h); 
    }//debug
    
// find solution
    h.Invert();

  if(debug>30) {
 cout << " TIB h after invertion----  "  << endl; 
   printh(h); 
  }//debug

 
  // fill results
 double* ps=new double[npar];
 
    for(int ip = 0; ip < npar; ++ip){
      ps[ip] = 0.;
      for(int ix = 0; ix < npar; ++ix){
	ps[ip] = ps[ip] + h(ip,ix)*bs[ix];
      }//ix
    }//ip


// calculate chi2
 int  ndf = 0;
  double  chisq = 0.;
    for(int bn = 0; bn < nbmax;  ++bn){
      for(int dn = 0; dn < ndmax;  ++dn){
	if(rms[bn][dn]==0.) {// should not occur
	  if(debug>10) cout<<" ***rms are 0  put them to 1 "<<rms[bn][dn]<<endl;
	  rms[bn][dn]=mcsig_mean;
	}
	double er=rms[bn][dn]*rms[bn][dn]/rad_tob/rad_tob; // take error of mean

	if( mean[bn][dn] > Cmin && mean[bn][dn] < Cmax   && fabs(res[bn][dn]) <  Cdiffmax ){ 
	  double rr=res[bn][dn]/rad_tob - (cos(phib[bn])*ps[1]/rad_tob - sin(phib[bn])*ps[0]/rad_tob + ps[4]) + ps[2]*izp_tob[dn]/rad_tob*cos(phib[bn]) + ps[3]*izp_tob[dn]/rad_tob*sin(phib[bn]);	  
	  chisq += rr*rr/er;	  
	  ndf++;
	}//if
      }//dn
    }//bn
    ndf-=npar;
    if(fabs(chisq)>1000.) chisq=10000.;
    vch2_tb.push_back(chisq);
    vndf_tb.push_back((ndf));

	if(bhist) {           
	  if(ndf==0) ndf=-1;
	  double ch2n= chisq/ndf;
	  thehisto->h1_tobch2->Fill(ch2n);
	}
    // scale to microns ans microrads 
  for(int ip = 0; ip < npar; ++ip){
       ps[ip]*=1000.;
    }
    ps[2]*=1000.;
    ps[3]*=1000.;
    ps[4]*=1000.;

  if(debug>1) {
      cout<<irun <<" "<<vruns[irun]<<" analTIB  final  ps =  "; 
      printa(ps,npar);
    }  
 
 // fill resutls
 vpar_tb.push_back(ps);
 
    if(debug>1) {
      cout<<irun <<" "<<vruns[irun]<<" analTIB  final  ps =  "; 
      printa(ps,npar);
    }  

  }//irun

  if(debug>1) {
    cout<<" ====all TOB params "<<vpar_tb.size()<<" "<<sfiles.size()<<endl;

    for(unsigned int irun=0;irun<nruns;irun++) {
      cout<<sfiles[irun]<<" "<<vruns[irun]<<" time="<<vtime[irun]<<endl;
	cout<<" dx= "<<vpar_tb.at(irun)[0]<<" dy= "<<vpar_tb.at(irun)[1]<<" Rx= "<<vpar_tb.at(irun)[2]<<" Ry= "<<vpar_tb.at(irun)[3]<<" Rz= "<<vpar_tb.at(irun)[4]<<endl;
      }

  }// debug

  return 0;
}



// algn TEC's +1, -1
int LAS_stabil::alignTEC(int sgn) {
  if(debug) cout<<"    alignTEC>>>"<<sgn<<endl;
  if(sgn==0) return 10;
  sgn/=abs(sgn);
  //tec config

int nbmax=8;
int ndmax=4;
int npar=3;//dx,dy,rz

int  nav[nbmax][ndmax];
double mean0[nbmax][ndmax];
double rms0[nbmax][ndmax];
double mean[nbmax][ndmax];
double rms[nbmax][ndmax];

 if(refrunn<0){
 for(int i=0;i<nbmax;i++)  
   for(int j=0;j<ndmax;j++) {
    mean0[i][j]=0; 
    nav[i][j]=0; 
       }
 }



 LASGlobalDataLoop* loop=0;
// get reference 
for(int irun = 0; irun < nruns;  ++irun){
   
    TFile f_res(sfiles[irun].c_str());
if(debug>2)  cout<<" read file ="<<sfiles[irun]<<endl;
 
     LASGlobalData<double> *mean_ptr=0;   
    f_res.GetObject("positions_strips", mean_ptr);   
    LASGlobalData<double> *rms_ptr=0;   
    f_res.GetObject("rms_av_0", rms_ptr);   

 if(sgn>0) 
  loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TEC_PLUS_AT);
 if(sgn<0) 
   loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS_AT);

   do{
      int bn = loop->get_beam();
      int dn = loop->get_zpos();
  if(refrunn<0) {
   double m= loop->GetEntry<double>( *mean_ptr );
    double r=loop->GetEntry<double>( *rms_ptr );
    	if(m > Cmin && m < Cmax   && r <Crms){  
mean0[bn][dn] =mean0[bn][dn]+m;
nav[bn][dn]++;
	}
  }

  if(vruns[irun]==refrunn&&refrunn>0) {
       // reference 
    mean0[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms0[bn][dn] = loop->GetEntry<double>( *rms_ptr );
     }
    }while(loop->next());   
   
delete loop;
  } // irun


 if(refrunn<0) {
 for(int bn=0; bn<nbmax;bn++) 
    for(int dn=0; dn<ndmax;dn++) 
       if(nav[bn][dn]!=0) mean0[bn][dn]=mean0[bn][dn]/nav[bn][dn]; 
}

double bs[npar];
double res[nbmax][ndmax];
TMatrixT<double> h(npar,npar);
// clean
 if(sgn>0){
clearv(&vpar_tecp);
vch2_tecp.clear();
vndf_tecp.clear();
 }
 if(sgn<0){
clearv(&vpar_tecm);

 vch2_tecm.clear();
 vndf_tecm.clear();
 }


// runs loop  
  for(int irun = 0; irun < nruns;  ++irun){ 

TFile f_res(sfiles[irun].c_str());
if(debug>2)  cout<<" analyze file TEC"<<sgn<<" ="<<sfiles[irun]<<endl;
// clear each run
 cleara(bs,npar);
 clearh(&h);


LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);   
LASGlobalData<int> *norm=0;   
f_res.GetObject("norm_0", norm);   



 if(sgn>0) 
  loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TEC_PLUS_AT);
 if(sgn<0) 
   loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS_AT);

    do{
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
     // boundaries
     if(bn>nbmax||bn<0) continue;
     if(dn>ndmax||dn<0) continue;

    mean[bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms[bn][dn] = loop->GetEntry<double>( *rms_ptr );
     int  anorm = loop->GetEntry<int>( *norm );
    if(debug>10) cout<<" read "<<bn<<dn<<" mean="<<mean[bn][dn]<<" "<<rms[bn][dn]<<endl;
  

  //hist
    if(bhist) {
      if(sgn>0) {
if(debug>11) cout<<" fill histoTECP "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[bn][dn]<<" rms="<<rms[bn][dn]<<endl;
	    thehisto->h2_tecpmean[bn]->Fill(double(dn),mean[bn][dn]);
	    thehisto->h2_tecprms[bn]->Fill(double(dn),rms[bn][dn]);
	    thehisto->h2_tecpnorm[bn]->Fill(double(dn),anorm);
      }
 if(sgn<0) {
if(debug>11) cout<<" fill histoTECM "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[bn][dn]<<" rms="<<rms[bn][dn]<<endl;
        thehisto->h2_tecmmean[bn]->Fill(double(dn),mean[bn][dn]);
        thehisto->h2_tecmrms[bn]->Fill(double(dn),rms[bn][dn]);
        thehisto->h2_tecmnorm[bn]->Fill(double(dn),anorm);
      }
      }// hist


      // new beam(assuming that the disks are stored inside the beams
  
//residuals loop over dn
	if(rms[bn][dn]==0.) {
	  if(debug>10) cout<<" ***rms are 0  put them to 1 "<<rms[bn][dn]<<endl;
          rms[bn][dn]=mcsig_mean;
	}
      if(anorm==0) anorm=1;
       if(anorm<10) 
           rms[bn][dn]=maxrms;

     rms[bn][dn]/=sqrt(anorm);


	// smear if MC
	if(bmc_mean) mcsmear(&mean[bn][dn],rms[bn][dn]);
	double er=  rms[bn][dn]* rms[bn][dn];

	double diff = sgn*(mean[bn][dn]-mean0[bn][dn])*pch_tec;     
        double iz=-1;
	if(sgn>0) iz=izp_tecp[dn];
	if(sgn<0) iz=izp_tecm[dn];

        if(bcor_at)
	res[bn][dn]= diff - varo.at(irun)[bn]*iz - vbro.at(irun)[bn];
	else
	  res[bn][dn]= diff;

	if(debug>8) cout<<" align : TECbn="<<bn<<"  dn="<<dn<<" diff="<<diff<<" res="<<res[bn][dn]<<endl;
  
	if( mean[bn][dn] > Cmin && mean[bn][dn] < Cmax   && fabs(res[bn][dn]) <  Cdiffmax ){ 
 
if(bhist) {
  if(sgn>0){
     thehisto->h1_tecpstat[bn]->Fill(dn);
     thehisto->h2_tecpstat->Fill(dn,bn);
     thehisto->h2_tecpres[bn]->Fill(dn,diff);
  } 
  if(sgn<0) {
    thehisto->h1_tecmstat[bn]->Fill(dn);
    thehisto->h2_tecmstat->Fill(dn,bn);
    thehisto->h2_tecmres[bn]->Fill(dn,diff);
  } 
 }


	  double dphi = res[bn][dn];
	  bs[0] = bs[0] + dphi*sin(phib[bn])/er;
	  bs[1] = bs[1] + dphi*cos(phib[bn])/er;
	  bs[2] = bs[2] + dphi/er;
// fill matrix
	  if(debug>10) cout<<bn<<dn<<" fill h "<<er<<endl;
	  if(debug>20) printh(h);

	  h(0,0) = h(0,0) - sin(phib[bn])*sin(phib[bn])/er;
	  h(0,1) = h(0,1) + sin(phib[bn])*cos(phib[bn])/er;
	  h(0,2) = h(0,2) + sin(phib[bn])/er;

	  h(1,0) = h(1,0) - sin(phib[bn])*cos(phib[bn])/er;
	  h(1,1) = h(1,1) + cos(phib[bn])*cos(phib[bn])/er;
	  h(1,2) = h(1,2) + cos(phib[bn])/er;

	  h(2,0) = h(2,0) - sin(phib[bn])/er;
	  h(2,1) = h(2,1) + cos(phib[bn])/er;
	  h(2,2) = h(2,2) + 1./er;	  

	}// ifselect
}while(loop->next()); 


    if(debug>5) {
      cout<<irun<<" analTEC"<<sgn<<" bs= ";   
      printa(bs,npar);
    }      

    if(debug>30) {
     cout << " TEC"<<sgn<<" h before invertion----  "  << endl; 
     printh(h); 
    }//debug
    
// find solution
    h.Invert();

  if(debug>30) {
 cout << "  TEC"<<sgn<<"h after invertion----  "  << endl; 
   printh(h); 
  }//debug

  // fill results
double* ps=new double[npar];
 

    for(int ip = 0; ip < npar; ++ip){
      ps[ip] = 0.;
      for(int ix = 0; ix < npar; ++ix){  
	ps[ip] = ps[ip] + h(ip,ix)*bs[ix];
       
	//	cout<<ip<<" TEC "<<ps[ip]<<" file "<<sfiles[irun]<<endl;
      }//ix
    }//ip

   ps[2] = ps[2]/rad_tec;
// calculate chi2
 int  ndf = 0;
  double  chisq = 0.;
    for(int bn = 0; bn < nbmax;  ++bn){
      for(int dn = 0; dn < ndmax;  ++dn){
	if(rms[bn][dn]==0.) {// should not occur
	  if(debug>10) cout<<" ***rms are 0  put them to 1 "<<rms[bn][dn]<<endl;
	  rms[bn][dn]=mcsig_mean;
	}
		double er=rms[bn][dn]*rms[bn][dn];

	if( mean[bn][dn] > Cmin && mean[bn][dn] < Cmax   && fabs(res[bn][dn]) <  Cdiffmax ){ 

	  double rr=res[bn][dn] - (cos(phib[bn])*(ps[1]) - sin(phib[bn])*(ps[0]) + ps[2]*rad_tec);
	  chisq += rr*rr/er;	  
	  ndf++;
	}//if
      }//dn
    }//bn
    ndf-=npar;
    if(fabs(chisq)>1000.) chisq=10000.;

      if(sgn>0) {
    vch2_tecp.push_back(chisq);
    vndf_tecp.push_back((ndf));
	if(bhist) {
	  if(ndf==0) ndf=-1;
	  double ch2n= chisq/ndf;
	  thehisto->h1_tecpch2->Fill(ch2n);
	}
      }
   if(sgn<0) {
    vch2_tecm.push_back(chisq);
    vndf_tecm.push_back((ndf));
	if(bhist) {
	  if(ndf==0) ndf=-1;
	  double ch2n= chisq/ndf;
	  thehisto->h1_tecmch2->Fill(ch2n);
	}
      }

for(int ip = 0; ip < npar; ++ip){
  if(fabs(ps[ip])<1.e+9)
       ps[ip]*=1000.;
      else 
	ps[ip]=0.;
    }
 ps[2] = ps[2]*1000;

    if(sgn>0) {
   vpar_tecp.push_back(ps);
    }
    if(sgn<0) {
   vpar_tecm.push_back(ps);
    }

    if(debug>1) {
      cout<<irun <<" "<<vruns[irun]<<" analTEC"<<sgn<<"  final  ps =  "; 
      printa(ps,npar);
    }  


  }//irun

  //for(unsigned int irun=0;irun<nruns;irun++) 
  //  cout<<"TEC  file"<<sfiles[irun]<<" dx= "<<vpar_tecp.at(irun)[0]<<" dy= "<<vpar_tecp.at(irun)[1]<<" Rz= "<<vpar_tecp.at(irun)[2]<<endl;

  if(debug>1) {
    cout<<" ====all TIB params "<<" "<<sfiles.size()<<endl;

    for(unsigned int irun=0;irun<nruns;irun++) {
      cout<<sfiles[irun]<<" "<<vruns[irun]<<" time="<<vtime[irun]<<endl;
  if(sgn>0)
      cout<<" dx= "<<vpar_tecp.at(irun)[0]<<" dy= "<<vpar_tecp.at(irun)[1]<<" Rz= "<<vpar_tecp.at(irun)[2]<<endl;
  if(sgn<0)
      cout<<" dx= "<<vpar_tecm.at(irun)[0]<<" dy= "<<vpar_tecm.at(irun)[1]<<" Rz= "<<vpar_tecm.at(irun)[2]<<endl;
      }

  }// debug

  return 0;
}

// internal alignement of TEC+,-
int LAS_stabil::alignTECint(int sgn) {
  if(debug) cout<<"alignintTEC>>>"<<sgn<<endl;
  if(sgn==0) return 10;
  sgn/=abs(sgn);
  //

int nbmax=8;
int ndmax=9;
int nrmax=2;
int nparmax=27;//dx,dy,rz for 9 rings
int npar=3;//dx,dy,rz 
int  nav[nrmax][nbmax][ndmax];
double mean0[nrmax][nbmax][ndmax];
double rms0[nrmax][nbmax][ndmax];
double mean[nrmax][nbmax][ndmax];
double rms[nrmax][nbmax][ndmax];

 if(refrunn<0){
   for(int k=0;k<nrmax;k++) {
   for(int i=0;i<nbmax;i++) {
   for(int j=0;j<ndmax;j++) {
    mean0[k][i][j]=0; 
    nav[k][i][j]=0; 
       }
    }
   }
 } //if

 LASGlobalDataLoop* loop=0;
// get reference 
for(int irun = 0; irun < nruns;  ++irun){
   
    TFile f_res(sfiles[irun].c_str());
if(debug>2)  cout<<" read file ="<<sfiles[irun]<<endl;
 
     LASGlobalData<double> *mean_ptr=0;   
    f_res.GetObject("positions_strips", mean_ptr);   
    LASGlobalData<double> *rms_ptr=0;   
    f_res.GetObject("rms_av_0", rms_ptr);   

 if(sgn>0) 
  loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TEC_PLUS);
 if(sgn<0) 
   loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS);

   do{
      int rn = loop->get_ring();
      int bn = loop->get_beam();
      int dn = loop->get_zpos();
     if(rn>nrmax||rn<0) continue;
     if(bn>nbmax||bn<0) continue;
     if(dn>ndmax||dn<0) continue;


  if(refrunn<0) {
   double m= loop->GetEntry<double>( *mean_ptr );
   double r=loop->GetEntry<double>( *rms_ptr );
    	if(m > Cmin && m < Cmax   && r <Crms){  
mean0[rn][bn][dn] =mean0[rn][bn][dn]+m;
nav[rn][bn][dn]++;
	}
  }

  if(vruns[irun]==refrunn&&refrunn>0) {
       // reference 
    mean0[rn][bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms0[rn][bn][dn] = loop->GetEntry<double>( *rms_ptr );
     }
    }while(loop->next());   
   
delete loop;
  } // irun



 if(refrunn<0) {
for(int rn=0; rn<nrmax;rn++)
  for(int bn=0; bn<nbmax;bn++) 
     for(int dn=0; dn<ndmax;dn++) 
        if(nav[rn][bn][dn]!=0) mean0[rn][bn][dn]=mean0[rn][bn][dn]/nav[rn][bn][dn]; 
}




double bs[npar];
double res[nrmax][nbmax][ndmax];
TMatrixT<double> h(npar,npar);
// clean
 if(sgn>0){
clearv(&vpar_tecpint);
vch2_tecpint.clear();
vndf_tecpint.clear();
 }
 if(sgn<0){
clearv(&vpar_tecmint);
 vch2_tecmint.clear();
 vndf_tecmint.clear();
 }



// runs loop  
  for(int irun = 0; irun < nruns;  ++irun){ 

TFile f_res(sfiles[irun].c_str());
if(debug>2)  cout<<" analyze file TECint"<<sgn<<" ="<<sfiles[irun]<<endl;
// clear each run
 cleara(bs,npar);
 clearh(&h);


LASGlobalData<double> *mean_ptr=0;      
f_res.GetObject("positions_strips", mean_ptr);   
LASGlobalData<double> *rms_ptr=0;   
f_res.GetObject("rms_av_0", rms_ptr);   
LASGlobalData<int> *norm=0;   
f_res.GetObject("norm_0", norm);   



 if(sgn>0) 
  loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TEC_PLUS);
 if(sgn<0) 
   loop=new  LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS);

    do{
     int rn = loop->get_ring();
     int  bn = loop->get_beam();
     int  dn = loop->get_zpos();
     // boundaries
    if(rn>nbmax||rn<0) continue;
     if(bn>nbmax||bn<0) continue;
     if(dn>ndmax||dn<0) continue;

    mean[rn][bn][dn] = loop->GetEntry<double>( *mean_ptr );
    rms[rn][bn][dn] = loop->GetEntry<double>( *rms_ptr );
     int  anorm = loop->GetEntry<int>( *norm );
    if(debug>10) cout<<" read "<<bn<<dn<<" mean="<<mean[rn][bn][dn]<<" "<<rms[rn][bn][dn]<<endl;
  
//hist
   if(bhist) {
      if(sgn>0) {
//if(debug>11) cout<<" fill histoTECP "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[rn][bn][dn]<<" rms="<<rms[rn][bn][dn]<<endl;
 	    thehisto->h2_tecpmeanint[rn][bn]->Fill(double(dn),mean[rn][bn][dn]);
// 	    thehisto->h2_tecprmsint[rn][bn]->Fill(double(dn),rms[rn][bn][dn]);
// 	    thehisto->h2_tecpnormint[rn][bn]->Fill(double(dn),anorm);
     }
     if(sgn<0) {
    //if(debug>11) cout<<" fill histoTECM "<<bn<<" dn= "<<double(dn)<<"  mean="<<mean[rn][bn][dn]<<" rms="<<rms[rn][bn][dn]<<endl;
         thehisto->h2_tecmmeanint[rn][bn]->Fill(double(dn),mean[rn][bn][dn]);
//         thehisto->h2_tecmrmsint[rn][bn]->Fill(double(dn),rms[rn][bn][dn]);
//         thehisto->h2_tecmnormint[rn][bn]->Fill(double(dn),anorm);
        }
         }// hist

 
	if(rms[rn][bn][dn]==0.) {
	  if(debug>10) cout<<" ***rms are 0  put them to 1 "<<rms[rn][bn][dn]<<endl;
          rms[rn][bn][dn]=mcsig_mean;
	}
      if(anorm==0) anorm=1;
       if(anorm<10) 
           rms[rn][bn][dn]=maxrms;
       // rms of mean
     rms[rn][bn][dn]/=sqrt(anorm);

    if(bmc_mean) mcsmear(&mean[rn][bn][dn],rms[rn][bn][dn]);

}while(loop->next()); 

    // emulate beam orientation errors 
    double ar[nrmax][nbmax];
    double br[nrmax][nbmax];
 for(int rn = 0; rn < nrmax;  ++rn){
         for(int bn = 0; bn < nbmax;  ++bn){
	   ar[rn][bn]=0.;
           br[rn][bn]=0.;
	 }
 }

 if(bmc_at) {
      if(debug>5) cout<<"smear TECint beam params "<<endl;
    
      //  cout<<"smear TECint beam params "<<sgn<<endl;
      for(int rn = 0; rn < nrmax;  ++rn){
         for(int bn = 0; bn < nbmax;  ++bn){

	mcsmear(&ar[rn][bn], mcsig_atec);
	mcsmear(&br[rn][bn], mcsig_btec);
	 }
      }
    }//bmc

  
 // calculate parameters in a sepaate loop
    //  dx0,dy0,rz0,dx1,dy1,rz1,....
double* ps = new double[ndmax*npar];

  for(int dn = 0; dn < ndmax;  ++dn){
// clear for each disk
 cleara(bs,npar);
 clearh(&h);
     for(int rn = 0; rn < nrmax;  ++rn){
         for(int bn = 0; bn < nbmax;  ++bn){
	double er=  rms[rn][bn][dn]* rms[rn][bn][dn];
	double diff = sgn*(mean[rn][bn][dn]-mean0[rn][bn][dn])*pch_tecint[rn];     
        double iz=-1;
	if(sgn>0) iz=izp_tecp[dn];
	if(sgn<0) iz=izp_tecm[dn];

	//  beam correction
	res[rn][bn][dn]= diff-ar[rn][bn]*iz-br[rn][bn];

	if(debug>8) cout<<" align : TECbn="<<bn<<"  dn="<<dn<<" diff="<<diff<<" res="<<res[rn][bn][dn]<<endl;

	//	cout<<" align : "<<sgn<<" TECbn= "<<bn<<"  dn="<<dn<<" diff="<<diff<<" res="<<res[rn][bn][dn]<<endl;

	if( mean[rn][bn][dn] > Cmin && mean[rn][bn][dn] < Cmax   && fabs(res[rn][bn][dn]) <  Cdiffmax ){ 

if(bhist) {
  if(sgn>0){
   //   thehisto->h1_tecpstatint[rn][bn]->Fill(dn);
//      thehisto->h2_tecpstatint[rn]->Fill(dn,bn);
//      thehisto->h2_tecpresint[rn][bn]->Fill(dn,diff);
  } 
  if(sgn<0) {
  //   thehisto->h1_tecmstatint[rn][bn]->Fill(dn);
//    thehisto->h2_tecmstatint[rn]->Fill(dn,bn);
//     thehisto->h2_tecmresint[rn][bn]->Fill(dn,diff);
  } 
 }	


	  double dphi = res[rn][bn][dn]/rad_tecint[rn];
            bs[0] = bs[0] + dphi*sin(phib[bn])*rad_tecint[rn]/er;
  	    bs[1] = bs[1] + dphi*cos(phib[bn])*rad_tecint[rn]/er;
  	    bs[2] = bs[2] + dphi/er/rad_tecint[rn];
 
 if(debug>10) cout<<bn<<dn<<" fill h "<<er<<endl;

	    h(0,0) = h(0,0) + sin(phib[bn])*sin(phib[bn])/er;
    	    h(0,1) = h(0,1) - cos(phib[bn])*sin(phib[bn])/er;
    	    h(0,2) = h(0,2) + sin(phib[bn])*rad_tecint[rn]/er;

    	    h(1,0) = h(1,0) + cos(phib[bn])*sin(phib[bn])/er;
    	    h(1,1) = h(1,1) - cos(phib[bn])*cos(phib[bn])/er;
    	    h(1,2) = h(1,2) + cos(phib[bn])*rad_tecint[rn]/er;

    	    h(2,0) = h(2,0) + sin(phib[bn])/rad_tecint[rn]/er;
    	    h(2,1) = h(2,1) - cos(phib[bn])/rad_tecint[rn]/er;
    	    h(2,2) = h(2,2) + 1./er;

      } //if res
    }//bn
  }// rn
     // solution for each disk, 3 params
  h.Invert();

  for(int ip = 0; ip < npar; ++ip){
      ps[dn*npar+ip] = 0.;
      for(int ix = 0; ix < npar; ++ix){
    	ps[dn*npar+ip] = ps[dn*npar+ip]+ h(ip,ix)*bs[ix];
      }
      ps[dn*npar+ip] = ps[dn*npar+ip]*1000;
    }
    ps[dn*npar+2] = ps[dn*npar+2]*1000;

    if(debug>5) cout<<" TECint  disk="<<dn<<" ps0= "<<ps[dn*npar+0]<<" ps1= "<<ps[dn*npar+1]<<"  ps2= "<<ps[dn*npar+2]<<endl;

  }//dn

 if(sgn>0) {
   vpar_tecpint.push_back(ps);
    }
    if(sgn<0) {
   vpar_tecmint.push_back(ps);
    }

    if(debug>1) {
      cout<<irun <<" "<<vruns[irun]<<" analintTEC"<<sgn<<"  final  ps =  "; 
      printa(ps,npar*ndmax);
    }  

  } //irun



 if(debug) cout<<">>>>alignintTEC done"<<sgn<<endl;
  return 0;
}

// extract TEC int from RDC
int LAS_stabil::alignTECintRDC() {

  return 0;
}
 
int LAS_stabil::alignAll() {
 if(debug) cout<<" alignALL>>>>"<<endl;
 
 if(btecint==1) {
if(alignTECint(+1)) cout<<"**algnTECint+ failed"<<endl; 
if(alignTECint(-1)) cout<<"**algnTECint- failed"<<endl; 
 }

if(btecint==2) {
  alignTECintRDC();
 }


 if(battob) {
 if(alignATtob())   cout<<"**algnATtob failed"<<endl; 
 if(alignTIB())   cout<<"**algnTIB failed"<<endl; 
 } else {
  
 if(alignATtib())   cout<<"**algnATtib failed"<<endl; 
 if(alignTOB())   cout<<"**algnT0B failed"<<endl; 
 }

 if(alignTEC(+1)) cout<<"**algnTEC+ failed"<<endl; 
 if(alignTEC(-1)) cout<<"**algnTEC+ failed"<<endl;

 if(sumch2())   cout<<"**sumch2  failed"<<endl;


if(debug) cout<<" >>>alignALL done"<<endl;
  return 0;
}

//get ch2
 //
int LAS_stabil::sumch2() {
if(debug) cout<<" sumch2 >>>"<<endl;
// calculate chi2
  vch2.clear();
if(int(vch2_at.size())!=nruns) return 1;
if(int(vch2_tb.size())!=nruns) return 1;
if(int(vch2_tecp.size())!=nruns) return 1;
if(int(vch2_tecm.size())!=nruns) return 1;

 int nbmax=8;
for(int irun=0; irun<nruns;irun++) {
  if(debug) cout<<" irun="<<irun<<" "<<vruns.at(irun)<<endl;
  // AT
  double Sch2=0;
  int Sndf=0;
  for(int k=0; k<nbmax;k++) {
    if(vndf_at.at(irun)[k]>0&&vch2_at.at(irun)[k]<1000.) {
    Sch2+=  vch2_at.at(irun)[k];
    Sndf+=  vndf_at.at(irun)[k];
    }
  }
  if(debug>15) cout<<" ch2at "<< Sch2<<" "<<Sndf<<endl;
  if(vndf_tb.at(irun)>0&&vch2_tb.at(irun)<1000. ) {
  Sch2+=vch2_tb.at(irun);
  Sndf+=vndf_tb.at(irun);
  }
 if(debug>15) cout<<" ch2at+tib "<< Sch2<<" "<<Sndf<<endl;
 if(vndf_tecp.at(irun)>0&&vch2_tecp.at(irun)<1000. ) {
  Sch2+=vch2_tecp.at(irun);
  Sndf+=vndf_tecp.at(irun);
 }
 if(debug>15) cout<<" ch2at+tib+tecp "<< Sch2<<" "<<Sndf<<endl;
 if(vndf_tecm.at(irun)>0&&vch2_tecm.at(irun)<1000.) {
  Sch2+=vch2_tecm.at(irun);
  Sndf+=vndf_tecm.at(irun);
 }
if(debug>15) cout<<" ch2at+tib+tecp+tecm "<< Sch2<<" "<<Sndf<<endl;
  //
 if(Sndf!=0&&Sndf <1.e+6) {
    Sch2=Sch2/Sndf;
    // Sch2=Sch2;
 }
  else Sch2=999.;

  vch2.push_back(Sch2);
  if(debug>10) cout<<" fill sumch2= "<< Sch2<<endl;


  if(debug>2) cout<<" run "<<vruns.at(irun)<<" ch2tot="<<Sch2<<endl;
 }//irun


if(debug) cout<<" >>>sumch2 done"<<endl;
  return 0;
}


void LAS_stabil::calculatesign() {
if(debug) cout<<" calcsign >>>"<<endl;
  // calculte  total significance Si*sigma1+S2*sigma2+.../sqrt(sigma^2+..) for all parameters
  // better sum/sqrt(n)
 vsign.clear();
for(int irun=0; irun<nruns;irun++) {
  double sign=0;
  int ntp=0;
 for(int k=0; k<5;k++) {
   if(vparer_tb.at(irun)[k]!=0) {
     sign+=vpar_tb.at(irun)[k]/vparer_tb.at(irun)[k];
     ntp++;
   } 
 }//k
for(int k=0; k<3;k++) {
  if(vparer_tecp.at(irun)[k]!=0) {
     sign+=vpar_tecp.at(irun)[k]/vparer_tecp.at(irun)[k];
     ntp++;
  } 
  if(vparer_tecm.at(irun)[k]!=0){
     sign+=vpar_tecm.at(irun)[k]/vparer_tecm.at(irun)[k];
     ntp++;   
  }
 } // k

 if(ntp!=0) sign/=sqrt(ntp);
 vsign.push_back(sign);
 }
 if(debug) cout<<" >>>calcsign "<<vsign.size()<<endl;
 return;
}

 // clear vector 
void LAS_stabil::clearv(vector<double*>* vv) {
   for(unsigned int  i=0; i<vv->size();i++) {
     delete [] vv->at(i);
   }
   vv->clear();
 }

// clear vector 
void LAS_stabil::clearv(vector<int*>* vv) {
   for(unsigned int  i=0; i<vv->size();i++) {
     delete [] vv->at(i);
   }
   vv->clear();
 }

//print array
void LAS_stabil::printa(double* a, int n) {
  cout<<" print array a("<<n<<")"<<endl;
  for(int i=0;i<n;i++) { 
      cout<<" a("<<i<<") =" << a[i]; 
  }
 cout<<endl;
  return;
}
//clear array
void LAS_stabil::cleara(double* a, int n) {
  for(int i=0;i<n;i++) { 
    a[i]=0; 
  }
  return;
}

//clear array
void LAS_stabil::cleara(int* a, int n) {
  for(int i=0;i<n;i++) { 
    a[i]=0; 
  }
  return;
}

//print h
void LAS_stabil::printh(TMatrixT<double>& h) {
  cout<<"print matrix  h("<<h.GetNcols()<<","<< h.GetNrows()<<")"<<endl;
  for(int i=0;i<h.GetNcols();i++) {
    for(int j=0;j<h.GetNrows();j++) {
      cout<<" h("<<i<<","<<j<<") = " << h(i,j); 
    }
    cout<<endl;
  }
  return;
}
//clear h
void LAS_stabil::clearh(TMatrixT<double>* h) {
  for(int i=0;i<h->GetNcols();i++) {
    for(int j=0;j<h->GetNrows();j++) {
      (*h)(i,j)=0.;
    }
  }
    return;
  }


//store in file
int LAS_stabil::printres() {
  if(debug) cout<<" ===printres "<<endl;
 
if(int(vpar_tb.size())!=nruns) return 1;
 
if(int(vpar_tecp.size())!=nruns) return 1;
 
if(int(vpar_tecm.size())!=nruns) return 1;
 
if(int(vch2_at.size())!=nruns) return 1;
 
if(int(vch2_tb.size())!=nruns) return 1;
 
if(int(vch2_tecp.size())!=nruns) return 1;
 
if(int(vch2_tecm.size())!=nruns) return 1;
 

   if(bmc==0){
cout<<"       time "
<<"  run   "
<<"   dx_tb "
<<"   dy_tb "
<<"   rx_tb "
<<"   ry_tb "
<<"   rz_tb "
<<"   dx_tecp "
<<"   dy_tecp "
<<"   rz_tecp "
<<"   dx_tecp "
<<"   dy_tecp "
<<"   rz_tecp "
<<"   ch2     "
<<endl;
 for(int irun=0; irun<nruns; irun++) {
    cout<< vtime.at(irun)
    <<" "<<vruns.at(irun)
    <<" "<<vpar_tb.at(irun)[0]
    <<" "<<vpar_tb.at(irun)[1]
    <<" "<<vpar_tb.at(irun)[2]
    <<" "<<vpar_tb.at(irun)[3]
    <<" "<<vpar_tb.at(irun)[4]
  
    <<" "<<vpar_tecp.at(irun)[0]
    <<" "<<vpar_tecp.at(irun)[1]
    <<" "<<vpar_tecp.at(irun)[2]
   
    <<" "<<vpar_tecm.at(irun)[0]
    <<" "<<vpar_tecm.at(irun)[1]
    <<" "<<vpar_tecm.at(irun)[2]
	<<" "<<vch2.at(irun)
<<endl;
 }
   } else {

cout<<"       time "
<<"  run   "
<<"  run   "
<<"   dx_tb "
<<"   dxer_tb "
<<"   dy_tb "
<<"   dyer_tb "
<<"   rx_tb "
<<"   rxer_tb "
<<"   ry_tb "
<<"   ryer_tb "
<<"   rz_tb "
<<"   rzer_tb "
<<"   dx_tecp "
<<"   dxer_tecp "
<<"   dy_tecp "
<<"   dyer_tecp "
<<"   rz_tecp "
<<"   rzer_tecp "
<<"   dx_tecp "
<<"   dxer_tecp "
<<"   dy_tecp "
<<"   dyer_tecp "
<<"   rz_tecp "
<<"   rzer_tecp "
<<"   ch2     "
<<endl;
 for(int irun=0; irun<nruns; irun++) {

    cout<< vtime.at(irun)
    <<" "<<vruns.at(irun)
    <<" "<<vpar_tb.at(irun)[0]
    <<" "<<vparer_tb.at(irun)[0]
    <<" "<<vpar_tb.at(irun)[1]
    <<" "<<vparer_tb.at(irun)[1]
    <<" "<<vpar_tb.at(irun)[2]
    <<" "<<vparer_tb.at(irun)[2]
    <<" "<<vpar_tb.at(irun)[3]
    <<" "<<vparer_tb.at(irun)[3]
    <<" "<<vpar_tb.at(irun)[4]
    <<" "<<vparer_tb.at(irun)[4]
  
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
<<endl;
 }

   }// else

   return 0;
}
//offset all params to mean. keeping errors
void  LAS_stabil::offsetpar() {
 if(debug ) cout<<" >>>> offset params  "<<endl; 

 int ntp=17; // 
 double avm[ntp];
 int avi[ntp];
 for(int i=0;i<ntp;i++) {
   avm[i]=0;
   avi[i]=0;
 }


 for(int irun=0; irun<nruns; irun++) {
   
   if(fabs(vpar_tb.at(irun)[0])<100)  { avm[0]+= vpar_tb.at(irun)[0]; avi[0]++;} 
   if(fabs(vpar_tb.at(irun)[1])<100)  {  avm[1]+= vpar_tb.at(irun)[1];avi[1]++;} 
   if(fabs(vpar_tb.at(irun)[2])<100)  {  avm[2]+= vpar_tb.at(irun)[2];avi[2]++;} 
   if(fabs(vpar_tb.at(irun)[3])<100)  {  avm[3]+= vpar_tb.at(irun)[3];avi[3]++;} 
   if(fabs(vpar_tb.at(irun)[4])<100)  { avm[4]+= vpar_tb.at(irun)[4];avi[4]++;} 

   if(fabs(vpar_tecp.at(irun)[0])<100)   { avm[5]+= vpar_tecp.at(irun)[0];avi[5]++;} 
   if(fabs(vpar_tecp.at(irun)[1])<100)   { avm[6]+= vpar_tecp.at(irun)[1];avi[6]++;} 
   if(fabs(vpar_tecp.at(irun)[2])<100)   { avm[7]+= vpar_tecp.at(irun)[2];avi[7]++;} 
 
   if(fabs(vpar_tecm.at(irun)[0])<100)   { avm[8]+= vpar_tecm.at(irun)[0];avi[8]++;} 
   if(fabs(vpar_tecm.at(irun)[1])<100)   { avm[9]+= vpar_tecm.at(irun)[1];avi[9]++;} 
   if(fabs(vpar_tecm.at(irun)[2])<100)   { avm[10]+= vpar_tecm.at(irun)[2];avi[10]++;} 


   if(fabs(vpar_tecpint.at(irun)[0])<100)   { avm[11]+= vpar_tecpint.at(irun)[0];avi[11]++;} 
   if(fabs(vpar_tecpint.at(irun)[1])<100)   { avm[12]+= vpar_tecpint.at(irun)[1];avi[12]++;} 
   if(fabs(vpar_tecpint.at(irun)[2])<100)   { avm[13]+= vpar_tecpint.at(irun)[2];avi[13]++;} 
 
   if(fabs(vpar_tecmint.at(irun)[0])<100)   { avm[14]+= vpar_tecmint.at(irun)[0];avi[14]++;} 
   if(fabs(vpar_tecmint.at(irun)[1])<100)   { avm[15]+= vpar_tecmint.at(irun)[1];avi[15]++;} 
   if(fabs(vpar_tecmint.at(irun)[2])<100)   { avm[16]+= vpar_tecmint.at(irun)[2];avi[16]++;} 
 }// irun

 //
 for(int i=0;i<ntp;i++) {
   if(avi[i]!=0)
   avm[i]= avm[i]/avi[i];
   else
    avm[i]=0;
 }



 // offset
 for(int irun=0; irun<nruns; irun++) {
   
   vpar_tb.at(irun)[0]=vpar_tb.at(irun)[0]+avm[0];
   vpar_tb.at(irun)[1]=vpar_tb.at(irun)[1]+avm[1];
   vpar_tb.at(irun)[2]=vpar_tb.at(irun)[2]+avm[2];
   vpar_tb.at(irun)[3]=vpar_tb.at(irun)[3]+avm[3];
   vpar_tb.at(irun)[4]=vpar_tb.at(irun)[4]+avm[4];



 vpar_tecp.at(irun)[0]=vpar_tecp.at(irun)[0]+avm[5];
 vpar_tecp.at(irun)[1]=vpar_tecp.at(irun)[1]+avm[6];
 vpar_tecp.at(irun)[2]=vpar_tecp.at(irun)[2]+avm[7];
 vpar_tecm.at(irun)[0]=vpar_tecm.at(irun)[0]+avm[8];
 vpar_tecm.at(irun)[1]=vpar_tecm.at(irun)[1]+avm[9];
 vpar_tecm.at(irun)[2]=vpar_tecm.at(irun)[2]+avm[10];


  vpar_tecpint.at(irun)[0]=vpar_tecpint.at(irun)[0]+avm[11];
  vpar_tecpint.at(irun)[1]=vpar_tecpint.at(irun)[1]+avm[12];
  vpar_tecpint.at(irun)[2]=vpar_tecpint.at(irun)[2]+avm[13];
  vpar_tecmint.at(irun)[0]=vpar_tecmint.at(irun)[0]+avm[14];
  vpar_tecmint.at(irun)[1]=vpar_tecmint.at(irun)[1]+avm[15];
  vpar_tecmint.at(irun)[2]=vpar_tecmint.at(irun)[2]+avm[16];
 }// irun


 return;
}

 //store in file
int  LAS_stabil::writeres(const string& fout) {
  if(debug ) cout<<" >>>> write to file "<<fout<<" "<<refrunn<<endl; 
 
if(int(vpar_tb.size())!=nruns) return 1;
if(int(vpar_tecp.size())!=nruns) return 1;
if(int(vpar_tecm.size())!=nruns) return 1;

if(int(vch2_at.size())!=nruns) return 1;
if(int(vch2_tb.size())!=nruns) return 1;
if(int(vch2_tecp.size())!=nruns) return 1;
if(int(vch2_tecm.size())!=nruns) return 1;


 if(debug>5) cout<<" sizes "<<vpar_tb.size()<<" "<<vpar_tecp.size()<<" "<<vpar_tecm.size()
		 <<" "<<vpar_tecpint.size()<<" "<<vpar_tecmint.size()<<" "<<vch2.size()<<" "<<vsign.size()<<" "<<vpar_tecpint.size()<<" "<<vpar_tecmint.size()
                 <<" "<<vparer_tb.size()<<" "<<vparer_tecpint.size()<<endl;

  // write by hand
fstream file(fout.c_str(), std::ios::out);
 if(!file) {  cout<<"****cant create fout "<<endl;return 1;}
  file << setprecision(6);
  // header
  bool bhead=false;

  if(bmc) {

 if(bhead&&!btecint) {
 file<<"       time "
<<"  run   "
<<"   dx_tb "
<<"   dxer_tb "
<<"   dy_tb "
<<"   dyer_tb "
<<"   rx_tb "
<<"   rxer_tb "
<<"   ry_tb "
<<"   ryer_tb "
<<"   rz_tb "
<<"   rzer_tb "
<<"   dx_tecp "
<<"   dxer_tecp "
<<"   dy_tecp "
<<"   dyer_tecp "
<<"   rz_tecp "
<<"   rzer_tecp "
<<"   dx_tecp "
<<"   dxer_tecp "
<<"   dy_tecp "
<<"   dyer_tecp "
<<"   rz_tecp "
<<"   rzer_tecp "
<<"   ch2        "
<<endl;
  }

 if(bhead&&btecint) {
 file<<"       time "
<<"  run   "
<<"   dx_tb "
<<"   dxer_tb "
<<"   dy_tb "
<<"   dyer_tb "
<<"   rx_tb "
<<"   rxer_tb "
<<"   ry_tb "
<<"   ryer_tb "
<<"   rz_tb "
<<"   rzer_tb "
<<"   dx_tecp "
<<"   dxer_tecp "
<<"   dy_tecp "
<<"   dyer_tecp "
<<"   rz_tecp "
<<"   rzer_tecp "
<<"   dx_tecp "
<<"   dxer_tecp "
<<"   dy_tecp "
<<"   dyer_tecp "
<<"   rz_tecp "
<<"   rzer_tecp "
<<"   ch2        "
<<"   sig        "
<<"   dx_tecpint "
<<"   dxer_tecpint "
<<"   dy_tecpint "
<<"   dyer_tecpint "
<<"   rz_tecpint "
<<"   rzer_tecpint "
<<"   dx_tecpint "
<<"   dxer_tecpint "
<<"   dy_tecpint "
<<"   dyer_tecpint "
<<"   rz_tecpint "
<<"   rzer_tecpint "


<<endl;


 }



 for(int irun=0; irun<nruns; irun++) {

// check for nan

 for(int i=0; i<5;i++) {
   if(vpar_tb.at(irun)[i]!=vpar_tb.at(irun)[i]) vpar_tb.at(irun)[i]=-999;
   if(vparer_tb.at(irun)[i]!=vparer_tb.at(irun)[i]) vparer_tb.at(irun)[i]=-999;
 }
for(int i=0; i<3;i++) {
   if(vpar_tecp.at(irun)[i]!=vpar_tecp.at(irun)[i]) vpar_tecp.at(irun)[i]=-999;
   if(vparer_tecp.at(irun)[i]!=vparer_tecp.at(irun)[i]) vpar_tecp.at(irun)[i]=-999;
   if(vpar_tecm.at(irun)[i]!=vpar_tecm.at(irun)[i]) vpar_tecm.at(irun)[i]=-999;
   if(vparer_tecm.at(irun)[i]!=vparer_tecm.at(irun)[i]) vpar_tecm.at(irun)[i]=-999;

  if(vpar_tecpint.at(irun)[i]!=vpar_tecpint.at(irun)[i]) vpar_tecpint.at(irun)[i]=-999;
   if(vparer_tecpint.at(irun)[i]!=vparer_tecpint.at(irun)[i]) vpar_tecpint.at(irun)[i]=-999;
   if(vpar_tecmint.at(irun)[i]!=vpar_tecmint.at(irun)[i]) vpar_tecmint.at(irun)[i]=-999;
   if(vparer_tecmint.at(irun)[i]!=vparer_tecmint.at(irun)[i]) vpar_tecmint.at(irun)[i]=-999;

 }




   // store only nonzero results
double sum=vpar_tb.at(irun)[0]
    +vpar_tb.at(irun)[1]
    +vpar_tb.at(irun)[2]
    +vpar_tb.at(irun)[3]
    +vpar_tb.at(irun)[4]
    +vpar_tecp.at(irun)[0]
    +vpar_tecp.at(irun)[1]
    +vpar_tecp.at(irun)[2]
    +vpar_tecm.at(irun)[0]
    +vpar_tecm.at(irun)[1]
  +vpar_tecm.at(irun)[2];

double sume=vparer_tb.at(irun)[0]
    +vparer_tb.at(irun)[1]
    +vparer_tb.at(irun)[2]
    +vparer_tb.at(irun)[3]
    +vparer_tb.at(irun)[4]
    +vparer_tecp.at(irun)[0]
    +vparer_tecp.at(irun)[1]
    +vparer_tecp.at(irun)[2]
    +vparer_tecm.at(irun)[0]
    +vparer_tecm.at(irun)[1]
  +vparer_tecm.at(irun)[2];


 if( (vruns.at(irun)!=refrunn&&sum==0)||(fabs(sum)>1000.)||fabs(sume)>1000.) continue;

 if(vch2.at(irun)<0.001) continue;

 

 if(!btecint) {
    file<< vtime.at(irun)
    <<" "<<vruns.at(irun)
    <<" "<<vpar_tb.at(irun)[0]
    <<" "<<vparer_tb.at(irun)[0]
    <<" "<<vpar_tb.at(irun)[1]
    <<" "<<vparer_tb.at(irun)[1]
    <<" "<<vpar_tb.at(irun)[2]
    <<" "<<vparer_tb.at(irun)[2]
    <<" "<<vpar_tb.at(irun)[3]
    <<" "<<vparer_tb.at(irun)[3]
    <<" "<<vpar_tb.at(irun)[4]
    <<" "<<vparer_tb.at(irun)[4]
  
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
<<endl;
 }

if(btecint) {
    file<< vtime.at(irun)
    <<" "<<vruns.at(irun)
    <<" "<<vpar_tb.at(irun)[0]
    <<" "<<vparer_tb.at(irun)[0]
    <<" "<<vpar_tb.at(irun)[1]
    <<" "<<vparer_tb.at(irun)[1]
    <<" "<<vpar_tb.at(irun)[2]
    <<" "<<vparer_tb.at(irun)[2]
    <<" "<<vpar_tb.at(irun)[3]
    <<" "<<vparer_tb.at(irun)[3]
    <<" "<<vpar_tb.at(irun)[4]
    <<" "<<vparer_tb.at(irun)[4]
  
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
 }



 } // irun
  } // ifmc
else {
  if(bhead) {
 file<<"       time "
<<"  run   "
<<"   dx_tb "
<<"   dy_tb "
<<"   rx_tb "
<<"   ry_tb "
<<"   rz_tb "
<<"   dx_tecp "
<<"   dy_tecp "
<<"   rz_tecp "
<<"   dx_tecp "
<<"   dy_tecp "
<<"   rz_tecp "
<<"   ch2     "
<<endl;
  }

 for(int irun=0; irun<nruns; irun++) {
   // store only nonzero results
double sum=
     vpar_tb.at(irun)[0]
    +vpar_tb.at(irun)[1]
    +vpar_tb.at(irun)[2]
    +vpar_tb.at(irun)[3]
    +vpar_tb.at(irun)[4]
    +vpar_tecp.at(irun)[0]
    +vpar_tecp.at(irun)[1]
    +vpar_tecp.at(irun)[2]
    +vpar_tecm.at(irun)[0]
    +vpar_tecm.at(irun)[1]
    +vpar_tecm.at(irun)[2];

  if(vruns.at(irun)!=refrunn&&sum==0) continue;

    file<< vtime.at(irun)
    <<" "<<vruns.at(irun)
    <<" "<<vpar_tb.at(irun)[0]
    <<" "<<vpar_tb.at(irun)[1]
    <<" "<<vpar_tb.at(irun)[2]
    <<" "<<vpar_tb.at(irun)[3]
    <<" "<<vpar_tb.at(irun)[4]
  
    <<" "<<vpar_tecp.at(irun)[0]
    <<" "<<vpar_tecp.at(irun)[1]
    <<" "<<vpar_tecp.at(irun)[2]
   
    <<" "<<vpar_tecm.at(irun)[0]
    <<" "<<vpar_tecm.at(irun)[1]
    <<" "<<vpar_tecm.at(irun)[2]
	<<" "<<vch2.at(irun)
<<endl;
 }


    }

  file.close();
  
  if(debug ) cout<<" done write to file>>>> "<<fout<<endl;
  return 0;
}

// enable histograms
void LAS_stabil:: histenable(int op) {
  if(debug) cout<<" enbable histo "<<endl;
  bhist=op;
  thehisto=new LAS_histo();
  thehisto->book();
  thehisto->setdebug(debug);
}

void LAS_stabil:: mcsmear(double*  d, double s) {
  mcsig_mean=1;
  *d=thernd->Gaus(*d,s);
  return;
}
//
void LAS_stabil:: mcoffset(double*  d, double s) {
  *d+=s;
}

//calculate running mean and rms for MC 
void LAS_stabil:: updateresmc(int s) {
  if(debug) cout<<" updateresmc>>" <<s<<endl;
updaterestb(s);
updaterestecp(s);
updaterestecm(s);
 if(btecint) {
updaterestecmint(s);
updaterestecpint(s);
}
 return;
}
//calculate final mean and rms for MC 
void LAS_stabil:: calculatermsmc() {
  if(debug) cout<<" calculate final rms>>" <<endl;
  if(utb<2||utecp<2||utecm<2) return;
  for (int ir=0; ir<nruns;ir++) {
    // tib
    int npar=5;
for(int k=0; k<npar; k++) {
  vparer_tb.at(ir)[k]=sqrt(vparer_tb.at(ir)[k]/(utb-1));
 }// tib

   npar=3;
for(int k=0; k<npar; k++) {
  vparer_tecp.at(ir)[k]=sqrt(vparer_tecp.at(ir)[k]/(utecp-1));
 }// tib
for(int k=0; k<npar; k++) {
  vparer_tecm.at(ir)[k]=sqrt(vparer_tecm.at(ir)[k]/(utecm-1));
 }// tib
  
 if(btecint) {
 npar=27;
for(int k=0; k<npar; k++) {
  vparer_tecpint.at(ir)[k]=sqrt(vparer_tecpint.at(ir)[k]/(utecpint-1));
 }// tib
for(int k=0; k<npar; k++) {
  vparer_tecmint.at(ir)[k]=sqrt(vparer_tecmint.at(ir)[k]/(utecmint-1));
 }// tib
 }// btec

}//ir


 return;
}


//calculate running mean and rms for MC 
void LAS_stabil:: updaterestb(int s) {
  if(debug) cout<<" updaterestb>>" <<s<<endl;

 

  int npar=5;
  if(s==0) {
    // init, set to zero (should be already)
    utb=0;
clearv(&vparer_tb);
clearv(&vparm_tb);
double* perr=0;
double* pm=0;

for (int ir=0; ir<nruns;ir++) {
perr=new double[npar];
pm=new double[npar];
  for(int ip = 0; ip < npar; ++ip){
      perr[ip]=0.;  
      pm[ip]=0.;
 }
 vparer_tb.push_back(perr);
 vparm_tb.push_back(pm);
 }//ir

  
  } else {
    // update
    utb++;
   for (int ir=0; ir<nruns;ir++) {
     if(debug>5)  cout<<" ---update run "<<vruns[ir]<<endl;
      for(int k=0; k<npar; k++) {
	// order is important!
        vparer_tb.at(ir)[k]=vparer_tb.at(ir)[k]+ (utb-1.)/utb *(vpar_tb.at(ir)[k]-vparm_tb.at(ir)[k])*(vpar_tb.at(ir)[k]-vparm_tb.at(ir)[k]);
	vparm_tb.at(ir)[k]=vparm_tb.at(ir)[k]+ (vpar_tb.at(ir)[k]-vparm_tb.at(ir)[k])/utb;
	if(debug>5) cout<<ir<<" "<<k<<"updaterestib "<<vparer_tb.at(ir)[k]<<" "<<vparm_tb.at(ir)[k]<<endl;
      }//k
   }// ir
  }// else
  if(debug)  cout<<"updatetbdone "<<endl;
  return;
}

void LAS_stabil:: updaterestecp(int s) {
 if(debug) cout<<" updaterestecp>>" <<s<<endl;
  int npar=3;
  if(s==0) {
    // init, set to zero (should be already)
    utecp=0;
clearv(&vparer_tecp);
clearv(&vparm_tecp);
double* perr=0;
double* pm=0;

for (int ir=0; ir<nruns;ir++) {
perr=new double[npar];
pm=new double[npar];
  for(int ip = 0; ip < npar; ++ip){
      perr[ip]=0.;  
      pm[ip]=0.;
 }
 vparer_tecp.push_back(perr);
 vparm_tecp.push_back(pm);
   }// ir

  } else {
    // update
    utecp++;
   for (unsigned int ir=0; ir<nruns;ir++) {
      for(int k=0; k<npar; k++) {
	// order is important!
        vparer_tecp.at(ir)[k]=vparer_tecp.at(ir)[k]+ (utecp-1.)/utecp *(vpar_tecp.at(ir)[k]-vparm_tecp.at(ir)[k])*(vpar_tecp.at(ir)[k]-vparm_tecp.at(ir)[k]);
	vparm_tecp.at(ir)[k]=vparm_tecp.at(ir)[k]+ (vpar_tecp.at(ir)[k]-vparm_tecp.at(ir)[k])/utecp;
      }//k
   }// ir
  }// else
  return;
}

void LAS_stabil:: updaterestecm(int s) {
 if(debug) cout<<" updaterestecm>>" <<s<<endl;
  int npar=3;
  if(s==0) {
    // init, set to zero (should be already)
    utecm=0;
clearv(&vparer_tecm);
clearv(&vparm_tecm);
double* perr=0;
double* pm=0;

for (int ir=0; ir<nruns;ir++) {
perr=new double[npar];
pm=new double[npar];
  for(int ip = 0; ip < npar; ++ip){
      perr[ip]=0.;  
      pm[ip]=0.;
 }
 vparer_tecm.push_back(perr);
 vparm_tecm.push_back(pm);
   }// ir

  } else {
    // update
    utecm++;
   for (unsigned int ir=0; ir<nruns;ir++) {
      for(int k=0; k<npar; k++) {
	// order is important!
        vparer_tecm.at(ir)[k]=vparer_tecm.at(ir)[k]+ (utecm-1.)/utecm *(vpar_tecm.at(ir)[k]-vparm_tecm.at(ir)[k])*(vpar_tecm.at(ir)[k]-vparm_tecm.at(ir)[k]);
	vparm_tecm.at(ir)[k]=vparm_tecm.at(ir)[k]+ (vpar_tecm.at(ir)[k]-vparm_tecm.at(ir)[k])/utecm;
      }//k
   }// ir
  }// else
  return;
}
//
void LAS_stabil:: updaterestecmint(int s) {
 if(debug) cout<<" updaterestecmint>>" <<s<<endl;
  int npar=27;
  if(s==0) {
    // init, set to zero (should be already)
    utecmint=0;
clearv(&vparer_tecmint);
clearv(&vparm_tecmint);
double* perr=0;
double* pm=0;

for (int ir=0; ir<nruns;ir++) {
perr=new double[npar];
pm=new double[npar];
  for(int ip = 0; ip < npar; ++ip){
      perr[ip]=0.;  
      pm[ip]=0.;
 }
 vparer_tecmint.push_back(perr);
 vparm_tecmint.push_back(pm);
   }// ir

  } else {
    // update
    utecmint++;
   for (unsigned int ir=0; ir<nruns;ir++) {
      for(int k=0; k<npar; k++) {
	// order is important!
        vparer_tecmint.at(ir)[k]=vparer_tecmint.at(ir)[k]+ (utecmint-1.)/utecmint *(vpar_tecmint.at(ir)[k]-vparm_tecmint.at(ir)[k])*(vpar_tecmint.at(ir)[k]-vparm_tecmint.at(ir)[k]);
	vparm_tecmint.at(ir)[k]=vparm_tecmint.at(ir)[k]+ (vpar_tecmint.at(ir)[k]-vparm_tecmint.at(ir)[k])/utecmint;
 
	//	cout<<ir<<" k="<<k<<" KKKKK updatetecmint "<<utecmint<<"  err="<< vparer_tecmint.at(ir)[k]<<" mean ="<<vparm_tecmint.at(ir)[k]<<endl;
     }//k
   }// ir
  }// else
  return;
}
//
void LAS_stabil::updaterestecpint(int s) {
 if(debug) cout<<" updaterestecpint>>" <<s<<endl;
  int npar=27;
  if(s==0) {
    // init, set to zero (should be already)
    utecpint=0;
clearv(&vparer_tecpint);
clearv(&vparm_tecpint);
double* perr=0;
double* pm=0;

for (int ir=0; ir<nruns;ir++) {
perr=new double[npar];
pm=new double[npar];
  for(int ip = 0; ip < npar; ++ip){
      perr[ip]=0.;  
      pm[ip]=0.;
 }
 vparer_tecpint.push_back(perr);
 vparm_tecpint.push_back(pm);
   }// ir

  } else {
    // update
    utecpint++;
   for (unsigned int ir=0; ir<nruns;ir++) {
      for(int k=0; k<npar; k++) {
	// order is important!
        vparer_tecpint.at(ir)[k]=vparer_tecpint.at(ir)[k]+ (utecpint-1.)/utecpint *(vpar_tecpint.at(ir)[k]-vparm_tecpint.at(ir)[k])*(vpar_tecpint.at(ir)[k]-vparm_tecpint.at(ir)[k]);
	vparm_tecpint.at(ir)[k]=vparm_tecpint.at(ir)[k]+ (vpar_tecpint.at(ir)[k]-vparm_tecpint.at(ir)[k])/utecpint;

	//	cout<<ir<<" k="<<k<<" KKKKK updatetecpint "<<utecpint<<"  err="<< vparer_tecpint.at(ir)[k]<<" mean ="<<vparm_tecpint.at(ir)[k]<<endl;

      }//k
   }// ir
  }// else
  return;
}

// summurize TECin alignment from 9 disks to 3 parameters
// for all runs
void LAS_stabil::tecintsummary() {
  if(debug) cout<<" Las_stabil::tecintsummary>>>>"<<endl;
 
clearv(&vpar_tecpinta);
clearv(&vpar_tecminta);
clearv(&vparer_tecpinta);
clearv(&vparer_tecminta);

 int ndmax=9;
 int npar=3;
 

 // cout<<" size "<<vpar_tecpint.size()<<endl;
 
 for (unsigned int run=0; run<nruns;run++) {

   //   cout<<" KKKK"<< run<<endl;
   double* pstecp=new double[npar];
   double* psertecp=new double[npar];
   double* pstecm=new double[npar];
   double* psertecm=new double[npar];
   cleara(pstecp,npar);
   cleara(psertecp,npar);
   cleara(pstecm,npar);
   cleara(psertecm,npar);   

  // average out all disks
     for(int dn=0; dn<ndmax;dn++) {
       for(int ip=0; ip<npar; ip++) {
	 // cout<<run<<" KKKKK dn="<<dn<<" ip="<<ip<<" tecpintper="<<vparer_tecpint.at(run)[dn*npar+ip]<<" tecmintper="<<vparer_tecmint.at(run)[dn*npar+ip]<<endl;

	 pstecp[ip]=pstecp[ip]+vpar_tecpint.at(run)[dn*npar+ip];
         psertecp[ip]=psertecp[ip]+vparer_tecpint.at(run)[dn*npar+ip]*vparer_tecpint.at(run)[dn*npar+ip];

	 pstecm[ip]=pstecm[ip]+vpar_tecmint.at(run)[dn*npar+ip];
         psertecm[ip]=psertecm[ip]+vparer_tecmint.at(run)[dn*npar+ip]*vparer_tecmint.at(run)[dn*npar+ip];
    }
  }// dn
 for(int ip=0; ip<npar; ip++) {
   pstecp[ip]/=ndmax;
   psertecp[ip]=sqrt(psertecp[ip])/ndmax;
   pstecm[ip]/=ndmax;
   psertecm[ip]= sqrt(psertecm[ip])/ndmax;
 }
  vpar_tecpinta.push_back(pstecp);
  vparer_tecpinta.push_back(psertecp);
  vpar_tecminta.push_back(pstecm);
  vparer_tecminta.push_back(psertecm);
	   } //irun

	   if(debug) cout<<" >>>>tecintsummary"<<endl;
	   return;
}
