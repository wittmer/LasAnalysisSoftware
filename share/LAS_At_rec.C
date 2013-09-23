#include <sstream>
#include <iostream>
#include <ctime>


#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>


#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMatrixT.h"
#include "TMultiGraph.h"
#include "TGraph.h"

#include "LASGlobalData.h"
#include "LASGlobalDataLoop.h"
#define __AVECROOT__
#include "Avec.h"
#include "Avec2D.h"
#include "LAS_vectorfloat_tools.h"
#include "LAS_basic_tools.h"
#include "LAS_globaldata_tools.h"
#include "LAS_RDC_tools.h"
#include "LAS_alpar.h"
#include "LAS_Tec_Reconstruction.h"
 
LASGlobalData<double> AT_spot_rec(const AtPar& par) //calculation of beam spot positions at TIB&TOB modules
{
  LASGlobalData<double> result;

  LASGlobalData<double> nzpos = zpos() - LAS::z_at_bs;

  LASGlobalDataLoop tib_loop(LASGlobalDataLoop::TIB);

  do{
    int bn = tib_loop.get_beam();
    int dn = tib_loop.get_zpos();
    double zv = nzpos.GetEntry(2,-1,bn,dn);
    double thv = theta().GetEntry(2,-1,bn,dn);
    double rv = r0().GetEntry(2,-1,bn,dn);
    tib_loop.GetEntry(result) = 
      par.beam_a[bn]*zv + par.beam_b[bn] 
      - par.tib_dx*sin(thv)/rv + par.tib_dy*cos(thv)/rv 
      + par.tib_rx*cos(thv)*zv/rv + par.tib_ry*sin(thv)*zv/rv 
      + par.tib_rz+ par.tib_tz*zv; 
  }while(tib_loop.next());
  
  LASGlobalDataLoop tecp_loop(LASGlobalDataLoop::TEC_PLUS_AT);

  do{
    tecp_loop.GetEntry(result) = 0;
//     int bn = tecp_loop.get_beam();
//     int dn = tecp_loop.get_zpos();
//     tecp_loop.GetEntry(result) = par.beam_a[bn]*zpos().GetEntry(0,-1,bn,dn) + par.beam_b[bn] - par.tecp_dx*sin(theta().GetEntry(0,-1,bn,dn))/r0().GetEntry(0,-1,bn,dn) + par.tecp_dy*cos(theta().GetEntry(0,-1,bn,dn))/r0().GetEntry(0,-1,bn,dn) + par.tecp_rx*cos(theta().GetEntry(0,-1,bn,dn))*zpos().GetEntry(0,-1,bn,dn)/r0().GetEntry(0,-1,bn,dn) + par.tecp_ry*sin(theta().GetEntry(0,-1,bn,dn))*zpos().GetEntry(0,-1,bn,dn)/r0().GetEntry(0,-1,bn,dn) + par.tecp_rz+ par.tecp_tz*zpos().GetEntry(0,-1,bn,dn); 
  }while(tecp_loop.next());
  
  LASGlobalDataLoop tecm_loop(LASGlobalDataLoop::TEC_MINUS_AT);

  do{
    tecm_loop.GetEntry(result) = 0;
//     int bn = tecm_loop.get_beam();
//     int dn = tecm_loop.get_zpos();
//     tecm_loop.GetEntry(result) = par.beam_a[bn]*zpos().GetEntry(1,-1,bn,dn) + par.beam_b[bn] - par.tecm_dx*sin(theta().GetEntry(1,-1,bn,dn))/r0().GetEntry(1,-1,bn,dn) + par.tecm_dy*cos(theta().GetEntry(1,-1,bn,dn))/r0().GetEntry(1,-1,bn,dn) + par.tecm_rx*cos(theta().GetEntry(1,-1,bn,dn))*zpos().GetEntry(1,-1,bn,dn)/r0().GetEntry(1,-1,bn,dn) + par.tecm_ry*sin(theta().GetEntry(1,-1,bn,dn))*zpos().GetEntry(1,-1,bn,dn)/r0().GetEntry(1,-1,bn,dn) + par.tecm_rz+ par.tecm_tz*zpos().GetEntry(1,-1,bn,dn); 
  }while(tecm_loop.next());
  
  LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);

  do{
    int bn = tob_loop.get_beam();
    int dn = tob_loop.get_zpos();
    double zv = nzpos.GetEntry(3,-1,bn,dn);
    tob_loop.GetEntry(result) = +par.beam_a[bn]*zv + par.beam_b[bn]; 
  }while(tob_loop.next());

  return result;
}

void AT_chi2(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,  AtPar& par) //calculation of chi2/ndf at TIB&TOB
{

  LASGlobalData<double> recpos = AT_spot_rec(par);
  LASGlobalData<double> resid = recpos - dif;
  LASGlobalData<double> res2 = pow(resid/err2,2);
  LASGlobalData<double> dif2 = pow(dif/err2,2);

  int nchi_tib,nchi_tob,nchi,nchib[8];
  
  nchi = 0;
  nchi_tib = 0;
  par.tib_chi = 0.;
  par.tib_chi0 = 0.;
  par.AT_chi = 0.;
  par.AT_chi0 = 0.;
  for(int bn=0;bn<8;bn++){
    nchib[bn] = 0;
    par.b_chi[bn] = 0.;
    par.b_chi0[bn] = 0.;
  }

  //std::cout << "tib res2 : \n";
  LASGlobalDataLoop tib_loop(LASGlobalDataLoop::TIB);
  do{
    int bn = tib_loop.get_beam();
    //nchib[bn] = 0;
    //par.b_chi[bn] = 0.;
    //par.b_chi0[bn] = 0.;
    //int dn = tib_loop.get_zpos();
    //std::cout << tib_loop.GetEntry(res2) << "  ";
    nchi = nchi + 1;
    nchib[bn] = nchib[bn] + 1;
    nchi_tib = nchi_tib + 1;
    par.AT_chi += tib_loop.GetEntry(res2);
    par.AT_chi0 += tib_loop.GetEntry(dif2);
    par.tib_chi += tib_loop.GetEntry(res2);
    par.tib_chi0 += tib_loop.GetEntry(dif2);
    par.b_chi[bn] += tib_loop.GetEntry(res2);
    par.b_chi0[bn] += tib_loop.GetEntry(dif2);
    }while(tib_loop.next());
  //std::cout << "\ntib_chi: " << par.tib_chi << std::endl;
  par.tib_chi = sqrt(par.tib_chi/nchi_tib);
  par.tib_chi0 = sqrt(par.tib_chi0/nchi_tib);
  //cout<<"irun: "<<irun<<"  nchi_tib: "<<nchi_tib<<"  chi_tib: "<<chi_tib[irun]<<endl;   
 
  nchi_tob = 0;
  par.tob_chi = 0.;
  par.tob_chi0 = 0.;

  LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);
  do{
    int bn = tob_loop.get_beam();
    //int dn = tob_loop.get_zpos();
    nchi = nchi + 1;
    nchib[bn] = nchib[bn] + 1;
    nchi_tob = nchi_tob + 1;
    par.AT_chi += tob_loop.GetEntry(res2);
    par.AT_chi0 += tob_loop.GetEntry(dif2);
    par.tob_chi += tob_loop.GetEntry(res2);
    par.tob_chi0 += tob_loop.GetEntry(dif2);
    par.b_chi[bn] += tob_loop.GetEntry(res2);
    par.b_chi0[bn] += tob_loop.GetEntry(dif2);
    }while(tob_loop.next());

  par.tob_chi = sqrt(par.tob_chi/nchi_tob);
  par.tob_chi0 = sqrt(par.tob_chi0/nchi_tob);
  //cout<<"irun: "<<irun<<"  nchi_tob: "<<nchi_tob<<"  chi_tob: "<<chi_tob[irun]<<endl;  
 
  par.AT_chi = sqrt(par.AT_chi/(nchi-22));
  par.AT_chi0 = sqrt(par.AT_chi0/nchi);

 
  for(int bnn = 0; bnn < 8;  ++bnn){
    par.b_chi[bnn] = sqrt(par.b_chi[bnn]/(nchib[bnn]-2));
    par.b_chi0[bnn] = sqrt(par.b_chi0[bnn]/nchib[bnn]);
  }
  
}

void DrawAtPar(std::vector<AtPar>& parlist, Avec& time)
{
  Avec tib_dx;
  Avec tib_dy;
  Avec tib_rx;
  Avec tib_ry;
  Avec tib_rz;
  Avec tib_tz;
  Avec tib_chi;
  Avec tib_chi0;
  Avec tob_chi;
  Avec tob_chi0;
  Avec2D beam_a;
  Avec2D beam_b;

  for(unsigned int i = 0; i < parlist.size(); i++){
    tib_dx.push_back(parlist[i].tib_dx);
    tib_dy.push_back(parlist[i].tib_dy);
    tib_rx.push_back(parlist[i].tib_rx);
    tib_ry.push_back(parlist[i].tib_ry);
    tib_rz.push_back(parlist[i].tib_rz);
    tib_tz.push_back(parlist[i].tib_tz);
    tib_chi.push_back(parlist[i].tib_chi);
    tib_chi0.push_back(parlist[i].tib_chi0);
    tob_chi.push_back(parlist[i].tob_chi);
    tob_chi0.push_back(parlist[i].tob_chi0);
    beam_a.push_back(parlist[i].beam_a);
    beam_b.push_back(parlist[i].beam_b);
  }
  beam_a = vtrans(beam_a);
  beam_b = vtrans(beam_b);

  // Move z origin to beamsplitter position
  //beam_b += beam_a * LAS::z_at_bs;
  //tib_dx -= tib_ry * LAS::z_at_bs;
  //tib_dy += tib_rx * LAS::z_at_bs;
  //tib_rz += tib_tz * LAS::z_at_bs;

  // Rescale to units for display
  tib_dx *= 1e3;
  tib_dy *= 1e3;
  tib_rx *= 1e6;
  tib_ry *= 1e6;
  tib_rz *= 1e6;
  tib_tz *= 1e9;
  beam_a *= 1e6;

  new TCanvas("TIB_AT","AT"); 
  TMultiGraph* gr = avec_draw(time, (tib_dx | tib_dy | tib_rx| tib_ry | tib_rz | tib_tz ),"TIB vs TOB (16 + 6p fit)","Date","#Deltax, #Deltay [#mum], Rx,Ry,Rz [#murad], Tz[#murad/m] ","AP");
  //TMultiGraph* gr = avec_draw(time, (dx | dy ),"TEC_AT wrt TOB Stability over Run (11 par)","Block Nr","#Deltax, #Deltay [#mum], Rz [#murad] ","AP");
  //gr->GetYaxis()->SetRangeUser(-5,5);
  gr->GetXaxis()->SetNdivisions(505, kTRUE);
  gr->GetXaxis()->SetTimeDisplay(1);
  gr->GetXaxis()->SetTimeFormat("%H %d/%m");
  //avec_draw(time, dy,"dy vs time","time","dy","AP");
  //avec_draw(time, dphi,"dx vs time","time","dphi","AP");
  //avec_draw(dx, dy,"dx vs dy","dx","dy","AP");
 
  std::vector<std::string> Legend;
  Legend.push_back("#Delta x");
  Legend.push_back("#Delta y");
  Legend.push_back("Rx");
  Legend.push_back("Ry");
  Legend.push_back("Rz");
  Legend.push_back("Tz");
  // Legend.push_back("#Deltax");
  // Legend.push_back("#Deltay");
  // Legend.push_back("Rx");
  // Legend.push_back("Ry");
  // Legend.push_back("Rz");
  AddLegend(gr, Legend);
  //TLine *L1=new TLine(time.front(),1, time.back(), 1);
  //TLine *L2=new TLine(time.front(), -1, time.back(), -1);
  //L1->Draw();
  //L2->Draw();

  //avec_draw(xvals, tib_tz,"TIB_Tz","Block Nr","Tz [micron]","AP");
  //avec_draw(tib_tz,"TIB_Tz","Block Nr","Tz [micron]","AP");

  //std::cout << "tib_chi:\n" << tib_chi << std::endl;
  TCanvas* cv_beampar = new TCanvas("beampar","beampar");
  cv_beampar->Divide(1,2);
  cv_beampar->cd(1);
  avec_draw(time, beam_a, "Beam Slope", "Date", "[#murad/m]", "AP");
  cv_beampar->cd(2);
  avec_draw(time, beam_b, "Beam Offset", "Date", "", "AP");

  new TCanvas("tibchi","tibchi");
  avec_draw(time, (tib_chi | tib_chi0 | tob_chi | tob_chi0), "TIB chi2", "Date", "", "AP");
}

// reconstruction of TIB wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TIB2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results)
{
  // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

  //LASGlobalData<double> nzpos = zpos();
  LASGlobalData<double> nzpos = zpos() - LAS::z_at_bs;

  double bs[22];
  double db[22][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TIB, TOB); k=0..7 (beam number); l=0..5 (module number)

  int npb = 16;// number of beam parameters
 
  LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
  LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
  LASGlobalData<double> dev2 = cos(theta()) * r0() * nzpos / err2;
  LASGlobalData<double> dev3 = sin(theta()) * r0() * nzpos / err2;
  LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
  LASGlobalData<double> dev5 = nzpos / err2;

//   LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
//   LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
//   LASGlobalData<double> par2 = dif * cos(theta()) * r0() * nzpos / err2;
//   LASGlobalData<double> par3 = dif * sin(theta()) * r0() * nzpos / err2;
//   LASGlobalData<double> par4 = dif / err2;
//   LASGlobalData<double> par5 = dif * nzpos/ err2;

  LASGlobalData<double> par0 = dif * dev0;
  LASGlobalData<double> par1 = dif * dev1;
  LASGlobalData<double> par2 = dif * dev2;
  LASGlobalData<double> par3 = dif * dev3;
  LASGlobalData<double> par4 = dif * dev4;
  LASGlobalData<double> par5 = dif * dev5;


  for(int ip = 0; ip < 22;  ++ip){ 
    bs[ip] = 0.;
  }
  for(int bn = 0; bn < 8;  ++bn){
    for(int dn = 0; dn < 6;  ++dn){
      for(int ip = 0; ip < 22;  ++ip){ 
	db[ip][0][bn][dn] = 0.;
	db[ip][1][bn][dn] = 0.;
      }
      if(err2.GetEntry(2,-1,bn,dn) < 1000.){
	    
	bs[2*bn] += par5.GetEntry(2,-1,bn,dn);
	bs[2*bn+1] += par4.GetEntry(2,-1,bn,dn);
	
	bs[npb+0] += par0.GetEntry(2,-1,bn,dn);
	bs[npb+1] += par1.GetEntry(2,-1,bn,dn);
	bs[npb+2] += par2.GetEntry(2,-1,bn,dn);
	bs[npb+3] += par3.GetEntry(2,-1,bn,dn);
	bs[npb+4] += par4.GetEntry(2,-1,bn,dn);
	bs[npb+5] += par5.GetEntry(2,-1,bn,dn);
        
	db[2*bn][0][bn][dn] = dev5.GetEntry(2,-1,bn,dn);
	db[2*bn+1][0][bn][dn] = dev4.GetEntry(2,-1,bn,dn);
	
	db[npb+0][0][bn][dn] = dev0.GetEntry(2,-1,bn,dn);
	db[npb+1][0][bn][dn] = dev1.GetEntry(2,-1,bn,dn);
	db[npb+2][0][bn][dn] = dev2.GetEntry(2,-1,bn,dn);
	db[npb+3][0][bn][dn] = dev3.GetEntry(2,-1,bn,dn);
	db[npb+4][0][bn][dn] = dev4.GetEntry(2,-1,bn,dn);
	db[npb+5][0][bn][dn] = dev5.GetEntry(2,-1,bn,dn);
      }
      if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
 	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
 	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
 	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
 	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);

	//bs[2*bn] -= par5.GetEntry(3,-1,bn,dn);
	//bs[2*bn+1] -= par4.GetEntry(3,-1,bn,dn);
	
	//db[2*bn][1][bn][dn] = -dev5.GetEntry(3,-1,bn,dn);
	//db[2*bn+1][1][bn][dn] = -dev4.GetEntry(3,-1,bn,dn);
      }
    }
  }
   
  TMatrixT<double> h(22,22);
  
  for(int ix = 0; ix < 22; ++ix){
    for(int iy = 0; iy < 22; ++iy){
      h(ix,iy) = 0.;
    }
  }
  
  LASGlobalData<double> haa =   nzpos  * nzpos  / err2;
  LASGlobalData<double> hab =   nzpos  / err2;
  LASGlobalData<double> hba =   nzpos  / err2;
  LASGlobalData<double> hbb =   LASGlobalData<double>(1.) / err2;


  LASGlobalData<double> h00 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) / err2;
  LASGlobalData<double> h01 =   sin(theta()) * cos(theta()) / err2;
  LASGlobalData<double> h02 =   sin(theta()) * cos(theta()) * nzpos  / err2;
  LASGlobalData<double> h03 =   sin(theta()) * sin(theta()) * nzpos  / err2;
  LASGlobalData<double> h04 =   sin(theta()) / err2;
  LASGlobalData<double> h05 =   sin(theta()) * nzpos  / err2;

  LASGlobalData<double> h10 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) / err2;
  LASGlobalData<double> h11 =   cos(theta()) * cos(theta()) / err2;
  LASGlobalData<double> h12 =   cos(theta()) * cos(theta()) * nzpos  / err2;
  LASGlobalData<double> h13 =   cos(theta()) * sin(theta()) * nzpos  / err2;
  LASGlobalData<double> h14 =   cos(theta()) / err2;
  LASGlobalData<double> h15 =   cos(theta()) * nzpos  / err2;

  LASGlobalData<double> h20 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) * nzpos  / err2;
  LASGlobalData<double> h21 =   cos(theta()) * cos(theta()) * nzpos  / err2;
  LASGlobalData<double> h22 =   cos(theta()) * cos(theta()) * nzpos  * nzpos  / err2;
  LASGlobalData<double> h23 =   cos(theta()) * sin(theta()) * nzpos  * nzpos  / err2;
  LASGlobalData<double> h24 =   cos(theta()) * nzpos / err2;
  LASGlobalData<double> h25 =   cos(theta()) * nzpos  * nzpos  / err2;

  LASGlobalData<double> h30 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) * nzpos  / err2;
  LASGlobalData<double> h31 =   sin(theta()) * cos(theta()) * nzpos  / err2;
  LASGlobalData<double> h32 =   sin(theta()) * cos(theta()) * nzpos  * nzpos  / err2;
  LASGlobalData<double> h33 =   sin(theta()) * sin(theta()) * nzpos  * nzpos  / err2;
  LASGlobalData<double> h34 =   sin(theta()) * nzpos / err2;
  LASGlobalData<double> h35 =   sin(theta()) * nzpos  * nzpos  / err2;

  LASGlobalData<double> h40 = LASGlobalData<double>(0.) - sin(theta()) / r0() / err2;
  LASGlobalData<double> h41 =   cos(theta()) / r0() / err2;
  LASGlobalData<double> h42 =   cos(theta()) * nzpos  / r0() / err2;
  LASGlobalData<double> h43 =   sin(theta()) * nzpos  / r0() / err2;
  LASGlobalData<double> h44 =   LASGlobalData<double>(1.) / err2;
  LASGlobalData<double> h45 =   nzpos  / err2;

  LASGlobalData<double> h50 = LASGlobalData<double>(0.) - sin(theta()) * nzpos  / r0() / err2;
  LASGlobalData<double> h51 =   cos(theta()) * nzpos  / r0() / err2;
  LASGlobalData<double> h52 =   cos(theta()) * nzpos  * nzpos  / r0() / err2;
  LASGlobalData<double> h53 =   sin(theta()) * nzpos  * nzpos  / r0() / err2;
  LASGlobalData<double> h54 =   nzpos  / err2;
  LASGlobalData<double> h55 =   nzpos  * nzpos  / err2;


//   LASGlobalData<double> haa =   zpos() * zpos() / err2;
//   LASGlobalData<double> hab =   zpos() / err2;
//   LASGlobalData<double> hba =   zpos() / err2;
//   LASGlobalData<double> hbb =   LASGlobalData<double>(1.) / err2;


//   LASGlobalData<double> h00 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) / err2;
//   LASGlobalData<double> h01 =   sin(theta()) * cos(theta()) / err2;
//   LASGlobalData<double> h02 =   sin(theta()) * cos(theta()) * zpos() / err2;
//   LASGlobalData<double> h03 =   sin(theta()) * sin(theta()) * zpos() / err2;
//   LASGlobalData<double> h04 =   sin(theta()) / err2;
//   LASGlobalData<double> h05 =   sin(theta()) * zpos() / err2;

//   LASGlobalData<double> h10 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) / err2;
//   LASGlobalData<double> h11 =   cos(theta()) * cos(theta()) / err2;
//   LASGlobalData<double> h12 =   cos(theta()) * cos(theta()) * zpos() / err2;
//   LASGlobalData<double> h13 =   cos(theta()) * sin(theta()) * zpos() / err2;
//   LASGlobalData<double> h14 =   cos(theta()) / err2;
//   LASGlobalData<double> h15 =   cos(theta()) * zpos() / err2;

//   LASGlobalData<double> h20 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) * zpos() / err2;
//   LASGlobalData<double> h21 =   cos(theta()) * cos(theta()) * zpos() / err2;
//   LASGlobalData<double> h22 =   cos(theta()) * cos(theta()) * zpos() * zpos() / err2;
//   LASGlobalData<double> h23 =   cos(theta()) * sin(theta()) * zpos() * zpos() / err2;
//   LASGlobalData<double> h24 =   cos(theta()) * zpos()/ err2;
//   LASGlobalData<double> h25 =   cos(theta()) * zpos() * zpos() / err2;

//   LASGlobalData<double> h30 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) * zpos() / err2;
//   LASGlobalData<double> h31 =   sin(theta()) * cos(theta()) * zpos() / err2;
//   LASGlobalData<double> h32 =   sin(theta()) * cos(theta()) * zpos() * zpos() / err2;
//   LASGlobalData<double> h33 =   sin(theta()) * sin(theta()) * zpos() * zpos() / err2;
//   LASGlobalData<double> h34 =   sin(theta()) * zpos()/ err2;
//   LASGlobalData<double> h35 =   sin(theta()) * zpos() * zpos() / err2;

//   LASGlobalData<double> h40 = LASGlobalData<double>(0.) - sin(theta()) / r0() / err2;
//   LASGlobalData<double> h41 =   cos(theta()) / r0() / err2;
//   LASGlobalData<double> h42 =   cos(theta()) * zpos() / r0() / err2;
//   LASGlobalData<double> h43 =   sin(theta()) * zpos() / r0() / err2;
//   LASGlobalData<double> h44 =   LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> h45 =   zpos() / err2;

//   LASGlobalData<double> h50 = LASGlobalData<double>(0.) - sin(theta()) * zpos() / r0() / err2;
//   LASGlobalData<double> h51 =   cos(theta()) * zpos() / r0() / err2;
//   LASGlobalData<double> h52 =   cos(theta()) * zpos() * zpos() / r0() / err2;
//   LASGlobalData<double> h53 =   sin(theta()) * zpos() * zpos() / r0() / err2;
//   LASGlobalData<double> h54 =   zpos() / err2;
//   LASGlobalData<double> h55 =   zpos() * zpos() / err2;


  for(int bn = 0; bn < 8;  ++bn){
    for(int dn = 0; dn < 6;  ++dn){
      
      if(err2.GetEntry(2,-1,bn,dn) < 1000.){
	    
	h(2*bn,2*bn)     += haa.GetEntry(2,-1,bn,dn);
	h(2*bn,2*bn+1)   += hab.GetEntry(2,-1,bn,dn);
	h(2*bn+1,2*bn)   += hba.GetEntry(2,-1,bn,dn);
	h(2*bn+1,2*bn+1) += hbb.GetEntry(2,-1,bn,dn);
	
	h(npb+0,npb+0)   += h00.GetEntry(2,-1,bn,dn);
	h(npb+0,npb+1)   += h01.GetEntry(2,-1,bn,dn);
	h(npb+0,npb+2)   += h02.GetEntry(2,-1,bn,dn);
	h(npb+0,npb+3)   += h03.GetEntry(2,-1,bn,dn);
	h(npb+0,npb+4)   += h04.GetEntry(2,-1,bn,dn);
	h(npb+0,npb+5)   += h05.GetEntry(2,-1,bn,dn);
	h(npb+0,2*bn)    += h05.GetEntry(2,-1,bn,dn);
	h(npb+0,2*bn+1)  += h04.GetEntry(2,-1,bn,dn);
	
	h(npb+1,npb+0)   += h10.GetEntry(2,-1,bn,dn);
	h(npb+1,npb+1)   += h11.GetEntry(2,-1,bn,dn);
	h(npb+1,npb+2)   += h12.GetEntry(2,-1,bn,dn);
	h(npb+1,npb+3)   += h13.GetEntry(2,-1,bn,dn);
	h(npb+1,npb+4)   += h14.GetEntry(2,-1,bn,dn);
	h(npb+1,npb+5)   += h15.GetEntry(2,-1,bn,dn);
	h(npb+1,2*bn)    += h15.GetEntry(2,-1,bn,dn);
	h(npb+1,2*bn+1)  += h14.GetEntry(2,-1,bn,dn);
	
	h(npb+2,npb+0)   += h20.GetEntry(2,-1,bn,dn);
	h(npb+2,npb+1)   += h21.GetEntry(2,-1,bn,dn);
	h(npb+2,npb+2)   += h22.GetEntry(2,-1,bn,dn);
	h(npb+2,npb+3)   += h23.GetEntry(2,-1,bn,dn);
	h(npb+2,npb+4)   += h24.GetEntry(2,-1,bn,dn);
	h(npb+2,npb+5)   += h25.GetEntry(2,-1,bn,dn);
	h(npb+2,2*bn)    += h25.GetEntry(2,-1,bn,dn);
	h(npb+2,2*bn+1)  += h24.GetEntry(2,-1,bn,dn);
	
	h(npb+3,npb+0)   += h30.GetEntry(2,-1,bn,dn);
	h(npb+3,npb+1)   += h31.GetEntry(2,-1,bn,dn);
	h(npb+3,npb+2)   += h32.GetEntry(2,-1,bn,dn);
	h(npb+3,npb+3)   += h33.GetEntry(2,-1,bn,dn);
	h(npb+3,npb+4)   += h34.GetEntry(2,-1,bn,dn);
	h(npb+3,npb+5)   += h35.GetEntry(2,-1,bn,dn);
	h(npb+3,2*bn)    += h35.GetEntry(2,-1,bn,dn);
	h(npb+3,2*bn+1)  += h34.GetEntry(2,-1,bn,dn);
	
	h(npb+4,npb+0)   += h40.GetEntry(2,-1,bn,dn);
	h(npb+4,npb+1)   += h41.GetEntry(2,-1,bn,dn);
	h(npb+4,npb+2)   += h42.GetEntry(2,-1,bn,dn);
	h(npb+4,npb+3)   += h43.GetEntry(2,-1,bn,dn);
	h(npb+4,npb+4)   += h44.GetEntry(2,-1,bn,dn);
	h(npb+4,npb+5)   += h45.GetEntry(2,-1,bn,dn);
	h(npb+4,2*bn)    += h45.GetEntry(2,-1,bn,dn);
	h(npb+4,2*bn+1)  += h44.GetEntry(2,-1,bn,dn);
	
	h(npb+5,npb+0)   += h50.GetEntry(2,-1,bn,dn);
	h(npb+5,npb+1)   += h51.GetEntry(2,-1,bn,dn);
	h(npb+5,npb+2)   += h52.GetEntry(2,-1,bn,dn);
	h(npb+5,npb+3)   += h53.GetEntry(2,-1,bn,dn);
	h(npb+5,npb+4)   += h54.GetEntry(2,-1,bn,dn);
	h(npb+5,npb+5)   += h55.GetEntry(2,-1,bn,dn);
	h(npb+5,2*bn)    += h55.GetEntry(2,-1,bn,dn);
	h(npb+5,2*bn+1)  += h54.GetEntry(2,-1,bn,dn);
	
	
	h(2*bn,npb+0)   += h50.GetEntry(2,-1,bn,dn);
	h(2*bn,npb+1)   += h51.GetEntry(2,-1,bn,dn);
	h(2*bn,npb+2)   += h52.GetEntry(2,-1,bn,dn);
	h(2*bn,npb+3)   += h53.GetEntry(2,-1,bn,dn);
	h(2*bn,npb+4)   += h54.GetEntry(2,-1,bn,dn);
	h(2*bn,npb+5)   += h55.GetEntry(2,-1,bn,dn);
	
	h(2*bn+1,npb+0) += h40.GetEntry(2,-1,bn,dn);
	h(2*bn+1,npb+1) += h41.GetEntry(2,-1,bn,dn);
	h(2*bn+1,npb+2) += h42.GetEntry(2,-1,bn,dn);
	h(2*bn+1,npb+3) += h43.GetEntry(2,-1,bn,dn);
	h(2*bn+1,npb+4) += h44.GetEntry(2,-1,bn,dn);
	h(2*bn+1,npb+5) += h45.GetEntry(2,-1,bn,dn);
      }
	
      if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
      }
    }
  }  

  
  h.Invert();

  double cp[22][2][8][6];
  double up[22][22];
  double ps[22];

  for(int ip = 0; ip < 22; ++ip){
    for(int ib = 0; ib < 2; ++ib){
      for(int bn = 0; bn < 8;  ++bn){
	for(int dn = 0; dn < 6;  ++dn){
	  cp[ip][ib][bn][dn] = 0.;
	  for(int ix = 0; ix < 22; ++ix){
	    cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
	    }
	}
      }
    }
  }

  for(int ip = 0; ip < 22; ++ip){
    for(int ix = 0; ix < 22; ++ix){
      up[ip][ix] = 0.;
      for(int ib = 0; ib < 2; ++ib){
	for(int bn = 0; bn < 8;  ++bn){
	  for(int dn = 0; dn < 6;  ++dn){
	    if( ib == 0){
	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(2,-1,bn,dn);
	    }
	    if( ib == 1 ){
	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
	    }
	  }
	}
      }
    }
  }
 
  for(int bn = 0; bn < 8; ++bn){
    results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
    results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
  }

  results.er_tib_dx = sqrt(up[npb+0][npb+0]);
  results.er_tib_dy = sqrt(up[npb+1][npb+1]);
  results.er_tib_rx = sqrt(up[npb+2][npb+2]);
  results.er_tib_ry = sqrt(up[npb+3][npb+3]);
  results.er_tib_rz = sqrt(up[npb+4][npb+4]);
  results.er_tib_tz = sqrt(up[npb+5][npb+5]);
  results.tib_dxdy  = up[npb+0][npb+1]/results.er_tib_dx/results.er_tib_dy;
  results.tib_dxrx  = up[npb+0][npb+2]/results.er_tib_dx/results.er_tib_rx;
  results.tib_dxry  = up[npb+0][npb+3]/results.er_tib_dx/results.er_tib_ry;
  results.tib_dxrz  = up[npb+0][npb+4]/results.er_tib_dx/results.er_tib_rz;
  results.tib_dxtz  = up[npb+0][npb+5]/results.er_tib_dx/results.er_tib_tz;
  results.tib_dyrx  = up[npb+1][npb+2]/results.er_tib_dy/results.er_tib_rx;
  results.tib_dyry  = up[npb+1][npb+3]/results.er_tib_dy/results.er_tib_ry;
  results.tib_dyrz  = up[npb+1][npb+4]/results.er_tib_dy/results.er_tib_rz;
  results.tib_dytz  = up[npb+1][npb+5]/results.er_tib_dy/results.er_tib_tz;
  results.tib_rxry  = up[npb+2][npb+3]/results.er_tib_rx/results.er_tib_ry;
  results.tib_rxrz  = up[npb+2][npb+4]/results.er_tib_rx/results.er_tib_rz;
  results.tib_rxtz  = up[npb+2][npb+5]/results.er_tib_rx/results.er_tib_tz;
  results.tib_ryrz  = up[npb+3][npb+4]/results.er_tib_ry/results.er_tib_rz;
  results.tib_rytz  = up[npb+3][npb+5]/results.er_tib_ry/results.er_tib_tz;
  results.tib_rztz  = up[npb+4][npb+5]/results.er_tib_rz/results.er_tib_tz;

  for(int ip = 0; ip < 22; ++ip){
    ps[ip] = 0.;
    for(int ix = 0; ix < 22; ++ix){
      ps[ip] += h(ip,ix)*bs[ix];
    }
  }

  for(int bn = 0; bn < 8; ++bn){
    results.beam_a[bn] = ps[2*bn];
    results.beam_b[bn] = ps[2*bn+1];
  }

  results.tib_dx = ps[npb+0];
  results.tib_dy = ps[npb+1];
  results.tib_rx = ps[npb+2];
  results.tib_ry = ps[npb+3];
  results.tib_rz = ps[npb+4];
  results.tib_tz = ps[npb+5];

}

// reconstruction of TEC+ wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TECP2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results)
{

  // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

  double bs[22];
  double db[22][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TEC+, TOB); k=0..7 (beam number); l=0..5 (module number)

  int npb = 16;// number of beam parameters
 
  LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
  LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
  LASGlobalData<double> par2 = dif * cos(theta()) * r0() * zpos()/ err2;
  LASGlobalData<double> par3 = dif * sin(theta()) * r0() * zpos()/ err2;
  LASGlobalData<double> par4 = dif / err2;
  LASGlobalData<double> par5 = dif * zpos()/ err2;


  LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
  LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
  LASGlobalData<double> dev2 = cos(theta()) * r0() * zpos()/ err2;
  LASGlobalData<double> dev3 = sin(theta()) * r0() * zpos()/ err2;
  LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
  LASGlobalData<double> dev5 = zpos()/ err2;

  for(int ip = 0; ip < 22;  ++ip){ 
    bs[ip] = 0.;
  }
  for(int bn = 0; bn < 8;  ++bn){
    for(int dn = 0; dn < 4;  ++dn){
      for(int ip = 0; ip < 22;  ++ip){ 
	db[ip][0][bn][dn] = 0.;
	db[ip][1][bn][dn] = 0.;
      }
      if(err2.GetEntry(0,-1,bn,dn) < 1000.){
	    
	bs[2*bn] += par5.GetEntry(0,-1,bn,dn);
	bs[2*bn+1] += par4.GetEntry(0,-1,bn,dn);
	
	bs[npb+0] += par0.GetEntry(0,-1,bn,dn);
	bs[npb+1] += par1.GetEntry(0,-1,bn,dn);
	bs[npb+2] += par2.GetEntry(0,-1,bn,dn);
	bs[npb+3] += par3.GetEntry(0,-1,bn,dn);
	bs[npb+4] += par4.GetEntry(0,-1,bn,dn);
	bs[npb+5] += par5.GetEntry(0,-1,bn,dn);
        
	db[2*bn][0][bn][dn] = dev5.GetEntry(0,-1,bn,dn);
	db[2*bn+1][0][bn][dn] = dev4.GetEntry(0,-1,bn,dn);
	
	db[npb+0][0][bn][dn] = dev0.GetEntry(0,-1,bn,dn);
	db[npb+1][0][bn][dn] = dev1.GetEntry(0,-1,bn,dn);
	db[npb+2][0][bn][dn] = dev2.GetEntry(0,-1,bn,dn);
	db[npb+3][0][bn][dn] = dev3.GetEntry(0,-1,bn,dn);
	db[npb+4][0][bn][dn] = dev4.GetEntry(0,-1,bn,dn);
	db[npb+5][0][bn][dn] = dev5.GetEntry(0,-1,bn,dn);
      }
      if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);
      }
    }
  }
   
  TMatrixT<double> h(22,22);
  
  for(int ix = 0; ix < 22; ++ix){
    for(int iy = 0; iy < 22; ++iy){
      h(ix,iy) = 0.;
    }
  }
  
  LASGlobalData<double> haa =   zpos() * zpos() / err2;
  LASGlobalData<double> hab =   zpos() / err2;
  LASGlobalData<double> hba =   zpos() / err2;
  LASGlobalData<double> hbb =   LASGlobalData<double>(1.) / err2;


  LASGlobalData<double> h00 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) / err2;
  LASGlobalData<double> h01 =   sin(theta()) * cos(theta()) / err2;
  LASGlobalData<double> h02 =   sin(theta()) * cos(theta()) * zpos() / err2;
  LASGlobalData<double> h03 =   sin(theta()) * sin(theta()) * zpos() / err2;
  LASGlobalData<double> h04 =   sin(theta()) / err2;
  LASGlobalData<double> h05 =   sin(theta()) * zpos() / err2;

  LASGlobalData<double> h10 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) / err2;
  LASGlobalData<double> h11 =   cos(theta()) * cos(theta()) / err2;
  LASGlobalData<double> h12 =   cos(theta()) * cos(theta()) * zpos() / err2;
  LASGlobalData<double> h13 =   cos(theta()) * sin(theta()) * zpos() / err2;
  LASGlobalData<double> h14 =   cos(theta()) / err2;
  LASGlobalData<double> h15 =   cos(theta()) * zpos() / err2;

  LASGlobalData<double> h20 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) * zpos() / err2;
  LASGlobalData<double> h21 =   cos(theta()) * cos(theta()) * zpos() / err2;
  LASGlobalData<double> h22 =   cos(theta()) * cos(theta()) * zpos() * zpos() / err2;
  LASGlobalData<double> h23 =   cos(theta()) * sin(theta()) * zpos() * zpos() / err2;
  LASGlobalData<double> h24 =   cos(theta()) * zpos()/ err2;
  LASGlobalData<double> h25 =   cos(theta()) * zpos() * zpos() / err2;

  LASGlobalData<double> h30 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) * zpos() / err2;
  LASGlobalData<double> h31 =   sin(theta()) * cos(theta()) * zpos() / err2;
  LASGlobalData<double> h32 =   sin(theta()) * cos(theta()) * zpos() * zpos() / err2;
  LASGlobalData<double> h33 =   sin(theta()) * sin(theta()) * zpos() * zpos() / err2;
  LASGlobalData<double> h34 =   sin(theta()) * zpos()/ err2;
  LASGlobalData<double> h35 =   sin(theta()) * zpos() * zpos() / err2;

  LASGlobalData<double> h40 = LASGlobalData<double>(0.) - sin(theta()) / r0() / err2;
  LASGlobalData<double> h41 =   cos(theta()) / r0() / err2;
  LASGlobalData<double> h42 =   cos(theta()) * zpos() / r0() / err2;
  LASGlobalData<double> h43 =   sin(theta()) * zpos() / r0() / err2;
  LASGlobalData<double> h44 =   LASGlobalData<double>(1.) / err2;
  LASGlobalData<double> h45 =   zpos() / err2;

  LASGlobalData<double> h50 = LASGlobalData<double>(0.) - sin(theta()) * zpos() / r0() / err2;
  LASGlobalData<double> h51 =   cos(theta()) * zpos() / r0() / err2;
  LASGlobalData<double> h52 =   cos(theta()) * zpos() * zpos() / r0() / err2;
  LASGlobalData<double> h53 =   sin(theta()) * zpos() * zpos() / r0() / err2;
  LASGlobalData<double> h54 =   zpos() / err2;
  LASGlobalData<double> h55 =   zpos() * zpos() / err2;


  for(int bn = 0; bn < 8;  ++bn){
    for(int dn = 0; dn < 4;  ++dn){
      
      if(err2.GetEntry(0,-1,bn,dn) < 1000.){
	    
	h(2*bn,2*bn)     += haa.GetEntry(0,-1,bn,dn);
	h(2*bn,2*bn+1)   += hab.GetEntry(0,-1,bn,dn);
	h(2*bn+1,2*bn)   += hba.GetEntry(0,-1,bn,dn);
	h(2*bn+1,2*bn+1) += hbb.GetEntry(0,-1,bn,dn);
	
	h(npb+0,npb+0)   += h00.GetEntry(0,-1,bn,dn);
	h(npb+0,npb+1)   += h01.GetEntry(0,-1,bn,dn);
	h(npb+0,npb+2)   += h02.GetEntry(0,-1,bn,dn);
	h(npb+0,npb+3)   += h03.GetEntry(0,-1,bn,dn);
	h(npb+0,npb+4)   += h04.GetEntry(0,-1,bn,dn);
	h(npb+0,npb+5)   += h05.GetEntry(0,-1,bn,dn);
	h(npb+0,2*bn)    += h05.GetEntry(0,-1,bn,dn);
	h(npb+0,2*bn+1)  += h04.GetEntry(0,-1,bn,dn);
	
	h(npb+1,npb+0)   += h10.GetEntry(0,-1,bn,dn);
	h(npb+1,npb+1)   += h11.GetEntry(0,-1,bn,dn);
	h(npb+1,npb+2)   += h12.GetEntry(0,-1,bn,dn);
	h(npb+1,npb+3)   += h13.GetEntry(0,-1,bn,dn);
	h(npb+1,npb+4)   += h14.GetEntry(0,-1,bn,dn);
	h(npb+1,npb+5)   += h15.GetEntry(0,-1,bn,dn);
	h(npb+1,2*bn)    += h15.GetEntry(0,-1,bn,dn);
	h(npb+1,2*bn+1)  += h14.GetEntry(0,-1,bn,dn);
	
	h(npb+2,npb+0)   += h20.GetEntry(0,-1,bn,dn);
	h(npb+2,npb+1)   += h21.GetEntry(0,-1,bn,dn);
	h(npb+2,npb+2)   += h22.GetEntry(0,-1,bn,dn);
	h(npb+2,npb+3)   += h23.GetEntry(0,-1,bn,dn);
	h(npb+2,npb+4)   += h24.GetEntry(0,-1,bn,dn);
	h(npb+2,npb+5)   += h25.GetEntry(0,-1,bn,dn);
	h(npb+2,2*bn)    += h25.GetEntry(0,-1,bn,dn);
	h(npb+2,2*bn+1)  += h24.GetEntry(0,-1,bn,dn);
	
	h(npb+3,npb+0)   += h30.GetEntry(0,-1,bn,dn);
	h(npb+3,npb+1)   += h31.GetEntry(0,-1,bn,dn);
	h(npb+3,npb+2)   += h32.GetEntry(0,-1,bn,dn);
	h(npb+3,npb+3)   += h33.GetEntry(0,-1,bn,dn);
	h(npb+3,npb+4)   += h34.GetEntry(0,-1,bn,dn);
	h(npb+3,npb+5)   += h35.GetEntry(0,-1,bn,dn);
	h(npb+3,2*bn)    += h35.GetEntry(0,-1,bn,dn);
	h(npb+3,2*bn+1)  += h34.GetEntry(0,-1,bn,dn);
	
	h(npb+4,npb+0)   += h40.GetEntry(0,-1,bn,dn);
	h(npb+4,npb+1)   += h41.GetEntry(0,-1,bn,dn);
	h(npb+4,npb+2)   += h42.GetEntry(0,-1,bn,dn);
	h(npb+4,npb+3)   += h43.GetEntry(0,-1,bn,dn);
	h(npb+4,npb+4)   += h44.GetEntry(0,-1,bn,dn);
	h(npb+4,npb+5)   += h45.GetEntry(0,-1,bn,dn);
	h(npb+4,2*bn)    += h45.GetEntry(0,-1,bn,dn);
	h(npb+4,2*bn+1)  += h44.GetEntry(0,-1,bn,dn);
	
	h(npb+5,npb+0)   += h50.GetEntry(0,-1,bn,dn);
	h(npb+5,npb+1)   += h51.GetEntry(0,-1,bn,dn);
	h(npb+5,npb+2)   += h52.GetEntry(0,-1,bn,dn);
	h(npb+5,npb+3)   += h53.GetEntry(0,-1,bn,dn);
	h(npb+5,npb+4)   += h54.GetEntry(0,-1,bn,dn);
	h(npb+5,npb+5)   += h55.GetEntry(0,-1,bn,dn);
	h(npb+5,2*bn)    += h55.GetEntry(0,-1,bn,dn);
	h(npb+5,2*bn+1)  += h54.GetEntry(0,-1,bn,dn);
	
	
	h(2*bn,npb+0)   += h50.GetEntry(0,-1,bn,dn);
	h(2*bn,npb+1)   += h51.GetEntry(0,-1,bn,dn);
	h(2*bn,npb+2)   += h52.GetEntry(0,-1,bn,dn);
	h(2*bn,npb+3)   += h53.GetEntry(0,-1,bn,dn);
	h(2*bn,npb+4)   += h54.GetEntry(0,-1,bn,dn);
	h(2*bn,npb+5)   += h55.GetEntry(0,-1,bn,dn);
	
	h(2*bn+1,npb+0) += h40.GetEntry(0,-1,bn,dn);
	h(2*bn+1,npb+1) += h41.GetEntry(0,-1,bn,dn);
	h(2*bn+1,npb+2) += h42.GetEntry(0,-1,bn,dn);
	h(2*bn+1,npb+3) += h43.GetEntry(0,-1,bn,dn);
	h(2*bn+1,npb+4) += h44.GetEntry(0,-1,bn,dn);
	h(2*bn+1,npb+5) += h45.GetEntry(0,-1,bn,dn);
      }
	
      if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
      }
    }
  }  

  
  h.Invert();

  double cp[22][2][8][6];
  double up[22][22];
  double ps[22];

  for(int ip = 0; ip < 22; ++ip){
    for(int ib = 0; ib < 2; ++ib){
      for(int bn = 0; bn < 8;  ++bn){
	for(int dn = 0; dn < 4;  ++dn){
	  cp[ip][ib][bn][dn] = 0.;
	  for(int ix = 0; ix < 22; ++ix){
	    cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
	    }
	}
      }
    }
  }

  for(int ip = 0; ip < 22; ++ip){
    for(int ix = 0; ix < 22; ++ix){
      up[ip][ix] = 0.;
      for(int ib = 0; ib < 2; ++ib){
	for(int bn = 0; bn < 8;  ++bn){
	  for(int dn = 0; dn < 4;  ++dn){
	    if( ib == 0){
	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(0,-1,bn,dn);
	    }
	    if( ib == 1 ){
	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
	    }
	  }
	}
      }
    }
  }
 
  for(int bn = 0; bn < 8; ++bn){
    results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
    results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
  }

  results.er_tecp_dx = sqrt(up[npb+0][npb+0]);
  results.er_tecp_dy = sqrt(up[npb+1][npb+1]);
  results.er_tecp_rx = sqrt(up[npb+2][npb+2]);
  results.er_tecp_ry = sqrt(up[npb+3][npb+3]);
  results.er_tecp_rz = sqrt(up[npb+4][npb+4]);
  results.er_tecp_tz = sqrt(up[npb+5][npb+5]);
  results.tecp_dxdy  = up[npb+0][npb+1]/results.er_tecp_dx/results.er_tecp_dy;
  results.tecp_dxrx  = up[npb+0][npb+2]/results.er_tecp_dx/results.er_tecp_rx;
  results.tecp_dxry  = up[npb+0][npb+3]/results.er_tecp_dx/results.er_tecp_ry;
  results.tecp_dxrz  = up[npb+0][npb+4]/results.er_tecp_dx/results.er_tecp_rz;
  results.tecp_dxtz  = up[npb+0][npb+5]/results.er_tecp_dx/results.er_tecp_tz;
  results.tecp_dyrx  = up[npb+1][npb+2]/results.er_tecp_dy/results.er_tecp_rx;
  results.tecp_dyry  = up[npb+1][npb+3]/results.er_tecp_dy/results.er_tecp_ry;
  results.tecp_dyrz  = up[npb+1][npb+4]/results.er_tecp_dy/results.er_tecp_rz;
  results.tecp_dytz  = up[npb+1][npb+5]/results.er_tecp_dy/results.er_tecp_tz;
  results.tecp_rxry  = up[npb+2][npb+3]/results.er_tecp_rx/results.er_tecp_ry;
  results.tecp_rxrz  = up[npb+2][npb+4]/results.er_tecp_rx/results.er_tecp_rz;
  results.tecp_rxtz  = up[npb+2][npb+5]/results.er_tecp_rx/results.er_tecp_tz;
  results.tecp_ryrz  = up[npb+3][npb+4]/results.er_tecp_ry/results.er_tecp_rz;
  results.tecp_rytz  = up[npb+3][npb+5]/results.er_tecp_ry/results.er_tecp_tz;
  results.tecp_rztz  = up[npb+4][npb+5]/results.er_tecp_rz/results.er_tecp_tz;

  for(int ip = 0; ip < 22; ++ip){
    ps[ip] = 0.;
    for(int ix = 0; ix < 22; ++ix){
      ps[ip] += h(ip,ix)*bs[ix];
    }
  }

  for(int bn = 0; bn < 8; ++bn){
    results.beam_a[bn] = ps[2*bn];
    results.beam_b[bn] = ps[2*bn+1];
  }

  results.tecp_dx = ps[npb+0];
  results.tecp_dy = ps[npb+1];
  results.tecp_rx = ps[npb+2];
  results.tecp_ry = ps[npb+3];
  results.tecp_rz = ps[npb+4];
  results.tecp_tz = ps[npb+5];

}
  

// reconstruction of TEC- wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TECM2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results)
{

  // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

  double bs[22];
  double db[22][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TEC-, TOB); k=0..7 (beam number); l=0..5 (module number)

  int npb = 16;// number of beam parameters
 
  LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
  LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
  LASGlobalData<double> par2 = dif * cos(theta()) * r0() * zpos()/ err2;
  LASGlobalData<double> par3 = dif * sin(theta()) * r0() * zpos()/ err2;
  LASGlobalData<double> par4 = dif / err2;
  LASGlobalData<double> par5 = dif * zpos()/ err2;


  LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
  LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
  LASGlobalData<double> dev2 = cos(theta()) * r0() * zpos()/ err2;
  LASGlobalData<double> dev3 = sin(theta()) * r0() * zpos()/ err2;
  LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
  LASGlobalData<double> dev5 = zpos()/ err2;

  for(int ip = 0; ip < 22;  ++ip){ 
    bs[ip] = 0.;
  }
  for(int bn = 0; bn < 8;  ++bn){
    for(int dn = 0; dn < 4;  ++dn){
      for(int ip = 0; ip < 22;  ++ip){ 
	db[ip][0][bn][dn] = 0.;
	db[ip][1][bn][dn] = 0.;
      }
      if(err2.GetEntry(1,-1,bn,dn) < 1000.){
	    
	bs[2*bn] += par5.GetEntry(1,-1,bn,dn);
	bs[2*bn+1] += par4.GetEntry(1,-1,bn,dn);
	
	bs[npb+0] += par0.GetEntry(1,-1,bn,dn);
	bs[npb+1] += par1.GetEntry(1,-1,bn,dn);
	bs[npb+2] += par2.GetEntry(1,-1,bn,dn);
	bs[npb+3] += par3.GetEntry(1,-1,bn,dn);
	bs[npb+4] += par4.GetEntry(1,-1,bn,dn);
	bs[npb+5] += par5.GetEntry(1,-1,bn,dn);
        
	db[2*bn][0][bn][dn] = dev5.GetEntry(1,-1,bn,dn);
	db[2*bn+1][0][bn][dn] = dev4.GetEntry(1,-1,bn,dn);
	
	db[npb+0][0][bn][dn] = dev0.GetEntry(1,-1,bn,dn);
	db[npb+1][0][bn][dn] = dev1.GetEntry(1,-1,bn,dn);
	db[npb+2][0][bn][dn] = dev2.GetEntry(1,-1,bn,dn);
	db[npb+3][0][bn][dn] = dev3.GetEntry(1,-1,bn,dn);
	db[npb+4][0][bn][dn] = dev4.GetEntry(1,-1,bn,dn);
	db[npb+5][0][bn][dn] = dev5.GetEntry(1,-1,bn,dn);
      }
      if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);
      }
    }
  }
   
  TMatrixT<double> h(22,22);
  
  for(int ix = 0; ix < 22; ++ix){
    for(int iy = 0; iy < 22; ++iy){
      h(ix,iy) = 0.;
    }
  }
  
  LASGlobalData<double> haa =   zpos() * zpos() / err2;
  LASGlobalData<double> hab =   zpos() / err2;
  LASGlobalData<double> hba =   zpos() / err2;
  LASGlobalData<double> hbb =   LASGlobalData<double>(1.) / err2;


  LASGlobalData<double> h00 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) / err2;
  LASGlobalData<double> h01 =   sin(theta()) * cos(theta()) / err2;
  LASGlobalData<double> h02 =   sin(theta()) * cos(theta()) * zpos() / err2;
  LASGlobalData<double> h03 =   sin(theta()) * sin(theta()) * zpos() / err2;
  LASGlobalData<double> h04 =   sin(theta()) / err2;
  LASGlobalData<double> h05 =   sin(theta()) * zpos() / err2;

  LASGlobalData<double> h10 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) / err2;
  LASGlobalData<double> h11 =   cos(theta()) * cos(theta()) / err2;
  LASGlobalData<double> h12 =   cos(theta()) * cos(theta()) * zpos() / err2;
  LASGlobalData<double> h13 =   cos(theta()) * sin(theta()) * zpos() / err2;
  LASGlobalData<double> h14 =   cos(theta()) / err2;
  LASGlobalData<double> h15 =   cos(theta()) * zpos() / err2;

  LASGlobalData<double> h20 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) * zpos() / err2;
  LASGlobalData<double> h21 =   cos(theta()) * cos(theta()) * zpos() / err2;
  LASGlobalData<double> h22 =   cos(theta()) * cos(theta()) * zpos() * zpos() / err2;
  LASGlobalData<double> h23 =   cos(theta()) * sin(theta()) * zpos() * zpos() / err2;
  LASGlobalData<double> h24 =   cos(theta()) * zpos()/ err2;
  LASGlobalData<double> h25 =   cos(theta()) * zpos() * zpos() / err2;

  LASGlobalData<double> h30 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) * zpos() / err2;
  LASGlobalData<double> h31 =   sin(theta()) * cos(theta()) * zpos() / err2;
  LASGlobalData<double> h32 =   sin(theta()) * cos(theta()) * zpos() * zpos() / err2;
  LASGlobalData<double> h33 =   sin(theta()) * sin(theta()) * zpos() * zpos() / err2;
  LASGlobalData<double> h34 =   sin(theta()) * zpos()/ err2;
  LASGlobalData<double> h35 =   sin(theta()) * zpos() * zpos() / err2;

  LASGlobalData<double> h40 = LASGlobalData<double>(0.) - sin(theta()) / r0() / err2;
  LASGlobalData<double> h41 =   cos(theta()) / r0() / err2;
  LASGlobalData<double> h42 =   cos(theta()) * zpos() / r0() / err2;
  LASGlobalData<double> h43 =   sin(theta()) * zpos() / r0() / err2;
  LASGlobalData<double> h44 =   LASGlobalData<double>(1.) / err2;
  LASGlobalData<double> h45 =   zpos() / err2;

  LASGlobalData<double> h50 = LASGlobalData<double>(0.) - sin(theta()) * zpos() / r0() / err2;
  LASGlobalData<double> h51 =   cos(theta()) * zpos() / r0() / err2;
  LASGlobalData<double> h52 =   cos(theta()) * zpos() * zpos() / r0() / err2;
  LASGlobalData<double> h53 =   sin(theta()) * zpos() * zpos() / r0() / err2;
  LASGlobalData<double> h54 =   zpos() / err2;
  LASGlobalData<double> h55 =   zpos() * zpos() / err2;


  for(int bn = 0; bn < 8;  ++bn){
    for(int dn = 0; dn < 4;  ++dn){
      
      if(err2.GetEntry(1,-1,bn,dn) < 1000.){
	    
	h(2*bn,2*bn)     += haa.GetEntry(1,-1,bn,dn);
	h(2*bn,2*bn+1)   += hab.GetEntry(1,-1,bn,dn);
	h(2*bn+1,2*bn)   += hba.GetEntry(1,-1,bn,dn);
	h(2*bn+1,2*bn+1) += hbb.GetEntry(1,-1,bn,dn);
	
	h(npb+0,npb+0)   += h00.GetEntry(1,-1,bn,dn);
	h(npb+0,npb+1)   += h01.GetEntry(1,-1,bn,dn);
	h(npb+0,npb+2)   += h02.GetEntry(1,-1,bn,dn);
	h(npb+0,npb+3)   += h03.GetEntry(1,-1,bn,dn);
	h(npb+0,npb+4)   += h04.GetEntry(1,-1,bn,dn);
	h(npb+0,npb+5)   += h05.GetEntry(1,-1,bn,dn);
	h(npb+0,2*bn)    += h05.GetEntry(1,-1,bn,dn);
	h(npb+0,2*bn+1)  += h04.GetEntry(1,-1,bn,dn);
	
	h(npb+1,npb+0)   += h10.GetEntry(1,-1,bn,dn);
	h(npb+1,npb+1)   += h11.GetEntry(1,-1,bn,dn);
	h(npb+1,npb+2)   += h12.GetEntry(1,-1,bn,dn);
	h(npb+1,npb+3)   += h13.GetEntry(1,-1,bn,dn);
	h(npb+1,npb+4)   += h14.GetEntry(1,-1,bn,dn);
	h(npb+1,npb+5)   += h15.GetEntry(1,-1,bn,dn);
	h(npb+1,2*bn)    += h15.GetEntry(1,-1,bn,dn);
	h(npb+1,2*bn+1)  += h14.GetEntry(1,-1,bn,dn);
	
	h(npb+2,npb+0)   += h20.GetEntry(1,-1,bn,dn);
	h(npb+2,npb+1)   += h21.GetEntry(1,-1,bn,dn);
	h(npb+2,npb+2)   += h22.GetEntry(1,-1,bn,dn);
	h(npb+2,npb+3)   += h23.GetEntry(1,-1,bn,dn);
	h(npb+2,npb+4)   += h24.GetEntry(1,-1,bn,dn);
	h(npb+2,npb+5)   += h25.GetEntry(1,-1,bn,dn);
	h(npb+2,2*bn)    += h25.GetEntry(1,-1,bn,dn);
	h(npb+2,2*bn+1)  += h24.GetEntry(1,-1,bn,dn);
	
	h(npb+3,npb+0)   += h30.GetEntry(1,-1,bn,dn);
	h(npb+3,npb+1)   += h31.GetEntry(1,-1,bn,dn);
	h(npb+3,npb+2)   += h32.GetEntry(1,-1,bn,dn);
	h(npb+3,npb+3)   += h33.GetEntry(1,-1,bn,dn);
	h(npb+3,npb+4)   += h34.GetEntry(1,-1,bn,dn);
	h(npb+3,npb+5)   += h35.GetEntry(1,-1,bn,dn);
	h(npb+3,2*bn)    += h35.GetEntry(1,-1,bn,dn);
	h(npb+3,2*bn+1)  += h34.GetEntry(1,-1,bn,dn);
	
	h(npb+4,npb+0)   += h40.GetEntry(1,-1,bn,dn);
	h(npb+4,npb+1)   += h41.GetEntry(1,-1,bn,dn);
	h(npb+4,npb+2)   += h42.GetEntry(1,-1,bn,dn);
	h(npb+4,npb+3)   += h43.GetEntry(1,-1,bn,dn);
	h(npb+4,npb+4)   += h44.GetEntry(1,-1,bn,dn);
	h(npb+4,npb+5)   += h45.GetEntry(1,-1,bn,dn);
	h(npb+4,2*bn)    += h45.GetEntry(1,-1,bn,dn);
	h(npb+4,2*bn+1)  += h44.GetEntry(1,-1,bn,dn);
	
	h(npb+5,npb+0)   += h50.GetEntry(1,-1,bn,dn);
	h(npb+5,npb+1)   += h51.GetEntry(1,-1,bn,dn);
	h(npb+5,npb+2)   += h52.GetEntry(1,-1,bn,dn);
	h(npb+5,npb+3)   += h53.GetEntry(1,-1,bn,dn);
	h(npb+5,npb+4)   += h54.GetEntry(1,-1,bn,dn);
	h(npb+5,npb+5)   += h55.GetEntry(1,-1,bn,dn);
	h(npb+5,2*bn)    += h55.GetEntry(1,-1,bn,dn);
	h(npb+5,2*bn+1)  += h54.GetEntry(1,-1,bn,dn);
	
	
	h(2*bn,npb+0)   += h50.GetEntry(1,-1,bn,dn);
	h(2*bn,npb+1)   += h51.GetEntry(1,-1,bn,dn);
	h(2*bn,npb+2)   += h52.GetEntry(1,-1,bn,dn);
	h(2*bn,npb+3)   += h53.GetEntry(1,-1,bn,dn);
	h(2*bn,npb+4)   += h54.GetEntry(1,-1,bn,dn);
	h(2*bn,npb+5)   += h55.GetEntry(1,-1,bn,dn);
	
	h(2*bn+1,npb+0) += h40.GetEntry(1,-1,bn,dn);
	h(2*bn+1,npb+1) += h41.GetEntry(1,-1,bn,dn);
	h(2*bn+1,npb+2) += h42.GetEntry(1,-1,bn,dn);
	h(2*bn+1,npb+3) += h43.GetEntry(1,-1,bn,dn);
	h(2*bn+1,npb+4) += h44.GetEntry(1,-1,bn,dn);
	h(2*bn+1,npb+5) += h45.GetEntry(1,-1,bn,dn);
      }
	
      if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
      }
    }
  }  

  
  h.Invert();

  double cp[22][2][8][6];
  double up[22][22];
  double ps[22];

  for(int ip = 0; ip < 22; ++ip){
    for(int ib = 0; ib < 2; ++ib){
      for(int bn = 0; bn < 8;  ++bn){
	for(int dn = 0; dn < 4;  ++dn){
	  cp[ip][ib][bn][dn] = 0.;
	  for(int ix = 0; ix < 22; ++ix){
	    cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
	    }
	}
      }
    }
  }

  for(int ip = 0; ip < 22; ++ip){
    for(int ix = 0; ix < 22; ++ix){
      up[ip][ix] = 0.;
      for(int ib = 0; ib < 2; ++ib){
	for(int bn = 0; bn < 8;  ++bn){
	  for(int dn = 0; dn < 4;  ++dn){
	    if( ib == 0){
	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(1,-1,bn,dn);
	    }
	    if( ib == 1 ){
	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
	    }
	  }
	}
      }
    }
  }
 
  for(int bn = 0; bn < 8; ++bn){
    results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
    results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
  }

  results.er_tecm_dx = sqrt(up[npb+0][npb+0]);
  results.er_tecm_dy = sqrt(up[npb+1][npb+1]);
  results.er_tecm_rx = sqrt(up[npb+2][npb+2]);
  results.er_tecm_ry = sqrt(up[npb+3][npb+3]);
  results.er_tecm_rz = sqrt(up[npb+4][npb+4]);
  results.er_tecm_tz = sqrt(up[npb+5][npb+5]);
  results.tecm_dxdy  = up[npb+0][npb+1]/results.er_tecm_dx/results.er_tecm_dy;
  results.tecm_dxrx  = up[npb+0][npb+2]/results.er_tecm_dx/results.er_tecm_rx;
  results.tecm_dxry  = up[npb+0][npb+3]/results.er_tecm_dx/results.er_tecm_ry;
  results.tecm_dxrz  = up[npb+0][npb+4]/results.er_tecm_dx/results.er_tecm_rz;
  results.tecm_dxtz  = up[npb+0][npb+5]/results.er_tecm_dx/results.er_tecm_tz;
  results.tecm_dyrx  = up[npb+1][npb+2]/results.er_tecm_dy/results.er_tecm_rx;
  results.tecm_dyry  = up[npb+1][npb+3]/results.er_tecm_dy/results.er_tecm_ry;
  results.tecm_dyrz  = up[npb+1][npb+4]/results.er_tecm_dy/results.er_tecm_rz;
  results.tecm_dytz  = up[npb+1][npb+5]/results.er_tecm_dy/results.er_tecm_tz;
  results.tecm_rxry  = up[npb+2][npb+3]/results.er_tecm_rx/results.er_tecm_ry;
  results.tecm_rxrz  = up[npb+2][npb+4]/results.er_tecm_rx/results.er_tecm_rz;
  results.tecm_rxtz  = up[npb+2][npb+5]/results.er_tecm_rx/results.er_tecm_tz;
  results.tecm_ryrz  = up[npb+3][npb+4]/results.er_tecm_ry/results.er_tecm_rz;
  results.tecm_rytz  = up[npb+3][npb+5]/results.er_tecm_ry/results.er_tecm_tz;
  results.tecm_rztz  = up[npb+4][npb+5]/results.er_tecm_rz/results.er_tecm_tz;

  for(int ip = 0; ip < 22; ++ip){
    ps[ip] = 0.;
    for(int ix = 0; ix < 22; ++ix){
      ps[ip] += h(ip,ix)*bs[ix];
    }
  }

  for(int bn = 0; bn < 8; ++bn){
    results.beam_a[bn] = ps[2*bn];
    results.beam_b[bn] = ps[2*bn+1];
  }

  results.tecm_dx = ps[npb+0];
  results.tecm_dy = ps[npb+1];
  results.tecm_rx = ps[npb+2];
  results.tecm_ry = ps[npb+3];
  results.tecm_rz = ps[npb+4];
  results.tecm_tz = ps[npb+5];

}
  
