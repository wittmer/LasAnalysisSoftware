#include "LAS_Tec_Reconstruction.h"

#include "TFile.h"
#include "TMath.h"
#include "LAS_globaldata_tools.h"
#include "Avec2D.h"

int count_masked(LASGlobalData<int>& mask, LAS::beam_group bg)
{
  int ctr = 0;

  LASGlobalDataLoop loop(convert_to_loop_type(bg));
  do{
    if(loop.GetEntry<int>(mask) != 1) ctr++;
  }while(loop.next());
  return ctr;
}

void TEC_Parameter_Reconstruction(const std::string& result_file, const LASGlobalData<double>& ref_pos, const LASGlobalData<int>& ref_mask, const std::string& dir_name)
{
  std::cout << "Starting TEC_Parameter_REconstruction" << std::endl;

  bool process_rings = false;
  const double rms_cut = 2.0;

  TECReconstructor TR;
  TR.Init_common();

  TFile rf(result_file.c_str(),"UPDATE");
  Avec timev = avec_get("block_timestamp", rf);
  Avec good_block = avec_get("good_block", rf);
  if (timev.empty()){
    std::cerr << "Could not find timestamps in file " << result_file << std::endl;
    return;
  }
  if (good_block.empty()){
    std::cerr << "Could not find Avec good_block in file " << result_file << std::endl;
    return;
  }
  if(good_block.size() != timev.size()){
    std::cerr << "Error, good_block.size() = " << good_block.size() << "  timev.size() = " << timev.size() << "   mismatch" << std::endl;
    return;
  }
  Avec good_fit(good_block.size(), 0.0);
  Avec ssq_tecp_before(good_block.size(), 0.0);
  Avec ssq_tecm_before(good_block.size(), 0.0);
  Avec ssq_tecp_after(good_block.size(), 0.0);
  Avec ssq_tecm_after(good_block.size(), 0.0);
  Avec nr_iter_tecp(good_block.size(), -1.0);
  Avec nr_iter_tecm(good_block.size(), -1.0);
  Avec nr_masked_tecp(good_block.size(), 0.0);
  Avec nr_masked_tecm(good_block.size(), 0.0);

  rf.mkdir(dir_name.c_str());
  for(Avec::size_type f=0; f < timev.size(); f++){
    //for(Avec::size_type f=0; f < 1; f++){
    if(good_block[f] != 1)continue;
    std::ostringstream pos_name;
    std::ostringstream rms_name;
    std::ostringstream mask_name;
    pos_name << "positions_" << f;
    rms_name << "rms_av_" << f;
    mask_name << "positions_mask_" << f;
    
    LASGlobalData<double> diff_rad = (ref_pos - global_data_get<double>(pos_name.str(),rf)) * pitch() / r0();
    correct_signs(diff_rad);
    LASGlobalData<double> diff_mm = (ref_pos - global_data_get<double>(pos_name.str(),rf)) * pitch();
    correct_signs(diff_mm);
    //print_global_data(diff_mm);
    //print_global_data(ref_mask);
    LASGlobalData<double> rms = global_data_get<double>(rms_name.str(),rf);
    LASGlobalData<int> mask = global_data_get<int>(mask_name.str(), rf) && ref_mask && (rms < rms_cut);

    nr_masked_tecp[f] = count_masked(mask, LAS::TEC_PLUS);
    nr_masked_tecm[f] = count_masked(mask, LAS::TEC_MINUS);
    if(nr_masked_tecp[f] > 10 || nr_masked_tecm[f] > 10)continue;

    LasAlPar alpar;
    if(process_rings){
      TR.Init_onering();
      TR.reconstruct_tec_ring(diff_rad, mask, alpar.tecp_r4, LAS::TEC_PLUS_R4);
      TR.reconstruct_tec_ring(diff_rad, mask, alpar.tecp_r6, LAS::TEC_PLUS_R6);
      TR.reconstruct_tec_ring(diff_rad, mask, alpar.tecm_r4, LAS::TEC_MINUS_R4);
      TR.reconstruct_tec_ring(diff_rad, mask, alpar.tecm_r6, LAS::TEC_MINUS_R6);
    }
    
    TR.Init_full();
    get_sum_sq(diff_mm, mask, ssq_tecp_before[f], ssq_tecm_before[f]);

    if(TR.reconstruct_tec_full(diff_mm, mask, alpar)){
      good_fit[f] = 1;
      LASGlobalData<double> residuals;
      TR.compute_residuals(diff_mm, alpar, residuals);

      std::ostringstream alpar_name;
      alpar_name << "alpar_" << f;
      std::ostringstream residual_name;
      residual_name << "residual_" << f;

      rf.cd(dir_name.c_str());
      alpar.Write(alpar_name.str().c_str(), TObject::kWriteDelete);
      residuals.Write(residual_name.str().c_str(), TObject::kWriteDelete);
      rf.cd();

      get_sum_sq(diff_mm - TR.dx_rec_tec(alpar), mask, ssq_tecp_after[f], ssq_tecm_after[f]);
      nr_iter_tecp[f] = TR.nr_iter_tecp;
      nr_iter_tecm[f] = TR.nr_iter_tecm;
    }
  }
  Int_t wflag = TObject::kWriteDelete;
  rf.cd(dir_name.c_str());
  good_fit.Write("good_fit", wflag);
  ssq_tecp_before.Write("ssq_tecp_before", wflag);
  ssq_tecp_after.Write("ssq_tecp_after", wflag);
  ssq_tecm_before.Write("ssq_tecm_before", wflag);
  ssq_tecm_after.Write("ssq_tecm_after", wflag);
  nr_iter_tecp.Write("nr_iter_tecp", wflag);
  nr_iter_tecm.Write("nr_iter_tecm", wflag);
  nr_masked_tecp.Write("nr_masked_tecp", wflag);
  nr_masked_tecm.Write("nr_masked_tecm", wflag);
  rf.Close();
}

void check_alpar(const TecPar& alpar)
{
  TECReconstructor TR;
  TR.Init_common();
  TR.Init_full();

  Avec zdat = TECReconstructor::make_zdat();
  std::cout << "vsum(Dphik): " << vsum(alpar.Dphik) << std::endl;
  std::cout << "vsum(Dxk): " << vsum(alpar.Dxk) << std::endl;
  std::cout << "vsum(Dyk): " << vsum(alpar.Dyk) << std::endl;
  std::cout << "vsum(Dphik*zdat): " << vsum(alpar.Dphik*zdat) << std::endl;
  std::cout << "vsum(Dxk*zdat): " << vsum(alpar.Dxk*zdat) << std::endl;
  std::cout << "vsum(Dyk*zdat): " << vsum(alpar.Dyk*zdat) << std::endl;

  std::cout << "vsum(DthetaA_r4): " << vsum(alpar.DthetaA_r4) << std::endl;
  std::cout << "vsum(DthetaA_r6): " << vsum(alpar.DthetaA_r6) << std::endl;
  std::cout << "vsum(DthetaB_r4): " << vsum(alpar.DthetaB_r4) << std::endl;
  std::cout << "vsum(DthetaB_r6): " << vsum(alpar.DthetaB_r6) << std::endl;

//   std::cout << "vsum(DthetaA_r4)*TR.r4: " << vsum(alpar.DthetaA_r4)*TR.r4 << std::endl;
//   std::cout << "vsum(DthetaA_r6)*TR.r6: " << vsum(alpar.DthetaA_r6)*TR.r6 << std::endl;
//   std::cout << "vsum(DthetaA_r4)*TR.r4*TR.r4: " << vsum(alpar.DthetaA_r4)*TR.r4*TR.r4 << std::endl;
//   std::cout << "vsum(DthetaA_r6)*TR.r6*TR.r6: " << vsum(alpar.DthetaA_r6)*TR.r6*TR.r6 << std::endl;

  std::cout << "vsum(DthetaA_r4)*TR.r4*TR.r4 + vsum(DthetaA_r6)*TR.r6*TR.r6: " << vsum(alpar.DthetaA_r4)*TR.r4*TR.r4 + vsum(alpar.DthetaA_r6)*TR.r6*TR.r6 << std::endl;
  std::cout << "vsum(DthetaB_r4)*TR.r4*TR.r4 + vsum(DthetaB_r6)*TR.r6*TR.r6: " << vsum(alpar.DthetaB_r4)*TR.r4*TR.r4 + vsum(alpar.DthetaB_r6)*TR.r6*TR.r6 << std::endl;
//   std::cout << "vsum(DthetaA_r6)*TR.r6*TR.r6: " << vsum(alpar.DthetaA_r6)*TR.r6*TR.r6 << std::endl;

//   std::cout << "vsum(DthetaB_r4)*TR.r4: " << vsum(alpar.DthetaB_r4)*TR.r4 << std::endl;
//   std::cout << "vsum(DthetaB_r6)*TR.r6: " << vsum(alpar.DthetaB_r6)*TR.r6 << std::endl;
//   std::cout << "vsum(DthetaB_r4)*TR.r4*TR.r4: " << vsum(alpar.DthetaB_r4)*TR.r4*TR.r4 << std::endl;
//   std::cout << "vsum(DthetaB_r6)*TR.r6*TR.r6: " << vsum(alpar.DthetaB_r6)*TR.r6*TR.r6 << std::endl;


}

void get_sum_sq(const LASGlobalData<double>& data, const LASGlobalData<int>& mask, double& ssq_tecp, double& ssq_tecm)
{
  ssq_tecp = 0;
  int norm_tecp = 0;

  double ssq_tecp_r4 = 0;
  int norm_tecp_r4 = 0;

  double ssq_tecp_r6 = 0;
  int norm_tecp_r6 = 0;

  ssq_tecm = 0;
  int norm_tecm = 0;

  double ssq_tecm_r4 = 0;
  int norm_tecm_r4 = 0;

  double ssq_tecm_r6 = 0;
  int norm_tecm_r6 = 0;

  LASGlobalDataLoop l1(LASGlobalDataLoop::TEC_PLUS_R4);
  do{
    if(l1.GetEntry<int>(mask)==0)continue;
    ssq_tecp_r4 += pow(l1.GetEntry<double>(data), 2);
    norm_tecp_r4 ++;
    ssq_tecp += pow(l1.GetEntry<double>(data), 2);
    norm_tecp ++;
  }while(l1.next());

  LASGlobalDataLoop l2(LASGlobalDataLoop::TEC_PLUS_R6);
  do{
    if(l2.GetEntry<int>(mask)==0)continue;
    ssq_tecp_r6 += pow(l2.GetEntry<double>(data), 2);
    norm_tecp_r6 ++;
    ssq_tecp += pow(l2.GetEntry<double>(data), 2);
    norm_tecp ++;
  }while(l2.next());

  LASGlobalDataLoop l3(LASGlobalDataLoop::TEC_MINUS_R4);
  do{
    if(l3.GetEntry<int>(mask)==0)continue;
    ssq_tecm_r4 += pow(l3.GetEntry<double>(data), 2);
    norm_tecm_r4 ++;
    ssq_tecm += pow(l3.GetEntry<double>(data), 2);
    norm_tecm ++;
  }while(l3.next());

  LASGlobalDataLoop l4(LASGlobalDataLoop::TEC_MINUS_R6);
  do{
    if(l4.GetEntry<int>(mask)==0)continue;
    ssq_tecm_r6 += pow(l4.GetEntry<double>(data), 2);
    norm_tecm_r6 ++;
    ssq_tecm += pow(l4.GetEntry<double>(data), 2);
    norm_tecm ++;
  }while(l4.next());

  ssq_tecp /= norm_tecp;
  ssq_tecm /= norm_tecm;
}

void print_sum_sq(const LASGlobalData<double>& data, const LASGlobalData<int>& mask)
{
  double ssq_tecp = 0;
  int norm_tecp = 0;

  double ssq_tecp_r4 = 0;
  int norm_tecp_r4 = 0;

  double ssq_tecp_r6 = 0;
  int norm_tecp_r6 = 0;

  double ssq_tecm = 0;
  int norm_tecm = 0;

  double ssq_tecm_r4 = 0;
  int norm_tecm_r4 = 0;

  double ssq_tecm_r6 = 0;
  int norm_tecm_r6 = 0;

  LASGlobalDataLoop l1(LASGlobalDataLoop::TEC_PLUS_R4);
  do{
    if(l1.GetEntry<int>(mask)==0)continue;
    ssq_tecp_r4 += pow(l1.GetEntry<double>(data), 2);
    norm_tecp_r4 ++;
    ssq_tecp += pow(l1.GetEntry<double>(data), 2);
    norm_tecp ++;
  }while(l1.next());

  LASGlobalDataLoop l2(LASGlobalDataLoop::TEC_PLUS_R6);
  do{
    if(l2.GetEntry<int>(mask)==0)continue;
    ssq_tecp_r6 += pow(l2.GetEntry<double>(data), 2);
    norm_tecp_r6 ++;
    ssq_tecp += pow(l2.GetEntry<double>(data), 2);
    norm_tecp ++;
  }while(l2.next());

  LASGlobalDataLoop l3(LASGlobalDataLoop::TEC_MINUS_R4);
  do{
    if(l3.GetEntry<int>(mask)==0)continue;
    ssq_tecm_r4 += pow(l3.GetEntry<double>(data), 2);
    norm_tecm_r4 ++;
    ssq_tecm += pow(l3.GetEntry<double>(data), 2);
    norm_tecm ++;
  }while(l3.next());

  LASGlobalDataLoop l4(LASGlobalDataLoop::TEC_MINUS_R6);
  do{
    if(l4.GetEntry<int>(mask)==0)continue;
    ssq_tecm_r6 += pow(l4.GetEntry<double>(data), 2);
    norm_tecm_r6 ++;
    ssq_tecm += pow(l4.GetEntry<double>(data), 2);
    norm_tecm ++;
  }while(l4.next());

  std::cout << "ssq TEC+ R4: " << sqrt(ssq_tecp_r4/norm_tecp_r4) << std::endl;
  std::cout << "ssq TEC+ R6: " << sqrt(ssq_tecp_r6/norm_tecp_r6) << std::endl;
  std::cout << "ssq TEC- R4: " << sqrt(ssq_tecm_r4/norm_tecm_r4) << std::endl;
  std::cout << "ssq TEC- R6: " << sqrt(ssq_tecm_r6/norm_tecm_r6) << std::endl;
 
  std::cout << "ssq TEC+: " << sqrt(ssq_tecp/norm_tecp) << std::endl;
  std::cout << "ssq TEC-: " << sqrt(ssq_tecm/norm_tecm) << std::endl;

  std::cout << std::endl;
}


//! Invert the sign for TOB modules with inverese strip orientation
void correct_signs(LASGlobalData<double>& data)
{

  LASGlobalDataLoop loop(LASGlobalDataLoop::TOB);
  do{
    double sign = 1.0;

    switch( loop.get_zpos() ){
    case 0:
    case 3:
    case 4:
      sign *= -1;
    }
    loop.GetEntry<double>(data) *= sign;
  }while(loop.next());

  loop = LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS);
  do{
    loop.GetEntry<double>(data) *= -1.0;
  }while(loop.next());

  //std::cout << "Not Inverting TEC- AT" << std::endl;
  loop = LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS_AT);
  do{
    loop.GetEntry<double>(data) *= -1.0;
  }while(loop.next());

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TECReconstructor::TECReconstructor():
  r4(LAS::r_tec_r4),
  r6(LAS::r_tec_r6),
  nbeam(8),
  ndisc(9),
  nring(2),
  pitch_r4(0.126029),
  pitch_r6(0.189),
  theta( Avec(nbeam) ),
  zdat( Avec(ndisc) ),
  //beam_mask( Avec(nbeam, 1) ),
  //beam_mask_r4( Avec(nbeam, 1) ),
  //beam_mask_r6( Avec(nbeam, 1) ),
  //disc_mask( Avec(ndisc, 1) ),
  L(),
  Sz(0),
  Sz2(0),
  St(0),
  St2(0),
  Ssin(0),
  Scos(0),
  Ssin2(0),
  Scos2(0),
  Scossin(0),
  A(0),
  B(0),
  C(0),
  D(0),
  E(0),
  F(0),
  denom1(0),
  term1(0),
  term2(0),
  nenn2(0),
  znorm(0),
  verbose_flag(1)
{
  Init_common();
}



//! TEC reconstruction

//! phi-position of TEC beams
Avec TECReconstructor::make_theta()
{
  Avec v1(9,0,TMath::Pi()*2);
  return v1(0,7)+0.3927;
}


void TECReconstructor::Init_common()
{
  theta = TECReconstructor::make_theta(); // Phi-positions of beams
}

// void TECReconstructor::Init_beam_mask(int b1, int b2, int b3, int b4, int b5, int b6, int b7, int b8, int ring)
// {
//   Avec bm(8,1);
//   bm[0]=b1;
//   bm[1]=b2;
//   bm[2]=b3;
//   bm[3]=b4;
//   bm[4]=b5;
//   bm[5]=b6;
//   bm[6]=b7;
//   bm[7]=b8;
//   switch(ring){
//   case 0: //Ring 4
//     beam_mask_r4=bm;
//     break;
//   case 1:
//     beam_mask_r6=bm;
//     break;
//   }
//   beam_mask=bm;
//   return;
// }

//! TEC z-positions
Avec TECReconstructor::make_zdat(Int_t bs_frame)
{
  Avec bzdat(9,0.0);
  bzdat[0] = 0;
  bzdat[1] = 140;
  bzdat[2] = 280;
  bzdat[3] = 420;
  bzdat[4] = 560;
  bzdat[5] = 735;
  bzdat[6] = 925;
  bzdat[7] = 1130;
  bzdat[8] = 1345;
  
  switch(bs_frame){
  case 0:
    bzdat+=1250;
    break;
  case 1:
    bzdat-=705;  // Beamsplitter at z=0
    break;
  }
  return bzdat;
}

//! Initialize for one ring reconstruction with symmetric (all) beams
/*! Discs can be masked out with disc_mask
 */
void TECReconstructor::Init_onering_symm(){
  if(verbose_flag > 1){
    std::cout << "**************************************\n";
    std::cout << "* Initializing TEC onering symmetric *\n";
    std::cout << "**************************************\n" << std::endl;
  }

  //using namespace TEC;
  theta = make_theta();
  zdat  = make_zdat();
  nbeam = 8;
  
  L = zdat[8] - zdat[0];
  
  // Constant parameters
  Sz=vsum(zdat);
  Sz2=vsum(zdat*zdat);
  St=vsum(theta);
  St2=vsum(theta*theta);
  // Sums over sin/cos of beam positions (using symmetry considerations)
  Ssin=0;
  Scos=0;
  Ssin2=nbeam/2;
  Scos2=nbeam/2;
  Scossin=0;
  
  A = 0;
  B = 0;
  C = 0;
  D = 0;
  E = 0;
  F = 0;
  denom1=nbeam*(Sz*Sz-ndisc*Sz2);
  term1 = 0;
  //nenn2=nbeam*St2-St*St;
  //znorm=ndisc*Sz2-Sz*Sz;
  
  return;
}

// Initialize for one ring reconstruction
// Some beams can be masked out using Avec beam_mask
void TECReconstructor::Init_onering(){
  if(verbose_flag > 2){
    std::cout << "****************************\n";
    std::cout << "* Initializing TEC onering *\n";
    std::cout << "****************************\n" << std::endl;
  }

  //using namespace TEC;
  //std::cout << "beam_mask:\n" << beam_mask << std::endl;
  theta = make_theta();
  zdat  = make_zdat(2);
  nbeam = 8;
  
  L=zdat[8]-zdat[0];
  // Calculate some constant parameters
  Sz=vsum(zdat);
  Sz2=vsum(zdat*zdat);
  St=vsum(theta);
  St2=vsum(theta*theta);
  
  // Sums over sin/cos of beam positions
  Ssin=vsum(sin(theta));
  Scos=vsum(cos(theta));
  Ssin2=vsum(pow(sin(theta),2));
  Scos2=vsum(pow(cos(theta),2));
  Scossin=vsum(cos(theta)*sin(theta));
  
  // Common terms of alignment parameters
  A = Scos2*Ssin - Scos*Scossin;
  B = Scos*Ssin2 - Scossin*Ssin;
  C = Scossin*Scossin - Scos2*Ssin2;
  D = Scos*Scos - nbeam*Scos2;
  E = nbeam*Scossin-Scos*Ssin;
  F = Ssin*Ssin - nbeam*Ssin2;
  term1 = Sz*Sz-ndisc*Sz2;
  denom1 = term1*(nbeam*C+(Scos*Scos*Ssin2+Scos2*Ssin*Ssin-2*Scos*Scossin*Ssin));
  term2=(nbeam*(Scos2*Ssin2-Scossin*Scossin)-Scos*Scos*Ssin2-Scos2*Ssin*Ssin+2*Scos*Scossin*Ssin);
  nenn2=nbeam*St2-St*St;
  znorm=ndisc*Sz2-Sz*Sz;

  //denom2= nbeam*Scos2*Ssin2-nbeam*Scossin*Scossin-Scos*Scos*Ssin2-Scos2*Ssin*Ssin+2*Scos*Scossin*Ssin;
  
  //nbeam*Scos2*Ssin2-nbeam*Scossin*Scossin-Scos*Scos*Ssin2-Scos2*Ssin*Ssin+2*Scos*Scossin*Ssin
  return;
}

//! Set masked modules to fit result
int TECReconstructor::set_masked_modules_tec(LASGlobalData<double>& data, const LASGlobalData<int>& modmask, const TecRingPar& alpar, LAS::beam_group bg, int print)
{
  if (bg != LAS::TEC_PLUS_R4 && bg != LAS::TEC_PLUS_R6 && bg != LAS::TEC_MINUS_R4 && bg != LAS::TEC_MINUS_R6 && bg != LAS::TEC_MINUS && bg != LAS::TEC_PLUS){
    std::cerr << "Error in set_masked_modules_tec: bg has invalid value " << bg << std::endl;
    return 0;
  }

  double epsilon = 1e-9;

  //using namespace TEC;
  int flag=0; // Flag to determine end of iteration

  LASGlobalDataLoop loop(convert_to_loop_type(bg));
  do{
    if( loop.GetEntry<int>(modmask) ) continue;
    
    if(print>1){
      std::cout << "Masked module " << GetModuleName(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos()) << "  data: " << loop.GetEntry<double>(data) << std::endl;
    }
    
    double dphi_rec = dphi_rec_tec(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos(), alpar, 0);
    double diff = loop.GetEntry<double>(data) - dphi_rec;
    if(print) std::cout << "phi_rec: " << dphi_rec << "   diff: " << diff << std::endl;
    if(fabs(diff) > epsilon) flag=1;
    loop.GetEntry<double>(data) = dphi_rec;
  }while( loop.next() );
  
  return flag;
}

//! Calculate TEC Alignment Parameters for one ring
/*! Calculate disc rotations and displacements for TEC and its beams\n
  Only one ring is used and beams can be masked with mask_beam \n
  TEC::Init_onering() should have been called first
*/
void TECReconstructor::reconstruct_tec_set_onering_weighted(const LASGlobalData<double>& data, const LASGlobalData<double>& weights, TecRingPar& alpar, LAS::beam_group bg, int print)
{
  // Radius and Phi-Positions of beams
  double rbeam;
  switch (bg){
  case LAS::TEC_PLUS_R4:
    rbeam = LAS::r_tec_r4; // radius of beams R4
    break;
  case LAS::TEC_PLUS_R6:
    rbeam = LAS::r_tec_r6; // radius of beams R6
    break;
  case LAS::TEC_MINUS_R4:
    rbeam = LAS::r_tec_r4; // radius of beams R4
    break;
  case LAS::TEC_MINUS_R6:
    rbeam = LAS::r_tec_r6; // radius of beams R6
    break;
  default:
    std::cerr << "Error in reconstruct_tec_ring: bg has invalid value " << bg << std::endl;
    return;
  }

  double Sy=0;
  double Syz=0;
  double Sysin=0;
  double Sycos=0;
  double Syzsin=0;
  double Syzcos=0;
  Avec Syi(ndisc);
  Avec Sysini(ndisc);
  Avec Sycosi(ndisc);
  Avec Syk(nbeam);
  Avec Syzk(nbeam);

  double Swy;
  Avec Swyi(ndisc);
  Avec Swyk(nbeam);


  // Loop over data and build coefficients
  LASGlobalDataLoop loop(convert_to_loop_type(bg));
  do{
    int i = loop.get_beam(); // beam index
    int k = loop.get_zpos(); // disc index
    //if(beam_mask[i]==0) continue ;
    
    double phi = loop.GetEntry<double>(data); // measured position (phi)
    double w = loop.GetEntry<double>(weights); // weight of this measurement
    
    double y = phi*rbeam * w; // measured position in mm (weighted)
    
    Sy      += y;
    Syz     += y * zdat[k];
    Sysin   += y * sin(theta[i]);
    Sycos   += y * cos(theta[i]);
    Syzsin  += y * zdat[k] * sin(theta[i]);
    Syzcos  += y * zdat[k] * cos(theta[i]);
    Syi[k]    += y;
    Sysini[k] += y * sin(theta[i]);
    Sycosi[k] += y * cos(theta[i]);
    Syk[i]    += y;
    Syzk[i]   += y * zdat[k];

    Swy += w;
    Swyi[k] += w;
    Swyk[i] += w;
  } while (loop.next());

  // Coefficients
  term1 = Sz*Sz-ndisc*Sz2;
  denom1 = term1*(nbeam*C+(Scos*Scos*Ssin2+Scos2*Ssin*Ssin-2*Scos*Scossin*Ssin));

  // Global Parameters
  alpar.Dphi0= ((A*Sysin+B*Sycos+C*Sy)*Sz2 - (A*Syzsin+B*Syzcos+C*Syz)*Sz) / denom1 / rbeam;
  alpar.Dx0=(Sz*(Syzcos*E+Syzsin*D+Syz*A)-Sz2*(Sycos*E+Sysin*D+Sy*A))/denom1;
  alpar.Dy0=(Sz2*(Sycos*F+Sysin*E+Sy*B)-Sz*(Syzcos*F+Syzsin*E+Syz*B))/denom1;

  alpar.Dphit=(ndisc*(Syz*C+Syzcos*B+Syzsin*A)-Sz*(Sy*C+Sycos*B+Sysin*A))*L/denom1/rbeam;
  alpar.Dxt=(Sz*(Sycos*E+Sysin*D+Sy*A)-ndisc*(Syzcos*E+Syzsin*D+Syz*A))*L/denom1;
  alpar.Dyt=(ndisc*(Syzcos*F+Syzsin*E+Syz*B)-Sz*(Sycos*F+Sysin*E+Sy*B))*L/denom1;

  // Disc Parameters
  alpar.Dphik=(
	       (-Syi*term1-Sy*Sz2+(Syz+zdat*Sy)*Sz-zdat*ndisc*Syz)*C
	       -(Sycosi*term1+Sycos*Sz2-(Syzcos+zdat*Sycos)*Sz+zdat*ndisc*Syzcos)*B
	       -(Sysini*term1+Sysin*Sz2-(Syzsin+zdat*Sysin)*Sz+zdat*ndisc*Syzsin)*A
	       )/denom1/rbeam;
  
  alpar.Dxk=(
	     (Syi*term1+Sy*Sz2-Syz*Sz-zdat*Sy*Sz+zdat*ndisc*Syz)*A
	     +(Sycosi*term1+Sycos*Sz2-Syzcos*Sz-zdat*Sycos*Sz+zdat*ndisc*Syzcos)*E
	     +(Sysini*term1+Sysin*Sz2-Syzsin*Sz-zdat*Sysin*Sz+zdat*ndisc*Syzsin)*D
	     )/denom1;

  alpar.Dyk=-(
	      (Sycosi*term1+Sycos*Sz2-Syzcos*Sz-zdat*Sycos*Sz+zdat*ndisc*Syzcos)*F
	      +(Sysini*term1+Sysin*Sz2-Syzsin*Sz-zdat*Sysin*Sz+zdat*ndisc*Syzsin)*E
	      +(Syi*term1+Sy*Sz2-Syz*Sz-zdat*Sy*Sz+zdat*ndisc*Syz)*B
	      )/denom1;

  // Beam Parameters
  alpar.DthetaA = (Syzk*Sz - Syk*Sz2)/term1/rbeam + alpar.Dphi0 - alpar.Dx0*sin(theta)/rbeam + alpar.Dy0*cos(theta)/rbeam;
  
  alpar.DthetaB = (Syzk*ndisc-Syk*Sz)*L/term1/rbeam + alpar.Dxt*sin(theta)/rbeam - alpar.Dyt*cos(theta)/rbeam - alpar.Dphit - alpar.DthetaA;
  
  if(print)alpar.print();
  
  return;
}

//! Calculate TEC Alignment Parameters for one ring
/*! Calculate disc rotations and displacements for TEC and its beams\n
  Only one ring is used and beams can be masked with mask_beam \n
  TEC::Init_onering() should have been called first
*/
void TECReconstructor::reconstruct_tec_set_onering(const LASGlobalData<double>& data, TecRingPar& alpar, LAS::beam_group bg, int print)
{
  
  // Radius and Phi-Positions of beams
  double rbeam;
  switch (bg){
  case LAS::TEC_PLUS_R4:
    rbeam = LAS::r_tec_r4; // radius of beams R4
    break;
  case LAS::TEC_PLUS_R6:
    rbeam = LAS::r_tec_r6; // radius of beams R6
    break;
  case LAS::TEC_MINUS_R4:
    rbeam = LAS::r_tec_r4; // radius of beams R4
    break;
  case LAS::TEC_MINUS_R6:
    rbeam = LAS::r_tec_r6; // radius of beams R6
    break;
  default:
    std::cerr << "Error in reconstruct_tec_ring: bg has invalid value " << bg << std::endl;
    return;
  }

  double Sy=0;
  double Syz=0;
  double Sysin=0;
  double Sycos=0;
  double Syzsin=0;
  double Syzcos=0;
  Avec Syi(ndisc);
  Avec Sysini(ndisc);
  Avec Sycosi(ndisc);
  Avec Syk(nbeam);
  Avec Syzk(nbeam);

  // Loop over data and build coefficients
  LASGlobalDataLoop loop(convert_to_loop_type(bg));
  do{
    int i = loop.get_beam(); // beam index
    int k = loop.get_zpos(); // disc index
    //if(beam_mask[i]==0) continue ;
    
    double phi = loop.GetEntry<double>(data); // measured position (phi)
    
    double y = phi*rbeam; // measured position in mm
    
    Sy      += y;
    Syz     += y * zdat[k];
    Sysin   += y * sin(theta[i]);
    Sycos   += y * cos(theta[i]);
    Syzsin  += y * zdat[k] * sin(theta[i]);
    Syzcos  += y * zdat[k] * cos(theta[i]);
    Syi[k]    += y;
    Sysini[k] += y * sin(theta[i]);
    Sycosi[k] += y * cos(theta[i]);
    Syk[i]    += y;
    Syzk[i]   += y * zdat[k];
  } while (loop.next());
  
  // Global Parameters
  alpar.Dphi0= ((A*Sysin + B*Sycos + C*Sy)*Sz2 - (A*Syzsin + B*Syzcos + C*Syz)*Sz) / denom1 / rbeam;
  alpar.Dx0=(Sz*(Syzcos*E+Syzsin*D+Syz*A)-Sz2*(Sycos*E+Sysin*D+Sy*A))/denom1;
  alpar.Dy0=(Sz2*(Sycos*F+Sysin*E+Sy*B)-Sz*(Syzcos*F+Syzsin*E+Syz*B))/denom1;

  alpar.Dphit=(ndisc*(Syz*C + Syzcos*B + Syzsin*A)-Sz*(Sy*C + Sycos*B + Sysin*A)) * L / denom1 / rbeam;
  alpar.Dxt=(Sz*(Sycos*E+Sysin*D+Sy*A)-ndisc*(Syzcos*E+Syzsin*D+Syz*A))*L/denom1;
  alpar.Dyt=(ndisc*(Syzcos*F+Syzsin*E+Syz*B)-Sz*(Sycos*F+Sysin*E+Sy*B))*L/denom1;

  // Disc Parameters
//   alpar.Dphik=(
// 	       (-Syi*term1-Sy*Sz2+(Syz+zdat*Sy)*Sz-zdat*ndisc*Syz)*C
// 	       -(Sycosi*term1+Sycos*Sz2-(Syzcos+zdat*Sycos)*Sz+zdat*ndisc*Syzcos)*B
// 	       -(Sysini*term1+Sysin*Sz2-(Syzsin+zdat*Sysin)*Sz+zdat*ndisc*Syzsin)*A
// 	       )/denom1/rbeam;

  alpar.Dphik = (A*Sysini + B*Sycosi + C*Syi)/term2/rbeam - alpar.Dphi0 - alpar.Dphit*zdat/L;
  
  alpar.Dxk=(
	     (Syi*term1+Sy*Sz2-Syz*Sz-zdat*Sy*Sz+zdat*ndisc*Syz)*A
	     +(Sycosi*term1+Sycos*Sz2-Syzcos*Sz-zdat*Sycos*Sz+zdat*ndisc*Syzcos)*E
	     +(Sysini*term1+Sysin*Sz2-Syzsin*Sz-zdat*Sysin*Sz+zdat*ndisc*Syzsin)*D
	     )/denom1;

  alpar.Dyk=-(
	      (Sycosi*term1+Sycos*Sz2-Syzcos*Sz-zdat*Sycos*Sz+zdat*ndisc*Syzcos)*F
	      +(Sysini*term1+Sysin*Sz2-Syzsin*Sz-zdat*Sysin*Sz+zdat*ndisc*Syzsin)*E
	      +(Syi*term1+Sy*Sz2-Syz*Sz-zdat*Sy*Sz+zdat*ndisc*Syz)*B
	      )/denom1;

  // Beam Parameters
  alpar.DthetaA = (Syzk*Sz - Syk*Sz2)/term1/rbeam + alpar.Dphi0 - alpar.Dx0*sin(theta)/rbeam + alpar.Dy0*cos(theta)/rbeam;
  
  alpar.DthetaB = (Syzk*ndisc-Syk*Sz)*L/term1/rbeam + alpar.Dxt*sin(theta)/rbeam - alpar.Dyt*cos(theta)/rbeam - alpar.Dphit - alpar.DthetaA;
  
  if(print)alpar.print();
  
  return;
}


//! Calculate reconstrucred phi-position for a given TEC module
double TECReconstructor::dphi_rec_tec(int subdet, int ring, int beam, int zpos, const TecRingPar& alpar, int print)
{
  if(subdet != 0 && subdet != 1)return 0.0;
  int i=beam;
  int k=zpos;
  double r0 =  (ring==0 || ring==-1) ? r4 : r6;
  
  double thet = 0;
  thet = theta[i];
  
  if(print>0){
    std::cout << "r0: " << r0 << std::endl;
    std::cout << "thet : " << thet << std::endl;
  }  

  double phi_reco = 
    - alpar.Dphi0 + sin(thet)*alpar.Dx0/r0 - cos(thet)*alpar.Dy0/r0
    - zdat[k]/L*(alpar.Dphit - sin(thet)*alpar.Dxt/r0 + cos(thet)*alpar.Dyt/r0)
    - alpar.Dphik[k] + sin(thet)*alpar.Dxk[k]/r0 - cos(thet)*alpar.Dyk[k]/r0
    + alpar.DthetaA[i] + zdat[k]/L*alpar.DthetaB[i];
  
  return phi_reco;
}
 

//! Calculate reconstrucred phi-position for a given TEC module
double TECReconstructor::dphi_rec_tec(int subdet, int ring, int beam, int zpos, const TecPar& alpar, int print)
{
  if(subdet != 0 && subdet != 1)return 0.0;
  int i=beam;
  int k=zpos;
  double r0 =  (ring==0 || ring==-1) ? r4 : r6;
  double DthetaA, DthetaB;
  if(ring == 0){
    DthetaA = alpar.DthetaA_r4[i];
    DthetaB = alpar.DthetaB_r4[i];
  }
  else{
    DthetaA = alpar.DthetaA_r6[i];
    DthetaB = alpar.DthetaB_r6[i];
  }

  double thet = 0;
  thet = theta[i];
  
  if(print>0){
    std::cout << "r0: " << r0 << std::endl;
    std::cout << "thet : " << thet << std::endl;
  }

  double phi_reco = 
    - alpar.Dphi0 + sin(thet)*alpar.Dx0/r0 - cos(thet)*alpar.Dy0/r0
    - zdat[k]/L*(alpar.Dphit - sin(thet)*alpar.Dxt/r0 + cos(thet)*alpar.Dyt/r0)
    - alpar.Dphik[k] + sin(thet)*alpar.Dxk[k]/r0 - cos(thet)*alpar.Dyk[k]/r0
    + DthetaA + zdat[k]/L*DthetaB;
  
  return phi_reco;
}
 
//! Calculate reconstrucred position in mm for a given TEC module
double TECReconstructor::dx_rec_tec(int subdet, int ring, int beam, int zpos, const TecRingPar& alpar, int print)
{
  if(subdet != 0 && subdet != 1)return 0.0;
  int i=beam;
  int k=zpos;
  double r0 =  (ring==0 || ring==-1) ? r4 : r6;

  double thet = 0;
  thet = theta[i];
  
  if(print>0){
    std::cout << "r0: " << r0 << std::endl;
    std::cout << "thet : " << thet << std::endl;
  }

  double dx_reco = 
    - alpar.Dphi0*r0 + sin(thet)*alpar.Dx0 - cos(thet)*alpar.Dy0
    - zdat[k]/L*(alpar.Dphit*r0 - sin(thet)*alpar.Dxt + cos(thet)*alpar.Dyt)
    - alpar.Dphik[k]*r0 + sin(thet)*alpar.Dxk[k] - cos(thet)*alpar.Dyk[k]
    + alpar.DthetaA[i]*r0 + zdat[k]/L*alpar.DthetaB[i]*r0;
  
  return dx_reco;
}
 
//! Calculate reconstrucred position in mm for a given TEC module
double TECReconstructor::dx_rec_tec(int subdet, int ring, int beam, int zpos, const TecPar& alpar, int print)
{
  if(subdet != 0 && subdet != 1)return 0.0;

  TecRingPar alpar2;
  alpar2.Dphi0 = alpar.Dphi0;
  alpar2.Dx0 = alpar.Dx0;
  alpar2.Dy0 = alpar.Dy0;

  alpar2.Dphit = alpar.Dphit;
  alpar2.Dxt = alpar.Dxt;
  alpar2.Dyt = alpar.Dyt;

  alpar2.Dphik = alpar.Dphik;
  alpar2.Dxk = alpar.Dxk;
  alpar2.Dyk = alpar.Dyk;

  if(ring == 0){
    alpar2.DthetaA = alpar.DthetaA_r4;
    alpar2.DthetaB = alpar.DthetaB_r4;
  }
  else{
    alpar2.DthetaA = alpar.DthetaA_r6;
    alpar2.DthetaB = alpar.DthetaB_r6;
  }

  return dx_rec_tec(subdet, ring, beam, zpos, alpar2, print);
}
 
//! Calculate reconstrucred x-position for both endcaps
LASGlobalData<double> TECReconstructor::dx_rec_tec(const LasAlPar& alpar, int print)
{

  LASGlobalData<double> retval(0.0);
  LASGlobalDataLoop loop(LASGlobalDataLoop::TEC);
  do{
    loop.GetEntry<double>(retval) = dx_rec_tec(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos(), loop.get_det() == 1 ? alpar.tecm : alpar.tecp, print);
  } while(loop.next());

  return retval;
}
 
//! Calculate reconstrucred x-position for one endcap
LASGlobalData<double> TECReconstructor::dx_rec_tec(const TecPar& alpar, LAS::beam_group bg, int print)
{
  if(bg != LAS::TEC_PLUS && bg != LAS::TEC_MINUS) {
    std::cerr << "Error in TECReconstructor::dx_rec_tec(const TecPar& alpar, LAS::beam_group bg, int print): beam group " << bg <<  " not implemented yet" << std::endl;
    return LASGlobalData<double>(0.0);
  }

  LASGlobalData<double> retval(0.0);
  LASGlobalDataLoop loop = convert_to_loop_type(bg);

  do{
    loop.GetEntry<double>(retval) = dx_rec_tec(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos(), alpar, print);
  } while(loop.next());

  return retval;
}
 
//! Calculate reconstrucred x-position for one endcap
LASGlobalData<double> TECReconstructor::dx_rec_tec(const TecRingPar& alpar, LAS::beam_group bg, int print)
{
  if(bg != LAS::TEC_PLUS_R4 && bg != LAS::TEC_PLUS_R6 && bg != LAS::TEC_MINUS_R4 && bg != LAS::TEC_MINUS_R6) {
    std::cerr << "Error in TECReconstructor::dx_rec_tec(const TecRingPar& alpar, LAS::beam_group bg, int print): beam group " << bg <<  " not implemented yet" << std::endl;
    return LASGlobalData<double>(0.0);
  }

  LASGlobalData<double> retval(0.0);
  LASGlobalDataLoop loop = convert_to_loop_type(bg);

  do{
    loop.GetEntry<double>(retval) = dx_rec_tec(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos(), alpar, print);
  } while(loop.next());

  return retval;
}
 
//! Calculate reconstrucred phi-position for a given ring
LASGlobalData<double> TECReconstructor::dphi_rec_tec(const TecRingPar& alpar, LAS::beam_group bg, int print)
{
  if(bg != LAS::TEC_PLUS_R4 && bg != LAS::TEC_PLUS_R6 && bg != LAS::TEC_MINUS_R4 && bg != LAS::TEC_MINUS_R6) {
    std::cerr << "Error in TECReconstructor::dphi_rec_tec(const TecRingPar& alpar, LAS::beam_group bg, int print): beam group " << bg <<  " not implemented yet" << std::endl;
    return LASGlobalData<double>(0.0);
  }

  LASGlobalData<double> retval(0.0);
  LASGlobalDataLoop loop = convert_to_loop_type(bg);

  do{
    loop.GetEntry<double>(retval) = dphi_rec_tec(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos(), alpar, print);
  } while(loop.next());

  return retval;
}
 

//! Calculate reconstrucred phi-position for a given beam group
LASGlobalData<double> TECReconstructor::dphi_rec_tec(const TecPar& alpar, LAS::beam_group bg, int print)
{
  if(bg != LAS::TEC_PLUS && bg != LAS::TEC_MINUS) {
    std::cerr << "Error in TECReconstructor::dphi_rec_tec(const TecPar& alpar, LAS::beam_group bg, int print): beam group " << bg <<  " not implemented yet" << std::endl;
    return LASGlobalData<double>(0.0);
  }

  LASGlobalData<double> retval(0.0);
  LASGlobalDataLoop loop = convert_to_loop_type(bg);

  do{
    loop.GetEntry<double>(retval) = dphi_rec_tec(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos(), alpar, print);
  } while(loop.next());

  return retval;
}
 

//! Calculate reconstrucred phi-position for both endcaps
LASGlobalData<double> TECReconstructor::dphi_rec_tec(const LasAlPar& alpar, int print)
{

  LASGlobalData<double> retval(0.0);
  LASGlobalDataLoop loop(LASGlobalDataLoop::TEC);
  do{
    loop.GetEntry<double>(retval) = dphi_rec_tec(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos(), loop.get_det() == 0 ? alpar.tecm : alpar.tecp, print);
  } while(loop.next());

  return retval;
}
 

//! Set masked modules to fit result
int TECReconstructor::set_masked_modules_tec(LASGlobalData<double>& data, const LASGlobalData<int>& modmask, const TecPar& alpar, LAS::beam_group bg, int print)
{
  if (bg != LAS::TEC_MINUS && bg != LAS::TEC_PLUS){
    std::cerr << "Error in set_masked_modules_tec: bg has invalid value " << bg << std::endl;
    return 0;
  }

  double epsilon = 1e-9;

  //using namespace TEC;
  int flag=0; // Flag to determine end of iteration

  LASGlobalDataLoop loop(convert_to_loop_type(bg));
  do{
    if( loop.GetEntry<int>(modmask) ) continue;
    
    if(print>1){
      std::cout << "Masked module " << GetModuleName(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos()) << "  data: " << loop.GetEntry<double>(data) << std::endl;
    }
    
    double dx_rec = dx_rec_tec(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos(), alpar, 0);
    double diff = loop.GetEntry<double>(data) - dx_rec;
    if(print) std::cout << "dx_rec: " << dx_rec << "[mm]   diff: " << diff << " [mm]" << std::endl;
    if(fabs(diff) > epsilon) flag=1;
    loop.GetEntry<double>(data) = dx_rec;
  }while( loop.next() );
  
  return flag;
}



// Initialize for full tEC reconstruction
// Some beams can be masked out using Avec beam_mask
void TECReconstructor::Init_full(){
  if(verbose_flag > 1){
    std::cout << "*************************\n";
    std::cout << "* Initializing TEC full *\n";
    std::cout << "*************************\n" << std::endl;
  }

//   std::cout << "beam_mask:\n" << beam_mask << std::endl;
  theta = make_theta();
  zdat  = make_zdat(2);
  nbeam = 8;
  ndisc = 9;
  nring = 2;
  
  L=zdat[8]-zdat[0];
  // Calculate some constant parameters
  Sz=vsum(zdat);
  Sz2=vsum(zdat*zdat);
//   St=vsum(theta*beam_mask);
//   St2=vsum(theta*theta*beam_mask);
  
  rbeam.clear();
  rbeam.push_back(r4);
  rbeam.push_back(r6);
  Sr = vsum(rbeam);
  Sr2 = vsum(rbeam * rbeam);

  // Sums over sin/cos of beam positions
  Ssin=vsum(sin(theta));
  Scos=vsum(cos(theta));
  Ssin2=vsum(pow(sin(theta),2));
  Scos2=vsum(pow(cos(theta),2));
  Scossin=vsum(cos(theta)*sin(theta));

  // Common terms of alignment parameters
  A = Sr*(Scos2*Ssin - Scos*Scossin);
  B = Sr*(Scos*Ssin2 - Scossin*Ssin);
  C = Scossin*Scossin - Scos2*Ssin2;

  D = nbeam*nring*Scossin*Sr2 - Scos*Ssin*Sr*Sr;
  E = Scos*Scos*Sr*Sr - nbeam*nring*Scos2*Sr2;
  F = Ssin*Ssin*Sr*Sr - nbeam*nring*Ssin2*Sr2;

  //D = nring*Sr2*(Scos*Scos - nbeam*Scos2);
  //E = Sr*Sr*(nbeam*Scossin-Scos*Ssin);
  //F = Sr2*(Ssin*Ssin - nbeam*Ssin2);
//   term1 = Sz*Sz-ndisc*Sz2;
//   denom1 = term1*(nbeam*C+(Scos*Scos*Ssin2+Scos2*Ssin*Ssin-2*Scos*Scossin*Ssin));
//   term2=(nbeam*(Scos2*Ssin2-Scossin*Scossin)-Scos*Scos*Ssin2-Scos2*Ssin*Ssin+2*Scos*Scossin*Ssin);
//   nenn2=nbeam*St2-St*St;
  znorm=ndisc*Sz2-Sz*Sz;

//   //denom2= nbeam*Scos2*Ssin2-nbeam*Scossin*Scossin-Scos*Scos*Ssin2-Scos2*Ssin*Ssin+2*Scos*Scossin*Ssin;
  
  //nbeam*Scos2*Ssin2-nbeam*Scossin*Scossin-Scos*Scos*Ssin2-Scos2*Ssin*Ssin+2*Scos*Scossin*Ssin
  return;
}

//! Calculate TEC Alignment Parameters for two rings
/*! Calculate disc rotations and displacements for TEC and its beams\n
  TECReconstructor::Init_full() should have been called first
*/
void TECReconstructor::reconstruct_tec_set_full(const LASGlobalData<double>& data, TecPar& alpar, LAS::beam_group bg, int print)
{
  if(bg != LAS::TEC_MINUS && bg != LAS::TEC_PLUS){
    std::cerr << "Error in reconstruct_tec_set_full: bg has invalid value " << bg << std::endl;
    return;
  }

  double Syr = 0;
  double Syzr = 0;
  double Sysin=0;
  double Sycos=0;
  double Syzsin=0;
  double Syzcos=0;
  Avec Syir(ndisc, 0.0);
  Avec Sysini(ndisc, 0.0);
  Avec Sycosi(ndisc,0.0);
  Avec2D Syk(nring, nbeam, 0);
  Avec2D Syzk(nring, nbeam, 0);
  
  //const double rnorm = nring*Sr2-Sr*Sr;
  const double denom2 = nbeam*nring*Scos2*Sr2*Ssin2-Scos*Scos*Sr*Sr*Ssin2-Scos2*Sr*Sr*Ssin*Ssin+2*Scos*Scossin*Sr*Sr*Ssin-nbeam*nring*Scossin*Scossin*Sr2;
  const double denom3 = denom2 * znorm;
  
  // Loop over data and build coefficients
  LASGlobalDataLoop loop(convert_to_loop_type(bg));
  do{
    int i = loop.get_beam(); // beam index
    int k = loop.get_zpos(); // disc index
    int j = loop.get_ring(); // ring index
    
    //double phi = loop.GetEntry<double>(data); // measured position (phi)
    //double y = phi*rbeam[j]; // measured position in mm
    double y = loop.GetEntry<double>(data); // measured position (mm)
    
    Syr += y * rbeam[j]; 
    Syzr += rbeam[j] * y * zdat[k];
    //     Sy      += y;
    //     Syz     += y * zdat[k];
    
    Sysin   += y * sin(theta[i]);
    Sycos   += y * cos(theta[i]);
    Syzsin  += y * zdat[k] * sin(theta[i]);
    Syzcos  += y * zdat[k] * cos(theta[i]);
    
    Syir[k]    += y * rbeam[j];
    Sysini[k] += y * sin(theta[i]);
    Sycosi[k] += y * cos(theta[i]);
    Syk[j][i]    += y;
    Syzk[j][i]   += y * zdat[k];
  } while (loop.next());

  const double s01 = Syr*Sz2 - Syzr*Sz;
  const double s02 = Sycos*Sz2 - Syzcos*Sz;
  const double s03 = Sysin*Sz2 - Syzsin*Sz;

  const double st1 = Syr*Sz - ndisc*Syzr;
  const double st2 = Sycos*Sz - ndisc*Syzcos;
  const double st3 = Sysin*Sz - ndisc*Syzsin;

  // Global Parameters
  alpar.Dphi0 =  (nring*s01*C + s02*B + s03*A) / denom3;
  alpar.Dx0   = -(nring*s01*A + s02*D + s03*E) / (denom3*nring);
  alpar.Dy0 =    (nring*s01*B + s02*F + s03*D) / (denom3*nring);

  alpar.Dphit = -(nring*st1*C + st2*B + st3*A) * L / denom3;
  alpar.Dxt   =  (nring*st1*A + st2*D + st3*E) * L / (denom3*nring);
  alpar.Dyt   = -(nring*st1*B + st2*F + st3*D) * L / (denom3*nring);

  // Disc Parameters
  alpar.Dphik = ( A*Sysini + B*Sycosi + C*nring*Syir )/denom2 - alpar.Dphi0 - alpar.Dphit*zdat/L;

  alpar.Dxk=-(
	      ((alpar.Dx0*nbeam*nring*nring*Scos2*Sr2-alpar.Dx0*nring*Scos*Scos*Sr*Sr)*Ssin2-alpar.Dx0*nring*Scos2*Sr*Sr*Ssin*Ssin+((2*alpar.Dx0*nring*Scos*Scossin-Sycosi*Scos)*Sr*Sr+Syir*nring*Scos2*Sr)*Ssin+(-alpar.Dx0*nbeam*nring*nring*Scossin*Scossin+Sycosi*nbeam*nring*Scossin-Sysini*nbeam*nring*Scos2)*Sr2+Sysini*Scos*Scos*Sr*Sr-Syir*nring*Scos*Scossin*Sr)*L+(alpar.Dxt*zdat*nbeam*nring*nring*Scos2*Sr2-alpar.Dxt*zdat*nring*Scos*Scos*Sr*Sr)*Ssin2-alpar.Dxt*zdat*nring*Scos2*Sr*Sr*Ssin*Ssin+2*alpar.Dxt*zdat*nring*Scos*Scossin*Sr*Sr*Ssin-alpar.Dxt*zdat*nbeam*nring*nring*Scossin*Scossin*Sr2
	      )/(
		 ((nbeam*nring*nring*Scos2*Sr2-nring*Scos*Scos*Sr*Sr)*Ssin2-nring*Scos2*Sr*Sr*Ssin*Ssin+2*nring*Scos*Scossin*Sr*Sr*Ssin-nbeam*nring*nring*Scossin*Scossin*Sr2)*L
		 );

  alpar.Dyk = -(
		Sycosi*(nbeam*nring*Sr2*Ssin2-Sr*Sr*Ssin*Ssin)*L-Syir*nring*Sr*(Scos*Ssin2-Scossin*Ssin)*L+Sysini*(Scos*Sr*Sr*Ssin-nbeam*nring*Scossin*Sr2)*L+denom2*alpar.Dy0*nring*L+denom2*alpar.Dyt*zdat*nring)/(denom2*nring*L);
  //   alpar.Dyk=-(
  // 	      (Sycosi*term1+Sycos*Sz2-Syzcos*Sz-zdat*Sycos*Sz+zdat*ndisc*Syzcos)*F
  // 	      +(Sysini*term1+Sysin*Sz2-Syzsin*Sz-zdat*Sysin*Sz+zdat*ndisc*Syzsin)*E
  // 	      +(Syi*term1+Sy*Sz2-Syz*Sz-zdat*Sy*Sz+zdat*ndisc*Syz)*B
  // 	      )/denom1;

  
  //   // Beam Parameters
 alpar.DthetaA_r4 = ( Syk[0]*Sz2 - Syzk[0]*Sz )/(rbeam[0]*znorm) + alpar.Dphi0 - alpar.Dx0*sin(theta)/rbeam[0] + alpar.Dy0*cos(theta)/rbeam[0];
 alpar.DthetaA_r6 = ( Syk[1]*Sz2 - Syzk[1]*Sz )/(rbeam[1]*znorm) + alpar.Dphi0 - alpar.Dx0*sin(theta)/rbeam[1] + alpar.Dy0*cos(theta)/rbeam[1];
  
 alpar.DthetaB_r4 = -(Syk[0]*Sz-Syzk[0]*ndisc)*L/rbeam[0]/znorm + alpar.Dphit - (alpar.Dxt*sin(theta)-alpar.Dyt*cos(theta))/rbeam[0];
 alpar.DthetaB_r6 = -(Syk[1]*Sz-Syzk[1]*ndisc)*L/rbeam[1]/znorm + alpar.Dphit - (alpar.Dxt*sin(theta)-alpar.Dyt*cos(theta))/rbeam[1];

 //alpar.DthetaB_r4 = (Syk[0]*Sz-Syzk[0]*ndisc)*L/rbeam[0]/znorm - alpar.Dphit + (alpar.Dxt*sin(theta)-alpar.Dyt*cos(theta))/rbeam[0];
 //alpar.DthetaB_r6 = (Syk[1]*Sz-Syzk[1]*ndisc)*L/rbeam[1]/znorm - alpar.Dphit + (alpar.Dxt*sin(theta)-alpar.Dyt*cos(theta))/rbeam[1];



 //DthetaB[i,j]=-(Syk[i,j]*Sz*L)/(R[j]*znorm)+(Syzk[i,j]*m*L)/(R[j]*znorm)  -  (Dxt*sin(theta[i]))/R[j]  +  (Dyt*cos(theta[i]))/R[j]  +  Dphit


//  alpar.DthetaB_r4 = 
//    ( Syk[0]*(Sz*L - Sz2) + Syzk[0]*(Sz - ndisc*L) )/(rbeam[0]*znorm)
//    - alpar.Dphit - alpar.Dphi0 + ( (alpar.Dxt + alpar.Dx0)*sin(theta) - (alpar.Dyt + alpar.Dy0)*cos(theta) )/rbeam[0];
//  alpar.DthetaB_r6 = 
//    ( Syk[1]*(Sz*L - Sz2) + Syzk[1]*(Sz - ndisc*L) )/(rbeam[1]*znorm)
//    - alpar.Dphit - alpar.Dphi0 + ( (alpar.Dxt + alpar.Dx0)*sin(theta) - (alpar.Dyt + alpar.Dy0)*cos(theta) )/rbeam[1];

//   alpar.DthetaB = (Syzk*ndisc-Syk*Sz)*L/term1/rbeam + alpar.Dxt*sin(theta)/rbeam - alpar.Dyt*cos(theta)/rbeam - alpar.Dphit - alpar.DthetaA;
  
   if(print > 3)alpar.print();
  
  return;
}

  //! Reconstruction of parameters for 1 TEC ring, that allows to mask out modules
void TECReconstructor::reconstruct_tec_ring(const LASGlobalData<double>& dphi, const LASGlobalData<int>& mask, TecRingPar& alpar, LAS::beam_group bg, int print)
{
  if (bg != LAS::TEC_PLUS_R4 && bg != LAS::TEC_PLUS_R6 && bg != LAS::TEC_MINUS_R4 && bg != LAS::TEC_MINUS_R6){
    std::cerr << "Error in reconstruct_tec_ring: beam group has invalid value. bg: " << bg << std::endl;
    return;
  }
  
  LASGlobalData<double> dphi2 = dphi; // A copy that can be modified
  
  set_masked_modules_tec(dphi2, mask, alpar, bg, print);
  // Reconstruct TEC parameters iteratively until the masked modules happen to lie on the fit
  int n_iter = 0;
  do{
    if(print > 1)std::cout << "Iteration " << n_iter << std::endl;
    reconstruct_tec_set_onering(dphi2, alpar, bg, print);
  }while(n_iter++ < 50 && set_masked_modules_tec(dphi2, mask, alpar, bg, print) != 0);
  
  if(print > 0) std::cout << "Finished after " << (n_iter) << " iterations" << std:: endl;
  if(print > 1) alpar.print();
}

//! Reconstruction of parameters for TEC, that allows to mask out modules
// The data is supposed to be given in mm
bool TECReconstructor::reconstruct_tec_full(const LASGlobalData<double>& data, const LASGlobalData<int>& mask, TecPar& alpar, LAS::beam_group bg, int print)
{
  if (bg != LAS::TEC_PLUS && bg != LAS::TEC_MINUS){
    std::cerr << "Error in reconstruct_tec_full: beam group has invalid value. bg: " << bg << std::endl;
    return false;
  }
  
  alpar *= 0;

  LASGlobalData<double> data2 = data; // A copy that can be modified
  
  set_masked_modules_tec(data2, mask, alpar, bg, print);
  // Reconstruct TEC parameters iteratively until the masked modules happen to lie on the fit
  int& n_iter =  (bg == LAS::TEC_PLUS ? nr_iter_tecp : nr_iter_tecm);

  n_iter = 0;
  int max_iter = 200;
  do{
    if(print > 1)std::cout << "Iteration " << n_iter << std::endl;
    reconstruct_tec_set_full(data2, alpar, bg, print);
  }while(n_iter++ < max_iter && set_masked_modules_tec(data2, mask, alpar, bg, print) != 0);
  if(n_iter > max_iter/2 && print > 1) std::cerr << "Warning, number of iterations in full reconstruction is " << n_iter << std::endl;  

  // If fit did not converge, set parameters to 0
  if(n_iter > max_iter){
    alpar *= 0;
    return false;
  }

  if(print > 0) std::cout << "Finished after " << (n_iter) << " iterations" << std:: endl;
  if(print > 1) alpar.print();
  return true;
}

//! Reconstruction of parameters for both TECs, that allows to mask out modules
bool TECReconstructor::reconstruct_tec_full(const LASGlobalData<double>& data, const LASGlobalData<int>& mask, LasAlPar& alpar, int print)
{
  if(print > 1){
    std::cout << "Sum of squares before fit:" << std::endl;
    print_sum_sq(data, mask);
  }


  bool retflag = reconstruct_tec_full(data, mask, alpar.tecp, LAS::TEC_PLUS, print);
  retflag = retflag && reconstruct_tec_full(data, mask, alpar.tecm, LAS::TEC_MINUS, print);

  if(print > 1){
    std::cout << "Sum of squares before fit:" << std::endl;
    print_sum_sq(data, mask);
    std::cout << "Sum of squares after fit:" << std::endl;
    print_sum_sq(data - dx_rec_tec(alpar), mask);
  }
  return retflag;
}

void TECReconstructor::compute_residuals(const LASGlobalData<double>& data, const LasAlPar& alpar, LASGlobalData<double>& residuals)
{
  LASGlobalDataLoop loop_tecp(LASGlobalDataLoop::TEC_PLUS);
  do{
    double reco = dx_rec_tec(loop_tecp.get_det(), loop_tecp.get_ring(), loop_tecp.get_beam(), loop_tecp.get_zpos(), alpar.tecp, 0);
    //std::cout <<  GetModuleName(loop_tecp_r4.get_det(), loop_tecp_r4.get_ring(), loop_tecp_r4.get_beam(), loop_tecp_r4.get_zpos(), 0, 1) << ": " << reco << " - " << loop_tecp_r4.GetEntry<double>(data) << " = " << (reco - loop_tecp_r4.GetEntry<double>(data)) << std::endl;
    loop_tecp.GetEntry<double>(residuals) = (reco - loop_tecp.GetEntry<double>(data));
  }while( loop_tecp.next());
  

  LASGlobalDataLoop loop_tecm(LASGlobalDataLoop::TEC_MINUS);
  do{
    double reco = dx_rec_tec(loop_tecm.get_det(), loop_tecm.get_ring(), loop_tecm.get_beam(), loop_tecm.get_zpos(), alpar.tecm, 0);
    //std::cout <<  GetModuleName(loop_tecp_r4.get_det(), loop_tecp_r4.get_ring(), loop_tecp_r4.get_beam(), loop_tecp_r4.get_zpos(), 0, 1) << ": " << reco << " - " << loop_tecp_r4.GetEntry<double>(data) << " = " << (reco - loop_tecp_r4.GetEntry<double>(data)) << std::endl;
    loop_tecm.GetEntry<double>(residuals) = (reco - loop_tecm.GetEntry<double>(data));
  }while( loop_tecm.next());
  

//   LASGlobalDataLoop loop_tecp_r4(LASGlobalDataLoop::TEC_PLUS_R4);
//   do{
//     double reco = dphi_rec_tec(loop_tecp_r4.get_det(), loop_tecp_r4.get_ring(), loop_tecp_r4.get_beam(), loop_tecp_r4.get_zpos(), alpar.tecp_r4, 0);
//     //std::cout <<  GetModuleName(loop_tecp_r4.get_det(), loop_tecp_r4.get_ring(), loop_tecp_r4.get_beam(), loop_tecp_r4.get_zpos(), 0, 1) << ": " << reco << " - " << loop_tecp_r4.GetEntry<double>(data) << " = " << (reco - loop_tecp_r4.GetEntry<double>(data)) << std::endl;
//     loop_tecp_r4.GetEntry<double>(residuals) = (reco - loop_tecp_r4.GetEntry<double>(data));
//   }while( loop_tecp_r4.next());
  
//   LASGlobalDataLoop loop_tecp_r6(LASGlobalDataLoop::TEC_PLUS_R6);
//   do{
//     double reco = dphi_rec_tec(loop_tecp_r6.get_det(), loop_tecp_r6.get_ring(), loop_tecp_r6.get_beam(), loop_tecp_r6.get_zpos(), alpar.tecp_r6, 0);
//     //std::cout <<  GetModuleName(loop_tecp_r6.get_det(), loop_tecp_r6.get_ring(), loop_tecp_r6.get_beam(), loop_tecp_r6.get_zpos(), 0, 1) << ": " << reco << " - " << loop_tecp_r6.GetEntry<double>(data) << " = " << (reco - loop_tecp_r6.GetEntry<double>(data)) << std::endl;
//     loop_tecp_r6.GetEntry<double>(residuals) = (reco - loop_tecp_r6.GetEntry<double>(data));
//   }while( loop_tecp_r6.next());
  
//   LASGlobalDataLoop loop_tecm_r4(LASGlobalDataLoop::TEC_MINUS_R4);
//   do{
//     double reco = dphi_rec_tec(loop_tecm_r4.get_det(), loop_tecm_r4.get_ring(), loop_tecm_r4.get_beam(), loop_tecm_r4.get_zpos(), alpar.tecm_r4, 0);
//     //std::cout <<  GetModuleName(loop_tecm_r4.get_det(), loop_tecm_r4.get_ring(), loop_tecm_r4.get_beam(), loop_tecm_r4.get_zpos(), 0, 1) << ": " << reco << " - " << loop_tecm_r4.GetEntry<double>(data) << " = " << (reco - loop_tecm_r4.GetEntry<double>(data)) << std::endl;
//     loop_tecm_r4.GetEntry<double>(residuals) = (reco - loop_tecm_r4.GetEntry<double>(data));
//   }while( loop_tecm_r4.next());
  
//   LASGlobalDataLoop loop_tecm_r6(LASGlobalDataLoop::TEC_MINUS_R6);
//   do{
//     double reco = dphi_rec_tec(loop_tecm_r6.get_det(), loop_tecm_r6.get_ring(), loop_tecm_r6.get_beam(), loop_tecm_r6.get_zpos(), alpar.tecm_r6, 0);
//     //std::cout <<  GetModuleName(loop_tecm_r6.get_det(), loop_tecm_r6.get_ring(), loop_tecm_r6.get_beam(), loop_tecm_r6.get_zpos(), 0, 1) << ": " << reco << " - " << loop_tecm_r6.GetEntry<double>(data) << " = " << (reco - loop_tecm_r6.GetEntry<double>(data)) << std::endl;
//     loop_tecm_r6.GetEntry<double>(residuals) = (reco - loop_tecm_r6.GetEntry<double>(data));
//   }while( loop_tecm_r6.next());
  
}

