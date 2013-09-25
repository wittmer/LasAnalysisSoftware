#include <sstream>
#include <iostream>
#include <ctime>
#include <iterator>
#include <algorithm>

#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>


#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TROOT.h"

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
#include "LAS_control_plots.h"
#include "LAS_Tec_Reconstruction.h"
#include "LAS_At_rec.h"

double global_data_sum(const LASGlobalData<double>& gd, LAS::beam_group beam_group)
{
  double retval = 0;
  LASGlobalDataLoop lp(convert_to_loop_type(beam_group));
  do{retval += lp.GetEntry(gd);}while(lp.next());
  return retval;
}


// Generalized function to perform the parameter fit
// Supply positions, weights, mask and a vector with coefficients to construct a matrix and a vector
// The matrix is inverted and multiplied with the vector, yielding a vector with the alignment parameters.
void LAS_Alpar_fit(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, const std::vector< LASGlobalData<double> >& coeff_list, std::vector<double>& par_list, TMatrixT<double>& U)
{
  unsigned int npar = coeff_list.size();
  if(npar == 0) throw LAS::Exception("Error in LAS_AlPar_fit: The list of coefficients is empty!");

  LASGlobalData<double> weight = LASGlobalData<double>(1.0) * mask / err2;
  remove_nan(weight);

  //// Weight for covariance matrix
  LASGlobalData<double> moderr2 = LASGlobalData<double>(1.0) * mask * err2;

  // Construct the matrix and vector
  TMatrixT<double> h(npar, npar);
  Avec beta;
  
  for(unsigned int ix = 0; ix < npar; ++ix){
    LASGlobalData<double> term1 = coeff_list[ix] * weight;
    LASGlobalData<double> vector_entry = dif * term1;
    beta.push_back(global_data_sum(vector_entry, LAS::ALL));  // Construct Vector
    for(unsigned int iy = 0; iy < npar; ++iy){
      LASGlobalData<double> matrix_entry = coeff_list[iy] * term1;
      h(ix,iy) = global_data_sum(matrix_entry);
    }
  }

  h.Invert();

  // Multiply the inverted matrix with the vector
  par_list.clear();
  for(unsigned int ip = 0; ip < npar; ++ip){
    double alpar = 0;
    for(unsigned int ix = 0; ix < npar; ++ix){
      alpar += h(ip,ix) * beta[ix];
    }
    par_list.push_back(alpar);
  }
  // Start to calculate the covariance matrix

  U.ResizeTo(npar, npar);

  std::vector<LASGlobalData<double> > cp(npar);

  for(unsigned int ip = 0; ip < npar; ++ip){
    cp[ip] = LASGlobalData<double>(0.);
    for(unsigned int ix = 0; ix < npar; ++ix){
      cp[ip] += coeff_list[ix] * h(ip,ix);
    }
  }

  for(unsigned int ip = 0; ip < npar; ++ip){
    for(unsigned int ix = 0; ix < npar; ++ix){
      //LASGlobalData<double> Umatrix_entry = cp[ip] * cp[ix] * moderr2;
      LASGlobalData<double> Umatrix_entry = cp[ip] * cp[ix] * weight;
      U(ip,ix) = global_data_sum(Umatrix_entry);
    }
  }

  // Finish to calculate the covariance matrix


}

// Create Coefficients for TIB beam parameters
void create_coefficients_TIB_notorsion( std::vector< LASGlobalData<double> >& coeff_list)
{ 
  // TIB translations
  coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TIB));
  coeff_list.push_back(cos(theta()) / r0() * group_mask(LAS::TIB));

  // TIB rotations
  coeff_list.push_back(cos(theta()) / r0() * zpos() * group_mask(LAS::TIB));
  coeff_list.push_back(sin(theta()) / r0() * zpos() * group_mask(LAS::TIB));
  coeff_list.push_back(group_mask(LAS::TIB));

}

// Create Coefficients for TIB beam parameters
void create_coefficients_TIB( std::vector< LASGlobalData<double> >& coeff_list)
{
  static bool initialized = false;
  static std::vector< LASGlobalData<double> > local_coeff_list;

  // Initialize the static list at the first call
  if( ! initialized){

    // TIB translations in mm
    local_coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TIB));
    local_coeff_list.push_back( cos(theta()) / r0() * group_mask(LAS::TIB));
    
    // TIB rotations in rad
    local_coeff_list.push_back(cos(theta()) / r0() * zpos() * group_mask(LAS::TIB));
    local_coeff_list.push_back(sin(theta()) / r0() * zpos() * group_mask(LAS::TIB));
    local_coeff_list.push_back(group_mask(LAS::TIB) );
    
    // TIB torsion in rad / mm 
    local_coeff_list.push_back(zpos() * group_mask(LAS::TIB));

    initialized = true;
  }
  std::back_insert_iterator< std::vector< LASGlobalData<double> > > ii(coeff_list);
  std::copy(local_coeff_list.begin(), local_coeff_list.end(), ii);
}

// Create Coefficients for AT beam parameters
void create_coefficients_AT_beams_pointing( std::vector< LASGlobalData<double> >& coeff_list, LAS::beam_group beam_group = LAS::AT)
{
  for(int beam_nr = 0; beam_nr < 8; ++beam_nr){
    LASGlobalData<double> at_beam_mask = group_mask(beam_group); // beam_group entries are 1 , the rest is 0
    // Mask everything but one beam
    LASGlobalDataLoop lp(LASGlobalDataLoop::AT);
    do{if(lp.get_beam() != beam_nr)lp.GetEntry(at_beam_mask)=0; }while(lp.next());

    coeff_list.push_back(at_beam_mask * (zpos() - LAS::z_at_bs)); // Coefficient for the slope
    coeff_list.push_back(at_beam_mask); // Coefficient for the offset
  }
}

// Create Coefficients for AT beam parameters
void create_coefficients_AT_beams( std::vector< LASGlobalData<double> >& coeff_list, LAS::beam_group beam_group)
{
  static bool initialized = false;
  static std::vector< LASGlobalData<double> > local_coeff_list;

  // Initialize the static list at the first call
  if( ! initialized){
    for(int beam_nr = 0; beam_nr < 8; ++beam_nr){
      LASGlobalData<double> at_beam_mask = group_mask(LAS::AT); // beam_group entries are 1 , the rest is 0
      // Mask everything but one beam
      LASGlobalDataLoop lp(LASGlobalDataLoop::AT);
      do{if(lp.get_beam() != beam_nr)lp.GetEntry(at_beam_mask)=0; }while(lp.next());
      
      local_coeff_list.push_back(at_beam_mask * (zpos() - LAS::z_at_bs) * LAS::r_at_tec / r0()); // Coefficient for the slope
      local_coeff_list.push_back(at_beam_mask * LAS::r_at_tec / r0()); // Coefficient for the offset
    }
    initialized = true;
  }

  std::back_insert_iterator< std::vector< LASGlobalData<double> > > ii(coeff_list);
  std::copy(local_coeff_list.begin(), local_coeff_list.end(), ii);
}


// Create Coefficients for TEC disc parameters
void create_coefficients_TEC_discs( std::vector< LASGlobalData<double> >& coeff_list, LAS::beam_group beam_group)
{
  static bool initialized = false;

  static std::vector< LASGlobalData<double> > tecm_coeff_list;
  static std::vector< LASGlobalData<double> > tecp_coeff_list;
  // Initialize the static lists at the first call
  if( ! initialized){
    // Disc parameters TEC
    for(int disc_nr = 0; disc_nr < 9; ++disc_nr){
      LASGlobalData<double> disc_mask(0);
      // Mask everything but one disc
      LASGlobalDataLoop lp1(LASGlobalDataLoop::TEC);
      do{if(lp1.get_zpos() == disc_nr)lp1.GetEntry(disc_mask) = 1; }while(lp1.next());
      LASGlobalDataLoop lp2(LASGlobalDataLoop::TEC_AT);
      do{if(lp2.get_zpos() == disc_nr)lp2.GetEntry(disc_mask) = 1; }while(lp2.next());

      // TEC-
      tecm_coeff_list.push_back( disc_mask * ( group_mask(LAS::TEC_MINUS) || group_mask(LAS::TEC_MINUS_AT) ) ); //Disc rotation
      tecm_coeff_list.push_back( -sin(theta()) / r0() * disc_mask * ( group_mask(LAS::TEC_MINUS) || group_mask(LAS::TEC_MINUS_AT) ) ); // Disc translation x
      tecm_coeff_list.push_back(  cos(theta()) / r0() * disc_mask * ( group_mask(LAS::TEC_MINUS) || group_mask(LAS::TEC_MINUS_AT) ) ); // Disc translation y

      // TEC+
      tecp_coeff_list.push_back( disc_mask * ( group_mask(LAS::TEC_PLUS) || group_mask(LAS::TEC_PLUS_AT) ) ); //Disc rotation
      tecp_coeff_list.push_back( -sin(theta()) / r0() * disc_mask * ( group_mask(LAS::TEC_PLUS) || group_mask(LAS::TEC_PLUS_AT) ) ); // Disc translation x
      tecp_coeff_list.push_back(  cos(theta()) / r0() * disc_mask * ( group_mask(LAS::TEC_PLUS) || group_mask(LAS::TEC_PLUS_AT) ) ); // Disc translation y

    }
    initialized = true;
  }

  //Fill the coefficients into the list, according to beam_group

  // Get an insertion iterator to fill the coefficient list (add to the back)
  std::back_insert_iterator< std::vector< LASGlobalData<double> > > ii(coeff_list);

  // TEC-
  if( beam_group == LAS::TEC_MINUS_R4 || beam_group == LAS::TEC_MINUS_R6 || beam_group == LAS::TEC_MINUS  || beam_group == LAS::TEC || beam_group == LAS::ALL)
    std::copy(tecm_coeff_list.begin(), tecm_coeff_list.end(), ii);

  // TEC+
  if( beam_group == LAS::TEC_PLUS_R4 || beam_group == LAS::TEC_PLUS_R6 || beam_group == LAS::TEC_PLUS  || beam_group == LAS::TEC || beam_group == LAS::ALL)
    std::copy(tecp_coeff_list.begin(), tecp_coeff_list.end(), ii);
}

// Create Coefficients for TEC beam parameters
void create_coefficients_TEC_beams( std::vector< LASGlobalData<double> >& coeff_list, LAS::beam_group beam_group)
{
  static bool initialized = false;

  static std::vector< LASGlobalData<double> > tecm_r4_coeff_list;
  static std::vector< LASGlobalData<double> > tecm_r6_coeff_list;
  static std::vector< LASGlobalData<double> > tecp_r4_coeff_list;
  static std::vector< LASGlobalData<double> > tecp_r6_coeff_list;

  // Initialize the static lists at the first call
  if( ! initialized){
    for(int beam_nr = 0; beam_nr < 8; ++beam_nr){

      LASGlobalData<double> beam_mask(0);
      // Mask everything but one beam
      LASGlobalDataLoop lp(LASGlobalDataLoop::TEC);
      do{if(lp.get_beam() == beam_nr)lp.GetEntry(beam_mask) = 1; }while(lp.next());

      // TEC- Ring 4
      tecm_r4_coeff_list.push_back(beam_mask * group_mask(LAS::TEC_MINUS_R4) * (zpos() + LAS::z_tec_bs + LAS::z_disc0)); // Coefficient for the slope
      tecm_r4_coeff_list.push_back(beam_mask * group_mask(LAS::TEC_MINUS_R4)); // Coefficient for the offset
      // TEC- Ring 6
      tecm_r6_coeff_list.push_back(beam_mask * group_mask(LAS::TEC_MINUS_R6) * (zpos() + LAS::z_tec_bs + LAS::z_disc0)); // Coefficient for the slope
      tecm_r6_coeff_list.push_back(beam_mask * group_mask(LAS::TEC_MINUS_R6)); // Coefficient for the offset

      // TEC+ Ring 4
      tecp_r4_coeff_list.push_back(beam_mask * group_mask(LAS::TEC_PLUS_R4) * (zpos() - LAS::z_tec_bs - LAS::z_disc0)); // Coefficient for the slope
      tecp_r4_coeff_list.push_back(beam_mask * group_mask(LAS::TEC_PLUS_R4)); // Coefficient for the offset
      // TEC+ Ring 6
      tecp_r6_coeff_list.push_back(beam_mask * group_mask(LAS::TEC_PLUS_R6) * (zpos() - LAS::z_tec_bs - LAS::z_disc0)); // Coefficient for the slope
      tecp_r6_coeff_list.push_back(beam_mask * group_mask(LAS::TEC_PLUS_R6)); // Coefficient for the offset
    }
    initialized = true;
  }

  //Fill the coefficients into the list, according to beam_group

  // Get an insertion iterator to fill the coefficient list (add to the back)
  std::back_insert_iterator< std::vector< LASGlobalData<double> > > ii(coeff_list);

  // TEC- Ring 4
  if( beam_group == LAS::TEC_MINUS_R4 || beam_group == LAS::TEC_MINUS  || beam_group == LAS::TEC || beam_group == LAS::ALL)
    std::copy(tecm_r4_coeff_list.begin(), tecm_r4_coeff_list.end(), ii);

  // TEC- Ring 6
  if( beam_group == LAS::TEC_MINUS_R6 || beam_group == LAS::TEC_MINUS  || beam_group == LAS::TEC || beam_group == LAS::ALL)
    std::copy(tecm_r6_coeff_list.begin(), tecm_r6_coeff_list.end(), ii);

  // TEC+ Ring 4
  if( beam_group == LAS::TEC_PLUS_R4 || beam_group == LAS::TEC_PLUS  || beam_group == LAS::TEC || beam_group == LAS::ALL)
    std::copy(tecp_r4_coeff_list.begin(), tecp_r4_coeff_list.end(), ii);

  // TEC+ Ring 6
  if( beam_group == LAS::TEC_PLUS_R6 || beam_group == LAS::TEC_PLUS  || beam_group == LAS::TEC || beam_group == LAS::ALL)
    std::copy(tecp_r6_coeff_list.begin(), tecp_r6_coeff_list.end(), ii);

}


void Calc_chi2(const LASGlobalData<double>& data, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, int ndof)
{
  // Calculate chi2
  LASGlobalData<double> chi = data * data / err2 * mask;
  remove_nan(chi);

  double tob_chi  = global_data_sum(chi * mask, LAS::TOB) / global_data_sum(mask && group_mask(LAS::TOB));
  std::cout << "TOB Chi2/ntob: " << tob_chi << std::endl;

  double tib_chi  = global_data_sum(chi * mask, LAS::TIB) / global_data_sum(mask && group_mask(LAS::TIB));
  std::cout << "TIB Chi2/ntib: " << tib_chi << std::endl;

  double tecm_chi  = global_data_sum(chi * mask, LAS::TEC_MINUS_AT) / global_data_sum(mask && group_mask(LAS::TEC_MINUS_AT));
  std::cout << "TEC- (AT) Chi2/ntecm: " << tecm_chi << std::endl;

  double tecp_chi  = global_data_sum(chi * mask, LAS::TEC_PLUS_AT) / global_data_sum(mask && group_mask(LAS::TEC_PLUS_AT));
  std::cout << "TEC+ (AT) Chi2/ntecp: " << tecp_chi << std::endl;

  double at_chi  = global_data_sum(chi * (mask && (group_mask(LAS::TIB) || group_mask(LAS::TOB) || group_mask(LAS::TEC_PLUS_AT) || group_mask(LAS::TEC_MINUS_AT))), LAS::AT) / ndof;
  std::cout << "AT Chi2/ndof: " << at_chi << std::endl;
}

unsigned int copy_at_beampar( LasAlPar& results, const std::vector<double>& alpar, const TMatrixT<double>& covariance_matrix, unsigned int offset = 0 )
{
  for(int bn = 0; bn < 8; ++bn){
    results.at.er_beam_a[bn]  = sqrt( covariance_matrix(offset, offset));
    results.at.beam_a[bn] = alpar[offset++];

    results.at.er_beam_b[bn]  = sqrt( covariance_matrix(offset, offset));
    results.at.beam_b[bn] = alpar[offset++];
  }
  return offset;
}

unsigned int copy_tib_alpar( LasAlPar& results, const std::vector<double>& alpar, const TMatrixT<double>& covariance_matrix, unsigned int offset = 0 )
{
  results.at.tib_dx  = alpar[offset+0];
  results.at.tib_dy  = alpar[offset+1];
  results.at.tib_rx  = alpar[offset+2];
  results.at.tib_ry  = alpar[offset+3];
  results.at.tib_rz  = alpar[offset+4];
  results.at.tib_tz  = alpar[offset+5];
  results.at.er_tib_dx  = sqrt( covariance_matrix(offset + 0, offset + 0));
  results.at.er_tib_dy  = sqrt( covariance_matrix(offset + 1, offset + 1));
  results.at.er_tib_rx  = sqrt( covariance_matrix(offset + 2, offset + 2));
  results.at.er_tib_ry  = sqrt( covariance_matrix(offset + 3, offset + 3));
  results.at.er_tib_rz  = sqrt( covariance_matrix(offset + 4, offset + 4));
  results.at.er_tib_tz  = sqrt( covariance_matrix(offset + 5, offset + 5));
  return offset + 6;
}

unsigned int copy_tecp_at_alpar( LasAlPar& results, const std::vector<double>& alpar, const TMatrixT<double>& covariance_matrix, unsigned int offset = 0 )
{
  results.at.tecp_dx  = alpar[offset+0];
  results.at.tecp_dy  = alpar[offset+1];
  results.at.tecp_rx  = alpar[offset+2];
  results.at.tecp_ry  = alpar[offset+3];
  results.at.tecp_rz  = alpar[offset+4];
  results.at.tecp_tz  = alpar[offset+5];
  results.at.er_tecp_dx  = sqrt( covariance_matrix(offset + 0, offset + 0));
  results.at.er_tecp_dy  = sqrt( covariance_matrix(offset + 1, offset + 1));
  results.at.er_tecp_rx  = sqrt( covariance_matrix(offset + 2, offset + 2));
  results.at.er_tecp_ry  = sqrt( covariance_matrix(offset + 3, offset + 3));
  results.at.er_tecp_rz  = sqrt( covariance_matrix(offset + 4, offset + 4));
  results.at.er_tecp_tz  = sqrt( covariance_matrix(offset + 5, offset + 5));
  return offset + 6;
}

unsigned int copy_tecm_at_alpar( LasAlPar& results, const std::vector<double>& alpar, const TMatrixT<double>& covariance_matrix, unsigned int offset = 0 )
{
  results.at.tecm_dx  = alpar[offset+0];
  results.at.tecm_dy  = alpar[offset+1];
  results.at.tecm_rx  = alpar[offset+2];
  results.at.tecm_ry  = alpar[offset+3];
  results.at.tecm_rz  = alpar[offset+4];
  results.at.tecm_tz  = alpar[offset+5];
  results.at.er_tecm_dx  = sqrt( covariance_matrix(offset + 0, offset + 0));
  results.at.er_tecm_dy  = sqrt( covariance_matrix(offset + 1, offset + 1));
  results.at.er_tecm_rx  = sqrt( covariance_matrix(offset + 2, offset + 2));
  results.at.er_tecm_ry  = sqrt( covariance_matrix(offset + 3, offset + 3));
  results.at.er_tecm_rz  = sqrt( covariance_matrix(offset + 4, offset + 4));
  results.at.er_tecm_tz  = sqrt( covariance_matrix(offset + 5, offset + 5));
  return offset + 6;
}

unsigned int copy_tec_beampar( LasAlPar& results, const std::vector<double>& alpar, const TMatrixT<double>& covariance_matrix, unsigned int offset = 0, LAS::beam_group beam_group = LAS::ALL )
{
  // TEC- Ring 4
  if( beam_group == LAS::TEC_MINUS_R4 || beam_group == LAS::TEC_MINUS  || beam_group == LAS::TEC || beam_group == LAS::ALL){
    for(int bn = 0; bn < 8; ++bn){
      results.tecm.DthetaA_r4_er[bn]  = sqrt( covariance_matrix(offset, offset));
      results.tecm.DthetaA_r4[bn] = alpar[offset++];
      results.tecm.DthetaB_r4_er[bn]  = sqrt( covariance_matrix(offset, offset));
      results.tecm.DthetaB_r4[bn] = alpar[offset++];
    }
  }
  // TEC- Ring 6
  if( beam_group == LAS::TEC_MINUS_R6 || beam_group == LAS::TEC_MINUS  || beam_group == LAS::TEC || beam_group == LAS::ALL){
    for(int bn = 0; bn < 8; ++bn){
      results.tecm.DthetaA_r6_er[bn]  = sqrt( covariance_matrix(offset, offset));
      results.tecm.DthetaA_r6[bn] = alpar[offset++];
      results.tecm.DthetaB_r6_er[bn]  = sqrt( covariance_matrix(offset, offset));
      results.tecm.DthetaB_r6[bn] = alpar[offset++];
    }
  }
  
  // TEC+ Ring 4
  if( beam_group == LAS::TEC_PLUS_R4 || beam_group == LAS::TEC_PLUS  || beam_group == LAS::TEC || beam_group == LAS::ALL){
    for(int bn = 0; bn < 8; ++bn){
      results.tecp.DthetaA_r4_er[bn]  = sqrt( covariance_matrix(offset, offset));
      results.tecp.DthetaA_r4[bn] = alpar[offset++];
      results.tecp.DthetaB_r4_er[bn]  = sqrt( covariance_matrix(offset, offset));
      results.tecp.DthetaB_r4[bn] = alpar[offset++];
    }
  }
  // TEC+ Ring 6
  if( beam_group == LAS::TEC_PLUS_R6 || beam_group == LAS::TEC_PLUS  || beam_group == LAS::TEC || beam_group == LAS::ALL){
    for(int bn = 0; bn < 8; ++bn){
      results.tecp.DthetaA_r6_er[bn]  = sqrt( covariance_matrix(offset, offset));
      results.tecp.DthetaA_r6[bn] = alpar[offset++];
      results.tecp.DthetaB_r6_er[bn]  = sqrt( covariance_matrix(offset, offset));
      results.tecp.DthetaB_r6[bn] = alpar[offset++];
    }
  }
  return offset;
}

unsigned int copy_tec_discpar( LasAlPar& results, const std::vector<double>& alpar, const TMatrixT<double>& covariance_matrix, unsigned int offset = 0, LAS::beam_group beam_group = LAS::ALL )
{

  // TEC-
  if( beam_group == LAS::TEC_MINUS_R4 || beam_group == LAS::TEC_MINUS_R6 || beam_group == LAS::TEC_MINUS  || beam_group == LAS::TEC || beam_group == LAS::ALL){
    for(int disc_nr = 0; disc_nr < 9; ++disc_nr){
      results.tecm.Dphik_er[disc_nr]  = sqrt( covariance_matrix(offset, offset));
      results.tecm.Dphik[disc_nr]  = alpar[offset++];
      results.tecm.Dxk_er[disc_nr]  = sqrt( covariance_matrix(offset, offset));
      results.tecm.Dxk[disc_nr]  = alpar[offset++];
      results.tecm.Dyk_er[disc_nr]  = sqrt( covariance_matrix(offset, offset));
      results.tecm.Dyk[disc_nr]  = alpar[offset++];
    }
  }

  // TEC+
  if( beam_group == LAS::TEC_PLUS_R4 || beam_group == LAS::TEC_PLUS_R6 || beam_group == LAS::TEC_PLUS  || beam_group == LAS::TEC || beam_group == LAS::ALL){
    for(int disc_nr = 0; disc_nr < 9; ++disc_nr){
      results.tecp.Dphik_er[disc_nr]  = sqrt( covariance_matrix(offset, offset));
      results.tecp.Dphik[disc_nr]  = alpar[offset++];
      results.tecp.Dxk_er[disc_nr]  = sqrt( covariance_matrix(offset, offset));
      results.tecp.Dxk[disc_nr]  = alpar[offset++];
      results.tecp.Dyk_er[disc_nr]  = sqrt( covariance_matrix(offset, offset));
      results.tecp.Dyk[disc_nr]  = alpar[offset++];
    }
  }
  return offset;
}


void Init_fit_type( std::vector< LASGlobalData<double> >& coeff_list, LASGlobalData<int>& gr_mask, LAS::FIT_TYPE::fit_type fit_type)
{
  switch(fit_type){
  case LAS::FIT_TYPE::ALL :
    create_coefficients_AT_beams(coeff_list);
    create_coefficients_TIB(coeff_list);
    create_coefficients_TEC_beams(coeff_list);
    create_coefficients_TEC_discs(coeff_list);
    break;
  case LAS::FIT_TYPE::TEC_MINUS :
    create_coefficients_TEC_beams(coeff_list, LAS::TEC_MINUS);
    gr_mask = group_mask(LAS::TEC_MINUS);
    break;
  case LAS::FIT_TYPE::TEC_MINUS_R4 :
    create_coefficients_TEC_beams(coeff_list, LAS::TEC_MINUS_R4);
    gr_mask = group_mask(LAS::TEC_MINUS_R4);
    break;
  case LAS::FIT_TYPE::TEC_MINUS_R6 :
    create_coefficients_TEC_beams(coeff_list, LAS::TEC_MINUS_R6);
    gr_mask = group_mask(LAS::TEC_MINUS_R6);
    break;
  case LAS::FIT_TYPE::TEC_PLUS :
    create_coefficients_TEC_beams(coeff_list, LAS::TEC_PLUS);
    gr_mask = group_mask(LAS::TEC_PLUS);
    break;
  case LAS::FIT_TYPE::TEC_PLUS_R4 :
    create_coefficients_TEC_beams(coeff_list, LAS::TEC_PLUS_R4);
    gr_mask = group_mask(LAS::TEC_PLUS_R4);
    break;
  case LAS::FIT_TYPE::TEC_PLUS_R6 :
    create_coefficients_TEC_beams(coeff_list, LAS::TEC_PLUS_R6);
    gr_mask = group_mask(LAS::TEC_PLUS_R6);
    break;
  case LAS::FIT_TYPE::TOB :
    gr_mask = group_mask(LAS::TOB);
    create_coefficients_AT_beams(coeff_list);
    break;
  case LAS::FIT_TYPE::TIB :
    gr_mask = group_mask(LAS::TOB) || group_mask(LAS::TIB);
    create_coefficients_AT_beams(coeff_list);
    create_coefficients_TIB(coeff_list);
    break;
  case LAS::FIT_TYPE::TEC_PLUS_AT :
    gr_mask = group_mask(LAS::TOB) || group_mask(LAS::TEC_PLUS_AT);
    create_coefficients_AT_beams(coeff_list);
    // AT TEC+ translations
    coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TEC_PLUS_AT));
    coeff_list.push_back( cos(theta()) / r0() * group_mask(LAS::TEC_PLUS_AT));

    // TEC+ rotation
    coeff_list.push_back(cos(theta()) / r0() * zpos() * group_mask(LAS::TEC_PLUS_AT));
    coeff_list.push_back(sin(theta()) / r0() * zpos() * group_mask(LAS::TEC_PLUS_AT));
    coeff_list.push_back(group_mask(LAS::TEC_PLUS_AT));

    // TEC+ torsion
    coeff_list.push_back((zpos() - LAS::z_disc0) * group_mask(LAS::TEC_PLUS_AT));

    break;
  case LAS::FIT_TYPE::TEC_MINUS_AT :
    gr_mask = group_mask(LAS::TOB) || group_mask(LAS::TEC_MINUS_AT);
    create_coefficients_AT_beams(coeff_list);
    // AT TEC- translations
    coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TEC_MINUS_AT));
    coeff_list.push_back( cos(theta()) / r0() * group_mask(LAS::TEC_MINUS_AT));

    // TEC- rotation
    coeff_list.push_back(cos(theta()) / r0() * zpos() * group_mask(LAS::TEC_MINUS_AT));
    coeff_list.push_back(sin(theta()) / r0() * zpos() * group_mask(LAS::TEC_MINUS_AT));
    coeff_list.push_back(group_mask(LAS::TEC_MINUS_AT));

    // TEC- torsion
    coeff_list.push_back(zpos() * group_mask(LAS::TEC_MINUS_AT));

    break;
  case LAS::FIT_TYPE::AT :
    gr_mask = group_mask(LAS::TOB) || group_mask(LAS::TIB) || group_mask(LAS::TEC_MINUS_AT) || group_mask(LAS::TEC_PLUS_AT);
    create_coefficients_AT_beams(coeff_list);
    create_coefficients_TIB(coeff_list);

    // AT TEC- translations
    coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TEC_MINUS_AT));
    coeff_list.push_back( cos(theta()) / r0() * group_mask(LAS::TEC_MINUS_AT));

    // TEC- rotation
    coeff_list.push_back(cos(theta()) / r0() * (zpos() + LAS::z_disc0) * group_mask(LAS::TEC_MINUS_AT));
    coeff_list.push_back(sin(theta()) / r0() * (zpos() + LAS::z_disc0) * group_mask(LAS::TEC_MINUS_AT));
    coeff_list.push_back(group_mask(LAS::TEC_MINUS_AT));

    // TEC- torsion
    coeff_list.push_back((zpos() + LAS::z_disc0) * group_mask(LAS::TEC_MINUS_AT));

    // AT TEC+ translations
    coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TEC_PLUS_AT));
    coeff_list.push_back( cos(theta()) / r0() * group_mask(LAS::TEC_PLUS_AT));

    // // TEC+ rotation
    coeff_list.push_back(cos(theta()) / r0() * (zpos() - LAS::z_disc0) * group_mask(LAS::TEC_PLUS_AT));
    coeff_list.push_back(sin(theta()) / r0() * (zpos() - LAS::z_disc0) * group_mask(LAS::TEC_PLUS_AT));
    coeff_list.push_back(group_mask(LAS::TEC_PLUS_AT));

    // // TEC+ torsion
    coeff_list.push_back((zpos() - LAS::z_disc0) * group_mask(LAS::TEC_PLUS_AT));

    break;
  default:
    throw LAS::Exception("Error in LAS_fit: this fit type is not implemented yet!");
  }
}


void copy_fit_results( const std::vector<double>& alpar, const TMatrixT<double>& covariance_matrix, LasAlPar& results, LAS::FIT_TYPE::fit_type fit_type)
{
  // Copy the alignment parameters to their locations
  switch(fit_type){
  case LAS::FIT_TYPE::TEC_MINUS :
    for(int bn = 0; bn < 8; ++bn){
      results.tecm.DthetaA_r4[bn] = alpar[4 * bn];
      results.tecm.DthetaB_r4[bn] = alpar[4 * bn + 1];
      results.tecm.DthetaA_r6[bn] = alpar[4 * bn + 2];
      results.tecm.DthetaB_r6[bn] = alpar[4 * bn + 3];
    }
    break;
  case LAS::FIT_TYPE::TEC_MINUS_R4 :
    for(int bn = 0; bn < 8; ++bn){
      results.tecm.DthetaA_r4[bn] = alpar[2*bn];
      results.tecm.DthetaB_r4[bn] = alpar[2*bn+1];
    }
    break;
  case LAS::FIT_TYPE::TEC_MINUS_R6 :
    for(int bn = 0; bn < 8; ++bn){
      results.tecm.DthetaA_r6[bn] = alpar[2*bn];
      results.tecm.DthetaB_r6[bn] = alpar[2*bn+1];
    }
    break;
  case LAS::FIT_TYPE::TEC_PLUS :
    for(int bn = 0; bn < 8; ++bn){
      results.tecp.DthetaA_r4[bn] = alpar[4 * bn];
      results.tecp.DthetaB_r4[bn] = alpar[4 * bn + 1];
      results.tecp.DthetaA_r6[bn] = alpar[4 * bn + 2];
      results.tecp.DthetaB_r6[bn] = alpar[4 * bn + 3];
    }
    break;
  case LAS::FIT_TYPE::TEC_PLUS_R4 :
    for(int bn = 0; bn < 8; ++bn){
      results.tecp.DthetaA_r4[bn] = alpar[2*bn];
      results.tecp.DthetaB_r4[bn] = alpar[2*bn+1];
    }
    break;
  case LAS::FIT_TYPE::TEC_PLUS_R6 :
    for(int bn = 0; bn < 8; ++bn){
      results.tecp.DthetaA_r6[bn] = alpar[2*bn];
      results.tecp.DthetaB_r6[bn] = alpar[2*bn+1];
    }
    break;
  case LAS::FIT_TYPE::TOB :
    copy_at_beampar(results, alpar, covariance_matrix, 0);
    break;
  case LAS::FIT_TYPE::TIB :
    {unsigned int offset = 0;
    offset = copy_at_beampar(results, alpar, covariance_matrix, offset);
    copy_tib_alpar(results, alpar, covariance_matrix, offset);
    break;}
  case LAS::FIT_TYPE::TEC_PLUS_AT :
    {unsigned int offset = 0;
    offset = copy_at_beampar(results, alpar, covariance_matrix, offset);
    copy_tecp_at_alpar(results, alpar, covariance_matrix, offset);
    break;}
  case LAS::FIT_TYPE::TEC_MINUS_AT :
    {unsigned int offset = 0;
    offset = copy_at_beampar(results, alpar, covariance_matrix, offset);
    copy_tecm_at_alpar(results, alpar, covariance_matrix, offset);
    break;}
  case LAS::FIT_TYPE::AT :
    {unsigned int offset = 0;
    copy_at_beampar(results, alpar, covariance_matrix, offset);
    offset = copy_tib_alpar(results, alpar, covariance_matrix, offset);
    offset = copy_tecm_at_alpar(results, alpar, covariance_matrix, offset);
    offset = copy_tecp_at_alpar(results, alpar, covariance_matrix, offset);
    break;}
  case LAS::FIT_TYPE::ALL :
    {unsigned int offset = 0;
    offset = copy_at_beampar(results, alpar, covariance_matrix, offset);
    offset = copy_tib_alpar(results, alpar, covariance_matrix, offset);
    offset = copy_tec_beampar(results, alpar, covariance_matrix, offset);
    copy_tec_discpar(results, alpar, covariance_matrix, offset);
    break;
    }
  default:
    break;
  }
}


void calc_chi2_at(LASGlobalData<double>& chi, LASGlobalData<int>& glob_mask, AtChi& results, LAS::FIT_TYPE::fit_type fit_type, int ndof)
{
  double tob_chi  = global_data_sum(chi * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
  double tib_chi  = global_data_sum(chi * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
  double tecp_chi  = global_data_sum(chi * glob_mask, LAS::TEC_PLUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_PLUS_AT));
  double tecm_chi  = global_data_sum(chi * glob_mask, LAS::TEC_MINUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_MINUS_AT));
  double AT_chi   = global_data_sum(chi * (glob_mask && group_mask(LAS::AT)), LAS::AT) / ndof;
  results.tib  = tib_chi;
  results.tob  = tob_chi;
  results.tecp = tecp_chi;
  results.tecm = tecm_chi;
  results.AT   = AT_chi;
}

// Fit a set of alignment parameters
void LAS_fit(const LASGlobalData<double>& dif_rad, const LASGlobalData<double>& err2_rad2, const LASGlobalData<int>& mask, LasAlPar& results, LASGlobalData<double>& rec_pos, LAS::FIT_TYPE::fit_type fit_type, bool verbose)
{
  std::vector< LASGlobalData<double> > coeff_list;
  LASGlobalData<int> gr_mask(1);

  Init_fit_type(coeff_list, gr_mask, fit_type);

  LASGlobalData<int> glob_mask = gr_mask && mask ;

  // Perform the minimization
  std::vector<double> alpar;
  TMatrixT<double> covariance_matrix;
  LAS_Alpar_fit(dif_rad, err2_rad2, glob_mask, coeff_list, alpar, covariance_matrix);

  // Calculate the reconstructed positions
  rec_pos = LASGlobalData<double>(0);
  for(unsigned int ip = 0; ip < alpar.size(); ++ip){
    rec_pos += coeff_list[ip] * alpar[ip];
  }

  copy_fit_results( alpar, covariance_matrix, results, fit_type);

  // Number of degrees of freedom
  int ndof = global_data_sum(glob_mask) - coeff_list.size();

  // Calculate chi2 before fit
  LASGlobalData<double> chi0 = dif_rad * dif_rad / err2_rad2 * glob_mask;
  remove_nan(chi0 );
  calc_chi2_at(chi0, glob_mask, results.at.chi0, fit_type, ndof);

  // Calculate chi2 after fit
  LASGlobalData<double> chi = (rec_pos - dif_rad) * (rec_pos - dif_rad) / err2_rad2 * glob_mask;
  remove_nan(chi );
  calc_chi2_at(chi, glob_mask, results.at.chi, fit_type, ndof);
}


void AT_Parameter_Reconstruction(const std::string& run_list_name, const std::string& path_prefix, const std::string& outfile, LAS::FIT_TYPE::fit_type fit_type, bool quick_draw)
{
  std::time_t t1 = std::time(NULL);

  // Read in the run list (first entry is reference file) 
  std::string ref_filename;
  std::vector<std::string> file_list;
  get_file_list(run_list_name, ref_filename, file_list, path_prefix);

  // Get the reference positions
  LASGlobalData<double> ref_pos;
  LASGlobalData<double> ref_err;
  LASGlobalData<int> ref_mask;
  get_ref_pos(ref_filename, ref_pos, ref_err, ref_mask);
 
  // Create vectors for resultes
  std::vector<LasAlPar> par_list;
  std::vector<AtPar> par_list_at;
  Avec time_axis;

  // Initialize the Fit coefficients and the module mask
  std::vector< LASGlobalData<double> > coeff_list;
  LASGlobalData<int> glob_mask = ref_mask;
  Init_fit_type(coeff_list, glob_mask, fit_type);
  mask_known_bad_modules(glob_mask);

  // Loop over all runs  
  for(unsigned int i = 0; i < file_list.size(); i++){
    std::cout << "\nProcessing file " << file_list[i] << std::endl;
    TFile file((file_list[i]).c_str(), "READ");

    // Check if the run is labled as good
    Avec good_run = avec_get("good_run", file);
    if(good_run.empty()){
      std::cerr << "Warning: No good_run quality flag in this file, will skip it." << std::endl;
      continue;
    }
    if( good_run[0] == 0){
      std::cerr << "Run is flagged as bad, but I will not skip it!" << std::endl;
    }

    // Get blocks flag and timestamp
    Avec good_block = avec_get("good_block",file);
    Avec block_timestamp = avec_get("block_timestamp", file);
    std::cout << "Nr. of good blocks: " << vsum(good_block) << " out of " << good_block.size() << ", ratio is " << double(vsum(good_block))/double(good_block.size()) << std::endl;
    if(vsum(good_block) < 1)continue;

    // Loop over all blocks in one run
    for(unsigned int blnr = 0; blnr < good_block.size(); blnr++){
      if(good_block[blnr] == 0){
	std::cout << "Bad block Nr." << blnr << std::endl;
	continue;
      }
      std::ostringstream block_nr;
      block_nr << blnr;

      LASGlobalData<double> pos = global_data_get<double>("positions_" + block_nr.str(), file);
      LASGlobalData<double> rms = global_data_get<double>("pos_error_" + block_nr.str(),file);
      LASGlobalData<int> norm = global_data_get<int>("norm_" + block_nr.str(),file);
      LASGlobalData<int> pos_mask = global_data_get<int>("positions_mask_" + block_nr.str(), file);
      rms /= sqrt(norm);
      remove_nan(rms);

      LASGlobalData<double> dif_rad = (pos - ref_pos)*pitch()/r0();
      correct_signs(dif_rad);
      LASGlobalData<double> error_rad = sqrt(rms*rms + ref_err*ref_err)*pitch()/r0();
      LASGlobalData<int> mask = glob_mask && pos_mask && fabs(dif_rad) < 1e-3;

      LasAlPar results;
      LASGlobalData<double> rec_pos_rad;
      LASGlobalData<double> err2_rad2 = error_rad*error_rad;
      
      // Perform the minimization
      std::vector<double> alpar;
      TMatrixT<double> covariance_matrix;
      LAS_Alpar_fit(dif_rad, err2_rad2, mask, coeff_list, alpar, covariance_matrix);
      
      copy_fit_results( alpar, covariance_matrix, results, fit_type);
      
      // Calculate the reconstructed positions
      rec_pos_rad = LASGlobalData<double>(0);
      for(unsigned int ip = 0; ip < alpar.size(); ++ip) rec_pos_rad += coeff_list[ip] * alpar[ip];
      
      // Number of degrees of freedom
      int ndof = global_data_sum(mask) - coeff_list.size();
      
      // Calculate chi2 before fit
      LASGlobalData<double> chi0 = dif_rad * dif_rad / err2_rad2 * mask;
      remove_nan(chi0 );
      calc_chi2_at(chi0, mask, results.at.chi0, fit_type, ndof);
      
      // Calculate chi2 after fit
      LASGlobalData<double> chi = (rec_pos_rad - dif_rad) * (rec_pos_rad - dif_rad) / err2_rad2 * mask;
      remove_nan(chi );
      calc_chi2_at(chi, mask, results.at.chi, fit_type, ndof);

      par_list.push_back(results);
      par_list_at.push_back(results.at);
      time_axis.push_back( block_timestamp[blnr] );
    }
    file.Close();
  } 

  Avec::VERBOSE_FLAG = 0;

  if(fit_type == LAS::FIT_TYPE::TEC_PLUS_R4 || fit_type == LAS::FIT_TYPE::ALL){
    alpar_history_draw(par_list, time_axis , quick_draw);
  }

  if(fit_type == LAS::FIT_TYPE::TOB){
    DrawATPar (par_list_at, time_axis , quick_draw, outfile, LAS::TOB);
  }

  if(fit_type == LAS::FIT_TYPE::TIB){
    DrawATPar (par_list_at, time_axis , quick_draw, outfile, LAS::TIB);
  }

  if(fit_type == LAS::FIT_TYPE::AT || fit_type == LAS::FIT_TYPE::ALL){
    DrawATPar (par_list_at, time_axis , quick_draw, outfile, LAS::AT);
  }

  if(fit_type == LAS::FIT_TYPE::TIB_OLD){
    DrawATPar (par_list_at, time_axis , quick_draw, outfile, LAS::AT);
  }

  if(fit_type == LAS::FIT_TYPE::TEC_MINUS_AT || fit_type == LAS::FIT_TYPE::TEC_PLUS_AT || fit_type == LAS::FIT_TYPE::TEC_AT || fit_type == LAS::FIT_TYPE::ALL){
    DrawATPar(par_list_at, time_axis, quick_draw, outfile, LAS::AT);
  }

  std::time_t t2 = std::time(NULL);
  std::time_t t3 = t2 - t1;
  std::cout << "AT_Parameter_Reconstruction took " << t3 << " seconds" << std::endl;

  Avec::VERBOSE_FLAG = 1;
}





// Fit only the AT beams to TOB
void AT_TOB_16par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos)
{
  // Mask containing only the modules that are to be used for the fit
  LASGlobalData<int> glob_mask = group_mask(LAS::TOB) && mask ;

  std::cout << "Performing fit AT_TOB_16par, number of modules: " << global_data_sum(glob_mask) << std::endl; 

  // Create the list of coefficients
  std::vector< LASGlobalData<double> > coeff_list;
  create_coefficients_AT_beams(coeff_list, LAS::TOB);

  // Perform the minimization
  //LASGlobalData<double> err2 = err * err ;
  std::vector<double> alpar;
  TMatrixT<double> covariance_matrix;
  LAS_Alpar_fit(dif, err2, glob_mask, coeff_list, alpar, covariance_matrix);

  // Calculate the reconstructed positions
  rec_pos = LASGlobalData<double>(0);
  for(unsigned int ip = 0; ip < alpar.size(); ++ip){
    rec_pos += coeff_list[ip] * alpar[ip];
  }

  // Calculate chi2 before and after the fit procedure
  LASGlobalData<double> chi0 = dif*dif / err2 * glob_mask;
  LASGlobalData<double> chi = (rec_pos - dif) * (rec_pos - dif) / err2 * glob_mask;
  remove_nan(chi0);
  remove_nan(chi);

  int ndof = global_data_sum(glob_mask) - coeff_list.size();

  double the_chi0 = global_data_sum(chi0, LAS::TOB) / ndof;
  double the_chi  = global_data_sum(chi,  LAS::TOB) / ndof;
  std::cout << "Chi2/ndof before: " << the_chi0 << std::endl;
  Calc_chi2(dif, err2, glob_mask, ndof);
  std::cout << "Chi2/ndof after:  " << the_chi  << std::endl;
  Calc_chi2(rec_pos - dif, err2, glob_mask, ndof);


  // Copy the alignment parameters to their locations
  for(int bn = 0; bn < 8; ++bn){
    results.beam_a[bn] = alpar[2*bn];
    results.beam_b[bn] = alpar[2*bn+1] + alpar[2*bn] * LAS::z_at_bs ;
  }
}

// // Fit only TIB and AT beams
// void AT_TIB_22par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose)
// {
//   LASGlobalData<int> glob_mask = (group_mask(LAS::TOB) || group_mask(LAS::TIB)) && mask ;

//   // Create the list of coefficients
//   std::vector< LASGlobalData<double> > coeff_list;
//   create_coefficients_AT_beams(coeff_list);
//   create_coefficients_TIB(coeff_list);

//   // Perform the minimization
//   std::vector<double> alpar;
//   TMatrixT<double> covariance_matrix;
//   LAS_Alpar_fit(dif, err2, glob_mask, coeff_list, alpar, covariance_matrix);

//   // Copy the alignment parameters to their locations
//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = alpar[2*bn];
//     results.beam_b[bn] = alpar[2*bn+1] + alpar[2*bn] * LAS::z_at_bs ;
//   }

//   results.tib_dx = alpar[16+0];
//   results.tib_dy = alpar[16+1];
//   results.tib_rx = alpar[16+2];
//   results.tib_ry = alpar[16+3];
//   results.tib_rz = alpar[16+4];
//   results.tib_tz = alpar[16+5];

//   // Calculate the reconstructed positions
//   rec_pos = LASGlobalData<double>(0);
//   for(unsigned int ip = 0; ip < alpar.size(); ++ip){
//     rec_pos += coeff_list[ip] * alpar[ip];
//   }

//   int ndof = global_data_sum(glob_mask) - coeff_list.size();

//   // Calculate chi2 before and after the fit procedure
//   if(verbose){
//     std::cout << "############## Chi2 calculation ################## " << std::endl;
//     std::cout << "Before fit: " << std::endl;
//     Calc_chi2(dif, err2, glob_mask, global_data_sum(glob_mask));
//     std::cout << "After fit: " << std::endl;
//     Calc_chi2(rec_pos - dif, err2, glob_mask, global_data_sum(glob_mask) - coeff_list.size());
//     //print_global_data(dif * 1e6, LAS::TIB);
//     //print_global_data(rec_pos * 1e6, LAS::TIB);
//     //print_global_data(TIB_spot_rec(results) * 1e6, LAS::TIB);
//   }

//   // Calculate chi2
//   LASGlobalData<double> chi0 = dif * dif / err2 * glob_mask;
//   remove_nan(chi0 );
//   double tob_chi0  = global_data_sum(chi0 * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi0  = global_data_sum(chi0 * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double AT_chi0   = global_data_sum(chi0 * (glob_mask && (group_mask(LAS::TIB) || group_mask(LAS::TOB))), LAS::AT) / ndof;
//   results.tib_chi0 = tib_chi0;
//   results.tob_chi0 = tob_chi0;
//   results.AT_chi0 = AT_chi0;

//   LASGlobalData<double> chi = (rec_pos - dif) * (rec_pos - dif) / err2 * glob_mask;
//   remove_nan(chi );
//   double tob_chi  = global_data_sum(chi * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi  = global_data_sum(chi * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double AT_chi   = global_data_sum(chi * (glob_mask && (group_mask(LAS::TIB) || group_mask(LAS::TOB))), LAS::AT) / ndof;
//   results.tib_chi = tib_chi;
//   results.tob_chi = tob_chi;
//   results.AT_chi = AT_chi;
// }


// // Fit only TIB and AT beams, no torsion
// void AT_TIB_21par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose)
// {
//   LASGlobalData<int> glob_mask = (group_mask(LAS::TOB) || group_mask(LAS::TIB)) && mask ;

//   // Create the list of coefficients
//   std::vector< LASGlobalData<double> > coeff_list;
//   create_coefficients_AT_beams(coeff_list);
//   create_coefficients_TIB_notorsion(coeff_list);

//   // Perform the minimization
//   std::vector<double> alpar;
//   TMatrixT<double> covariance_matrix;
//   LAS_Alpar_fit(dif, err2, glob_mask, coeff_list, alpar, covariance_matrix);

//   // Copy the alignment parameters to their locations
//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = alpar[2*bn];
//     results.beam_b[bn] = alpar[2*bn+1] + alpar[2*bn] * LAS::z_at_bs ;
//   }

//   results.tib_dx = alpar[16+0];
//   results.tib_dy = alpar[16+1];
//   results.tib_rx = alpar[16+2];
//   results.tib_ry = alpar[16+3];
//   results.tib_rz = alpar[16+4];

//   // Calculate the reconstructed positions
//   rec_pos = LASGlobalData<double>(0);
//   for(unsigned int ip = 0; ip < alpar.size(); ++ip){
//     rec_pos += coeff_list[ip] * alpar[ip];
//   }

//   int ndof = global_data_sum(glob_mask) - coeff_list.size();

//   // Calculate chi2 before and after the fit procedure
//   if(verbose){
//     std::cout << "############## Chi2 calculation ################## " << std::endl;
//     std::cout << "Before fit: " << std::endl;
//     Calc_chi2(dif, err2, glob_mask, global_data_sum(glob_mask));
//     std::cout << "After fit: " << std::endl;
//     Calc_chi2(rec_pos - dif, err2, glob_mask, global_data_sum(glob_mask) - coeff_list.size());
//     print_global_data(dif * 1e6, LAS::TIB);
//     print_global_data(rec_pos * 1e6, LAS::TIB);
//     print_global_data(TIB_spot_rec(results) * 1e6, LAS::TIB);
//   }

//   // Calculate chi2
//   LASGlobalData<double> chi0 = dif * dif / err2 * glob_mask;
//   remove_nan(chi0 );
//   double tob_chi0  = global_data_sum(chi0 * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi0  = global_data_sum(chi0 * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double AT_chi0   = global_data_sum(chi0 * (glob_mask && (group_mask(LAS::TIB) || group_mask(LAS::TOB))), LAS::AT) / ndof;
//   results.tib_chi0 = tib_chi0;
//   results.tob_chi0 = tob_chi0;
//   results.AT_chi0 = AT_chi0;

//   LASGlobalData<double> chi = (rec_pos - dif) * (rec_pos - dif) / err2 * glob_mask;
//   remove_nan(chi );
//   double tob_chi  = global_data_sum(chi * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi  = global_data_sum(chi * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double AT_chi   = global_data_sum(chi * (glob_mask && (group_mask(LAS::TIB) || group_mask(LAS::TOB))), LAS::AT) / ndof;
//   results.tib_chi = tib_chi;
//   results.tob_chi = tob_chi;
//   results.AT_chi = AT_chi;


// }


// // Fit only TIB, TEC- and AT beams
// void AT_TIB_TECM_25par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose)
// {
//   LASGlobalData<int> glob_mask = (group_mask(LAS::TOB) || group_mask(LAS::TIB) || group_mask(LAS::TEC_MINUS_AT)) && mask ;

//   // Create the list of coefficients
//   std::vector< LASGlobalData<double> > coeff_list;
//   create_coefficients_AT_beams(coeff_list);
//   create_coefficients_TIB(coeff_list);

//   // AT TEC- translations
//   coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TEC_MINUS_AT));
//   coeff_list.push_back( cos(theta()) / r0() * group_mask(LAS::TEC_MINUS_AT));

//   // TEC- rotation
//   coeff_list.push_back(group_mask(LAS::TEC_MINUS_AT));

//   coeff_list.push_back(cos(theta()) / r0() * zpos() * group_mask(LAS::TEC_MINUS_AT));
//   coeff_list.push_back(sin(theta()) / r0() * zpos() * group_mask(LAS::TEC_MINUS_AT));

//   // TEC- torsion
//   coeff_list.push_back(zpos() * group_mask(LAS::TEC_MINUS_AT));


//   // Perform the minimization
//   std::vector<double> alpar;
//   TMatrixT<double> covariance_matrix;
//   LAS_Alpar_fit(dif, err2, glob_mask, coeff_list, alpar, covariance_matrix);

//   // Copy the alignment parameters to their locations
//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = alpar[2*bn];
//     results.beam_b[bn] = alpar[2*bn+1] + alpar[2*bn] * LAS::z_at_bs ;
//   }

//   results.tib_dx  = alpar[16+0];
//   results.tib_dy  = alpar[16+1];
//   results.tib_rx  = alpar[16+2];
//   results.tib_ry  = alpar[16+3];
//   results.tib_rz  = alpar[16+4];
//   results.tib_tz  = alpar[16+5];
//   results.tecm_dx = alpar[16+6];
//   results.tecm_dy = alpar[16+7];
//   results.tecm_rz = alpar[16+8];

//   // Calculate the reconstructed positions
//   rec_pos = LASGlobalData<double>(0);
//   for(unsigned int ip = 0; ip < alpar.size(); ++ip){
//     rec_pos += coeff_list[ip] * alpar[ip];
//   }

//   int ndof = global_data_sum(glob_mask) - coeff_list.size();

//   // Calculate chi2 before and after the fit procedure
//   if(verbose){
//     std::cout << "############## Chi2 calculation ################## " << std::endl;
//     std::cout << "Before fit: " << std::endl;
//     Calc_chi2(dif, err2, glob_mask, global_data_sum(glob_mask));
//     std::cout << "After fit: " << std::endl;
//     Calc_chi2(rec_pos - dif, err2, glob_mask, global_data_sum(glob_mask) - coeff_list.size());
//     //print_global_data(dif * 1e6, LAS::TIB);
//     //print_global_data(rec_pos * 1e6, LAS::TIB);
//     //print_global_data(TIB_spot_rec(results) * 1e6, LAS::TIB);
//   }

//   // Calculate chi2
//   LASGlobalData<double> chi0 = dif * dif / err2 * glob_mask;
//   remove_nan(chi0 );
//   double tob_chi0  = global_data_sum(chi0 * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi0  = global_data_sum(chi0 * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double tecm_chi0  = global_data_sum(chi0 * glob_mask, LAS::TEC_MINUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_MINUS_AT));
//   double AT_chi0   = global_data_sum(chi0 * (glob_mask && (group_mask(LAS::TIB) || group_mask(LAS::TOB) || group_mask(LAS::TEC_MINUS_AT))), LAS::AT) / ndof;
//   results.tib_chi0 = tib_chi0;
//   results.tob_chi0 = tob_chi0;
//   results.tecm_chi0 = tecm_chi0;
//   results.AT_chi0 = AT_chi0;

//   LASGlobalData<double> chi = (rec_pos - dif) * (rec_pos - dif) / err2 * glob_mask;
//   remove_nan(chi );
//   double tob_chi  = global_data_sum(chi * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi  = global_data_sum(chi * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double tecm_chi  = global_data_sum(chi * glob_mask, LAS::TEC_MINUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_MINUS_AT));
//   double AT_chi   = global_data_sum(chi * (glob_mask && (group_mask(LAS::TIB) || group_mask(LAS::TOB) || group_mask(LAS::TEC_MINUS_AT))), LAS::AT) / ndof;
//   results.tib_chi = tib_chi;
//   results.tob_chi = tob_chi;
//   results.tecm_chi = tecm_chi;
//   results.AT_chi = AT_chi;
// }

// // Fit only TIB, TEC+ and AT beams
// void AT_TIB_TECP(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose)
// {
//   LASGlobalData<int> glob_mask = (group_mask(LAS::TOB) || group_mask(LAS::TIB) || group_mask(LAS::TEC_PLUS_AT)) && mask ;

//   // Create the list of coefficients
//   std::vector< LASGlobalData<double> > coeff_list;
//   create_coefficients_AT_beams(coeff_list);
//   create_coefficients_TIB(coeff_list);

//   // AT TEC+ translations
//   coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TEC_PLUS_AT));
//   coeff_list.push_back( cos(theta()) / r0() * group_mask(LAS::TEC_PLUS_AT));

//   // TEC+ rotation
//   coeff_list.push_back(group_mask(LAS::TEC_PLUS_AT));

//   coeff_list.push_back(cos(theta()) / r0() * zpos() * group_mask(LAS::TEC_PLUS_AT));
//   coeff_list.push_back(sin(theta()) / r0() * zpos() * group_mask(LAS::TEC_PLUS_AT));

//   // TEC+ torsion
//   coeff_list.push_back(zpos() * group_mask(LAS::TEC_PLUS_AT));


//   // Perform the minimization
//   std::vector<double> alpar;
//   TMatrixT<double> covariance_matrix;
//   LAS_Alpar_fit(dif, err2, glob_mask, coeff_list, alpar, covariance_matrix);

//   // Copy the alignment parameters to their locations
//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = alpar[2*bn];
//     results.beam_b[bn] = alpar[2*bn+1] + alpar[2*bn] * LAS::z_at_bs ;
//   }

//   results.tib_dx  = alpar[16+0];
//   results.tib_dy  = alpar[16+1];
//   results.tib_rx  = alpar[16+2];
//   results.tib_ry  = alpar[16+3];
//   results.tib_rz  = alpar[16+4];
//   results.tib_tz  = alpar[16+5];
//   results.tecp_dx = alpar[16+6];
//   results.tecp_dy = alpar[16+7];
//   results.tecp_rz = alpar[16+8];

//   // Calculate the reconstructed positions
//   rec_pos = LASGlobalData<double>(0);
//   for(unsigned int ip = 0; ip < alpar.size(); ++ip){
//     rec_pos += coeff_list[ip] * alpar[ip];
//   }

//   int ndof = global_data_sum(glob_mask) - coeff_list.size();

//   // Calculate chi2 before and after the fit procedure
//   if(verbose){
//     std::cout << "############## Chi2 calculation ################## " << std::endl;
//     std::cout << "Before fit: " << std::endl;
//     Calc_chi2(dif, err2, glob_mask, global_data_sum(glob_mask));
//     std::cout << "After fit: " << std::endl;
//     Calc_chi2(rec_pos - dif, err2, glob_mask, global_data_sum(glob_mask) - coeff_list.size());
//     //print_global_data(dif * 1e6, LAS::TIB);
//     //print_global_data(rec_pos * 1e6, LAS::TIB);
//     //print_global_data(TIB_spot_rec(results) * 1e6, LAS::TIB);
//   }

//   // Calculate chi2
//   LASGlobalData<double> chi0 = dif * dif / err2 * glob_mask;
//   remove_nan(chi0 );
//   double tob_chi0  = global_data_sum(chi0 * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi0  = global_data_sum(chi0 * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double tecp_chi0  = global_data_sum(chi0 * glob_mask, LAS::TEC_PLUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_PLUS_AT));
//   double AT_chi0   = global_data_sum(chi0 * (glob_mask && (group_mask(LAS::TIB) || group_mask(LAS::TOB) || group_mask(LAS::TEC_PLUS_AT))), LAS::AT) / ndof;
//   results.tib_chi0 = tib_chi0;
//   results.tob_chi0 = tob_chi0;
//   results.tecp_chi0 = tecp_chi0;
//   results.AT_chi0 = AT_chi0;

//   LASGlobalData<double> chi = (rec_pos - dif) * (rec_pos - dif) / err2 * glob_mask;
//   remove_nan(chi );
//   double tob_chi  = global_data_sum(chi * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi  = global_data_sum(chi * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double tecp_chi  = global_data_sum(chi * glob_mask, LAS::TEC_PLUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_PLUS_AT));
//   double AT_chi   = global_data_sum(chi * (glob_mask && (group_mask(LAS::TIB) || group_mask(LAS::TOB) || group_mask(LAS::TEC_PLUS_AT))), LAS::AT) / ndof;
//   results.tib_chi = tib_chi;
//   results.tob_chi = tob_chi;
//   results.tecp_chi = tecp_chi;
//   results.AT_chi = AT_chi;
// }



// // Fit only TIB, TEC+, TEC- and AT beams
// void AT_full(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose)
// {
//   LASGlobalData<int> glob_mask = group_mask(LAS::AT) && mask ;

//   // Create the list of coefficients
//   std::vector< LASGlobalData<double> > coeff_list;
//   create_coefficients_AT_beams(coeff_list);
//   create_coefficients_TIB(coeff_list);

//   // AT TEC+ translations
//   coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TEC_PLUS_AT));
//   coeff_list.push_back( cos(theta()) / r0() * group_mask(LAS::TEC_PLUS_AT));

//   // TEC+ rotation
//   coeff_list.push_back(group_mask(LAS::TEC_PLUS_AT));

//   //coeff_list.push_back(cos(theta()) / r0() * zpos() * group_mask(LAS::TEC_PLUS_AT));
//   //coeff_list.push_back(sin(theta()) / r0() * zpos() * group_mask(LAS::TEC_PLUS_AT));

//   // TEC+ torsion
//   //coeff_list.push_back(zpos() * group_mask(LAS::TEC_PLUS_AT));


//   // AT TEC- translations
//   coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TEC_MINUS_AT));
//   coeff_list.push_back( cos(theta()) / r0() * group_mask(LAS::TEC_MINUS_AT));

//   // TEC- rotation
//   coeff_list.push_back(group_mask(LAS::TEC_MINUS_AT));

//   //coeff_list.push_back(cos(theta()) / r0() * zpos() * group_mask(LAS::TEC_MINUS_AT));
//   //coeff_list.push_back(sin(theta()) / r0() * zpos() * group_mask(LAS::TEC_MINUS_AT));

//   // TEC- torsion
//   //coeff_list.push_back(zpos() * group_mask(LAS::TEC_MINUS_AT));


//   // Perform the minimization
//   std::vector<double> alpar;
//   TMatrixT<double> covariance_matrix;
//   LAS_Alpar_fit(dif, err2, glob_mask, coeff_list, alpar, covariance_matrix);

//   // Copy the alignment parameters to their locations
//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = alpar[2*bn];
//     results.beam_b[bn] = alpar[2*bn+1] + alpar[2*bn] * LAS::z_at_bs ;
//   }

//   results.tib_dx  = alpar[16+0];
//   results.tib_dy  = alpar[16+1];
//   results.tib_rx  = alpar[16+2];
//   results.tib_ry  = alpar[16+3];
//   results.tib_rz  = alpar[16+4];
//   results.tib_tz  = alpar[16+5];

//   results.tecp_dx = alpar[16+6];
//   results.tecp_dy = alpar[16+7];
//   results.tecp_rz = alpar[16+8];

//   results.tecm_dx = alpar[16+9];
//   results.tecm_dy = alpar[16+10];
//   results.tecm_rz = alpar[16+11];

//   //results.tecm_dx = alpar[16+12];
//   //results.tecm_dy = alpar[16+13];
//   //results.tecm_rz = alpar[16+14];

//   // Calculate the reconstructed positions
//   rec_pos = LASGlobalData<double>(0);
//   for(unsigned int ip = 0; ip < alpar.size(); ++ip){
//     rec_pos += coeff_list[ip] * alpar[ip];
//   }

//   int ndof = global_data_sum(glob_mask) - coeff_list.size();

//   // Calculate chi2 before and after the fit procedure
//   if(verbose){
//     std::cout << "############## Chi2 calculation ################## " << std::endl;
//     std::cout << "Before fit: " << std::endl;
//     Calc_chi2(dif, err2, glob_mask, global_data_sum(glob_mask));
//     std::cout << "After fit: " << std::endl;
//     Calc_chi2(rec_pos - dif, err2, glob_mask, global_data_sum(glob_mask) - coeff_list.size());
//     //print_global_data(dif * 1e6, LAS::TIB);
//     //print_global_data(rec_pos * 1e6, LAS::TIB);
//     //print_global_data(TIB_spot_rec(results) * 1e6, LAS::TIB);
//   }

//   // Calculate chi2
//   LASGlobalData<double> chi0 = dif * dif / err2 * glob_mask;
//   remove_nan(chi0 );
//   double tob_chi0  = global_data_sum(chi0 * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi0  = global_data_sum(chi0 * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double tecp_chi0  = global_data_sum(chi0 * glob_mask, LAS::TEC_PLUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_PLUS_AT));
//   double tecm_chi0  = global_data_sum(chi0 * glob_mask, LAS::TEC_MINUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_MINUS_AT));
//   double AT_chi0   = global_data_sum(chi0 * (glob_mask && group_mask(LAS::AT)), LAS::AT) / ndof;
//   results.tib_chi0 = tib_chi0;
//   results.tob_chi0 = tob_chi0;
//   results.tecp_chi0 = tecp_chi0;
//   results.tecm_chi0 = tecm_chi0;
//   results.AT_chi0 = AT_chi0;

//   LASGlobalData<double> chi = (rec_pos - dif) * (rec_pos - dif) / err2 * glob_mask;
//   remove_nan(chi );
//   double tob_chi  = global_data_sum(chi * glob_mask, LAS::TOB) / global_data_sum(glob_mask && group_mask(LAS::TOB));
//   double tib_chi  = global_data_sum(chi * glob_mask, LAS::TIB) / global_data_sum(glob_mask && group_mask(LAS::TIB));
//   double tecp_chi  = global_data_sum(chi * glob_mask, LAS::TEC_PLUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_PLUS_AT));
//   double tecm_chi  = global_data_sum(chi * glob_mask, LAS::TEC_MINUS_AT) / global_data_sum(glob_mask && group_mask(LAS::TEC_MINUS_AT));
//   double AT_chi   = global_data_sum(chi * (glob_mask && group_mask(LAS::AT)), LAS::AT) / ndof;
//   results.tib_chi = tib_chi;
//   results.tob_chi = tob_chi;
//   results.tecp_chi = tecp_chi;
//   results.tecm_chi = tecm_chi;
//   results.AT_chi = AT_chi;
// }



// // reconstruction of alignment parameters in a generalized way
// // Under construction
// // Goal is to supply vector with coefficients, weights and spot positions to construct a matrix and a vector
// // The matrix is inverted and multiplied with the vector, yielding a vector with the alignment parameters.
// void AlPar_fit_old(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, bool swap_flag, bool verbose_flag)
// {
//   std::cout << "Calling AlPar_fit" << std::endl;

//   // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

//   double bs[22];
//   double db[22][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TIB, TOB); k=0..7 (beam number); l=0..5 (module number)

//   int npb = 16;// number of beam parameters
 
//   LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
//   LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
//   LASGlobalData<double> par2 = dif * cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par3 = dif * sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par4 = dif / err2;
//   LASGlobalData<double> par5 = dif * zpos()/ err2;


//   LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
//   LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
//   LASGlobalData<double> dev2 = cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev3 = sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> dev5 = zpos()/ err2;


//   // Create the list of coefficients
//   std::vector< LASGlobalData<double> > coeff_list;

//   // Coefficients for beam parameters
//   for(int beam_nr = 0; beam_nr < 8; ++beam_nr){
//     LASGlobalData<double> at_beam_mask = group_mask(LAS::TOB) || group_mask(LAS::TIB); // TOB and TIB are 1 , the rest is 0
//     // Mask everything but one beam
//     LASGlobalDataLoop lp(LASGlobalDataLoop::AT);
//     do{if(lp.get_beam() != beam_nr)lp.GetEntry(at_beam_mask)=0; }while(lp.next());

//     coeff_list.push_back(at_beam_mask * zpos()); // Coefficient for the slope
//     coeff_list.push_back(at_beam_mask); // Coefficient for the offset
//   }

//   // TIB translations
//   coeff_list.push_back(-sin(theta()) / r0() * group_mask(LAS::TIB));
//   coeff_list.push_back(cos(theta()) / r0() * group_mask(LAS::TIB));

//   // TIB rotations
//   coeff_list.push_back(cos(theta()) / r0() * zpos() * group_mask(LAS::TIB));
//   coeff_list.push_back(sin(theta()) / r0() * zpos() * group_mask(LAS::TIB));
//   coeff_list.push_back(group_mask(LAS::TIB));

//   // TIB torsion
//   coeff_list.push_back(zpos() * group_mask(LAS::TIB));

//   //std::cout << "coeff_list.size(): " << coeff_list.size() << std::endl;


//   //////////////////////////////////////////////////////////////////////////////////////////////////////////
//   // From here on the function should be generic (independent of the fit type)
//   // Construct Vector
//   Avec beta;
//   for(unsigned int idx = 0; idx < coeff_list.size(); idx++){
//     LASGlobalData<double> vector_entry = dif * coeff_list[idx] * mask / err2;
//     remove_nan(vector_entry);
//     beta.push_back(global_data_sum(vector_entry, LAS::ALL));
//   }

//   // Construct the matrix
//   TMatrixT<double> h2(coeff_list.size(), coeff_list.size());
  
//   for(unsigned int ix = 0; ix < coeff_list.size(); ++ix){
//     for(unsigned int iy = 0; iy < coeff_list.size(); ++iy){
//       LASGlobalData<double> matrix_entry = coeff_list[ix] * coeff_list[iy] * mask / err2;
//       h2(ix,iy) = global_data_sum(matrix_entry, LAS::AT);
//     }
//   }

//   //////////////////////////////////////////////////////////////////////////////////////////////////////////

//   for(int ip = 0; ip < 22;  ++ip){ 
//     bs[ip] = 0.;
//   }
//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 6;  ++dn){
//       for(int ip = 0; ip < 22;  ++ip){ 
// 	db[ip][0][bn][dn] = 0.;
// 	db[ip][1][bn][dn] = 0.;
//       }
//       if(err2.GetEntry(2,-1,bn,dn) < 1000. && mask.GetEntry(2, -1, bn, dn)){
	    
// 	bs[2*bn] += par5.GetEntry(2,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(2,-1,bn,dn);
	
// 	bs[npb+0] += par0.GetEntry(2,-1,bn,dn);
// 	bs[npb+1] += par1.GetEntry(2,-1,bn,dn);
// 	bs[npb+2] += par2.GetEntry(2,-1,bn,dn);
// 	bs[npb+3] += par3.GetEntry(2,-1,bn,dn);
// 	bs[npb+4] += par4.GetEntry(2,-1,bn,dn);
// 	bs[npb+5] += par5.GetEntry(2,-1,bn,dn);
        
// 	db[2*bn][0][bn][dn] = dev5.GetEntry(2,-1,bn,dn);
// 	db[2*bn+1][0][bn][dn] = dev4.GetEntry(2,-1,bn,dn);
	
// 	db[npb+0][0][bn][dn] = dev0.GetEntry(2,-1,bn,dn);
// 	db[npb+1][0][bn][dn] = dev1.GetEntry(2,-1,bn,dn);
// 	db[npb+2][0][bn][dn] = dev2.GetEntry(2,-1,bn,dn);
// 	db[npb+3][0][bn][dn] = dev3.GetEntry(2,-1,bn,dn);
// 	db[npb+4][0][bn][dn] = dev4.GetEntry(2,-1,bn,dn);
// 	db[npb+5][0][bn][dn] = dev5.GetEntry(2,-1,bn,dn);
//       }
//       if(err2.GetEntry(3,-1,bn,dn) < 1000. && mask.GetEntry(3, -1, bn, dn)){
	    
// 	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
// 	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
// 	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);
//       }
//     }
//   }
   
//   TMatrixT<double> h(22,22);
  
//   for(int ix = 0; ix < 22; ++ix){
//     for(int iy = 0; iy < 22; ++iy){
//       h(ix,iy) = 0.;
//     }
//   }
  
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


//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 6;  ++dn){
      
//       if(err2.GetEntry(2,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(2,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(2,-1,bn,dn);
	
// 	h(npb+0,npb+0)   += h00.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+1)   += h01.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+2)   += h02.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+3)   += h03.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+4)   += h04.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+5)   += h05.GetEntry(2,-1,bn,dn);
// 	h(npb+0,2*bn)    += h05.GetEntry(2,-1,bn,dn);
// 	h(npb+0,2*bn+1)  += h04.GetEntry(2,-1,bn,dn);
	
// 	h(npb+1,npb+0)   += h10.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+1)   += h11.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+2)   += h12.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+3)   += h13.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+4)   += h14.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+5)   += h15.GetEntry(2,-1,bn,dn);
// 	h(npb+1,2*bn)    += h15.GetEntry(2,-1,bn,dn);
// 	h(npb+1,2*bn+1)  += h14.GetEntry(2,-1,bn,dn);
	
// 	h(npb+2,npb+0)   += h20.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+1)   += h21.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+2)   += h22.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+3)   += h23.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+4)   += h24.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+5)   += h25.GetEntry(2,-1,bn,dn);
// 	h(npb+2,2*bn)    += h25.GetEntry(2,-1,bn,dn);
// 	h(npb+2,2*bn+1)  += h24.GetEntry(2,-1,bn,dn);
	
// 	h(npb+3,npb+0)   += h30.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+1)   += h31.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+2)   += h32.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+3)   += h33.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+4)   += h34.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+5)   += h35.GetEntry(2,-1,bn,dn);
// 	h(npb+3,2*bn)    += h35.GetEntry(2,-1,bn,dn);
// 	h(npb+3,2*bn+1)  += h34.GetEntry(2,-1,bn,dn);
	
// 	h(npb+4,npb+0)   += h40.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+1)   += h41.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+2)   += h42.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+3)   += h43.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+4)   += h44.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+5)   += h45.GetEntry(2,-1,bn,dn);
// 	h(npb+4,2*bn)    += h45.GetEntry(2,-1,bn,dn);
// 	h(npb+4,2*bn+1)  += h44.GetEntry(2,-1,bn,dn);
	
// 	h(npb+5,npb+0)   += h50.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+1)   += h51.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+2)   += h52.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+3)   += h53.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+4)   += h54.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+5)   += h55.GetEntry(2,-1,bn,dn);
// 	h(npb+5,2*bn)    += h55.GetEntry(2,-1,bn,dn);
// 	h(npb+5,2*bn+1)  += h54.GetEntry(2,-1,bn,dn);
	
	
// 	h(2*bn,npb+0)   += h50.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+1)   += h51.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+2)   += h52.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+3)   += h53.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+4)   += h54.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+5)   += h55.GetEntry(2,-1,bn,dn);
	
// 	h(2*bn+1,npb+0) += h40.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+1) += h41.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+2) += h42.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+3) += h43.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+4) += h44.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+5) += h45.GetEntry(2,-1,bn,dn);
//       }
	
//       if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
//       }
//     }
//   }  

//   if(verbose_flag){
//     for(int ix = 0; ix < 22; ++ix){
//       for(int iy = 0; iy < 22; ++iy){
// 	if(1 - fabs(h(ix, iy) / h2(ix, iy)) > 1e-3) std::cout << "--------------->";
// 	std::cout << "ix: " << ix << "  iy: " << iy << "  h(ix,iy): " << h(ix,iy) << "  h2(ix,iy): " << h2(ix,iy) << std::endl;
//       }
//     }
//   }

//   if(swap_flag){
//     std::cout << "Using new matrix" << std::endl;
//     h=h2;
//   }
//   h.Invert();

//   double cp[22][2][8][6];
//   double up[22][22];
//   double ps[22];

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ib = 0; ib < 2; ++ib){
//       for(int bn = 0; bn < 8;  ++bn){
// 	for(int dn = 0; dn < 6;  ++dn){
// 	  cp[ip][ib][bn][dn] = 0.;
// 	  for(int ix = 0; ix < 22; ++ix){
// 	    cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
// 	    }
// 	}
//       }
//     }
//   }

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ix = 0; ix < 22; ++ix){
//       up[ip][ix] = 0.;
//       for(int ib = 0; ib < 2; ++ib){
// 	for(int bn = 0; bn < 8;  ++bn){
// 	  for(int dn = 0; dn < 6;  ++dn){
// 	    if( ib == 0){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(2,-1,bn,dn);
// 	    }
// 	    if( ib == 1 ){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
 
//   for(int bn = 0; bn < 8; ++bn){
//     results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
//     results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
//   }

//   results.er_tib_dx = sqrt(up[npb+0][npb+0]);
//   results.er_tib_dy = sqrt(up[npb+1][npb+1]);
//   results.er_tib_rx = sqrt(up[npb+2][npb+2]);
//   results.er_tib_ry = sqrt(up[npb+3][npb+3]);
//   results.er_tib_rz = sqrt(up[npb+4][npb+4]);
//   results.er_tib_tz = sqrt(up[npb+5][npb+5]);
//   results.tib_dxdy  = up[npb+0][npb+1]/results.er_tib_dx/results.er_tib_dy;
//   results.tib_dxrx  = up[npb+0][npb+2]/results.er_tib_dx/results.er_tib_rx;
//   results.tib_dxry  = up[npb+0][npb+3]/results.er_tib_dx/results.er_tib_ry;
//   results.tib_dxrz  = up[npb+0][npb+4]/results.er_tib_dx/results.er_tib_rz;
//   results.tib_dxtz  = up[npb+0][npb+5]/results.er_tib_dx/results.er_tib_tz;
//   results.tib_dyrx  = up[npb+1][npb+2]/results.er_tib_dy/results.er_tib_rx;
//   results.tib_dyry  = up[npb+1][npb+3]/results.er_tib_dy/results.er_tib_ry;
//   results.tib_dyrz  = up[npb+1][npb+4]/results.er_tib_dy/results.er_tib_rz;
//   results.tib_dytz  = up[npb+1][npb+5]/results.er_tib_dy/results.er_tib_tz;
//   results.tib_rxry  = up[npb+2][npb+3]/results.er_tib_rx/results.er_tib_ry;
//   results.tib_rxrz  = up[npb+2][npb+4]/results.er_tib_rx/results.er_tib_rz;
//   results.tib_rxtz  = up[npb+2][npb+5]/results.er_tib_rx/results.er_tib_tz;
//   results.tib_ryrz  = up[npb+3][npb+4]/results.er_tib_ry/results.er_tib_rz;
//   results.tib_rytz  = up[npb+3][npb+5]/results.er_tib_ry/results.er_tib_tz;
//   results.tib_rztz  = up[npb+4][npb+5]/results.er_tib_rz/results.er_tib_tz;

//   LASGlobalData<double> rec_pos(0);
//   for(int ip = 0; ip < 22; ++ip){
//     if(verbose_flag){
//       std::cout << "idx: " << ip << "  bs[idx]: " << bs[ip] << "  beta[idx]: " << beta[ip] << std::endl;
//     }
//     ps[ip] = 0.;
//     for(int ix = 0; ix < 22; ++ix){
//       ps[ip] += h(ip,ix)*bs[ix];
//     }
//     rec_pos += coeff_list[ip] * ps[ip];// Calculate the reconstructed positions
//   }

//   LASGlobalData<double> chi0 = fabs(dif) / err2;
//   LASGlobalData<double> chi = fabs(rec_pos) / err2;
//   remove_nan(chi0);
//   remove_nan(chi);


//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = ps[2*bn];
//     results.beam_b[bn] = ps[2*bn+1];
//   }

//   results.tib_dx = ps[npb+0];
//   results.tib_dy = ps[npb+1];
//   results.tib_rx = ps[npb+2];
//   results.tib_ry = ps[npb+3];
//   results.tib_rz = ps[npb+4];
//   results.tib_tz = ps[npb+5];

//   results.tib_chi0 = sqrt( global_data_sum(pow(chi0, 2), LAS::TIB) / global_data_sum(mask, LAS::TIB));
//   results.tib_chi  = sqrt( global_data_sum(pow(chi,  2), LAS::TIB) / global_data_sum(mask, LAS::TIB));
//   std::cout << "Number of TIB modules: " << global_data_sum(mask, LAS::TIB) << std::endl;

//   std::cout << "----------------TIB_spot_rec: -----------------" << std::endl;
//   print_global_data(TIB_spot_rec(results), LAS::TIB);
//   std::cout << "----------------spot_rec: ---------------------" << std::endl;
//   print_global_data(rec_pos, LAS::TIB);
//   std::cout << "----------------dif: --------------------------" << std::endl;
//   print_global_data(dif, LAS::TIB);
//   std::cout << "----------------err2: -------------------------" << std::endl;
//   print_global_data(err2, LAS::TIB);

// }

LASGlobalData<double> AT_spot_rec(const AtPar& par) //calculation of beam spot positions at TIB&TOB modules
{
  LASGlobalData<double> result;
  LASGlobalDataLoop tib_loop(LASGlobalDataLoop::TIB);

  do{
    int bn = tib_loop.get_beam();
    int dn = tib_loop.get_zpos();
    tib_loop.GetEntry(result) = par.beam_a[bn]*zpos().GetEntry(2,-1,bn,dn) + par.beam_b[bn] - par.tib_dx*sin(theta().GetEntry(2,-1,bn,dn))/r0().GetEntry(2,-1,bn,dn) + par.tib_dy*cos(theta().GetEntry(2,-1,bn,dn))/r0().GetEntry(2,-1,bn,dn) + par.tib_rx*cos(theta().GetEntry(2,-1,bn,dn))*zpos().GetEntry(2,-1,bn,dn)/r0().GetEntry(2,-1,bn,dn) + par.tib_ry*sin(theta().GetEntry(2,-1,bn,dn))*zpos().GetEntry(2,-1,bn,dn)/r0().GetEntry(2,-1,bn,dn) + par.tib_rz+ par.tib_tz*zpos().GetEntry(2,-1,bn,dn); 
  }while(tib_loop.next());
  
  LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);

  do{
    int bn = tob_loop.get_beam();
    int dn = tob_loop.get_zpos();
    tob_loop.GetEntry(result) = par.beam_a[bn]*zpos().GetEntry(3,-1,bn,dn) + par.beam_b[bn]; 
  }while(tob_loop.next());

  LASGlobalDataLoop tecp_loop(LASGlobalDataLoop::TEC_PLUS_AT);

  do{
    int bn = tecp_loop.get_beam();
    int dn = tecp_loop.get_zpos();
    tecp_loop.GetEntry(result) = par.beam_a[bn]*zpos().GetEntry(0,-1,bn,dn) + par.beam_b[bn] - par.tecp_dx*sin(theta().GetEntry(0,-1,bn,dn))/r0().GetEntry(0,-1,bn,dn) + par.tecp_dy*cos(theta().GetEntry(0,-1,bn,dn))/r0().GetEntry(0,-1,bn,dn) + par.tecp_rz; 
    //cout<<"  bn= "<<bn<<"  dn= "<<dn<<"  recpos: "<<tecp_loop.GetEntry(result)<<endl;  
  }while(tecp_loop.next());

  LASGlobalDataLoop tecm_loop(LASGlobalDataLoop::TEC_MINUS_AT);

  do{
    int bn = tecm_loop.get_beam();
    int dn = tecm_loop.get_zpos();
    tecm_loop.GetEntry(result) = par.beam_a[bn]*zpos().GetEntry(1,-1,bn,dn) + par.beam_b[bn] - par.tecm_dx*sin(theta().GetEntry(1,-1,bn,dn))/r0().GetEntry(1,-1,bn,dn) + par.tecm_dy*cos(theta().GetEntry(1,-1,bn,dn))/r0().GetEntry(1,-1,bn,dn) + par.tecm_rz; 
    //cout<<"  bn= "<<bn<<"  dn= "<<dn<<"  recpos: "<<tecm_loop.GetEntry(result)<<endl;  
  }while(tecm_loop.next());

  return result;
}

LASGlobalData<double> TIB_spot_rec(const AtPar& par) //calculation of beam spot positions at TIB&TOB modules
{
  LASGlobalData<double> result;

  const LASGlobalData<double> coeff_dx = -sin(theta()) / r0();
  const LASGlobalData<double> coeff_dy =  cos(theta()) / r0();
  const LASGlobalData<double> coeff_rx =  cos(theta()) / r0() * zpos();
  const LASGlobalData<double> coeff_ry =  sin(theta()) / r0() * zpos();

  LASGlobalDataLoop tib_loop(LASGlobalDataLoop::TIB);
  do{
    int bn = tib_loop.get_beam();
    tib_loop.GetEntry(result) = 
      par.beam_a[bn] * tib_loop.GetEntry(zpos()) 
      + par.beam_b[bn] 
      + par.tib_dx * tib_loop.GetEntry( coeff_dx ) 
      + par.tib_dy * tib_loop.GetEntry( coeff_dy ) 
      + par.tib_rx * tib_loop.GetEntry( coeff_rx ) 
      + par.tib_ry * tib_loop.GetEntry( coeff_ry )
      + par.tib_rz 
      + par.tib_tz * tib_loop.GetEntry( zpos() ); 
  }while(tib_loop.next());
  
  LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);
  do{
    int bn = tob_loop.get_beam();
    tob_loop.GetEntry(result) = 
      par.beam_a[bn] * tob_loop.GetEntry(zpos())
      + par.beam_b[bn]; 
  }while(tob_loop.next());
  return result;
}

LASGlobalData<double> TECP_spot_rec(const AtPar& par) //calculation of beam spot positions at TEC+&TOB modules
{
  LASGlobalData<double> result;

  const LASGlobalData<double> coeff_dx = -sin(theta()) / r0();
  const LASGlobalData<double> coeff_dy =  cos(theta()) / r0();

  LASGlobalDataLoop tecp_loop(LASGlobalDataLoop::TEC_PLUS_AT);
  do{
    int bn = tecp_loop.get_beam();
    tecp_loop.GetEntry(result) = 
      par.beam_a[bn] * tecp_loop.GetEntry( zpos() )
      + par.beam_b[bn] 
      + par.tecp_dx * tecp_loop.GetEntry( coeff_dx )
      + par.tecp_dy * tecp_loop.GetEntry( coeff_dy )
      + par.tecp_rz; 
  }while(tecp_loop.next());
   
  LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);
  do{
    int bn = tob_loop.get_beam();
    tob_loop.GetEntry(result) = 
      par.beam_a[bn] * tob_loop.GetEntry(zpos())
      + par.beam_b[bn]; 
  }while(tob_loop.next());

  return result;
}

LASGlobalData<double> TECM_spot_rec(const AtPar& par) //calculation of beam spot positions at TEC-&TOB modules
{
  LASGlobalData<double> result;

  const LASGlobalData<double> coeff_dx = -sin(theta()) / r0();
  const LASGlobalData<double> coeff_dy =  cos(theta()) / r0();

  LASGlobalDataLoop tecm_loop(LASGlobalDataLoop::TEC_MINUS_AT);
  do{
    int bn = tecm_loop.get_beam();
    tecm_loop.GetEntry(result) = 
      par.beam_a[bn] * tecm_loop.GetEntry( zpos() )
      + par.beam_b[bn] 
      + par.tecm_dx * tecm_loop.GetEntry( coeff_dx ) 
      + par.tecm_dy * tecm_loop.GetEntry( coeff_dy )
      + par.tecm_rz; 
  }while(tecm_loop.next());
   
  LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);
  do{
    int bn = tob_loop.get_beam();
    tob_loop.GetEntry(result) = 
      par.beam_a[bn] * tob_loop.GetEntry(zpos())
      + par.beam_b[bn]; 
  }while(tob_loop.next());

  return result;
}




LASGlobalData<double> AT_spot_rec_old(const AtPar& par) //calculation of beam spot positions at TIB&TOB modules
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

// void AT_chi2_old(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,  AtPar& par) //calculation of chi2/ndf at TIB&TOB
// {

//   LASGlobalData<double> recpos = AT_spot_rec(par);
//   LASGlobalData<double> resid = recpos - dif;
//   LASGlobalData<double> res2 = pow(resid/err2,2);
//   LASGlobalData<double> dif2 = pow(dif/err2,2);

//   int nchi_tib,nchi_tob,nchi,nchib[8];
  
//   nchi = 0;
//   nchi_tib = 0;
//   par.tib_chi = 0.;
//   par.tib_chi0 = 0.;
//   par.AT_chi = 0.;
//   par.AT_chi0 = 0.;
//   for(int bn=0;bn<8;bn++){
//     nchib[bn] = 0;
//     par.b_chi[bn] = 0.;
//     par.b_chi0[bn] = 0.;
//   }

//   //std::cout << "tib res2 : \n";
//   LASGlobalDataLoop tib_loop(LASGlobalDataLoop::TIB);
//   do{
//     int bn = tib_loop.get_beam();
//     //nchib[bn] = 0;
//     //par.b_chi[bn] = 0.;
//     //par.b_chi0[bn] = 0.;
//     //int dn = tib_loop.get_zpos();
//     //std::cout << tib_loop.GetEntry(res2) << "  ";
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tib = nchi_tib + 1;
//     par.AT_chi += tib_loop.GetEntry(res2);
//     par.AT_chi0 += tib_loop.GetEntry(dif2);
//     par.tib_chi += tib_loop.GetEntry(res2);
//     par.tib_chi0 += tib_loop.GetEntry(dif2);
//     par.b_chi[bn] += tib_loop.GetEntry(res2);
//     par.b_chi0[bn] += tib_loop.GetEntry(dif2);
//     }while(tib_loop.next());
//   //std::cout << "\ntib_chi: " << par.tib_chi << std::endl;
//   par.tib_chi = sqrt(par.tib_chi/nchi_tib);
//   par.tib_chi0 = sqrt(par.tib_chi0/nchi_tib);
//   //cout<<"irun: "<<irun<<"  nchi_tib: "<<nchi_tib<<"  chi_tib: "<<chi_tib[irun]<<endl;   
 
//   nchi_tob = 0;
//   par.tob_chi = 0.;
//   par.tob_chi0 = 0.;

//   LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);
//   do{
//     int bn = tob_loop.get_beam();
//     //int dn = tob_loop.get_zpos();
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tob = nchi_tob + 1;
//     par.AT_chi += tob_loop.GetEntry(res2);
//     par.AT_chi0 += tob_loop.GetEntry(dif2);
//     par.tob_chi += tob_loop.GetEntry(res2);
//     par.tob_chi0 += tob_loop.GetEntry(dif2);
//     par.b_chi[bn] += tob_loop.GetEntry(res2);
//     par.b_chi0[bn] += tob_loop.GetEntry(dif2);
//     }while(tob_loop.next());

//   par.tob_chi = sqrt(par.tob_chi/nchi_tob);
//   par.tob_chi0 = sqrt(par.tob_chi0/nchi_tob);
//   //cout<<"irun: "<<irun<<"  nchi_tob: "<<nchi_tob<<"  chi_tob: "<<chi_tob[irun]<<endl;  
 
//   par.AT_chi = sqrt(par.AT_chi/(nchi-22));
//   par.AT_chi0 = sqrt(par.AT_chi0/nchi);

 
//   for(int bnn = 0; bnn < 8;  ++bnn){
//     par.b_chi[bnn] = sqrt(par.b_chi[bnn]/(nchib[bnn]-2));
//     par.b_chi0[bnn] = sqrt(par.b_chi0[bnn]/nchib[bnn]);
//   }
  
// }

LASGlobalData<int> group_mask(LAS::beam_group beam_group)
{
  LASGlobalData<int> result(0);
  LASGlobalDataLoop loop(convert_to_loop_type(beam_group));
  do{ loop.GetEntry( result ) = 1;}while(loop.next());
  return result;
}

void tweak_tecm_data(LASGlobalData<double>& data, double tweak_tecm_tor, double tweak_tecm_shx, double tweak_tecm_shy)
{
  LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_MINUS_AT);
  do{
    loop.GetEntry(data) += loop.GetEntry(zpos()) * tweak_tecm_tor;
    loop.GetEntry(data) += loop.GetEntry(zpos()) * cos(loop.GetEntry(theta())) * tweak_tecm_shx;
    loop.GetEntry(data) -= loop.GetEntry(zpos()) * sin(loop.GetEntry(theta())) * tweak_tecm_shy;
  }while(loop.next());
}

void mask_tecm_disc(LASGlobalData<int>& mask, int disc)
{
  LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_MINUS_AT);
  do{
    if(loop.get_zpos() == disc) loop.GetEntry(mask) =0;
  }while(loop.next());
}

void AT_residuals(const std::string& ref_filename, const std::string& filename, unsigned int blnr, LAS::FIT_TYPE::fit_type fit_type, double tweak_tecm_tor, double tweak_tecm_shx, double tweak_tecm_shy)
{
  // Get the reference positions
  std::cout << "Reference file is: " << ref_filename << std::endl;
  LASGlobalData<double> ref_pos;
  LASGlobalData<double> ref_err;
  LASGlobalData<int> ref_mask;
  get_ref_pos(ref_filename, ref_pos, ref_err, ref_mask);
  
  Avec good_block = avec_get("good_block",filename);
  
  if( good_block[blnr] == 0.0 || blnr > good_block.size() - 1){
    std::cerr << "Block is marked as not good or out of range" << std::endl ;
    std::cerr << "Selected block: " << blnr << "  Total number of blocks: " << good_block.size() << std::endl;
    return;
  }

  std::ostringstream blk_nr;
  blk_nr << blnr;
  std::string blk_str = blk_nr.str();

  LASGlobalData<double> pos = global_data_get<double>("positions_" + blk_str, filename);
  LASGlobalData<double> rms = global_data_get<double>("pos_error_" + blk_str, filename);
  LASGlobalData<int> norm = global_data_get<int>("norm_" + blk_str, filename);
  LASGlobalData<double> dif = (pos - ref_pos)*pitch()/r0();
  correct_signs(dif);
  LASGlobalData<double> error = sqrt(rms*rms/norm + ref_err*ref_err)*pitch()/r0();
  LASGlobalData<int> mask = fabs(dif) < 1e-3;;

  mask_known_bad_modules(mask);

  LasAlPar results;
  LASGlobalData<double> rec_pos;
  switch(fit_type){
  case LAS::FIT_TYPE::TEC_PLUS_AT:
    AT_TECP2TOB_3par(dif, error, ref_mask && mask, results.at);
    rec_pos = TECP_spot_rec(results.at);
    mask &= group_mask(LAS::TOB) || group_mask(LAS::TEC_PLUS_AT); 
    break;
  case LAS::FIT_TYPE::TEC_MINUS_AT:
    //tweak_tecm_data(dif, tweak_tecm_tor, tweak_tecm_shx, tweak_tecm_shy);
    //mask_tecm_disc( mask, 3);
    //mask_tecm_disc( mask, 2);
    //mask_tecm_disc( mask, 1);
    AT_TECM2TOB_3par(dif, error, ref_mask && mask, results.at);
    rec_pos = TECM_spot_rec(results.at);
    mask &= group_mask(LAS::TOB) || group_mask(LAS::TEC_MINUS_AT); 
    break;
  case LAS::FIT_TYPE::TIB_OLD:
    AT_TIB2TOB_6par(dif, error, ref_mask && mask, results.at);
    //AlPar_fit(dif, error, ref_mask && mask, results.at, true, true);
    TIB_chi2( dif, error, ref_mask && mask, results.at );
    rec_pos = TIB_spot_rec(results.at);
    mask &= group_mask(LAS::TOB) || group_mask(LAS::TIB); 
    break;
  case LAS::FIT_TYPE::ALL:
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true); 
    break;
  case LAS::FIT_TYPE::AT:
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true);
    mask &= group_mask(LAS::AT); 
    break;
  case LAS::FIT_TYPE::TOB:
    //AT_TOB_16par(dif, error*error, ref_mask && mask, results.at, rec_pos);
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true);
    mask &= group_mask(LAS::TOB); 
    break;
  case LAS::FIT_TYPE::TIB:
    //AT_TIB_22par(dif, error*error, ref_mask && mask, results.at, rec_pos, true);
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true);
    mask &= group_mask(LAS::TOB) || group_mask(LAS::TIB); 
    break;
  case LAS::FIT_TYPE::TEC_MINUS:
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true);
    mask &= group_mask(LAS::TEC_MINUS); 
    break;
  case LAS::FIT_TYPE::TEC_MINUS_R4:
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true);
    mask &= group_mask(LAS::TEC_MINUS_R4); 
    break;
  case LAS::FIT_TYPE::TEC_MINUS_R6:
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true);
    mask &= group_mask(LAS::TEC_MINUS_R6); 
    break;
  case LAS::FIT_TYPE::TEC_PLUS:
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true);
    mask &= group_mask(LAS::TEC_PLUS); 
    break;
  case LAS::FIT_TYPE::TEC_PLUS_R4:
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true);
    mask &= group_mask(LAS::TEC_PLUS_R4); 
    break;
  case LAS::FIT_TYPE::TEC_PLUS_R6:
    LAS_fit(dif, error*error, ref_mask && mask, results, rec_pos, fit_type, true);
    mask &= group_mask(LAS::TEC_PLUS_R6); 
    break;
  default:
    AT_28par(dif, error, ref_mask && mask, results.at);
    rec_pos = AT_spot_rec(results.at);
    mask &= group_mask(LAS::AT);
  }
  
  LASGlobalData<double> residuals = rec_pos - dif;
  if(fit_type == LAS::FIT_TYPE::ALL || fit_type == LAS::FIT_TYPE::AT || fit_type == LAS::FIT_TYPE::TEC_PLUS_AT || fit_type == LAS::FIT_TYPE::TEC_MINUS_AT || fit_type == LAS::FIT_TYPE::TIB || fit_type == LAS::FIT_TYPE::TOB || fit_type == LAS::FIT_TYPE::TIB_OLD ){
    draw_global_data_z(ref_pos, error, ref_mask && mask, LAS::ALL, "cv_ref_pos");
    draw_global_data_z(pos, error, ref_mask && mask, LAS::ALL, "cv_abs_pos");
    draw_global_data_z(dif, error, ref_mask && mask, LAS::ALL, "cv_dif");
    draw_global_data_z(residuals, error, ref_mask && mask, LAS::ALL, "cv_residuals");
    results.at.print();
  }
  if(fit_type == LAS::FIT_TYPE::TEC_PLUS || fit_type == LAS::FIT_TYPE::TEC_PLUS_R4 || fit_type == LAS::FIT_TYPE::TEC_PLUS_R6 || fit_type == LAS::FIT_TYPE::TEC_MINUS || fit_type == LAS::FIT_TYPE::TEC_MINUS_R4 || fit_type == LAS::FIT_TYPE::TEC_MINUS_R6){
    draw_global_data_z(dif, error, ref_mask && mask, LAS::TEC, "cv_positions");
    draw_global_data_z(residuals, error, ref_mask && mask, LAS::TEC, "cv_residuals");
    results.at.print();
  }
}

// void AT_chi2(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask,  AtPar& par) //calculation of chi2/ndf at TIB&TOB
// {

//   LASGlobalData<double> recpos = AT_spot_rec(par);
//   LASGlobalData<double> resid = recpos - dif;
//   LASGlobalData<double> res2 = pow(resid/(err2),2);
//   LASGlobalData<double> dif2 = pow(dif/(err2),2);

//   int nchi_tecp,nchi_tecm,nchi_tib,nchi_tob,nchi,nchib[8];
 
//   nchi = 0;
//   par.AT_chi = 0.;
//   par.AT_chi0 = 0.;

//   for(int bnn = 0; bnn < 8;  ++bnn){
//     nchib[bnn] = 0; 
//     par.b_chi[bnn] = 0.;
//     par.b_chi0[bnn] = 0;
//   } 

//   nchi_tecp = 0;
//   nchi_tecm = 0;
//   nchi_tib = 0;
//   nchi_tob = 0;
//   par.tecp_chi = 0.;
//   par.tecp_chi0 = 0.;
//   par.tecm_chi = 0.;
//   par.tecm_chi0 = 0.;
//   par.tib_chi = 0.;
//   par.tib_chi0 = 0.;
//   par.tob_chi = 0.;
//   par.tob_chi0 = 0.;

//   LASGlobalDataLoop tecp_loop(LASGlobalDataLoop::TEC_PLUS_AT);
//   do{
//     int bn = tecp_loop.get_beam();
//     //int dn = tecp_loop.get_zpos();
//     //cout<<"  bn= "<<bn<<"  dn= "<<dn<<"  res: "<<tecp_loop.GetEntry(resid)<<"  dif: "<<tecp_loop.GetEntry(dif)<<endl;  
//     if(tecp_loop.GetEntry(mask) != 1 ) continue;
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tecp = nchi_tecp + 1;
//     par.AT_chi += tecp_loop.GetEntry(res2);
//     par.AT_chi0 += tecp_loop.GetEntry(dif2);
//     par.tecp_chi += tecp_loop.GetEntry(res2);
//     par.tecp_chi0 += tecp_loop.GetEntry(dif2);
//     par.b_chi[bn] += tecp_loop.GetEntry(res2);
//     par.b_chi0[bn] += tecp_loop.GetEntry(dif2);
//   }while(tecp_loop.next());
 
  
//   if(nchi_tecp > 0){
//     par.tecp_chi = sqrt(par.tecp_chi/nchi_tecp);
//     par.tecp_chi0 = sqrt(par.tecp_chi0/nchi_tecp);
//   }
//   //cout<<"  nchi_tecp: "<<nchi_tecp<<"  tecp_chi: "<<par.tecp_chi<<" tecp_chi0: "<<par.tecp_chi0<<endl;  

//   LASGlobalDataLoop tecm_loop(LASGlobalDataLoop::TEC_MINUS_AT);
//   do{
//     //int dn = tecm_loop.get_zpos();
//     //cout<<"  bn= "<<bn<<"  dn= "<<dn<<"  res: "<<tecm_loop.GetEntry(resid)<<"  dif: "<<tecm_loop.GetEntry(dif)<<endl;  
//     if(tecm_loop.GetEntry(mask) != 1 ) continue;
//     int bn = tecm_loop.get_beam();
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tecm = nchi_tecm + 1;
//     par.AT_chi += tecm_loop.GetEntry(res2);
//     par.AT_chi0 += tecm_loop.GetEntry(dif2);
//     par.tecm_chi += tecm_loop.GetEntry(res2);
//     par.tecm_chi0 += tecm_loop.GetEntry(dif2);
//     par.b_chi[bn] += tecm_loop.GetEntry(res2);
//     par.b_chi0[bn] += tecm_loop.GetEntry(dif2);
    
//   }while(tecm_loop.next());
 
  
//   if(nchi_tecm > 0){
//     par.tecm_chi = sqrt(par.tecm_chi/nchi_tecm);
//     par.tecm_chi0 = sqrt(par.tecm_chi0/nchi_tecm);
//   }
//   //cout<<"  nchi_tecm: "<<nchi_tecm<<"  tecm_chi: "<<par.tecm_chi<<" tecm_chi0: "<<par.tecm_chi0<<endl;  

//   LASGlobalDataLoop tib_loop(LASGlobalDataLoop::TIB);
//   do{
//     if( tib_loop.GetEntry(mask) != 1 ) continue;
//     int bn = tib_loop.get_beam();
// //     int dn = tecp_loop.get_zpos();
// //     cout<<"  bn= "<<bn<<"  dn= "<<dn<<"  res: "<<tib_loop.GetEntry(resid)<<"  dif: "<<tib_loop.GetEntry(dif)<<endl;  
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tib = nchi_tib + 1;
//     par.AT_chi += tib_loop.GetEntry(res2);
//     par.AT_chi0 += tib_loop.GetEntry(dif2);
//     par.tib_chi += tib_loop.GetEntry(res2);
//     par.tib_chi0 += tib_loop.GetEntry(dif2);
//     par.b_chi[bn] += tib_loop.GetEntry(res2);
//     par.b_chi0[bn] += tib_loop.GetEntry(dif2);
//   }while(tib_loop.next());

//   if(nchi_tib > 0){
//     par.tib_chi = sqrt(par.tib_chi/nchi_tib);
//     par.tib_chi0 = sqrt(par.tib_chi0/nchi_tib);
//   }
//   //cout<<"  nchi_tib: "<<nchi_tib<<"  tib_chi: "<<par.tib_chi<<"  tib_chi0: "<<par.tib_chi0<<endl;  
 
//   LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);

//   do{
//     if( tob_loop.GetEntry(mask) != 1 ) continue;
//     int bn = tob_loop.get_beam();
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tob = nchi_tob + 1;
//     par.AT_chi += tob_loop.GetEntry(res2);
//     par.AT_chi0 += tob_loop.GetEntry(dif2);
//     par.tob_chi += tob_loop.GetEntry(res2);
//     par.tob_chi0 += tob_loop.GetEntry(dif2);
//     par.b_chi[bn] += tob_loop.GetEntry(res2);
//     par.b_chi0[bn] += tob_loop.GetEntry(dif2);
//   }while(tob_loop.next());

//   if(nchi_tob > 0){
//     par.tob_chi = sqrt(par.tob_chi/nchi_tob);
//     par.tob_chi0 = sqrt(par.tob_chi0/nchi_tob);
//   }
//   //cout<<"  nchi_tob: "<<nchi_tob<<"  tob_chi: "<<par.tob_chi<<"  tob_chi0: "<<par.tob_chi0<<endl;  
 
//   for(int bnn = 0; bnn < 8;  ++bnn){
//     if(nchib[bnn] > 2){
//       par.b_chi[bnn] = sqrt(par.b_chi[bnn]/(nchib[bnn]-2));
//       par.b_chi0[bnn] = sqrt(par.b_chi0[bnn]/nchib[bnn]);
//     }
//     //cout<<" bn= "<<bnn<<"  nchi_b: "<<nchib[bnn]<<"  b_chi: "<<par.b_chi[bnn]<<" b_chi0: "<<par.b_chi0[bnn]<<endl;  
//   }
  
//   if(nchi > 28){
//     par.AT_chi = sqrt(par.AT_chi/(nchi-28));
//     par.AT_chi0 = sqrt(par.AT_chi0/nchi);
//   }
//   //cout<<"  nchi_AT: "<<nchi<<"  AT_chi: "<<par.AT_chi<<"  AT_chi0: "<<par.AT_chi0<<endl;  
  
// }

// // reconstruction of TEC+, TEC-, TIB wrt TOB with 28 alignment parameters( 16 beam parameters: (a&b)*8; 3 TEC+ parameters: dxp, dyp, rzp; 
// //3 TEC- parameters: dxm, dym, rzm; 6 TIB parameters:dx, dy, rx, ry, rz, tz;)
// void AT_28par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results)
// {

//   // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

//   double bs[28];
//   double db[28][4][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..27 (parameter number); j=0,3 (TEC+ ,TEC-, TIB, TOB); k=0..7 (beam number); l=0..5 (module number)

//   int npb[3];
//   npb[0] = 16;
//   npb[1] = 19;
//   npb[2] = 22;

//   int mdn[4],md;
//   mdn[0] = 4; // number of TEC+ discs
//   mdn[1] = 4; // number of TEC1 discs
//   mdn[2] = 6; // number of TIB discs
//   mdn[3] = 6; // number of TOB discs

 
//   LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
//   LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
//   LASGlobalData<double> par2 = dif * cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par3 = dif * sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par4 = dif / err2;
//   LASGlobalData<double> par5 = dif * zpos()/ err2;


//   LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
//   LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
//   LASGlobalData<double> dev2 = cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev3 = sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> dev5 = zpos()/ err2;

//   for(int ip = 0; ip < 28;  ++ip){ 
//     bs[ip] = 0.;
//   }
 
//   for(int bn = 0; bn < 8;  ++bn){
//     for(int id = 0; id < 4;  ++id){ 
//       for(int dn = 0; dn < 6;  ++dn){
// 	for(int ip = 0; ip < 28;  ++ip){ 
// 	  db[ip][id][bn][dn] = 0.;
// 	}
//       }
//     }
//   }

//   for(int bn = 0; bn < 8;  ++bn){
//     for(int id = 0; id < 4;  ++id){ 
//       md = mdn[id];
//       for(int dn = 0; dn < md;  ++dn){
	
// 	if(id < 2){
// 	  if(err2.GetEntry(id,-1,bn,dn) < 1000. && mask.GetEntry(id, -1, bn, dn)){
	    
// 	    bs[2*bn] += par5.GetEntry(id,-1,bn,dn);
// 	    bs[2*bn+1] += par4.GetEntry(id,-1,bn,dn);
	    
// 	    bs[npb[id]+0] += par0.GetEntry(id,-1,bn,dn);
// 	    bs[npb[id]+1] += par1.GetEntry(id,-1,bn,dn);
// 	    bs[npb[id]+2] += par4.GetEntry(id,-1,bn,dn);
	    
// 	    db[2*bn][id][bn][dn] = dev5.GetEntry(id,-1,bn,dn);
// 	    db[2*bn+1][id][bn][dn] = dev4.GetEntry(id,-1,bn,dn);
	    
// 	    db[npb[id]+0][id][bn][dn] = dev0.GetEntry(id,-1,bn,dn);
// 	    db[npb[id]+1][id][bn][dn] = dev1.GetEntry(id,-1,bn,dn);
// 	    db[npb[id]+2][id][bn][dn] = dev4.GetEntry(id,-1,bn,dn);
// 	  }
// 	}
// 	if(id == 2){
// 	  if(err2.GetEntry(id,-1,bn,dn) < 1000. && mask.GetEntry(id, -1, bn, dn)){
	    
// 	    bs[2*bn] += par5.GetEntry(id,-1,bn,dn);
// 	    bs[2*bn+1] += par4.GetEntry(id,-1,bn,dn);
	    
// 	    bs[npb[id]+0] += par0.GetEntry(id,-1,bn,dn);
// 	    bs[npb[id]+1] += par1.GetEntry(id,-1,bn,dn);
// 	    bs[npb[id]+2] += par2.GetEntry(id,-1,bn,dn);
// 	    bs[npb[id]+3] += par3.GetEntry(id,-1,bn,dn);
// 	    bs[npb[id]+4] += par4.GetEntry(id,-1,bn,dn);
// 	    bs[npb[id]+5] += par5.GetEntry(id,-1,bn,dn);
	    
// 	    db[2*bn][id][bn][dn] = dev5.GetEntry(id,-1,bn,dn);
// 	    db[2*bn+1][id][bn][dn] = dev4.GetEntry(id,-1,bn,dn);
	    
// 	    db[npb[id]+0][id][bn][dn] = dev0.GetEntry(id,-1,bn,dn);
// 	    db[npb[id]+1][id][bn][dn] = dev1.GetEntry(id,-1,bn,dn);
// 	    db[npb[id]+2][id][bn][dn] = dev2.GetEntry(id,-1,bn,dn);
// 	    db[npb[id]+3][id][bn][dn] = dev3.GetEntry(id,-1,bn,dn);
// 	    db[npb[id]+4][id][bn][dn] = dev4.GetEntry(id,-1,bn,dn);
// 	    db[npb[id]+5][id][bn][dn] = dev5.GetEntry(id,-1,bn,dn);
// 	  }
// 	}
      
// 	if(id == 3){
// 	  if(err2.GetEntry(id,-1,bn,dn) < 1000. && mask.GetEntry(id, -1, bn, dn)){
	    
// 	    bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
// 	    bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	    
// 	    db[2*bn][3][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
// 	    db[2*bn+1][3][bn][dn] = dev4.GetEntry(3,-1,bn,dn);
// 	  }
// 	}
//       }
//     }
//   }
   
//   TMatrixT<double> h(28,28);
  
//   for(int ix = 0; ix < 28; ++ix){
//     for(int iy = 0; iy < 28; ++iy){
//       h(ix,iy) = 0.;
//     }
//   }
  
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


//   for(int bn = 0; bn < 8;  ++bn){
//     for(int id = 0; id < 4;  ++id){ 
//       md = mdn[id];
//       for(int dn = 0; dn < md;  ++dn){
	
// 	if(id < 2){
// 	  if(err2.GetEntry(id,-1,bn,dn) < 1000. && mask.GetEntry(id, -1, bn, dn)){
	    
// 	    h(2*bn,2*bn)             += haa.GetEntry(id,-1,bn,dn);
// 	    h(2*bn,2*bn+1)           += hab.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,2*bn)           += hba.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,2*bn+1)         += hbb.GetEntry(id,-1,bn,dn);
	    
// 	    h(npb[id]+0,npb[id]+0)   += h00.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,npb[id]+1)   += h01.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,npb[id]+2)   += h04.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,2*bn)        += h05.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,2*bn+1)      += h04.GetEntry(id,-1,bn,dn);
	    
// 	    h(npb[id]+1,npb[id]+0)   += h10.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,npb[id]+1)   += h11.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,npb[id]+2)   += h14.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,2*bn)        += h15.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,2*bn+1)      += h14.GetEntry(id,-1,bn,dn);
	    
// 	    h(npb[id]+2,npb[id]+0)   += h40.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,npb[id]+1)   += h41.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,npb[id]+2)   += h44.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,2*bn)        += h45.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,2*bn+1)      += h44.GetEntry(id,-1,bn,dn);
	    
// 	    h(2*bn,npb[id]+0)        += h50.GetEntry(id,-1,bn,dn);
// 	    h(2*bn,npb[id]+1)        += h51.GetEntry(id,-1,bn,dn);
// 	    h(2*bn,npb[id]+2)        += h54.GetEntry(id,-1,bn,dn);
	    
// 	    h(2*bn+1,npb[id]+0)      += h40.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,npb[id]+1)      += h41.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,npb[id]+2)      += h44.GetEntry(id,-1,bn,dn);
// 	  }
// 	}
	
// 	if(id == 2){
// 	  if(err2.GetEntry(id,-1,bn,dn) < 1000. && mask.GetEntry(id, -1, bn, dn)){
	    
// 	    h(2*bn,2*bn)             += haa.GetEntry(id,-1,bn,dn);
// 	    h(2*bn,2*bn+1)           += hab.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,2*bn)           += hba.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,2*bn+1)         += hbb.GetEntry(id,-1,bn,dn);
	    
// 	    h(npb[id]+0,npb[id]+0)   += h00.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,npb[id]+1)   += h01.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,npb[id]+2)   += h02.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,npb[id]+3)   += h03.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,npb[id]+4)   += h04.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,npb[id]+5)   += h05.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,2*bn)        += h05.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+0,2*bn+1)      += h04.GetEntry(id,-1,bn,dn);
	    
// 	    h(npb[id]+1,npb[id]+0)   += h10.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,npb[id]+1)   += h11.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,npb[id]+2)   += h12.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,npb[id]+3)   += h13.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,npb[id]+4)   += h14.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,npb[id]+5)   += h15.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,2*bn)        += h15.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+1,2*bn+1)      += h14.GetEntry(id,-1,bn,dn);
	    
// 	    h(npb[id]+2,npb[id]+0)   += h20.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,npb[id]+1)   += h21.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,npb[id]+2)   += h22.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,npb[id]+3)   += h23.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,npb[id]+4)   += h24.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,npb[id]+5)   += h25.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,2*bn)        += h25.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+2,2*bn+1)      += h24.GetEntry(id,-1,bn,dn);
	    
// 	    h(npb[id]+3,npb[id]+0)   += h30.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+3,npb[id]+1)   += h31.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+3,npb[id]+2)   += h32.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+3,npb[id]+3)   += h33.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+3,npb[id]+4)   += h34.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+3,npb[id]+5)   += h35.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+3,2*bn)        += h35.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+3,2*bn+1)      += h34.GetEntry(id,-1,bn,dn);
	    
// 	    h(npb[id]+4,npb[id]+0)   += h40.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+4,npb[id]+1)   += h41.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+4,npb[id]+2)   += h42.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+4,npb[id]+3)   += h43.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+4,npb[id]+4)   += h44.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+4,npb[id]+5)   += h45.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+4,2*bn)        += h45.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+4,2*bn+1)      += h44.GetEntry(id,-1,bn,dn);
	    
// 	    h(npb[id]+5,npb[id]+0)   += h50.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+5,npb[id]+1)   += h51.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+5,npb[id]+2)   += h52.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+5,npb[id]+3)   += h53.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+5,npb[id]+4)   += h54.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+5,npb[id]+5)   += h55.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+5,2*bn)        += h55.GetEntry(id,-1,bn,dn);
// 	    h(npb[id]+5,2*bn+1)      += h54.GetEntry(id,-1,bn,dn);
	    
	    
// 	    h(2*bn,npb[id]+0)        += h50.GetEntry(id,-1,bn,dn);
// 	    h(2*bn,npb[id]+1)        += h51.GetEntry(id,-1,bn,dn);
// 	    h(2*bn,npb[id]+2)        += h52.GetEntry(id,-1,bn,dn);
// 	    h(2*bn,npb[id]+3)        += h53.GetEntry(id,-1,bn,dn);
// 	    h(2*bn,npb[id]+4)        += h54.GetEntry(id,-1,bn,dn);
// 	    h(2*bn,npb[id]+5)        += h55.GetEntry(id,-1,bn,dn);
	    
// 	    h(2*bn+1,npb[id]+0)      += h40.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,npb[id]+1)      += h41.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,npb[id]+2)      += h42.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,npb[id]+3)      += h43.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,npb[id]+4)      += h44.GetEntry(id,-1,bn,dn);
// 	    h(2*bn+1,npb[id]+5)      += h45.GetEntry(id,-1,bn,dn);
// 	  }
// 	}
	
// 	if(id == 3){
// 	  if(err2.GetEntry(id,-1,bn,dn) < 1000. && mask.GetEntry(id, -1, bn, dn)){
	    
// 	    h(2*bn,2*bn)             += haa.GetEntry(3,-1,bn,dn);
// 	    h(2*bn,2*bn+1)           += hab.GetEntry(3,-1,bn,dn);
// 	    h(2*bn+1,2*bn)           += hba.GetEntry(3,-1,bn,dn);
// 	    h(2*bn+1,2*bn+1)         += hbb.GetEntry(3,-1,bn,dn);	
// 	  }
// 	}
//       }
//     }
//   }  

  
//   h.Invert();

//   double cp[28][4][8][6];
//   double up[28][28];
//   double ps[28];

//   for(int ip = 0; ip < 28; ++ip){
//     for(int id = 0; id < 4; ++id){
//       for(int bn = 0; bn < 8;  ++bn){
// 	for(int dn = 0; dn < mdn[id];  ++dn){
// 	  cp[ip][id][bn][dn] = 0.;
// 	  for(int ix = 0; ix < 28; ++ix){
// 	    cp[ip][id][bn][dn] = cp[ip][id][bn][dn] + h[ip][ix]*db[ix][id][bn][dn];
// 	    //std::cout <<"ip="<< ip <<" id="<< id <<" bn="<< bn <<" dn="<< dn <<" cp="<< cp[ip][id][bn][dn]  << std::endl;
// 	  }
// 	}
//       }
//     }
//   }

//   for(int ip = 0; ip < 28; ++ip){
//     for(int ix = 0; ix < 28; ++ix){
//       up[ip][ix] = 0.;
//       for(int id = 0; id < 4; ++id){
// 	for(int bn = 0; bn < 8;  ++bn){
// 	  for(int dn = 0; dn < mdn[id];  ++dn){
// 	    up[ip][ix] = up[ip][ix] + cp[ip][id][bn][dn]*cp[ix][id][bn][dn]*err2.GetEntry(id,-1,bn,dn);
// 	  }
// 	}
//       }
//     }
//   }
 
//   for(int bn = 0; bn < 8; ++bn){
//     results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
//     results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
//   }

//   int id = 0;
//   results.er_tecp_dx = sqrt(up[npb[id]+0][npb[id]+0]);
//   results.er_tecp_dy = sqrt(up[npb[id]+1][npb[id]+1]);
//   results.er_tecp_rz = sqrt(up[npb[id]+2][npb[id]+2]);
// //     std::cout <<"err tecp_dx "<< up[npb[id]+0][npb[id]+0] << std::endl;
// //     std::cout <<"err tecp_dy "<< up[npb[id]+1][npb[id]+1] << std::endl;
// //     std::cout <<"err tecp_rz "<< up[npb[id]+2][npb[id]+2] << std::endl;

//   id = 1;
//   results.er_tecm_dx = sqrt(up[npb[id]+0][npb[id]+0]);
//   results.er_tecm_dy = sqrt(up[npb[id]+1][npb[id]+1]);
//   results.er_tecm_rz = sqrt(up[npb[id]+2][npb[id]+2]);

//   id = 2;
//   results.er_tib_dx = sqrt(up[npb[id]+0][npb[id]+0]);
//   results.er_tib_dy = sqrt(up[npb[id]+1][npb[id]+1]);
//   results.er_tib_rx = sqrt(up[npb[id]+2][npb[id]+2]);
//   results.er_tib_ry = sqrt(up[npb[id]+3][npb[id]+3]);
//   results.er_tib_rz = sqrt(up[npb[id]+4][npb[id]+4]);
//   results.er_tib_tz = sqrt(up[npb[id]+5][npb[id]+5]);
//   results.tib_dxdy  = up[npb[id]+0][npb[id]+1]/results.er_tib_dx/results.er_tib_dy;
//   results.tib_dxrx  = up[npb[id]+0][npb[id]+2]/results.er_tib_dx/results.er_tib_rx;
//   results.tib_dxry  = up[npb[id]+0][npb[id]+3]/results.er_tib_dx/results.er_tib_ry;
//   results.tib_dxrz  = up[npb[id]+0][npb[id]+4]/results.er_tib_dx/results.er_tib_rz;
//   results.tib_dxtz  = up[npb[id]+0][npb[id]+5]/results.er_tib_dx/results.er_tib_tz;
//   results.tib_dyrx  = up[npb[id]+1][npb[id]+2]/results.er_tib_dy/results.er_tib_rx;
//   results.tib_dyry  = up[npb[id]+1][npb[id]+3]/results.er_tib_dy/results.er_tib_ry;
//   results.tib_dyrz  = up[npb[id]+1][npb[id]+4]/results.er_tib_dy/results.er_tib_rz;
//   results.tib_dytz  = up[npb[id]+1][npb[id]+5]/results.er_tib_dy/results.er_tib_tz;
//   results.tib_rxry  = up[npb[id]+2][npb[id]+3]/results.er_tib_rx/results.er_tib_ry;
//   results.tib_rxrz  = up[npb[id]+2][npb[id]+4]/results.er_tib_rx/results.er_tib_rz;
//   results.tib_rxtz  = up[npb[id]+2][npb[id]+5]/results.er_tib_rx/results.er_tib_tz;
//   results.tib_ryrz  = up[npb[id]+3][npb[id]+4]/results.er_tib_ry/results.er_tib_rz;
//   results.tib_rytz  = up[npb[id]+3][npb[id]+5]/results.er_tib_ry/results.er_tib_tz;
//   results.tib_rztz  = up[npb[id]+4][npb[id]+5]/results.er_tib_rz/results.er_tib_tz;

//   for(int ip = 0; ip < 28; ++ip){
//     ps[ip] = 0.;
//     for(int ix = 0; ix < 28; ++ix){
//       ps[ip] += h(ip,ix)*bs[ix];
//     }
//     //std::cout <<"ip= "<<ip<< "  ps =  " << ps[ip] << std::endl;
//   }

//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = ps[2*bn];
//     results.beam_b[bn] = ps[2*bn+1];
//   }

//   id = 0;
//   results.tecp_dx = ps[npb[id]+0];
//   results.tecp_dy = ps[npb[id]+1];
//   results.tecp_rz = ps[npb[id]+2];

//   id = 1;
//   results.tecm_dx = ps[npb[id]+0];
//   results.tecm_dy = ps[npb[id]+1];
//   results.tecm_rz = ps[npb[id]+2];

//   id = 2;
//   results.tib_dx = ps[npb[id]+0];
//   results.tib_dy = ps[npb[id]+1];
//   results.tib_rx = ps[npb[id]+2];
//   results.tib_ry = ps[npb[id]+3];
//   results.tib_rz = ps[npb[id]+4];
//   results.tib_tz = ps[npb[id]+5];

// }

// // reconstruction of TIB wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
// void AT_TIB2TOB_6par_old(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results)
// {
//   // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

//   //LASGlobalData<double> nzpos = zpos();
//   LASGlobalData<double> nzpos = zpos() - LAS::z_at_bs;

//   double bs[22];
//   double db[22][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TIB, TOB); k=0..7 (beam number); l=0..5 (module number)

//   int npb = 16;// number of beam parameters
 
//   LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
//   LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
//   LASGlobalData<double> dev2 = cos(theta()) * r0() * nzpos / err2;
//   LASGlobalData<double> dev3 = sin(theta()) * r0() * nzpos / err2;
//   LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> dev5 = nzpos / err2;

// //   LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
// //   LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
// //   LASGlobalData<double> par2 = dif * cos(theta()) * r0() * nzpos / err2;
// //   LASGlobalData<double> par3 = dif * sin(theta()) * r0() * nzpos / err2;
// //   LASGlobalData<double> par4 = dif / err2;
// //   LASGlobalData<double> par5 = dif * nzpos/ err2;

//   LASGlobalData<double> par0 = dif * dev0;
//   LASGlobalData<double> par1 = dif * dev1;
//   LASGlobalData<double> par2 = dif * dev2;
//   LASGlobalData<double> par3 = dif * dev3;
//   LASGlobalData<double> par4 = dif * dev4;
//   LASGlobalData<double> par5 = dif * dev5;


//   for(int ip = 0; ip < 22;  ++ip){ 
//     bs[ip] = 0.;
//   }
//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 6;  ++dn){
//       for(int ip = 0; ip < 22;  ++ip){ 
// 	db[ip][0][bn][dn] = 0.;
// 	db[ip][1][bn][dn] = 0.;
//       }
//       if(err2.GetEntry(2,-1,bn,dn) < 1000.){
	    
// 	bs[2*bn] += par5.GetEntry(2,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(2,-1,bn,dn);
	
// 	bs[npb+0] += par0.GetEntry(2,-1,bn,dn);
// 	bs[npb+1] += par1.GetEntry(2,-1,bn,dn);
// 	bs[npb+2] += par2.GetEntry(2,-1,bn,dn);
// 	bs[npb+3] += par3.GetEntry(2,-1,bn,dn);
// 	bs[npb+4] += par4.GetEntry(2,-1,bn,dn);
// 	bs[npb+5] += par5.GetEntry(2,-1,bn,dn);
        
// 	db[2*bn][0][bn][dn] = dev5.GetEntry(2,-1,bn,dn);
// 	db[2*bn+1][0][bn][dn] = dev4.GetEntry(2,-1,bn,dn);
	
// 	db[npb+0][0][bn][dn] = dev0.GetEntry(2,-1,bn,dn);
// 	db[npb+1][0][bn][dn] = dev1.GetEntry(2,-1,bn,dn);
// 	db[npb+2][0][bn][dn] = dev2.GetEntry(2,-1,bn,dn);
// 	db[npb+3][0][bn][dn] = dev3.GetEntry(2,-1,bn,dn);
// 	db[npb+4][0][bn][dn] = dev4.GetEntry(2,-1,bn,dn);
// 	db[npb+5][0][bn][dn] = dev5.GetEntry(2,-1,bn,dn);
//       }
//       if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
//  	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
//  	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
//  	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
//  	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);

// 	//bs[2*bn] -= par5.GetEntry(3,-1,bn,dn);
// 	//bs[2*bn+1] -= par4.GetEntry(3,-1,bn,dn);
	
// 	//db[2*bn][1][bn][dn] = -dev5.GetEntry(3,-1,bn,dn);
// 	//db[2*bn+1][1][bn][dn] = -dev4.GetEntry(3,-1,bn,dn);
//       }
//     }
//   }
   
//   TMatrixT<double> h(22,22);
  
//   for(int ix = 0; ix < 22; ++ix){
//     for(int iy = 0; iy < 22; ++iy){
//       h(ix,iy) = 0.;
//     }
//   }
  
//   LASGlobalData<double> haa =   nzpos  * nzpos  / err2;
//   LASGlobalData<double> hab =   nzpos  / err2;
//   LASGlobalData<double> hba =   nzpos  / err2;
//   LASGlobalData<double> hbb =   LASGlobalData<double>(1.) / err2;


//   LASGlobalData<double> h00 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) / err2;
//   LASGlobalData<double> h01 =   sin(theta()) * cos(theta()) / err2;
//   LASGlobalData<double> h02 =   sin(theta()) * cos(theta()) * nzpos  / err2;
//   LASGlobalData<double> h03 =   sin(theta()) * sin(theta()) * nzpos  / err2;
//   LASGlobalData<double> h04 =   sin(theta()) / err2;
//   LASGlobalData<double> h05 =   sin(theta()) * nzpos  / err2;

//   LASGlobalData<double> h10 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) / err2;
//   LASGlobalData<double> h11 =   cos(theta()) * cos(theta()) / err2;
//   LASGlobalData<double> h12 =   cos(theta()) * cos(theta()) * nzpos  / err2;
//   LASGlobalData<double> h13 =   cos(theta()) * sin(theta()) * nzpos  / err2;
//   LASGlobalData<double> h14 =   cos(theta()) / err2;
//   LASGlobalData<double> h15 =   cos(theta()) * nzpos  / err2;

//   LASGlobalData<double> h20 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) * nzpos  / err2;
//   LASGlobalData<double> h21 =   cos(theta()) * cos(theta()) * nzpos  / err2;
//   LASGlobalData<double> h22 =   cos(theta()) * cos(theta()) * nzpos  * nzpos  / err2;
//   LASGlobalData<double> h23 =   cos(theta()) * sin(theta()) * nzpos  * nzpos  / err2;
//   LASGlobalData<double> h24 =   cos(theta()) * nzpos / err2;
//   LASGlobalData<double> h25 =   cos(theta()) * nzpos  * nzpos  / err2;

//   LASGlobalData<double> h30 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) * nzpos  / err2;
//   LASGlobalData<double> h31 =   sin(theta()) * cos(theta()) * nzpos  / err2;
//   LASGlobalData<double> h32 =   sin(theta()) * cos(theta()) * nzpos  * nzpos  / err2;
//   LASGlobalData<double> h33 =   sin(theta()) * sin(theta()) * nzpos  * nzpos  / err2;
//   LASGlobalData<double> h34 =   sin(theta()) * nzpos / err2;
//   LASGlobalData<double> h35 =   sin(theta()) * nzpos  * nzpos  / err2;

//   LASGlobalData<double> h40 = LASGlobalData<double>(0.) - sin(theta()) / r0() / err2;
//   LASGlobalData<double> h41 =   cos(theta()) / r0() / err2;
//   LASGlobalData<double> h42 =   cos(theta()) * nzpos  / r0() / err2;
//   LASGlobalData<double> h43 =   sin(theta()) * nzpos  / r0() / err2;
//   LASGlobalData<double> h44 =   LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> h45 =   nzpos  / err2;

//   LASGlobalData<double> h50 = LASGlobalData<double>(0.) - sin(theta()) * nzpos  / r0() / err2;
//   LASGlobalData<double> h51 =   cos(theta()) * nzpos  / r0() / err2;
//   LASGlobalData<double> h52 =   cos(theta()) * nzpos  * nzpos  / r0() / err2;
//   LASGlobalData<double> h53 =   sin(theta()) * nzpos  * nzpos  / r0() / err2;
//   LASGlobalData<double> h54 =   nzpos  / err2;
//   LASGlobalData<double> h55 =   nzpos  * nzpos  / err2;


// //   LASGlobalData<double> haa =   zpos() * zpos() / err2;
// //   LASGlobalData<double> hab =   zpos() / err2;
// //   LASGlobalData<double> hba =   zpos() / err2;
// //   LASGlobalData<double> hbb =   LASGlobalData<double>(1.) / err2;


// //   LASGlobalData<double> h00 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) / err2;
// //   LASGlobalData<double> h01 =   sin(theta()) * cos(theta()) / err2;
// //   LASGlobalData<double> h02 =   sin(theta()) * cos(theta()) * zpos() / err2;
// //   LASGlobalData<double> h03 =   sin(theta()) * sin(theta()) * zpos() / err2;
// //   LASGlobalData<double> h04 =   sin(theta()) / err2;
// //   LASGlobalData<double> h05 =   sin(theta()) * zpos() / err2;

// //   LASGlobalData<double> h10 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) / err2;
// //   LASGlobalData<double> h11 =   cos(theta()) * cos(theta()) / err2;
// //   LASGlobalData<double> h12 =   cos(theta()) * cos(theta()) * zpos() / err2;
// //   LASGlobalData<double> h13 =   cos(theta()) * sin(theta()) * zpos() / err2;
// //   LASGlobalData<double> h14 =   cos(theta()) / err2;
// //   LASGlobalData<double> h15 =   cos(theta()) * zpos() / err2;

// //   LASGlobalData<double> h20 = LASGlobalData<double>(0.) - cos(theta()) * sin(theta()) * zpos() / err2;
// //   LASGlobalData<double> h21 =   cos(theta()) * cos(theta()) * zpos() / err2;
// //   LASGlobalData<double> h22 =   cos(theta()) * cos(theta()) * zpos() * zpos() / err2;
// //   LASGlobalData<double> h23 =   cos(theta()) * sin(theta()) * zpos() * zpos() / err2;
// //   LASGlobalData<double> h24 =   cos(theta()) * zpos()/ err2;
// //   LASGlobalData<double> h25 =   cos(theta()) * zpos() * zpos() / err2;

// //   LASGlobalData<double> h30 = LASGlobalData<double>(0.) - sin(theta()) * sin(theta()) * zpos() / err2;
// //   LASGlobalData<double> h31 =   sin(theta()) * cos(theta()) * zpos() / err2;
// //   LASGlobalData<double> h32 =   sin(theta()) * cos(theta()) * zpos() * zpos() / err2;
// //   LASGlobalData<double> h33 =   sin(theta()) * sin(theta()) * zpos() * zpos() / err2;
// //   LASGlobalData<double> h34 =   sin(theta()) * zpos()/ err2;
// //   LASGlobalData<double> h35 =   sin(theta()) * zpos() * zpos() / err2;

// //   LASGlobalData<double> h40 = LASGlobalData<double>(0.) - sin(theta()) / r0() / err2;
// //   LASGlobalData<double> h41 =   cos(theta()) / r0() / err2;
// //   LASGlobalData<double> h42 =   cos(theta()) * zpos() / r0() / err2;
// //   LASGlobalData<double> h43 =   sin(theta()) * zpos() / r0() / err2;
// //   LASGlobalData<double> h44 =   LASGlobalData<double>(1.) / err2;
// //   LASGlobalData<double> h45 =   zpos() / err2;

// //   LASGlobalData<double> h50 = LASGlobalData<double>(0.) - sin(theta()) * zpos() / r0() / err2;
// //   LASGlobalData<double> h51 =   cos(theta()) * zpos() / r0() / err2;
// //   LASGlobalData<double> h52 =   cos(theta()) * zpos() * zpos() / r0() / err2;
// //   LASGlobalData<double> h53 =   sin(theta()) * zpos() * zpos() / r0() / err2;
// //   LASGlobalData<double> h54 =   zpos() / err2;
// //   LASGlobalData<double> h55 =   zpos() * zpos() / err2;


//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 6;  ++dn){
      
//       if(err2.GetEntry(2,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(2,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(2,-1,bn,dn);
	
// 	h(npb+0,npb+0)   += h00.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+1)   += h01.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+2)   += h02.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+3)   += h03.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+4)   += h04.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+5)   += h05.GetEntry(2,-1,bn,dn);
// 	h(npb+0,2*bn)    += h05.GetEntry(2,-1,bn,dn);
// 	h(npb+0,2*bn+1)  += h04.GetEntry(2,-1,bn,dn);
	
// 	h(npb+1,npb+0)   += h10.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+1)   += h11.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+2)   += h12.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+3)   += h13.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+4)   += h14.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+5)   += h15.GetEntry(2,-1,bn,dn);
// 	h(npb+1,2*bn)    += h15.GetEntry(2,-1,bn,dn);
// 	h(npb+1,2*bn+1)  += h14.GetEntry(2,-1,bn,dn);
	
// 	h(npb+2,npb+0)   += h20.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+1)   += h21.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+2)   += h22.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+3)   += h23.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+4)   += h24.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+5)   += h25.GetEntry(2,-1,bn,dn);
// 	h(npb+2,2*bn)    += h25.GetEntry(2,-1,bn,dn);
// 	h(npb+2,2*bn+1)  += h24.GetEntry(2,-1,bn,dn);
	
// 	h(npb+3,npb+0)   += h30.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+1)   += h31.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+2)   += h32.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+3)   += h33.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+4)   += h34.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+5)   += h35.GetEntry(2,-1,bn,dn);
// 	h(npb+3,2*bn)    += h35.GetEntry(2,-1,bn,dn);
// 	h(npb+3,2*bn+1)  += h34.GetEntry(2,-1,bn,dn);
	
// 	h(npb+4,npb+0)   += h40.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+1)   += h41.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+2)   += h42.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+3)   += h43.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+4)   += h44.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+5)   += h45.GetEntry(2,-1,bn,dn);
// 	h(npb+4,2*bn)    += h45.GetEntry(2,-1,bn,dn);
// 	h(npb+4,2*bn+1)  += h44.GetEntry(2,-1,bn,dn);
	
// 	h(npb+5,npb+0)   += h50.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+1)   += h51.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+2)   += h52.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+3)   += h53.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+4)   += h54.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+5)   += h55.GetEntry(2,-1,bn,dn);
// 	h(npb+5,2*bn)    += h55.GetEntry(2,-1,bn,dn);
// 	h(npb+5,2*bn+1)  += h54.GetEntry(2,-1,bn,dn);
	
	
// 	h(2*bn,npb+0)   += h50.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+1)   += h51.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+2)   += h52.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+3)   += h53.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+4)   += h54.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+5)   += h55.GetEntry(2,-1,bn,dn);
	
// 	h(2*bn+1,npb+0) += h40.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+1) += h41.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+2) += h42.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+3) += h43.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+4) += h44.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+5) += h45.GetEntry(2,-1,bn,dn);
//       }
	
//       if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
//       }
//     }
//   }  

  
//   h.Invert();

//   double cp[22][2][8][6];
//   double up[22][22];
//   double ps[22];

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ib = 0; ib < 2; ++ib){
//       for(int bn = 0; bn < 8;  ++bn){
// 	for(int dn = 0; dn < 6;  ++dn){
// 	  cp[ip][ib][bn][dn] = 0.;
// 	  for(int ix = 0; ix < 22; ++ix){
// 	    cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
// 	    }
// 	}
//       }
//     }
//   }

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ix = 0; ix < 22; ++ix){
//       up[ip][ix] = 0.;
//       for(int ib = 0; ib < 2; ++ib){
// 	for(int bn = 0; bn < 8;  ++bn){
// 	  for(int dn = 0; dn < 6;  ++dn){
// 	    if( ib == 0){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(2,-1,bn,dn);
// 	    }
// 	    if( ib == 1 ){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
 
//   for(int bn = 0; bn < 8; ++bn){
//     results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
//     results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
//   }

//   results.er_tib_dx = sqrt(up[npb+0][npb+0]);
//   results.er_tib_dy = sqrt(up[npb+1][npb+1]);
//   results.er_tib_rx = sqrt(up[npb+2][npb+2]);
//   results.er_tib_ry = sqrt(up[npb+3][npb+3]);
//   results.er_tib_rz = sqrt(up[npb+4][npb+4]);
//   results.er_tib_tz = sqrt(up[npb+5][npb+5]);
//   results.tib_dxdy  = up[npb+0][npb+1]/results.er_tib_dx/results.er_tib_dy;
//   results.tib_dxrx  = up[npb+0][npb+2]/results.er_tib_dx/results.er_tib_rx;
//   results.tib_dxry  = up[npb+0][npb+3]/results.er_tib_dx/results.er_tib_ry;
//   results.tib_dxrz  = up[npb+0][npb+4]/results.er_tib_dx/results.er_tib_rz;
//   results.tib_dxtz  = up[npb+0][npb+5]/results.er_tib_dx/results.er_tib_tz;
//   results.tib_dyrx  = up[npb+1][npb+2]/results.er_tib_dy/results.er_tib_rx;
//   results.tib_dyry  = up[npb+1][npb+3]/results.er_tib_dy/results.er_tib_ry;
//   results.tib_dyrz  = up[npb+1][npb+4]/results.er_tib_dy/results.er_tib_rz;
//   results.tib_dytz  = up[npb+1][npb+5]/results.er_tib_dy/results.er_tib_tz;
//   results.tib_rxry  = up[npb+2][npb+3]/results.er_tib_rx/results.er_tib_ry;
//   results.tib_rxrz  = up[npb+2][npb+4]/results.er_tib_rx/results.er_tib_rz;
//   results.tib_rxtz  = up[npb+2][npb+5]/results.er_tib_rx/results.er_tib_tz;
//   results.tib_ryrz  = up[npb+3][npb+4]/results.er_tib_ry/results.er_tib_rz;
//   results.tib_rytz  = up[npb+3][npb+5]/results.er_tib_ry/results.er_tib_tz;
//   results.tib_rztz  = up[npb+4][npb+5]/results.er_tib_rz/results.er_tib_tz;

//   for(int ip = 0; ip < 22; ++ip){
//     ps[ip] = 0.;
//     for(int ix = 0; ix < 22; ++ix){
//       ps[ip] += h(ip,ix)*bs[ix];
//     }
//   }

//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = ps[2*bn];
//     results.beam_b[bn] = ps[2*bn+1];
//   }

//   results.tib_dx = ps[npb+0];
//   results.tib_dy = ps[npb+1];
//   results.tib_rx = ps[npb+2];
//   results.tib_ry = ps[npb+3];
//   results.tib_rz = ps[npb+4];
//   results.tib_tz = ps[npb+5];

// }


// void TIB_chi2(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask,  AtPar& par) //calculation of chi2/ndf at TIB&TOB
// {
//   LASGlobalData<double> recpos = TIB_spot_rec(par);
//   LASGlobalData<double> resid = recpos - dif;
//   LASGlobalData<double> res2 = pow(resid/(err2),2);
//   LASGlobalData<double> dif2 = pow(dif/(err2),2);

//   int nchi_tib,nchi_tob,nchi,nchib[8];
 
//   nchi = 0;
//   par.AT_chi = 0.;
//   par.AT_chi0 = 0.;

//   nchi_tib = 0;
//   nchi_tob = 0;
//   par.tib_chi = 0.;
//   par.tib_chi0 = 0.;
//   par.tob_chi = 0.;
//   par.tob_chi0 = 0.;
//   for(int bnn = 0; bnn < 8;  ++bnn){
//     nchib[bnn] = 0; 
//     par.b_chi[bnn] = 0.;
//     par.b_chi0[bnn] = 0;
//   } 

//   LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);
//   do{
//     if( tob_loop.GetEntry(mask) != 1 ) continue;
//     int bn = tob_loop.get_beam();
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tob = nchi_tob + 1;
//     par.AT_chi += tob_loop.GetEntry(res2);
//     par.AT_chi0 += tob_loop.GetEntry(dif2);
//     par.tob_chi += tob_loop.GetEntry(res2);
//     par.tob_chi0 += tob_loop.GetEntry(dif2);
//     par.b_chi[bn] += tob_loop.GetEntry(res2);
//     par.b_chi0[bn] += tob_loop.GetEntry(dif2);
//   }while(tob_loop.next());

//   if(nchi_tob > 0){
//     par.tob_chi = sqrt(par.tob_chi/nchi_tob);
//     par.tob_chi0 = sqrt(par.tob_chi0/nchi_tob);
//   }
//   //cout<<"  nchi_tob: "<<nchi_tob<<"  tob_chi: "<<par.tob_chi<<"  tob_chi0: "<<par.tob_chi0<<endl;  
 
//   LASGlobalDataLoop tib_loop(LASGlobalDataLoop::TIB);
//   do{
//     if( tib_loop.GetEntry(mask) != 1 ) continue;
//     int bn = tib_loop.get_beam();
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tib = nchi_tib + 1;
//     par.AT_chi += tib_loop.GetEntry(res2);
//     par.AT_chi0 += tib_loop.GetEntry(dif2);
//     par.tib_chi += tib_loop.GetEntry(res2);
//     par.tib_chi0 += tib_loop.GetEntry(dif2);
//     par.b_chi[bn] += tib_loop.GetEntry(res2);
//     par.b_chi0[bn] += tib_loop.GetEntry(dif2);
//   }while(tib_loop.next());

//   if(nchi_tib > 0){
//     par.tib_chi = sqrt(par.tib_chi/nchi_tib);
//     par.tib_chi0 = sqrt(par.tib_chi0/nchi_tib);
//   }
//   //cout<<"  nchi_tib: "<<nchi_tib<<"  tib_chi: "<<par.tib_chi<<"  tib_chi0: "<<par.tib_chi0<<endl;  
 
//   for(int bnn = 0; bnn < 8;  ++bnn){
//     par.b_chi[bnn] = sqrt(par.b_chi[bnn]/(nchib[bnn]-2));
//     par.b_chi0[bnn] = sqrt(par.b_chi0[bnn]/nchib[bnn]);
//   }

//   for(int bnn = 0; bnn < 8;  ++bnn){
//     if(nchib[bnn] > 2){
//       par.b_chi[bnn] = sqrt(par.b_chi[bnn]/(nchib[bnn]-2));
//       par.b_chi0[bnn] = sqrt(par.b_chi0[bnn]/nchib[bnn]);
//     }
//     //cout<<" bn= "<<bnn<<"  nchi_b: "<<nchib[bnn]<<"  b_chi: "<<par.b_chi[bnn]<<" b_chi0: "<<par.b_chi0[bnn]<<endl;  
//   }
  
//   if(nchi > 22){
//     par.AT_chi = sqrt(par.AT_chi/(nchi-22));
//     par.AT_chi0 = sqrt(par.AT_chi0/nchi);
//   }
//   //cout<<"  nchi_AT: "<<nchi<<"  AT_chi: "<<par.AT_chi<<"  AT_chi0: "<<par.AT_chi0<<endl;  
  
// }

// void TECP_chi2(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, const LASGlobalData<int>& mask,  AtPar& par) //calculation of chi2/ndf at TOB&TEC+
// {
//   LASGlobalData<double> recpos = TECP_spot_rec(par);
//   LASGlobalData<double> resid = recpos - dif;
//   LASGlobalData<double> res2 = pow(resid/err2,2);
//   LASGlobalData<double> dif2 = pow(dif/err2,2);

//   int nchi_tob,nchi_tecp,nchi,nchib[8];
 
//   nchi = 0;
//   par.AT_chi = 0.;
//   par.AT_chi0 = 0.;

//   nchi_tob = 0;
//   nchi_tecp = 0;
//   par.tob_chi = 0.;
//   par.tob_chi0 = 0.;
//   par.tecp_chi = 0.;
//   par.tecp_chi0 = 0.;
//   for(int bnn = 0; bnn < 8;  ++bnn){
//     nchib[bnn] = 0; 
//     par.b_chi[bnn] = 0.;
//     par.b_chi0[bnn] = 0;
//   } 
 
//   LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);
//   do{
//     if( tob_loop.GetEntry(mask) != 1 ) continue;
//     int bn = tob_loop.get_beam();
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tob = nchi_tob + 1;
//     par.AT_chi += tob_loop.GetEntry(res2);
//     par.AT_chi0 += tob_loop.GetEntry(dif2);
//     par.tob_chi += tob_loop.GetEntry(res2);
//     par.tob_chi0 += tob_loop.GetEntry(dif2);
//     par.b_chi[bn] += tob_loop.GetEntry(res2);
//     par.b_chi0[bn] += tob_loop.GetEntry(dif2);
//   }while(tob_loop.next());

//   if(nchi_tob > 0){
//     par.tob_chi = sqrt(par.tob_chi/nchi_tob);
//     par.tob_chi0 = sqrt(par.tob_chi0/nchi_tob);
//   }
//   //cout<<"  nchi_tob: "<<nchi_tob<<"  tob_chi: "<<par.tob_chi<<"  tob_chi0: "<<par.tob_chi0<<endl;  
 

//   LASGlobalDataLoop tecp_loop(LASGlobalDataLoop::TEC_PLUS_AT);
//   do{
//     int bn = tecp_loop.get_beam();
//     //int dn = tecp_loop.get_zpos();
//     //cout<<"  bn= "<<bn<<"  dn= "<<dn<<"  res: "<<tecp_loop.GetEntry(resid)<<"  dif: "<<tecp_loop.GetEntry(dif)<<endl;  
//     if( tecp_loop.GetEntry(mask) != 1 ) continue;
//       nchi = nchi + 1;
//       nchib[bn] = nchib[bn] + 1;
//       nchi_tecp = nchi_tecp + 1;
//       par.AT_chi += tecp_loop.GetEntry(res2);
//       par.AT_chi0 += tecp_loop.GetEntry(dif2);
//       par.tecp_chi += tecp_loop.GetEntry(res2);
//       par.tecp_chi0 += tecp_loop.GetEntry(dif2);
//       par.b_chi[bn] += tecp_loop.GetEntry(res2);
//       par.b_chi0[bn] += tecp_loop.GetEntry(dif2);
//       //}
//   }while(tecp_loop.next());
 
  
//   if(nchi_tecp > 0){
//     par.tecp_chi = sqrt(par.tecp_chi/nchi_tecp);
//     par.tecp_chi0 = sqrt(par.tecp_chi0/nchi_tecp);
//   }
//   //cout<<"  nchi_tecp: "<<nchi_tecp<<"  tecp_chi: "<<par.tecp_chi<<" tecp_chi0: "<<par.tecp_chi0<<endl;  
//   for(int bnn = 0; bnn < 8;  ++bnn){
//     if(nchib[bnn] > 2){
//       par.b_chi[bnn] = sqrt(par.b_chi[bnn]/(nchib[bnn]-2));
//       par.b_chi0[bnn] = sqrt(par.b_chi0[bnn]/nchib[bnn]);
//     }
//     //cout<<" bn= "<<bnn<<"  nchi_b: "<<nchib[bnn]<<"  b_chi: "<<par.b_chi[bnn]<<" b_chi0: "<<par.b_chi0[bnn]<<endl;  
//   }
  
//   if(nchi > 19){
//     par.AT_chi = sqrt(par.AT_chi/(nchi-19));
//     par.AT_chi0 = sqrt(par.AT_chi0/nchi);
//   }
//   //cout<<"  nchi_AT: "<<nchi<<"  AT_chi: "<<par.AT_chi<<"  AT_chi0: "<<par.AT_chi0<<endl;  
  
// }

// void TECM_chi2(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, const LASGlobalData<int>& mask,  AtPar& par) //calculation of chi2/ndf at TOB&TEC-
// {
//   LASGlobalData<double> recpos = TECM_spot_rec(par);
//   LASGlobalData<double> resid = recpos - dif;
//   LASGlobalData<double> res2 = pow(resid/err2,2);
//   LASGlobalData<double> dif2 = pow(dif/err2,2);

//   int nchi_tob,nchi_tecm,nchi,nchib[8];
 
//   nchi = 0;
//   par.AT_chi = 0.;
//   par.AT_chi0 = 0.;

//   nchi_tob = 0;
//   nchi_tecm = 0;
//   par.tob_chi = 0.;
//   par.tob_chi0 = 0.;
//   par.tecm_chi = 0.;
//   par.tecm_chi0 = 0.;
//   for(int bnn = 0; bnn < 8;  ++bnn){
//     nchib[bnn] = 0; 
//     par.b_chi[bnn] = 0.;
//     par.b_chi0[bnn] = 0;
//   } 
 
//   LASGlobalDataLoop tob_loop(LASGlobalDataLoop::TOB);
//   do{
//     if( tob_loop.GetEntry(mask) != 1 ) continue;
//     int bn = tob_loop.get_beam();
//     nchi = nchi + 1;
//     nchib[bn] = nchib[bn] + 1;
//     nchi_tob = nchi_tob + 1;
//     par.AT_chi += tob_loop.GetEntry(res2);
//     par.AT_chi0 += tob_loop.GetEntry(dif2);
//     par.tob_chi += tob_loop.GetEntry(res2);
//     par.tob_chi0 += tob_loop.GetEntry(dif2);
//     par.b_chi[bn] += tob_loop.GetEntry(res2);
//     par.b_chi0[bn] += tob_loop.GetEntry(dif2);
//   }while(tob_loop.next());

//   if(nchi_tob > 0){
//     par.tob_chi = sqrt(par.tob_chi/nchi_tob);
//     par.tob_chi0 = sqrt(par.tob_chi0/nchi_tob);
//   }
//   //cout<<"  nchi_tob: "<<nchi_tob<<"  tob_chi: "<<par.tob_chi<<"  tob_chi0: "<<par.tob_chi0<<endl;  
 

//   LASGlobalDataLoop tecm_loop(LASGlobalDataLoop::TEC_MINUS_AT);
//   do{
//     int bn = tecm_loop.get_beam();
//     //int dn = tecm_loop.get_zpos();
//     //cout<<"  bn= "<<bn<<"  dn= "<<dn<<"  res: "<<tecm_loop.GetEntry(resid)<<"  dif: "<<tecm_loop.GetEntry(dif)<<endl;  
//     if( tecm_loop.GetEntry(mask) != 1 ) continue;
//       nchi = nchi + 1;
//       nchib[bn] = nchib[bn] + 1;
//       nchi_tecm = nchi_tecm + 1;
//       par.AT_chi += tecm_loop.GetEntry(res2);
//       par.AT_chi0 += tecm_loop.GetEntry(dif2);
//       par.tecm_chi += tecm_loop.GetEntry(res2);
//       par.tecm_chi0 += tecm_loop.GetEntry(dif2);
//       par.b_chi[bn] += tecm_loop.GetEntry(res2);
//       par.b_chi0[bn] += tecm_loop.GetEntry(dif2);
//       //}
//   }while(tecm_loop.next());
 
  
//   if(nchi_tecm > 0){
//     par.tecm_chi = sqrt(par.tecm_chi/nchi_tecm);
//     par.tecm_chi0 = sqrt(par.tecm_chi0/nchi_tecm);
//   }
//   //cout<<"  nchi_tecm: "<<nchi_tecm<<"  tecm_chi: "<<par.tecm_chi<<" tecm_chi0: "<<par.tecm_chi0<<endl;  
//   for(int bnn = 0; bnn < 8;  ++bnn){
//     if(nchib[bnn] > 2){
//       par.b_chi[bnn] = sqrt(par.b_chi[bnn]/(nchib[bnn]-2));
//       par.b_chi0[bnn] = sqrt(par.b_chi0[bnn]/nchib[bnn]);
//     }
//     //cout<<" bn= "<<bnn<<"  nchi_b: "<<nchib[bnn]<<"  b_chi: "<<par.b_chi[bnn]<<" b_chi0: "<<par.b_chi0[bnn]<<endl;  
//   }
  
//   if(nchi > 19){
//     par.AT_chi = sqrt(par.AT_chi/(nchi-19));
//     par.AT_chi0 = sqrt(par.AT_chi0/nchi);
//   }
//   //cout<<"  nchi_AT: "<<nchi<<"  AT_chi: "<<par.AT_chi<<"  AT_chi0: "<<par.AT_chi0<<endl;  
  
// }

// // reconstruction of TIB wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
// void AT_TIB2TOB_6par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results)
// {

//   // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

//   double bs[22];
//   double db[22][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TIB, TOB); k=0..7 (beam number); l=0..5 (module number)

//   int npb = 16;// number of beam parameters
 
//   LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
//   LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
//   LASGlobalData<double> par2 = dif * cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par3 = dif * sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par4 = dif / err2;
//   LASGlobalData<double> par5 = dif * zpos()/ err2;


//   LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
//   LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
//   LASGlobalData<double> dev2 = cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev3 = sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> dev5 = zpos()/ err2;

//   for(int ip = 0; ip < 22;  ++ip){ 
//     bs[ip] = 0.;
//   }
//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 6;  ++dn){
//       for(int ip = 0; ip < 22;  ++ip){ 
// 	db[ip][0][bn][dn] = 0.;
// 	db[ip][1][bn][dn] = 0.;
//       }
//       if(err2.GetEntry(2,-1,bn,dn) < 1000. && mask.GetEntry(2, -1, bn, dn)){
	    
// 	bs[2*bn] += par5.GetEntry(2,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(2,-1,bn,dn);
	
// 	bs[npb+0] += par0.GetEntry(2,-1,bn,dn);
// 	bs[npb+1] += par1.GetEntry(2,-1,bn,dn);
// 	bs[npb+2] += par2.GetEntry(2,-1,bn,dn);
// 	bs[npb+3] += par3.GetEntry(2,-1,bn,dn);
// 	bs[npb+4] += par4.GetEntry(2,-1,bn,dn);
// 	bs[npb+5] += par5.GetEntry(2,-1,bn,dn);
        
// 	db[2*bn][0][bn][dn] = dev5.GetEntry(2,-1,bn,dn);
// 	db[2*bn+1][0][bn][dn] = dev4.GetEntry(2,-1,bn,dn);
	
// 	db[npb+0][0][bn][dn] = dev0.GetEntry(2,-1,bn,dn);
// 	db[npb+1][0][bn][dn] = dev1.GetEntry(2,-1,bn,dn);
// 	db[npb+2][0][bn][dn] = dev2.GetEntry(2,-1,bn,dn);
// 	db[npb+3][0][bn][dn] = dev3.GetEntry(2,-1,bn,dn);
// 	db[npb+4][0][bn][dn] = dev4.GetEntry(2,-1,bn,dn);
// 	db[npb+5][0][bn][dn] = dev5.GetEntry(2,-1,bn,dn);
//       }
//       if(err2.GetEntry(3,-1,bn,dn) < 1000. && mask.GetEntry(3, -1, bn, dn)){
	    
// 	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
// 	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
// 	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);
//       }
//     }
//   }
   
//   TMatrixT<double> h(22,22);
  
//   for(int ix = 0; ix < 22; ++ix){
//     for(int iy = 0; iy < 22; ++iy){
//       h(ix,iy) = 0.;
//     }
//   }
  
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


//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 6;  ++dn){
      
//       if(err2.GetEntry(2,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(2,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(2,-1,bn,dn);
	
// 	h(npb+0,npb+0)   += h00.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+1)   += h01.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+2)   += h02.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+3)   += h03.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+4)   += h04.GetEntry(2,-1,bn,dn);
// 	h(npb+0,npb+5)   += h05.GetEntry(2,-1,bn,dn);
// 	h(npb+0,2*bn)    += h05.GetEntry(2,-1,bn,dn);
// 	h(npb+0,2*bn+1)  += h04.GetEntry(2,-1,bn,dn);
	
// 	h(npb+1,npb+0)   += h10.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+1)   += h11.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+2)   += h12.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+3)   += h13.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+4)   += h14.GetEntry(2,-1,bn,dn);
// 	h(npb+1,npb+5)   += h15.GetEntry(2,-1,bn,dn);
// 	h(npb+1,2*bn)    += h15.GetEntry(2,-1,bn,dn);
// 	h(npb+1,2*bn+1)  += h14.GetEntry(2,-1,bn,dn);
	
// 	h(npb+2,npb+0)   += h20.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+1)   += h21.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+2)   += h22.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+3)   += h23.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+4)   += h24.GetEntry(2,-1,bn,dn);
// 	h(npb+2,npb+5)   += h25.GetEntry(2,-1,bn,dn);
// 	h(npb+2,2*bn)    += h25.GetEntry(2,-1,bn,dn);
// 	h(npb+2,2*bn+1)  += h24.GetEntry(2,-1,bn,dn);
	
// 	h(npb+3,npb+0)   += h30.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+1)   += h31.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+2)   += h32.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+3)   += h33.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+4)   += h34.GetEntry(2,-1,bn,dn);
// 	h(npb+3,npb+5)   += h35.GetEntry(2,-1,bn,dn);
// 	h(npb+3,2*bn)    += h35.GetEntry(2,-1,bn,dn);
// 	h(npb+3,2*bn+1)  += h34.GetEntry(2,-1,bn,dn);
	
// 	h(npb+4,npb+0)   += h40.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+1)   += h41.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+2)   += h42.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+3)   += h43.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+4)   += h44.GetEntry(2,-1,bn,dn);
// 	h(npb+4,npb+5)   += h45.GetEntry(2,-1,bn,dn);
// 	h(npb+4,2*bn)    += h45.GetEntry(2,-1,bn,dn);
// 	h(npb+4,2*bn+1)  += h44.GetEntry(2,-1,bn,dn);
	
// 	h(npb+5,npb+0)   += h50.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+1)   += h51.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+2)   += h52.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+3)   += h53.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+4)   += h54.GetEntry(2,-1,bn,dn);
// 	h(npb+5,npb+5)   += h55.GetEntry(2,-1,bn,dn);
// 	h(npb+5,2*bn)    += h55.GetEntry(2,-1,bn,dn);
// 	h(npb+5,2*bn+1)  += h54.GetEntry(2,-1,bn,dn);
	
	
// 	h(2*bn,npb+0)   += h50.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+1)   += h51.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+2)   += h52.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+3)   += h53.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+4)   += h54.GetEntry(2,-1,bn,dn);
// 	h(2*bn,npb+5)   += h55.GetEntry(2,-1,bn,dn);
	
// 	h(2*bn+1,npb+0) += h40.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+1) += h41.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+2) += h42.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+3) += h43.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+4) += h44.GetEntry(2,-1,bn,dn);
// 	h(2*bn+1,npb+5) += h45.GetEntry(2,-1,bn,dn);
//       }
	
//       if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
//       }
//     }
//   }  

  
//   h.Invert();

//   double cp[22][2][8][6];
//   double up[22][22];
//   double ps[22];

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ib = 0; ib < 2; ++ib){
//       for(int bn = 0; bn < 8;  ++bn){
// 	for(int dn = 0; dn < 6;  ++dn){
// 	  cp[ip][ib][bn][dn] = 0.;
// 	  for(int ix = 0; ix < 22; ++ix){
// 	    cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
// 	    }
// 	}
//       }
//     }
//   }

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ix = 0; ix < 22; ++ix){
//       up[ip][ix] = 0.;
//       for(int ib = 0; ib < 2; ++ib){
// 	for(int bn = 0; bn < 8;  ++bn){
// 	  for(int dn = 0; dn < 6;  ++dn){
// 	    if( ib == 0){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(2,-1,bn,dn);
// 	    }
// 	    if( ib == 1 ){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
 
//   for(int bn = 0; bn < 8; ++bn){
//     results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
//     results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
//   }

//   results.er_tib_dx = sqrt(up[npb+0][npb+0]);
//   results.er_tib_dy = sqrt(up[npb+1][npb+1]);
//   results.er_tib_rx = sqrt(up[npb+2][npb+2]);
//   results.er_tib_ry = sqrt(up[npb+3][npb+3]);
//   results.er_tib_rz = sqrt(up[npb+4][npb+4]);
//   results.er_tib_tz = sqrt(up[npb+5][npb+5]);
//   results.tib_dxdy  = up[npb+0][npb+1]/results.er_tib_dx/results.er_tib_dy;
//   results.tib_dxrx  = up[npb+0][npb+2]/results.er_tib_dx/results.er_tib_rx;
//   results.tib_dxry  = up[npb+0][npb+3]/results.er_tib_dx/results.er_tib_ry;
//   results.tib_dxrz  = up[npb+0][npb+4]/results.er_tib_dx/results.er_tib_rz;
//   results.tib_dxtz  = up[npb+0][npb+5]/results.er_tib_dx/results.er_tib_tz;
//   results.tib_dyrx  = up[npb+1][npb+2]/results.er_tib_dy/results.er_tib_rx;
//   results.tib_dyry  = up[npb+1][npb+3]/results.er_tib_dy/results.er_tib_ry;
//   results.tib_dyrz  = up[npb+1][npb+4]/results.er_tib_dy/results.er_tib_rz;
//   results.tib_dytz  = up[npb+1][npb+5]/results.er_tib_dy/results.er_tib_tz;
//   results.tib_rxry  = up[npb+2][npb+3]/results.er_tib_rx/results.er_tib_ry;
//   results.tib_rxrz  = up[npb+2][npb+4]/results.er_tib_rx/results.er_tib_rz;
//   results.tib_rxtz  = up[npb+2][npb+5]/results.er_tib_rx/results.er_tib_tz;
//   results.tib_ryrz  = up[npb+3][npb+4]/results.er_tib_ry/results.er_tib_rz;
//   results.tib_rytz  = up[npb+3][npb+5]/results.er_tib_ry/results.er_tib_tz;
//   results.tib_rztz  = up[npb+4][npb+5]/results.er_tib_rz/results.er_tib_tz;

//   for(int ip = 0; ip < 22; ++ip){
//     ps[ip] = 0.;
//     for(int ix = 0; ix < 22; ++ix){
//       ps[ip] += h(ip,ix)*bs[ix];
//     }
//   }

//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = ps[2*bn];
//     results.beam_b[bn] = ps[2*bn+1];
//   }

//   results.tib_dx = ps[npb+0];
//   results.tib_dy = ps[npb+1];
//   results.tib_rx = ps[npb+2];
//   results.tib_ry = ps[npb+3];
//   results.tib_rz = ps[npb+4];
//   results.tib_tz = ps[npb+5];

// }

// // reconstruction of TEC+ wrt TOB 3 alignment parameters: dx, dy, rz and 16 beam parameters: (a&b)*8
// void AT_TECP2TOB_3par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,const LASGlobalData<int>& mask, AtPar& results)
// {

//   // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

//   double bs[19];
//   double db[19][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TEC+, TOB); k=0..7 (beam number); l=0..5 (module number)


//   int npb = 16;// number of beam parameters
//   int mdn = 5;// number of TEC discs
 
//   LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
//   LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
//   LASGlobalData<double> par2 = dif * cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par3 = dif * sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par4 = dif / err2;
//   LASGlobalData<double> par5 = dif * zpos()/ err2;


//   LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
//   LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
//   LASGlobalData<double> dev2 = cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev3 = sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> dev5 = zpos()/ err2;

//   for(int ip = 0; ip < 19;  ++ip){ 
//     bs[ip] = 0.;
//   }
//   for(int bn = 0; bn < 8;  ++bn){
//     //for(int dn = 0; dn < mdn;  ++dn){
//     for(int dn = 0; dn < 6;  ++dn){
//       for(int ip = 0; ip < 19;  ++ip){ 
// 	db[ip][0][bn][dn] = 0.;
// 	db[ip][1][bn][dn] = 0.;
//       }
//       if(dn < mdn && err2.GetEntry(0,-1,bn,dn) < 1000. && mask.GetEntry(0, -1, bn, dn)){
// 	bs[2*bn] += par5.GetEntry(0,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(0,-1,bn,dn);
	
// 	bs[npb+0] += par0.GetEntry(0,-1,bn,dn);
// 	bs[npb+1] += par1.GetEntry(0,-1,bn,dn);
// 	bs[npb+2] += par4.GetEntry(0,-1,bn,dn);
        
// 	db[2*bn][0][bn][dn] = dev5.GetEntry(0,-1,bn,dn);
// 	db[2*bn+1][0][bn][dn] = dev4.GetEntry(0,-1,bn,dn);
	
// 	db[npb+0][0][bn][dn] = dev0.GetEntry(0,-1,bn,dn);
// 	db[npb+1][0][bn][dn] = dev1.GetEntry(0,-1,bn,dn);
// 	db[npb+2][0][bn][dn] = dev4.GetEntry(0,-1,bn,dn);
//       }
//       if(err2.GetEntry(3,-1,bn,dn) < 1000. && mask.GetEntry(3, -1, bn, dn)){
	    
// 	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
// 	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
// 	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);
//       }
//     }
//   }

// //   std::cout << "bs(0) =  " << bs[0] << std::endl;
// //   std::cout << "bs(1) =  " << bs[1] << std::endl;
// //   std::cout << "bs(2) =  " << bs[2] << std::endl;
// //   std::cout << "bs(3) =  " << bs[3] << std::endl;
// //   std::cout << "bs(4) =  " << bs[4] << std::endl;
// //   std::cout << "bs(5) =  " << bs[5] << std::endl;
// //   std::cout << "bs(6) =  " << bs[6] << std::endl;
// //   std::cout << "bs(7) =  " << bs[7] << std::endl;
// //   std::cout << "bs(8) =  " << bs[8] << std::endl;
// //   std::cout << "bs(9) =  " << bs[9] << std::endl;
// //   std::cout << "bs(10) =  " << bs[10] << std::endl;
// //   std::cout << "bs(11) =  " << bs[11] << std::endl;
// //   std::cout << "bs(12) =  " << bs[12] << std::endl;
// //   std::cout << "bs(13) =  " << bs[13] << std::endl;
// //   std::cout << "bs(14) =  " << bs[14] << std::endl;
// //   std::cout << "bs(15) =  " << bs[15] << std::endl;
// //   std::cout << "bs(16) =  " << bs[16] << std::endl;
// //   std::cout << "bs(17) =  " << bs[17] << std::endl;

   
//   TMatrixT<double> h(19,19);
//   TMatrixT<double> h1(19,19);
//   TMatrixT<double> h2(19,19);
  
//   for(int ix = 0; ix < 19; ++ix){
//     for(int iy = 0; iy < 19; ++iy){
//       h(ix,iy) = 0.;
//     }
//   }
  
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


//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 6;  ++dn){

//       if(dn < mdn && err2.GetEntry(0,-1,bn,dn) < 1000. && mask.GetEntry(0, -1, bn, dn)){      
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(0,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(0,-1,bn,dn);
	
// 	h(npb+0,npb+0)   += h00.GetEntry(0,-1,bn,dn);
// 	h(npb+0,npb+1)   += h01.GetEntry(0,-1,bn,dn);
// 	h(npb+0,npb+2)   += h04.GetEntry(0,-1,bn,dn);
// 	h(npb+0,2*bn)    += h05.GetEntry(0,-1,bn,dn);
// 	h(npb+0,2*bn+1)  += h04.GetEntry(0,-1,bn,dn);
	
// 	h(npb+1,npb+0)   += h10.GetEntry(0,-1,bn,dn);
// 	h(npb+1,npb+1)   += h11.GetEntry(0,-1,bn,dn);
// 	h(npb+1,npb+2)   += h14.GetEntry(0,-1,bn,dn);
// 	h(npb+1,2*bn)    += h15.GetEntry(0,-1,bn,dn);
// 	h(npb+1,2*bn+1)  += h14.GetEntry(0,-1,bn,dn);
	
// 	h(npb+2,npb+0)   += h40.GetEntry(0,-1,bn,dn);
// 	h(npb+2,npb+1)   += h41.GetEntry(0,-1,bn,dn);
// 	h(npb+2,npb+2)   += h44.GetEntry(0,-1,bn,dn);
// 	h(npb+2,2*bn)    += h45.GetEntry(0,-1,bn,dn);
// 	h(npb+2,2*bn+1)  += h44.GetEntry(0,-1,bn,dn);
	
// 	h(2*bn,npb+0)   += h50.GetEntry(0,-1,bn,dn);
// 	h(2*bn,npb+1)   += h51.GetEntry(0,-1,bn,dn);
// 	h(2*bn,npb+2)   += h54.GetEntry(0,-1,bn,dn);
	
// 	h(2*bn+1,npb+0) += h40.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,npb+1) += h41.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,npb+2) += h44.GetEntry(0,-1,bn,dn);
//       }
	
//       if(err2.GetEntry(3,-1,bn,dn) < 1000. && mask.GetEntry(3, -1, bn, dn)){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
//       }
//     }
//   }  

  
//   for(int ix = 0; ix < 19; ++ix){
//     for(int iy = 0; iy < 19; ++iy){
//       h1(ix,iy) = h(ix,iy);
//     }
//   }
  
// //   std::cout << "h(0,0) =  " << h(0,0) << std::endl;
// //   std::cout << "h(1,1) =  " << h(1,1) << std::endl;
// //   std::cout << "h(2,2) =  " << h(2,2) << std::endl;
// //   std::cout << "h(3,3) =  " << h(3,3) << std::endl;
// //   std::cout << "h(4,4) =  " << h(4,4) << std::endl;
// //   std::cout << "h(5,5) =  " << h(5,5) << std::endl;
// //   std::cout << "h(6,6) =  " << h(6,6) << std::endl;
// //   std::cout << "h(7,7) =  " << h(7,7) << std::endl;
// //   std::cout << "h(8,8) =  " << h(8,8) << std::endl;
// //   std::cout << "h(9,9) =  " << h(9,9) << std::endl;
// //   std::cout << "h(10,10) =  " << h(10,10) << std::endl;
// //   std::cout << "h(11,11) =  " << h(11,11) << std::endl;
// //   std::cout << "h(12,12) =  " << h(12,12) << std::endl;
// //   std::cout << "h(13,13) =  " << h(13,13) << std::endl;
// //   std::cout << "h(14,14) =  " << h(14,14) << std::endl;
// //   std::cout << "h(15,15) =  " << h(15,15) << std::endl;
// //   std::cout << "h(16,16) =  " << h(16,16) << std::endl;
// //   std::cout << "h(17,17) =  " << h(17,17) << std::endl;
// //   std::cout << "h(18,18) =  " << h(18,18) << std::endl;
  
// //   std::cout << "h0,1) =  " << h(0,1) << std::endl;
// //   std::cout << "h(1,2) =  " << h(1,2) << std::endl;
// //   std::cout << "h(2,3) =  " << h(2,3) << std::endl;
// //   std::cout << "h(3,4) =  " << h(3,4) << std::endl;
// //   std::cout << "h(4,5) =  " << h(4,5) << std::endl;
// //   std::cout << "h(5,6) =  " << h(5,6) << std::endl;
// //   std::cout << "h(6,7) =  " << h(6,7) << std::endl;
// //   std::cout << "h(7,8) =  " << h(7,8) << std::endl;
// //   std::cout << "h(8,9) =  " << h(8,9) << std::endl;
// //   std::cout << "h(9,8) =  " << h(9,8) << std::endl;
// //   std::cout << "h(10,11) =  " << h(10,11) << std::endl;
// //   std::cout << "h(11,12) =  " << h(11,12) << std::endl;
// //   std::cout << "h(12,13) =  " << h(12,13) << std::endl;
// //   std::cout << "h(13,14) =  " << h(13,14) << std::endl;
// //   std::cout << "h(14,15) =  " << h(14,15) << std::endl;
// //   std::cout << "h(15,16) =  " << h(15,16) << std::endl;
// //   std::cout << "h(16,17) =  " << h(16,17) << std::endl;
// //   std::cout << "h(17,18) =  " << h(17,18) << std::endl;
// //   std::cout << "h(18,17) =  " << h(18,17) << std::endl;

//   h.Invert();

//   for(int ix = 0; ix < 19; ++ix){
//     for(int iy = 0; iy < 19; ++iy){
//       h2(ix,iy) = 0.;
//       for(int ixx = 0; ixx < 19; ++ixx){
// 	h2(ix,iy) = h2(ix,iy) + h(ix,ixx)*h1(ixx,iy);
//       }
//     }
//   }
  
// //   std::cout << "h2(0,0) =  " << h2(0,0) << std::endl;
// //   std::cout << "h2(1,1) =  " << h2(1,1) << std::endl;
// //   std::cout << "h2(2,2) =  " << h2(2,2) << std::endl;
// //   std::cout << "h2(3,3) =  " << h2(3,3) << std::endl;
// //   std::cout << "h2(4,4) =  " << h2(4,4) << std::endl;
// //   std::cout << "h2(5,5) =  " << h2(5,5) << std::endl;
// //   std::cout << "h2(6,6) =  " << h2(6,6) << std::endl;
// //   std::cout << "h2(7,7) =  " << h2(7,7) << std::endl;
// //   std::cout << "h2(8,8) =  " << h2(8,8) << std::endl;
// //   std::cout << "h2(9,9) =  " << h2(9,9) << std::endl;
// //   std::cout << "h2(10,10) =  " << h2(10,10) << std::endl;
// //   std::cout << "h2(11,11) =  " << h2(11,11) << std::endl;
// //   std::cout << "h2(12,12) =  " << h2(12,12) << std::endl;
// //   std::cout << "h2(13,13) =  " << h2(13,13) << std::endl;
// //   std::cout << "h2(14,14) =  " << h2(14,14) << std::endl;
// //   std::cout << "h2(15,15) =  " << h2(15,15) << std::endl;
// //   std::cout << "h2(16,16) =  " << h2(16,16) << std::endl;
// //   std::cout << "h2(17,17) =  " << h2(17,17) << std::endl;
// //   std::cout << "h2(18,18) =  " << h2(18,18) << std::endl;
  
// //   std::cout << "h2(0,1) =  " << h2(0,1) << std::endl;
// //   std::cout << "h2(1,2) =  " << h2(1,2) << std::endl;
// //   std::cout << "h2(2,3) =  " << h2(2,3) << std::endl;
// //   std::cout << "h2(3,4) =  " << h2(3,4) << std::endl;
// //   std::cout << "h2(4,5) =  " << h2(4,5) << std::endl;
// //   std::cout << "h2(5,6) =  " << h2(5,6) << std::endl;
// //   std::cout << "h2(6,7) =  " << h2(6,7) << std::endl;
// //   std::cout << "h2(7,8) =  " << h2(7,8) << std::endl;
// //   std::cout << "h2(8,9) =  " << h2(8,9) << std::endl;
// //   std::cout << "h2(9,8) =  " << h2(9,8) << std::endl;
// //   std::cout << "h2(10,11) =  " << h2(10,11) << std::endl;
// //   std::cout << "h2(11,12) =  " << h2(11,12) << std::endl;
// //   std::cout << "h2(12,13) =  " << h2(12,13) << std::endl;
// //   std::cout << "h2(13,14) =  " << h2(13,14) << std::endl;
// //   std::cout << "h2(14,15) =  " << h2(14,15) << std::endl;
// //   std::cout << "h2(15,16) =  " << h2(15,16) << std::endl;
// //   std::cout << "h2(16,17) =  " << h2(16,17) << std::endl;
// //   std::cout << "h2(17,18) =  " << h2(17,18) << std::endl;
// //   std::cout << "h2(18,17) =  " << h2(18,17) << std::endl;

//   double cp[19][2][8][6];
//   double up[19][19];
//   double ps[19];

//   for(int ip = 0; ip < 19; ++ip){
//     for(int ib = 0; ib < 2; ++ib){
//       for(int bn = 0; bn < 8;  ++bn){
// 	for(int dn = 0; dn < 6;  ++dn){
// 	  cp[ip][ib][bn][dn] = 0.;
// 	  for(int ix = 0; ix < 19; ++ix){
// 	    if(ib == 0 && dn < mdn && err2.GetEntry(0,-1,bn,dn) != 0.){
// 	      cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
// 	    }
// 	    if(ib == 1 && err2.GetEntry(3,-1,bn,dn) != 0.){
// 		cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
// 	    }
// 	  }
// 	}
//       }
//     }
//   }

//   for(int ip = 0; ip < 19; ++ip){
//     for(int ix = 0; ix < 19; ++ix){
//       up[ip][ix] = 0.;
//       for(int ib = 0; ib < 2; ++ib){
// 	for(int bn = 0; bn < 8;  ++bn){
// 	  for(int dn = 0; dn < 6;  ++dn){
// 	    if(ib == 0 && dn < mdn && err2.GetEntry(0,-1,bn,dn) != 0.){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(0,-1,bn,dn);
// 	    }
// 	    if( ib == 1 && err2.GetEntry(3,-1,bn,dn) != 0.){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
 
//   for(int bn = 0; bn < 8; ++bn){
//     results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
//     results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
//   }

//   results.er_tecp_dx = sqrt(up[npb+0][npb+0]);
//   results.er_tecp_dy = sqrt(up[npb+1][npb+1]);
//   results.er_tecp_rz = sqrt(up[npb+2][npb+2]);
//   results.tecp_dxdy  = up[npb+0][npb+1]/results.er_tecp_dx/results.er_tecp_dy;
//   results.tecp_dxrz  = up[npb+0][npb+2]/results.er_tecp_dx/results.er_tecp_rz;

//   for(int ip = 0; ip < 19; ++ip){
//     ps[ip] = 0.;
//     for(int ix = 0; ix < 19; ++ix){
//       ps[ip] += h(ip,ix)*bs[ix];
//     }
//     //std::cout <<"ip= "<<ip<< "  ps =  " << ps[ip] << std::endl;
//   }

//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = ps[2*bn];
//     results.beam_b[bn] = ps[2*bn+1];
//   }

//   results.tecp_dx = ps[npb+0];
//   results.tecp_dy = ps[npb+1];
//   results.tecp_rz = ps[npb+2];

// }
  
// // reconstruction of TEC- wrt TOB 3 alignment parameters: dx, dy, rz and 16 beam parameters: (a&b)*8
// void AT_TECM2TOB_3par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,const LASGlobalData<int>& mask, AtPar& results)
// {

//   // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

//   double bs[19];
//   double db[19][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TEC-, TOB); k=0..7 (beam number); l=0..5 (module number)

//   int npb = 16;// number of beam parameters
//   int mdn = 5;// number of TEC discs
 
//   LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
//   LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
//   LASGlobalData<double> par2 = dif * cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par3 = dif * sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par4 = dif / err2;
//   LASGlobalData<double> par5 = dif * zpos()/ err2;

//   LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
//   LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
//   LASGlobalData<double> dev2 = cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev3 = sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> dev5 = zpos()/ err2;

//   for(int ip = 0; ip < 19;  ++ip){ 
//     bs[ip] = 0.;
//   }
//   for(int bn = 0; bn < 8;  ++bn){
//     //for(int dn = 0; dn < mdn;  ++dn){
//     for(int dn = 0; dn < 6;  ++dn){
//       for(int ip = 0; ip < 19;  ++ip){ 
// 	db[ip][0][bn][dn] = 0.;
// 	db[ip][1][bn][dn] = 0.;
//       }
//       if(dn < mdn)if(err2.GetEntry(1,-1,bn,dn) < 1000. && mask.GetEntry(1, -1, bn, dn)){
// 	bs[2*bn] += par5.GetEntry(1,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(1,-1,bn,dn);
	
// 	bs[npb+0] += par0.GetEntry(1,-1,bn,dn);
// 	bs[npb+1] += par1.GetEntry(1,-1,bn,dn);
// 	bs[npb+2] += par4.GetEntry(1,-1,bn,dn);
        
// 	db[2*bn][0][bn][dn] = dev5.GetEntry(1,-1,bn,dn);
// 	db[2*bn+1][0][bn][dn] = dev4.GetEntry(1,-1,bn,dn);
	
// 	db[npb+0][0][bn][dn] = dev0.GetEntry(1,-1,bn,dn);
// 	db[npb+1][0][bn][dn] = dev1.GetEntry(1,-1,bn,dn);
// 	db[npb+2][0][bn][dn] = dev4.GetEntry(1,-1,bn,dn);
//       }
//       if(err2.GetEntry(3,-1,bn,dn) < 1000. && mask.GetEntry(3, -1, bn, dn)){
	    
// 	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
// 	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
// 	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);
//       }
//     }
//   }

//   TMatrixT<double> h(19,19);
//   TMatrixT<double> h1(19,19);
//   TMatrixT<double> h2(19,19);

//   for(int ix = 0; ix < 19; ++ix){
//     for(int iy = 0; iy < 19; ++iy){
//       h(ix,iy) = 0.;
//     }
//   }
  
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


//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 6;  ++dn){

//       if(dn < mdn && err2.GetEntry(1,-1,bn,dn) < 1000. && mask.GetEntry(1, -1, bn, dn)){      

// 	h(2*bn,2*bn)     += haa.GetEntry(1,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(1,-1,bn,dn);
	
// 	h(npb+0,npb+0)   += h00.GetEntry(1,-1,bn,dn);
// 	h(npb+0,npb+1)   += h01.GetEntry(1,-1,bn,dn);
// 	h(npb+0,npb+2)   += h04.GetEntry(1,-1,bn,dn);
// 	h(npb+0,2*bn)    += h05.GetEntry(1,-1,bn,dn);
// 	h(npb+0,2*bn+1)  += h04.GetEntry(1,-1,bn,dn);
	
// 	h(npb+1,npb+0)   += h10.GetEntry(1,-1,bn,dn);
// 	h(npb+1,npb+1)   += h11.GetEntry(1,-1,bn,dn);
// 	h(npb+1,npb+2)   += h14.GetEntry(1,-1,bn,dn);
// 	h(npb+1,2*bn)    += h15.GetEntry(1,-1,bn,dn);
// 	h(npb+1,2*bn+1)  += h14.GetEntry(1,-1,bn,dn);
	
// 	h(npb+2,npb+0)   += h40.GetEntry(1,-1,bn,dn);
// 	h(npb+2,npb+1)   += h41.GetEntry(1,-1,bn,dn);
// 	h(npb+2,npb+2)   += h44.GetEntry(1,-1,bn,dn);
// 	h(npb+2,2*bn)    += h45.GetEntry(1,-1,bn,dn);
// 	h(npb+2,2*bn+1)  += h44.GetEntry(1,-1,bn,dn);
	
// 	h(2*bn,npb+0)   += h50.GetEntry(1,-1,bn,dn);
// 	h(2*bn,npb+1)   += h51.GetEntry(1,-1,bn,dn);
// 	h(2*bn,npb+2)   += h54.GetEntry(1,-1,bn,dn);
	
// 	h(2*bn+1,npb+0) += h40.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,npb+1) += h41.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,npb+2) += h44.GetEntry(1,-1,bn,dn);
//       }
	
//       if(err2.GetEntry(3,-1,bn,dn) < 1000. && mask.GetEntry(3, -1, bn, dn)){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
//       }
//     }
//   }  

//   for(int ix = 0; ix < 19; ++ix){
//     for(int iy = 0; iy < 19; ++iy){
//       h1(ix,iy) = h(ix,iy);
//     }
//   }

//   h.Invert();

//   for(int ix = 0; ix < 19; ++ix){
//     for(int iy = 0; iy < 19; ++iy){
//       h2(ix,iy) = 0.;
//       for(int ixx = 0; ixx < 19; ++ixx){
// 	h2(ix,iy) = h2(ix,iy) + h(ix,ixx)*h1(ixx,iy);
//       }
//     }
//   }  

//   double cp[19][2][8][6];
//   double up[19][19];
//   double ps[19];

//   for(int ip = 0; ip < 19; ++ip){
//     for(int ib = 0; ib < 2; ++ib){
//       for(int bn = 0; bn < 8;  ++bn){
// 	for(int dn = 0; dn < 6;  ++dn){
// 	  cp[ip][ib][bn][dn] = 0.;
// 	  for(int ix = 0; ix < 19; ++ix){
// 	    if(ib == 0 && dn < mdn && err2.GetEntry(1,-1,bn,dn) != 0.){
// 	      cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
// 	    }
// 	    if(ib == 1 && err2.GetEntry(3,-1,bn,dn) != 0.){
// 	      cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
// 	    }
// 	  }
// 	}
//       }
//     }
//   }

//   for(int ip = 0; ip < 19; ++ip){
//     for(int ix = 0; ix < 19; ++ix){
//       up[ip][ix] = 0.;
//       for(int ib = 0; ib < 2; ++ib){
// 	for(int bn = 0; bn < 8;  ++bn){
// 	  for(int dn = 0; dn < 6;  ++dn){
// 	    if(ib == 0 && dn < mdn && err2.GetEntry(1,-1,bn,dn) != 0.){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(1,-1,bn,dn);
// 	    }
// 	    if( ib == 1 && err2.GetEntry(3,-1,bn,dn) != 0.){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
 
//   for(int bn = 0; bn < 8; ++bn){
//     results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
//     results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
//   }

//   results.er_tecm_dx = sqrt(up[npb+0][npb+0]);
//   results.er_tecm_dy = sqrt(up[npb+1][npb+1]);
//   results.er_tecm_rz = sqrt(up[npb+2][npb+2]);
//   results.tecm_dxdy  = up[npb+0][npb+1]/results.er_tecm_dx/results.er_tecm_dy;
//   results.tecm_dxrz  = up[npb+0][npb+2]/results.er_tecm_dx/results.er_tecm_rz;

//   for(int ip = 0; ip < 19; ++ip){
//     ps[ip] = 0.;
//     for(int ix = 0; ix < 19; ++ix){
//       ps[ip] += h(ip,ix)*bs[ix];
//     }
//   }

//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = ps[2*bn];
//     results.beam_b[bn] = ps[2*bn+1];
//   }

//   results.tecm_dx = ps[npb+0];
//   results.tecm_dy = ps[npb+1];
//   results.tecm_rz = ps[npb+2];
// }

// // reconstruction of TEC+ wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
// void AT_TECP2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results)
// {

//   // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

//   double bs[22];
//   double db[22][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TEC+, TOB); k=0..7 (beam number); l=0..5 (module number)

//   int npb = 16;// number of beam parameters
 
//   LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
//   LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
//   LASGlobalData<double> par2 = dif * cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par3 = dif * sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par4 = dif / err2;
//   LASGlobalData<double> par5 = dif * zpos()/ err2;


//   LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
//   LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
//   LASGlobalData<double> dev2 = cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev3 = sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> dev5 = zpos()/ err2;

//   for(int ip = 0; ip < 22;  ++ip){ 
//     bs[ip] = 0.;
//   }
//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 4;  ++dn){
//       for(int ip = 0; ip < 22;  ++ip){ 
// 	db[ip][0][bn][dn] = 0.;
// 	db[ip][1][bn][dn] = 0.;
//       }
//       if(err2.GetEntry(0,-1,bn,dn) < 1000.){
	    
// 	bs[2*bn] += par5.GetEntry(0,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(0,-1,bn,dn);
	
// 	bs[npb+0] += par0.GetEntry(0,-1,bn,dn);
// 	bs[npb+1] += par1.GetEntry(0,-1,bn,dn);
// 	bs[npb+2] += par2.GetEntry(0,-1,bn,dn);
// 	bs[npb+3] += par3.GetEntry(0,-1,bn,dn);
// 	bs[npb+4] += par4.GetEntry(0,-1,bn,dn);
// 	bs[npb+5] += par5.GetEntry(0,-1,bn,dn);
        
// 	db[2*bn][0][bn][dn] = dev5.GetEntry(0,-1,bn,dn);
// 	db[2*bn+1][0][bn][dn] = dev4.GetEntry(0,-1,bn,dn);
	
// 	db[npb+0][0][bn][dn] = dev0.GetEntry(0,-1,bn,dn);
// 	db[npb+1][0][bn][dn] = dev1.GetEntry(0,-1,bn,dn);
// 	db[npb+2][0][bn][dn] = dev2.GetEntry(0,-1,bn,dn);
// 	db[npb+3][0][bn][dn] = dev3.GetEntry(0,-1,bn,dn);
// 	db[npb+4][0][bn][dn] = dev4.GetEntry(0,-1,bn,dn);
// 	db[npb+5][0][bn][dn] = dev5.GetEntry(0,-1,bn,dn);
//       }
//       if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
// 	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
// 	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
// 	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);
//       }
//     }
//   }
   
//   TMatrixT<double> h(22,22);
  
//   for(int ix = 0; ix < 22; ++ix){
//     for(int iy = 0; iy < 22; ++iy){
//       h(ix,iy) = 0.;
//     }
//   }
  
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


//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 4;  ++dn){
      
//       if(err2.GetEntry(0,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(0,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(0,-1,bn,dn);
	
// 	h(npb+0,npb+0)   += h00.GetEntry(0,-1,bn,dn);
// 	h(npb+0,npb+1)   += h01.GetEntry(0,-1,bn,dn);
// 	h(npb+0,npb+2)   += h02.GetEntry(0,-1,bn,dn);
// 	h(npb+0,npb+3)   += h03.GetEntry(0,-1,bn,dn);
// 	h(npb+0,npb+4)   += h04.GetEntry(0,-1,bn,dn);
// 	h(npb+0,npb+5)   += h05.GetEntry(0,-1,bn,dn);
// 	h(npb+0,2*bn)    += h05.GetEntry(0,-1,bn,dn);
// 	h(npb+0,2*bn+1)  += h04.GetEntry(0,-1,bn,dn);
	
// 	h(npb+1,npb+0)   += h10.GetEntry(0,-1,bn,dn);
// 	h(npb+1,npb+1)   += h11.GetEntry(0,-1,bn,dn);
// 	h(npb+1,npb+2)   += h12.GetEntry(0,-1,bn,dn);
// 	h(npb+1,npb+3)   += h13.GetEntry(0,-1,bn,dn);
// 	h(npb+1,npb+4)   += h14.GetEntry(0,-1,bn,dn);
// 	h(npb+1,npb+5)   += h15.GetEntry(0,-1,bn,dn);
// 	h(npb+1,2*bn)    += h15.GetEntry(0,-1,bn,dn);
// 	h(npb+1,2*bn+1)  += h14.GetEntry(0,-1,bn,dn);
	
// 	h(npb+2,npb+0)   += h20.GetEntry(0,-1,bn,dn);
// 	h(npb+2,npb+1)   += h21.GetEntry(0,-1,bn,dn);
// 	h(npb+2,npb+2)   += h22.GetEntry(0,-1,bn,dn);
// 	h(npb+2,npb+3)   += h23.GetEntry(0,-1,bn,dn);
// 	h(npb+2,npb+4)   += h24.GetEntry(0,-1,bn,dn);
// 	h(npb+2,npb+5)   += h25.GetEntry(0,-1,bn,dn);
// 	h(npb+2,2*bn)    += h25.GetEntry(0,-1,bn,dn);
// 	h(npb+2,2*bn+1)  += h24.GetEntry(0,-1,bn,dn);
	
// 	h(npb+3,npb+0)   += h30.GetEntry(0,-1,bn,dn);
// 	h(npb+3,npb+1)   += h31.GetEntry(0,-1,bn,dn);
// 	h(npb+3,npb+2)   += h32.GetEntry(0,-1,bn,dn);
// 	h(npb+3,npb+3)   += h33.GetEntry(0,-1,bn,dn);
// 	h(npb+3,npb+4)   += h34.GetEntry(0,-1,bn,dn);
// 	h(npb+3,npb+5)   += h35.GetEntry(0,-1,bn,dn);
// 	h(npb+3,2*bn)    += h35.GetEntry(0,-1,bn,dn);
// 	h(npb+3,2*bn+1)  += h34.GetEntry(0,-1,bn,dn);
	
// 	h(npb+4,npb+0)   += h40.GetEntry(0,-1,bn,dn);
// 	h(npb+4,npb+1)   += h41.GetEntry(0,-1,bn,dn);
// 	h(npb+4,npb+2)   += h42.GetEntry(0,-1,bn,dn);
// 	h(npb+4,npb+3)   += h43.GetEntry(0,-1,bn,dn);
// 	h(npb+4,npb+4)   += h44.GetEntry(0,-1,bn,dn);
// 	h(npb+4,npb+5)   += h45.GetEntry(0,-1,bn,dn);
// 	h(npb+4,2*bn)    += h45.GetEntry(0,-1,bn,dn);
// 	h(npb+4,2*bn+1)  += h44.GetEntry(0,-1,bn,dn);
	
// 	h(npb+5,npb+0)   += h50.GetEntry(0,-1,bn,dn);
// 	h(npb+5,npb+1)   += h51.GetEntry(0,-1,bn,dn);
// 	h(npb+5,npb+2)   += h52.GetEntry(0,-1,bn,dn);
// 	h(npb+5,npb+3)   += h53.GetEntry(0,-1,bn,dn);
// 	h(npb+5,npb+4)   += h54.GetEntry(0,-1,bn,dn);
// 	h(npb+5,npb+5)   += h55.GetEntry(0,-1,bn,dn);
// 	h(npb+5,2*bn)    += h55.GetEntry(0,-1,bn,dn);
// 	h(npb+5,2*bn+1)  += h54.GetEntry(0,-1,bn,dn);
	
	
// 	h(2*bn,npb+0)   += h50.GetEntry(0,-1,bn,dn);
// 	h(2*bn,npb+1)   += h51.GetEntry(0,-1,bn,dn);
// 	h(2*bn,npb+2)   += h52.GetEntry(0,-1,bn,dn);
// 	h(2*bn,npb+3)   += h53.GetEntry(0,-1,bn,dn);
// 	h(2*bn,npb+4)   += h54.GetEntry(0,-1,bn,dn);
// 	h(2*bn,npb+5)   += h55.GetEntry(0,-1,bn,dn);
	
// 	h(2*bn+1,npb+0) += h40.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,npb+1) += h41.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,npb+2) += h42.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,npb+3) += h43.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,npb+4) += h44.GetEntry(0,-1,bn,dn);
// 	h(2*bn+1,npb+5) += h45.GetEntry(0,-1,bn,dn);
//       }
	
//       if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
//       }
//     }
//   }  

  
//   h.Invert();

//   double cp[22][2][8][6];
//   double up[22][22];
//   double ps[22];

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ib = 0; ib < 2; ++ib){
//       for(int bn = 0; bn < 8;  ++bn){
// 	for(int dn = 0; dn < 4;  ++dn){
// 	  cp[ip][ib][bn][dn] = 0.;
// 	  for(int ix = 0; ix < 22; ++ix){
// 	    cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
// 	    }
// 	}
//       }
//     }
//   }

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ix = 0; ix < 22; ++ix){
//       up[ip][ix] = 0.;
//       for(int ib = 0; ib < 2; ++ib){
// 	for(int bn = 0; bn < 8;  ++bn){
// 	  for(int dn = 0; dn < 4;  ++dn){
// 	    if( ib == 0){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(0,-1,bn,dn);
// 	    }
// 	    if( ib == 1 ){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
 
//   for(int bn = 0; bn < 8; ++bn){
//     results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
//     results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
//   }

//   results.er_tecp_dx = sqrt(up[npb+0][npb+0]);
//   results.er_tecp_dy = sqrt(up[npb+1][npb+1]);
//   results.er_tecp_rx = sqrt(up[npb+2][npb+2]);
//   results.er_tecp_ry = sqrt(up[npb+3][npb+3]);
//   results.er_tecp_rz = sqrt(up[npb+4][npb+4]);
//   results.er_tecp_tz = sqrt(up[npb+5][npb+5]);
//   results.tecp_dxdy  = up[npb+0][npb+1]/results.er_tecp_dx/results.er_tecp_dy;
//   results.tecp_dxrx  = up[npb+0][npb+2]/results.er_tecp_dx/results.er_tecp_rx;
//   results.tecp_dxry  = up[npb+0][npb+3]/results.er_tecp_dx/results.er_tecp_ry;
//   results.tecp_dxrz  = up[npb+0][npb+4]/results.er_tecp_dx/results.er_tecp_rz;
//   results.tecp_dxtz  = up[npb+0][npb+5]/results.er_tecp_dx/results.er_tecp_tz;
//   results.tecp_dyrx  = up[npb+1][npb+2]/results.er_tecp_dy/results.er_tecp_rx;
//   results.tecp_dyry  = up[npb+1][npb+3]/results.er_tecp_dy/results.er_tecp_ry;
//   results.tecp_dyrz  = up[npb+1][npb+4]/results.er_tecp_dy/results.er_tecp_rz;
//   results.tecp_dytz  = up[npb+1][npb+5]/results.er_tecp_dy/results.er_tecp_tz;
//   results.tecp_rxry  = up[npb+2][npb+3]/results.er_tecp_rx/results.er_tecp_ry;
//   results.tecp_rxrz  = up[npb+2][npb+4]/results.er_tecp_rx/results.er_tecp_rz;
//   results.tecp_rxtz  = up[npb+2][npb+5]/results.er_tecp_rx/results.er_tecp_tz;
//   results.tecp_ryrz  = up[npb+3][npb+4]/results.er_tecp_ry/results.er_tecp_rz;
//   results.tecp_rytz  = up[npb+3][npb+5]/results.er_tecp_ry/results.er_tecp_tz;
//   results.tecp_rztz  = up[npb+4][npb+5]/results.er_tecp_rz/results.er_tecp_tz;

//   for(int ip = 0; ip < 22; ++ip){
//     ps[ip] = 0.;
//     for(int ix = 0; ix < 22; ++ix){
//       ps[ip] += h(ip,ix)*bs[ix];
//     }
//   }

//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = ps[2*bn];
//     results.beam_b[bn] = ps[2*bn+1];
//   }

//   results.tecp_dx = ps[npb+0];
//   results.tecp_dy = ps[npb+1];
//   results.tecp_rx = ps[npb+2];
//   results.tecp_ry = ps[npb+3];
//   results.tecp_rz = ps[npb+4];
//   results.tecp_tz = ps[npb+5];

// }
  

// // reconstruction of TEC- wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
// void AT_TECM2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results)
// {

//   // vector of parameters (ps) is defined by linear equation: ps = h^{-1}*bs, where h - matrix, and bs - vector

//   double bs[22];
//   double db[22][2][8][6]; // matrix of derivatives d(ps_i)/d(mes_j_k_l), i=0..21 (parameter number); j=0,1 (TEC-, TOB); k=0..7 (beam number); l=0..5 (module number)

//   int npb = 16;// number of beam parameters
 
//   LASGlobalData<double> par0 = dif * sin(theta()) * r0() / err2;
//   LASGlobalData<double> par1 = dif * cos(theta()) * r0() / err2;
//   LASGlobalData<double> par2 = dif * cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par3 = dif * sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> par4 = dif / err2;
//   LASGlobalData<double> par5 = dif * zpos()/ err2;


//   LASGlobalData<double> dev0 = sin(theta()) * r0() / err2;
//   LASGlobalData<double> dev1 = cos(theta()) * r0() / err2;
//   LASGlobalData<double> dev2 = cos(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev3 = sin(theta()) * r0() * zpos()/ err2;
//   LASGlobalData<double> dev4 = LASGlobalData<double>(1.) / err2;
//   LASGlobalData<double> dev5 = zpos()/ err2;

//   for(int ip = 0; ip < 22;  ++ip){ 
//     bs[ip] = 0.;
//   }
//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 4;  ++dn){
//       for(int ip = 0; ip < 22;  ++ip){ 
// 	db[ip][0][bn][dn] = 0.;
// 	db[ip][1][bn][dn] = 0.;
//       }
//       if(err2.GetEntry(1,-1,bn,dn) < 1000.){
	    
// 	bs[2*bn] += par5.GetEntry(1,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(1,-1,bn,dn);
	
// 	bs[npb+0] += par0.GetEntry(1,-1,bn,dn);
// 	bs[npb+1] += par1.GetEntry(1,-1,bn,dn);
// 	bs[npb+2] += par2.GetEntry(1,-1,bn,dn);
// 	bs[npb+3] += par3.GetEntry(1,-1,bn,dn);
// 	bs[npb+4] += par4.GetEntry(1,-1,bn,dn);
// 	bs[npb+5] += par5.GetEntry(1,-1,bn,dn);
        
// 	db[2*bn][0][bn][dn] = dev5.GetEntry(1,-1,bn,dn);
// 	db[2*bn+1][0][bn][dn] = dev4.GetEntry(1,-1,bn,dn);
	
// 	db[npb+0][0][bn][dn] = dev0.GetEntry(1,-1,bn,dn);
// 	db[npb+1][0][bn][dn] = dev1.GetEntry(1,-1,bn,dn);
// 	db[npb+2][0][bn][dn] = dev2.GetEntry(1,-1,bn,dn);
// 	db[npb+3][0][bn][dn] = dev3.GetEntry(1,-1,bn,dn);
// 	db[npb+4][0][bn][dn] = dev4.GetEntry(1,-1,bn,dn);
// 	db[npb+5][0][bn][dn] = dev5.GetEntry(1,-1,bn,dn);
//       }
//       if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
// 	bs[2*bn] += par5.GetEntry(3,-1,bn,dn);
// 	bs[2*bn+1] += par4.GetEntry(3,-1,bn,dn);
	
// 	db[2*bn][1][bn][dn] = dev5.GetEntry(3,-1,bn,dn);
// 	db[2*bn+1][1][bn][dn] = dev4.GetEntry(3,-1,bn,dn);
//       }
//     }
//   }
   
//   TMatrixT<double> h(22,22);
  
//   for(int ix = 0; ix < 22; ++ix){
//     for(int iy = 0; iy < 22; ++iy){
//       h(ix,iy) = 0.;
//     }
//   }
  
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


//   for(int bn = 0; bn < 8;  ++bn){
//     for(int dn = 0; dn < 4;  ++dn){
      
//       if(err2.GetEntry(1,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(1,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(1,-1,bn,dn);
	
// 	h(npb+0,npb+0)   += h00.GetEntry(1,-1,bn,dn);
// 	h(npb+0,npb+1)   += h01.GetEntry(1,-1,bn,dn);
// 	h(npb+0,npb+2)   += h02.GetEntry(1,-1,bn,dn);
// 	h(npb+0,npb+3)   += h03.GetEntry(1,-1,bn,dn);
// 	h(npb+0,npb+4)   += h04.GetEntry(1,-1,bn,dn);
// 	h(npb+0,npb+5)   += h05.GetEntry(1,-1,bn,dn);
// 	h(npb+0,2*bn)    += h05.GetEntry(1,-1,bn,dn);
// 	h(npb+0,2*bn+1)  += h04.GetEntry(1,-1,bn,dn);
	
// 	h(npb+1,npb+0)   += h10.GetEntry(1,-1,bn,dn);
// 	h(npb+1,npb+1)   += h11.GetEntry(1,-1,bn,dn);
// 	h(npb+1,npb+2)   += h12.GetEntry(1,-1,bn,dn);
// 	h(npb+1,npb+3)   += h13.GetEntry(1,-1,bn,dn);
// 	h(npb+1,npb+4)   += h14.GetEntry(1,-1,bn,dn);
// 	h(npb+1,npb+5)   += h15.GetEntry(1,-1,bn,dn);
// 	h(npb+1,2*bn)    += h15.GetEntry(1,-1,bn,dn);
// 	h(npb+1,2*bn+1)  += h14.GetEntry(1,-1,bn,dn);
	
// 	h(npb+2,npb+0)   += h20.GetEntry(1,-1,bn,dn);
// 	h(npb+2,npb+1)   += h21.GetEntry(1,-1,bn,dn);
// 	h(npb+2,npb+2)   += h22.GetEntry(1,-1,bn,dn);
// 	h(npb+2,npb+3)   += h23.GetEntry(1,-1,bn,dn);
// 	h(npb+2,npb+4)   += h24.GetEntry(1,-1,bn,dn);
// 	h(npb+2,npb+5)   += h25.GetEntry(1,-1,bn,dn);
// 	h(npb+2,2*bn)    += h25.GetEntry(1,-1,bn,dn);
// 	h(npb+2,2*bn+1)  += h24.GetEntry(1,-1,bn,dn);
	
// 	h(npb+3,npb+0)   += h30.GetEntry(1,-1,bn,dn);
// 	h(npb+3,npb+1)   += h31.GetEntry(1,-1,bn,dn);
// 	h(npb+3,npb+2)   += h32.GetEntry(1,-1,bn,dn);
// 	h(npb+3,npb+3)   += h33.GetEntry(1,-1,bn,dn);
// 	h(npb+3,npb+4)   += h34.GetEntry(1,-1,bn,dn);
// 	h(npb+3,npb+5)   += h35.GetEntry(1,-1,bn,dn);
// 	h(npb+3,2*bn)    += h35.GetEntry(1,-1,bn,dn);
// 	h(npb+3,2*bn+1)  += h34.GetEntry(1,-1,bn,dn);
	
// 	h(npb+4,npb+0)   += h40.GetEntry(1,-1,bn,dn);
// 	h(npb+4,npb+1)   += h41.GetEntry(1,-1,bn,dn);
// 	h(npb+4,npb+2)   += h42.GetEntry(1,-1,bn,dn);
// 	h(npb+4,npb+3)   += h43.GetEntry(1,-1,bn,dn);
// 	h(npb+4,npb+4)   += h44.GetEntry(1,-1,bn,dn);
// 	h(npb+4,npb+5)   += h45.GetEntry(1,-1,bn,dn);
// 	h(npb+4,2*bn)    += h45.GetEntry(1,-1,bn,dn);
// 	h(npb+4,2*bn+1)  += h44.GetEntry(1,-1,bn,dn);
	
// 	h(npb+5,npb+0)   += h50.GetEntry(1,-1,bn,dn);
// 	h(npb+5,npb+1)   += h51.GetEntry(1,-1,bn,dn);
// 	h(npb+5,npb+2)   += h52.GetEntry(1,-1,bn,dn);
// 	h(npb+5,npb+3)   += h53.GetEntry(1,-1,bn,dn);
// 	h(npb+5,npb+4)   += h54.GetEntry(1,-1,bn,dn);
// 	h(npb+5,npb+5)   += h55.GetEntry(1,-1,bn,dn);
// 	h(npb+5,2*bn)    += h55.GetEntry(1,-1,bn,dn);
// 	h(npb+5,2*bn+1)  += h54.GetEntry(1,-1,bn,dn);
	
	
// 	h(2*bn,npb+0)   += h50.GetEntry(1,-1,bn,dn);
// 	h(2*bn,npb+1)   += h51.GetEntry(1,-1,bn,dn);
// 	h(2*bn,npb+2)   += h52.GetEntry(1,-1,bn,dn);
// 	h(2*bn,npb+3)   += h53.GetEntry(1,-1,bn,dn);
// 	h(2*bn,npb+4)   += h54.GetEntry(1,-1,bn,dn);
// 	h(2*bn,npb+5)   += h55.GetEntry(1,-1,bn,dn);
	
// 	h(2*bn+1,npb+0) += h40.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,npb+1) += h41.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,npb+2) += h42.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,npb+3) += h43.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,npb+4) += h44.GetEntry(1,-1,bn,dn);
// 	h(2*bn+1,npb+5) += h45.GetEntry(1,-1,bn,dn);
//       }
	
//       if(err2.GetEntry(3,-1,bn,dn) < 1000.){
	    
// 	h(2*bn,2*bn)     += haa.GetEntry(3,-1,bn,dn);
// 	h(2*bn,2*bn+1)   += hab.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn)   += hba.GetEntry(3,-1,bn,dn);
// 	h(2*bn+1,2*bn+1) += hbb.GetEntry(3,-1,bn,dn);	
	
//       }
//     }
//   }  

  
//   h.Invert();

//   double cp[22][2][8][6];
//   double up[22][22];
//   double ps[22];

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ib = 0; ib < 2; ++ib){
//       for(int bn = 0; bn < 8;  ++bn){
// 	for(int dn = 0; dn < 4;  ++dn){
// 	  cp[ip][ib][bn][dn] = 0.;
// 	  for(int ix = 0; ix < 22; ++ix){
// 	    cp[ip][ib][bn][dn] = cp[ip][ib][bn][dn] + h[ip][ix]*db[ix][ib][bn][dn];
// 	    }
// 	}
//       }
//     }
//   }

//   for(int ip = 0; ip < 22; ++ip){
//     for(int ix = 0; ix < 22; ++ix){
//       up[ip][ix] = 0.;
//       for(int ib = 0; ib < 2; ++ib){
// 	for(int bn = 0; bn < 8;  ++bn){
// 	  for(int dn = 0; dn < 4;  ++dn){
// 	    if( ib == 0){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][0][bn][dn]*cp[ix][0][bn][dn]*err2.GetEntry(1,-1,bn,dn);
// 	    }
// 	    if( ib == 1 ){
// 	      up[ip][ix] = up[ip][ix] + cp[ip][1][bn][dn]*cp[ix][1][bn][dn]*err2.GetEntry(3,-1,bn,dn);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
 
//   for(int bn = 0; bn < 8; ++bn){
//     results.er_beam_a[bn] = sqrt(up[2*bn][2*bn]);
//     results.er_beam_b[bn] = sqrt(up[2*bn+1][2*bn+1]);
//   }

//   results.er_tecm_dx = sqrt(up[npb+0][npb+0]);
//   results.er_tecm_dy = sqrt(up[npb+1][npb+1]);
//   results.er_tecm_rx = sqrt(up[npb+2][npb+2]);
//   results.er_tecm_ry = sqrt(up[npb+3][npb+3]);
//   results.er_tecm_rz = sqrt(up[npb+4][npb+4]);
//   results.er_tecm_tz = sqrt(up[npb+5][npb+5]);
//   results.tecm_dxdy  = up[npb+0][npb+1]/results.er_tecm_dx/results.er_tecm_dy;
//   results.tecm_dxrx  = up[npb+0][npb+2]/results.er_tecm_dx/results.er_tecm_rx;
//   results.tecm_dxry  = up[npb+0][npb+3]/results.er_tecm_dx/results.er_tecm_ry;
//   results.tecm_dxrz  = up[npb+0][npb+4]/results.er_tecm_dx/results.er_tecm_rz;
//   results.tecm_dxtz  = up[npb+0][npb+5]/results.er_tecm_dx/results.er_tecm_tz;
//   results.tecm_dyrx  = up[npb+1][npb+2]/results.er_tecm_dy/results.er_tecm_rx;
//   results.tecm_dyry  = up[npb+1][npb+3]/results.er_tecm_dy/results.er_tecm_ry;
//   results.tecm_dyrz  = up[npb+1][npb+4]/results.er_tecm_dy/results.er_tecm_rz;
//   results.tecm_dytz  = up[npb+1][npb+5]/results.er_tecm_dy/results.er_tecm_tz;
//   results.tecm_rxry  = up[npb+2][npb+3]/results.er_tecm_rx/results.er_tecm_ry;
//   results.tecm_rxrz  = up[npb+2][npb+4]/results.er_tecm_rx/results.er_tecm_rz;
//   results.tecm_rxtz  = up[npb+2][npb+5]/results.er_tecm_rx/results.er_tecm_tz;
//   results.tecm_ryrz  = up[npb+3][npb+4]/results.er_tecm_ry/results.er_tecm_rz;
//   results.tecm_rytz  = up[npb+3][npb+5]/results.er_tecm_ry/results.er_tecm_tz;
//   results.tecm_rztz  = up[npb+4][npb+5]/results.er_tecm_rz/results.er_tecm_tz;

//   for(int ip = 0; ip < 22; ++ip){
//     ps[ip] = 0.;
//     for(int ix = 0; ix < 22; ++ix){
//       ps[ip] += h(ip,ix)*bs[ix];
//     }
//   }

//   for(int bn = 0; bn < 8; ++bn){
//     results.beam_a[bn] = ps[2*bn];
//     results.beam_b[bn] = ps[2*bn+1];
//   }

//   results.tecm_dx = ps[npb+0];
//   results.tecm_dy = ps[npb+1];
//   results.tecm_rx = ps[npb+2];
//   results.tecm_ry = ps[npb+3];
//   results.tecm_rz = ps[npb+4];
//   results.tecm_tz = ps[npb+5];

// }



// ////////////////////////////////////////////////////////////////////////////////////////////
// // Draw Alignment Tube Parameters

// void DrawAtPar_old(std::vector<AtPar>& parlist, Avec& time)
// {
//   Avec tib_dx;
//   Avec tib_dy;
//   Avec tib_rx;
//   Avec tib_ry;
//   Avec tib_rz;
//   Avec tib_tz;
//   Avec tib_chi;
//   Avec tib_chi0;
//   Avec tob_chi;
//   Avec tob_chi0;
//   Avec2D beam_a;
//   Avec2D beam_b;

//   for(unsigned int i = 0; i < parlist.size(); i++){
//     tib_dx.push_back(parlist[i].tib_dx);
//     tib_dy.push_back(parlist[i].tib_dy);
//     tib_rx.push_back(parlist[i].tib_rx);
//     tib_ry.push_back(parlist[i].tib_ry);
//     tib_rz.push_back(parlist[i].tib_rz);
//     tib_tz.push_back(parlist[i].tib_tz);
//     tib_chi.push_back(parlist[i].tib_chi);
//     tib_chi0.push_back(parlist[i].tib_chi0);
//     tob_chi.push_back(parlist[i].tob_chi);
//     tob_chi0.push_back(parlist[i].tob_chi0);
//     beam_a.push_back(parlist[i].beam_a);
//     beam_b.push_back(parlist[i].beam_b);
//   }
//   beam_a = vtrans(beam_a);
//   beam_b = vtrans(beam_b);

//   // Move z origin to beamsplitter position
//   //beam_b += beam_a * LAS::z_at_bs;
//   //tib_dx -= tib_ry * LAS::z_at_bs;
//   //tib_dy += tib_rx * LAS::z_at_bs;
//   //tib_rz += tib_tz * LAS::z_at_bs;

//   // Rescale to units for display
//   tib_dx *= 1e3;
//   tib_dy *= 1e3;
//   tib_rx *= 1e6;
//   tib_ry *= 1e6;
//   tib_rz *= 1e6;
//   tib_tz *= 1e9;
//   beam_a *= 1e6;

//   new TCanvas("TIB_AT","AT"); 
//   TMultiGraph* gr = avec_draw(time, (tib_dx | tib_dy | tib_rx| tib_ry | tib_rz | tib_tz ),"TIB vs TOB (16 + 6p fit)","Date","#Deltax, #Deltay [#mum], Rx,Ry,Rz [#murad], Tz[#murad/m] ","AP");
//   //TMultiGraph* gr = avec_draw(time, (dx | dy ),"TEC_AT wrt TOB Stability over Run (11 par)","Block Nr","#Deltax, #Deltay [#mum], Rz [#murad] ","AP");
//   //gr->GetYaxis()->SetRangeUser(-5,5);
//   gr->GetXaxis()->SetNdivisions(505, kTRUE);
//   gr->GetXaxis()->SetTimeDisplay(1);
//   gr->GetXaxis()->SetTimeFormat("%H %d/%m");
//   //avec_draw(time, dy,"dy vs time","time","dy","AP");
//   //avec_draw(time, dphi,"dx vs time","time","dphi","AP");
//   //avec_draw(dx, dy,"dx vs dy","dx","dy","AP");
 
//   std::vector<std::string> Legend;
//   Legend.push_back("#Delta x");
//   Legend.push_back("#Delta y");
//   Legend.push_back("Rx");
//   Legend.push_back("Ry");
//   Legend.push_back("Rz");
//   Legend.push_back("Tz");
//   // Legend.push_back("#Deltax");
//   // Legend.push_back("#Deltay");
//   // Legend.push_back("Rx");
//   // Legend.push_back("Ry");
//   // Legend.push_back("Rz");
//   AddLegend(gr, Legend);
//   //TLine *L1=new TLine(time.front(),1, time.back(), 1);
//   //TLine *L2=new TLine(time.front(), -1, time.back(), -1);
//   //L1->Draw();
//   //L2->Draw();

//   //avec_draw(xvals, tib_tz,"TIB_Tz","Block Nr","Tz [micron]","AP");
//   //avec_draw(tib_tz,"TIB_Tz","Block Nr","Tz [micron]","AP");

//   //std::cout << "tib_chi:\n" << tib_chi << std::endl;
//   TCanvas* cv_beampar = new TCanvas("beampar","beampar");
//   cv_beampar->Divide(1,2);
//   cv_beampar->cd(1);
//   avec_draw(time, beam_a, "Beam Slope", "Date", "[#murad/m]", "AP");
//   cv_beampar->cd(2);
//   avec_draw(time, beam_b, "Beam Offset", "Date", "", "AP");

//   new TCanvas("tibchi","tibchi");
//   avec_draw(time, (tib_chi | tib_chi0 | tob_chi | tob_chi0), "TIB chi2", "Date", "", "AP");
// }

void DrawATPar(std::vector<AtPar>& parlist, Avec& time, Bool_t quick_draw, const std::string& output_file, LAS::beam_group beam_group)
{
  Avec tecp_dx;
  Avec tecp_dy;
  Avec tecp_rx;
  Avec tecp_ry;
  Avec tecp_rz;
  Avec tecp_tz;
  Avec tecm_dx;
  Avec tecm_dy;
  Avec tecm_rx;
  Avec tecm_ry;
  Avec tecm_rz;
  Avec tecm_tz;
  Avec tib_dx;
  Avec tib_dy;
  Avec tib_rx;
  Avec tib_ry;
  Avec tib_rz;
  Avec tib_tz;

  Avec tecp_er_dx;
  Avec tecp_er_dy;
  Avec tecp_er_rx;
  Avec tecp_er_ry;
  Avec tecp_er_rz;
  Avec tecp_er_tz;
  Avec tecm_er_dx;
  Avec tecm_er_dy;
  Avec tecm_er_rx;
  Avec tecm_er_ry;
  Avec tecm_er_rz;
  Avec tecm_er_tz;
  Avec tib_er_dx;
  Avec tib_er_dy;
  Avec tib_er_rx;
  Avec tib_er_ry;
  Avec tib_er_rz;
  Avec tib_er_tz;

  Avec tecp_chi;
  Avec tecp_chi0;
  Avec tecm_chi;
  Avec tecm_chi0;
  Avec tib_chi;
  Avec tib_chi0;
  Avec tob_chi;
  Avec tob_chi0;
  Avec AT_chi;
  Avec AT_chi0;

  Avec beam_slp[8];
  Avec beam_int[8];

  Avec beam_er_slp[8];
  Avec beam_er_int[8];

  Avec b_chi[8];
  Avec b_chi0[8];

  for(unsigned int i = 0; i < parlist.size(); i++){
    tecp_dx.push_back(parlist[i].tecp_dx * 1e3);
    tecp_dy.push_back(parlist[i].tecp_dy * 1e3);
    tecp_rx.push_back(parlist[i].tecp_rx * 1e6);
    tecp_ry.push_back(parlist[i].tecp_ry * 1e6);
    tecp_rz.push_back(parlist[i].tecp_rz * 1e6);
    tecp_tz.push_back(parlist[i].tecp_tz * 1e9);

    tecp_er_dx.push_back(parlist[i].er_tecp_dx * 1e3);
    tecp_er_dy.push_back(parlist[i].er_tecp_dy * 1e3);
    tecp_er_rx.push_back(parlist[i].er_tecp_rx * 1e6);
    tecp_er_ry.push_back(parlist[i].er_tecp_ry * 1e6);
    tecp_er_rz.push_back(parlist[i].er_tecp_rz * 1e6);
    tecp_er_tz.push_back(parlist[i].er_tecp_tz * 1e9);

    tecp_chi.push_back(parlist[i].chi.tecp);
    tecp_chi0.push_back(parlist[i].chi0.tecp);
  }

  for(unsigned int i = 0; i < parlist.size(); i++){
    tecm_dx.push_back(parlist[i].tecm_dx * 1e3);
    tecm_dy.push_back(parlist[i].tecm_dy * 1e3);
    tecm_rx.push_back(parlist[i].tecm_rx * 1e6);
    tecm_ry.push_back(parlist[i].tecm_ry * 1e6);
    tecm_rz.push_back(parlist[i].tecm_rz * 1e6);
    tecm_tz.push_back(parlist[i].tecm_tz * 1e9);

    tecm_er_dx.push_back(parlist[i].er_tecm_dx * 1e3);
    tecm_er_dy.push_back(parlist[i].er_tecm_dy * 1e3);
    tecm_er_rx.push_back(parlist[i].er_tecm_rx * 1e6);
    tecm_er_ry.push_back(parlist[i].er_tecm_ry * 1e6);
    tecm_er_rz.push_back(parlist[i].er_tecm_rz * 1e6);
    tecm_er_tz.push_back(parlist[i].er_tecm_tz * 1e9);

    tecm_chi.push_back(parlist[i].chi.tecm);
    tecm_chi0.push_back(parlist[i].chi0.tecm);
  }

  for(unsigned int i = 0; i < parlist.size(); i++){
    //if (parlist[i].tib_chi > 10.0) continue;
    tib_dx.push_back(parlist[i].tib_dx * 1e3);
    tib_dy.push_back(parlist[i].tib_dy * 1e3);
    tib_rx.push_back(parlist[i].tib_rx * 1e6);
    tib_ry.push_back(parlist[i].tib_ry * 1e6);
    tib_rz.push_back(parlist[i].tib_rz * 1e6);
    tib_tz.push_back(parlist[i].tib_tz * 1e9);

    tib_er_dx.push_back(parlist[i].er_tib_dx * 1e3);
    tib_er_dy.push_back(parlist[i].er_tib_dy * 1e3);
    tib_er_rx.push_back(parlist[i].er_tib_rx * 1e6);
    tib_er_ry.push_back(parlist[i].er_tib_ry * 1e6);
    tib_er_rz.push_back(parlist[i].er_tib_rz * 1e6);
    tib_er_tz.push_back(parlist[i].er_tib_tz * 1e9);

    tib_chi.push_back(parlist[i].chi.tib);
    tib_chi0.push_back(parlist[i].chi0.tib);
    tob_chi.push_back(parlist[i].chi.tob);
    tob_chi0.push_back(parlist[i].chi0.tob);
    AT_chi.push_back(parlist[i].chi.AT);
    AT_chi0.push_back(parlist[i].chi0.AT);
    for(int nb = 0; nb < 8; nb++){
      beam_slp[nb].push_back( parlist[i].beam_a[nb] * LAS::r_at_tec * 1000000);
      beam_int[nb].push_back( parlist[i].beam_b[nb] * LAS::r_at_tec * 1000 );

      beam_er_slp[nb].push_back( parlist[i].er_beam_a[nb] * LAS::r_at_tec * 1000000);
      beam_er_int[nb].push_back( parlist[i].er_beam_b[nb] * LAS::r_at_tec * 1000 );

      b_chi[nb].push_back(parlist[i].chi.beam[nb]);
      b_chi0[nb].push_back(parlist[i].chi0.beam[nb]);
    }
  }

  gROOT->SetStyle("Plain");

  TFile of(output_file.c_str(), "RECREATE");
  TDirectory* details_at = of.mkdir("Details_AT");

  TCanvas* cv = 0;
  TMultiGraph* gr = 0;

  // Canvas with TIB parameters
  if(beam_group == LAS::TIB || beam_group == LAS::AT || beam_group == LAS::ALL){
    cv = new TCanvas("TIB_AT","TIB_AT");
    cv->SetGridy(); 
    gr = avec_draw(time, (tib_dx | tib_dy | tib_rx| tib_ry | tib_rz | tib_tz ), (tib_er_dx | tib_er_dy | tib_er_rx| tib_er_ry | tib_er_rz | tib_er_tz ),"TIB vs TOB","Date","#Deltax, #Deltay [#mum], Rx,Ry,Rz [#murad], Tz[#murad/m] ","AP", 1, "mgr_tib_par");
    gr->GetXaxis()->SetNdivisions(505, kTRUE);
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
    std::vector<std::string> Legend;
    Legend.push_back("#Delta x");
    Legend.push_back("#Delta y");
    Legend.push_back("Rx");
    Legend.push_back("Ry");
    Legend.push_back("Rz");
    Legend.push_back("Tz");
    AddLegend(gr, Legend, 1, 0.9, 0.6, 1.0, 1.0);
    cv->Write();
    gr->Write("mgr_tib_par");
  }

  // Canvas with TEC+ parameters
  if(beam_group == LAS::TEC_PLUS_AT || beam_group == LAS::TEC_AT || beam_group == LAS::AT || beam_group == LAS::ALL){
    cv = new TCanvas("TECP_AT","TECP_AT");
    cv->SetGridy(); 
    TMultiGraph* grp = avec_draw(time, (tecp_dx | tecp_dy | tecp_rx | tecp_ry | tecp_rz | tecp_tz ), (tecp_er_dx | tecp_er_dy | tecp_er_rx | tecp_er_ry | tecp_er_rz | tecp_er_tz ),"TEC+ vs TOB","Date","#Deltax, #Deltay [#mum],Rz [#murad] ","AP");
    grp->GetXaxis()->SetNdivisions(505, kTRUE);
    grp->GetXaxis()->SetTimeDisplay(1);
    grp->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
    std::vector<std::string> Legendp;
    Legendp.push_back("#Delta x");
    Legendp.push_back("#Delta y");
    Legendp.push_back("Rx");
    Legendp.push_back("Ry");
    Legendp.push_back("Rz");
    Legendp.push_back("Tz");
    AddLegend(gr, Legendp, 1, 0.9, 0.6, 1.0, 1.0);
    cv->Write();
  }

  // Canvas with TEC- parameters
  if(beam_group == LAS::TEC_MINUS_AT || beam_group == LAS::TEC_AT || beam_group == LAS::AT || beam_group == LAS::ALL){
    cv = new TCanvas("TECM_AT","TECM_AT");
    cv->SetGridy(); 
    TMultiGraph* grm = avec_draw(time, (tecm_dx | tecm_dy | tecm_rx | tecm_ry | tecm_rz | tecm_tz ), (tecm_er_dx | tecm_er_dy | tecm_er_rx | tecm_er_ry | tecm_er_rz | tecm_er_tz ),"TEC+ vs TOB","Date","#Deltax, #Deltay [#mum],Rz [#murad] ","AP");
    grm->GetXaxis()->SetNdivisions(505, kTRUE);
    grm->GetXaxis()->SetTimeDisplay(1);
    grm->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
    std::vector<std::string> Legendm;
    Legendm.push_back("#Delta x");
    Legendm.push_back("#Delta y");
    Legendm.push_back("Rx");
    Legendm.push_back("Ry");
    Legendm.push_back("Rz");
    Legendm.push_back("Tz");
    AddLegend(gr, Legendm, 1, 0.9, 0.6, 1.0, 1.0);
    cv->Write();
  }

  ////////////////////////////////////////
  // AT Details
  /////////////////

  // If desired, switch temporarily to batch mode, to suppress graphical output
  Bool_t batch_mode = gROOT->IsBatch();
  gROOT->SetBatch(quick_draw);

  of.cd("Details_AT");

  // Common Details that are always there (TOB, Beams, Total Chi2)
  cv = new TCanvas("CH_TOB","CH_TOB"); 
  TMultiGraph* grco = avec_draw(time, (tob_chi | tob_chi0),"TIB vs TOB (16 + 6p fit)","Date","#chi2/ndf","AP");
  grco->GetXaxis()->SetNdivisions(505, kTRUE);
  grco->GetXaxis()->SetTimeDisplay(1);
  grco->GetXaxis()->SetTimeFormat("%d/%m");

  grco->GetYaxis()->SetRangeUser(0.1, 1e4);

  std::vector<std::string> Legendco;
  Legendco.push_back("TOB after fit");
  Legendco.push_back("TOB before fit");
  AddLegend(grco, Legendco);
  cv->SetLogy();
  cv->SetGridy();
  cv->Write();

  cv = new TCanvas("CH_AT","CH_AT"); 
  TMultiGraph* grca = avec_draw(time, (AT_chi | AT_chi0),"TIB vs TOB (16 + 6p fit)","Date","#chi2/ndf","AP");
  grca->GetXaxis()->SetNdivisions(505, kTRUE);
  grca->GetXaxis()->SetTimeDisplay(1);
  grca->GetXaxis()->SetTimeFormat("%d/%m");

  grca->GetYaxis()->SetRangeUser(0.1, 1e4);
 
  std::vector<std::string> Legendca;
  Legendca.push_back("AT after fit");
  Legendca.push_back("AT before fit");
  AddLegend(grca, Legendca);
  cv->SetLogy();
  cv->SetGridy();
  cv->Write();

  cv = new TCanvas("BEAM_A","BEAM_A"); 
  TMultiGraph* gra = avec_draw(time, (beam_slp[0] | beam_slp[1] | beam_slp[2] | beam_slp[3] | beam_slp[4] | beam_slp[5] | beam_slp[6] | beam_slp[7]), (beam_er_slp[0] | beam_er_slp[1] | beam_er_slp[2] | beam_er_slp[3] | beam_er_slp[4] | beam_er_slp[5] | beam_er_slp[6] | beam_er_slp[7]),"TIB vs TOB (16 + 6p fit)","Date","Beam slopes[#murad] ","AP");
  gra->GetXaxis()->SetNdivisions(505, kTRUE);
  gra->GetXaxis()->SetTimeDisplay(1);
  gra->GetXaxis()->SetTimeFormat("%d/%m");
  cv->Write();

  cv = new TCanvas("BEAM_B","BEAM_B"); 
  TMultiGraph* grb = avec_draw(time, (beam_int[0] | beam_int[1] | beam_int[2] | beam_int[3] | beam_int[4] | beam_int[5] | beam_int[6] | beam_int[7]), (beam_er_int[0] | beam_er_int[1] | beam_er_int[2] | beam_er_int[3] | beam_er_int[4] | beam_er_int[5] | beam_er_int[6] | beam_er_int[7]),"TIB vs TOB (16 + 6p fit)","Date","Beam intersepts[#mum] ","AP");
  grb->GetXaxis()->SetNdivisions(505, kTRUE);
  grb->GetXaxis()->SetTimeDisplay(1);
  grb->GetXaxis()->SetTimeFormat("%d/%m");
  cv->Write();
 
  // TIB Details 
  if(beam_group == LAS::TIB || beam_group == LAS::AT || beam_group == LAS::ALL){
    cv = new TCanvas("TIB_ER","TIB_ER"); 
    TMultiGraph* gre = avec_draw(time, (tib_er_dx | tib_er_dy | tib_er_rx| tib_er_ry | tib_er_rz | tib_er_tz ),"TIB vs TOB (16 + 6p fit)","Date","#sigma#Deltax, #sigma#Deltay [#mum], #sigmaRx,#sigmaRy,#sigmaRz [#murad], #sigmaTz[#murad/m] ","AP");
    gre->GetXaxis()->SetNdivisions(505, kTRUE);
    gre->GetXaxis()->SetTimeDisplay(1);
    gre->GetXaxis()->SetTimeFormat("%d/%m");
 
    std::vector<std::string> Legende;
    Legende.push_back("#sigma#Delta x");
    Legende.push_back("#sigma#Delta y");
    Legende.push_back("#sigmaRx");
    Legende.push_back("#sigmaRy");
    Legende.push_back("#sigmaRz");
    Legende.push_back("#sigmaTz");
    AddLegend(gre, Legende);
    cv->Write();

    cv = new TCanvas("CH_TIB","CH_TIB"); 
    TMultiGraph* grcp = avec_draw(time, (tib_chi | tib_chi0),"TIB vs TOB (16 + 6p fit)","Date","#chi2/ndf","AP");
    grcp->GetXaxis()->SetNdivisions(505, kTRUE);
    grcp->GetXaxis()->SetTimeDisplay(1);
    grcp->GetXaxis()->SetTimeFormat("%d/%m");

    grcp->GetYaxis()->SetRangeUser(0.1, 1e4);
    
    std::vector<std::string> Legendcp;
    Legendcp.push_back("TIB after fit");
    Legendcp.push_back("TIB before fit");
    AddLegend(grcp, Legendcp);
    cv->SetLogy();
    cv->SetGridy();
    cv->Write();
  }

  // TEC+ Details
  if(beam_group == LAS::TEC_PLUS_AT || beam_group == LAS::TEC_AT || beam_group == LAS::AT || beam_group == LAS::ALL){
    cv = new TCanvas("TECP_ER","TECP_ER"); 
    TMultiGraph* grpe = avec_draw(time, (tecp_er_dx | tecp_er_dy | tecp_er_rz),"TEC+ vs TOB (16 + 3p fit)","Date","#sigma#Deltax, #sigma#Deltay [#mum], #sigmaRz [#murad] ","AP");
    grpe->GetXaxis()->SetNdivisions(505, kTRUE);
    grpe->GetXaxis()->SetTimeDisplay(1);
    grpe->GetXaxis()->SetTimeFormat("%d/%m");
 
    std::vector<std::string> Legendpe;
    Legendpe.push_back("#sigma#Delta x");
    Legendpe.push_back("#sigma#Delta y");
    Legendpe.push_back("#sigmaRz");
    AddLegend(grpe, Legendpe);
    cv->Write();

    cv = new TCanvas("CH_TECP","CH_TECP"); 
    TMultiGraph* grpcp = avec_draw(time, (tecp_chi | tecp_chi0),"TEC+ vs TOB (16 + 3p fit)","Date","#chi2/ndf","AP");
    grpcp->GetXaxis()->SetNdivisions(505, kTRUE);
    grpcp->GetXaxis()->SetTimeDisplay(1);
    grpcp->GetXaxis()->SetTimeFormat("%d/%m");

    grpcp->GetYaxis()->SetRangeUser(0.1, 1e4);
 
    std::vector<std::string> Legendpcp;
    Legendpcp.push_back("TEC+ after fit");
    Legendpcp.push_back("TEC+ before fit");
    AddLegend(grpcp, Legendpcp);
    cv->SetLogy();
    cv->SetGridy();

    cv->Write();
  }

  // TEC- Details
  if(beam_group == LAS::TEC_MINUS_AT || beam_group == LAS::TEC_AT || beam_group == LAS::AT || beam_group == LAS::ALL){
    cv = new TCanvas("TECM_ER","TECM_ER"); 
    TMultiGraph* grme = avec_draw(time, (tecm_er_dx | tecm_er_dy | tecm_er_rz),"TEC- vs TOB (16 + 3p fit)","Date","#sigma#Deltax, #sigma#Deltay [#mum], #sigmaRz [#murad] ","AP");
    grme->GetXaxis()->SetNdivisions(505, kTRUE);
    grme->GetXaxis()->SetTimeDisplay(1);
    grme->GetXaxis()->SetTimeFormat("%d/%m");
 
    std::vector<std::string> Legendme;
    Legendme.push_back("#sigma#Delta x");
    Legendme.push_back("#sigma#Delta y");
    Legendme.push_back("#sigmaRz");
    AddLegend(grme, Legendme);
    cv->Write();


    cv = new TCanvas("CH_TECM","CH_TECM"); 
    TMultiGraph* grmcp = avec_draw(time, (tecm_chi | tecm_chi0),"TEC- vs TOB (16 + 3p fit)","Date","#chi2/ndf","AP");
    grmcp->GetXaxis()->SetNdivisions(505, kTRUE);
    grmcp->GetXaxis()->SetTimeDisplay(1);
    grmcp->GetXaxis()->SetTimeFormat("%d/%m");

    grmcp->GetYaxis()->SetRangeUser(0.1, 1e4);
 
    std::vector<std::string> Legendmcp;
    Legendmcp.push_back("TEC- after fit");
    Legendmcp.push_back("TEC- before fit");
    AddLegend(grmcp, Legendmcp);
    cv->SetLogy();
    cv->SetGridy();
    cv->Write();
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Save Avecs with patrameters, errors and time
  details_at->mkdir("data");
  details_at->cd("data");

  // Common Data
  time.Write("time");

  tob_chi.Write("tob_chi");
  tob_chi0.Write("tob_chi0");
  AT_chi.Write("AT_chi");
  AT_chi0.Write("AT_chi0");

  beam_slp[0].Write("beam_slp_0");
  beam_slp[1].Write("beam_slp_1");
  beam_slp[2].Write("beam_slp_2");
  beam_slp[3].Write("beam_slp_3");
  beam_slp[4].Write("beam_slp_4");
  beam_slp[5].Write("beam_slp_5");
  beam_slp[6].Write("beam_slp_6");
  beam_slp[7].Write("beam_slp_7");

  beam_int[0].Write("beam_int_0");
  beam_int[1].Write("beam_int_1");
  beam_int[2].Write("beam_int_2");
  beam_int[3].Write("beam_int_3");
  beam_int[4].Write("beam_int_4");
  beam_int[5].Write("beam_int_5");
  beam_int[6].Write("beam_int_6");
  beam_int[7].Write("beam_int_7");

  b_chi[0].Write("b_chi_0");
  b_chi[1].Write("b_chi_1");
  b_chi[2].Write("b_chi_2");
  b_chi[3].Write("b_chi_3");
  b_chi[4].Write("b_chi_4");
  b_chi[5].Write("b_chi_5");
  b_chi[6].Write("b_chi_6");
  b_chi[7].Write("b_chi_7");

  b_chi0[0].Write("b_chi0_0");
  b_chi0[1].Write("b_chi0_1");
  b_chi0[2].Write("b_chi0_2");
  b_chi0[3].Write("b_chi0_3");
  b_chi0[4].Write("b_chi0_4");
  b_chi0[5].Write("b_chi0_5");
  b_chi0[6].Write("b_chi0_6");
  b_chi0[7].Write("b_chi0_7");

  // TIB Data
  if(beam_group == LAS::TIB || beam_group == LAS::AT || beam_group == LAS::ALL){
    tib_dx.Write("tib_dx ");
    tib_dy.Write("tib_dy ");
    tib_rx.Write("tib_rx");
    tib_ry.Write("tib_ry ");
    tib_rz.Write("tib_rz ");
    tib_tz.Write("tib_tz ");
    tib_er_dx.Write("tib_er_dx ");
    tib_er_dy.Write("tib_er_dy ");
    tib_er_rx.Write("tib_er_rx");
    tib_er_ry.Write("tib_er_ry ");
    tib_er_rz.Write("tib_er_rz ");
    tib_er_tz.Write("tib_er_tz");
    tib_chi.Write("tib_chi");
    tib_chi0.Write("tib_chi0");
  }

  // TEC+ Data
  if(beam_group == LAS::TEC_PLUS_AT || beam_group == LAS::TEC_AT || beam_group == LAS::AT || beam_group == LAS::ALL){
    tecp_dx.Write("tecp_dx ");
    tecp_dy.Write("tecp_dy ");
    tecp_rz.Write("tecp_rz ");
    tecp_er_dx.Write("tecp_er_dx ");
    tecp_er_dy.Write("tecp_er_dy ");
    tecp_er_rz.Write("tecp_er_rz ");
    tecp_chi.Write("tecp_chi");
    tecp_chi0.Write("tecp_chi0");
  }

  // TEC- Data
  if(beam_group == LAS::TEC_MINUS_AT || beam_group == LAS::TEC_AT || beam_group == LAS::AT || beam_group == LAS::ALL){
    tecm_dx.Write("tecm_dx ");
    tecm_dy.Write("tecm_dy ");
    tecm_rz.Write("tecm_rz ");
    tecm_er_dx.Write("tecm_er_dx ");
    tecm_er_dy.Write("tecm_er_dy ");
    tecm_er_rz.Write("tecm_er_rz ");
    tecm_chi.Write("tecm_chi");
    tecm_chi0.Write("tecm_chi0");
  }

  of.Close();

  gROOT->SetBatch(batch_mode); // Return to original batch mode
}


// void DrawTIBPar(std::vector<AtPar>& parlist, Avec& time, bool quick_draw){
//   Avec tib_dx;
//   Avec tib_dy;
//   Avec tib_rx;
//   Avec tib_ry;
//   Avec tib_rz;
//   Avec tib_tz;
//   Avec tib_er_dx;
//   Avec tib_er_dy;
//   Avec tib_er_rx;
//   Avec tib_er_ry;
//   Avec tib_er_rz;
//   Avec tib_er_tz;
//   Avec tib_chi;
//   Avec tib_chi0;
//   Avec tob_chi;
//   Avec tob_chi0;
//   Avec AT_chi;
//   Avec AT_chi0;
//   Avec beam_slp[8];
//   Avec beam_int[8];
//   Avec b_chi[8];
//   Avec b_chi0[8];
//   for(unsigned int i = 0; i < parlist.size(); i++){
//     tib_dx.push_back(parlist[i].tib_dx*1000);
//     tib_dy.push_back(parlist[i].tib_dy*1000);
//     tib_rx.push_back(parlist[i].tib_rx*1000000);
//     tib_ry.push_back(parlist[i].tib_ry*1000000);
//     tib_rz.push_back(parlist[i].tib_rz*1000000);
//     tib_tz.push_back(parlist[i].tib_tz*1000000000);

//     // tib_er_dx.push_back(parlist[i].er_tib_dx*1000);
//     // tib_er_dy.push_back(parlist[i].er_tib_dy*1000);
//     // tib_er_rx.push_back(parlist[i].er_tib_rx*1000000);
//     // tib_er_ry.push_back(parlist[i].er_tib_ry*1000000);
//     // tib_er_rz.push_back(parlist[i].er_tib_rz*1000000);
//     // tib_er_tz.push_back(parlist[i].er_tib_tz*1000000000);

//     tib_er_dx.push_back(parlist[i].er_tib_dx);
//     tib_er_dy.push_back(parlist[i].er_tib_dy);
//     tib_er_rx.push_back(parlist[i].er_tib_rx*1000);
//     tib_er_ry.push_back(parlist[i].er_tib_ry*1000);
//     tib_er_rz.push_back(parlist[i].er_tib_rz*1000);
//     tib_er_tz.push_back(parlist[i].er_tib_tz*1000000);

//     tib_chi.push_back(parlist[i].tib_chi);
//     tib_chi0.push_back(parlist[i].tib_chi0);
//     tob_chi.push_back(parlist[i].tob_chi);
//     tob_chi0.push_back(parlist[i].tob_chi0);
//     AT_chi.push_back(parlist[i].AT_chi);
//     AT_chi0.push_back(parlist[i].AT_chi0);
//     for(int nb = 0; nb < 8; nb++){
//       beam_slp[nb].push_back(parlist[i].beam_a[nb]*1000000*564.);
//       beam_int[nb].push_back(parlist[i].beam_b[nb]*1000*564. + (LAS::z_at_bs + 0.0 ) * parlist[i].beam_a[nb]*1000*564. );
//       b_chi[nb].push_back(parlist[i].b_chi[nb]);
//       b_chi0[nb].push_back(parlist[i].b_chi0[nb]);
//     }
//   }

//   TFile f_tib_at("TIB_AT.root","RECREATE");
//   time -= 14400;
//   (new TCanvas("TIB_AT","TIB_AT"))->SetGridy(); 
//   TMultiGraph* gr = avec_draw(time, (tib_dx | tib_dy | tib_rx| tib_ry | tib_rz | tib_tz ), (tib_er_dx | tib_er_dy | tib_er_rx| tib_er_ry | tib_er_rz | tib_er_tz ),"TIB vs TOB (16 + 6p fit)","Date","#Deltax, #Deltay [#mum], Rx,Ry,Rz [#murad], Tz[#murad/m] ","AP");
//   gr->GetXaxis()->SetNdivisions(505, kTRUE);
//   gr->GetXaxis()->SetTimeDisplay(1);
//   gr->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
//   std::vector<std::string> Legend;
//   Legend.push_back("#Delta x");
//   Legend.push_back("#Delta y");
//   Legend.push_back("Rx");
//   Legend.push_back("Ry");
//   Legend.push_back("Rz");
//   Legend.push_back("Tz");
//   AddLegend(gr, Legend, 1, 0.9, 0.6, 1.0, 1.0);
//   f_tib_at.Close();

//   TFile f("TIB_AT_parameters.root", "RECREATE");
//   tib_dx.Write("tib_dx ");
//   tib_dy.Write("tib_dy ");
//   tib_rx.Write("tib_rx");
//   tib_ry.Write("tib_ry ");
//   tib_rz.Write("tib_rz ");
//   tib_tz.Write("tib_tz ");
//   tib_er_dx.Write("tib_er_dx ");
//   tib_er_dy.Write("tib_er_dy ");
//   tib_er_rx.Write("tib_er_rx");
//   tib_er_ry.Write("tib_er_ry ");
//   tib_er_rz.Write("tib_er_rz ");
//   tib_er_tz.Write("tib_er_tz");
//   time.Write("time");
//   f.Close();

//   if(quick_draw) return;

//   new TCanvas("TIB_ER","TIB_ER"); 
//   TMultiGraph* gre = avec_draw(time, (tib_er_dx | tib_er_dy | tib_er_rx| tib_er_ry | tib_er_rz | tib_er_tz ),"TIB vs TOB (16 + 6p fit)","Date","#sigma#Deltax, #sigma#Deltay [#mum], #sigmaRx,#sigmaRy,#sigmaRz [#murad], #sigmaTz[#murad/m] ","AP");
//   gre->GetXaxis()->SetNdivisions(505, kTRUE);
//   gre->GetXaxis()->SetTimeDisplay(1);
//   gre->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legende;
//   Legende.push_back("#sigma#Delta x");
//   Legende.push_back("#sigma#Delta y");
//   Legende.push_back("#sigmaRx");
//   Legende.push_back("#sigmaRy");
//   Legende.push_back("#sigmaRz");
//   Legende.push_back("#sigmaTz");
//   AddLegend(gre, Legende);


//   new TCanvas("CH_TOB","CH_TOB"); 
//   TMultiGraph* grco = avec_draw(time, (tob_chi | tob_chi0),"TIB vs TOB (16 + 6p fit)","Date","#chi2/ndf","AP");
//   grco->GetXaxis()->SetNdivisions(505, kTRUE);
//   grco->GetXaxis()->SetTimeDisplay(1);
//   grco->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legendco;
//   Legendco.push_back("TOB after fit");
//   Legendco.push_back("TOB before fit");
//   AddLegend(grco, Legendco);

//   new TCanvas("CH_TIB","CH_TIB"); 
//   TMultiGraph* grcp = avec_draw(time, (tib_chi | tib_chi0),"TIB vs TOB (16 + 6p fit)","Date","#chi2/ndf","AP");
//   grcp->GetXaxis()->SetNdivisions(505, kTRUE);
//   grcp->GetXaxis()->SetTimeDisplay(1);
//   grcp->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legendcp;
//   Legendcp.push_back("TIB after fit");
//   Legendcp.push_back("TIB before fit");
//   AddLegend(grcp, Legendcp);

//   new TCanvas("CH_AT","CH_AT"); 
//   TMultiGraph* grca = avec_draw(time, (AT_chi | AT_chi0),"TIB vs TOB (16 + 6p fit)","Date","#chi2/ndf","AP");
//   grca->GetXaxis()->SetNdivisions(505, kTRUE);
//   grca->GetXaxis()->SetTimeDisplay(1);
//   grca->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legendca;
//   Legendca.push_back("AT after fit");
//   Legendca.push_back("AT before fit");
//   AddLegend(grca, Legendca);

//   new TCanvas("TIB_A","TIB_A"); 
//   TMultiGraph* gra = avec_draw(time, (beam_slp[0] | beam_slp[1] | beam_slp[2] | beam_slp[3] | beam_slp[4] | beam_slp[5] | beam_slp[6] | beam_slp[7]),"TIB vs TOB (16 + 6p fit)","Date","Beam slopes[#murad] ","AP");
//   gra->GetXaxis()->SetNdivisions(505, kTRUE);
//   gra->GetXaxis()->SetTimeDisplay(1);
//   gra->GetXaxis()->SetTimeFormat("%d/%m");
 
//   new TCanvas("TIB_B","TIB_B"); 
//   TMultiGraph* grb = avec_draw(time, (beam_int[0] | beam_int[1] | beam_int[2] | beam_int[3] | beam_int[4] | beam_int[5] | beam_int[6] | beam_int[7]),"TIB vs TOB (16 + 6p fit)","Date","Beam intersepts[#mum] ","AP");
//   grb->GetXaxis()->SetNdivisions(505, kTRUE);
//   grb->GetXaxis()->SetTimeDisplay(1);
//   grb->GetXaxis()->SetTimeFormat("%d/%m");
 
//   TCanvas * cvx=new TCanvas("TIB #Deltax","TIB Dx ");
//   cvx->cd(1);
//   avec_plot(tib_dx,200,"TIB #Deltax","TIB #Deltax [#mum]");
//   TCanvas * cvy=new TCanvas("TIB #Deltay","TIB Dy ");
//   cvy->cd(1);
//   avec_plot(tib_dy,200,"TIB #Deltay","TIB #Deltay [#mum]");
//   TCanvas * cvrx=new TCanvas("TIB Rx","TIB Rx ");
//   cvrx->cd(1);
//   avec_plot(tib_rx,200,"TIB Rx","Rx [#murad]");
//   TCanvas * cvry=new TCanvas("TIB Ry","TIB Ry ");
//   cvry->cd(1);
//   avec_plot(tib_ry,200,"TIB Ry","Ry [#murad]");
//   TCanvas * cvrz=new TCanvas("TIB Rz","TIB Rz ");
//   cvrz->cd(1);
//   avec_plot(tib_rz,200,"TIB Rz","Rz [#murad]");
//   TCanvas * cvtz=new TCanvas("TIB Tz","TIB Tz ");
//   cvtz->cd(1);
//   avec_plot(tib_tz,200,"TIB Tz","Tz [#murad]");

// }

// void DrawTECPPar(std::vector<AtPar>& parlist, Avec& time, bool quick_draw){
//   Avec tecp_dx;
//   Avec tecp_dy;
//   Avec tecp_rz;
//   Avec tecp_chi;
//   Avec tecp_chi0;
//   Avec tob_chi;
//   Avec tob_chi0;
//   Avec AT_chi;
//   Avec AT_chi0;
//   Avec beam_slp[8];
//   Avec beam_int[8];
//   Avec b_chi[8];
//   Avec b_chi0[8];
//   Avec tecp_er_dx;
//   Avec tecp_er_dy;
//   Avec tecp_er_rz;
//   for(unsigned int i = 0; i < parlist.size(); i++){
//     tecp_dx.push_back(parlist[i].tecp_dx*1000);
//     tecp_dy.push_back(parlist[i].tecp_dy*1000);
//     tecp_rz.push_back(parlist[i].tecp_rz*1000000);

//     tecp_er_dx.push_back(parlist[i].er_tecp_dx);
//     tecp_er_dy.push_back(parlist[i].er_tecp_dy);
//     tecp_er_rz.push_back(parlist[i].er_tecp_rz*1000);

//     tecp_chi.push_back(parlist[i].tecp_chi);
//     tecp_chi0.push_back(parlist[i].tecp_chi0);
//     tob_chi.push_back(parlist[i].tob_chi);
//     tob_chi0.push_back(parlist[i].tob_chi0);
//     AT_chi.push_back(parlist[i].AT_chi);
//     AT_chi0.push_back(parlist[i].AT_chi0);
//     for(int nb = 0; nb < 8; nb++){
//       beam_slp[nb].push_back(parlist[i].beam_a[nb]*1000000*564.);
//       beam_int[nb].push_back(parlist[i].beam_b[nb]*1000*564.);
//       b_chi[nb].push_back(parlist[i].b_chi[nb]);
//       b_chi0[nb].push_back(parlist[i].b_chi0[nb]);
//     }
//   }

//   TFile f_tecp_at("TECP_AT.root","RECREATE");
//   time -= 14400;

//   (new TCanvas("TECP_AT","TECP_AT"))->SetGridy(); 
//   TMultiGraph* gr = avec_draw(time, (tecp_dx | tecp_dy | tecp_rz), (tecp_er_dx | tecp_er_dy | tecp_er_rz),"TEC+ vs TOB (16 + 3p fit)","Date","#Deltax, #Deltay [#mum], Rz [#murad] ","AP");
//   gr->GetXaxis()->SetNdivisions(505, kTRUE);
//   gr->GetXaxis()->SetTimeDisplay(1);
//   gr->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legend;
//   Legend.push_back("#Delta x");
//   Legend.push_back("#Delta y");
//   Legend.push_back("Rz");
//   AddLegend(gr, Legend, 1, 0.9, 0.6, 1.0, 1.0);

//   f_tecp_at.Close();

//   TFile f("TECP_AT_parameters.root", "RECREATE");

//   tecp_dx.Write("tecp_dx ");
//   tecp_dy.Write("tecp_dy ");
//   tecp_rz.Write("tecp_rz ");
//   tecp_er_dx.Write("tecp_er_dx ");
//   tecp_er_dy.Write("tecp_er_dy ");
//   tecp_er_rz.Write("tecp_er_rz ");

//   time.Write("time");
//   f.Close();

//   if(quick_draw) return;

//   new TCanvas("TECP_ER","TECP_ER"); 
//   TMultiGraph* gre = avec_draw(time, (tecp_er_dx | tecp_er_dy | tecp_er_rz),"TEC+ vs TOB (16 + 3p fit)","Date","#sigma#Deltax, #sigma#Deltay [#mum], #sigmaRz [#murad] ","AP");
//   gre->GetXaxis()->SetNdivisions(505, kTRUE);
//   gre->GetXaxis()->SetTimeDisplay(1);
//   gre->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legende;
//   Legende.push_back("#sigma#Delta x");
//   Legende.push_back("#sigma#Delta y");
//   Legende.push_back("#sigmaRz");
//   AddLegend(gre, Legende);


//   new TCanvas("CH_TOB","CH_TOB"); 
//   TMultiGraph* grco = avec_draw(time, (tob_chi | tob_chi0),"TEC+ vs TOB (16 + 3p fit)","Date","#chi2/ndf","AP");
//   grco->GetXaxis()->SetNdivisions(505, kTRUE);
//   grco->GetXaxis()->SetTimeDisplay(1);
//   grco->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legendco;
//   Legendco.push_back("TOB after fit");
//   Legendco.push_back("TOB before fit");
//   AddLegend(grco, Legendco);

//   new TCanvas("CH_TECP","CH_TECP"); 
//   TMultiGraph* grcp = avec_draw(time, (tecp_chi | tecp_chi0),"TEC+ vs TOB (16 + 3p fit)","Date","#chi2/ndf","AP");
//   grcp->GetXaxis()->SetNdivisions(505, kTRUE);
//   grcp->GetXaxis()->SetTimeDisplay(1);
//   grcp->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legendcp;
//   Legendcp.push_back("TEC+ after fit");
//   Legendcp.push_back("TEC+ before fit");
//   AddLegend(grcp, Legendcp);

//   new TCanvas("CH_AT","CH_AT"); 
//   TMultiGraph* grca = avec_draw(time, (AT_chi | AT_chi0),"TEC+ vs TOB (16 + 3p fit)","Date","#chi2/ndf","AP");
//   grca->GetXaxis()->SetNdivisions(505, kTRUE);
//   grca->GetXaxis()->SetTimeDisplay(1);
//   grca->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legendca;
//   Legendca.push_back("AT after fit");
//   Legendca.push_back("AT before fit");
//   AddLegend(grca, Legendca);

//   new TCanvas("TECP_A","TECP_A"); 
//   TMultiGraph* gra = avec_draw(time, (beam_slp[0] | beam_slp[1] | beam_slp[2] | beam_slp[3] | beam_slp[4] | beam_slp[5] | beam_slp[6] | beam_slp[7]),"TEC+ vs TOB (16 + 3p fit)","Date","Beam slopes[#murad] ","AP");
//   gra->GetXaxis()->SetNdivisions(505, kTRUE);
//   gra->GetXaxis()->SetTimeDisplay(1);
//   gra->GetXaxis()->SetTimeFormat("%d/%m");
 
//   new TCanvas("TECP_B","TECP_B"); 
//   TMultiGraph* grb = avec_draw(time, (beam_int[0] | beam_int[1] | beam_int[2] | beam_int[3] | beam_int[4] | beam_int[5] | beam_int[6] | beam_int[7]),"TEC+ vs TOB (16 + 3p fit)","Date","Beam intersepts[#mum] ","AP");
//   grb->GetXaxis()->SetNdivisions(505, kTRUE);
//   grb->GetXaxis()->SetTimeDisplay(1);
//   grb->GetXaxis()->SetTimeFormat("%d/%m");
 
//   TCanvas * cvx=new TCanvas("TEC+ #Deltax","TEC+ Dx ");
//   cvx->cd(1);
//   avec_plot(tecp_dx,200,"TEC+ #Deltax","TEC+ #Deltax [#mum]");
//   TCanvas * cvy=new TCanvas("TEC+ #Deltay","TEC+ Dy ");
//   cvy->cd(1);
//   avec_plot(tecp_dy,200,"TEC+ #Deltay","TEC+ #Deltay [#mum]");
//   TCanvas * cvrz=new TCanvas("TEC+ Rz","TEC+ Rz ");
//   cvrz->cd(1);
//   avec_plot(tecp_rz,200,"TEC+ Rz","TEC+ Rz [#murad]");

// }

// void DrawTECMPar(std::vector<AtPar>& parlist, Avec& time, bool quick_draw)
// {
//   Avec tecm_dx;
//   Avec tecm_dy;
//   Avec tecm_rz;
//   Avec tecm_chi;
//   Avec tecm_chi0;
//   Avec tob_chi;
//   Avec tob_chi0;
//   Avec AT_chi;
//   Avec AT_chi0;
//   Avec beam_slp[8];
//   Avec beam_int[8];
//   Avec b_chi[8];
//   Avec b_chi0[8];
//   Avec tecm_er_dx;
//   Avec tecm_er_dy;
//   Avec tecm_er_rz;
//   for(unsigned int i = 0; i < parlist.size(); i++){
//     tecm_dx.push_back(parlist[i].tecm_dx*1000);
//     tecm_dy.push_back(parlist[i].tecm_dy*1000);
//     tecm_rz.push_back(parlist[i].tecm_rz*1000000);

//     tecm_er_dx.push_back(parlist[i].er_tecm_dx);
//     tecm_er_dy.push_back(parlist[i].er_tecm_dy);
//     tecm_er_rz.push_back(parlist[i].er_tecm_rz*1000);

//     tecm_chi.push_back(parlist[i].tecm_chi);
//     tecm_chi0.push_back(parlist[i].tecm_chi0);
//     tob_chi.push_back(parlist[i].tob_chi);
//     tob_chi0.push_back(parlist[i].tob_chi0);
//     AT_chi.push_back(parlist[i].AT_chi);
//     AT_chi0.push_back(parlist[i].AT_chi0);
//     for(int nb = 0; nb < 8; nb++){
//       beam_slp[nb].push_back(parlist[i].beam_a[nb]*1000000*564.);
//       beam_int[nb].push_back(parlist[i].beam_b[nb]*1000*564.);
//       b_chi[nb].push_back(parlist[i].b_chi[nb]);
//       b_chi0[nb].push_back(parlist[i].b_chi0[nb]);
//     }
//   }

//   TFile f_tecm_at("TECM_AT.root","RECREATE");
//   time -= 14400;

//   (new TCanvas("TECM_AT","TECM_AT"))->SetGridy(); 
//   TMultiGraph* gr = avec_draw(time, (tecm_dx | tecm_dy | tecm_rz), (tecm_er_dx | tecm_er_dy | tecm_er_rz),"TEC- vs TOB (16 + 3p fit)","Date","#Deltax, #Deltay [#mum], Rz [#murad] ","AP");
//   gr->GetXaxis()->SetNdivisions(505, kTRUE);
//   gr->GetXaxis()->SetTimeDisplay(1);
//   gr->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legend;
//   Legend.push_back("#Delta x");
//   Legend.push_back("#Delta y");
//   Legend.push_back("Rz");
//   AddLegend(gr, Legend, 1, 0.9, 0.6, 1.0, 1.0);

//   f_tecm_at.Close();

//   TFile f("TECM_AT_parameters.root", "RECREATE");

//   tecm_dx.Write("tecm_dx ");
//   tecm_dy.Write("tecm_dy ");
//   tecm_rz.Write("tecm_rz ");
//   tecm_er_dx.Write("tecm_er_dx ");
//   tecm_er_dy.Write("tecm_er_dy ");
//   tecm_er_rz.Write("tecm_er_rz ");

//   time.Write("time");
//   f.Close();

//   if(quick_draw) return;

//   new TCanvas("TECM_ER","TECM_ER"); 
//   TMultiGraph* gre = avec_draw(time, (tecm_er_dx | tecm_er_dy | tecm_er_rz),"TEC- vs TOB (16 + 3p fit)","Date","#sigma#Deltax, #sigma#Deltay [#mum], #sigmaRz [#murad] ","AP");
//   gre->GetXaxis()->SetNdivisions(505, kTRUE);
//   gre->GetXaxis()->SetTimeDisplay(1);
//   gre->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legende;
//   Legende.push_back("#sigma#Delta x");
//   Legende.push_back("#sigma#Delta y");
//   Legende.push_back("#sigmaRz");
//   AddLegend(gre, Legende);


//   new TCanvas("CH_TOB","CH_TOB"); 
//   TMultiGraph* grco = avec_draw(time, (tob_chi | tob_chi0),"TEC- vs TOB (16 + 3p fit)","Date","#chi2/ndf","AP");
//   grco->GetXaxis()->SetNdivisions(505, kTRUE);
//   grco->GetXaxis()->SetTimeDisplay(1);
//   grco->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legendco;
//   Legendco.push_back("TOB after fit");
//   Legendco.push_back("TOB before fit");
//   AddLegend(grco, Legendco);

//   new TCanvas("CH_TECM","CH_TECM"); 
//   TMultiGraph* grcp = avec_draw(time, (tecm_chi | tecm_chi0),"TEC- vs TOB (16 + 3p fit)","Date","#chi2/ndf","AP");
//   grcp->GetXaxis()->SetNdivisions(505, kTRUE);
//   grcp->GetXaxis()->SetTimeDisplay(1);
//   grcp->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legendcp;
//   Legendcp.push_back("TEC- after fit");
//   Legendcp.push_back("TEC- before fit");
//   AddLegend(grcp, Legendcp);

//   new TCanvas("CH_AT","CH_AT"); 
//   TMultiGraph* grca = avec_draw(time, (AT_chi | AT_chi0),"TEC- vs TOB (16 + 3p fit)","Date","#chi2/ndf","AP");
//   grca->GetXaxis()->SetNdivisions(505, kTRUE);
//   grca->GetXaxis()->SetTimeDisplay(1);
//   grca->GetXaxis()->SetTimeFormat("%d/%m");
 
//   std::vector<std::string> Legendca;
//   Legendca.push_back("AT after fit");
//   Legendca.push_back("AT before fit");
//   AddLegend(grca, Legendca);

//   new TCanvas("TECM_A","TECM_A"); 
//   TMultiGraph* gra = avec_draw(time, (beam_slp[0] | beam_slp[1] | beam_slp[2] | beam_slp[3] | beam_slp[4] | beam_slp[5] | beam_slp[6] | beam_slp[7]),"TEC- vs TOB (16 + 3p fit)","Date","Beam slopes[#murad] ","AP");
//   gra->GetXaxis()->SetNdivisions(505, kTRUE);
//   gra->GetXaxis()->SetTimeDisplay(1);
//   gra->GetXaxis()->SetTimeFormat("%d/%m");
 
//   new TCanvas("TECM_B","TECM_B"); 
//   TMultiGraph* grb = avec_draw(time, (beam_int[0] | beam_int[1] | beam_int[2] | beam_int[3] | beam_int[4] | beam_int[5] | beam_int[6] | beam_int[7]),"TEC- vs TOB (16 + 3p fit)","Date","Beam intersepts[#mum] ","AP");
//   grb->GetXaxis()->SetNdivisions(505, kTRUE);
//   grb->GetXaxis()->SetTimeDisplay(1);
//   grb->GetXaxis()->SetTimeFormat("%d/%m");
 
//   TCanvas * cvx=new TCanvas("TEC- #Deltax","TEC- Dx ");
//   cvx->cd(1);
//   avec_plot(tecm_dx,200,"TEC- #Deltax","TEC- #Deltax [#mum]");
//   TCanvas * cvy=new TCanvas("TEC- #Deltay","TEC- Dy ");
//   cvy->cd(1);
//   avec_plot(tecm_dy,200,"TEC- #Deltay","TEC- #Deltay [#mum]");
//   TCanvas * cvrz=new TCanvas("TEC- Rz","TEC- Rz ");
//   cvrz->cd(1);
//   avec_plot(tecm_rz,200,"TEC- Rz","TEC- Rz [#murad]");

// }

