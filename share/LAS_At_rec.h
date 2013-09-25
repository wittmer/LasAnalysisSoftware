#ifndef __INCLUDE_LAS_AT_REC_H__
#define __INCLUDE_LAS_AT_REC_H__

#include "TMatrixT.h"

LASGlobalData<int> group_mask(LAS::beam_group beam_group);

void AT_Parameter_Reconstruction(const std::string& run_list_name, const std::string& path_prefix="", const std::string& outfile = "AT_latest_results.root", LAS::FIT_TYPE::fit_type fit_type = LAS::FIT_TYPE::AT, bool quick_draw = true);

// Create Coefficients for Alignment parameters
void create_coefficients_TIB( std::vector< LASGlobalData<double> >& coeff_list);
void create_coefficients_AT_beams( std::vector< LASGlobalData<double> >& coeff_list, LAS::beam_group beam_group = LAS::AT);
void create_coefficients_TEC_beams( std::vector< LASGlobalData<double> >& coeff_list, LAS::beam_group beam_group = LAS::TEC);
void create_coefficients_TEC_discs( std::vector< LASGlobalData<double> >& coeff_list, LAS::beam_group beam_group = LAS::TEC);

double global_data_sum(const LASGlobalData<double>& gd, LAS::beam_group beam_group = LAS::ALL);
void Calc_chi2(const LASGlobalData<double>& data, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, int ndof);

void LAS_Alpar_fit(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, const std::vector< LASGlobalData<double> >& coeff_list, std::vector<double>& par_list, TMatrixT<double>& U);

void AlPar_fit_old(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, bool swap_flag = false, bool verbose_flag = false);
void AT_residuals(const std::string& ref_filename, const std::string& filename, unsigned int blnr = 0, LAS::FIT_TYPE::fit_type fit_type = LAS::FIT_TYPE::AT, double tweak_tecm_tor = 0, double tweak_tecm_shx = 0, double tweak_tecm_shy = 0);

// Fit a set of alignment parameters
void LAS_fit(const LASGlobalData<double>& dif_rad, const LASGlobalData<double>& err2_rad2, const LASGlobalData<int>& mask, LasAlPar& results, LASGlobalData<double>& rec_pos, LAS::FIT_TYPE::fit_type fit_type = LAS::FIT_TYPE::AT, bool verbose = false);


// Draw the AT alignment parameters
void DrawATPar(std::vector<AtPar>& parlist, Avec& time, Bool_t quick_draw = kTRUE, const std::string& output_file = "AT_parameters.root", LAS::beam_group beam_group = LAS::AT);
void DrawTIBPar(std::vector<AtPar>& parlist, Avec& time, bool quick_draw = false);
void DrawTECPPar(std::vector<AtPar>& parlist, Avec& time, bool quick_draw = false);
void DrawTECMPar(std::vector<AtPar>& parlist, Avec& time, bool quick_draw = false);

// Draw the AT alignment parameters
void DrawAtPar_old(std::vector<AtPar>& parlist, Avec& time);

void AT_TIB_22par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose = false);
void AT_TIB_21par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose = false);
void AT_TIB_TECM_25par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose = false);
void AT_TIB_TECP(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose = false);
void AT_full(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos, bool verbose = false);
// Fit only the AT beams to TOB
void AT_TOB_16par(const LASGlobalData<double>& dif, const LASGlobalData<double>& err, const LASGlobalData<int>& mask, AtPar& results, LASGlobalData<double>& rec_pos);
// reconstruction of TIB wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TIB2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,const LASGlobalData<int>& mask, AtPar& results);
void AT_TECP2TOB_3par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,const LASGlobalData<int>& mask, AtPar& results);
void AT_TECM2TOB_3par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,const LASGlobalData<int>& mask, AtPar& results);
void AT_28par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,const LASGlobalData<int>& mask, AtPar& results);

LASGlobalData<double> TIB_spot_rec(const AtPar& par); //calculation of beam spot positions at TIB&TOB modules
LASGlobalData<double> TECP_spot_rec(const AtPar& par); //calculation of beam spot positions at TEC+&TOB modules
LASGlobalData<double> TECM_spot_rec(const AtPar& par); //calculation of beam spot positions at TEC-&TOB modules
LASGlobalData<double> AT_spot_rec(const AtPar& par); // calculation of beam spot positions for AT modules

// Calculation of chi2/ndf for AT reconstruction
void   AT_chi2(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask,  AtPar& par);
void  TIB_chi2(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask,  AtPar& par); //calculation of chi2/ndf at TIB&TOB
void TECP_chi2(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask,  AtPar& par); //calculation of chi2/ndf at TOB&TEC+
void TECM_chi2(const LASGlobalData<double>& dif, const LASGlobalData<double>& err2, const LASGlobalData<int>& mask,  AtPar& par); //calculation of chi2/ndf at TOB&TEC-

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculation of beam spot positions for AT modules
LASGlobalData<double> AT_spot_rec_old(const AtPar& par);

// Calculation of chi2/ndf for AT reconstruction
void AT_chi2_old(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,  AtPar& par);

// reconstruction of TIB wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TIB2TOB_6par_old(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results);

// reconstruction of TEC+ wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TECP2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results);

// reconstruction of TEC- wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TECM2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results);


#endif
