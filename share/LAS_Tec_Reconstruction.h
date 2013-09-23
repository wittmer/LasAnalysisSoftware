#ifndef LAS_TEC_RECONTRUCTION_H
#define LAS_TEC_RECONTRUCTION_H

#include "LAS_basic_tools.h"
#include "LASGlobalData.h"
#include "LAS_alpar.h"

void TEC_Parameter_Reconstruction(const std::string& result_file, const LASGlobalData<double>& ref_pos, const LASGlobalData<int>& ref_mask, const std::string& dir_name="tec_rec");
void check_alpar(const TecPar& alpar);
void print_sum_sq(const LASGlobalData<double>& data, const LASGlobalData<int>& mask);
void correct_signs(LASGlobalData<double>& data);
void get_sum_sq(const LASGlobalData<double>& data, const LASGlobalData<int>& mask, double& ssq_tecp, double& ssq_tecm);

//! TEC reconstruction
class TECReconstructor
{
 public:
  TECReconstructor();
  Avec make_theta();
  static Avec make_zdat(Int_t bs_frame=0);
  void Init_onering_symm();
  void Init_onering();
  void Init_full();
  void Init_common();

  double dphi_rec_tec(int subdet, int ring, int beam, int zpos, const TecRingPar& alpar, int print=0);
  double dphi_rec_tec(int subdet, int ring, int beam, int zpos, const TecPar& alpar, int print=0);
  LASGlobalData<double> dphi_rec_tec(const TecRingPar& alpar, LAS::beam_group bg, int print=0);
  LASGlobalData<double> dphi_rec_tec(const TecPar& alpar, LAS::beam_group bg, int print=0);
  LASGlobalData<double> dphi_rec_tec(const LasAlPar& alpar, int print=0);

  // Calculation of beam spot deviation for given alignment parameters
  double dx_rec_tec(int subdet, int ring, int beam, int zpos, const TecRingPar& alpar, int print=0);
  double dx_rec_tec(int subdet, int ring, int beam, int zpos, const TecPar& alpar, int print=0);

  LASGlobalData<double> dx_rec_tec(const TecRingPar& alpar, LAS::beam_group bg, int print=0);
  LASGlobalData<double> dx_rec_tec(const TecPar& alpar, LAS::beam_group bg, int print=0);
  LASGlobalData<double> dx_rec_tec(const LasAlPar& alpar, int print=0);

  int set_masked_modules_tec(LASGlobalData<double>& data, const LASGlobalData<int>& modmask, const TecRingPar& alpar, LAS::beam_group bg, int print=0);

  void reconstruct_tec_set_onering_weighted(const LASGlobalData<double>& data, const LASGlobalData<double>& weights, TecRingPar& alpar, LAS::beam_group bg, int print);
  void reconstruct_tec_set_onering(const LASGlobalData<double>& data, TecRingPar& alpar, LAS::beam_group bg, int print=0);
  void reconstruct_tec_ring(const LASGlobalData<double>& dphi, const LASGlobalData<int>& mask, TecRingPar& alpar, LAS::beam_group bg, int print=0);

  int set_masked_modules_tec(LASGlobalData<double>& data, const LASGlobalData<int>& modmask, const TecPar& alpar, LAS::beam_group bg, int print=0);
  void reconstruct_tec_set_full(const LASGlobalData<double>& data, TecPar& alpar, LAS::beam_group bg, int print=0);
  bool reconstruct_tec_full(const LASGlobalData<double>& dphi, const LASGlobalData<int>& mask, TecPar& alpar, LAS::beam_group bg, int print=0);
  bool reconstruct_tec_full(const LASGlobalData<double>& dphi, const LASGlobalData<int>& mask, LasAlPar& alpar, int print=0);

  void compute_residuals(const LASGlobalData<double>& data, const LasAlPar& alpar, LASGlobalData<double>& residuals);

 public:
  double r4; // Radius Ring 4
  double r6; // Radius Ring 6
  int nr_iter_tecp; // Nr. of iterations needed for alpar reconstruction in TEC+
  int nr_iter_tecm; // Nr. of iterations needed for alpar reconstruction in TEC-

 private:
  int nbeam; // Number of laser beams
  int ndisc; // Number of endcap discs
  int nring; // Number of rings
  double pitch_r4;  // Pitch Ring 4
  double pitch_r6;  // Pitch Ring 6

  Avec theta; // Phi-positions of beams
  Avec zdat; // z-position of discs
  double L; // Endcap length
  Avec rbeam; // Radius of beams Ring 4/ Ring6

  // Sums over r
  double Sr;
  double Sr2;

  // Sums over z-positions
  double Sz;
  double Sz2;
  
  // Sums over beam positions and their sin/cos
  double St;
  double St2;
  double Ssin;
  double Scos;
  double Ssin2;
  double Scos2;
  double Scossin;

  // Common terms of alignment parameters
  double A;
  double B;
  double C;
  double D;
  double E;
  double F;
  double denom1;
  double denom2;
  double term1;
  double term2;
  double nenn2;
  double znorm;

 public:
  int verbose_flag;
};



#endif
