#ifndef __LAS_BASIC_TOOLS_H__
#define __LAS_BASIC_TOOLS_H__

#include <ctime>
#include <vector>
#include <string>
#include <iostream>
#include "TMultiGraph.h"
#include "TMath.h"
#define __AVECROOT__
#include "Avec.h"
#include "LASGlobalDataLoop.h"

#include <stdexcept>

namespace LAS{
  class Exception : public std::runtime_error
  {
  public:
    Exception(const std::string& message);
    virtual ~Exception() throw();
  };

  enum beam_group{ALL, TEC_PLUS, TEC_MINUS, TEC, AT, TIB, TOB, TEC_PLUS_AT, TEC_MINUS_AT, TEC_AT, TEC_PLUS_R4, TEC_PLUS_R6, TEC_MINUS_R4, TEC_MINUS_R6};
  namespace FIT_TYPE{
    enum fit_type{ALL, TEC_PLUS, TEC_MINUS, TEC, AT, TIB, TIB_OLD, TOB, TEC_PLUS_AT, TEC_MINUS_AT, TEC_AT, TEC_PLUS_R4, TEC_PLUS_R6, TEC_MINUS_R4, TEC_MINUS_R6};
  };

  const double r_tec_r4 = 564.0;
  const double r_tec_r6 = 840.0;
  const double r_at_tec = r_tec_r4;
  const double r_at_tib = 514.0;
  const double r_at_tob = 600.0;

  const double tec_zvalues[9] = {0.0, 140.0, 280.0, 420.0, 560.0, 735.0, 925.0, 1130.0, 1345.0};
  const double z_disc0 = 1250.0;
  const double z_tec_bs = 705.0;  // z-position of tec beamsplitters (relative to disc 0, i.e. disc 0 at z = 0)
  const double L_tec = 1345.0;
  const Avec zdat_tec(9, tec_zvalues);

  const Avec theta_tec = 0.3927 + Avec(8,0,7) * TMath::Pi() / 4.0;

  const double dr_tib = r_at_tec - r_at_tib;
  const double tib_zvalues[6] = {620.0, 380.0, 180.0, -100.0, -340.0, -540.0};
  const Avec zdat_tib(6, tib_zvalues);
  const Avec zdat_tib_eff = Avec(6, tib_zvalues) + dr_tib;

  const double dr_tob = r_at_tob - r_at_tec;
  const double tob_zvalues[6] = {1040.0, 580.0, 220.0, -140.0, -500.0, -860.0};
  const Avec zdat_tob(6, tob_zvalues);
  const Avec zdat_tob_eff = Avec(6, tob_zvalues) + dr_tob;

  const double z_at_bs = 1080.0;  // z-position of AT beamsplitters

  //! phi-position of AT beams
  const double at_theta_values[8] = { 0.392699, 1.289799, 1.851794, 2.748894, 3.645995, 4.319690, 5.216791, 5.778784 };
  const Avec theta_at(8, at_theta_values);

  const double cos9deg = 0.987688341;
  const double pitch_tec_r4 = 0.126029;  // Pitch Ring 4
  const double pitch_tec_r6 = 0.189;     // Pitch Ring 6
  const double pitch_at_tec = 0.126029;
  const double pitch_at_tib = 0.120 * cos9deg;
  const double pitch_at_tob = 0.183;

};


//class TMultiGraph;
class TLegend;
class TTree;

void mask_known_bad_modules(LASGlobalData<int>& mask);

std::string get_next_line(std::istream& input);

time_t make_root_time(double year, double month, double day, double hour, double min, double sec);
time_t make_root_time_offset();

void TreeMerge(const std::vector<std::string>& file_list, const std::string& tree_name, const std::string& outfile);
void sort_tree(const std::string& input_filename, const std::string& treename, const std::string& output_filename);
void sort_tree(TTree* tree, const std::string& output_filename);

TLegend* AddLegend(TMultiGraph* mgr, std::vector<std::string>& leg_tit, int ncol=1, float xlo=0.6, float ylo=0.5, float xhi=0.9, float yhi=0.9, bool plot = true);
void SetLineWidth(TMultiGraph* mgr, Width_t lwidth);
void SetXRangeUser(TMultiGraph* mgr, Double_t ufirst, Double_t ulast);

std::string  GetModuleName(const LASGlobalDataLoop& loop, int underscore=0, int prefix=0);
std::string  GetModuleName( int subdetector, int ring, int beam, int zpos , int underscore=0, int prefix=0);
std::string  GetTECName( int subdetector, int ring, int beam, int disc , int underscore=0, int prefix=0);
std::string  GetTIBTOBName( int subdetector, int beam, int pos, int underscore=0 , int prefix=0);
std::string  GetTEC2TECName( int subdetector, int beam, int disc, int underscore=0 , int prefix=0);
int invalid_indices(int subdet, int beam, int zpos, int ring=-1, int type=0, int verbose=0);

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ namespace LAS;
#pragma link C++ namespace LAS::FIT_TYPE;

#pragma link C++ function make_root_time(double, double, double, double, double, double);
#pragma link C++ function sort_tree(const std::string&, const std::string&, const std::string&);
#pragma link C++ function sort_tree(TTree*, const std::string&);
#endif

#endif
