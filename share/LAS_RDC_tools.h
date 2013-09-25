#ifndef __LAS_RDC_TOOLS_H__
#define __LAS_RDC_TOOLS_H__

#include "LASGlobalData.h"
#include "LAS_basic_tools.h"
#include "LAS_alpar.h"
#include <algorithm>


// Class to define all sort of parameters that are used in the reconstruction
namespace LAS{
  class AnalysisParameters
  {
  public:
    AnalysisParameters() : 
      noise_signal_thresh(6000),
      bad_ratio_cut(0.05), 
      noise_algo_type(0),
      edge_threshold(0.5)
	{;}
  public:
    double noise_signal_thresh; 
    double bad_ratio_cut; 
    int noise_algo_type;
    double edge_threshold; 
  };
  
};

double calc_chi2(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, bool verbose = false);
void check_run(const std::string& result_file, const std::string& pos_prefix = "positions_");
void check_run_list(const std::string& run_list_name, const std::string& pos_prefix = "positions_");
void check_run_list(const std::string& ref_filename, const std::vector<std::string>& file_list, const std::string& pos_prefix = "positions_");

void get_file_list(const std::string& filename, std::vector<std::string>& file_list, const std::string& path_prefix = "");
void get_file_list(const std::string& filename, std::string& ref_file, std::vector<std::string>& file_list, const std::string& path_prefix = "");


void check_ref_pos(const std::string& ref_file, int block_nr= -1, double rms_cut = 0.5, const std::string& pos_prefix = "positions_");
bool get_ref_pos(const std::string& ref_file, LASGlobalData<double>& ref_pos, LASGlobalData<double>& ref_err, LASGlobalData<int>& ref_mask, int block_nr=-1, double rms_cut=0.5, const std::string& pos_prefix = "positions_");
void fill_global_tec_data(const LASGlobalData<double>& theData, const LASGlobalData<int>& mask, Avec& tecp_r4, Avec& tecp_r6, Avec& tecm_r4, Avec& tecm_r6);

void Calculate_Difference(const std::string& filename, LASGlobalData<double>& ref_pos, LASGlobalData<int>& ref_mask);
float calc_pos(const std::vector<float>& buffer, std::vector<float>::size_type max_idx, float ratio = 0.5);
void calc_pos(LASGlobalData<Avec>& profiles, LASGlobalData<double>& positions, LASGlobalData<int>& results_mask,   double ratio = 0.5);


void alpar_history(const std::vector<std::string>& file_list, const std::string& ref_file);

void get_pos_mask(const std::string& filename, int block_nr, LASGlobalData<double>& pos, LASGlobalData<int>& mask, int flag = 0, double rms_cut = 1.0);
void at_beam_draw(const std::string& data_file, int block_nr, const std::string& ref_file, int ref_block = -1, int flag = 0, double rms_cut = 1.0);

void Run_History(const std::string& filename, bool plot = false);
void Position_Calculation(const std::string& filename, const LAS::AnalysisParameters& par = LAS::AnalysisParameters());

void Signal_Integration(const std::string& data_file, const std::string& label_file);
void Signal_Integration_old(const std::string& data_file, const std::string& label_file, double edge_threshold = 0.5);

void replace_step_mask(const std::string& mask_filename, const std::string& result_filename);
void Block_Label_Check(const std::string& result_file, Avec::size_type expected_block_size = 2000);
void Event_Labels(const std::string& filename, bool plot = false);
void Block_Slicer(const std::string& output_filename, bool plot=false, double factor = 20, double min_block_size = 1000);
//void Noise_Filter(const std::string& filename, double signal_thresh = 6000, double bad_ratio_cut = 0.05, int algo_type = 0, const LAS::AnalysisParameters& par = LAS::default_parameters);
void Noise_Filter(const std::string& filename, const LAS::AnalysisParameters& par = LAS::AnalysisParameters());
void get_max_signal(const std::string& filename, const std::string& output_rootfilename, bool plot=false, bool full_output = false, Long64_t nentries = -1);

// Directives for creation of ROOT Dictionary
#ifdef __CINT__

#pragma link C++ function check_run(const std::string& result_file);
#pragma link C++ function calc_chi2(const LASGlobalData<double>&, const LASGlobalData<double>&, bool);
#pragma link C++ function check_run_list(const std::string&);
#pragma link C++ function get_file_list(const std::string&, std::vector<std::string>&);
#pragma link C++ function get_file_list(const std::string&, std::string&, std::vector<std::string>&);
#pragma link C++ function check_ref_pos(const std::string&, int, double, const std::string& );
#pragma link C++ function get_ref_pos(const std::string&, LASGlobalData<double>&, LASGlobalData<double>&, LASGlobalData<int>&, int, double, const std::string&);
#pragma link C++ function fill_global_tec_data(const LASGlobalData<double>&, const LASGlobalData<int>&, Avec&, Avec&, Avec&, Avec&);

#pragma link C++ function Calculate_Difference(const std::string&, LASGlobalData<double>&, LASGlobalData<int>&);
#pragma link C++ function calc_pos(LASGlobalData<Avec>&, LASGlobalData<double>&, LASGlobalData<int>&,   double);

#pragma link C++ function alpar_history(const std::vector<std::string>&, const std::string&);

#pragma link C++ function get_pos_mask(const std::string&, int, LASGlobalData<double>&, LASGlobalData<int>&, int, double);
#pragma link C++ function at_beam_draw(const std::string&, int, const std::string&, int, int, double);

#pragma link C++ function Run_History(const std::string&, bool);
#pragma link C++ function Position_Calculation(const std::string&, const LAS::AnalysisParameters&);
#pragma link C++ function Signal_Integration(const std::string&, const std::string&);
#pragma link C++ function replace_step_mask(const std::string&, const std::string&);
#pragma link C++ function Block_Label_Check(const std::string&, Avec::size_type);
#pragma link C++ function Event_Labels(const std::string&, bool);
#pragma link C++ function Block_Slicer(const std::string&, bool, double, double);
#pragma link C++ function Noise_Filter(const std::string&, const LAS::AnalysisParameters&);
#pragma link C++ function get_max_signal(const std::string&, const std::string&, bool, bool, Long64_t);

#endif


template <typename T>
double compute_position(T it_first, T it_last, double ratio = 0.5)
{
  T max_it = std::max_element(it_first, it_last);
  
  if(*max_it == 0) return -1;
  
  double threshold = *max_it * ratio;
  
  T  left_it = max_it;
  T right_it = max_it;
  
  while( *left_it > threshold &&  left_it != it_first   )  left_it--;
  while(*right_it > threshold && right_it != it_last - 1 ) right_it++;
  
//   //class T::size_type  left_idx =  left_it - data.begin();
//   //class T::size_type right_idx = right_it - data.begin();
  
  double  left_xt = (threshold - *( left_it  )) / (*( left_it + 1) - *( left_it    )) + (double)( left_it - it_first)    ;
  double right_xt = (threshold - *(right_it-1)) / (*(right_it    ) - *(right_it - 1)) + (double)(right_it - it_first - 1);
     
  return (left_xt + right_xt)/2;
}

template <class T>
double compute_position(const T& data, double ratio = 0.5)
{
  return compute_position(data.begin(), data.end(), ratio);
/*   class T::const_iterator max_it = std::max_element(data.begin(), data.end()); */
  
/*   if(*max_it == 0) return -1; */
  
/*   double threshold = *max_it * ratio; */
  
/*   class T::const_iterator left_it = max_it; */
/*   class T::const_iterator right_it = max_it; */
  
/*   while( *left_it > threshold &&  left_it != data.begin()   )  left_it--; */
/*   while(*right_it > threshold && right_it != data.end() - 1 ) right_it++; */
  
/*   class T::size_type  left_idx =  left_it - data.begin(); */
/*   class T::size_type right_idx = right_it - data.begin(); */
  
/*   double  left_xt = (threshold - *( left_it  )) / (*( left_it + 1) - *( left_it    )) +  left_idx    ; */
/*   double right_xt = (threshold - *(right_it-1)) / (*(right_it    ) - *(right_it - 1)) + right_idx - 1; */
     
/*   return (left_xt + right_xt)/2; */
}

#endif
