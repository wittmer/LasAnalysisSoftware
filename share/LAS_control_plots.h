#ifndef __LAS_CONTROL_PLOTS_H__
#define __LAS_CONTROL_PLOTS_H__

#include <string>
#include <vector>

#include "TMultiGraph.h"
#include "LAS_alpar.h"

void control_profile_maxima_long(const std::string& run_list_name, const std::string& path_prefix = "", int subdet= 3, int beam = 0, int ring = -1, int zpos = 0, int algo_type = 0);
void control_pos_int_correlation(const std::string& run_list_name, const std::string& path_prefix="", LAS::beam_group group = LAS::ALL, int algo_type = 0);
void control_profile_maxima(const std::string& run_list_name, const std::string& path_prefix="", LAS::beam_group group = LAS::ALL, int algo_type = 0);

bool run_selector(const std::string& result_file, double& bad_ratio, Avec::size_type& nr_blocks, double& block_ratio, bool print = false);

double SetDiscParAxes(TMultiGraph* mgr, double min_yrange = 10.0);
void alpar_history_fine(const std::vector<std::string>& file_list);
void alpar_history_draw(const std::vector<LasAlPar>& par_list, const Avec& xvals , bool bs_frame = true, bool process_rings = false);
//void alpar_draw(const LasAlPar& alpar );

void control_pos_av(const std::string& filename);
void control_rms_av(const std::string& filename);
void control_positions(const std::string& results_file);
void control_positions(const std::string& results_file, const std::string& ref_file);
void control_pos_err(const std::string& results_file, const std::string& ref_file);
void control_integrated_pos(const std::string& filename);
void control_profiles(const std::string& results_file, int block_nr = -1);
void control_striprms(const std::string& results_file);
void control_history(const std::string& results_file);
void control_intensities(const std::string& result_filename);
void control_blocks(const std::string& output_filename, bool save_plots);
void control_max_signal(const std::string& max_file, bool full_output = false);

/* // Directives for creation of ROOT Dictionary */
/* #ifdef __CINT__ */


/* #pragma link C++ function control_pos_av(const std::string&); */
/* #pragma link C++ function control_rms_av(const std::string&); */
/* #pragma link C++ function control_positions(const std::string&, const std::string&); */
/* #pragma link C++ function control_positions(const std::string&); */
/* #pragma link C++ function control_integrated_pos(const std::string&); */
/* #pragma link C++ function control_profiles(const std::string&, int); */

/* #pragma link C++ function control_striprms(const std::string&); */
/* #pragma link C++ function control_history(const std::string&); */
/* #pragma link C++ function control_intensities(const std::string&); */
/* #pragma link C++ function control_blocks(const std::string&, bool); */
/* #pragma link C++ function control_max_signal(const std::string&, bool); */

/* #endif */

#endif
