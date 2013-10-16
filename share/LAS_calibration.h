#ifndef __INCLUDE_LAS_CALIBRATION_H__
#define __INCLUDE_LAS_CALIBRATION_H__

#include <string>
#include <vector>

#include <TGraphErrors.h>

#define __AVECROOT__
#include "Avec.h"

#include "LASGlobalData.h"


enum fibre_swap{ NONE, TECP_R6, TECP_R4, AT, ALL};

void lsboard_scan_analysis(const std::string& data_filename, const std::string& xdaq_log_filename, const std::string& results_filename, const std::string& xml_template_file, const std::string& xml_output_file);

void generate_xml_config(const std::string& results_filename, const std::string& xml_template_file, const std::string& xml_output_file);
void find_intensity_settings(const std::string& output_filename);
void fit_intensity_ramps(const std::string& results_filename);
void lsboard_intensity_ramp(const std::string& results_file);
void lsboard_fit_pulseshapes(const std::string& shapes_file);
void lsboard_create_pulseshapes(const std::string& results_file);
void label_events(const std::string& results_filename);
void make_lsboard_level_labels(const std::string& result_file, const std::string& directory, int board_id, int UTC_hour_offset = 0);
void decode_counter(const std::string& input_file, const std::string& directory="", int subdet=0, int ring=0, fibre_swap swapped_fibres_pattern = TECP_R6, bool plot = false, double threshold = 400);
void slice_run(const std::string& input_file);
void find_scan_boundaries(const std::string& output_filename);
void control_scan_signal(const std::string& resultfile);
void get_scan_signal(const std::string& filename, const std::string& output_rootfilename, Bool_t plot = kFALSE);
void get_eventnr_and_time(const std::string& filename, const std::string& output_rootfilename, Bool_t plot = kFALSE);
void lsboard_get_log_data(const std::string& filename, const std::string& output_filename);

void avec_linreg(const Avec& x, const Avec& y, double& a, double& b);

TGraphErrors* pulseshape_fit(const Avec& xval, const Avec& data, double& xmax, double& ymax, int fit_type = 1, bool output = true);

void cali_pos_control(const std::string& raw_data_filename, const std::string& results_filename, unsigned int board = 2, int det = 3, int ring = -1, int beam = 0, int zpos = 0);

void create_directory_list(const Avec& board_list, std::vector<std::string>& directory_list);

#endif
