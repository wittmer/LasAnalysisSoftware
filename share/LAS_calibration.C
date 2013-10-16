#include "LAS_calibration.h"

#include <iomanip>

#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TAxis.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>

#include "LASGlobalDataLoop.h"
#include "LAS_basic_tools.h"
#include "LAS_globaldata_tools.h"
#include "LAS_vectorfloat_tools.h"
#include "LAS_RDC_tools.h"

#include "Avec2D.h"

void lsboard_scan_analysis(const std::string& data_filename, const std::string& xdaq_log_filename, const std::string& results_filename, const std::string& xml_template_file, const std::string& xml_output_file)
{
  gSystem->Exec( ("../scripts/extract_lsboard_scan.sh " + xdaq_log_filename + " " + "xdaq_log").c_str() );

  Avec::VERBOSE_FLAG = 0;

  TFile outf(results_filename.c_str(), "RECREATE");
  outf.mkdir("Control_Plots");
  outf.mkdir("Details");
  outf.Close();

  lsboard_get_log_data("xdaq_log_values.dat", results_filename);
  
  get_eventnr_and_time(data_filename, results_filename);
  
  find_scan_boundaries(results_filename);

  get_scan_signal(data_filename, results_filename);
  
  //slice_run(results_filename);
  
  label_events(results_filename);
  
  lsboard_create_pulseshapes(results_filename);
  
  lsboard_fit_pulseshapes(results_filename);
  
  lsboard_intensity_ramp(results_filename);

  find_intensity_settings(results_filename);

  generate_xml_config(results_filename, xml_template_file, xml_output_file);

}

void generate_xml_config(const std::string& results_filename, const std::string& xml_template_file, const std::string& xml_output_file)
{
  TFile f(results_filename.c_str(), "READ");

  Avec2D settings_board_1 = avec2d_get("intensity_setting_board_1", f);
  Avec2D settings_board_2 = avec2d_get("intensity_setting_board_2", f);
  Avec2D settings_board_3 = avec2d_get("intensity_setting_board_3", f);
  Avec2D settings_board_4 = avec2d_get("intensity_setting_board_4", f);
  Avec2D settings_board_5 = avec2d_get("intensity_setting_board_5", f);

  Avec2D delay_board_1 = avec2d_get("delay_setting_board_1", f);
  Avec2D delay_board_2 = avec2d_get("delay_setting_board_2", f);
  Avec2D delay_board_3 = avec2d_get("delay_setting_board_3", f);
  Avec2D delay_board_4 = avec2d_get("delay_setting_board_4", f);
  Avec2D delay_board_5 = avec2d_get("delay_setting_board_5", f);

  Avec2D norm_board_1 = avec2d_get("norm_board_1", f);
  Avec2D norm_board_2 = avec2d_get("norm_board_2", f);
  Avec2D norm_board_3 = avec2d_get("norm_board_3", f);
  Avec2D norm_board_4 = avec2d_get("norm_board_4", f);
  Avec2D norm_board_5 = avec2d_get("norm_board_5", f);

  std::ofstream sf("sed_file.txt");
  for(unsigned int diode = 0; diode < 8; diode++){
    for(unsigned int level = 0; level < 5; level++){
      sf << "s/__INT_1_" << diode+1 << "_" << level + 1 << "__/" << settings_board_1[level][diode] << "/\n";
      sf << "s/__INT_2_" << diode+1 << "_" << level + 1 << "__/" << settings_board_2[level][diode] << "/\n";
      sf << "s/__INT_3_" << diode+1 << "_" << level + 1 << "__/" << settings_board_3[level][diode] << "/\n";
      sf << "s/__INT_4_" << diode+1 << "_" << level + 1 << "__/" << settings_board_4[level][diode] << "/\n";
      sf << "s/__INT_5_" << diode+1 << "_" << level + 1 << "__/" << settings_board_5[level][diode] << "/\n";

      sf << "s/__DEL_1_" << diode+1 << "_" << level + 1 << "__/" <<  (norm_board_1[level][diode] ? (int)delay_board_1[level][diode] : 0) << "/\n"; 
      sf << "s/__DEL_2_" << diode+1 << "_" << level + 1 << "__/" <<  (norm_board_2[level][diode] ? (int)delay_board_2[level][diode] : 0) << "/\n"; 
      sf << "s/__DEL_3_" << diode+1 << "_" << level + 1 << "__/" <<  (norm_board_3[level][diode] ? (int)delay_board_3[level][diode] : 0) << "/\n"; 
      sf << "s/__DEL_4_" << diode+1 << "_" << level + 1 << "__/" <<  (norm_board_4[level][diode] ? (int)delay_board_4[level][diode] : 0) << "/\n"; 
      sf << "s/__DEL_5_" << diode+1 << "_" << level + 1 << "__/" <<  (norm_board_5[level][diode] ? (int)delay_board_5[level][diode] : 0) << "/\n"; 
    }
  }
  sf.close();

  gSystem->Exec( ("sed -f sed_file.txt <" + xml_template_file + " >" + xml_output_file).c_str() );

}

void find_intensity_settings(const std::string& output_filename)
{

  TFile f(output_filename.c_str(), "UPDATE");

  LASGlobalData<Avec>& x_vals = global_data_get<Avec>("x_vals", f);
  LASGlobalData<Avec>& y_vals = global_data_get<Avec>("y_vals", f);

  LASGlobalData<int>& bad_modules = global_data_get<int>("bad_modules", f);
  LASGlobalData<double>& optimum_delay = global_data_get<double>("optimum_delay", f);

  Avec2D settings_board_1(5,0);
  Avec2D settings_board_2(5,0);
  Avec2D settings_board_3(5,0);
  Avec2D settings_board_4(5,0);
  Avec2D settings_board_5(5,0);

  Avec2D delay_board_1(5,8,0);
  Avec2D delay_board_2(5,8,0);
  Avec2D delay_board_3(5,8,0);
  Avec2D delay_board_4(5,8,0);
  Avec2D delay_board_5(5,8,0);

  Avec2D norm_board_1(5,8,0);
  Avec2D norm_board_2(5,8,0);
  Avec2D norm_board_3(5,8,0);
  Avec2D norm_board_4(5,8,0);
  Avec2D norm_board_5(5,8,0);

  // Find upper and lower limit of intensity setting for each module
  LASGlobalData<double> lo_val;
  LASGlobalData<double> hi_val;
  double lower_limit = 40;
  double upper_limit = 200;

  LASGlobalDataLoop lp;
  do{
    if(lp.GetEntry(bad_modules) != 0) continue;
    Avec& x = lp.GetEntry(x_vals);
    Avec& y = lp.GetEntry(y_vals);
    Avec local_mask = y >= lower_limit && y <= upper_limit;
    Avec range = avec_mask(x, local_mask);
    if(range.empty()){
      lp.GetEntry(bad_modules) = 1;
      continue;
    }
    lp.GetEntry(lo_val) = vmin(range);
    lp.GetEntry(hi_val) = vmax(range);
  }while(lp.next());

  std::cout << "Bad Modules:" << std::endl;
  print_global_data(bad_modules, LAS::AT);
  std::cout << "Low Values:" << std::endl;
  print_global_data(lo_val, LAS::AT);
  std::cout << "High Values:" << std::endl;
  print_global_data(hi_val, LAS::AT);


  LASGlobalData<int> flagged = bad_modules;
  LASGlobalData<double> av_delay;
  LASGlobalData<int> step_mask(0);

  for(int level = 5; level > 0; level--){
    std::cout << "Level " << level << std::endl;

    /////////////
    // TEC- R4 //
    /////////////

    Avec largest_lo(8, 0.0);
    lp = LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS_R4);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(largest_lo[beam] < lp.GetEntry(lo_val)) largest_lo[beam] = lp.GetEntry(lo_val); 
    }while(lp.next());

    std::cout << "Largest lo_vals for TEC- R4 beams: \n" << largest_lo << std::endl;

    lp = LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS_R4);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(lp.GetEntry(hi_val) >= largest_lo[beam] - 0.5){
	lp.GetEntry(flagged) = 1;
	delay_board_1[level-1][beam] += lp.GetEntry(optimum_delay);
	norm_board_1[level-1][beam]++;
	lp.GetEntry(step_mask) = level;
	//std::cout << "beam " << beam << " det " << lp.get_det() << " zpos " << lp.get_zpos() << std::endl;
      }
    }while(lp.next());
    settings_board_1[level-1] = largest_lo;

    ////////////////////////////////////////////////////////////////////////////////////

    /////////////
    // TEC- R6 //
    /////////////

    largest_lo = Avec(8, 0.0);
    lp = LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS_R6);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(largest_lo[beam] < lp.GetEntry(lo_val)) largest_lo[beam] = lp.GetEntry(lo_val); 
    }while(lp.next());

    std::cout << "Largest lo_vals for TEC- R6 beams: \n" << largest_lo << std::endl;

    lp = LASGlobalDataLoop(LASGlobalDataLoop::TEC_MINUS_R6);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(lp.GetEntry(hi_val) >= largest_lo[beam] - 0.5){
	lp.GetEntry(flagged) = 1;
	delay_board_2[level-1][beam] += lp.GetEntry(optimum_delay);
	norm_board_2[level-1][beam]++;
	lp.GetEntry(step_mask) = level;
	//std::cout << "beam " << beam << " det " << lp.get_det() << " zpos " << lp.get_zpos() << std::endl;
      }
    }while(lp.next());
    settings_board_2[level-1] = largest_lo;

    ////////////////////////////////////////////////////////////////////////////////////

    ////////
    // AT //
    ////////

    largest_lo = Avec(8, 0.0);
    lp = LASGlobalDataLoop(LASGlobalDataLoop::AT);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(largest_lo[beam] < lp.GetEntry(lo_val)) largest_lo[beam] = lp.GetEntry(lo_val); 
    }while(lp.next());

    std::cout << "Largest lo_vals for AT beams: \n" << largest_lo << std::endl;

    lp = LASGlobalDataLoop(LASGlobalDataLoop::AT);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(lp.GetEntry(hi_val) >= largest_lo[beam] - 0.5){
	lp.GetEntry(flagged) = 1;
	delay_board_3[level-1][beam] += lp.GetEntry(optimum_delay);
	norm_board_3[level-1][beam]++;
	lp.GetEntry(step_mask) = level;
	//std::cout << "beam " << beam << " det " << lp.get_det() << " zpos " << lp.get_zpos() << std::endl;
      }
    }while(lp.next());
    settings_board_3[level-1] = largest_lo;

    ////////////////////////////////////////////////////////////////////////////////////

    /////////////
    // TEC+ R4 //
    /////////////

    largest_lo = Avec(8, 0.0);
    lp = LASGlobalDataLoop(LASGlobalDataLoop::TEC_PLUS_R4);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(largest_lo[beam] < lp.GetEntry(lo_val)) largest_lo[beam] = lp.GetEntry(lo_val); 
    }while(lp.next());

    std::cout << "Largest lo_vals for TEC+ R4 beams: \n" << largest_lo << std::endl;

    lp = LASGlobalDataLoop(LASGlobalDataLoop::TEC_PLUS_R4);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(lp.GetEntry(hi_val) >= largest_lo[beam] - 0.5){
	lp.GetEntry(flagged) = 1;
	delay_board_4[level-1][beam] += lp.GetEntry(optimum_delay);
	norm_board_4[level-1][beam]++;
	lp.GetEntry(step_mask) = level;
	//std::cout << "beam " << beam << " det " << lp.get_det() << " zpos " << lp.get_zpos() << std::endl;
      }
    }while(lp.next());
    settings_board_4[level-1] = largest_lo;

    ////////////////////////////////////////////////////////////////////////////////////

    /////////////
    // TEC+ R6 //
    /////////////

    largest_lo = Avec(8, 0.0);
    lp = LASGlobalDataLoop(LASGlobalDataLoop::TEC_PLUS_R6);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(largest_lo[beam] < lp.GetEntry(lo_val)) largest_lo[beam] = lp.GetEntry(lo_val); 
    }while(lp.next());

    std::cout << "Largest lo_vals for TEC+ R6 beams: \n" << largest_lo << std::endl;

    lp = LASGlobalDataLoop(LASGlobalDataLoop::TEC_PLUS_R6);
    do{
      if(lp.GetEntry(flagged)) continue;
      int beam = lp.get_beam();
      if(lp.GetEntry(hi_val) >= largest_lo[beam] - 0.5){
	lp.GetEntry(flagged) = 1;
	delay_board_5[level-1][beam] += lp.GetEntry(optimum_delay);
	norm_board_5[level-1][beam]++;
	lp.GetEntry(step_mask) = level;
	//std::cout << "beam " << beam << " det " << lp.get_det() << " zpos " << lp.get_zpos() << std::endl;
      }
    }while(lp.next());
    settings_board_5[level-1] = largest_lo;

  }

  // Average over all delays for one setting
  delay_board_1 /= norm_board_1;
  delay_board_2 /= norm_board_2;
  delay_board_3 /= norm_board_3;
  delay_board_4 /= norm_board_4;
  delay_board_5 /= norm_board_5;

  avec_remove_nan(delay_board_1);
  avec_remove_nan(delay_board_2);
  avec_remove_nan(delay_board_3);
  avec_remove_nan(delay_board_4);
  avec_remove_nan(delay_board_5);

  settings_board_1.Write("intensity_setting_board_1");
  settings_board_2.Write("intensity_setting_board_2");
  settings_board_3.Write("intensity_setting_board_3");
  settings_board_4.Write("intensity_setting_board_4");
  settings_board_5.Write("intensity_setting_board_5");

  delay_board_1.Write("delay_setting_board_1");
  delay_board_2.Write("delay_setting_board_2");
  delay_board_3.Write("delay_setting_board_3");
  delay_board_4.Write("delay_setting_board_4");
  delay_board_5.Write("delay_setting_board_5");

  norm_board_1.Write("norm_board_1");
  norm_board_2.Write("norm_board_2");
  norm_board_3.Write("norm_board_3");
  norm_board_4.Write("norm_board_4");
  norm_board_5.Write("norm_board_5");

  step_mask.Write("step_mask");

  f.Close();

  std::cout << "Flagged or bad modules:\n";
  print_global_data(flagged);

  std::cout << "Step Mask:\n";
  print_global_data(step_mask);

  delete &bad_modules;
  delete &optimum_delay;
  delete &x_vals;
  delete &y_vals;
}

// Linear regression through points (x,y), no errors, no weighting
void avec_linreg(const Avec& x, const Avec& y, double& a, double& b)
{
  double N = x.size();
  a = (vsum(x)*vsum(y)/N - vsum(x*y))/(vsum(x)*vsum(x)/N - vsum(x*x));
  b = (vsum(y) - a * vsum(x))/N;
}

//! Fit straight lines through the intensity ramps
void fit_intensity_ramps(const std::string& results_filename)
{
  TFile f(results_filename.c_str(), "UPDATE");

  LASGlobalData<Avec>& x_vals = global_data_get<Avec>("x_vals", f);
  LASGlobalData<Avec>& y_vals = global_data_get<Avec>("y_vals", f);

  LASGlobalData<int>& bad_modules = global_data_get<int>("bad_modules", f);

  // Accepted signal range
  double lower_limit = 20;
  double upper_limit = 200;

  LASGlobalData<double> slope(0);
  LASGlobalData<double> offset(0);

  LASGlobalDataLoop lp;
  do{
    if(lp.GetEntry(bad_modules) != 0) continue;
    Avec& x = lp.GetEntry(x_vals);
    Avec& y = lp.GetEntry(y_vals);
    Avec local_mask = y >= lower_limit && y <= upper_limit;
    Avec x_range = avec_mask(x, local_mask);
    Avec y_range = avec_mask(y, local_mask);
    if(x_range.empty()){
      lp.GetEntry(bad_modules) = 1;
      continue;
    }
    avec_linreg(y_range, x_range, lp.GetEntry(slope), lp.GetEntry(offset));

  }while(lp.next());

  offset.Write("intensity_ramp_offset");
  slope.Write("intensity_ramp_slope");
  f.Close();
}

//! The intensity labels in the scan correspond to specific register settings of the board
Avec convert_label_to_intensity(const Avec& data)
{
  // The first 40 intensity steps range from 900 to 1290 in steps of 10
  Avec fine(40, 900, 1290);
  // The last 27 intensity steps range from 1300 to 2600 in steps of 50
  Avec coarse(27, 1300, 2600);
  Avec intensity_val = Avec(1,0.0) & fine & coarse;
  //std::cout << "intensity_val of size " << intensity_val.size() << std::endl << intensity_val << std::endl;

  // Now convert the input Avec with intensity labels to register settings
  Avec retval(data.size());
  for(unsigned int i = 0; i < data.size(); i++){
    unsigned int idx = data[i];
    if(idx < intensity_val.size()) retval[i] = intensity_val[idx];
  }
  return retval;
}


//! Map the 5 Laser board Ids to the 5 beam groups
LASGlobalDataLoop::loop_type convert_boardid_to_loop( int boardid)
{
  LASGlobalDataLoop::loop_type loop_t;
  switch(boardid){
  case 1:
    loop_t = LASGlobalDataLoop::TEC_MINUS_R4;
    break;
  case 2:
    loop_t = LASGlobalDataLoop::TEC_MINUS_R6;
    break;
  case 3:
    loop_t = LASGlobalDataLoop::AT;
    break;
  case 4:
    loop_t = LASGlobalDataLoop::TEC_PLUS_R4;
    break;
  case 5:
    loop_t = LASGlobalDataLoop::TEC_PLUS_R6;
    break;
  default:
    std::ostringstream err_msg;
    err_msg << "Unrecognized board id (" << boardid << ")";
    throw LAS::Exception(err_msg.str());
  }
  return loop_t;
}

// Sort x in ascending order and put y in the same order
// Implemented with a bubble-sort algorithm
void sort_avecs(Avec& x, Avec& y)
{
  if(x.size() != y.size()) throw LAS::Exception("Error in sort_avecs: x and y do not have the same size");
  if( x.size() == 1 ) return;

  bool flag = false;

  do{
    flag = false;
    for(Avec::size_type idx = 0; idx < x.size()-1; idx++){
      if( x[idx] > x[idx+1]){
	std::swap( x[idx], x[idx+1]);
	std::swap( y[idx], y[idx+1]);
	flag = true;
      }
    }
  }while (flag);

}

//! Use the best delay from lsboard_fit_pulseshapes to determine the intensity ramp
void lsboard_intensity_ramp(const std::string& results_file)
{
  try{
    TFile f(results_file.c_str(), "UPDATE");

    LASGlobalData<Avec> x_vals; 
    LASGlobalData<Avec> y_vals;

    Avec board_list = avec_get("board_list", f);
    std::vector<std::string> directory_list;
    create_directory_list(board_list, directory_list);


    for(unsigned int i = 0; i < board_list.size(); i++){ // loop over all boards
      int board_id = board_list[i];
      std::string& directory = directory_list[i];
      std::cout << "Creating intensity ramp for board " << board_id << std::endl;

      Avec label_board = avec_get(directory + "/label_board", f);
      Avec label_level = avec_get(directory + "/label_level", f);
      Avec label_fine = avec_get(directory + "/fine_label", f);
      if(label_board.empty() || label_level.empty() || label_fine.empty())
	throw LAS::Exception("Could not find all Avecs with labels");
      
      LASGlobalData<double>& optimum_delay = global_data_get<double>("optimum_delay", f);
      LASGlobalData<int>& bad_modules = global_data_get<int>("bad_modules", f);

      // Object with signal maxima
      LASGlobalData<Avec>* max_signal=0;

      f.GetObject((directory + "/max_signal").c_str(), max_signal);
      if(!max_signal) throw LAS::Exception("Did not find max_signal object");
  
      LASGlobalDataLoop::loop_type loop_t = convert_boardid_to_loop(board_id);
      LASGlobalDataLoop loop(loop_t);
      for(;!loop.finished();loop.next()){ // loop over all modules of the board
	if(loop.GetEntry(bad_modules) != 0) continue;

	double delay = loop.GetEntry(optimum_delay);

	Avec& buffer = loop.GetEntry<Avec>(*max_signal);

	Avec mask = label_board == board_id && label_fine <= delay + 0.5 && label_fine > delay - 0.5;

	Avec intensity_values = convert_label_to_intensity(avec_mask(label_level, mask));
	Avec signal_values = avec_mask(buffer, mask);
	sort_avecs(intensity_values, signal_values);
	Avec::iterator it_x = intensity_values.begin();
	Avec::iterator it_y = signal_values.begin();
	for(;it_x != intensity_values.end(); it_x++, it_y++)if(*it_y >= 220)break;
	intensity_values.erase(it_x, intensity_values.end());
	signal_values.erase(it_y, signal_values.end());
	loop.GetEntry(x_vals) = intensity_values;
	loop.GetEntry(y_vals) = signal_values;
      }
    }

    x_vals.Write("x_vals");
    y_vals.Write("y_vals");
    
    f.Close();

  }
  catch(LAS::Exception e){
    std::cerr << "Caught LAS::Exception saying: " << e.what() << std::endl;
  }

}

Double_t pulse_function(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  if(xx > par[2]) return 0;
  Double_t f = par[0]*(par[2]-xx)*exp(par[1]*(xx-par[2]));
  return f;
}

TF1* pulse_function()
{
  TF1 *f1 = new TF1("pulse_func",pulse_function,0,100,3);
  f1->SetParameters(1, 1, 1);
  f1->SetParNames("constant","exponent","start");
  return f1;
}


TGraphErrors* pulseshape_fit(const Avec& xval, const Avec& data, double& xmax, double& ymax, int fit_type, bool output)
{
  if(xval.empty() || data.empty())return 0;

  Avec::size_type idx = lvmax(data);
  ymax = data[idx];
  xmax = xval[idx];

//   TF1 *f1 = new TF1("pulseshape","[0]*([2]-x)*exp([1]*(x-[2]))",0,37);
//   f1->SetParameter(0, 1e3);
//   f1->SetParameter(1, 0.2);
//   f1->SetParameter(2, 37);

  TF1* f1 = pulse_function();
  double est_max = vmax(data);
  //double est_sta = vmax(xval);
  //double est_sta = (vmax(xval) + xval[lvmax(data)])/2;
  double est_sta = xval[lvmax(data)] + 2;

  double est_exp = 1/((est_sta - xval[lvmax(data)]));
  if(output){
    std::cout << "est_max = " << est_max << std::endl;
    std::cout << "est_exp = " << est_exp << std::endl;
    std::cout << "est_sta = " << est_sta << std::endl;
  }
  f1->SetParameter(0,est_max);
  f1->SetParameter(1,est_exp);
  f1->SetParameter(2,est_sta);

  TGraphErrors* gre;
  double standard_noise = 0.04;

  Avec::size_type xsize = xval.size();
  Avec::VERBOSE_FLAG = 0;
  gre = avec_draw(xval, data, Avec(xsize), Avec(xsize,standard_noise), "Pulsefit", "Entry", "max_signal", "AP", (output ? 1 : 0) );
  Avec::VERBOSE_FLAG = 1;
  double hymax = vmax(data)*1.05;
  //gre->GetYaxis()->SetRangeUser(0,hymax);

  Int_t fit_result = -1;
  if(fit_type == 2) fit_result = gre->Fit("gaus", "Q0");
  else fit_result = gre->Fit("pulse_func", "Q0");
  
  double p0 = f1->GetParameter(0);
  double p1 = f1->GetParameter(1);
  double p2 = f1->GetParameter(2);
  
  if(fit_type == 2){
    TF1 *fit_dat = gre->GetFunction("gaus");
    
    xmax = fit_dat->GetParameter(1);
    ymax = fit_dat->GetParameter(0);
  }
  else{
    xmax = p2-1/p1;
    ymax = f1->Eval(xmax, 0, 0);
  }
  if(output){
    std::cout << "p0: " << p0 << std::endl;
    std::cout << "p1: " << p1 << std::endl;
    std::cout << "p2: " << p2 << std::endl;
    std::cout << "p2 - xmax: " << (p2 - xmax) << std::endl;
    std::cout << "chi2: " << f1->GetChisquare() << std::endl;
    std::cout << "fit_result: " << fit_result << std::endl;
    std::cout << "xmax: " << xmax << "  ymax: " << ymax << std::endl;
    TLine* l = new TLine(xmax,0,xmax,hymax);
    l->Draw();
  }
  return gre;
}

//! Compute the average delay of TEC modules as function of z-position
Avec average_delay(const LASGlobalData<double>& delay)
{
  Avec average_delay(9,0.);
  Avec norm(9, 0.);
  LASGlobalDataLoop loop(LASGlobalDataLoop::TEC);

  do{
    if(loop.GetEntry(delay) > 20){
      average_delay[loop.get_zpos()] += loop.GetEntry(delay);
      norm[loop.get_zpos()] ++;
    }
  }while(loop.next());
  average_delay /= norm;
  avec_remove_nan(average_delay);

  // new TCanvas("av_tec_delay", "average TEC delay");
  //avec_draw(average_delay, "Average TEC delay", "Disc Nr.", "delay [ns]", "AP");

  return average_delay;
}

//! Fit the pulseshapes created with lsboard_create_pulseshapes
void lsboard_fit_pulseshapes(const std::string& results_file)
{
  const unsigned int shape_size = 100;
  LASGlobalData<Avec>& shapes = global_data_get<Avec>("pulse_shapes", results_file);

  unsigned int range = 15;

  LASGlobalData<double> delay(0);
  LASGlobalData<int> bad_modules(0);

  LASGlobalDataLoop loop;
  do{
    Avec& buffer = loop.GetEntry(shapes);
    if(buffer.size() != shape_size){
      loop.GetEntry(bad_modules) = 1;
      continue;
    }
    unsigned int max_pos = lvmax(buffer);
    if(max_pos == 0 ){
      continue;
      loop.GetEntry(bad_modules) = 1;
    }
    unsigned int low_pos = (max_pos < range ? 0 : max_pos - range);
    unsigned int hig_pos = (max_pos + range >= buffer.size() ? buffer.size() - 1 : max_pos + range);
    //avec_draw(Avec(2*range+1, low_pos, hig_pos), buffer(low_pos, hig_pos), "", "", "", "AP");
    double xmax, ymax;
    pulseshape_fit(Avec(hig_pos - low_pos + 1, low_pos, hig_pos), buffer(low_pos, hig_pos), xmax, ymax, 2, false);
    if(xmax < 1.0){
      loop.GetEntry(bad_modules) = 1;
      continue;
    }
    loop.GetEntry(delay) = xmax;

  }while(loop.next());
  delete &shapes;

  TFile of(results_file.c_str(), "UPDATE");

  Bool_t batch_mode = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  bad_modules.Write("bad_modules");
  delay.Write("optimum_delay");

  // draw_global_data(delay, bad_modules * -1 + 1);

  average_delay(delay).Write("average_TEC_delay");
  // std::cout << "List of modules with no signal:" << std::endl;
  // print_global_data(bad_modules);

  gROOT->SetBatch(batch_mode);

  of.Close();

}

//! Create pulseshapes for all modules
void lsboard_create_pulseshapes(const std::string& results_file)
{
  TFile f(results_file.c_str(),"UPDATE");
  
  Avec board_list = avec_get("board_list", f);
  std::vector<std::string> directory_list;
  create_directory_list(board_list, directory_list);
  
  LASGlobalData<Avec> shapes_y;
  LASGlobalData<Avec> shapes_x;
  LASGlobalData<int> shape_norm(0);
  
  for(unsigned int i = 0; i < board_list.size(); i++){
    int board_id = board_list[i];
    std::string& directory = directory_list[i];
    
    std::cout << "Creating pulseshape for board " << board_id << std::endl;
    LASGlobalDataLoop::loop_type loop_t = convert_boardid_to_loop(board_id);
    
    Avec label_board = avec_get(directory + "/label_board", f);
    Avec label_level = avec_get(directory + "/label_level", f);
    Avec label_fine = avec_get(directory + "/fine_label", f);
    if(label_board.empty() || label_level.empty() || label_fine.empty())
      throw LAS::Exception("Could not find all Avecs with labels");
    
    // Object with signal maxima all modules, all settings
    LASGlobalData<Avec>* max_signal=0;
    f.GetObject( (directory + "/max_signal").c_str(), max_signal);
    if(!max_signal) throw LAS::Exception("Did not find max_signal object");

    int min_level = vmin(label_level);
    int max_level = vmax(label_level);
    
    double thresh_max_signal = 254;
    double thresh_min_signal = 20;
    
    // Loop over all modules
    LASGlobalDataLoop loop(loop_t);
    for(;!loop.finished();loop.next()){
      //std::cout << "Next module" << std::endl;
      Avec& buffer = loop.GetEntry<Avec>(*max_signal); // signal maxima for one module, all settings
      Avec local_shape;
      int norm = 0;

      for(int level = min_level; level <= max_level; level++){ // loop over all intensity settings
	//std::cout << "Level " << level << std::endl;	
	Avec mask = label_board == board_id && label_level == level;
	double maximum = vmax( mask * buffer );
	if( maximum >= thresh_max_signal || maximum <= thresh_min_signal) continue;

	Avec extracted_y;
	Avec extracted_x;
	for(Avec::size_type i = 0 ; i < mask.size(); i++){
	  if(mask[i]){
	    extracted_y.push_back( buffer[i] );
	    extracted_x.push_back( label_fine[i] );
	  }
	}
	sort_avecs(extracted_x, extracted_y);
	if( local_shape.empty()) local_shape = extracted_y;
	else local_shape += extracted_y;
	norm ++;
      }

      // Divide by the number of curves that were added up
      if (norm > 0){
	local_shape /= norm;
	if(vmax(local_shape) > 10){
	  loop.GetEntry<Avec>(shapes_y) = local_shape;
	  loop.GetEntry<int>(shape_norm) = norm;
	}
      }
    }
  }

  shapes_y.Write("pulse_shapes");
  f.Close();
}

//! Put intensity and delay labels on each event
void label_events(const std::string& results_filename)
{

  Avec board_list = avec_get("board_list", results_filename);
  std::vector<std::string> directory_list;
  create_directory_list(board_list, directory_list);


  for(unsigned int i = 0; i < board_list.size(); i++){
    int board_id = board_list[i];
    std::cout << "Creating labels for board " << board_id << std::endl;
    int subdet = 0; int ring = 0; fibre_swap swapped_fibres_pattern = TECP_R6;
    switch(board_id){
    case 3:
    case 1:
    case 4:
      subdet = 1; ring = 1; swapped_fibres_pattern = TECP_R6;
      break;
    case 2:
    case 5:
      subdet = 3; ring = -1; swapped_fibres_pattern = AT;
      break;
    default:
      std::cerr << "Error in lsboard_scan_analysis: Unknown board id " << board_id;
      return;
    }
    std::cout << "Creating delay labels with subdet = " << subdet << "  ring = " << ring << std::endl;
    decode_counter(results_filename, directory_list[i], subdet, ring, swapped_fibres_pattern);
    std::cout << "Creating level labels" << std::endl;
    make_lsboard_level_labels(results_filename, directory_list[i], board_list[i], 0);
  }

}

void control_lsboard_labels(const std::string& max_file, const std::string& label_file)
{
  Avec label_board, label_level;
  vecread( ( label_board | label_level), label_file);

  Avec eventnumber, unixTime, signals_tecplus, signals_tecminus, signals_at;
  vecread((eventnumber | unixTime | signals_tecplus | signals_tecminus | signals_at), max_file);

  new TCanvas("control_lsboard_labels","LsBoard Labels");
  avec_draw(eventnumber, (label_board | label_level), "LsBoard Labels", "event nr.", "board id", "AP");

}

// Find the index in Avec 'time' that is compatible with the timestamp
Avec::size_type find_setting_start(const Avec& time, const Avec& level, double timestamp)
{
  for(Avec::size_type f = time.size(); f != 0; f--) if(time[f-1] <= timestamp && level[f-1] >= 0) return f-1;
  throw LAS::Exception("Error in find_setting_start: Could not find setting");
  return 0;
}

//! Create labels with the board inensity settings
void make_lsboard_level_labels(const std::string& result_file, const std::string& directory, int board_id, int UTC_hour_offset)
{
  TFile f(result_file.c_str(), "UPDATE");

  if(!f.cd(directory.c_str()))
    throw LAS::Exception( "Error in make_lsboard_level_labels, file " + result_file + " has no directory named " + directory );

  Avec fine_label = avec_get(directory + "/fine_label", f);
  Avec unixTime = avec_get(directory + "/unixTime", f);
  Avec eventnumber = avec_get(directory + "/eventnumber", f);
  Avec settings_board = avec_get("lsboard_boardid", f);
  Avec settings_level = avec_get("lsboard_level", f);
  Avec settings_time = avec_get("lsboard_set_time", f);
  if(fine_label.empty() || unixTime.empty() || eventnumber.empty() || settings_board.empty() || settings_level.empty() || settings_time.empty() )
    throw LAS::Exception( "Error in  make_lsboard_level_labels: Not all required Avecs were found" );

  for(Avec::size_type idx = 1; idx < settings_time.size(); idx++){
    if( settings_time[idx] == settings_time[idx-1] && settings_level[idx] > 0 && settings_level[idx-1] > 0 ){
      std::cout << "Warning in make_lsboard_level_labels: There are two levels withthe same rimestamp:" << std::endl;
      std::cout << "board " << settings_board[idx-1] << ", level " << settings_level[idx-1] << ", timestamp " << std::setprecision(12) << settings_time[idx-1] << std::endl;
      std::cout << "board " << settings_board[idx] << ", level " << settings_level[idx] << ", timestamp " << std::setprecision(12) << settings_time[idx] << std::endl;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Avec t1 = settings_time;
  unixTime += UTC_hour_offset * 3600;

  Bool_t batch_mode = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  TCanvas* cv_t = new TCanvas("time_test","Time Test");
  cv_t->Divide(2,1);
  cv_t->cd(1);
  TGraph* gr = avec_draw(settings_time - make_root_time_offset(), "","", "settings_time","AP");
  gr->GetYaxis()->SetNdivisions(505, kTRUE);
  gr->GetYaxis()->SetTimeDisplay(1);
  gr->GetYaxis()->SetTimeFormat("%d/%m/%y %H:%M:%S"); 

  cv_t->cd(2);
  gr = avec_draw(unixTime - make_root_time_offset(), "","", "unixTime","AP");
  gr->GetYaxis()->SetNdivisions(505, kTRUE);
  gr->GetYaxis()->SetTimeDisplay(1);
  gr->GetYaxis()->SetTimeFormat("%d/%m/%y %H:%M:%S"); 

  cv_t->Write();

  gROOT->SetBatch(batch_mode);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double fine_min = vmin(fine_label);
  double fine_max = vmax(fine_label);
  std::cout << "Fine Label min: " << fine_min << "  max: " << fine_max << std::endl;

  Avec::size_type size = fine_label.size();
  Avec label_board(size);
  Avec label_level(size);

  Avec::size_type set_size = settings_time.size();
  Avec set_start(set_size, -1);
  Avec set_end(set_size, -1);

  int start_ctr = 0;
  int end_ctr = 0;

  std::cout << "Finding start and end in " << size << " events" << std::endl;
  // Find first event of each intensity step
  for(Avec::size_type i = 0; i < size; i++){ // loop over all events in this scan
    if(fine_label[i] == fine_min){
      start_ctr++;
      Avec::size_type midx = find_setting_start(settings_time, settings_level, unixTime[i]);
      if( set_start[midx] >= 0 ) throw LAS::Exception("Error in make_lsboard_level_labels: Found another start event for same level");
      set_start[midx] = eventnumber[i];
    }
  }

  // Find last event of each intensity step
  for(Avec::size_type i = 0; i < size; i++){ // loop over all events in this scan
    if(fine_label[i] == fine_max){
      end_ctr++;
      // From all start events that are smaller than the current eventnumber, take the one that is largest.
      Avec::size_type midx = lvmax(set_start * (set_start < eventnumber[i]));
      if( set_end[midx] >= 0 ) throw LAS::Exception("Error in make_lsboard_level_labels: Found another end event for same level");
      set_end[midx] = eventnumber[i];
    }
  }

  std::cout << "Found " << start_ctr << " start labels and " << end_ctr << " end labels" << std::endl;

  // Generate the labels for the events
  for(Avec::size_type is = 0; is < set_size; is++){
    if(settings_board[is] != board_id) continue;
    if(settings_level[is] < 0)continue;

    if(set_start[is] < 0){
       std::cerr << "Did not find settings start for level " << settings_level[is] << std::endl;
       continue;
    }
    if(set_end[is] < 0){
      std::cerr << "Did not find settings end for level " << settings_level[is] << std::endl;
      continue;
    }
    Avec mask = (eventnumber >= set_start[is]) && (eventnumber <= set_end[is]);
    std::cout << "Level " << settings_level[is] << " has " << vsum(mask) << " events, first is " << std::setprecision(12) << set_start[is] << ", last is " << std::setprecision(12) << set_end[is] << std::endl;
    label_level += mask * settings_level[is];
    label_board += mask * board_id;
  }

  label_board.Write("label_board");
  label_level.Write("label_level");
  f.Close();
}

//! Convert the signals from the beam group that was used as event counter into labels for each event
void decode_counter(const std::string& input_file, const std::string& directory, int subdet, int ring, fibre_swap swapped_fibres_pattern, bool plot, double threshold)
{
  TFile f(input_file.c_str(), "UPDATE");
  if(!f.cd(directory.c_str()))
    throw LAS::Exception( "Error in decode_counter: File " + input_file + " has no directory named " + directory );
    
    // Object with signal maxima
  std::string full_obj_name = directory + "/max_signal";
  LASGlobalData<Avec>* max_signal=0;
  f.GetObject(full_obj_name.c_str(),max_signal);
  if(!max_signal)
    throw LAS::Exception( "Error in decode_counter: Did not find max_signal object" );

  Avec eventnumber = avec_get(directory + "/eventnumber", f);
  Avec unixTime = avec_get(directory + "/unixTime", f);
  Avec signal_tecplus = avec_get(directory + "/signal_tecplus", f);
  Avec signal_tecminus = avec_get(directory + "/signal_tecminus", f);
  Avec signal_at = avec_get(directory + "/signal_at", f);
  if(signal_tecplus.empty() || signal_tecminus.empty() || signal_at.empty() || eventnumber.empty() || unixTime.empty())
    throw LAS::Exception( "Error in decode_counter: Could not find all Avecs" );

  // Convert signals to digital 1 and 0
  Avec sigsum_all;
  std::vector<Avec> bits(8);

  int max_zpos = ((subdet == 2 || subdet == 3) ? 6 : 9);
  std::cout << "looping over " << eventnumber.size() << " events" << std::endl;
  for(Avec::size_type i = 0; i < eventnumber.size(); i++){
    for(int beam = 0; beam < 8; beam++){
      double sigsum = 0;
      for(int zpos = 0; zpos < max_zpos; zpos++){
	if(subdet < 2 && zpos == 0)continue;
	sigsum += max_signal->GetEntry(subdet, ring, beam, zpos)[i];
      }
      sigsum_all.push_back(sigsum);
      bits[beam].push_back(sigsum > threshold ? 1 : 0);
    }
  }

  Bool_t batch_mode = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  TCanvas* cv = new TCanvas("binary_counter_distribution", "Binary Counter Distribution");
  cv->SetLogy();
  TH1* hist = avec_plot(sigsum_all);
  TLine thr_line(threshold, 0, threshold, hist->GetMaximum());
  thr_line.SetLineColor(2);
  thr_line.SetLineWidth(3);
  thr_line.SetLineStyle(2);
  thr_line.Draw();
  cv->Write();
  gROOT->SetBatch(batch_mode);

  // Convert digital 1 and 0 to decimal number
  Avec ctr;
  ctr = bits[0] + bits[1]*2 + bits[2]*4 + bits[3]*8 + bits[4]*16 + bits[5]*32 + bits[6]*64 + bits[7]*128;
  switch(swapped_fibres_pattern){
  case TECP_R6:
    if(subdet == 0 && ring == 1) ctr = bits[0] + bits[1]*2 + bits[2]*4 + bits[3]*8 + bits[7]*16 + bits[6]*32 + bits[5]*64 + bits[4]*128;
    break;
  case AT:
    if( (subdet == 2 || subdet == 3) && ring == -1) ctr = bits[1] + bits[0]*2 + bits[3]*4 + bits[2]*8 + bits[7]*16 + bits[5]*32 + bits[6]*64 + bits[4]*128;
    break;
  case NONE:
    break;
  default:
    std::cerr << "Error: Unknown fibre swap pattern." << std::endl;
  }

  ctr.Write("fine_label");

  std::cout << "sigsum_all.size(): " << sigsum_all.size() << std::endl;

  if(plot){
    new TCanvas("cv_sigsum_all","Digital Signals");
    avec_plot(sigsum_all, 1000);

    new TCanvas("cv_ctr","Digital Counter");
    //avec_draw(ctr,"","","","AP");
    avec_draw(eventnumber,ctr,"","","","AP");
  }
  f.Close();
}

//! Slice run into parts for each board or scan type
void slice_run(const std::string& output_filename)
{
  Avec signal_threshold;
  TFile f(output_filename.c_str(), "UPDATE");

  Avec board_list = avec_get("board_list", f);
  Avec ev1 = avec_get("ev_start", f);
  Avec ev2 = avec_get("ev_end", f);

  std::vector<std::string> directory_list;
  create_directory_list(board_list, directory_list);

  if(directory_list.size() != ev1.size() || ev1.size() != ev2.size()){
    std::cerr << "Inconsistent output settings" << std::endl;
    return;
  }

  if (signal_threshold.empty()) signal_threshold = Avec(ev1.size(), 0.0);

  // Object with signal maxima
  LASGlobalData<Avec>* max_signal = 0;
  f.GetObject("max_signal",max_signal);
  if(!max_signal){
    std::cerr << "Did not find max_signal object" << std::endl;
    return;
  }

  // Object with single positions
  LASGlobalData<Avec>* single_pos = 0;
  f.GetObject( "single_pos", single_pos );
  if( !single_pos ){
    std::cerr << "Did not find single_pos object" << std::endl;
    return;
  }

  Avec eventnumber = avec_get("eventnumber", f);
  Avec unixTime = avec_get("unixTime", f);
  Avec signal_tecp_r4 = avec_get("Details/signal_tecp_r4", f);
  Avec signal_tecp_r6 = avec_get("Details/signal_tecp_r6", f);
  Avec signal_tecm_r4 = avec_get("Details/signal_tecm_r4", f);
  Avec signal_tecm_r6 = avec_get("Details/signal_tecm_r6", f);
  Avec signal_tecplus = avec_get("Details/signal_tecplus", f);
  Avec signal_tecminus = avec_get("Details/signal_tecminus", f);
  Avec signal_at = avec_get("Details/signal_at", f);
  Avec signal_all = avec_get("signal_all", f);
  if(eventnumber.empty() || 
     unixTime.empty() || 
     signal_tecp_r4.empty() || 
     signal_tecp_r6.empty() || 
     signal_tecm_r4.empty() || 
     signal_tecm_r6.empty() || 
     signal_tecplus.empty() || 
     signal_tecminus.empty() || 
     signal_at.empty() || 
     signal_all.empty()){
    std::cerr << "Could not find all Avecs" << std::endl;
    return;
  }
  
  Avec::size_type nentries = eventnumber.size();
  std::string pwd = gDirectory->GetName();

  for(Avec::size_type out_idx = 0; out_idx < ev1.size(); out_idx++){
    std::cout << "Processing " << directory_list[out_idx] << std::endl;
    LASGlobalData<Avec> max_signal_board;
    LASGlobalData<Avec> single_pos_board;
    Avec eventnumber_board, unixTime_board, signal_tecplus_board, signal_tecminus_board, signal_at_board, signal_all_board;
    Avec signal_tecp_r4_board, signal_tecp_r6_board, signal_tecm_r4_board, signal_tecm_r6_board;

    for(Avec::size_type i = 0; i < nentries; i++){
      if(eventnumber[i] >= ev1[out_idx] && eventnumber[i] <= ev2[out_idx] && signal_all[i] > signal_threshold[out_idx]){
	eventnumber_board.push_back(eventnumber[i]);
	unixTime_board.push_back(unixTime[i]);
	signal_tecp_r4_board.push_back(signal_tecp_r4[i]);
	signal_tecp_r6_board.push_back(signal_tecp_r6[i]);
	signal_tecm_r4_board.push_back(signal_tecm_r4[i]);
	signal_tecm_r6_board.push_back(signal_tecm_r6[i]);
	signal_tecplus_board.push_back(signal_tecplus[i]);
	signal_tecminus_board.push_back(signal_tecminus[i]);
	signal_at_board.push_back(signal_at[i]);
	signal_all_board.push_back(signal_all[i]);
	LASGlobalDataLoop loop;
	do{
	  loop.GetEntry<Avec>(max_signal_board).push_back(loop.GetEntry<Avec>(*max_signal)[i]);
	  loop.GetEntry<Avec>(single_pos_board).push_back(loop.GetEntry<Avec>(*single_pos)[i]);
	}while(loop.next());
      }
    }

    std::cout << "Accumulated " << eventnumber_board.size() << " entries" << std::endl;

    if(!f.cd(directory_list[out_idx].c_str())){
      f.mkdir(directory_list[out_idx].c_str());
      f.cd(directory_list[out_idx].c_str());
    }

    max_signal_board.Write("max_signal");
    single_pos_board.Write("single_pos");
    eventnumber_board.Write("eventnumber");
    unixTime_board.Write("unixTime");
    signal_tecp_r4_board.Write("signal_tecp_r4");
    signal_tecp_r6_board.Write("signal_tecp_r6");
    signal_tecm_r4_board.Write("signal_tecm_r4");
    signal_tecm_r6_board.Write("signal_tecm_r6");
    signal_tecplus_board.Write("signal_tecplus");
    signal_tecminus_board.Write("signal_tecminus");
    signal_at_board.Write("signal_at");
    signal_all_board.Write("signal_all");
    f.cd(pwd.c_str());
  }
  f.Close();
  max_signal->Delete();
  eventnumber.Delete();
  unixTime.Delete();
  signal_tecplus.Delete();
  signal_tecminus.Delete();
  signal_at.Delete();
  signal_all.Delete();
}

//! Verify the results from get_scan_signal
void control_scan_signal(const std::string& resultfile)
{
  std::cout << "Drawing control plot" << std::endl;

  TFile f(resultfile.c_str(), "READ");
  if( ! f.IsOpen() ){
    std::cerr << "Could not open file " << resultfile << std::endl;
    return;
  }

  Avec unixTime = avec_get("unixTime", f);
  Avec eventnumber = avec_get("eventnumber", f);
  Avec signal_tecp_r4 = avec_get("signal_tecp_r4", f);
  Avec signal_tecp_r6 = avec_get("signal_tecp_r6", f);
  Avec signal_tecm_r4 = avec_get("signal_tecm_r4", f);
  Avec signal_tecm_r6 = avec_get("signal_tecm_r6", f);
  Avec signal_at = avec_get("signal_at", f);
  Avec signal_all = avec_get("signal_all", f);
  f.Close();

  if( unixTime.empty() || eventnumber.empty() || signal_tecp_r4.empty() || signal_tecp_r6.empty() || signal_tecm_r4.empty() || signal_tecm_r6.empty() || signal_at.empty() || signal_all.empty()){
    std::cerr << "Error in control_scan_signal: Some of the Avecs are empty" << std::endl;
    return;
  }

  std::vector<std::string> legend;
  legend.push_back("Endcap + Ring 4");
  legend.push_back("Endcap + Ring 6");
  legend.push_back("Endcap - Ring 4");
  legend.push_back("Endcap - Ring 6");
  legend.push_back("Alignment Tubes");

  new TCanvas("cv_evtnr", "Scan Signals");
  TMultiGraph* mgr = avec_draw(eventnumber, (signal_tecp_r4 | signal_tecp_r6 | signal_tecm_r4 | signal_tecm_r6 | signal_at), "Max Signal", "Event Nr.", "signal", "AP");
  //TMultiGraph* mgr = avec_draw(unixTime, (signal_tecp_r4 | signal_tecp_r6 | signal_tecm_r4 | signal_tecm_r6 | signal_at), "Max Signal", "evt ntr.", "signal", "AP");
  AddLegend(mgr, legend);

  new TCanvas("cv_signal_all", "Signal Sum");
  avec_plot(signal_all, 10000);

  new TCanvas("cv_evtnr_vs_time", "Event Numbers");
  TGraph* gr = avec_draw(unixTime - make_root_time_offset(), eventnumber, "EvtNr vs. time", "time", "Event Nr.");
  gr->GetXaxis()->SetNdivisions(505);
  gr->GetXaxis()->SetTimeDisplay(kTRUE);
  gr->GetXaxis()->SetTimeFormat("#splitline{%d.%m.%Y}{%H:%M}");

}

bool check_range(double event_nr, const Avec& ev_start, const Avec& ev_end)
{
  if(ev_start.size() != ev_end.size()) throw LAS::Exception("Error in check_range: ev_start and ev_end do not have the same size");

  for(Avec::size_type idx = 0; idx < ev_start.size(); idx++)
    if( event_nr >= ev_start[idx] && event_nr <= ev_end[idx] ) return true;

  return false;
}

Avec::size_type check_range_idx(double event_nr, const Avec& ev_start, const Avec& ev_end)
{
  if(ev_start.size() != ev_end.size()) throw LAS::Exception("Error in check_range_idx: ev_start and ev_end do not have the same size");

  for( Avec::size_type idx = 0; idx < ev_start.size(); idx++ )
    if( event_nr >= ev_start[idx] && event_nr <= ev_end[idx] ) return idx;

  return ev_start.size();
}

//! Retrieve the signal maxima for all modules and all events
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_scan_signal(const std::string& filename, const std::string& output_rootfilename, Bool_t plot)
{
  std::cout << "Running get_scan_signal" << std::endl;

  float good_profile_thresh_lo = 30;
  float good_profile_thresh_hi = 180;


  // Retrieve information from previous analysis steps
  TFile outf(output_rootfilename.c_str(),"UPDATE");
  if( !outf.IsOpen() ) throw LAS::Exception( "Error in get_scan_signal: Could not open file " + output_rootfilename);

  Avec eventnumber = avec_get("eventnumber", outf);
  Avec board_list = avec_get("board_list", outf);
  Avec ev_start = avec_get("ev_start", outf);
  Avec ev_end = avec_get("ev_end", outf);

  if( eventnumber.empty() ) 
    throw LAS::Exception("Error in get_scan_signal: Could not find Avec eventnumber");

  if( board_list.empty() || ev_start.empty() || ev_end.empty()) 
    throw LAS::Exception("Error in get_scan_signal: Could not find all required Avecs with scan boundaries (board_list, ev_start, ev_end)");

  Avec::size_type nr_of_scans = board_list.size();;
  std::vector<std::string> directory_list;
  create_directory_list(board_list, directory_list);
  if(directory_list.size() != nr_of_scans || ev_start.size() != nr_of_scans || ev_end.size() != nr_of_scans )
    throw LAS::Exception( "Error in get_scan_signal: Inconsistent information on scan boundaries");


  // Open the data file
  TFile f(filename.c_str());
  if(!f.IsOpen()) throw LAS::Exception( "Error in get_scan_signal: Could not open file " + filename);

  // Find the Tree
  TTree* tree = 0;
  f.GetObject("lasRawDataTree",tree);
  if(!tree) throw LAS::Exception( "Error in get_scan_signal: Did not find rawDataTree" );
  Int_t nentries = (Int_t)tree->GetEntries();
  std::cout << "The Tree has " << nentries << " Entries" << std::endl;
  if((unsigned int)nentries != eventnumber.size())
    throw LAS::Exception( "Error in get_scan_signal: Avec 'eventnumber' and Tree do not have the same size" );

  // Create container for Branch content
  std::vector<float> empty_buffer(512,0);
  LASGlobalData<std::vector<float> > theData( empty_buffer );
  LASGlobalData<std::vector<float> > *ptr = &theData;
  tree->SetBranchAddress("lasRawData",&ptr);
  Int_t evt_nr;
  tree->SetBranchAddress("eventnumber",&evt_nr);
  UInt_t unixTime;
  tree->SetBranchAddress("unixTime",&unixTime);

  // Define a strip range in which to look for the maximum
  // In this case the 60 central strips
  LASGlobalData<int> range_low, range_high;
  fill_strip_ranges(range_low, range_high);


  // Object for results
  std::vector< LASGlobalData<Avec> > max_signal( nr_of_scans );
  std::vector< LASGlobalData<Avec> > single_pos( nr_of_scans );
  std::vector< Avec > signal_all( nr_of_scans );
  std::vector< Avec::size_type > scan_ctr(nr_of_scans, 0);
  std::vector< Avec > unixTime_scan( nr_of_scans );
  std::vector< Avec > eventnumber_scan( nr_of_scans );

  // Loop over the tree
  for (Int_t entry_nr = 0; entry_nr < nentries; entry_nr++) {
    if(! (entry_nr%1000))std::cout << "Processing entry nr. " << entry_nr << "  Event Nr. " << std::setprecision(12) << eventnumber[entry_nr] << std::endl; 

    LASGlobalDataLoop loop;

    // Check to which scan this event belongs
    Avec::size_type scan_idx = check_range_idx( eventnumber[entry_nr], ev_start, ev_end );
    if( scan_idx >= nr_of_scans ) continue; // skip eventnumbers that are not within the range of any scan
    scan_ctr[scan_idx] ++;

    // Get an Entry from the tree (one event)
    tree->GetEntry(entry_nr);

    // Copy the data from the tree to the vectors of the corresponding scan
    unixTime_scan[scan_idx].push_back(unixTime);
    eventnumber_scan[scan_idx].push_back(evt_nr);
    double sig_all = 0;
    do{
      // Get the data for one specific module
      std::vector<float>& buff = loop.GetEntry<std::vector<float> >(theData);

      // Find the index with the highest entry      
      std::vector<float>::size_type max_idx = vector_max_idx(buff, loop.GetEntry<int>(range_low), loop.GetEntry<int>(range_high));
      float max =  buff.at(max_idx);

      loop.GetEntry<Avec>( max_signal[scan_idx] ).push_back( max );
      sig_all += max;
      loop.GetEntry<Avec>( single_pos[scan_idx] ).push_back( ( (max <= good_profile_thresh_hi) && (max >= good_profile_thresh_lo) ) ? calc_pos( buff, max_idx) : 0 );
    } while(loop.next());
    signal_all[scan_idx].push_back(sig_all);
  }

  f.Close();

  std::cout << "Computing sums and writing to file" << std::endl;
  outf.cd();
  std::string pwd = gDirectory->GetName();

  for(Avec::size_type scan_idx = 0; scan_idx < nr_of_scans; scan_idx++){
    std::cout << "Scan of board " << board_list[scan_idx] << " contains " << scan_ctr[scan_idx] << " entries" << std::endl;

    Avec signal_tecp_r6(scan_ctr[scan_idx]);
    Avec signal_tecp_r4(scan_ctr[scan_idx]);
    Avec signal_tecm_r6(scan_ctr[scan_idx]);
    Avec signal_tecm_r4(scan_ctr[scan_idx]);

    Avec signal_tib(scan_ctr[scan_idx]);
    Avec signal_tob(scan_ctr[scan_idx]);
    Avec signal_tecp_at(scan_ctr[scan_idx]);
    Avec signal_tecm_at(scan_ctr[scan_idx]);

    LASGlobalDataLoop loop1(LASGlobalDataLoop::TEC_PLUS_R4);
    do{
      signal_tecp_r4 += loop1.GetEntry<Avec>(max_signal[scan_idx]);
    }while(loop1.next());
    
    LASGlobalDataLoop loop2(LASGlobalDataLoop::TEC_PLUS_R6);
    do{
      signal_tecp_r6 += loop2.GetEntry<Avec>(max_signal[scan_idx]);
    }while(loop2.next());
    
    LASGlobalDataLoop loop3(LASGlobalDataLoop::TEC_MINUS_R4);
    do{
      signal_tecm_r4 += loop3.GetEntry<Avec>(max_signal[scan_idx]);
    }while(loop3.next());
    
    LASGlobalDataLoop loop4(LASGlobalDataLoop::TEC_MINUS_R6);
    do{
      signal_tecm_r6 += loop4.GetEntry<Avec>(max_signal[scan_idx]);
    }while(loop4.next());
    
    LASGlobalDataLoop loop5(LASGlobalDataLoop::TIB);
    do{
      signal_tib += loop5.GetEntry<Avec>(max_signal[scan_idx]);
    }while(loop5.next());

    LASGlobalDataLoop loop6(LASGlobalDataLoop::TOB);
    do{
      signal_tob += loop6.GetEntry<Avec>(max_signal[scan_idx]);
    }while(loop6.next());

    LASGlobalDataLoop loop7(LASGlobalDataLoop::TEC_PLUS_AT);
    do{
      signal_tecp_at += loop7.GetEntry<Avec>(max_signal[scan_idx]);
    }while(loop7.next());

    LASGlobalDataLoop loop8(LASGlobalDataLoop::TEC_MINUS_AT);
    do{
      signal_tecm_at += loop8.GetEntry<Avec>(max_signal[scan_idx]);
    }while(loop8.next());
    
    Avec signal_tecplus = signal_tecp_r4 + signal_tecp_r6;
    Avec signal_tecminus = signal_tecm_r4 + signal_tecm_r6;
    Avec signal_at = signal_tib + signal_tob + signal_tecp_at + signal_tecm_at;

    if(!outf.cd(directory_list[scan_idx].c_str())){
      outf.mkdir(directory_list[scan_idx].c_str());
      outf.cd(directory_list[scan_idx].c_str());
    }

    max_signal[scan_idx].Write("max_signal");
    single_pos[scan_idx].Write("single_pos");
    eventnumber_scan[scan_idx].Write("eventnumber");
    unixTime_scan[scan_idx].Write("unixTime");
    signal_all[scan_idx].Write("signal_all");
    signal_tecp_r4.Write("signal_tecp_r4");
    signal_tecp_r6.Write("signal_tecp_r6");
    signal_tecm_r4.Write("signal_tecm_r4");
    signal_tecm_r6.Write("signal_tecm_r6");
    signal_tecplus.Write("signal_tecplus");
    signal_tecminus.Write("signal_tecminus");
    signal_at.Write("signal_at");

    outf.cd(pwd.c_str());
  }

  outf.Close();

  // // Switch temporarily to batch mode, to suppress graphical output if not desired
  // Bool_t batch_mode = gROOT->IsBatch();
  // gROOT->SetBatch(!plot);

  // std::vector<std::string> legend;
  // legend.push_back("Endcap + Ring 4");
  // legend.push_back("Endcap + Ring 6");
  // legend.push_back("Endcap - Ring 4");
  // legend.push_back("Endcap - Ring 6");
  // legend.push_back("Alignment Tubes");

  // TCanvas* cv = new TCanvas("cv_evtnr", "Scan Signals");
  // TMultiGraph* mgr = avec_draw(eventnumber, (signal_tecp_r4 | signal_tecp_r6 | signal_tecm_r4 | signal_tecm_r6 | signal_at), "Max Signal", "evt ntr.", "signal", "AP");
  // AddLegend(mgr, legend);

  // // Create the subdirectory if needed
  // if(outf.GetDirectory("Control_Plots") == 0) outf.mkdir("Control_Plots");
  // outf.cd("Control_Plots");

  // cv->Write();

  // gROOT->SetBatch(batch_mode);
}

void create_directory_list(const Avec& board_list, std::vector<std::string>& directory_list)
{
  for(unsigned int i = 0; i < board_list.size(); i++){
    std::ostringstream of;
    of << "Scan_LsBoard_" << board_list[i]; 
    directory_list.push_back(of.str());
  }
}

//! Find first and last event for the scans of each board
void find_scan_boundaries(const std::string& output_filename)
{
  TFile f(output_filename.c_str(), "UPDATE");

  Avec lsboard_boardid = avec_get("lsboard_boardid", f);
  Avec lsboard_level = avec_get("lsboard_level", f);
  Avec lsboard_set_time = avec_get("lsboard_set_time", f);
  Avec unixTime = avec_get("unixTime", f);
  Avec eventnumber = avec_get("eventnumber", f);

  if(lsboard_boardid.empty() || lsboard_level.empty() || lsboard_set_time.empty() || unixTime.empty() || eventnumber.empty()){
    std::cerr << "Error in lsboard_scan_analysis: Could not find all required Avecs" << std::endl;
    return;
  }

  time_t utstart = unixTime.front();
  time_t utend = unixTime.back();
  std::cout << "unixTime start: " << ctime(&utstart) << "  unixTime end: " << ctime(&utend) << std::endl;

  // Find beginning and end of scans for each board
  Avec ev_start, ev_end, board_list;
  int current_board  = -1;
  unsigned int current_start = 0;

  const double max_evt_nr = vmax(eventnumber);

  for(unsigned int idx = 0; idx <= lsboard_set_time.size(); idx++){
    time_t set_time = lsboard_set_time[idx];
    if(lsboard_level[idx] == -1){
      Avec m1 = unixTime >= set_time; // All events that come after the start
      unsigned int utidx = lvmin( max_evt_nr * (1 - m1) + m1 * eventnumber ); // Index of event with smallest eventnumber and coming after the start

      current_board = lsboard_boardid[idx];
      current_start = eventnumber[utidx];
      std::cout << "Found Scan start for board " << current_board 
		<< " at time " << ctime(&set_time) << "(" << set_time << ")" 
		<< " Event Nr. " << std::setprecision(12) << current_start << std::endl;
    }
    if(lsboard_level[idx] < -1){
      if(current_board != lsboard_boardid[idx]){
	std::cerr << "Found end of scan tag for a different board id!\n"
		  << "Start board ID: " << current_board << "  End board ID: " << lsboard_boardid[idx] << std::endl;
	continue;
      }
      Avec m1 = unixTime <= set_time; // All events that come before the end
      unsigned int utidx = lvmax( m1 * eventnumber ); // Index of event with largest eventnumber and coming before the start

      std::cout << "Found Scan end for board " << lsboard_boardid[idx] 
		<< " at time " << ctime(&set_time) << "(" << set_time << ")" 
		<< " Event Nr. " << eventnumber[utidx] << std::endl;
      ev_start.push_back(current_start);
      ev_end.push_back(eventnumber[utidx]);
      board_list.push_back(current_board);
      Avec scan_mask = (eventnumber >= current_start) && (eventnumber <= eventnumber[utidx]);
      std::cout << "This scan has " << vsum(scan_mask) << " entries" << std::endl;
      current_board = -1;
      current_start = 0;
    }
  }
  if(board_list.empty()){
    std::cerr << "Error in lsboard_scan_analysis, found no lsboard scans" << std::endl;
    return;
  }

  board_list.Write("board_list");
  ev_start.Write("ev_start");
  ev_end.Write("ev_end");
  f.Close();
}


//! Retrieve the event numbers and timestamps
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_eventnr_and_time(const std::string& filename, const std::string& output_rootfilename, Bool_t plot)
{
  std::cout << "Running get_event_and_time" << std::endl;

  // Open the data file
  TFile f(filename.c_str());
  if(!f.IsOpen()){
    std::cerr << "Could not open file " << filename << std::endl;
    return;
  }

  // Find the Tree
  TTree* tree = 0;
  f.GetObject("lasRawDataTree",tree);
  if(!tree){
    std::cerr << "Did not find rawDataTree" << std::endl;
    return;
  }

  // Create container for Branch content
  Int_t evt_nr;
  tree->SetBranchAddress("eventnumber",&evt_nr);
  // Create container for Branch content
  UInt_t unixTime;
  tree->SetBranchAddress("unixTime",&unixTime);

  tree->SetBranchStatus("*",0); //disable all branches
  tree->SetBranchStatus("eventnumber",1);
  tree->SetBranchStatus("unixTime",1);

  // Object for results
  Avec v_unixTime;
  Avec eventnumber;

  Int_t nentries = (Int_t)tree->GetEntries();
  std::cout << "The Tree has " << nentries << " Entries" << std::endl;

  // Loop over the tree
  for (Int_t entry_nr=0; entry_nr<nentries; entry_nr++) {
    // Get an Entry from the tree (one event)
    tree->GetEntry(entry_nr);

    v_unixTime.push_back(unixTime);
    eventnumber.push_back(evt_nr);

  }
  f.Close();

  TFile outf(output_rootfilename.c_str(),"UPDATE");

  eventnumber.Write("eventnumber");
  v_unixTime.Write("unixTime");

  // Switch temporarily to batch mode, to suppress graphical output if not desired
  Bool_t batch_mode = gROOT->IsBatch();
  gROOT->SetBatch(!plot);

  TCanvas* cv = new TCanvas("cv_eventnuber_vs_unixTime", "Event Numbers");
  TGraph* gr = avec_draw(eventnumber, v_unixTime, "Event Numbers", "time", "evt ntr.", "AP");
  gr->GetXaxis()->SetNdivisions(505);
  gr->GetXaxis()->SetTimeDisplay(kTRUE);
  gr->GetXaxis()->SetTimeFormat("#splitline{%d.%m.%Y}{%H:%M}");

  // Create the subdirectory if needed
  if(outf.GetDirectory("Control_Plots") == 0) outf.mkdir("Control_Plots");
  outf.cd("Control_Plots");

  cv->Write();

  gROOT->SetBatch(batch_mode);

  outf.Close();
}

//! Convert the text file with Laser Board settings into root format
void lsboard_get_log_data(const std::string& filename, const std::string& output_filename)
{
  const bool debug_flag = false;
  if(debug_flag) std::cout << "Debug flag is " << (debug_flag ? "true" : "false") << std::endl;
  Avec day, mon, year, hour, min, sec, board, level, delay1, delay2;
  vecread( (day | mon | year | hour | min | sec | board | level | delay1 | delay2), filename);

  // Needed to make mktime behave properly. Best would be to use <chrono> or some boost library for time manipulations
  setenv("TZ", "", 1);
  tzset();

  Avec t1;
  double current_board = -1;
  bool scan_started = false;
  Avec::size_type last_start_idx = 0;
  for(Avec::size_type i = 0; i < day.size(); i++){
    tm mydate;
    mydate.tm_sec  = (int)sec[i];
    mydate.tm_min  = (int)min[i];
    mydate.tm_hour = (int)hour[i];
    mydate.tm_mday = (int)day[i];
    mydate.tm_mon  = (int)mon[i]-1;
    mydate.tm_year = (int)year[i]-1900;
    time_t mytime = mktime(&mydate);

    if(debug_flag && (i<=1 || i >= day.size() - 2))
      std::cout << day[i] << "."<< mon[i] << "."<< year[i] << " "<< hour[i] << ":"<< min[i] << ":"<< sec[i] 
     		<< " gives " << asctime(&mydate) << " and " << mytime << " " << ctime(&mytime) << std::endl;

    t1.push_back(mytime);

    if(scan_started){
      switch ((int)board[i]){
      case -2:
	board[i] = current_board;
	board[last_start_idx] = current_board;
	current_board = -1;
	scan_started = false;
	break;
      case -1:
	std::cerr << "Warning in lsboard_get_log_data: Found start of scan tag, but current scan has not ended yet." << std::endl;
	break;
      default:
	if(current_board != -1 && board[i] != current_board) std::cerr << "Warning in lsboard_get_log_data: Board Id changed during a scan" << std::endl;
	current_board = board[i];
      }
    }
    else{
      switch ((int)board[i]){
      case -2:
	std::cerr << "Warning in lsboard_get_log_data: Found end of scan tag, but no scan has started yet" << std::endl;
	break;
      case -1:
	scan_started = true;
	last_start_idx = i;
	break;
      default:
	std::cerr << "Warning in lsboard_get_log_data: Found scan entry, but scan has not started yet" << std::endl;
      }
    }
  }

  TFile of(output_filename.c_str(),"UPDATE");
  board.Write("lsboard_boardid");
  level.Write("lsboard_level");
  t1.Write("lsboard_set_time");

  // Switch temporarily to batch mode, to suppress graphical output
  Bool_t batch_mode = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  // Create the subdirectory if needed
  if(of.GetDirectory("Control_Plots") == 0) of.mkdir("Control_Plots");
  of.cd("Control_Plots");

  // Canvas with board Id as function of time
  TCanvas* cv = new TCanvas("cv_lsboard_boardid","LsBoard Board Id");
  TGraph* gr = avec_draw(t1 - make_root_time_offset(), board, "Laser Board Id", "Settings Time", "Id", "AP");
  gr->GetXaxis()->SetTimeDisplay(1);
  gr->GetXaxis()->SetTimeFormat("#splitline{%d.%m.%Y}{%H:%M:%S}"); 
  gr->GetXaxis()->SetNdivisions(505, kTRUE);
  cv->Write();

  // Canvas with intensity setting as function of time
  cv = new TCanvas("cv_lsboard_level","LsBoard Intensity Level Settings");
  gr = avec_draw(t1 - make_root_time_offset(), level, "Laser Board Intensity Level", "Settings Time", "level", "AP");
  gr->GetXaxis()->SetTimeDisplay(1);
  gr->GetXaxis()->SetTimeFormat("#splitline{%d.%m.%Y}{%H:%M:%S}"); 
  gr->GetXaxis()->SetNdivisions(505, kTRUE);
  cv->Write();

  gROOT->SetBatch(batch_mode);

  of.Close();
}



void cali_pos_control(const std::string& raw_data_filename, const std::string& results_filename, unsigned int board, int det, int ring, int beam, int zpos)
{
  try{
    std::cout << "Running cali_pos_control" << std::endl;

    // Open the data file
    TFile f(raw_data_filename.c_str());
    if(!f.IsOpen()){
      std::cerr << "Could not open file " << raw_data_filename << std::endl;
      return;
    }
    
    Avec eventnumber_all = avec_get("eventnumber", results_filename);
    if(eventnumber_all.empty()) throw LAS::Exception("Did not find Avec eventnumber");

    // Find the Tree
    TTree* tree = 0;
    f.GetObject("lasRawDataTree",tree);
    if(!tree) throw LAS::Exception("Did not find rawDataTree");

    // Create container for Branch content
    std::vector<float> empty_buffer(512,0);
    LASGlobalData<std::vector<float> > theData( empty_buffer );
    LASGlobalData<std::vector<float> > *ptr = &theData;
    tree->SetBranchAddress("lasRawData",&ptr);

    Int_t evt_nr;
    Int_t retval = tree->SetBranchAddress("eventnumber", &evt_nr);
    if( retval < 0 ) throw LAS::Exception("Failed to set branch address for eventnumber");
    Int_t orb_nr;
    retval = tree->SetBranchAddress("orbitNumber", &orb_nr);
    if( retval < 0 ) throw LAS::Exception("Failed to set branch address for orbitNumber");
    Int_t lumi_bl;
    retval = tree->SetBranchAddress("lumiblock", &lumi_bl);
    if( retval < 0 ) throw LAS::Exception("Failed to set branch address for lumiblock");
    UInt_t unix_time;
    retval = tree->SetBranchAddress("unixTime", &unix_time);
    if( retval < 0 ) throw LAS::Exception("Failed to set branch address for unixTime");
    UInt_t microsec_off;
    retval = tree->SetBranchAddress("microsecondOffset", &microsec_off);
    if( retval < 0 ) throw LAS::Exception("Failed to set branch address for microsecondOffset");


    Avec board_list = avec_get("board_list", results_filename);
    std::vector<std::string> directory_list;
    create_directory_list(board_list, directory_list);

    int board_id = board_list[board];
    std::string& directory = directory_list[board];

    std::cout << "Creating profiles for board " << board_id << std::endl;
    //LASGlobalDataLoop::loop_type loop_t = convert_boardid_to_loop(board_id);

    // Open the results file
    TFile res_f(results_filename.c_str());
    if(!res_f.IsOpen()){
      std::cerr << "Could not open file " << results_filename << std::endl;
      return;
    }
    
    Avec eventnumber = avec_get(directory + "/eventnumber", res_f);
    Avec label_board = avec_get(directory + "/label_board", res_f);
    Avec label_level = avec_get(directory + "/label_level", res_f);
    Avec label_fine = avec_get(directory + "/fine_label", res_f);
    if(label_board.empty() || label_level.empty() || label_fine.empty())
      throw LAS::Exception("Could not find all Avecs with labels");
    

    // Object with signal maxima
    LASGlobalData<Avec>* max_signal=0;
    res_f.GetObject( (directory + "/max_signal").c_str(), max_signal);
    if(!max_signal) throw LAS::Exception("Did not find max_signal object");

    // Object with single positions
    LASGlobalData<Avec>* single_pos=0;
    res_f.GetObject( (directory + "/single_pos").c_str(), single_pos);
    if(!single_pos) throw LAS::Exception("Did not find single_pos object");

    // Define a strip range in which to look for the maximum
    // In this case the 60 central strips
    LASGlobalData<int> range_low, range_high;
    fill_strip_ranges(range_low, range_high);

    float good_profile_thresh_lo = 30;
    float good_profile_thresh_hi = 200;

    int ctr = 0;
    std::vector<Avec> profile_list;
    Avec propos;
    Avec propos2;
    Avec delay_vals;
    for( unsigned int i = 0; i < label_board.size(); i++){
      if( label_board[i] != board_id || label_level[i] != 30 ) continue;
      if( (single_pos->GetEntry( det, ring, beam, zpos ))[i] < 1 ) continue;
      Avec::size_type li;
      for( li = 0; li < eventnumber_all.size(); li++) if(eventnumber_all[li] == eventnumber[i]) break;
      if(li == eventnumber_all.size()) throw LAS::Exception("Mismatch in eventnumbers");
      tree->GetEntry(li);
      std::vector<float>& buff = theData.GetEntry( det, ring, beam, zpos );

      std::vector<float>::size_type max_idx = vector_max_idx(buff, range_low.GetEntry( det, ring, beam, zpos ), range_high.GetEntry( det, ring, beam, zpos ) );
      float max =  buff.at(max_idx);
      float pos2 = calc_pos( buff, max_idx);
      std::cout << "eventnumber: " << std::setprecision(12) << eventnumber[i] << "  li: " << li;
      std::cout << "   max: " << max;
      std::cout << "   pos: " << ( ( (max <= good_profile_thresh_hi) && (max >= good_profile_thresh_lo) ) ? pos2 : 0 );
      std::cout << "   single_pos: " << (single_pos->GetEntry( det, ring, beam, zpos ))[i] << std::endl;

      profile_list.push_back( Avec( buff ) );
      propos.push_back( (single_pos->GetEntry( det, ring, beam, zpos ))[i] );
      propos2.push_back(pos2);
      delay_vals.push_back(label_fine[i]);
      ctr ++;
    }
    std::cout << "Found " << ctr << " entries for this setting" << std::endl;

    res_f.Close();
    f.Close();

    if(!profile_list.empty()){
      TCanvas* cv1 = new TCanvas("profiles", "profiles");
      avec_draw(Avec(512, 1, 512), profile_list);
    }
    if( ! propos.empty() ){
      TCanvas* cv2 = new TCanvas("positions", "positions");
      //avec_draw( Avec(propos.size(), 0, propos.size()), (propos | propos2), "Positions", "delay", "signal", "AP" );
      avec_draw( delay_vals, (propos | propos2), "Positions", "delay", "signal", "AP" );
    }
  }
  catch (LAS::Exception& e)  {std::cerr << "Caught LAS::Exception, saying: " << e.what() << std::endl;}
  
  return;

}

