#define __AVECROOT__
#include "Avec.h"
#include "Avec2D.h"
#include "LASGlobalData.h"
#include "LAS_basic_tools.h"
#include "LAS_globaldata_tools.h"
#include "LAS_vectorfloat_tools.h"
#include "LAS_Tec_Reconstruction.h"
#include "LAS_control_plots.h"
#include "LAS_RDC_tools.h"


#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TTree.h"
#include "TLine.h"
#include "THStack.h"
#include "TPaveLabel.h"
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <ctime>

double calc_chi2(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, bool verbose)
{
  double chi2 = 0;
  double max = 0;
  int ctr = 0;

  LASGlobalDataLoop lp(LASGlobalDataLoop::AT);
  do{
    if(lp.GetEntry(err2) > 1e-10){
      ctr ++;
      chi2 += fabs(lp.GetEntry(dif)) / lp.GetEntry(err2);
      max = std::max(fabs(lp.GetEntry(dif)) , max);
    }
  } while (lp.next());
  
  if(verbose){
    std::cout << "ctr: " << ctr << "  max: " << max << std::endl;
  }

  return (ctr > 0 ? chi2 / ctr : 0);
}

void check_ref_pos(const std::string& ref_file, int block_nr, double rms_cut, const std::string& pos_prefix)
{
  LASGlobalData<double> ref_pos;
  LASGlobalData<double> ref_err;
  LASGlobalData<int> ref_mask;

  get_ref_pos(ref_file, ref_pos, ref_err, ref_mask, block_nr, rms_cut, pos_prefix);

  //print_global_data(ref_pos);
  draw_global_data(ref_pos, ref_err, ref_mask);

}

void check_run(const std::string& result_file, const std::string& pos_prefix)
{
  std::cout << "Calling check_run( "<< result_file << ", " << pos_prefix << " );" << std::endl;
  std::vector<std::string> file_list(1, result_file);
  check_run_list(result_file, file_list, pos_prefix);

}

void check_run_list(const std::string& run_list_name, const std::string& pos_prefix)
{
  std::cout << "Calling check_run_list( " << run_list_name << ", " << pos_prefix << " );" << std::endl;

  std::string ref_filename;
  std::vector<std::string> file_list;

  get_file_list(run_list_name, ref_filename, file_list);

  check_run_list(ref_filename, file_list, pos_prefix);

}


void fill_globaldata_avec(const LASGlobalData<double>& data, LASGlobalData<Avec>& history)
{
  LASGlobalDataLoop lp;
  do{
    lp.GetEntry(history).push_back(lp.GetEntry(data));
  }while (lp.next());
}

void check_run_list(const std::string& ref_filename, const std::vector<std::string>& file_list, const std::string& pos_prefix)
{
  std::cout << "Calling check_run_list( " << ref_filename << ", (...), " << pos_prefix << " );" << std::endl;

  Avec time_axis;
  LASGlobalData<Avec> diff_hist;
  LASGlobalData<Avec> err_hist;

  LASGlobalData<double> ref_pos;
  LASGlobalData<double> ref_err;
  LASGlobalData<int> ref_mask;
  get_ref_pos(ref_filename, ref_pos, ref_err, ref_mask, -1, 0.5, pos_prefix);
  
  for(unsigned int i = 0; i < file_list.size(); i++){
    std::cout << "Processing file " << file_list[i] << "   ";

    Avec good_run = avec_get("good_run", file_list[i]);
    if(good_run.empty()){
      std::cout << "Warning: No good_run quality flag in this file " << std::endl;
    }
    else
      std::cout << "Run is flagged as " << (good_run[0] == 0 ? "bad,  " : "good, ");

    Avec good_block = avec_get("good_block",file_list[i]);
    std::cout << "Nr. of good blocks: " << vsum(good_block) << "\n" << std::endl;
    Avec block_timestamp = avec_get("block_timestamp", file_list[i]);

    for(unsigned int blnr = 0; blnr < good_block.size(); blnr++){
      if(good_block[blnr] == 0)continue;
      std::ostringstream blkmean;
      std::ostringstream blkrms;
      std::ostringstream blknorm;
      blkmean << pos_prefix <<blnr;
      //blkrms  <<    "rms_av_" <<blnr;
      blkrms  <<    "pos_error_" <<blnr;
      blknorm <<      "norm_" <<blnr;

      LASGlobalData<double> pos = global_data_get<double>(blkmean.str(),file_list[i]);
      LASGlobalData<int> norm = global_data_get<int>(blknorm.str(),file_list[i]);
      LASGlobalData<double> rms = global_data_get<double>(blkrms.str(),file_list[i]);

      LASGlobalData<double> dif = (pos - ref_pos)*pitch()/r0();
      correct_signs(dif);
      LASGlobalData<double> error = sqrt((rms*rms + ref_err*ref_err)/norm)*pitch()/r0();

      fill_globaldata_avec(dif, diff_hist);
      fill_globaldata_avec(error, err_hist);
      time_axis.push_back( block_timestamp[blnr] );
    }

  }

  diff_hist *= 1e6;
  err_hist *= 1e6;
  LASGlobalData<Avec> xval(time_axis);

  int flag = Avec::VERBOSE_FLAG;
  Avec::VERBOSE_FLAG = 0;
  draw_global_data(diff_hist, err_hist, xval, "AP", "time");
  Avec::VERBOSE_FLAG = flag;

}


void get_file_list(const std::string& filename, std::vector<std::string>& file_list, const std::string& path_prefix )
{
  std::ifstream f(filename.c_str());
  
  while(f){
    std::string buffer = get_next_line(f);
    if(buffer != ""){
      file_list.push_back( (path_prefix != "" ? path_prefix + "/" : "" ) + buffer );
    }
  }
}

void get_file_list(const std::string& filename, std::string& ref_file, std::vector<std::string>& file_list, const std::string& path_prefix )
{
  std::ifstream f(filename.c_str());
  
  ref_file = (path_prefix != "" ? path_prefix + "/" : "" ) + get_next_line(f);
  std::cout << "Reference file: " << ref_file << std::endl;

  while(f){
    std::string buffer = get_next_line(f);
    if(buffer != ""){
      file_list.push_back( (path_prefix != "" ? path_prefix + "/" : "" ) + buffer );
    }
  }
}

bool get_ref_pos_average(const std::string& ref_file, LASGlobalData<double>& ref_pos, LASGlobalData<double>& ref_err, LASGlobalData<int>& ref_mask, double rms_cut, const std::string& pos_prefix)
{
  TFile f(ref_file.c_str(), "READ");
  Avec good_block = avec_get("good_block", f);
  if(vsum(good_block) == 0){
    std::cerr << "No good block in reference file " << ref_file << std::endl;
    return false;
  }

  std::cout << vsum(good_block) << " good blocks out of " << good_block.size() << std::endl;

  ref_err = LASGlobalData<double>(0.);

  LASGlobalData<int> norm(0);

  for(Avec::size_type idx = 0; idx < good_block.size(); idx++){
    if(good_block[idx] != 1) continue;

    std::ostringstream block_nr_stream;
    block_nr_stream << idx;
    std::string block_nr = block_nr_stream.str();

    //LASGlobalData<int> norm;
    LASGlobalData<int> block_norm = global_data_get<int>("norm_" + block_nr, f);
    LASGlobalData<double> pos = global_data_get<double>(pos_prefix + block_nr, f);
    LASGlobalData<double> rms = global_data_get<double>("pos_error_" + block_nr, f) / sqrt(block_norm);
    LASGlobalData<int> block_mask = global_data_get<int>("positions_mask_" + block_nr, f) && (rms < rms_cut) && !(rms != rms);

    remove_nan(rms);
    pos *= block_mask;
    rms *= block_mask;
    norm += block_mask;

    //ref_pos +=  global_data_get<double>(pos_prefix + block_nr, f) * block_mask;
    ref_pos +=  pos * block_mask;
    ref_err +=  rms * rms;
  }
  f.Close();
  ref_pos = ref_pos / norm;
  ref_err = sqrt(ref_err / norm);
  remove_nan(ref_pos);
  remove_nan(ref_err);
  ref_mask = norm > 0;
  return true;
}

bool get_ref_pos(const std::string& ref_file, LASGlobalData<double>& ref_pos, LASGlobalData<double>& ref_err, LASGlobalData<int>& ref_mask, int block_nr, double rms_cut, const std::string& pos_prefix)
{
  bool retval = true;

  Avec good_run = avec_get("good_run", ref_file);
  if(good_run.empty()){
    std::cerr << "Warning in get_ref_pos: No good_run quality flag for reference file " << std::endl;
  }
  else{
    if(good_run[0] == 0){
      std::cerr << "Error in get_ref_pos: Reference run is flagged as not good" << std::endl;
      retval = false;
    }
  }

  Avec good_block = avec_get("good_block", ref_file);
  if(vsum(good_block) == 0){
    std::cerr << "No good block in reference file " << ref_file << std::endl;
    return false;
  }

  // If block_nr == -1 take average of all blocks
  if(block_nr == -1){
    std::cout << "Getting reference position from average of file " << ref_file << std::endl;
    return get_ref_pos_average(ref_file, ref_pos, ref_err, ref_mask, rms_cut, pos_prefix) ? retval : false;
  }

  // If block_nr == -2, use last block of file
  if(block_nr == -2) block_nr = good_block.size() - 1;

  if(block_nr < 0 || (unsigned int)block_nr >= good_block.size()){
    std:: cerr << "Error in get_ref_pos: Invalid block Nr.: " << block_nr << "  (highest block nr. is " << good_block.size() << std::endl;
    return false;
  }

  if(good_block[block_nr] != 1){
    std::cerr << "Error in get_ref_pos: Block Nr. " << block_nr << " is marked as 'not good'" << std::endl;
    return false;
  }

  ref_mask = LASGlobalData<int>(1);

  std::ostringstream blnr;
  blnr << block_nr;

  ref_pos =  global_data_get<double>(pos_prefix + blnr.str(), ref_file);

  LASGlobalData<double> rms = global_data_get<double>("pos_error_" + blnr.str(), ref_file);
  LASGlobalData<double> rms_mask = rms < rms_cut;
  ref_mask = rms_mask && global_data_get<int>("positions_mask_" + blnr.str(), ref_file);

  return retval;
}

void fill_global_tec_data(const LASGlobalData<double>& theData, const LASGlobalData<int>& mask, Avec& tecp_r4, Avec& tecp_r6, Avec& tecm_r4, Avec& tecm_r6)
{
  LASGlobalDataLoop l1(LASGlobalDataLoop::TEC_PLUS_R4);
  do{
    if(l1.GetEntry<int>(mask) != 0) tecp_r4.push_back(l1.GetEntry<double>(theData));
  }while(l1.next());
  
  LASGlobalDataLoop l2(LASGlobalDataLoop::TEC_PLUS_R6);
  do{
    if(l2.GetEntry<int>(mask) != 0) tecp_r6.push_back(l2.GetEntry<double>(theData));
  }while(l2.next());
  
  LASGlobalDataLoop l3(LASGlobalDataLoop::TEC_MINUS_R4);
  do{
    if(l3.GetEntry<int>(mask) != 0) tecm_r4.push_back(l3.GetEntry<double>(theData));
  }while(l3.next());
  
  LASGlobalDataLoop l4(LASGlobalDataLoop::TEC_MINUS_R6);
  do{
    if(l4.GetEntry<int>(mask) != 0) tecm_r6.push_back(l4.GetEntry<double>(theData));
  }while(l4.next());

}


//! Convert alignment parameter index to module acces index //
void apidx2modidx(int subdet, int zpos, int&subdet2, int& zpos2)
{
  subdet2=subdet;
  zpos2=zpos;
  switch(subdet){
  case 0: // TEC+
    break;
  case 1: // TEC-
    break;
  case 2: // TIB+
    break;
  case 3: // TIB-
    subdet2=2;
    zpos2+=3;
    break;
  case 4: // TOB+
    subdet2=3;
    break;
  case 5: // TOB-
    subdet2=3;
    zpos2+=3;
    break;
  }
  return;
}

//! Convert module acces index to alignment parameter index //
void modidx2apidx(int subdet, int zpos, int&subdet2, int& zpos2)
{
  subdet2=subdet;
  zpos2=zpos;
  switch(subdet){
  case 0: // TEC+
    break;
  case 1: // TEC-
    break;
  case 2: // TIB
    if(zpos>2){
      subdet2=3;
      zpos2-=3;
    }
    break;
  case 3: // TOB
    if(zpos<3) subdet2=4;
    else{
      subdet2=5;
      zpos2-=3;
    }
    break;
  }
  return;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 



// get a set of positions and a mask for processing
// if block_nr == -1, average over the entire file
// flag == 0: use pos_av_xx and rms_av_xx (use rms_cut for creating the mask)
// flag == 1: use positions_xx and positions_mask_xx
void get_pos_mask(const std::string& filename, int block_nr, LASGlobalData<double>& pos, LASGlobalData<int>& mask, int flag, double rms_cut)
{
  std::ostringstream pos_name;
  std::ostringstream rms_name;
  std::ostringstream mask_name;


  if(block_nr == -1){
    switch(flag){
    case 0:{
      mask *= 0;
      mask += 1;
      pos *= 0;
      Avec time = avec_get("block_timestamp", filename);
      Avec::size_type size = time.size();
      for(Avec::size_type idx = 0; idx < size; idx++){
	pos_name.clear();
	pos_name.str("");
	rms_name.clear();
	rms_name.str("");
	pos_name << "pos_av_" << idx;
	rms_name << "rms_av_" << idx;
	pos +=  global_data_get<double>(pos_name.str(), filename);
	mask = mask && ( global_data_get<double>(rms_name.str(), filename) < rms_cut);
      }
      pos /= size;
      break;}
    case 1:
      pos = global_data_get<double>("positions_strips", filename);
      mask = global_data_get<int>("results_mask", filename);
      break;
    default:
      std::cerr << "Error in get_pos_mask: flag " << flag << " has no implementation for block_nr == -1" << std::endl;
      return;
    }
  }
  else{
    switch(flag){
    case 0:
      pos_name << "pos_av_" << block_nr;
      rms_name << "rms_av_" << block_nr;
      pos = global_data_get<double>(pos_name.str(), filename);
      mask *= 0;
      mask += 1;
      mask = mask && (global_data_get<double>(rms_name.str(), filename) < rms_cut);
      break;
    case 1:
      pos_name << "positions_" << block_nr;
      mask_name << "positions_mask_" << block_nr;
      pos = global_data_get<double>(pos_name.str(), filename);
      mask = global_data_get<int>(mask_name.str(), filename);
      break;
    default:
      std::cerr << "Error in get_pos_mask: flag " << flag << " has no implementation for block_nr != -1" << std::endl;
      return;
    }
  }
}

void at_beam_draw(const std::string& data_file, int block_nr, const std::string& ref_file, int ref_block, int flag, double rms_cut)
{
  std::cout << "Reference file: " << ref_file << std::endl;

  LASGlobalData<double> ref_pos;
  LASGlobalData<int> ref_mask;
  get_pos_mask(ref_file, ref_block, ref_pos, ref_mask, flag, rms_cut);

  LASGlobalData<double> pos;
  LASGlobalData<int> mask;
  get_pos_mask(data_file, block_nr, pos, mask, flag, rms_cut);

  LASGlobalData<double> diff =  (ref_pos - pos)*pitch()/r0();
  mask = ref_mask && mask;
  correct_signs(diff);

  print_global_data<int>(mask);

  Avec2D diff_tob(8);
  Avec2D xval_tob(8);

  LASGlobalDataLoop loop(LASGlobalDataLoop::TOB);
  do{
    if(loop.GetEntry<int>(mask)){
      int beam = loop.get_beam();
      diff_tob[beam].push_back(loop.GetEntry<double>(diff));
      xval_tob[beam].push_back(loop.GetEntry<double>(zpos()));
    }
  } while (loop.next());

  Avec2D diff_tib(8);
  Avec2D xval_tib(8);
  loop = LASGlobalDataLoop(LASGlobalDataLoop::TIB);
  do{
    if(loop.GetEntry<int>(mask)){
      int beam = loop.get_beam();
      diff_tib[beam].push_back(loop.GetEntry<double>(diff));
      xval_tib[beam].push_back(loop.GetEntry<double>(zpos()));
    }
  } while (loop.next());

  Avec2D diff_tec(8);
  Avec2D xval_tec(8);
  loop = LASGlobalDataLoop(LASGlobalDataLoop::TEC_AT);
  do{
    if(loop.GetEntry<int>(mask)){
      int beam = loop.get_beam();
      diff_tec[beam].push_back(loop.GetEntry<double>(diff));
      xval_tec[beam].push_back(loop.GetEntry<double>(zpos()));
    }
  } while (loop.next());

  TCanvas* cv = new TCanvas("at_beam_draw","AT beams");
  TMultiGraph* mgr;
  cv->Divide(3,3);
  for(int i = 0; i<8; i++){
    cv->cd(i + 1);
    mgr = avec_draw((xval_tob[i] | xval_tib[i] | xval_tec[i]) , (diff_tob[i] | diff_tib[i] | diff_tec[i]), "", "zpos [mm]", "diff","AP");
  }

  std::vector<std::string> leg;
  leg.push_back("TOB");
  leg.push_back("TIB");
  leg.push_back("TEC");
  cv->cd(9);
  AddLegend(mgr, leg);
}



// Helping functions

//! Calculate the positions given the profiles
float calc_pos(const std::vector<float>& buffer, std::vector<float>::size_type max_idx, float ratio)
{
  float threshold = buffer[max_idx] * ratio;

  std::vector<float>::size_type left_idx = max_idx;
  std::vector<float>::size_type right_idx = max_idx;

  do{left_idx--;}while(buffer[left_idx] > threshold && left_idx > 0);
  do{right_idx++;}while(buffer[right_idx] > threshold && right_idx < buffer.size());

  float left_xt = (threshold - buffer[left_idx])/(buffer[left_idx+1] - buffer[left_idx]) + left_idx;
  float right_xt = (threshold - buffer[right_idx-1])/(buffer[right_idx] - buffer[right_idx-1]) + right_idx -1;

  return (left_xt + right_xt)/2;
}

//! Calculate the positions given the profiles
void calc_pos(LASGlobalData<Avec>& profiles, LASGlobalData<double>& positions, LASGlobalData<int>& results_mask,   double ratio)
{
  LASGlobalDataLoop loop;

  do{
    Avec& buffer = loop.GetEntry<Avec>(profiles);

    Avec::size_type max_idx = lvmax(buffer);

    if(max_idx == 0)continue;

    double threshold = buffer[max_idx]*ratio;

    Avec::size_type left_idx = max_idx;
    Avec::size_type right_idx = max_idx;

    do{left_idx--;}while(buffer[left_idx] > threshold && left_idx > 0);
    do{right_idx++;}while(buffer[right_idx] > threshold && right_idx < buffer.size());

    double left_xt = (threshold - buffer[left_idx])/(buffer[left_idx+1] - buffer[left_idx]) + left_idx;
    double right_xt = (threshold - buffer[right_idx-1])/(buffer[right_idx] - buffer[right_idx-1]) + right_idx -1;

    loop.GetEntry<double>(positions) = (left_xt + right_xt)/2;
    loop.GetEntry<int>(results_mask) = 1;

  } while(loop.next());

}

//! Calculate the positions given the profiles
void calc_pos_2(LASGlobalData<Avec>& profiles, LASGlobalData<Avec>& striprms, LASGlobalData<double>& positions, LASGlobalData<double>& pos_error, LASGlobalData<int>& results_mask,   double ratio)
{
  LASGlobalDataLoop loop;

  do{
    Avec& buffer = loop.GetEntry<Avec>(profiles);
    Avec& noise_buffer = loop.GetEntry<Avec>(striprms);

    Avec::size_type max_idx = lvmax(buffer);

    if(max_idx == 0)continue;

    double threshold = buffer[max_idx]*ratio;

    Avec::size_type left_idx = max_idx;
    Avec::size_type right_idx = max_idx;

    do{left_idx--;}while(buffer[left_idx] > threshold && left_idx > 0);
    do{right_idx++;}while(buffer[right_idx] > threshold && right_idx < buffer.size());

    double left_xt = (threshold - buffer[left_idx])/(buffer[left_idx+1] - buffer[left_idx]) + left_idx;
    double right_xt = (threshold - buffer[right_idx-1])/(buffer[right_idx] - buffer[right_idx-1]) + right_idx -1;

    loop.GetEntry<double>(positions) = (left_xt + right_xt)/2;
    loop.GetEntry<int>(results_mask) = 1;

    // compute the error on the position, given the error of the profiles
    double denom_left = (buffer[left_idx] - buffer[left_idx + 1]) * (buffer[left_idx] - buffer[left_idx + 1]);
    double dxl_dsl1 = (threshold - buffer[left_idx + 1]) / denom_left;
    double dxl_dsl2 = (buffer[left_idx] - threshold) / denom_left;

    double denom_right = (buffer[right_idx] - buffer[right_idx - 1]) * (buffer[right_idx] - buffer[right_idx - 1]);
    double dxr_dsr1 = (threshold - buffer[right_idx]) / denom_right;
    double dxr_dsr2 = (buffer[right_idx - 1] - threshold ) / denom_right;

    double contr_sl1 = pow( dxl_dsl1 * noise_buffer[left_idx]    , 2);
    double contr_sl2 = pow( dxl_dsl2 * noise_buffer[left_idx + 1], 2);
    double contr_sr1 = pow( dxr_dsr1 * noise_buffer[right_idx - 1], 2);
    double contr_sr2 = pow( dxr_dsr2 * noise_buffer[right_idx]    , 2);

    loop.GetEntry(pos_error) = sqrt(contr_sr1 + contr_sr2 + contr_sl1 + contr_sl2) / 2;
    //std::cout << loop.GetEntry(pos_error) << "  ";

  } while(loop.next());

  //std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Calculate cut to separate AT from TEC events
double AtTecSeparationFactor(const Avec& signal_at, const Avec& signal_tec, bool verbose = false)
{
  double norm = vmax(signal_tec) / vmax(signal_at);

  double min_phi = vmin(atan2(signal_tec , signal_at * norm));
  double max_phi = vmax(atan2(signal_tec , signal_at * norm));
  double step_phi = (max_phi - min_phi) / 99;

  std::cout << "min_phi: " << min_phi << "  max_phi: " << max_phi << std::endl;

  std::map<int, int> map_ctr;
  std::map<int, double> map_fct;

  for(double phi = min_phi; phi <= max_phi; phi += step_phi){
    double factor = tan(phi);
    //std::cout << "Factor: " << factor << ": " << AT_count << " AT events and " << vsum(TEC_event) << " TEC events" << std::endl;
    int AT_count = (int)vsum(signal_at * norm * factor > signal_tec);

    map_ctr[AT_count] += 1;
    map_fct[AT_count] += factor;
  }

  int max_val = 0;
  int max_ctr = 0;
  for(std::map<int,int>::iterator i = map_ctr.begin(); i != map_ctr.end(); i++){
    //std::cout << "entry " << i->first << " has counted " << i->second << std::endl;
    if(i->second > max_ctr){
      max_ctr = i->second;
      max_val = i->first;
    }
  }

  double best_factor = norm * (map_fct[max_val] / map_ctr[max_val]);

  if(verbose){
    std::cout << "AtTecSeparationFactor" << std::endl;
    std::cout << "\n plateau is at " << max_val << " AT events" << std::endl;
    std::cout << "norm: " << norm << " * factorsum: " << map_fct[max_val] << " / counter: " << map_ctr[max_val] << " = " << best_factor << std::endl;
  }

  return best_factor;
}

void group_separation(const Avec& data, const Avec& mask, Avec& gaps, int nr_groups = 5)
{
  double sep_factor = 10000;
  Avec data_copy1 = avec_mask(data, mask);
  std::sort(data_copy1.begin(), data_copy1.end());
  Avec diff1 = vdiff(data_copy1);
  Avec data_copy2 = avec_mask( data_copy1, diff1 < vmax(diff1)/sep_factor );
  std::sort(data_copy2.begin(), data_copy2.end());
  Avec diff2 = vdiff(data_copy2);

  for(int i = 0; i < (nr_groups - 1); i++){
    Avec::size_type mt1 = lvmax(diff2);
    std::cout << "vmax(diff2): " << diff2[mt1] << "  from " << data_copy2[mt1-1] << " to " << data_copy2[mt1] << std::endl;
    double gap_pos = (data_copy2[mt1-1] + data_copy2[mt1])/2;
    diff2[mt1]=0;
    gaps.push_back(gap_pos);
  }
  std::sort(gaps.begin(), gaps.end());
}


/////////////////////////////////////////////////////////
// Candidates for Production
/////////////////////////////////////////////////////////

//! Calculate differences of positions 
// reference position is supplied in LASGlobalData<double> object together with a mask
// Positions are fetched from the file, the object called 'position_strips' has to be found in there 
void Calculate_Difference(const std::string& filename, LASGlobalData<double>& ref_pos, LASGlobalData<int>& ref_mask)
{

  LASGlobalData<double>& positions_strips = global_data_get<double>("positions_strips", filename);
  LASGlobalData<int>& results_mask = global_data_get<int>("results_mask", filename);

  //LASGlobalData<double> pitch;
  //make_pitch(pitch);

  LASGlobalData<int> mask = ref_mask * results_mask;
  LASGlobalData<double> diff_strips = (ref_pos - positions_strips);
  correct_signs(diff_strips);
  LASGlobalData<double> diff_mu = diff_strips * pitch() * 1000;

  TFile f(filename.c_str(), "UPDATE");

  diff_strips.Write("diff_strips");
  diff_mu.Write("diff_mu");
  mask.Write("diff_mask");

  f.Close();
}


void alpar_history(const std::vector<std::string>& file_list, const std::string& ref_file)
{
  TECReconstructor TR;
  TR.Init_common();
  TR.Init_onering();

  LASGlobalData<double>& ref_pos = global_data_get<double>("positions_strips", ref_file);
  LASGlobalData<int>& ref_mask = global_data_get<int>("results_mask", ref_file);

  std::vector<LasAlPar> alpar_list;
  Avec xvals;
  for(std::vector<std::string>::size_type i = 0; i < file_list.size(); i++){
    LASGlobalData<double> diff_rad = (ref_pos - global_data_get<double>("positions_strips",file_list[i])) * pitch() / r0();
    correct_signs(diff_rad);
    LASGlobalData<int> mask = global_data_get<int>("results_mask",file_list[i]) && ref_mask;
    Avec timev = avec_get("block_timestamp", file_list[i]);
    if (timev.empty()){
      std::cerr << "Could not find timestamps in file " << file_list[i] << std::endl;
      return;
    }
    LasAlPar alpar;
    TR.reconstruct_tec_ring(diff_rad, mask, alpar.tecp_r4, LAS::TEC_PLUS_R4);
    TR.reconstruct_tec_ring(diff_rad, mask, alpar.tecp_r6, LAS::TEC_PLUS_R6);
    TR.reconstruct_tec_ring(diff_rad, mask, alpar.tecm_r4, LAS::TEC_MINUS_R4);
    TR.reconstruct_tec_ring(diff_rad, mask, alpar.tecm_r6, LAS::TEC_MINUS_R6);

    alpar_list.push_back(alpar);
    xvals.push_back(timev[0]);
  }
  
  //int npoints = file_list.size();
  //alpar_history_draw(alpar_list, Avec(npoints,1,npoints));

  alpar_history_draw(alpar_list, xvals);

}


/////////////////////////
// Main Analysis steps //
/////////////////////////


//! Calculate evolution of positions within one run
void Run_History(const std::string& filename, bool plot)
{
  TECReconstructor TR;
  TR.Init_common();
  TR.Init_onering();

  TFile f(filename.c_str(),"UPDATE");
  Avec event_block = avec_get("event_block", f);
  Avec block_nr = avec_get("block_nr", f);
  Avec block_timestamp = avec_get("block_timestamp", f);
  unsigned int nblocks = (unsigned int)vmax(event_block) + 1;

  std::cout << "There are " << nblocks << " blocks in this run (" << block_nr.size() << ")" << std::endl;

  if(nblocks == 0) return;

  //LASGlobalData<Avec> ref_profiles = global_data_get<Avec>("profiles_0", f);
  LASGlobalData<Avec> ref_profiles = global_data_get<Avec>("profiles", f);
  LASGlobalData<double> ref_pos;
  LASGlobalData<int> ref_mask;
  calc_pos(ref_profiles, ref_pos, ref_mask);

  LASGlobalData<Avec> pos_history;
  LASGlobalData<Avec> diff_history;

  for(unsigned int i = 0; i< nblocks; i++){
    std::cout << "Processing block " << i << std::endl;
    std::ostringstream prof_name;
    std::ostringstream pos_name;
    std::ostringstream mask_name;
    std::ostringstream diff_name;
    std::ostringstream alpar_name;
    prof_name << "profiles_" << i;
    pos_name << "positions_" << i;
    mask_name << "positions_mask_" << i;
    diff_name << "diff_mm_" << i;
    alpar_name << "alpar_" << i;

    LASGlobalData<Avec>& profiles = global_data_get<Avec>(prof_name.str(), f);
    LASGlobalData<double> positions;
    LASGlobalData<int> mask;
    calc_pos(profiles, positions, mask);
    LASGlobalData<double> diff_mm = (ref_pos - positions) * pitch();
    LASGlobalData<int> diff_mask = ref_mask && mask;
    correct_signs(diff_mm);

    LASGlobalData<double> diff_rad = diff_mm / r0();
    LasAlPar alpar;
    TR.reconstruct_tec_ring(diff_rad, diff_mask, alpar.tecp_r4, LAS::TEC_PLUS_R4);
    TR.reconstruct_tec_ring(diff_rad, diff_mask, alpar.tecp_r6, LAS::TEC_PLUS_R6);
    TR.reconstruct_tec_ring(diff_rad, diff_mask, alpar.tecm_r4, LAS::TEC_MINUS_R4);
    TR.reconstruct_tec_ring(diff_rad, diff_mask, alpar.tecm_r6, LAS::TEC_MINUS_R6);


    positions.Write(pos_name.str().c_str());
    mask.Write(mask_name.str().c_str());
    alpar.Write(alpar_name.str().c_str());

    LASGlobalDataLoop loop;
    do{
      loop.GetEntry<Avec>(pos_history).push_back(loop.GetEntry<double>(positions));
      loop.GetEntry<Avec>(diff_history).push_back(loop.GetEntry<double>(diff_mm));
    }while(loop.next());

  }
  pos_history.Write("pos_history");
  diff_history.Write("diff_history");
  f.Close();

  if(plot) control_history(filename);
}

//! Calculate the positions
void Position_Calculation(const std::string& filename, const LAS::AnalysisParameters& par)
{
  std::cout << "Position Calculation" << std::endl;
  std::cout << "edge_threshold=" << par.edge_threshold << std::endl;

  // Open the file
  TFile f(filename.c_str(),"UPDATE");
  if(!f.IsOpen()) throw LAS::Exception( "Error in Position_Calculation: Could not open file " + filename );

  // Blockwise
  Avec block_nr = avec_get("block_nr", f);
  Avec::size_type nblocks = block_nr.size();
  std::cout << "There are " << nblocks << " blocks in this run" << std::endl;
  if(nblocks == 0) throw LAS::Exception( "Error in Position_Calculation: There are no blocks in this run");

  for(unsigned int i = 0; i< nblocks; i++){
    //std::cout << "Processing block " << i << std::endl;
    std::ostringstream prof_name;
    std::ostringstream srms_name;
    std::ostringstream pos_name;
    std::ostringstream err_name;
    std::ostringstream mask_name;
    prof_name << "profiles_" << i;
    srms_name << "striprms_" << i;
    pos_name << "positions_" << i;
    err_name << "pos_error_" << i;
    mask_name << "positions_mask_" << i;

    LASGlobalData<Avec>& profiles = global_data_get<Avec>(prof_name.str(), f);
    LASGlobalData<Avec>& striprms = global_data_get<Avec>(srms_name.str(), f);
    LASGlobalData<double> positions;
    LASGlobalData<double> pos_error;
    LASGlobalData<int> mask;
    //calc_pos(profiles, positions, mask, edge_threshold);
    calc_pos_2(profiles, striprms, positions, pos_error, mask, par.edge_threshold);

    positions.Write(pos_name.str().c_str());
    pos_error.Write(err_name.str().c_str());
    mask.Write(mask_name.str().c_str());
  }

  // Integrated over whole run
  LASGlobalData<Avec>* profiles = 0;
  f.GetObject("profiles",profiles);
  if(!profiles){
    f.Close();
    throw LAS::Exception( "Error in Position_Calculation: Did not find integrated profiles" );
  }
  LASGlobalData<double> positions(0);
  LASGlobalData<int> results_mask(0);

  calc_pos(*profiles, positions, results_mask, par.edge_threshold);

  positions.Write("positions_strips");
  results_mask.Write("results_mask");

  f.Close();
}

void process_candidates(const Avec2D& prof_cands, const Avec& sum_cand, Avec& theprofile, Avec& striprms, int& norm, double& optimal_threshold, double& pos_av, double& therms)
{
  if(prof_cands.empty())return;
  Avec sum_prof = sum_cand / prof_cands.size();

  // Construct a mask that surrounds the center of the laser profile
  Avec mask2(512);
  Avec::size_type idx = lvmax(sum_prof);
  double maxval = sum_prof[idx];
  double mask_threshold= 0.09;

  do{
    mask2[idx] = 1;
  }while (sum_prof[idx] >= mask_threshold * maxval && --idx > 0);
  idx = lvmax(sum_prof);
  do{
    mask2[idx] = 1;
  }while (sum_prof[idx] >= mask_threshold * maxval && ++idx < sum_prof.size());
  

  // Remove candidates with no signal in the profile range
  Avec2D  good_profs;
  Avec sum_prof2(512, 0.0);
  Avec ssq_prof2(512, 0.0);
  for(Avec::size_type i = 0; i < prof_cands.size(); i++){
    Avec profile = prof_cands[i] * mask2;
    if(vsum(profile) > 0){
      good_profs.push_back(profile);
      sum_prof2 += profile;
      ssq_prof2 += profile * profile;
    }
  }

  // Normalize the integrated profile
  Avec2D::size_type npr = good_profs.size();
  if(npr == 0){
    //std::cerr << "No good profiles left" << std::endl;
    return;
  }
  //ssq_prof2 = sqrt((ssq_prof2 - sum_prof2 * sum_prof2 / npr) / (npr - 1)/ npr);
  ssq_prof2 = sqrt((ssq_prof2 - sum_prof2 * sum_prof2 / npr) / (npr - 1));
  sum_prof2 /= npr;

  // Now perform the edge threshold optimization
  optimal_threshold = 0;
  double smallest_rms = 1e6;

  for(double edg_thr = 0.4; edg_thr <= 0.8; edg_thr += 0.05){
    double sum = 0;
    double ssq = 0;
    for(Avec2D::size_type i = 0; i < npr; i++){
      double pos = compute_position(good_profs[i].begin(), good_profs[i].end(), edg_thr);
      sum += pos;
      ssq += pos*pos;
    }
    double rms = sqrt((ssq - sum * sum/npr)/(npr-1));
    if(rms < smallest_rms){
      optimal_threshold = edg_thr;
      smallest_rms = rms;
      pos_av = sum / npr;
    }
  }
  //std::cout << "Optimal threshold is " << optimal_threshold << " with rms " << smallest_rms << std::endl;
  therms = smallest_rms < 100 ? smallest_rms : 0.0;
  theprofile = sum_prof2;
  striprms = ssq_prof2;
  norm = good_profs.size();
}

void process_block(
		   const LASGlobalData<Avec2D>& profile_candidates,
		   LASGlobalData<Avec>& sum_candidates,
		   LASGlobalData<Avec>& profiles,
		   LASGlobalData<Avec>& striprms,
		   LASGlobalData<int>& norm,
		   LASGlobalData<double>& edge_thresholds,
		   LASGlobalData<double>& pos_av,
		   LASGlobalData<double>& rms_av,
		   const std::string& output_file,
		   unsigned int current_block
		   )
{
  LASGlobalData<double> positions;

  LASGlobalDataLoop lp;
  do{
    Avec&   prof    = lp.GetEntry(profiles);
    Avec&   srms    = lp.GetEntry(striprms);
    int&    thenorm = lp.GetEntry(norm);
    double& edg_thr = lp.GetEntry(edge_thresholds);
    double& pos     = lp.GetEntry(pos_av);
    double& rms     = lp.GetEntry(rms_av);
    process_candidates(lp.GetEntry(profile_candidates), lp.GetEntry(sum_candidates), prof, srms, thenorm, edg_thr, pos, rms);
    lp.GetEntry(positions) = compute_position(prof, edg_thr);
  } while(lp.next());


  // Write results to file
  TFile of(output_file.c_str(),"UPDATE");

  std::ostringstream prof_name;
  std::ostringstream srms_name;
  std::ostringstream norm_name;
  std::ostringstream pos_name;
  std::ostringstream pos_av_name;
  std::ostringstream rms_av_name;
  std::ostringstream edg_thr_name;

  prof_name << "profiles_" << current_block;
  srms_name << "striprms_" << current_block;
  norm_name << "norm_" << current_block;
  pos_name << "positions_" << current_block;
  pos_av_name << "pos_av_" << current_block;
  rms_av_name << "rms_av_" << current_block;
  edg_thr_name << "edg_thr_" << current_block;

  profiles.Write(prof_name.str().c_str());
  striprms.Write(srms_name.str().c_str());
  norm.Write(norm_name.str().c_str());
  pos_av.Write(pos_av_name.str().c_str());
  positions.Write(pos_name.str().c_str());
  rms_av.Write(rms_av_name.str().c_str());
  edge_thresholds.Write(edg_thr_name.str().c_str());

  of.Close();

  return;
}


void process_block_old(
		   unsigned int current_block,
		   const std::vector< LASGlobalData<Avec> >& profile_candidates,
		   const LASGlobalData<Avec>& flag_candidates, 
		   LASGlobalData<Avec>& sum_candidates,
		   LASGlobalData<Avec>& profiles,
		   LASGlobalData<int>& norm,
		   std::vector<LASGlobalData<Avec> >& block_profiles,
		   std::vector<LASGlobalData<int> >& block_norm,
		   std::vector<LASGlobalData<double> >& block_sum_pos, 
		   std::vector<LASGlobalData<double> >& block_ssq_pos,
		   double edge_threshold = 0.5
		   )
{
  // Construct a mask that surrounds the center of the laser profile
  LASGlobalData<Avec> narr_mask( Avec(512, 0.0) );
  LASGlobalDataLoop maskloop;
  const double mask_threshold = 0.1;

  // Estimate good edge thresholds
  LASGlobalData<double> edg_thr(0.0);
  double amp1 = 10;
  double amp2 = 255;
  double th1 = 0.7;
  double th2 = 0.2;
  double slope = (th2 - th1) / (amp2 - amp1);

  do{
    // Normalize the integrated candidates
    Avec& sum_cand = maskloop.GetEntry(sum_candidates);
    double norm_cand = vsum(maskloop.GetEntry(flag_candidates));
    if(norm_cand != 0) sum_cand /= norm_cand;
    
    Avec& mask2 = maskloop.GetEntry(narr_mask);
    Avec::size_type idx = lvmax(sum_cand);
    double maxval = sum_cand[idx];
    do{
      mask2[idx] = 1;
    }while (sum_cand[idx] > mask_threshold * maxval && --idx > 0);
    idx = lvmax(sum_cand);
    do{
      mask2[idx] = 1;
    }while (sum_cand[idx] > mask_threshold * maxval && ++idx < sum_cand.size());

    maskloop.GetEntry(edg_thr) = slope * (maxval - amp1) + th1;
  }while (maskloop.next());

  // Go through candidates and keep only good ones
  for(unsigned i = 0; i < profile_candidates.size(); i++){
    LASGlobalDataLoop selectloop;
    do{
      // Check if this entry was flagged as candidate
      if(selectloop.GetEntry(flag_candidates)[i] != 1) continue;
      
      const Avec prof_cand = selectloop.GetEntry(profile_candidates[i]) * selectloop.GetEntry(narr_mask);
      
      // Check if there is signal within the restricted mask
      if(vsum( prof_cand ) > 0){

	// Integrate average profile of run
	selectloop.GetEntry<Avec>(profiles) += prof_cand;
	selectloop.GetEntry<int>(norm) ++;
	
	// Integrate profile for blocks
	selectloop.GetEntry<Avec>(block_profiles[current_block]) += prof_cand;
	selectloop.GetEntry<int>(block_norm[current_block]) ++;	
	
	// Single event statistics
	double the_edge_threshold = edge_threshold < 0 ? selectloop.GetEntry(edg_thr) : edge_threshold;
	double pos = compute_position( prof_cand.begin(), prof_cand.end(), the_edge_threshold);
	selectloop.GetEntry<double>(block_sum_pos.at(current_block)) += pos;
	selectloop.GetEntry<double>(block_ssq_pos.at(current_block)) += pos*pos;      
      }
    }while (selectloop.next());
  }
}

//! Integrate the events for the different modules
void Signal_Integration(const std::string& data_file, const std::string& label_file)
{
  std::cout << "Running Signal_Integration" << std::endl;
  //std::cout << "edge_threshold=" << edge_threshold << std::endl;

  std::time_t t1 = std::time(NULL);

  double compl_signal_cut = 100.0;

  // Get labels and step mask from the result file
  Avec label_step = avec_get("label_step", label_file);
  Avec block_nr = avec_get("event_block", label_file);
  LASGlobalData<int> step_mask = global_data_get<int>("step_mask", label_file);

  if(label_step.empty() || block_nr.empty()) throw LAS::Exception( "Error in Signal_Integration: could not find all Avecs (label_step, event_block)");

  Avec::size_type nr_entries = label_step.size();
  Avec::size_type nblocks = (unsigned int)vmax(block_nr) + 1; // Number of Blocks (Block numbering starts from 0)

  // Open the data file
  TFile df(data_file.c_str(),"READ");
  if(!df.IsOpen()) throw LAS::Exception( "Error in Signal_Integration: Could not open file " + data_file );

  // Find the Raw data Tree
  TTree* tree = 0;
  df.GetObject("lasRawDataTree",tree);
  if(!tree) throw LAS::Exception( "Error in Signal_Integration: Did not find lasRawDataTree" );

  // Create container for Branch content
  std::vector<float> empty_buffer(512,0);
  LASGlobalData<std::vector<float> > theData( empty_buffer );
  LASGlobalData<std::vector<float> > *ptr = &theData;
  if( tree->SetBranchAddress("lasRawData",&ptr) < 0) throw LAS::Exception( "Error in Signal_Integration: could not set branch address for lasRawData");

  /////////////////////////
  // Objects for results //
  /////////////////////////

  // Average profiles for entire run
  LASGlobalData<Avec> profiles(Avec(512,0.0));
  LASGlobalData<int> norm(0);

  // Define a strip range in which to look for the maximum
  // In this case the 100 central strips
  LASGlobalData<int> range_low, range_high;
  fill_strip_ranges(range_low, range_high);
  // Define masks for each module, which are 0 within the valid range of strips and 1 otherwise
  LASGlobalData<Avec> compl_mask = LASGlobalData<Avec>(Avec(512, 0, 511)) < range_low || LASGlobalData<Avec>(Avec(512, 0, 511)) > range_high;

  // Define masks for each module, which are 1 within the valid range of strips and 0 otherwise
  LASGlobalData<Avec> strip_mask = LASGlobalData<Avec>(Avec(512, 0, 511)) >= range_low && LASGlobalData<Avec>(Avec(512, 0, 511)) <= range_high;
  LASGlobalData<std::vector<int> > strip_mask_int(std::vector<int>(512));
  LASGlobalDataLoop cvlp;
  do{
    for(Avec::size_type idx = 0; idx < 512; idx++) cvlp.GetEntry(strip_mask_int)[idx] = (cvlp.GetEntry(strip_mask)[idx] > 0.0 ? 1 : 0);
  }while(cvlp.next());


  Avec comp_sig_tec;
  Avec comp_sig_tob;
  Avec comp_sig_tib;

  LASGlobalData<Avec2D> profile_candidates;
  LASGlobalData<Avec> sum_candidates(Avec(512,0.0));

  // Profiles per block
  std::vector<LASGlobalData<Avec> > block_profiles(nblocks, LASGlobalData<Avec>(Avec(512,0.0)));
  std::vector<LASGlobalData<Avec> > block_striprms(nblocks, LASGlobalData<Avec>(Avec(512,0.0)));
  std::vector<LASGlobalData<int> > block_norm(nblocks, 0);
  std::vector<LASGlobalData<int> > block_mask(nblocks, 0);
  std::vector<LASGlobalData<double> > block_pos_av(nblocks, 0.0);
  std::vector<LASGlobalData<double> > block_rms_av(nblocks, 0.0);
  std::vector<LASGlobalData<double> > block_edg_th(nblocks, 0.0);
  std::vector<LASGlobalData<double> > block_position(nblocks, 0.0);

  unsigned int current_block = 0;
 
  // Loop over all Entries
  for(Avec::size_type i = 0; i < nr_entries; i++){
    if(!(i%5000)) std::cout << "Processing entry " << i  << " of " << nr_entries << std::endl;

    if(block_nr[i] < 0) continue;

    unsigned int bl_nr = (unsigned int)block_nr[i];

    // Check if a new block has started
    if(bl_nr != current_block){
      process_block(
		   profile_candidates, 
		   sum_candidates,
		   block_profiles[current_block],
		   block_striprms[current_block],
		   block_norm[current_block],
		   block_edg_th[current_block],
		   block_pos_av[current_block],
		   block_rms_av[current_block],
		   label_file,
		   current_block
		   );
      profiles += block_profiles[current_block] * block_norm[current_block];
      norm += block_norm[current_block];
      std::cout << "Finished processing block " << current_block << std::endl;

      // Swap to new block
      current_block = bl_nr;
      profile_candidates.Reset();
      sum_candidates.Reset(Avec(512, 0.0));
    }

    tree->GetEntry(i);

    // Check for every module if it is a candidate for integration
    LASGlobalDataLoop loop;
    do{
      // Check if laser settings are the ones for this module
      if(loop.GetEntry<int>(step_mask) != label_step[i]) continue;

      Avec buffer(512);
      double compl_sig=0;
      // Skip empty entries
      std::vector<float>& theProfile = loop.GetEntry<std::vector<float> >(theData);

      switch(0){
      case 1:
	{std::vector<int>& mask  = loop.GetEntry(strip_mask_int);
	std::vector<float>::size_type idx;
	float maxv = 0;
	for(idx = 0; idx < theProfile.size(); idx++){
	  switch(mask[idx]){
	  case 1:
	    maxv = std::max(maxv, theProfile[idx]);
	    buffer[idx] = theProfile[idx];
	    break;
	  case 0:
	    compl_sig += theProfile[idx];
	    break;
	  }
	}
	if(maxv == 0) continue;
	break;}
      default:
	if( vector_max( theProfile , loop.GetEntry(range_low), loop.GetEntry(range_high) ) == 0 ) continue;
	buffer = loop.GetEntry(theData);
	// Check the signal outside valid range, if too high, this profile is considered bad
	compl_sig = vsum(loop.GetEntry<Avec>(compl_mask) * buffer );
      }

      if(compl_sig > compl_signal_cut){
	switch(loop.get_det()){
	case 0:
	case 1:
	  comp_sig_tec.push_back( compl_sig );
	  break;
	case 2:
	  comp_sig_tib.push_back( compl_sig );
	  break;
	case 3:
	  comp_sig_tob.push_back( compl_sig );
	  break;
	}
	continue;
      }

      loop.GetEntry(sum_candidates) += buffer;
      loop.GetEntry(profile_candidates).push_back(buffer);

    }while (loop.next());

  } // End of loop over all entries

  // Here we process the last block
  process_block(
		profile_candidates, 
		sum_candidates,
		block_profiles[current_block],
		block_striprms[current_block],
		block_norm[current_block],
		block_edg_th[current_block],
		block_pos_av[current_block],
		block_rms_av[current_block],
		label_file,
		current_block
		);
  profiles += block_profiles[current_block];
  norm += block_norm[current_block];
  std::cout << "Finished processing block " << current_block << std::endl;
  //profile_candidates.Reset();
  //sum_candidates.Reset(Avec(512, 0.0));

  df.Close();

  // Normalize the integrated profile
  LASGlobalDataLoop loop;
  do{
    if(loop.GetEntry<int>(norm) != 0) loop.GetEntry<Avec>(profiles) /= loop.GetEntry<int>(norm);
  }while (loop.next());


  // Write results to file
  TFile of(label_file.c_str(),"UPDATE");
  profiles.Write("profiles");
  comp_sig_tec.Write("comp_sig_tec");
  comp_sig_tib.Write("comp_sig_tib");
  comp_sig_tob.Write("comp_sig_tob");

//   for(unsigned int i = 0; i < nblocks; i++){

//     std::ostringstream prof_name;
//     std::ostringstream norm_name;
//     std::ostringstream pos_av_name;
//     std::ostringstream rms_av_name;

//     prof_name << "profiles_" << i;
//     norm_name << "norm_" << i;
//     pos_av_name << "pos_av_" << i;
//     rms_av_name << "rms_av_" << i;

//     block_profiles[i].Write(prof_name.str().c_str());
//     block_norm[i].Write(norm_name.str().c_str());
//     block_pos_av[i].Write(pos_av_name.str().c_str());
//     block_rms_av[i].Write(rms_av_name.str().c_str());
//   }
  of.Close();

  std::time_t t2 = std::time(NULL);
  std::time_t t3 = t2 - t1;
  std::cout << "SignalIntegration took " << t3 << " seconds" << std::endl;

}

//! Integrate the events for the different modules
void Signal_Integration_old(const std::string& data_file, const std::string& label_file, double edge_threshold)
{
  std::cout << "Running Signal_Integration" << std::endl;
  std::cout << "edge_threshold=" << edge_threshold << std::endl;

  std::time_t t1 = std::time(NULL);

  double compl_signal_cut = 0.0;

  // Get labels and step mask from the result file
  Avec label_step = avec_get("label_step", label_file);
  Avec block_nr = avec_get("event_block", label_file);
  LASGlobalData<int> step_mask = global_data_get<int>("step_mask", label_file);

  if(label_step.empty() || block_nr.empty()) throw LAS::Exception("Error in SignalIntegration: could not find all Avecs (label_step, event_block)");

  Avec::size_type nr_entries = label_step.size();
  Avec::size_type nblocks = (unsigned int)vmax(block_nr) + 1; // Number of Blocks (Block numbering starts from 0)

  // Open the data file
  TFile df(data_file.c_str(),"READ");
  if(!df.IsOpen()) throw LAS::Exception( "Could not open file " + data_file );

  // Find the Raw data Tree
  TTree* tree = 0;
  df.GetObject("lasRawDataTree",tree);
  if(!tree) throw LAS::Exception( "Error in Signal_Integration: Did not find lasRawDataTree" );

  // Create container for Branch content
  std::vector<float> empty_buffer(512,0);
  LASGlobalData<std::vector<float> > theData( empty_buffer );
  LASGlobalData<std::vector<float> > *ptr = &theData;
  Int_t retval = tree->SetBranchAddress("lasRawData",&ptr);
  if(retval < 0) throw LAS::Exception( "Error in SignalIntegration: Could not set branch address for lasRawData");

  /////////////////////////
  // Objects for results //
  /////////////////////////

  // Average profiles for entire run
  LASGlobalData<Avec> profiles(Avec(512,0.0));
  LASGlobalData<int> norm(0);

  // Profiles per block
  std::vector<LASGlobalData<Avec> > block_profiles(nblocks, LASGlobalData<Avec>(Avec(512,0.0)));
  std::vector<LASGlobalData<int> > block_norm(nblocks, 0);

  std::vector<LASGlobalData<double> > block_sum_pos(nblocks, 0.0); // sum of positions from all events corresponding to a block 
  std::vector<LASGlobalData<double> > block_ssq_pos(nblocks, 0.0); // sum of squares of positions from all events corresponding to a block 
  std::vector<LASGlobalData<double> > block_rms_pos(nblocks, 0.0); // rms of positions corresponding to a block 
  std::vector<LASGlobalData<Avec> > block_norm_strip_sum(nblocks, LASGlobalData<Avec>(Avec(512,0.0)));
  std::vector<LASGlobalData<Avec> > block_norm_strip_ssq(nblocks, LASGlobalData<Avec>(Avec(512,0.0)));

  // Define a strip range in which to look for the maximum
  // In this case the 100 central strips
  LASGlobalData<int> range_low, range_high;
  fill_strip_ranges(range_low, range_high);
  // Define masks for each module, which are 0 within the valid range of strips and 1 otherwise
  LASGlobalData<Avec> compl_mask = LASGlobalData<Avec>(Avec(512, 0, 511)) < range_low || LASGlobalData<Avec>(Avec(512, 0, 511)) > range_high;

  // Define masks for each module, which are 1 within the valid range of strips and 0 otherwise
  LASGlobalData<Avec> strip_mask = LASGlobalData<Avec>(Avec(512, 0, 511)) >= range_low && LASGlobalData<Avec>(Avec(512, 0, 511)) <= range_high;

  Avec comp_sig_tec;
  Avec comp_sig_tob;
  Avec comp_sig_tib;

  std::vector< LASGlobalData<Avec> > profile_candidates;
  LASGlobalData<Avec> flag_candidates;
  LASGlobalData<Avec> sum_candidates(Avec(512,0.0));
  //LASGlobalData<int> ctr_candidates(0);
  unsigned int current_block = 0;
 
  // Loop over all Entries
  for(Avec::size_type i = 0; i < nr_entries; i++){
    if(!(i%5000)) std::cout << "Processing entry " << i  << " of " << nr_entries << std::endl;

    if(block_nr[i] < 0) continue;

    unsigned int bl_nr = (unsigned int)block_nr[i];

    // Check if a new block has started
    if(bl_nr != current_block){
      process_block_old(
		   current_block,
		   profile_candidates,
		   flag_candidates, 
		   sum_candidates,
		   profiles,
		   norm,
		   block_profiles,
		   block_norm,
		   block_sum_pos, 
		   block_ssq_pos,
		   edge_threshold
		   );
      std::cout << "Finished processing block " << current_block << std::endl;

      // Swap to new block
      current_block = bl_nr;
      profile_candidates.clear();
      flag_candidates.Reset();
      sum_candidates.Reset(Avec(512, 0.0));
    }

    if(profile_candidates.size() >= 3100) continue; // Upper limit for block size to avoid memory issues
    tree->GetEntry(i);
    profile_candidates.push_back( LASGlobalData<Avec>()); // Empty Entry to hold the profiles
    LASGlobalData<Avec>& prof_cand = profile_candidates.back();

    // Check for every module if it is a candidate for integration
    LASGlobalDataLoop loop;
    do{
      loop.GetEntry( flag_candidates ).push_back(0); // By default this is no candidate

      // Check if laser settings are the ones for this module
      if(loop.GetEntry<int>(step_mask) != label_step[i]) continue;
      // Skip empty entries
      std::vector<float>& theProfile = loop.GetEntry<std::vector<float> >(theData);
      if( vector_max( theProfile , loop.GetEntry(range_low), loop.GetEntry(range_high) ) == 0 ) continue;

      loop.GetEntry(prof_cand) = loop.GetEntry<std::vector<float> >(theData); // Copy the profile
      Avec& buffer = loop.GetEntry(prof_cand);

      // Check the signal outside valid range, if too high, this profile is considered bad	
      double compl_sig = vsum(loop.GetEntry<Avec>(compl_mask) * buffer );

      if(compl_sig > compl_signal_cut){
	switch(loop.get_det()){
	case 0:
	case 1:
	  comp_sig_tec.push_back( compl_sig );
	  break;
	case 2:
	  comp_sig_tib.push_back( compl_sig );
	  break;
	case 3:
	  comp_sig_tob.push_back( compl_sig );
	  break;
	}
	continue;
      }

      loop.GetEntry(flag_candidates).back() = 1;
      loop.GetEntry(sum_candidates) += buffer;

    }while (loop.next());

  } // End of loop over all entries

  // Here we process the last block
  process_block_old(
		current_block,
		profile_candidates,
		flag_candidates, 
		sum_candidates,
		profiles,
		norm,
		block_profiles,
		block_norm,
		block_sum_pos, 
		block_ssq_pos,
		edge_threshold
		);
  std::cout << "Finished processing block " << current_block << std::endl;
  profile_candidates.clear();
  flag_candidates.Reset();
  sum_candidates.Reset(Avec(512, 0.0));

  df.Close();

  // Normalize the intgrated profile
  LASGlobalDataLoop loop;
  do{
    if(loop.GetEntry<int>(norm) != 0) loop.GetEntry<Avec>(profiles) /= loop.GetEntry<int>(norm);
  }while (loop.next());


  // Write results to file

  TFile of(label_file.c_str(),"UPDATE");
  profiles.Write("profiles");
  comp_sig_tec.Write("comp_sig_tec");
  comp_sig_tib.Write("comp_sig_tib");
  comp_sig_tob.Write("comp_sig_tob");

  for(unsigned int i = 0; i < nblocks; i++){

    // Normalize the intgrated block profiles and compute positions and rms for each block
    LASGlobalDataLoop loop;
    do{
      int bl_norm =  loop.GetEntry<int>(block_norm[i]);
      if(bl_norm > 1){
	loop.GetEntry<Avec>(block_profiles[i]) /= bl_norm;
	loop.GetEntry<double>(block_rms_pos[i]) = sqrt((loop.GetEntry<double>(block_ssq_pos[i]) - (loop.GetEntry<double>(block_sum_pos[i]) * loop.GetEntry<double>(block_sum_pos[i])) / bl_norm) / (bl_norm -1));
	loop.GetEntry<double>(block_sum_pos[i]) /= bl_norm;
      }
    }while (loop.next());

    std::ostringstream prof_name;
    std::ostringstream norm_name;
    std::ostringstream pos_av_name;
    std::ostringstream rms_av_name;

    prof_name << "profiles_" << i;
    norm_name << "norm_" << i;
    pos_av_name << "pos_av_" << i;
    rms_av_name << "rms_av_" << i;

    block_profiles[i].Write(prof_name.str().c_str());
    block_norm[i].Write(norm_name.str().c_str());
    block_sum_pos[i].Write(pos_av_name.str().c_str());
    block_rms_pos[i].Write(rms_av_name.str().c_str());
  }
  of.Close();

  std::time_t t2 = std::time(NULL);
  std::time_t t3 = t2 - t1;
  std::cout << "SignalIntegration took " << t3 << " seconds" << std::endl;

}

//! Integrate the events for the different modules (obsolete version)
void Signal_Integration_very_old(const std::string& data_file, const std::string& label_file)
{
  std::cout << "Running Signal_Integration" << std::endl;

  double compl_signal_cut = 1.0;

  // Get labels and step mask from the result file
  Avec label_step = avec_get("label_step", label_file);
  Avec block_nr = avec_get("event_block", label_file);
  Avec unixTime = avec_get("unixTime", label_file);
  LASGlobalData<int> step_mask = global_data_get<int>("step_mask", label_file);

  if(label_step.empty() || block_nr.empty() || unixTime.empty()){
    std::cerr << "Error, could not find all Avecs (label_step, block_nr and unixTime)" << std::endl;
    return;
  }
  Avec::size_type nr_entries = label_step.size();
  Avec::size_type nblocks = (unsigned int)vmax(block_nr) + 1; // Number of Blocks (Block numbering starts from 0)

  // Open the data file
  TFile df(data_file.c_str(),"READ");
  if(!df.IsOpen()){
    std::cerr << "Could not open file " << data_file << std::endl;
    return;
  }

  // Find the Raw data Tree
  TTree* tree = 0;
  df.GetObject("lasRawDataTree",tree);
  if(!tree){
    std::cerr << "Did not find lasRawDataTree" << std::endl;
    return;
  }

  // Create container for Branch content
  std::vector<float> empty_buffer(512,0);
  LASGlobalData<std::vector<float> > theData( empty_buffer );
  LASGlobalData<std::vector<float> > *ptr = &theData;
  tree->SetBranchAddress("lasRawData",&ptr);

  /////////////////////////
  // Objects for results //
  /////////////////////////

  // Average profiles for entire run
  LASGlobalData<Avec> profiles(Avec(512,0.0));
  LASGlobalData<int> norm(0);

  // Profiles per block
  std::vector<LASGlobalData<Avec> > block_profiles(nblocks, LASGlobalData<Avec>(Avec(512,0.0)));
  std::vector<LASGlobalData<int> > block_norm(nblocks, 0);

  std::vector<LASGlobalData<double> > block_sum_pos(nblocks, 0.0); // sum of positions from all events corresponding to a block 
  std::vector<LASGlobalData<double> > block_ssq_pos(nblocks, 0.0); // sum of squares of positions from all events corresponding to a block 
  std::vector<LASGlobalData<double> > block_rms_pos(nblocks, 0.0); // rms of positions corresponding to a block 
  std::vector<LASGlobalData<Avec> > block_norm_strip_sum(nblocks, LASGlobalData<Avec>(Avec(512,0.0)));
  std::vector<LASGlobalData<Avec> > block_norm_strip_ssq(nblocks, LASGlobalData<Avec>(Avec(512,0.0)));

  // Define a strip range in which to look for the maximum
  // In this case the 60 central strips
  LASGlobalData<int> range_low, range_high;
  fill_strip_ranges(range_low, range_high);
  // Define masks for each module, which are 0 within the valid range of strips and 1 otherwise
  LASGlobalData<Avec> compl_mask = LASGlobalData<Avec>(Avec(512, 0, 511)) < range_low || LASGlobalData<Avec>(Avec(512, 0, 511)) > range_high;

  // Define masks for each module, which are 1 within the valid range of strips and 0 otherwise
  LASGlobalData<Avec> strip_mask = LASGlobalData<Avec>(Avec(512, 0, 511)) >= range_low && LASGlobalData<Avec>(Avec(512, 0, 511)) <= range_high;

  Avec comp_sig_acc_tec;
  Avec comp_sig_rej_tec;
  Avec comp_sig_acc_tob;
  Avec comp_sig_rej_tob;
  Avec comp_sig_acc_tib;
  Avec comp_sig_rej_tib;

  //nr_entries = 10000;
  for(Avec::size_type i = 0; i < nr_entries; i++){
    if(!(i%1000)) std::cout << "Processing entry " << i  << " of " << nr_entries << std::endl;

    if(block_nr[i] < 0) continue;
    unsigned int bl_nr = (unsigned int)block_nr[i];

    tree->GetEntry(i);

    LASGlobalDataLoop loop;
    do{
      if(loop.GetEntry<int>(step_mask) == label_step[i]){
	//unsigned int max_idx = vector_max_idx(loop.GetEntry<std::vector<float> >(theData), loop.GetEntry<int>(range_low), loop.GetEntry<int>(range_high));
	Avec buffer = loop.GetEntry<std::vector<float> >(theData);

	// Check the signal outside valid range, if too high, this profile is considered bad	
	double compl_sig = vsum(loop.GetEntry<Avec>(compl_mask) * buffer );
	if(compl_sig > compl_signal_cut){
	  switch(loop.get_det()){
	  case 0:
	  case 1:
	    comp_sig_rej_tec.push_back( compl_sig );
	    break;
	  case 2:
	    comp_sig_rej_tib.push_back( compl_sig );
	    break;
	  case 3:
	    comp_sig_rej_tob.push_back( compl_sig );
	    break;
	  }
	  continue;
	}
	
	switch(loop.get_det()){
	case 0:
	case 1:
	  comp_sig_acc_tec.push_back( compl_sig );
	  break;
	case 2:
	  comp_sig_acc_tib.push_back( compl_sig );
	  break;
	case 3:
	  comp_sig_acc_tob.push_back( compl_sig );
	  break;
	}

	buffer *= loop.GetEntry(strip_mask);
	// Integrate average profile of run
	loop.GetEntry<Avec>(profiles) += buffer;
	loop.GetEntry<int>(norm) ++;

	// Integrate profile for blocks
	loop.GetEntry<Avec>(block_profiles[bl_nr]) += buffer;
	loop.GetEntry<int>(block_norm[bl_nr]) ++;

	//std::vector<float>::iterator lo_it = buffer.begin() + loop.GetEntry<int>(range_low);
	//std::vector<float>::iterator hi_it = buffer.begin() + loop.GetEntry<int>(range_high);
	//double pos = compute_position(lo_it, hi_it, 0.7);

	double pos = compute_position(buffer.begin(), buffer.end(), 0.7);
	loop.GetEntry<double>(block_sum_pos.at(bl_nr)) += pos;
	loop.GetEntry<double>(block_ssq_pos.at(bl_nr)) += pos*pos;

      }
      
    }while (loop.next());
  } // End of loop over all entries
  df.Close();

  // Normalize the intgrated profile
  LASGlobalDataLoop loop;
  do{
    if(loop.GetEntry<int>(norm) != 0) loop.GetEntry<Avec>(profiles) /= loop.GetEntry<int>(norm);
  }while (loop.next());

  TFile of(label_file.c_str(),"UPDATE");
  profiles.Write("profiles");
  comp_sig_acc_tec.Write("comp_sig_acc_tec");
  comp_sig_rej_tec.Write("comp_sig_rej_tec");
  comp_sig_acc_tib.Write("comp_sig_acc_tib");
  comp_sig_rej_tib.Write("comp_sig_rej_tib");
  comp_sig_acc_tob.Write("comp_sig_acc_tob");
  comp_sig_rej_tob.Write("comp_sig_rej_tob");

  for(unsigned int i = 0; i < nblocks; i++){

    //LASGlobalData<double> pos_prof; // Positions calculated from block profiles

    // Normalize the intgrated block profiles and compute positions and rms for each block
    LASGlobalDataLoop loop;
    do{
      //Avec& profile = loop.GetEntry(block_profiles);
      //loop.GetEntry(pos_prof) = compute_position(profile.begin(), profile.end(), 0.7);

      int bl_norm =  loop.GetEntry<int>(block_norm[i]);
      if(bl_norm > 1){
	loop.GetEntry<Avec>(block_profiles[i]) /= bl_norm;
	loop.GetEntry<double>(block_rms_pos[i]) = sqrt((loop.GetEntry<double>(block_ssq_pos[i]) - (loop.GetEntry<double>(block_sum_pos[i]) * loop.GetEntry<double>(block_sum_pos[i])) / bl_norm) / (bl_norm -1));
	loop.GetEntry<double>(block_sum_pos[i]) /= bl_norm;
      }
    }while (loop.next());
    std::ostringstream prof_name;
    std::ostringstream norm_name;
    std::ostringstream pos_av_name;
    std::ostringstream rms_av_name;
    //std::ostringstream pos_prof_name;
    prof_name << "profiles_" << i;
    norm_name << "norm_" << i;
    pos_av_name << "pos_av_" << i;
    rms_av_name << "rms_av_" << i;
    //pos_prof_name << "pos_prof_" << i;
    block_profiles[i].Write(prof_name.str().c_str());
    block_norm[i].Write(norm_name.str().c_str());
    block_sum_pos[i].Write(pos_av_name.str().c_str());
    block_rms_pos[i].Write(rms_av_name.str().c_str());
    //pos_prof.Write(pos_prof_name.str().c_str());
  }
  of.Close();

}

//! Dump the step_mask from an ASCII file into the root file as LASGlobalData<int> object
void replace_step_mask(const std::string& mask_filename, const std::string& result_filename)
{
  bool headers = false;
  LASGlobalData<int> read;

  std::ifstream in(mask_filename.c_str());
  if(!in) throw LAS::Exception( "Error in replace_step_mask: Could not open file " + mask_filename );
  read_global_data(read, LAS::ALL, headers, in);

  TFile f(result_filename.c_str(),"UPDATE");
  if( !f.IsOpen()) throw LAS::Exception( "Error in replace_step_mask: Could not open output file " + result_filename );
  read.Write("step_mask");
  f.Close();
}


//////////////////////////////////////////////////////////////////////////////////
//! Check the blocs and intensity labels and mark blocks, that are not complete or have wrong labelling
void Block_Label_Check(const std::string& result_file, Avec::size_type expected_block_size)
{
  Avec good_event = avec_get("good_event", result_file);
  Avec good_run = avec_get("good_run", result_file);

  if(good_event.empty()) throw LAS::Exception( "Error in Block_Label_Check, could not find Avec 'good_event'" );
  if(good_run.empty()) throw LAS::Exception( "Error in Block_Label_Check, could not find Avec 'good_run'" );

  Avec::size_type nr_total = good_event.size();

  Avec block_nr = avec_get("block_nr", result_file);
  Avec block_timestamp = avec_get("block_timestamp", result_file);
  Avec block_size = avec_get("block_size", result_file);

  if(block_nr.empty() || block_timestamp.empty() || block_size.empty()) throw LAS::Exception( "Error in Block_Label_Check, could not find all Avecs" );

  Avec label_step = avec_get("label_step", result_file);
  Avec event_block = avec_get("event_block", result_file);
  if(label_step.empty() || event_block.empty()) throw LAS::Exception( "Error in Block_Label_Check, could not find all Avecs" );

  Avec::size_type nr_blocks = block_nr.size();
  Avec good_block(nr_blocks, 1);

  for(Avec::size_type i = 0; i < nr_blocks; i++){
    if(block_size[i] != expected_block_size){
      std::cout << "block " << block_nr[i] << " has " << block_size[i] << " events instead of " << expected_block_size << std::endl;
      good_block[i] = 0;
    }
    else{
      Avec labels_block;
      for(Avec::size_type g = 0; g < nr_total; g++){
	if(event_block[g] == block_nr[i]){
	  labels_block.push_back(label_step[g]);
	}
      }
      if(labels_block.size() != expected_block_size){
	std::cerr << "Error in Block_Label_Check: Inconsistency between block_size and number of labelled events in event_block" << std::endl;
	good_block[i] = 0;
      }
      else{
	for(Avec::size_type h = 0; h < expected_block_size; h++){
	  if(labels_block[h] != (h%10) + 1){
	    std::cerr << "Error in Block_Label_Check: Incorrect label sequence" << std::endl;
	    good_block[i] = 0;
	  }
	}
      }
    }
  }

  bool retval = true;
  std::cout << vsum(good_block) << " good blocks out of " << nr_blocks << std::endl;
  if(vsum(good_block) == 0){
    std::cerr << "Warning: No good blocks in this run" << std::endl;
    retval = false;
  }

  good_run *= (retval ? 1 : 0);
  TFile f(result_file.c_str(), "UPDATE");
  good_block.Write("good_block");
  good_run.Write("good_run");
  f.Close();

  return;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Create Labels for the events of a run
void Event_Labels(const std::string& filename, bool plot)
{
  std::cout << "Labelling the Events" << std::endl;

  Avec good_event = avec_get("good_event",filename);
  Avec signal_tec_sep = avec_get("signal_tec_r4_nooverlap",filename);
  Avec signal_at_sep = avec_get("signal_at_nooverlap",filename);


  Avec signal_at = avec_get("signal_at",filename);
  Avec signal_tec_r4 = avec_get("signal_tec_r4",filename);
  Avec signal_tec_r6 = avec_get("signal_tec_r6",filename);
  Avec signal_tec_disc1 = avec_get("signal_tec_disc1",filename);
  Avec signal_tec_far = avec_get("signal_tec_far",filename);


  if(good_event.empty() || 
     signal_tec_sep.empty() ||
     signal_at_sep.empty() || 
     signal_at.empty() || 
     signal_tec_r4.empty() || 
     signal_tec_r6.empty() || 
     signal_tec_disc1.empty() || 
     signal_tec_far.empty()
     ) throw LAS::Exception("Could not get all Avecs");


  // AT / TEC Separation
  signal_tec_sep *= good_event;
  signal_at_sep *= good_event;

  double best_factor = AtTecSeparationFactor(signal_at_sep, signal_tec_sep, true);

  Avec AT_event = signal_at_sep * best_factor > signal_tec_sep;
  Avec TEC_event = -AT_event + 1;

  AT_event *= good_event;
  TEC_event *= good_event;

  // Report the results of labelling
  std::cout << vsum(AT_event) << " AT events and " << vsum(TEC_event) << " TEC events" << std::endl;

  if(plot){
    new TCanvas("AT_TEC_separation","AT / TEC Separation");
    TGraph * gr = avec_draw(signal_at_sep, signal_tec_sep, "Event Type Separation", "AT signal", "TEC signal", "AP");
    //gr->GetXaxis()->SetRangeUser(0, 1.1);
    //gr->GetYaxis()->SetRangeUser(0, 1.1);
    Double_t max = std::min(gr->GetXaxis()->GetXmax(), gr->GetYaxis()->GetXmax() / best_factor );
    TLine* l1 = new TLine(0, 0, max, max * best_factor );
    l1->SetLineColor(2);
    l1->SetLineWidth(3);
    l1->SetLineStyle(2);
    l1->Draw();
  }

  //////////////////////
  // Intensity Labels //
  //////////////////////

  Avec::size_type nr_entries = AT_event.size();

  // 1. Alignment Tube Events

  // Normalize
  Avec s1 = signal_at / vmax(signal_at * AT_event);
  Avec s2 = signal_tec_r6 / vmax(signal_tec_r6 * AT_event);
  Avec s3 = signal_tec_disc1 / vmax(signal_tec_disc1 * AT_event);
  // Find cuts between intensity groups
  Avec dist_at = sqrt(s1*s1 + s2*s2 + s3*s3);
  Avec cuts_at;
  group_separation(dist_at, AT_event, cuts_at, 5);

  Avec label_step(nr_entries, 0.0);
  std::vector<Avec> masks;
  double prev_cut = 0;
  for(Avec::size_type ci = 0; ci < cuts_at.size(); ci++){
    Avec mask = dist_at > prev_cut && dist_at < cuts_at[ci] && AT_event != 0;
    masks.push_back(mask);
    label_step += mask * (ci+1);
    prev_cut = cuts_at[ci];
  }
  Avec mask = dist_at > prev_cut && AT_event != 0;
  masks.push_back(mask);
  label_step += mask * (cuts_at.size() + 1);

  label_step += AT_event * (cuts_at.size() + 1);

  // 2. TEC Events

  // Normalize
  Avec t1 = signal_tec_r4 / vmax(signal_tec_r4 * TEC_event);
  Avec t2 = signal_tec_r6 / vmax(signal_tec_r6 * TEC_event);
  Avec t3 = signal_tec_far / vmax(signal_tec_far * TEC_event);
  // Find cuts between intensity groups
  Avec dist_tec = sqrt(t1*t1 + t2*t2 + t3*t3);
  Avec cuts_tec;
  group_separation(dist_tec, TEC_event, cuts_tec, 5);

  prev_cut = 0;
  for(Avec::size_type ci = 0; ci < cuts_tec.size(); ci++){
    Avec mask = dist_tec > prev_cut && dist_tec < cuts_tec[ci] && TEC_event != 0;
    masks.push_back(mask);
    label_step += mask * (ci+1);
    prev_cut = cuts_tec[ci];
  }
  Avec mask2 = dist_tec > prev_cut && TEC_event != 0;
  masks.push_back(mask2);
  label_step += mask2 * (cuts_tec.size() + 1);

//   Avec label_step(nr_entries, 0.0);
//   std::vector<Avec> masks;
//   for(Avec::size_type i = 0; i < nr_levels; i++){
//     Avec inset =          signal_at > at_lo[i]    &&         signal_at < at_hi[i] 
//                 &&   signal_tecplus > tecp_lo[i]  &&    signal_tecplus < tecp_hi[i] 
//                 &&  signal_tecminus > tecm_lo[i]  &&   signal_tecminus < tecm_hi[i]
//                 && signal_tec_disc1 > tecd1_lo[i] &&  signal_tec_disc1 < tecd1_hi[i];
//     std::cout << "Level " << i << "  inset: " << vsum(inset) << " events" << std::endl;
//     masks.push_back(inset);
//     label_step += inset * (i+1);
//   }

//   // Report the results of labelling
//   std::cout << vsum(AT_event) << " AT events and " << vsum(TEC_event) << " TEC events" << std::endl;

  Avec overlap(signal_at.size(), 0.0);
  for(Avec::size_type i = 0; i < masks.size(); i++) overlap += masks[i];

  std::cout << "overlap: " << vsum(overlap > 1) << "  unassigned: " << vsum(overlap == 0) << std::endl;

  Avec overlap_TEC = overlap * TEC_event;
  std::cout << "TEC assigned: " << vsum(overlap_TEC > 0) << "  unassigned: " << vsum(overlap_TEC == 0) << std::endl;

  Avec overlap_AT = overlap * AT_event;
  std::cout << "AT assigned: " << vsum(overlap_AT > 0) << "  unassigned: " << vsum(overlap_AT == 0) << std::endl;



  // Save the results
  TFile f(filename.c_str(), "UPDATE");
  AT_event.Write("AT_event");
  TEC_event.Write("TEC_event");
  dist_at.Write("dist_at");
  cuts_at.Write("cuts_at");
  dist_tec.Write("dist_tec");
  cuts_tec.Write("cuts_tec");

  label_step.Write("label_step");
  for(Avec::size_type i = 0; i < masks.size(); i++){
    std::ostringstream vecname;
    vecname << "step_" << (i+1) ;
    masks[i].Write(vecname.str().c_str());
  }
  f.Close();

}

//! Determine the blocks of 2000 events that form 1 snapshot
// (This needs to be revisited)
void Block_Slicer(const std::string& output_filename, bool plot, double factor, double min_block_size)
{
  std::cout << "Running Block Slicer" << std::endl;

  Avec good_event = avec_get("good_event", output_filename);
  Avec eventnumber = avec_get("eventnumber", output_filename);
  Avec lumiblock = avec_get("lumiblock", output_filename);
  Avec unixTime = avec_get("unixTime", output_filename);
  if(good_event.empty() || eventnumber.empty() || unixTime.empty()) throw LAS::Exception("Error in Block_Slicer: could not find all Avecs");

  // Vector for event differences
  Avec event_diff;
  double last_event = 0;
  double diff;

  // First loop to generate event nr. differences
  //std::cout << "Generating Event Nr. differences" << std::endl;
  Avec::size_type size = eventnumber.size();
  for(Avec::size_type i = 0; i < size; i++){
    //if(i%10000 == 0) std::cout << "processing entry Nr. " << i << std::endl;

    diff = 0;
    if(good_event[i] > 0){
      diff = eventnumber[i] - last_event;
      last_event = eventnumber[i];
    }
    event_diff.push_back(diff);
  }
  if(! event_diff.empty())event_diff[0]=0;
  //std::cout << "Computed " << event_diff.size() << " differences" << std::endl;

  // Estimate the typical difference between event numbers and define a threshold for block changes
  Avec sorted_diff = avec_mask(event_diff, good_event);
  //std::cout << "Sorting " << event_diff.size() << " differences" << std::endl;

  std::sort(sorted_diff.begin(), sorted_diff.end());
  double estimated_diff = sorted_diff.at(sorted_diff.size()/2);
  int threshold = (int)(factor * estimated_diff);
  std::cout << "Estimated event diff: " << estimated_diff << "   threshold: " << threshold << std::endl;
  //threshold = 1e6;

  // Assign block numbers
  std::cout << "Generating block numbers" << std::endl;
  Int_t block_ctr = 0;
  Avec event_block(size, -1);
  Avec::size_type idx;
  for(idx = 0 ; idx< size - 1; idx++){
    if(good_event[idx] > 0){
      if(event_diff.at(idx) < threshold){
	event_block[idx] = block_ctr;
      }
      else{
	if(event_diff.at(idx+1) < threshold){
	  block_ctr++;
	  event_block[idx] = block_ctr;
	}
      }
    }
  }
  // Process last entry separately
  if(event_diff.back() < threshold && good_event[idx]){
    event_block.back() = block_ctr;
  }

  // Eliminate blocks with small size;
  //std::cout << "Eliminating small blocks" << std::endl;
  Avec correction(size);
  for(int bc = 0; bc <= block_ctr; bc++){
    double block_size = vsum(event_block == bc);
    if(block_size < min_block_size){
      // Ugly vector arithmetics ...
      event_block -= (event_block == bc)*(bc+1);
      correction += (event_block > bc);
    }
  }  
  event_block -= correction;

  // Generate block overview
  Avec block_nr, block_timestamp, block_size, block_firstlumi, block_lastlumi;
  // Print block sizes
  std::cout << std::endl << "Block Sizes: " << std::endl;
  for(int bc = -1; bc <= block_ctr; bc++){
    double bsize = vsum(event_block == bc);
    std::cout << "Block Nr. " << bc << " has " << bsize << " entries, timestamp: " << vmax(unixTime * (event_block == bc)) << std::endl;
    if(bc >= 0 && vsum(event_block == bc) > 0){
      block_nr.push_back(bc);
      block_timestamp.push_back(vmax(unixTime * (event_block == bc)));
      block_firstlumi.push_back( vmin( avec_mask(lumiblock, (event_block == bc)) ) );
      block_lastlumi.push_back( vmax( avec_mask(lumiblock, (event_block == bc)) ) );
      block_size.push_back( bsize );
    }
  }

  // Write results to file
  TFile outfile(output_filename.c_str(),"UPDATE");
  event_diff.Write("event_diff");
  event_block.Write("event_block");
  block_nr.Write("block_nr");
  block_timestamp.Write("block_timestamp");
  block_size.Write("block_size");
  block_firstlumi.Write("block_firstlumi");
  block_lastlumi.Write("block_lastlumi");
  Avec(1, threshold).Write("block_threshold");
  outfile.Close();

  // Plot the data
  if(plot)control_blocks(output_filename, true);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Try to spot events that are only noise (this version works for good LAS settings since run 161439)
void Noise_Filter(const std::string& filename, const LAS::AnalysisParameters& par )
{
  std::cout << "Running Noise Filter" << std::endl;

  TFile f(filename.c_str(),"UPDATE");
  if(!f.IsOpen()) throw LAS::Exception("Could not open file " + filename);

  Avec signal_all = avec_get("signal_all",f);
  Avec eventnumber = avec_get("eventnumber",f);
  Avec signal_tec_r4_nooverlap = avec_get("signal_tec_r4_nooverlap", f);
  Avec signal_at_nooverlap = avec_get("signal_at_nooverlap", f);

  if( signal_all.empty() ) throw LAS::Exception("Error in Noise_Filter: Could not find Avec signal_all");
  //if( eventnumber.empty() ) throw LAS::Exception("Error in Noise_Filter: Could not find Avec eventnumber");

  Avec::size_type nr_total = signal_all.size();

  Avec good_event;
  switch(par.noise_algo_type){
  case 2:
    good_event = fabs(signal_tec_r4_nooverlap - signal_at_nooverlap/3*2) > par.noise_signal_thresh;
    break;
  case 3:
    {Avec preselect =  signal_all > par.noise_signal_thresh * 10.0;
    good_event = fabs(signal_tec_r4_nooverlap - signal_at_nooverlap/3*2) * preselect > par.noise_signal_thresh;
    break;}
  default:
    good_event = signal_all > par.noise_signal_thresh;
  }
  //good_event = good_event && (eventnumber > 500000);

  int nr_good = (int)vsum(good_event);
  int nr_bad  = nr_total - nr_good;
  double bad_ratio = (double)nr_bad/nr_good;

  bool good_run = bad_ratio < par.bad_ratio_cut;

  std::cout << nr_good << " good events, " << nr_bad << " bad events, " << nr_total << " total" << std::endl;
  std::cout << "bad events / good events = " << bad_ratio << std::endl;

  good_event.Write("good_event");
  Avec(1, good_run ? 1 : 0 ).Write("good_run");
  Avec(1, par.noise_signal_thresh).Write("signal_thresh");
  Avec(1, par.bad_ratio_cut).Write("bad_ratio_cut");
  Avec(1, bad_ratio).Write("bad_ratio");

  f.Close();
  if( nr_good < 2000) throw  LAS::Exception("Exception in noise_Filter: There are less than 2000 good events in this run");
}

//! Retrieve the signal maxima for all modules and all events
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_max_signal(const std::string& filename, const std::string& output_rootfilename, bool plot, bool full_output, Long64_t nentries)
{
  std::cout << "Running get_max_signal" << std::endl;
  std::time_t start_time = std::time(NULL);

  // Open the data file
  TFile f(filename.c_str());
  if(!f.IsOpen()) throw LAS::Exception("Could not open file " + filename);

  // Find the Tree
  TTree* tree = 0;
  f.GetObject("lasRawDataTree",tree);
  if( !tree ) throw LAS::Exception("Did not find rawDataTree in " + filename);    

  // Create container for Branch content
  std::vector<float> empty_buffer(512,0);
  LASGlobalData<std::vector<float> > theData( empty_buffer );
  LASGlobalData<std::vector<float> > *ptr = &theData;
  Int_t retval = tree->SetBranchAddress("lasRawData", &ptr);
  if( retval < 0 ) throw LAS::Exception("Failed to set branch address for lasRawData");
  //UInt_t eventnumber;
  Int_t eventnumber;
  retval = tree->SetBranchAddress("eventnumber", &eventnumber);
  if( retval < 0 ) throw LAS::Exception("Failed to set branch address for eventnumber");
  Int_t lumiblock;
  retval = tree->SetBranchAddress("lumiblock", &lumiblock);
  if( retval < 0 ) throw LAS::Exception("Failed to set branch address for lumiblock");
  UInt_t unixTime;
  retval = tree->SetBranchAddress("unixTime", &unixTime);
  if( retval < 0 ) throw LAS::Exception("Failed to set branch address for unixTime");

  // Define a strip range in which to look for the maximum
  // In this case the 60 central strips
  LASGlobalData<int> range_low, range_high;; 
  fill_strip_ranges(range_low, range_high);

  // Object for results
  LASGlobalData<Avec> max_signal;
  Avec v_unixTime;
  Avec v_eventnumber;
  Avec v_lumiblock;

  //Int_t nentries = (Int_t)tree->GetEntries();
  if(nentries <= 0 || nentries > tree->GetEntries()) nentries = tree->GetEntries();
  std::cout << "Process " << nentries << " out of " << tree->GetEntries() << " Entries" << std::endl;

  Avec signal_all(nentries);

  // Loop over the tree
  for (Long64_t entry_nr = 0; entry_nr < nentries; entry_nr++) {
    // Get an Entry from the tree (one event)
    tree->GetEntry(entry_nr);
    if(! (entry_nr % 5000))std::cout << "Processing entry nr. " << entry_nr << "  Event Nr. " << eventnumber << std::endl; 

    v_unixTime.push_back(unixTime);
    v_eventnumber.push_back(eventnumber);
    v_lumiblock.push_back(lumiblock);

    LASGlobalDataLoop loop;

    do{
      // Get the data for one specific module
      std::vector<float>& buff = loop.GetEntry<std::vector<float> >(theData);
      
      float max =  vector_max(buff, loop.GetEntry<int>(range_low), loop.GetEntry<int>(range_high));
      loop.GetEntry<Avec>(max_signal).push_back(max);
      signal_all.at(entry_nr) += max;
    } while(loop.next());
  }

  // Quantities for AT / TEC separation
  Avec signal_tec_r4_nooverlap(nentries);
  Avec signal_at_nooverlap(nentries);

  LASGlobalDataLoop loop1(LASGlobalDataLoop::TEC);
  do{
    if( loop1.get_ring() == 1 || ((loop1.get_beam() == 0 || loop1.get_beam() == 3 || loop1.get_beam() == 5) && loop1.get_ring() == 0) ) continue;
    signal_tec_r4_nooverlap += loop1.GetEntry<Avec>(max_signal);
  }while(loop1.next());

  LASGlobalDataLoop loop2(LASGlobalDataLoop::AT);
  do{
    if( loop2.get_beam() != 0 && loop2.get_beam() != 3 && loop2.get_beam() != 5 )
      signal_at_nooverlap += loop2.GetEntry<Avec>(max_signal);
  }while(loop2.next());


  // Quantities for intensity labelling
  Avec signal_tecplus(nentries);
  Avec signal_tecminus(nentries);

  Avec signal_tec_r6(nentries);
  Avec signal_tec_r4(nentries);

  Avec signal_tec_close(nentries);
  Avec signal_tec_far(nentries);

  Avec signal_at(nentries);
  Avec signal_tec_disc1(nentries);

  LASGlobalDataLoop loop3(LASGlobalDataLoop::TEC_PLUS);
  do{
    signal_tecplus += loop3.GetEntry<Avec>(max_signal);
  }while(loop3.next());

  LASGlobalDataLoop loop4(LASGlobalDataLoop::TEC_MINUS);
  do{
    signal_tecminus += loop4.GetEntry<Avec>(max_signal);
  }while(loop4.next());

  LASGlobalDataLoop loop5(LASGlobalDataLoop::TEC);
  do{
    if(loop5.get_ring() == 1) signal_tec_r6 += loop5.GetEntry<Avec>(max_signal);
    if(loop5.get_ring() == 0) signal_tec_r4 += loop5.GetEntry<Avec>(max_signal);
    if(loop5.get_zpos() == 0)signal_tec_disc1 += loop5.GetEntry<Avec>(max_signal);
    if(loop5.get_zpos() <  2 || loop5.get_zpos() >  7)signal_tec_far += loop5.GetEntry<Avec>(max_signal);
    if(loop5.get_zpos() >= 2 && loop5.get_zpos() <= 7)signal_tec_close += loop5.GetEntry<Avec>(max_signal);
  }while(loop5.next());

  LASGlobalDataLoop loop6(LASGlobalDataLoop::AT);
  do{
    signal_at += loop6.GetEntry<Avec>(max_signal);
  }while(loop6.next());


  v_unixTime -= make_root_time_offset(); // Subtract the standard root time offset

//   Avec::OUTPUT_PRECISION = 10;
//   vecwrite((v_eventnumber | v_unixTime | signal_tecplus | signal_tecminus | signal_at), output_filename);

  TFile outf(output_rootfilename.c_str(),"RECREATE");
  if(full_output) max_signal.Write("max_signal");
  v_eventnumber.Write("eventnumber");
  v_unixTime.Write("unixTime");
  v_lumiblock.Write("lumiblock");

  signal_all.Write("signal_all");

  signal_tec_r4_nooverlap.Write("signal_tec_r4_nooverlap");
  signal_at_nooverlap.Write("signal_at_nooverlap");

  signal_tec_close.Write("signal_tec_close");
  signal_tec_far.Write("signal_tec_far");

  signal_tec_r6.Write("signal_tec_r6");
  signal_tec_r4.Write("signal_tec_r4");

  signal_tecplus.Write("signal_tecplus");
  signal_tecminus.Write("signal_tecminus");
  signal_tec_disc1.Write("signal_tec_disc1");
  signal_at.Write("signal_at");

  outf.Close();

  if(plot)  control_max_signal(output_rootfilename, false);

  std::time_t stop_time = std::time(NULL);
  std::cout << "Wall time passed: " << std::difftime(stop_time, start_time) << " s.\n";
}

