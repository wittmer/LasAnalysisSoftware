//#include "LASGlobalData.h"
//#include "LASGlobalDataLoop.h"

//#include <boost/program_options.hpp>

#define __AVECROOT__
#include "Avec.h"

#include "LAS_basic_tools.h"
#include "LAS_globaldata_tools.h"
#include "LAS_Tec_Reconstruction.h"
#include "LAS_control_plots.h"
#include "LAS_RDC_tools.h"
#include "LAS_At_rec.h"

#include "TROOT.h"
#include "TCanvas.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Main analysis function
void main_analysis(const std::string& filename, const std::string& output_filename, Long64_t nentries = -1, bool max_signal = true, bool noise_filter = true, bool block_slicer = true, bool labels = true, bool profiles = true, bool positions = true, bool tec_reco = true)
{
  Avec::VERBOSE_FLAG = 0;

  LAS::AnalysisParameters par;

  // New Noise Filter Algorithm
  par.noise_algo_type = 2;
  par.noise_signal_thresh = 600.0;

  //par.noise_algo_type = 0;
  //par.noise_signal_thresh = 6000.0;


  try{
    if(max_signal)   get_max_signal(filename, output_filename, false, false, nentries);
    if(noise_filter) Noise_Filter  (output_filename, par);
    if(block_slicer) Block_Slicer  (output_filename, false);
    if(labels)       Event_Labels  (output_filename);
    if(labels)       Block_Label_Check(output_filename);
    replace_step_mask("step_mask_1.txt", output_filename);
    //replace_step_mask("step_mask_improved_disc6.txt", output_filename);
    if(profiles)     Signal_Integration(filename, output_filename);
    if(positions)    Position_Calculation(output_filename, par);

    if(tec_reco){
      std::cout << "Computing TEC Alignment Parameters (within the run)" << std::endl;
      LASGlobalData<Avec> ref_profiles = global_data_get<Avec>("profiles", output_filename);
      LASGlobalData<double> ref_pos;
      LASGlobalData<int> ref_mask;
      calc_pos(ref_profiles, ref_pos, ref_mask, par.edge_threshold);
      TEC_Parameter_Reconstruction(output_filename, ref_pos, ref_mask, "tec_rec_run");
    }
  }
  catch(LAS::Exception& e){
    std::cerr << "LAS::Exception when processing file " << filename << ":\nMessage: "<< e.what() << std::endl;
  }
}


//! Draw a set of Canvasses with TEC alignment parameter history
// The canvasses are saved in PNG format and into a root file
void history_plots(const std::string& listfile, const std::string& outfile, Bool_t suppress_graph = kTRUE)
{
  gROOT->SetStyle("Plain");
  // If desired, switch temporarily to batch mode, to suppress graphical output
  Bool_t batch_mode = gROOT->IsBatch();
  gROOT->SetBatch(suppress_graph);

  std::ifstream f(listfile.c_str());
  if(!f){
    std::cerr << "Could not open file " << listfile << std::endl;
    return;
  }

  std::string ref_file = get_next_line(f);
  std::cout << "Reference file: " << ref_file << std::endl;

  std::vector<std::string> file_list;
  while(f){
    std::string buffer = get_next_line(f);
    if(buffer != "")file_list.push_back(buffer);
  }

  std::cout << "Found " << file_list.size() << " lines in file " << listfile << std::endl;
  if(file_list.empty()){
    std::cerr << "These are too few for a history plot!!" << std::endl;
    return;
  }

  alpar_history_fine(file_list);

  //Write results to file
  TFile of(outfile.c_str(), "RECREATE");
  of.mkdir("Details_AT");
  of.mkdir("Details_TEC_internal");
  of.cd("Details_TEC_internal");
  TCanvas* cv = (TCanvas*)(gROOT->FindObject("beampar_tecp_full"));
  if(cv) cv->Write("beampar_tecp");
  cv = (TCanvas*)(gROOT->FindObject("beampar_tecm_full"));
  if(cv) cv->Write("beampar_tecm");
  cv = (TCanvas*)(gROOT->FindObject("cv_tecp"));
  if(cv) cv->Write("discpar_tecp");
  cv = (TCanvas*)(gROOT->FindObject("cv_tecm"));
  if(cv) cv->Write("discpar_tecm");
  cv = (TCanvas*)(gROOT->FindObject("cv_globpar_full"));
  if(cv) cv->Write("globpar");
  of.Close();

  gROOT->SetBatch(batch_mode);
}

// Perform the AT reconstruction for all files listed in 'listfile'
// The refernce position is the average of all positions from the first entry found in 'listfile'
void rec_AT(const std::string& listfile)
{
  // List and vector for the results
  std::vector<AtPar> par_list;
  Avec time_axis;

  // Get the list of files and the name of the reference file
  std::string ref_filename;
  std::vector<std::string> file_list;
  get_file_list(listfile, ref_filename, file_list);


  //Get the reference position
  std::cout << "Reference file is: " << ref_filename << std::endl;
  LASGlobalData<double> ref_pos;
  LASGlobalData<double> ref_err;
  LASGlobalData<int> ref_mask;
  get_ref_pos( ref_filename, ref_pos, ref_err, ref_mask);
  
  // Loop over all files in the list
  for(unsigned int i = 0; i < file_list.size(); i++){
    std::cout << "Processing file " << file_list[i] << std::endl;

    //Avec block_num = avec_get("block_nr",file_list.back());
    Avec good_block = avec_get("good_block",file_list[i]);
    Avec block_timestamp = avec_get("block_timestamp", file_list[i]);

    std::cout << vsum(good_block) << " good blocks out of " << good_block.size() << std::endl;

    for(unsigned int blnr = 0; blnr < good_block.size(); blnr++){
      if(good_block[blnr] == 0)continue;
      std::ostringstream blkmean;
      std::ostringstream blkrms;
      std::ostringstream blknorm;
      blkmean<<"positions_"<<blnr;
      blkrms<<"rms_av_"<<blnr;
      blknorm<<"norm_"<<blnr;
      LASGlobalData<double> pos = global_data_get<double>(blkmean.str(),file_list[i]);
      LASGlobalData<double> rms = global_data_get<double>(blkrms.str(),file_list[i]);
      LASGlobalData<int> norm = global_data_get<int>(blknorm.str(),file_list[i]);
      LASGlobalData<double> dif = (pos - ref_pos)*pitch()/r0();
      correct_signs(dif);
      LASGlobalData<double> error = rms/sqrt(norm)*pitch()/r0();
      AtPar results;
      AT_TIB2TOB_6par(dif,error,results);
      AT_chi2(dif,error,results);
      par_list.push_back(results);
      time_axis.push_back( block_timestamp[blnr] );
    }
  }

  DrawAtPar(par_list, time_axis);
}
 

// Do the AT reconstruction for one specific block
// Reference is the average position of all blocks in ref_filename
// Results are printed and residuals are drawn
void rec_AT_block(const std::string& filename, unsigned int block_nr, const std::string& ref_filename)
{

  std::cout << "Reference file is: " << ref_filename << std::endl;
  LASGlobalData<double> ref_pos;
  LASGlobalData<double> ref_err;
  LASGlobalData<int> ref_mask;
  get_ref_pos( ref_filename, ref_pos, ref_err, ref_mask);
  
  std::cout << "Processing file " << filename << std::endl;
  
  TFile f(filename.c_str(), "READ");
  if(!f.IsOpen()){
    std::cerr << "Could not open file " << filename << std::endl;
    return;
  }

  Avec good_block = avec_get("good_block",f);

  if(block_nr >= good_block.size()){
    std::cerr << "Error in rec_AT_TIB: block number " << block_nr << " is larger than good_block.size(): " << good_block.size() << std::endl;
    return;
  }
    
  if(good_block[block_nr] == 0){
    std::cerr << "Error in rec_AT_TIB: block nr. " << block_nr << " is 'not good'" << std::endl;
    return;
  }

  std::ostringstream blkmean;
  std::ostringstream blkrms;
  std::ostringstream blknorm;
  std::ostringstream blkmask;

  blkmean << "positions_"<< block_nr;
  blkrms  << "rms_av_"   << block_nr;
  blknorm << "norm_"     << block_nr;
  blkmask << "positions_mask_"<< block_nr;

  LASGlobalData<double> pos = global_data_get<double>(blkmean.str(), f);
  LASGlobalData<double> rms = global_data_get<double>(blkrms.str(), f);
  LASGlobalData<int> norm = global_data_get<int>(blknorm.str(), f);
  LASGlobalData<int> mask = global_data_get<int>(blkmask.str(), f);

  LASGlobalData<double> dif = (pos - ref_pos)*pitch()/r0();
  correct_signs(dif);
  LASGlobalData<double> error = rms/sqrt(norm)*pitch()/r0();

  AtPar results;
  AT_TIB2TOB_6par(dif,error,results);
  AT_chi2(dif,error,results);
  std::cout << results << std::endl;

  LASGlobalData<double> resid = AT_spot_rec(results) - dif;

  draw_global_data_z(resid, mask, LAS::TIB, "TIB");
  draw_global_data_z(resid, mask, LAS::TOB, "TOB");
  draw_global_data_z(resid, mask, LAS::AT, "AT");

  LASGlobalData<double> xpro = resid * cos(theta());
  LASGlobalData<double> ypro = resid * sin(theta());
  draw_global_data_z(xpro, mask, LAS::AT, "AT_X");
  draw_global_data_z(ypro, mask, LAS::AT, "AT_Y");
}



//! Make a MC experiment for tEC internal alignment
void tec_MC()
{
  TECReconstructor TR;
  TR.Init_common();
  TR.Init_full();

  // Generate random alignment parameters
  // Rotations are in rad, displacements are in mm
  double sigma_dphi = 1e-3; 
  double sigma_dx = 1;
  double sigma_dy = 1;
  LasAlPar alpar;
  alpar_rand_gauss(alpar, sigma_dphi, sigma_dx, sigma_dy);
  separate_correlations_tec(alpar.tecp);
  separate_correlations_tec(alpar.tecm);
  //alpar.print();

  // Calculate Laser spot positions (in mm)
  LASGlobalData<double> diff = TR.dx_rec_tec(alpar);

  // Smear laser spot positions
  LASGlobalData<double> sigma_pos(0.005); // This can be specific for every module
  diff += global_data_random_gauss() * sigma_pos;

  // Possibility to mask out some modules
  LASGlobalData<int> mask(1);
  diff *= mask;

  // Reconstruct parameters
  LasAlPar alpar2;
  TR.reconstruct_tec_full(diff, mask, alpar2);
  //alpar2.print();

  std::cout << "Difference between reconstructed and simulated parameters:" << std::endl;
  (alpar2 - alpar).print();


  LASGlobalData<double> residuals;
  TR.compute_residuals(diff, alpar2, residuals);
  plot_global_data(residuals, mask, LAS::TEC);
}



void rec_LAS_TEC(const std::string& run_list_name, const std::string& path_prefix = "")
{
 
  std::string ref_filename;
  std::vector<std::string> file_list;
  get_file_list(run_list_name, ref_filename, file_list);

  if (path_prefix != "") ref_filename = path_prefix + "/" + ref_filename;

  // Get the reference positions
  std::cout << "Reference file is: " << ref_filename << std::endl;
  LASGlobalData<double> ref_pos;
  LASGlobalData<double> ref_err;
  LASGlobalData<int> ref_mask;
  get_ref_pos(ref_filename, ref_pos, ref_err, ref_mask);
  
  // Loop over all runs  
  for(unsigned int i = 0; i < file_list.size(); i++){
    std::string filename =   (path_prefix != "" ? path_prefix + "/" : "" ) + file_list[i];
    std::cout << "\nProcessing file " << filename << std::endl;
    TEC_Parameter_Reconstruction(file_list[i], ref_pos, ref_mask);
  }

}
