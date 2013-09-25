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
#include "TAxis.h"
#include "TMatrixT.h"

#include <ctime>

//! Draw a set of Canvasses with TEC alignment parameter history
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Main analysis function
void main_analysis(const std::string& filename, const std::string& output_filename, Long64_t nentries, bool max_signal, bool noise_filter, bool block_slicer, bool labels, bool profiles, bool positions, bool tec_reco)
{
  Avec::VERBOSE_FLAG = 0;

  LAS::AnalysisParameters par;

  // New Noise Filter Algorithm
  //par.noise_algo_type = 2;
  //par.noise_signal_thresh = 600.0;

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


// The canvasses are saved in PNG format and into a root file
void history_plots(const std::string& listfile, const std::string& outfile, Bool_t suppress_graph)
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

// Check if parts of the results exceed thresholds
bool bad_AT_results(const AtPar& results)
{
  double tib_par_threshold = 1000.;

  if(fabs(results.tib_dx) > tib_par_threshold) return true; //TIB x-displacement
  if(fabs(results.tib_dy) > tib_par_threshold) return true; //TIB y-displacement
  if(fabs(results.tib_rx) > tib_par_threshold) return true; //TIB x-rotation
  if(fabs(results.tib_ry) > tib_par_threshold) return true; //TIB y-rotation
  if(fabs(results.tib_rz) > tib_par_threshold) return true; //TIB z-rotation
  if(fabs(results.tib_tz) > tib_par_threshold) return true; //TIB z-torsion

  return false;
}

void compress_AT_plots(const std::string& parameters_file, Avec::size_type n, const std::string& directory, LAS::beam_group beam_group)
{
  TFile f(parameters_file.c_str(), "READ");

  Avec time = avec_get(directory + "time", f);
  if(time.empty()) throw LAS::Exception("Error in compress_AT_plots: vector with timestamps is empty");


  Avec beam_slp_0 = avec_get(directory + "beam_slp_0", f);
  Avec beam_slp_1 = avec_get(directory + "beam_slp_1", f);
  Avec beam_slp_2 = avec_get(directory + "beam_slp_2", f);
  Avec beam_slp_3 = avec_get(directory + "beam_slp_3", f);
  Avec beam_slp_4 = avec_get(directory + "beam_slp_4", f);
  Avec beam_slp_5 = avec_get(directory + "beam_slp_5", f);
  Avec beam_slp_6 = avec_get(directory + "beam_slp_6", f);
  Avec beam_slp_7 = avec_get(directory + "beam_slp_7", f);

  Avec beam_int_0 = avec_get(directory + "beam_int_0", f);
  Avec beam_int_1 = avec_get(directory + "beam_int_1", f);
  Avec beam_int_2 = avec_get(directory + "beam_int_2", f);
  Avec beam_int_3 = avec_get(directory + "beam_int_3", f);
  Avec beam_int_4 = avec_get(directory + "beam_int_4", f);
  Avec beam_int_5 = avec_get(directory + "beam_int_5", f);
  Avec beam_int_6 = avec_get(directory + "beam_int_6", f);
  Avec beam_int_7 = avec_get(directory + "beam_int_7", f);

  Avec er_beam_slp_0(beam_slp_0.size());
  Avec er_beam_slp_1(beam_slp_1.size());
  Avec er_beam_slp_2(beam_slp_2.size());
  Avec er_beam_slp_3(beam_slp_3.size());
  Avec er_beam_slp_4(beam_slp_4.size());
  Avec er_beam_slp_5(beam_slp_5.size());
  Avec er_beam_slp_6(beam_slp_6.size());
  Avec er_beam_slp_7(beam_slp_7.size());

  Avec er_beam_int_0(beam_int_0.size());
  Avec er_beam_int_1(beam_int_1.size());
  Avec er_beam_int_2(beam_int_2.size());
  Avec er_beam_int_3(beam_int_3.size());
  Avec er_beam_int_4(beam_int_4.size());
  Avec er_beam_int_5(beam_int_5.size());
  Avec er_beam_int_6(beam_int_6.size());
  Avec er_beam_int_7(beam_int_7.size());


  Avec time_c = time;
  avec_compress(time_c, 
		(beam_slp_0   | beam_slp_1   | beam_slp_2   | beam_slp_3   | beam_slp_4   | beam_slp_5   | beam_slp_6   | beam_slp_7),
		(er_beam_slp_0   | er_beam_slp_1   | er_beam_slp_2   | er_beam_slp_3   | er_beam_slp_4   | er_beam_slp_5   | er_beam_slp_6   | er_beam_slp_7), 
		n, true);

  (new TCanvas("BEAM_A_compressed","BEAM_A_compressed"))->SetGridy(); 
  TMultiGraph* gra = avec_draw(time_c, (beam_slp_0 | beam_slp_1 | beam_slp_2 | beam_slp_3 | beam_slp_4 | beam_slp_5 | beam_slp_6 | beam_slp_7), (er_beam_slp_0 | er_beam_slp_1 | er_beam_slp_2 | er_beam_slp_3 | er_beam_slp_4 | er_beam_slp_5 | er_beam_slp_6 | er_beam_slp_7),"Beam Slopes","Date","Beam slopes[#murad] ","AP");
  gra->GetXaxis()->SetNdivisions(505, kTRUE);
  gra->GetXaxis()->SetTimeDisplay(1);
  gra->GetXaxis()->SetTimeFormat("%d/%m");

  time_c = time;
  avec_compress(time_c, 
		(beam_int_0   | beam_int_1   | beam_int_2   | beam_int_3   | beam_int_4   | beam_int_5   | beam_int_6   | beam_int_7),
		(er_beam_int_0   | er_beam_int_1   | er_beam_int_2   | er_beam_int_3   | er_beam_int_4   | er_beam_int_5   | er_beam_int_6   | er_beam_int_7), 
		n, true);

  (new TCanvas("BEAM_B_compressed","BEAM_B_compressed"))->SetGridy(); 
  TMultiGraph* grb = avec_draw(time_c, (beam_int_0 | beam_int_1 | beam_int_2 | beam_int_3 | beam_int_4 | beam_int_5 | beam_int_6 | beam_int_7), (er_beam_int_0 | er_beam_int_1 | er_beam_int_2 | er_beam_int_3 | er_beam_int_4 | er_beam_int_5 | er_beam_int_6 | er_beam_int_7),"Beam Offsets","Date","Beam intersepts[#mum] ","AP");
  grb->GetXaxis()->SetNdivisions(505, kTRUE);
  grb->GetXaxis()->SetTimeDisplay(1);
  grb->GetXaxis()->SetTimeFormat("%d/%m");
 

  Avec AT_chi = avec_get(directory + "AT_chi", f);
  Avec AT_chi0 = avec_get(directory + "AT_chi0", f);
  time_c = time;
  avec_compress(time_c, (AT_chi | AT_chi0), n, true);

  Avec tob_chi = avec_get(directory + "tob_chi", f);
  Avec tob_chi0 = avec_get(directory + "tob_chi0", f);
  time_c = time;
  avec_compress(time_c, (tob_chi | tob_chi0), n, true);

  Avec tib_chi = avec_get(directory + "tib_chi", f);
  Avec tib_chi0 = avec_get(directory + "tib_chi0", f);
  time_c = time;
  avec_compress(time_c, (tib_chi | tib_chi0), n, true);

  TCanvas* cvca = new TCanvas("CH_AT_compressed","CH_AT_compressed");
  cvca->SetLogy();
  cvca->SetGridy();
  TMultiGraph* grca = avec_draw(time_c, (AT_chi | AT_chi0),"#chi2/ndf for AT modules","Date","#chi2/ndf","AP");
  grca->GetXaxis()->SetNdivisions(505, kTRUE);
  grca->GetXaxis()->SetTimeDisplay(1);
  grca->GetXaxis()->SetTimeFormat("%d/%m");
  grca->GetXaxis()->SetRangeUser(547.1e6, 572.5e6);

  grca->GetYaxis()->SetRangeUser(0.1, 1e4);
 
  std::vector<std::string> Legendca;
  Legendca.push_back("AT after fit");
  Legendca.push_back("AT before fit");
  AddLegend(grca, Legendca, 1, 0.65, 0.7);
  
  TCanvas* cvcb = new TCanvas("CH_TIB_compressed","CH_TIB_compressed");
  cvcb->SetLogy();
  cvcb->SetGridy();
  TMultiGraph* grcb = avec_draw(time_c, (tib_chi | tib_chi0),"#chi2/n for TIB modules","Date","#chi2/n","AP");
  grcb->GetXaxis()->SetNdivisions(505, kTRUE);
  grcb->GetXaxis()->SetTimeDisplay(1);
  grcb->GetXaxis()->SetTimeFormat("%d/%m");
  grcb->GetXaxis()->SetRangeUser(547.1e6, 572.5e6);

  grcb->GetYaxis()->SetRangeUser(0.1, 1e4);
 
  std::vector<std::string> Legendcb;
  Legendcb.push_back("TIB after fit");
  Legendcb.push_back("TIB before fit");
  AddLegend(grcb, Legendcb, 1, 0.65, 0.7);
  
  TCanvas* cvcc = new TCanvas("CH_TOB_compressed","CH_TOB_compressed");
  cvcc->SetLogy();
  cvcc->SetGridy();
  TMultiGraph* grcc = avec_draw(time_c, (tob_chi | tob_chi0),"#chi2/n for TOB modules","Date","#chi2/n","AP");
  grcc->GetXaxis()->SetNdivisions(505, kTRUE);
  grcc->GetXaxis()->SetTimeDisplay(1);
  grcc->GetXaxis()->SetTimeFormat("%d/%m");
  grcc->GetXaxis()->SetRangeUser(547.1e6, 572.5e6);

  grcc->GetYaxis()->SetRangeUser(0.1, 1e4);
 
  std::vector<std::string> Legendcc;
  Legendcc.push_back("TOB after fit");
  Legendcc.push_back("TOB before fit");
  AddLegend(grcc, Legendcc, 1, 0.65, 0.7);
  
  /////////
  // TIB //
  /////////
  if(beam_group == LAS::TIB || beam_group == LAS::AT || beam_group == LAS::ALL){
    Avec tib_dx = avec_get(directory + "tib_dx", f);
    Avec tib_dy = avec_get(directory + "tib_dy", f);
    Avec tib_rx = avec_get(directory + "tib_rx", f);
    Avec tib_ry = avec_get(directory + "tib_ry", f);
    Avec tib_rz = avec_get(directory + "tib_rz", f);
    Avec tib_tz = avec_get(directory + "tib_tz", f);
    Avec tib_er_dx = avec_get(directory + "tib_er_dx", f);
    Avec tib_er_dy = avec_get(directory + "tib_er_dy", f);
    Avec tib_er_rx = avec_get(directory + "tib_er_rx", f);
    Avec tib_er_ry = avec_get(directory + "tib_er_ry", f);
    Avec tib_er_rz = avec_get(directory + "tib_er_rz", f);
    Avec tib_er_tz = avec_get(directory + "tib_er_tz", f);

    (new TCanvas("TIB_AT_compressed", "TIB_AT_compressed"))->SetGridy();

    time_c = time;
    avec_compress(time_c, (tib_dx | tib_dy | tib_rx | tib_ry | tib_rz | tib_tz), (tib_er_dx | tib_er_dy | tib_er_rx | tib_er_ry | tib_er_rz | tib_er_tz), n, true);
    TMultiGraph* gr = avec_draw(time_c, (tib_dx | tib_dy | tib_rx | tib_ry | tib_rz | tib_tz), (tib_er_dx | tib_er_dy | tib_er_rx | tib_er_ry | tib_er_rz | tib_er_tz), "TIB vs TOB", "time", "#Deltax, #Deltay [#mum], Rx,Ry,Rz [#murad], Tz[#murad/m] ", "AP");
    gr->GetXaxis()->SetNdivisions(505, kTRUE);
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat("%d/%m");

    std::vector<std::string> Legend;
    Legend.push_back("#Delta x");
    Legend.push_back("#Delta y");
    Legend.push_back("Rx");
    Legend.push_back("Ry");
    Legend.push_back("Rz");
    Legend.push_back("Tz");
    AddLegend(gr, Legend, 1, 0.9, 0.6, 1.0, 1.0);
  }

  /////////////
  // TEC+ AT //
  /////////////
  if(beam_group == LAS::TEC_PLUS_AT || beam_group == LAS::TEC_AT || beam_group == LAS::AT || beam_group == LAS::ALL){
    Avec tecp_dx = avec_get(directory + "tecp_dx", f);
    Avec tecp_dy = avec_get(directory + "tecp_dy", f);
    Avec tecp_rz = avec_get(directory + "tecp_rz", f);
    Avec tecp_er_dx = avec_get(directory + "tecp_er_dx", f);
    Avec tecp_er_dy = avec_get(directory + "tecp_er_dy", f);
    Avec tecp_er_rz = avec_get(directory + "tecp_er_rz", f);

    (new TCanvas("TECP_AT_compressed", "TECP_AT_compressed"))->SetGridy();

    time_c = time;
    avec_compress(time_c, (tecp_dx | tecp_dy | tecp_rz), (tecp_er_dx | tecp_er_dy | tecp_er_rz), n, true);
    TMultiGraph* gr = avec_draw(time_c, (tecp_dx | tecp_dy | tecp_rz), (tecp_er_dx | tecp_er_dy | tecp_er_rz), "TECP vs TOB", "time", "#Deltax, #Deltay [#mum], Rz [#murad]", "AP");
    gr->GetXaxis()->SetNdivisions(505, kTRUE);
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat("%d/%m");

    std::vector<std::string> Legendp;
    Legendp.push_back("#Delta x");
    Legendp.push_back("#Delta y");
    Legendp.push_back("Rz");
    AddLegend(gr, Legendp, 1, 0.9, 0.6, 1.0, 1.0);

  }

  /////////////
  // TEC- AT //
  /////////////
  if(beam_group == LAS::TEC_MINUS_AT || beam_group == LAS::TEC_AT || beam_group == LAS::AT || beam_group == LAS::ALL){
    Avec tecm_dx = avec_get(directory + "tecm_dx", f);
    Avec tecm_dy = avec_get(directory + "tecm_dy", f);
    Avec tecm_rz = avec_get(directory + "tecm_rz", f);
    Avec tecm_er_dx = avec_get(directory + "tecm_er_dx", f);
    Avec tecm_er_dy = avec_get(directory + "tecm_er_dy", f);
    Avec tecm_er_rz = avec_get(directory + "tecm_er_rz", f);

    (new TCanvas("TECM_AT_compressed", "TECM_AT_compressed"))->SetGridy();

    time_c = time;
    avec_compress(time_c, (tecm_dx | tecm_dy | tecm_rz), (tecm_er_dx | tecm_er_dy | tecm_er_rz), n, true);
    TMultiGraph* gr = avec_draw(time_c, (tecm_dx | tecm_dy | tecm_rz), (tecm_er_dx | tecm_er_dy | tecm_er_rz), "TECM vs TOB", "time", "#Deltax, #Deltay [#mum], Rz [#murad]", "AP");

    gr->GetXaxis()->SetNdivisions(505, kTRUE);
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat("%d/%m");

    std::vector<std::string> Legendp;
    Legendp.push_back("#Delta x");
    Legendp.push_back("#Delta y");
    Legendp.push_back("Rz");
    AddLegend(gr, Legendp, 1, 0.9, 0.6, 1.0, 1.0);
  }

  f.Close();
}


// Perform the AT reconstruction for all files listed in 'listfile'
// The refernce position is the average of all positions from the first entry found in 'listfile'
// void rec_AT(const std::string& listfile)
// {
//   // List and vector for the results
//   std::vector<AtPar> par_list;
//   Avec time_axis;

//   // Get the list of files and the name of the reference file
//   std::string ref_filename;
//   std::vector<std::string> file_list;
//   get_file_list(listfile, ref_filename, file_list);


//   //Get the reference position
//   std::cout << "Reference file is: " << ref_filename << std::endl;
//   LASGlobalData<double> ref_pos;
//   LASGlobalData<double> ref_err;
//   LASGlobalData<int> ref_mask;
//   get_ref_pos( ref_filename, ref_pos, ref_err, ref_mask);
  
//   // Loop over all files in the list
//   for(unsigned int i = 0; i < file_list.size(); i++){
//     std::cout << "Processing file " << file_list[i] << std::endl;

//     //Avec block_num = avec_get("block_nr",file_list.back());
//     Avec good_block = avec_get("good_block",file_list[i]);
//     Avec block_timestamp = avec_get("block_timestamp", file_list[i]);

//     std::cout << vsum(good_block) << " good blocks out of " << good_block.size() << std::endl;

//     for(unsigned int blnr = 0; blnr < good_block.size(); blnr++){
//       if(good_block[blnr] == 0)continue;
//       std::ostringstream blkmean;
//       std::ostringstream blkrms;
//       std::ostringstream blknorm;
//       blkmean<<"positions_"<<blnr;
//       blkrms<<"rms_av_"<<blnr;
//       blknorm<<"norm_"<<blnr;
//       LASGlobalData<double> pos = global_data_get<double>(blkmean.str(),file_list[i]);
//       LASGlobalData<double> rms = global_data_get<double>(blkrms.str(),file_list[i]);
//       LASGlobalData<int> norm = global_data_get<int>(blknorm.str(),file_list[i]);
//       LASGlobalData<double> dif = (pos - ref_pos)*pitch()/r0();
//       correct_signs(dif);
//       LASGlobalData<double> error = rms/sqrt(norm)*pitch()/r0();
//       AtPar results;
//       AT_TIB2TOB_6par_old(dif,error,results);
//       AT_chi2_old(dif,error,results);
//       par_list.push_back(results);
//       time_axis.push_back( block_timestamp[blnr] );
//     }
//   }

//   DrawAtPar_old(par_list, time_axis);
// }
 

// // Do the AT reconstruction for one specific block
// // Reference is the average position of all blocks in ref_filename
// // Results are printed and residuals are drawn
// void rec_AT_block(const std::string& filename, unsigned int block_nr, const std::string& ref_filename)
// {

//   std::cout << "Reference file is: " << ref_filename << std::endl;
//   LASGlobalData<double> ref_pos;
//   LASGlobalData<double> ref_err;
//   LASGlobalData<int> ref_mask;
//   get_ref_pos( ref_filename, ref_pos, ref_err, ref_mask);
  
//   std::cout << "Processing file " << filename << std::endl;
  
//   TFile f(filename.c_str(), "READ");
//   if(!f.IsOpen()){
//     std::cerr << "Could not open file " << filename << std::endl;
//     return;
//   }

//   Avec good_block = avec_get("good_block",f);

//   if(block_nr >= good_block.size()){
//     std::cerr << "Error in rec_AT_TIB: block number " << block_nr << " is larger than good_block.size(): " << good_block.size() << std::endl;
//     return;
//   }
    
//   if(good_block[block_nr] == 0){
//     std::cerr << "Error in rec_AT_TIB: block nr. " << block_nr << " is 'not good'" << std::endl;
//     return;
//   }

//   std::ostringstream blkmean;
//   std::ostringstream blkrms;
//   std::ostringstream blknorm;
//   std::ostringstream blkmask;

//   blkmean << "positions_"<< block_nr;
//   blkrms  << "rms_av_"   << block_nr;
//   blknorm << "norm_"     << block_nr;
//   blkmask << "positions_mask_"<< block_nr;

//   LASGlobalData<double> pos = global_data_get<double>(blkmean.str(), f);
//   LASGlobalData<double> rms = global_data_get<double>(blkrms.str(), f);
//   LASGlobalData<int> norm = global_data_get<int>(blknorm.str(), f);
//   LASGlobalData<int> mask = global_data_get<int>(blkmask.str(), f);

//   LASGlobalData<double> dif = (pos - ref_pos)*pitch()/r0();
//   correct_signs(dif);
//   LASGlobalData<double> error = rms/sqrt(norm)*pitch()/r0();

//   AtPar results;
//   AT_TIB2TOB_6par_old(dif,error,results);
//   AT_chi2_old(dif,error,results);
//   std::cout << results << std::endl;

//   LASGlobalData<double> resid = AT_spot_rec(results) - dif;

//   draw_global_data_z(resid, mask, LAS::TIB, "TIB");
//   draw_global_data_z(resid, mask, LAS::TOB, "TOB");
//   draw_global_data_z(resid, mask, LAS::AT, "AT");

//   LASGlobalData<double> xpro = resid * cos(theta());
//   LASGlobalData<double> ypro = resid * sin(theta());
//   draw_global_data_z(xpro, mask, LAS::AT, "AT_X");
//   draw_global_data_z(ypro, mask, LAS::AT, "AT_Y");
// }


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
