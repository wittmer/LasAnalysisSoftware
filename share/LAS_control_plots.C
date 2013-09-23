#include "LAS_control_plots.h"

#include <iomanip>

#define __AVECROOT__
#include "Avec.h"
#include "Avec2D.h"

#include "LASGlobalData.h"
#include "LASGlobalDataLoop.h"

#include "LAS_basic_tools.h"
#include "LAS_globaldata_tools.h"

#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TPaveLabel.h"
#include "TLine.h"
#include "TH1.h"
#include "THStack.h"


bool run_selector(const std::string& result_file, double& bad_ratio, Avec::size_type& nr_blocks, double& block_ratio, bool print)
{
  Avec::size_type expected_block_size = 2000;
  double bad_ratio_cut = 0.15;
  bool retval = true;

  Avec good_event = avec_get("good_event", result_file);
  bad_ratio = 1;

  if(good_event.empty()){
    if(print) std::cout << "Error in run_selector, could not find Avec 'good_event'" << std::endl;
    retval =  false;
  }

  Avec::size_type nr_total = good_event.size();
  int nr_good = (int)vsum(good_event);
  int nr_bad  = nr_total - nr_good;
  bad_ratio = (double)nr_bad/nr_good;

  //std::cout << nr_good << " good events, " << nr_bad << " bad events, " << nr_total << " total" << std::endl;
  //std::cout << "bad events / all events = " << bad_ratio << std::endl;

  if (bad_ratio >= bad_ratio_cut) retval = false;

  Avec block_nr = avec_get("block_nr", result_file);
  Avec block_timestamp = avec_get("block_timestamp", result_file);
  Avec block_size = avec_get("block_size", result_file);

  if(block_nr.empty() || block_timestamp.empty() || block_size.empty()){
    if(print) std::cout << "Error in run_selector, could not find all Avecs" << std::endl;
    retval = false;
  }

  nr_blocks = block_nr.size();
  Avec label_step = avec_get("label_step", result_file);
  Avec event_block = avec_get("event_block", result_file);
  if(label_step.empty() || event_block.empty()){
    if(print) std::cout << "Error in run_selector, could not find all Avecs" << std::endl;
    retval = false;
  }

  Avec good_block(nr_blocks, 1);


  for(Avec::size_type i = 0; i < nr_blocks; i++){
    if(block_size[i] != expected_block_size){
      if(print) std::cout << "block " << block_nr[i] << " has " << block_size[i] << " events instead of " << expected_block_size << std::endl;
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
	if(print) std::cerr << "Error in run_selector: Inconsistency between block_size and number of labelled events in event_block" << std::endl;
	good_block[i] = 0;
      }
      else{
	for(Avec::size_type h = 0; h < expected_block_size; h++){
	  if(labels_block[h] != (h%10) + 1){
	    if(print) std::cerr << "Error in run_selector: Incorrect label sequence" << std::endl;
	    good_block[i] = 0;
	  }
	}
      }
    }
  }

  block_ratio = vsum(good_block) / nr_blocks;
  if(print) std::cout << vsum(good_block) << " good blocks out of " << nr_blocks << " (" << block_ratio*100 << "%)" << std::endl;
  if(vsum(good_block) == 0){
    if(print) std::cerr << "No good blocks in this run" << std::endl;
    retval = false;
  }

  if(block_ratio < 0.5){
    if(print) std::cerr << "Too few good blocks in this run" << std::endl;
    retval = false;
  }

  //TFile f(result_file.c_str(), "UPDATE");
  //good_block.Write("good_block");
  //Avec(1, retval ? 1 : 0 ).Write("good_run");
  //f.Close();

  return retval;
}

/////////////////////////////////////////////////////////////////////
// Control plots for the main analysis steps
/////////////////////////////////////////////////////////////////////

double SetDiscParAxes(TMultiGraph* mgr, double min_yrange)
{
  if(!mgr)return 0;
  double yrange = std::max(min_yrange, std::max( fabs(mgr->GetYaxis()->GetXmax()), fabs(mgr->GetYaxis()->GetXmin())));
  mgr->GetYaxis()->SetRangeUser(-yrange, yrange);
  mgr->GetXaxis()->SetNdivisions(505);
  mgr->GetXaxis()->SetTimeDisplay(1);
  mgr->GetXaxis()->SetLabelOffset(0.02);
  mgr->GetXaxis()->SetTimeFormat("#splitline{%d.%m.%y}{%H:%M}");
  return yrange;
}

double calc_rms(const Avec2D& alpar)
{
  Avec sumvl = vsumr(alpar);
  Avec sumsq = vsumr(alpar * alpar);
  Avec rms = sqrt((sumsq - sumvl * sumvl/alpar.size())/(alpar.size()-1));
  //std::cout << "sumsq: " << sumsq << std::endl;
  //std::cout << "sumvl: " << sumvl << std::endl;
  //std::cout << "RMS: " << rms << std::endl;
  return vsum(rms)/rms.size();
}

void draw_rms_label( double value, const std::string& unit, double xlo, double ylo, double xhi, double yhi)
{
  //double ymax = 10;
  std::ostringstream label;
  label << "RMS: " << std::setprecision(3) << value << " " << unit;
  //double ylo=ymax*4/7;
  //double yhi=ymax*6/7;
  TPaveLabel* tpl=new TPaveLabel(xlo,ylo,xhi,yhi,label.str().c_str());
  tpl->SetBorderSize(2);
  tpl->SetTextSize(0.5);
  tpl->Draw();

}


//! Control plots for output of Run_History
void alpar_history_draw(const std::vector<LasAlPar>& par_list, const Avec& xvals, bool bs_frame, bool process_rings )
{
  std::cout << "Running alpar_history_draw with bs_frame=" << (bs_frame ? "true" : "false") << std::endl;
  const char options[] = "AP";

  std::vector<LasAlPar>::size_type nblocks = par_list.size();
  if(nblocks != xvals.size()){
    std::cerr << "Error in alpar_history_draw: par_list and xvals have different size (" << nblocks << " and " << xvals.size() << ")" << std::endl;
    return;
  }
  
  std::cout << "There are " << nblocks << " blocks in this run" << std::endl;

  Avec2D tecp_r4_dphik;
  Avec2D tecp_r4_dxk;
  Avec2D tecp_r4_dyk;
  Avec2D tecp_r6_dphik;
  Avec2D tecp_r6_dxk;
  Avec2D tecp_r6_dyk;

  Avec2D tecm_r4_dphik;
  Avec2D tecm_r4_dxk;
  Avec2D tecm_r4_dyk;
  Avec2D tecm_r6_dphik;
  Avec2D tecm_r6_dxk;
  Avec2D tecm_r6_dyk;

  Avec2D tecp_dphik;
  Avec2D tecp_dxk;
  Avec2D tecp_dyk;

  Avec2D tecm_dphik;
  Avec2D tecm_dxk;
  Avec2D tecm_dyk;

  Avec2D globpar_tecp_r4;
  Avec2D globpar_tecp_r6;
  Avec2D globpar_tecm_r4;
  Avec2D globpar_tecm_r6;

  Avec2D globpar_tecp;
  Avec2D globpar_tecm;

  Avec2D beampar_tecp_r4;
  Avec2D beampar_tecp_r6;
  Avec2D beampar_tecm_r4;
  Avec2D beampar_tecm_r6;

  Avec2D beamparA_full_tecp_r4;
  Avec2D beamparA_full_tecp_r6;
  Avec2D beamparA_full_tecm_r4;
  Avec2D beamparA_full_tecm_r6;

  Avec2D beamparB_full_tecp_r4;
  Avec2D beamparB_full_tecp_r6;
  Avec2D beamparB_full_tecm_r4;
  Avec2D beamparB_full_tecm_r6;

  for(std::vector<LasAlPar>::size_type i = 0; i < nblocks; i++){
    LasAlPar alpar = par_list[i];

    alpar.RescaleRot(1e6);
    alpar.RescaleTrans(1e3);


    // Disc Parameters
    tecp_r4_dphik.push_back(alpar.tecp_r4.Dphik);
    tecp_r4_dxk.push_back(alpar.tecp_r4.Dxk);
    tecp_r4_dyk.push_back(alpar.tecp_r4.Dyk);
    tecp_r6_dphik.push_back(alpar.tecp_r6.Dphik);
    tecp_r6_dxk.push_back(alpar.tecp_r6.Dxk);
    tecp_r6_dyk.push_back(alpar.tecp_r6.Dyk);

    tecm_r4_dphik.push_back(alpar.tecm_r4.Dphik);
    tecm_r4_dxk.push_back(alpar.tecm_r4.Dxk);
    tecm_r4_dyk.push_back(alpar.tecm_r4.Dyk);
    tecm_r6_dphik.push_back(alpar.tecm_r6.Dphik);
    tecm_r6_dxk.push_back(alpar.tecm_r6.Dxk);
    tecm_r6_dyk.push_back(alpar.tecm_r6.Dyk);

    tecp_dphik.push_back(alpar.tecp.Dphik);
    tecp_dxk.push_back(alpar.tecp.Dxk);
    tecp_dyk.push_back(alpar.tecp.Dyk);

    tecm_dphik.push_back(alpar.tecm.Dphik);
    tecm_dxk.push_back(alpar.tecm.Dxk);
    tecm_dyk.push_back(alpar.tecm.Dyk);

    // Global Parameters
    Avec globpar;
    globpar.push_back(alpar.tecp_r4.Dphi0);
    globpar.push_back(alpar.tecp_r4.Dphit);
    globpar.push_back(alpar.tecp_r4.Dx0);
    globpar.push_back(alpar.tecp_r4.Dxt);
    globpar.push_back(alpar.tecp_r4.Dy0);
    globpar.push_back(alpar.tecp_r4.Dyt);
    globpar_tecp_r4.push_back(globpar);
    globpar.clear();
    globpar.push_back(alpar.tecp_r6.Dphi0);
    globpar.push_back(alpar.tecp_r6.Dphit);
    globpar.push_back(alpar.tecp_r6.Dx0);
    globpar.push_back(alpar.tecp_r6.Dxt);
    globpar.push_back(alpar.tecp_r6.Dy0);
    globpar.push_back(alpar.tecp_r6.Dyt);
    globpar_tecp_r6.push_back(globpar);

    globpar.clear();
    globpar.push_back(alpar.tecm_r4.Dphi0);
    globpar.push_back(alpar.tecm_r4.Dphit);
    globpar.push_back(alpar.tecm_r4.Dx0);
    globpar.push_back(alpar.tecm_r4.Dxt);
    globpar.push_back(alpar.tecm_r4.Dy0);
    globpar.push_back(alpar.tecm_r4.Dyt);
    globpar_tecm_r4.push_back(globpar);
    globpar.clear();
    globpar.push_back(alpar.tecm_r6.Dphi0);
    globpar.push_back(alpar.tecm_r6.Dphit);
    globpar.push_back(alpar.tecm_r6.Dx0);
    globpar.push_back(alpar.tecm_r6.Dxt);
    globpar.push_back(alpar.tecm_r6.Dy0);
    globpar.push_back(alpar.tecm_r6.Dyt);
    globpar_tecm_r6.push_back(globpar);

    globpar.clear();
    globpar.push_back(alpar.tecp.Dphi0);
    globpar.push_back(alpar.tecp.Dphit);
    globpar.push_back(alpar.tecp.Dx0);
    globpar.push_back(alpar.tecp.Dxt);
    globpar.push_back(alpar.tecp.Dy0);
    globpar.push_back(alpar.tecp.Dyt);
    globpar_tecp.push_back(globpar);

    globpar.clear();
    globpar.push_back(alpar.tecm.Dphi0);
    globpar.push_back(alpar.tecm.Dphit);
    globpar.push_back(alpar.tecm.Dx0);
    globpar.push_back(alpar.tecm.Dxt);
    globpar.push_back(alpar.tecm.Dy0);
    globpar.push_back(alpar.tecm.Dyt);
    globpar_tecm.push_back(globpar);

    // Beam Parameters
    beampar_tecp_r4.push_back(alpar.tecp_r4.DthetaA & alpar.tecp_r4.DthetaB);
    beampar_tecp_r6.push_back(alpar.tecp_r6.DthetaA & alpar.tecp_r6.DthetaB);
    beampar_tecm_r4.push_back(alpar.tecm_r4.DthetaA & alpar.tecm_r4.DthetaB);
    beampar_tecm_r6.push_back(alpar.tecm_r6.DthetaA & alpar.tecm_r6.DthetaB);

    beamparA_full_tecp_r4.push_back(alpar.tecp.DthetaA_r4);
    beamparA_full_tecp_r6.push_back(alpar.tecp.DthetaA_r6);
    beamparA_full_tecm_r4.push_back(alpar.tecm.DthetaA_r4);
    beamparA_full_tecm_r6.push_back(alpar.tecm.DthetaA_r6);
    
    beamparB_full_tecp_r4.push_back(alpar.tecp.DthetaB_r4);
    beamparB_full_tecp_r6.push_back(alpar.tecp.DthetaB_r6);
    beamparB_full_tecm_r4.push_back(alpar.tecm.DthetaB_r4);
    beamparB_full_tecm_r6.push_back(alpar.tecm.DthetaB_r6);

  }

  globpar_tecp = vtrans(globpar_tecp);
  globpar_tecm = vtrans(globpar_tecm);
  beamparA_full_tecp_r4 = vtrans(beamparA_full_tecp_r4);
  beamparA_full_tecp_r6 = vtrans(beamparA_full_tecp_r6);
  beamparA_full_tecm_r4 = vtrans(beamparA_full_tecm_r4);
  beamparA_full_tecm_r6 = vtrans(beamparA_full_tecm_r6);
  beamparB_full_tecp_r4 = vtrans(beamparB_full_tecp_r4);
  beamparB_full_tecp_r6 = vtrans(beamparB_full_tecp_r6);
  beamparB_full_tecm_r4 = vtrans(beamparB_full_tecm_r4);
  beamparB_full_tecm_r6 = vtrans(beamparB_full_tecm_r6);

  //Avec err_tecp = vsumr(fabs(beamparB_full_tecp_r4 + beamparB_full_tecp_r6))/16;
  //Avec err_tecm = vsumr(fabs(beamparB_full_tecm_r4 + beamparB_full_tecm_r6))/16;
  Avec sum_tecp = vsumr(beamparB_full_tecp_r4 + beamparB_full_tecp_r6);
  Avec sum_tecm = vsumr(beamparB_full_tecm_r4 + beamparB_full_tecm_r6);

  Avec sumsq_tecp = vsumr(beamparB_full_tecp_r4 * beamparB_full_tecp_r4 + beamparB_full_tecp_r6 * beamparB_full_tecp_r6);
  Avec sumsq_tecm = vsumr(beamparB_full_tecm_r4 * beamparB_full_tecm_r4 + beamparB_full_tecm_r6 * beamparB_full_tecm_r6);

  Avec err_tecp = sqrt((sumsq_tecp - sum_tecp*sum_tecp/16)/15);
  Avec err_tecm = sqrt((sumsq_tecm - sum_tecm*sum_tecm/16)/15);

  Avec::size_type size = err_tecp.size();
  Avec2D glob_err_tecp;
  glob_err_tecp.push_back(Avec(size));
  //glob_err_tecp.push_back(err_tecp);
  glob_err_tecp.push_back(Avec(size));
  glob_err_tecp.push_back(Avec(size));
  glob_err_tecp.push_back(Avec(size));
  glob_err_tecp.push_back(Avec(size));
  glob_err_tecp.push_back(Avec(size));
  Avec2D glob_err_tecm;
  glob_err_tecm.push_back(Avec(size));
  //glob_err_tecm.push_back(err_tecm);
  glob_err_tecm.push_back(Avec(size));
  glob_err_tecm.push_back(Avec(size));
  glob_err_tecm.push_back(Avec(size));
  glob_err_tecm.push_back(Avec(size));
  glob_err_tecm.push_back(Avec(size));

  if(bs_frame){
    //std::cout << "Converting global parameters" << std::endl;
    globpar_tecp[0] += globpar_tecp[1] * LAS::z_tec_bs / LAS::L_tec;
    globpar_tecp[2] += globpar_tecp[3] * LAS::z_tec_bs / LAS::L_tec;
    globpar_tecp[4] += globpar_tecp[5] * LAS::z_tec_bs / LAS::L_tec;
    globpar_tecm[0] += globpar_tecm[1] * LAS::z_tec_bs / LAS::L_tec;
    globpar_tecm[2] += globpar_tecm[3] * LAS::z_tec_bs / LAS::L_tec;
    globpar_tecm[4] += globpar_tecm[5] * LAS::z_tec_bs / LAS::L_tec;

    beamparA_full_tecp_r4 += beamparB_full_tecp_r4 * LAS::z_tec_bs / LAS::L_tec;
    beamparA_full_tecp_r6 += beamparB_full_tecp_r6 * LAS::z_tec_bs / LAS::L_tec;
    beamparA_full_tecm_r4 += beamparB_full_tecm_r4 * LAS::z_tec_bs / LAS::L_tec;
    beamparA_full_tecm_r6 += beamparB_full_tecm_r6 * LAS::z_tec_bs / LAS::L_tec;
  }

  //double xlo = xvals.front() * (1 - 0.6) + xvals.back()*0.6;
  //double xhi = xvals.back();
  double xlo = vmin(xvals) * (1 - 0.6) + vmax(xvals)*0.6;
  double xhi = vmax(xvals);
  double ymax = 0;

  TMultiGraph* mgr = 0;

  if(process_rings){
    TCanvas* cv_beampar = new TCanvas("beampar","TEC Beam Parameters", 1500, 1000);
    cv_beampar->Divide(2,2);
    cv_beampar->cd(1);
    mgr = avec_draw(xvals, vtrans(beampar_tecp_r4), "TEC+ Ring 4","","#Delta#theta [#murad]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(beampar_tecp_r4), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);
    
    cv_beampar->cd(2);
    mgr = avec_draw(xvals, vtrans(beampar_tecp_r6), "TEC+ Ring 6","","#Delta#theta [#murad]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(beampar_tecp_r6), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_beampar->cd(3);
    mgr = avec_draw(xvals, vtrans(beampar_tecm_r4), "TEC- Ring 4","","#Delta#theta [#murad]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(beampar_tecm_r4), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_beampar->cd(4);
    mgr = avec_draw(xvals, vtrans(beampar_tecm_r6), "TEC- Ring 6","","#Delta#theta [#murad]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(beampar_tecm_r6), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_beampar->Print("beampar.png");
  }

  TCanvas* cv_beampar_tecp_full = new TCanvas("beampar_tecp_full","TEC+ Beam Parameters", 1500, 1000);
  cv_beampar_tecp_full->Divide(2,2);

  cv_beampar_tecp_full->cd(1);
  mgr = avec_draw(xvals, beamparA_full_tecp_r4, "TEC+ Ring 4 offset","","#Delta#theta_{A} [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(beamparA_full_tecp_r4), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);
		  
  cv_beampar_tecp_full->cd(2);
  mgr = avec_draw(xvals, beamparA_full_tecp_r6, "TEC+ Ring 6 offset","","#Delta#theta_{A} [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(beamparA_full_tecp_r6), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_beampar_tecp_full->cd(3);
  mgr = avec_draw(xvals, beamparB_full_tecp_r4, "TEC+ Ring 4 slope","","#Delta#theta_{B} [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(beamparB_full_tecp_r4), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_beampar_tecp_full->cd(4);
  mgr = avec_draw(xvals, beamparB_full_tecp_r6, "TEC+ Ring 6 slope","","#Delta#theta_{B} [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(beamparB_full_tecp_r6), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_beampar_tecp_full->Print("beampar_tecp_full.png");

  TCanvas* cv_beampar_tecm_full = new TCanvas("beampar_tecm_full","TEC- Beam Parameters", 1500, 1000);
  cv_beampar_tecm_full->Divide(2,2);

  cv_beampar_tecm_full->cd(1);
  mgr = avec_draw(xvals, beamparA_full_tecm_r4, "TEC- Ring 4 offset","","#Delta#theta_{A} [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(beamparA_full_tecm_r4), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_beampar_tecm_full->cd(2);
  mgr = avec_draw(xvals, beamparA_full_tecm_r6, "TEC- Ring 6 offset","","#Delta#theta_{A} [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(beamparA_full_tecm_r6), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_beampar_tecm_full->cd(3);
  mgr = avec_draw(xvals, beamparB_full_tecm_r4, "TEC- Ring 4 slope","","#Delta#theta_{B} [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(beamparB_full_tecm_r4), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_beampar_tecm_full->cd(4);
  mgr = avec_draw(xvals, beamparB_full_tecm_r6, "TEC- Ring 6 slope","","#Delta#theta_{B} [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(beamparB_full_tecm_r6), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_beampar_tecm_full->Print("beampar_tecm_full.png");


  std::vector<std::string> glob_legend;
  glob_legend.push_back("#Delta#phi_{0}");
  glob_legend.push_back("#Delta#phi_{t}");
  glob_legend.push_back("#Deltax_{0}");
  glob_legend.push_back("#Deltax_{t}");
  glob_legend.push_back("#Deltay_{0}");
  glob_legend.push_back("#Deltay_{t}");

  if(process_rings){
    TCanvas* cv_globpar = new TCanvas("cv_globpar","Global TEC Parameters", 1500, 1000);
    cv_globpar->Divide(2,2);
    cv_globpar->cd(1);
    mgr = avec_draw(xvals, vtrans(globpar_tecp_r4), "TEC+ Ring 4","","#Delta#phi/#Deltaxy [#murad/#mum]",options);
    ymax = SetDiscParAxes(mgr);
    AddLegend(mgr, glob_legend, 2, 0.7, 0.7, 0.9, 0.9);
    cv_globpar->cd(2);
    mgr = avec_draw(xvals, vtrans(globpar_tecp_r6), "TEC+ Ring 6","","#Delta#phi/#Deltaxy [#murad/#mum]",options);
    ymax = SetDiscParAxes(mgr);
    AddLegend(mgr, glob_legend, 2, 0.7, 0.7, 0.9, 0.9);
    cv_globpar->cd(3);
    mgr = avec_draw(xvals, vtrans(globpar_tecm_r4), "TEC- Ring 4","","#Delta#phi/#Deltaxy [#murad/#mum]",options);
    ymax = SetDiscParAxes(mgr);
    AddLegend(mgr, glob_legend, 2, 0.7, 0.7, 0.9, 0.9);
    cv_globpar->cd(4);
    mgr = avec_draw(xvals, vtrans(globpar_tecm_r6), "TEC- Ring 6","","#Delta#phi/#Deltaxy [#murad/#mum]",options);
    ymax = SetDiscParAxes(mgr);
    AddLegend(mgr, glob_legend, 2, 0.7, 0.7, 0.9, 0.9);
    cv_globpar->Print("globpar.png");
  }

  //glob_err_tecp.print();
  //glob_err_tecm.print();

  std::vector<std::string> legend;
  legend.push_back("Disc 1");
  legend.push_back("Disc 2");
  legend.push_back("Disc 3");
  legend.push_back("Disc 4");
  legend.push_back("Disc 5");
  legend.push_back("Disc 6");
  legend.push_back("Disc 7");
  legend.push_back("Disc 8");
  legend.push_back("Disc 9");

  if(process_rings){
    TCanvas* cv_tecp_r4 = new TCanvas("cv_tecp_r4","Disc Parameters TEC+ Ring 4", 1500, 1000);
    cv_tecp_r4->Divide(2,2);
    
    cv_tecp_r4->cd(1);
    mgr = avec_draw(xvals, vtrans(tecp_r4_dphik), "TEC+ Ring 4 #Delta#phi","","#Delta#phi [#murad]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecp_r4_dphik), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_tecp_r4->cd(2);
    mgr = avec_draw(xvals, vtrans(tecp_r4_dxk), "TEC+ Ring 4 #Deltax ","","#Deltax [#mum]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecp_r4_dxk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_tecp_r4->cd(3);
    mgr = avec_draw(xvals, vtrans(tecp_r4_dyk), "TEC+ Ring 4 #Deltay ","","#Deltay [#mum]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecp_r4_dyk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_tecp_r4->cd(4);
    AddLegend(mgr, legend, 2, 0.3, 0.2, 0.8, 0.8);
    cv_tecp_r4->Print("discpar_tecp_r4.png");
    
    TCanvas* cv_tecp_r6 = new TCanvas("cv_tecp_r6","Disc Parameters TEC+ Ring 6", 1500, 1000);
    cv_tecp_r6->Divide(2,2);
    
    cv_tecp_r6->cd(1);
    mgr = avec_draw(xvals, vtrans(tecp_r6_dphik), "TEC+ Ring 6 #Delta#phi","","#Delta#phi [#murad]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecp_r6_dphik), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_tecp_r6->cd(2);
    mgr = avec_draw(xvals, vtrans(tecp_r6_dxk), "TEC+ Ring 6 #Deltax ","","#Deltax [#mum]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecp_r6_dxk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_tecp_r6->cd(3);
    mgr = avec_draw(xvals, vtrans(tecp_r6_dyk), "TEC+ Ring 6 #Deltay ","","#Deltay [#mum]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecp_r6_dyk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_tecp_r6->cd(4);
    AddLegend(mgr, legend, 2, 0.3, 0.2, 0.8, 0.8);
    cv_tecp_r6->Print("discpar_tecp_r6.png");
    
    
    TCanvas* cv_tecm_r4 = new TCanvas("cv_tecm_r4","Disc Parameters TEC- Ring 4", 1500, 1000);
    cv_tecm_r4->Divide(2,2);
    
    cv_tecm_r4->cd(1);
    mgr = avec_draw(xvals, vtrans(tecm_r4_dphik), "TEC- Ring 4 #Delta#phi","","#Delta#phi [#murad]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecm_r4_dphik), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_tecm_r4->cd(2);
    mgr = avec_draw(xvals, vtrans(tecm_r4_dxk), "TEC- Ring 4 #Deltax ","","#Deltax [#mum]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecm_r4_dxk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_tecm_r4->cd(3);
    mgr = avec_draw(xvals, vtrans(tecm_r4_dyk), "TEC- Ring 4 #Deltay ","","#Deltay [#mum]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecm_r4_dyk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);
    cv_tecm_r4->cd(4);
    AddLegend(mgr, legend, 2, 0.3, 0.2, 0.8, 0.8);
    cv_tecm_r4->Print("discpar_tecm_r4.png");

    
    std::cout << "Average RMS Dphik: " << calc_rms(tecm_r6_dphik) << std::endl;
    std::cout << "Average RMS   Dxk: " << calc_rms(tecm_r6_dxk) << std::endl;
    std::cout << "Average RMS   Dyk: " << calc_rms(tecm_r6_dyk) << std::endl;
    
    TCanvas* cv_tecm_r6 = new TCanvas("cv_tecm_r6","Disc Parameters TEC- Ring 6", 1500, 1000);
    cv_tecm_r6->Divide(2,2);
    
    cv_tecm_r6->cd(1);
    mgr = avec_draw(xvals, vtrans(tecm_r6_dphik), "TEC- Ring 6 #Delta#phi","","#Delta#phi [#murad]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecm_r6_dphik), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);
    
    cv_tecm_r6->cd(2);
    mgr = avec_draw(xvals, vtrans(tecm_r6_dxk), "TEC- Ring 6 #Deltax ","","#Deltax [#mum]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecm_r6_dxk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);
    
    cv_tecm_r6->cd(3);
    mgr = avec_draw(xvals, vtrans(tecm_r6_dyk), "TEC- Ring 6 #Deltay ","","#Deltay [#mum]",options);
    ymax = SetDiscParAxes(mgr);
    draw_rms_label( calc_rms(tecm_r6_dyk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);
    
    cv_tecm_r6->cd(4);
    AddLegend(mgr, legend, 2, 0.3, 0.2, 0.8, 0.8);
    cv_tecm_r6->Print("discpar_tecm_r6.png");
  }

  TCanvas* cv_tecp = new TCanvas("cv_tecp","Disc Parameters TEC+", 1500, 1000);
  cv_tecp->Divide(2,2);

  cv_tecp->cd(1);
  mgr = avec_draw(xvals, vtrans(tecp_dphik), "TEC+ #Delta#phi","","#Delta#phi [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(tecp_dphik), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_tecp->cd(2);
  mgr = avec_draw(xvals, vtrans(tecp_dxk), "TEC+ #Deltax ","","#Deltax [#mum]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(tecp_dxk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_tecp->cd(3);
  mgr = avec_draw(xvals, vtrans(tecp_dyk), "TEC+ #Deltay ","","#Deltay [#mum]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(tecp_dyk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_tecp->cd(4);
  AddLegend(mgr, legend, 2, 0.3, 0.2, 0.8, 0.8);
  cv_tecp->Print("discpar_tecp.png");

  TCanvas* cv_tecm = new TCanvas("cv_tecm","Disc Parameters TEC-", 1500, 1000);
  cv_tecm->Divide(2,2);

  cv_tecm->cd(1);
  mgr = avec_draw(xvals, vtrans(tecm_dphik), "TEC- #Delta#phi","","#Delta#phi [#murad]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(tecm_dphik), "#murad", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_tecm->cd(2);
  mgr = avec_draw(xvals, vtrans(tecm_dxk), "TEC- #Deltax ","","#Deltax [#mum]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(tecm_dxk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_tecm->cd(3);
  mgr = avec_draw(xvals, vtrans(tecm_dyk), "TEC- #Deltay ","","#Deltay [#mum]",options);
  ymax = SetDiscParAxes(mgr);
  draw_rms_label( calc_rms(tecm_dyk), "#mum", xlo, ymax*4/7, xhi, ymax*6/7);

  cv_tecm->cd(4);
  AddLegend(mgr, legend, 2, 0.3, 0.2, 0.8, 0.8);
  cv_tecm->Print("discpar_tecm.png");


  TCanvas* cv_globpar_full = new TCanvas("cv_globpar_full","Global TEC Parameters (full reconstruction)", 700, 1000);
  cv_globpar_full->Divide(1,2);
  cv_globpar_full->cd(1);
  //mgr = avec_draw(xvals, globpar_tecp, "TEC+ (full)","","#Delta#phi/#Deltaxy [#murad/#mum]",options);
  mgr = avec_draw(xvals, globpar_tecp, glob_err_tecp, "TEC+ (full)","","#Delta#phi/#Deltaxy [#murad/#mum]",options, 1, "mgr_tec_globpar");
  ymax = SetDiscParAxes(mgr);
  AddLegend(mgr, glob_legend, 2, 0.7, 0.7, 0.9, 0.9);
  cv_globpar_full->cd(2);
  mgr = avec_draw(xvals, globpar_tecm, glob_err_tecm, "TEC- (full)","","#Delta#phi/#Deltaxy [#murad/#mum]",options);
  ymax = SetDiscParAxes(mgr);
  AddLegend(mgr, glob_legend, 2, 0.7, 0.7, 0.9, 0.9);
  cv_globpar_full->Print("globpar_full.png");




}

void alpar_history_fine(const std::vector<std::string>& file_list)
{
  if(file_list.empty()){
    std::cerr << "Error in alpar_history_fine: file_list is empty!" << std::endl;
    return;
  }

  bool process_rings = false;
  bool bs_frame = true;

  std::vector<LasAlPar> alpar_list;
  Avec xvals;

  double b1,b2;
  Avec::size_type b3;
  // Loop over all files
  for(std::vector<std::string>::size_type i = 0; i < file_list.size(); i++){
    std::cout << "Processing file " << file_list[i] << std::endl;
    if(!run_selector(file_list[i], b1, b3, b2))continue;

    Avec timev = avec_get("block_timestamp", file_list[i]);
    Avec good_block = avec_get("good_block", file_list[i]);
    Avec good_fit = avec_get("tec_rec/good_fit", file_list[i]);

    if (timev.empty()){
      std::cerr << "Could not find Avec block_timestamp in file " << file_list[i] << std::endl;
      continue;
    }
    if (good_block.empty()){
      std::cerr << "Could not find Avec good_block in file " << file_list[i] << std::endl;
      continue;
    }
    if (good_fit.empty()){
      std::cerr << "Could not find Avec tec_rec/good_fit in file " << file_list[i] << std::endl;
      continue;
    }
    if(good_block.size() != timev.size() || good_block.size() != good_fit.size()){
      std::cerr << "Error, good_block.size() = " << good_block.size() << "  timev.size() = " << timev.size() << "   mismatch" << std::endl;
      continue;
    }

    // Loop over all blocks
    for(Avec::size_type f=0; f < timev.size(); f++){
      if(good_block[f] != 1 || good_fit[f] != 1)continue;
      std::ostringstream alpar_name;
      alpar_name << "tec_rec/alpar_" << f;

      LasAlPar alpar = alpar_get(alpar_name.str(), file_list[i]);
      alpar_list.push_back(alpar);
      xvals.push_back(timev[f]);
    }
  }
  
  alpar_history_draw(alpar_list, xvals, bs_frame, process_rings);

}

void control_positions(const std::string& filename)
{
  TFile f(filename.c_str(), "READ");
  Avec block_nr = avec_get("block_nr", f);
  Avec block_timestamp = avec_get("block_timestamp", f);
  if(block_nr.empty() || block_timestamp.empty()){
    std::cerr << "Error in control_positions: not all objects were found" << std::endl;
    return;
  }

  LASGlobalData<Avec> xval;
  LASGlobalData<Avec> pos_history;
  LASGlobalData<Avec> err_history;
  LASGlobalData<Avec> rms_history;
  //LASGlobalData<Avec> mask_history;

  for(Avec::size_type i = 0; i < block_nr.size(); i++){
    std::cout << "Processing block " << i << std::endl;

    std::ostringstream norm_name;
    std::ostringstream pos_name;
    std::ostringstream rms_av_name;
    std::ostringstream mask_name;
    norm_name << "norm_" << i;
    pos_name << "positions_" << i;
    //rms_av_name << "rms_av_" << i;
    rms_av_name << "pos_error_" << i;
    mask_name << "positions_mask_" << i;
    LASGlobalData<int> norm = global_data_get<int>(norm_name.str(), f);
    LASGlobalData<double> pos = global_data_get<double>(pos_name.str(), f);
    LASGlobalData<double> rms_av = global_data_get<double>(rms_av_name.str(), f);
    LASGlobalData<int> mask = global_data_get<int>(mask_name.str(), f);
    LASGlobalData<double> err = sqrt(rms_av*rms_av/norm);

    LASGlobalDataLoop loop;
    do{
      if(loop.GetEntry(mask)){
	loop.GetEntry(xval).push_back(block_timestamp[i]);
	loop.GetEntry(pos_history).push_back(loop.GetEntry<double>(pos));
	loop.GetEntry(rms_history).push_back(loop.GetEntry<double>(rms_av));
	loop.GetEntry(err_history).push_back(loop.GetEntry<double>(err));
	//loop.GetEntry(mask_history).push_back(loop.GetEntry(mask));
      }
    }while(loop.next());
  }
  f.Close();

  Avec v_avg_err;
  Avec v_run_rms;
  LASGlobalDataLoop checkloop;
  do{
    std::cout << GetModuleName(checkloop);
    Avec& err = checkloop.GetEntry(err_history);
    Avec& pos = checkloop.GetEntry(pos_history);
    if(! err.empty()){
      //pos -= vsum(pos)/pos.size();
      double average_err = (vsum(err)/err.size());
      double run_rms = (sqrt((vsum(pos*pos) - vsum(pos) * vsum(pos)/pos.size())/(pos.size()-1)));
      std::cout << ": has average err of " << average_err;
      std::cout << "  and run rms of " << run_rms;
      v_avg_err.push_back(average_err);
      v_run_rms.push_back(run_rms);
    }
    std::cout << std::endl;
  }while (checkloop.next());

  Avec::VERBOSE_FLAG = 0;
  draw_global_data(pos_history, err_history, xval, "AP", "time");
  Avec::VERBOSE_FLAG = 1;
  new TCanvas("err_corr", "Correlations of Fluctuations");
  //avec_draw(v_avg_err, v_run_rms, "", "Average block error", "Run rms");
  avec_plot(v_avg_err, v_run_rms);
}

void control_pos_av(const std::string& filename)
{
  TFile f(filename.c_str(), "READ");
  Avec block_nr = avec_get("block_nr", f);
  Avec block_timestamp = avec_get("block_timestamp", f);
  if(block_nr.empty() || block_timestamp.empty()){
    std::cerr << "Error in control_pos_av: not all objects were found" << std::endl;
    return;
  }

  LASGlobalData<Avec> pos_history;
  LASGlobalData<Avec> err_history;
  LASGlobalData<Avec> rms_history;

  for(Avec::size_type i = 0; i < block_nr.size(); i++){
    std::cout << "Processing block " << i << std::endl;

    std::ostringstream norm_name;
    std::ostringstream pos_av_name;
    std::ostringstream rms_av_name;
    norm_name << "norm_" << i;
    pos_av_name << "pos_av_" << i;
    rms_av_name << "rms_av_" << i;
    LASGlobalData<int> norm = global_data_get<int>(norm_name.str(), f);
    LASGlobalData<double> pos_av = global_data_get<double>(pos_av_name.str(), f);
    LASGlobalData<double> rms_av = global_data_get<double>(rms_av_name.str(), f);
    LASGlobalData<double> err = sqrt(rms_av*rms_av/norm);

    LASGlobalDataLoop loop;
    do{
      loop.GetEntry<Avec>(pos_history).push_back(loop.GetEntry<double>(pos_av));
      loop.GetEntry<Avec>(rms_history).push_back(loop.GetEntry<double>(rms_av));
      loop.GetEntry<Avec>(err_history).push_back(loop.GetEntry<double>(err));
    }while(loop.next());
  }
  f.Close();

  draw_global_data(pos_history, err_history, block_timestamp, "AP", "time");
}

void control_rms_av(const std::string& filename)
{
  TFile f(filename.c_str(), "READ");
  Avec block_nr = avec_get("block_nr", f);
  Avec block_timestamp = avec_get("block_timestamp", f);
  if(block_nr.empty() || block_timestamp.empty()){
    std::cerr << "Error in control_rms_av: not all objects were found" << std::endl;
    return;
  }

  LASGlobalData<Avec> pos_history;
  LASGlobalData<Avec> rms_history;

  for(Avec::size_type i = 0; i < block_nr.size(); i++){
    std::cout << "Processing block " << i << std::endl;

    std::ostringstream norm_name;
    std::ostringstream pos_av_name;
    std::ostringstream rms_av_name;
    norm_name << "norm_" << i;
    pos_av_name << "pos_av_" << i;
    rms_av_name << "rms_av_" << i;
    LASGlobalData<double> pos_av = global_data_get<double>(pos_av_name.str(), f);
    LASGlobalData<double> rms_av = global_data_get<double>(rms_av_name.str(), f);
    LASGlobalDataLoop loop;
    do{
      loop.GetEntry<Avec>(pos_history).push_back(loop.GetEntry<double>(pos_av));
      loop.GetEntry<Avec>(rms_history).push_back(loop.GetEntry<double>(rms_av));
    }while(loop.next());
  }
  f.Close();

  LASGlobalData<Avec> xval(block_timestamp);
  //draw_global_data(pos_history, xval, "AP");
  draw_global_data(rms_history, xval, "AP");
}


//! Control plots for output of Run_History
void control_history(const std::string& results_file)
{
  Avec block_timestamp = avec_get("block_timestamp", results_file);
  LASGlobalData<Avec> diff_history = global_data_get<Avec>("diff_history", results_file);

  int nblocks = block_timestamp.size();

  std::cout << "There are " << nblocks << " blocks in this run" << std::endl;

  Avec2D tecp_r4_dphik;
  Avec2D tecp_r4_dxk;
  Avec2D tecp_r4_dyk;
  Avec2D tecp_r6_dphik;
  Avec2D tecp_r6_dxk;
  Avec2D tecp_r6_dyk;
  Avec2D tecm_r4_dphik;
  Avec2D tecm_r4_dxk;
  Avec2D tecm_r4_dyk;
  Avec2D tecm_r6_dphik;
  Avec2D tecm_r6_dxk;
  Avec2D tecm_r6_dyk;

  Avec2D globpar_tecp_r4;
  Avec2D globpar_tecp_r6;
  Avec2D globpar_tecm_r4;
  Avec2D globpar_tecm_r6;

  Avec2D beampar_tecp_r4;
  Avec2D beampar_tecp_r6;
  Avec2D beampar_tecm_r4;
  Avec2D beampar_tecm_r6;

  TFile f(results_file.c_str(),"READ");
  std::vector<LasAlPar> par_list;

  f.cd("tec_rec_run");
  for(int i = 0; i < nblocks; i++){
    std::ostringstream alpar_name;
    alpar_name << "tec_rec_run/alpar_" << i;

    LasAlPar* alpar_ptr = 0;
    f.GetObject(alpar_name.str().c_str(), alpar_ptr);
    if(!alpar_ptr){
      std::cerr << "Did not find " << alpar_name.str() << std::endl;
      par_list.push_back(LasAlPar());
      continue;
      //return;
    }
    par_list.push_back(*alpar_ptr);

  }
  f.Close();


  alpar_history_draw(par_list, block_timestamp);

  return;

  LASGlobalData<Avec> xval(block_timestamp);
  //draw_global_data(pos_history, xval, "AP");
  draw_global_data(diff_history, xval, "AP");

}

void control_compl_sig(const std::string& filename){
  TFile f(filename.c_str(), "READ");
  if( ! f.IsOpen() ){
    std::cerr << "Could not open file " << filename << std::endl;
    return;
  }

  Avec comp_sig_acc_tec = avec_get("comp_sig_acc_tec", f);
  Avec comp_sig_rej_tec = avec_get("comp_sig_rej_tec", f);
  Avec comp_sig_acc_tib = avec_get("comp_sig_acc_tib", f);
  Avec comp_sig_rej_tib = avec_get("comp_sig_rej_tib", f);
  Avec comp_sig_acc_tob = avec_get("comp_sig_acc_tob", f);
  Avec comp_sig_rej_tob = avec_get("comp_sig_rej_tob", f);

  f.Close();

  TCanvas* cv = new TCanvas("cv_compl_sig","Complementary Signal");
  cv->Divide(3,2);
  cv->cd(1);
  avec_plot(comp_sig_acc_tec, 200, "TEC Accepted");
  cv->cd(4);
  avec_plot(comp_sig_rej_tec, 200, "TEC Rejected");
  cv->cd(2);
  avec_plot(comp_sig_acc_tib, 200, "TIB Accepted");
  cv->cd(5);
  avec_plot(comp_sig_rej_tib, 200, "TIB Rejected");
  cv->cd(3);
  avec_plot(comp_sig_acc_tob, 200, "TOB Accepted");
  cv->cd(6);
  avec_plot(comp_sig_rej_tob, 200, "TOB Rejected");

}


void control_intensities(const std::string& result_filename)
{
  Avec AT_event = avec_get("AT_event", result_filename);
  Avec TEC_event = avec_get("TEC_event", result_filename);

  Avec signal_at = avec_get("signal_at",result_filename);
  Avec signal_tec_r4 = avec_get("signal_tec_r4",result_filename);
  Avec signal_tec_r6 = avec_get("signal_tec_r6",result_filename);
  Avec signal_tec_disc1 = avec_get("signal_tec_disc1",result_filename);
  Avec signal_tec_far = avec_get("signal_tec_far",result_filename);

  Avec dist_at = avec_get("dist_at", result_filename);
  Avec cuts_at = avec_get("cuts_at", result_filename);
  Avec dist_tec = avec_get("dist_tec", result_filename);
  Avec cuts_tec = avec_get("cuts_tec", result_filename);

  if(AT_event.empty() || 
     TEC_event.empty() || 
     signal_at.empty() || 
     signal_tec_r4.empty() || 
     signal_tec_r6.empty() || 
     signal_tec_disc1.empty() || 
     signal_tec_far.empty() ||
     dist_at.empty() ||
     cuts_at.empty() ||
     dist_tec.empty() ||
     cuts_tec.empty()
     ){
    std::cerr << "Could not get all Avecs" << std::endl;
    return;
  }

  TCanvas* cv4 = new TCanvas("cv_int1","Intensities correlations AT");
  cv4->Divide(2,2);
  cv4->cd(1);
  avec_plot_mask(signal_at, signal_tec_r6, AT_event, 100, 100, "TEC R6 % AT");
  cv4->cd(2);
  avec_plot_mask(signal_tec_r6, signal_tec_disc1, AT_event, 100, 100, "TEC Disc 1 % TEC_R6");
  cv4->cd(3);
  avec_plot_mask(signal_tec_disc1, signal_at, AT_event, 100, 100, "AT % TEC Disc 1");
  cv4->cd(4)->SetLogy();
  avec_plot_mask(dist_at, AT_event, 1000);
  for(Avec::size_type i = 0; i < cuts_at.size(); i++){
    TLine* l = new TLine(cuts_at[i], 0, cuts_at[i], 100);
    l->SetLineColor(2);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->Draw();
  }

  TCanvas* cv5 = new TCanvas("cv_insignal_tec_r6","Intensities correlations TEC");
  cv5->Divide(2,2);
  cv5->cd(1);
  avec_plot_mask(signal_tec_r4, signal_tec_r6, TEC_event, 100, 100, "TEC R6 % TEC R4");
  cv5->cd(2);
  avec_plot_mask(signal_tec_r6, signal_tec_far, TEC_event, 100, 100, "TEC FAR % TEC R6");
  cv5->cd(3);
  avec_plot_mask(signal_tec_far, signal_tec_r4, TEC_event, 100, 100, "TEC R4 % TEC FAR");
  cv5->cd(4)->SetLogy();
  avec_plot_mask(dist_tec, TEC_event, 1000);
  for(Avec::size_type  i = 0; i < cuts_tec.size(); i++){
    TLine* l = new TLine(cuts_tec[i], 0, cuts_tec[i], 100);
    l->SetLineColor(2);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->Draw();
  }

  return;
}

void control_blocks(const std::string& output_filename, bool save_plots)
{
  Avec event_diff = avec_get("event_diff", output_filename);
  Avec event_block = avec_get("event_block", output_filename);
  Avec block_nr = avec_get("block_nr", output_filename);
  Avec block_timestamp = avec_get("block_timestamp", output_filename);
  Avec block_size = avec_get("block_size", output_filename);
  Avec block_threshold = avec_get("block_threshold", output_filename);

  if(event_diff.empty() || event_block.empty() || block_nr.empty() || block_size.empty() || block_threshold.empty()){
    std::cerr << "Error in control_blocks: Could not find all Avecs" << std::endl;
    return;
  }

  Avec::size_type size = event_diff.size();
  double threshold = block_threshold[0];

  TCanvas* cv = new TCanvas("cv_block_slicer","Canvas Block Slicer");
  cv->Divide(2,2);
  cv->cd(1)->SetLogy();
  TGraph* gr = avec_draw(event_diff, "Event Nr. Differences", "entry nr.", "evt. nr. diff", "AP", 1, "event_diff");
  double max_yaxis = std::max( pow(10, ceil( log10( vmax(event_diff) ) ) ), threshold*1.1);
  //std::cout << "maxy: " << max_yaxis << std::endl;
  gr->GetYaxis()->SetRangeUser(1,max_yaxis);
  gr->GetXaxis()->SetNdivisions( 505, kTRUE);

  TLine* l3 = new TLine(0, threshold, size, threshold);
  l3->SetLineColor(2);
  l3->SetLineWidth(3);
  l3->SetLineStyle(2);
  l3->Draw();
    
  cv->cd(2);
  event_block = (event_block > -1)*(event_block +1) - 1.0;
  gr = avec_draw(event_block, "Block Assignment", "entry nr.", "Block Nr.", "AP");
  gr->GetXaxis()->SetNdivisions( 505, kTRUE);

  cv->cd(3)->SetLogy();
  //avec_plot(log10(event_diff)(1,size));

  cv->cd(4);
  gr = avec_draw(block_timestamp, block_size, "Blocks", "time", "block size", "AP");
  gr->GetXaxis()->SetTimeDisplay(1);
  gr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 

  if(save_plots){
    TFile f(output_filename.c_str(), "UPDATE");
    cv->Write();
    f.Close();
  }

  return;
}

//! Control plots for output of get_max_signal
void control_max_signal(const std::string& max_file, bool full_output)
{
  // Object with signal maxima
  LASGlobalData<Avec>* max_signal=0;
  TFile f(max_file.c_str(),"READ");

  if(full_output){
    f.GetObject("max_signal",max_signal);
    if(!max_signal) throw LAS::Exception("Did not find max_signal object");
  }

  Avec eventnumber = avec_get("eventnumber", f);
  Avec unixTime = avec_get("unixTime", f);
  Avec signal_tecplus = avec_get("signal_tecplus", f);
  Avec signal_tecminus = avec_get("signal_tecminus", f);
  Avec signal_tec_disc1 = avec_get("signal_tec_disc1", f);
  Avec signal_at = avec_get("signal_at", f);
  Avec signal_all = avec_get("signal_all", f);
  Avec signal_tec_r4_nooverlap = avec_get("signal_tec_r4_nooverlap", f);
  Avec signal_at_nooverlap = avec_get("signal_at_nooverlap", f);
  Avec signal_tec_close = avec_get("signal_tec_close", f);
  Avec signal_tec_far = avec_get("signal_tec_far", f);
  Avec signal_tec_r6 = avec_get("signal_tec_r6", f);
  Avec signal_tec_r4 = avec_get("signal_tec_r4", f);

  if(
     eventnumber.empty() || 
     unixTime.empty() ||
     signal_tecplus.empty() || 
     signal_tecminus.empty() || 
     signal_tec_disc1.empty() || 
     signal_at.empty() || 
     signal_all.empty() ||
     signal_tec_r4_nooverlap.empty() ||
     signal_at_nooverlap.empty() ||
     signal_tec_close.empty() ||
     signal_tec_far.empty() ||
     signal_tec_r6.empty() ||
     signal_tec_r4.empty()
     ){
    throw LAS::Exception("Could not find all Avecs");
  }

  f.Close();

  if(0){
    std::cout << "unixTime.size(): " << unixTime.size() << std::endl;    
    Avec mask = fabs(signal_tec_r4_nooverlap - signal_at_nooverlap/3*2) > 600;
    unixTime = avec_mask(unixTime, mask);
    eventnumber = avec_mask(eventnumber, mask);
    signal_tecplus = avec_mask(signal_tecplus, mask);
    signal_tecminus = avec_mask(signal_tecminus, mask);
    signal_at = avec_mask(signal_at, mask);
    //signal_at_nooverlap = avec_mask(signal_at_nooverlap, mask);
    //signal_tec_r4_nooverlap = avec_mask(signal_tec_r4_nooverlap, mask);
    
    //avec_mask(, mask);
    std::cout << "unixTime.size(): " << unixTime.size() << std::endl;
  }

  TMultiGraph* mgr;
  std::vector<std::string> legend;
  legend.push_back("Endcap +");
  legend.push_back("Endcap -");
  legend.push_back("Alignment Tubes");

  new TCanvas("cv_time", "Signals Time");
  mgr = avec_draw(unixTime, (signal_tecplus | signal_tecminus | signal_at), "Max Signal", "time", "signal", "AP");
  mgr->GetXaxis()->SetNdivisions(505, kTRUE);
  mgr->GetXaxis()->SetTimeDisplay(1);
  mgr->GetXaxis()->SetTimeFormat("%H:%M"); 
  AddLegend(mgr, legend);

  new TCanvas("cv_evtnr", "Signals Evtnr");
  avec_draw(eventnumber, (signal_tecplus | signal_tecminus | signal_at), "Max Signal", "evt ntr.", "signal", "AP");
  AddLegend(mgr, legend);

  TCanvas* cv_tmp = new TCanvas("cv_tmp","Temporary Canvas");
  int nbins = 500;
  THStack* stack = new THStack("thestack", "thestack");
  TH1* h;
  h = avec_plot(signal_all, nbins);
  h->SetLineColor(1);
  stack->Add(h);
  h = avec_plot(signal_at, nbins);
  h->SetLineColor(2);
  stack->Add(h);
  h = avec_plot(signal_tecplus, nbins);
  h->SetLineColor(3);
  stack->Add(h);
  h = avec_plot(signal_tecminus, nbins);
  h->SetLineColor(4);
  stack->Add(h);
  h = avec_plot(signal_tec_disc1, nbins);
  h->SetLineColor(6);
  stack->Add(h);
  h = avec_plot(signal_tec_r4, nbins);
  h->SetLineColor(7);
  stack->Add(h);
  h = avec_plot(signal_tec_r6, nbins);
  h->SetLineColor(8);
  stack->Add(h);
  h = avec_plot(signal_tec_close, nbins);
  h->SetLineColor(9);
  stack->Add(h);
  h = avec_plot(signal_tec_far, nbins);
  h->SetLineColor(11);
  stack->Add(h);
  h = avec_plot(signal_tec_r4_nooverlap, nbins);
  h->SetLineColor(41);
  stack->Add(h);
  h = avec_plot(signal_at_nooverlap, nbins);
  h->SetLineColor(46);
  stack->Add(h);
  delete cv_tmp;

  (new TCanvas("cv_sigsum","Signal Sum"))->SetLogy();
  //stack->Draw("nostack");
  stack->Draw();

  TCanvas* cv_dist = new TCanvas("cv_dist", "Dist");
  cv_dist->Divide(2,2);
  cv_dist->cd(1);
  avec_plot(signal_at_nooverlap);
  cv_dist->cd(2);
  avec_plot(signal_tec_r4_nooverlap);
  cv_dist->cd(3);
  avec_plot(signal_at_nooverlap, signal_tec_r4_nooverlap);
  cv_dist->cd(4);
  avec_plot(signal_tec_r4_nooverlap - signal_at_nooverlap/3*2, 400);


  if(full_output) draw_global_data(*max_signal,"AP"); // This will be too much for large runs (nr. of events > 2000)

}


void control_integrated_pos(const std::string& filename)
{

  // Open the file
  TFile f(filename.c_str(),"UPDATE");
  if(!f.IsOpen()){
    std::cerr << "Could not open file " << filename << std::endl;
    return;
  }

  LASGlobalData<double>* positions_strips = 0;
  f.GetObject("positions_strips",positions_strips);
  if(!positions_strips){
    std::cerr << "Did not find positions_strips" << std::endl;
    return;
  }

  LASGlobalData<int>* results_mask = 0;
  f.GetObject("results_mask",results_mask);
  if(!results_mask){
    std::cerr << "Did not find results_mask" << std::endl;
    return;
  }

  draw_global_data(*positions_strips, *results_mask);

}

// void control_profiles(const std::string& filename, int block_nr)
// {
//   // Open the file
//   TFile f(filename.c_str(),"READ");
//   if(!f.IsOpen()){
//     std::cerr << "Could not open file " << filename << std::endl;
//     return;
//   }


//   LASGlobalData<Avec>* profiles = 0;
//   f.GetObject(prof_name.str(), profiles);
//   if(!profiles){
//     std::cerr << "Did not find "  << prof_name.str() << std::endl;
//     return;
//   }

//   draw_global_data(*profiles,"AL");

//   return;
// }

void control_profiles(const std::string& filename, int block_nr)
{
  //if(block_nr < 0) return control_profiles(filename);

  std::ostringstream prof_name("profiles");
  if( block_nr >= 0) prof_name << "profiles_" << block_nr;

  Avec xval(512, 0, 511);

  //std::ostringstream prof_name;
  //prof_name << "profiles_" << block_nr; 
  int flag = Avec::VERBOSE_FLAG;
  Avec::VERBOSE_FLAG = 0;
  draw_global_data(global_data_get<Avec>( prof_name.str(), filename ), global_data_get<Avec>("striprms_0", filename) , xval, "AL");
  Avec::VERBOSE_FLAG = flag;

  return;
}


void control_striprms(const std::string& filename)
{
  // Open the file
  TFile f(filename.c_str(),"READ");
  if(!f.IsOpen()){
    std::cerr << "Could not open file " << filename << std::endl;
    return;
  }

  LASGlobalData<Avec>* striprms = 0;
  f.GetObject("striprms_0", striprms);
  if(!striprms){
    std::cerr << "Did not find striprms" << std::endl;
    return;
  }

  draw_global_data(*striprms,"AL");

  return;
}

void control_step_mask(const std::string& filename)
{
  LASGlobalData<int>* step_mask=0;
  TFile f(filename.c_str(),"READ");
  f.GetObject("step_mask",step_mask);
  if(!step_mask){
    std::cerr << "Did not find step_mask object" << std::endl;
    return;
  }
  f.Close();

  print_global_data(*step_mask);

}

void control_step(const std::string& filename, int step)
{
  Avec label_step = avec_get("label_step", filename);
  Avec::size_type nr_entries = label_step.size();
  if(nr_entries == 0){
    std::cerr << "Error, could not find label_step" << std::endl;
    return;
  }

  // Object with signal maxima
  LASGlobalData<Avec>* max_signal=0;
  TFile f(filename.c_str(),"READ");
  f.GetObject("max_signal",max_signal);
  if(!max_signal){
    std::cerr << "Did not find max_signal object" << std::endl;
    return;
  }
  LASGlobalData<Avec> work_copy = *max_signal;
  f.Close();


  Avec mask = label_step == step;

  std::cout << vsum(mask) << " entries for step " << step << std::endl;

  LAS::beam_group group = (step <= 5 ? LAS::TEC : LAS::AT);


  LASGlobalDataLoop loop(convert_to_loop_type(group));
  do{
    loop.GetEntry<Avec>(work_copy) *= mask;
  }while (loop.next());

  draw_global_data(work_copy, "AP", group);
}
