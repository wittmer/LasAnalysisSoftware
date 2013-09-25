#include "LAS_globaldata_tools.h"

#define __AVECROOT__
#include "Avec.h"

#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TChain.h"
#include "TROOT.h"

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"

// Specializtion for drawing Avecs
template <>
void global_data_draw<Avec>(const LASGlobalData<Avec>& data)
{
  std::cout << "global_data_draw<Avec>" << std::endl;
  draw_global_data(data);
}

template void global_data_draw<Avec>(const LASGlobalData<Avec>&);

// Specializtion for drawing doubles
template <>
void global_data_draw<double>(const LASGlobalData<double>& data)
{
  std::cout << "global_data_draw<double>" << std::endl;
  draw_global_data(data);
}

template void global_data_draw<double>(const LASGlobalData<double>&);

// Specializtion for drawing floats
template <>
void global_data_draw<float>(const LASGlobalData<float>& data)
{
  std::cout << "global_data_draw<float>" << std::endl;
  draw_global_data(data, "");
}

template void global_data_draw<float>(const LASGlobalData<float>&);

// Specializtion for drawing ints
template <>
void global_data_draw<int>(const LASGlobalData<int>& data)
{
  std::cout << "global_data_draw<int>" << std::endl;
  draw_global_data<int>(data);
}

template void global_data_draw<int>(const LASGlobalData<int>&);

//! Generate random values following a gaussian distribution with mean = 0 and sigma = 1
LASGlobalData<double> global_data_random_gauss()
{
  // So far only TEC internal values are implemented
  Avec random = grand(288, 1.0); // There are 288 TEC modules
  Avec::iterator it = random.begin();
  LASGlobalData<double> random_values;
  LASGlobalDataLoop loop(LASGlobalDataLoop::TEC);
  do{
    loop.GetEntry(random_values) = *it;
  } while (loop.next() && ++it != random.end());
  return random_values;
}

//! Take global data object with Avecs and compute mean, rms and number of entries for each Avec
void global_data_stat(const LASGlobalData<Avec>& data, LASGlobalData<double>& mean, LASGlobalData<double>& rms, LASGlobalData<int>& count)
{
  LASGlobalDataLoop loop;
  do{
    const Avec& buffer = loop.GetEntry<Avec>(data);
    Avec::size_type size = buffer.size();
    double sum = vsum(buffer);
    double sumsq = vsum(buffer * buffer);
    loop.GetEntry<double>(mean) = sum/size;
    loop.GetEntry<double>(rms) = sqrt((sumsq - sum * sum / size)/(size - 1));
    loop.GetEntry<int>(count) = (int)size;
  }while(loop.next());
}

//! Draw global data values against z-positions
TMultiGraph* draw_global_data_z(const LASGlobalData<double>& data, const LASGlobalData<double>& error, const LASGlobalData<int>& mask, LAS::beam_group group, const std::string& canvas_name)
{
  TMultiGraph* retval = 0;

  std::vector<std::string> legend;
  legend.push_back("beam 0");
  legend.push_back("beam 1");
  legend.push_back("beam 2");
  legend.push_back("beam 3");
  legend.push_back("beam 4");
  legend.push_back("beam 5");
  legend.push_back("beam 6");
  legend.push_back("beam 7");


  if(group != LAS::TEC){
    std::vector<Avec> xvals(8);
    std::vector<Avec> yvals(8);
    std::vector<Avec> err(8);
    
    LASGlobalDataLoop loop(LASGlobalDataLoop::AT);
    //LASGlobalDataLoop loop(convert_to_loop_type(group));
    do{
      if(loop.GetEntry(mask) != 1) continue;
      xvals[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
      yvals[loop.get_beam()].push_back( loop.GetEntry(data) );
      err[loop.get_beam()].push_back( loop.GetEntry(error) );
      
    }while (loop.next());

    std::string canvas_name_at = canvas_name + "_AT";
    new TCanvas(canvas_name_at.c_str(), canvas_name_at.c_str());
    
    TMultiGraph* mgr = avec_draw(xvals, yvals, err, "Alignment Tubes","zpos [mm]","","AP");
    AddLegend(mgr, legend, 2, 0.6, 0.8, 1, 1);
    retval = mgr;
  }

  if(group == LAS::TEC || group == LAS::ALL){
    std::vector<Avec> xvals_tecp_r4(8);
    std::vector<Avec> yvals_tecp_r4(8);
    std::vector<Avec> xvals_tecp_r6(8);
    std::vector<Avec> yvals_tecp_r6(8);
    std::vector<Avec> xvals_tecm_r4(8);
    std::vector<Avec> yvals_tecm_r4(8);
    std::vector<Avec> xvals_tecm_r6(8);
    std::vector<Avec> yvals_tecm_r6(8);
    
    LASGlobalDataLoop loop(LASGlobalDataLoop::TEC);
    do{
      if(loop.GetEntry(mask) == 0) continue;
      switch(loop.get_det()){
      case 0:
	switch(loop.get_ring()){
	case 0:
	  xvals_tecp_r4[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
	  yvals_tecp_r4[loop.get_beam()].push_back( loop.GetEntry(data) );
	  break;
	case 1:
	  xvals_tecp_r6[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
	  yvals_tecp_r6[loop.get_beam()].push_back( loop.GetEntry(data) );
	  break;
	}
	break;
      case 1:
	switch(loop.get_ring()){
	case 0:
	  xvals_tecm_r4[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
	  yvals_tecm_r4[loop.get_beam()].push_back( loop.GetEntry(data) );
	  break;
	case 1:
	  xvals_tecm_r6[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
	  yvals_tecm_r6[loop.get_beam()].push_back( loop.GetEntry(data) );
	  break;
	}
	break;
      }
      
    }while (loop.next());

    TCanvas* cv_tec = new TCanvas(canvas_name.c_str(), canvas_name.c_str());
    cv_tec->Divide(2,2);
    TMultiGraph* mgr; 
    cv_tec->cd(1);
    mgr = avec_draw(xvals_tecp_r4, yvals_tecp_r4, "TEC+ Ring 4","zpos [mm]","","AP");
    cv_tec->cd(2);
    mgr = avec_draw(xvals_tecp_r6, yvals_tecp_r6, "TEC+ Ring 6","zpos [mm]","","AP");
    AddLegend(mgr, legend, 2, 0.6, 0.8, 1, 1);
    cv_tec->cd(3);
    mgr = avec_draw(xvals_tecm_r4, yvals_tecm_r4, "TEC- Ring 4","zpos [mm]","","AP");
    cv_tec->cd(4);
    mgr = avec_draw(xvals_tecm_r6, yvals_tecm_r6, "TEC- Ring 6","zpos [mm]","","AP");
    retval = mgr;
  }

  return retval;
}


//! Draw global data values against z-positions
TMultiGraph* draw_global_data_z(const LASGlobalData<double>& data, const LASGlobalData<int>& mask, LAS::beam_group group, const std::string& canvas_name)
{
  TMultiGraph* retval = 0;

  std::vector<std::string> legend;
  legend.push_back("beam 0");
  legend.push_back("beam 1");
  legend.push_back("beam 2");
  legend.push_back("beam 3");
  legend.push_back("beam 4");
  legend.push_back("beam 5");
  legend.push_back("beam 6");
  legend.push_back("beam 7");


  if(group != LAS::TEC){
    std::vector<Avec> xvals(8);
    std::vector<Avec> yvals(8);
    
    //LASGlobalDataLoop loop(LASGlobalDataLoop::AT);
    LASGlobalDataLoop loop(convert_to_loop_type(group));
    do{
      if(loop.GetEntry(mask) != 1) continue;
      xvals[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
      yvals[loop.get_beam()].push_back( loop.GetEntry(data) );
      
    }while (loop.next());
    
    new TCanvas(canvas_name.c_str(), canvas_name.c_str());
    
    TMultiGraph* mgr = avec_draw(xvals, yvals, "Alignment Tubes","zpos [mm]","","AP");
    AddLegend(mgr, legend, 2);
    retval = mgr;
  }

  if(group == LAS::TEC){
    std::vector<Avec> xvals_tecp_r4(8);
    std::vector<Avec> yvals_tecp_r4(8);
    std::vector<Avec> xvals_tecp_r6(8);
    std::vector<Avec> yvals_tecp_r6(8);
    std::vector<Avec> xvals_tecm_r4(8);
    std::vector<Avec> yvals_tecm_r4(8);
    std::vector<Avec> xvals_tecm_r6(8);
    std::vector<Avec> yvals_tecm_r6(8);
    
    LASGlobalDataLoop loop(LASGlobalDataLoop::TEC);
    do{
      if(loop.GetEntry(mask) == 0) continue;
      switch(loop.get_det()){
      case 0:
	switch(loop.get_ring()){
	case 0:
	  xvals_tecp_r4[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
	  yvals_tecp_r4[loop.get_beam()].push_back( loop.GetEntry(data) );
	  break;
	case 1:
	  xvals_tecp_r6[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
	  yvals_tecp_r6[loop.get_beam()].push_back( loop.GetEntry(data) );
	  break;
	}
	break;
      case 1:
	switch(loop.get_ring()){
	case 0:
	  xvals_tecm_r4[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
	  yvals_tecm_r4[loop.get_beam()].push_back( loop.GetEntry(data) );
	  break;
	case 1:
	  xvals_tecm_r6[loop.get_beam()].push_back( loop.GetEntry(zpos()) );
	  yvals_tecm_r6[loop.get_beam()].push_back( loop.GetEntry(data) );
	  break;
	}
	break;
      }
      
    }while (loop.next());

    TCanvas* cv_tec = new TCanvas("TEC", "TEC Internal");
    cv_tec->Divide(2,2);
    TMultiGraph* mgr; 
    cv_tec->cd(1);
    mgr = avec_draw(xvals_tecp_r4, yvals_tecp_r4, "TEC+ Ring 4","zpos [mm]","","AP");
    cv_tec->cd(2);
    mgr = avec_draw(xvals_tecp_r6, yvals_tecp_r6, "TEC+ Ring 6","zpos [mm]","","AP");
    cv_tec->cd(3);
    mgr = avec_draw(xvals_tecm_r4, yvals_tecm_r4, "TEC- Ring 4","zpos [mm]","","AP");
    cv_tec->cd(4);
    mgr = avec_draw(xvals_tecm_r6, yvals_tecm_r6, "TEC- Ring 6","zpos [mm]","","AP");
    AddLegend(mgr, legend, 2);
    retval = mgr;
  }

  return retval;
}


void global_data_save_canvas(const std::string& filename, const std::string& prefix)
{
  TFile f(filename.c_str(), "UPDATE");
  if(!f.IsOpen()) throw LAS::Exception("Error in global_data_save_canvas, could not open file " + filename);

  global_data_save_canvas(prefix);
  return;
}

void global_data_save_canvas(TFile& f, const std::string& prefix)
{
  f.cd();
  global_data_save_canvas(prefix);
  return;
}

void global_data_save_canvas(const std::string& prefix)
{
  TCanvas* cv = 0;

  cv = (TCanvas*)(gROOT->FindObject("TECplus_R4"));
  if(cv) cv->Write((std::string(prefix + "_tecp_r4")).c_str());
  cv = (TCanvas*)(gROOT->FindObject("TECplus_R6"));
  if(cv) cv->Write((std::string(prefix + "_tecp_r6")).c_str());
  cv = (TCanvas*)(gROOT->FindObject("TECminus_R4"));
  if(cv) cv->Write((std::string(prefix + "_tecm_r4")).c_str());
  cv = (TCanvas*)(gROOT->FindObject("TECminus_R6"));
  if(cv) cv->Write((std::string(prefix + "_tecm_r6")).c_str());
  cv = (TCanvas*)(gROOT->FindObject("TIB"));
  if(cv) cv->Write((std::string(prefix + "_tib")).c_str());
  cv = (TCanvas*)(gROOT->FindObject("TOB"));
  if(cv) cv->Write((std::string(prefix + "_tob")).c_str());
  cv = (TCanvas*)(gROOT->FindObject("TECplus_AT"));
  if(cv) cv->Write((std::string(prefix + "_tecp_at")).c_str());
  cv = (TCanvas*)(gROOT->FindObject("TECminus_AT"));
  if(cv) cv->Write((std::string(prefix + "_tecm_at")).c_str());

  return;
}


//! Plot the doubles that are contained in the LASGlobalData object with the corresponding errors using only non-masked entries
void draw_global_data(const LASGlobalData<double>& theData, const LASGlobalData<double>& theErrors, const LASGlobalData<int>& mask, LAS::beam_group group)
{
  group = LAS::ALL;


  TCanvas* cv_tec = new TCanvas("TEC", "TEC internal");
  cv_tec->Divide(2,2);

  std::vector<Avec> list_tecp_r4;
  std::vector<Avec> list_tecp_r6;
  std::vector<Avec> list_tecm_r4;
  std::vector<Avec> list_tecm_r6;
  
  std::vector<Avec> error_list_tecp_r4;
  std::vector<Avec> error_list_tecp_r6;
  std::vector<Avec> error_list_tecm_r4;
  std::vector<Avec> error_list_tecm_r6;
  
  std::vector<Avec> xlist_tecp_r4;
  std::vector<Avec> xlist_tecp_r6;
  std::vector<Avec> xlist_tecm_r4;
  std::vector<Avec> xlist_tecm_r6;

  for(int beam = 0; beam < 8; beam++){
    Avec data_tecp_r4;
    Avec data_tecp_r6;
    Avec data_tecm_r4;
    Avec data_tecm_r6;

    Avec errors_tecp_r4;
    Avec errors_tecp_r6;
    Avec errors_tecm_r4;
    Avec errors_tecm_r6;

    Avec xval_tecp_r4;
    Avec xval_tecp_r6;
    Avec xval_tecm_r4;
    Avec xval_tecm_r6;

    for(int disc = 0; disc < 9; disc++){

      if(mask.GetEntry(0, 0, beam, disc)){
	data_tecp_r4.push_back(theData.GetEntry(0, 0, beam, disc));
	errors_tecp_r4.push_back(theErrors.GetEntry(0, 0, beam, disc));
	xval_tecp_r4.push_back(disc +1);
      }
      if(mask.GetEntry(0, 1, beam, disc)){
	data_tecp_r6.push_back(theData.GetEntry(0, 1, beam, disc));
	errors_tecp_r6.push_back(theErrors.GetEntry(0, 1, beam, disc));
	xval_tecp_r6.push_back(disc +1);
      }
      if(mask.GetEntry(1, 0, beam, disc)){
	data_tecm_r4.push_back(theData.GetEntry(1, 0, beam, disc));
	errors_tecm_r4.push_back(theErrors.GetEntry(1, 0, beam, disc));
	xval_tecm_r4.push_back(disc +1);
      }
      if(mask.GetEntry(1, 1, beam, disc)){
	data_tecm_r6.push_back(theData.GetEntry(1, 1, beam, disc));
	errors_tecm_r6.push_back(theErrors.GetEntry(1, 1, beam, disc));
	xval_tecm_r6.push_back(disc +1);
      }

    }

    list_tecp_r4.push_back(data_tecp_r4);
    list_tecp_r6.push_back(data_tecp_r6);
    list_tecm_r4.push_back(data_tecm_r4);
    list_tecm_r6.push_back(data_tecm_r6);
    
    error_list_tecp_r4.push_back(errors_tecp_r4);
    error_list_tecp_r6.push_back(errors_tecp_r6);
    error_list_tecm_r4.push_back(errors_tecm_r4);
    error_list_tecm_r6.push_back(errors_tecm_r6);
    
    xlist_tecp_r4.push_back(xval_tecp_r4);
    xlist_tecp_r6.push_back(xval_tecp_r6);
    xlist_tecm_r4.push_back(xval_tecm_r4);
    xlist_tecm_r6.push_back(xval_tecm_r6);

    std::ostringstream title;
    title << "Beam " << beam;

  }

  cv_tec->cd(1);
  avec_draw(xlist_tecp_r4, list_tecp_r4, error_list_tecp_r4, "TEC plus Ring 4","Disc","","AP");

  cv_tec->cd(2);
  avec_draw(xlist_tecp_r6, list_tecp_r6, error_list_tecp_r6, "TEC plus Ring 6","Disc","","AP");

  cv_tec->cd(3);
  avec_draw(xlist_tecm_r4, list_tecm_r4, error_list_tecm_r4, "TEC minus Ring 4","Disc","","AP");

  cv_tec->cd(4);
  avec_draw(xlist_tecm_r6, list_tecm_r6, error_list_tecm_r6, "TEC minus Ring 6","Disc","","AP");


  TCanvas* cv_AT = new TCanvas("AT", "Alignment Tubes");
  cv_AT->Divide(2,2);


  std::vector<Avec> list_tib;
  std::vector<Avec> list_tob;
  std::vector<Avec> list_tecp_at;
  std::vector<Avec> list_tecm_at;
  
  std::vector<Avec> error_list_tib;
  std::vector<Avec> error_list_tob;
  std::vector<Avec> error_list_tecp_at;
  std::vector<Avec> error_list_tecm_at;
  
  std::vector<Avec> xlist_tib;
  std::vector<Avec> xlist_tob;
  std::vector<Avec> xlist_tecp_at;
  std::vector<Avec> xlist_tecm_at;

  for(int beam = 0; beam < 8; beam++){

    Avec data_tib;
    Avec data_tob;
    Avec data_tecp_at;
    Avec data_tecm_at;

    Avec errors_tib;
    Avec errors_tob;
    Avec errors_tecp_at;
    Avec errors_tecm_at;

    Avec xval_tib;
    Avec xval_tob;
    Avec xval_tecp_at;
    Avec xval_tecm_at;

    for(int zpos = 0; zpos < 6; zpos++){

      if(mask.GetEntry(2, -1, beam, zpos)){
	data_tib.push_back(theData.GetEntry(2, -1, beam, zpos));
	errors_tib.push_back(theErrors.GetEntry(2, -1, beam, zpos));
	xval_tib.push_back(zpos +1);
      }

      if(mask.GetEntry(3, -1, beam, zpos)){
	data_tob.push_back(theData.GetEntry(3, -1, beam, zpos));
	errors_tob.push_back(theErrors.GetEntry(3, -1, beam, zpos));
	xval_tob.push_back(zpos +1);
      }
    }

    for(int zpos = 0; zpos < 5; zpos++){
      if(mask.GetEntry(0, -1, beam, zpos)){
	data_tecp_at.push_back(theData.GetEntry(0, -1, beam, zpos));
	errors_tecp_at.push_back(theErrors.GetEntry(0, -1, beam, zpos));
	xval_tecp_at.push_back(zpos +1);
      }
      if(mask.GetEntry(1, -1, beam, zpos)){
	data_tecm_at.push_back(theData.GetEntry(1, -1, beam, zpos));
	errors_tecm_at.push_back(theErrors.GetEntry(1, -1, beam, zpos));
	xval_tecm_at.push_back(zpos +1);
      }
    }


    list_tib.push_back(data_tib);
    list_tob.push_back(data_tob);
    list_tecp_at.push_back(data_tecp_at);
    list_tecm_at.push_back(data_tecm_at);
    
    error_list_tib.push_back(errors_tib);
    error_list_tob.push_back(errors_tob);
    error_list_tecp_at.push_back(errors_tecp_at);
    error_list_tecm_at.push_back(errors_tecm_at);
    
    xlist_tib.push_back(xval_tib);
    xlist_tob.push_back(xval_tob);
    xlist_tecp_at.push_back(xval_tecp_at);
    xlist_tecm_at.push_back(xval_tecm_at);

    std::ostringstream title;
    title << "Beam " << beam;

  }

  cv_AT->cd(1);
  avec_draw(xlist_tib, list_tib, error_list_tib, "TIB","zpos","","AP");
  cv_AT->cd(2);
  avec_draw(xlist_tob, list_tob, error_list_tob, "TOB","zpos","","AP");
  cv_AT->cd(3);
  avec_draw(xlist_tecp_at, list_tecp_at, error_list_tecp_at, "TEC+ AT","Disc","","AP");
  cv_AT->cd(4);
  avec_draw(xlist_tecm_at, list_tecm_at, error_list_tecm_at, "TEC- AT","Disc","","AP");

}


//! Plot the doubles that are contained in the LASGlobalData object using only non-masked entries
void draw_global_data(const LASGlobalData<double>& theData, const LASGlobalData<int>& mask, LAS::beam_group group)
{
  group = LAS::ALL;

  TCanvas* cv_tec = new TCanvas("TEC", "TEC internal");
  cv_tec->Divide(2,2);

  std::vector<Avec> list_tecp_r4;
  std::vector<Avec> list_tecp_r6;
  std::vector<Avec> list_tecm_r4;
  std::vector<Avec> list_tecm_r6;
  
  std::vector<Avec> xlist_tecp_r4;
  std::vector<Avec> xlist_tecp_r6;
  std::vector<Avec> xlist_tecm_r4;
  std::vector<Avec> xlist_tecm_r6;

  for(int beam = 0; beam < 8; beam++){
    Avec data_tecp_r4;
    Avec data_tecp_r6;
    Avec data_tecm_r4;
    Avec data_tecm_r6;

    Avec xval_tecp_r4;
    Avec xval_tecp_r6;
    Avec xval_tecm_r4;
    Avec xval_tecm_r6;

    for(int disc = 0; disc < 9; disc++){

      if(mask.GetEntry(0, 0, beam, disc)){
	data_tecp_r4.push_back(theData.GetEntry(0, 0, beam, disc));
	xval_tecp_r4.push_back(disc +1);
      }
      if(mask.GetEntry(0, 1, beam, disc)){
	data_tecp_r6.push_back(theData.GetEntry(0, 1, beam, disc));
	xval_tecp_r6.push_back(disc +1);
      }
      if(mask.GetEntry(1, 0, beam, disc)){
	data_tecm_r4.push_back(theData.GetEntry(1, 0, beam, disc));
	xval_tecm_r4.push_back(disc +1);
      }
      if(mask.GetEntry(1, 1, beam, disc)){
	data_tecm_r6.push_back(theData.GetEntry(1, 1, beam, disc));
	xval_tecm_r6.push_back(disc +1);
      }

    }

    list_tecp_r4.push_back(data_tecp_r4);
    list_tecp_r6.push_back(data_tecp_r6);
    list_tecm_r4.push_back(data_tecm_r4);
    list_tecm_r6.push_back(data_tecm_r6);
    
    xlist_tecp_r4.push_back(xval_tecp_r4);
    xlist_tecp_r6.push_back(xval_tecp_r6);
    xlist_tecm_r4.push_back(xval_tecm_r4);
    xlist_tecm_r6.push_back(xval_tecm_r6);

    std::ostringstream title;
    title << "Beam " << beam;

  }

  cv_tec->cd(1);
  avec_draw(xlist_tecp_r4, list_tecp_r4, "TEC plus Ring 4","Disc","","AP");

  cv_tec->cd(2);
  avec_draw(xlist_tecp_r6, list_tecp_r6, "TEC plus Ring 6","Disc","","AP");

  cv_tec->cd(3);
  avec_draw(xlist_tecm_r4, list_tecm_r4, "TEC minus Ring 4","Disc","","AP");

  cv_tec->cd(4);
  avec_draw(xlist_tecm_r6, list_tecm_r6, "TEC minus Ring 6","Disc","","AP");


  TCanvas* cv_AT = new TCanvas("AT", "Alignment Tubes");
  cv_AT->Divide(2,2);


  std::vector<Avec> list_tib;
  std::vector<Avec> list_tob;
  std::vector<Avec> list_tecp_at;
  std::vector<Avec> list_tecm_at;
  
  std::vector<Avec> xlist_tib;
  std::vector<Avec> xlist_tob;
  std::vector<Avec> xlist_tecp_at;
  std::vector<Avec> xlist_tecm_at;

  for(int beam = 0; beam < 8; beam++){

    Avec data_tib;
    Avec data_tob;
    Avec data_tecp_at;
    Avec data_tecm_at;

    Avec xval_tib;
    Avec xval_tob;
    Avec xval_tecp_at;
    Avec xval_tecm_at;

    for(int zpos = 0; zpos < 6; zpos++){

      if(mask.GetEntry(2, -1, beam, zpos)){
	data_tib.push_back(theData.GetEntry(2, -1, beam, zpos));
	xval_tib.push_back(zpos +1);
      }

      if(mask.GetEntry(3, -1, beam, zpos)){
	data_tob.push_back(theData.GetEntry(3, -1, beam, zpos));
	xval_tob.push_back(zpos +1);
      }
    }

    for(int zpos = 0; zpos < 5; zpos++){
      if(mask.GetEntry(0, -1, beam, zpos)){
	data_tecp_at.push_back(theData.GetEntry(0, -1, beam, zpos));
	xval_tecp_at.push_back(zpos +1);
      }
      if(mask.GetEntry(1, -1, beam, zpos)){
	data_tecm_at.push_back(theData.GetEntry(1, -1, beam, zpos));
	xval_tecm_at.push_back(zpos +1);
      }
    }


    list_tib.push_back(data_tib);
    list_tob.push_back(data_tob);
    list_tecp_at.push_back(data_tecp_at);
    list_tecm_at.push_back(data_tecm_at);
    
    xlist_tib.push_back(xval_tib);
    xlist_tob.push_back(xval_tob);
    xlist_tecp_at.push_back(xval_tecp_at);
    xlist_tecm_at.push_back(xval_tecm_at);

    std::ostringstream title;
    title << "Beam " << beam;

  }

  cv_AT->cd(1);
  avec_draw(xlist_tib, list_tib, "TIB","zpos","","AP");
  cv_AT->cd(2);
  avec_draw(xlist_tob, list_tob, "TOB","zpos","","AP");
  cv_AT->cd(3);
  avec_draw(xlist_tecp_at, list_tecp_at, "TEC+ AT","Disc","","AP");
  cv_AT->cd(4);
  avec_draw(xlist_tecm_at, list_tecm_at, "TEC- AT","Disc","","AP");

}


//! Plot the doubles that are contained in the LASGlobalData object
void draw_global_data(const LASGlobalData<double>& theData, LAS::beam_group group)
{
  group = LAS::ALL;

  TCanvas* cv_tec = new TCanvas("TEC", "TEC");
  cv_tec->Divide(2,2);

  std::vector<Avec> dlist_tecp_r4;
  std::vector<Avec> dlist_tecp_r6;
  std::vector<Avec> dlist_tecm_r4;
  std::vector<Avec> dlist_tecm_r6;

  std::vector<Avec> xlist_tecp_r4;
  std::vector<Avec> xlist_tecp_r6;
  std::vector<Avec> xlist_tecm_r4;
  std::vector<Avec> xlist_tecm_r6;

  for(int beam = 0; beam < 8; beam++){
    Avec data_tecp_r4;
    Avec data_tecp_r6;
    Avec data_tecm_r4;
    Avec data_tecm_r6;

    Avec xval_tecp_r4(9,1,9);
    Avec xval_tecp_r6(9,1,9);
    Avec xval_tecm_r4(9,1,9);
    Avec xval_tecm_r6(9,1,9);

    for(int disc = 0; disc < 9; disc++){

      data_tecp_r4.push_back(theData.GetEntry(0, 0, beam, disc));
      data_tecp_r6.push_back(theData.GetEntry(0, 1, beam, disc));
      data_tecm_r4.push_back(theData.GetEntry(1, 0, beam, disc));
      data_tecm_r6.push_back(theData.GetEntry(1, 1, beam, disc));

    }
    dlist_tecp_r4.push_back(data_tecp_r4);
    dlist_tecp_r6.push_back(data_tecp_r6);
    dlist_tecm_r4.push_back(data_tecm_r4);
    dlist_tecm_r6.push_back(data_tecm_r6);

    xlist_tecp_r4.push_back(xval_tecp_r4);
    xlist_tecp_r6.push_back(xval_tecp_r6);
    xlist_tecm_r4.push_back(xval_tecm_r4);
    xlist_tecm_r6.push_back(xval_tecm_r6);
  }


  cv_tec->cd(1);
  avec_draw(xlist_tecp_r4, dlist_tecp_r4,"TEC+ R4","Disc","","AP");
  cv_tec->cd(2);
  avec_draw(xlist_tecp_r6, dlist_tecp_r6,"TEC+ R6","Disc","","AP");
  cv_tec->cd(3);
  avec_draw(xlist_tecm_r4, dlist_tecm_r4,"TEC- R4","Disc","","AP");
  cv_tec->cd(4);
  avec_draw(xlist_tecm_r6, dlist_tecm_r6,"TEC- R6","Disc","","AP");


  TCanvas* cv_AT = new TCanvas("AT", "Alignment Tubes");
  cv_AT->Divide(2,2);


  std::vector<Avec> dlist_tib;
  std::vector<Avec> dlist_tob;
  std::vector<Avec> dlist_tecp_at;
  std::vector<Avec> dlist_tecm_at;

  std::vector<Avec> xlist_tib;
  std::vector<Avec> xlist_tob;
  std::vector<Avec> xlist_tecp_at;
  std::vector<Avec> xlist_tecm_at;

  for(int beam = 0; beam < 8; beam++){

    Avec data_tib;
    Avec data_tob;
    Avec data_tecp_at;
    Avec data_tecm_at;

    Avec xval_tib(6,1,6);
    Avec xval_tob(6,1,6);
    Avec xval_tecp_at(5,1,5);
    Avec xval_tecm_at(5,1,5);

    for(int zpos = 0; zpos < 6; zpos++){
      data_tib.push_back(theData.GetEntry(2, -1, beam, zpos));
      data_tob.push_back(theData.GetEntry(3, -1, beam, zpos));
    }

    for(int zpos = 0; zpos < 5; zpos++){
      data_tecp_at.push_back(theData.GetEntry(0, -1, beam, zpos));
      data_tecm_at.push_back(theData.GetEntry(1, -1, beam, zpos));
    }

    dlist_tib.push_back(data_tib);
    dlist_tob.push_back(data_tob);
    dlist_tecp_at.push_back(data_tecp_at);
    dlist_tecm_at.push_back(data_tecm_at);

    xlist_tib.push_back(xval_tib);
    xlist_tob.push_back(xval_tob);
    xlist_tecp_at.push_back(xval_tecp_at);
    xlist_tecm_at.push_back(xval_tecm_at);

  }

  cv_AT->cd(1);
  avec_draw(xlist_tib, dlist_tib, "TIB","zpos","","AP");
  cv_AT->cd(2);
  avec_draw(xlist_tob, dlist_tob, "TOB","zpos","","AP");
  cv_AT->cd(3);
  avec_draw(xlist_tecp_at, dlist_tecp_at, "TEC+ AT","Disc","","AP");
  cv_AT->cd(4);
  avec_draw(xlist_tecm_at, dlist_tecm_at, "TEC- AT","Disc","","AP");

}

//! Plot the Avecs that are contained in the LASGlobalData object, a total of 8 Canvasses is created
void draw_global_data(const LASGlobalData<Avec>& theData, const std::string& options, LAS::beam_group group)
{

  if(group == LAS::ALL || group == LAS::TEC || group == LAS::TEC_PLUS || group == LAS::TEC_MINUS || group == LAS::TEC_PLUS_R4 || group == LAS::TEC_PLUS_R6 || group == LAS::TEC_MINUS_R4 || group == LAS::TEC_MINUS_R6){

  TCanvas* cv_tecp_r6 = new TCanvas("TECplus_R6", "TEC+ Ring 6");
  cv_tecp_r6->Divide(3,3);
  TCanvas* cv_tecp_r4 = new TCanvas("TECplus_R4", "TEC+ Ring 4");
  cv_tecp_r4->Divide(3,3);
  TCanvas* cv_tecm_r6 = new TCanvas("TECminus_R6", "TEC- Ring 6");
  cv_tecm_r6->Divide(3,3);
  TCanvas* cv_tecm_r4 = new TCanvas("TECminus_R4", "TEC- Ring 4");
  cv_tecm_r4->Divide(3,3);
  TMultiGraph* mgr=0;


  for(int beam = 0; beam < 8; beam++){
    std::vector<Avec> data_tecp_r4;
    std::vector<Avec> data_tecp_r6;
    std::vector<Avec> data_tecm_r4;
    std::vector<Avec> data_tecm_r6;

    std::vector<Avec> xval_tecp_r4;
    std::vector<Avec> xval_tecp_r6;
    std::vector<Avec> xval_tecm_r4;
    std::vector<Avec> xval_tecm_r6;

    for(int disc = 0; disc < 9; disc++){

      Avec data;

      data = theData.GetEntry(0, 0, beam, disc);
      data_tecp_r4.push_back(data);
      xval_tecp_r4.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(0, 1, beam, disc);
      data_tecp_r6.push_back(data);
      xval_tecp_r6.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(1, 0, beam, disc);
      data_tecm_r4.push_back(data);
      xval_tecm_r4.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(1, 1, beam, disc);
      data_tecm_r6.push_back(data);
      xval_tecm_r6.push_back(Avec(data.size(), 0, data.size() -1));

    }

    std::ostringstream title;
    title << "Beam " << beam;

    if(group == LAS::ALL || group == LAS::TEC || group == LAS::TEC_PLUS || group == LAS::TEC_PLUS_R4){
      //std::cout << "Drawing TEC+ Ring 4" << std::endl;
      cv_tecp_r4->cd(beam + 1);
      avec_draw(xval_tecp_r4, data_tecp_r4, title.str(), "", "", options);
    }

    if(group == LAS::ALL || group == LAS::TEC || group == LAS::TEC_PLUS || group == LAS::TEC_PLUS_R6){
      //std::cout << "Drawing TEC+ R6 beam " << beam << " xval_tecp_r6.size(): " << xval_tecp_r6.size() << std::endl;
      cv_tecp_r6->cd(beam + 1);
      avec_draw(xval_tecp_r6, data_tecp_r6, title.str(), "", "", options, 1, title.str());
    }

    if(group == LAS::ALL || group == LAS::TEC || group == LAS::TEC_MINUS || group == LAS::TEC_MINUS_R4){
      cv_tecm_r4->cd(beam + 1);
      avec_draw(xval_tecm_r4, data_tecm_r4, title.str(), "", "", options);
    }

    if(group == LAS::ALL || group == LAS::TEC || group == LAS::TEC_MINUS || group == LAS::TEC_MINUS_R6){
      cv_tecm_r6->cd(beam + 1);
      mgr = avec_draw(xval_tecm_r6, data_tecm_r6, title.str(), "", "", options);
    }
  }
  std::vector<std::string> leg_tit;
  leg_tit.push_back("Disc 1");
  leg_tit.push_back("Disc 2");
  leg_tit.push_back("Disc 3");
  leg_tit.push_back("Disc 4");
  leg_tit.push_back("Disc 5");
  leg_tit.push_back("Disc 6");
  leg_tit.push_back("Disc 7");
  leg_tit.push_back("Disc 8");
  leg_tit.push_back("Disc 9");

  cv_tecp_r4->cd(9);
  AddLegend(mgr, leg_tit);
  cv_tecp_r6->cd(9);
  AddLegend(mgr, leg_tit);
  cv_tecm_r4->cd(9);
  AddLegend(mgr, leg_tit);
  cv_tecm_r6->cd(9);
  AddLegend(mgr, leg_tit);
  }


  if(group == LAS::ALL || group == LAS::AT || group == LAS::TIB || group == LAS::TOB){

    TCanvas* cv_tib = 0;
    if(group != LAS::TOB){
      cv_tib = new TCanvas("TIB", "TIB");
      cv_tib->Divide(3,3);
    }
    TCanvas* cv_tob = 0;
    if(group != LAS::TIB){
      cv_tob = new TCanvas("TOB", "TOB");
      cv_tob->Divide(3,3);
    }
    TMultiGraph* mgr = 0;

    for(int beam = 0; beam < 8; beam++){

      std::vector<Avec> data_tib;
      std::vector<Avec> data_tob;

      std::vector<Avec> xval_tib;
      std::vector<Avec> xval_tob;

      for(int zpos = 0; zpos < 6; zpos++){

	Avec data;

	if(group != LAS::TOB){
	  data = theData.GetEntry(2, -1, beam, zpos);
	  data_tib.push_back(data);
	  xval_tib.push_back(Avec(data.size(), 0, data.size() -1));
	}

	if(group != LAS::TIB){
	  data = theData.GetEntry(3, -1, beam, zpos);
	  data_tob.push_back(data);
	  xval_tob.push_back(Avec(data.size(), 0, data.size() -1));
	}
      }

      std::ostringstream title;
      title << "Beam " << beam;
      
      if(group != LAS::TOB){
	cv_tib->cd(beam + 1);
	avec_draw(xval_tib, data_tib, title.str(), "", "", options);
      }
      if(group != LAS::TIB){
	cv_tob->cd(beam + 1);
	mgr = avec_draw(xval_tob, data_tob, title.str(), "", "", options);
      }
    }
    std::vector<std::string> leg_tit;
    leg_tit.push_back("Pos 1");
    leg_tit.push_back("Pos 2");
    leg_tit.push_back("Pos 3");
    leg_tit.push_back("Pos 4");
    leg_tit.push_back("Pos 5");
    leg_tit.push_back("Pos 6");
    
    if(group != LAS::TOB){
      cv_tib->cd(9);
      AddLegend(mgr, leg_tit);
    }
    if(group != LAS::TIB){
      cv_tob->cd(9);
      AddLegend(mgr, leg_tit);
    }
  }
  
  
  if(group == LAS::ALL || group == LAS::AT || group == LAS::TEC_PLUS_AT || group == LAS::TEC_MINUS_AT){
  TCanvas* cv_tecp_at = new TCanvas("TECplus_AT", "TEC+ AT");
  cv_tecp_at->Divide(3,3);
  TCanvas* cv_tecm_at = new TCanvas("TECminus_AT", "TEC- AT");
  cv_tecm_at->Divide(3,3);
  TMultiGraph* mgr=0;

  for(int beam = 0; beam < 8; beam++){
    std::vector<Avec> data_tecp_at;
    std::vector<Avec> data_tecm_at;

    std::vector<Avec> xval_tecp_at;
    std::vector<Avec> xval_tecm_at;

    for(int zpos = 0; zpos < 5; zpos++){

      Avec data;

      data = theData.GetEntry(0, -1, beam, zpos);
      data_tecp_at.push_back(data);
      xval_tecp_at.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(1, -1, beam, zpos);
      data_tecm_at.push_back(data);
      xval_tecm_at.push_back(Avec(data.size(), 0, data.size() -1));

    }

    std::ostringstream title;
    title << "Beam " << beam;

    //std::cout << "Title: " << title.str() << std::endl;

    cv_tecp_at->cd(beam + 1);
    avec_draw(xval_tecp_at, data_tecp_at, title.str(), "", "", options);

    cv_tecm_at->cd(beam + 1);
    mgr = avec_draw(xval_tecm_at, data_tecm_at, title.str(), "", "", options);
  }
  std::vector<std::string> leg_tit;
  leg_tit.push_back("Disc 1");
  leg_tit.push_back("Disc 2");
  leg_tit.push_back("Disc 3");
  leg_tit.push_back("Disc 4");
  leg_tit.push_back("Disc 5");

  cv_tecp_at->cd(9);
  AddLegend(mgr, leg_tit);
  cv_tecm_at->cd(9);
  AddLegend(mgr, leg_tit);
  }

  return;
}

//! Plot the Avecs that are contained in the LASGlobalData object, a total of 8 Canvasses is created
void draw_global_data(const LASGlobalData<Avec>& theData, const Avec& xval, const std::string& options, const std::string& axis, LAS::beam_group group)
{
  draw_global_data(theData, LASGlobalData<Avec>(xval), options, axis, group);
}

//! Plot the Avecs that are contained in the LASGlobalData object, a total of 8 Canvasses is created
void draw_global_data(const LASGlobalData<Avec>& theData, const LASGlobalData<Avec>& theErrors, const Avec& xval, const std::string& options, const std::string& axis, LAS::beam_group group)
{
  draw_global_data(theData, theErrors, LASGlobalData<Avec>(xval), options, axis, group);
}

//! Plot the Avecs that are contained in the LASGlobalData object, a total of 8 Canvasses is created
void draw_global_data(const LASGlobalData<Avec>& theData, const LASGlobalData<Avec>& theErrors, const LASGlobalData<Avec>& xval, const std::string& options, const std::string& axis, LAS::beam_group group)
{
  group = LAS::ALL;

  TCanvas* cv_tecp_r6 = new TCanvas("TECplus_R6", "TEC+ Ring 6");
  cv_tecp_r6->Divide(3,3);
  TCanvas* cv_tecp_r4 = new TCanvas("TECplus_R4", "TEC+ Ring 4");
  cv_tecp_r4->Divide(3,3);
  TCanvas* cv_tecm_r6 = new TCanvas("TECminus_R6", "TEC- Ring 6");
  cv_tecm_r6->Divide(3,3);
  TCanvas* cv_tecm_r4 = new TCanvas("TECminus_R4", "TEC- Ring 4");
  cv_tecm_r4->Divide(3,3);

  TMultiGraph* mgr=0;
  
  for(int beam = 0; beam < 8; beam++){
    std::vector<Avec> data_tecp_r4;
    std::vector<Avec> data_tecp_r6;
    std::vector<Avec> data_tecm_r4;
    std::vector<Avec> data_tecm_r6;

    std::vector<Avec> err_tecp_r4;
    std::vector<Avec> err_tecp_r6;
    std::vector<Avec> err_tecm_r4;
    std::vector<Avec> err_tecm_r6;

    std::vector<Avec> xval_tecp_r4;
    std::vector<Avec> xval_tecp_r6;
    std::vector<Avec> xval_tecm_r4;
    std::vector<Avec> xval_tecm_r6;

    for(int disc = 0; disc < 9; disc++){

      data_tecp_r4.push_back(   theData.GetEntry(0, 0, beam, disc) );
       err_tecp_r4.push_back( theErrors.GetEntry(0, 0, beam, disc) );
      xval_tecp_r4.push_back(      xval.GetEntry(0, 0, beam, disc) );

      data_tecp_r6.push_back(   theData.GetEntry(0, 1, beam, disc) );
       err_tecp_r6.push_back( theErrors.GetEntry(0, 1, beam, disc) );
      xval_tecp_r6.push_back(      xval.GetEntry(0, 1, beam, disc) );

      data_tecm_r4.push_back(   theData.GetEntry(1, 0, beam, disc) );
       err_tecm_r4.push_back( theErrors.GetEntry(1, 0, beam, disc) );
      xval_tecm_r4.push_back(      xval.GetEntry(1, 0, beam, disc) );

      data_tecm_r6.push_back(   theData.GetEntry(1, 1, beam, disc) );
       err_tecm_r6.push_back( theErrors.GetEntry(1, 1, beam, disc) );
      xval_tecm_r6.push_back(      xval.GetEntry(1, 1, beam, disc) );
    }

    std::ostringstream title;
    title << "Beam " << beam;

    cv_tecp_r4->cd(beam + 1);
    mgr = avec_draw(xval_tecp_r4, data_tecp_r4, err_tecp_r4, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tecp_r6->cd(beam + 1);
    mgr = avec_draw(xval_tecp_r6, data_tecp_r6, err_tecp_r6, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tecm_r4->cd(beam + 1);
    mgr = avec_draw(xval_tecm_r4, data_tecm_r4, err_tecm_r4, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tecm_r6->cd(beam + 1);
    mgr = avec_draw(xval_tecm_r6, data_tecm_r6, err_tecm_r6, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }
  }

  std::vector<std::string> tec_legend;
  tec_legend.push_back("Disc 1");
  tec_legend.push_back("Disc 2");
  tec_legend.push_back("Disc 3");
  tec_legend.push_back("Disc 4");
  tec_legend.push_back("Disc 5");
  tec_legend.push_back("Disc 6");
  tec_legend.push_back("Disc 7");
  tec_legend.push_back("Disc 8");
  tec_legend.push_back("Disc 9");
  cv_tecp_r4->cd(9);
  AddLegend(mgr, tec_legend, 1, 0.1, 0, 0.5, 1);
  cv_tecp_r6->cd(9);
  AddLegend(mgr, tec_legend, 1, 0.1, 0, 0.5, 1);
  cv_tecm_r4->cd(9);
  AddLegend(mgr, tec_legend, 1, 0.1, 0, 0.5, 1);
  cv_tecm_r6->cd(9);
  AddLegend(mgr, tec_legend, 1, 0.1, 0, 0.5, 1);



  TCanvas* cv_tib = new TCanvas("TIB", "TIB");
  cv_tib->Divide(3,3);
  TCanvas* cv_tob = new TCanvas("TOB", "TOB");
  cv_tob->Divide(3,3);

  for(int beam = 0; beam < 8; beam++){

    std::vector<Avec> data_tib;
    std::vector<Avec> data_tob;

    std::vector<Avec> err_tib;
    std::vector<Avec> err_tob;

    std::vector<Avec> xval_tib;
    std::vector<Avec> xval_tob;

    for(int zpos = 0; zpos < 6; zpos++){

      data_tib.push_back(   theData.GetEntry(2, -1, beam, zpos));
       err_tib.push_back( theErrors.GetEntry(2, -1, beam, zpos));
      xval_tib.push_back(      xval.GetEntry(2, -1, beam, zpos));

      data_tob.push_back(   theData.GetEntry(3, -1, beam, zpos));
       err_tob.push_back( theErrors.GetEntry(3, -1, beam, zpos));
      xval_tob.push_back(      xval.GetEntry(3, -1, beam, zpos));
    }

    std::ostringstream title;
    title << "Beam " << beam;

    TMultiGraph* mgr;
    cv_tib->cd(beam + 1);
    mgr = avec_draw(xval_tib, data_tib, err_tib, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tob->cd(beam + 1);
    mgr = avec_draw(xval_tob, data_tob, err_tob, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }
  }
  std::vector<std::string> leg_tit;
  leg_tit.push_back("Pos 1");
  leg_tit.push_back("Pos 2");
  leg_tit.push_back("Pos 3");
  leg_tit.push_back("Pos 4");
  leg_tit.push_back("Pos 5");
  leg_tit.push_back("Pos 6");

  cv_tib->cd(9);
  AddLegend(mgr, leg_tit);
  cv_tob->cd(9);
  AddLegend(mgr, leg_tit);

  TCanvas* cv_tecp_at = new TCanvas("TECplus_AT", "TEC+ AT");
  cv_tecp_at->Divide(3,3);
  TCanvas* cv_tecm_at = new TCanvas("TECminus_AT", "TEC- AT");
  cv_tecm_at->Divide(3,3);

  for(int beam = 0; beam < 8; beam++){
    std::vector<Avec> data_tecp_at;
    std::vector<Avec> data_tecm_at;

    std::vector<Avec> err_tecp_at;
    std::vector<Avec> err_tecm_at;

    std::vector<Avec> xval_tecp_at;
    std::vector<Avec> xval_tecm_at;

    for(int zpos = 0; zpos < 5; zpos++){

      data_tecp_at.push_back(   theData.GetEntry(0, -1, beam, zpos));
       err_tecp_at.push_back( theErrors.GetEntry(0, -1, beam, zpos));
      xval_tecp_at.push_back(      xval.GetEntry(0, -1, beam, zpos));

      data_tecm_at.push_back(   theData.GetEntry(1, -1, beam, zpos));
       err_tecm_at.push_back( theErrors.GetEntry(1, -1, beam, zpos));
      xval_tecm_at.push_back(      xval.GetEntry(1, -1, beam, zpos));
    }

    std::ostringstream title;
    title << "Beam " << beam;

    //std::cout << "Title: " << title.str() << std::endl;

    TMultiGraph* mgr;
    cv_tecp_at->cd(beam + 1);
    mgr = avec_draw(xval_tecp_at, data_tecp_at, err_tecp_at, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tecm_at->cd(beam + 1);
    mgr = avec_draw(xval_tecm_at, data_tecm_at, err_tecm_at, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }
  }

  return;
}


//! Plot the Avecs that are contained in the LASGlobalData object, a total of 8 Canvasses is created
void draw_global_data(const LASGlobalData<Avec>& theData, const LASGlobalData<Avec>& xval, const std::string& options, const std::string& axis, LAS::beam_group group)
{
  group = LAS::ALL;

  TCanvas* cv_tecp_r6 = new TCanvas("TECplus_R6", "TEC+ Ring 6");
  cv_tecp_r6->Divide(3,3);
  TCanvas* cv_tecp_r4 = new TCanvas("TECplus_R4", "TEC+ Ring 4");
  cv_tecp_r4->Divide(3,3);
  TCanvas* cv_tecm_r6 = new TCanvas("TECminus_R6", "TEC- Ring 6");
  cv_tecm_r6->Divide(3,3);
  TCanvas* cv_tecm_r4 = new TCanvas("TECminus_R4", "TEC- Ring 4");
  cv_tecm_r4->Divide(3,3);

  TMultiGraph* mgr=0;
  
  for(int beam = 0; beam < 8; beam++){
    std::vector<Avec> data_tecp_r4;
    std::vector<Avec> data_tecp_r6;
    std::vector<Avec> data_tecm_r4;
    std::vector<Avec> data_tecm_r6;

    std::vector<Avec> xval_tecp_r4;
    std::vector<Avec> xval_tecp_r6;
    std::vector<Avec> xval_tecm_r4;
    std::vector<Avec> xval_tecm_r6;

    for(int disc = 0; disc < 9; disc++){

      Avec data;

      data = theData.GetEntry(0, 0, beam, disc);
      data_tecp_r4.push_back(data);
      xval_tecp_r4.push_back(xval.GetEntry(0, 0, beam, disc));
      //       xval_tecp_r4.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(0, 1, beam, disc);
      data_tecp_r6.push_back(data);
      xval_tecp_r6.push_back(xval.GetEntry(0, 1, beam, disc));
      //xval_tecp_r6.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(1, 0, beam, disc);
      data_tecm_r4.push_back(data);
      xval_tecm_r4.push_back(xval.GetEntry(1, 0, beam, disc));
      //xval_tecm_r4.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(1, 1, beam, disc);
      data_tecm_r6.push_back(data);
      xval_tecm_r6.push_back(xval.GetEntry(1, 1, beam, disc));
      //xval_tecm_r6.push_back(Avec(data.size(), 0, data.size() -1));

    }

    std::ostringstream title;
    title << "Beam " << beam;

    cv_tecp_r4->cd(beam + 1);
    mgr = avec_draw(xval_tecp_r4, data_tecp_r4, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tecp_r6->cd(beam + 1);
    mgr = avec_draw(xval_tecp_r6, data_tecp_r6, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tecm_r4->cd(beam + 1);
    mgr = avec_draw(xval_tecm_r4, data_tecm_r4, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tecm_r6->cd(beam + 1);
    mgr = avec_draw(xval_tecm_r6, data_tecm_r6, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }
  }

  std::vector<std::string> tec_legend;
  tec_legend.push_back("Disc 1");
  tec_legend.push_back("Disc 2");
  tec_legend.push_back("Disc 3");
  tec_legend.push_back("Disc 4");
  tec_legend.push_back("Disc 5");
  tec_legend.push_back("Disc 6");
  tec_legend.push_back("Disc 7");
  tec_legend.push_back("Disc 8");
  tec_legend.push_back("Disc 9");
  cv_tecp_r4->cd(9);
  AddLegend(mgr, tec_legend, 1, 0.1, 0, 0.5, 1);
  cv_tecp_r6->cd(9);
  AddLegend(mgr, tec_legend, 1, 0.1, 0, 0.5, 1);
  cv_tecm_r4->cd(9);
  AddLegend(mgr, tec_legend, 1, 0.1, 0, 0.5, 1);
  cv_tecm_r6->cd(9);
  AddLegend(mgr, tec_legend, 1, 0.1, 0, 0.5, 1);



  TCanvas* cv_tib = new TCanvas("TIB", "TIB");
  cv_tib->Divide(3,3);
  TCanvas* cv_tob = new TCanvas("TOB", "TOB");
  cv_tob->Divide(3,3);

  for(int beam = 0; beam < 8; beam++){

    std::vector<Avec> data_tib;
    std::vector<Avec> data_tob;

    std::vector<Avec> xval_tib;
    std::vector<Avec> xval_tob;

    for(int zpos = 0; zpos < 6; zpos++){

      Avec data;

      data = theData.GetEntry(2, -1, beam, zpos);
      data_tib.push_back(data);
      xval_tib.push_back(xval.GetEntry(2, -1, beam, zpos));
      //xval_tib.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(3, -1, beam, zpos);
      data_tob.push_back(data);
      xval_tob.push_back(xval.GetEntry(3, -1, beam, zpos));
      //xval_tob.push_back(Avec(data.size(), 0, data.size() -1));
    }


    std::ostringstream title;
    title << "Beam " << beam;

    TMultiGraph* mgr;
    cv_tib->cd(beam + 1);
    mgr = avec_draw(xval_tib, data_tib, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tob->cd(beam + 1);
    mgr = avec_draw(xval_tob, data_tob, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }
  }
  std::vector<std::string>  leg_tit;
  leg_tit.push_back("Pos 1");
  leg_tit.push_back("Pos 2");
  leg_tit.push_back("Pos 3");
  leg_tit.push_back("Pos 4");
  leg_tit.push_back("Pos 5");
  leg_tit.push_back("Pos 6");

  cv_tib->cd(9);
  AddLegend(mgr, leg_tit);
  cv_tob->cd(9);
  AddLegend(mgr, leg_tit);


  TCanvas* cv_tecp_at = new TCanvas("TECplus_AT", "TEC+ AT");
  cv_tecp_at->Divide(3,3);
  TCanvas* cv_tecm_at = new TCanvas("TECminus_AT", "TEC- AT");
  cv_tecm_at->Divide(3,3);

  for(int beam = 0; beam < 8; beam++){
    std::vector<Avec> data_tecp_at;
    std::vector<Avec> data_tecm_at;

    std::vector<Avec> xval_tecp_at;
    std::vector<Avec> xval_tecm_at;

    for(int zpos = 0; zpos < 5; zpos++){

      Avec data;

      data = theData.GetEntry(0, -1, beam, zpos);
      data_tecp_at.push_back(data);
      xval_tecp_at.push_back(xval.GetEntry(0, -1, beam, zpos));
      //xval_tecp_at.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(1, -1, beam, zpos);
      data_tecm_at.push_back(data);
      xval_tecm_at.push_back(xval.GetEntry(1, -1, beam, zpos));
      //xval_tecm_at.push_back(Avec(data.size(), 0, data.size() -1));

    }

    std::ostringstream title;
    title << "Beam " << beam;

    TMultiGraph* mgr;
    cv_tecp_at->cd(beam + 1);
    mgr = avec_draw(xval_tecp_at, data_tecp_at, title.str(),axis,"",options);
    if(axis == "time"){
      std::cout << "Time axis" << std::endl;
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }

    cv_tecm_at->cd(beam + 1);
    mgr = avec_draw(xval_tecm_at, data_tecm_at, title.str(),axis,"",options);
    if(axis == "time"){
      mgr->GetXaxis()->SetNdivisions(505, kTRUE);
      mgr->GetXaxis()->SetTimeDisplay(1);
      //mgr->GetXaxis()->SetTimeFormat("%d\/%m\/%y %H:%M");
      mgr->GetXaxis()->SetTimeFormat("%H:%M:%S"); 
    }
  }

  leg_tit.clear();
  leg_tit.push_back("Disc 1");
  leg_tit.push_back("Disc 2");
  leg_tit.push_back("Disc 3");
  leg_tit.push_back("Disc 4");
  leg_tit.push_back("Disc 5");

  cv_tecp_at->cd(9);
  AddLegend(mgr, leg_tit);
  cv_tecm_at->cd(9);
  AddLegend(mgr, leg_tit);

  return;
}

//! Plot the vectors contained in a LASGlobalData object
void draw_global_data(LASGlobalData<std::vector<float> >& theData, LAS::beam_group group)
{
  group = LAS::ALL;

  TMultiGraph* mgr = 0;

  TCanvas* cv_tecp_r6 = new TCanvas("TECplus_R6", "TEC+ Ring 6");
  cv_tecp_r6->Divide(3,3);
  TCanvas* cv_tecp_r4 = new TCanvas("TECplus_R4", "TEC+ Ring 4");
  cv_tecp_r4->Divide(3,3);
  TCanvas* cv_tecm_r6 = new TCanvas("TECminus_R6", "TEC- Ring 6");
  cv_tecm_r6->Divide(3,3);
  TCanvas* cv_tecm_r4 = new TCanvas("TECminus_R4", "TEC- Ring 4");
  cv_tecm_r4->Divide(3,3);

  for(int beam = 0; beam < 8; beam++){
    std::vector<Avec> data_tecp_r4;
    std::vector<Avec> data_tecp_r6;
    std::vector<Avec> data_tecm_r4;
    std::vector<Avec> data_tecm_r6;

    std::vector<Avec> xval_tecp_r4;
    std::vector<Avec> xval_tecp_r6;
    std::vector<Avec> xval_tecm_r4;
    std::vector<Avec> xval_tecm_r6;

    for(int disc = 0; disc < 9; disc++){

      Avec data;

      data = theData.GetEntry(0, 0, beam, disc);
      data_tecp_r4.push_back(data);
      xval_tecp_r4.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(0, 1, beam, disc);
      data_tecp_r6.push_back(data);
      xval_tecp_r6.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(1, 0, beam, disc);
      data_tecm_r4.push_back(data);
      xval_tecm_r4.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(1, 1, beam, disc);
      data_tecm_r6.push_back(data);
      xval_tecm_r6.push_back(Avec(data.size(), 0, data.size() -1));

    }
    cv_tecp_r4->cd(beam + 1);
    mgr = avec_draw(xval_tecp_r4, data_tecp_r4);
    SetLineWidth(mgr, 3);

    cv_tecp_r6->cd(beam + 1);
    mgr = avec_draw(xval_tecp_r6, data_tecp_r6);
    SetLineWidth(mgr, 3);

    cv_tecm_r4->cd(beam + 1);
    mgr = avec_draw(xval_tecm_r4, data_tecm_r4);
    SetLineWidth(mgr, 3);

    cv_tecm_r6->cd(beam + 1);
    mgr = avec_draw(xval_tecm_r6, data_tecm_r6);
    SetLineWidth(mgr, 3);
  }

  std::vector<std::string> leg_titles;
  leg_titles.push_back("disc 1");
  leg_titles.push_back("disc 2");
  leg_titles.push_back("disc 3");
  leg_titles.push_back("disc 4");
  leg_titles.push_back("disc 5");
  leg_titles.push_back("disc 6");
  leg_titles.push_back("disc 7");
  leg_titles.push_back("disc 8");
  leg_titles.push_back("disc 9");

  cv_tecp_r4->cd(9);
  TLegend* leg = AddLegend(mgr, leg_titles);
  cv_tecp_r6->cd(9);
  leg->Draw("ALP");
  cv_tecm_r4->cd(9);
  leg->Draw("ALP");
  cv_tecm_r6->cd(9);
  leg->Draw("ALP");

  TCanvas* cv_tecp_at = new TCanvas("TECplus_AT", "TEC+ AT");
  cv_tecp_at->Divide(3,3);
  TCanvas* cv_tecm_at = new TCanvas("TECminus_AT", "TEC- AT");
  cv_tecm_at->Divide(3,3);

  for(int beam = 0; beam < 8; beam++){
    std::vector<Avec> data_tecp_at;
    std::vector<Avec> data_tecm_at;

    std::vector<Avec> xval_tecp_at;
    std::vector<Avec> xval_tecm_at;

    for(int zpos = 0; zpos < 5; zpos++){

      Avec data;

      data = theData.GetEntry(0, -1, beam, zpos);
      data_tecp_at.push_back(data);
      xval_tecp_at.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(1, -1, beam, zpos);
      data_tecm_at.push_back(data);
      xval_tecm_at.push_back(Avec(data.size(), 0, data.size() -1));

    }

    std::ostringstream title;
    title << "Beam " << beam + 1;

    cv_tecp_at->cd(beam + 1);
    mgr = avec_draw(xval_tecp_at, data_tecp_at, title.str());
    SetLineWidth(mgr, 3);

    cv_tecm_at->cd(beam + 1);
    mgr = avec_draw(xval_tecm_at, data_tecm_at, title.str());
    SetLineWidth(mgr, 3);
  }

  leg_titles.clear();
  leg_titles.push_back("zpos 1");
  leg_titles.push_back("zpos 2");
  leg_titles.push_back("zpos 3");
  leg_titles.push_back("zpos 4");
  leg_titles.push_back("zpos 5");

  cv_tecm_at->cd(9);
  leg = AddLegend(mgr, leg_titles);
  cv_tecp_at->cd(9);
  leg->Draw("ALP");


  TCanvas* cv_tib = new TCanvas("TIB", "TIB");
  cv_tib->Divide(3,3);
  TCanvas* cv_tob = new TCanvas("TOB", "TOB");
  cv_tob->Divide(3,3);

  for(int beam = 0; beam < 8; beam++){
    std::vector<Avec> data_tib;
    std::vector<Avec> data_tob;

    std::vector<Avec> xval_tib;
    std::vector<Avec> xval_tob;

    for(int zpos = 0; zpos < 6; zpos++){

      Avec data;

      data = theData.GetEntry(2, -1, beam, zpos);
      data_tib.push_back(data);
      xval_tib.push_back(Avec(data.size(), 0, data.size() -1));

      data = theData.GetEntry(3, -1, beam, zpos);
      data_tob.push_back(data);
      xval_tob.push_back(Avec(data.size(), 0, data.size() -1));
    }

    cv_tib->cd(beam + 1);
    mgr = avec_draw(xval_tib, data_tib);
    SetLineWidth(mgr, 3);

    cv_tob->cd(beam + 1);
    mgr = avec_draw(xval_tob, data_tob);
    SetLineWidth(mgr, 3);
  }

  leg_titles.push_back("zpos 6");
  cv_tib->cd(9);
  leg = AddLegend(mgr, leg_titles);
  cv_tob->cd(9);
  leg->Draw("ALP");


  return;
}



//! Plot the doubles that are contained in the LASGlobalData object using only non-masked entries
void plot_global_data(const LASGlobalData<double>& theData, const LASGlobalData<int>& mask, LAS::beam_group group)
{
  //group = LAS::ALL;

  TCanvas* cv_tec = new TCanvas("TEC", "TEC internal");
  cv_tec->Divide(2,2);

  Avec tecp_r4;
  Avec tecp_r6;
  Avec tecm_r4;
  Avec tecm_r6;

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
  
  cv_tec->cd(1);
  avec_plot(tecp_r4, 30, "TEC plus Ring 4");

  cv_tec->cd(2);
  avec_plot(tecp_r6, 30, "TEC plus Ring 6");

  cv_tec->cd(3);
  avec_plot(tecm_r4, 30, "TEC minus Ring 4");

  cv_tec->cd(4);
  avec_plot(tecm_r6, 30, "TEC minus Ring 6");


  if(group == LAS::ALL || group == LAS::AT){
    TCanvas* cv_AT = new TCanvas("AT", "Alignment Tubes");
    cv_AT->Divide(2,2);
    
    Avec tib;
    Avec tob;
    Avec tecp_at;
    Avec tecm_at;
    
    LASGlobalDataLoop l5(LASGlobalDataLoop::TIB);
    do{
      if(l5.GetEntry<int>(mask) != 0) tib.push_back(l5.GetEntry<double>(theData));
    }while(l5.next());
    
    LASGlobalDataLoop l6(LASGlobalDataLoop::TOB);
    do{
      if(l6.GetEntry<int>(mask) != 0) tob.push_back(l6.GetEntry<double>(theData));
    }while(l6.next());
    
    LASGlobalDataLoop l7(LASGlobalDataLoop::TEC_PLUS_AT);
    do{
      if(l7.GetEntry<int>(mask) != 0) tecp_at.push_back(l7.GetEntry<double>(theData));
    }while(l7.next());
    
    LASGlobalDataLoop l8(LASGlobalDataLoop::TEC_MINUS_AT);
    do{
      if(l8.GetEntry<int>(mask) != 0) tecm_at.push_back(l8.GetEntry<double>(theData));
    }while(l8.next());
    
    cv_AT->cd(1);
    avec_plot(tib, 30, "TIB");
    cv_AT->cd(2);
    avec_plot(tob, 30, "TOB");
    cv_AT->cd(3);
    avec_plot(tecp_at, 30, "TEC+ AT");
    cv_AT->cd(4);
    avec_plot(tecm_at, 30, "TEC- AT");
  }
}

//! Fill an object with the ranges on where to look for the signal
void fill_strip_ranges(LASGlobalData<int>& range_low, LASGlobalData<int>& range_high)
{
  // Define a strip range in which to look for the maximum

  const int tob_left_low  =  90;
  const int tob_left_high = 189;

  const int tob_right_low  = 330;
  const int tob_right_high = 429;

  // Default is the 100 central strips
  range_low  = LASGlobalData<int>(206); 
  range_high = LASGlobalData<int>(315); 


  // Some of the TOB modules get the beam far from the center
  // Beam 2
  range_low.GetTIBTOBEntry( 3,2,0) = tob_right_low;
  range_high.GetTIBTOBEntry(3,2,0) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,2,1) = tob_left_low;
  range_high.GetTIBTOBEntry(3,2,1) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,2,2) = tob_left_low;
  range_high.GetTIBTOBEntry(3,2,2) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,2,3) = tob_right_low;
  range_high.GetTIBTOBEntry(3,2,3) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,2,4) = tob_right_low;
  range_high.GetTIBTOBEntry(3,2,4) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,2,5) = tob_left_low;
  range_high.GetTIBTOBEntry(3,2,5) = tob_left_high;

  // Beam 3
  range_low.GetTIBTOBEntry( 3,3,0) = tob_right_low;
  range_high.GetTIBTOBEntry(3,3,0) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,3,1) = tob_left_low;
  range_high.GetTIBTOBEntry(3,3,1) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,3,2) = tob_left_low;
  range_high.GetTIBTOBEntry(3,3,2) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,3,3) = tob_right_low;
  range_high.GetTIBTOBEntry(3,3,3) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,3,4) = tob_right_low;
  range_high.GetTIBTOBEntry(3,3,4) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,3,5) = tob_left_low;
  range_high.GetTIBTOBEntry(3,3,5) = tob_left_high;

  // Beam 4
  range_low.GetTIBTOBEntry( 3,4,0) = tob_right_low;
  range_high.GetTIBTOBEntry(3,4,0) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,4,1) = tob_left_low;
  range_high.GetTIBTOBEntry(3,4,1) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,4,2) = tob_left_low;
  range_high.GetTIBTOBEntry(3,4,2) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,4,3) = tob_right_low;
  range_high.GetTIBTOBEntry(3,4,3) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,4,4) = tob_right_low;
  range_high.GetTIBTOBEntry(3,4,4) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,4,5) = tob_left_low;
  range_high.GetTIBTOBEntry(3,4,5) = tob_left_high;

  // Beam 5
  range_low.GetTIBTOBEntry( 3,5,0) = tob_left_low;
  range_high.GetTIBTOBEntry(3,5,0) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,5,1) = tob_right_low;
  range_high.GetTIBTOBEntry(3,5,1) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,5,2) = tob_right_low;
  range_high.GetTIBTOBEntry(3,5,2) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,5,3) = tob_left_low;
  range_high.GetTIBTOBEntry(3,5,3) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,5,4) = tob_left_low;
  range_high.GetTIBTOBEntry(3,5,4) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,5,5) = tob_right_low;
  range_high.GetTIBTOBEntry(3,5,5) = tob_right_high;

  // Beam 6
  range_low.GetTIBTOBEntry( 3,6,0) = tob_left_low;
  range_high.GetTIBTOBEntry(3,6,0) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,6,1) = tob_right_low;
  range_high.GetTIBTOBEntry(3,6,1) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,6,2) = tob_right_low;
  range_high.GetTIBTOBEntry(3,6,2) = tob_right_high;

  range_low.GetTIBTOBEntry( 3,6,3) = tob_left_low;
  range_high.GetTIBTOBEntry(3,6,3) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,6,4) = tob_left_low;
  range_high.GetTIBTOBEntry(3,6,4) = tob_left_high;

  range_low.GetTIBTOBEntry( 3,6,5) = tob_right_low;
  range_high.GetTIBTOBEntry(3,6,5) = tob_right_high;

}


LASGlobalDataLoop::loop_type convert_to_loop_type(LAS::beam_group group)
{
  switch(group){
  case LAS::TEC_PLUS :
    return LASGlobalDataLoop::TEC_PLUS ;
  case LAS::TEC_MINUS :
    return LASGlobalDataLoop::TEC_MINUS ;
  case LAS::TEC :
    return LASGlobalDataLoop::TEC ;
  case LAS::AT :
    return LASGlobalDataLoop::AT ;
  case LAS::TIB :
    return LASGlobalDataLoop::TIB ;
  case LAS::TOB :
    return LASGlobalDataLoop::TOB ;
  case LAS::TEC_PLUS_AT :
    return LASGlobalDataLoop::TEC_PLUS_AT ;
  case LAS::TEC_MINUS_AT :
    return LASGlobalDataLoop::TEC_MINUS_AT ;
  case LAS::TEC_AT :
    return LASGlobalDataLoop::TEC_AT ;
  case LAS::TEC_PLUS_R4 :
    return LASGlobalDataLoop::TEC_PLUS_R4 ;
  case LAS::TEC_PLUS_R6 :
    return LASGlobalDataLoop::TEC_PLUS_R6 ;
  case LAS::TEC_MINUS_R4 :
    return LASGlobalDataLoop::TEC_MINUS_R4 ;
  case LAS::TEC_MINUS_R6 :
    return LASGlobalDataLoop::TEC_MINUS_R6 ;
  default:
    return LASGlobalDataLoop::ALL;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// LASGlobalData Objects with module properties:
// pitch
// zpos
// r0
// theta

const LASGlobalData<double>& pitch()
{
  static LASGlobalData<double> pitch;
  static bool not_initialized = true;

  if(not_initialized){
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_AT);       !loop.finished(); loop.next() ) loop.GetEntry(pitch) = LAS::pitch_at_tec;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TIB);          !loop.finished(); loop.next() ) loop.GetEntry(pitch) = LAS::pitch_at_tib;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TOB);          !loop.finished(); loop.next() ) loop.GetEntry(pitch) = LAS::pitch_at_tob;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_PLUS_R4) ; !loop.finished(); loop.next() ) loop.GetEntry(pitch) = LAS::pitch_tec_r4;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_PLUS_R6) ; !loop.finished(); loop.next() ) loop.GetEntry(pitch) = LAS::pitch_tec_r6;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_MINUS_R4); !loop.finished(); loop.next() ) loop.GetEntry(pitch) = LAS::pitch_tec_r4;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_MINUS_R6); !loop.finished(); loop.next() ) loop.GetEntry(pitch) = LAS::pitch_tec_r6;

    not_initialized = false;
  }

  return pitch;
}

const LASGlobalData<double>& zpos()
{
  static LASGlobalData<double> zpos;
  static bool not_initialized = true;

  if(not_initialized){
    Avec tecz = LAS::zdat_tec + LAS::z_disc0;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC)   ; !loop.finished(); loop.next() ) loop.GetEntry<double>(zpos) = tecz[loop.get_zpos()] * (1 - 2 * loop.get_det());
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TIB)   ; !loop.finished(); loop.next() ) loop.GetEntry<double>(zpos) = LAS::zdat_tib[loop.get_zpos()];
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TOB)   ; !loop.finished(); loop.next() ) loop.GetEntry<double>(zpos) = LAS::zdat_tob[loop.get_zpos()];
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_AT); !loop.finished(); loop.next() ) loop.GetEntry<double>(zpos) = tecz[loop.get_zpos()] * (1 - 2 * loop.get_det());
    not_initialized = false;
  }
  return zpos;
}

const LASGlobalData<double>& zpos_eff()
{
  static LASGlobalData<double> zpos;
  static bool not_initialized = true;

  if(not_initialized){
    Avec tecz = LAS::zdat_tec + LAS::z_disc0;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC)   ; !loop.finished(); loop.next() ) loop.GetEntry<double>(zpos) = tecz[loop.get_zpos()] * (1 - 2 * loop.get_det());
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TIB)   ; !loop.finished(); loop.next() ) loop.GetEntry<double>(zpos) = LAS::zdat_tib_eff[loop.get_zpos()];
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TOB)   ; !loop.finished(); loop.next() ) loop.GetEntry<double>(zpos) = LAS::zdat_tob_eff[loop.get_zpos()];
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_AT); !loop.finished(); loop.next() ) loop.GetEntry<double>(zpos) = tecz[loop.get_zpos()] * (1 - 2 * loop.get_det());
    not_initialized = false;
  }
  return zpos;
}

const LASGlobalData<double>& r0()
{
  static LASGlobalData<double> r0;
  static bool not_initialized = true;

  if(not_initialized){
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_AT);       !loop.finished(); loop.next() ) loop.GetEntry<double>(r0) = LAS::r_at_tec;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TIB);          !loop.finished(); loop.next() ) loop.GetEntry<double>(r0) = LAS::r_at_tib;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TOB);          !loop.finished(); loop.next() ) loop.GetEntry<double>(r0) = LAS::r_at_tob;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_PLUS_R4) ; !loop.finished(); loop.next() ) loop.GetEntry<double>(r0) = LAS::r_tec_r4;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_PLUS_R6) ; !loop.finished(); loop.next() ) loop.GetEntry<double>(r0) = LAS::r_tec_r6;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_MINUS_R4); !loop.finished(); loop.next() ) loop.GetEntry<double>(r0) = LAS::r_tec_r4;
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC_MINUS_R6); !loop.finished(); loop.next() ) loop.GetEntry<double>(r0) = LAS::r_tec_r6;

    not_initialized = false;
  }
  return r0;
}

const LASGlobalData<double>& theta()
{
  static LASGlobalData<double> theta;
  static bool not_initialized = true;

  if(not_initialized){
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::TEC); !loop.finished(); loop.next() ) loop.GetEntry<double>(theta) = LAS::theta_tec.at(loop.get_beam());
    for(LASGlobalDataLoop loop(LASGlobalDataLoop::AT ); !loop.finished(); loop.next() ) loop.GetEntry<double>(theta) = LAS::theta_at.at(loop.get_beam());

    not_initialized = false;
  }
  return theta;
}
