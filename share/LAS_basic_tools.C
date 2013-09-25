#include "LAS_basic_tools.h"
#include "LAS_vectorfloat_tools.h"

#include "LASGlobalDataLoop.h"
#include <iostream>
#include <sstream>
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TChain.h"
#include "Avec2D.h"

LAS::Exception::Exception(const std::string& message) : std::runtime_error(message)
{
}

LAS::Exception::~Exception() throw()
{
}

//! Reads a line from a stream, empty lines and lines beginning with '#' are skipped 
std::string get_next_line(std::istream& input)
{
  while(input){
    std::string the_line;
    std::getline(input, the_line);
    if(the_line != "" && the_line[0] != '#') return the_line;
  }
  return "";
}

void mask_known_bad_modules(LASGlobalData<int>& mask)
{
  // TOB
  mask.GetEntry(3, -1, 6, 0) = 0;
  mask.GetEntry(3, -1, 7, 5) = 0;
  mask.GetEntry(3, -1, 2, 2) = 0;
  mask.GetEntry(3, -1, 5, 3) = 0;
  // mask.GetEntry(3, -1, 7, 4) = 0;

  // TIB
  mask.GetEntry(2, -1, 2, 0) = 0;

  // TEC- AT
  mask.GetEntry(1, -1, 6, 2) = 0;
  //mask.GetEntry(1, -1, 2, 0) = 0;
  //mask.GetEntry(1, -1, 2, 1) = 0;
  mask.GetEntry(1, -1, 3, 0) = 0;
  mask.GetEntry(1, -1, 3, 1) = 0;
  mask.GetEntry(1, -1, 3, 2) = 0;
  mask.GetEntry(1, -1, 3, 3) = 0;
  mask.GetEntry(1, -1, 3, 4) = 0;

  // TEC+ AT
  mask.GetEntry(0, -1, 3, 0) = 0;
  mask.GetEntry(0, -1, 2, 0) = 0;
  mask.GetEntry(0, -1, 2, 1) = 0;
  mask.GetEntry(0, -1, 2, 2) = 0;

  // TEC+ R4
  mask.GetEntry(0, 0, 7, 0) = 0;


}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Function to create a time structure for rROOT format
time_t make_root_time(double year, double month, double day, double hour, double min, double sec)
{
  // Needed to make mktime behave properly. Best would be to use <chrono> or some boost library for time manipulations
  setenv("TZ", "", 1);
  tzset();


  tm date;
  date.tm_sec  = (int)sec;
  date.tm_min  = (int)min;
  date.tm_hour = (int)hour;
  date.tm_mday = (int)day;
  date.tm_mon = (int)month - 1;
  date.tm_year = (int)year;
  return mktime(&date);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Function to compute the difference between root time and unix time
time_t make_root_time_offset()
{
  // Root defines the standard time offset to 1. Jan. 1995 00:00:00
  // But my system gives 01:00:00 for the hour of the offset, which I think is related to the daylight saving shift of one hour
  // I am not sure if I am doing things right, so there may be 1 hour shift w.r.t. the real time...
  // Note that the valid month range is 0-11, while the valid days run from 1-31 ...

  // Needed to make mktime behave properly. Best would be to use <chrono> or some boost library for time manipulations
  setenv("TZ", "", 1);
  tzset();

  tm mydate;
  mydate.tm_sec  = 0;
  mydate.tm_min  = 0;
  mydate.tm_hour = 0;
  mydate.tm_mday = 1;
  mydate.tm_mon = 0;
  mydate.tm_year = 95;
  time_t root_offset = mktime(&mydate);
  //std::cout << "root_offset: " << root_offset << "  =  " << ctime(&root_offset) << std::endl;
  return root_offset;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Sort the Tree by ascending event number
void sort_tree(const std::string& input_filename, const std::string& treename, const std::string& output_filename)
{
  // Open the data file
  TFile f(input_filename.c_str(),"READ");
  if(!f.IsOpen()){
    std::cerr << "Could not open file " << input_filename << std::endl;
    return;
  }

  // Find the Tree
  TTree* tree = 0;
  f.GetObject(treename.c_str(),tree);
  if(!tree){
    std::cerr << "Did not find Tree " << treename << std::endl;
    return;
  }
  sort_tree(tree, output_filename);
}

//! Merge the trees in the files supplied by the list into one single tree that is stored in outfile
void TreeMerge(const std::vector<std::string>& file_list, const std::string& tree_name, const std::string& outfile)
{
  // First open output file
  std::cout << "Opening output file: " << outfile << std::endl;
  TFile output(outfile.c_str(), "RECREATE");
  if(!output.IsOpen()){
    std::cerr << "Could not open file " << outfile << std::endl;
    return;
  }

  TChain chain(tree_name.c_str());

  // Loop over all files in the list, look for the tree and add it to the chain
  std::vector<std::string>::const_iterator i;
  for(i = file_list.begin(); i != file_list.end(); i++){
    std::cout << "Processing " << *i << std::endl;

    // Open the data file, skip if not possible
    TFile f(i->c_str(),"READ");
    if(!f.IsOpen()){
      std::cerr << "Could not open file " << *i << std::endl;
      continue;
    }

    // Find the Tree, skip if not found
    TTree* tree = 0;
    //f.GetObject("makeroottree/AC1B",tree);
    f.GetObject(tree_name.c_str(),tree);
    if(!tree){
      std::cerr << "Did not find Tree" << std::endl;
      continue;
    }

    // Add to chain
    chain.Add(i->c_str());
  }

  std::cout << "Cloning Tree, please be patient..." << std::endl;
  output.cd();
  TTree* tree =  chain.CloneTree();
  tree->Write();

  std::cout << "Total Events: " << tree->GetEntries() << std::endl;
  output.Close();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Sort the Tree by ascending event number
void sort_tree(TTree* tree, const std::string& output_filename)
{

  if(!tree){
    std::cerr << "Error in sort_tree: supplied tree pointer is NULL" << std::endl;
    return;
  }

  Int_t nentries = (Int_t)tree->GetEntries();

  //Drawing variable eventnumber with no graphics option.
  //variable eventnumber stored in array fV1 (see TTree::Draw)
  tree->Draw("eventnumber","","goff");
  Int_t *index = new Int_t[nentries];
  //sort array containing eventnumber in decreasing order
  //The array index contains the entry numbers in decreasing order of eventnumber
  TMath::Sort(nentries,tree->GetV1(),index, kFALSE);
	
  //open new file to store the sorted Tree
  TFile f2(output_filename.c_str(), "recreate");
  //Create an empty clone of the original tree
  TTree *tsorted = (TTree*)tree->CloneTree(0);
  for (Int_t i=0;i<nentries;i++) {
    tree->GetEntry(index[i]);
    tsorted->Fill();
  }
  tsorted->Write();
  delete [] index;
}


//! Returns the Module name
std::string  GetModuleName(const LASGlobalDataLoop& loop, int underscore, int prefix)
{
  return GetModuleName(loop.get_det(), loop.get_ring(), loop.get_beam(), loop.get_zpos(), underscore, prefix);
}

//! Returns the Module name
std::string  GetModuleName( int subdetector, int ring, int beam, int zpos , int underscore, int prefix)
{
  if(invalid_indices(subdetector, beam, zpos, ring))return "";

  switch(subdetector){
  case 0:
  case 1:
    switch(ring){
    case -1:
      return GetTEC2TECName(subdetector, beam, zpos, underscore, prefix);
    case 0:
    case 1:
      return GetTECName(subdetector, ring, beam, zpos, underscore, prefix);
    }
  case 2:
  case 3:
    return GetTIBTOBName(subdetector, beam, zpos, underscore, prefix);
  }
  return "";
}

//! Returns the TEC module name
std::string  GetTECName( int subdetector, int ring, int beam, int disc , int underscore, int prefix)
{
  if(invalid_indices(subdetector, beam, disc, ring))return "";

  std::string sep=" ";
  if(underscore)sep="_";

  std::ostringstream name;
  if(prefix){
    switch(subdetector){
    case 0: 
      name << "TECP" << sep;
      break;
    case 1:
      name << "TECM" << sep;
      break;
    }
  }
  name << "Ring" << sep << ring*2+4 << sep << "Beam" << sep << beam+1 << sep << "Disc" << sep << disc+1;
  return name.str();
}

//! returns the TIB or TOB module name
std::string  GetTIBTOBName( int subdetector, int beam, int pos, int underscore , int prefix)
{
  if(invalid_indices(subdetector, beam, pos))return "";
  std::string sep=" ";
  if(underscore)sep="_";

  std::ostringstream name;
  if(prefix){
    switch(subdetector){
    case 2: 
      name << "TIB" << sep;
      break;
    case 3:
      name << "TOB" << sep;
      break;
    }
  }
  name << "Beam" << sep << beam+1 << sep << "Pos" << sep << pos+1;
  return name.str();
}

//! Returns  TEC Module name for alignment tubes
std::string  GetTEC2TECName( int subdetector, int beam, int disc, int underscore , int prefix)
{
  if(invalid_indices(subdetector, beam, disc))return "";
  std::string sep=" ";
  if(underscore)sep="_";

  std::ostringstream name;
  if(prefix){
    switch(subdetector){
    case 0: 
      name << "TECP" << sep;
      break;
    case 1:
      name << "TECM" << sep;
      break;
    }
  }
  name << "Beam" << sep << beam+1 << sep << "Disc" << sep << disc+1;
  return name.str();
}

/// Helping function to validate set of indices
//! Check if the set of indices correspond to an existing module
/*! type specifies the kind of check to do:\n
    type = 0 : any LAS module\n
    type = 1 : any LAS AT module\n
    type = 2 : any LAS TEC internal module\n
    Supplying only the first 3 arguments is equivalent to choosing type = 1\n
    return values:\n
    0 indices are a valid combination\n
    1 invalid subdet\n
    2 invalid beam\n
    3 invalid zpos\n
    4 invalid ring\n
    If verbose is nonzero, error messages are printed
*/
int invalid_indices(int subdet, int beam, int zpos, int ring, int type, int verbose)
{
  // Check if beam is valid
  if( beam < 0 || beam > 7){
    if(verbose) std::cerr << "Invalid indices subdet" << subdet << " ->beam " << beam << "<- zpos: " << zpos << " ring " << ring << std::endl;
    return 2;
  }

  switch(subdet){
  case 0: // TEC+
  case 1: // TEC-
    switch(ring){
    case -1: //AT
      if( type == 2 ){ // only TEC internal modules were asked
	if(verbose) std::cerr << "Invalid indices subdet" << subdet << " beam " << beam << " zpos: " << zpos << " ->ring " << ring << "<-" << std::endl;
	return 4;
      }
      if( zpos < 0 || zpos > 4 ){ // AT TEC has only 5 discs
	if(verbose) std::cerr << "Invalid indices subdet" << subdet << " beam " << beam << " ->zpos: " << zpos << "<- ring " << ring << std::endl;
	return 3;
      }
      break;
    case 0: // Ring 4
    case 1: // Ring 6
      if( type == 1 ){ // only AT modules were asked
	if(verbose) std::cerr << "Invalid indices subdet" << subdet << " beam " << beam << " zpos: " << zpos << " ->ring " << ring << "<-" << std::endl;
	return 4;
      }
      if( zpos < 0 || zpos > 8 ){ // TEC internal has 9 discs
	if(verbose) std::cerr << "Invalid indices subdet" << subdet << " beam " << beam << " ->zpos: " << zpos << "<- ring " << ring << std::endl;
	return 3;
      }
      break;
    default: // Non-valid ring number
      if(verbose) std::cerr << "Invalid indices subdet" << subdet << " beam " << beam << " zpos: " << zpos << " ->ring " << ring << "<-" << std::endl;
      return 4;
    }
    break;
  case 2: // TIB
  case 3: // TOB
    if( type == 2 ){ // only TEC internal modules were asked
      if(verbose) std::cerr << "Invalid indices ->subdet" << subdet << "<- beam " << beam << " zpos: " << zpos << " ring " << ring << std::endl;
      return 1;
    }
    if( ring != -1 ){ // AT modules have always ring -1
      if(verbose) std::cerr << "Invalid indices subdet" << subdet << " beam " << beam << " zpos: " << zpos << " ->ring " << ring << "<-" << std::endl;
      return 4;
    }
    if( zpos < 0 || zpos > 5 ){ // 6 z-positions for TIB and TOB
      if(verbose) std::cerr << "Invalid indices subdet" << subdet << " beam " << beam << " ->zpos: " << zpos << "<- ring " << ring << std::endl;
      return 3;
    }
    break;
  default: // Invalid subdet Nr.
    if(verbose) std::cerr << "Invalid indices ->subdet" << subdet << "<- beam " << beam << " zpos: " << zpos << " ring " << ring << std::endl;
    return 1;
  }

  return 0;
}


///////////////////////////////////////
// Functions for handling Root stuff //
///////////////////////////////////////

//! Create a Legend from a TMultiGraph and draw it if plot==true
/*!  If leg_tit has a differnt number of entries than there are Graphs in mgr, the smaller number of entries is relevant
 */
TLegend* AddLegend(TMultiGraph* mgr, std::vector<std::string>& leg_tit, int ncol, float xlo, float ylo, float xhi, float yhi, bool plot)
{
  TLegend * leg = new TLegend(xlo, ylo, xhi, yhi);
  leg->SetNColumns(ncol);
  TList* gr_list = mgr->GetListOfGraphs();
  TIter next(gr_list);
  TObject* obj;
  std::vector<std::string>::size_type i=0;
  while ((obj = next()) && (i < leg_tit.size()) ){
    //std::cout << "Add Entry" << std::endl;
    TLegendEntry* tle=leg->AddEntry(obj,leg_tit[i].c_str(),"P");
    tle->SetMarkerSize(2);
    i++;
  }
  if(plot)leg->Draw("ALP");
  return leg;
}

//! set the Line Width for all Graphs of the MultiGraph
void SetLineWidth(TMultiGraph* mgr, Width_t lwidth)
{
 TList* gr_list = mgr->GetListOfGraphs();
 TIter next(gr_list);
 TGraph* graph;
 while ((graph = (TGraph*)next())) graph->SetLineWidth(lwidth);
 return;
} 

//! set the X-range for all Graphs of the MultiGraph
void SetXRangeUser(TMultiGraph* mgr, Double_t ufirst, Double_t ulast)
{
 TList* gr_list = mgr->GetListOfGraphs();
 TIter next(gr_list);
 TGraph* graph;
 while ((graph = (TGraph*)next())) graph->GetXaxis()->SetRangeUser(ufirst, ulast);
 return;
} 

