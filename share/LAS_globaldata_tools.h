#ifndef __LAS_GLOBALDATA_TOOLS_H__
#define __LAS_GLOBALDATA_TOOLS_H__

class Avec;

#include "LAS_basic_tools.h"
#include "LASGlobalData.h"
#include "LASGlobalDataLoop.h"
#include "TFile.h"
#include "TCanvas.h"
#include <vector>


const LASGlobalData<double>& pitch();
const LASGlobalData<double>& zpos();
const LASGlobalData<double>& zpos_eff();
const LASGlobalData<double>& r0();
const LASGlobalData<double>& theta();


LASGlobalData<double> global_data_random_gauss();

void global_data_stat(const LASGlobalData<Avec>& data, LASGlobalData<double>& mean, LASGlobalData<double>& rms, LASGlobalData<int>& count);

//! Draw global data values against z-positions
//void draw_global_data_z(const LASGlobalData<double>& data, const LASGlobalData<double>& error, const LASGlobalData<int>& mask, LAS::beam_group group, const std::string& canvas_name);
TMultiGraph* draw_global_data_z(const LASGlobalData<double>& data, const LASGlobalData<double>& error, const LASGlobalData<int>& mask = LASGlobalData<int>(1), LAS::beam_group group = LAS::ALL, const std::string& canvas_name = "cv_global_data_z");
TMultiGraph* draw_global_data_z(const LASGlobalData<double>& data, const LASGlobalData<int>& mask = LASGlobalData<int>(1), LAS::beam_group group = LAS::ALL, const std::string& canvas_name = "cv_global_data_z");

void global_data_save_canvas(const std::string& filename, const std::string& prefix);
void global_data_save_canvas(TFile& f, const std::string& prefix);
void global_data_save_canvas(const std::string& prefix);

void draw_global_data(const LASGlobalData<double>& theData, LAS::beam_group group = LAS::ALL);
void draw_global_data(const LASGlobalData<double>& theData, const LASGlobalData<int>& mask, LAS::beam_group group = LAS::ALL);
void draw_global_data(const LASGlobalData<double>& theData, const LASGlobalData<double>& theErrors, const LASGlobalData<int>& mask, LAS::beam_group group = LAS::ALL);

//! Draw the Avecs that are contained in the LASGlobalData object, a total of 8 Canvasses is created
// Without errors
void draw_global_data(const LASGlobalData<Avec>& theData, const std::string& options="AP", LAS::beam_group group = LAS::ALL);
void draw_global_data(const LASGlobalData<Avec>& theData, const Avec& xval, const std::string& options="AP", const std::string& axis="", LAS::beam_group group = LAS::ALL);
void draw_global_data(const LASGlobalData<Avec>& theData, const LASGlobalData<Avec>& xval, const std::string& options="AP", const std::string& axis="", LAS::beam_group group = LAS::ALL);

// With errors
void draw_global_data(const LASGlobalData<Avec>& theData, const LASGlobalData<Avec>& theErrors, const Avec& xval, const std::string& options="AP", const std::string& axis="", LAS::beam_group group = LAS::ALL);
void draw_global_data(const LASGlobalData<Avec>& theData, const LASGlobalData<Avec>& theErrors, const LASGlobalData<Avec>& xval, const std::string& options="AP", const std::string& axis="", LAS::beam_group group = LAS::ALL);

//! Draw the std::vector<float> arrays
void draw_global_data(LASGlobalData<std::vector<float> >& theData, LAS::beam_group group = LAS::ALL);


void plot_global_data(const LASGlobalData<double>& theData, const LASGlobalData<int>& mask, LAS::beam_group group = LAS::ALL);

LASGlobalDataLoop::loop_type convert_to_loop_type(LAS::beam_group group);

//! Fill an object with the ranges on where to look for the signal
void fill_strip_ranges(LASGlobalData<int>& range_low, LASGlobalData<int>& range_high);

///////////////////////////////////
///          Templates          ///
///////////////////////////////////

template <class T>
//! Remove nan entries and replace them by 'replacement'
void remove_nan(LASGlobalData<T>& data, T replacement = 0)
{
  LASGlobalDataLoop lp;
  do{
    if(lp.GetEntry(data) != lp.GetEntry(data))lp.GetEntry(data) = replacement;
  }while ( lp.next() );
}

template <class T>
//! Retrieve a LASGlobalData Object from a root file through TFile reference
LASGlobalData<T>& global_data_get(const std::string& name, TFile& file)
{
  static LASGlobalData<T> empty;

  LASGlobalData<T>* obj_ptr = 0;
  file.GetObject(name.c_str(), obj_ptr);
  if(!obj_ptr){
    std::cerr << "Did not find " << name << std::endl;
    return empty;
  }
  return *obj_ptr;
}

template <class T>
//! Retrieve a LASGlobalData Object from a root file through file name
LASGlobalData<T>& global_data_get(const std::string& objname, const std::string& filename)
{
  static LASGlobalData<T> empty;
  TFile f(filename.c_str(),"READ");
  if(!f.IsOpen()){
    std::cerr << "Could not open file " << filename << std::endl;
    return empty;
  }
  return global_data_get<T>(objname, f);
}

/* template <class T> */
/* void print_global_data(const LASGlobalData<T>& theData, LAS::beam_group group = LAS::ALL, bool headers = true, std::ostream& out = std::cout); */

template <class T>
//! print the global data to a stream
// If 'headers' is true the beam groups are labelled
void print_global_data(const LASGlobalData<T>& theData, LAS::beam_group group = LAS::ALL, bool headers = true, std::ostream& out = std::cout)
{

  // TEC+ Ring 4
  if(group == LAS::ALL || group == LAS::TEC || group == LAS::TEC_PLUS || group == LAS::TEC_PLUS_R4){
    if(headers) out << "TEC+ Ring 4" << std::endl;
    for(int beam = -1; beam < 8; beam ++){
      if(headers){
	if(beam == -1) out << "Disc";
	else out << "Beam " << beam+1;
      }
      else if(beam == -1)continue;
      for(int disc = 0; disc < 9; disc++){
	if(headers && beam == -1) out << "\t" << disc;
	else out << "\t" << theData.GetEntry(0, 0, beam, disc);
      }
      out << std::endl;
    }
  }

  // TEC+ Ring 6
  if(group == LAS::ALL || group == LAS::TEC || group == LAS::TEC_PLUS || group == LAS::TEC_PLUS_R6){
    if(headers) out << "TEC+ Ring 6" << std::endl;
    for(int beam = -1; beam < 8; beam ++){
      if(headers){
	if(beam == -1) out << "Disc";
	else out << "Beam " << beam+1;
      } 
      else if(beam == -1)continue;
      for(int disc = 0; disc < 9; disc++){
	if(headers && beam == -1) out << "\t" << disc;
	else out << "\t" << theData.GetEntry(0, 1, beam, disc);
      }
      out << std::endl;
    }
  }


  // TEC- Ring 4
  if(group == LAS::ALL || group == LAS::TEC || group == LAS::TEC_MINUS || group == LAS::TEC_MINUS_R4){
    if(headers) out << "TEC- Ring 4" << std::endl;
    for(int beam = -1; beam < 8; beam ++){
      if(headers){
	if(beam == -1) out << "Disc";
	else out << "Beam " << beam+1;
      } 
      else if(beam == -1)continue;
      for(int disc = 0; disc < 9; disc++){
	if(headers && beam == -1) out << "\t" << disc;
	else out << "\t" << theData.GetEntry(1, 0, beam, disc);
      }
      out << std::endl;
    }
  }

  // TEC- Ring 6
  if(group == LAS::ALL || group == LAS::TEC || group == LAS::TEC_MINUS || group == LAS::TEC_MINUS_R6){
  if(headers) out << "TEC- Ring 6" << std::endl;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) out << "Disc";
      else out << "Beam " << beam+1;
    } 
    else if(beam == -1)continue;
    for(int disc = 0; disc < 9; disc++){
      if(headers && beam == -1) out << "\t" << disc;
      else out << "\t" << theData.GetEntry(1, 1, beam, disc);
    }
    out << std::endl;
  }
  }

  // TIB
  if(group == LAS::ALL || group == LAS::AT || group == LAS::TIB){
  if(headers) out << "TIB" << std::endl;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) out << "Zpos";
      else out << "Beam " << beam+1;
    } 
    else if(beam == -1)continue;
    for(int zpos = 0; zpos < 6; zpos++){
      if(headers && beam == -1) out << "\t" << zpos;
      else out << "\t" << theData.GetEntry(2, -1, beam, zpos);
    }
    out << std::endl;
  }
  }

  // TOB
  if(group == LAS::ALL || group == LAS::AT || group == LAS::TOB){
  if(headers) out << "TOB" << std::endl;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) out << "Zpos";
      else out << "Beam " << beam+1;
    } 
    else if(beam == -1)continue;
    for(int zpos = 0; zpos < 6; zpos++){
      if(headers && beam == -1) out << "\t" << zpos;
      else out << "\t" << theData.GetEntry(3, -1, beam, zpos);
    }
    out << std::endl;
  }
  }

  // TEC+ AT
  if(group == LAS::ALL || group == LAS::AT || group == LAS::TEC_AT || group == LAS::TEC_PLUS_AT){
  if(headers) out << "TEC+ AT" << std::endl;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) out << "Disc";
      else out << "Beam " << beam+1;
    } 
    else if(beam == -1)continue;
    for(int zpos = 0; zpos < 5; zpos++){
      if(headers && beam == -1) out << "\t" << zpos;
      else out << "\t" << theData.GetEntry(0, -1, beam, zpos);
    }
    out << std::endl;
  }
  }

  // TEC- AT
  if(group == LAS::ALL || group == LAS::AT || group == LAS::TEC_AT || group == LAS::TEC_MINUS_AT){
  if(headers) out << "TEC- AT" << std::endl;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) out << "Disc";
      else out << "Beam " << beam+1;
    } 
    else if(beam == -1)continue;
    for(int zpos = 0; zpos < 5; zpos++){
      if(headers && beam == -1) out << "\t" << zpos;
      else out << "\t" << theData.GetEntry(1, -1, beam, zpos);
    }
    out << std::endl;
  }
  }

}

template <class T>
//! Read the values from a stream and write them into the LASGlobalData object. The order is the same as print_global_data
void read_global_data(LASGlobalData<T>& theData, LAS::beam_group group = LAS::ALL, bool headers = true, std::istream& in = std::cin)
{
  group = LAS::ALL; 

  std::string dummy;
  T phony;

  // TEC+ Ring 4
  if(headers) in >> dummy;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) in >> dummy;
      else {in >> dummy; in >> phony;}
    }
    else if(beam == -1)continue;
    for(int disc = 0; disc < 9; disc++){
      if(headers && beam == -1) in >> phony;
      else in >> theData.GetEntry(0, 0, beam, disc);
    }
  }

  // TEC+ Ring 6
  if(headers) in >> dummy;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) in >> dummy;
      else {in >> dummy; in >> phony;}
    }
    else if(beam == -1)continue;
    for(int disc = 0; disc < 9; disc++){
      if(headers && beam == -1) in >> phony;
      else in >> theData.GetEntry(0, 1, beam, disc);
    }
  }

  // TEC- Ring 4
  if(headers) in >> dummy;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) in >> dummy;
      else {in >> dummy; in >> phony;}
    }
    else if(beam == -1)continue;
    for(int disc = 0; disc < 9; disc++){
      if(headers && beam == -1) in >> phony;
      else in >> theData.GetEntry(1, 0, beam, disc);
    }
  }

  // TEC- Ring 6
  if(headers) in >> dummy;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) in >> dummy;
      else {in >> dummy; in >> phony;}
    }
    else if(beam == -1)continue;
    for(int disc = 0; disc < 9; disc++){
      if(headers && beam == -1) in >> phony;
      else in >> theData.GetEntry(1, 1, beam, disc);
    }
  }

  // TIB
  if(headers) in >> dummy;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) in >> dummy;
      else {in >> dummy; in >> phony;}
    }
    else if(beam == -1)continue;
    for(int zpos = 0; zpos < 6; zpos++){
      if(headers && beam == -1) in >> phony;
      else in >> theData.GetEntry(2, -1, beam, zpos);
    }
  }

  // TOB
  if(headers) in >> dummy;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) in >> dummy;
      else {in >> dummy; in >> phony;}
    }
    else if(beam == -1)continue;
    for(int zpos = 0; zpos < 6; zpos++){
      if(headers && beam == -1) in >> phony;
      else in >> theData.GetEntry(3, -1, beam, zpos);
    }
  }


  // TEC+ AT
  if(headers) in >> dummy;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) in >> dummy;
      else {in >> dummy; in >> phony;}
    }
    else if(beam == -1)continue;
    for(int disc = 0; disc < 5; disc++){
      if(headers && beam == -1) in >> phony;
      else in >> theData.GetEntry(0, -1, beam, disc);
    }
  }

  // TEC- AT
  if(headers) in >> dummy;
  for(int beam = -1; beam < 8; beam ++){
    if(headers){
      if(beam == -1) in >> dummy;
      else {in >> dummy; in >> phony;}
    }
    else if(beam == -1)continue;
    for(int disc = 0; disc < 5; disc++){
      if(headers && beam == -1) in >> phony;
      else in >> theData.GetEntry(1, -1, beam, disc);
    }
  }
  return;
}

////////////////////////////////////////////////////
// Drawing of global data


template <class T>
//! Plot the simple types that are contained in the LASGlobalData object
void draw_global_data(const LASGlobalData<T>& theData, LAS::beam_group group = LAS::ALL)
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


/////////////////////////////////////////////
// Boolean Operators

template <class T>
LASGlobalData<T> operator&&(const LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = loop.GetEntry<T>(gd1) && loop.GetEntry<T>(gd2);} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator&&(const LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) && loop.GetEntry<T2>(gd2);} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> operator||(const LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = loop.GetEntry<T>(gd1) || loop.GetEntry<T>(gd2);} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> operator!(const LASGlobalData<T>& gd)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = !(loop.GetEntry<T>(gd));} while (loop.next() );
  return retval;
}

/////////////////////////////////////////////
// Comparison Operators

template <class T>
LASGlobalData<T> operator==(const LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = loop.GetEntry<T>(gd1) == loop.GetEntry<T>(gd2);} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> operator!=(const LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = loop.GetEntry<T>(gd1) != loop.GetEntry<T>(gd2);} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> operator<(const LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = loop.GetEntry<T>(gd1) < loop.GetEntry<T>(gd2);} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> operator>(const LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = loop.GetEntry<T>(gd1) > loop.GetEntry<T>(gd2);} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> operator<=(const LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = loop.GetEntry<T>(gd1) <= loop.GetEntry<T>(gd2);} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> operator>=(const LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = loop.GetEntry<T>(gd1) >= loop.GetEntry<T>(gd2);} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator>(const LASGlobalData<T1>& gd1, T2 comparison)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) > comparison;} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator<(const LASGlobalData<T1>& gd1, T2 comparison)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) < comparison;} while (loop.next() );
  return retval;
}


//////////////////////////
// Arithmetic Operators //
//////////////////////////

// Two LASGlobalData objects holding different types

/* template <class T> */
/* LASGlobalData<T> operator+(const LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2) */
/* { */
/*   LASGlobalData<T> retval; */
/*   LASGlobalDataLoop loop; */
/*   do{ loop.GetEntry<T>(retval) = loop.GetEntry<T>(gd1) + loop.GetEntry<T>(gd2);} while (loop.next() ); */
/*   return retval; */
/* } */

template <class T1, class T2>
LASGlobalData<T1> operator+(const LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) + loop.GetEntry<T2>(gd2);} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator-(const LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) - loop.GetEntry<T2>(gd2);} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator*(const LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) * loop.GetEntry<T2>(gd2);} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator/(const LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = ( loop.GetEntry<T2>(gd2) != 0 ? loop.GetEntry<T1>(gd1) / loop.GetEntry<T2>(gd2) : 0 );} while (loop.next() );
  return retval;
}


template <class T>
LASGlobalData<T> operator-(const LASGlobalData<T>& gd)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = - loop.GetEntry<T>(gd);} while (loop.next() );
  return retval;
}


// Inplace operators
template <class T>
void operator+=(LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(gd1) += loop.GetEntry<T>(gd2);} while (loop.next() );
}

template <class T1, class T2>
void operator+=(LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(gd1) += loop.GetEntry<T2>(gd2);} while (loop.next() );
}

template <class T>
void operator-=(LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(gd1) -= loop.GetEntry<T>(gd2);} while (loop.next() );
}

template <class T1, class T2>
void operator-=(LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(gd1) -= loop.GetEntry<T2>(gd2);} while (loop.next() );
}

template <class T1, class T2>
void operator*=(LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(gd1) *= loop.GetEntry<T2>(gd2);} while (loop.next() );
}

template <class T>
void operator/=(LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(gd1) /= (loop.GetEntry<T>(gd2) != 0 ? loop.GetEntry<T>(gd2) : 1.0 );} while (loop.next() );
}

template <class T1, class T2>
void operator/=(LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(gd1) /= loop.GetEntry<T2>(gd2);} while (loop.next() );
}

template <class T>
void operator&=(LASGlobalData<T>& gd1, const LASGlobalData<T>& gd2)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(gd1) &= loop.GetEntry<T>(gd2);} while (loop.next() );
}

template <class T1, class T2>
void operator&=(LASGlobalData<T1>& gd1, const LASGlobalData<T2>& gd2)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(gd1) &= (T1)loop.GetEntry<T2>(gd2);} while (loop.next() );
}


/////////////////////////////////////////////////////////////////////////////////
// Arithmetic with simple types


template <class T1, class T2>
LASGlobalData<T1> operator+(const LASGlobalData<T1>& gd1, T2 offset)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) + offset;} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator-(const LASGlobalData<T1>& gd1, T2 offset)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) - offset;} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator*(const LASGlobalData<T1>& gd1, T2 factor)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) * factor;} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator/(const LASGlobalData<T1>& gd1, T2 factor)
{
  LASGlobalData<T1> retval;
  if(factor != 0){
    LASGlobalDataLoop loop;
    do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) / factor;} while (loop.next() );
  }
  return retval;
}

template <class T1, class T2>
void operator+=(LASGlobalData<T1>& gd1, T2 offset)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(gd1) += offset;} while (loop.next() );
  return;
}

template <class T1, class T2>
void operator-=(LASGlobalData<T1>& gd1, T2 offset)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(gd1) -= offset;} while (loop.next() );
  return;
}

template <class T1, class T2>
void operator*=(LASGlobalData<T1>& gd1, T2 factor)
{
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(gd1) *= factor;} while (loop.next() );
  return;
}

template <class T1, class T2>
void operator/=(LASGlobalData<T1>& gd1, T2 factor)
{
  if(factor == 0) return;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(gd1) /= factor;} while (loop.next() );
  return;
}


////////////////////////////////////////////////////////
// Compare two different types

template <class T1, class T2>
LASGlobalData<T1> operator<(const LASGlobalData<T1>& gd1, const LASGlobalData<T2>& comparison)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) < loop.GetEntry<T2>(comparison);} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator>(const LASGlobalData<T1>& gd1, const LASGlobalData<T2>& comparison)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) > loop.GetEntry<T2>(comparison);} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator<=(const LASGlobalData<T1>& gd1, const LASGlobalData<T2>& comparison)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) <= loop.GetEntry<T2>(comparison);} while (loop.next() );
  return retval;
}

template <class T1, class T2>
LASGlobalData<T1> operator>=(const LASGlobalData<T1>& gd1, const LASGlobalData<T2>& comparison)
{
  LASGlobalData<T1> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T1>(retval) = loop.GetEntry<T1>(gd1) >= loop.GetEntry<T2>(comparison);} while (loop.next() );
  return retval;
}




/////////////////////////////////////////////////////////////////////////////////////////////
//// math.h functions

template <class T>
LASGlobalData<double> sqrt(const LASGlobalData<T>& gdata)
{
  LASGlobalData<double> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<double>(retval) = sqrt(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> sin(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = sin(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> cos(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = cos(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> tan(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = tan(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> asin(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = asin(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> acos(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = acos(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> atan(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = atan(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> sinh(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = sinh(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> cosh(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = cosh(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> tanh(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = tanh(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> ceil(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = ceil(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> floor(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = floor(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> fabs(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = fabs(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> exp(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = exp(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> log(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = log(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> log10(const LASGlobalData<T>& gdata)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = log10(loop.GetEntry<T>(gdata));} while (loop.next() );
  return retval;
}

template <class T>
LASGlobalData<T> pow(const LASGlobalData<T>& gdata, double expo)
{
  LASGlobalData<T> retval;
  LASGlobalDataLoop loop;
  do{ loop.GetEntry<T>(retval) = pow(loop.GetEntry<T>(gdata), expo);} while (loop.next() );
  return retval;
}


#endif
