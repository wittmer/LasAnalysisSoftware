#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#include "LASGlobalData.h"
#include "LASGlobalDataLoop.h"
#include "LAS_basic_tools.h"
#include "LAS_globaldata_tools.h"

#define __AVECROOT__
#include "Avec.h"
#include "Avec2D.h"

#include "LAS_alpar.h"
#include "LAS_RDC_tools.h"
#include "LAS_control_plots.h"
#include "LAS_At_rec.h"

#pragma link C++ class LASGlobalDataLoop+;
#pragma link C++ function LASGlobalDataLoop::GetEntry<float>(LASGlobalData<float>&)+;
#pragma link C++ function LASGlobalDataLoop::GetEntry<double>(LASGlobalData<double>&)+;
#pragma link C++ function LASGlobalDataLoop::GetEntry<int>(LASGlobalData<int>&)+;
#pragma link C++ function LASGlobalDataLoop::GetEntry<Avec>(LASGlobalData<Avec>&)+;
#pragma link C++ function LASGlobalDataLoop::GetEntry< vector<float> >(LASGlobalData< vector<float> >&)+;
#pragma link C++ function LASGlobalDataLoop::GetEntry< vector<int> >(LASGlobalData< vector<int> >&)+;


#pragma link C++ class LASGlobalData< std::vector<int> >+;
#pragma link C++ class LASGlobalData< std::vector<float> >+;
#pragma link C++ class LASGlobalData< std::vector<double> >+;
#pragma link C++ class LASGlobalData<float>+;
#pragma link C++ class LASGlobalData<double>+;
#pragma link C++ class LASGlobalData<int>+;
#pragma link C++ class LASGlobalData<Avec>+;
#pragma link C++ class LASGlobalData<Avec2D>+;

#pragma link C++ function print_global_data<int>(const LASGlobalData<int>&)+;
#pragma link C++ function print_global_data<double>(const LASGlobalData<double>&)+;

#pragma link C++ function global_data_get<int>(const std::string&, TFile&)+;
#pragma link C++ function global_data_get<int>(const std::string&, const std::string&)+;

#pragma link C++ function global_data_get<double>(const std::string&, TFile&)+;
#pragma link C++ function global_data_get<double>(const std::string&, const std::string&)+;

#pragma link C++ function global_data_get<Avec>(const std::string&, TFile&)+;
#pragma link C++ function global_data_get<Avec>(const std::string&, const std::string&)+;

#pragma link C++ function operator+(const LASGlobalData<double>&, double)+;


#pragma link C++ class AtPar+;
#pragma link C++ class TecPar+;
#pragma link C++ class TecRingPar+;
#pragma link C++ class LasAlPar+;
#pragma link C++ class LAS::AnalysisParameters+;

#pragma link C++ namespace LAS;

#pragma link C++ defined_in "LAS_globaldata_tools.h";
#pragma link C++ defined_in "LAS_basic_tools.h";
#pragma link C++ defined_in "LAS_alpar.h";
#pragma link C++ defined_in "LAS_At_rec.h";
#pragma link C++ defined_in "LAS_control_plots.h";
#pragma link C++ defined_in "LAS_calibration.h";
#pragma link C++ defined_in "LAS_data_processing.h";

#endif
