#ifndef __INCLUDE_LAS_AT_REC_H__
#define __INCLUDE_LAS_AT_REC_H__

// Calculation of beam spot positions for AT modules
LASGlobalData<double> AT_spot_rec(const AtPar& par);

// Calculation of chi2/ndf for AT reconstruction
void AT_chi2(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2,  AtPar& par);

// Draw the AT alignment parameters
void DrawAtPar(std::vector<AtPar>& parlist, Avec& time);

// reconstruction of TIB wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TIB2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results);

// reconstruction of TEC+ wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TECP2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results);

// reconstruction of TEC- wrt TOB 6 alignment parameters: dx,dy, rx, ry, rz, tz and 16 beam parameters: (a&b)*8
void AT_TECM2TOB_6par(const LASGlobalData<double>& dif,const LASGlobalData<double>& err2, AtPar& results);


#endif
