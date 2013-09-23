#ifndef __LAS_ALPAR_H__
#define __LAS_ALPAR_H__

#include <iostream>
#define __AVECROOT__
#include "Avec.h"

class LasAlPar;
class TecPar;


LasAlPar& alpar_get(const std::string& name, TFile& file);
LasAlPar& alpar_get(const std::string& objname, const std::string& filename);

void separate_correlations_tec(Avec& v, double& offset, double& slope);
void separate_correlations_tec(TecPar& alpar, bool keep_flag = false);
void alpar_rand_gauss(TecPar& alpar, double sigma_dphik, double sigma_dxk, double sigma_dyk);
void alpar_rand_gauss(LasAlPar& alpar, double sigma_dphik, double sigma_dxk, double sigma_dyk);

class AtPar : public TNamed
{
public:
  AtPar();
public:
  double AT_chi; //chi^{2}/ndf for AT
  double AT_chi0; //chi^{2}/ndf for AT before fit

  Avec beam_a; //beam slopes
  Avec beam_b; //beam offsets
  Avec er_beam_a; //errors of beam slopes 
  Avec er_beam_b; //errors of beam offsets
  Avec b_chi; // chi2/ndf of beams
  Avec b_chi0; // chi2/ndf of beams before fit

  double tob_dx; //TOB x-displacement
  double tob_dy; //TOB y-displacement
  double tob_rx; //TOB x-rotation
  double tob_ry; //TOB y-rotation
  double tob_rz; //TOB z-rotation
  double tob_tz; //TOB z-torsion
  double er_tob_dx; //error of TOB x-displacement
  double er_tob_dy; //error of TOB y-displacement
  double er_tob_rx; //error of TOB x-rotation
  double er_tob_ry; //error of TOB y-rotation
  double er_tob_rz; //error of TOB z-rotation
  double er_tob_tz; //error of TOB z-torsion
  double tob_dxdy; //correlation TOB dx&dy
  double tob_dxrx; //correlation TOB dx&rx
  double tob_dxry; //correlation TOB dx&ry
  double tob_dxrz; //correlation TOB dx&rz
  double tob_dxtz; //correlation TOB dx&tz
  double tob_dyrx; //correlation TOB dy&rx
  double tob_dyry; //correlation TOB dy&ry
  double tob_dyrz; //correlation TOB dy&rz
  double tob_dytz; //correlation TOB dy&tz
  double tob_rxry; //correlation TOB rx&ry
  double tob_rxrz; //correlation TOB rx&rz
  double tob_rxtz; //correlation TOB rx&tz
  double tob_ryrz; //correlation TOB ry&rz
  double tob_rytz; //correlation TOB ry&tz
  double tob_rztz; //correlation TOB rz&tz
  double tob_chi; //chi^{2}/ndf for TOB
  double tob_chi0; //chi^{2}/ndf for TOB before fit

  double tib_dx; //TIB x-displacement
  double tib_dy; //TIB y-displacement
  double tib_rx; //TIB x-rotation
  double tib_ry; //TIB y-rotation
  double tib_rz; //TIB z-rotation
  double tib_tz; //TIB z-torsion
  double er_tib_dx; //error of TIB x-displacement
  double er_tib_dy; //error of TIB y-displacement
  double er_tib_rx; //error of TIB x-rotation
  double er_tib_ry; //error of TIB y-rotation
  double er_tib_rz; //error of TIB z-rotation
  double er_tib_tz; //error of TIB z-torsion
  double tib_dxdy; //correlation TIB dx&dy
  double tib_dxrx; //correlation TIB dx&rx
  double tib_dxry; //correlation TIB dx&ry
  double tib_dxrz; //correlation TIB dx&rz
  double tib_dxtz; //correlation TIB dx&tz
  double tib_dyrx; //correlation TIB dy&rx
  double tib_dyry; //correlation TIB dy&ry
  double tib_dyrz; //correlation TIB dy&rz
  double tib_dytz; //correlation TIB dy&tz
  double tib_rxry; //correlation TIB rx&ry
  double tib_rxrz; //correlation TIB rx&rz
  double tib_rxtz; //correlation TIB rx&tz
  double tib_ryrz; //correlation TIB ry&rz
  double tib_rytz; //correlation TIB ry&tz
  double tib_rztz; //correlation TIB rz&tz
  double tib_chi; //chi^{2}/ndf for TIB
  double tib_chi0; //chi^{2}/ndf for TIB before fit

  double tecp_dx; //TEC+ x-displacement
  double tecp_dy; //TEC+ y-displacement
  double tecp_rx; //TEC+ x-rotation
  double tecp_ry; //TEC+ y-rotation
  double tecp_rz; //TEC+ z-rotation
  double tecp_tz; //TEC+ z-torsion
  double er_tecp_dx; //error of TEC+ x-displacement
  double er_tecp_dy; //error of TEC+ y-displacement
  double er_tecp_rx; //error of TEC+ x-rotation
  double er_tecp_ry; //error of TEC+ y-rotation
  double er_tecp_rz; //error of TEC+ z-rotation
  double er_tecp_tz; //error of TEC+ z-torsion
  double tecp_dxdy; //correlation TEC+ dx&dy
  double tecp_dxrx; //correlation TEC+ dx&rx
  double tecp_dxry; //correlation TEC+ dx&ry
  double tecp_dxrz; //correlation TEC+ dx&rz
  double tecp_dxtz; //correlation TEC+ dx&tz
  double tecp_dyrx; //correlation TEC+ dy&rx
  double tecp_dyry; //correlation TEC+ dy&ry
  double tecp_dyrz; //correlation TEC+ dy&rz
  double tecp_dytz; //correlation TEC+ dy&tz
  double tecp_rxry; //correlation TEC+ rx&ry
  double tecp_rxrz; //correlation TEC+ rx&rz
  double tecp_rxtz; //correlation TEC+ rx&tz
  double tecp_ryrz; //correlation TEC+ ry&rz
  double tecp_rytz; //correlation TEC+ ry&tz
  double tecp_rztz; //correlation TEC+ rz&tz
  double tecp_chi; //chi^{2}/ndf for TEC+ 
  double tecp_chi0; //chi^{2}/ndf for TEC+ before fit

  double tecm_dx; //TEC- x-displacement
  double tecm_dy; //TEC- y-displacement
  double tecm_rx; //TEC- x-rotation
  double tecm_ry; //TEC- y-rotation
  double tecm_rz; //TEC- z-rotation
  double tecm_tz; //TEC- z-torsion
  double er_tecm_dx; //error of TEC- x-displacement
  double er_tecm_dy; //error of TEC- y-displacement
  double er_tecm_rx; //error of TEC- x-rotation
  double er_tecm_ry; //error of TEC- y-rotation
  double er_tecm_rz; //error of TEC- z-rotation
  double er_tecm_tz; //error of TEC- z-torsion
  double tecm_dxdy; //correlation TEC- dx&dy
  double tecm_dxrx; //correlation TEC- dx&rx
  double tecm_dxry; //correlation TEC- dx&ry
  double tecm_dxrz; //correlation TEC- dx&rz
  double tecm_dxtz; //correlation TEC- dx&tz
  double tecm_dyrx; //correlation TEC- dy&rx
  double tecm_dyry; //correlation TEC- dy&ry
  double tecm_dyrz; //correlation TEC- dy&rz
  double tecm_dytz; //correlation TEC- dy&tz
  double tecm_rxry; //correlation TEC- rx&ry
  double tecm_rxrz; //correlation TEC- rx&rz
  double tecm_rxtz; //correlation TEC- rx&tz
  double tecm_ryrz; //correlation TEC- ry&rz
  double tecm_rytz; //correlation TEC- ry&tz
  double tecm_rztz; //correlation TEC- rz&tz
  double tecm_chi; //chi^{2}/ndf for TEC- to TOB
  double tecm_chi0; //chi^{2}/ndf for TEC- before fit

  ClassDef( AtPar, 1 );
};

std::ostream& operator<<(std::ostream& out, const AtPar& ap);

class TecRingPar : public TObject
{
public:
  TecRingPar();
  void print();

  void operator+=(double offset);
  void operator-=(double offset);
  void operator*=(double factor); 
  void operator/=(double factor);

  void operator+=(const TecRingPar& p2);
  void operator-=(const TecRingPar& p2);
  void operator*=(const TecRingPar& p2);
  void operator/=(const TecRingPar& p2);

  TecRingPar operator+(double offset);
  TecRingPar operator-(double offset);
  TecRingPar operator*(double factor);
  TecRingPar operator/(double factor);

  TecRingPar operator+(const TecRingPar& p2);
  TecRingPar operator-(const TecRingPar& p2);
  TecRingPar operator*(const TecRingPar& p2);
  TecRingPar operator/(const TecRingPar& p2);

  double Dphi0, Dx0, Dy0, Dphit, Dxt, Dyt;
  Avec Dphik, Dxk, Dyk;
  Avec DthetaA, DthetaB;

  ClassDef( TecRingPar, 1 );
};

class TecPar : public TObject
{
public:
  TecPar();
  void print();

  void operator+=(double offset);
  void operator-=(double offset);
  void operator*=(double factor); 
  void operator/=(double factor);

  void operator+=(const TecPar& p2);
  void operator-=(const TecPar& p2);
  void operator*=(const TecPar& p2);
  void operator/=(const TecPar& p2);

  TecPar operator+(double offset);
  TecPar operator-(double offset);
  TecPar operator*(double factor);
  TecPar operator/(double factor);

  TecPar operator+(const TecPar& p2);
  TecPar operator-(const TecPar& p2);
  TecPar operator*(const TecPar& p2);
  TecPar operator/(const TecPar& p2);

  double Dphi0, Dx0, Dy0, Dphit, Dxt, Dyt;
  Avec Dphik, Dxk, Dyk;
  Avec DthetaA_r4, DthetaB_r4;
  Avec DthetaA_r6, DthetaB_r6;

  ClassDef( TecPar, 1 );
};



class LasAlPar : public TNamed
{
public:
  void print();

  void operator+=(double offset);
  void operator-=(double offset);
  void operator*=(double factor); 
  void operator/=(double factor);

  void operator+=(const LasAlPar& p2);
  void operator-=(const LasAlPar& p2);
  void operator*=(const LasAlPar& p2);
  void operator/=(const LasAlPar& p2);

  LasAlPar operator+(double offset);
  LasAlPar operator-(double offset);
  LasAlPar operator*(double factor);
  LasAlPar operator/(double factor);

  LasAlPar operator+(const LasAlPar& p2);
  LasAlPar operator-(const LasAlPar& p2);
  LasAlPar operator*(const LasAlPar& p2);
  LasAlPar operator/(const LasAlPar& p2);


  void RescaleRot(double factor);
  void RescaleTrans(double factor);

  TecRingPar tecp_r4;
  TecRingPar tecp_r6;
  TecRingPar tecm_r4;
  TecRingPar tecm_r6;

  TecPar tecm;
  TecPar tecp;

  AtPar at;

  ClassDef( LasAlPar, 3 );
};




#endif
