#include "LAS_alpar.h"
#include "LAS_basic_tools.h"

#include "TFile.h"


LasAlPar& alpar_get(const std::string& name, TFile& file)
{
  static LasAlPar empty;

  LasAlPar* alpar_ptr = 0;
  file.GetObject(name.c_str(), alpar_ptr);
  if(!alpar_ptr){
    std::cerr << "Did not find LasAlPar object " << name << std::endl;
    return empty;
  }
  return *alpar_ptr;
}

LasAlPar& alpar_get(const std::string& objname, const std::string& filename)
{
  static LasAlPar empty;
  TFile f(filename.c_str(),"READ");
  if(!f.IsOpen()){
    std::cerr << "Could not open file " << filename << std::endl;
    return empty;
  }
  return alpar_get(objname, f);
}


void separate_correlations_tec(Avec& v, double& offset, double& slope)
{
  const Avec& zdat = LAS::zdat_tec;
  //zdat.print();
  double A = vsum(v);
  double B = v.size();
  double C = vsum(zdat);
  double D = vsum(v*zdat);
  double E = vsum(zdat);
  double F = vsum(zdat*zdat);

  double denom = B*F - C*E;

  offset = (C*D - A*F) / denom;
  slope  = (A*E - B*D) / denom;
  v += slope*zdat + offset;
}

void separate_correlations_tec(TecPar& alpar, bool keep_flag)
{
  double offset, slope;
  separate_correlations_tec(alpar.Dphik, offset, slope);
  if(keep_flag){
    alpar.Dphi0 += offset;
    alpar.Dphit += slope;
  }
  separate_correlations_tec(alpar.Dxk, offset, slope);
  if(keep_flag){
    alpar.Dx0 += offset;
    alpar.Dxt += slope;
  }
  separate_correlations_tec(alpar.Dyk, offset, slope);
  if(keep_flag){
    alpar.Dy0 += offset;
    alpar.Dyt += slope;
  }
}


void alpar_rand_gauss(TecPar& alpar, double sigma_dphik, double sigma_dxk, double sigma_dyk)
{
  alpar.Dphik = grand(9, sigma_dphik);
  alpar.Dxk = grand(9, sigma_dxk);
  alpar.Dyk = grand(9, sigma_dyk);

  Avec rand_phi = grand(2, sigma_dphik);
  Avec rand_x = grand(2, sigma_dxk);
  Avec rand_y = grand(2, sigma_dyk);

  alpar.Dphi0 = rand_phi[0];
  alpar.Dx0 = rand_x[0];
  alpar.Dy0 = rand_y[0];

  alpar.Dphit = rand_phi[1];
  alpar.Dxt = rand_x[1];
  alpar.Dyt = rand_y[1];
}

void alpar_rand_gauss(LasAlPar& alpar, double sigma_dphik, double sigma_dxk, double sigma_dyk)
{
  alpar_rand_gauss(alpar.tecp, sigma_dphik, sigma_dxk, sigma_dyk);
  alpar_rand_gauss(alpar.tecm, sigma_dphik, sigma_dxk, sigma_dyk);

}





///////////////////////////////////////////////////////////////////////////////////
// AtPar implementation

ClassImp(AtPar)

  Avec beam_a; //beam slopes
  Avec beam_b; //beam offsets
  Avec er_beam_a; //errors of beam slopes 
  Avec er_beam_b; //errors of beam offsets
  Avec b_chi; // chi2/ndf of beams
  Avec b_chi0; // chi2/ndf of beams before fit


AtPar::AtPar():
  beam_a(Avec(8)),
  beam_b(Avec(8)),
  er_beam_a(Avec(8)),
  er_beam_b(Avec(8)),
  b_chi(Avec(8)),
  b_chi0(Avec(8))
{

}

std::ostream& operator<<(std::ostream& out, const AtPar& ap)
{
  out << "TIB x-displacement: " << ap.tib_dx*1e3 << " +- " << ap.er_tib_dx << " micron" << std::endl; //
  out << "TIB y-displacement: " << ap.tib_dy*1e3 << " +- " << ap.er_tib_dy << " micron" << std::endl; //
  out << "TIB x-rotation: " << ap.tib_rx*1e6 << " +- " << ap.er_tib_rx << " microrad" << std::endl; //
  out << "TIB y-rotation: " << ap.tib_ry*1e6 << " +- " << ap.er_tib_ry << " microrad" << std::endl; //
  out << "TIB z-rotation: " << ap.tib_rz*1e6 << " +- " << ap.er_tib_rz << " microrad" << std::endl; //
  out << "TIB z-torsion: " << ap.tib_tz*1e9 << " +- " << ap.er_tib_tz << " microrad / m" << std::endl; //
  out << "TIB chi2: " << ap.tib_chi << " was " << ap.tib_chi0 << " before the fit" << std::endl;
  return out;
}

///////////////////////////////////////////////////////////////////////////////////
// TecPar implementation

ClassImp(TecPar)

TecPar::TecPar() :
  Dphi0(0), 
  Dx0(0), 
  Dy0(0), 
  Dphit(0), 
  Dxt(0), 
  Dyt(0),
  Dphik(Avec(9,0.0)),
  Dxk(Avec(9,0.0)),
  Dyk(Avec(9,0.0)),
  DthetaA_r4(Avec(8,0.0)),
  DthetaB_r4(Avec(8,0.0)),
  DthetaA_r6(Avec(8,0.0)),
  DthetaB_r6(Avec(8,0.0))
{;}

void TecPar::print()
{
  std::cout << "Dphi0 [mrad]: " << Dphi0*1000 << std::endl;
  std::cout << "  Dx0   [mm]: " << Dx0        << std::endl;
  std::cout << "  Dy0   [mm]: " << Dy0        << std::endl;
  
  std::cout << "Dphit [mrad]: " << Dphit*1000 << std::endl;
  std::cout << "  Dxt   [mm]: " << Dxt        << std::endl;
  std::cout << "  Dyt   [mm]: " << Dyt        << std::endl;
  
  std::cout << "Dphik [mrad]\n" << Dphik*1000 << std::endl;
  std::cout << "  Dxk   [mm]\n" <<   Dxk << std::endl;
  std::cout << "  Dyk   [mm]\n" <<   Dyk << std::endl;
  
  std::cout << "  DthetaA_r4 [mrad]\n" << DthetaA_r4*1000 << std::endl;
  std::cout << "  DthetaB_r4 [mrad]\n" << DthetaB_r4*1000 << std::endl;
  std::cout << "  DthetaA_r6 [mrad]\n" << DthetaA_r6*1000 << std::endl;
  std::cout << "  DthetaB_r6 [mrad]\n" << DthetaB_r6*1000 << std::endl;
}

//////////////////////////
// Arithmetic operators //
//////////////////////////
void TecPar::operator+=(double offset)
{
  Dphi0 += offset;
  Dx0 += offset;
  Dy0 += offset;
  Dphit += offset;
  Dxt += offset;
  Dyt += offset;
  Dphik += offset;
  Dxk += offset;
  Dyk += offset;
  DthetaA_r4 += offset;
  DthetaB_r4 += offset;
  DthetaA_r6 += offset;
  DthetaB_r6 += offset;
}

void TecPar::operator-=(double offset)
{
  Dphi0 -= offset;
  Dx0 -= offset;
  Dy0 -= offset;
  Dphit -= offset;
  Dxt -= offset;
  Dyt -= offset;
  Dphik -= offset;
  Dxk -= offset;
  Dyk -= offset;
  DthetaA_r4 -= offset;
  DthetaB_r4 -= offset;
  DthetaA_r6 -= offset;
  DthetaB_r6 -= offset;
}

void TecPar::operator*=(double factor)
{
  Dphi0 *= factor;
  Dx0 *= factor;
  Dy0 *= factor;
  Dphit *= factor;
  Dxt *= factor;
  Dyt *= factor;
  Dphik *= factor;
  Dxk *= factor;
  Dyk *= factor;
  DthetaA_r4 *= factor;
  DthetaB_r4 *= factor;
  DthetaA_r6 *= factor;
  DthetaB_r6 *= factor;
}

void TecPar::operator/=(double factor)
{
  Dphi0 /= factor;
  Dx0 /= factor;
  Dy0 /= factor;
  Dphit /= factor;
  Dxt /= factor;
  Dyt /= factor;
  Dphik /= factor;
  Dxk /= factor;
  Dyk /= factor;
  DthetaA_r4 /= factor;
  DthetaB_r4 /= factor;
  DthetaA_r6 /= factor;
  DthetaB_r6 /= factor;
}

void TecPar::operator+=(const TecPar& p2)
{
  Dphi0 += p2.Dphi0;
  Dx0 += p2.Dx0;
  Dy0 += p2.Dy0;
  Dphit += p2.Dphit;
  Dxt += p2.Dxt;
  Dyt += p2.Dyt;
  Dphik += p2.Dphik;
  Dxk += p2.Dxk;
  Dyk += p2.Dyk;
  DthetaA_r4 += p2.DthetaA_r4;
  DthetaB_r4 += p2.DthetaB_r4;
  DthetaA_r6 += p2.DthetaA_r6;
  DthetaB_r6 += p2.DthetaB_r6;
}

void TecPar::operator-=(const TecPar& p2)
{
  Dphi0 -= p2.Dphi0;
  Dx0 -= p2.Dx0;
  Dy0 -= p2.Dy0;
  Dphit -= p2.Dphit;
  Dxt -= p2.Dxt;
  Dyt -= p2.Dyt;
  Dphik -= p2.Dphik;
  Dxk -= p2.Dxk;
  Dyk -= p2.Dyk;
  DthetaA_r4 -= p2.DthetaA_r4;
  DthetaB_r4 -= p2.DthetaB_r4;
  DthetaA_r6 -= p2.DthetaA_r6;
  DthetaB_r6 -= p2.DthetaB_r6;
}

void TecPar::operator*=(const TecPar& p2)
{
  Dphi0 *= p2.Dphi0;
  Dx0 *= p2.Dx0;
  Dy0 *= p2.Dy0;
  Dphit *= p2.Dphit;
  Dxt *= p2.Dxt;
  Dyt *= p2.Dyt;
  Dphik *= p2.Dphik;
  Dxk *= p2.Dxk;
  Dyk *= p2.Dyk;
  DthetaA_r4 *= p2.DthetaA_r4;
  DthetaB_r4 *= p2.DthetaB_r4;
  DthetaA_r6 *= p2.DthetaA_r6;
  DthetaB_r6 *= p2.DthetaB_r6;
}

void TecPar::operator/=(const TecPar& p2)
{
  Dphi0 /= p2.Dphi0;
  Dx0 /= p2.Dx0;
  Dy0 /= p2.Dy0;
  Dphit /= p2.Dphit;
  Dxt /= p2.Dxt;
  Dyt /= p2.Dyt;
  Dphik /= p2.Dphik;
  Dxk /= p2.Dxk;
  Dyk /= p2.Dyk;
  DthetaA_r4 /= p2.DthetaA_r4;
  DthetaB_r4 /= p2.DthetaB_r4;
  DthetaA_r6 /= p2.DthetaA_r6;
  DthetaB_r6 /= p2.DthetaB_r6;
}

TecPar TecPar::operator+(double offset)
{
  TecPar trp(*this);
  trp += offset;
  return trp;
}

TecPar TecPar::operator-(double offset)
{
  TecPar trp(*this);
  trp -= offset;
  return trp;
}

TecPar TecPar::operator*(double factor)
{
  TecPar trp(*this);
  trp *= factor;
  return trp;
}

TecPar TecPar::operator/(double factor)
{
  TecPar trp(*this);
  trp /= factor;
  return trp;
}

TecPar TecPar::operator+(const TecPar& offset)
{
  TecPar trp(*this);
  trp += offset;
  return trp;
}

TecPar TecPar::operator-(const TecPar& offset)
{
  TecPar trp(*this);
  trp -= offset;
  return trp;
}

TecPar TecPar::operator*(const TecPar& factor)
{
  TecPar trp(*this);
  trp *= factor;
  return trp;
}

TecPar TecPar::operator/(const TecPar& factor)
{
  TecPar trp(*this);
  trp /= factor;
  return trp;
}



///////////////////////////////////////////////////////////////////////////////////
// TecRingPar implementation
TecRingPar::TecRingPar() :
  Dphi0(0),
  Dx0(0), 
  Dy0(0), 
  Dphit(0), 
  Dxt(0), 
  Dyt(0),
  Dphik(Avec(9,0.0)),
  Dxk(Avec(9,0.0)),
  Dyk(Avec(9,0.0)),
  DthetaA(Avec(8,0.0)),
  DthetaB(Avec(8,0.0))
{;}

void TecRingPar::print()
{
  std::cout << "Dphi0 [mrad]: " << Dphi0*1000 << std::endl;
  std::cout << "  Dx0   [mm]: " << Dx0        << std::endl;
  std::cout << "  Dy0   [mm]: " << Dy0        << std::endl;
  
  std::cout << "Dphit [mrad]: " << Dphit*1000 << std::endl;
  std::cout << "  Dxt   [mm]: " << Dxt        << std::endl;
  std::cout << "  Dyt   [mm]: " << Dyt        << std::endl;
  
  std::cout << "Dphik [mrad]\n" << Dphik*1000 << std::endl;
  std::cout << "  Dxk   [mm]\n" <<   Dxk << std::endl;
  std::cout << "  Dyk   [mm]\n" <<   Dyk << std::endl;
  
  std::cout << "  DthetaA [mrad]\n" << DthetaA*1000 << std::endl;
  std::cout << "  DthetaB [mrad]\n" << DthetaB*1000 << std::endl;
}

//////////////////////////
// Arithmetic operators //
//////////////////////////
void TecRingPar::operator+=(double offset)
{
  Dphi0 += offset;
  Dx0 += offset;
  Dy0 += offset;
  Dphit += offset;
  Dxt += offset;
  Dyt += offset;
  Dphik += offset;
  Dxk += offset;
  Dyk += offset;
  DthetaA += offset;
  DthetaB += offset;
}

void TecRingPar::operator-=(double offset)
{
  Dphi0 -= offset;
  Dx0 -= offset;
  Dy0 -= offset;
  Dphit -= offset;
  Dxt -= offset;
  Dyt -= offset;
  Dphik -= offset;
  Dxk -= offset;
  Dyk -= offset;
  DthetaA -= offset;
  DthetaB -= offset;
}

void TecRingPar::operator*=(double factor)
{
  Dphi0 *= factor;
  Dx0 *= factor;
  Dy0 *= factor;
  Dphit *= factor;
  Dxt *= factor;
  Dyt *= factor;
  Dphik *= factor;
  Dxk *= factor;
  Dyk *= factor;
  DthetaA *= factor;
  DthetaB *= factor;
}

void TecRingPar::operator/=(double factor)
{
  Dphi0 /= factor;
  Dx0 /= factor;
  Dy0 /= factor;
  Dphit /= factor;
  Dxt /= factor;
  Dyt /= factor;
  Dphik /= factor;
  Dxk /= factor;
  Dyk /= factor;
  DthetaA /= factor;
  DthetaB /= factor;
}

void TecRingPar::operator+=(const TecRingPar& p2)
{
  Dphi0 += p2.Dphi0;
  Dx0 += p2.Dx0;
  Dy0 += p2.Dy0;
  Dphit += p2.Dphit;
  Dxt += p2.Dxt;
  Dyt += p2.Dyt;
  Dphik += p2.Dphik;
  Dxk += p2.Dxk;
  Dyk += p2.Dyk;
  DthetaA += p2.DthetaA;
  DthetaB += p2.DthetaB;
}

void TecRingPar::operator-=(const TecRingPar& p2)
{
  Dphi0 -= p2.Dphi0;
  Dx0 -= p2.Dx0;
  Dy0 -= p2.Dy0;
  Dphit -= p2.Dphit;
  Dxt -= p2.Dxt;
  Dyt -= p2.Dyt;
  Dphik -= p2.Dphik;
  Dxk -= p2.Dxk;
  Dyk -= p2.Dyk;
  DthetaA -= p2.DthetaA;
  DthetaB -= p2.DthetaB;
}

void TecRingPar::operator*=(const TecRingPar& p2)
{
  Dphi0 *= p2.Dphi0;
  Dx0 *= p2.Dx0;
  Dy0 *= p2.Dy0;
  Dphit *= p2.Dphit;
  Dxt *= p2.Dxt;
  Dyt *= p2.Dyt;
  Dphik *= p2.Dphik;
  Dxk *= p2.Dxk;
  Dyk *= p2.Dyk;
  DthetaA *= p2.DthetaA;
  DthetaB *= p2.DthetaB;
}

void TecRingPar::operator/=(const TecRingPar& p2)
{
  Dphi0 /= p2.Dphi0;
  Dx0 /= p2.Dx0;
  Dy0 /= p2.Dy0;
  Dphit /= p2.Dphit;
  Dxt /= p2.Dxt;
  Dyt /= p2.Dyt;
  Dphik /= p2.Dphik;
  Dxk /= p2.Dxk;
  Dyk /= p2.Dyk;
  DthetaA /= p2.DthetaA;
  DthetaB /= p2.DthetaB;
}

TecRingPar TecRingPar::operator+(double offset)
{
  TecRingPar trp(*this);
  trp += offset;
  return trp;
}

TecRingPar TecRingPar::operator-(double offset)
{
  TecRingPar trp(*this);
  trp -= offset;
  return trp;
}

TecRingPar TecRingPar::operator*(double factor)
{
  TecRingPar trp(*this);
  trp *= factor;
  return trp;
}

TecRingPar TecRingPar::operator/(double factor)
{
  TecRingPar trp(*this);
  trp /= factor;
  return trp;
}

TecRingPar TecRingPar::operator+(const TecRingPar& offset)
{
  TecRingPar trp(*this);
  trp += offset;
  return trp;
}

TecRingPar TecRingPar::operator-(const TecRingPar& offset)
{
  TecRingPar trp(*this);
  trp -= offset;
  return trp;
}

TecRingPar TecRingPar::operator*(const TecRingPar& factor)
{
  TecRingPar trp(*this);
  trp *= factor;
  return trp;
}

TecRingPar TecRingPar::operator/(const TecRingPar& factor)
{
  TecRingPar trp(*this);
  trp /= factor;
  return trp;
}



///////////////////////////////////////////////////////////////////////////////////
// LasAlPar implementation

ClassImp(LasAlPar)

void LasAlPar::print()
{
  std::cout << "\n TEC+ full" << std::endl;
  tecp.print();
  std::cout << "\n TEC- full" << std::endl;
  tecm.print();

  std::cout << "\n TEC+ Ring4" << std::endl;
  tecp_r4.print();
  std::cout << "\n TEC+ Ring6" << std::endl;
  tecp_r6.print();
  std::cout << "\n TEC- Ring4" << std::endl;
  tecm_r4.print();
  std::cout << "\n TEC- Ring6" << std::endl;
  tecm_r6.print();
}

//////////////////////////
// Arithmetic operators //
//////////////////////////
void LasAlPar::operator+=(double offset)
{
  tecp_r4 += offset;
  tecp_r6 += offset;
  tecm_r4 += offset;
  tecm_r6 += offset;

  tecp += offset;
  tecm += offset;
}

void LasAlPar::operator-=(double offset)
{
  tecp_r4 -= offset;
  tecp_r6 -= offset;
  tecm_r4 -= offset;
  tecm_r6 -= offset;

  tecp -= offset;
  tecm -= offset;
}

void LasAlPar::operator*=(double factor)
{
  tecp_r4 *= factor;
  tecp_r6 *= factor;
  tecm_r4 *= factor;
  tecm_r6 *= factor;

  tecp *= factor;
  tecm *= factor;
}

void LasAlPar::operator/=(double factor)
{
  tecp_r4 /= factor;
  tecp_r6 /= factor;
  tecm_r4 /= factor;
  tecm_r6 /= factor;

  tecp /= factor;
  tecm /= factor;
}

void LasAlPar::operator+=(const LasAlPar& offset)
{
  tecp_r4 += offset.tecp_r4;
  tecp_r6 += offset.tecp_r6;
  tecm_r4 += offset.tecm_r4;
  tecm_r6 += offset.tecm_r6;

  tecp += offset.tecp;
  tecm += offset.tecm;
}

void LasAlPar::operator-=(const LasAlPar& offset)
{
  tecp_r4 -= offset.tecp_r4;
  tecp_r6 -= offset.tecp_r6;
  tecm_r4 -= offset.tecm_r4;
  tecm_r6 -= offset.tecm_r6;

  tecp -= offset.tecp;
  tecm -= offset.tecm;
}

void LasAlPar::operator*=(const LasAlPar& factor)
{
  tecp_r4 *= factor.tecp_r4;
  tecp_r6 *= factor.tecp_r6;
  tecm_r4 *= factor.tecm_r4;
  tecm_r6 *= factor.tecm_r6;

  tecp *= factor.tecp;
  tecm *= factor.tecm;
}

void LasAlPar::operator/=(const LasAlPar& factor)
{
  tecp_r4 /= factor.tecp_r4;
  tecp_r6 /= factor.tecp_r6;
  tecm_r4 /= factor.tecm_r4;
  tecm_r6 /= factor.tecm_r6;

  tecp /= factor.tecp;
  tecm /= factor.tecm;
}

LasAlPar LasAlPar::operator+(double offset)
{
  LasAlPar trp(*this);
  trp += offset;
  return trp;
}

LasAlPar LasAlPar::operator-(double offset)
{
  LasAlPar trp(*this);
  trp -= offset;
  return trp;
}

LasAlPar LasAlPar::operator*(double factor)
{
  LasAlPar trp(*this);
  trp *= factor;
  return trp;
}

LasAlPar LasAlPar::operator/(double factor)
{
  LasAlPar trp(*this);
  trp /= factor;
  return trp;
}

LasAlPar LasAlPar::operator+(const LasAlPar& offset)
{
  LasAlPar trp(*this);
  trp += offset;
  return trp;
}

LasAlPar LasAlPar::operator-(const LasAlPar& offset)
{
  LasAlPar trp(*this);
  trp -= offset;
  return trp;
}

LasAlPar LasAlPar::operator*(const LasAlPar& factor)
{
  LasAlPar trp(*this);
  trp *= factor;
  return trp;
}

LasAlPar LasAlPar::operator/(const LasAlPar& factor)
{
  LasAlPar trp(*this);
  trp /= factor;
  return trp;
}


void LasAlPar::RescaleRot(double factor)
{
  tecp_r4.Dphi0 *= factor;
  tecp_r4.Dphit *= factor;
  tecp_r4.Dphik *= factor;
  tecp_r6.Dphi0 *= factor;
  tecp_r6.Dphit *= factor;
  tecp_r6.Dphik *= factor;
  tecm_r4.Dphi0 *= factor;
  tecm_r4.Dphit *= factor;
  tecm_r4.Dphik *= factor;
  tecm_r6.Dphi0 *= factor;
  tecm_r6.Dphit *= factor;
  tecm_r6.Dphik *= factor;

  tecp.Dphi0 *= factor;
  tecp.Dphit *= factor;
  tecp.Dphik *= factor;
  tecm.Dphi0 *= factor;
  tecm.Dphit *= factor;
  tecm.Dphik *= factor;

  tecp_r4.DthetaA *= factor;
  tecp_r4.DthetaB *= factor;
  tecp_r6.DthetaA *= factor;
  tecp_r6.DthetaB *= factor;
  tecm_r4.DthetaA *= factor;
  tecm_r4.DthetaB *= factor;
  tecm_r6.DthetaA *= factor;
  tecm_r6.DthetaB *= factor;

  tecp.DthetaA_r4 *= factor;
  tecp.DthetaB_r4 *= factor;
  tecp.DthetaA_r6 *= factor;
  tecp.DthetaB_r6 *= factor;
  tecm.DthetaA_r4 *= factor;
  tecm.DthetaB_r4 *= factor;
  tecm.DthetaA_r6 *= factor;
  tecm.DthetaB_r6 *= factor;

}

void LasAlPar::RescaleTrans(double factor)
{
  tecp_r4.Dx0 *= factor;
  tecp_r4.Dy0 *= factor;
  tecp_r4.Dxt *= factor;
  tecp_r4.Dyt *= factor;
  tecp_r4.Dxk *= factor;
  tecp_r4.Dyk *= factor;

  tecp_r6.Dx0 *= factor;
  tecp_r6.Dy0 *= factor;
  tecp_r6.Dxt *= factor;
  tecp_r6.Dyt *= factor;
  tecp_r6.Dxk *= factor;
  tecp_r6.Dyk *= factor;

  tecm_r4.Dx0 *= factor;
  tecm_r4.Dy0 *= factor;
  tecm_r4.Dxt *= factor;
  tecm_r4.Dyt *= factor;
  tecm_r4.Dxk *= factor;
  tecm_r4.Dyk *= factor;

  tecm_r6.Dx0 *= factor;
  tecm_r6.Dy0 *= factor;
  tecm_r6.Dxt *= factor;
  tecm_r6.Dyt *= factor;
  tecm_r6.Dxk *= factor;
  tecm_r6.Dyk *= factor;

  tecp.Dx0 *= factor;
  tecp.Dy0 *= factor;
  tecp.Dxt *= factor;
  tecp.Dyt *= factor;
  tecp.Dxk *= factor;
  tecp.Dyk *= factor;

  tecm.Dx0 *= factor;
  tecm.Dy0 *= factor;
  tecm.Dxt *= factor;
  tecm.Dyt *= factor;
  tecm.Dxk *= factor;
  tecm.Dyk *= factor;

}

