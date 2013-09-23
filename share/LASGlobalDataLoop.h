#ifndef __LASGLOBALDATALOOP_H
#define __LASGLOBALDATALOOP_H

#include <stdexcept>

#include "LASGlobalData.h"
//#include <iostream>

///
/// helper class for looping over LASGlobalData objects
/// Use like this, where T can be any class or type supported by LASGlobalData:
/// \code
/// LASGlobalData<T> mydata;
/// // Fill mydata with something...
///
/// LASGlobalDataLoop theLoop();
/// do
/// {
///   T& entry_ref = theLoop<T>.GetEntry(mydata);
///   // Now entry_ref is refering to a valid entry
/// }while ( theLoop.next() );
/// 
/// // Alternative:
/// for( LASGlobalDataLoop theLoop(); ! theLoop.finished(); theLoop.next()){
///   T& entry_ref = theLoop.GetEntry<T>(mydata);
///   // Now entry_ref is refering to a valid entry
/// }
/// \endcode
///

class LASGlobalDataLoop {
 public:
  enum loop_type{ALL, TEC_PLUS, TEC_MINUS, TEC, AT, TIB, TOB, TEC_PLUS_AT, TEC_MINUS_AT, TEC_AT, TEC_PLUS_R4, TEC_PLUS_R6, TEC_MINUS_R4, TEC_MINUS_R6};
  LASGlobalDataLoop(loop_type lp_tp = ALL);
  bool next();
  bool finished(){return loop_finished;}
  template <class T> T& GetEntry(LASGlobalData<T>& data){return data.GetEntry(det, ring, beam, zpos);}
  template <class T> const T& GetEntry(const LASGlobalData<T>& data) {return data.GetEntry(det, ring, beam, zpos);}
  void inspect(std::ostream & out = std::cout);
  int get_det() const {return det;}
  int get_ring() const {return ring;}
  int get_beam() const {return beam;}
  int get_zpos() const {return zpos;}
 private:
  loop_type the_loop_type;
  int det;
  int beam;
  int ring;
  int zpos;
  bool loop_finished;
  int max_det;
  int max_beam;
  int max_ring;
  int max_zpos;
};



#endif
