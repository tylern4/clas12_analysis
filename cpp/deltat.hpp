/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DT_H_GUARD
#define DT_H_GUARD
#include "constants.hpp"

class DeltaT {
 private:
  std::vector<double> masses = {MASS_E, MASS_P, MASS_PIP};
  double vertex = 0.0;
  double dt_E = 0.0;
  double dt_P = 0.0;
  double dt_Pi = 0.0;

  double vertex_time(double sc_time, double sc_pathlength,
                     double relatavistic_beta);

 public:
  DeltaT(double sc_time, double sc_pathlength);
  ~DeltaT();

  void deltat(double momentum, double sc_t, double sc_r);
  double Get_dt_E();
  double Get_dt_P();
  double Get_dt_Pi();
};

#endif
