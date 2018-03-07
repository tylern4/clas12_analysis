/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "deltat.hpp"
#include "constants.hpp"

Delta_T::Delta_T(double sc_time, double sc_pathlength) { vertex = vertex_time(sc_time, sc_pathlength, 1.0); }

Delta_T::~Delta_T() {}

double Delta_T::vertex_time(double sc_time, double sc_pathlength, double relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

void Delta_T::deltat(double momentum, double sc_t, double sc_r) {
  double beta = 0.0;
  beta = 1.0 / sqrt(1.0 + (masses.at(0) / momentum) * (masses.at(0) / momentum));
  dt_E = vertex - vertex_time(sc_t, sc_r, beta);

  beta = 1.0 / sqrt(1.0 + (masses.at(1) / momentum) * (masses.at(1) / momentum));
  dt_P = vertex - vertex_time(sc_t, sc_r, beta);

  beta = 1.0 / sqrt(1.0 + (masses.at(2) / momentum) * (masses.at(2) / momentum));
  dt_Pi = vertex - vertex_time(sc_t, sc_r, beta);

  beta = 1.0 / sqrt(1.0 + (masses.at(3) / momentum) * (masses.at(3) / momentum));
  dt_K = vertex - vertex_time(sc_t, sc_r, beta);
}

double Delta_T::Get_dt_E() { return dt_E; }
double Delta_T::Get_dt_P() { return dt_P; }
double Delta_T::Get_dt_Pi() { return dt_Pi; }
double Delta_T::Get_dt_K() { return dt_K; }
