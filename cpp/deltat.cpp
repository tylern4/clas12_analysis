/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "deltat.hpp"
#include "constants.hpp"

Delta_T::Delta_T(double sc_time, double sc_pathlength) {
  this->vertex = vertex_time(sc_time, sc_pathlength, 1.0);
}

Delta_T::~Delta_T() {}

double Delta_T::vertex_time(double sc_time, double sc_pathlength,
                           double relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

void Delta_T::deltat(double momentum, double sc_t, double sc_r) {
  double beta = 0.0;
  beta = 1.0 / sqrt(1.0 +
                    (this->masses.at(0) / momentum) *
                        (this->masses.at(0) / momentum));
  this->dt_E = this->vertex - vertex_time(sc_t, sc_r, beta);

  beta = 1.0 / sqrt(1.0 +
                    (this->masses.at(1) / momentum) *
                        (this->masses.at(1) / momentum));
  this->dt_P = this->vertex - vertex_time(sc_t, sc_r, beta);

  beta = 1.0 / sqrt(1.0 +
                    (this->masses.at(2) / momentum) *
                        (this->masses.at(2) / momentum));
  this->dt_Pi = this->vertex - vertex_time(sc_t, sc_r, beta);
}

double Delta_T::Get_dt_E() { return this->dt_E; }
double Delta_T::Get_dt_P() { return this->dt_P; }
double Delta_T::Get_dt_Pi() { return this->dt_Pi; }
