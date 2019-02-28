/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "deltat.hpp"

Delta_T::Delta_T(std::shared_ptr<Clas12Branches> data) {
  _data = data;
  _time_1b_vert = data->sc_ftof_1b_time(0);
  _path_1b_vert = data->sc_ftof_1b_path(0);
  _time_1a_vert = data->sc_ftof_1a_time(0);
  _path_1a_vert = data->sc_ftof_1a_path(0);
  _time_2_vert = data->sc_ftof_2_time(0);
  _path_2_vert = data->sc_ftof_2_path(0);
  if (!std::isnan(_time_1b_vert)) {
    _sc_t_v = _time_1b_vert;
    _sc_r_v = _path_1b_vert;
  } else if (!std::isnan(_time_1a_vert)) {
    _sc_t_v = _time_1a_vert;
    _sc_r_v = _path_1a_vert;
  } else if (!std::isnan(_time_2_vert)) {
    _sc_t_v = _time_2_vert;
    _sc_r_v = _path_2_vert;
  }

  _vertex = _vertex_time(_sc_t_v, _sc_r_v, 1.0);
}

Delta_T::~Delta_T() {}

float Delta_T::_vertex_time(float sc_time, float sc_pathlength, float relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

float Delta_T::_deltat(int pid) {
  _beta = 1.0 / sqrt(1.0 + (_mass_map[pid] / _momentum) * (_mass_map[pid] / _momentum));
  if (_sc_t == _sc_t && _sc_r == _sc_r) {
    return _vertex - _vertex_time(_sc_t, _sc_r, _beta);
  } else {
    return NAN;
  }
}

void Delta_T::get_det_info(int i) {
  _momentum = _data->p(i);
  float _time_1b = _data->sc_ftof_1b_time(i);
  float _path_1b = _data->sc_ftof_1b_path(i);
  float _time_1a = _data->sc_ftof_1a_time(i);
  float _path_1a = _data->sc_ftof_1a_path(i);
  float _time_2 = _data->sc_ftof_2_time(i);
  float _path_2 = _data->sc_ftof_2_path(i);
  float _time_ctof = _data->sc_ctof_time(i);
  float _path_ctof = _data->sc_ctof_path(i);

  if (!std::isnan(_time_1b)) {
    _sc_t = _time_1b;
    _sc_r = _path_1b;
  } else if (!std::isnan(_time_1a)) {
    _sc_t = _time_1a;
    _sc_r = _path_1a;
  } else if (!std::isnan(_time_ctof)) {
    _sc_t = _time_ctof;
    _sc_r = _path_ctof;
  } else if (!std::isnan(_time_2)) {
    _sc_t = _time_2;
    _sc_r = _path_2;
  }
}

float Delta_T::dt_E(int i) {
  get_det_info(i);
  return _deltat(ELECTRON);
}
float Delta_T::dt_P(int i) {
  get_det_info(i);
  return _deltat(PROTON);
}
float Delta_T::dt_Pi(int i) {
  get_det_info(i);
  return _deltat(PIP);
}
float Delta_T::dt_K(int i) {
  get_det_info(i);
  return _deltat(KP);
}
float Delta_T::dt(int i, int pid) {
  get_det_info(i);
  return _deltat(pid);
}
float Delta_T::momentum(int i) { return _data->p(i); }
