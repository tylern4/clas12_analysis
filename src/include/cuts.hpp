/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "branches.hpp"
#include "constants.hpp"
#include "deltat.hpp"
#include "physics.hpp"

class Cuts {
 protected:
  std::shared_ptr<Branches12> _data = nullptr;
  std::shared_ptr<Delta_T> _dt = nullptr;

 public:
  Cuts(const std::shared_ptr<Branches12>& data);
  Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt);
  ~Cuts();

  bool ElectronCuts();
  bool FiducialCuts();
  bool IsPip(int i);
  bool IsProton(int i);
  bool IsPim(int i);
};

class rga_Cuts : public Cuts {
  rga_Cuts(const std::shared_ptr<Branches12>& data) : Cuts(data) {}
  rga_Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt) : Cuts(data, dt){};
};

class rgf_Cuts : public Cuts {
  rgf_Cuts(const std::shared_ptr<Branches12>& data) : Cuts(data) {}
  rgf_Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt) : Cuts(data, dt){};
};

#endif
