from __future__ import print_function
import os

get_id = {'PROTON': 2212, 'NEUTRON': 2112, 'PIP': 211, 'PIM': -211,
          'PI0': 111, 'KP': 321, 'KM': -321, 'PHOTON': 22, 'ELECTRON': 11}

masses = {11: 0.000511, 211: 0.13957, -211: 0.13957, 2212: 0.93827,
          2112: 0.939565, 321: 0.493667, -321: 0.493667, 22: 0}

CLAS12_E0 = 10.6 if not os.getenv('BEAM_E') else float(os.getenv('BEAM_E'))

MAX_PARTS = 100
PI = 3.14159
SOL = 29.9792458
FSC = 0.00729735253
NA = 6.02214129E23
QE = 1.60217646E-19

MASS_P = 0.93827203
MASS_N = 0.93956556
MASS_E = 0.000511
MASS_PIP = 0.13957018
MASS_PIM = 0.13957018
MASS_PI0 = 0.1349766
MASS_KP = 0.493677
MASS_KM = 0.493677
MASS_G = 0.0
MASS_OMEGA = 0.78265


def get_mass(x): return masses[get_id[x]]
def getM2(x): return [masses[pid]**2 for pid in x]
def DegToRad(x): return x * (PI / 180.0)
def RadToDeg(x): return x * (180.0 / PI)
