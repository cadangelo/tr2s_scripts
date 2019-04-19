#/usr/env python
from pyne.material import Material, MaterialLibrary, MultiMaterial
import os

matlib = MaterialLibrary()


Lead = Material({'Pb': 1.00}) 
Lead = Lead.expand_elements()
Lead.name = "Lead"
Lead.density = 11.36
matlib["Lead"] = Lead
print Lead

Iron = Material({'Fe': 1.00})
Iron = Iron.expand_elements()
Iron.name = "Iron"
Iron.density = 7.87
matlib["Iron"] = Iron
print Iron

#espi metals
ss316 = Material({'Fe': 65.34, 'Cr': 17, 'Ni': 12, 'Mo': 2.5, 'Mn': 2.0, 'Si':1.00, 'C12':0.08, 'P':0.05, 'S': 0.03 })
ss316 = ss316.expand_elements()
#ss316.name = "ss316"
ss316.name = "m1"
ss316.density = 8.03
#matlib["ss316"] = ss316
matlib["m1"] = ss316
print ss316

#SS316 = Material({'Fe': 65.34, 'Cr': 17, 'Ni': 12, 'Mo': 2.5, 'Mn': 2.0, 'Si':1.00, 'C12':0.08, 'P':0.05, 'S': 0.03 })
#SS316 = SS316.expand_elements()
#SS316.name = "SS316"
#SS316.density = 8.03
#matlib["SS316"] = SS316
#print SS316

m2 = Material({'Fe': 65.34, 'Cr': 17, 'Ni': 12, 'Mo': 2.5, 'Mn': 2.0, 'Si':1.00, 'C12':0.08, 'P':0.05, 'S': 0.03 })
m2 = m2.expand_elements()
m2.name = "mat:m2/rho:8.0E+00"
m2.density = 8.03
#matlib["mat:m2/rho:8.0E+00"] = m2
print m2

air = Material({'C12': 0.000124, 'N':0.755268, 'O':0.231781, 'Ar':0.012827})
air = air.expand_elements()
air.name = "air"
air.density = 0.001205
matlib["air"] = air
print air

air_comp = Material({'C12': 0.000124, 'N':0.755268, 'O':0.231781, 'Ar':0.012827})
air_comp = air_comp.expand_elements()
air_comp.name = "air_comp"
air_comp.density = 0.001205
matlib["air_comp"] = air_comp
print air_comp

ic_comp = Material({'C12': 0.000124, 'N':0.755268, 'O':0.231781, 'Ar':0.012827})
ic_comp = ic_comp.expand_elements()
ic_comp.name = "ic_comp"
ic_comp.density = 0.001205
matlib["ic_comp"] = ic_comp
matlib["m9"] = ic_comp
print ic_comp

ic = Material({'C12': 0.000124, 'N':0.755268, 'O':0.231781, 'Ar':0.012827})
ic = ic.expand_elements()
ic.name = "ic"
ic.density = 0.001205
matlib["ic"] = ic
print ic

m3 = Material({'C12': 0.000124, 'N':0.755268, 'O':0.231781, 'Ar':0.012827})
m3 = m3.expand_elements()
m3.name = "mat:m3/rho:1.2E-03"
m3.density = 0.001205
#matlib["mat:m3/rho:1.2E-03"] = m3 
print m3

detector = Material({'C12': 0.915, 'H': 0.085})
detector = detector.expand_elements()
detector.name = "detector"
detector.density = 1.037 
matlib["detector"] = detector
matlib["m2"] = detector
print detector

#m1 = Material({'C12': 0.915, 'H': 0.085})
#m1 = m1.expand_elements()
#m1.name = "mat:m1/rho:1.0E+00"
#m1.density = 1.037 
#matlib["mat:m1/rho:1.0E+00"] = m1
#print m1


he = Material({'He': 1.0})
he = he.expand_elements()
he.name = "he"
he.density = 0.000166
matlib["he"] = he
matlib["m3"] = he
print he

#m4 = Material({'He': 1.0})
#m4 = m4.expand_elements()
#m4.name = "mat:m4/rho:1.7E-4"
#m4.density = 0.000166
#matlib["m4"] = m4
#print m4

#steel_he = MultiMaterial({he:0.75, ss316:0.25}).mix_by_mass()
#steel_he.name = "steel_he"
#matlib["steel_he"] = steel_he
m4 = MultiMaterial({he:0.75, ss316:0.25}).mix_by_volume()
m4.name = "m4"
matlib["m4"] = m4
print m4


matlib.write_hdf5("basic_mat_lib_3.h5")
