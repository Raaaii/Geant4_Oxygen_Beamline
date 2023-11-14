"""
    For  additional questions or comments, please contact the authors:
    p.r.rajka@gmail.com (Rajka Pejanovic) and 
    pablo.cirrone@cern.ch (Pablo Cirrone)
"""

# ========================

from scipy import stats
import scipy.stats.distributions
import matplotlib.pyplot as plt

import numpy
import numpy as np
import os
import math

import pandas as pd

import matplotlib.mlab as mlab
import matplotlib
matplotlib.get_backend()
import numpy.polynomial.polynomial as poly

import scipy.integrate as integrate

from usefulFunctions import getColumns

voxelSizeX = 0.1  # mm
voxelSizeY = 0.1  # mm
voxelSizeZ = 0.1  # mm

file_LET_Binary = open("Let1m.out", 'r')
cols_LET_Binary, indexToName_LET_Binary = getColumns(file_LET_Binary)
file_LET_Binary.close()



numberOfColumns_LET_Binary = len(indexToName_LET_Binary)

columnsHeaderFull_LET_Binary = {}
columnsContentFull_LET_Binary = {}

for i in range(0, numberOfColumns_LET_Binary):
    header_LET_Binary = indexToName_LET_Binary[i]
    columnsHeaderFull_LET_Binary[i] = header_LET_Binary

    
    columnsContent_LET_Binary = cols_LET_Binary[columnsHeaderFull_LET_Binary[i]]
    columnsContentArray_LET_Binary = np.array(columnsContent_LET_Binary)
    columnsContentArrayNumber_LET_Binary = columnsContentArray_LET_Binary.astype(np.float)

    columnsContentFull_LET_Binary[i] = columnsContentArrayNumber_LET_Binary

depthX_Binary = columnsContentFull_LET_Binary[0] * voxelSizeX
depthY_Binary = columnsContentFull_LET_Binary[1] * voxelSizeY
depthZ_Binary = columnsContentFull_LET_Binary[2] * voxelSizeZ



LET_Dose_Total_Binary = columnsContentFull_LET_Binary[3]
LET_Track_Total_Binary = columnsContentFull_LET_Binary[4]


Hydrogen_list = ['proton_D','deuteron_D','triton_D','H4_D','H5_D','H6_D','H7_D','H8_D','H9_D','H10_D','H11_D']
print(columnsHeaderFull_LET_Binary.values())
Helium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "He" in s and "_D" in s]
Helium_list += ['alpha_D']
Lithium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Li" in s and "_D" in s]
Beryllium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Be" in s and "_D" in s]
Borum_list = [s for s in columnsHeaderFull_LET_Binary.values() if "B" in s and "_D" in s and "Be" not in s]
Carbon_list = [s for s in columnsHeaderFull_LET_Binary.values() if "C" in s and "_D" in s]
Nitrogen_list = [s for s in columnsHeaderFull_LET_Binary.values() if "N" in s and "_D" in s and "Ne" not in s and "Na" not in s]
Oxygen_list = [s for s in columnsHeaderFull_LET_Binary.values() if "O" in s and "_D" in s]
Fluorine_list = [s for s in columnsHeaderFull_LET_Binary.values() if "F" in s and "_D" in s]
Neon_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Ne" in s and "_D" in s]
Sodium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Na" in s and "_D" in s]
Magnesium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Mg" in s and "_D" in s]
Aluminium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Al" in s and "_D" in s]
Silicon_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Si" in s and "_D" in s]
Sulfur_list = [s for s in columnsHeaderFull_LET_Binary.values() if "S" in s and "_D" in s and "Si" not in s]

print(Hydrogen_list)
print(Helium_list)
print(Lithium_list)
print(Beryllium_list)
print(Borum_list)
print(Carbon_list)
print(Nitrogen_list)
print(Oxygen_list)
print(Fluorine_list)
print(Neon_list)
print(Sodium_list)
print(Magnesium_list)
print(Aluminium_list)
#print(Silicon_list)
#print(Sulfur_list)

Hydrogen_Binary = [0]*len(depthX_Binary)
Helium_Binary = [0]*len(depthX_Binary)
Lithium_Binary = [0]*len(depthX_Binary)
Beryllium_Binary = [0]*len(depthX_Binary)
Borum_Binary = [0]*len(depthX_Binary)
Carbon_Binary = [0]*len(depthX_Binary)
Nitrogen_Binary = [0]*len(depthX_Binary)
Oxygen_Binary = [0]*len(depthX_Binary)
Fluorine_Binary = [0]*len(depthX_Binary)
Neon_Binary = [0]*len(depthX_Binary)
Sodium_Binary = [0]*len(depthX_Binary)
Magnesium_Binary = [0]*len(depthX_Binary)
Aluminium_Binary = [0]*len(depthX_Binary)
Silicon_Binary = [0]*len(depthX_Binary)
Sulfur_Binary = [0]*len(depthX_Binary)
for i in range(1,numberOfColumns_LET_Binary,2):
    if columnsHeaderFull_LET_Binary[i] in Hydrogen_list:
        Hydrogen_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Helium_list:
        Helium_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Lithium_list:
        Lithium_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Beryllium_list:
        Beryllium_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Borum_list:
        Borum_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Carbon_list:
        Carbon_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Nitrogen_list:
        Nitrogen_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Oxygen_list:
        Oxygen_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Fluorine_list:
        Fluorine_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Neon_list:
        Neon_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Sodium_list:
        Sodium_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Magnesium_list:
        Magnesium_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Aluminium_list:
        Aluminium_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Silicon_list:
        Silicon_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Sulfur_list:
        Sulfur_Binary += columnsContentFull_LET_Binary[i]

figure51, ax = plt.subplots(figsize=(14, 10))

ax.scatter(depthX_Binary, Hydrogen_Binary,label='Binary_Hydrogen', color='C0')
ax.scatter(depthX_Binary, Helium_Binary,label='Binary_Helium', color='C1')
ax.scatter(depthX_Binary, Lithium_Binary,label='Binary_Lithium',  color='C2')
ax.scatter(depthX_Binary, Beryllium_Binary,label='Binary_Beryllium',  color='C3')
ax.scatter(depthX_Binary, Borum_Binary,label='Binary_Borum',  color='C4')
ax.scatter(depthX_Binary, Carbon_Binary,label='Binary_Carbon',  color='C5')
ax.scatter(depthX_Binary, Nitrogen_Binary,label='Binary_Nitrogen', color='C6')
ax.scatter(depthX_Binary, Oxygen_Binary,label='Binary_Oxygen',  color='C7')
ax.scatter(depthX_Binary, Fluorine_Binary,label='Binary_Fluorine', color='C8')
ax.scatter(depthX_Binary, Neon_Binary,label='Binary_Neon', color='C9')

plt.title('LET Dose of secondaries (no primaries) grouped by element')
plt.xlabel('Depth in water [ mm ]')
plt.ylabel('LET [ keV/um ]')
plt.ylim(0,8000)
plt.legend(fontsize = 13, loc='upper left',bbox_to_anchor=(1, 1.05), prop={'size':6})
plt.draw()




Hydrogen_list = ['proton_T','deuteron_T','triton_T','H4_T','H5_T','H6_T','H7_T','H8_T','H9_T','H10_T','H11_T']
Helium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "He" in s and "_T" in s]
Helium_list += ['alpha_T']
Lithium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Li" in s and "_T" in s]
Beryllium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Be" in s and "_T" in s]
Borum_list = [s for s in columnsHeaderFull_LET_Binary.values() if "B" in s and "_T" in s and "Be" not in s]
Carbon_list = [s for s in columnsHeaderFull_LET_Binary.values() if "C" in s and "_T" in s]
Nitrogen_list = [s for s in columnsHeaderFull_LET_Binary.values() if "N" in s and "_T" in s and "Ne" not in s and "Na" not in s]
Oxygen_list = [s for s in columnsHeaderFull_LET_Binary.values() if "O" in s and "_T" in s]
Fluorine_list = [s for s in columnsHeaderFull_LET_Binary.values() if "F" in s and "_T" in s]
Neon_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Ne" in s and "_T" in s]
Sodium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Na" in s and "_T" in s]
Magnesium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Mg" in s and "_T" in s]
Aluminium_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Al" in s and "_T" in s]
Silicon_list = [s for s in columnsHeaderFull_LET_Binary.values() if "Si" in s and "_T" in s]
Sulfur_list = [s for s in columnsHeaderFull_LET_Binary.values() if "S" in s and "_T" in s and "Si" not in s]

Hydrogen_T_Binary = [0]*len(depthX_Binary)
Helium_T_Binary = [0]*len(depthX_Binary)
Lithium_T_Binary = [0]*len(depthX_Binary)
Beryllium_T_Binary = [0]*len(depthX_Binary)
Borum_T_Binary = [0]*len(depthX_Binary)
Carbon_T_Binary = [0]*len(depthX_Binary)
Nitrogen_T_Binary = [0]*len(depthX_Binary)
Oxygen_T_Binary = [0]*len(depthX_Binary)
Fluorine_T_Binary = [0]*len(depthX_Binary)
Neon_T_Binary = [0]*len(depthX_Binary)
Sodium_T_Binary = [0]*len(depthX_Binary)
Magnesium_T_Binary = [0]*len(depthX_Binary)
Aluminium_T_Binary = [0]*len(depthX_Binary)
Silicon_T_Binary = [0]*len(depthX_Binary)
Sulfur_T_Binary = [0]*len(depthX_Binary)
for i in range(2,numberOfColumns_LET_Binary,2):
    if columnsHeaderFull_LET_Binary[i] in Hydrogen_list:
        Hydrogen_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Helium_list:
        Helium_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Lithium_list:
        Lithium_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Beryllium_list:
        Beryllium_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Borum_list:
        Borum_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Carbon_list:
        Carbon_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Nitrogen_list:
        Nitrogen_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Oxygen_list:
        Oxygen_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Fluorine_list:
        Fluorine_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Neon_list:
        Neon_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Sodium_list:
        Sodium_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Magnesium_list:
        Magnesium_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Aluminium_list:
        Aluminium_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Silicon_list:
        Silicon_T_Binary += columnsContentFull_LET_Binary[i]
    if columnsHeaderFull_LET_Binary[i] in Sulfur_list:
        Sulfur_T_Binary += columnsContentFull_LET_Binary[i]

figure60, ax = plt.subplots(figsize=(14, 10))

ax.scatter(depthX_Binary, Hydrogen_T_Binary,label='Binary_Hydrogen', color='C0')
ax.scatter(depthX_Binary, Helium_T_Binary,label='Binary_Helium', color='C1')
ax.scatter(depthX_Binary, Lithium_T_Binary,label='Binary_Lithium',  color='C2')
ax.scatter(depthX_Binary, Beryllium_T_Binary,label='Binary_Beryllium',  color='C3')
ax.scatter(depthX_Binary, Borum_T_Binary,label='Binary_Borum', color='C4')
ax.scatter(depthX_Binary, Carbon_T_Binary,label='Binary_Carbon',  color='C5')
ax.scatter(depthX_Binary, Nitrogen_T_Binary,label='Binary_Nitrogen', color='C6')
ax.scatter(depthX_Binary, Oxygen_T_Binary,label='Binary_Oxygen', color='C7')
ax.scatter(depthX_Binary, Fluorine_T_Binary,label='Binary_Fluorine',  color='C8')
ax.scatter(depthX_Binary, Neon_T_Binary,label='Binary_Neon',  color='C9')

plt.title('LET Total track of secondaries (no primaries) grouped by element')
plt.xlabel('Depth in water [ mm ]')
plt.ylabel('LET [ keV/um ]')
#limit y to 0-100
plt.ylim(0,8000)
plt.legend(fontsize = 13, loc='upper left',bbox_to_anchor=(1, 1.05), prop={'size':6})
plt.draw()


plt.show()

