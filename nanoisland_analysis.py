"""
Nanoisland Analysis Demo
========================

Post-processing workflow for ImageJ Analyze Particles CSV output.

The script reads an ImageJ results table and calculates:
- equivalent island diameter
- area statistics
- edge-to-edge spacing statistics
- circularity
- roundness
- island density

Expected CSV columns:
    Area, X, Y, Circ., Round
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

cm = 1 / 2.54

# =========================
# User input
# =========================
directory = "sample_data/"
filename = "Results.csv"

sizeTh = 50       # particle size threshold in nm
image_size_x = 3  # um
image_size_y = 2.05  # um

# =========================
# Load ImageJ CSV
# =========================
filepath = directory + filename

data_col_raw = pd.read_csv(filepath, sep=',')
df = pd.DataFrame(data_col_raw)

area_1 = df['Area']
X_1 = df['X']
Y_1 = df['Y']
circ_1 = df['Circ.']
Round_1 = df['Round']

# =========================
# Remove small particles
# =========================
diameter_nm = ((area_1 / 3.14) ** 0.5) * 2000

sort_diameter = []
area_2 = []
circ = []
Round = []
X = []
Y = []

for n in range(len(diameter_nm)):
    if diameter_nm[n] > sizeTh:
        sort_diameter.append(diameter_nm[n])
        area_2.append(area_1[n])
        circ.append(circ_1[n])
        Round.append(Round_1[n])
        X.append(X_1[n])
        Y.append(Y_1[n])

area = area_2
area_3 = [x * 1000000 for x in area_2]  # intended conversion to nm^2
ave_diameter = np.sum(sort_diameter) / len(sort_diameter)

# =========================
# Area statistics
# =========================
ave_area = np.sum(area_3) / len(area_3)

area_d = []
for n in range(len(area_3)):
    area_d.append((area_3[n] - ave_area) ** 2)

STD_area = (sum(area_d) / len(area_3)) ** 0.5
SE_area = STD_area / math.sqrt(len(area_3))

fig, ax = plt.subplots(figsize=(9 * cm, 8 * cm))
plt.hist(area_3, edgecolor="red", bins=10)
ax.tick_params(axis="both", which='both', direction="in", labelsize=15)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.xlabel('Area (nm$^2$)', size=15)
plt.ylabel('Counts', size=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.tight_layout()
plt.savefig('area-histogram.png', dpi=300)
plt.show()

# =========================
# Edge-to-edge spacing calculation
# =========================
r_NW = []
for i in range(len(area)):
    r_NW.append((area[i] / 3.14) ** 0.5)

All_D = []
for r in range(len(X)):
    D = []
    for n in range(len(X)):
        if r != n:
            distance = (((X[r] - X[n]) ** 2 + (Y[r] - Y[n]) ** 2) ** 0.5) - r_NW[r] - r_NW[n]
            D.append(distance)
    All_D.append(D)

# =========================
# Alternative center-to-center neighbor method (optional)
# =========================
coords = np.column_stack((X, Y))
k = 5
all_nearest_distances = []

for i in range(len(coords)):
    distances = np.sqrt(np.sum((coords - coords[i]) ** 2, axis=1))
    distances = np.sort(distances[distances > 0])[:k]
    all_nearest_distances.extend(distances)

mean_dist = np.mean(all_nearest_distances)
median_dist = np.median(all_nearest_distances)

print("Alternative center-to-center nearest-neighbor method")
print("Mean distance:", mean_dist)
print("Median distance:", median_dist)
print()

# =========================
# Nearest-neighbor extraction from edge-to-edge spacing
# =========================
D_min_all = []
n_NW_sum_mix = []
n_NW_sum = []

for i in range(len(All_D)):
    sort_All_D = sorted(All_D[i])

    NND = np.amin(All_D[i])
    D_min_all.append(NND)

    n_NW_1 = sort_All_D[0]
    n_NW_2 = sort_All_D[1]
    n_NW_3 = sort_All_D[2]
    n_NW_4 = sort_All_D[3]
    n_NW_5 = sort_All_D[4]

    n_NW_sum = [n_NW_1, n_NW_2, n_NW_3, n_NW_4, n_NW_5]

    n_NW_sum_mix.append(n_NW_1)
    n_NW_sum_mix.append(n_NW_2)
    n_NW_sum_mix.append(n_NW_3)
    n_NW_sum_mix.append(n_NW_4)
    n_NW_sum_mix.append(n_NW_5)

# =========================
# Diameter statistics and histogram
# =========================
STD_diameter = np.std(sort_diameter)
SE_diameter = STD_diameter / math.sqrt(len(sort_diameter))

density_NW = len(sort_diameter) / (image_size_x * image_size_y)
area_fraction = sum(area) / (image_size_x * image_size_y)

fig, ax = plt.subplots(figsize=(9 * cm, 8 * cm), layout="constrained")
plt.hist(sort_diameter, edgecolor="red", bins=10)
ax.tick_params(axis="both", which='both', direction="in", labelsize=15)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.xlabel('Diameter (nm)', size=15)
plt.ylabel('Counts', size=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig('diameter-histogram.png', dpi=300)
plt.show()

# =========================
# Minimum spacing statistics and histogram
# =========================
D_min_nm = [1000 * x_1 for x_1 in D_min_all]
ave_D_min_nm = np.sum(D_min_nm) / len(D_min_nm)

D_min_d = []
for n_1 in range(len(D_min_all) - 1):
    D_min_d.append((D_min_nm[n_1] - ave_D_min_nm) ** 2)

STD_min_d = (sum(D_min_d) / len(D_min_all)) ** 0.5
SE_min_d = STD_min_d / math.sqrt(len(D_min_all))

fig, ax = plt.subplots(figsize=(9 * cm, 8 * cm))
plt.hist(D_min_nm, edgecolor="red", bins='auto')
ax.tick_params(axis="both", which='both', direction="in", labelsize=15)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.xlabel('Distance (nm)', size=15)
plt.ylabel('Counts', size=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig('minimum-distance-histogram.png', dpi=300)
plt.show()

# =========================
# Pooled neighbor spacing statistics and histogram
# =========================
D_all_nm = [1000 * x_1 for x_1 in n_NW_sum_mix]
ave_D_all_nm = np.sum(D_all_nm) / len(D_all_nm)

D_all_d = []
for n_2 in range(len(D_all_nm) - 1):
    D_all_d.append((D_all_nm[n_2] - ave_D_all_nm) ** 2)

STD_all_d = (sum(D_all_d) / len(D_all_nm)) ** 0.5
SE_all_d = STD_all_d / math.sqrt(len(D_all_nm))

fig, ax = plt.subplots(figsize=(9 * cm, 8 * cm), layout="constrained")
plt.hist(D_all_nm, edgecolor="red", bins=12)
ax.tick_params(axis="both", which='both', direction="in", labelsize=15)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.xlabel('Spacing (nm)', size=15)
plt.ylabel('Count', size=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig('spacing-histogram.png', dpi=300)
plt.show()

# =========================
# Circularity statistics and histogram
# =========================
ave_circ = np.sum(circ) / len(circ)

circ_d = []
for n_3 in range(len(circ) - 1):
    circ_d.append((circ[n_3] - ave_circ) ** 2)

STD_circ = (sum(circ_d) / len(circ)) ** 0.5
SE_circ = STD_circ / math.sqrt(len(circ))

fig, ax = plt.subplots(figsize=(9 * cm, 8 * cm), layout="constrained")
plt.hist(circ, edgecolor="red", bins=10)
ax.tick_params(axis="both", which='both', direction="in", labelsize=15)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.xlabel('Circularity', size=15)
plt.ylabel('Count', size=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig('circularity-histogram.png', dpi=300)
plt.show()

# =========================
# Roundness statistics and histogram
# =========================
ave_round = np.sum(Round) / len(Round)

round_d = []
for n_4 in range(len(Round) - 1):
    round_d.append((Round[n_4] - ave_round) ** 2)

STD_round = (sum(round_d) / len(Round)) ** 0.5
SE_round = STD_round / math.sqrt(len(Round))

fig, ax = plt.subplots(figsize=(9 * cm, 8 * cm), layout="constrained")
plt.hist(Round, edgecolor="red", bins=10)
ax.tick_params(axis="both", which='both', direction="in", labelsize=15)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.xlabel('Roundness', size=15)
plt.ylabel('Count', size=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig('roundness-histogram.png', dpi=300)
plt.show()

# =========================
# Print summary
# =========================
print("Island diameter mean =", round(ave_diameter, 3), "nm")
print("Island diameter STD / SE =", round(STD_diameter, 3), "/", round(SE_diameter, 3), "nm")
print()

print("Minimum neighbor distance mean / STD / SE =", round(ave_D_min_nm, 3), "/", round(STD_min_d, 3), "/", round(SE_min_d, 3), "nm")
print("Neighbor distance mean / STD / SE =", round(ave_D_all_nm, 3), "/", round(STD_all_d, 3), "/", round(SE_all_d, 3), "nm")
print()

print("Circularity mean / STD / SE =", round(ave_circ, 3), "/", round(STD_circ, 3), "/", round(SE_circ, 3))
print("Roundness mean / STD / SE =", round(ave_round, 3), "/", round(STD_round, 3), "/", round(SE_round, 3))
print()

print("Island density =", round(density_NW, 3), "counts/µm2")
print("Nanowire area fraction =", round(area_fraction, 3), "%")
print("Image size =", image_size_x, "x", image_size_y, "µm²")
print("Area mean / STD / SE =", round(ave_area, 3), "/", round(STD_area, 3), "/", round(SE_area, 3), "nm²")

# =========================
# Export processed data
# =========================
Path_3 = directory + "diameter and circularity_" + filename
Path_4 = directory + "spacing_" + filename

combined_df = pd.DataFrame({
    'diameter_nm': sort_diameter,
    'circularity': circ
})

combined_df_2 = pd.DataFrame({
    'spacing_nm': D_all_nm
})

combined_df.to_csv(Path_3, index=False)
combined_df_2.to_csv(Path_4, index=False)

print()
print("Saved:", Path_3)
print("Saved:", Path_4)
