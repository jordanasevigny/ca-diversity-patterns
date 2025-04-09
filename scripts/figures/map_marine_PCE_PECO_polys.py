#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 11:12:40 2025

@author: jsevigny
"""

# Load libraries
import netCDF4 as nc
import pandas as pd
import zipfile
import matplotlib.pyplot as plt
import shapely as sh
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D
import re
import numpy as np

#%%

# LOAD ROMS GRID DATA
file_path_grid = '/home/jsevigny/grid_data/wc15n_grd.nc'
grid_ds = nc.Dataset(file_path_grid, 'r')
# mask_rho dataset
mask_rho = grid_ds.variables['mask_rho'][:]
# lat and lon rho from grid
lat_rho = grid_ds.variables['lat_rho'][:]
lon_rho = grid_ds.variables['lon_rho'][:]

# BL ORIGIN POLYGONS
points_ex_data = pd.read_csv('/home/jsevigny/git/CA_diversity_patterns/data/coast_polygons/wc15n_300km2_settling_polygons.txt', delimiter=',')
# Make the vertices into polygons
# Function that groups the lat and lon of the verteces then zips them into polygons
def create_polygon(group):
    return sh.geometry.Polygon(zip(group[' lon'], group[' lat']))
polygons = (
    points_ex_data.groupby('cell #')[[' lon', ' lat']]
    .apply(create_polygon))
polygon_df = polygons.reset_index(name='polygon')


#%%
# Connectivity polygon plot

fig, ax = plt.subplots(figsize=(12, 8))

# Plot land mass
plt.imshow(mask_rho, extent=[lon_rho.min(), lon_rho.max(), lat_rho.min(), lat_rho.max()], 
           cmap='pink', origin='lower', alpha=0.15)

# Loop through all rows in the polygon DataFrame
for idx, row in polygon_df.iterrows():
    polygon = row['polygon']
    poly_id = row['cell #']
    x, y = polygon.exterior.xy
    ax.plot(x, y, alpha=0.9, color='black')


# Adjust plot limits to fit all polygons
ax.autoscale()  # Automatically scale the axes to fit the polygons
ax.set_aspect('equal', adjustable='box')  # Keep aspect ratio equal for accurate representation

plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Bruce's Polygons")
plt.show()

#%%
# MARINe Sites

# Load the CSV file into a DataFrame
marine = pd.read_csv('/home/jsevigny/git/CA_diversity_patterns/data/processed_data/biodiversity/marine_species_20241025_dates_20241028_merged.csv')

# Remove duplicates based on latitude and longitude
unique_coords_marine = marine[['latitude', 'longitude']].drop_duplicates()

# Plot the unique coordinates (MARINe sites)
plt.figure(figsize=(10, 6))
plt.scatter(unique_coords_marine['longitude'], unique_coords_marine['latitude'], c='#00BFC4', marker='*')

# Add labels and title
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('MARINe Sites')

# Show the plot
plt.show()


# MARINE Sites Extended
# Specify the ZIP file and the CSV file inside it

zf = zipfile.ZipFile('/home/jsevigny/git/CA_diversity_patterns/data/marine_data/gbif_marine_biodiversity_20250409.zip') 
marine_ex = pd.read_csv(zf.open('0014067-250402121839773.csv'), sep='\t')

# Remove duplicates based on latitude and longitude
unique_coords_marine_ex = marine_ex[['decimalLatitude', 'decimalLongitude']].drop_duplicates()

# Plot the unique coordinates (MARINe sites)
plt.figure(figsize=(10, 6))
plt.scatter(unique_coords_marine_ex['decimalLongitude'], unique_coords_marine_ex['decimalLatitude'], c='#00BFC4', marker='*')

# Add labels and title
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('MARINe Sites - Extended')

# Show the plot
plt.show()

#%%
# PCE (PC-eDNA) Sites

# Load the CSV file into a DataFrame
PCE = pd.read_csv('/home/jsevigny/git/CA_diversity_patterns/data/edna_jv_data/summer_2024/Summer_Exped_And_Central_Coast_All_Metadata_2024.csv')

# Drop odd sample sites (inland)
PCE = PCE.dropna(subset=['WP'])

# Remove duplicates based on latitude and longitude
unique_coords_PCE = PCE[['Latitude', 'Longitude']].drop_duplicates()

# Plot the unique coordinates (MARINe sites)
plt.figure(figsize=(10, 6))
plt.scatter(unique_coords_PCE['Longitude'], unique_coords_PCE['Latitude'], c='#F8766D', marker='*')

# Add labels and title
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('PC-eDNA Sites')

# Show the plot
plt.show()

#%%
# PECO Sites

# Load the CSV file into a DataFrame
peco = pd.read_csv('/home/jsevigny/git/CA_diversity_patterns/data/peco_data/peco_sites_manual.csv')

# Remove duplicates based on latitude and longitude
unique_coords_peco = peco[['Latitude', 'Longitude']].drop_duplicates()

# Plot the unique coordinates (MARINe sites)
plt.figure(figsize=(10, 6))
plt.scatter(unique_coords_peco['Longitude'], unique_coords_peco['Latitude'], c='#C77CFF', marker='*')

# Add labels and title
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('PECO Sites')

# Show the plot
plt.show()

#%%
# Map with combined polygons and sites

fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# Add coastlines to the map
#ax.coastlines(resolution='110m', color='black', linewidth=0.2)
# Optional: Add additional geographical features (e.g., land, rivers)
ax.add_feature(cfeature.LAND, edgecolor='black', linewidth=0.2)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.1, edgecolor='gray')
ax.add_feature(cfeature.STATES, linestyle='-', linewidth=0.1, edgecolor='gray')

# Add sample sites
ax.scatter(unique_coords_marine_ex['decimalLongitude'], unique_coords_marine_ex['decimalLatitude'], c='#00BFC4', marker='*')
ax.scatter(unique_coords_peco['Longitude'], unique_coords_peco['Latitude'], c='#C77CFF', marker='*')
ax.scatter(unique_coords_PCE['Longitude'], unique_coords_PCE['Latitude'], c='#F8766D', marker='*')

# Add settlement polygons
# Loop through all rows in the polygon DataFrame
for idx, row in polygon_df.iterrows():
    polygon = row['polygon']
    
    # Limit settlement polygons to max lat of california
    if polygon.bounds[3] < 42.0:  # bounds = (minx, miny, maxx, maxy)
        x, y = polygon.exterior.xy
        ax.plot(x, y, alpha=0.9, color='black', linewidth=0.5)
    
    # x, y = polygon.exterior.xy
    # ax.plot(x, y, alpha=0.9, color='black', linewidth=0.5)

# Adjust plot limits to fit all polygons
ax.autoscale()  # Automatically scale the axes to fit the polygons
ax.set_aspect('equal', adjustable='box')  # Keep aspect ratio equal for accurate representation

# Custom legend handles
legend_elements = [
    Line2D([0], [0], marker='*', color='w', label='MARINe Sites', markerfacecolor='#00BFC4', markersize=10),
    Line2D([0], [0], marker='*', color='w', label='PCE Sites', markerfacecolor='#F8766D', markersize=10),
    Line2D([0], [0], marker='*', color='w', label='PECO Sites', markerfacecolor='#C77CFF', markersize=10),
    Line2D([0], [0], color='black', lw=1, label='Settlement Cells')
]

# Add the legend to the plot
ax.legend(handles=legend_elements, loc='lower left')  # or another location like 'upper right'

# Add gridlines for latitude and longitude
ax.grid(True, which='both', linestyle=':', color='gray', linewidth=0.5)

# Customize latitude and longitude ticks and labels
ax.set_xticks(range(-138, -110, 3))  # Longitude range (modify as needed)
ax.set_yticks(range(25, 60, 3))  # Latitude range (modify as needed)

# Optional: Customize tick labels
ax.set_xticklabels([f"{i}°" for i in range(-138, -110, 3)])  # Longitude labels
ax.set_yticklabels([f"{i}°" for i in range(25, 60, 3)])  # Latitude labels

plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Settlement Cells & Sample Sites")

plt.savefig("/home/jsevigny/git/CA_diversity_patterns/figures/biodiversity/map_marine_PCE_PECO_polys.pdf", bbox_inches='tight')
plt.show()

