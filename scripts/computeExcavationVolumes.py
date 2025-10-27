import geopandas as gpd
import rasterio
import rasterio.mask
import numpy as np

# Charger MNT et emprises depuis le répertoire racine
mnt = rasterio.open("swissALTI3D.tif")
buildings = gpd.read_file("NouveauxBatiments.fg.json", layer="Building")

# Ajoute un tampon de 2m d'excavation
buildings['geometry'] = buildings.geometry.buffer(2)

volumes = []

for idx, b in buildings.iterrows():
    # Masque raster sous l'emprise
    geom = [b['geometry']]
    try:
        out_image, out_transform = rasterio.mask.mask(mnt, geom, crop=True, filled=True)
        data = out_image[0]
        # Ignore les valeurs nodata
        data = np.where(data == mnt.nodata, np.nan, data)
        zmin = np.nanmin(data)
        depth = 3
        excav = np.where(data < zmin + depth, depth - (data - zmin), 0)
        pixel_area = abs(mnt.res[0] * mnt.res[1])
        volume = round(np.nansum(excav) * pixel_area, 2)  # m² * m = m³, arrondi à 2 décimales
        print(f"Bâtiment {idx}: zmin={zmin:.2f}, volume={volume:.2f} m³")
        volumes.append(volume)
    except ValueError:
        print(f"Bâtiment {idx}: hors raster, ignoré.")
        volumes.append(np.nan)

if np.isnan(volumes).any():
    print("Attention: certains volumes sont NaN (bâtiment hors raster ou géométrie non recouvrante)")

total_volume = round(np.nansum(volumes), 2)
print("Volume total d'excavation (m³):", total_volume)
