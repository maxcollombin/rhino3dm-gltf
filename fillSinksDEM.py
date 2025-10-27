import whitebox
import rasterio
import numpy as np
import os

# === PARAMÈTRES ===
input_dem = os.path.abspath("mntR21.tif")
filled_dem_path = os.path.abspath("mntR21_comble.tif")
diff_path = os.path.abspath("mntR21_diff.tif")

# === CONFIGURATION WHITEBOX ===
wbt = whitebox.WhiteboxTools()
wbt.verbose = True  # Active l'affichage des messages de whitebox dans la console

print("WhiteboxTools version:", wbt.version())
print("WhiteboxTools exe path:", wbt.exe_path)

if not os.path.exists(input_dem):
    raise FileNotFoundError(f"Le fichier DEM {input_dem} n'existe pas !")

# === COMBLER LES DÉPRESSIONS ===
print(">> Exécution de FillDepressions...")
success = wbt.fill_depressions(
    dem=input_dem,
    output=filled_dem_path
)
print("Résultat WhiteboxTools:", success)

if not os.path.exists(filled_dem_path):
    raise RuntimeError(f"Le fichier {filled_dem_path} n'a pas été créé. "
                       "Vérifie les messages ci-dessus pour comprendre pourquoi.")

# === CALCUL DE LA DIFFÉRENCE ===
with rasterio.open(input_dem) as src:
    dem_array = src.read(1)
    profile = src.profile

with rasterio.open(filled_dem_path) as src_filled:
    filled_array = src_filled.read(1)

diff = np.array(filled_array - dem_array, dtype=np.float32)

# === DÉTECTION DES CELLULES AVEC DIFFÉRENCE >= 0.5m ===
threshold = 0.5
significant_diff_mask = diff >= threshold
significant_diff_count = np.sum(significant_diff_mask)
total_cells = diff.size
percentage = (significant_diff_count / total_cells) * 100

print(f"\n=== ANALYSE DES DIFFÉRENCES ===")
print(f"Seuil de détection : {threshold}m")
print(f"Cellules avec différence >= {threshold}m : {significant_diff_count}")
print(f"Total de cellules : {total_cells}")
print(f"Pourcentage de cellules affectées : {percentage:.2f}%")

if significant_diff_count > 0:
    print(f"Différence minimale détectée : {np.min(diff[significant_diff_mask]):.3f}m")
    print(f"Différence maximale détectée : {np.max(diff[significant_diff_mask]):.3f}m")
    print(f"Différence moyenne (cellules >= {threshold}m) : {np.mean(diff[significant_diff_mask]):.3f}m")

# === CRÉER UN MASQUE BINAIRE DES ZONES SIGNIFICATIVES ===
mask_path = os.path.abspath("mntR21_mask_05m.tif")
mask_array = significant_diff_mask.astype(np.uint8)

profile_mask = profile.copy()
profile_mask.update(dtype=rasterio.uint8, nodata=0)
with rasterio.open(mask_path, "w", **profile_mask) as dst:
    dst.write(mask_array, 1)

profile.update(dtype=rasterio.float32, nodata=None)
with rasterio.open(diff_path, "w", **profile) as dst:
    dst.write(diff, 1)

print(f"\nDEM comblé sauvegardé : {filled_dem_path}")
print(f"Raster de différence sauvegardé : {diff_path}")
print(f"Masque des zones >= {threshold}m sauvegardé : {mask_path}")
