# rhino3dm-gltf

Conversion de géodatabases GDB vers CityJSON avec validation val3dity et regroupement sémantique intégré.

## Activation de l'environnement

`. start.sh`

## Script principal

### Conversion GDB vers CityJSON avec regroupement sémantique intégré
```bash
# Conversion directe avec regroupement sémantique
python3 scripts/gdb_to_cityjson.py chemin/vers/fichier.gdb -o output/buildings_semantic.city.json

# Avec système de coordonnées spécifique
python3 scripts/gdb_to_cityjson.py swissbuildings3d.gdb --epsg 2056 -o output/swissbuildings3d.city.json

# Lister toutes les couches disponibles
python3 scripts/gdb_to_cityjson.py swissbuildings3d.gdb --list-layers

# Traiter des couches spécifiques
python3 scripts/gdb_to_cityjson.py swissbuildings3d.gdb --layers nom_couche1 nom_couche2 -o output/specific_layers.city.json
```

### Autres scripts de conversion
- `fix_closed_rings.py` - Corrige les rings fermés pour compatibilité val3dity (si nécessaire)
- `convert_to_v11.py` - Convertit CityJSON 2.0 vers 1.1
- `check_semantics.py` - Vérifie et analyse la sémantique des surfaces dans un fichier CityJSON

### Vérification de la sémantique
```bash
# Analyser la sémantique d'un fichier CityJSON
python3 scripts/check_semantics.py output/buildings_semantic.city.json
```

## Structure sémantique des bâtiments

Le fichier de sortie contient des bâtiments avec regroupement sémantique automatique :

- **Bâtiments** de type `Building` avec surfaces regroupées
- **Surfaces sémantiques** :
  - `WallSurface` : surfaces de mur (détection automatique par nom de couche)
  - `RoofSurface` : surfaces de toit (keywords: roof, toit, toiture, top)
  - `GroundSurface` : surfaces de sol (keywords: floor, sol, ground, base)

**Extraction automatique d'ID de bâtiment** :
- Recherche dans les attributs : `building_id`, `buildingid`, `id`, `name`, `gid`, `ewid`
- Génération automatique basée sur FID ou coordonnées si pas trouvé

**Sélection de couches** :
- `--list-layers` : Affiche toutes les couches disponibles avec leurs statistiques
- `--layers nom1 nom2` : Traite uniquement les couches spécifiées
- Par défaut : traite toutes les couches de la géodatabase
- Vérification automatique de l'existence des couches demandées

**Filtrage automatique par OBJEKTART** :
- Applique un filtre SQL au niveau de chaque couche : `OBJEKTART NOT IN ('Flugdach', 'Offenes Gebaeude')`
- Filtrage efficace une seule fois en amont, pas de vérification feature par feature
- Affichage des statistiques de filtrage dans les logs

**Classification automatique des surfaces** :
- Analyse des noms de couches pour déterminer le type de surface
- Support des keywords français et anglais
- Fallback intelligent sur WallSurface par défaut

## Accomplishments

- ✅ **Sélection de couches** : Support des noms de couches spécifiques avec --layers et listing avec --list-layers
- ✅ **Filtrage OBJEKTART optimisé** : Filtre SQL unique au niveau des couches pour exclure Flugdach et Offenes Gebaeude
- ✅ **Résolution des erreurs val3dity** : Correction des erreurs CONSECUTIVE_POINTS_SAME avec génération de rings ouverts
- ✅ **Regroupement sémantique intégré** : Classification automatique Wall/Roof/Floor par bâtiment directement lors de la conversion
- ✅ **Script unifié** : Fusion de `group_semantic_buildings.py` dans `gdb_to_cityjson.py` pour traitement direct
- ✅ **Classification intelligente** : Détection automatique du type de surface basée sur les noms de couches
- ✅ **Extraction d'ID robuste** : Recherche étendue d'identifiants de bâtiment dans les attributs

## TODO

- [x] Filtrer tous les bâtiments autres que Type IN ('Bâtiment ouvert','Toit flottant') en recherchant la correspondance allemande → **Terminé** : Filtrage automatique de Flugdach et Offenes Gebaeude


```python
python3 scripts/gdb_to_cityjson.py swissbuilding3d.gdb -o output/swissbuildings3d.city.json --layers Buildings --include-attrs UUID  STATUT --epsg 2056
```

Post-traitement des bâtiments swisstopo:

python scripts/postprocess_swissbuildings.py --input output/swissbuildings3d.city.json --output output/swissbuildings3d.city.json


fusion de fichiers avec cijo: (il ne semble apparemment pas possible de fusionner plus de 2 fichiers à la fois)

cjio output/Ronquoz21.city.json merge output/swissALTI3D.city.json save output/Ronquoz21Bis.city.json


cjio output/Ronquoz21.city.json merge output/swissALTI3D.city.json save output/Ronquoz21bis.city.json


## Conversion GeoTIFF TINRelief (CityJSON 2.0)

**Étape préalable recommandée : Simplification du raster (optionnelle)**

```bash
# Pour raster très haute résolution (< 1m) - Recommandé
gdalwarp -tr 1 1 -r average swissALTI3D.tif swissALTI3D_1m.tif

# Pour zones étendues - Réduction plus aggressive
gdalwarp -tr 2 2 -r average swissALTI3D.tif swissALTI3D_2m.tif

# Découpage d'une zone d'intérêt (recommandé)
gdalwarp -te 2592000 1119000 2594000 1121000 swissALTI3D.tif zone_interesse.tif
```

**Workflow optimal :**

```bash
# 1. Découpage + simplification (si nécessaire)
gdalwarp -te 2592000 1119000 2594000 1121000 -tr 1 1 -r average \
  swissALTI3D.tif swissALTI3D_zone_1m.tif

# 2. Conversion avec optimisation intelligente
python3 scripts/GeoTIFF2TINRelief.py -i swissALTI3D_zone_1m.tif \
  --optimization importance --tolerance 1.0 --materials \
  -o output/terrain_optimise.city.json

# 3. Alternative sans pré-traitement (utilise l'optimisation du script)
python3 scripts/GeoTIFF2TINRelief.py -i swissALTI3D.tif \
  --optimization adaptive --tolerance 2.0 --materials \
  -o output/terrain_direct.city.json
```

```bash
# Optimisation intelligente du terrain (recommandé pour raster non pré-traité)
python3 scripts/GeoTIFF2TINRelief.py -i mnt_5m.tif --optimization adaptive --tolerance 2.0

# Optimisation par importance géomorphologique (pics, vallées)
python3 scripts/GeoTIFF2TINRelief.py -i mnt_5m.tif --optimization importance --tolerance 1.5

# Simplification Douglas-Peucker 3D
python3 scripts/GeoTIFF2TINRelief.py -i mnt_5m.tif --optimization douglas-peucker --tolerance 1.0

# Avec point de référence personnalisé pour coordonnées locales
python3 scripts/GeoTIFF2TINRelief.py -i mnt_5m.tif --reference-point 2592000 1119000 500

# Configuration complète avec raster pré-simplifié
python3 scripts/GeoTIFF2TINRelief.py -i swissALTI3D_1m.tif \
  --optimization importance --tolerance 0.5 --materials \
  --object-id "terrain_final"

# Sans optimisation (raster déjà optimisé par gdalwarp)
python3 scripts/GeoTIFF2TINRelief.py -i swissALTI3D_1m.tif --materials

# Aide
python3 scripts/GeoTIFF2TINRelief.py --help
```

**Stratégies d'optimisation recommandées :**

| **Cas d'usage** | **Pré-traitement GDAL** | **Optimisation script** | **Résultat** |
|----------------|------------------------|------------------------|--------------|
| **Raster 0.25m, zone étendue** | `gdalwarp -tr 1 1` | `--optimization importance` | Optimal |
| **Raster 1m, zone modérée** | Optionnel | `--optimization adaptive` | Bon compromis |
| **Raster 5m, petit territoire** | Non nécessaire | `--optimization none` ou `-s 0.5` | Direct |
| **Performance maximale** | `gdalwarp -tr 2 2` | `--optimization douglas-peucker` | Ultra-léger |

**Avantages du pré-traitement GDAL :**
- ✅ Réduction massive à la source (plus efficace)
- ✅ Lissage naturel avec méthode `average`
- ✅ Découpage précis de zone d'intérêt
- ✅ Contrôle de la résolution finale

**Avantages de l'optimisation script :**
- ✅ Préservation intelligente des détails
- ✅ Algorithmes géomorphologiques
- ✅ Contrôle fin de la tolérance
- ✅ Pas de perte de données source

python3 scripts/GeoTIFF2TINRelief.py -i swissALTI3D_1m.tif -o output/swissALTI3D.city.json \
  --materials --object-id "swissALTI3D"

avec simplification:

python3 scripts/GeoTIFF2TINRelief.py -i swissALTI3D_1m.tif -o output/swissALTI3Dopt.city.json \
  --optimization importance --tolerance 1.0 --materials --object-id "swissALTI3D"

python3 scripts/GeoTIFF2TINRelief.py -i swissALTI3D_1m.tif -o output/swissALTI3Dopt.city.json \
  --optimization adaptive --tolerance 2.0 --materials --object-id "swissALTI3D"

extraction d'un subset du modèle

cjio output/LoD3_Railway.city.json subset --id UUID_d96effed-08fe-4f74-b134-05b194aa3cff save output/TINReliefExample.city.json

Point d'amélioration: créer un fichier CityJSONSeq (.jsonl) 
Commande finalement fonctionnelle:

python3 scripts/GeoTIFF2TINRelief.py --input swissALTI3D_1m.tif --output output/swissALTI3D.city.json --output-format cityjson


cjio output/Ronquoz21.city.json merge output/secteurs3d.city.json save output/Ronquoz21Restrictions.city.json


calcule des volumes 3d

python3 scripts/computeBoundingVolumes.py

comblement des dépressions:

python3 scripts/computeLandfillVolumes.py BiensFondsRemblaiements.gpkg -o output/LandfillVolumes.city.json


cjio output/swissALTI3D.city.json merge output/LandfillVolumes.city.json save output/LanfdillTerrain.city.json

Conversion JSON-FG
cityjson2jsonfg <input.city.json> <output.fg.json>

cjio output/NouveauxBatiments.city.json merge output/secteurs3d.city.json save output/BatimentsSecteurs.city.json