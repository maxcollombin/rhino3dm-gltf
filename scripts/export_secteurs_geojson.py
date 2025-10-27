import os
import json
import rhino3dm
import numpy as np

# ------------------------------
# Configuration
# ------------------------------
input_file = "./R21_3d.3dm"
output_geojson = "./output/secteurs.geojson"

# Secteurs à extraire (basé sur votre liste)
secteurs_to_extract = [
    "S-SECTEURS",
    "R21.1", "R21.2", "R21.3", "R21.4", "R21.5", "R21.6", "R21.7", "R21.8", "R21.9",
    "R21.10", "R21.11", "R21.12", "R21.13", "R21.14", "R21.15", "R21.16", "R21.17", 
    "R21.18", "R21.19", "R21.20", "R21.21"
]

# ------------------------------
# Chargement du modèle
# ------------------------------
model = rhino3dm.File3dm.Read(input_file)
if not model:
    raise ValueError(f"Erreur : Impossible de charger {input_file}")
print(f"Modèle chargé : {input_file}")

# ------------------------------
# Configuration transformation (même que CityJSON)
# ------------------------------
SCALE = [1.0, 1.0, 1.0]
TRANSLATE = [2592980.685, 1119281.703, 0.0]

# ------------------------------
# Fonctions utilitaires
# ------------------------------
def rhino_point_to_coords(point):
    """Convertit un point Rhino3D en coordonnées [x, y] avec transformation"""
    # Appliquer la transformation comme pour CityJSON
    x_transformed = point.X * SCALE[0] + TRANSLATE[0]
    y_transformed = point.Y * SCALE[1] + TRANSLATE[1]
    return [x_transformed, y_transformed]

def extract_curve_coordinates(geometry):
    """Extrait les coordonnées d'une courbe Rhino3D"""
    coordinates = []
    
    try:
        # Essayer d'obtenir une polyligne
        if hasattr(geometry, 'TryGetPolyline'):
            try:
                result = geometry.TryGetPolyline()
                if isinstance(result, tuple) and len(result) >= 2:
                    success, polyline = result[0], result[1]
                    if success and polyline:
                        coordinates = [rhino_point_to_coords(pt) for pt in polyline]
                        # Fermer le polygone si nécessaire
                        if len(coordinates) > 2 and coordinates[0] != coordinates[-1]:
                            coordinates.append(coordinates[0])
                        return [coordinates]  # GeoJSON Polygon format
                elif hasattr(result, '__len__') and len(result) > 0:
                    # Résultat direct de polyline
                    coordinates = [rhino_point_to_coords(pt) for pt in result]
                    if len(coordinates) > 2 and coordinates[0] != coordinates[-1]:
                        coordinates.append(coordinates[0])
                    return [coordinates]
            except:
                pass
        
        # Pour les courbes complexes, échantillonner des points
        if hasattr(geometry, 'DivideByCount'):
            try:
                # Diviser la courbe en segments
                points = geometry.DivideByCount(50, True)  # 50 points
                if points:
                    coordinates = [rhino_point_to_coords(pt) for pt in points]
                    # Fermer le polygone si nécessaire
                    if len(coordinates) > 2 and coordinates[0] != coordinates[-1]:
                        coordinates.append(coordinates[0])
                    return [coordinates]  # GeoJSON Polygon format
            except:
                pass
        
        # Fallback: essayer d'obtenir des points de contrôle
        if hasattr(geometry, 'PointAtStart') and hasattr(geometry, 'PointAtEnd'):
            try:
                start = geometry.PointAtStart
                end = geometry.PointAtEnd
                # Échantillonner la courbe
                t_values = [i/20.0 for i in range(21)]  # 21 points
                points = []
                for t in t_values:
                    pt = geometry.PointAt(t)
                    points.append(rhino_point_to_coords(pt))
                
                if len(points) > 2:
                    # Fermer si nécessaire
                    if points[0] != points[-1]:
                        points.append(points[0])
                    return [points]
            except:
                pass
    
    except Exception as e:
        print(f"  Erreur extraction courbe: {e}")
    
    return None

def convert_geometry_to_geojson(geometry, layer_name):
    """Convertit une géométrie Rhino3D en feature GeoJSON"""
    feature = {
        "type": "Feature",
        "properties": {
            "layer": layer_name,
            "type": "Secteur"
        },
        "geometry": None
    }
    
    # Traitement selon le type de géométrie
    if isinstance(geometry, rhino3dm.Curve):
        coordinates = extract_curve_coordinates(geometry)
        if coordinates:
            feature["geometry"] = {
                "type": "Polygon",
                "coordinates": coordinates
            }
    
    elif isinstance(geometry, rhino3dm.Point):
        feature["geometry"] = {
            "type": "Point",
            "coordinates": rhino_point_to_coords(geometry.Location)
        }
    
    elif isinstance(geometry, rhino3dm.Mesh):
        # Pour les mesh, créer un polygone à partir des vertices
        vertices = []
        for vertex in geometry.Vertices:
            vertices.append(rhino_point_to_coords(vertex))
        
        if len(vertices) >= 3:
            # Fermer le polygone
            if vertices[0] != vertices[-1]:
                vertices.append(vertices[0])
            
            feature["geometry"] = {
                "type": "Polygon",
                "coordinates": [vertices]
            }
    
    elif isinstance(geometry, rhino3dm.Brep):
        # Pour les Brep, essayer d'extraire les contours
        try:
            # Obtenir les courbes de bordure
            edges = geometry.Edges
            if edges:
                all_coords = []
                for edge in edges:
                    curve = edge.DuplicateCurve()
                    if curve:
                        coords = extract_curve_coordinates(curve)
                        if coords:
                            all_coords.extend(coords[0])
                
                if all_coords:
                    # Fermer le polygone
                    if all_coords[0] != all_coords[-1]:
                        all_coords.append(all_coords[0])
                    
                    feature["geometry"] = {
                        "type": "Polygon",
                        "coordinates": [all_coords]
                    }
        except:
            pass
    
    return feature if feature["geometry"] else None

# ------------------------------
# Extraction des secteurs
# ------------------------------
features = []
secteurs_trouves = set()

print("🔍 Recherche des secteurs...")

for obj in model.Objects:
    layer = model.Layers[obj.Attributes.LayerIndex]
    layer_name = layer.Name
    full_path = layer.FullPath
    
    # Vérifier si ce calque correspond à un secteur recherché
    secteur_match = None
    for secteur in secteurs_to_extract:
        if layer_name == secteur or full_path.endswith(secteur) or secteur in full_path:
            secteur_match = secteur
            break
    
    if secteur_match:
        secteurs_trouves.add(secteur_match)
        
        # Convertir la géométrie en feature GeoJSON
        feature = convert_geometry_to_geojson(obj.Geometry, secteur_match)
        
        if feature:
            # Ajouter des propriétés supplémentaires
            feature["properties"].update({
                "secteur_id": secteur_match,
                "full_path": full_path,
                "layer_name": layer_name
            })
            
            features.append(feature)
            print(f"   {secteur_match}: {len(features)} géométries trouvées")

# ------------------------------
# Création du GeoJSON
# ------------------------------
geojson = {
    "type": "FeatureCollection",
    "crs": {
        "type": "name",
        "properties": {
            "name": "urn:ogc:def:crs:EPSG::2056"
        }
    },
    "features": features,
    "metadata": {
        "title": "Secteurs R21 - Extraction depuis Rhino3dm",
        "description": "Secteurs extraits du modèle R21_3d.3dm avec transformation EPSG:2056",
        "source": input_file,
        "coordinate_system": "EPSG:2056 - CH1903+ / LV95",
        "transformation": {
            "scale": SCALE,
            "translate": TRANSLATE,
            "description": "Coordonnées transformées depuis le système local Rhino vers EPSG:2056"
        },
        "secteurs_recherches": secteurs_to_extract,
        "secteurs_trouves": list(secteurs_trouves),
        "secteurs_manquants": list(set(secteurs_to_extract) - secteurs_trouves),
        "total_features": len(features)
    }
}

# ------------------------------
# Sauvegarde
# ------------------------------
os.makedirs(os.path.dirname(output_geojson), exist_ok=True)

with open(output_geojson, "w", encoding="utf-8") as f:
    json.dump(geojson, f, indent=2, ensure_ascii=False)

# ------------------------------
# Rapport final
# ------------------------------
print(f"\ Extraction terminée !")
print(f" Fichier créé : {output_geojson}")
print(f" Système de coordonnées : EPSG:2056 (CH1903+ / LV95)")
print(f" Transformation appliquée : translate={TRANSLATE}")
print(f" Statistiques :")
print(f"  • Features créées : {len(features)}")
print(f"  • Secteurs trouvés : {len(secteurs_trouves)}/{len(secteurs_to_extract)}")
print(f"  • Secteurs trouvés : {', '.join(sorted(secteurs_trouves))}")

if secteurs_trouves != set(secteurs_to_extract):
    manquants = set(secteurs_to_extract) - secteurs_trouves
    print(f"   Secteurs non trouvés : {', '.join(sorted(manquants))}")

print(f"\n GeoJSON créé avec succès avec coordonnées EPSG:2056 !")
