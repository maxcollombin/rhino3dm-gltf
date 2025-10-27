#!/usr/bin/env python3
"""
Génère des volumes 3D avec comblements (fills) pour les parcelles du fichier BienFondsRemblaiements.gpkg
Inspiré de fillSinksDEM.py mais adapté pour la création de maillages 3D volumétriques.
"""

import geopandas as gpd
import rasterio
import rasterio.features
import numpy as np
from shapely.geometry import mapping, Polygon, MultiPolygon
import warnings
import json
import uuid
from pathlib import Path
try:
    from scipy.spatial import Delaunay
    from scipy.ndimage import gaussian_filter
    from scipy.interpolate import griddata
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("⚠️ SciPy non disponible - utilisation d'algorithmes simplifiés")

warnings.filterwarnings("ignore", category=UserWarning)

class FilledVolumeGenerator:
    """
    Générateur de volumes 3D comblés pour les parcelles
    """
    
    def __init__(self, parcelles_path, mnt_path, output_path, min_fill_height=0.5):
        """
        Initialise le générateur de volumes comblés
        
        Args:
            parcelles_path: Chemin vers BienFondsRemblaiements.gpkg
            mnt_path: Chemin vers le MNT (swissALTI3D.tif)
            output_path: Chemin de sortie pour le CityJSON
            min_fill_height: Hauteur minimale de comblement en mètres (défaut: 0.5m)
        """
        self.parcelles_path = Path(parcelles_path)
        self.mnt_path = Path(mnt_path)
        self.output_path = Path(output_path)
        self.min_fill_height = min_fill_height
        
        # Variables de traitement
        self.parcelles = None
        self.mnt = None
        self.vertices = []
        self.vertex_index_map = {}
        self.city_objects = {}
        
        # Point de référence pour transformation locale
        self.ref_point = (2592980.685, 1119281.703, 483.8)
        
        # Statistiques
        self.stats = {
            'parcelles_processed': 0,
            'volumes_created': 0,
            'total_fill_volume': 0.0,
            'vertices_total': 0,
            'avg_fill_height': 0.0,
            'max_fill_height': 0.0,
            'min_calculated_height': float('inf')
        }
    
    def load_data(self):
        """Charge les données parcelles et MNT"""
        print(f"🔄 Chargement des données...")
        
        # Charger les parcelles
        print(f"📁 Chargement des parcelles: {self.parcelles_path}")
        self.parcelles = gpd.read_file(self.parcelles_path)
        print(f"📊 {len(self.parcelles)} parcelles chargées")
        
        # Afficher les colonnes disponibles
        print(f"📋 Colonnes disponibles: {list(self.parcelles.columns)}")
        print(f"📍 CRS parcelles: {self.parcelles.crs}")
        
        # Charger le MNT
        print(f"🗻 Chargement du MNT: {self.mnt_path}")
        self.mnt = rasterio.open(self.mnt_path)
        print(f"📐 Résolution MNT: {self.mnt.res}")
        print(f"📍 CRS MNT: {self.mnt.crs}")
        
        # Reprojeter les parcelles si nécessaire
        if self.parcelles.crs != self.mnt.crs:
            print(f"🔄 Reprojection des parcelles vers {self.mnt.crs}")
            self.parcelles = self.parcelles.to_crs(self.mnt.crs)
        
        return True
    
    def add_vertex(self, x, y, z):
        """Ajoute un vertex avec déduplication"""
        ref_x, ref_y, ref_z = self.ref_point
        
        # Transformation en coordonnées locales
        x_local = x - ref_x
        y_local = y - ref_y
        z_local = z - ref_z
        
        # Clé pour éviter les doublons
        key = (round(x_local, 3), round(y_local, 3), round(z_local, 3))
        
        if key in self.vertex_index_map:
            return self.vertex_index_map[key]
        
        idx = len(self.vertices)
        self.vertices.append([x_local, y_local, z_local])
        self.vertex_index_map[key] = idx
        return idx
    
    def sample_mnt_for_polygon(self, polygon, resolution=1.0):
        """
        Échantillonne le MNT dans un polygone avec une résolution donnée
        
        Args:
            polygon: Polygone Shapely
            resolution: Résolution d'échantillonnage en mètres
            
        Returns:
            Liste de points (x, y, z)
        """
        # Obtenir les limites du polygone
        minx, miny, maxx, maxy = polygon.bounds
        
        # Créer une grille de points
        x_coords = np.arange(minx, maxx, resolution)
        y_coords = np.arange(miny, maxy, resolution)
        
        points = []
        
        for x in x_coords:
            for y in y_coords:
                from shapely.geometry import Point
                point = Point(x, y)
                
                # Vérifier si le point est dans le polygone
                if polygon.contains(point):
                    # Interpoler l'altitude depuis le MNT
                    try:
                        row, col = self.mnt.index(x, y)
                        if 0 <= row < self.mnt.height and 0 <= col < self.mnt.width:
                            altitude = self.mnt.read(1)[row, col]
                            if not np.isnan(altitude) and altitude != self.mnt.nodata:
                                points.append((x, y, float(altitude)))
                    except:
                        continue
        
        return points
    
    def calculate_dynamic_fill_height(self, points, polygon):
        """
        Calcule la hauteur de comblement dynamique basée sur l'analyse du terrain
        
        Args:
            points: Liste de points (x, y, z)
            polygon: Polygone de la parcelle
            
        Returns:
            Dictionnaire avec les informations de comblement calculées
        """
        if len(points) < 4:
            return None
        
        # Convertir en arrays numpy
        coords = np.array(points)
        x_coords = coords[:, 0]
        y_coords = coords[:, 1]
        z_coords = coords[:, 2]
        
        # Statistiques de base du terrain
        min_altitude = np.min(z_coords)
        max_altitude = np.max(z_coords)
        mean_altitude = np.mean(z_coords)
        std_altitude = np.std(z_coords)
        
        # Appliquer un filtre gaussien pour lisser et identifier les tendances
        # Créer une grille régulière pour l'interpolation
        minx, maxx = x_coords.min(), x_coords.max()
        miny, maxy = y_coords.min(), y_coords.max()
        
        # Grille pour l'analyse
        grid_resolution = min(50, len(points) // 4)  # Résolution adaptative
        if grid_resolution < 10:
            grid_resolution = 10
        
        xi = np.linspace(minx, maxx, grid_resolution)
        yi = np.linspace(miny, maxy, grid_resolution)
        XI, YI = np.meshgrid(xi, yi)
        
        # Interpolation des altitudes sur la grille
        if SCIPY_AVAILABLE:
            ZI = griddata((x_coords, y_coords), z_coords, (XI, YI), method='cubic', fill_value=np.nan)
            # Appliquer un filtre gaussien pour identifier les tendances de surface
            ZI_smooth = gaussian_filter(ZI, sigma=2.0)
        else:
            # Version simplifiée sans SciPy
            ZI = np.full_like(XI, mean_altitude)
            ZI_smooth = ZI.copy()
        
        # Calculer les dépressions (différence entre surface lissée et originale)
        depressions = ZI_smooth - ZI
        depressions = np.where(depressions > 0, depressions, 0)  # Seules les dépressions positives
        
        # Calculer les statistiques des dépressions
        total_depression_volume = np.nansum(depressions) * (xi[1] - xi[0]) * (yi[1] - yi[0])
        max_depression = np.nanmax(depressions) if np.any(~np.isnan(depressions)) else 0
        mean_depression = np.nanmean(depressions[depressions > 0]) if np.any(depressions > 0) else 0
        
        # CALCUL DYNAMIQUE DE LA HAUTEUR DE COMBLEMENT
        # Stratégie : combiner plusieurs facteurs pour déterminer la hauteur optimale
        
        # 1. Hauteur basée sur les dépressions détectées
        depression_fill_height = max_depression * 1.2  # 20% de marge
        
        # 2. Hauteur basée sur la variabilité du terrain
        variability_fill_height = std_altitude * 0.8  # Proportionnel à la variabilité
        
        # 3. Hauteur basée sur la pente générale
        # Calculer la pente moyenne du terrain
        slope_magnitude = 0
        if len(points) > 10:
            # Régression linéaire simple pour estimer la pente
            A = np.column_stack([x_coords, y_coords, np.ones(len(x_coords))])
            try:
                coeffs, _, _, _ = np.linalg.lstsq(A, z_coords, rcond=None)
                slope_magnitude = np.sqrt(coeffs[0]**2 + coeffs[1]**2)
                slope_fill_height = slope_magnitude * 2.0  # Facteur basé sur la pente
            except:
                slope_fill_height = 0
        else:
            slope_fill_height = 0
        
        # 4. Hauteur basée sur la surface de la parcelle (plus grande parcelle = plus de comblement)
        area_factor = min(polygon.area / 1000.0, 5.0)  # Maximum 5m pour très grandes parcelles
        area_fill_height = area_factor * 0.5
        
        # Combiner les différents facteurs
        calculated_fill_height = max(
            depression_fill_height,
            variability_fill_height,
            slope_fill_height,
            area_fill_height,
            self.min_fill_height  # Garantir une hauteur minimale
        )
        
        # Limiter la hauteur maximale pour éviter des valeurs aberrantes
        max_allowed_fill = min(10.0, max_altitude - min_altitude + 2.0)
        calculated_fill_height = min(calculated_fill_height, max_allowed_fill)
        
        # Altitude de comblement
        base_altitude = mean_altitude  # Utiliser l'altitude moyenne comme référence
        fill_altitude = base_altitude + calculated_fill_height
        
        return {
            'base_altitude': base_altitude,
            'fill_altitude': fill_altitude,
            'calculated_fill_height': calculated_fill_height,
            'max_depression': max_depression,
            'mean_depression': mean_depression,
            'natural_depression_volume': total_depression_volume,
            'fill_volume': polygon.area * calculated_fill_height,
            'terrain_stats': {
                'min_altitude': min_altitude,
                'max_altitude': max_altitude,
                'mean_altitude': mean_altitude,
                'std_altitude': std_altitude,
                'slope_magnitude': slope_magnitude
            },
            'fill_factors': {
                'depression_factor': depression_fill_height,
                'variability_factor': variability_fill_height,
                'slope_factor': slope_fill_height,
                'area_factor': area_fill_height
            },
            'grid_x': XI,
            'grid_y': YI,
            'grid_z_original': ZI,
            'grid_z_smooth': ZI_smooth,
            'depressions': depressions
        }
    
    def create_filled_volume_geometry(self, polygon, fill_data):
        """
        Crée la géométrie 3D du volume comblé
        
        Args:
            polygon: Polygone de la parcelle
            fill_data: Données de comblement
            
        Returns:
            Dictionnaire avec vertices et boundaries pour CityJSON
        """
        if isinstance(polygon, MultiPolygon):
            # Prendre le plus grand polygone
            polygon = max(polygon.geoms, key=lambda p: p.area)
        
        if not hasattr(polygon, 'exterior'):
            return None
        
        # Extraire les coordonnées du contour
        exterior_coords = list(polygon.exterior.coords[:-1])  # Enlever le dernier point dupliqué
        
        if len(exterior_coords) < 3:
            return None
        
        # Créer les vertices pour le sol et le plafond
        bottom_indices = []
        top_indices = []
        
        base_altitude = fill_data['base_altitude']
        fill_altitude = fill_data['fill_altitude']
        
        for x, y in exterior_coords:
            # Interpoler l'altitude du sol depuis le MNT
            try:
                row, col = self.mnt.index(x, y)
                if 0 <= row < self.mnt.height and 0 <= col < self.mnt.width:
                    ground_altitude = self.mnt.read(1)[row, col]
                    if np.isnan(ground_altitude) or ground_altitude == self.mnt.nodata:
                        ground_altitude = base_altitude
                else:
                    ground_altitude = base_altitude
            except:
                ground_altitude = base_altitude
            
            # Ajouter les vertices du sol et du plafond
            bottom_idx = self.add_vertex(x, y, ground_altitude)
            top_idx = self.add_vertex(x, y, fill_altitude)
            
            bottom_indices.append(bottom_idx)
            top_indices.append(top_idx)
        
        # Construire les faces du solide
        boundaries = []
        
        # Face du sol (orientée vers le bas)
        bottom_face = list(reversed(bottom_indices))
        boundaries.append([bottom_face])
        
        # Face du plafond (orientée vers le haut)
        top_face = top_indices.copy()
        boundaries.append([top_face])
        
        # Faces latérales
        n_points = len(bottom_indices)
        for i in range(n_points):
            next_i = (i + 1) % n_points
            
            # Face latérale (quadrilatère)
            wall_face = [
                bottom_indices[i],
                bottom_indices[next_i],
                top_indices[next_i],
                top_indices[i]
            ]
            boundaries.append([wall_face])
        
        return {
            'boundaries': boundaries,
            'fill_data': fill_data
        }
    
    def create_cityjson_object(self, parcelle_row, geometry_3d):
        """
        Crée un objet CityJSON pour une parcelle comblée
        
        Args:
            parcelle_row: Ligne de la parcelle depuis le GeoDataFrame
            geometry_3d: Géométrie 3D générée
            
        Returns:
            Tuple (obj_id, city_object)
        """
        obj_id = str(uuid.uuid4())
        
        # Extraire les attributs de la parcelle
        attributes = {}
        for col in parcelle_row.index:
            if col != 'geometry' and parcelle_row[col] is not None:
                value = parcelle_row[col]
                # Convertir les types NumPy en types Python natifs pour la sérialisation JSON
                if isinstance(value, (np.integer, np.int64, np.int32)):
                    attributes[col] = int(value)
                elif isinstance(value, (np.floating, np.float64, np.float32)):
                    attributes[col] = float(value)
                else:
                    attributes[col] = value
        
        # Ajouter les informations de comblement
        fill_data = geometry_3d['fill_data']
        attributes.update({
            'fillHeight': round(fill_data['calculated_fill_height'], 2),
            'fillAltitude': round(fill_data['fill_altitude'], 2),
            'fillVolume': round(fill_data['fill_volume'], 2)
        })
        
        city_object = {
            "type": "GenericCityObject",
            "attributes": attributes,
            "geometry": [{
                "type": "Solid",
                "lod": "3.0",
                "boundaries": [geometry_3d['boundaries']]  # Un shell
            }]
        }
        
        return obj_id, city_object
    
    def process_parcelles(self):
        """Traite toutes les parcelles pour créer les volumes comblés"""
        print(f"\n Traitement des parcelles...")
        
        fill_heights = []  # Pour calculer les statistiques
        
        for idx, parcelle in self.parcelles.iterrows():
            self.stats['parcelles_processed'] += 1
            
            # Identifier la parcelle
            parcelle_id = parcelle.get('id', parcelle.get('ID', parcelle.get('fid', f'parcelle_{idx}')))
            print(f"📍 Parcelle {parcelle_id} ({self.stats['parcelles_processed']}/{len(self.parcelles)})")
            
            geometry = parcelle.geometry
            if not geometry or geometry.is_empty:
                print(f"   ⚠️ Géométrie vide, ignorée")
                continue
            
            # Échantillonner le MNT dans la parcelle
            print(f"   🗻 Échantillonnage du MNT (résolution: 0.5m)")
            points = self.sample_mnt_for_polygon(geometry, resolution=0.5)
            
            if len(points) < 4:
                print(f"   ⚠️ Pas assez de points MNT ({len(points)}), ignorée")
                continue
            
            print(f"   📊 {len(points)} points échantillonnés")
            
            # Calculer la hauteur de comblement dynamique
            print(f"   🧮 Calcul dynamique de la hauteur de comblement")
            fill_data = self.calculate_dynamic_fill_height(points, geometry)
            
            if not fill_data:
                print(f"   ❌ Impossible de calculer le comblement")
                continue
            
            calculated_height = fill_data['calculated_fill_height']
            fill_heights.append(calculated_height)
            
            print(f"   📏 Altitude de base: {fill_data['base_altitude']:.2f}m")
            print(f"   🏗️ Hauteur de comblement calculée: {calculated_height:.2f}m")
            print(f"   🏁 Altitude de comblement: {fill_data['fill_altitude']:.2f}m")
            print(f"   📦 Volume de comblement: {fill_data['fill_volume']:.1f} m³")
            print(f"   🕳️ Volume naturel des dépressions: {fill_data['natural_depression_volume']:.1f} m³")
            print(f"   📉 Dépression maximale: {fill_data['max_depression']:.2f}m")
            print(f"   📊 Variabilité du terrain: ±{fill_data['terrain_stats']['std_altitude']:.2f}m")
            
            # Créer la géométrie 3D
            geometry_3d = self.create_filled_volume_geometry(geometry, fill_data)
            
            if not geometry_3d:
                print(f"   ❌ Impossible de créer la géométrie 3D")
                continue
            
            # Créer l'objet CityJSON
            obj_id, city_object = self.create_cityjson_object(parcelle, geometry_3d)
            
            self.city_objects[obj_id] = city_object
            self.stats['volumes_created'] += 1
            self.stats['total_fill_volume'] += fill_data['fill_volume']
            
            print(f"   ✅ Volume 3D créé: {len(geometry_3d['boundaries'])} faces")
            
            # Affichage périodique du progrès
            if self.stats['parcelles_processed'] % 10 == 0:
                print(f"\n📊 Progression: {self.stats['parcelles_processed']}/{len(self.parcelles)} parcelles")
                print(f"   • Volumes créés: {self.stats['volumes_created']}")
                print(f"   • Volume total de comblement: {self.stats['total_fill_volume']:.0f} m³")
        
        # Calculer les statistiques des hauteurs de comblement
        if fill_heights:
            self.stats['avg_fill_height'] = np.mean(fill_heights)
            self.stats['max_fill_height'] = np.max(fill_heights)
            self.stats['min_calculated_height'] = np.min(fill_heights)
    
    def convert_numpy_types(self, obj):
        """Convertit récursivement les types NumPy en types Python natifs pour JSON"""
        if isinstance(obj, dict):
            return {k: self.convert_numpy_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self.convert_numpy_types(item) for item in obj]
        elif isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return [self.convert_numpy_types(item) for item in obj.tolist()]
        else:
            return obj
    
    def create_cityjson_document(self):
        """Crée le document CityJSON final"""
        print(f"\n Création du document CityJSON...")
        
        ref_x, ref_y, ref_z = self.ref_point
        
        # Calculer l'étendue géographique en coordonnées locales
        if self.vertices:
            xs, ys, zs = zip(*self.vertices)
            geographical_extent = [
                min(xs), min(ys), min(zs),  # min x, y, z
                max(xs), max(ys), max(zs)   # max x, y, z
            ]
        else:
            geographical_extent = [0, 0, 0, 0, 0, 0]
        
        cityjson = {
            "type": "CityJSON",
            "version": "2.0",
            "CityObjects": self.city_objects,
            "vertices": self.vertices,
            "metadata": {
                "referenceSystem": "https://www.opengis.net/def/crs/EPSG/0/2056",
                "geographicalExtent": geographical_extent,
                "presentLoDs": ["2.0"],
                "datasetTitle": f"Backfill volumes"
            },
            "transform": {
                "scale": [1.0, 1.0, 1.0],
                "translate": [ref_x, ref_y, ref_z]
            }
        }
        
        # Convertir tous les types NumPy en types Python natifs
        cityjson = self.convert_numpy_types(cityjson)
        
        # Sauvegarder le fichier
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(self.output_path, 'w', encoding='utf-8') as f:
            json.dump(cityjson, f, indent=2, ensure_ascii=False)
        
        self.stats['vertices_total'] = len(self.vertices)
        
        return True
    
    def generate_filled_volumes(self):
        """Génère tous les volumes comblés"""
        print(f"🏗️ Génération des volumes comblés pour les parcelles")
        print(f"📁 Source: {self.parcelles_path}")
        print(f"🗻 MNT: {self.mnt_path}")
        print(f"📄 Sortie: {self.output_path}")
        print(f"📏 Hauteur minimale de comblement: {self.min_fill_height}m")
        print(f"🧮 Mode: Calcul dynamique des hauteurs de comblement")
        
        # Charger les données
        if not self.load_data():
            return False
        
        # Traiter les parcelles
        self.process_parcelles()
        
        # Créer le document CityJSON
        if not self.create_cityjson_document():
            return False
        
        # Afficher les résultats
        print(f"\n✅ Génération terminée avec succès !")
        print(f"📄 Fichier CityJSON créé: {self.output_path}")
        print(f"📊 Statistiques finales:")
        print(f"   • Parcelles traitées: {self.stats['parcelles_processed']}")
        print(f"   • Volumes 3D créés: {self.stats['volumes_created']}")
        print(f"   • Volume total de comblement: {self.stats['total_fill_volume']:.0f} m³")
        print(f"   • Vertices totaux: {self.stats['vertices_total']}")
        
        if self.stats['volumes_created'] > 0:
            print(f"   • Hauteur de comblement moyenne: {self.stats['avg_fill_height']:.2f}m")
            print(f"   • Hauteur de comblement min/max: {self.stats['min_calculated_height']:.2f}m / {self.stats['max_fill_height']:.2f}m")
        
        efficiency = (self.stats['volumes_created'] / self.stats['parcelles_processed'] * 100) if self.stats['parcelles_processed'] > 0 else 0
        print(f"   • Taux de réussite: {efficiency:.1f}%")
        
        return True


def main():
    """Fonction principale"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Génère des volumes 3D comblés pour les parcelles avec calcul dynamique des hauteurs")
    parser.add_argument("parcelles", help="Chemin vers BienFondsRemblaiements.gpkg")
    parser.add_argument("-m", "--mnt", default="swissALTI3D.tif", help="Chemin vers le MNT (défaut: swissALTI3D.tif)")
    parser.add_argument("-o", "--output", default="output/LandfillVolumes.city.json", help="Fichier de sortie CityJSON")
    parser.add_argument("--min-fill-height", type=float, default=0.5, help="Hauteur minimale de comblement en mètres (défaut: 0.5)")
    parser.add_argument("--list-columns", action="store_true", help="Afficher les colonnes disponibles dans le fichier parcelles")
    
    args = parser.parse_args()
    
    # Vérifier que les fichiers existent
    parcelles_path = Path(args.parcelles)
    mnt_path = Path(args.mnt)
    
    if not parcelles_path.exists():
        print(f"❌ Fichier parcelles non trouvé: {args.parcelles}")
        return 1
    
    if not mnt_path.exists():
        print(f"❌ Fichier MNT non trouvé: {args.mnt}")
        return 1
    
    # Si demandé, afficher les colonnes et quitter
    if args.list_columns:
        print(f"📋 Analyse du fichier: {args.parcelles}")
        try:
            gdf = gpd.read_file(args.parcelles)
            print(f"📊 {len(gdf)} parcelles trouvées")
            print(f"📍 CRS: {gdf.crs}")
            print(f"📋 Colonnes disponibles:")
            for col in gdf.columns:
                if col != 'geometry':
                    print(f"   • {col}")
                    # Afficher quelques exemples de valeurs
                    sample_values = gdf[col].dropna().unique()[:5]
                    if len(sample_values) > 0:
                        print(f"     Exemples: {list(sample_values)}")
        except Exception as e:
            print(f"❌ Erreur lors de la lecture du fichier: {e}")
            return 1
        return 0
    
    # Créer le générateur et lancer la génération
    generator = FilledVolumeGenerator(
        parcelles_path=args.parcelles,
        mnt_path=args.mnt,
        output_path=args.output,
        min_fill_height=args.min_fill_height
    )
    
    success = generator.generate_filled_volumes()
    return 0 if success else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
