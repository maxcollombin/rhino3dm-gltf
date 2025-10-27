#!/usr/bin/env python3
"""
Script de conversion GDB/GeoPackage vers CityJSON 2.0 (optimisé)
Convertit des géométries 3D depuis une géodatabase (.gdb) ou GeoPackage (.gpkg) vers le format CityJSON
avec regroupement sémantique des surfaces par bâtiment.

Gestion UUID conforme à exportRhinoLocalRef2CityJSON.py:
- Génération de nouveaux UUIDs aléatoires pour tous les objets (Building et BuildingPart)
- Suppression des attributs UUID redondants dans les objets
- Relations parent-enfant explicites selon la spécification CityJSON
- Pas de pattern descriptif dans les UUIDs
- Support des formats GDB et GeoPackage
"""

import json
import uuid
import os
import sys
import fnmatch
from pathlib import Path
from osgeo import ogr, osr
import numpy as np

class GeoDataToCityJSON:
    def __init__(self, data_path, output_path, crs_epsg=2056, layer_filter=None, attribute_filter=None):
        """
        Initialise le convertisseur GDB/GeoPackage vers CityJSON
        
        Args:
            data_path: Chemin vers la géodatabase (.gdb) ou GeoPackage (.gpkg)
            output_path: Chemin de sortie pour le fichier CityJSON
            crs_epsg: Code EPSG du système de coordonnées de sortie (défaut: 2056 pour CH1903+ / LV95)
            layer_filter: Liste des noms de couches à traiter (None = toutes les couches)
            attribute_filter: Dictionnaire de filtrage des attributs {
                'include': ['attr1', 'attr2'],  # Attributs à inclure (si spécifié, seuls ces attributs seront inclus)
                'exclude': ['attr3', 'attr4'],  # Attributs à exclure
                'patterns': ['*ID*', 'UUID*']   # Patterns à inclure (wildcards supportés)
            }
        """
        self.data_path = Path(data_path)
        self.output_path = Path(output_path)
        self.crs_epsg = crs_epsg
        self.layer_filter = set(layer_filter) if layer_filter else None
        self.attribute_filter = attribute_filter or {}
        
        # Détecter le type de fichier
        self.data_format = self.detect_data_format()
        
        # Initialiser les variables de transformation
        self.transform = None
        self.vertices = []
        self.vertex_index_map = {}
        self.city_objects = {}
        
        # Dictionnaire pour regrouper les surfaces par bâtiment
        self.buildings_surfaces = {}  # building_id -> {"walls": [], "roofs": [], "floors": []}
        
        # Statistiques
        self.stats = {
            'layers_processed': 0,
            'features_processed': 0,
            'features_filtered': 0,
            'geometries_converted': 0,
            'vertices_total': 0,
            'vertices_cleaned': 0,
            'invalid_geometries_fixed': 0,
            'buildings_created': 0
        }
    
    def detect_data_format(self):
        """Détecte le format du fichier de données"""
        if self.data_path.suffix.lower() == '.gpkg':
            return 'gpkg'
        elif self.data_path.suffix.lower() == '.gdb' or self.data_path.is_dir():
            return 'gdb'
        else:
            raise ValueError(f"Format de fichier non supporté: {self.data_path.suffix}")
    
    def open_datasource(self):
        """Ouvre la source de données selon le format détecté"""
        if self.data_format == 'gpkg':
            driver = ogr.GetDriverByName("GPKG")
            if not driver:
                raise RuntimeError("Driver GeoPackage (GPKG) non disponible")
            
            datasource = driver.Open(str(self.data_path), 0)  # 0 = lecture seule
            if not datasource:
                raise RuntimeError(f"Impossible d'ouvrir le GeoPackage: {self.data_path}")
            
            print(f"📦 GeoPackage ouvert: {self.data_path}")
            return datasource
            
        elif self.data_format == 'gdb':
            # Essayer d'abord OpenFileGDB, puis FileGDB
            driver = ogr.GetDriverByName("OpenFileGDB")
            if not driver:
                driver = ogr.GetDriverByName("FileGDB")
            
            if not driver:
                raise RuntimeError("Driver GDB non disponible. Installez FileGDB ou utilisez OpenFileGDB")
            
            datasource = driver.Open(str(self.data_path), 0)  # 0 = lecture seule
            if not datasource:
                raise RuntimeError(f"Impossible d'ouvrir la géodatabase: {self.data_path}")
            
            print(f"🗃️ Géodatabase ouverte: {self.data_path}")
            return datasource
        
        else:
            raise ValueError(f"Format non supporté: {self.data_format}")
    
    def setup_coordinate_transformation(self, source_srs):
        """Configure la transformation de coordonnées"""
        target_srs = osr.SpatialReference()
        target_srs.ImportFromEPSG(self.crs_epsg)
        
        if source_srs and not source_srs.IsSame(target_srs):
            self.transform = osr.CoordinateTransformation(source_srs, target_srs)
            print(f" Transformation configurée vers EPSG:{self.crs_epsg}")
        else:
            print(f" Coordonnées déjà en EPSG:{self.crs_epsg}")
    
    def add_vertex(self, x, y, z):
        """Ajoute un vertex et retourne son index, avec déduplication ABSOLUE"""
        # Conversion en float pour éviter les problèmes de précision
        x, y, z = float(x), float(y), float(z)
        
        # Clé avec arrondi très précis
        key = (round(x, 3), round(y, 3), round(z, 3))
        
        if key in self.vertex_index_map:
            return self.vertex_index_map[key]
        
        # Double vérification : chercher des coordonnées identiques existantes
        exact_coords = [x, y, z]
        for existing_idx, existing_coords in enumerate(self.vertices):
            # Comparaison avec tolérance très stricte
            if (abs(existing_coords[0] - x) < 1e-6 and 
                abs(existing_coords[1] - y) < 1e-6 and 
                abs(existing_coords[2] - z) < 1e-6):
                self.vertex_index_map[key] = existing_idx
                return existing_idx
        
        idx = len(self.vertices)
        self.vertices.append(exact_coords)
        self.vertex_index_map[key] = idx
        return idx
    
    def transform_geometry(self, geometry):
        """Transforme une géométrie OGR en vertices CityJSON"""
        if not geometry:
            return None
        
        # Appliquer la transformation de coordonnées si nécessaire
        if self.transform:
            geometry.Transform(self.transform)
        
        return self.extract_vertices_from_geometry(geometry)
    
    def extract_vertices_from_geometry(self, geometry):
        """Extrait les vertices d'une géométrie OGR"""
        if not geometry:
            return None
            
        geom_type = geometry.GetGeometryType()
        
        try:
            if geom_type == ogr.wkbPoint25D or geom_type == ogr.wkbPoint:
                x, y, z = geometry.GetX(), geometry.GetY(), geometry.GetZ()
                return [[self.add_vertex(x, y, z)]]  # Point comme surface simple
            
            elif geom_type == ogr.wkbLineString25D or geom_type == ogr.wkbLineString:
                vertices = []
                for i in range(geometry.GetPointCount()):
                    x, y, z = geometry.GetX(i), geometry.GetY(i), geometry.GetZ(i)
                    vertices.append(self.add_vertex(x, y, z))
                if len(vertices) >= 3:
                    # S'assurer que le LineString est fermé pour former une surface valide
                    cleaned_vertices = self.clean_consecutive_vertices(vertices)
                    if cleaned_vertices and len(cleaned_vertices) >= 3:  # 3 points minimum pour ring ouvert
                        return [cleaned_vertices]  # LineString comme surface
                return None
            
            elif geom_type == ogr.wkbPolygon25D or geom_type == ogr.wkbPolygon:
                return self.extract_polygon_vertices(geometry)
            
            elif geom_type == ogr.wkbTriangle or geom_type == 1017:  # Triangle
                return self.extract_triangle_vertices(geometry)
            
            elif geom_type == ogr.wkbMultiPolygon25D or geom_type == ogr.wkbMultiPolygon:
                surfaces = []
                for i in range(geometry.GetGeometryCount()):
                    poly = geometry.GetGeometryRef(i)
                    if poly:
                        poly_surfaces = self.extract_polygon_vertices(poly)
                        if poly_surfaces:
                            surfaces.extend(poly_surfaces)
                return surfaces if surfaces else None
            
            elif geom_type == ogr.wkbPolyhedralSurface or geom_type == ogr.wkbTIN:
                return self.extract_polyhedral_vertices(geometry)
            
            elif geom_type == 1016:  # Type spécifique TIN (Triangulated Irregular Network)
                return self.extract_tin_vertices(geometry)
            
            elif geom_type == -2147483642:  # Type spécifique négatif (souvent TIN avec extension)
                return self.extract_polyhedral_vertices(geometry)
            
            elif geom_type == ogr.wkbGeometryCollection:
                # Collection de géométries - traiter récursivement
                surfaces = []
                for i in range(geometry.GetGeometryCount()):
                    sub_geom = geometry.GetGeometryRef(i)
                    if sub_geom:
                        sub_surfaces = self.extract_vertices_from_geometry(sub_geom)
                        if sub_surfaces:
                            surfaces.extend(sub_surfaces)
                return surfaces if surfaces else None
            
            else:
                # Pour les types non reconnus, essayer d'extraire comme polyédrique
                if geometry.GetGeometryCount() > 0:
                    return self.extract_polyhedral_vertices(geometry)
                print(f" Type de géométrie non supporté : {geom_type} (nom: {geometry.GetGeometryName()})")
                return None
                
        except Exception as e:
            print(f" Erreur extraction géométrie type {geom_type}: {e}")
            return None
    
    def clean_consecutive_vertices(self, vertex_indices):
        """
        Supprime ABSOLUMENT TOUS les vertices consécutifs identiques
        
        IMPORTANT: Cette fonction génère des rings OUVERTS (non fermés explicitement)
        car val3dity considère les rings fermés comme une erreur CONSECUTIVE_POINTS_SAME.
        Le premier et dernier point NE SONT PAS identiques dans le résultat.
        """
        if not vertex_indices or len(vertex_indices) < 3:
            return vertex_indices
        
        # ÉTAPE 0: Supprimer la fermeture explicite si elle existe
        # Si le premier et dernier point sont identiques, supprimer le dernier
        if len(vertex_indices) > 3 and vertex_indices[0] == vertex_indices[-1]:
            vertex_indices = vertex_indices[:-1]
            self.stats['vertices_cleaned'] += 1
        
        original_count = len(vertex_indices)
        
        # Phase 1: Nettoyer tous les doublons consécutifs de manière ultra-agressive
        cleaned = []
        
        for i, vertex_idx in enumerate(vertex_indices):
            if vertex_idx >= len(self.vertices):
                continue  # Index invalide
                
            curr_coords = self.vertices[vertex_idx]
            
            # Ajouter le vertex seulement s'il est différent de TOUS les précédents récents
            should_add = True
            
            if len(cleaned) > 0:
                # Vérifier contre le précédent
                prev_coords = self.vertices[cleaned[-1]]
                if curr_coords == prev_coords:
                    should_add = False
                    self.stats['vertices_cleaned'] += 1
            
            if should_add and len(cleaned) > 1:
                # Vérifier aussi contre l'avant-précédent pour éviter A-B-A patterns
                prev_prev_coords = self.vertices[cleaned[-2]]
                if curr_coords == prev_prev_coords:
                    should_add = False
                    self.stats['vertices_cleaned'] += 1
            
            if should_add:
                cleaned.append(vertex_idx)
        
        # Phase 2: Nettoyage final - vérifier tous les points consécutifs encore une fois
        # Note: Ne PAS fermer explicitement les rings car val3dity considère cela comme une erreur
        if len(cleaned) >= 3:
            final_cleaned = []
            for i, vertex_idx in enumerate(cleaned):
                curr_coords = self.vertices[vertex_idx]
                
                # Pour tous sauf le dernier, vérifier contre le suivant
                if i < len(cleaned) - 1:
                    next_coords = self.vertices[cleaned[i + 1]]
                    if curr_coords == next_coords:
                        # Point identique au suivant - ignorer ce point
                        self.stats['vertices_cleaned'] += 1
                        continue
                
                final_cleaned.append(vertex_idx)
            
            # Vérifier une dernière fois que nous avons assez de points
            # Note: val3dity attend des rings ouverts (pas de fermeture explicite)
            if len(final_cleaned) >= 3:
                return final_cleaned
        
        return None

    def extract_triangle_vertices(self, triangle):
        """Extrait les vertices d'un triangle et s'assure qu'il est correctement fermé"""
        if not triangle:
            return None
        
        try:
            vertices = []
            
            # Si le triangle a des sous-géométries (comme un LINEARRING)
            if triangle.GetGeometryCount() > 0:
                # Prendre la première sous-géométrie (généralement le LINEARRING)
                ring = triangle.GetGeometryRef(0)
                if ring and ring.GetPointCount() >= 3:
                    # Extraire tous les points du ring
                    for i in range(ring.GetPointCount()):
                        x, y, z = ring.GetX(i), ring.GetY(i), ring.GetZ(i)
                        vertices.append(self.add_vertex(x, y, z))
            
            # Si le triangle a des points directement
            elif triangle.GetPointCount() >= 3:
                for i in range(triangle.GetPointCount()):
                    x, y, z = triangle.GetX(i), triangle.GetY(i), triangle.GetZ(i)
                    vertices.append(self.add_vertex(x, y, z))
            
            # Nettoyer et fermer correctement le ring
            if len(vertices) >= 3:
                cleaned_vertices = self.clean_consecutive_vertices(vertices)
                if cleaned_vertices and len(cleaned_vertices) >= 3:  # 3 points minimum pour ring ouvert
                    return [[cleaned_vertices]]
            
        except Exception as e:
            print(f"     Erreur extraction triangle: {e}")
        
        return None

    def extract_tin_vertices(self, tin_geometry):
        """Extrait les vertices d'un TIN (Triangulated Irregular Network)"""
        if not tin_geometry:
            return None
        
        surfaces = []
        
        try:
            # Un TIN contient de nombreux triangles
            for i in range(tin_geometry.GetGeometryCount()):
                triangle = tin_geometry.GetGeometryRef(i)
                if triangle and triangle.GetGeometryType() in [ogr.wkbTriangle, 1017]:
                    triangle_surfaces = self.extract_triangle_vertices(triangle)
                    if triangle_surfaces:
                        surfaces.extend(triangle_surfaces)
            
            return surfaces if surfaces else None
            
        except Exception as e:
            print(f"     Erreur extraction TIN: {e}")
            return None

    def extract_polygon_vertices(self, polygon):
        """Extrait les vertices d'un polygone (exterior + holes) avec fermeture correcte"""
        if not polygon or polygon.GetGeometryCount() == 0:
            return None
            
        surfaces = []
        
        try:
            # Ring extérieur
            exterior_ring = polygon.GetGeometryRef(0)
            if exterior_ring and exterior_ring.GetPointCount() >= 3:
                exterior_vertices = []
                for i in range(exterior_ring.GetPointCount()):
                    x, y, z = exterior_ring.GetX(i), exterior_ring.GetY(i), exterior_ring.GetZ(i)
                    exterior_vertices.append(self.add_vertex(x, y, z))
                
                # Nettoyer et s'assurer de la fermeture
                cleaned_exterior = self.clean_consecutive_vertices(exterior_vertices)
                
                if cleaned_exterior and len(cleaned_exterior) >= 3:  # 3 points minimum pour ring ouvert
                    current_surface = [cleaned_exterior]
                    
                    # Rings intérieurs (trous)
                    for i in range(1, polygon.GetGeometryCount()):
                        interior_ring = polygon.GetGeometryRef(i)
                        if interior_ring and interior_ring.GetPointCount() >= 3:
                            interior_vertices = []
                            for j in range(interior_ring.GetPointCount()):
                                x, y, z = interior_ring.GetX(j), interior_ring.GetY(j), interior_ring.GetZ(j)
                                interior_vertices.append(self.add_vertex(x, y, z))
                            
                            # Nettoyer et s'assurer de la fermeture
                            cleaned_interior = self.clean_consecutive_vertices(interior_vertices)
                            
                            if cleaned_interior and len(cleaned_interior) >= 3:  # 3 points minimum pour ring ouvert
                                current_surface.append(cleaned_interior)
                    
                    surfaces.append(current_surface)
                    
        except Exception as e:
            print(f"     Erreur extraction polygone: {e}")
            return None
        
        return surfaces if surfaces else None
    
    def extract_polyhedral_vertices(self, geometry):
        """Extrait les vertices d'une surface polyédrique"""
        if not geometry:
            return None
            
        surfaces = []
        
        try:
            # Si c'est directement un polygone
            if geometry.GetGeometryType() in [ogr.wkbPolygon, ogr.wkbPolygon25D]:
                return self.extract_polygon_vertices(geometry)
            
            # Si c'est directement un triangle
            elif geometry.GetGeometryType() in [ogr.wkbTriangle, 1017]:
                return self.extract_triangle_vertices(geometry)
            
            # Traiter les sous-géométries
            for i in range(geometry.GetGeometryCount()):
                surface = geometry.GetGeometryRef(i)
                if surface:
                    surface_type = surface.GetGeometryType()
                    
                    if surface_type in [ogr.wkbPolygon25D, ogr.wkbPolygon]:
                        poly_surfaces = self.extract_polygon_vertices(surface)
                        if poly_surfaces:
                            surfaces.extend(poly_surfaces)
                            
                    elif surface_type in [ogr.wkbTriangle, 1017]:
                        # Triangle - traiter comme un polygone simple
                        triangle_surfaces = self.extract_triangle_vertices(surface)
                        if triangle_surfaces:
                            surfaces.extend(triangle_surfaces)
                            
                    elif surface_type in [ogr.wkbLineString, ogr.wkbLineString25D]:
                        # Traiter comme un contour fermé
                        if surface.GetPointCount() >= 3:
                            vertices = []
                            for j in range(surface.GetPointCount()):
                                x, y, z = surface.GetX(j), surface.GetY(j), surface.GetZ(j)
                                vertices.append(self.add_vertex(x, y, z))
                            
                            # Nettoyer et s'assurer de la fermeture
                            cleaned_vertices = self.clean_consecutive_vertices(vertices)
                            
                            if cleaned_vertices and len(cleaned_vertices) >= 3:  # 3 points minimum pour ring ouvert
                                surfaces.append([cleaned_vertices])
                    
                    else:
                        # Traiter récursivement d'autres types
                        sub_surfaces = self.extract_vertices_from_geometry(surface)
                        if sub_surfaces:
                            surfaces.extend(sub_surfaces)
                            
        except Exception as e:
            print(f"     Erreur extraction polyédrique : {e}")
            return None
        
        return surfaces if surfaces else None
    
    def filter_attributes(self, attributes):
        """
        Filtre les attributs selon les critères définis
        
        Args:
            attributes: Dictionnaire des attributs à filtrer
            
        Returns:
            Dictionnaire des attributs filtrés
        """
        if not self.attribute_filter:
            return attributes
        
        import fnmatch
        
        filtered_attributes = {}
        include_list = self.attribute_filter.get('include', [])
        exclude_list = self.attribute_filter.get('exclude', [])
        pattern_list = self.attribute_filter.get('patterns', [])
        
        for attr_name, attr_value in attributes.items():
            should_include = True
            
            # Si une liste d'inclusion est spécifiée, seuls ces attributs sont autorisés
            if include_list:
                should_include = attr_name in include_list
            
            # Vérifier les patterns d'inclusion (wildcards)
            if not should_include and pattern_list:
                for pattern in pattern_list:
                    if fnmatch.fnmatch(attr_name, pattern):
                        should_include = True
                        break
            
            # Exclure les attributs dans la liste d'exclusion
            if should_include and exclude_list:
                if attr_name in exclude_list:
                    should_include = False
            
            # Exclure systématiquement les attributs UUID/ID redondants (conforme à exportRhinoLocalRef2CityJSON.py)
            if attr_name.lower() in ['buildingid', 'uuid', 'id', 'objectid', 'gid', 'fid']:
                should_include = False
            
            if should_include:
                filtered_attributes[attr_name] = attr_value
        
        return filtered_attributes

    def extract_building_id(self, feature, layer_name):
        """Extrait l'ID du bâtiment depuis les attributs de la feature avec logique améliorée"""
        # Récupérer les attributs de la feature
        feature_defn = feature.GetDefnRef()
        building_id = None
        
        # Chercher des champs qui pourraient contenir l'ID du bâtiment (élargi)
        possible_fields = [
            'building_id', 'buildingid', 'id', 'name', 'nom', 'bat_id', 'building',
            'gid', 'objectid', 'fid', 'uid', 'uuid', 'egid', 'ewid'
        ]
        
        # Premier passage: recherche exacte
        for i in range(feature_defn.GetFieldCount()):
            field_name = feature_defn.GetFieldDefn(i).GetName().lower()
            field_value = feature.GetField(i)
            
            if field_value is not None and field_name in possible_fields:
                building_id = str(field_value)
                break
        
        # Deuxième passage: recherche partielle si pas trouvé
        if not building_id:
            for i in range(feature_defn.GetFieldCount()):
                field_name = feature_defn.GetFieldDefn(i).GetName().lower()
                field_value = feature.GetField(i)
                
                if field_value is not None and any(pf in field_name for pf in possible_fields):
                    building_id = str(field_value)
                    break
        
        # Si pas d'ID trouvé, générer un ID basé sur la position ou les coordonnées
        if not building_id:
            fid = feature.GetFID()
            if fid >= 0:
                building_id = f"{layer_name}_{fid}"
            else:
                # Utiliser les coordonnées de la géométrie pour créer un ID unique
                geometry = feature.GetGeometryRef()
                if geometry:
                    envelope = geometry.GetEnvelope()  # (minX, maxX, minY, maxY)
                    coord_hash = hash((round(envelope[0], 1), round(envelope[2], 1))) % 10000
                    building_id = f"{layer_name}_coord_{coord_hash}"
                else:
                    building_id = f"{layer_name}_unknown_{len(self.buildings_surfaces)}"
        
        return building_id
    
    def classify_geometry(self, feature, layer_name):
        """Classifie une géométrie pour déterminer le type d'objet CityJSON"""
        # Classification basée sur les attributs ou le nom de couche
        attributes = {}
        
        # Récupérer les attributs de la feature
        feature_defn = feature.GetDefnRef()
        for i in range(feature_defn.GetFieldCount()):
            field_name = feature_defn.GetFieldDefn(i).GetName()
            field_value = feature.GetField(i)
            if field_value is not None:
                attributes[field_name] = field_value
        
        # Classification basique (à adapter selon vos données)
        city_object_type = "GenericCityObject"
        function = "unknown"
        
        # Exemples de classification (à adapter)
        layer_lower = layer_name.lower()
        if "building" in layer_lower or "batiment" in layer_lower:
            city_object_type = "Building"
            function = "residential"
        elif "road" in layer_lower or "route" in layer_lower or "rue" in layer_lower:
            city_object_type = "Road"
            function = "traffic"
        elif "bridge" in layer_lower or "pont" in layer_lower:
            city_object_type = "Bridge"
            function = "transport"
        elif "water" in layer_lower or "eau" in layer_lower:
            city_object_type = "WaterBody"
            function = "hydrography"
        elif "vegetation" in layer_lower or "plant" in layer_lower:
            city_object_type = "PlantCover"
            function = "vegetation"
        
        return city_object_type, function, attributes
    
    def validate_and_clean_boundaries(self, boundaries):
        """Valide et nettoie les boundaries en supprimant les vertices consécutifs identiques"""
        if not boundaries:
            return boundaries
            
        def clean_ring(ring):
            if not ring or len(ring) < 3:  # 3 points minimum pour ring ouvert
                return None
            
            # Les ring sont des listes d'indices de vertices
            cleaned_ring = self.clean_consecutive_vertices(ring)
            
            if cleaned_ring and len(cleaned_ring) >= 3:  # 3 points minimum pour ring ouvert
                return cleaned_ring
            return None
        
        # Nettoyer selon la structure des boundaries
        # boundaries est généralement une liste de surfaces, chaque surface étant une liste de rings
        cleaned_boundaries = []
        
        for surface in boundaries:
            if not surface:
                continue
                
            cleaned_surface = []
            
            # Chaque surface est une liste de rings (exterior + holes)
            for ring in surface:
                if isinstance(ring, list) and len(ring) >= 3:  # 3 points minimum pour ring ouvert
                    cleaned_ring = clean_ring(ring)
                    if cleaned_ring:
                        cleaned_surface.append(cleaned_ring)
            
            if cleaned_surface:
                cleaned_boundaries.append(cleaned_surface)
        
        return cleaned_boundaries if cleaned_boundaries else None

    def create_city_object(self, feature, layer_name, boundaries):
        """Crée un objet CityJSON à partir d'une feature"""
        # Valider et nettoyer les boundaries
        cleaned_boundaries = self.validate_and_clean_boundaries(boundaries)
        if not cleaned_boundaries:
            return None, None
        
        city_object_type, function, attributes = self.classify_geometry(feature, layer_name)
        
        obj_id = str(uuid.uuid4())
        
        # Analyser la structure des boundaries pour choisir le bon type de géométrie
        total_surfaces = len(cleaned_boundaries)
        
        if total_surfaces == 1:
            # Une seule surface
            surface = cleaned_boundaries[0]
            if isinstance(surface, list) and len(surface) == 1:
                # Surface simple avec un seul ring
                geometry = {
                    "type": "MultiSurface",
                    "lod": "2.0",
                    "boundaries": cleaned_boundaries
                }
            else:
                # Surface avec trous ou géométrie complexe
                geometry = {
                    "type": "MultiSurface",
                    "lod": "2.0",
                    "boundaries": cleaned_boundaries
                }
        elif total_surfaces <= 6:
            # Nombre de surfaces compatible avec un Solid (6 faces max pour un cube)
            geometry = {
                "type": "Solid",
                "lod": "2.0",
                "boundaries": [cleaned_boundaries]  # Un seul shell
            }
        else:
            # Beaucoup de surfaces - utiliser MultiSurface
            geometry = {
                "type": "MultiSurface",
                "lod": "2.0",
                "boundaries": cleaned_boundaries
            }
        
        city_object = {
            "type": city_object_type,
            "attributes": {
                "function": function,
                "layer": layer_name,
                "surfaces_count": total_surfaces,
                **attributes  # Ajouter tous les attributs de la feature
            },
            "geometry": [geometry]
        }
        
        return obj_id, city_object
    
    def analyze_objektart_values(self, layer):
        """Analyse les valeurs d'OBJEKTART présentes dans la couche pour débogage"""
        objektart_values = set()
        objektart_index = -1
        
        # Vérifier si le champ existe
        feature_defn = layer.GetLayerDefn()
        for i in range(feature_defn.GetFieldCount()):
            if feature_defn.GetFieldDefn(i).GetName() == 'OBJEKTART':
                objektart_index = i
                break
        
        if objektart_index < 0:
            print(f"     Aucun champ OBJEKTART trouvé dans la couche {layer.GetName()}")
            return objektart_values
        
        # Parcourir les features pour collecter les valeurs
        layer.ResetReading()
        feature = layer.GetNextFeature()
        sample_count = 0
        
        while feature and sample_count < 1000:  # Limiter à 1000 samples pour performance
            objektart_value = feature.GetField('OBJEKTART')
            if objektart_value is not None:
                objektart_values.add(str(objektart_value).strip())
            sample_count += 1
            feature = layer.GetNextFeature()
        
        # Afficher les valeurs trouvées
        if objektart_values:
            print(f"     Valeurs OBJEKTART trouvées dans {layer.GetName()} (échantillon de {sample_count}):")
            for value in sorted(objektart_values):
                print(f"       • '{value}'")
        
        layer.ResetReading()  # Reset pour le traitement principal
        return objektart_values
    
    def apply_objektart_filter(self, layer):
        """Applique un filtre SQL pour exclure Flugdach et Offenes Gebaeude au niveau de la couche"""
        # Vérifier si le champ OBJEKTART existe
        feature_defn = layer.GetLayerDefn()
        has_objektart = False
        
        for i in range(feature_defn.GetFieldCount()):
            if feature_defn.GetFieldDefn(i).GetName() == 'OBJEKTART':
                has_objektart = True
                break
        
        if not has_objektart:
            print(f"     Aucun champ OBJEKTART trouvé - aucun filtrage appliqué")
            return
        
        # Compter les features avant filtrage
        total_before = layer.GetFeatureCount()
        
        # Appliquer le filtre SQL pour exclure les valeurs non désirées
        filter_sql = "OBJEKTART NOT IN ('Flugdach', 'Offenes Gebaeude')"
        layer.SetAttributeFilter(filter_sql)
        
        # Compter les features après filtrage
        total_after = layer.GetFeatureCount()
        filtered_count = total_before - total_after
        
        print(f"     Filtre OBJEKTART appliqué: {filtered_count} features exclues ({total_after} restantes sur {total_before})")
        
        return filtered_count
    
    def process_layer(self, layer):
        """Traite une couche de la géodatabase et regroupe par bâtiment"""
        layer_name = layer.GetName()
        total_features_original = layer.GetFeatureCount()
        
        print(f" Traitement de la couche : {layer_name} ({total_features_original} features)")
        
        # Analyser les valeurs OBJEKTART pour débogage (avant filtrage)
        objektart_values = self.analyze_objektart_values(layer)
        
        # Appliquer le filtre OBJEKTART au niveau de la couche
        filtered_count = self.apply_objektart_filter(layer)
        
        # Nouveau total après filtrage
        total_features = layer.GetFeatureCount()
        
        # Configuration de la transformation de coordonnées
        source_srs = layer.GetSpatialRef()
        if not self.transform:
            self.setup_coordinate_transformation(source_srs)
        
        # Statistiques détaillées
        feature_count = 0
        converted_count = 0
        error_count = 0
        empty_geometry_count = 0
        
        # Traiter chaque feature (déjà filtrée)
        layer.ResetReading()
        feature = layer.GetNextFeature()
        
        while feature:
            feature_count += 1
            try:
                geometry = feature.GetGeometryRef()
                if geometry:
                    # Extraire l'ID du bâtiment
                    building_id = self.extract_building_id(feature, layer_name)
                    
                    # Déterminer le type de surface pour cette feature spécifique
                    surface_type = self.determine_surface_type_from_feature(feature, layer_name)
                    
                    # Cloner la géométrie pour éviter les modifications
                    geom_clone = geometry.Clone()
                    geom_type = geom_clone.GetGeometryType()
                    
                    # Diagnostic de la géométrie
                    if feature_count <= 5:  # Debug pour les premières features
                        objektart = feature.GetField('OBJEKTART') if feature.GetFieldIndex('OBJEKTART') >= 0 else 'N/A'
                        surface_type_attr = feature.GetField('SurfaceType') if feature.GetFieldIndex('SurfaceType') >= 0 else 'N/A'
                        print(f"     Feature {feature_count}: type géométrie = {geom_type}, building_id = {building_id}, OBJEKTART = '{objektart}', SurfaceType = '{surface_type_attr}', surface_type = '{surface_type}'")
                    
                    # Transformer la géométrie
                    boundaries = self.transform_geometry(geom_clone)
                    
                    if boundaries and len(boundaries) > 0:
                        # Vérifier que les boundaries ne sont pas vides
                        valid_boundaries = [b for b in boundaries if b and len(b) > 0]
                        
                        if valid_boundaries:
                            # Ajouter les surfaces au bâtiment correspondant
                            self.add_surfaces_to_building(building_id, surface_type, valid_boundaries, feature, layer_name)
                            converted_count += 1
                        else:
                            empty_geometry_count += 1
                            if feature_count <= 3:
                                print(f"     Feature {feature_count}: boundaries vides après filtrage")
                    else:
                        empty_geometry_count += 1
                        if feature_count <= 3:
                            print(f"     Feature {feature_count}: aucune boundary générée")
                else:
                    empty_geometry_count += 1
                    if feature_count <= 3:
                        print(f"     Feature {feature_count}: géométrie nulle")
                
                # Affichage périodique du progrès
                if feature_count % 50 == 0:
                    print(f"   Progression: {feature_count}/{total_features} features ({converted_count} converties, {error_count} erreurs, {empty_geometry_count} vides)")
                    
            except Exception as e:
                error_count += 1
                if error_count <= 5:  # Afficher les premières erreurs
                    print(f"   Erreur feature {feature_count}: {e}")
            
            feature = layer.GetNextFeature()
        
        # Résumé détaillé
        print(f"   Couche {layer_name} terminée:")
        print(f"    • Features originales: {total_features_original}")
        print(f"    • Features filtrées (OBJEKTART): {filtered_count if filtered_count else 0}")
        print(f"    • Features restantes: {total_features}")
        print(f"    • Features traitées: {feature_count}")
        print(f"    • Converties avec succès: {converted_count}")
        print(f"    • Géométries vides/nulles: {empty_geometry_count}")
        print(f"    • Erreurs: {error_count}")
        
        self.stats['features_processed'] += feature_count
        self.stats['features_filtered'] += (filtered_count if filtered_count else 0)
        self.stats['layers_processed'] += 1
    
    def convert(self):
        """Effectue la conversion complète"""
        print(f"🔄 Ouverture de la source de données : {self.data_path}")
        
        # Ouvrir la source de données selon le format
        try:
            datasource = self.open_datasource()
        except (RuntimeError, ValueError) as e:
            print(f"❌ Erreur: {e}")
            return False
        
        total_layers = datasource.GetLayerCount()
        print(f"📊 Source de données ouverte : {total_layers} couches trouvées")
        
        # Lister toutes les couches disponibles
        available_layers = []
        for i in range(total_layers):
            layer = datasource.GetLayer(i)
            available_layers.append(layer.GetName())
        
        # Filtrer les couches si nécessaire
        layers_to_process = []
        if self.layer_filter:
            print(f" Filtrage des couches : {', '.join(self.layer_filter)}")
            for i in range(total_layers):
                layer = datasource.GetLayer(i)
                layer_name = layer.GetName()
                if layer_name in self.layer_filter:
                    layers_to_process.append((i, layer_name))
            
            # Vérifier que toutes les couches demandées existent
            missing_layers = self.layer_filter - {name for _, name in layers_to_process}
            if missing_layers:
                print(f" ⚠️ Couches non trouvées : {', '.join(missing_layers)}")
                print(f" Couches disponibles : {', '.join(available_layers)}")
        else:
            # Traiter toutes les couches
            for i in range(total_layers):
                layer = datasource.GetLayer(i)
                layers_to_process.append((i, layer.GetName()))
        
        print(f" Couches à traiter : {len(layers_to_process)} sur {total_layers}")
        
        # Traiter chaque couche sélectionnée
        for layer_idx, layer_name in layers_to_process:
            layer = datasource.GetLayer(layer_idx)
            self.process_layer(layer)
        
        # Créer les bâtiments à partir des surfaces regroupées
        self.create_buildings_from_surfaces()
        
        # Créer le document CityJSON
        self.create_cityjson_document()
        
        # Nettoyer
        datasource = None
        
        return True
    
    def create_cityjson_document(self):
        """Crée le document CityJSON final"""
        print(f" Création du document CityJSON...")
        
        # Point de référence pour la transformation locale (même que NouveauxBatiments.city.json)
        translate_x = 2592980.685
        translate_y = 1119281.703
        translate_z = 483.8  # Altitude de référence
        
        # Appliquer la transformation locale aux vertices
        local_vertices = []
        for vertex in self.vertices:
            x_local = vertex[0] - translate_x
            y_local = vertex[1] - translate_y
            z_local = vertex[2] - translate_z
            local_vertices.append([x_local, y_local, z_local])
        
        # Calculer l'étendue géographique en coordonnées locales
        if local_vertices:
            xs, ys, zs = zip(*local_vertices)
            geographical_extent = [
                min(xs),  # minx
                min(ys),  # miny
                min(zs),  # minz
                max(xs),  # maxx
                max(ys),  # maxy
                max(zs)   # maxz
            ]
        else:
            geographical_extent = [0, 0, 0, 0, 0, 0]
        
        # Créer le document CityJSON
        cityjson = {
            "type": "CityJSON",
            "version": "2.0",
            "CityObjects": self.city_objects,
            "vertices": local_vertices,
            "metadata": {
                "referenceSystem": f"https://www.opengis.net/def/crs/EPSG/0/{self.crs_epsg}",
                "geographicalExtent": geographical_extent,
                "presentLoDs": ["2.0"],
                "datasetTitle": "swissBUILDINGS3D"
            },
            "transform": {
                "scale": [1.0, 1.0, 1.0],
                "translate": [translate_x, translate_y, translate_z]
            }
        }
        
        # Sauvegarder
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(self.output_path, 'w', encoding='utf-8') as f:
            json.dump(cityjson, f, indent=2, ensure_ascii=False)
        
        # Statistiques finales
        self.stats['vertices_total'] = len(self.vertices)
        self.stats['geometries_converted'] = len(self.city_objects)
        
        # Compter les différents types d'objets
        buildings_count = sum(1 for obj in self.city_objects.values() if obj['type'] == 'Building')
        building_parts_count = sum(1 for obj in self.city_objects.values() if obj['type'] == 'BuildingPart')
        
        print(f" Conversion terminée !")
        print(f" Fichier créé : {self.output_path}")
        print(f" Système de coordonnées : EPSG:{self.crs_epsg}")
        print(f" Transformation locale appliquée : translate=[{translate_x}, {translate_y}, {translate_z}]")
        print(f" Étendue géographique locale : {geographical_extent}")
        print(f" Statistiques :")
        print(f"  • Couches traitées : {self.stats['layers_processed']}")
        print(f"  • Features examinées : {self.stats['features_processed']}")
        print(f"  • Features filtrées (OBJEKTART) : {self.stats['features_filtered']}")
        print(f"  • Objets CityJSON créés : {self.stats['geometries_converted']}")
        print(f"    - Buildings : {buildings_count}")
        print(f"    - BuildingParts : {building_parts_count}")
        print(f"  • Vertices totaux : {self.stats['vertices_total']}")
        print(f"  • Vertices nettoyés (doublons supprimés) : {self.stats['vertices_cleaned']}")
        print(" CityJSON créé avec succès !")
        print(" Validation recommandée avec http://geovalidation.bk.tudelft.nl/val3dity/")

    def determine_surface_type_from_feature(self, feature, layer_name):
        """Détermine le type de surface depuis les attributs de la feature, avec fallback sur le nom de couche"""
        # D'abord, essayer d'utiliser l'attribut SurfaceType de la feature
        surface_type_field_idx = feature.GetFieldIndex('SurfaceType')
        if surface_type_field_idx >= 0:
            surface_type_value = feature.GetField('SurfaceType')
            if surface_type_value:
                surface_type_lower = str(surface_type_value).lower()
                # Mapper les valeurs connues
                if 'roof' in surface_type_lower or 'toit' in surface_type_lower:
                    return 'RoofSurface'
                elif 'wall' in surface_type_lower or 'mur' in surface_type_lower:
                    return 'WallSurface'
                elif 'floor' in surface_type_lower or 'ground' in surface_type_lower or 'sol' in surface_type_lower:
                    return 'GroundSurface'
        
        # Si pas d'attribut SurfaceType ou valeur non reconnue, utiliser la logique basée sur le nom de couche
        return self.determine_surface_type(layer_name)

    def determine_surface_type(self, layer_name):
        """Détermine le type de surface selon le nom de la couche avec logique enrichie"""
        layer_lower = layer_name.lower()
        
        # Classification des murs
        wall_keywords = ['wall', 'mur', 'facade', 'façade', 'side', 'côté', 'paroi', 'vertical']
        if any(keyword in layer_lower for keyword in wall_keywords):
            return 'WallSurface'
        
        # Classification des toits
        roof_keywords = ['roof', 'toit', 'toiture', 'couverture', 'top', 'sommet', 'haut']
        if any(keyword in layer_lower for keyword in roof_keywords):
            return 'RoofSurface'
        
        # Classification des sols/planchers
        floor_keywords = ['floor', 'sol', 'plancher', 'ground', 'base', 'bottom', 'bas', 'fondation']
        if any(keyword in layer_lower for keyword in floor_keywords):
            return 'GroundSurface'
        
        # Classification par position dans le nom de la couche
        if any(pos in layer_lower for pos in ['_top_', '_roof_', '_r_']):
            return 'RoofSurface'
        elif any(pos in layer_lower for pos in ['_bottom_', '_floor_', '_f_', '_ground_']):
            return 'GroundSurface'
        elif any(pos in layer_lower for pos in ['_wall_', '_w_', '_side_']):
            return 'WallSurface'
        
        # Par défaut, considérer comme un mur (surface verticale la plus commune)
        return 'WallSurface'
    
    def add_surfaces_to_building(self, building_id, surface_type, boundaries, feature, layer_name):
        """Ajoute des surfaces à un bâtiment"""
        if building_id not in self.buildings_surfaces:
            self.buildings_surfaces[building_id] = {
                'walls': [],
                'roofs': [],
                'floors': [],
                'attributes': self.extract_building_attributes(feature, layer_name)
            }
        
        # Ajouter les surfaces selon leur type
        if surface_type == 'WallSurface':
            self.buildings_surfaces[building_id]['walls'].extend(boundaries)
        elif surface_type == 'RoofSurface':
            self.buildings_surfaces[building_id]['roofs'].extend(boundaries)
        elif surface_type == 'GroundSurface':
            self.buildings_surfaces[building_id]['floors'].extend(boundaries)
    
    def extract_building_attributes(self, feature, layer_name):
        """Extrait les attributs d'un bâtiment depuis une feature et les filtre (optimisé CityJSON 2.0)"""
        attributes = {}
        
        # Récupérer les attributs de la feature
        feature_defn = feature.GetDefnRef()
        for i in range(feature_defn.GetFieldCount()):
            field_name = feature_defn.GetFieldDefn(i).GetName()
            field_value = feature.GetField(i)
            if field_value is not None:
                attributes[field_name] = field_value
        
        # Appliquer le filtrage des attributs
        # Note: Plus besoin de gérer buildingID ici car il sera supprimé lors de la création des objets
        return self.filter_attributes(attributes)
    
    def create_buildings_from_surfaces(self):
        """Crée les objets Building et BuildingPart à partir des surfaces regroupées avec sémantique optimisée selon CityJSON 2.0"""
        print(f" Création de {len(self.buildings_surfaces)} bâtiments avec BuildingParts...")
        
        # Statistiques de regroupement
        total_wall_surfaces = 0
        total_roof_surfaces = 0
        total_floor_surfaces = 0
        
        for building_id, building_data in self.buildings_surfaces.items():
            # --- Créer l'ID principal du Building (génération de nouveaux UUIDs sans pattern descriptif) ---
            # Génération de nouveaux UUIDs pour tous les objets comme dans exportRhinoLocalRef2CityJSON.py
            building_uuid = str(uuid.uuid4())
            
            child_uuids = []
            
            # Combiner toutes les surfaces pour le Building
            all_boundaries = []
            all_semantics_values = []
            all_semantics_faces = []
            face_index = 0
            
            # Ajouter les murs
            for wall_boundary in building_data['walls']:
                all_boundaries.append(wall_boundary)
                all_semantics_values.append("WallSurface")
                all_semantics_faces.append(face_index)
                face_index += 1
                total_wall_surfaces += 1
            
            # Ajouter les toits
            for roof_boundary in building_data['roofs']:
                all_boundaries.append(roof_boundary)
                all_semantics_values.append("RoofSurface")
                all_semantics_faces.append(face_index)
                face_index += 1
                total_roof_surfaces += 1
            
            # Ajouter les sols
            for floor_boundary in building_data['floors']:
                all_boundaries.append(floor_boundary)
                all_semantics_values.append("GroundSurface")
                all_semantics_faces.append(face_index)
                face_index += 1
                total_floor_surfaces += 1
            
            if all_boundaries:
                # Nettoyer les boundaries
                cleaned_boundaries = self.validate_and_clean_boundaries(all_boundaries)
                
                if cleaned_boundaries:
                    # Préparer les attributs (suppression de tous les attributs UUID)
                    attributes = building_data['attributes'].copy()
                    # Supprimer tous les attributs UUID redondants car l'ID est déjà la clé dans CityObjects
                    attributes.pop('buildingID', None)
                    attributes.pop('UUID', None)
                    attributes.pop('uuid', None)
                    
                    # Ajuster la sémantique pour correspondre aux boundaries nettoyées
                    final_semantics_values = all_semantics_values[:len(cleaned_boundaries)]
                    final_semantics_faces = list(range(len(cleaned_boundaries)))
                    
                    # Créer l'objet Building principal
                    self.city_objects[building_uuid] = {
                        "type": "Building",
                        "attributes": attributes,
                        "geometry": [{
                            "type": "MultiSurface",
                            "lod": "2.0",
                            "boundaries": cleaned_boundaries,
                            "semantics": {
                                "surfaces": [{"type": s} for s in final_semantics_values],
                                "values": final_semantics_faces
                            }
                        }],
                        "children": []
                    }
                    
                    # --- Créer les BuildingParts par type de surface (génération de nouveaux UUIDs) ---
                    
                    # BuildingPart pour les murs
                    if building_data['walls']:
                        wall_boundaries = self.validate_and_clean_boundaries(building_data['walls'])
                        if wall_boundaries:
                            # Génération de nouveaux UUIDs pour les BuildingParts (pas de pattern descriptif)
                            wall_part_uuid = str(uuid.uuid4())
                            child_uuids.append(wall_part_uuid)
                            
                            self.city_objects[wall_part_uuid] = {
                                "type": "BuildingPart",
                                "attributes": {
                                    "surfaceType": "WallSurface"
                                },
                                "geometry": [{
                                    "type": "MultiSurface",
                                    "lod": "2.0",
                                    "boundaries": wall_boundaries,
                                    "semantics": {
                                        "surfaces": [{"type": "WallSurface"} for _ in wall_boundaries],
                                        "values": list(range(len(wall_boundaries)))
                                    }
                                }],
                                "parents": [building_uuid]
                            }
                    
                    # BuildingPart pour les toits
                    if building_data['roofs']:
                        roof_boundaries = self.validate_and_clean_boundaries(building_data['roofs'])
                        if roof_boundaries:
                            # Génération de nouveaux UUIDs pour les BuildingParts (pas de pattern descriptif)
                            roof_part_uuid = str(uuid.uuid4())
                            child_uuids.append(roof_part_uuid)
                            
                            self.city_objects[roof_part_uuid] = {
                                "type": "BuildingPart",
                                "attributes": {
                                    "surfaceType": "RoofSurface"
                                },
                                "geometry": [{
                                    "type": "MultiSurface",
                                    "lod": "2.0",
                                    "boundaries": roof_boundaries,
                                    "semantics": {
                                        "surfaces": [{"type": "RoofSurface"} for _ in roof_boundaries],
                                        "values": list(range(len(roof_boundaries)))
                                    }
                                }],
                                "parents": [building_uuid]
                            }
                    
                    # BuildingPart pour les sols
                    if building_data['floors']:
                        floor_boundaries = self.validate_and_clean_boundaries(building_data['floors'])
                        if floor_boundaries:
                            # Génération de nouveaux UUIDs pour les BuildingParts (pas de pattern descriptif)
                            floor_part_uuid = str(uuid.uuid4())
                            child_uuids.append(floor_part_uuid)
                            
                            self.city_objects[floor_part_uuid] = {
                                "type": "BuildingPart",
                                "attributes": {
                                    "surfaceType": "GroundSurface"
                                },
                                "geometry": [{
                                    "type": "MultiSurface",
                                    "lod": "2.0",
                                    "boundaries": floor_boundaries,
                                    "semantics": {
                                        "surfaces": [{"type": "GroundSurface"} for _ in floor_boundaries],
                                        "values": list(range(len(floor_boundaries)))
                                    }
                                }],
                                "parents": [building_uuid]
                            }
                    
                    # Mettre à jour la liste des enfants du Building
                    self.city_objects[building_uuid]["children"] = child_uuids
                    self.stats['buildings_created'] += 1
        
        print(f" {self.stats['buildings_created']} bâtiments créés avec BuildingParts (UUID aléatoires)")
        print(f" Surfaces regroupées: {total_wall_surfaces} murs, {total_roof_surfaces} toits, {total_floor_surfaces} sols")
        print(f" Total surfaces: {total_wall_surfaces + total_roof_surfaces + total_floor_surfaces}")
        
        # Vérifier la conservation de la sémantique
        total_semantic_objects = 0
        buildings_with_semantics = 0
        building_parts_with_semantics = 0
        
        for obj in self.city_objects.values():
            if obj['type'] == 'Building' and 'geometry' in obj:
                if any('semantics' in geom for geom in obj['geometry']):
                    buildings_with_semantics += 1
                    total_semantic_objects += 1
            elif obj['type'] == 'BuildingPart' and 'geometry' in obj:
                if any('semantics' in geom for geom in obj['geometry']):
                    building_parts_with_semantics += 1
                    total_semantic_objects += 1
        
        print(f" ✅ Sémantique conservée :")
        print(f"    • Buildings avec sémantique: {buildings_with_semantics}")
        print(f"    • BuildingParts avec sémantique: {building_parts_with_semantics}")
        print(f"    • Total objets avec sémantique: {total_semantic_objects}")
        print(f" ✅ Gestion UUID conforme à exportRhinoLocalRef2CityJSON.py :")
        print(f"    • Génération de nouveaux UUIDs pour tous les objets")
        print(f"    • Suppression des attributs UUID redondants")
        print(f"    • Relations parent-enfant avec UUIDs aléatoires")


def list_available_layers(data_path):
    """Liste toutes les couches disponibles dans une géodatabase ou GeoPackage"""
    data_path = Path(data_path)
    print(f"📋 Analyse des couches dans : {data_path}")
    
    # Détecter le format
    if data_path.suffix.lower() == '.gpkg':
        driver = ogr.GetDriverByName("GPKG")
        format_name = "GeoPackage"
    elif data_path.suffix.lower() == '.gdb' or data_path.is_dir():
        driver = ogr.GetDriverByName("OpenFileGDB")
        if not driver:
            driver = ogr.GetDriverByName("FileGDB")
        format_name = "Géodatabase"
    else:
        print(f"❌ Format non supporté: {data_path.suffix}")
        return
    
    if not driver:
        print(f"❌ Driver {format_name} non disponible")
        return
    
    datasource = driver.Open(str(data_path), 0)
    if not datasource:
        print(f"❌ Impossible d'ouvrir {data_path}")
        return
    
    total_layers = datasource.GetLayerCount()
    print(f"📁 {total_layers} couches trouvées dans le {format_name} :")
    print()
    
    for i in range(total_layers):
        layer = datasource.GetLayer(i)
        layer_name = layer.GetName()
        feature_count = layer.GetFeatureCount()
        
        # Analyser la géométrie
        geom_types = set()
        layer.ResetReading()
        sample_count = 0
        
        feature = layer.GetNextFeature()
        while feature and sample_count < 10:  # Échantillon de 10 features
            geometry = feature.GetGeometryRef()
            if geometry:
                geom_types.add(geometry.GetGeometryName())
            sample_count += 1
            feature = layer.GetNextFeature()
        
        geom_type_str = ", ".join(geom_types) if geom_types else "Aucune géométrie"
        
        print(f"   {i+1:2d}. {layer_name}")
        print(f"       • Features: {feature_count}")
        print(f"       • Géométrie: {geom_type_str}")
        
        # Vérifier la présence du champ OBJEKTART
        feature_defn = layer.GetLayerDefn()
        has_objektart = False
        for j in range(feature_defn.GetFieldCount()):
            if feature_defn.GetFieldDefn(j).GetName() == 'OBJEKTART':
                has_objektart = True
                break
        
        if has_objektart:
            print(f"       • ✅ Champ OBJEKTART présent (filtrage disponible)")
        else:
            print(f"       • ⚠️ Pas de champ OBJEKTART")
        
        print()
    
    datasource = None
    
    print("💡 Usage :")
    print(f"   # Traiter toutes les couches :")
    print(f"   python3 gdb_to_cityjson.py {data_path}")
    print()
    print(f"   # Traiter des couches spécifiques :")
    print(f"   python3 gdb_to_cityjson.py {data_path} --layers nom_couche1 nom_couche2")


def list_available_attributes(data_path, layer_filter=None):
    """Liste tous les attributs disponibles dans les couches d'une géodatabase ou GeoPackage"""
    data_path = Path(data_path)
    print(f"📋 Analyse des attributs dans : {data_path}")
    
    # Détecter le format et ouvrir
    if data_path.suffix.lower() == '.gpkg':
        driver = ogr.GetDriverByName("GPKG")
        format_name = "GeoPackage"
    elif data_path.suffix.lower() == '.gdb' or data_path.is_dir():
        driver = ogr.GetDriverByName("OpenFileGDB")
        if not driver:
            driver = ogr.GetDriverByName("FileGDB")
        format_name = "Géodatabase"
    else:
        print(f"❌ Format non supporté: {data_path.suffix}")
        return
    
    if not driver:
        print(f"❌ Driver {format_name} non disponible")
        return
    
    datasource = driver.Open(str(data_path), 0)
    if not datasource:
        print(f"❌ Impossible d'ouvrir {data_path}")
        return
    
    total_layers = datasource.GetLayerCount()
    layers_to_check = []
    
    # Filtrer les couches si nécessaire
    if layer_filter:
        layer_filter_set = set(layer_filter)
        for i in range(total_layers):
            layer = datasource.GetLayer(i)
            if layer.GetName() in layer_filter_set:
                layers_to_check.append(layer)
    else:
        # Prendre juste la première couche avec des features
        for i in range(total_layers):
            layer = datasource.GetLayer(i)
            if layer.GetFeatureCount() > 0:
                layers_to_check.append(layer)
                break
    
    if not layers_to_check:
        print(" ❌ Aucune couche avec des features trouvée")
        return
    
    all_attributes = set()
    
    for layer in layers_to_check:
        layer_name = layer.GetName()
        feature_count = layer.GetFeatureCount()
        
        print(f"\n 📁 Couche : {layer_name} ({feature_count} features)")
        
        # Analyser les attributs de la couche
        feature_defn = layer.GetLayerDefn()
        layer_attributes = []
        
        for i in range(feature_defn.GetFieldCount()):
            field_defn = feature_defn.GetFieldDefn(i)
            field_name = field_defn.GetName()
            field_type = field_defn.GetTypeName()
            layer_attributes.append((field_name, field_type))
            all_attributes.add(field_name)
        
        print(f"   Attributs ({len(layer_attributes)}) :")
        for field_name, field_type in sorted(layer_attributes):
            print(f"     • {field_name} ({field_type})")
        
        # Analyser quelques valeurs d'exemple
        layer.ResetReading()
        feature = layer.GetNextFeature()
        if feature:
            print(f"\n   Exemple de valeurs (première feature) :")
            for field_name, field_type in sorted(layer_attributes)[:10]:  # Limiter à 10 pour la lisibilité
                field_value = feature.GetField(field_name)
                if field_value is not None:
                    # Tronquer les valeurs trop longues
                    str_value = str(field_value)
                    if len(str_value) > 50:
                        str_value = str_value[:47] + "..."
                    print(f"     • {field_name}: {str_value}")
    
    # Résumé global
    print(f"\n 📊 Résumé global :")
    print(f"   • Couches analysées : {len(layers_to_check)}")
    print(f"   • Attributs uniques : {len(all_attributes)}")
    
    print(f"\n 💡 Exemples d'utilisation du filtrage :")
    print(f"   # Inclure uniquement certains attributs :")
    print(f"   --include-attrs UUID OBJEKTART SurfaceType")
    print(f"   ")
    print(f"   # Exclure des attributs spécifiques :")
    print(f"   --exclude-attrs HERKUNFT_MONAT REVISION_MONAT")
    print(f"   ")
    print(f"   # Utiliser des patterns (wildcards) :")
    print(f"   --attr-patterns '*ID*' 'UUID*' '*TYPE*'")
    print(f"   ")
    print(f"   # Combiner inclusion et exclusion :")
    print(f"   --include-attrs UUID OBJEKTART --exclude-attrs TEMP_FIELD")
    
    datasource = None


def main():
    """Fonction principale - Conversion directe GDB/GeoPackage vers CityJSON avec regroupement sémantique, BuildingParts et gestion optimisée des ID"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Convertit une géodatabase (.gdb) ou GeoPackage (.gpkg) en CityJSON avec regroupement sémantique des surfaces par bâtiment, création de BuildingParts et gestion optimisée des ID conforme CityJSON 2.0")
    parser.add_argument("data_path", help="Chemin vers la géodatabase (.gdb) ou GeoPackage (.gpkg)")
    parser.add_argument("-o", "--output", help="Fichier de sortie CityJSON", 
                       default="./output/buildings_semantic.city.json")
    parser.add_argument("--epsg", type=int, default=2056, 
                       help="Code EPSG du système de coordonnées de sortie (défaut: 2056)")
    parser.add_argument("--layers", nargs="*", help="Noms spécifiques des couches à traiter (par défaut: toutes)")
    parser.add_argument("--list-layers", action="store_true", help="Afficher la liste des couches disponibles et quitter")
    
    # Filtrage des attributs
    attr_group = parser.add_argument_group('Filtrage des attributs')
    attr_group.add_argument("--include-attrs", nargs="*", 
                           help="Attributs à inclure uniquement (si spécifié, seuls ces attributs seront inclus)")
    attr_group.add_argument("--exclude-attrs", nargs="*", 
                           help="Attributs à exclure")
    attr_group.add_argument("--attr-patterns", nargs="*", 
                           help="Patterns d'attributs à inclure (wildcards supportés, ex: '*ID*', 'UUID*')")
    attr_group.add_argument("--no-attrs", action="store_true",
                           help="Ne conserver aucun attribut (seulement la géométrie)")
    attr_group.add_argument("--list-attrs", action="store_true", 
                           help="Lister tous les attributs disponibles dans la première couche et quitter")
    
    args = parser.parse_args()
    
    # Vérifier que le fichier existe
    data_path = Path(args.data_path)
    if not data_path.exists():
        print(f"❌ Fichier non trouvé : {args.data_path}")
        sys.exit(1)
    
    # Vérifier que c'est un format supporté
    if not (data_path.suffix.lower() in ['.gdb', '.gpkg'] or data_path.is_dir()):
        print(f"❌ Format non supporté. Le fichier doit être une géodatabase (.gdb) ou un GeoPackage (.gpkg) : {data_path.suffix}")
        sys.exit(1)
    
    # Si demandé, lister les couches et quitter
    if args.list_layers:
        list_available_layers(args.data_path)
        sys.exit(0)
    
    # Si demandé, lister les attributs et quitter
    if args.list_attrs:
        list_available_attributes(args.data_path, args.layers)
        sys.exit(0)
    
    # Construire le filtre d'attributs
    attribute_filter = {}
    
    # Si --no-attrs est spécifié, ne conserver aucun attribut
    if args.no_attrs:
        attribute_filter['include'] = []  # Liste vide = aucun attribut
        print("🚫 Mode sans attributs activé")
    else:
        if args.include_attrs:
            attribute_filter['include'] = args.include_attrs
        if args.exclude_attrs:
            attribute_filter['exclude'] = args.exclude_attrs
        if args.attr_patterns:
            attribute_filter['patterns'] = args.attr_patterns
    
    print(f"🔄 Conversion GDB/GeoPackage vers CityJSON avec regroupement sémantique, BuildingParts et ID optimisés")
    print(f"📁 Source: {args.data_path}")
    print(f"📄 Destination: {args.output}")
    print(f"🌐 Système de coordonnées: EPSG:{args.epsg}")
    
    if args.layers:
        print(f"🎯 Couches spécifiées: {', '.join(args.layers)}")
    else:
        print(f"📊 Traitement de toutes les couches disponibles")
    
    if attribute_filter:
        print(f"🔍 Filtrage des attributs activé:")
        if 'include' in attribute_filter:
            print(f"   • Inclure uniquement: {', '.join(attribute_filter['include'])}")
        if 'exclude' in attribute_filter:
            print(f"   • Exclure: {', '.join(attribute_filter['exclude'])}")
        if 'patterns' in attribute_filter:
            print(f"   • Patterns: {', '.join(attribute_filter['patterns'])}")
    
    print(f"⚙️ Optimisations CityJSON 2.0: ID Building = clé CityObjects, BuildingPart = UUID aléatoire")
    
    # Effectuer la conversion avec regroupement intégré
    converter = GeoDataToCityJSON(args.data_path, args.output, args.epsg, 
                                 layer_filter=args.layers, attribute_filter=attribute_filter)
    success = converter.convert()
    
    if success:
        print(f"\n✅ Conversion terminée avec succès !")
        print(f"📄 Fichier CityJSON créé: {args.output}")
        print(f"🏗️ Bâtiments regroupés avec classification sémantique des surfaces")
        print(f"🧱 BuildingParts créés par type de surface (murs, toits, sols)")
        print(f"🔗 Relations parent-enfant établies entre Buildings et BuildingParts")
        print(f"🆔 Gestion des ID optimisée conforme CityJSON 2.0 (suppression des redondances)")
        if attribute_filter:
            print(f"🔍 Filtrage des attributs appliqué")
    else:
        print(f"\n❌ Échec de la conversion")
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
