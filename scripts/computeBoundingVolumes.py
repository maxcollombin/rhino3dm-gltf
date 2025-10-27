import geopandas as gpd
import rasterio
import rasterio.features
import numpy as np
from shapely.geometry import mapping, Polygon
import warnings
import json
import uuid
from pathlib import Path

warnings.filterwarnings("ignore", category=UserWarning)

# --- charger données ---
secteurs = gpd.read_file("secteurs.gpkg", layer="secteurs")

# Charger toutes les surfaces de limitation Takeoff_climb_surface
# (pour avoir plus de chances d'intersection avec les secteurs)
filter_query = "\"SurfaceType\" = 'Takeoff_climb_surface'"

print(f"🔍 Chargement des surfaces avec filtre: {filter_query}")

plancher = gpd.read_file("hindernisbegrenzungsflaechen-kataster_2056.gpkg", 
                        layer="CadastreOfObstacleLimitationSurfaces_V2_CadastreOfObstacleLimitationSurfaces_WithLatestModification_Ols",
                        where=filter_query)

print(f"📄 Surfaces de limitation chargées: {len(plancher)}")
print(f"📍 Zone des secteurs: {secteurs.total_bounds}")
print(f"📍 Zone des surfaces: {plancher.total_bounds}")

mnt = rasterio.open("swissALTI3D.tif")

# reprojecter si nécessaire
if secteurs.crs != mnt.crs:
    secteurs = secteurs.to_crs(mnt.crs)
if plancher.crs != mnt.crs:
    plancher = plancher.to_crs(mnt.crs)

def interpolate_altitude_on_triangle(point_x, point_y, triangle_coords_3d):
    """
    Interpolation bilinéaire de l'altitude d'un point sur un triangle 3D
    en utilisant les coordonnées barycentriques.
    
    Args:
        point_x, point_y: Coordonnées du point à interpoler
        triangle_coords_3d: Liste de 3 tuples (x, y, z) définissant le triangle
        
    Returns:
        float: Altitude interpolée
    """
    if len(triangle_coords_3d) < 3:
        return None
    
    # Points du triangle
    x1, y1, z1 = triangle_coords_3d[0]
    x2, y2, z2 = triangle_coords_3d[1] 
    x3, y3, z3 = triangle_coords_3d[2]
    
    # Calculer les coordonnées barycentriques
    denom = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    
    if abs(denom) < 1e-10:  # Triangle dégénéré
        return (z1 + z2 + z3) / 3  # Moyenne des altitudes
    
    a = ((y2 - y3) * (point_x - x3) + (x3 - x2) * (point_y - y3)) / denom
    b = ((y3 - y1) * (point_x - x3) + (x1 - x3) * (point_y - y3)) / denom
    c = 1 - a - b
    
    # Interpolation
    interpolated_z = a * z1 + b * z2 + c * z3
    
    return interpolated_z

# utilitaire : extraire altitude minimale d'une géométrie 3D pour un secteur donné
def get_minimum_ceiling_altitude(sector_geom, ceiling_surfaces):
    """
    Calcule l'altitude minimale du plafond en interpolant l'altitude de chaque point
    du secteur sur les surfaces de limitation d'obstacles aéronautiques.
    
    Args:
        sector_geom: Géométrie du secteur (Polygon ou MultiPolygon)
        ceiling_surfaces: GeoDataFrame des surfaces de plafond intersectantes
    
    Returns:
        float: Altitude minimale du plafond pour ce secteur
    """
    interpolated_altitudes = []
    
    # Gérer les MultiPolygon en prenant le plus grand
    from shapely.geometry import MultiPolygon
    if isinstance(sector_geom, MultiPolygon):
        sector_geom = max(sector_geom.geoms, key=lambda p: p.area)
    
    if not hasattr(sector_geom, 'exterior'):
        return None
    
    # Extraire tous les points du contour du secteur
    sector_points = list(sector_geom.exterior.coords[:-1])  # Enlever le dernier point dupliqué
    
    print(f"   🔍 Interpolation sur {len(sector_points)} points du secteur avec {len(ceiling_surfaces)} surface(s)")
    
    # Pour chaque point du secteur, calculer l'altitude interpolée
    for i, (x, y) in enumerate(sector_points):
        # Chercher dans quelle surface ce point se trouve
        from shapely.geometry import Point
        point_2d = Point(x, y)
        
        point_altitude = None
        for _, surface in ceiling_surfaces.iterrows():
            surface_geom = surface.geometry
            
            if surface_geom.contains(point_2d) or surface_geom.intersects(point_2d):
                # Récupérer les coordonnées 3D de la surface
                surface_coords_3d = []
                if hasattr(surface_geom, 'exterior') and surface_geom.exterior is not None:
                    for coord in surface_geom.exterior.coords[:-1]:  # Enlever dernier point dupliqué
                        if len(coord) >= 3:
                            surface_coords_3d.append(coord)
                
                if len(surface_coords_3d) >= 3:
                    # Interpoler l'altitude
                    interpolated_alt = interpolate_altitude_on_triangle(x, y, surface_coords_3d)
                    point_altitude = interpolated_alt
                    break
        
        if point_altitude is not None:
            interpolated_altitudes.append(point_altitude)
            if i < 3:  # Afficher quelques exemples
                print(f"   📏 Point {i+1} ({x:.1f}, {y:.1f}): altitude interpolée = {point_altitude:.2f}m")
    
    if interpolated_altitudes:
        global_min = min(interpolated_altitudes)
        global_max = max(interpolated_altitudes)
        print(f"   🎯 Altitude du plafond: {global_min:.2f}m à {global_max:.2f}m → MIN = {global_min:.2f}m (sur {len(interpolated_altitudes)} points)")
        return {"min": global_min, "max": global_max, "all_altitudes": interpolated_altitudes}
    else:
        print(f"   ⚠️ Aucune altitude interpolée trouvée pour les points du secteur")
        return None

# calculer pixel area (m^2) à partir de l'affine transform
transform = mnt.transform
pixel_area = abs(transform.a * transform.e)  # a = xres, e = yres (négatif souvent)

# helper : extraire valeurs MNT (masque) pour une géométrie (retourne 1D array valid)
def sample_mnt_values_for_geom(raster, geom):
    # masque True à l'extérieur, invert True => mask is area of interest
    out_shape = (raster.height, raster.width)
    mask = rasterio.features.geometry_mask([mapping(geom)],
                                           out_shape=out_shape,
                                           transform=raster.transform,
                                           invert=True)
    arr = raster.read(1, masked=True)
    # mask True means inside geometry; arr is masked array so we index with mask
    vals = arr[mask]
    # retirer nodata
    if hasattr(arr, "mask"):
        # vals est un numpy.ma.MaskedArray already; convert to normal array removing masked
        if isinstance(vals, np.ma.MaskedArray):
            vals = vals.compressed()
    return np.array(vals)

def interpolate_mnt_altitude(x, y, mnt_raster):
    """
    Interpoler l'altitude du MNT pour un point donné
    """
    try:
        # Convertir les coordonnées du monde en indices de pixel
        row, col = mnt_raster.index(x, y)
        
        # Vérifier que les indices sont dans les limites
        if 0 <= row < mnt_raster.height and 0 <= col < mnt_raster.width:
            altitude = mnt_raster.read(1)[row, col]
            if not np.isnan(altitude) and altitude != mnt_raster.nodata:
                return float(altitude)
    except:
        pass
    return None

def create_3d_volume_geometry(geom_2d, z_bottom_ref, z_top, ref_point=(2592980.685, 1119281.703, 483.8), mnt_raster=None, ceiling_surfaces=None):
    """
    Crée une géométrie 3D Solid CityJSON à partir d'un polygone 2D et de hauteurs
    
    Args:
        geom_2d: Polygone 2D Shapely (Polygon ou MultiPolygon)
        z_bottom_ref: Altitude de référence du sol (float) ou None pour utiliser le MNT
        z_top: Altitude du plafond (float ou array pour variation)
        ref_point: Point de référence pour transformation locale
        mnt_raster: Raster MNT pour interpolation des altitudes du sol
        ceiling_surfaces: GeoDataFrame des surfaces de plafond pour interpolation
    
    Returns:
        dict: Géométrie CityJSON avec vertices et boundaries
    """
    from shapely.geometry import MultiPolygon
    
    # Gérer les MultiPolygon en prenant le plus grand
    if isinstance(geom_2d, MultiPolygon):
        # Prendre le polygone avec la plus grande aire
        largest_poly = max(geom_2d.geoms, key=lambda p: p.area)
        geom_2d = largest_poly
        print(f"   📐 MultiPolygon détecté, utilisation du plus grand polygone (aire: {geom_2d.area:.0f} m²)")
    
    if not hasattr(geom_2d, 'exterior'):
        print(f"   ❌ Géométrie non supportée: {type(geom_2d)}")
        return None
    
    ref_x, ref_y, ref_z = ref_point
    
    # Extraire les coordonnées du contour extérieur
    exterior_coords = list(geom_2d.exterior.coords[:-1])  # Enlever le dernier point dupliqué
    
    if len(exterior_coords) < 3:
        return None
    
    # Créer les vertices en coordonnées locales
    vertices = []
    vertex_index_map = {}
    
    def add_vertex(x, y, z):
        # Transformation en coordonnées locales
        x_local = x - ref_x
        y_local = y - ref_y
        z_local = z - ref_z
        
        # Clé pour éviter les doublons
        key = (round(x_local, 3), round(y_local, 3), round(z_local, 3))
        
        if key in vertex_index_map:
            return vertex_index_map[key]
        
        idx = len(vertices)
        vertices.append([x_local, y_local, z_local])
        vertex_index_map[key] = idx
        return idx
    
    # Créer les vertices pour le sol et le plafond
    bottom_indices = []
    top_indices = []
    
    for i, (x, y) in enumerate(exterior_coords):
        # Sol - interpoler depuis le MNT si disponible
        if mnt_raster is not None:
            z_bot = interpolate_mnt_altitude(x, y, mnt_raster)
            if z_bot is None:
                z_bot = z_bottom_ref  # Fallback sur l'altitude de référence
        else:
            if isinstance(z_bottom_ref, (list, np.ndarray)):
                z_bot = z_bottom_ref[min(i, len(z_bottom_ref)-1)]
            else:
                z_bot = z_bottom_ref
        
        # Plafond - interpoler depuis les surfaces de plafond si disponibles
        if ceiling_surfaces is not None and len(ceiling_surfaces) > 0:
            # Interpoler l'altitude du plafond pour ce point spécifique
            from shapely.geometry import Point
            point_2d = Point(x, y)
            z_tp = None
            
            for _, surface in ceiling_surfaces.iterrows():
                surface_geom = surface.geometry
                if surface_geom.contains(point_2d) or surface_geom.intersects(point_2d):
                    # Récupérer les coordonnées 3D de la surface
                    surface_coords_3d = []
                    if hasattr(surface_geom, 'exterior') and surface_geom.exterior is not None:
                        for coord in surface_geom.exterior.coords[:-1]:
                            if len(coord) >= 3:
                                surface_coords_3d.append(coord)
                    
                    if len(surface_coords_3d) >= 3:
                        # Interpoler l'altitude
                        z_tp = interpolate_altitude_on_triangle(x, y, surface_coords_3d)
                        break
            
            # Fallback si pas d'interpolation possible
            if z_tp is None:
                if isinstance(z_top, (list, np.ndarray)):
                    z_tp = z_top[min(i, len(z_top)-1)]
                else:
                    z_tp = z_top
        else:
            if isinstance(z_top, (list, np.ndarray)):
                z_tp = z_top[min(i, len(z_top)-1)]
            else:
                z_tp = z_top
        
        bottom_idx = add_vertex(x, y, z_bot)
        top_idx = add_vertex(x, y, z_tp)
        
        bottom_indices.append(bottom_idx)
        top_indices.append(top_idx)
    
    # Construire les faces du solide
    boundaries = []
    
    # Face du sol (orientée vers le bas - normale négative Z)
    bottom_face = list(reversed(bottom_indices))  # Inverser pour normale correcte
    boundaries.append([bottom_face])
    
    # Face du plafond (orientée vers le haut - normale positive Z)
    top_face = top_indices.copy()
    boundaries.append([top_face])
    
    # Faces latérales
    n_points = len(bottom_indices)
    for i in range(n_points):
        next_i = (i + 1) % n_points
        
        # Face latérale (quadrilatère)
        # Ordre des points pour normale extérieure
        wall_face = [
            bottom_indices[i],
            bottom_indices[next_i], 
            top_indices[next_i],
            top_indices[i]
        ]
        boundaries.append([wall_face])
    
    return {
        'vertices': vertices,
        'boundaries': boundaries
    }

def create_cityjson_object(sector_geom, z_bottom_ref, z_top, sector_id, object_type="GenericCityObject", z_top_max=None, mnt_raster=None, ceiling_surfaces=None):
    """
    Crée un objet CityJSON complet pour un secteur
    """
    geom_3d = create_3d_volume_geometry(sector_geom, z_bottom_ref, z_top, mnt_raster=mnt_raster, ceiling_surfaces=ceiling_surfaces)
    
    if not geom_3d:
        return None, None, None
    
    obj_id = str(uuid.uuid4())
    
    # Utiliser z_top_max si fourni, sinon z_top
    z_top_value = z_top_max if z_top_max is not None else z_top
    
    # Calculer les valeurs de sol et plafond correctement
    # Les vertices sont organisés par paires : sol puis plafond pour chaque point du contour
    ref_x, ref_y, ref_z = 2592980.685, 1119281.703, 483.8
    
    # Identifier les vertices du sol et du plafond
    # Les faces du sol et plafond sont les 2 premières dans boundaries
    bottom_face_indices = geom_3d['boundaries'][0][0]  # Face du sol
    top_face_indices = geom_3d['boundaries'][1][0]     # Face du plafond
    
    # Extraire les altitudes du sol
    bottom_altitudes = []
    for idx in bottom_face_indices:
        if idx < len(geom_3d['vertices']):
            vertex = geom_3d['vertices'][idx]
            global_z = vertex[2] + ref_z
            bottom_altitudes.append(global_z)
    
    # Extraire les altitudes du plafond  
    top_altitudes = []
    for idx in top_face_indices:
        if idx < len(geom_3d['vertices']):
            vertex = geom_3d['vertices'][idx]
            global_z = vertex[2] + ref_z
            top_altitudes.append(global_z)
    
    # Calculer les valeurs min/max
    if bottom_altitudes:
        z_bottom_min = float(min(bottom_altitudes))
        z_bottom_max = float(max(bottom_altitudes))
    else:
        z_bottom_min = float(z_bottom_ref) if isinstance(z_bottom_ref, (int, float)) else 0.0
        z_bottom_max = z_bottom_min
        
    if top_altitudes:
        z_top_min_actual = float(min(top_altitudes))
        z_top_max_actual = float(max(top_altitudes))
    else:
        z_top_min_actual = float(z_top)
        z_top_max_actual = float(z_top_value)
    
    # Calculer la hauteur minimale (sécurité aéronautique)
    maxHeight = z_top_min_actual - z_bottom_min
    height_min = z_top_min_actual - z_bottom_max  # Hauteur minimale = plafond min - sol max
    height_max = z_top_max_actual - z_bottom_min  # Hauteur maximale = plafond max - sol min
    
    city_object = {
        "type": object_type,
        "attributes": {
            "id": sector_id,
            "terrainMin": round(z_bottom_min, 2),
            "terrainMax": round(z_bottom_max, 2),
            "topBoundaryMin": round(z_top_min_actual, 2),
            "topBoundaryMax": round(z_top_max_actual, 2),
            "maxHeight": round(maxHeight, 0),
            "area": round(float(sector_geom.area), 2)
        },
        "geometry": [{
            "type": "Solid",
            "lod": "2.0",
            "boundaries": [geom_3d['boundaries']]  # Un shell
        }]
    }
    
    return obj_id, city_object, geom_3d['vertices']

results = []
city_objects = {}
all_vertices = []
vertex_offset = 0

print(f"🔄 Traitement de {len(secteurs)} secteurs...")

for idx, row in secteurs.iterrows():
    geom = row.geometry
    # Utiliser l'attribut 'name' du GeoPackage comme ID, avec fallback sur sector_idx
    sector_id = row.get("name", row.get("id", f"sector_{idx}"))
    
    print(f"📍 Secteur {sector_id}: {geom.area:.0f} m²")
    
    # 1) Extraire z plafond : intersection avec plancher
    planchers_locaux = plancher[plancher.intersects(geom)]
    if planchers_locaux.empty:
        print(f"   ⚠️ Pas de plafond trouvé, ignoré")
        continue
    
    # Unionner et intersecter pour la portion au-dessus du secteur
    try:
        plafond_inter = planchers_locaux.union_all().intersection(geom)
    except AttributeError:
        # Fallback pour versions plus anciennes de geopandas
        plafond_inter = planchers_locaux.unary_union.intersection(geom)
    
    # Calculer l'altitude minimale du plafond pour ce secteur
    altitude_data = get_minimum_ceiling_altitude(geom, planchers_locaux)
    if altitude_data is None:
        print(f"   ⚠️ Impossible d'extraire l'altitude du plafond, ignoré")
        continue
    
    z_plafond_min = altitude_data["min"]
    z_plafond_max = altitude_data["max"]
    
    # 2) Échantillonnage du MNT sous la géométrie pour altitude du sol
    vals = sample_mnt_values_for_geom(mnt, geom)
    if vals.size == 0:
        print(f"   ⚠️ Aucune donnée MNT trouvée, ignoré")
        continue
    
    # Utiliser z_mean comme altitude de référence du sol
    z_sol = float(np.mean(vals))
    z_min = float(np.min(vals))
    z_max = float(np.max(vals))
    
    print(f"   📏 Sol: {z_min:.1f}m - {z_max:.1f}m (moy: {z_sol:.1f}m)")
    print(f"   🏢 Plafond: {z_plafond_min:.2f}m - {z_plafond_max:.2f}m")
    print(f"   📐 Hauteur: {z_plafond_min - z_sol:.1f}m (min)")
    
    # 3) Créer l'objet CityJSON 3D (utiliser l'altitude minimale pour la sécurité)
    obj_id, city_object, vertices = create_cityjson_object(
        geom, z_sol, z_plafond_min, sector_id, "GenericCityObject", z_plafond_max, mnt_raster=mnt, ceiling_surfaces=planchers_locaux
    )
    
    if obj_id and city_object:
        # Ajuster les indices des vertices dans la géométrie
        adjusted_boundaries = []
        for shell in city_object['geometry'][0]['boundaries']:
            adjusted_shell = []
            for face in shell:
                adjusted_face = []
                for ring in face:
                    adjusted_ring = [idx + vertex_offset for idx in ring]
                    adjusted_face.append(adjusted_ring)
                adjusted_shell.append(adjusted_face)
            adjusted_boundaries.append(adjusted_shell)
        
        city_object['geometry'][0]['boundaries'] = adjusted_boundaries
        city_objects[obj_id] = city_object
        all_vertices.extend(vertices)
        vertex_offset += len(vertices)
        
        # Statistiques pour résultats - calcul correct du volume et des hauteurs
        z_bottom_min = float(np.min(vals)) if vals.size > 0 else z_sol
        z_bottom_max = float(np.max(vals)) if vals.size > 0 else z_sol
        
        # Hauteur minimale pour sécurité aéronautique (plafond min - sol max)
        hauteur_min = z_plafond_min - z_bottom_max
        hauteur_max = z_plafond_max - z_bottom_min
        
        # Volume approximatif basé sur hauteur minimale (conservateur)
        volume_m3 = hauteur_min * geom.area
        
        results.append({
            "id": sector_id,
            "object_id": obj_id,
            "area_m2": round(float(geom.area), 2),
            "z_sol_min": round(z_bottom_min, 2),
            "z_sol_max": round(z_bottom_max, 2),
            "z_sol_moy": round(z_sol, 2),
            "z_top_min": round(z_plafond_min, 2),
            "z_top_max": round(z_plafond_max, 2),
            "hauteur_min": round(hauteur_min, 2),
            "hauteur_max": round(hauteur_max, 2),
            "volume_m3": round(volume_m3, 0),
            "vertices_count": len(vertices),
            "faces_count": len(city_object['geometry'][0]['boundaries'][0])
        })
        
        print(f"   ✅ Volume créé: {volume_m3:.0f} m³ (hauteur: {hauteur_min:.1f}m-{hauteur_max:.1f}m, {len(vertices)} vertices)")
    else:
        print(f"   ❌ Échec création géométrie 3D")

print(f"\n📊 Résumé: {len(city_objects)} volumes 3D créés sur {len(secteurs)} secteurs")

print(f"\n📊 Résumé: {len(city_objects)} volumes 3D créés sur {len(secteurs)} secteurs")

# Créer le document CityJSON
if city_objects:
    # Point de référence (même que le script gdb_to_cityjson.py)
    ref_x, ref_y, ref_z = 2592980.685, 1119281.703, 483.8
    
    # Calculer l'étendue géographique en coordonnées locales
    if all_vertices:
        xs, ys, zs = zip(*all_vertices)
        geographical_extent = [
            min(xs), min(ys), min(zs),  # min x, y, z
            max(xs), max(ys), max(zs)   # max x, y, z
        ]
    else:
        geographical_extent = [0, 0, 0, 0, 0, 0]
    
    cityjson = {
        "type": "CityJSON",
        "version": "2.0",
        "CityObjects": city_objects,
        "vertices": all_vertices,
        "metadata": {
            "referenceSystem": "https://www.opengis.net/def/crs/EPSG/0/2056",
            "geographicalExtent": geographical_extent,
            "presentLoDs": ["2.0"],
            "datasetTitle": "Volumes de limitation d'obstacles aéronautiques"
        },
        "transform": {
            "scale": [1.0, 1.0, 1.0],
            "translate": [ref_x, ref_y, ref_z]
        }
    }
    
    # Sauvegarder le fichier CityJSON
    output_file = "output/secteurs3d.city.json"
    Path("output").mkdir(exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(cityjson, f, indent=2, ensure_ascii=False)
    
    print(f"\n✅ Fichier CityJSON créé: {output_file}")
    print(f"📐 Transformation locale: translate=[{ref_x}, {ref_y}, {ref_z}]")
    print(f"📏 Étendue géographique: {geographical_extent}")
    print(f"🔢 Statistiques:")
    print(f"   • Objets CityJSON: {len(city_objects)}")
    print(f"   • Vertices totaux: {len(all_vertices)}")
    
    # Statistiques détaillées
    total_volume = sum(r['volume_m3'] for r in results)
    total_area = sum(r['area_m2'] for r in results)
    avg_height_min = sum(r['hauteur_min'] for r in results) / len(results) if results else 0
    avg_height_max = sum(r['hauteur_max'] for r in results) / len(results) if results else 0
    
    print(f"   • Volume total: {total_volume:.0f} m³")
    print(f"   • Surface totale: {total_area:.0f} m²")
    print(f"   • Hauteur moyenne: {avg_height_min:.1f}m - {avg_height_max:.1f}m")
    
    # Afficher quelques exemples
    print(f"\n📋 Exemples de volumes créés:")
    for r in results[:3]:
        print(f"   • Secteur {r['id']}: {r['volume_m3']:.0f} m³ (hauteur: {r['hauteur_min']:.1f}m-{r['hauteur_max']:.1f}m)")

else:
    print("❌ Aucun volume 3D créé")

# Afficher les résultats détaillés si demandé
if results:
    print(f"\n📊 Résultats détaillés:")
    for r in results:
        print(f"   {r}")
