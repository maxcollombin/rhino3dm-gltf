import rasterio
import numpy as np
from scipy.spatial import Delaunay
from scipy.ndimage import gaussian_filter
import json
import argparse
import os
import uuid

def create_spatial_tiles(points, tile_size):
    """Diviser les points en tuiles spatiales pour CityJSONSeq"""
    print(f"🗂️ Division en tuiles spatiales ({tile_size}m x {tile_size}m)")
    
    # Calculer les limites
    min_x, min_y = np.min(points[:, 0]), np.min(points[:, 1])
    max_x, max_y = np.max(points[:, 0]), np.max(points[:, 1])
    
    # Calculer le nombre de tuiles
    n_tiles_x = int(np.ceil((max_x - min_x) / tile_size))
    n_tiles_y = int(np.ceil((max_y - min_y) / tile_size))
    
    print(f"   Grille de tuiles: {n_tiles_x} x {n_tiles_y} = {n_tiles_x * n_tiles_y} tuiles")
    
    tiles = {}
    
    for i, (x, y, z) in enumerate(points):
        # Déterminer la tuile
        tile_x = int((x - min_x) / tile_size)
        tile_y = int((y - min_y) / tile_size)
        
        # S'assurer que les indices sont dans les limites
        tile_x = min(tile_x, n_tiles_x - 1)
        tile_y = min(tile_y, n_tiles_y - 1)
        
        tile_key = (tile_x, tile_y)
        
        if tile_key not in tiles:
            tiles[tile_key] = []
        tiles[tile_key].append(i)
    
    print(f"   Tuiles non-vides: {len(tiles)}")
    
    return tiles, (n_tiles_x, n_tiles_y), (min_x, min_y)

def create_cityjsonseq_tile(points_subset, tile_id, ref_point, boundaries_subset):
    """Créer une tuile CityJSONSeq"""
    ref_x, ref_y, ref_z = ref_point
    
    # Transformation en coordonnées locales
    vertices_local = []
    for point in points_subset:
        x_local = point[0] - ref_x
        y_local = point[1] - ref_y  
        z_local = point[2] - ref_z
        vertices_local.append([x_local, y_local, z_local])
    
    # Calculer l'étendue géographique
    if vertices_local:
        xs_local, ys_local, zs_local = zip(*vertices_local)
        geographical_extent = [min(xs_local), min(ys_local), min(zs_local),
                              max(xs_local), max(ys_local), max(zs_local)]
    else:
        geographical_extent = [0, 0, 0, 0, 0, 0]
    
    # Créer l'objet CityJSON pour cette tuile
    tile_cityjson = {
        "type": "CityJSON",
        "version": "2.0",
        "transform": {
            "scale": [1.0, 1.0, 1.0],
            "translate": [ref_x, ref_y, ref_z]
        },
        "CityObjects": {
            f"TINRelief_tile_{tile_id}": {
                "type": "TINRelief",
                "geometry": [{
                    "type": "CompositeSurface",
                    "lod": "3.0",
                    "boundaries": boundaries_subset
                }]
            }
        },
        "vertices": vertices_local,
        "metadata": {
            "referenceSystem": "https://www.opengis.net/def/crs/EPSG/0/2056",
            "geographicalExtent": geographical_extent
        }
    }
    
    return tile_cityjson

def main():
    # --- Configuration des arguments ---
    parser = argparse.ArgumentParser(description='Convertir un GeoTIFF en TINRelief CityJSON')
    parser.add_argument('--input', '-i', 
                       required=True,
                       help='Fichier GeoTIFF d\'entrée (ex: mnt_5m.tif)')
    parser.add_argument('--output', '-o',
                       help='Fichier CityJSON de sortie (optionnel, généré automatiquement si non spécifié)')
    parser.add_argument('--output-format',
                       choices=['cityjson', 'cityjsonseq'],
                       default='cityjson',
                       help='Format de sortie: cityjson (fichier unique) ou cityjsonseq (streaming .jsonl)')
    parser.add_argument('--tile-size',
                       type=float,
                       default=1000.0,
                       help='Taille des tuiles en mètres pour CityJSONSeq (défaut: 1000m)')
    parser.add_argument('--object-id',
                       default=None,
                       help='ID personnalisé pour l\'objet TINRelief (génération UUID par défaut)')
    parser.add_argument('--reference-point',
                       nargs=3,
                       type=float,
                       default=[2592980.685, 1119281.703, 483.8],
                       metavar=('X', 'Y', 'Z'),
                       help='Point de référence pour la transformation locale (défaut: 2592980.685 1119281.703 483.8)')
    parser.add_argument('--no-transform',
                       action='store_true',
                       help='Désactiver la transformation en coordonnées locales')
    
    args = parser.parse_args()
    
    # --- Validation des paramètres ---
    input_tif = args.input
    if not os.path.exists(input_tif):
        print(f"❌ Erreur: Le fichier {input_tif} n'existe pas")
        return
    
    # Génération automatique du nom de sortie si non spécifié
    if args.output:
        output_file = args.output
    else:
        base_name = os.path.splitext(os.path.basename(input_tif))[0]
        if args.output_format == 'cityjsonseq':
            output_file = f"{base_name}_tinrelief.jsonl"
        else:
            output_file = f"{base_name}_tinrelief.city.json"
    
    # Génération de l'ID d'objet
    object_id = args.object_id if args.object_id else f"UUID_{uuid.uuid4()}"
    
    # Point de référence pour transformation
    ref_x, ref_y, ref_z = args.reference_point
    
    print(f"🔄 Traitement de {input_tif}")
    print(f" Format de sortie: {args.output_format}")
    if args.output_format == 'cityjsonseq':
        print(f"🗂️ Taille des tuiles: {args.tile_size}m")
    print(f" Sortie: {output_file}")
    print(f"🆔 ID objet: {object_id}")
    if not args.no_transform:
        print(f"📍 Point de référence: ({ref_x:.3f}, {ref_y:.3f}, {ref_z:.1f})")
    else:
        print(f"📍 Coordonnées absolues (pas de transformation)")

    # --- lecture du raster ---
    print(f"📖 Lecture du raster...")
    with rasterio.open(input_tif) as src:
        print(f"   Dimensions: {src.width} x {src.height} pixels")
        print(f"   Résolution: {src.res[0]:.2f} x {src.res[1]:.2f} m/pixel")
        print(f"   Système de coordonnées: {src.crs}")
        
        # Stocker les métadonnées pour la fin
        src_bounds = src.bounds
        src_crs = src.crs
        
        z = src.read(1)
        mask = z != src.nodata
        
        print(f"🔍 Création de la grille de points...")
        xs, ys = np.meshgrid(
            np.linspace(src.bounds.left, src.bounds.right, src.width),
            np.linspace(src.bounds.top, src.bounds.bottom, src.height)
        )
        xs, ys, zs = xs[mask], ys[mask], z[mask]
        
        print(f"   Points valides extraits: {len(xs):,} / {src.width * src.height:,}")
        print(f"   Altitude min/max: {np.min(zs):.1f}m / {np.max(zs):.1f}m")

    # --- Traitement direct sans sous-échantillonnage ---
    print(f"📍 Utilisation de tous les points")
    points = np.c_[xs, ys, zs]
    print(f"📍 Points traités: {len(points):,}")

    # --- triangulation ---
    print("🔺 Triangulation Delaunay en cours...")
    print(f"   Points d'entrée: {len(points):,}")
    tri = Delaunay(points[:, :2])
    print(f"   Triangles générés: {len(tri.simplices):,}")

    # --- construction CityJSON ---
    print("🏗️ Construction du CityJSON...")
    vertices = points.tolist()
    boundaries = tri.simplices.tolist()

    print(f"📐 Conversion des boundaries au format CityJSON...")
    # Convertir les boundaries en format CityJSON (triangles avec indices fermés)
    boundaries_cityjson = [[[triangle[0], triangle[1], triangle[2]]] for triangle in boundaries]
    print(f"   Boundaries converties: {len(boundaries_cityjson):,}")

    # Transformation en coordonnées locales si demandée
    if not args.no_transform:
        print(f"🌐 Transformation en coordonnées locales...")
        vertices_local = []
        for vertex in vertices:
            x_local = vertex[0] - ref_x
            y_local = vertex[1] - ref_y
            z_local = vertex[2] - ref_z
            vertices_local.append([x_local, y_local, z_local])
        
        # Calculer l'étendue géographique en coordonnées locales
        xs_local, ys_local, zs_local = zip(*vertices_local)
        geographical_extent = [min(xs_local), min(ys_local), min(zs_local), 
                             max(xs_local), max(ys_local), max(zs_local)]
        
        transform = {
            "scale": [1.0, 1.0, 1.0],
            "translate": [ref_x, ref_y, ref_z]
        }
        vertices_final = vertices_local
    else:
        # Coordonnées absolues
        xs, ys, zs = zip(*vertices)
        geographical_extent = [min(xs), min(ys), min(zs), 
                             max(xs), max(ys), max(zs)]
        
        transform = {
            "scale": [1.0, 1.0, 1.0],
            "translate": [0.0, 0.0, 0.0]
        }
        vertices_final = vertices

    cityjson = {
        "type": "CityJSON",
        "version": "2.0",
        "transform": transform,
        "CityObjects": {
            object_id: {
                "type": "TINRelief",
                "geometry": [{
                    "type": "CompositeSurface",
                    "lod": "3.0",
                    "boundaries": boundaries_cityjson
                }]
            }
        },
        "vertices": vertices_final,
        "metadata": {
            "referenceSystem": f"https://www.opengis.net/def/crs/EPSG/0/2056",
            "geographicalExtent": geographical_extent
        }
    }

    # --- sauvegarde selon le format choisi ---
    if args.output_format == 'cityjsonseq':
        print(f"�️ Génération CityJSONSeq par tuiles...")
        
        # Diviser en tuiles spatiales
        tiles, (n_tiles_x, n_tiles_y), (min_x, min_y) = create_spatial_tiles(points, args.tile_size)
        
        print(f"�💾 Sauvegarde vers {output_file} (format CityJSONSeq)...")
        os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
        
        with open(output_file, "w") as f:
            total_tiles = len(tiles)
            processed_tiles = 0
            
            for (tile_x, tile_y), point_indices in tiles.items():
                processed_tiles += 1
                
                if processed_tiles % max(1, total_tiles // 10) == 0:
                    print(f"   Tuiles traitées: {processed_tiles}/{total_tiles} ({processed_tiles/total_tiles*100:.1f}%)")
                
                # Points de cette tuile
                tile_points = points[point_indices]
                
                # Triangulation pour cette tuile
                if len(tile_points) >= 3:
                    try:
                        tile_tri = Delaunay(tile_points[:, :2])
                        # Convertir les indices en entiers Python pour la sérialisation JSON
                        tile_boundaries = [[[int(tri[0]), int(tri[1]), int(tri[2])]] for tri in tile_tri.simplices]
                        
                        # Créer l'objet CityJSON pour cette tuile
                        tile_id = f"{tile_x}_{tile_y}"
                        tile_cityjson = create_cityjsonseq_tile(
                            tile_points, tile_id, (ref_x, ref_y, ref_z), tile_boundaries
                        )
                        
                        # Écrire une ligne CityJSONSeq
                        f.write(json.dumps(tile_cityjson, separators=(',', ':')) + '\n')
                    except Exception as e:
                        print(f"   ⚠️ Erreur tuile {tile_id}: {e}")
        
        print(f"✅ TINRelief CityJSONSeq créé : {output_file}")
        print(f"📊 Statistiques finales:")
        print(f"   - Tuiles générées: {len(tiles)}")
        print(f"   - Points totaux: {len(points):,}")
        
    else:
        # Format CityJSON classique
        print(f"💾 Sauvegarde vers {output_file} (format CityJSON)...")
        os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
        with open(output_file, "w") as f:
            json.dump(cityjson, f, indent=2)
        
        print(f"✅ TINRelief CityJSON créé : {output_file}")
        print(f"📊 Statistiques finales:")
        print(f"   - Triangles: {len(boundaries):,}")
        print(f"   - Vertices: {len(vertices_final):,}")
    print(f"   - Altitude min/max: {np.min(zs):.1f}m / {np.max(zs):.1f}m")
    print(f"   - Étendue: {src_bounds.left:.1f}, {src_bounds.bottom:.1f} → {src_bounds.right:.1f}, {src_bounds.top:.1f}")
    if not args.no_transform:
        print(f"   - Coordonnées: Locales (transformées)")
        print(f"   - Référence: ({ref_x:.3f}, {ref_y:.3f}, {ref_z:.1f})")
    else:
        print(f"   - Coordonnées: Absolues (non transformées)")
    print(f"   - CityJSON: version 2.0")
    print(f"   - Géométrie: CompositeSurface (conforme CityJSON 2.0)")

if __name__ == "__main__":
    main()
