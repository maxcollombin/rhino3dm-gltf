import os
import json
import uuid
import rhino3dm
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN

input_file = "./R21_3d.3dm"
output_cityjson = "./output/NouveauxBatiments.city.json"

model = rhino3dm.File3dm.Read(input_file)
if not model:
    raise ValueError(f"Erreur : Impossible de charger {input_file}")
print(f" Modèle chargé : {input_file}")

# ------------------------------
# Fonctions utilitaires
# ------------------------------
def remove_consecutive_duplicates(ring):
    if len(ring) < 3:
        return []
    cleaned = []
    for i, vertex in enumerate(ring):
        next_vertex = ring[(i + 1) % len(ring)]
        if vertex != next_vertex:
            cleaned.append(vertex)
    if len(cleaned) >= 3 and cleaned[0] == cleaned[-1]:
        cleaned = cleaned[:-1]
    return cleaned if len(cleaned) >= 3 else []

def validate_face(ring):
    if len(ring) < 3:
        return False
    for i in range(len(ring)):
        if ring[i] == ring[(i + 1) % len(ring)]:
            return False
    return True

def project_to_plane(vertices_coords):
    if len(vertices_coords) < 4:
        return vertices_coords
    points = np.array(vertices_coords)
    centroid = np.mean(points, axis=0)
    centered_points = points - centroid
    pca = PCA(n_components=3)
    pca.fit(centered_points)
    normal = pca.components_[2]
    projected_centered = centered_points - np.outer(np.dot(centered_points, normal), normal)
    return (projected_centered + centroid).tolist()

def convert_to_mesh(geometry):
    if isinstance(geometry, rhino3dm.Mesh):
        return geometry
    if isinstance(geometry, rhino3dm.Extrusion):
        for mesh_type in (rhino3dm.MeshType.Render, rhino3dm.MeshType.Default):
            mesh = geometry.GetMesh(mesh_type)
            if mesh and len(mesh.Vertices) and len(mesh.Faces):
                mesh.Normals.ComputeNormals()
                mesh.Compact()
                return mesh
        brep = geometry.ToBrep(True) or geometry.ToBrep(False)
        if brep:
            return convert_to_mesh(brep)
        return None
    if isinstance(geometry, rhino3dm.Brep):
        meshes = []
        for face in geometry.Faces:
            mesh = face.GetMesh(rhino3dm.MeshType.Default)
            if mesh and len(mesh.Vertices) and len(mesh.Faces):
                meshes.append(mesh)
        if meshes:
            combined = rhino3dm.Mesh()
            for m in meshes:
                combined.Append(m)
            combined.Normals.ComputeNormals()
            combined.Compact()
            return combined
    return None

def classify_surface(vertices_coords):
    points = np.array(vertices_coords)
    if len(points) < 3:
        return "WallSurface"  # Fallback conforme à CityJSON 2.0
    centroid = np.mean(points, axis=0)
    centered_points = points - centroid
    pca = PCA(n_components=3)
    pca.fit(centered_points)
    normal = pca.components_[2]
    z_axis = np.array([0,0,1])
    vertical = np.abs(np.dot(normal, z_axis)) < 0.5
    horizontal = np.abs(np.dot(normal, z_axis)) > 0.9
    if horizontal:
        z_mean = np.mean(points[:,2])
        if z_mean < centroid[2]:
            return "GroundSurface"
        else:
            return "RoofSurface"
    elif vertical:
        return "WallSurface"
    else:
        return "ClosureSurface"  # Conforme à CityJSON 2.0

# ------------------------------
# Gestion des sommets et bâtiments
# ------------------------------
vertices_global = []
vertex_index_map = {}

def add_vertex(v):
    key = (round(v.X,6), round(v.Y,6), round(v.Z,6))
    if key in vertex_index_map:
        return vertex_index_map[key]
    idx = len(vertices_global)
    vertices_global.append([v.X, v.Y, v.Z])
    vertex_index_map[key] = idx
    return idx

# Collecte des meshes et calcul des barycentres
meshes_all = []
barycentres = []
for obj in model.Objects:
    layer = model.Layers[obj.Attributes.LayerIndex]
    full_path = layer.FullPath
    if not full_path.startswith("Bâtiments Nouveaux"):
        continue
    if "EXISTANT" in full_path or "DEMOLITION" in full_path:
        continue
    mesh = convert_to_mesh(obj.Geometry)
    if mesh:
        meshes_all.append(mesh)
        # Calcul barycentre
        verts = np.array([[v.X, v.Y, v.Z] for v in mesh.Vertices])
        bary = np.mean(verts, axis=0)
        barycentres.append(bary)

print(f"{len(meshes_all)} meshes collectés.")

# Clustering spatial des barycentres
if len(barycentres) == 0:
    raise ValueError("Aucun mesh détecté.")

barycentres_np = np.array(barycentres)
clustering = DBSCAN(eps=10.0, min_samples=1).fit(barycentres_np)  # eps=10m, à ajuster selon la densité
labels = clustering.labels_
print(f"{len(set(labels))} bâtiments détectés par clustering spatial.")

# Regroupement des meshes par cluster
buildings = {}
for mesh, label in zip(meshes_all, labels):
    buildings.setdefault(label, []).append(mesh)

print(f" {len(buildings)} bâtiments détectés.")

# ------------------------------
# Fonction de traitement des faces
# ------------------------------
def process_mesh_faces(mesh):
    shells = []
    semantics_values = []
    semantics_faces = []
    mesh.Normals.ComputeNormals()
    for face in mesh.Faces:
        face_vertices = list(face) if len(face)==4 else list(face[:3])
        ring = [add_vertex(mesh.Vertices[j]) for j in face_vertices]
        cleaned_ring = remove_consecutive_duplicates(ring)
        if not validate_face(cleaned_ring):
            continue
        vertices_coords = [vertices_global[idx] for idx in cleaned_ring]
        projected_coords = project_to_plane(vertices_coords)
        # classification
        v1, v2, v3 = np.array(projected_coords[:3])
        normal = np.cross(v2 - v1, v3 - v1)
        normal = normal / (np.linalg.norm(normal)+1e-9)
        nz = normal[2]
        if abs(nz) > 0.9:
            face_type = "RoofSurface" if nz>0 else "GroundSurface"
        else:
            face_type = "WallSurface"
        semantics_values.append(face_type)
        semantics_faces.append(len(shells))
        shells.append([[add_vertex(rhino3dm.Point3d(*c)) for c in projected_coords]])
    return shells, semantics_values, semantics_faces

# ------------------------------
# Construction du CityJSON
# ------------------------------
cityobjects = {}
for building_label, meshes in buildings.items():
    building_id = str(uuid.uuid4())
    print(f"Création du bâtiment cluster {building_label} avec UUID {building_id}")
    child_ids = []
    # --- Building ---
    all_shells = []
    all_semantics_values = []
    all_semantics_faces = []
    for mesh in meshes:
        shells, s_values, s_faces = process_mesh_faces(mesh)
        all_shells.extend(shells)
        all_semantics_values.extend(s_values)
        all_semantics_faces.extend(s_faces)
    # --- BuildingParts ---
    mesh_data = []
    for i, mesh in enumerate(meshes):
        z_coords = [vertex.Z for vertex in mesh.Vertices]
        avg_z = sum(z_coords) / len(z_coords) if z_coords else 0
        mesh_data.append((mesh, avg_z, i))
    mesh_data.sort(key=lambda x: x[1])
    total_floors = len(mesh_data)
    all_z_coords = []
    for mesh, _, _ in mesh_data:
        all_z_coords.extend([vertex.Z for vertex in mesh.Vertices])
    min_z = min(all_z_coords) if all_z_coords else 0
    max_z = max(all_z_coords) if all_z_coords else 0
    building_height = max_z - min_z
    cityobjects[building_id] = {
        "type": "Building",
        "attributes": {
            "buildingStatus": "planned",
            "developmentPhase": 0,
            "storeysAboveGround": total_floors,
            "computedHeight": int(round(building_height))
        },
        "geometry": [{
            "type": "MultiSurface",
            "lod": "2.0",
            "boundaries": all_shells,
            "semantics": {
                "surfaces": [{"type": s} for s in all_semantics_values],
                "values": all_semantics_faces
            }
        }],
        "children": []
    }
    for order_index, (mesh, avg_z, original_index) in enumerate(mesh_data):
        part_id = str(uuid.uuid4())
        print(f"  Création étage {order_index+1} pour bâtiment cluster {building_label} (parent: {building_id})")
        child_ids.append(part_id)
        shells, s_values, s_faces = process_mesh_faces(mesh)
        cityobjects[part_id] = {
            "type": "BuildingPart",
            "attributes": {
                "floorOrder": order_index + 1
            },
            "geometry": [{
                "type": "MultiSurface",
                "lod": "2.0",
                "boundaries": shells,
                "semantics": {
                    "surfaces": [{"type": s} for s in s_values],
                    "values": s_faces
                }
            }],
            "parents": [building_id]
        }
    cityobjects[building_id]["children"] = child_ids

# ------------------------------
# Nettoyage des vertices non utilisés
# ------------------------------
def clean_unused_vertices(cityobjects, vertices_global):
    """Supprime les vertices non utilisés et réindexe les géométries"""
    used_vertices = set()
    
    # Collecter tous les indices de vertices utilisés
    for city_obj in cityobjects.values():
        if "geometry" in city_obj:
            for geom in city_obj["geometry"]:
                if "boundaries" in geom:
                    for shell in geom["boundaries"]:
                        for ring in shell:
                            for vertex_idx in ring:
                                used_vertices.add(vertex_idx)
    
    # Créer la nouvelle liste de vertices et le mapping
    old_to_new = {}
    new_vertices = []
    new_idx = 0
    
    for old_idx in range(len(vertices_global)):
        if old_idx in used_vertices:
            old_to_new[old_idx] = new_idx
            new_vertices.append(vertices_global[old_idx])
            new_idx += 1
    
    # Réindexer toutes les géométries
    for city_obj in cityobjects.values():
        if "geometry" in city_obj:
            for geom in city_obj["geometry"]:
                if "boundaries" in geom:
                    for shell in geom["boundaries"]:
                        for ring in shell:
                            for i in range(len(ring)):
                                ring[i] = old_to_new[ring[i]]
    
    return new_vertices

# Nettoyer les vertices non utilisés
vertices_cleaned = clean_unused_vertices(cityobjects, vertices_global)

# Appliquer la transformation locale aux vertices nettoyés
translate_x = 2592980.685
translate_y = 1119281.703
translate_z = 483.8

vertices_local = []
for vertex in vertices_cleaned:
    # Les coordonnées Rhino sont déjà locales, nous les gardons telles quelles
    # et ajustons seulement l'altitude par rapport à la référence
    x_local = vertex[0]  # Garder tel quel
    y_local = vertex[1]  # Garder tel quel  
    z_local = vertex[2] - translate_z  # Soustraire la référence d'altitude
    vertices_local.append([x_local, y_local, z_local])

# ------------------------------
# Écriture CityJSON
# ------------------------------
xs, ys, zs = zip(*vertices_local)
geographical_extent = [min(xs), min(ys), min(zs), max(xs), max(ys), max(zs)]

cityjson = {
    "type": "CityJSON",
    "version": "2.0",
    "CityObjects": cityobjects,
    "vertices": vertices_local,
    "metadata": {
        "referenceSystem": "https://www.opengis.net/def/crs/EPSG/0/2056",
        "geographicalExtent": geographical_extent,
        "presentLoDs": ["2.0"],
        "datasetTitle": "Planned buildings"
    },
    "transform": {
        "scale":[1.0,1.0,1.0],
        "translate":[2592980.685,1119281.703,483.8]
    }
}

os.makedirs(os.path.dirname(output_cityjson), exist_ok=True)
with open(output_cityjson, "w", encoding="utf-8") as f:
    json.dump(cityjson, f, indent=2)

print(f"CityJSON exporté avec {len(cityobjects)} objets et {len(vertices_local)} sommets (nettoyés et transformés localement).")
