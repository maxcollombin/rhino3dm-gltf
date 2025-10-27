import os
import json
import uuid
import rhino3dm
import numpy as np
from sklearn.decomposition import PCA

# ------------------------------
# 1. Fichiers
# ------------------------------
input_file = "./R21_3d.3dm"
output_cityjson = "./output/NouveauxBatiments_FME.city.json"

# ------------------------------
# 2. Charger le modèle
# ------------------------------
model = rhino3dm.File3dm.Read(input_file)
if not model:
    raise ValueError(f"Erreur : Impossible de charger {input_file}")
print(f" Modèle chargé : {input_file}")

# ------------------------------
# 3. Fonctions utilitaires
# ------------------------------
def remove_consecutive_duplicates(ring):
    if len(ring) < 3: return []
    cleaned = []
    for i, vertex in enumerate(ring):
        next_vertex = ring[(i + 1) % len(ring)]
        if vertex != next_vertex:
            cleaned.append(vertex)
    if len(cleaned) >= 3 and cleaned[0] == cleaned[-1]:
        cleaned = cleaned[:-1]
    return cleaned if len(cleaned) >= 3 else []

def validate_face(ring):
    if len(ring) < 3: return False
    for i in range(len(ring)):
        if ring[i] == ring[(i + 1) % len(ring)]:
            return False
    return True

def project_to_plane(vertices_coords):
    if len(vertices_coords) < 4: return vertices_coords
    points = np.array(vertices_coords)
    centroid = np.mean(points, axis=0)
    centered = points - centroid
    pca = PCA(n_components=3)
    pca.fit(centered)
    normal = pca.components_[2]
    projected_centered = centered - np.outer(np.dot(centered, normal), normal)
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
    if len(points) < 3: return "UndefinedSurface"
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
        return "GroundSurface" if z_mean < centroid[2] else "RoofSurface"
    elif vertical:
        return "WallSurface"
    else:
        return "ClosureSurface"

# ------------------------------
# 4. Collecte des bâtiments
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

buildings = {}
for obj in model.Objects:
    layer = model.Layers[obj.Attributes.LayerIndex]
    full_path = layer.FullPath
    if not full_path.startswith("Bâtiments Nouveaux"): continue
    if "EXISTANT" in full_path or "DEMOLITION" in full_path: continue
    parts = full_path.split("::")
    if len(parts) < 5: continue
    building_name = parts[-1]
    mesh = convert_to_mesh(obj.Geometry)
    if mesh:
        buildings.setdefault(building_name, []).append(mesh)

print(f" {len(buildings)} bâtiments détectés.")

# ------------------------------
# 5. Construction CityJSON FME-ready
# ------------------------------
cityobjects = {}

for building_name, meshes in buildings.items():
    building_uuid = str(uuid.uuid4())
    all_shells = []
    semantics_values = []
    semantics_faces = []

    for mesh in meshes:
        mesh.Normals.ComputeNormals()
        for face in mesh.Faces:
            face_vertices = list(face) if len(face)==4 else list(face[:3])
            ring = [add_vertex(mesh.Vertices[j]) for j in face_vertices]
            cleaned = remove_consecutive_duplicates(ring)
            if not validate_face(cleaned): continue
            coords = [vertices_global[idx] for idx in cleaned]
            projected = project_to_plane(coords)

            # Classification
            v1,v2,v3 = np.array(projected[:3])
            normal = np.cross(v2-v1, v3-v1)
            normal /= (np.linalg.norm(normal)+1e-9)
            nz = normal[2]
            face_type = "RoofSurface" if nz>0.9 else ("GroundSurface" if nz<-0.9 else "WallSurface")

            semantics_values.append(face_type)
            semantics_faces.append(len(all_shells))

            all_shells.append([[add_vertex(rhino3dm.Point3d(*c)) for c in projected]])

    if all_shells:
        cityobjects[building_uuid] = {
            "type": "Building",
            "attributes": {"buildingID": building_name},
            "geometry": [{
                "type": "MultiSurface",
                "lod": 2.0,
                "boundaries": all_shells,
                "semantics": {
                    "surfaces":[{"type":s} for s in semantics_values],
                    "values": semantics_faces
                }
            }]
        }

# ------------------------------
# 6. Écriture CityJSON
# ------------------------------
scale = [1.0,1.0,1.0]
translate = [2592980.685,1119281.703,0.0]

xs, ys, zs = zip(*vertices_global)
xs_t = [x*scale[0]+translate[0] for x in xs]
ys_t = [y*scale[1]+translate[1] for y in ys]
zs_t = [z*scale[2]+translate[2] for z in zs]
geographical_extent = [min(xs_t), min(ys_t), min(zs_t), max(xs_t), max(ys_t), max(zs_t)]

cityjson = {
    "type": "CityJSON",
    "version": "2.0",
    "CityObjects": cityobjects,
    "vertices": vertices_global,
    "metadata": {
        "referenceSystem": "EPSG:2056",
        "geographicalExtent": geographical_extent,
        "presentLoDs": [2.0],
        "datasetTitle": "Nouveaux bâtiments exportés depuis Rhino3dm pour FME"
    },
    "transform": {"scale": scale, "translate": translate}
}

os.makedirs(os.path.dirname(output_cityjson), exist_ok=True)
with open(output_cityjson,"w",encoding="utf-8") as f:
    json.dump(cityjson,f,indent=2)

print(f" CityJSON FME-ready exporté avec {len(cityobjects)} bâtiments et {len(vertices_global)} sommets.")
