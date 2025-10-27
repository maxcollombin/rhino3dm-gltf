import json
import os
import math
import numpy as np

# -----------------------------
# Paramètres procéduraux
# -----------------------------
FEN_WIDTH = 1.0
FEN_HEIGHT = 1.5
FEN_SPACING_X = 0.5
FLOOR_HEIGHT = 3.0
DOOR_WIDTH = 1.0
DOOR_HEIGHT = 2.0
BALCON_DEPTH = 0.8
BALCON_WIDTH_NB_WINDOWS = 2
BEAM_WIDTH = 0.2

# -----------------------------
# Fonctions utilitaires
# -----------------------------
def add_vertex(vertices, vertex, vertex_map):
    key = tuple(round(c,5) for c in vertex)
    if key in vertex_map:
        return vertex_map[key]
    idx = len(vertices)
    vertices.append(list(key))
    vertex_map[key] = idx
    return idx

def create_face(indices):
    return indices  # pas besoin de répéter le premier sommet pour CityJSON

def create_multisurface(coords_list, vertices, vertex_map):
    indices = [add_vertex(vertices, tuple(v), vertex_map) for v in coords_list]
    return {"type": "MultiSurface", "boundaries": [[create_face(indices)]]}

def create_solid(faces_coords_list, vertices, vertex_map):
    faces=[]
    for face in faces_coords_list:
        idx = [add_vertex(vertices, tuple(v), vertex_map) for v in face]
        faces.append(create_face(idx))
    return {"type": "Solid", "boundaries": [[faces]]}

def create_box_solid(corners, vertices, vertex_map):
    # corners: list of 8 np.array points (box corners)
    idx = [add_vertex(vertices, v, vertex_map) for v in corners]
    faces = [
        [idx[0], idx[1], idx[2], idx[3]], # bottom
        [idx[4], idx[5], idx[6], idx[7]], # top
        [idx[0], idx[1], idx[5], idx[4]], # side 1
        [idx[1], idx[2], idx[6], idx[5]], # side 2
        [idx[2], idx[3], idx[7], idx[6]], # side 3
        [idx[3], idx[0], idx[4], idx[7]]  # side 4
    ]
    return {"type": "Solid", "boundaries": [[faces]]}

def calculate_bbox(surface):
    xs = [p[0] for p in surface]
    ys = [p[1] for p in surface]
    return min(xs), min(ys), max(xs), max(ys)

def get_wall_frame(surface_coords):
    p0 = np.array(surface_coords[0])
    p1 = np.array(surface_coords[1])
    p2 = np.array(surface_coords[2])
    x_axis = p1 - p0
    x_axis /= np.linalg.norm(x_axis)
    v = p2 - p0
    y_axis = v - np.dot(v, x_axis) * x_axis
    y_axis /= np.linalg.norm(y_axis)
    z_axis = np.cross(x_axis, y_axis)
    z_axis /= np.linalg.norm(z_axis)
    return p0, x_axis, y_axis, z_axis

# -----------------------------
# Charger CityJSON LOD2
# -----------------------------
with open("batiment_unique.city.json", "r") as f:
    cj = json.load(f)

cj.setdefault("vertices", [])
cj.setdefault("transform", {"scale": [1,1,1], "translate": [0,0,0]})
vertex_map = {}
proc_obj_count = 0

# -----------------------------
# Traitement 1er bâtiment
# -----------------------------
for bldg_id, bldg in cj["CityObjects"].items():
    if bldg.get("type") != "Building":
        continue

    print(f"Traitement bâtiment {bldg_id}...")

    for child_id in bldg.get("children", []):
        child = cj["CityObjects"][child_id]
        part_id = f"{child_id}_part"
        cj["CityObjects"][part_id] = {"type":"BuildingPart","parents":[child_id],"geometry":[],"children":[]}
        child.setdefault("children",[]).append(part_id)

        fen_geoms=[]
        balcon_geoms=[]
        door_geom_added=False

        for geom in child.get("geometry",[]):
            if geom.get("type") != "MultiSurface":
                continue
            boundaries = geom.get("boundaries",[])
            semantics = geom.get("semantics",{})
            surfaces_types = semantics.get("surfaces",[])
            surfaces_values = semantics.get("values",[])

            for i, surface in enumerate(boundaries):
                surface_type_idx = surfaces_values[i] if i < len(surfaces_values) else None
                surface_type = (
                    surfaces_types[surface_type_idx]["type"]
                    if surface_type_idx is not None and surface_type_idx < len(surfaces_types)
                    else None
                )
                if surface_type != "WallSurface":
                    continue

                surface_indices = surface[0] if isinstance(surface[0], list) else surface
                if len(surface_indices) < 3:
                    continue

                surface_coords = [cj["vertices"][idx] for idx in surface_indices]
                width = np.linalg.norm(np.array(surface_coords[1]) - np.array(surface_coords[0]))
                height = np.linalg.norm(np.array(surface_coords[2]) - np.array(surface_coords[1]))
                p0, x_axis, y_axis, z_axis = get_wall_frame(surface_coords)

                num_windows_x = max(1, math.floor(width / (FEN_WIDTH + FEN_SPACING_X)))
                step = max(1, BALCON_WIDTH_NB_WINDOWS)

                # Fenêtres
                for i_win in range(num_windows_x):
                    x0 = i_win * (FEN_WIDTH + FEN_SPACING_X)
                    x1 = x0 + FEN_WIDTH
                    y0 = 0
                    y1 = FEN_HEIGHT
                    rect = create_multisurface(
                        [p0 + x_axis*x0 + y_axis*y0,
                         p0 + x_axis*x1 + y_axis*y0,
                         p0 + x_axis*x1 + y_axis*y1,
                         p0 + x_axis*x0 + y_axis*y1],
                        cj["vertices"], vertex_map)
                    fen_geoms.append(rect)

                # Porte seulement premier étage
                if not door_geom_added:
                    door_x0 = width/2 - DOOR_WIDTH/2
                    door_x1 = door_x0 + DOOR_WIDTH
                    door_y0 = 0
                    door_y1 = DOOR_HEIGHT
                    door_rect = create_multisurface(
                        [p0 + x_axis*door_x0 + y_axis*door_y0,
                         p0 + x_axis*door_x1 + y_axis*door_y0,
                         p0 + x_axis*door_x1 + y_axis*door_y1,
                         p0 + x_axis*door_x0 + y_axis*door_y1],
                        cj["vertices"], vertex_map)
                    door_id = "door"
                    cj["CityObjects"][door_id] = {"type":"BuildingInstallation",
                                                  "geometry":[door_rect],
                                                  "parents":[part_id],
                                                  "attributes":{"function":"Door"}}
                    cj["CityObjects"][part_id]["children"].append(door_id)
                    proc_obj_count+=1
                    door_geom_added=True

                # Balcons + poutres correctement placés
                y_pos = FEN_HEIGHT
                # Calcul du centre vertical du mur
                mur_base = p0
                mur_top = p0 + y_axis*height
                for i_win in range(0, num_windows_x, step):
                    x0 = i_win * (FEN_WIDTH + FEN_SPACING_X)
                    x1 = x0 + step*FEN_WIDTH + (step-1)*FEN_SPACING_X
                    # Plancher balcon centré devant les fenêtres
                    bal_y0 = y_pos
                    bal_y1 = y_pos + BALCON_DEPTH
                    # Rectangle sur la façade
                    base0 = p0 + x_axis*x0 + y_axis*bal_y0
                    base1 = p0 + x_axis*x1 + y_axis*bal_y0
                    base2 = p0 + x_axis*x1 + y_axis*bal_y1
                    base3 = p0 + x_axis*x0 + y_axis*bal_y1
                    # Extrusion perpendiculaire
                    top0 = base0 + z_axis*BALCON_DEPTH
                    top1 = base1 + z_axis*BALCON_DEPTH
                    top2 = base2 + z_axis*BALCON_DEPTH
                    top3 = base3 + z_axis*BALCON_DEPTH
                    corners = [base0, base1, base2, base3, top0, top1, top2, top3]
                    box = create_box_solid(corners, cj["vertices"], vertex_map)
                    balcon_geoms.append(box)
                    # Poutres sur les bords verticaux du mur
                    # Poutre gauche
                    b0 = p0 + x_axis*x0
                    b1 = b0 + y_axis*0
                    b2 = b0 + y_axis*height
                    b3 = b0 + y_axis*height
                    t0 = b0 + z_axis*BALCON_DEPTH
                    t1 = b1 + z_axis*BALCON_DEPTH
                    t2 = b2 + z_axis*BALCON_DEPTH
                    t3 = b3 + z_axis*BALCON_DEPTH
                    corners1 = [b0, b1, b2, b3, t0, t1, t2, t3]
                    box1 = create_box_solid(corners1, cj["vertices"], vertex_map)
                    # Poutre droite
                    b0 = p0 + x_axis*(x1-BEAM_WIDTH)
                    b1 = b0 + y_axis*0
                    b2 = b0 + y_axis*height
                    b3 = b0 + y_axis*height
                    t0 = b0 + z_axis*BALCON_DEPTH
                    t1 = b1 + z_axis*BALCON_DEPTH
                    t2 = b2 + z_axis*BALCON_DEPTH
                    t3 = b3 + z_axis*BALCON_DEPTH
                    corners2 = [b0, b1, b2, b3, t0, t1, t2, t3]
                    box2 = create_box_solid(corners2, cj["vertices"], vertex_map)
                    balcon_geoms.extend([box1, box2])

        # Ajouter fenêtres
        if fen_geoms:
            cj["CityObjects"]["windows"] = {"type":"BuildingInstallation",
                                            "geometry": fen_geoms,
                                            "parents":[part_id],
                                            "attributes":{"function":"Windows"}}
            cj["CityObjects"][part_id]["children"].append("windows")
            proc_obj_count += 1

        # Ajouter balcons
        if balcon_geoms:
            cj["CityObjects"]["balcon_beams"] = {"type":"BuildingInstallation",
                                                 "geometry": balcon_geoms,
                                                 "parents":[part_id],
                                                 "attributes":{"function":"Balcony+Beams"}}
            cj["CityObjects"][part_id]["children"].append("balcon_beams")
            proc_obj_count += 1

    break  # seulement premier bâtiment

# -----------------------------
# Sauvegarde CityJSON
# -----------------------------
os.makedirs("output", exist_ok=True)
output_file = "output/test.city.json"
with open(output_file, "w") as f:
    json.dump(cj, f, indent=2)

print(f"\n✅ LOD3 test généré pour 1 bâtiment avec {proc_obj_count} installations.")
print(f"📄 Fichier : {output_file}")
