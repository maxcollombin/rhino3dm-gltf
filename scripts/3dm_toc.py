import rhino3dm
from collections import defaultdict

# Charger le modèle 3DM
# input_file = "./R21_3d_trsl.3dm"
input_file = "./R21_3d.3dm"
output_file = "./output/toc1.txt"
model = rhino3dm.File3dm.Read(input_file)

# Vérifier si le modèle a été chargé
if not model:
    raise ValueError(f"Erreur : Impossible de charger {input_file}")

# Construire une hiérarchie imbriquée des calques
def build_hierarchy(layers):
    hierarchy = defaultdict(list)
    for layer in layers:
        parts = layer.FullPath.split("::")
        for i in range(len(parts) - 1):
            parent = "::".join(parts[:i + 1])
            child = "::".join(parts[:i + 2])
            hierarchy[parent].append(child)
    return hierarchy

# Fonction récursive pour écrire la hiérarchie dans le fichier
def write_hierarchy(f, hierarchy, parent, level=0):
    children = sorted(set(hierarchy.get(parent, [])))
    for child in children:
        name = child.split("::")[-1]
        f.write(f"{'  ' * level}- {name}\n")
        write_hierarchy(f, hierarchy, child, level + 1)

# Collecter les calques et construire la hiérarchie
layers = model.Layers
hierarchy = build_hierarchy(layers)

# Écrire la hiérarchie dans un fichier
with open(output_file, "w", encoding="utf-8") as f:
    f.write("Hiérarchie des calques regroupés :\n")
    root_layers = sorted(set(layer.FullPath.split("::")[0] for layer in layers))
    for root in root_layers:
        f.write(f"- {root}\n")
        write_hierarchy(f, hierarchy, root)

print(f"✅ Hiérarchie des calques regroupés écrite dans le fichier : {output_file}")