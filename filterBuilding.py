from cjio import cityjson
import json

input_file = "NouveauxBatiments.city.json"
output_file = "batiment_unique.city.json"
building_id = "3475bdd2-f071-410e-96c5-a9af168f1f44"

with open(input_file, "r", encoding="utf-8") as f:
    cj = cityjson.reader(file=f)

# Work directly with raw CityObjects
cityobjects = cj.j["CityObjects"]

if building_id not in cityobjects:
    print(f"Error: Building ID {building_id} not found in raw CityObjects.")
    exit(1)

# Extraire le bâtiment et ses enfants
ids_to_keep = set([building_id])
ids_to_keep.update(cityobjects[building_id].get("children", []))

# Filtrer les CityObjects
filtered_cityobjects = {k: v for k, v in cityobjects.items() if k in ids_to_keep}

# Créer un nouveau CityJSON minimal
new_cityjson = cj.j.copy()
new_cityjson["CityObjects"] = filtered_cityobjects

# Optionally, remove unused vertices (not handled here, but can be added)

with open(output_file, "w", encoding="utf-8") as f:
    json.dump(new_cityjson, f, ensure_ascii=False, indent=2)

print(f"Bâtiment {building_id} exporté dans {output_file}")