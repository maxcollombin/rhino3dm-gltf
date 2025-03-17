import rhino3dm

def inspect_3dm(file_3dm):
    model = rhino3dm.File3dm.Read(file_3dm)
    print(f"📂 Fichier chargé : {file_3dm}\n")

    print(f"📌 Nombre total d'objets : {len(model.Objects)}")
    print(f"📦 Nombre total de groupes : {len(model.Groups)}")
    print(f"🗂️ Nombre total de couches (calques) : {len(model.Layers)}\n")

    # 🔹 Liste des groupes avec le nombre d'objets associés
    group_dict = {i: {"name": group.Name, "count": 0} for i, group in enumerate(model.Groups)}
    for obj in model.Objects:
        attr = obj.Attributes
        if attr.GroupCount > 0:
            group_ids = attr.GetGroupList()
            for g in group_ids:
                if g in group_dict:
                    group_dict[g]["count"] += 1

    print("📦 Groupes d'objets :")
    if group_dict:
        for g_id, data in group_dict.items():
            print(f"   - {data['name']} (ID {g_id}) : {data['count']} objets")
    else:
        print("   ❌ Aucun groupe trouvé.")

    # 🔹 Liste des couches (calques) avec leurs noms et indices
    print("\n🗂️ Couches disponibles :")
    for layer in model.Layers:
        print(f"   - {layer.Index} : {layer.Name}")

    # 🔹 Informations détaillées sur chaque objet
    for i, obj in enumerate(model.Objects):
        geom = obj.Geometry
        attr = obj.Attributes

        obj_name = attr.Name if attr.Name else f"Objet_{i}"
        obj_type = type(geom).__name__

        # Groupes associés
        group_info = "Aucun"
        if attr.GroupCount > 0:
            group_ids = attr.GetGroupList()
            group_names = [group_dict[g]["name"] for g in group_ids if g in group_dict]
            group_info = ", ".join(group_names) if group_names else "Groupes sans nom"

        # Couches associées
        layer_name = model.Layers[attr.LayerIndex].Name if 0 <= attr.LayerIndex < len(model.Layers) else "Inconnu"

        print(f"\n🔹 Objet {i + 1}:")
        print(f"   ├── Nom : {obj_name}")
        print(f"   ├── Type : {obj_type}")
        print(f"   ├── Groupe(s) : {group_info}")
        print(f"   ├── Calque : {layer_name}")

    print("\n✅ Inspection terminée.")

# 🔹 Utilisation :
inspect_3dm("input/R21_3d.3dm")
