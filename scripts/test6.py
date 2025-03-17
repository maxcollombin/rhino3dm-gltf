import rhino3dm

def inspect_3dm(file_3dm):
    model = rhino3dm.File3dm.Read(file_3dm)
    print(f"📂 Fichier chargé : {file_3dm}\n")

    print(f"📌 Nombre total d'objets : {len(model.Objects)}")
    print(f"📦 Nombre total de groupes d'objets: {len(model.Groups)}")
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

    # 🔹 Regroupement des objets par groupe de calque (groupe -> calque)
    group_layer_dict = {}
    for obj in model.Objects:
        attr = obj.Attributes
        group_ids = attr.GetGroupList()
        layer_name = model.Layers[attr.LayerIndex].Name if 0 <= attr.LayerIndex < len(model.Layers) else "Inconnu"

        for g in group_ids:
            group_name = group_dict[g]["name"] if g in group_dict else "Inconnu"
            if group_name not in group_layer_dict:
                group_layer_dict[group_name] = {}
            if layer_name not in group_layer_dict[group_name]:
                group_layer_dict[group_name][layer_name] = []
            group_layer_dict[group_name][layer_name].append(obj)

    # 🔹 Affichage du regroupement par groupe de calque
    print("\n📦 Regroupement des objets par Groupe de Calque :")
    for group_name, layer_dict in group_layer_dict.items():
        print(f"  Groupe : {group_name}")
        for layer_name, objs in layer_dict.items():
            print(f"   ├── Calque : {layer_name} ({len(objs)} objets)")
            for obj in objs:
                print(f"     ├── {obj.Attributes.Name if obj.Attributes.Name else 'Objet sans nom'}")

    # 🔹 Regroupement des objets par groupe de calques
    layer_group_dict = {}
    for obj in model.Objects:
        attr = obj.Attributes
        layer_name = model.Layers[attr.LayerIndex].Name if 0 <= attr.LayerIndex < len(model.Layers) else "Inconnu"
        group_ids = attr.GetGroupList()

        if layer_name not in layer_group_dict:
            layer_group_dict[layer_name] = {}

        for g in group_ids:
            group_name = group_dict[g]["name"] if g in group_dict else "Inconnu"
            if group_name not in layer_group_dict[layer_name]:
                layer_group_dict[layer_name][group_name] = []
            layer_group_dict[layer_name][group_name].append(obj)

    # 🔹 Affichage du regroupement par groupe de calques
    print("\n🗂️ Regroupement des objets par groupe de calques:")
    for layer_name, group_dict in layer_group_dict.items():
        print(f"  Calque : {layer_name}")
        for group_name, objs in group_dict.items():
            print(f"   ├── Groupe : {group_name} ({len(objs)} objets)")
            for obj in objs:
                print(f"     ├── {obj.Attributes.Name if obj.Attributes.Name else 'Objet sans nom'}")

    print("\n✅ Inspection terminée.")

# 🔹 Utilisation :
inspect_3dm("input/R21_3d.3dm")
