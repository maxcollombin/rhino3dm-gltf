import rhino3dm

def inspect_3dm(file_3dm):
    model = rhino3dm.File3dm.Read(file_3dm)
    print(f"📂 Fichier chargé : {file_3dm}")
    print(f"📌 Nombre total d'objets : {len(model.Objects)}")
    print(f"📦 Nombre de groupes : {len(model.Groups)}")

    # Récupérer les groupes sous forme de dictionnaire {index: nom}
    group_dict = {i: group.Name for i, group in enumerate(model.Groups)}

    for i, obj in enumerate(model.Objects):
        geom = obj.Geometry
        attr = obj.Attributes

        # Récupérer le nom de l'objet
        obj_name = attr.Name if attr.Name else f"Objet_{i}"

        # Déterminer le type de l'objet
        obj_type = type(geom).__name__

        # Vérifier s'il appartient à un ou plusieurs groupes
        group_info = "Aucun"
        if attr.GroupCount > 0:
            group_ids = attr.GetGroupList()  # Liste des indices des groupes
            group_names = [group_dict[g] for g in group_ids if g in group_dict]
            group_info = ", ".join(group_names) if group_names else "Groupes sans nom"

        print(f"\n🔹 Objet {i + 1}:")
        print(f"   ├── Nom : {obj_name}")
        print(f"   ├── Type : {obj_type}")
        print(f"   ├── Groupes : {group_info}")

        # Afficher les informations de géométrie si c'est un mesh
        if isinstance(geom, rhino3dm.Mesh):
            print(f"   ├── Nombre de sommets : {len(geom.Vertices)}")
            print(f"   ├── Nombre de faces : {geom.Faces.Count}")

    print("\n✅ Inspection terminée.")

# 🔹 Utilisation :
inspect_3dm("input/R21_3d.3dm")
