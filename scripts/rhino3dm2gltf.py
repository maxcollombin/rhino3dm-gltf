import rhino3dm
import os
import csv

# Chemin d'accès au fichier
file_path = "./input/R21_3d.3dm"
output_path = "./output/first_layer.gltf"

# Vérification de l'existence  du fichier
if not os.path.exists(file_path):
    print(f"Error: The file '{file_path}' does not exist.")
else:
    # Lecture du fichier .3dm
    model = rhino3dm.File3dm.Read(file_path)
    if model is None:
        print("Erreur lors du chargement du fichier.")
    else:
        print("Fichier correctement chargé.")

        # Affichage des informations générales du fichier
        print(f"Systèmes d'unité de mesures: {model.Settings.ModelUnitSystem}")
        print(f"Point de référence du modèle: {model.Settings.ModelBasePoint}")
        
        # Affichage des informations relatives aux couches
        layers = model.Layers
        print(f"Nombre de couches: {len(layers)}")
        
        # Définition du fichier en sortie
        csv_output_path = "./output/layers_info.csv"

        # Ecriture du fichier CSV avec les entêtes
        with open(csv_output_path, mode='w', newline='') as csv_file:
            fieldnames = ['Groupe de couches', 'Index de la couche', 'Nom de la couche', "Nombre d'objets"]
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()

            # Itération sur les couches
            for i, layer in enumerate(layers):
                layer_objects = [obj for obj in model.Objects if obj.Attributes.LayerIndex == layer.Index]
                writer.writerow({'Groupe de couches': layer.ParentLayerId, 'Index de la couche': layer.Index, 'Nom de la couche': layer.Name, "Nombre d'objets" : len(layer_objects)})

        print(f"Informations des couches écrites dans le fichier CSV : {csv_output_path}")
