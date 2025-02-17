import rhino3dm
import os

# Lecture du fichier .3dm
file_path = "./input/R21_3d.3dm"

# Vérifier si le fichier est accessible
if not os.path.exists(file_path):
    print(f"Erreur : le fichier '{file_path}' n'existe pas.")
else:
    model = rhino3dm.File3dm.Read(file_path)
    if model is None:
        print("Erreur : Impossible de charger le fichier .3dm.")
    else:
        print("Fichier chargé avec succès.")
        # Vérifier le système d'unités du modèle
        unit_system = model.Settings.ModelUnitSystem
        print(f"Système d'unités du modèle : {unit_system}")
        # Obtenir la BoundingBox du modèle
        bbox = model.Objects.GetBoundingBox()
        if bbox:
            print(f"BoundingBox: Min=({bbox.Min.X}, {bbox.Min.Y}, {bbox.Min.Z}), Max=({bbox.Max.X}, {bbox.Max.Y}, {bbox.Max.Z})")
        else:
            print("Erreur : Impossible d'obtenir la BoundingBox.")