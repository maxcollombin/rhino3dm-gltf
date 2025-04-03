import rhino3dm

# Charger le modèle 3DM
input_file = "R21_3d.3dm"   # Remplace par ton fichier source
output_file = "R21_3d_trsl.3dm"

# Définir l'origine en EPSG:2056 (ex: Berne ou autre point de référence)
translation_vector = (2592980.6849999987, 1119281.7030000016, 483.8)

# Créer un objet Vector3d pour la translation
translation = rhino3dm.Vector3d(*translation_vector)

# Charger le fichier Rhino
model = rhino3dm.File3dm.Read(input_file)

# Vérifier si le modèle a été chargé
if not model:
    raise ValueError(f"Erreur : Impossible de charger {input_file}")

# Appliquer la translation aux objets du modèle
for obj in model.Objects:
    geometry = obj.Geometry
    if isinstance(geometry, rhino3dm.Point3d):  # Pour les points
        geometry.X += translation.X
        geometry.Y += translation.Y
        geometry.Z += translation.Z
    elif hasattr(geometry, "Translate"):  # Pour courbes, surfaces, etc.
        geometry.Translate(translation)

# Sauvegarder le modèle transformé
success = model.Write(output_file, 6)  # Version 6 de Rhino ou adapter

if success:
    print(f"✅ Modèle converti et enregistré sous {output_file}")
else:
    print("❌ Erreur lors de l'enregistrement")