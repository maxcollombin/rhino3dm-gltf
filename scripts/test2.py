import rhino3dm
import os

def log_error(log_file, message):
    """Écrit un message d'erreur dans le fichier de log."""
    print(f"Logging error: {message}")  # Débogage
    with open(log_file, "a") as log:
        log.write(message + "\n")

def export_to_obj_by_layer(file_3dm, output_dir):
    model = rhino3dm.File3dm.Read(file_3dm)
    print(f"Objets dans le fichier: {len(model.Objects)}")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    log_file = os.path.join(output_dir, "errors.log")  # Fichier de log pour les erreurs

    # Initialiser le fichier de log
    with open(log_file, "w") as log:
        log.write("Log d'erreurs pour l'exportation des fichiers OBJ\n")

    print(f"Log file path: {log_file}")  # Débogage

    layers = model.Layers
    for layer in layers:
        layer_name = layer.Name
        layer_index = layer.Index
        layer_objects = [obj for obj in model.Objects if obj.Attributes.LayerIndex == layer_index]

        if not layer_objects:
            error_message = f"Aucune géométrie trouvée pour la couche {layer_name}."
            print(error_message)
            log_error(log_file, error_message)
            continue

        file_obj = os.path.join(output_dir, f"{layer_name}.obj")
        with open(file_obj, "w") as obj_file:
            obj_file.write("# Exported from Rhino3dm\n")

            vertex_offset = 1  # Indice des sommets pour l'OBJ
            has_valid_geometry = False  # Flag to check if there is any valid geometry

            for obj in layer_objects:
                print(f"Processing object with ID: {obj.Attributes.Id}")

                geom = obj.Geometry
                try:
                    if isinstance(geom, rhino3dm.Mesh):
                        print(f"Object {obj.Attributes.Id} is a Mesh with {len(geom.Vertices)} vertices and {geom.Faces.Count} faces.")
                        has_valid_geometry = True
                        # Écriture des sommets
                        for v in geom.Vertices:
                            obj_file.write(f"v {v.X} {v.Y} {v.Z}\n")

                        # Écriture des faces
                        for i in range(geom.Faces.Count):
                            face = geom.Faces[i]  # Accéder à la face directement
                            A, B, C, D = face

                            if len(face) == 4:  # Quad
                                obj_file.write(f"f {A+vertex_offset} {B+vertex_offset} {C+vertex_offset} {D+vertex_offset}\n")
                            else:  # Triangle
                                obj_file.write(f"f {A+vertex_offset} {B+vertex_offset} {C+vertex_offset}\n")

                        vertex_offset += len(geom.Vertices)
                    elif isinstance(geom, rhino3dm.Brep):
                        print(f"Object {obj.Attributes.Id} is a Brep.")
                        # Convert Brep to Mesh
                        for face in geom.Faces:
                            brep_mesh = face.GetMesh(rhino3dm.MeshType.Render)
                            if brep_mesh:
                                print(f"Converted Brep face to Mesh with {len(brep_mesh.Vertices)} vertices and {brep_mesh.Faces.Count} faces.")
                                has_valid_geometry = True
                                # Écriture des sommets
                                for v in brep_mesh.Vertices:
                                    obj_file.write(f"v {v.X} {v.Y} {v.Z}\n")

                                # Écriture des faces
                                for j in range(brep_mesh.Faces.Count):
                                    face = brep_mesh.Faces[j]  # Accéder à la face directement
                                    A, B, C, D = face

                                    if len(face) == 4:  # Quad
                                        obj_file.write(f"f {A+vertex_offset} {B+vertex_offset} {C+vertex_offset} {D+vertex_offset}\n")
                                    else:  # Triangle
                                        obj_file.write(f"f {A+vertex_offset} {B+vertex_offset} {C+vertex_offset}\n")

                                vertex_offset += len(brep_mesh.Vertices)
                            else:
                                error_message = f"Impossible de convertir le Brep en Mesh pour l'objet {obj.Attributes.Id}."
                                print(error_message)
                                log_error(log_file, error_message)
                    elif isinstance(geom, rhino3dm.Extrusion):
                        print(f"Object {obj.Attributes.Id} is an Extrusion.")
                        # Tentative de conversion de l'Extrusion en Brep
                        brep = geom.ToBrep(splitKinkyFaces=True)
                        if brep:
                            print(f"Extrusion convertie en Brep, traitement comme Brep.")
                            for edge in brep.Edges:
                                curve = edge.ToNurbsCurve()
                                if curve:
                                    print(f"Exporting curve from extrusion.")
                                    has_valid_geometry = True
                                    # Écriture des points de la courbe comme sommets
                                    for point in curve.Points:
                                        # Conversion des coordonnées homogènes en coordonnées cartésiennes
                                        x = point.X / point.W if point.W != 0 else point.X
                                        y = point.Y / point.W if point.W != 0 else point.Y
                                        z = point.Z / point.W if point.W != 0 else point.Z
                                        obj_file.write(f"v {x} {y} {z}\n")
                        else:
                            error_message = f"Échec de la conversion de l'Extrusion en Brep pour l'objet {obj.Attributes.Id}."
                            print(error_message)
                            log_error(log_file, error_message)
                    else:
                        error_message = f"Type de géométrie non pris en charge : {type(geom)} pour l'objet {obj.Attributes.Id}."
                        print(error_message)
                        log_error(log_file, error_message)
                except Exception as e:
                    error_message = f"Erreur lors du traitement de l'objet {obj.Attributes.Id} : {e}"
                    print(error_message)
                    log_error(log_file, error_message)

            if not has_valid_geometry:
                error_message = f"Aucune géométrie valide trouvée pour la couche {layer_name}, suppression du fichier vide."
                print(error_message)
                log_error(log_file, error_message)
                os.remove(file_obj)
            else:
                print(f"Export terminé pour la couche {layer_name} : {file_obj}")

# Utilisation :
export_to_obj_by_layer("input/R21_3d.3dm", "output/obj")