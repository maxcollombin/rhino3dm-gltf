import rhino3dm
import os

def export_to_obj_by_layer(file_3dm, output_dir):
    model = rhino3dm.File3dm.Read(file_3dm)
    print(f"Objets dans le fichier: {len(model.Objects)}")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    layers = model.Layers
    for layer in layers:
        layer_name = layer.Name
        layer_index = layer.Index
        layer_objects = [obj for obj in model.Objects if obj.Attributes.LayerIndex == layer_index]

        if not layer_objects:
            continue

        file_obj = os.path.join(output_dir, f"{layer_name}.obj")
        with open(file_obj, "w") as obj_file:
            obj_file.write("# Exported from Rhino3dm\n")

            vertex_offset = 1  # Indice des sommets pour l'OBJ
            has_valid_geometry = False  # Flag to check if there is any valid geometry

            for obj in layer_objects:
                print(f"Processing object with ID: {obj.Attributes.Id}")

                geom = obj.Geometry
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

            if not has_valid_geometry:
                print(f"Aucune géométrie valide trouvée pour la couche {layer_name}, suppression du fichier vide.")
                os.remove(file_obj)
            else:
                print(f"Export terminé pour la couche {layer_name} : {file_obj}")

# Utilisation :
export_to_obj_by_layer("input/R21_3d.3dm", "output/obj")