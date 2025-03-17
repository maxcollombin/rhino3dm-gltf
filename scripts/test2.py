import rhino3dm

def export_to_obj(file_3dm, file_obj):
    model = rhino3dm.File3dm.Read(file_3dm)
    print(f"Objets dans le fichier: {len(model.Objects)}")

    with open(file_obj, "w") as obj_file:
        obj_file.write("# Exported from Rhino3dm\n")

        vertex_offset = 1  # Indice des sommets pour l'OBJ
        for obj in model.Objects:
            print(f"Processing object with ID: {obj.Attributes.Id}")

            geom = obj.Geometry
            if isinstance(geom, rhino3dm.Mesh):
                print(f"Object {obj.Attributes.Id} is a Mesh with {len(geom.Vertices)} vertices and {geom.Faces.Count} faces.")
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

    print(f"Export terminé : {file_obj}")

# 🔹 Utilisation :
export_to_obj("input/R21_3d.3dm", "output/test2.obj")
