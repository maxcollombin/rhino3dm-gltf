import rhino3dm

def export_to_obj(file_3dm, file_obj):
    model = rhino3dm.File3dm.Read(file_3dm)
    print(f"Objets dans le fichier: {len(model.Objects)}")

    with open(file_obj, "w") as obj_file:
        obj_file.write("# Exported from Rhino3dm\n")

        vertex_offset = 1  # Indice des sommets pour l'OBJ
        for obj in model.Objects:
            geom = obj.Geometry
            if isinstance(geom, rhino3dm.Mesh):
                # Écriture des sommets
                for v in geom.Vertices:
                    obj_file.write(f"v {v.X} {v.Y} {v.Z}\n")

                # Écriture des faces
                for i in range(geom.Faces.Count):
                    face = geom.Faces[i]  # Accéder à la face directement

                    if len(face) == 4:  # Quad
                        obj_file.write(f"f {face[0]+vertex_offset} {face[1]+vertex_offset} {face[2]+vertex_offset} {face[3]+vertex_offset}\n")
                    else:  # Triangle
                        obj_file.write(f"f {face[0]+vertex_offset} {face[1]+vertex_offset} {face[2]+vertex_offset}\n")

                vertex_offset += len(geom.Vertices)
            elif isinstance(geom, rhino3dm.Brep):
                # Convert Brep to Mesh
                brep_mesh = geom.Faces[0].GetMesh(rhino3dm.MeshType.Render)
                if brep_mesh:
                    # Écriture des sommets
                    for v in brep_mesh.Vertices:
                        obj_file.write(f"v {v.X} {v.Y} {v.Z}\n")

                    # Écriture des faces
                    for i in range(brep_mesh.Faces.Count):
                        face = brep_mesh.Faces[i]  # Accéder à la face directement

                        if len(face) == 4:  # Quad
                            obj_file.write(f"f {face[0]+vertex_offset} {face[1]+vertex_offset} {face[2]+vertex_offset} {face[3]+vertex_offset}\n")
                        else:  # Triangle
                            obj_file.write(f"f {face[0]+vertex_offset} {face[1]+vertex_offset} {face[2]+vertex_offset}\n")

                    vertex_offset += len(brep_mesh.Vertices)

    print(f"Export terminé : {file_obj}")

# 🔹 Utilisation :
export_to_obj("input/R21_3d.3dm", "output/test1.obj")
