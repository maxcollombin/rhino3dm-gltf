import rhino3dm

model = rhino3dm.File3dm.Read("input/R21_3d.3dm")
print("Objets dans le fichier:", len(model.Objects))

# Exporter un objet (exemple: premier objet en JSON)
obj = model.Objects[0]
print(obj.Geometry)
