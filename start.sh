#!/bin/bash
set -eu

VENV_NAME="venv"

# Vérifier la présence de GDAL Python
if ! python3 -c "import osgeo" >/dev/null 2>&1; then
    sudo apt-get update -qq
    sudo apt-get install -y -qq python3-gdal gdal-bin libgdal-dev
fi

# Créer l'environnement virtuel si nécessaire
if [ ! -d "$VENV_NAME" ]; then
    python3 -m venv "$VENV_NAME"
fi

# Activer l'environnement virtuel
# (doit être sourcé pour persister dans le shell actuel)
. "$VENV_NAME/bin/activate"

# Installer les dépendances Python
pip install --upgrade pip -q
pip install -q numpy rhino3dm pygltflib scikit-learn
pip install -q cjio
pip install -q rasterio
pip install -q geopandas
pip install -q shapely
pip install -q scipy
pip install -q fiona
pip install -q pyproj
pip install -q cityjson2jsonfg

# Créer un lien symbolique vers GDAL système si absent
SITE_PACKAGES=$(python3 -c "import site; print(site.getsitepackages()[0])")
if [ ! -e "$SITE_PACKAGES/osgeo" ]; then
    ln -sf /usr/lib/python3/dist-packages/osgeo "$SITE_PACKAGES/"
fi
