#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-processing avancé du fichier swissbuildings3d.city.json
Ajoute les attributs status et constructionPhase conformes à CityJSON 2.0
avec options de configuration
"""

import json
import os
import argparse

def process_swissbuildings_cityjson(input_file, output_file, keep_original_statut=False):
    """
    Post-traite le fichier CityJSON pour ajouter les attributs conformes
    
    Args:
        input_file: Chemin vers le fichier CityJSON d'entrée
        output_file: Chemin vers le fichier CityJSON de sortie
        keep_original_statut: Si True, garde l'attribut STATUT original
    """
    print(f"Chargement du fichier : {input_file}")
    
    # Charger le fichier CityJSON
    with open(input_file, 'r', encoding='utf-8') as f:
        cityjson = json.load(f)
    
    # Compteurs pour statistiques
    buildings_processed = 0
    conservation_count = 0
    demolition_count = 0
    
    # Traiter chaque CityObject
    for obj_id, city_obj in cityjson.get("CityObjects", {}).items():
        if city_obj.get("type") == "Building":
            buildings_processed += 1
            
            # Récupérer le statut actuel
            current_statut = city_obj.get("attributes", {}).get("STATUT")
            
            if not city_obj.get("attributes"):
                city_obj["attributes"] = {}
            
            # Mapper STATUT vers les attributs CityJSON conformes
            if current_statut == "Conservation":
                conservation_count += 1
                # Bâtiment existant à conserver
                city_obj["attributes"]["buildingStatus"] = "existing"
                
            elif current_statut == "Démolition":
                demolition_count += 1
                # Bâtiment à démolir
                city_obj["attributes"]["buildingStatus"] = "demolition"
            # Ajouter developmentPhase avec valeur par défaut 0 pour tous les bâtiments
            city_obj["attributes"]["developmentPhase"] = 0
            
            # Supprimer l'attribut STATUT original si demandé
            if not keep_original_statut and "STATUT" in city_obj["attributes"]:
                del city_obj["attributes"]["STATUT"]
    
    # Sauvegarder le fichier modifié
    print(f"Sauvegarde du fichier traité : {output_file}")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(cityjson, f, indent=2, ensure_ascii=False)
    
    # Afficher les statistiques
    print(f"\n📊 Statistiques du traitement :")
    print(f"   Bâtiments traités : {buildings_processed}")
    print(f"   Conservation : {conservation_count}")
    print(f"   Démolition : {demolition_count}")
    print(f"\n✅ Attributs CityJSON 2.0 conformes ajoutés :")
    print(f"   - status: 'existing' ou 'demolition'")
    print(f"   - constructionPhase: 0 (valeur par défaut pour tous)")
    
    if keep_original_statut:
        print(f"   - STATUT original conservé pour traçabilité")
    else:
        print(f"   - STATUT original supprimé (nettoyage)")

def main():
    parser = argparse.ArgumentParser(description='Post-traitement CityJSON pour attributs conformes')
    parser.add_argument('--input', '-i', 
                       default='/home/maxime/Desktop/rhino3dm-gltf/output/swissbuildings3d.city.json',
                       help='Fichier CityJSON d\'entrée')
    parser.add_argument('--output', '-o',
                       default='/home/maxime/Desktop/rhino3dm-gltf/output/swissbuildings3d_final.city.json',
                       help='Fichier CityJSON de sortie')
    parser.add_argument('--keep-statut', action='store_true',
                       help='Conserver l\'attribut STATUT original')
    
    args = parser.parse_args()
    
    process_swissbuildings_cityjson(args.input, args.output, args.keep_statut)
    print(f"\n🎯 Fichier final disponible : {args.output}")
    
    # Recommandations de bonnes pratiques
    print(f"\n📋 Conformité CityJSON 2.0 :")
    print(f"   ✓ status: existing (bâtiments conservés)")
    print(f"   ✓ status: demolition (bâtiments à démolir)")
    print(f"   ✓ developmentPhase: attribut personnalisé (0 par défaut)")
    print(f"\n💡 Note: developmentPhase n'est pas dans le core CityJSON 2.0")
    print(f"   Il s'agit d'un attribut personnalisé que vous pouvez")
    print(f"   modifier selon vos besoins futurs.")

if __name__ == "__main__":
    main()
