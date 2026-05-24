#!/usr/bin/env python3
import os
import urllib.request
import urllib.error

def descargar_3d_pubchem(cid, nombre_salida):
    """
    Descarga la estructura tridimensional en formato SDF de un compuesto químico
    en PubChem empleando su identificador numérico CID (Compound ID).
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    carpeta_destino = os.path.join(base_dir, "data", "ligands")
    os.makedirs(carpeta_destino, exist_ok=True)
    
    ruta_completa = os.path.join(carpeta_destino, f"{nombre_salida}_3D.sdf")
    
    # Endpoint oficial de PubChem PUG-REST para conformaciones 3D
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    
    print(f"Intentando descargar {nombre_salida} (CID {cid}) desde PubChem...")
    
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read()
            
        with open(ruta_completa, "wb") as file:
            file.write(data)
            
        print(f"  [EXITO] Guardado en: {os.path.relpath(ruta_completa, base_dir)}")
        return ruta_completa
        
    except urllib.error.HTTPError as e:
        print(f"  [ERROR HTTP] PubChem devolvió código {e.code}. Verifique si el compuesto tiene confórmero 3D.")
    except Exception as e:
        print(f"  [ERROR CONEXION] No se pudo conectar a los servidores: {e}")
        
    return None

def descargar_receptor_pdb(pdb_id, nombre_salida):
    """
    Descarga una estructura macromolecular (molde biológico) en formato PDB
    desde la base de datos oficial RCSB Protein Data Bank.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    carpeta_destino = os.path.join(base_dir, "data", "receptors")
    os.makedirs(carpeta_destino, exist_ok=True)
    
    ruta_completa = os.path.join(carpeta_destino, f"{nombre_salida}.pdb")
    
    # Endpoint de descarga de archivos oficiales de RCSB
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    
    print(f"Intentando descargar molde PDB {pdb_id.upper()} desde RCSB...")
    
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read()
            
        with open(ruta_completa, "wb") as file:
            file.write(data)
            
        print(f"  [EXITO] Guardado en: {os.path.relpath(ruta_completa, base_dir)}")
        return ruta_completa
        
    except urllib.error.HTTPError as e:
        print(f"  [ERROR HTTP] RCSB devolvió código {e.code}. Verifique si el PDB {pdb_id} es correcto.")
    except Exception as e:
        print(f"  [ERROR CONEXION] No se pudo conectar a los servidores: {e}")
        
    return None

if __name__ == "__main__":
    # Descargas de ejemplo y preparación de controles
    print("=== Iniciando Descargas Programadas del Proyecto ===")
    
    # 1. Descargar la Quercetina (Control de MurG)
    descargar_3d_pubchem(5280343, "Quercetin")
    
    # 2. Descargar el Molde Cristalino de MurG (1F0K) para alineamientos futuros
    descargar_receptor_pdb("1F0K", "MurG_template_1F0K")
    
    print("\n=== Descargas Finalizadas ===")
