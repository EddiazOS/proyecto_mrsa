#!/usr/bin/env python3
"""
preparar_ligandos.py
====================
Prepara todos los ligandos del proyecto (compuestos HPLC y controles positivos)
en formato PDBQT usando RDKit y Meeko.

Pasos:
  1. Identifica todos los archivos SDF en data/ligands/
  2. Filtra múltiples fragmentos (se queda con el fragmento principal, p.ej. eliminando sales/iones/aguas)
  3. Asegura que existan hidrógenos explícitos en 3D
  4. Ejecuta mk_prepare_ligand.py de Meeko para generar los archivos .pdbqt

Uso:
  conda run -n docking_env python scripts/preparar_ligandos.py
"""

import os
import glob
import subprocess
from rdkit import Chem

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, os.pardir))
    ligands_dir = os.path.join(project_root, "data", "ligands")
    
    sdf_files = glob.glob(os.path.join(ligands_dir, "*.sdf"))
    
    print("====================================================================")
    print("  PREPARACIÓN AUTOMÁTICA DE LIGANDOS (Meeko + RDKit)")
    print("====================================================================")
    
    success_count = 0
    fail_count = 0
    
    for sdf_path in sdf_files:
        basename = os.path.basename(sdf_path)
        name = os.path.splitext(basename)[0]
        
        # Omitir archivos temporales si existen
        if name.startswith("tmp_"):
            continue
            
        out_pdbqt = os.path.join(ligands_dir, f"{name}.pdbqt")
        print(f"\n[*] Procesando: {basename}")
        
        try:
            # Leer molécula con hidrógenos explícitos del SDF
            suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
            mols = [m for m in suppl if m is not None]
            if not mols:
                print(f"  [ERROR] No se pudo leer la molécula de {sdf_path}")
                fail_count += 1
                continue
                
            mol = mols[0]
            
            # Filtro de fragmentos (quedarse con el fragmento de mayor número de átomos pesados)
            frags = Chem.GetMolFrags(mol, asMols=True)
            if len(frags) > 1:
                print(f"  [WARN] Múltiples fragmentos ({len(frags)}) detectados. Conservando el más grande...")
                mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
            else:
                mol = frags[0]
            
            # Agregar hidrógenos explícitos si faltan (manteniendo coordenadas 3D)
            mol = Chem.AddHs(mol, addCoords=True)
            
            # Escribir molécula limpia en un SDF temporal
            tmp_sdf = os.path.join(ligands_dir, f"tmp_{name}.sdf")
            writer = Chem.SDWriter(tmp_sdf)
            writer.write(mol)
            writer.close()
            
            # Ejecutar mk_prepare_ligand.py
            cmd = [
                "mk_prepare_ligand.py",
                "-i", tmp_sdf,
                "-o", out_pdbqt
            ]
            
            # Ejecutar en el mismo entorno o con el ejecutable en el path
            res = subprocess.run(cmd, capture_output=True, text=True)
            
            if res.returncode == 0:
                print(f"  [OK] Exitosamente preparado -> {os.path.basename(out_pdbqt)}")
                success_count += 1
            else:
                print(f"  [ERROR] Error en mk_prepare_ligand.py:")
                print(res.stderr.strip())
                fail_count += 1
                
            # Eliminar temporal
            if os.path.exists(tmp_sdf):
                os.remove(tmp_sdf)
                
        except Exception as e:
            print(f"  [ERROR] Error inesperado procesando {basename}: {e}")
            fail_count += 1
            
    print("\n====================================================================")
    print(f"  Resumen: {success_count} exitosos, {fail_count} fallidos.")
    print("====================================================================")

if __name__ == "__main__":
    main()
