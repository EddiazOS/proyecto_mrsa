#!/usr/bin/env python3
import os
import shutil

def main():
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print(f"Reorganizando espacio de trabajo en: {project_root}")
    
    # 1. Crear carpetas de destino
    folders = [
        "docs",
        "notebooks",
        "manuscript",
        "scripts/diagnostics",
        "docking/vina_output/GyrB",
        "docking/vina_output/MurG",
        "docking/vina_output/PBP2a"
    ]
    
    for folder in folders:
        path = os.path.join(project_root, folder)
        os.makedirs(path, exist_ok=True)
        print(f"  [OK] Carpeta creada o verificada: {folder}")
        
    # 2. Mover archivos de documentación del root a docs/
    docs_to_move = [
        "context.md",
        "reorganizacion_proyecto.md",
        "limpieza_receptores.md",
        "metodologia_gridbox_murg.md",
        "auditoria_gridbox_y_murg.md"
    ]
    
    print("\n--- Organizando Documentación ---")
    for doc in docs_to_move:
        src = os.path.join(project_root, doc)
        dst = os.path.join(project_root, "docs", doc)
        if os.path.exists(src):
            shutil.move(src, dst)
            print(f"  [MOVED] {doc} -> docs/")
            
    # 3. Mover notebooks a notebooks/
    notebooks_to_move = [
        ("scripts/dinamica_colab.ipynb", "notebooks/dinamica_colab.ipynb"),
        ("scripts/guia_docking.ipynb", "notebooks/guia_docking.ipynb")
    ]
    
    print("\n--- Organizando Notebooks ---")
    for src_rel, dst_rel in notebooks_to_move:
        src = os.path.join(project_root, src_rel)
        dst = os.path.join(project_root, dst_rel)
        if os.path.exists(src):
            shutil.move(src, dst)
            print(f"  [MOVED] {os.path.basename(src)} -> notebooks/")
            
    # 4. Mover manuscrito a manuscript/
    manuscript_to_move = [
        "manuscript_metodologia.typ",
        "manuscript_metodologia.pdf"
    ]
    
    print("\n--- Organizando Manuscrito ---")
    for ms in manuscript_to_move:
        src = os.path.join(project_root, ms)
        dst = os.path.join(project_root, "manuscript", ms)
        if os.path.exists(src):
            shutil.move(src, dst)
            print(f"  [MOVED] {ms} -> manuscript/")
            
    # 5. Mover scripts de diagnóstico a scripts/diagnostics/
    diagnostics_to_move = [
        "check_lit_residues.py",
        "find_close_residues.py",
        "analyze_afzelin_contacts.py",
        "analyze_quercetin_contacts.py"
    ]
    
    print("\n--- Organizando Scripts de Diagnóstico ---")
    for scr in diagnostics_to_move:
        src = os.path.join(project_root, "scripts", scr)
        dst = os.path.join(project_root, "scripts", "diagnostics", scr)
        if os.path.exists(src):
            shutil.move(src, dst)
            print(f"  [MOVED] {scr} -> scripts/diagnostics/")
            
    # 6. Organizar docking/vina_output/ por diana
    print("\n--- Organizando Resultados de Docking ---")
    output_dir = os.path.join(project_root, "docking", "vina_output")
    if os.path.exists(output_dir):
        files = os.listdir(output_dir)
        for f in files:
            src = os.path.join(output_dir, f)
            if os.path.isdir(src):
                continue
                
            # Determinar diana por prefijo
            diana_dest = None
            if f.startswith("GyrB_"):
                diana_dest = "GyrB"
            elif f.startswith("MurG_"):
                diana_dest = "MurG"
            elif f.startswith("PBP2a_"):
                diana_dest = "PBP2a"
                
            if diana_dest:
                dst = os.path.join(output_dir, diana_dest, f)
                shutil.move(src, dst)
                print(f"  [MOVED] {f} -> docking/vina_output/{diana_dest}/")
                
    print("\n[ÉXITO] Espacio de trabajo reorganizado limpiamente.")

if __name__ == "__main__":
    main()
