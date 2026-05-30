#!/usr/bin/env python3
import json
import os

def main():
    notebook_path = "notebooks/dinamica_colab.ipynb"
    if not os.path.exists(notebook_path):
        print(f"Error: {notebook_path} not found.")
        return
        
    with open(notebook_path, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    updated = False
    for cell in data.get("cells", []):
        if cell.get("cell_type") == "code":
            source = cell.get("source", [])
            # Search for the cell with condacolab.check()
            if any("condacolab.check()" in line for line in source):
                new_source = []
                for line in source:
                    new_source.append(line)
                    if "condacolab.check()" in line:
                        # Append the fix lines
                        new_source.append("\n")
                        new_source.append("# Eliminar el archivo de anclaje de conda para evitar errores de incompatibilidad con Colab\n")
                        new_source.append("!rm -f /usr/local/conda-meta/pinned\n")
                cell["source"] = new_source
                updated = True
                print("  [OK] Se actualizó con éxito la celda de condacolab en el cuaderno.")
                break
                
    if updated:
        with open(notebook_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
            # Add a trailing newline as per standard json formatting
            f.write("\n")
        print("  [SUCCESS] Archivo notebooks/dinamica_colab.ipynb modificado de forma segura.")
    else:
        print("  [WARNING] No se encontró la celda objetivo en el cuaderno.")

if __name__ == "__main__":
    main()
