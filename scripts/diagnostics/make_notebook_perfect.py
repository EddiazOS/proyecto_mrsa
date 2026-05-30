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
        
    new_cells = []
    
    # Vamos a procesar cada celda y automatizarla
    for cell in data.get("cells", []):
        cell_type = cell.get("cell_type")
        source = cell.get("source", [])
        
        # 1. Inyectar la celda de configuración al principio del Paso 1
        if cell_type == "markdown" and any("## Paso 1: Configuración del Entorno" in line for line in source):
            new_cells.append(cell)
            
            # Crear celda de configuración con dropdown para Colab
            config_cell = {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# =====================================================================\n",
                    "# CONFIGURACIÓN DEL COMPLEJO A SIMULAR (Google Colab Form)\n",
                    "# =====================================================================\n",
                    "# Selecciona el complejo de interés de la lista desplegable:\n",
                    "COMPLEJO = \"MurG_Afzelin\"  # @param [\"MurG_Afzelin\", \"MurG_Quercetin\", \"PBP2a_Afzelin\", \"PBP2a_Ceftaroline\"]\n",
                    "\n",
                    "print(f\"[*] Sistema seleccionado: {COMPLEJO}\")\n",
                    "\n",
                    "# Definir rutas y variables dinámicas según la selección\n",
                    "if COMPLEJO == \"MurG_Afzelin\":\n",
                    "    protein_file = \"data/receptors/MurG_AF-Q6GGZ0-F1.pdb\"\n",
                    "    ligand_file = \"data/ligands/Afzelin_3D.sdf\"\n",
                    "    ligand_name = \"Afzelin_3D\"\n",
                    "    charge = 0\n",
                    "    resname = \"AFZ\"\n",
                    "elif COMPLEJO == \"MurG_Quercetin\":\n",
                    "    protein_file = \"data/receptors/MurG_AF-Q6GGZ0-F1.pdb\"\n",
                    "    ligand_file = \"data/ligands/Quercetin_3D.sdf\"\n",
                    "    ligand_name = \"Quercetin_3D\"\n",
                    "    charge = 0\n",
                    "    resname = \"QUE\"\n",
                    "elif COMPLEJO == \"PBP2a_Afzelin\":\n",
                    "    protein_file = \"data/receptors/PBP2a_3ZG0_clean.pdb\"\n",
                    "    ligand_file = \"data/ligands/Afzelin_3D.sdf\"\n",
                    "    ligand_name = \"Afzelin_3D\"\n",
                    "    charge = 0\n",
                    "    resname = \"AFZ\"\n",
                    "elif COMPLEJO == \"PBP2a_Ceftaroline\":\n",
                    "    protein_file = \"data/receptors/PBP2a_3ZG0_clean.pdb\"\n",
                    "    ligand_file = \"data/ligands/Ceftaroline_3D.sdf\"\n",
                    "    ligand_name = \"Ceftaroline_3D\"\n",
                    "    charge = 0\n",
                    "    resname = \"CEF\"\n",
                    "\n",
                    "print(f\"  Proteína : {protein_file}\")\n",
                    "print(f\"  Ligando  : {ligand_file}\")\n",
                    "print(f\"  Residuo  : {resname}\")\n"
                ]
            }
            new_cells.append(config_cell)
            continue
            
        # 2. Automatizar Paso 2 (acpype)
        if cell_type == "code" and any("print(\"=== GENERANDO TOPOLOGÍA PARA AFZELINA ===\")" in line for line in source):
            # Reemplazar por el comando automatizado usando las variables de configuración
            automated_acpype = {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Ejecutar ACPYPE dinámicamente según el complejo seleccionado\n",
                    "import subprocess\n",
                    "print(f\"[*] Ejecutando ACPYPE para: {ligand_name} (Carga: {charge})...\")\n",
                    "\n",
                    "# Comando acpype automatizado\n",
                    "cmd = f\"acpype -i {ligand_file} -c bcc -n {charge} -m amber -f\"\n",
                    "subprocess.run(cmd, shell=True)\n",
                    "print(\"[*] Topología del ligando generada exitosamente.\")\n"
                ]
            }
            new_cells.append(automated_acpype)
            continue
            
        # 3. Automatizar Paso 3 (pdb2gmx)
        if cell_type == "code" and any("!gmx pdb2gmx" in line for line in source):
            automated_pdb2gmx = {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Generar topología de la proteína dinámicamente en GROMACS\n",
                    "print(f\"[*] Ejecutando pdb2gmx sobre: {protein_file}...\")\n",
                    "!gmx pdb2gmx -f {protein_file} -o protein_processed.gro -water tip3p -ignh -p topol.top\n"
                ]
            }
            new_cells.append(automated_pdb2gmx)
            continue
            
        # 4. Automatizar Paso 4 (fusionar y actualizar topol)
        if cell_type == "code" and any("# Ejemplo de uso de la fusión en Colab:" in line for line in source):
            automated_fusion = {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# Ejecutar la fusión de la estructura y la topología dinámicamente\n",
                    "print(f\"[*] Fusionando coordenadas de proteína y ligando ({resname})...\")\n",
                    "protein_gro_input = \"protein_processed.gro\"\n",
                    "ligand_gro_input = f\"{ligand_name}.acpype/{ligand_name}_GMX.gro\"\n",
                    "ligand_itp_input = f\"{ligand_name}.acpype/{ligand_name}_GMX.itp\"\n",
                    "\n",
                    "fusionar_gro(protein_gro_input, ligand_gro_input, \"complex.gro\")\n",
                    "actualizar_topol(\"topol.top\", ligand_itp_input, resname)\n",
                    "print(\"[*] Complejo complex.gro e itp acoplados con éxito.\")\n"
                ]
            }
            new_cells.append(automated_fusion)
            continue
            
        # Mantener las demás celdas tal cual
        new_cells.append(cell)
        
    data["cells"] = new_cells
    
    with open(notebook_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")
        
    print("  [SUCCESS] El cuaderno notebooks/dinamica_colab.ipynb ha sido automatizado al 100% con un selector de complejos interactivo.")

if __name__ == "__main__":
    main()
