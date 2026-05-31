#!/usr/bin/env python3
"""
ejecutar_dinamica.py
====================
Script de automatización para ejecutar dinámicas moleculares de GROMACS en Colab.
Permite correr todo el pipeline científico en un solo comando:
  python scripts/ejecutar_dinamica.py <COMPLEJO> <TIEMPO_NS>

Sistemas soportados:
  - MurG_Afzelin
  - MurG_Quercetin
  - PBP2a_Afzelin
  - PBP2a_Ceftaroline

Uso:
  python scripts/ejecutar_dinamica.py MurG_Afzelin 10
  (ejecuta 10 ns de dinámica de producción)
"""

import os
import sys
import shutil
import subprocess

# =====================================================================
# Definición de Funciones de Fusión Estructural y Topológica
# =====================================================================

def fusionar_gro(protein_gro, ligand_gro, output_gro):
    """
    Combina físicamente el archivo .gro de la proteína con el .gro del ligando.
    Ajusta el recuento total de átomos en la línea 2.
    """
    if not os.path.exists(protein_gro):
        print(f"Error: No se encontró {protein_gro}")
        sys.exit(1)
    if not os.path.exists(ligand_gro):
        print(f"Error: No se encontró {ligand_gro}")
        sys.exit(1)
        
    with open(protein_gro, 'r') as f:
        p_lines = f.readlines()
    with open(ligand_gro, 'r') as f:
        l_lines = f.readlines()
        
    # Extraer recuentos de átomos (línea 2, indexada en 1)
    p_atoms = int(p_lines[1].strip())
    l_atoms = int(l_lines[1].strip())
    total_atoms = p_atoms + l_atoms
    
    new_lines = []
    # Línea 1: Título
    new_lines.append(f"Complejo Proteina-Ligando Fusionado por Semillero\n")
    # Línea 2: Número de átomos total
    new_lines.append(f"{total_atoms:>5}\n")
    
    # Añadir los átomos de la proteína (desde la línea 3 hasta la penúltima línea)
    for i in range(2, len(p_lines) - 1):
        new_lines.append(p_lines[i])
        
    # Añadir los átomos del ligando (desde la línea 3 hasta la penúltima línea)
    # Ajustando la numeración atómica si es necesario
    for i in range(2, len(l_lines) - 1):
        new_lines.append(l_lines[i])
        
    # Añadir la última línea (dimensiones de la caja de la proteína)
    new_lines.append(p_lines[-1])
    
    with open(output_gro, 'w') as f:
        f.writelines(new_lines)
    print(f"  [OK] Fusión estructural completada en: {output_gro} (Total Átomos: {total_atoms})")

def actualizar_topol(topol_file, itp_file, ligand_resname):
    """
    Inserta la directiva #include de la topología del ligando en topol.top
    y registra la molécula al final en la sección [ molecules ].
    """
    if not os.path.exists(topol_file):
        print(f"Error: {topol_file} no existe.")
        sys.exit(1)
        
    with open(topol_file, 'r') as f:
        lines = f.readlines()
        
    new_lines = []
    inserted_itp = False
    
    for line in lines:
        new_lines.append(line)
        # Insertar el include del ligando inmediatamente después de la topología de agua/iones
        if 'forcefield.itp' in line and not inserted_itp:
            new_lines.append(f'#include "{itp_file}"\n')
            inserted_itp = True
            
    # Añadir el ligando en la sección final [ molecules ]
    # Verificar si ya existe al final para no duplicar
    has_ligand = any(ligand_resname in line for line in lines[-5:])
    if not has_ligand:
        new_lines.append(f"{ligand_resname:<20} 1\n")
    
    with open(topol_file, 'w') as f:
        f.writelines(new_lines)
    print(f"  [OK] Archivo topol.top actualizado lógicamente con el ligando {ligand_resname}.")

# =====================================================================
# Definición de Parámetros .mdp de GROMACS
# =====================================================================

MINIM_MDP = """
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 1
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
"""

NVT_MDP = """
define      = -DPOSRES
integrator  = md
nsteps      = 50000     ; 100 ps
dt          = 0.002
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
continuation = no
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 10
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
pme_order   = 4
fourierspacing = 0.16
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pcoupl      = no
pbc         = xyz
DispCorr    = EnerPres
velocity_generating = yes
gen_temp    = 300
gen_seed    = -1
"""

NPT_MDP = """
define      = -DPOSRES
integrator  = md
nsteps      = 50000     ; 100 ps
dt          = 0.002
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 10
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
pme_order   = 4
fourierspacing = 0.16
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pcoupl      = C-rescale
pcoupltype  = isotropic
tau_p       = 2.0
compressibility = 4.5e-5
ref_p       = 1.0
pbc         = xyz
DispCorr    = EnerPres
"""

def generar_md_mdp(nsteps):
    return f"""
integrator  = md
nsteps      = {nsteps}
dt          = 0.002
nstxout     = 5000
nstvout     = 5000
nstenergy   = 5000
nstlog      = 5000
nstxout-compressed = 5000
compressed-x-grps = System
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 10
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
pme_order   = 4
fourierspacing = 0.16
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = 300     300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
compressibility = 4.5e-5
ref_p       = 1.0
pbc         = xyz
DispCorr    = EnerPres
"""

def run_cmd(cmd_str, error_msg, silence=True):
    res = subprocess.run(cmd_str, shell=True, capture_output=silence, text=True)
    if res.returncode != 0:
        print(f"\n[ERROR] {error_msg}")
        if silence:
            print("=== STDOUT ===")
            print(res.stdout)
            print("=== STDERR ===")
            print(res.stderr)
        sys.exit(1)
    return res

# =====================================================================
# Orquestación de Ejecución Principal
# =====================================================================

def main():
    if len(sys.argv) < 2:
        print("\n[ERROR] Debes especificar el complejo como argumento.")
        print("Uso:")
        print("  python scripts/ejecutar_dinamica.py <COMPLEJO> [TIEMPO_NS]")
        print("Ejemplos:")
        print("  python scripts/ejecutar_dinamica.py MurG_Afzelin 10")
        print("  python scripts/ejecutar_dinamica.py PBP2a_Afzelin 100")
        sys.exit(1)
        
    complejo = sys.argv[1]
    
    # Tiempo de producción por defecto: 10 ns (5,000,000 pasos de 2 fs)
    tiempo_ns = 10.0
    if len(sys.argv) >= 3:
        try:
            tiempo_ns = float(sys.argv[2])
        except ValueError:
            print("[WARNING] Tiempo de simulación inválido, usando por defecto 10 ns.")
            
    nsteps_prod = int((tiempo_ns * 1000 * 1000) / 2) # dt = 2 fs = 0.002 ps
    
    # 1. Mapeo de archivos y variables
    sistemas = {
        "MurG_Afzelin": {
            "protein": "data/receptors/MurG_AF-Q6GGZ0-F1.pdb",
            "ligand": "data/ligands/Afzelin_3D.sdf",
            "ligand_name": "Afzelin_3D",
            "charge": 0,
            "resname": "AFZ"
        },
        "MurG_Quercetin": {
            "protein": "data/receptors/MurG_AF-Q6GGZ0-F1.pdb",
            "ligand": "data/ligands/Quercetin_3D.sdf",
            "ligand_name": "Quercetin_3D",
            "charge": 0,
            "resname": "QUE"
        },
        "PBP2a_Afzelin": {
            "protein": "data/receptors/PBP2a_3ZG0_clean.pdb",
            "ligand": "data/ligands/Afzelin_3D.sdf",
            "ligand_name": "Afzelin_3D",
            "charge": 0,
            "resname": "AFZ"
        },
        "PBP2a_Ceftaroline": {
            "protein": "data/receptors/PBP2a_3ZG0_clean.pdb",
            "ligand": "data/ligands/Ceftaroline_3D.sdf",
            "ligand_name": "Ceftaroline_3D",
            "charge": 0,
            "resname": "CEF"
        }
    }
    
    if complejo not in sistemas:
        print(f"\n[ERROR] Complejo '{complejo}' no soportado.")
        print(f"Sistemas válidos: {', '.join(sistemas.keys())}")
        sys.exit(1)
        
    sys_info = sistemas[complejo]
    
    # Verificar ejecutabilidad de GROMACS
    if not shutil.which("gmx"):
        print("\n[ERROR] GROMACS ('gmx') no se encuentra en el PATH actual.")
        print("Asegúrate de ejecutar las celdas de instalación de dependencias en Colab.")
        sys.exit(1)
        
    # Verificar ejecutabilidad de ACPYPE
    if not shutil.which("acpype"):
        print("\n[ERROR] ACPYPE ('acpype') no se encuentra en el PATH actual.")
        sys.exit(1)
        
    print("\n" + "=" * 65)
    print(f"  INICIANDO PIPELINE DE DINÁMICA MOLECULAR")
    print(f"  Sistema  : {complejo}")
    print(f"  Tiempo   : {tiempo_ns} ns ({nsteps_prod} pasos)")
    print("=" * 65)
    
    # 2. Crear carpeta de simulación aislada
    run_dir = f"md_run_{complejo}"
    os.makedirs(run_dir, exist_ok=True)
    print(f"\n[*] Carpeta de trabajo creada: {run_dir}")
    
    # Copiar ligandos y proteínas al directorio
    shutil.copy(sys_info["protein"], os.path.join(run_dir, "receptor.pdb"))
    shutil.copy(sys_info["ligand"], os.path.join(run_dir, "ligand.sdf"))
    
    # Movernos al directorio de trabajo
    os.chdir(run_dir)
    
    # 3. Escribir archivos MDP
    with open("minim.mdp", "w") as f:
        f.write(MINIM_MDP)
    with open("nvt.mdp", "w") as f:
        f.write(NVT_MDP)
    with open("npt.mdp", "w") as f:
        f.write(NPT_MDP)
    with open("md.mdp", "w") as f:
        f.write(generar_md_mdp(nsteps_prod))
    print("[*] Archivos de parámetros (.mdp) creados con éxito.")
    
    # =====================================================================
    # PASO A: Parametrización del Ligando (ACPYPE)
    # =====================================================================
    print("\n--- PASO A: Parametrizando el ligando con ACPYPE (GAFF2/AM1-BCC) ---")
    acpype_cmd = f"acpype -i ligand.sdf -c bcc -n {sys_info['charge']} -f"
    res_acpype = subprocess.run(acpype_cmd, shell=True, capture_output=True, text=True)
    
    ligand_folder = f"ligand.acpype"
    ligand_gro = f"{ligand_folder}/ligand_GMX.gro"
    ligand_itp = f"{ligand_folder}/ligand_GMX.itp"
    
    if not os.path.exists(ligand_gro):
        print("[ERROR] Fallo en la parametrización de ACPYPE:")
        print(res_acpype.stderr)
        sys.exit(1)
    print("[*] Ligando parametrizado exitosamente.")
    
    # =====================================================================
    # PASO B: Topología de la Proteína (pdb2gmx)
    # =====================================================================
    print("\n--- PASO B: Generando topología de la proteína (AMBER99SB-ILDN) ---")
    pdb2gmx_cmd = "gmx pdb2gmx -f receptor.pdb -o protein_processed.gro -water tip3p -ignh -p topol.top <<EOF\n6\nEOF"
    # El '6' en la entrada estándar selecciona AMBER99SB-ILDN por defecto en la instalación estándar de conda
    res_pdb = subprocess.run(pdb2gmx_cmd, shell=True, capture_output=True, text=True)
    
    if not os.path.exists("protein_processed.gro"):
        print("[ERROR] Fallo en pdb2gmx de GROMACS:")
        print(res_pdb.stderr)
        sys.exit(1)
    print("[*] Topología y coordenadas de la proteína listas.")
    
    # =====================================================================
    # PASO C: Fusión del Complejo (🧩)
    # =====================================================================
    print("\n--- PASO C: Fusionando coordenadas y topología del complejo ---")
    fusionar_gro("protein_processed.gro", ligand_gro, "complex.gro")
    actualizar_topol("topol.top", ligand_itp, sys_info["resname"])
    print("[*] Complejo fusionado física y lógicamente.")
    
    # =====================================================================
    # PASO D: Solvatación y Neutralización
    # =====================================================================
    print("\n--- PASO D: Solvatando el sistema en caja dodecaédrica (1.2 nm buffer) ---")
    run_cmd("gmx editconf -f complex.gro -o complex_box.gro -c -d 1.2 -bt dodecahedron", "Fallo en gmx editconf")
    run_cmd("gmx solvate -cs spc216.gro -cp complex_box.gro -o complex_solv.gro -p topol.top", "Fallo en gmx solvate")
    print("[*] Sistema solvatado.")
    
    # Agregar iones para neutralizar y llegar a 0.15 M NaCl
    print("\n--- PASO E: Neutralizando el sistema con Na+/Cl- a 0.15 M ---")
    # 1. Compilar mdp temporal para genion
    run_cmd("gmx grompp -f minim.mdp -c complex_solv.gro -p topol.top -o ions.tpr -maxwarn 1", "Fallo al compilar mdp para genion")
    
    # 2. Agregar iones (seleccionando 'SOL' como grupo solvente, típicamente grupo 13 o 15)
    # Mandamos 'SOL' a genion
    genion_cmd = "echo 'SOL' | gmx genion -s ions.tpr -o complex_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15"
    run_cmd(genion_cmd, "Fallo al agregar iones con gmx genion")
    print("[*] Sistema neutralizado.")
    
    # =====================================================================
    # PASO F: Simulaciones (Minimización y Equilibración)
    # =====================================================================
    print("\n--- PASO F: Ejecutando Minimización de Energía (EM) ---")
    run_cmd("gmx grompp -f minim.mdp -c complex_solv_ions.gro -p topol.top -o em.tpr", "Fallo al compilar mdp para minimización")
    run_cmd("gmx mdrun -v -deffnm em", "Fallo al correr minimización de energía")
    print("[*] Minimización energética completada.")
    
    print("\n--- PASO G: Ejecutando Equilibración NVT (100 ps, 300 K) ---")
    run_cmd("gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr", "Fallo al compilar mdp para NVT")
    run_cmd("gmx mdrun -deffnm nvt", "Fallo al correr equilibración NVT")
    print("[*] Equilibración NVT completada.")
    
    print("\n--- PASO H: Ejecutando Equilibración NPT (100 ps, 1 bar) ---")
    run_cmd("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr", "Fallo al compilar mdp para NPT")
    run_cmd("gmx mdrun -deffnm npt", "Fallo al correr equilibración NPT")
    print("[*] Equilibración NPT completada.")
    
    # =====================================================================
    # PASO I: Producción MD
    # =====================================================================
    print(f"\n--- PASO I: Iniciando simulación de producción ({tiempo_ns} ns) ---")
    print("[INFO] Esta corrida puede tomar varias horas. GROMACS usará la GPU CUDA automáticamente.")
    run_cmd("gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_production.tpr", "Fallo al compilar mdp para producción")
    
    # Para la ejecución en Colab en segundo plano o interactiva, se puede correr directamente mdrun
    # Aquí iniciamos el comando mdrun. El usuario verá el archivo de log progresar.
    print("[*] Corriendo mdrun de producción...")
    res_md = subprocess.run("gmx mdrun -deffnm md_production", shell=True)
    
    if res_md.returncode == 0:
        print(f"\n[ÉXITO] ¡Simulación de {tiempo_ns} ns completada exitosamente!")
        print(f"Los resultados se encuentran archivados en la carpeta: md_run_{complejo}/")
    else:
        print("\n[WARNING] La simulación de producción se detuvo o fue interrumpida.")

if __name__ == "__main__":
    main()
