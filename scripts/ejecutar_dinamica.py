#!/usr/bin/env python3
"""
ejecutar_dinamica.py
====================
Script de automatización para ejecutar dinámicas moleculares de GROMACS.
Soporta ejecución con aceleración GPU completa (RunPod / Google Colab con GPU).

Uso:
  python scripts/ejecutar_dinamica.py <COMPLEJO> <TIEMPO_NS>

Sistemas soportados:
  - MurG_Afzelin
  - MurG_Quercetin
  - PBP2a_Afzelin
  - PBP2a_Ceftaroline

Ejemplos:
  python scripts/ejecutar_dinamica.py MurG_Afzelin 10
  python scripts/ejecutar_dinamica.py PBP2a_Afzelin 100

Notas GPU:
  Requiere GROMACS compilado con CUDA (imagen nvcr.io/hpc/gromacs:2024.1).
  Los flags GPU activan offloading de: non-bonded, PME, bonded, integrador.
"""

import os
import sys
import shutil
import subprocess

# =====================================================================
# Flags GPU para gmx mdrun
# -ntmpi 1     : 1 proceso MPI (1 GPU)
# -ntomp 8     : 8 threads OpenMP por proceso
# -nb gpu      : fuerzas no-enlazadas en GPU
# -pme gpu     : electrostáticas PME en GPU
# -bonded gpu  : fuerzas enlazadas en GPU
# -update gpu  : integrador de movimiento en GPU (solo integradores 'md')
# -gpu_id 0    : usar la primera GPU disponible
# =====================================================================
GPU_FLAGS_FULL  = "-ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu -gpu_id 0"
GPU_FLAGS_MINIM = "-ntmpi 1 -ntomp 8 -nb gpu -gpu_id 0"
# Nota: -update gpu y -bonded gpu NO aplican al integrador 'steep' (minimización)

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

    p_atoms = int(p_lines[1].strip())
    l_atoms = int(l_lines[1].strip())
    total_atoms = p_atoms + l_atoms

    new_lines = []
    new_lines.append(f"Complejo Proteina-Ligando Fusionado por Semillero\n")
    new_lines.append(f"{total_atoms:>5}\n")
    for i in range(2, len(p_lines) - 1):
        new_lines.append(p_lines[i])
    for i in range(2, len(l_lines) - 1):
        new_lines.append(l_lines[i])
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
        if 'forcefield.itp' in line and not inserted_itp:
            new_lines.append(f'#include "{itp_file}"\n')
            inserted_itp = True

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
gen_vel     = yes
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
refcoord-scaling = com
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

    tiempo_ns = 10.0
    if len(sys.argv) >= 3:
        try:
            tiempo_ns = float(sys.argv[2])
        except ValueError:
            print("[WARNING] Tiempo de simulación inválido, usando por defecto 10 ns.")

    nsteps_prod = int((tiempo_ns * 1000 * 1000) / 2)  # dt = 2 fs = 0.002 ps

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

    if not shutil.which("gmx"):
        print("\n[ERROR] GROMACS ('gmx') no se encuentra en el PATH actual.")
        print("Asegúrate de usar la imagen Docker nvcr.io/hpc/gromacs:2024.1 en RunPod.")
        sys.exit(1)

    if not shutil.which("acpype"):
        print("\n[ERROR] ACPYPE ('acpype') no se encuentra en el PATH actual.")
        sys.exit(1)

    print("\n" + "=" * 65)
    print(f"  INICIANDO PIPELINE DE DINÁMICA MOLECULAR (GPU-ENABLED)")
    print(f"  Sistema  : {complejo}")
    print(f"  Tiempo   : {tiempo_ns} ns ({nsteps_prod} pasos)")
    print(f"  GPU flags: {GPU_FLAGS_FULL}")
    print("=" * 65)

    run_dir = f"md_run_{complejo}"
    os.makedirs(run_dir, exist_ok=True)
    print(f"\n[*] Carpeta de trabajo creada: {run_dir}")

    resname = sys_info["resname"]
    shutil.copy(sys_info["protein"], os.path.join(run_dir, "receptor.pdb"))
    shutil.copy(sys_info["ligand"], os.path.join(run_dir, f"{resname}.sdf"))
    os.chdir(run_dir)

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
    acpype_cmd = f"acpype -i {resname}.sdf -c bcc -n {sys_info['charge']} -f"
    res_acpype = subprocess.run(acpype_cmd, shell=True, capture_output=True, text=True)

    ligand_folder = f"{resname}.acpype"
    ligand_gro = f"{ligand_folder}/{resname}_GMX.gro"
    ligand_itp = f"{ligand_folder}/{resname}_GMX.itp"

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
    res_pdb = subprocess.run(pdb2gmx_cmd, shell=True, capture_output=True, text=True)

    if not os.path.exists("protein_processed.gro"):
        print("[ERROR] Fallo en pdb2gmx de GROMACS:")
        print(res_pdb.stderr)
        sys.exit(1)
    print("[*] Topología y coordenadas de la proteína listas.")

    # =====================================================================
    # PASO C: Fusión del Complejo
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

    print("\n--- PASO E: Neutralizando el sistema con Na+/Cl- a 0.15 M ---")
    run_cmd("gmx grompp -f minim.mdp -c complex_solv.gro -p topol.top -o ions.tpr -maxwarn 1", "Fallo al compilar mdp para genion")
    genion_cmd = "echo 'SOL' | gmx genion -s ions.tpr -o complex_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15"
    run_cmd(genion_cmd, "Fallo al agregar iones con gmx genion")
    print("[*] Sistema neutralizado.")

    # =====================================================================
    # PASO F: Minimización de Energía (GPU: nb solamente)
    # Nota: integrador 'steep' no soporta -update gpu ni -bonded gpu
    # =====================================================================
    print("\n--- PASO F: Ejecutando Minimización de Energía (EM) con GPU ---")
    run_cmd("gmx grompp -f minim.mdp -c complex_solv_ions.gro -p topol.top -o em.tpr", "Fallo al compilar mdp para minimización")
    print("[INFO] Corriendo mdrun de minimización con aceleración GPU parcial...")
    run_cmd(f"gmx mdrun -v -deffnm em {GPU_FLAGS_MINIM}", "Fallo al correr minimización de energía", silence=False)
    print("[*] Minimización energética completada.")

    # =====================================================================
    # PASO G: Equilibración NVT (GPU completo)
    # =====================================================================
    print("\n--- PASO G: Ejecutando Equilibración NVT (100 ps, 300 K) con GPU ---")
    run_cmd("gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr", "Fallo al compilar mdp para NVT")
    print("[INFO] Corriendo mdrun NVT con offloading GPU completo...")
    run_cmd(f"gmx mdrun -deffnm nvt {GPU_FLAGS_FULL}", "Fallo al correr equilibración NVT", silence=False)
    print("[*] Equilibración NVT completada.")

    # =====================================================================
    # PASO H: Equilibración NPT (GPU completo)
    # =====================================================================
    print("\n--- PASO H: Ejecutando Equilibración NPT (100 ps, 1 bar) con GPU ---")
    run_cmd("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr", "Fallo al compilar mdp para NPT")
    print("[INFO] Corriendo mdrun NPT con offloading GPU completo...")
    run_cmd(f"gmx mdrun -deffnm npt {GPU_FLAGS_FULL}", "Fallo al correr equilibración NPT", silence=False)
    print("[*] Equilibración NPT completada.")

    # =====================================================================
    # PASO I: Producción MD (GPU completo)
    # =====================================================================
    print(f"\n--- PASO I: Iniciando simulación de producción ({tiempo_ns} ns) con GPU ---")
    print(f"[INFO] GPU flags activos: {GPU_FLAGS_FULL}")
    run_cmd("gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_production.tpr", "Fallo al compilar mdp para producción")
    print("[*] Corriendo mdrun de producción con offloading GPU completo...")
    res_md = subprocess.run(f"gmx mdrun -deffnm md_production {GPU_FLAGS_FULL}", shell=True)

    if res_md.returncode == 0:
        print(f"\n[ÉXITO] Simulación de {tiempo_ns} ns completada exitosamente!")
        print(f"Resultados en: md_run_{complejo}/")
    else:
        print("\n[WARNING] La simulación de producción se detuvo o fue interrumpida.")
        print("Revisa md_production.log para diagnóstico.")


if __name__ == "__main__":
    main()
