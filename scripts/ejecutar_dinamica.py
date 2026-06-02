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
  Requiere GROMACS compilado con CUDA (imagen nvcr.io/hpc/gromacs:2023.3).
  Los flags GPU activan offloading de: non-bonded, PME, bonded, integrador.
"""

import os
import sys
import shutil
import subprocess

# =====================================================================
# Flags GPU para gmx mdrun
# =====================================================================
GPU_FLAGS_FULL  = "-ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu -gpu_id 0"
GPU_FLAGS_MINIM = "-ntmpi 1 -ntomp 8 -nb gpu -gpu_id 0"


# =====================================================================
# Separar proteína y ligando desde el PDB del complejo
# =====================================================================

def separar_proteina_ligando(complex_pdb, protein_out, ligand_out, resname):
    """
    Lee el PDB del complejo y escribe dos archivos separados:
      - protein_out : solo líneas ATOM (proteína estándar)
      - ligand_out  : solo líneas HETATM del ligando con el resname dado
    pdb2gmx solo acepta residuos estándar; el ligando debe excluirse.
    """
    if not os.path.exists(complex_pdb):
        print(f"[ERROR] No se encontró el PDB del complejo: {complex_pdb}")
        sys.exit(1)

    protein_lines = []
    ligand_lines  = []

    with open(complex_pdb, 'r') as f:
        for line in f:
            rec = line[:6].strip()
            if rec == 'ATOM':
                protein_lines.append(line)
            elif rec == 'HETATM':
                # Extraer nombre de residuo (columnas 17-20)
                res = line[17:21].strip()
                if res == resname:
                    ligand_lines.append(line)
                # HOH/WAT del cristal se descarta (no se añaden)
            elif rec in ('TER', 'END'):
                protein_lines.append(line)

    if not protein_lines:
        print(f"[ERROR] No se encontraron líneas ATOM en {complex_pdb}")
        sys.exit(1)

    with open(protein_out, 'w') as f:
        f.writelines(protein_lines)
        if not protein_lines[-1].startswith('END'):
            f.write('END\n')

    with open(ligand_out, 'w') as f:
        f.writelines(ligand_lines)
        if ligand_lines:
            f.write('END\n')

    print(f"  [OK] Proteína extraida: {len(protein_lines)} líneas → {protein_out}")
    print(f"  [OK] Ligando extraido : {len(ligand_lines)} líneas → {ligand_out}")

    if not ligand_lines:
        print(f"  [WARNING] No se encontraron HETATM con resname '{resname}' en el PDB.")
        print(f"            Verifica que el nombre del residuo del ligando sea correcto.")


# =====================================================================
# Fusión de archivos .gro
# =====================================================================

def fusionar_gro(protein_gro, ligand_gro, output_gro):
    if not os.path.exists(protein_gro):
        print(f"[ERROR] No se encontró {protein_gro}")
        sys.exit(1)
    if not os.path.exists(ligand_gro):
        print(f"[ERROR] No se encontró {ligand_gro}")
        sys.exit(1)

    with open(protein_gro, 'r') as f:
        p_lines = f.readlines()
    with open(ligand_gro, 'r') as f:
        l_lines = f.readlines()

    p_atoms = int(p_lines[1].strip())
    l_atoms = int(l_lines[1].strip())
    total_atoms = p_atoms + l_atoms

    new_lines = []
    new_lines.append("Complejo Proteina-Ligando Fusionado por Semillero\n")
    new_lines.append(f"{total_atoms:>5}\n")
    for i in range(2, len(p_lines) - 1):
        new_lines.append(p_lines[i])
    for i in range(2, len(l_lines) - 1):
        new_lines.append(l_lines[i])
    new_lines.append(p_lines[-1])

    with open(output_gro, 'w') as f:
        f.writelines(new_lines)
    print(f"  [OK] Fusión .gro completada: {output_gro} ({total_atoms} átomos)")


# =====================================================================
# Actualizar topol.top con el ligando
# =====================================================================

def actualizar_topol(topol_file, itp_file, ligand_resname):
    if not os.path.exists(topol_file):
        print(f"[ERROR] {topol_file} no existe.")
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
    print(f"  [OK] topol.top actualizado con ligando {ligand_resname}")


# =====================================================================
# Parámetros .mdp
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
nsteps      = 50000
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
nsteps      = 50000
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
# Pipeline principal
# =====================================================================

def main():
    if len(sys.argv) < 2:
        print("\n[ERROR] Debes especificar el complejo como argumento.")
        print("Uso: python scripts/ejecutar_dinamica.py <COMPLEJO> [TIEMPO_NS]")
        sys.exit(1)

    complejo = sys.argv[1]
    tiempo_ns = 10.0
    if len(sys.argv) >= 3:
        try:
            tiempo_ns = float(sys.argv[2])
        except ValueError:
            print("[WARNING] Tiempo inválido, usando 10 ns.")

    nsteps_prod = int((tiempo_ns * 1000 * 1000) / 2)

    sistemas = {
        "MurG_Afzelin": {
            "complex_pdb": "data/complexes/complex_MurG_Afzelin.pdb",
            "ligand_sdf":  "data/ligands/Afzelin_3D.sdf",
            "charge": 0,
            "resname": "AFZ"
        },
        "MurG_Quercetin": {
            "complex_pdb": "data/complexes/complex_MurG_Quercetin.pdb",
            "ligand_sdf":  "data/ligands/Quercetin_3D.sdf",
            "charge": 0,
            "resname": "QUE"
        },
        "PBP2a_Afzelin": {
            "complex_pdb": "data/complexes/complex_PBP2a_Afzelin.pdb",
            "ligand_sdf":  "data/ligands/Afzelin_3D.sdf",
            "charge": 0,
            "resname": "AFZ"
        },
        "PBP2a_Ceftaroline": {
            "complex_pdb": "data/complexes/complex_PBP2a_Ceftaroline.pdb",
            "ligand_sdf":  "data/ligands/Ceftaroline_3D.sdf",
            "charge": 0,
            "resname": "CEF"
        }
    }

    if complejo not in sistemas:
        print(f"\n[ERROR] Complejo '{complejo}' no soportado.")
        print(f"Válidos: {', '.join(sistemas.keys())}")
        sys.exit(1)

    sys_info = sistemas[complejo]

    if not shutil.which("gmx"):
        print("\n[ERROR] GROMACS ('gmx') no encontrado en PATH.")
        sys.exit(1)
    if not shutil.which("acpype"):
        print("\n[ERROR] ACPYPE no encontrado en PATH.")
        sys.exit(1)

    for key, path in [("complex_pdb", sys_info["complex_pdb"]), ("ligand_sdf", sys_info["ligand_sdf"])]:
        if not os.path.exists(path):
            print(f"\n[ERROR] Archivo no encontrado: {path}")
            sys.exit(1)

    print("\n" + "=" * 65)
    print(f"  INICIANDO PIPELINE DE DINÁMICA MOLECULAR (GPU-ENABLED)")
    print(f"  Sistema  : {complejo}")
    print(f"  Tiempo   : {tiempo_ns} ns ({nsteps_prod} pasos)")
    print(f"  Complejo : {sys_info['complex_pdb']}")
    print(f"  Ligando  : {sys_info['ligand_sdf']}")
    print(f"  GPU flags: {GPU_FLAGS_FULL}")
    print("=" * 65)

    run_dir = f"md_run_{complejo}"
    os.makedirs(run_dir, exist_ok=True)
    print(f"\n[*] Carpeta de trabajo: {run_dir}")

    resname = sys_info["resname"]
    shutil.copy(sys_info["complex_pdb"], os.path.join(run_dir, "complex_input.pdb"))
    shutil.copy(sys_info["ligand_sdf"],  os.path.join(run_dir, f"{resname}.sdf"))
    os.chdir(run_dir)

    for fn, content in [("minim.mdp", MINIM_MDP), ("nvt.mdp", NVT_MDP), ("npt.mdp", NPT_MDP)]:
        with open(fn, "w") as f:
            f.write(content)
    with open("md.mdp", "w") as f:
        f.write(generar_md_mdp(nsteps_prod))
    print("[*] Archivos .mdp creados.")

    # =====================================================================
    # PASO A: Parametrización del ligando (ACPYPE)
    # =====================================================================
    print("\n--- PASO A: Parametrizando ligando con ACPYPE (GAFF2/AM1-BCC) ---")
    res_acpype = subprocess.run(
        f"acpype -i {resname}.sdf -c bcc -n {sys_info['charge']} -f",
        shell=True, capture_output=True, text=True
    )
    ligand_folder = f"{resname}.acpype"
    ligand_gro = f"{ligand_folder}/{resname}_GMX.gro"
    ligand_itp = f"{ligand_folder}/{resname}_GMX.itp"
    if not os.path.exists(ligand_gro):
        print("[ERROR] Fallo en ACPYPE:")
        print(res_acpype.stderr)
        sys.exit(1)
    print("[*] Ligando parametrizado.")

    # =====================================================================
    # PASO B: Separar proteína del PDB del complejo
    # pdb2gmx no reconoce HETATM de ligandos — deben excluirse
    # =====================================================================
    print("\n--- PASO B: Separando proteína y ligando del PDB del complejo ---")
    separar_proteina_ligando("complex_input.pdb", "protein_only.pdb", "ligand_from_complex.pdb", resname)

    # =====================================================================
    # PASO C: Topología de la proteína (pdb2gmx)
    # =====================================================================
    print("\n--- PASO C: Generando topología de la proteína (AMBER99SB-ILDN) ---")
    pdb2gmx_cmd = "gmx pdb2gmx -f protein_only.pdb -o protein_processed.gro -water tip3p -ignh -p topol.top <<EOF\n6\nEOF"
    res_pdb = subprocess.run(pdb2gmx_cmd, shell=True, capture_output=True, text=True)
    if not os.path.exists("protein_processed.gro"):
        print("[ERROR] Fallo en pdb2gmx:")
        print(res_pdb.stderr)
        sys.exit(1)
    print("[*] Topología de proteína lista.")

    # =====================================================================
    # PASO D: Fusión del complejo
    # =====================================================================
    print("\n--- PASO D: Fusionando coordenadas y topología del complejo ---")
    fusionar_gro("protein_processed.gro", ligand_gro, "complex.gro")
    actualizar_topol("topol.top", ligand_itp, resname)
    print("[*] Complejo fusionado.")

    # =====================================================================
    # PASO E: Solvatación
    # =====================================================================
    print("\n--- PASO E: Solvatando el sistema (caja dodecaédrica, 1.2 nm buffer) ---")
    run_cmd("gmx editconf -f complex.gro -o complex_box.gro -c -d 1.2 -bt dodecahedron", "Fallo en editconf")
    run_cmd("gmx solvate -cs spc216.gro -cp complex_box.gro -o complex_solv.gro -p topol.top", "Fallo en solvate")
    print("[*] Sistema solvatado.")

    # =====================================================================
    # PASO F: Neutralización
    # =====================================================================
    print("\n--- PASO F: Neutralizando con Na+/Cl- a 0.15 M ---")
    run_cmd("gmx grompp -f minim.mdp -c complex_solv.gro -p topol.top -o ions.tpr -maxwarn 1", "Fallo en grompp para genion")
    run_cmd("echo 'SOL' | gmx genion -s ions.tpr -o complex_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15", "Fallo en genion")
    print("[*] Sistema neutralizado.")

    # =====================================================================
    # PASO G: Minimización de energía
    # =====================================================================
    print("\n--- PASO G: Minimización de energía (GPU parcial) ---")
    run_cmd("gmx grompp -f minim.mdp -c complex_solv_ions.gro -p topol.top -o em.tpr", "Fallo en grompp minimización")
    run_cmd(f"gmx mdrun -v -deffnm em {GPU_FLAGS_MINIM}", "Fallo en minimización", silence=False)
    print("[*] Minimización completada.")

    # =====================================================================
    # PASO H: Equilibración NVT
    # =====================================================================
    print("\n--- PASO H: Equilibración NVT (100 ps, 300 K) ---")
    run_cmd("gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr", "Fallo en grompp NVT")
    run_cmd(f"gmx mdrun -deffnm nvt {GPU_FLAGS_FULL}", "Fallo en NVT", silence=False)
    print("[*] NVT completado.")

    # =====================================================================
    # PASO I: Equilibración NPT
    # =====================================================================
    print("\n--- PASO I: Equilibración NPT (100 ps, 1 bar) ---")
    run_cmd("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr", "Fallo en grompp NPT")
    run_cmd(f"gmx mdrun -deffnm npt {GPU_FLAGS_FULL}", "Fallo en NPT", silence=False)
    print("[*] NPT completado.")

    # =====================================================================
    # PASO J: Producción MD
    # =====================================================================
    print(f"\n--- PASO J: Producción MD ({tiempo_ns} ns) ---")
    run_cmd("gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_production.tpr", "Fallo en grompp producción")
    res_md = subprocess.run(f"gmx mdrun -deffnm md_production {GPU_FLAGS_FULL}", shell=True)

    if res_md.returncode == 0:
        print(f"\n[ÉXITO] Simulación de {tiempo_ns} ns completada!")
        print(f"Resultados en: md_run_{complejo}/")
    else:
        print("\n[WARNING] Producción interrumpida. Revisa md_production.log.")


if __name__ == "__main__":
    main()
