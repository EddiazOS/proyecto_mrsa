#!/usr/bin/env python3
import os

def extract_ligand_and_clean_receptor(pdb_path, residue_name, ligand_output_path, receptor_output_path):
    """
    Parses a PDB file, extracts atoms belonging to the specified residue_name (ligand)
    into ligand_output_path, and saves the protein-only atoms (ATOM/TER) into receptor_output_path.
    Also preserves relevant CONECT records for the ligand.
    """
    if not os.path.exists(pdb_path):
        print(f"Error: {pdb_path} no existe.")
        return False

    print(f"Procesando {pdb_path}...")
    ligand_lines = []
    receptor_lines = []
    ligand_atom_serials = set()

    with open(pdb_path, 'r') as f:
        lines = f.readlines()

    # Primera pasada: Extraer átomos
    for line in lines:
        if line.startswith("ATOM") or line.startswith("TER"):
            receptor_lines.append(line)
        elif line.startswith("HETATM"):
            # El nombre del residuo está en las columnas 17-20 (17, 18, 19, 20 en 1-indexed, es decir, cols 17-20)
            res_name_field = line[17:21].strip()
            if res_name_field == residue_name:
                ligand_lines.append(line)
                # El número de serie del átomo está en las columnas 6-11 (índice 6:11)
                try:
                    atom_serial = int(line[6:11].strip())
                    ligand_atom_serials.add(atom_serial)
                except ValueError:
                    pass
            # Las aguas e iones no se agregan a la proteína limpia

    # Segunda pasada: Extraer registros CONECT para el ligando
    conect_lines = []
    for line in lines:
        if line.startswith("CONECT"):
            parts = line.split()
            # Si todos los números en el registro CONECT corresponden a nuestros átomos de ligando, lo guardamos
            try:
                serials = [int(p) for p in parts[1:]]
                if any(s in ligand_atom_serials for s in serials):
                    # Solo nos quedamos con las conexiones que involucren átomos del ligando
                    conect_lines.append(line)
            except ValueError:
                pass

    # Escribir el archivo del ligando
    if ligand_lines:
        with open(ligand_output_path, 'w') as f_lig:
            f_lig.write(f"REMARK   6 EXTRAIDO AUTOMATICAMENTE POR ANTIGRAVITY\n")
            f_lig.write(f"REMARK   6 LIGANDO DE CONTROL: {residue_name}\n")
            for line in ligand_lines:
                f_lig.write(line)
            for line in conect_lines:
                f_lig.write(line)
            f_lig.write("END\n")
        print(f"  -> Ligando {residue_name} extraído en {ligand_output_path} ({len(ligand_lines)} átomos)")
    else:
        print(f"  [ADVERTENCIA] No se encontraron átomos para el ligando {residue_name}")

    # Escribir el archivo del receptor limpio
    if receptor_lines:
        with open(receptor_output_path, 'w') as f_rec:
            f_rec.write(f"REMARK   6 RECEPTOR LIMPIADO AUTOMATICAMENTE POR ANTIGRAVITY\n")
            f_rec.write(f"REMARK   6 SIN SOLVENTE NI DIANAS CO-CRISTALIZADAS DE HETATM\n")
            for line in receptor_lines:
                f_rec.write(line)
            f_rec.write("END\n")
        print(f"  -> Receptor purificado en {receptor_output_path} ({len(receptor_lines)} líneas)")
    else:
        print(f"  [ERROR] No se pudieron extraer coordenadas para el receptor limpio.")

    return True

if __name__ == "__main__":
    # Rutas por defecto asumiendo que el script se corre en la estructura organizada
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # 1. PBP2a (3ZG0) -> Ligando AI8 (Ceftarolina)
    pbp2a_in = os.path.join(base_dir, "data", "receptors", "PBP2a_3ZG0.pdb")
    pbp2a_lig_out = os.path.join(base_dir, "data", "ligands", "Ceftaroline_3D.pdb")
    pbp2a_rec_out = os.path.join(base_dir, "data", "receptors", "PBP2a_3ZG0_clean.pdb")
    extract_ligand_and_clean_receptor(pbp2a_in, "AI8", pbp2a_lig_out, pbp2a_rec_out)

    # 2. GyrB (3TTZ) -> Ligando 07N (Pirrolamida 07N)
    gyrb_in = os.path.join(base_dir, "data", "receptors", "GyrB_3TTZ.pdb")
    gyrb_lig_out = os.path.join(base_dir, "data", "ligands", "07N_GyrB_3D.pdb")
    gyrb_rec_out = os.path.join(base_dir, "data", "receptors", "GyrB_3TTZ_clean.pdb")
    extract_ligand_and_clean_receptor(gyrb_in, "07N", gyrb_lig_out, gyrb_rec_out)
