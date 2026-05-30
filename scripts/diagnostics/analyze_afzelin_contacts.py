#!/usr/bin/env python3
import os
import math

def parse_pdbqt_ligand_pose1(ligand_path):
    coords = []
    if not os.path.exists(ligand_path):
        print(f"Error: {ligand_path} not found.")
        return None
        
    with open(ligand_path, "r") as f:
        in_model_1 = False
        for line in f:
            if line.startswith("MODEL 1"):
                in_model_1 = True
                continue
            if line.startswith("ENDMDL"):
                if in_model_1:
                    break
            if in_model_1 and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue
    return coords

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, os.pardir))
    
    receptor_path = os.path.join(project_root, "data", "receptors", "MurG_AF-Q6GGZ0-F1.pdb")
    afzelin_path = os.path.join(project_root, "docking", "vina_output", "MurG_Afzelin_3D_out.pdbqt")
    
    # 1. Parse docked Afzelin pose 1 coords
    lig_coords = parse_pdbqt_ligand_pose1(afzelin_path)
    if not lig_coords:
        print("Could not parse Afzelin pose 1.")
        return
        
    print(f"Parsed {len(lig_coords)} atoms from Afzelin's best docked pose.")
    
    # 2. Find close residues in receptor
    contact_residues = {}
    
    with open(receptor_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    atom_name = line[12:16].strip()
                    
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    # Check distance to any ligand atom
                    for lx, ly, lz in lig_coords:
                        dist = math.sqrt((x - lx)**2 + (y - ly)**2 + (z - lz)**2)
                        if dist <= 4.5:
                            if res_num not in contact_residues:
                                contact_residues[res_num] = {
                                    "name": res_name,
                                    "min_dist": dist,
                                    "closest_atom": atom_name,
                                    "contacts_count": 1
                                }
                            else:
                                contact_residues[res_num]["contacts_count"] += 1
                                if dist < contact_residues[res_num]["min_dist"]:
                                    contact_residues[res_num]["min_dist"] = dist
                                    contact_residues[res_num]["closest_atom"] = atom_name
                except (ValueError, IndexError):
                    continue
                    
    print("\n--- Residues in MurG forming contacts (<= 4.5 A) with our Afzelin (Best Pose) ---")
    print(f"{'Residue':<10} {'Name':<6} {'Min Distance (A)':<18} {'Closest Atom':<15} {'Total Contacts'}")
    print("-" * 65)
    for res_num in sorted(contact_residues.keys()):
        info = contact_residues[res_num]
        print(f"{res_num:<10} {info['name']:<6} {info['min_dist']:<18.3f} {info['closest_atom']:<15} {info['contacts_count']}")

if __name__ == "__main__":
    main()
