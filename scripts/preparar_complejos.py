#!/usr/bin/env python3
import os
import sys

def get_element(atom_name, pdbqt_type):
    # Try to map from the pdbqt_type at the end
    t = pdbqt_type.strip().upper()
    if t in ['A', 'C', 'CA']:
        return 'C'
    elif t in ['OA', 'O']:
        return 'O'
    elif t in ['HD', 'H']:
        return 'H'
    elif t in ['N', 'NA']:
        return 'N'
    elif t in ['SA', 'S']:
        return 'S'
    elif t in ['CL']:
        return 'CL'
    elif t in ['F']:
        return 'F'
    elif t in ['P']:
        return 'P'
    
    # Fallback to the atom name
    name = atom_name.strip()
    if not name:
        return 'C'
    # Keep only letters at start
    elem = ""
    for char in name:
        if char.isalpha():
            elem += char
          
        else:
            break
    elem = elem.upper()
    if elem in ['CL', 'BR', 'ZN', 'FE', 'MG', 'CA']:
        return elem
    elif elem:
        return elem[0]
    return 'C'

def format_atom_name(name, element):
    # Standard PDB: 4-character atom name (columns 13-16).
    if len(name) >= 4:
        return name[:4]
    
    if len(element) == 1:
        return f" {name:<3}"[:4]
    else:
        return f"{name:<4}"[:4]

def make_pdb_line(record_type, serial, atom_name, res_name, chain, res_seq, x, y, z, occupancy, temp_factor, element):
    record_type = f"{record_type:<6}"[:6]
    serial_str = f"{serial:>5}"[:5]
    res_name = f"{res_name:>3}"[:3]
    chain = f"{chain:>1}"[:1] if chain else " "
    res_seq_str = f"{res_seq:>4}"[:4]
    
    x_str = f"{x:>8.3f}"[:8]
    y_str = f"{y:>8.3f}"[:8]
    z_str = f"{z:>8.3f}"[:8]
    
    occ_str = f"{occupancy:>6.2f}"[:6]
    temp_str = f"{temp_factor:>6.2f}"[:6]
    element_str = f"{element:>2}"[:2]
    
    line = (
        f"{record_type}"
        f"{serial_str}"
        f" "
        f"{atom_name}"
        f" "
        f"{res_name}"
        f" "
        f"{chain}"
        f"{res_seq_str}"
        f" "
        f"   "
        f"{x_str}"
        f"{y_str}"
        f"{z_str}"
        f"{occ_str}"
        f"{temp_str}"
        f"          "
        f"{element_str}"
        f"  "
    )
    return line

def parse_pdbqt_model1(pdbqt_path):
    atoms = []
    if not os.path.exists(pdbqt_path):
        print(f"Error: Docking file not found: {pdbqt_path}")
        sys.exit(1)
        
    with open(pdbqt_path, 'r') as f:
        in_model1 = False
        for line in f:
            if line.startswith('MODEL 1'):
                in_model1 = True
                continue
            elif line.startswith('ENDMDL') and in_model1:
                break
            if in_model1:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    parts = line.split()
                    if len(parts) < 8:
                        continue
                    
                    # Columns in standard PDBQT lines
                    name_str = line[12:16].strip()
                    x_str = line[30:38].strip()
                    y_str = line[38:46].strip()
                    z_str = line[46:54].strip()
                    
                    pdbqt_type = parts[-1]
                    
                    x = float(x_str)
                    y = float(y_str)
                    z = float(z_str)
                    
                    elem = get_element(name_str, pdbqt_type)
                    atoms.append({
                        'element': elem,
                        'x': x,
                        'y': y,
                        'z': z,
                        'original_name': name_str
                    })
                    
    # Assign unique names
    elem_counts = {}
    unique_atoms = []
    for atom in atoms:
        elem = atom['element']
        elem_counts[elem] = elem_counts.get(elem, 0) + 1
        num = elem_counts[elem]
        atom['unique_name'] = f"{elem}{num}"
        unique_atoms.append(atom)
        
    return unique_atoms

def build_complex(protein_pdb, ligand_pdbqt, ligand_resname, output_pdb):
    print(f"\nBuilding complex: {output_pdb}")
    print(f"  Protein: {protein_pdb}")
    print(f"  Ligand: {ligand_pdbqt} (Residue: {ligand_resname})")
    
    if not os.path.exists(protein_pdb):
        print(f"Error: Protein file not found: {protein_pdb}")
        sys.exit(1)
        
    # 1. Parse ligand coordinates (Model 1)
    ligand_atoms = parse_pdbqt_model1(ligand_pdbqt)
    print(f"  Parsed {len(ligand_atoms)} ligand atoms from Model 1 of docking output.")
    
    # 2. Read protein records
    protein_lines = []
    last_serial = 0
    with open(protein_pdb, 'r') as f:
        for line in f:
            # Check if line is ATOM/HETATM to update last_serial
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    serial = int(line[6:11].strip())
                    if serial > last_serial:
                        last_serial = serial
                except ValueError:
                    pass
                protein_lines.append(line.rstrip('\n'))
            elif line.startswith('TER'):
                protein_lines.append(line.rstrip('\n'))
            elif line.startswith('ANISOU'):
                protein_lines.append(line.rstrip('\n'))
            # Ignore END, ENDMDL, or other footers so we can append ligand
            
    print(f"  Read protein structure. Last atom serial: {last_serial}")
    
    # 3. Create complex records
    complex_lines = []
    complex_lines.append(f"REMARK   4 PROTEIN-LIGAND COMPLEX GENERATED BY PREPARAR_COMPLEJOS.PY")
    complex_lines.append(f"REMARK   4 PROTEIN: {os.path.basename(protein_pdb)}")
    complex_lines.append(f"REMARK   4 LIGAND: {os.path.basename(ligand_pdbqt)} (RESTORING AS {ligand_resname})")
    
    # Add protein lines
    complex_lines.extend(protein_lines)
    
    # Check if protein lines already ended with TER, if not add one
    if protein_lines and not protein_lines[-1].startswith('TER'):
        last_ter_serial = last_serial + 1
        # Extract chain from last ATOM line to match chain
        chain_id = 'A'
        for line in reversed(protein_lines):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21]
                break
        complex_lines.append(f"TER   {last_ter_serial:>5}      UNK {chain_id} 999")
        last_serial = last_ter_serial
        
    # Append ligand atoms as HETATM
    ligand_serial = last_serial + 1
    ligand_resseq = 900 # A high residue sequence number to avoid collision
    ligand_chain = 'A' # Standard chain assignment
    
    for atom in ligand_atoms:
        atom_name_formatted = format_atom_name(atom['unique_name'], atom['element'])
        pdb_line = make_pdb_line(
            record_type='HETATM',
            serial=ligand_serial,
            atom_name=atom_name_formatted,
            res_name=ligand_resname,
            chain=ligand_chain,
            res_seq=ligand_resseq,
            x=atom['x'],
            y=atom['y'],
            z=atom['z'],
            occupancy=1.00,
            temp_factor=0.00,
            element=atom['element']
        )
        complex_lines.append(pdb_line)
        ligand_serial += 1
        
    # Add TER and END at the very end
    complex_lines.append(f"TER   {ligand_serial:>5}      {ligand_resname} {ligand_chain} {ligand_resseq}")
    complex_lines.append("END")
    
    # 4. Write complex PDB file
    os.makedirs(os.path.dirname(output_pdb), exist_ok=True)
    with open(output_pdb, 'w') as f:
        for line in complex_lines:
            f.write(line + '\n')
            
    print(f"  Successfully wrote complex file: {output_pdb}")

def main():
    base_dir = "/home/echoes/projects/cursos-semillero/proyecto_mrsa"
    
    # 4 complexes to generate
    jobs = [
        {
            'protein': os.path.join(base_dir, 'data/receptors/MurG_AF-Q6GGZ0-F1.pdb'),
            'ligand': os.path.join(base_dir, 'docking/vina_output/MurG_Afzelin_3D_out.pdbqt'),
            'resname': 'AFZ',
            'output': os.path.join(base_dir, 'data/complexes/complex_MurG_Afzelin.pdb')
        },
        {
            'protein': os.path.join(base_dir, 'data/receptors/MurG_AF-Q6GGZ0-F1.pdb'),
            'ligand': os.path.join(base_dir, 'docking/vina_output/MurG_Quercetin_3D_out.pdbqt'),
            'resname': 'QUE',
            'output': os.path.join(base_dir, 'data/complexes/complex_MurG_Quercetin.pdb')
        },
        {
            'protein': os.path.join(base_dir, 'data/receptors/PBP2a_3ZG0_clean.pdb'),
            'ligand': os.path.join(base_dir, 'docking/vina_output/PBP2a_Afzelin_3D_out.pdbqt'),
            'resname': 'AFZ',
            'output': os.path.join(base_dir, 'data/complexes/complex_PBP2a_Afzelin.pdb')
        },
        {
            'protein': os.path.join(base_dir, 'data/receptors/PBP2a_3ZG0_clean.pdb'),
            'ligand': os.path.join(base_dir, 'docking/vina_output/PBP2a_Ceftaroline_3D_out.pdbqt'),
            'resname': 'CEF',
            'output': os.path.join(base_dir, 'data/complexes/complex_PBP2a_Ceftaroline.pdb')
        }
    ]
    
    for job in jobs:
        build_complex(job['protein'], job['ligand'], job['resname'], job['output'])
        
    print("\nAll 4 complexes have been generated successfully!")

if __name__ == '__main__':
    main()
