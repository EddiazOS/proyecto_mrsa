#!/usr/bin/env python3
import os
import math

def main():
    cx, cy, cz = 4.729, 1.061, -3.361
    pdb_path = "data/receptors/MurG_AF-Q6GGZ0-F1.pdb"
    
    if not os.path.exists(pdb_path):
        print(f"Error: {pdb_path} not found.")
        return
        
    close_residues = {}
    
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    chain_id = line[21].strip()
                    
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    dist = math.sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
                    
                    if dist <= 8.0:
                        if res_num not in close_residues:
                            close_residues[res_num] = {
                                "name": res_name,
                                "min_dist": dist,
                                "atoms": 1
                            }
                        else:
                            close_residues[res_num]["atoms"] += 1
                            if dist < close_residues[res_num]["min_dist"]:
                                close_residues[res_num]["min_dist"] = dist
                except (ValueError, IndexError):
                    continue
                    
    print(f"Residues within 8.0 A of center ({cx}, {cy}, {cz}):")
    print(f"{'Residue':<10} {'Name':<6} {'Min Distance (A)':<18} {'Atoms in sphere'}")
    print("-" * 55)
    for res_num in sorted(close_residues.keys()):
        info = close_residues[res_num]
        print(f"{res_num:<10} {info['name']:<6} {info['min_dist']:<18.3f} {info['atoms']}")

if __name__ == "__main__":
    main()
