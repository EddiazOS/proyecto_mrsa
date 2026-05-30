#!/usr/bin/env python3
import os
import math

def main():
    cx, cy, cz = 4.729, 1.061, -3.361
    pdb_path = "data/receptors/MurG_AF-Q6GGZ0-F1.pdb"
    
    if not os.path.exists(pdb_path):
        print(f"Error: {pdb_path} not found.")
        return
        
    target_resnums = {14, 90, 125, 263, 284, 285, 287, 289}
    residues_data = {}
    
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    
                    if res_num in target_resnums:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        
                        dist = math.sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
                        
                        if res_num not in residues_data:
                            residues_data[res_num] = {
                                "name": res_name,
                                "min_dist": dist,
                                "coords": (x, y, z)
                            }
                        else:
                            if dist < residues_data[res_num]["min_dist"]:
                                residues_data[res_num]["min_dist"] = dist
                                residues_data[res_num]["coords"] = (x, y, z)
                except (ValueError, IndexError):
                    continue
                    
    print("Distance of literature-reported binding residues to our grid center (4.729, 1.061, -3.361):")
    print(f"{'Residue':<10} {'Name':<6} {'Min Distance (A)':<18} {'Coords (Closest Atom)'}")
    print("-" * 65)
    for num in sorted(target_resnums):
        if num in residues_data:
            info = residues_data[num]
            print(f"{num:<10} {info['name']:<6} {info['min_dist']:<18.3f} ({info['coords'][0]:.2f}, {info['coords'][1]:.2f}, {info['coords'][2]:.2f})")
        else:
            print(f"{num:<10} {'?':<6} {'Not found in PDB':<18}")

if __name__ == "__main__":
    main()
