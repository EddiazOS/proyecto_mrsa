#!/usr/bin/env python3
import os
import numpy as np

def find_coordinates(path_file):
    """
    Lee un archivo PDB y extrae el centro geométrico (promedio de coordenadas)
    de todas las líneas ATOM y HETATM.
    """
    if not os.path.exists(path_file):
        print(f"Error: El archivo {path_file} no existe.")
        return None
    
    atoms_coordinates = []

    with open(path_file, "r") as file:
        for line in file:
            # Evaluamos si la línea representa coordenadas atómicas
            if line.startswith("HETATM") or line.startswith("ATOM"):
                # Rango exacto del estándar PDB (columnas 31-38, 39-46, 47-54)
                x = line[30:38].strip()
                y = line[38:46].strip()
                z = line[46:54].strip()

                try:
                    atoms_coordinates.append([float(x), float(y), float(z)])
                except ValueError:
                    # En caso de filas mal formateadas o cabeceras
                    continue
    
    if not atoms_coordinates:
        print(f"Advertencia: No se encontraron coordenadas válidas en {path_file}.")
        return None
    
    # Convertimos a matriz NumPy una sola vez para eficiencia
    coordinates_matrix = np.array(atoms_coordinates)
    
    # Calculamos el promedio sobre el eje 0 (columnas) para obtener el centro geométrico [X, Y, Z]
    return coordinates_matrix.mean(axis=0)

def generar_config_vina(receptor_nombre, centro, tamaño_caja, archivo_salida):
    """
    Genera automáticamente el archivo de configuración conf.txt para AutoDock Vina
    empleando los parámetros del protocolo de investigación.
    """
    # La ruta de salida será docking/vina_input/ dentro de la raíz del proyecto
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    carpeta_salida = os.path.join(base_dir, "docking", "vina_input")
    os.makedirs(carpeta_salida, exist_ok=True)

    ruta_completa = os.path.join(carpeta_salida, archivo_salida)

    if centro is None:
        print(f"Error: No se puede escribir la configuración de {receptor_nombre} porque el centro es nulo.")
        return

    try:
        with open(ruta_completa, "w") as f:
            f.write(f"# Configuración de AutoDock Vina para {receptor_nombre}\n")
            f.write("# Generado automáticamente por Andrés Pérez Palacios\n\n")

            # Rutas de entrada relativas (asumiendo que Vina se corre desde la raíz)
            f.write(f"receptor       = data/receptors/{receptor_nombre}_clean.pdbqt\n")
            f.write(f"ligand         = data/ligands/cambiar_por_ligando.pdbqt\n\n")

            # Centro de la Rejilla (Grid Box Center)
            f.write(f"center_x       = {centro[0]:.3f}\n")
            f.write(f"center_y       = {centro[1]:.3f}\n")
            f.write(f"center_z       = {centro[2]:.3f}\n\n")

            # Tamaño de la Rejilla (Grid Box Size)
            f.write(f"size_x         = {tamaño_caja:.1f}\n")
            f.write(f"size_y         = {tamaño_caja:.1f}\n")
            f.write(f"size_z         = {tamaño_caja:.1f}\n\n")

            # Parámetros metodológicos fijos para publicación
            f.write("exhaustiveness = 32\n")
            f.write("num_modes      = 9\n")
            f.write("energy_range   = 3\n")

        print(f"¡Configuración guardada con éxito en: {os.path.relpath(ruta_completa, base_dir)}!")

    except Exception as e:
        print(f"Error al escribir el archivo de configuración: {e}")

if __name__ == "__main__":
    # Definición de rutas relativas basadas en la estructura del proyecto
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    
    ligand_pbp2a = os.path.join(base_dir, "data", "ligands", "Ceftaroline_3D.pdb")
    ligand_gyrb = os.path.join(base_dir, "data", "ligands", "07N_GyrB_3D.pdb")

    # 1. Procesar Centro e Iniciar configuraciones para PBP2a
    print("Procesando datos para PBP2a...")
    centro_pbp2a = find_coordinates(ligand_pbp2a)
    if centro_pbp2a is not None:
        print(f"  Centro calculado -> X: {centro_pbp2a[0]:.3f}, Y: {centro_pbp2a[1]:.3f}, Z: {centro_pbp2a[2]:.3f}")
        generar_config_vina("PBP2a_3ZG0", centro_pbp2a, 22.0, "conf_pbp2a.txt")

    # 2. Procesar Centro e Iniciar configuraciones para GyrB
    print("\nProcesando datos para GyrB...")
    centro_gyrb = find_coordinates(ligand_gyrb)
    if centro_gyrb is not None:
        print(f"  Centro calculado -> X: {centro_gyrb[0]:.3f}, Y: {centro_gyrb[1]:.3f}, Z: {centro_gyrb[2]:.3f}")
        generar_config_vina("GyrB_3TTZ", centro_gyrb, 22.0, "conf_gyrb.txt")
