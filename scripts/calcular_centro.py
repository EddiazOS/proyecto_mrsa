#!/usr/bin/env python3
"""
Script para calcular el centro geométrico (centroide) de un ligando
a partir de su archivo PDB, y determinar el tamaño mínimo de caja
para AutoDock Vina.

Correcciones de auditoría aplicadas:
  - Filtro ALTLOC: solo conformación A o vacía (descarta B, C, etc.)
  - Filtro de solvente/iones: excluye HOH, SO4, EDO, GOL, PEG, etc.
  - Cálculo dinámico del tamaño de caja basado en las dimensiones del ligando

Uso:
    python scripts/calcular_centro.py                    # procesa PBP2a y GyrB
    python scripts/calcular_centro.py <archivo.pdb>      # procesa un ligando individual
"""

import os
import sys
import numpy as np


# Residuos de solvente, tampón e iones que NO son parte del ligando de interés.
# Se excluyen del cálculo del centroide para evitar desplazamientos espurios.
EXCLUIR_RESIDUOS = {
    'HOH', 'WAT',          # agua
    'SO4', 'PO4',          # tampones inorgánicos
    'EDO', 'GOL', 'PEG',   # crioprotectores / polietilenglicol
    'NA', 'CL', 'MG', 'ZN', 'CA', 'MN', 'FE', 'K',  # iones comunes
    'ACE', 'NH2',          # caps de modelado
}


def _extraer_coords_limpias(path_file, chain_id=None):
    """
    Lee un archivo PDB y retorna un array numpy (N×3) con las coordenadas
    de átomos ATOM/HETATM, aplicando:
      1. Filtro de cadena: si chain_id no es None, solo retiene esa cadena.
      2. Filtro ALTLOC: solo conformación ' ', 'A' o '' (columna 17 del PDB).
      3. Filtro de residuos: excluye solvente, iones y tampones.

    Args:
        path_file: ruta al archivo PDB
        chain_id:  cadena a retener ('A', 'B', etc.) o None para todas.
                   Para ligandos co-cristalizados extraídos de homodímeros,
                   usar 'A' para evitar promediar copias de la unidad asimétrica.

    Retorna None si no se encuentran coordenadas válidas.
    """
    if not os.path.exists(path_file):
        print(f"Error: El archivo {path_file} no existe.")
        return None

    coords = []

    with open(path_file, "r") as f:
        for line in f:
            if not (line.startswith("HETATM") or line.startswith("ATOM")):
                continue

            # --- Filtro de cadena (columna 22, índice 21) ---
            if chain_id is not None:
                pdb_chain = line[21] if len(line) > 21 else ' '
                if pdb_chain != chain_id:
                    continue

            # --- Filtro ALTLOC (columna 17, índice 16) ---
            altloc = line[16] if len(line) > 16 else ' '
            if altloc not in (' ', 'A', ''):
                continue

            # --- Filtro de residuos excluidos ---
            res_name = line[17:20].strip()
            if res_name in EXCLUIR_RESIDUOS:
                continue

            # --- Extracción de coordenadas (columnas 31-54 del estándar PDB) ---
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue

    if not coords:
        print(f"Advertencia: No se encontraron coordenadas válidas en {path_file}.")
        return None

    return np.array(coords)


def find_coordinates(path_file, chain_id=None):
    """
    Lee un archivo PDB y retorna el centro geométrico (centroide)
    como un array numpy [X, Y, Z].

    Aplica filtros de cadena, ALTLOC y solvente/iones antes del cálculo.
    NOTA: Esto es un centroide geométrico (promedio simple), NO un
    centro de masa ponderado por masa atómica. Para AutoDock Vina,
    el centroide geométrico es el método estándar.

    Args:
        path_file: ruta al archivo PDB
        chain_id:  cadena a retener ('A', etc.) o None para todas
    """
    coords = _extraer_coords_limpias(path_file, chain_id=chain_id)
    if coords is None:
        return None

    centroide = coords.mean(axis=0)
    return centroide


def calcular_tamaño_caja(path_file, margen=10.0, minimo=22.0, chain_id=None):
    """
    Calcula el tamaño mínimo de caja cúbica para contener el ligando
    en cualquier orientación, con el margen especificado.

    La fórmula es: tamaño = dimensión_máxima_del_ligando + 2 × margen
    Se garantiza un mínimo absoluto.

    Args:
        path_file: ruta al PDB del ligando
        margen:    espacio adicional por lado en Å (default: 10.0)
        minimo:    tamaño mínimo absoluto en Å (default: 22.0)
        chain_id:  cadena a retener ('A', etc.) o None para todas
    Returns:
        float: tamaño de caja en Å (misma dimensión para X, Y, Z)
    """
    coords = _extraer_coords_limpias(path_file, chain_id=chain_id)
    if coords is None:
        print(f"  Usando tamaño de caja mínimo ({minimo} Å) por falta de coordenadas.")
        return minimo

    dims = coords.max(axis=0) - coords.min(axis=0)
    tamaño_necesario = dims.max() + 2 * margen
    tamaño_final = round(max(tamaño_necesario, minimo), 1)

    return tamaño_final


def calcular_tamaño_caja_sdf(path_file, margen=10.0, minimo=22.0):
    """
    Versión para archivos SDF (V2000 MOL format).
    Extrae coordenadas de la sección de átomos del SDF.

    Args:
        path_file: ruta al archivo SDF
        margen:    espacio adicional por lado en Å (default: 10.0)
        minimo:    tamaño mínimo absoluto en Å (default: 22.0)
    Returns:
        float: tamaño de caja en Å
    """
    if not os.path.exists(path_file):
        print(f"Error: El archivo {path_file} no existe.")
        return minimo

    coords = []

    with open(path_file, "r") as f:
        lines = f.readlines()

    # Buscar la línea de conteo (línea 4 en formato V2000, índice 3)
    # Formato: "aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv"
    # aaa = número de átomos, bbb = número de enlaces
    n_atoms = None
    counts_line_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.endswith('V2000') or stripped.endswith('V3000'):
            try:
                n_atoms = int(stripped[:3])
                counts_line_idx = i
                break
            except (ValueError, IndexError):
                continue

    if n_atoms is None or counts_line_idx is None:
        print(f"Advertencia: No se pudo parsear el SDF {path_file}.")
        return minimo

    # Las líneas de átomos empiezan después de la línea de conteo
    for i in range(counts_line_idx + 1, counts_line_idx + 1 + n_atoms):
        if i >= len(lines):
            break
        line = lines[i]
        try:
            x = float(line[0:10].strip())
            y = float(line[10:20].strip())
            z = float(line[20:30].strip())
            coords.append([x, y, z])
        except (ValueError, IndexError):
            continue

    if not coords:
        print(f"Advertencia: No se encontraron coordenadas en {path_file}.")
        return minimo

    coords = np.array(coords)
    dims = coords.max(axis=0) - coords.min(axis=0)
    tamaño_necesario = dims.max() + 2 * margen
    return round(max(tamaño_necesario, minimo), 1)


def generar_config_vina(receptor_nombre, centro, tamaño_caja, archivo_salida,
                        comentario_centro=""):
    """
    Genera automáticamente el archivo de configuración conf.txt para AutoDock Vina
    empleando los parámetros del protocolo de investigación.

    Args:
        receptor_nombre: nombre base del receptor (sin extensión)
        centro:          array [X, Y, Z] del centroide
        tamaño_caja:     tamaño de la caja en Å (mismo para X, Y, Z)
        archivo_salida:  nombre del archivo de salida (ej. "conf_pbp2a.txt")
        comentario_centro: texto descriptivo sobre el origen de las coordenadas
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    carpeta_salida = os.path.join(base_dir, "docking", "vina_input")
    os.makedirs(carpeta_salida, exist_ok=True)

    ruta_completa = os.path.join(carpeta_salida, archivo_salida)

    if centro is None:
        print(f"Error: No se puede escribir la configuración de {receptor_nombre} "
              f"porque el centro es nulo.")
        return

    try:
        with open(ruta_completa, "w") as f:
            f.write(f"# Configuración de AutoDock Vina para {receptor_nombre}\n")
            f.write("# Generado automáticamente por Andrés Pérez Palacios\n")
            if comentario_centro:
                f.write(f"# Centro: {comentario_centro}\n")
            f.write(f"# Tamaño de caja: {tamaño_caja:.1f} Å "
                    f"(calculado dinámicamente, margen 10 Å)\n\n")

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

        print(f"  ¡Configuración guardada en: "
              f"{os.path.relpath(ruta_completa, base_dir)}!")

    except Exception as e:
        print(f"Error al escribir el archivo de configuración: {e}")


def _calcular_caja_max_ligandos(base_dir, archivos_ligandos, margen=10.0, minimo=22.0):
    """
    Calcula el tamaño máximo de caja necesario entre un conjunto de ligandos.
    Soporta archivos PDB y SDF.
    Para archivos PDB, usa chain_id='A' para evitar promediar copias
    de la unidad asimétrica en homodímeros.

    Returns:
        float: tamaño máximo de caja para el conjunto
    """
    tamaños = []
    for archivo in archivos_ligandos:
        ruta = os.path.join(base_dir, "data", "ligands", archivo)
        if not os.path.exists(ruta):
            continue

        if archivo.endswith('.sdf'):
            t = calcular_tamaño_caja_sdf(ruta, margen=margen, minimo=minimo)
        else:
            t = calcular_tamaño_caja(ruta, margen=margen, minimo=minimo,
                                     chain_id='A')
        tamaños.append(t)
        print(f"    {archivo}: {t:.1f} Å")

    if not tamaños:
        return minimo
    return max(tamaños)


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)

    # --- Modo individual: calcular centro de un solo archivo ---
    if len(sys.argv) >= 2:
        path_file = sys.argv[1]
        print(f"\n{'='*60}")
        print(f"Calculando centro geométrico para: {path_file}")
        print(f"{'='*60}")

        centro = find_coordinates(path_file)
        if centro is not None:
            print(f"\nCentro geométrico del ligando:")
            print(f"  center_x = {centro[0]:.3f}")
            print(f"  center_y = {centro[1]:.3f}")
            print(f"  center_z = {centro[2]:.3f}")

            tamaño = calcular_tamaño_caja(path_file)
            print(f"\nTamaño de caja sugerido: {tamaño:.1f} Å")
            print(f"  (dimensión_máx_ligando + 2 × 10.0 Å margen)")
        sys.exit(0)

    # --- Modo completo: procesar PBP2a y GyrB ---
    print("=" * 60)
    print("CÁLCULO DE CENTROS Y TAMAÑOS DE CAJA")
    print("Auditoría: filtros ALTLOC y solvente aplicados")
    print("=" * 60)

    # Ligandos co-cristalizados (para centros)
    ligand_pbp2a = os.path.join(base_dir, "data", "ligands", "Ceftaroline_3D.pdb")
    ligand_gyrb = os.path.join(base_dir, "data", "ligands", "07N_GyrB_3D.pdb")

    # Todos los ligandos del estudio (para tamaño de caja)
    todos_ligandos = [
        "Avocadenofuran_3D.sdf",
        "Avocadyne_acetate_3D.sdf",
        "Avocadene_acetate_3D.sdf",
        "Ceftaroline_3D.pdb",
        "07N_GyrB_3D.pdb",
        "Quercetin_3D.sdf",
    ]

    # ---- Calcular tamaño máximo de caja entre todos los ligandos ----
    print("\n--- Tamaño de caja por ligando (margen = 10 Å) ---")
    tamaño_max = _calcular_caja_max_ligandos(base_dir, todos_ligandos)
    print(f"\n  → Tamaño máximo de caja (para todas las dianas): {tamaño_max:.1f} Å")

    # ---- PBP2a ----
    print(f"\n{'='*60}")
    print("1. PBP2a (PDB: 3ZG0)")
    print(f"{'='*60}")
    centro_pbp2a = find_coordinates(ligand_pbp2a, chain_id='A')
    if centro_pbp2a is not None:
        print(f"  Centro geométrico (Ceftarolina AI8):")
        print(f"    center_x = {centro_pbp2a[0]:.3f}")
        print(f"    center_y = {centro_pbp2a[1]:.3f}")
        print(f"    center_z = {centro_pbp2a[2]:.3f}")
        print(f"  Tamaño de caja: {tamaño_max:.1f} Å")
        generar_config_vina(
            "PBP2a_3ZG0", centro_pbp2a, tamaño_max, "conf_pbp2a.txt",
            comentario_centro="centroide geométrico del ligando co-cristalizado (Ceftarolina AI8)"
        )

    # ---- GyrB ----
    print(f"\n{'='*60}")
    print("2. GyrB (PDB: 3TTZ)")
    print(f"{'='*60}")
    centro_gyrb = find_coordinates(ligand_gyrb, chain_id='A')
    if centro_gyrb is not None:
        print(f"  Centro geométrico (07N):")
        print(f"    center_x = {centro_gyrb[0]:.3f}")
        print(f"    center_y = {centro_gyrb[1]:.3f}")
        print(f"    center_z = {centro_gyrb[2]:.3f}")
        print(f"  Tamaño de caja: {tamaño_max:.1f} Å")
        generar_config_vina(
            "GyrB_3TTZ", centro_gyrb, tamaño_max, "conf_gyrb.txt",
            comentario_centro="centroide geométrico del ligando co-cristalizado (07N)"
        )

    # ---- MurG (no se procesa aquí, requiere validar_sitio_murg.py) ----
    print(f"\n{'='*60}")
    print("3. MurG (AlphaFold: AF-Q6GGZ0)")
    print(f"{'='*60}")
    print("  ⚠ MurG no tiene ligando co-cristalizado.")
    print("  → Ejecutar scripts/validar_sitio_murg.py para obtener")
    print("    las coordenadas del sitio activo por superposición con 1F0K.")
    print(f"  → Usar tamaño de caja: {tamaño_max:.1f} Å")
