#!/usr/bin/env python3
import os
import shutil
import subprocess

def setup():
    # El directorio base es el contenedor de 'scripts'
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print(f"=== Iniciando reorganización del proyecto en: {base_dir} ===")

    # 1. Definición de directorios a crear
    dirs_to_create = [
        "data/receptors",
        "data/ligands/raw_2d",
        "admet",
        "docking/vina_input",
        "docking/vina_output",
        "md/MurG",
        "md/PBP2a",
        "md/GyrB",
        "scripts"
    ]

    print("\n[1/4] Creando estructura de carpetas...")
    for d in dirs_to_create:
        d_path = os.path.join(base_dir, d)
        if not os.path.exists(d_path):
            os.makedirs(d_path)
            print(f"  + Creado: {d}")
        else:
            print(f"  . Ya existe: {d}")

    # 2. Definición de reubicaciones y renames
    # Llave: nombre original en la raíz, Valor: ruta destino relativa a base_dir
    file_mapping = {
        "3ZG0.pdb": "data/receptors/PBP2a_3ZG0.pdb",
        "3TTZ.pdb": "data/receptors/GyrB_3TTZ.pdb",
        "AF-Q6GGZ0-F1-model_v6.pdb": "data/receptors/MurG_AF-Q6GGZ0-F1.pdb",
        "4Hydroxybenzoic_acid_CID135_3D.sdf": "data/ligands/4Hydroxybenzoic_acid_3D.sdf",
        "Afzelin_CID5316673_3D.sdf": "data/ligands/Afzelin_3D.sdf",
        "Avocadenofuran_CID6857792_3D.sdf": "data/ligands/Avocadenofuran_3D.sdf",
        "Avocadyne_acetate_3D.sdf": "data/ligands/Avocadyne_acetate_3D.sdf",
        "Avocadene_acetate_3D.sdf": "data/ligands/Avocadene_acetate_3D.sdf",
        "Avocadyne_acetate_CID3952079_2D_ONLY.sdf": "data/ligands/raw_2d/Avocadyne_acetate_2D.sdf",
        "Avocadene_acetate_CID3624980_2D_ONLY.sdf": "data/ligands/raw_2d/Avocadene_acetate_2D.sdf"
    }

    print("\n[2/4] Moviendo y renombrando archivos estructurales...")
    for src_name, dst_rel in file_mapping.items():
        src_path = os.path.join(base_dir, src_name)
        dst_path = os.path.join(base_dir, dst_rel)

        if os.path.exists(src_path):
            # Si el archivo destino ya existe, lo sobrescribimos para asegurar consistencia
            shutil.move(src_path, dst_path)
            print(f"  -> Movido: {src_name} a {dst_rel}")
        else:
            # Si ya se había movido antes
            if os.path.exists(dst_path):
                print(f"  . Ya posicionado: {dst_rel}")
            else:
                print(f"  [ALERTA] Archivo no encontrado en la raíz: {src_name}")

    # 3. Eliminar archivos NTFS Zone.Identifier
    print("\n[3/4] Eliminando archivos temporales NTFS (*:Zone.Identifier)...")
    cleaned_count = 0
    # Buscar en todo el directorio base recursivamente
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if "Zone.Identifier" in file:
                file_path = os.path.join(root, file)
                try:
                    os.remove(file_path)
                    print(f"  x Eliminado: {os.path.relpath(file_path, base_dir)}")
                    cleaned_count += 1
                except Exception as e:
                    print(f"  [ERROR] No se pudo eliminar {file}: {e}")
    print(f"  Total de archivos de metadatos eliminados: {cleaned_count}")

    # 4. Ejecutar la extracción de controles positivos y purificación de receptores
    print("\n[4/4] Extrayendo controles positivos de ligandos y purificando receptores...")
    extract_script = os.path.join(base_dir, "scripts", "extract_controls.py")
    if os.path.exists(extract_script):
        try:
            # Ejecutamos el script de extracción de controles usando python3
            result = subprocess.run(["python3", extract_script], capture_output=True, text=True, check=True)
            print(result.stdout)
            if result.stderr:
                print(f"Errores en la extracción:\n{result.stderr}")
        except subprocess.CalledProcessError as e:
            print(f"  [ERROR] Falló la ejecución de extract_controls.py: {e}")
            print(f"Detalles:\n{e.output}\n{e.stderr}")
    else:
        print("  [ERROR] No se encontró el script extract_controls.py en scripts/")

    print("\n=== ¡Organización inicial completada con éxito! ===")

if __name__ == "__main__":
    setup()
