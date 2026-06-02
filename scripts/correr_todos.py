#!/usr/bin/env python3
"""
correr_todos.py
===============
Orquestador del pipeline de dinámica molecular para los 4 complejos MRSA.
Corre secuencialmente cada complejo y respalda los resultados en el
Network Volume de RunPod (/workspace/backups/).

Uso:
  python scripts/correr_todos.py [TIEMPO_NS]

Ejemplos:
  python scripts/correr_todos.py          # corre 100 ns por complejo (default)
  python scripts/correr_todos.py 10       # corre 10 ns por complejo (prueba rápida)
"""

import sys
import os
import shutil
import subprocess


def respaldar_resultados(local_folder, backup_base):
    """
    Copia la carpeta de resultados al directorio de backup en el
    Network Volume de RunPod (/workspace/backups/).
    Si ya existe un backup previo del mismo complejo, lo sobreescribe.
    """
    if not os.path.exists(local_folder):
        print(f"[Backup WARNING] La carpeta '{local_folder}' no existe. Omitiendo backup.")
        return False

    try:
        os.makedirs(backup_base, exist_ok=True)
        dest_path = os.path.join(backup_base, local_folder)

        if os.path.exists(dest_path):
            shutil.rmtree(dest_path)

        shutil.copytree(local_folder, dest_path)
        print(f"[Backup OK] Resultados guardados en: {dest_path}")
        return True
    except Exception as e:
        print(f"[Backup WARNING] No se pudo realizar el backup: {e}")
        return False


def main():
    # Tiempo de producción por defecto: 100 ns
    tiempo_ns = 100.0
    if len(sys.argv) >= 2:
        try:
            tiempo_ns = float(sys.argv[1])
        except ValueError:
            print("[WARNING] Tiempo de simulación inválido, usando 100.0 ns por defecto.")

    complex_list = [
        "MurG_Afzelin",
        "MurG_Quercetin",
        "PBP2a_Afzelin",
        "PBP2a_Ceftaroline"
    ]

    # Directorio de backup en el Network Volume de RunPod
    # Si no existe /workspace (local o Colab), usa el directorio actual
    if os.path.exists("/workspace"):
        backup_base = "/workspace/backups"
    else:
        backup_base = os.path.join(os.getcwd(), "backups")

    print("=" * 70)
    print("      INICIANDO PIPELINE DE SIMULACIÓN MULTI-COMPLEJO (GPU)")
    print(f"      Tiempo de producción por complejo : {tiempo_ns} ns")
    print(f"      Total de complejos a simular      : {len(complex_list)}")
    print(f"      Directorio de backup               : {backup_base}")
    print("=" * 70)

    resultados = {}

    for i, name in enumerate(complex_list, 1):
        print(f"\n[{i}/{len(complex_list)}] >>> INICIANDO SIMULACIÓN DE: {name} ({tiempo_ns} ns) <<<")

        cmd = f"python scripts/ejecutar_dinamica.py {name} {tiempo_ns}"
        res = subprocess.run(cmd, shell=True)

        if res.returncode == 0:
            print(f"\n[OK] Simulación de {name} finalizada exitosamente.")
            exito_backup = respaldar_resultados(f"md_run_{name}", backup_base)
            resultados[name] = "✓ Completado" + (" + Backup OK" if exito_backup else " (sin backup)")
        else:
            print(f"\n[ERROR] La simulación de {name} falló o fue interrumpida.")
            print("Continuando con el siguiente complejo...")
            resultados[name] = "✗ Falló"

    # Resumen final
    print("\n" + "=" * 70)
    print("      ¡PIPELINE MULTI-COMPLEJO FINALIZADO!")
    print("=" * 70)
    print("\nResumen de resultados:")
    for name, estado in resultados.items():
        print(f"  {name:<25} → {estado}")
    print(f"\nBackups disponibles en: {backup_base}")
    print("=" * 70)


if __name__ == "__main__":
    main()
