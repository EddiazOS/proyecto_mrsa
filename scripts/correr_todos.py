#!/usr/bin/env python3
import sys
import subprocess

def main():
    tiempo_ns = 10.0
    if len(sys.argv) >= 2:
        try:
            tiempo_ns = float(sys.argv[1])
        except ValueError:
            print("[WARNING] Tiempo de simulación inválido, usando 10.0 ns por defecto.")

    complex_list = [
        "MurG_Afzelin",
        "MurG_Quercetin",
        "PBP2a_Afzelin",
        "PBP2a_Ceftaroline"
    ]

    print("=" * 70)
    print("      INICIANDO PIPELINE DE SIMULACIÓN MULTI-COMPLEJO")
    print(f"      Tiempo de producción por complejo: {tiempo_ns} ns")
    print(f"      Total de complejos a simular     : {len(complex_list)}")
    print("=" * 70)

    for i, name in enumerate(complex_list, 1):
        print(f"\n[{i}/{len(complex_list)}] >>> INICIANDO SIMULACIÓN DE: {name} ({tiempo_ns} ns) <<<")
        # Se ejecuta secuencialmente la dinámica para cada complejo
        cmd = f"python scripts/ejecutar_dinamica.py {name} {tiempo_ns}"
        res = subprocess.run(cmd, shell=True)
        
        if res.returncode == 0:
            print(f"\n[OK] Simulación de {name} finalizada exitosamente.")
        else:
            print(f"\n[ERROR] La simulación de {name} falló o fue interrumpida.")
            print("Continuando con el siguiente complejo de la lista...")

    print("\n" + "=" * 70)
    print("      ¡PIPELINE MULTI-COMPLEJO FINALIZADO COMPLETAMENTE!")
    print("=" * 70)

if __name__ == "__main__":
    main()
