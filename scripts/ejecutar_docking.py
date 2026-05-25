#!/usr/bin/env python3
"""
ejecutar_docking.py
===================
Ejecuta el screening virtual completo para el proyecto de MRSA:
  3 receptores (PBP2a, GyrB, MurG) x (5 compuestos HPLC + 1 control c/u) = 18 experimentos.

El script:
  1. Identifica los archivos PDBQT de receptores y ligandos.
  2. Ejecuta AutoDock Vina 1.2 usando las cajas optimizadas en conf_*.txt.
  3. Organiza los resultados en docking/vina_output/.
  4. Parsea los archivos de registro (.log) para extraer las mejores afinidades (kcal/mol).
  5. Genera una tabla de resumen final y un reporte en Markdown.

Uso:
  conda run -n docking_env python scripts/ejecutar_docking.py
"""

import os
import sys
import subprocess
import time
import re

# =====================================================================
# Configuración del Screening
# =====================================================================
HPLC_LIGANDS = [
    "Afzelin_3D",
    "4Hydroxybenzoic_acid_3D",
    "Avocadenofuran_3D",
    "Avocadyne_acetate_3D",
    "Avocadene_acetate_3D"
]

DIANAS = {
    "PBP2a": {
        "config": "docking/vina_input/conf_pbp2a.txt",
        "control": "Ceftaroline_3D",
        "receptor": "data/receptors/PBP2a_3ZG0_clean.pdbqt",
        "desc": "Proteína de unión a penicilina 2a (3ZG0)"
    },
    "GyrB": {
        "config": "docking/vina_input/conf_gyrb.txt",
        "control": "07N_GyrB_3D",
        "receptor": "data/receptors/GyrB_3TTZ_clean.pdbqt",
        "desc": "DNA girasa subunidad B (3TTZ)"
    },
    "MurG": {
        "config": "docking/vina_input/conf_murg.txt",
        "control": "Quercetin_3D",
        "receptor": "data/receptors/MurG_AF-Q6GGZ0-F1_clean.pdbqt",
        "desc": "Glucosiltransferasa MurG (AlphaFold AF-Q6GGZ0-F1)"
    }
}

def parsear_energia_vina(log_path):
    """Parsea el log de Vina y extrae la afinidad de la mejor pose (modo 1)."""
    if not os.path.exists(log_path):
        return None
    
    with open(log_path, "r") as f:
        content = f.read()
        
    # Buscar la primera línea de la tabla de resultados
    # Ejemplo:
    #    1      -8.4      0.000      0.000
    match = re.search(r"\s+1\s+(-?\d+\.\d+)\s+", content)
    if match:
        return float(match.group(1))
    return None

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, os.pardir))
    
    output_dir = os.path.join(project_root, "docking", "vina_output")
    os.makedirs(output_dir, exist_ok=True)
    
    print("====================================================================")
    print("  SCREENING VIRTUAL DE DOCKING MOLECULAR (AutoDock Vina)")
    print("  Proyecto: Compuestos de Persea americana vs. MRSA")
    print("====================================================================")
    
    # 1. Validaciones previas de archivos
    for diana_name, info in DIANAS.items():
        conf_path = os.path.join(project_root, info["config"])
        rec_path = os.path.join(project_root, info["receptor"])
        
        if not os.path.exists(conf_path):
            print(f"[ERROR] No se encontró el archivo de configuración: {conf_path}")
            sys.exit(1)
        if not os.path.exists(rec_path):
            print(f"[ERROR] No se encontró el receptor PDBQT: {rec_path}")
            sys.exit(1)
            
    # Verificar existencia de Vina
    try:
        subprocess.run(["vina", "--version"], capture_output=True, text=True)
    except FileNotFoundError:
        print("[ERROR] AutoDock Vina ('vina') no se encuentra en el PATH actual.")
        print("Asegúrate de ejecutar este script dentro del entorno conda 'docking_env'.")
        sys.exit(1)

    # 2. Bucle de docking (18 corridas)
    resultados = []
    total_runs = len(DIANAS) * (len(HPLC_LIGANDS) + 1)
    run_idx = 1
    
    t_inicio_global = time.time()
    
    for diana_name, info in DIANAS.items():
        conf_path = os.path.join(project_root, info["config"])
        
        # Preparar lista de ligandos (Compuestos HPLC + su Control respectivo)
        ligandos_para_correr = HPLC_LIGANDS + [info["control"]]
        
        print(f"\n[DIANA] {diana_name} - {info['desc']}")
        print("-" * 68)
        
        for lig_name in ligandos_para_correr:
            lig_path = os.path.join(project_root, "data", "ligands", f"{lig_name}.pdbqt")
            
            # Verificar existencia del ligando pdbqt
            if not os.path.exists(lig_path):
                print(f"  [ERROR] No se encontró el ligando PDBQT: {lig_path}")
                continue
                
            out_pdbqt = os.path.join(output_dir, f"{diana_name}_{lig_name}_out.pdbqt")
            log_txt = os.path.join(output_dir, f"{diana_name}_{lig_name}_log.txt")
            
            is_control = "Sí" if lig_name == info["control"] else "No"
            
            print(f"  [{run_idx}/{total_runs}] Docking: {lig_name} (Control: {is_control})...")

            if os.path.exists(log_txt):
                run_idx += 1
                print(f"{lig_name} ya ha sido procesado y sus resultados están en {log_txt}|{out_pdbqt}")
                continue
            
            # Comando de Vina (sin parámetro --log descontinuado)
            cmd = [
                "vina",
                "--config", conf_path,
                "--ligand", lig_path,
                "--out", out_pdbqt
            ]

            t0 = time.time()
            res = subprocess.run(cmd, capture_output=True, text=True)
            elapsed = time.time() - t0

            # Guardamos la salida estándar (stdout) directamente en el archivo log_txt
            with open(log_txt, "w") as log_file:
                log_file.write(res.stdout)
                # Opcional: si hubo errores de consola, también los guardamos al final
                if res.stderr:
                    log_file.write("\n=== ERRORES (STDERR) ===\n")
                    log_file.write(res.stderr)

            if res.returncode == 0:
                # Ahora parseamos el archivo recién escrito. ¡Funcionará idénticamente!
                afinidad = parsear_energia_vina(log_txt)
                print(f"    -> Completado en {elapsed:.1f}s | Mejor Afinidad (ΔG): {afinidad} kcal/mol")
                resultados.append({
                    "Diana": diana_name,
                    "Ligando": lig_name,
                    "Control": is_control,
                    "Energia_kcal_mol": afinidad,
                    "Tiempo_s": round(elapsed, 1)
                })
                
            run_idx += 1
            
    t_fin_global = time.time()
    duracion_total = t_fin_global - t_inicio_global
    
    # 3. Guardar resultados y generar reporte
    generar_reporte(resultados, output_dir, duracion_total)

def generar_reporte(resultados, output_dir, duracion_total):
    """Guarda los resultados del docking en CSV y genera un reporte legible."""
    csv_path = os.path.join(output_dir, "resultados_docking.csv")
    md_path = os.path.join(output_dir, "reporte_docking.md")
    
    # Escribir CSV
    with open(csv_path, "w") as f:
        f.write("Diana,Ligando,Control,Energia_kcal_mol,Tiempo_s\n")
        for r in resultados:
            f.write(f"{r['Diana']},{r['Ligando']},{r['Control']},{r['Energia_kcal_mol']},{r['Tiempo_s']}\n")
            
    # Escribir Reporte Markdown
    with open(md_path, "w") as f:
        f.write("# Reporte del Screening Virtual de Docking Molecular\n\n")
        f.write(f"**Fecha y Hora:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"**Duración total del screening:** {duracion_total/60:.2f} minutos\n")
        f.write("**Algoritmo:** AutoDock Vina 1.2.0\n")
        f.write("**Exhaustiveness:** 32 (Exhaustividad alta para publicación)\n\n")
        f.write("Este reporte consolida los resultados del docking molecular de los 5 compuestos de *Persea americana* y los controles positivos correspondientes frente a las tres dianas de *S. aureus* MRSA.\n\n")
        
        # Agrupar resultados por diana
        dianas = sorted(list(set(r["Diana"] for r in resultados)))
        
        for d in dianas:
            f.write(f"## Resultados para {d}\n\n")
            f.write("| Posición | Ligando | ¿Control? | Energía de Unión (ΔG, kcal/mol) | Tiempo (s) |\n")
            f.write("| :---: | :--- | :---: | :---: | :---: |\n")
            
            d_res = [r for r in resultados if r["Diana"] == d]
            # Ordenar de mayor a menor afinidad (energías más negativas primero)
            # Manejar posibles valores None
            d_res_sorted = sorted(
                d_res,
                key=lambda x: x["Energia_kcal_mol"] if x["Energia_kcal_mol"] is not None else 999.0
            )
            
            for idx, r in enumerate(d_res_sorted):
                pos = idx + 1
                energia = f"{r['Energia_kcal_mol']:.1f}" if r["Energia_kcal_mol"] is not None else "Error"
                is_ctrl = "Sí" if r["Control"] == "Sí" else "No"
                
                # Resaltar si supera o iguala al control positivo
                ctrl_energy = next((x["Energia_kcal_mol"] for x in d_res if x["Control"] == "Sí"), None)
                if r["Control"] == "No" and ctrl_energy is not None and r["Energia_kcal_mol"] is not None:
                    if r["Energia_kcal_mol"] <= ctrl_energy:
                        energia = f"**{energia}** 🌟"
                        
                f.write(f"| {pos} | {r['Ligando']} | {is_ctrl} | {energia} | {r['Tiempo_s']} |\n")
                
            f.write("\n*Nota: Las energías más negativas indican mayor afinidad teórica. Las celdas en negrita y marcadas con 🌟 representan compuestos de interés que igualan o superan la afinidad del control positivo de esa diana.*\n\n")
            
    print("\n====================================================================")
    print("  SCREENING COMPLETADO CON ÉXITO")
    print(f"  Resultados guardados en CSV: {csv_path}")
    print(f"  Reporte generado en Markdown: {md_path}")
    print("====================================================================")

if __name__ == "__main__":
    main()
