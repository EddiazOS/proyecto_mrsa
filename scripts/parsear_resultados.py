#!/usr/bin/env python3
"""
parsear_resultados.py
=====================
Script independiente para escanear y parsear todos los archivos log de docking
(*_log.txt) en la carpeta docking/vina_output/, extraer las mejores afinidades,
recuperar los tiempos de ejecución previos del archivo CSV (si existe), y generar
los archivos consolidados:
  - resultados_docking.csv
  - reporte_docking.md

Uso:
  python scripts/parsear_resultados.py
"""

import os
import re
import sys
import time

# =====================================================================
# Configuración y Constantes
# =====================================================================
DIANAS = {
    "PBP2a": {
        "control": "Ceftaroline_3D",
        "desc": "Proteína de unión a penicilina 2a (3ZG0)"
    },
    "GyrB": {
        "control": "07N_GyrB_3D",
        "desc": "DNA girasa subunidad B (3TTZ)"
    },
    "MurG": {
        "control": "Quercetin_3D",
        "desc": "Glucosiltransferasa MurG (AlphaFold AF-Q6GGZ0-F1)"
    }
}

def parsear_energia_vina(log_path):
    """Parsea el log de Vina y extrae la afinidad de la mejor pose (modo 1)."""
    if not os.path.exists(log_path):
        return None
    
    with open(log_path, "r", encoding="utf-8") as f:
        content = f.read()
        
    # Buscar la primera línea de la tabla de resultados (modo 1)
    # Ejemplo:
    #    1      -8.4      0.000      0.000
    match = re.search(r"\s+1\s+(-?\d+\.\d+)\s+", content)
    if match:
        return float(match.group(1))
    return None

def cargar_tiempos_previos(csv_path):
    """Carga los tiempos de ejecución anteriores desde el CSV para no perderlos."""
    tiempos = {}
    if os.path.exists(csv_path):
        try:
            with open(csv_path, "r", encoding="utf-8") as f:
                lines = f.readlines()
                # Saltar la cabecera
                for line in lines[1:]:
                    parts = line.strip().split(",")
                    if len(parts) >= 5:
                        diana = parts[0]
                        ligando = parts[1]
                        try:
                            tiempo = float(parts[4])
                            tiempos[(diana, ligando)] = tiempo
                        except ValueError:
                            pass
        except Exception as e:
            print(f"[WARNING] No se pudieron leer algunos tiempos del CSV previo: {e}")
    return tiempos

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, os.pardir))
    output_dir = os.path.join(project_root, "docking", "vina_output")
    
    csv_path = os.path.join(output_dir, "resultados_docking.csv")
    md_path = os.path.join(output_dir, "reporte_docking.md")
    
    print("====================================================================")
    print("  PARSER INDEPENDIENTE DE RESULTADOS DE DOCKING (AutoDock Vina)")
    print("  Proyecto: Compuestos de Persea americana vs. MRSA")
    print("====================================================================")
    
    if not os.path.exists(output_dir):
        print(f"[ERROR] La carpeta de salidas no existe: {output_dir}")
        sys.exit(1)
        
    # 1. Cargar tiempos anteriores
    tiempos_previos = cargar_tiempos_previos(csv_path)
    if tiempos_previos:
        print(f"[INFO] Se cargaron {len(tiempos_previos)} registros de tiempos desde el CSV previo.")
    
    # 2. Escanear archivos *_log.txt
    archivos = os.listdir(output_dir)
    log_files = [f for f in archivos if f.endswith("_log.txt")]
    
    if not log_files:
        print("[WARNING] No se encontraron archivos de log (*_log.txt) en la carpeta de salida.")
        sys.exit(0)
        
    print(f"[INFO] Se encontraron {len(log_files)} archivos de log para procesar.")
    
    resultados = []
    
    # Definimos los nombres de las dianas conocidas
    dianas_keys = list(DIANAS.keys())
    
    for filename in sorted(log_files):
        log_path = os.path.join(output_dir, filename)
        
        # Determinar Diana y Ligando a partir del nombre del archivo
        diana_detectada = None
        for dk in dianas_keys:
            if filename.startswith(dk + "_"):
                diana_detectada = dk
                break
                
        if not diana_detectada:
            # Si no empieza con dianas conocidas, intentamos parsear por el primer guión bajo
            parts = filename.split("_")
            if len(parts) >= 3:
                diana_detectada = parts[0]
            else:
                print(f"[WARNING] No se pudo determinar la diana para el archivo: {filename}")
                continue
                
        # El ligando es todo lo que queda entre la diana y el sufijo _log.txt
        # Ejemplo: GyrB_07N_GyrB_3D_log.txt -> diana="GyrB", ligando="07N_GyrB_3D"
        prefix_len = len(diana_detectada) + 1
        ligando = filename[prefix_len:-8] # Quita {diana}_ del inicio y _log.txt del final
        
        # Parsear afinidad de energía
        afinidad = parsear_energia_vina(log_path)
        
        # Verificar si es un control
        is_control = "No"
        if diana_detectada in DIANAS:
            if ligando == DIANAS[diana_detectada]["control"]:
                is_control = "Sí"
                
        # Recuperar tiempo de ejecución anterior
        tiempo_s = tiempos_previos.get((diana_detectada, ligando), 0.0)
        
        resultados.append({
            "Diana": diana_detectada,
            "Ligando": ligando,
            "Control": is_control,
            "Energia_kcal_mol": afinidad,
            "Tiempo_s": tiempo_s
        })
        
        print(f"  - Parsed: {diana_detectada} | {ligando} | Control: {is_control} | ΔG: {afinidad} kcal/mol | Tiempo: {tiempo_s}s")
        
    if not resultados:
        print("[ERROR] No se pudo extraer ningún resultado válido.")
        sys.exit(1)
        
    # 3. Guardar resultados y generar reporte
    # Escribir CSV
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("Diana,Ligando,Control,Energia_kcal_mol,Tiempo_s\n")
        for r in resultados:
            # Si la energía es None, escribimos vacío o N/A
            e_str = f"{r['Energia_kcal_mol']}" if r['Energia_kcal_mol'] is not None else ""
            f.write(f"{r['Diana']},{r['Ligando']},{r['Control']},{e_str},{r['Tiempo_s']}\n")
            
    # Escribir Reporte Markdown
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# Reporte del Screening Virtual de Docking Molecular\n\n")
        f.write(f"**Fecha y Hora del Reporte:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("**Algoritmo:** AutoDock Vina 1.2.0\n")
        f.write("**Exhaustiveness:** 32 (Exhaustividad alta para publicación)\n\n")
        f.write("Este reporte consolida los resultados del docking molecular de los 5 compuestos de *Persea americana* y los controles positivos correspondientes frente a las tres dianas de *S. aureus* MRSA.\n\n")
        
        # Agrupar resultados por diana
        dianas_presentes = sorted(list(set(r["Diana"] for r in resultados)))
        
        for d in dianas_presentes:
            desc_diana = DIANAS[d]["desc"] if d in DIANAS else "Diana terapéutica"
            f.write(f"## Resultados para {d} ({desc_diana})\n\n")
            f.write("| Posición | Ligando | ¿Control? | Energía de Unión (ΔG, kcal/mol) | Tiempo (s) |\n")
            f.write("| :---: | :--- | :---: | :---: | :---: |\n")
            
            d_res = [r for r in resultados if r["Diana"] == d]
            # Ordenar de mayor a menor afinidad (energías más negativas primero)
            # Manejar posibles valores None asignando un valor positivo alto para que queden al final
            d_res_sorted = sorted(
                d_res,
                key=lambda x: x["Energia_kcal_mol"] if x["Energia_kcal_mol"] is not None else 999.0
            )
            
            # Buscar la energía del control para esta diana
            ctrl_energy = None
            for x in d_res:
                if x["Control"] == "Sí" and x["Energia_kcal_mol"] is not None:
                    ctrl_energy = x["Energia_kcal_mol"]
                    break
            
            for idx, r in enumerate(d_res_sorted):
                pos = idx + 1
                energia = f"{r['Energia_kcal_mol']:.1f}" if r["Energia_kcal_mol"] is not None else "Error"
                is_ctrl = "Sí" if r["Control"] == "Sí" else "No"
                tiempo = f"{r['Tiempo_s']:.1f}" if r["Tiempo_s"] > 0 else "N/A"
                
                # Resaltar si supera o iguala al control positivo
                if r["Control"] == "No" and ctrl_energy is not None and r["Energia_kcal_mol"] is not None:
                    if r["Energia_kcal_mol"] <= ctrl_energy:
                        energia = f"**{energia}** 🌟"
                        
                f.write(f"| {pos} | {r['Ligando']} | {is_ctrl} | {energia} | {tiempo} |\n")
                
            f.write("\n*Nota: Las energías más negativas indican mayor afinidad teórica. Las celdas en negrita y marcadas con 🌟 representan compuestos de interés que igualan o superan la afinidad del control positivo de esa diana.*\n\n")
            
    print("\n====================================================================")
    print("  REPORTE CONSOLIDADO CON ÉXITO")
    print(f"  Resultados consolidados en CSV: {csv_path}")
    print(f"  Reporte regenerado en Markdown: {md_path}")
    print("====================================================================")

if __name__ == "__main__":
    main()
