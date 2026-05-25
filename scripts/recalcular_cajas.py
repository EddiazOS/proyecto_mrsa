#!/usr/bin/env python3
"""
Script para recalcular los tamaños de caja de todas las dianas,
basándose en las dimensiones reales de los ligandos del estudio.

Para cada diana (PBP2a, GyrB, MurG), el tamaño de caja se define como:
  max(tamaño_necesario_por_ligando) para todos los ligandos que se
  acoplarán a esa diana.

Cada ligando individual requiere:
  tamaño = dimensión_máxima + 2 × margen (10 Å por defecto)

Este script NO modifica archivos de configuración directamente.
Solo reporta los tamaños recomendados.

Uso:
    python scripts/recalcular_cajas.py
"""

import os
import sys

# Agregar directorio de scripts al path para importar calcular_centro
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from calcular_centro import (
    calcular_tamaño_caja,
    calcular_tamaño_caja_sdf,
)


def obtener_tamaño_ligando(ruta, margen=10.0, minimo=22.0):
    """
    Determina el tamaño de caja necesario para un ligando,
    detectando automáticamente el formato (PDB o SDF).
    Para PDB, filtra por cadena A para evitar promediar copias del homodímero.
    """
    if ruta.endswith('.sdf'):
        return calcular_tamaño_caja_sdf(ruta, margen=margen, minimo=minimo)
    elif ruta.endswith('.pdb'):
        return calcular_tamaño_caja(ruta, margen=margen, minimo=minimo,
                                    chain_id='A')
    else:
        print(f"  ⚠ Formato no soportado: {ruta}")
        return None


def main():
    base_dir = os.path.dirname(script_dir)
    ligands_dir = os.path.join(base_dir, "data", "ligands")

    # ---- Definición de ligandos por diana ----
    # Los 5 compuestos de prueba se acoplan a las 3 dianas.
    # Cada diana tiene además su control positivo específico.
    compuestos_prueba = [
        "Avocadenofuran_3D.sdf",
        "Avocadyne_acetate_3D.sdf",
        "Avocadene_acetate_3D.sdf",
    ]

    # Otros ligandos disponibles que podrían incluirse
    otros_ligandos = [
        "Quercetin_3D.sdf",
        "Afzelin_3D.sdf",
        "4Hydroxybenzoic_acid_3D.sdf",
    ]

    controles_por_diana = {
        "PBP2a": "Ceftaroline_3D.pdb",
        "GyrB":  "07N_GyrB_3D.pdb",
        "MurG":  "Quercetin_3D.sdf",
    }

    print("=" * 70)
    print("RECÁLCULO DE TAMAÑOS DE CAJA POR DIANA")
    print(f"Margen: 10 Å por lado | Mínimo absoluto: 22 Å")
    print("=" * 70)

    # ---- Calcular tamaño individual de cada ligando ----
    print("\n--- Tamaños individuales por ligando ---\n")
    tamaños_individuales = {}

    todos_archivos = list(set(
        compuestos_prueba +
        list(controles_por_diana.values()) +
        otros_ligandos
    ))
    todos_archivos.sort()

    for archivo in todos_archivos:
        ruta = os.path.join(ligands_dir, archivo)
        if not os.path.exists(ruta):
            print(f"  ⚠ No encontrado: {archivo}")
            continue

        tamaño = obtener_tamaño_ligando(ruta)
        if tamaño is not None:
            tamaños_individuales[archivo] = tamaño
            print(f"  {archivo:<35s} → {tamaño:.1f} Å")

    # ---- Calcular tamaño por diana ----
    print(f"\n{'='*70}")
    print("--- Tamaño de caja recomendado por diana ---\n")

    resultados_diana = {}

    for diana, control in controles_por_diana.items():
        ligandos_diana = compuestos_prueba + [control]
        tamaños_diana = []

        for lig in ligandos_diana:
            if lig in tamaños_individuales:
                tamaños_diana.append(tamaños_individuales[lig])

        if tamaños_diana:
            caja_diana = max(tamaños_diana)
        else:
            caja_diana = 22.0  # fallback

        resultados_diana[diana] = caja_diana
        print(f"  {diana:<8s} → {caja_diana:.1f} Å  "
              f"(máx de {len(tamaños_diana)} ligandos)")

    # ---- Tamaño unificado (máximo global) ----
    if resultados_diana:
        caja_global = max(resultados_diana.values())
    else:
        caja_global = 22.0

    print(f"\n  {'GLOBAL':<8s} → {caja_global:.1f} Å  "
          f"(máximo entre todas las dianas)")

    # ---- Resumen final ----
    print(f"\n{'='*70}")
    print("RESUMEN: TAMAÑOS DE CAJA PARA ARCHIVOS DE CONFIGURACIÓN")
    print(f"{'='*70}")
    print(f"\n  Diana     Caja actual    Caja recalculada    ¿Cambio?")
    print(f"  {'─'*55}")

    cajas_actuales = {"PBP2a": 22.0, "GyrB": 22.0, "MurG": 24.0}

    for diana in ["PBP2a", "GyrB", "MurG"]:
        actual = cajas_actuales.get(diana, 22.0)
        nueva = resultados_diana.get(diana, caja_global)
        cambio = "✓ SÍ" if abs(nueva - actual) > 0.1 else "  NO"
        print(f"  {diana:<8s}   {actual:>6.1f} Å       {nueva:>8.1f} Å        {cambio}")

    print(f"\n  NOTA: Se recomienda usar el tamaño GLOBAL ({caja_global:.1f} Å) para")
    print(f"  las 3 dianas, garantizando consistencia metodológica.\n")

    # ---- Valores para copiar ----
    print("Valores para conf_*.txt:")
    print(f"  size_x = {caja_global:.1f}")
    print(f"  size_y = {caja_global:.1f}")
    print(f"  size_z = {caja_global:.1f}")


if __name__ == "__main__":
    main()
