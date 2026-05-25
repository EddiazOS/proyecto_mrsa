#!/usr/bin/env python3
"""
validar_sitio_murg.py
=====================
Validación del sitio activo de MurG (modelo AlphaFold, MRSA S. aureus)
mediante superposición estructural con el cristal 1F0K (E. coli).

Pasos:
  1. Carga de estructuras
  2. Alineamiento de secuencia → mapeo de residuos
  3. Superposición global (todos los CA alineados)
  4. Superposición local (residuos del sitio activo)
  5. Centroide del sitio activo en el marco de referencia de MurG_MRSA
  6. Verificación de pLDDT en residuos del sitio activo
  7. Resumen formateado

Uso:
  python scripts/validar_sitio_murg.py
  (ejecutar desde la raíz del proyecto)

Autor : Semillero Bioinformática
Fecha : 2026
"""

import os
import sys
import warnings
import copy

import numpy as np

from Bio.PDB import PDBParser, Superimposer
from Bio.PDB import PPBuilder
from Bio.Align import PairwiseAligner
from Bio import BiopythonWarning

# ── Silenciar warnings menores de BioPython ──────────────────────────
warnings.simplefilter("ignore", BiopythonWarning)

# =====================================================================
# Configuración
# =====================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir))

AF_PDB = os.path.join(PROJECT_ROOT, "data", "receptors", "MurG_AF-Q6GGZ0-F1.pdb")
TEMPLATE_PDB = os.path.join(PROJECT_ROOT, "data", "receptors", "MurG_template_1F0K.pdb")

# Residuos del sitio activo en la numeración de 1F0K (E. coli)
ACTIVE_SITE_RESIDS_1F0K = [14, 15, 169, 269, 302, 317]

# Umbrales de aceptación
RMSD_GLOBAL_THRESHOLD = 3.0   # Å
RMSD_LOCAL_THRESHOLD  = 3.5   # Å
PLDDT_THRESHOLD       = 70.0

# Nombres de HETATM a excluir al buscar ligandos
EXCLUDE_HETATM = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO",
                  "ACT", "NA", "CL", "MG", "CA", "ZN", "MN", "FE",
                  "K", "IOD", "BR", "NI", "CO", "CU"}

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


# =====================================================================
# Funciones auxiliares
# =====================================================================

def separador(titulo: str) -> None:
    """Imprime un separador con título."""
    ancho = 68
    print(f"\n{'=' * ancho}")
    print(f"  {titulo}")
    print(f"{'=' * ancho}")


def cargar_estructura(ruta: str, nombre: str):
    """Carga una estructura PDB y devuelve (structure, model)."""
    parser = PDBParser(QUIET=True)
    estructura = parser.get_structure(nombre, ruta)
    modelo = estructura[0]
    return estructura, modelo


def extraer_secuencia(modelo) -> tuple:
    """
    Extrae la secuencia de aminoácidos y la lista de residuos de la
    primera cadena proteica usando Polypeptides.

    Devuelve:
        (secuencia_str, lista_de_residuos)
    """
    ppbuilder = PPBuilder()
    # Tomar la primera cadena que tenga polipéptidos
    for cadena in modelo.get_chains():
        pps = ppbuilder.build_peptides(cadena)
        if pps:
            seq = ""
            residuos = []
            for pp in pps:
                seq += str(pp.get_sequence())
                residuos.extend(pp)
            return seq, residuos
    raise ValueError("No se encontraron polipéptidos en el modelo.")


def buscar_hetatm_ligando(modelo):
    """
    Busca residuos HETATM que podrían ser el ligando UDP-GlcNAc.
    Excluye agua e iones comunes.

    Devuelve: lista de (resname, chain_id, res_id, n_atoms)
    """
    ligandos = []
    for cadena in modelo.get_chains():
        for residuo in cadena.get_residues():
            hetflag = residuo.get_id()[0]
            if hetflag.startswith("H_") or hetflag == "W":
                resname = residuo.get_resname().strip()
                if resname not in EXCLUDE_HETATM:
                    n_atoms = len(list(residuo.get_atoms()))
                    ligandos.append((resname, cadena.id, residuo.get_id(), n_atoms))
    return ligandos


# =====================================================================
# PASO 1 — Carga de estructuras
# =====================================================================

def paso1_cargar():
    """Carga ambas estructuras PDB."""
    separador("PASO 1: Carga de estructuras")

    if not os.path.isfile(AF_PDB):
        print(f"  ✗ ERROR: No se encontró el modelo AlphaFold: {AF_PDB}")
        sys.exit(1)
    if not os.path.isfile(TEMPLATE_PDB):
        print(f"  ✗ ERROR: No se encontró el template 1F0K: {TEMPLATE_PDB}")
        sys.exit(1)

    struct_af, model_af = cargar_estructura(AF_PDB, "MurG_MRSA")
    struct_tmpl, model_tmpl = cargar_estructura(TEMPLATE_PDB, "1F0K")

    # Información básica
    for tag, mdl, ruta in [("AlphaFold (MRSA)", model_af, AF_PDB),
                           ("Template 1F0K (E. coli)", model_tmpl, TEMPLATE_PDB)]:
        cadenas = list(mdl.get_chains())
        n_res = sum(1 for r in mdl.get_residues()
                    if r.get_id()[0] == " ")
        print(f"\n  ▸ {tag}")
        print(f"    Archivo  : {os.path.basename(ruta)}")
        print(f"    Cadenas  : {', '.join(c.id for c in cadenas)}")
        print(f"    Residuos : {n_res}")

    return struct_af, model_af, struct_tmpl, model_tmpl


# =====================================================================
# PASO 2 — Alineamiento de secuencia y mapeo de residuos
# =====================================================================

def paso2_mapeo(model_af, model_tmpl):
    """
    Alinea las secuencias y mapea los residuos del sitio activo de
    1F0K a la numeración de MurG_MRSA.

    Devuelve:
        mapeo : dict {resid_1F0K: resid_MRSA}
        residuos_af : list de objetos Residue de AlphaFold
        residuos_tmpl : list de objetos Residue de 1F0K
    """
    separador("PASO 2: Alineamiento de secuencia y mapeo de residuos")

    seq_af, res_af = extraer_secuencia(model_af)
    seq_tmpl, res_tmpl = extraer_secuencia(model_tmpl)

    print(f"\n  Longitud secuencia MRSA  : {len(seq_af)} aa")
    print(f"  Longitud secuencia 1F0K : {len(seq_tmpl)} aa")

    # Alineamiento global con PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    alineamientos = aligner.align(seq_af, seq_tmpl)
    aln = alineamientos[0]

    # Calcular identidad
    aln_af_str, aln_tmpl_str = str(aln).split("\n")[0], str(aln).split("\n")[2]
    # Usar aligned property para obtener mapeo
    # Construir mapeo posición-a-posición desde el alineamiento
    # aln.aligned da bloques de coordenadas alineadas

    # Construir arrays de índice para cada posición alineada
    idx_af = []    # índices en seq_af
    idx_tmpl = []  # índices en seq_tmpl

    # Recorrer el alineamiento carácter por carácter
    aligned_seqs = aln.format().split("\n")
    # El formato puede variar; usemos indices directamente
    pos_af = 0
    pos_tmpl = 0
    pares = []  # (pos_af, pos_tmpl) para posiciones alineadas

    # Usar aln.indices que da un array 2×L con índices (-1 = gap)
    # Disponible en Bio.Align moderno
    try:
        indices = aln.indices  # shape (2, L)
        for col in range(indices.shape[1]):
            i_af = indices[0, col]
            i_tmpl = indices[1, col]
            if i_af >= 0 and i_tmpl >= 0:
                pares.append((int(i_af), int(i_tmpl)))
    except AttributeError:
        # Fallback: parsear el alineamiento formateado
        lines = aln.format().strip().split("\n")
        seq_top = ""
        seq_bot = ""
        for i, line in enumerate(lines):
            if i % 4 == 0:
                seq_top += line.split()[-1] if line.strip() else ""
            elif i % 4 == 2:
                seq_bot += line.split()[-1] if line.strip() else ""
        pos_af = 0
        pos_tmpl = 0
        for c_af, c_tmpl in zip(seq_top, seq_bot):
            if c_af != "-" and c_tmpl != "-":
                pares.append((pos_af, pos_tmpl))
            if c_af != "-":
                pos_af += 1
            if c_tmpl != "-":
                pos_tmpl += 1

    # Mapear residuos del sitio activo
    # Necesitamos ir de resid_1F0K → posición en seq_tmpl → posición en seq_af → resid_MRSA
    # Construir mapeo: posición_secuencia_tmpl → resid_1F0K
    pos_to_resid_tmpl = {}
    for i, res in enumerate(res_tmpl):
        pos_to_resid_tmpl[i] = res.get_id()[1]

    resid_to_pos_tmpl = {v: k for k, v in pos_to_resid_tmpl.items()}

    pos_to_resid_af = {}
    for i, res in enumerate(res_af):
        pos_to_resid_af[i] = res.get_id()[1]

    # Mapeo tmpl_pos → af_pos
    tmpl_to_af = {}
    for pos_af_i, pos_tmpl_i in pares:
        tmpl_to_af[pos_tmpl_i] = pos_af_i

    # Identidad de secuencia
    matches = sum(1 for pa, pt in pares if seq_af[pa] == seq_tmpl[pt])
    identidad = 100.0 * matches / max(len(pares), 1)
    print(f"  Posiciones alineadas     : {len(pares)}")
    print(f"  Identidad de secuencia   : {identidad:.1f}%")

    # Mapeo del sitio activo
    mapeo = {}
    print(f"\n  {'Residuo 1F0K':>15s}  →  {'Residuo MRSA':<15s}  {'Coincide?'}")
    print(f"  {'─' * 55}")

    for resid_tmpl in ACTIVE_SITE_RESIDS_1F0K:
        if resid_tmpl not in resid_to_pos_tmpl:
            print(f"  {'?':>15s}  →  {'?':<15s}  ✗ No encontrado en template")
            continue

        pos_t = resid_to_pos_tmpl[resid_tmpl]
        if pos_t not in tmpl_to_af:
            print(f"  {resid_tmpl:>15d}  →  {'gap':<15s}  ✗ Gap en alineamiento")
            continue

        pos_a = tmpl_to_af[pos_t]
        resid_af = pos_to_resid_af[pos_a]
        mapeo[resid_tmpl] = resid_af

        aa_tmpl = seq_tmpl[pos_t]
        aa_af = seq_af[pos_a]
        match = "✓" if aa_tmpl == aa_af else "≈"
        res_name_tmpl = THREE_TO_ONE.get(res_tmpl[pos_t].get_resname(), "?") if pos_t < len(res_tmpl) else "?"
        print(f"  {aa_tmpl}{resid_tmpl:>4d} (1F0K)    →  "
              f"{aa_af}{resid_af:>4d} (MRSA)    {match}")

    if len(mapeo) < len(ACTIVE_SITE_RESIDS_1F0K):
        print(f"\n  ⚠ ADVERTENCIA: No se pudieron mapear todos los residuos "
              f"del sitio activo ({len(mapeo)}/{len(ACTIVE_SITE_RESIDS_1F0K)})")

    return mapeo, res_af, res_tmpl, pares


# =====================================================================
# PASO 3 — Superposición global
# =====================================================================

def paso3_superposicion_global(model_af, model_tmpl, res_af, res_tmpl, pares):
    """
    Superposición global usando todos los CA alineados.
    Alinea 1F0K sobre MurG_MRSA (para que las coordenadas del ligando
    queden en el marco de referencia de MurG_MRSA).

    Devuelve:
        sup       : Superimposer ya ajustado
        rmsd_glob : RMSD global
    """
    separador("PASO 3: Superposición global (todos los CA alineados)")

    # Construir listas de átomos CA pareados
    ca_af = []
    ca_tmpl = []

    # Crear diccionarios de posición → residuo para acceso rápido
    res_af_dict = {i: r for i, r in enumerate(res_af)}
    res_tmpl_dict = {i: r for i, r in enumerate(res_tmpl)}

    for pos_a, pos_t in pares:
        r_af = res_af_dict.get(pos_a)
        r_tmpl = res_tmpl_dict.get(pos_t)
        if r_af is None or r_tmpl is None:
            continue
        if "CA" in r_af and "CA" in r_tmpl:
            ca_af.append(r_af["CA"])
            ca_tmpl.append(r_tmpl["CA"])

    print(f"\n  Pares de CA utilizados: {len(ca_af)}")

    if len(ca_af) < 10:
        print("  ✗ ERROR: Muy pocos CA pareados para superposición.")
        sys.exit(1)

    # Superimposer: fijo = MurG_MRSA, móvil = 1F0K
    # Queremos mover 1F0K al marco de referencia de MurG_MRSA
    sup = Superimposer()
    sup.set_atoms(ca_af, ca_tmpl)  # fijo, móvil

    rmsd_glob = sup.rms
    print(f"  RMSD global: {rmsd_glob:.3f} Å")

    if rmsd_glob < RMSD_GLOBAL_THRESHOLD:
        print(f"  ✓ RMSD global < {RMSD_GLOBAL_THRESHOLD} Å → ACEPTABLE")
    else:
        print(f"  ⚠ ADVERTENCIA: RMSD global ≥ {RMSD_GLOBAL_THRESHOLD} Å → "
              f"La superposición global es pobre")

    # Aplicar transformación a TODOS los átomos de 1F0K
    sup.apply(model_tmpl.get_atoms())
    print("  Transformación aplicada a todos los átomos de 1F0K.")

    return sup, rmsd_glob


# =====================================================================
# PASO 4 — Superposición local (sitio activo)
# =====================================================================

def paso4_superposicion_local(model_af, model_tmpl, mapeo):
    """
    Superposición local usando solo los CA del sitio activo.

    Devuelve:
        rmsd_local : RMSD local
    """
    separador("PASO 4: Superposición local (residuos del sitio activo)")

    # Seleccionar la primera cadena proteica de cada modelo
    cadena_af = next(model_af.get_chains())
    cadena_tmpl = next(model_tmpl.get_chains())

    ca_af = []
    ca_tmpl = []
    print()
    for resid_tmpl, resid_af in sorted(mapeo.items()):
        try:
            r_af = cadena_af[(" ", resid_af, " ")]
            r_tmpl = cadena_tmpl[(" ", resid_tmpl, " ")]
        except KeyError as e:
            print(f"  ⚠ No se encontró el residuo {e} — omitido")
            continue

        if "CA" in r_af and "CA" in r_tmpl:
            ca_af.append(r_af["CA"])
            ca_tmpl.append(r_tmpl["CA"])
            dist = r_af["CA"] - r_tmpl["CA"]
            print(f"  {r_tmpl.get_resname()}{resid_tmpl:>4d} (1F0K) ↔ "
                  f"{r_af.get_resname()}{resid_af:>4d} (MRSA)  "
                  f"dist CA = {dist:.2f} Å")

    if len(ca_af) < 3:
        print("  ✗ ERROR: Muy pocos residuos para superposición local.")
        return float("inf")

    sup_local = Superimposer()
    sup_local.set_atoms(ca_af, ca_tmpl)  # fijo=AF, móvil=tmpl (ya transformados)
    rmsd_local = sup_local.rms

    print(f"\n  RMSD local (sitio activo): {rmsd_local:.3f} Å")

    if rmsd_local < RMSD_LOCAL_THRESHOLD:
        print(f"  ✓ RMSD local < {RMSD_LOCAL_THRESHOLD} Å → ACEPTABLE")
    else:
        print(f"  ⚠ ADVERTENCIA: RMSD local ≥ {RMSD_LOCAL_THRESHOLD} Å")
        print(f"    Considere usar fpocket para identificar bolsillos de unión")
        print(f"    alternativos en MurG_MRSA como respaldo.")

    return rmsd_local


# =====================================================================
# PASO 5 — Centroide del sitio activo
# =====================================================================

def paso5_centroide(model_af, model_tmpl, mapeo):
    """
    Calcula el centroide del sitio activo en el marco de referencia
    de MurG_MRSA.

    Estrategia:
      1. Buscar ligando UDP-GlcNAc en 1F0K (ya superpuesto sobre MRSA).
         Si existe → centroide del ligando.
      2. Si no hay ligando → centroide de los CA del sitio activo en
         MurG_MRSA directamente.

    Devuelve:
        (cx, cy, cz) : coordenadas del centroide
    """
    separador("PASO 5: Centroide del sitio activo")

    # --- Intentar encontrar ligando en 1F0K (ya superpuesto) ---
    ligandos = buscar_hetatm_ligando(model_tmpl)

    centroide = None

    if ligandos:
        print("\n  Ligandos encontrados en 1F0K (ya superpuesto):")
        for resname, chain_id, res_id, n_atoms in ligandos:
            print(f"    • {resname} cadena {chain_id}, "
                  f"id={res_id}, átomos={n_atoms}")

        # Usar el ligando más grande (mayor número de átomos)
        best = max(ligandos, key=lambda x: x[3])
        resname, chain_id, res_id, _ = best
        print(f"\n  Usando ligando: {resname} (cadena {chain_id})")

        # Extraer coordenadas
        cadena = model_tmpl[chain_id]
        residuo_lig = cadena[res_id]
        coords = np.array([a.get_vector().get_array()
                           for a in residuo_lig.get_atoms()])
        centroide = coords.mean(axis=0)
        print(f"  Centroide calculado a partir de {len(coords)} átomos del ligando.")
    else:
        print("\n  ℹ No se encontró ligando UDP-GlcNAc en 1F0K.")
        print("    Se usará el centroide de los CA del sitio activo de MurG_MRSA.")

        # Centroide de los CA del sitio activo en MurG_MRSA
        cadena_af = next(model_af.get_chains())
        coords = []
        for resid_af in sorted(mapeo.values()):
            try:
                r = cadena_af[(" ", resid_af, " ")]
                if "CA" in r:
                    coords.append(r["CA"].get_vector().get_array())
            except KeyError:
                pass

        if not coords:
            print("  ✗ ERROR: No se pudieron obtener coordenadas del sitio activo.")
            sys.exit(1)

        coords = np.array(coords)
        centroide = coords.mean(axis=0)
        print(f"  Centroide calculado a partir de {len(coords)} CA del sitio activo.")

    cx, cy, cz = centroide
    print(f"\n  ── Coordenadas del centro del sitio activo (marco MurG_MRSA) ──")
    print(f"  center_x = {cx:.3f}")
    print(f"  center_y = {cy:.3f}")
    print(f"  center_z = {cz:.3f}")

    return cx, cy, cz


# =====================================================================
# PASO 6 — Verificación de pLDDT
# =====================================================================

def paso6_plddt(model_af, mapeo):
    """
    Extrae pLDDT (B-factor en modelos AlphaFold) de los residuos del
    sitio activo.

    Devuelve:
        plddt_dict : dict {resid_MRSA: pLDDT}
        mean_plddt : float
    """
    separador("PASO 6: Verificación de pLDDT (confianza AlphaFold)")

    cadena_af = next(model_af.get_chains())
    plddt_dict = {}
    hay_warning = False

    print(f"\n  {'Residuo':>12s}  {'pLDDT':>7s}  {'Estado'}")
    print(f"  {'─' * 40}")

    for resid_tmpl, resid_af in sorted(mapeo.items()):
        try:
            r = cadena_af[(" ", resid_af, " ")]
            if "CA" in r:
                plddt = r["CA"].get_bfactor()
                plddt_dict[resid_af] = plddt
                estado = "✓" if plddt >= PLDDT_THRESHOLD else "⚠ BAJA CONFIANZA"
                if plddt < PLDDT_THRESHOLD:
                    hay_warning = True
                resname = r.get_resname()
                print(f"  {resname}{resid_af:>4d} (MRSA)  {plddt:>7.2f}  {estado}")
        except KeyError:
            print(f"  Residuo {resid_af} no encontrado")

    if not plddt_dict:
        print("  ✗ ERROR: No se obtuvieron valores de pLDDT.")
        return {}, 0.0

    mean_plddt = np.mean(list(plddt_dict.values()))
    print(f"\n  pLDDT promedio (sitio activo): {mean_plddt:.2f}")

    if hay_warning:
        print(f"  ⚠ ADVERTENCIA: Algunos residuos del sitio activo tienen "
              f"pLDDT < {PLDDT_THRESHOLD}")
        print(f"    La predicción de AlphaFold puede ser menos confiable "
              f"en estas regiones.")
    else:
        print(f"  ✓ Todos los residuos del sitio activo tienen pLDDT ≥ "
              f"{PLDDT_THRESHOLD}")

    return plddt_dict, mean_plddt


# =====================================================================
# PASO 7 — Resumen
# =====================================================================

def paso7_resumen(rmsd_glob, rmsd_local, cx, cy, cz,
                  plddt_dict, mean_plddt, mapeo):
    """Imprime un resumen final formateado."""
    separador("PASO 7: Resumen de validación")

    ancho = 62
    print(f"\n  ┌{'─' * ancho}┐")
    print(f"  │{'VALIDACIÓN DEL SITIO ACTIVO DE MurG (MRSA)':^{ancho}s}│")
    print(f"  ├{'─' * ancho}┤")

    # RMSD global
    estado_glob = "✓ PASA" if rmsd_glob < RMSD_GLOBAL_THRESHOLD else "✗ FALLA"
    print(f"  │ {'RMSD global':.<35s} {rmsd_glob:>7.3f} Å  {estado_glob:>10s} │")

    # RMSD local
    estado_local = "✓ PASA" if rmsd_local < RMSD_LOCAL_THRESHOLD else "✗ FALLA"
    print(f"  │ {'RMSD local (sitio activo)':.<35s} {rmsd_local:>7.3f} Å  {estado_local:>10s} │")

    # pLDDT
    estado_plddt = "✓ PASA" if all(v >= PLDDT_THRESHOLD for v in plddt_dict.values()) else "⚠ REVISAR"
    print(f"  │ {'pLDDT promedio':.<35s} {mean_plddt:>7.2f}     {estado_plddt:>10s} │")

    print(f"  ├{'─' * ancho}┤")
    print(f"  │{'COORDENADAS DEL SITIO ACTIVO':^{ancho}s}│")
    print(f"  ├{'─' * ancho}┤")
    print(f"  │ {'center_x':.<20s} {cx:>10.3f}{' ' * 30}│")
    print(f"  │ {'center_y':.<20s} {cy:>10.3f}{' ' * 30}│")
    print(f"  │ {'center_z':.<20s} {cz:>10.3f}{' ' * 30}│")

    print(f"  ├{'─' * ancho}┤")
    print(f"  │{'pLDDT POR RESIDUO DEL SITIO ACTIVO':^{ancho}s}│")
    print(f"  ├{'─' * ancho}┤")
    for resid_tmpl, resid_af in sorted(mapeo.items()):
        if resid_af in plddt_dict:
            plddt = plddt_dict[resid_af]
            marca = "✓" if plddt >= PLDDT_THRESHOLD else "⚠"
            linea = f"  Res {resid_af:>4d} (MRSA) ← {resid_tmpl:>4d} (1F0K)  pLDDT={plddt:.1f} {marca}"
            print(f"  │ {linea:<{ancho - 2}s} │")
    print(f"  └{'─' * ancho}┘")

    # Formato para conf_murg.txt
    print(f"\n  ── Para pegar en conf_murg.txt ──────────────────────────────")
    print(f"  center_x = {cx:.3f}")
    print(f"  center_y = {cy:.3f}")
    print(f"  center_z = {cz:.3f}")
    print()

    # Resultado final
    todo_ok = (rmsd_glob < RMSD_GLOBAL_THRESHOLD and
               rmsd_local < RMSD_LOCAL_THRESHOLD and
               all(v >= PLDDT_THRESHOLD for v in plddt_dict.values()))

    if todo_ok:
        print("  ✅ VALIDACIÓN EXITOSA: El sitio activo de MurG_MRSA es "
              "estructuralmente compatible con el template 1F0K.")
    else:
        print("  ⚠  VALIDACIÓN CON OBSERVACIONES: Revise los puntos "
              "marcados arriba antes de proceder con docking.")


# =====================================================================
# MAIN
# =====================================================================

def main():
    print("\n" + "▓" * 68)
    print("  VALIDACIÓN ESTRUCTURAL DEL SITIO ACTIVO DE MurG")
    print("  Modelo: AlphaFold Q6GGZ0 (S. aureus MRSA)")
    print("  Template: 1F0K (E. coli, cristal 1.9 Å)")
    print("▓" * 68)

    # Paso 1
    struct_af, model_af, struct_tmpl, model_tmpl = paso1_cargar()

    # Paso 2
    mapeo, res_af, res_tmpl, pares = paso2_mapeo(model_af, model_tmpl)

    # Necesitamos una copia profunda del model_tmpl para el paso 5,
    # ya que aplicaremos la transformación global en el paso 3
    # y queremos las coordenadas del ligando transformadas.
    # (Aplicamos la transformación a model_tmpl, luego extraemos
    #  las coordenadas del ligando — ya en el marco de MurG_MRSA.)

    # Paso 3 — mover 1F0K al marco de MurG_MRSA
    sup, rmsd_glob = paso3_superposicion_global(
        model_af, model_tmpl, res_af, res_tmpl, pares
    )

    # Paso 4 — superposición local (usando coords ya transformadas)
    rmsd_local = paso4_superposicion_local(model_af, model_tmpl, mapeo)

    # Paso 5 — centroide (1F0K ya está en el marco de MurG_MRSA)
    cx, cy, cz = paso5_centroide(model_af, model_tmpl, mapeo)

    # Paso 6 — pLDDT
    plddt_dict, mean_plddt = paso6_plddt(model_af, mapeo)

    # Paso 7 — resumen
    paso7_resumen(rmsd_glob, rmsd_local, cx, cy, cz,
                  plddt_dict, mean_plddt, mapeo)


if __name__ == "__main__":
    main()
