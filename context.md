# Estructura general.

proyecto_mrsa/
  data/
    receptors/
      MurG_AF-Q6GGZ0-F1.pdb
      PBP2a_3ZG0.pdb
      GyrB_3TTZ.pdb       # opcional
    ligands/
      Afzelin_3D.sdf
      4Hydroxybenzoic_acid_3D.sdf
      Avocadenofuran_3D.sdf
      Avocadyne_acetate_3D.sdf
      Avocadene_acetate_3D.sdf
      Quercetin_3D.sdf            # control MurG
      Ceftaroline_3D.sdf          # control PBP2a
      07N_GyrB_3D.sdf             # control GyrB (desde 3TTZ)
  admet/
  docking/
    vina_input/
    vina_output/
  md/
    MurG/
    PBP2a/
    GyrB/

| Criterio                   | MurG (AF-Q6GGZ0)            | PBP2a (3ZG0)                       | GyrB (3TTZ)                               |
| -------------------------- | --------------------------- | ---------------------------------- | ----------------------------------------- |
| Organismo                  | S. aureus MRSA252           | MRSA                               | S. aureus                                 |
| Resolución                 | N/A (AlphaFold, pLDDT 93)   | 2.37 Å                             | 1.63 Å ← mejor                            |
| Ligando co-cristalizado    | Sin ligando                 | Ceftarolina                        | Pirrolamida 07N                           |
| Mecanismo                  | Biosíntesis pared celular   | Resistencia β-lactámica            | Replicación DNA                           |
| Precedente con flavonoides | Quercetina pmc.ncbi.nlm.nih | Kaempferol, Catequina ncbi.nlm.nih | Limitado en S. aureus                     |
| Relevancia para MRSA       | Alta (pared celular)        | Alta (diana definitoria)           | Moderada (conservada, no específica MRSA) |

Perfecto. Aquí está la metodología completa y ordenada, consolidando todo lo que hemos trabajado, con los parámetros específicos que la literatura reciente respalda.

***

## Diseño general del estudio computacional

El estudio es un análisis **in silico** de tipo estructura-actividad (SBDD), compuesto por cuatro etapas secuenciales: (1) preparación y validación de ligandos y receptores, (2) filtro de propiedades fisicoquímicas ADMET, (3) docking molecular con AutoDock Vina 1.2.0, y (4) simulaciones de dinámica molecular (MD) con cálculo de energía libre de unión MM-GBSA. El objetivo es identificar cuáles de los cinco compuestos de *Persea americana* presentan mayor potencial inhibitorio frente a las tres dianas de *S. aureus* MRSA seleccionadas. [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC10683950/)

***

## Fase 1 — Preparación de ligandos y receptores

### Ligandos

Los cinco compuestos identificados por HPLC-QTof (Afzelin, Ácido 4-Hidroxibenzoico, Avocadenofuran, Avocadyne acetate, Avocadene acetate) se obtuvieron desde PubChem en formato SDF. Las estructuras 2D se convirtieron a confórmeros 3D utilizando el algoritmo **ETKDGv3 de RDKit** seguido de minimización energética con campo de fuerzas **MMFF94**. Se añaden hidrógenos explícitos antes de la minimización. Los ligandos de control positivo — **Quercetina** (control para MurG ), **Ceftarolina** (control para PBP2a ), y el inhibidor **07N** co-cristalizado en 3TTZ (control para GyrB) — se preparan de forma idéntica. Finalmente, cada ligando se convierte a formato **.pdbqt** usando `prepare_ligand4.py` de MGLTools con el parámetro `-A checkhydrogens`. [downloads.hindawi](https://downloads.hindawi.com/journals/ijmicro/2022/9130700.pdf)

### Receptores

Se utilizan tres dianas estructurales:

| Diana | Fuente | Resolución / Score | Justificación |
|---|---|---|---|
| **MurG** (Q6GGZ0) | AlphaFold AF-Q6GGZ0-F1 | pLDDT 93.56 | Validado con Quercetina en literatura  [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC11873250/) |
| **PBP2a** | PDB 3ZG0 | 2.37 Å | Diana definitoria de resistencia MRSA  [downloads.hindawi](https://downloads.hindawi.com/journals/ijmicro/2022/9130700.pdf) |
| **GyrB** (*S. aureus*) | PDB 3TTZ | 1.63 Å | Estructura nativa de alta resolución; inhibidor co-cristalizado 07N |

Cada receptor se prepara eliminando moléculas de agua (excepto las del sitio activo si forman puentes clave) e iones no estructurales con PyMOL (`remove solvent`). La protonación se asigna a pH 7.4 fisiológico. La conversión a PDBQT se realiza con `prepare_receptor4.py -A checkhydrogens`. [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC10683950/)

***

## Fase 2 — Filtro ADMET

Antes del docking, los cinco compuestos se evalúan con **SwissADME** (swissadme.ch) y **pkCSM** para los siguientes parámetros, usando los criterios de Lipinski y Veber como filtros de drug-likeness:

| Parámetro | Umbral Lipinski/Veber |
|---|---|
| Peso molecular | ≤ 500 g/mol |
| LogP (lipofilia) | ≤ 5 |
| HBD (dadores H) | ≤ 5 |
| HBA (aceptores H) | ≤ 10 |
| Rot. bonds | ≤ 10 (Veber) |
| PSA | ≤ 140 Å² |

Los avocadenoides (14 enlaces rotables) se incluyen igualmente en el docking pero su violación de la regla de Veber se reporta explícitamente como limitación, dado que las reglas de drug-likeness fueron diseñadas para moléculas sintéticas de administración oral, no para compuestos naturales. Esta es una posición estándar en la literatura reciente de productos naturales. [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC12466889/)

***

## Fase 3 — Docking molecular con AutoDock Vina 1.2.0

### Validación del protocolo (paso obligatorio previo)

Antes de dockar los compuestos de interés, se realiza **redocking del ligando co-cristalizado** para cada diana experimental (Ceftarolina en 3ZG0, 07N en 3TTZ). El criterio de aceptación es un RMSD ≤ 2.0 Å entre la pose re-docked y la pose cristalográfica. Si el RMSD supera este umbral, se ajustan los parámetros de la caja antes de continuar. Para MurG (AlphaFold), la validación se realiza re-docking de Quercetina y comparando los residuos de interacción con los reportados por Kim et al. 2025. [bio-protocol](https://bio-protocol.org/exchange/minidetail?id=2531599&type=30)

### Definición de cajas de búsqueda (grid boxes)

Las coordenadas se definen así para cada diana:

**PBP2a (3ZG0):** Centro en el sitio alostérico definido por las coordenadas del ligando co-cristalizado (Ceftarolina). Tamaño de caja: 22 × 22 × 22 Å. [bio-protocol](https://bio-protocol.org/exchange/minidetail?id=2531599&type=30)

**GyrB (3TTZ):** Centro en el sitio de unión de ATP (dominio ATPasa N-terminal), definido por el inhibidor 07N co-cristalizado. Tamaño de caja: 22 × 22 × 22 Å.

**MurG (AF-Q6GGZ0-F1):** Centro en el bolsillo de unión de UDP-GlcNAc, identificado por superposición estructural con el template 1F0K y las coordenadas reportadas en Kim et al. 2025. Tamaño de caja: 24 × 24 × 24 Å (mayor por tratarse de un modelo computacional sin ligando co-cristalizado). [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC11873250/)

### Parámetros fijos de Vina

```toml
exhaustiveness = 32      # estándar para rigor en publicación
num_modes      = 9       # poses a retener por corrida
energy_range   = 3       # kcal/mol entre mejor y peor pose
```

Estos parámetros se fijan **antes de ver los resultados** y no se modifican después, para evitar sesgo de selección. [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC10683950/)

### Criterios de selección de poses

Se selecciona la pose de menor energía de unión (ΔG, kcal/mol) como representativa de cada complejo. Se reportan todos los modos y se analiza el RMSD entre las tres mejores poses para evaluar la consistencia conformacional. Para avanzar a MD, se seleccionan los **2–3 complejos con mejor ΔG por diana**, siempre que sean mejores que el control positivo o comparables a él.

***

## Fase 4 — Dinámica molecular (MD) con GROMACS

### Configuración del sistema

Los complejos seleccionados se simulan con **GROMACS 2023+** usando el campo de fuerzas **AMBER99SB-ILDN** para la proteína  y **GAFF2** para los ligandos (parámetros generados con `acpype` o `antechamber`). El sistema se solva en caja octaédrica truncada con agua explícita TIP3P, neutralizando con iones Na⁺/Cl⁻ a 0.15 M de concentración fisiológica. [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC12240381/)

### Protocolo de simulación

```
1. Minimización energética       → método steepest descent, convergencia < 1000 kJ/mol/nm
2. Equilibración NVT             → 100 ps, 300 K, termostato V-rescale
3. Equilibración NPT             → 100 ps, 1 bar, barostato Parrinello-Rahman
4. Producción MD                 → 100 ns, paso de integración 2 fs
5. Análisis de trayectoria       → RMSD, RMSF, Rg, SASA, número de puentes H
6. Energía libre de unión        → MM-GBSA sobre últimos 50 ns (estado estacionario)
```

La simulación de 100 ns es el mínimo aceptable para flavonoides en dianas bacterianas según la literatura actual. Si los recursos computacionales lo permiten, se amplía a 200 ns para los mejores candidatos. [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC12389164/)

### Parámetros de análisis de trayectoria

```bash
# RMSD del backbone (Cα) del receptor
gmx rms -s topol.tpr -f traj.xtc -o rmsd_receptor.xvg

# RMSF por residuo (flexibilidad)
gmx rmsf -s topol.tpr -f traj.xtc -o rmsf.xvg -res

# Radio de giro (compactación global)
gmx gyrate -s topol.tpr -f traj.xtc -o gyrate.xvg

# MM-GBSA con gmx_MMPBSA
gmx_MMPBSA -O -i mmgbsa.in -cs topol.tpr -ct traj.xtc -ci index.ndx
```

El criterio de convergencia es que el RMSD del backbone se estabilice antes de los 20 ns. Si no converge, la simulación no es reportable sin ajuste del protocolo de equilibración.

***

## Controles del experimento computacional

| Control | Tipo | Diana | Propósito |
|---|---|---|---|
| Quercetina | Positivo | MurG | ΔG de referencia ya publicado  [pmc.ncbi.nlm.nih](https://pmc.ncbi.nlm.nih.gov/articles/PMC11873250/) |
| Ceftarolina | Positivo | PBP2a | Antibiótico MRSA aprobado, co-cristalizado |
| 07N | Positivo | GyrB | Inhibidor nanomolar co-cristalizado en 3TTZ |
| Sin ligando (apo) | Negativo MD | Todas | Comparar estabilidad basal del receptor |

***

## Flujo de trabajo completo (resumen ejecutivo)

```
SDF PubChem → ETKDGv3 (RDKit) → MMFF94 → PDBQT (MGLTools)
                                              ↓
PDB/AlphaFold → Limpieza PyMOL → PDBQT (MGLTools)
                                              ↓
                        SwissADME / pkCSM (filtro ADMET)
                                              ↓
                   Redocking validación (RMSD ≤ 2.0 Å)
                                              ↓
                AutoDock Vina 1.2.0 (exhaustiveness=32)
                                              ↓
              Selección top 2-3 complejos por diana
                                              ↓
           GROMACS 100 ns · AMBER99SB-ILDN + GAFF2
                                              ↓
        RMSD / RMSF / Rg / SASA / MM-GBSA (50 ns finales)
```

¿Quieres que ahora desarrollemos el script de Colab para ejecutar Vina de forma automatizada para los 5 ligandos × 3 dianas (15 experimentos de docking), o prefieres que ajustemos primero la definición numérica de las cajas de cada diana?
