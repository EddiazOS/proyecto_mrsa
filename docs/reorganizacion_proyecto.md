# Reporte de Reorganización y Plan del Proyecto Virtual Screening Persea

Este documento sirve como constancia del plan implementado y de las acciones realizadas para la reestructuración del espacio de trabajo del proyecto *Virtual Screening Persea* contra dianas de *S. aureus* MRSA.

---

## 1. Diseño y Estructura del Workspace

El directorio raíz del proyecto (`/home/echoes/projects/cursos-semillero/VirtualScreening_Persea/`) ha sido estructurado siguiendo los mejores estándares de bioinformática estructural:

```text
VirtualScreening_Persea/
├── data/                            # Datos estructurales primarios
│   ├── receptors/                   # Receptores proteicos (originales y purificados)
│   │   ├── PBP2a_3ZG0.pdb           # Receptor PBP2a completo original
│   │   ├── PBP2a_3ZG0_clean.pdb     # Receptor PBP2a purificado (sin solvente ni ligandos)
│   │   ├── GyrB_3TTZ.pdb            # Receptor GyrB completo original
│   │   ├── GyrB_3TTZ_clean.pdb      # Receptor GyrB purificado (sin solvente ni ligandos)
│   │   └── MurG_AF-Q6GGZ0-F1.pdb    # Receptor MurG AlphaFold
│   └── ligands/                     # Ligandos en formato 3D minimizado
│       ├── Afzelin_3D.sdf           # Metabolito Persea (3D)
│       ├── 4Hydroxybenzoic_acid_3D.sdf
│       ├── Avocadenofuran_3D.sdf
│       ├── Avocadyne_acetate_3D.sdf
│       ├── Avocadene_acetate_3D.sdf
│       ├── Ceftaroline_3D.pdb       # Control positivo extraído de 3ZG0 (Ceftarolina)
│       ├── 07N_GyrB_3D.pdb          # Control positivo extraído de 3TTZ (Pirrolamida 07N)
│       └── raw_2d/                  # Estructuras 2D originales de PubChem (resguardo)
│           ├── Avocadyne_acetate_2D.sdf
│           └── Avocadene_acetate_2D.sdf
├── admet/                           # Espacio para análisis de drug-likeness (SwissADME)
├── docking/                         # Preparación y ejecución de docking molecular
│   ├── vina_input/                  # Archivos PDBQT, redes y parámetros (.txt)
│   └── vina_output/                 # Modos de pose resultantes de AutoDock Vina
├── md/                              # Simulaciones de Dinámica Molecular (GROMACS)
│   ├── MurG/
│   ├── PBP2a/
│   └── GyrB/
├── scripts/                         # Scripts de automatización y procesamiento
│   ├── setup_project.py             # Script de organización de archivos
│   └── extract_controls.py          # Script extractor de ligandos co-cristalizados
├── context.md                       # Documento metodológico de contexto original
├── README.txt                       # Información y procedencia de PubChem
└── reorganizacion_proyecto.md       # Este documento (constancia del proceso)
```

---

## 2. Resumen de Acciones Implementadas

### A. Creación y Clasificación de Datos
1. **Creación de Subdirectorios:** Se crearon de forma limpia las carpetas correspondientes a cada etapa del estudio (`data/receptors`, `data/ligands/raw_2d`, `admet`, `docking/vina_input`, `docking/vina_output`, `md/MurG`, `md/PBP2a`, `md/GyrB`).
2. **Reubicación de Receptores:** Se clasificaron y renombraron las dianas de estudio:
   - `3ZG0.pdb` → `data/receptors/PBP2a_3ZG0.pdb`
   - `3TTZ.pdb` → `data/receptors/GyrB_3TTZ.pdb`
   - `AF-Q6GGZ0-F1-model_v6.pdb` → `data/receptors/MurG_AF-Q6GGZ0-F1.pdb`
3. **Normalización de Ligandos:** Se renombraron los metabolitos de *Persea* a nombres legibles y estandarizados (eliminando CIDs y marcando claramente su formato 3D). Las versiones en 2D puro se trasladaron a `data/ligands/raw_2d/` para evitar confusiones en los dockings.

### B. Depuración de Workspace
Se eliminaron de forma recursiva **11 archivos de metadatos NTFS de Windows** (`*:Zone.Identifier`) que se habían generado automáticamente durante las descargas de los archivos de PubChem y PDB. Esto purga el espacio de trabajo de archivos redundantes e incompatibles con Git y Linux.

### C. Extracción Automatizada de Controles y Purificación (Valor Agregado)
Se diseñaron dos scripts en `scripts/` para procesar las estructuras biológicas:
1. **`extract_controls.py` (Procesamiento Molecular):**
   - **Extracción de Ceftarolina (`AI8`):** Extrajo el ligando co-cristalizado de PBP2a (`3ZG0`) generando el archivo `data/ligands/Ceftaroline_3D.pdb` (78 átomos) conservando los enlaces `CONECT` necesarios para docking y visualización.
   - **Extracción de Pirrolamida 07N (`07N`):** Extrajo el inhibidor de GyrB (`3TTZ`) generando el archivo `data/ligands/07N_GyrB_3D.pdb` (52 átomos).
   - **Purificación de Receptores:** Generó copias limpias de las proteínas (`PBP2a_3ZG0_clean.pdb` y `GyrB_3TTZ_clean.pdb`) removiendo solventes (aguas), iones cristalinos y la molécula de ligando original. Esto previene interferencias físicas durante los ensayos de docking.
2. **`setup_project.py` (Orquestador):**
   - Coordina de forma jerárquica la creación de las carpetas, la migración de los archivos, la limpieza de metadatos e invoca al script de procesamiento molecular.

---

## 3. Próximos Pasos Recomendados en el Protocolo

1. **Filtro ADMET (`admet/`):**
   - Usar las estructuras 3D optimizadas en `data/ligands/` para realizar la predicción de absorción, distribución, metabolismo y excreción en SwissADME (`swissadme.ch`). Guardar las hojas de cálculo resultantes en `admet/`.
2. **Preparación de Docking (`docking/vina_input/`):**
   - Convertir los receptores limpios (`*_clean.pdb` y MurG AF) a formato `.pdbqt` usando `prepare_receptor4.py`.
   - Convertir los metabolitos y controles de SDF/PDB a `.pdbqt` usando `prepare_ligand4.py`.
   - Configurar la caja de docking (`grid box`) de acuerdo a las coordenadas de los sitios activos detalladas en `context.md` y escribir los archivos de configuración `conf.txt`.
