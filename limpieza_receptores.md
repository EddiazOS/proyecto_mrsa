# Guía Metodológica: Purificación de Receptores y Extracción de Controles

Este documento detalla el procedimiento bioinformático implementado para la purificación de los receptores proteicos (**PBP2a** y **GyrB**) y la extracción de sus respectivos ligandos de control co-cristalizados (**Ceftarolina** y **Pirrolamida 07N**).

---

## 1. Fundamento de la Limpieza Molecular

En cristalografía de rayos X, las estructuras depositadas en el *Protein Data Bank* (PDB) contienen no solo la proteína de interés, sino también:
- **Solvente:** Moléculas de agua de cristalización necesarias para estabilizar la red cristalina.
- **Iones y Aditivos:** Compuestos químicos añadidos para facilitar la nucleación del cristal (ej. sales, agentes precipitantes).
- **Ligandos:** Moléculas orgánicas (inhibidores, co-factores) co-cristalizados en el sitio activo para caracterizar el complejo.

Para llevar a cabo un estudio de **docking molecular (cribado virtual)**, es indispensable purificar el receptor eliminando estos elementos. De lo contrario, se producirían colisiones estéricas artificiales y errores graves en el cálculo de mapas de afinidad de AutoDock Vina.

---

## 2. Anatomía de un Archivo PDB y Lógica de Parseo

Un archivo `.pdb` es un formato de texto plano con columnas fijas. El script automatizado `scripts/extract_controls.py` procesa los archivos estructurales línea por línea basándose en las especificaciones del formato oficial de PDB:

### Registro ATOM
Representa los átomos estándar que constituyen la cadena polipeptídica de la proteína (carbono alfa, nitrógeno de amina, oxígeno de carbonilo, etc. de los aminoácidos).
- **Lógica:** Se conservan **todos** los registros que inicien exactamente con `ATOM`.

### Registro TER
Indica el término (finalización) de una cadena de aminoácidos o nucleótidos.
- **Lógica:** Se conservan **todos** los registros que inicien con `TER`.

### Registro HETATM (Heteroátomos)
Representa átomos pertenecientes a moléculas no estándar (solventes, iones, cofactores, ligandos).
- **Lógica:** El script evalúa el campo del **Nombre del Residuo** (columnas 18 a 20, indexación 0-based: `[17:21]`):
  1. Si coincide con el ligando de control (`AI8` o `07N`), extrae la línea y la guarda en el archivo del ligando independiente.
  2. Si es agua (`HOH` o `WAT`) o iones no estructurales (ej. `CD`, `CL`), la línea es **descartada**, purificando el receptor.

### Registro CONECT
Define los enlaces covalentes explícitos para átomos que no pertenecen a residuos estándar (los heteroátomos del ligando).
- **Lógica:** El script recopila los números de serie de los átomos de ligando extraídos y filtra los registros `CONECT` al final del archivo original, conservando solo aquellos que involucran los átomos del ligando para garantizar que el PDB del control positivo conserve sus enlaces correctos al ser cargado en PyMOL o AutoDock.

---

## 3. Detalle por Diana de Estudio

### A. PBP2a (Origen: PDB `3ZG0`)
- **Resolución:** 2.37 Å.
- **Ligando de Control:** Ceftarolina (Residuo en PDB: **`AI8`**, Cadena A, número de residuo `1403`).
- **Componentes Eliminados del Receptor:**
  - **Solvente:** 602 moléculas de agua (`HOH`).
  - **Iones de Cristalización:** 4 iones de Cadmio (`CD`) y 2 iones de Cloro (`CL`).
- **Resultados:**
  - `data/receptors/PBP2a_3ZG0_clean.pdb`: Contiene únicamente la proteína (10,254 registros `ATOM` y `TER`).
  - `data/ligands/Ceftaroline_3D.pdb`: Contiene exclusivamente los 78 átomos de Ceftarolina con sus registros de conectividad molecular (`CONECT`).

### B. GyrB (Origen: PDB `3TTZ`)
- **Resolución:** 1.63 Å.
- **Ligando de Control:** Pirrolamida 07N (Residuo en PDB: **`07N`**, Cadena A y B).
- **Componentes Eliminados del Receptor:**
  - **Solvente:** Moléculas de agua (`HOH`).
  - **Iones:** 2 iones de Magnesio (`MG`) libres de la fase móvil de cristalización.
- **Resultados:**
  - `data/receptors/GyrB_3TTZ_clean.pdb`: Contiene la proteína de alta resolución limpia.
  - `data/ligands/07N_GyrB_3D.pdb`: Contiene los 52 átomos del inhibidor experimental.

### C. MurG (Origen: AlphaFold `AF-Q6GGZ0-F1`)
- **Resolución / Score:** pLDDT 93.56.
- **Características:** Al ser una estructura predicha computacionalmente (AlphaFold v6), el archivo PDB de origen **no contiene solvente, iones ni cofactores**.
- **Acción:** No requirió purificación algorítmica. Se posicionó directamente en `data/receptors/MurG_AF-Q6GGZ0-F1.pdb`.

---

## 4. Importancia Crítica en el Cribado Virtual

El docking molecular busca predecir el modo de unión termodinámicamente más estable de un ligando en el bolsillo activo de un receptor. La purificación previa implementada previene los siguientes errores catastróficos:

1. **Traslape de Átomos (Clashing):** Si el ligando co-cristalizado permaneciera en el receptor, el algoritmo de Vina detectaría colisión física insalvable y las energías de unión ($\Delta G$) calculadas para los avocadenoides o la afzelina serían erróneamente de $+10,000 \text{ kcal/mol}$ o simplemente no encontraría poses viables.
2. **Obstrucción por Solvente:** Las moléculas de agua cristalinas ocupan volumen físico en el sitio activo. En un sistema real de docking rígido o semi-flexible, impedirían que tu ligando explore las interacciones en el fondo del bolsillo (como los puentes de hidrógeno con residuos clave de unión).
3. **Redocking de Validación:** Al aislar el ligando en `Ceftaroline_3D.pdb` y `07N_GyrB_3D.pdb`, ahora puedes realizar el **redocking de validación** (re-acoplar el control en su propia proteína limpia) para medir el RMSD entre la pose experimental de rayos X y la predicha por AutoDock Vina. Un RMSD $\leq 2.0\text{ \AA}$ validará tu protocolo de docking antes de ensayar tus metabolitos.
