# Auditoría: Grid Box y Validación del Sitio Activo de MurG

> **Estado:** Pendiente de corrección antes de ejecutar cualquier corrida de docking.  
> **Afecta:** `scripts/calcular_centro.py` y todos los archivos `conf_*.txt` en `docking/vina_input/`.

---

## 1. Problemas identificados en `calcular_centro.py`

### 1.1 El script no filtra conformaciones alternativas (ALTLOC)

El parser actual incluye **todas** las líneas `ATOM` y `HETATM` sin discriminar por el campo ALTLOC (columna 17 del estándar PDB). Si el archivo PDB del ligando co-cristalizado tiene dos conformaciones alternativas (A y B), el centroide calculado será el promedio de ambas nubes de átomos y **no corresponderá a ninguna de las dos poses reales**. Antes de calcular el centro, se debe filtrar para retener solo la conformación `A` (o la que tenga mayor ocupancia).

**Corrección requerida en `find_coordinates()`:**

```python
altloc = line[16]  # columna 17, índice 16
if altloc not in (' ', 'A', ''):  # descartar conformaciones B, C, etc.
    continue
```

### 1.2 No se valida que el archivo contenga un único ligando

Si el PDB de entrada tiene moléculas de solvente residual (HOH, SO4, EDO, etc.) además del ligando de interés, el centroide se desplaza. El script debe filtrar por nombre de residuo del ligando objetivo o, como mínimo, excluir explícitamente registros de agua y iones comunes.

**Corrección requerida:**

```python
EXCLUIR_RESIDUOS = {'HOH', 'WAT', 'SO4', 'EDO', 'GOL', 'PEG', 'NA', 'CL', 'MG', 'ZN'}
res_name = line[17:20].strip()
if res_name in EXCLUIR_RESIDUOS:
    continue
```

### 1.3 El centro calculado no es "centro de masa"

El promedio simple de coordenadas es el **centroide geométrico**, no el centro de masa ponderado por masa atómica. Esto es aceptable para AutoDock Vina, pero debe reportarse correctamente en el manuscrito como "centro geométrico del ligando co-cristalizado". No usar el término "centro de masa" en la sección de métodos.

---

## 2. El tamaño de caja de 22 Å es insuficiente para ligandos grandes

### 2.1 Fundamento del problema

AutoDock Vina requiere que el ligando quepa **completamente** dentro de la caja en cualquier orientación posible. La dimensión mínima necesaria es:

```
tamaño_caja = dimensión_máxima_del_ligando + 2 × margen
```

El margen estándar para exploración libre de orientaciones es **8–10 Å por lado**. Los avocadenoides (Avocadenofuran, Avocadyne acetate, Avocadene acetate) tienen 14 enlaces rotables y en conformación extendida pueden superar los 18–20 Å de longitud. Con una caja de 22 Å y un ligando de 20 Å extendido, quedan solo **1 Å de margen** — las poses que requieran rotaciones libres serán truncadas silenciosamente sin generar error.

### 2.2 Corrección: función `calcular_tamaño_caja()` a agregar en el script

Agregar la siguiente función en `calcular_centro.py`:

```python
def calcular_tamaño_caja(path_file, margen=10.0, minimo=22.0):
    """
    Calcula el tamaño mínimo de caja cúbica para contener el ligando
    en cualquier orientación, con el margen especificado.
    
    Args:
        path_file: ruta al PDB del ligando
        margen:    espacio adicional por lado en Å (default: 10.0)
        minimo:    tamaño mínimo absoluto en Å (default: 22.0)
    Returns:
        float: tamaño de caja en Å (misma dimensión para X, Y, Z)
    """
    coords = _extraer_coords_limpias(path_file)  # usa la misma lógica filtrada
    if coords is None:
        return minimo
    dims = coords.max(axis=0) - coords.min(axis=0)
    tamaño_necesario = dims.max() + 2 * margen
    return round(max(tamaño_necesario, minimo), 1)
```

### 2.3 Tamaños de caja esperados por diana (estimación)

| Diana   | Ligando control | Ligandos grandes (avocadenoides) | Caja recomendada |
|---------|----------------|----------------------------------|------------------|
| PBP2a   | Ceftarolina     | Avocadenofuran / acetatos        | ≥ 30 Å           |
| GyrB    | 07N             | Avocadenofuran / acetatos        | ≥ 30 Å           |
| MurG    | Quercetina      | Avocadenofuran / acetatos        | ≥ 30 Å           |

> **Nota:** Calcular el tamaño real ejecutando `calcular_tamaño_caja()` sobre cada ligando SDF/PDB antes de fijar las cajas definitivas. Los valores de la tabla son estimaciones conservadoras. La caja debe ser **la misma para todos los ligandos en una misma diana** — se usa el valor máximo entre todos los ligandos del conjunto.

---

## 3. Validación del sitio activo de MurG (AlphaFold AF-Q6GGZ0)

MurG no tiene ligando co-cristalizado porque proviene de un modelo AlphaFold. Las coordenadas del sitio activo deben derivarse por **superposición estructural con una homóloga cristalográfica**, no por predicción de bolsillos como método primario.

### 3.1 Estructura de referencia recomendada

Usar **PDB 1F0K** — MurG de *E. coli* con UDP-GlcNAc co-cristalizado, resolución 1.90 Å. Esta es la estructura más citada para este propósito y está referenciada en `context.md`.

| Parámetro | Valor |
|-----------|-------|
| PDB ID    | 1F0K  |
| Organismo | *E. coli* |
| Ligando co-cristalizado | UDP-GlcNAc (sustrato natural) |
| Resolución | 1.90 Å |
| Identidad de secuencia con Q6GGZ0 | ~35–40% (suficiente para alineamiento estructural) |

### 3.2 Protocolo de superposición en PyMOL (paso a paso)

Ejecutar los siguientes comandos en la consola de PyMOL **en el orden indicado**:

```python
# Paso 1: cargar las dos estructuras
fetch 1F0K                                    # MurG E. coli con UDP-GlcNAc
load data/receptors/MurG_AF-Q6GGZ0-F1.pdb, MurG_MRSA   # modelo AlphaFold

# Paso 2: alineamiento global (secuencia + estructura)
align MurG_MRSA, 1F0K
# → registrar el RMSD global reportado. Criterio de aceptación: < 3.0 Å

# Paso 3: identificar residuos del sitio activo en 1F0K
# Los residuos que contactan UDP-GlcNAc a < 4 Å en 1F0K son:
# G14, T15, R169, E269, H302, N317 (numeración 1F0K)
select sitio_activo_1F0K, 1F0K and resi 14+15+169+269+302+317

# Paso 4: alineamiento LOCAL restringido al sitio activo
pair_fit MurG_MRSA and resi 14+15+169+269+302+317, sitio_activo_1F0K
# → registrar el RMSD LOCAL. Criterio de aceptación: < 3.5 Å
# Si RMSD local > 3.5 Å → las coordenadas trasladadas NO son confiables

# Paso 5: visualizar el ligando de referencia en el contexto de MurG_MRSA
select UDP_GlcNAc, 1F0K and resn UDP
show sticks, UDP_GlcNAc
color yellow, UDP_GlcNAc

# Paso 6: extraer las coordenadas del centroide del ligando
python
import pymol
model = cmd.get_model('UDP_GlcNAc')
coords = [atom.coord for atom in model.atom]
centro = [sum(c[i] for c in coords)/len(coords) for i in range(3)]
print(f'Centro MurG: X={centro[0]:.3f}  Y={centro[1]:.3f}  Z={centro[2]:.3f}')
python end
```

### 3.3 Criterios de validación del alineamiento

| Métrica | Criterio de aceptación | Acción si falla |
|---------|------------------------|------------------|
| RMSD global (align) | < 3.0 Å | Revisar si hay dominios móviles; intentar `cealign` |
| RMSD local sitio activo | < 3.5 Å | Coordenadas no trasladables; usar fpocket como alternativa |
| Residuos catalíticos alineados | ≥ 4 de 6 coinciden estructuralmente | Validar manualmente uno a uno |

### 3.4 Si el RMSD local supera 3.5 Å: protocolo alternativo con fpocket

Solo si el alineamiento local falla, ejecutar fpocket directamente sobre el modelo AlphaFold:

```bash
fpocket -f data/receptors/MurG_AF-Q6GGZ0-F1.pdb
```

fpocket genera un ranking de bolsillos por `druggability score`. Seleccionar el bolsillo con mayor score **que coincida topológicamente con el dominio de unión a UDP** (cara N-terminal del barrel TIM). **No seleccionar el bolsillo de mayor score si está en una región de pLDDT < 70** — verificar contra el gráfico de pLDDT por residuo antes de decidir.

### 3.5 Verificación del pLDDT en la región del bolsillo

Independientemente del método de localización del sitio activo, se debe generar el perfil de pLDDT por residuo en la zona del bolsillo. El valor de pLDDT 93.56 reportado es el **promedio global** del modelo; la confianza local en los loops del sitio de unión puede ser significativamente menor.

```python
# Script para extraer pLDDT (almacenado en columna B-factor del AlphaFold PDB)
import numpy as np

residues_plddt = {}
with open('data/receptors/MurG_AF-Q6GGZ0-F1.pdb') as f:
    for line in f:
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':  # solo Cα
            res_num = int(line[22:26].strip())
            plddt = float(line[60:66].strip())
            residues_plddt[res_num] = plddt

# Reportar los residuos del sitio activo (tras alineamiento con 1F0K)
# Sustituir los números con la numeración equivalente en Q6GGZ0
residues_sitio = [14, 15, 169, 269, 302, 317]  # actualizar tras alineamiento
for r in residues_sitio:
    if r in residues_plddt:
        print(f'Residuo {r}: pLDDT = {residues_plddt[r]:.1f}')
```

**Criterio:** Si algún residuo catalítico tiene pLDDT < 70, esta limitación debe reportarse explícitamente en la sección de métodos del manuscrito.

---

## 4. Resumen de acciones antes de correr Vina

1. **Corregir `calcular_centro.py`**: agregar filtro ALTLOC, filtro de solvente/iones, y función `calcular_tamaño_caja()`.
2. **Calcular tamaño de caja real**: ejecutar la función sobre los 8 ligandos (5 compuestos + 3 controles) y usar el valor máximo por diana.
3. **Ejecutar protocolo PyMOL para MurG**: alineamiento 1F0K → AF-Q6GGZ0, registrar RMSD global y local.
4. **Verificar pLDDT en residuos del sitio activo de MurG** con el script incluido arriba.
5. **Actualizar `conf_murg.txt`** con las coordenadas obtenidas del paso 3 y el tamaño de caja del paso 2.
6. **Ningún archivo en `docking/`** debe existir hasta completar los pasos 1–5 y registrar los RMSDs en este documento.
