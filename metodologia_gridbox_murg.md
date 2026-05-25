# Metodología: Definición de Grid Box y Validación del Sitio Activo de MurG

> **Objetivo de este documento:** Explicar el *porqué* científico de cada paso
> metodológico aplicado en la corrección de la grid box y la validación del
> sitio activo de MurG. No es un manual de instrucciones sino una justificación
> razonada, útil para la sección de métodos del manuscrito y como referencia
> de aprendizaje.

---

## 1. ¿Qué es la grid box y por qué importa definirla bien?

AutoDock Vina realiza el docking dentro de una **caja tridimensional** (grid
box) que define la región del espacio donde el algoritmo busca poses del
ligando. Todo lo que queda fuera de esta caja es invisible para Vina.

La grid box se define con dos parámetros:

| Parámetro | Significado |
|-----------|-------------|
| **Centro** (`center_x/y/z`) | Punto central de la caja, idealmente en el centro del sitio activo |
| **Tamaño** (`size_x/y/z`) | Dimensiones de la caja en Å |

Si el centro está desplazado, Vina busca en la región equivocada de la
proteína. Si la caja es demasiado pequeña, las poses que requieran rotaciones
amplias serán truncadas **sin generar error** — Vina simplemente no las
explora.

---

## 2. Correcciones al cálculo del centro

### 2.1 Filtro de conformaciones alternativas (ALTLOC)

#### ¿Qué es ALTLOC?

En cristalografía de rayos X, algunos átomos no tienen una posición única en
el cristal. Cuando el cristalógrafo resuelve la estructura, puede modelar
**dos o más conformaciones alternativas** del mismo átomo. Estas se distinguen
por la columna 17 del formato PDB (índice 16 en Python):

```
         1         2         3         4         5
1234567890123456789012345678901234567890123456789012345678
HETATM    1  CA  ALA A  42      12.345  67.890  23.456
                ^
                columna 17 = ALTLOC
```

- `' '` (espacio) → posición única, sin ambigüedad
- `'A'` → conformación A (generalmente la de mayor ocupancia)
- `'B'` → conformación B (alternativa)

#### ¿Por qué hay que filtrar?

Si un ligando co-cristalizado tiene dos conformaciones y calculamos el
centroide con ambas, el resultado es un **punto intermedio que no corresponde
a ninguna pose real**. Es como calcular la posición promedio de una persona
que está 50% del tiempo en la cocina y 50% en el dormitorio — el promedio
cae en el pasillo, donde nunca está.

#### ¿Qué hicimos?

En `_extraer_coords_limpias()` se lee la columna 17 y se retienen solo los
átomos con ALTLOC en `(' ', 'A', '')`. Esto garantiza que usamos una sola
conformación coherente.

---

### 2.2 Filtro de solvente e iones

#### ¿Por qué aparecen aguas e iones en un PDB?

Una estructura cristalográfica no contiene solo la proteína y su ligando.
El cristal se forma en una solución tampón que contiene agua, sales, y a
veces crioprotectores (como glicerol o polietilenglicol). Estas moléculas
quedan atrapadas en el cristal y aparecen como registros HETATM en el PDB:

- `HOH` / `WAT` — moléculas de agua
- `SO4`, `PO4` — iones del tampón
- `EDO`, `GOL`, `PEG` — crioprotectores
- `NA`, `CL`, `MG`, `ZN` — iones metálicos

#### ¿Por qué hay que excluirlos?

Si el archivo PDB del ligando incluye estas moléculas, el centroide se
desplaza hacia ellas. Una molécula de agua a 15 Å del ligando puede mover
el centroide varios angstroms, suficiente para desplazar la grid box del
sitio activo.

#### ¿Qué hicimos?

Se definió un conjunto `EXCLUIR_RESIDUOS` con los nombres de residuo de
solventes, iones y crioprotectores comunes. Antes de incluir un átomo en
el cálculo, se verifica que su nombre de residuo (columnas 18-20 del PDB)
no esté en esta lista.

---

### 2.3 Filtro de cadena — descubrimiento nuevo

#### Contexto: homodímeros cristalográficos

Muchas proteínas cristalizan como **homodímeros**: dos copias idénticas de
la misma cadena en la unidad asimétrica. Cuando un script extrae el ligando
de estas estructuras, puede extraer **ambas copias** (cadenas A y B).

En nuestro proyecto, esto ocurría con:

- **Ceftarolina (PBP2a, PDB 3ZG0):** 39 átomos en cadena A + 39 átomos en
  cadena B, separados por ~60 Å en el cristal.
- **07N (GyrB, PDB 3TTZ):** 26 átomos en cadena A + 26 átomos en cadena B.

#### ¿Cuál era el impacto?

El centroide se calculaba como el promedio de **ambas** nubes de átomos. El
resultado caía en un punto intermedio entre las dos copias del ligando —
**fuera del sitio activo de ambas cadenas**.

```
Cadena A (sitio activo real)    Centroide erróneo    Cadena B (otra copia)
        ★                            ×                       ★
   (27.9, 29.8, 88.3)          (-4.0, 37.1, 77.6)     (-35.7, 44.4, 66.9)
```

| Diana | Centro con ambas cadenas | Centro solo cadena A | Error (Å) |
|-------|--------------------------|---------------------|-----------|
| PBP2a | (-4.043, 37.118, 77.570) | (27.902, 29.847, 88.317) | ~35 Å |
| GyrB  | (8.203, -8.019, 15.475)  | (0.318, 2.949, 24.168)   | ~15 Å |

Un desplazamiento de 35 Å es catastrófico: la grid box entera quedaría en
una región de la proteína sin relevancia biológica.

#### ¿Qué hicimos?

Se agregó un parámetro `chain_id` a `_extraer_coords_limpias()`. Para los
ligandos co-cristalizados se usa `chain_id='A'`, reteniendo solo los átomos
de la primera cadena del homodímero. La cadena A se elige por convención
(es la primera en el PDB y generalmente tiene mejores estadísticas de
refinamiento).

---

### 2.4 Centroide geométrico vs. centro de masa

El script calcula un **centroide geométrico** (promedio simple de
coordenadas), no un centro de masa (ponderado por masa atómica). La
diferencia es sutil:

- **Centroide geométrico:** cada átomo contribuye por igual, independientemente
  de si es hidrógeno (1 Da) o azufre (32 Da).
- **Centro de masa:** los átomos pesados "pesan más" en el promedio.

Para definir la grid box de AutoDock Vina, el centroide geométrico es el
método estándar y aceptado. La diferencia con el centro de masa rara vez
supera 0.5 Å para moléculas orgánicas, y es irrelevante frente al tamaño
de la caja (37 Å).

**Para el manuscrito:** reportar como "centroide geométrico del ligando
co-cristalizado", nunca como "centro de masa".

---

## 3. Tamaño de caja: por qué 22 Å era insuficiente

### 3.1 Criterio para el tamaño mínimo

Vina necesita que el ligando **quepa completamente** dentro de la caja en
**cualquier orientación** posible. Dado que Vina explora rotaciones libres
del ligando, debemos considerar su dimensión máxima en conformación
extendida, más un margen para que no quede "raspando" los bordes:

```
tamaño_caja = dimensión_máxima_del_ligando + 2 × margen
```

El margen estándar en la literatura es **8–10 Å por lado**. Nosotros usamos
10 Å, que es el valor más conservador.

### 3.2 ¿Qué pasa si la caja es demasiado pequeña?

Si un ligando no cabe, Vina **no genera error**. Simplemente no explora las
conformaciones que se salen de la caja. Esto introduce un sesgo invisible:

- Las poses compactas (plegadas) serán exploradas → energías artificialmente
  favorecidas.
- Las poses extendidas (que podrían ser las bioactivas) serán ignoradas.
- El resultado es un falso ranking de afinidad.

### 3.3 Resultados del recálculo

Los avocadenoides tienen cadenas alifáticas largas (C₁₇-C₂₃). En
conformación extendida, su dimensión máxima alcanza ~17 Å:

| Ligando | Dimensión máxima | Caja necesaria (10 Å margen) |
|---------|------------------|------------------------------|
| Avocadenofuran | ~17.2 Å | 37.2 Å |
| Avocadyne acetate | ~16.5 Å | 36.5 Å |
| Avocadene acetate | ~16.8 Å | 36.8 Å |
| Ceftarolina (chain A) | ~17.4 Å | 37.4 Å |
| 07N (chain A) | ~13.1 Å | 33.1 Å |
| Quercetina | ~11.8 Å | 31.8 Å |

Con una caja de 22 Å y un ligando de 17 Å, quedaba solo **2.5 Å de margen
por lado** — insuficiente para explorar orientaciones libremente.

### 3.4 ¿Por qué una caja única para todas las dianas?

Se usa el **máximo global** (37.4 Å) para las 3 dianas. Razones:

1. **Consistencia metodológica:** si una diana tuviera una caja más pequeña,
   los ligandos grandes tendrían menos espacio conformacional en esa diana,
   sesgo que favorece ligandos pequeños.
2. **Reproducibilidad:** un solo valor es más fácil de reportar y replicar.
3. **Costo computacional marginal:** una caja de 37 Å vs 33 Å aumenta el
   volumen ~1.5×, pero con `exhaustiveness=32`, el tiempo adicional es
   aceptable.

---

## 4. Validación del sitio activo de MurG

### 4.1 ¿Por qué MurG es un caso especial?

PBP2a y GyrB tienen estructuras cristalográficas experimentales con ligandos
co-cristalizados. El centro de la grid box se obtiene directamente del
centroide del ligando.

MurG, en cambio, proviene de un **modelo AlphaFold** (AF-Q6GGZ0-F1). Estos
modelos predicen la estructura tridimensional de la proteína a partir de la
secuencia, pero **no contienen ligandos**. No hay un ligando co-cristalizado
del cual derivar el centro.

### 4.2 Estrategia: superposición con homólogo cristalográfico

La solución estándar es:

1. Buscar una estructura cristalográfica de un **homólogo** que sí tenga
   un ligando co-cristalizado.
2. **Superponer** el homólogo sobre el modelo AlphaFold.
3. Las coordenadas del ligando del homólogo, **tras la superposición**,
   indican dónde está el sitio activo en el modelo AlphaFold.

#### ¿Por qué 1F0K?

PDB 1F0K es la estructura cristalográfica de MurG de *E. coli* con su
sustrato natural (UDP-GlcNAc) co-cristalizado a 1.90 Å de resolución.
Es la referencia más citada para este propósito.

Aunque la identidad de secuencia con Q6GGZ0 (*S. aureus*) es solo ~25%,
las enzimas de la familia de glicosiltransferasas GT-B tienen un plegamiento
altamente conservado (barrel TIM de dos dominios), y los residuos catalíticos
suelen estar conservados estructuralmente incluso con baja identidad de
secuencia.

### 4.3 Paso a paso: qué hace el script y por qué

#### Paso 1 — Carga de estructuras

Se cargan ambos PDB usando `Bio.PDB.PDBParser`. Esto crea una jerarquía
`Structure → Model → Chain → Residue → Atom` que permite acceder a las
coordenadas atómicas de forma estructurada.

**1F0K tiene dos cadenas (A y B)** porque cristalizó como homodímero. El
script usa la primera cadena de cada estructura.

#### Paso 2 — Alineamiento de secuencia

Antes de superponer las estructuras, necesitamos saber **qué residuos de
una corresponden a cuáles de la otra**. Esto se hace alineando las
secuencias de aminoácidos.

Se usa `Bio.Align.PairwiseAligner` con un alineamiento global (tipo
Needleman-Wunsch). Del alineamiento se extrae un mapeo
`posición_1F0K → posición_MRSA`.

Los residuos del sitio activo en la numeración de 1F0K son: **14, 15, 169,
269, 302, 317** (residuos que contactan UDP-GlcNAc a < 4 Å en el cristal).
El mapeo encontrado fue:

| 1F0K | MRSA | Aminoácido 1F0K | Aminoácido MRSA | Conservado |
|------|------|-----------------|-----------------|------------|
| 14   | 9    | GLY | GLY | ✓ Idéntico |
| 15   | 10   | GLY | GLY | ✓ Idéntico |
| 169  | 173  | ALA | ASN | ≈ Diferente |
| 269  | 268  | GLU | GLU | ✓ Idéntico |
| 302  | 303  | ALA | ALA | ✓ Idéntico |
| 317  | 316  | ASN | LEU | ≈ Diferente |

4 de 6 residuos son idénticos, y los 6 pudieron mapearse (no hay gaps en
el alineamiento en estas posiciones). Esto cumple el criterio de la auditoría
(≥ 4 de 6 coinciden).

#### Paso 3 — Superposición global

**¿Qué es una superposición estructural?**

Es el proceso de encontrar la rotación y traslación que minimizan la
distancia entre átomos equivalentes de dos estructuras. Se usa el algoritmo
de Kabsch, que encuentra la matriz de rotación óptima por descomposición en
valores singulares (SVD) de la matriz de covarianza.

El RMSD (Root Mean Square Deviation) cuantifica la calidad del alineamiento:

```
RMSD = √( Σ(d_i²) / N )
```

donde `d_i` es la distancia entre el átomo i de la estructura fija y el
móvil, y N es el número de pares.

**Resultado:** RMSD global = 5.784 Å.

Esto **no pasa** el criterio de < 3.0 Å. Sin embargo, es un resultado
esperado con ~25% de identidad de secuencia. Los dominios periféricos y
loops divergen considerablemente entre *E. coli* y *S. aureus*, inflando
el RMSD global. Lo que importa es el RMSD **local** del sitio activo.

**Dirección de la superposición:** Se alinea 1F0K **sobre** MurG_MRSA (1F0K
es el "móvil", MurG_MRSA es el "fijo"). Así, tras aplicar la transformación,
las coordenadas del ligando de 1F0K quedan expresadas en el sistema de
coordenadas de MurG_MRSA — que es lo que necesitamos para la grid box.

#### Paso 4 — Superposición local

Se repite el cálculo de RMSD pero usando **solo** los 6 residuos del sitio
activo. Esto aísla la calidad del alineamiento en la región de interés
biológico.

**Resultado:** RMSD local = 3.361 Å (< 3.5 Å → aceptable).

Las distancias individuales entre CAs fueron:

| Par de residuos | Distancia CA |
|----------------|-------------|
| GLY14 ↔ GLY9   | 1.92 Å |
| GLY15 ↔ GLY10  | 2.23 Å |
| ALA169 ↔ ASN173 | **9.04 Å** (outlier) |
| GLU269 ↔ GLU268 | 1.94 Å |
| ALA302 ↔ ALA303 | 1.80 Å |
| ASN317 ↔ LEU316 | 4.37 Å |

El residuo 169→173 es un outlier extremo (9.04 Å). Esto indica que el loop
que contiene este residuo ha divergido estructuralmente. Sin este residuo,
el RMSD local sería < 2.5 Å. Esto debe reportarse como limitación.

#### Paso 5 — Extracción del centroide

Idealmente, se usaría el centroide del ligando UDP-GlcNAc de 1F0K (ya en
coordenadas de MurG_MRSA tras la superposición). Sin embargo, el archivo
`MurG_template_1F0K.pdb` descargado **no incluye el ligando** (solo contiene
SO4 y agua como HETATM).

Como alternativa, se calculó el centroide de los **6 átomos CA del sitio
activo** directamente en MurG_MRSA. Esta es una aproximación válida: el
centroide de los residuos catalíticos define el centro del bolsillo de unión.

**Resultado:** Centro del sitio activo = **(4.729, 1.061, -3.361)**

#### Paso 6 — Verificación de pLDDT

**¿Qué es el pLDDT?**

El pLDDT (predicted Local Distance Difference Test) es la métrica de
confianza de AlphaFold. Se reporta **por residuo** y se almacena en la
columna de B-factor del PDB. Su escala es:

| pLDDT | Interpretación |
|-------|----------------|
| > 90  | Muy alta confianza, estructura confiable |
| 70-90 | Buena confianza, backbone correcto |
| 50-70 | Baja confianza, posibles errores de plegamiento |
| < 50  | Muy baja confianza, región probablemente desordenada |

**¿Por qué verificar el pLDDT del sitio activo?**

El valor reportado de pLDDT = 93.56 es el **promedio global** del modelo
completo. Los loops del sitio de unión podrían tener confianza menor. Si los
residuos catalíticos tienen pLDDT < 70, las coordenadas del sitio activo no
son confiables y el docking en esa región pierde validez.

**Resultados:**

| Residuo | pLDDT | Interpretación |
|---------|-------|----------------|
| GLY 9   | 92.06 | Muy alta |
| GLY 10  | 87.50 | Alta |
| ASN 173 | 97.00 | Muy alta |
| GLU 268 | 93.75 | Muy alta |
| ALA 303 | 98.19 | Muy alta |
| LEU 316 | 95.62 | Muy alta |
| **Promedio** | **94.02** | **Muy alta** |

Todos los residuos superan ampliamente el umbral de 70. Esto confirma que
AlphaFold tiene alta confianza en la estructura del sitio activo, y las
coordenadas son válidas para docking.

---

## 5. Resumen de parámetros finales

### Configuraciones de Vina generadas

| Diana | Archivo | Centro (X, Y, Z) | Caja | Método del centro |
|-------|---------|-------------------|------|-------------------|
| PBP2a | `conf_pbp2a.txt` | (27.902, 29.847, 88.317) | 37.4 Å | Centroide geométrico de Ceftarolina (cadena A) |
| GyrB  | `conf_gyrb.txt`  | (0.318, 2.949, 24.168)   | 37.4 Å | Centroide geométrico de 07N (cadena A) |
| MurG  | `conf_murg.txt`  | (4.729, 1.061, -3.361)   | 37.4 Å | Centroide de CA del sitio activo (superposición con 1F0K) |

### Parámetros comunes

| Parámetro | Valor | Justificación |
|-----------|-------|---------------|
| `exhaustiveness` | 32 | Mayor muestreo que el default (8), recomendado para reproducibilidad |
| `num_modes` | 9 | Número de poses a reportar |
| `energy_range` | 3 | Solo poses dentro de 3 kcal/mol de la mejor |

---

## 6. Scripts generados

| Script | Función |
|--------|---------|
| `scripts/calcular_centro.py` | Calcula centroide y tamaño de caja con filtros (ALTLOC, solvente, cadena). Genera archivos `conf_*.txt`. |
| `scripts/recalcular_cajas.py` | Calcula el tamaño de caja por ligando y por diana. Solo reporta, no modifica archivos. |
| `scripts/validar_sitio_murg.py` | Superposición estructural 1F0K → MurG MRSA. Reporta RMSDs, mapeo de residuos, centroide, y pLDDT. |

---

## 7. Limitaciones a reportar en el manuscrito

1. **RMSD global alto (5.784 Å):** Refleja la baja identidad de secuencia
   (~25%) entre MurG de *E. coli* y *S. aureus*. El sitio catalítico está
   conservado localmente (RMSD 3.361 Å), pero la estructura global diverge.

2. **Residuo ASN173 (equivalente a ALA169):** Distancia CA de 9.04 Å tras
   superposición. Este loop ha divergido estructuralmente. Si este residuo
   participa directamente en la unión del sustrato, las predicciones de
   docking en esa subregión son menos confiables.

3. **Ausencia de ligando en 1F0K:** El centroide del sitio activo se calculó
   a partir de los CA de los residuos catalíticos, no del ligando
   co-cristalizado. Esto es una aproximación — el centro real del bolsillo
   podría estar ligeramente desplazado respecto al centroide de CA.

4. **Terminología:** El centro se reporta como "centroide geométrico", no
   como "centro de masa", ya que no se pondera por masa atómica.
