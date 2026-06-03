#set text(
  size: 11pt,
  lang: "es"
)

#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 3cm),
  header: align(right, text(size: 8pt, fill: gray)[
    Metodología Computacional y Validación - MRSA
  ]),
  footer: [
    #set text(size: 8pt, fill: gray)
    #grid(
      columns: (1fr, 1fr),
      align(left)[Sección de Metodología y Validación],
      align(right)[Pág. #context { counter(page).display() }]
    )
  ]
)

#set par(
  justify: true,
  leading: 0.65em
)

#show heading: set text(fill: rgb("#1a3a5f"))
#show heading.where(level: 1): it => {
  v(1.5em, weak: true)
  it
  v(0.8em, weak: true)
}
#show heading.where(level: 2): it => {
  v(1.2em, weak: true)
  it
  v(0.6em, weak: true)
}

#align(center)[
  #text(size: 18pt, weight: "bold", fill: rgb("#1a3a5f"))[
    Metodología Computacional y Validación del Sitio Activo en Dianas de _Staphylococcus aureus_ Resistente a Meticilina (MRSA)
  ]
  #v(1em)
  #text(size: 11pt, style: "italic", fill: gray)[
    Sección de Metodología y Validación para Publicación Científica
  ]
]

#v(2em)

= 1. Introducción y Objetivos del Estudio Computacional

El aumento global de cepas de _Staphylococcus aureus_ resistente a meticilina (MRSA, por sus siglas en inglés) representa una de las amenazas más críticas para la salud pública contemporánea. MRSA ha desarrollado resistencia a prácticamente todos los antibióticos $beta$-lactámicos tradicionales mediante la adquisición del gen _mecA_, el cual codifica para la proteína de unión a penicilina 2a (PBP2a), una transpeptidasa con baja afinidad por estos fármacos. Adicionalmente, enzimas esenciales para la replicación del ADN, como la ADN Girasa (específicamente la subunidad GyrB), y para la biosíntesis de la pared celular, como la glicosiltransferasa MurG, constituyen dianas terapéuticas de alto valor para el desarrollo de nuevos agentes antimicrobianos que evadan los mecanismos de resistencia existentes.

En este contexto, el presente estudio computacional tiene como objetivo evaluar el potencial inhibitorio y el mecanismo de interacción a nivel molecular de metabolitos secundarios bioactivos (avocadenoides: avocadenofurano, acetato de avocadina y acetato de avocadeno) frente a tres dianas moleculares clave de _S. aureus_: PBP2a, GyrB y MurG. Mediante un protocolo riguroso de acoplamiento molecular (_molecular docking_) y simulaciones de dinámica molecular (MD), se busca caracterizar la afinidad termodinámica, la estabilidad estructural de los complejos ligando-receptor y los perfiles de interacción a corto y largo plazo. Este enfoque _in silico_ permite dilucidar si los avocadenoides representan candidatos viables para el desarrollo de nuevas terapias dirigidas a mitigar la resistencia de MRSA.

= 2. Fase 1: Preparación del Receptor y de los Ligandos

La calidad y fiabilidad de los resultados de acoplamiento molecular dependen de manera crítica de la preparación geométrica y química de los modelos estructurales de partida. Por ello, se implementó un flujo de trabajo riguroso y automatizado para la adecuación de receptores y ligandos.

== 2.1 Generación de Conformómeros 3D y Minimización de Ligandos
Las estructuras bidimensionales de los compuestos de prueba (avocadenoides) y de los controles de referencia (Ceftarolina para PBP2a, Pirrolamida 07N para GyrB, y Quercetina para MurG) se convirtieron a coordenadas tridimensionales mediante el kit de herramientas de quimioinformática RDKit. Para la generación de conformómeros 3D se empleó el algoritmo *ETKDGv3* (_Experimental-Torsion Knowledge Distance Geometry, versión 3_), el cual combina las restricciones de distancia geométrica tradicionales con reglas empíricas de ángulos de torsión derivadas de estructuras cristalográficas de ligandos de alta resolución. Este método asegura que las conformaciones iniciales sean físicamente realistas y bioactivamente relevantes.

A todos los ligandos se les añadieron hidrógenos explícitos para representar adecuadamente sus estados de protonación a pH fisiológico (7.4). Posteriormente, se llevó a cabo una minimización de energía utilizando el campo de fuerza de valencia molecular *MMFF94* (_Merck Molecular Force Field 94_) con un umbral de convergencia de gradiente de energía de $10^(-6) " kcal mol"^(-1) "Å"^(-1)$. Las estructuras minimizadas se exportaron en formato `.pdbqt` empleando herramientas de preparación de ligandos para preservar las cargas atómicas parciales y la asignación de enlaces rotables.

== 2.2 Purificación Estructural de Receptores Cristalográficos
Las coordenadas estructurales de los receptores se obtuvieron del _Protein Data Bank_ (PDB):
- *PBP2a*: Resuelto mediante cristalografía de rayos X a una resolución de 2.37 Å (código PDB: *3ZG0*), co-cristalizado con Ceftarolina (código de residuo: `AI8`).
- *GyrB*: Resuelto por cristalografía de rayos X a una resolución de 1.63 Å (código PDB: *3TTZ*), co-cristalizado con el inhibidor experimental Pirrolamida 07N (código de residuo: `07N`).

El proceso de purificación de estas estructuras cristalográficas se diseñó para eliminar todos los componentes de la fase cristalina que no forman parte del marco estructural biológicamente activo, evitando así colisiones estéricas artificiales y sesgos en el cálculo de mapas de afinidad. El protocolo consistió en:
1. *Eliminación de solvente*: Se removieron todas las moléculas de agua de cristalización (`HOH` / `WAT`), dado que su inclusión en el acoplamiento rígido restringe indebidamente el volumen útil del bolsillo activo.
2. *Remoción de aditivos y amortiguadores*: Se excluyeron iones del tampón y crioprotectores como sulfatos ($"SO"_4^(2-)$), fosfatos ($"PO"_4^(3-)$), etilenglicol (`EDO`) y glicerol (`GOL`).
3. *Exclusión de metales no estructurales*: Se eliminaron los iones de cadmio ($"Cd"^(2+)$) y cloro ($"Cl"^(-)$) de la estructura 3ZG0, así como los iones de magnesio ($"Mg"^(2+)$) libres de la estructura 3TTZ, dado que representaban aditivos de la cristalización y no cofactores del sitio catalítico de unión.

== 2.3 Filtrado de Cadenas en Homodímeros Cristalográficos
Un hallazgo crítico durante la auditoría metodológica reveló que tanto PBP2a (3ZG0) como GyrB (3TTZ) cristalizan como homodímeros en la unidad asimétrica del cristal. En consecuencia, el archivo original contiene dos sitios de unión biológicamente equivalentes por cada proteína (Ceftarolina unida a las cadenas A y B en 3ZG0; 07N unido a las cadenas A y B en 3TTZ). 

Si la extracción del ligando de referencia y el posterior cálculo del centro de la caja se realizan sobre la estructura sin discriminar cadenas, el algoritmo calcula el centro geométrico promediando las coordenadas de ambos ligandos. Al estar separados por distancias macroscópicas en el cristal (~60 Å en PBP2a y ~15 Å en GyrB), el centroide intermedio resultante se ubica en un espacio vacío fuera de la proteína, lo que desplaza catastróficamente la caja de docking fuera de cualquier sitio activo real (un desplazamiento de ~35 Å para PBP2a y ~15 Å para GyrB).

Para corregir esta anomalía, se implementó un script de parseo estricto que retiene de manera exclusiva la *Cadena A* tanto del receptor como del ligando co-cristalizado. Esto asegura que la definición del bolsillo y el cálculo de la _grid box_ se centren de manera inequívoca en un único sitio activo biológicamente coherente.

= 3. Fase 2: Validación del Sitio Activo de MurG (AlphaFold AF-Q6GGZ0-F1)

A diferencia de PBP2a y GyrB, la proteína MurG de _S. aureus_ carece de una estructura cristalográfica experimental depositada en el PDB con un ligando co-cristalizado. En su lugar, se utilizó el modelo estructural predicho por *AlphaFold Database* (código de acceso: *AF-Q6GGZ0-F1*). Debido a la ausencia de un ligando de referencia en el modelo computacional, fue indispensable implementar un protocolo riguroso de validación y localización del sitio activo basado en homología estructural.

== 3.1 Alineamiento y Superposición Estructural con la Plantilla 1F0K
El sitio activo de MurG se cartografió mediante una superposición estructural con su homóloga cristalográfica de _Escherichia coli_ (código PDB: *1F0K*), la cual se encuentra resuelta a una resolución de 1.90 Å y co-cristalizada con su sustrato natural, la uridina difosfato _N_-acetilglucosamina (UDP-GlcNAc). A pesar de que la identidad de secuencia global entre la MurG de _E. coli_ y la de _S. aureus_ es moderada (~25%), el plegamiento tridimensional de la familia de glicosiltransferasas GT-B (compuesto por dos dominios tipo _TIM-barrel_ enfrentados) está altamente conservado.

El protocolo de alineamiento consistió en:
1. *Alineamiento de Secuencia*: Se realizó un alineamiento de secuencia global por el método de Needleman-Wunsch para establecer la correspondencia residuo a residuo. Los residuos del sitio de unión de UDP-GlcNAc en la plantilla 1F0K (definidos como aquellos a una distancia $< 4$ Å del sustrato co-cristalizado) se mapearon con sus equivalentes en la MurG de _S. aureus_ (AF-Q6GGZ0-F1). El mapeo resultante se detalla en la @residues-table.

#v(0.8em)
#align(center)[
  #figure(
    table(
      columns: (1.2fr, 1.25fr, 1.6fr, 1.6fr, 1.4fr, 1.35fr),
      align: center + horizon,
      stroke: (x, y) => if y == 0 { (bottom: 1.5pt + black) } else { 0.5pt + gray },
      fill: (x, y) => if y == 0 { rgb("#e6f0fa") } else { none },
      [*Residuo 1F0K*], [*Residuo MRSA*], [*Aminoácido 1F0K*], [*Aminoácido MRSA*], [*Distancia $C_alpha$ (Å)*], [*pLDDT (MRSA)*],
      [14], [9], [GLY], [GLY], [1.92], [92.06],
      [15], [10], [GLY], [GLY], [2.23], [87.50],
      [169], [173], [ALA], [ASN], [9.04], [97.00],
      [269], [268], [GLU], [GLU], [1.94], [93.75],
      [302], [303], [ALA], [ALA], [1.80], [98.19],
      [317], [316], [ASN], [LEU], [4.37], [95.62]
    ),
    caption: [Mapeo de residuos del sitio activo entre la plantilla 1F0K (E. coli) y el modelo AlphaFold AF-Q6GGZ0-F1 (S. aureus), incluyendo distancias interatómicas $C_alpha$ post-superposición y valores de confianza local pLDDT.]
  ) <residues-table>
]
#v(0.8em)

2. *Superposición Estructural Global*: Utilizando el algoritmo de Kabsch basado en descomposición en valores singulares (SVD) de la matriz de covarianza de coordenadas, se superpuso la plantilla 1F0K sobre el modelo de _S. aureus_. El alineamiento global reportó un *RMSD global de 5.784 Å*. Este valor, relativamente elevado, refleja la divergencia evolutiva y la flexibilidad intrínseca en los dominios y bucles periféricos de ambas enzimas, lo cual es normal para identidades de secuencia cercanas al 25%.
3. *Superposición Estructural Local*: Para evaluar si el sitio catalítico mantiene una arquitectura espacial idéntica a pesar de la divergencia global, se calculó el RMSD restringido únicamente a los átomos de carbono alfa ($C_alpha$) de los 6 residuos del sitio activo mapeados. El *RMSD local del sitio activo fue de 3.361 Å*, el cual cumple con éxito el criterio de aceptación metodológico preestablecido ($< 3.5$ Å).

Una inspección detallada de las distancias individuales de $C_alpha$ (véase la @residues-table) indica que cinco de los seis residuos presentan un acoplamiento espacial sobresaliente (distancias entre 1.80 Å y 4.37 Å). No obstante, se identificó al par *ALA169 ↔ ASN173* como un valor atípico (_outlier_) con una distancia de *9.04 Å*. Esto demuestra una divergencia conformacional o un desplazamiento local significativo del bucle que aloja este residuo en _S. aureus_. Esta observación se reporta como una limitación metodológica del modelo estático para el docking de esta subregión específica, la cual será evaluada bajo el régimen dinámico en las simulaciones de MD.

== 3.2 Evaluación de Confianza Local (pLDDT) del Modelo de MurG
Para asegurar la viabilidad del docking en un modelo computacional, es imperativo comprobar que la región del sitio activo posea una alta fiabilidad predictiva. El indicador de confianza local de AlphaFold, *pLDDT* (_predicted Local Distance Difference Test_), el cual fluctúa entre 0 y 100, se extrajo directamente de los factores B del archivo PDB.

A pesar de que el pLDDT global promedio de toda la proteína MurG de _S. aureus_ es de 93.56, se evaluó individualmente el pLDDT de cada residuo del sitio de unión para descartar la presencia de bucles desordenados de baja confianza (pLDDT $< 70$). Los resultados detallados en la @residues-table revelan una confianza local excepcional, con una media del sitio activo de *94.02* y todos los residuos individuales superando ampliamente el umbral de 70 (rango de 87.50 a 98.19). Esto valida inequívocamente la rigidez estructural y la precisión geométrica de la cadena principal en el sitio activo predicho, justificando su uso para estudios de cribado virtual.

== 3.3 Definición de las Coordenadas Finales de la Grid Box
Debido a la ausencia física del ligando UDP-GlcNAc en las coordenadas originales de la plantilla PDB de MurG, el centro del bolsillo catalítico se definió formalmente como el *centroide geométrico de los átomos $C_alpha$ de los seis residuos del sitio activo* mapeados de _S. aureus_. El cálculo arrojó las coordenadas tridimensionales de centrado: *(4.729, 1.061, -3.361)*, las cuales se establecieron de forma definitiva para la ejecución de AutoDock Vina.

= 4. Fase 3: Protocolo de Acoplamiento Molecular (Molecular Docking)

El acoplamiento molecular se llevó a cabo utilizando el motor de docking *AutoDock Vina versión 1.2.0*, ampliamente validado por su precisión y velocidad en el cribado virtual de bibliotecas de metabolitos.

== 4.1 Auditoría del Tamaño de la Grid Box: Justificación del Límite de 37.4 Å
Durante la fase de diseño experimental, se identificó que la dimensión inicial de la _grid box_ de 22 Å cúbicos resultaba gravemente deficiente para la evaluación de metabolitos de cadena alifática larga, como los avocadenoides. Los compuestos en estudio (avocadenofurano, acetato de avocadina y acetato de avocadeno) poseen una estructura altamente lineal y flexible, caracterizada por *14 enlaces rotables libres* en sus colas alifáticas.

En una conformación completamente extendida, los avocadenoides alcanzan una longitud máxima ($d_("max")$) de aproximadamente *16.5 a 17.2 Å*. El algoritmo de búsqueda de AutoDock Vina opera bajo la restricción estricta de que el ligando debe caber *en su totalidad* dentro del volumen delimitado por la caja en cualquier orientación tridimensional. Si el tamaño de la caja es menor que la longitud del ligando más un margen de amortiguación, el software trunca silenciosamente las poses extendidas sin emitir errores de ejecución. Esto introduce un sesgo metodológico severo: se favorecen artificialmente las conformaciones enrolladas de alta energía interna y se omiten las conformaciones extendidas de unión bioactiva, invalidando los rankings de energía libre de Gibbs ($Delta G$).

Para solucionar esta limitación, se implementó un algoritmo dinámico de cálculo de caja basado en la siguiente formulación matemática:

$ L_("box") = d_("max") + 2 times text("margen") $

Donde el margen representa la holgura espacial mínima requerida a cada lado para permitir giros libres de $360°$ del ligando sin rozar los límites de la caja. Se adoptó un margen conservador y estándar en la literatura de *10.0 Å* por lado. Aplicando esta lógica a los ligandos del conjunto de datos y los controles, se obtuvieron las siguientes necesidades individuales de caja:
- *Avocadenofurano*: $17.2 " Å" + 20 " Å" = 37.2 " Å"$
- *Acetato de avocadina*: $16.5 " Å" + 20 " Å" = 36.5 " Å"$
- *Acetato de avocadeno*: $16.8 " Å" + 20 " Å" = 36.8 " Å"$
- *Ceftarolina* (control de PBP2a, cadena A): $17.4 " Å" + 20 " Å" = 37.4 " Å"$
- *Pirrolamida 07N* (control de GyrB, cadena A): $13.1 " Å" + 20 " Å" = 33.1 " Å"$
- *Quercetina* (control de MurG): $11.8 " Å" + 20 " Å" = 31.8 " Å"$

Para mantener una rigurosa consistencia metodológica y evitar sesgos de volumen en la exploración conformacional entre diferentes sistemas, se fijó el tamaño máximo calculado de *37.4 Å cúbicos* de manera uniforme para todas las corridas de docking en las tres dianas. En la @gridbox-table se resumen los parámetros finales de centrado y dimensión del espacio de búsqueda.

#v(0.8em)
#align(center)[
  #figure(
    table(
      columns: (1fr, 2.2fr, 1.2fr, 2.6fr),
      align: center + horizon,
      stroke: (x, y) => if y == 0 { (bottom: 1.5pt + black) } else { 0.5pt + gray },
      fill: (x, y) => if y == 0 { rgb("#e6f0fa") } else { none },
      [*Diana*], [*Centro geométrico (X, Y, Z)*], [*Tamaño (Å)*], [*Origen de la Definición del Centro*],
      [PBP2a], [(27.902, 29.847, 88.317)], [37.4], [Centroide de Ceftarolina (Cadena A, PDB 3ZG0)],
      [GyrB], [(0.318, 2.949, 24.168)], [37.4], [Centroide de Pirrolamida 07N (Cadena A, PDB 3TTZ)],
      [MurG], [(4.729, 1.061, -3.361)], [37.4], [Centroide de $C_alpha$ del Sitio Activo (1F0K)]
    ),
    caption: [Parámetros definitivos de configuración de la grid box en AutoDock Vina 1.2.0 para cada una de las dianas biológicas de estudio.]
  ) <gridbox-table>
]
#v(0.8em)

== 4.2 Parámetros Fijos de Búsqueda Estocástica
Debido al aumento del volumen de búsqueda que representa una caja de 37.4 Å (aproximadamente un incremento de 5 veces en volumen comparado con una caja de 22 Å) y al elevado número de grados de libertad conformacionales de los avocadenoides, se ajustaron los parámetros heurísticos de AutoDock Vina para garantizar la reproducibilidad global y evitar quedar atrapados en mínimos locales de energía:
- *Exhaustividad del muestreo* (`exhaustiveness`): Fijado en *32* (cuatro veces superior al valor por defecto de 8). Este muestreo exhaustivo asegura que los algoritmos genéticos iteren el tiempo suficiente para cartografiar el espacio de fases tridimensional ampliado.
- *Número de modos de unión* (`num_modes`): Establecido en *9*, para retener un conjunto diverso de orientaciones y conformaciones bioactivas alternativas.
- *Rango energético* (`energy_range`): Fijado en *3.0 kcal/mol*, lo que descarta de forma automática poses de acoplamiento de baja estabilidad termodinámica y retiene únicamente aquellas dentro de una ventana termodinámicamente accesible desde la pose de menor energía libre de Gibbs ($Delta G$).

= 5. Fase 4: Plan de Simulación de Dinámica Molecular (MD)

Para superar las limitaciones inherentes del docking molecular rígido (como el desprecio de la flexibilidad inducida del receptor y los efectos del solvente explícito), se estructuró un plan riguroso para la ejecución de simulaciones de Dinámica Molecular (MD) en complejos de alta afinidad seleccionados. Las simulaciones se llevarán a cabo en entornos de computación de alto rendimiento acelerados por GPU empleando el motor *GROMACS (versión 2023 o superior)*.

== 5.1 Parámetros Físicos y Configuración del Sistema
El flujo de trabajo sistemático para la construcción y parametrización de los sistemas solutos se describe a continuación:
1. *Parametrización del Receptor*: Las estructuras proteicas limpias se modelarán bajo el campo de fuerza de todos los átomos *AMBER99SB-ILDN*, el cual incluye correcciones refinadas para los potenciales torsionales del esqueleto y cadenas laterales de aminoácidos, garantizando la estabilidad conformacional del receptor durante escalas de tiempo extendidas.
2. *Parametrización del Ligando*: Las topologías y parámetros de fuerza de los avocadenoides se construirán mediante el campo de fuerza general de Amber *GAFF2* (_General Amber Force Field, versión 2_). Los cargos atómicos parciales se calcularán bajo el esquema semi-empírico de cargas de enlace *AM1-BCC* empleando la herramienta `antechamber` a través de la interfaz de automatización Python `acpype`. Este nivel de teoría es el estándar de oro en dinámica molecular clásica para fármacos pequeños y metabolitos flexibles.
3. *Solvatación Explícita*: Los complejos ligando-receptor se insertarán en cajas de simulación cúbicas respetando una distancia amortiguadora mínima de *10.0 Å* (1.0 nm) entre el átomo de soluto más externo y las fronteras de la caja de simulación. Se empleará el modelo de agua explícita de tres puntos *TIP3P*.
4. *Neutralización e Ionicidad*: Para simular condiciones fisiológicas realistas y balancear cargas del sistema, se agregarán contraiones neutros de sodio ($"Na"^(+)$) y cloruro ($"Cl"^(-)$) hasta alcanzar una concentración salina final de *0.15 M*.
5. *Condiciones de Frontera*: Se aplicarán condiciones de frontera periódicas (PBC, _Periodic Boundary Conditions_) en las tres direcciones ortogonales del espacio para modelar un medio termodinámico continuo e infinito.

== 5.2 Protocolo de Relajación Termodinámica (Equilibración)
Antes de proceder con la recolección de trayectorias de producción, cada complejo solvatado se someterá a un protocolo de equilibración multietapa para eliminar tensiones estéricas e integrar el solvente:
- *Minimización de Energía*: Se ejecutará un algoritmo de descenso de gradiente más pronunciado (_steepest descent_) por un máximo de 50,000 pasos o hasta que la fuerza máxima del sistema sea inferior a $1000 " kJ mol"^(-1) "nm"^(-1)$, resolviendo cualquier traslape atómico inicial.
- *Equilibración NVT* (Temperatura Constante): Se llevará a cabo durante *100 ps* a una temperatura constante de *300 K*, utilizando el termostato de re-escalado de velocidad modificado (_V-rescale_) con un tiempo de acoplamiento térmico ($tau_t$) de 0.1 ps. Se aplicarán restricciones posicionales armónicas de $1000 " kJ mol"^(-1) "nm"^(-2)$ sobre todos los átomos pesados de la proteína y el ligando para estabilizar la fase de empaquetamiento del solvente.
- *Equilibración NPT* (Presión Constante): Se ejecutará durante *100 ps* a una presión de referencia constante de *1.0 bar* y temperatura de *300 K*, empleando el baróstato isotrópico C-rescale con un tiempo de acoplamiento de presión ($tau_p$) de 5.0 ps y una compresibilidad isotérmica de $4.5 times 10^(-5)$ a $4.5 times 10^(-6) " bar"^(-1)$ para garantizar la estabilidad termodinámica del sistema y mitigar oscilaciones extremas de presión.

== 5.3 Producción de Trayectorias y Parámetros Numéricos
Las simulaciones de producción libres se realizarán por un intervalo acumulado de *100 ns* por cada sistema de complejo seleccionado, prescindiendo de toda restricción posicional para permitir el ajuste flexible inducido y la exploración termodinámica del ligando. Los parámetros numéricos y de integración serán:
- *Paso de integración temporal*: 2.0 femtosegundos (fs).
- *Tratamiento de enlaces*: Todos los enlaces covalentes que involucren átomos de hidrógeno se restringirán rígidamente utilizando el algoritmo *LINCS* (_Linear Constraint Solver_), permitiendo el incremento del paso de integración a 2 fs sin inestabilidades numéricas.
- *Interacciones de Rango Largo*: Las fuerzas electrostáticas de largo alcance se evaluarán empleando el formalismo tridimensional *PME* (_Particle Mesh Ewald_) con un espaciado de cuadrícula de 0.12 nm y un límite de corte (_cutoff_) en espacio real de 1.0 nm.
- *Fuerzas de Van der Waals*: Se truncarán a una distancia de corte de 1.0 nm mediante una función de suavizado continuo de fuerzas.
- *Frecuencia de Registro*: Se guardarán las coordenadas atómicas y energías cada 10.0 picosegundos (ps), generando trayectorias individuales con un total de 10,000 marcos estructurales por cada corrida para su posterior análisis estadístico.

== 5.4 Flujo de Análisis de Trayectorias
Para cuantificar la estabilidad y el comportamiento termodinámico a lo largo de los 100 ns de simulación, se calcularán y graficarán las siguientes métricas bioinformáticas:
1. *Desviación Media Cuadrática (RMSD)*: Evaluada sobre los carbonos alfa de la proteína y sobre los átomos pesados del ligando para monitorizar la estabilidad estructural del complejo y determinar el tiempo de equilibrio.
2. *Fluctuación Media Cuadrática (RMSF)*: Calculada por residuo en la proteína para discernir qué regiones o bucles catalíticos adquieren mayor flexibilidad o rigidez debido a la presencia de los avocadenoides.
3. *Radio de Giro ($R_g$)*: Monitoreado para determinar si la proteína sufre transiciones conformacionales de colapso o apertura, sirviendo como medida indirecta de la compacidad de la diana.
4. *Análisis de Enlaces de Hidrógeno (H-Bonds)*: Cuantificación del número de puentes de hidrógeno directos establecidos entre el ligando y la proteína como función del tiempo, determinando su estabilidad temporal, distancias de enlace y ocupancias específicas.
