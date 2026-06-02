# Guía de Ejecución Local para Dinámica Molecular en GROMACS

Esta guía contiene las instrucciones detalladas para configurar el entorno y ejecutar el pipeline de dinámica molecular multi-complejo de forma local en los computadores del laboratorio, aprovechando la aceleración por hardware (GPU NVIDIA CUDA).

---

## 1. Configuración del Entorno de Conda

Hemos creado el archivo [environment_md.yml](file:///home/echoes/projects/cursos-semillero/proyecto_mrsa/environment_md.yml) que empaqueta todas las dependencias necesarias de Python, ACPYPE y GROMACS.

### Paso A: Clonar el Repositorio (Si es un PC nuevo en el laboratorio)
```bash
git clone https://github.com/EddiazOS/proyecto_mrsa.git
cd proyecto_mrsa
```

### Paso B: Crear el Entorno Virtual
Asegúrate de tener Anaconda o Miniconda instalado. Ejecuta en tu terminal:
```bash
conda env create -f environment_md.yml
```
*Esto creará un entorno aislado llamado `md_env` con GROMACS, ACPYPE, OpenBabel y las librerías de análisis de datos.*

### Paso C: Activar el Entorno
```bash
conda activate md_env
```

---

## 2. Aceleración por GPU (CUDA) en Computadores de Laboratorio

Para obtener el máximo rendimiento en un PC con tarjeta gráfica dedicada (NVIDIA GeForce RTX, Tesla, Quadro, etc.), es fundamental que GROMACS utilice aceleración CUDA:

### En Linux (Ubuntu / Debian / CentOS)
Si el computador corre Linux nativo y tiene controladores de NVIDIA instalados, el paquete de GROMACS instalado por Conda (`gromacs`) detectará y usará automáticamente CUDA. 
* Puedes verificar que GROMACS ve la GPU corriendo:
  ```bash
  gmx mdrun -h | grep -i cuda
  ```

### En Windows (WSL2 - Windows Subsystem for Linux)
En computadores del laboratorio que utilicen Windows, **no se recomienda ejecutar GROMACS nativamente en Windows**, ya que ACPYPE y GROMACS funcionan con mucha mayor estabilidad en entornos basados en Unix.
1. Abre una terminal de **WSL2 (Ubuntu)** en Windows.
2. Instala la versión de Conda para Linux dentro de WSL2.
3. Activa tu entorno virtual `md_env` e inicia la simulación. WSL2 pasará automáticamente las instrucciones de cálculo GPU al controlador de Windows (NVIDIA CUDA en WSL es soportado de forma nativa).

---

## 3. Cómo Lanzar y Monitorear la Simulación Secuencial

Una vez configurado y activo el entorno `md_env`:

### Lanzar todas las simulaciones (secuencial de los 4 complejos):
Puedes especificar el tiempo de simulación deseado en nanosegundos (por ejemplo, `10` ns para pruebas o `100` ns para datos finales de publicación):
```bash
# Sintaxis: python scripts/correr_todos.py [TIEMPO_NS]
python scripts/correr_todos.py 100
```

### Monitoreo en Tiempo Real:
Como la salida de GROMACS está activada en tiempo real (`silence=False`), verás el progreso exacto en tu terminal. Adicionalmente, puedes monitorear de forma externa con:
```bash
# Ver uso y temperatura de la GPU NVIDIA en tiempo real
watch -n 1 nvidia-smi
```

---

## 4. Respaldos y Resultados
* **Ubicación local:** Los resultados de cada complejo se almacenan en carpetas dedicadas del repositorio llamadas `md_run_MurG_Afzelin/`, `md_run_MurG_Quercetin/`, etc.
* **Google Drive:** El script tiene lógica para detectar si la ruta de Google Drive `/content/drive` está activa. En una máquina local del laboratorio, simplemente omitirá el respaldo a Drive e informará que la simulación local finalizó con éxito, manteniendo tus datos locales a salvo en el disco duro físico del computador.
