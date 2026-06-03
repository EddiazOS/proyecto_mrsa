#!/usr/bin/env bash
# =====================================================================
# correr_PBP2a_Ceftaroline.sh
# =====================================================================
# Script para reanudar NPT, correr producción y comprimir resultados
# diseñado para ejecutarse en el Pod de RunPod.
# =====================================================================

set -e # Terminar inmediatamente si algún comando falla

echo "=== [1/5] ACTIVANDO ENTORNO CONDA ==="
if [ -f /opt/miniconda/etc/profile.d/conda.sh ]; then
    source /opt/miniconda/etc/profile.d/conda.sh
    conda activate base
    echo "  [OK] Entorno conda activado desde /opt/miniconda."
elif [ -f /opt/conda/etc/profile.d/conda.sh ]; then
    source /opt/conda/etc/profile.d/conda.sh
    conda activate base
    echo "  [OK] Entorno conda activado desde /opt/conda."
else
    echo "  [WARNING] No se encontró conda.sh en rutas estándar. Continuando con PATH actual."
fi

# Verificar dependencias
if ! command -v gmx &> /dev/null; then
    echo "[ERROR] GROMACS ('gmx') no está en el PATH."
    exit 1
fi

RUN_DIR="md_run_PBP2a_Ceftaroline"
if [ ! -d "$RUN_DIR" ]; then
    echo "[ERROR] No se encontró el directorio: $RUN_DIR"
    exit 1
fi

cd "$RUN_DIR"
echo "  [OK] Directorio de trabajo: $(pwd)"

echo "=== [2/5] AJUSTANDO MDP DE NPT ==="
# Asegurar parámetros más estables
sed -i 's/tau_p *= *[0-9.]*/tau_p = 5.0/g' npt.mdp
sed -i 's/compressibility *= *[0-9.e-]*\b/compressibility = 4.5e-6/g' npt.mdp
echo "  [OK] Parámetros modificados en npt.mdp (tau_p = 5.0, compressibility = 4.5e-6)."

echo "=== [3/5] EJECUTANDO EQUILIBRACIÓN NPT ==="
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 5
# Usamos -update cpu para garantizar estabilidad en NPT
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update cpu -gpu_id 0
echo "  [OK] Equilibración NPT finalizada con éxito."

echo "=== [4/5] EJECUTANDO PRODUCCIÓN MD ==="
if [ ! -f "md.mdp" ]; then
    echo "[ERROR] No se encontró md.mdp. Asegúrate de que el script ejecutar_dinamica.py lo haya creado."
    exit 1
fi
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_production.tpr -maxwarn 5
# En producción podemos usar -update gpu para velocidad máxima
gmx mdrun -deffnm md_production -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu -gpu_id 0
echo "  [OK] Producción MD finalizada con éxito."

echo "=== [5/5] CREANDO BACKUP COMPRIMIDO ==="
cd ..
BACKUP_NAME="md_run_PBP2a_Ceftaroline_backup_$(date +%Y%m%d_%H%M%S).tar.gz"
echo "Comprimiendo carpeta $RUN_DIR en $BACKUP_NAME..."
tar -czf "$BACKUP_NAME" "$RUN_DIR"

echo "====================================================================="
# Calcular tamaño en MB/GB
FILE_SIZE=$(du -sh "$BACKUP_NAME" | cut -f1)
echo "  [ÉXITO] Proceso completado exitosamente."
echo "  [INFO] Respaldo guardado en: $BACKUP_NAME (Tamaño: $FILE_SIZE)"
echo "====================================================================="
