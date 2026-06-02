# Base: GROMACS 2023.3 compilado con soporte CUDA completo (NVIDIA NGC)
FROM nvcr.io/hpc/gromacs:2023.3

# Instalar dependencias del sistema
RUN apt-get update && apt-get install -y wget git curl && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Instalar Miniconda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

ENV PATH="/opt/conda/bin:$PATH"

# Instalar dependencias Python (SIN gromacs — ya está en la imagen base con CUDA)
RUN conda install -y -c conda-forge \
    python=3.10 \
    acpype \
    openbabel \
    ambertools \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    ipykernel \
    jupyterlab && \
    conda clean -afy

# Directorio de trabajo persistente (montar Network Volume aquí en RunPod)
WORKDIR /workspace

# Exponer JupyterLab
EXPOSE 8888

# Comando por defecto: JupyterLab sin autenticación (seguro en entorno privado de RunPod)
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", \
     "--allow-root", "--NotebookApp.token=''", "--NotebookApp.password=''"]
