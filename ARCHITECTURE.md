# XfCapture - Arquitectura Modular

## 📋 Resumen de la Transición

Se ha completado exitosamente la transición de un script monolítico (`xf_capture.py`) a una arquitectura modular profesional.

## 🏗️ Nueva Estructura

```
src/xf_capture/
├── __init__.py              # Módulo principal
├── cli.py                   # Interfaz de línea de comandos (CLI)
├── setup.py                 # Gestión de setup y configuración
├── runner.py                # Wrapper de Snakemake para ejecución
├── common.py                # Utilidades compartidas
├── workflows/
│   ├── Snakefile           # Pipeline de Snakemake
│   └── envs/               # Entornos Conda
├── resources/
│   └── ref_seqs.tar.gz     # Referencias empaquetadas
└── utils/
    ├── __init__.py
    ├── auxiliar_functions.py
    ├── plot_tree.R
    └── reconstruction_summary.R
```

## 📦 Módulos Principales

### 1. **cli.py** - Interfaz de Línea de Comandos
- **Responsabilidad**: Define la interfaz CLI con `argparse`
- **Comandos**:
  - `xf_capture setup --dir <path>`: Configura el entorno
  - `xf_capture run -i <input> -o <output>`: Ejecuta el pipeline
- **Características**:
  - Parsing de argumentos con grupos lógicos
  - Soporte para argumentos extra (pasados directamente a Snakemake)
  - Documentación integrada con `--help`

### 2. **setup.py** - Gestión de Setup
- **Responsabilidad**: Configuración inicial del entorno
- **Funciones principales**:
  - `setup_workflow()`: Orquesta el setup completo
  - `download_kraken2_database()`: Descarga y verifica base de datos Kraken2
  - `verify_md5()`: Verificación de integridad de archivos
  - `save_user_config()`: Guarda configuración persistente
  - `get_default_workflow_dir()`: Recupera configuración guardada
- **Características**:
  - Descarga automática de Kraken2 PlusPF 8GB
  - Verificación MD5 con reintentos
  - Extracción de secuencias de referencia
  - Configuración persistente en `~/.xf_capture/config.yaml`
  - Idempotencia (detecta recursos ya instalados)

### 3. **runner.py** - Wrapper de Snakemake
- **Responsabilidad**: Ejecución del pipeline de Snakemake
- **Funciones principales**:
  - `run_pipeline()`: Orquesta la ejecución en dos fases
  - `generate_config()`: Genera configuración dinámica para Snakemake
  - `check_snakemake()`: Verifica disponibilidad de Snakemake
  - `check_conda()`: Verifica disponibilidad de Conda/Mamba
- **Características**:
  - Ejecución bifásica:
    - **Fase 1**: QC → Taxonomy → Reconstruction → MLST
    - **Fase 2**: Phylogenetic analysis (solo muestras exitosas)
  - Modo automático vs. interactivo
  - Gestión de recursos (cores, jobs paralelos)
  - Integración con configuración persistente
  - Soporte para parámetros extra de Snakemake

### 4. **common.py** - Utilidades Compartidas
- **Responsabilidad**: Funciones compartidas entre módulos
- **Funciones**:
  - `check_conda()`: Detecta conda/mamba
  - `check_snakemake()`: Detecta snakemake
  - `print_header()`: Banner de la aplicación
  - `print_version()`: Información de versión
  - `format_bytes()`: Formateo de tamaños de archivo

## 🔄 Flujo de Ejecución

### Setup Workflow
```
xf_capture setup --dir /path/to/workflow
    ↓
cli.py (parse args)
    ↓
setup.py::setup_workflow()
    ↓
1. Crear estructura de directorios
2. Extraer secuencias de referencia
3. Descargar Kraken2 DB (con MD5)
4. Guardar configuración
5. Guardar config persistente
```

### Run Pipeline
```
xf_capture run -i data/ -o results/
    ↓
cli.py (parse args + extra_args)
    ↓
runner.py::run_pipeline()
    ↓
1. Verificar Snakemake instalado
2. Cargar workflow_dir persistente
3. Generar config.yaml dinámico
4. FASE 1: Ejecutar pipeline principal
5. Verificar muestras exitosas
6. FASE 2: Análisis filogenético (si hay éxito)
```

## 🎯 Mejoras Implementadas

### 1. **Separación de Responsabilidades**
- CLI separado de lógica de negocio
- Setup independiente de ejecución
- Utilidades compartidas reutilizables

### 2. **Configuración Persistente**
- Archivo de usuario en `~/.xf_capture/config.yaml`
- No necesidad de especificar `--workflow-dir` en cada run
- Override explícito posible

### 3. **Generación Dinámica de Configuración**
- Config YAML generado en tiempo de ejecución
- Integración con workflow_dir automática
- Override flexible de parámetros

### 4. **Gestión Robusta de Recursos**
- Verificación MD5 con reintentos
- Detección de recursos ya instalados
- Manejo de errores graceful

### 5. **Modo Bifásico Inteligente**
- Checkpoint entre fases
- Solo procesa muestras exitosas en Fase 2
- Confirmación interactiva vs. modo automático

### 6. **Extensibilidad**
- Argumentos extra pasados directamente a Snakemake
- Parámetros de recursos configurables
- Fácil adición de nuevas funcionalidades

## 🔧 Integración con Snakemake

### Ajustes en Snakefile
- Import de utilidades usando rutas relativas al paquete
- Scripts R accesibles desde `UTILS_DIR`
- Independiente de working directory

### Conda Environments
- Todos los entornos definidos en `workflows/envs/`
- Instalación automática en `workflow_dir/conda_envs/`
- Reutilización entre ejecuciones

## 📝 Configuración del Paquete

### pyproject.toml
```toml
[project.scripts]
xf_capture = "xf_capture.cli:main"

[tool.setuptools.package-data]
xf_capture = [
    "workflows/**/*",
    "resources/**/*",
    "utils/*.R",
    "utils/*.py"
]
```

## 🚀 Uso

### Setup inicial
```bash
xf_capture setup --dir ~/xf_workflow
```

### Ejecución básica
```bash
xf_capture run -i data/ -o results/
```

### Ejecución con parámetros personalizados
```bash
xf_capture run \
    -i data/ \
    -o results/ \
    --cores 16 \
    --kraken-jobs 2 \
    --alignment-jobs 6 \
    --iqtree-jobs 4 \
    --k2-mapping-memory \
    --no-auto
```

### Con argumentos extra de Snakemake
```bash
xf_capture run -i data/ -o results/ --dry-run
xf_capture run -i data/ -o results/ --forcerun verify_reconstruction
```

## 🔍 Verificación de la Implementación

### Tests realizados
- ✅ Instalación del paquete: `pip install -e .`
- ✅ Comando principal: `xf_capture --help`
- ✅ Comando setup: `xf_capture setup --help`
- ✅ Comando run: `xf_capture run --help`
- ✅ Imports correctos en todos los módulos
- ✅ Rutas de recursos configuradas
- ✅ Snakefile actualizado con imports correctos

## 📊 Ventajas de la Nueva Arquitectura

1. **Mantenibilidad**: Código organizado por responsabilidades
2. **Testabilidad**: Módulos independientes fáciles de testear
3. **Reusabilidad**: Funciones compartidas en `common.py`
4. **Escalabilidad**: Fácil agregar nuevos comandos o funcionalidades
5. **Profesionalismo**: Estructura estándar de proyectos Python
6. **Documentación**: Código autodocumentado con docstrings
7. **Configuración**: Sistema flexible de configuración persistente
8. **Robustez**: Manejo de errores y verificaciones

## 🎓 Principios Aplicados

- **Single Responsibility Principle**: Cada módulo tiene una responsabilidad clara
- **DRY (Don't Repeat Yourself)**: Funciones compartidas evitan duplicación
- **Separation of Concerns**: CLI, lógica y utilidades separados
- **Fail Fast**: Verificaciones tempranas de requisitos
- **Idempotencia**: Operaciones seguras para re-ejecución
