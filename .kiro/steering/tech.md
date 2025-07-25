# Technology Stack

## Build System
- **Build Backend**: setuptools with pyproject.toml configuration
- **Package Manager**: PDM (Python Dependency Management)
- **Python Versions**: 3.8+ (supports 3.7-3.12)

## Core Dependencies
- **NumPy** (>=1.21.0) - Array operations and linear algebra
- **SciPy** (>=1.7.0) - Scientific computing, eigenvalue solvers
- **Matplotlib** (>=3.4.0) - Plotting and visualization
- **ASE** (>=3.19) - Atomic Simulation Environment

## Optional Dependencies
- **sisl** (>=0.9.0) - SIESTA integration (optional-dependencies)
- **pytest** (>=6.2.0) - Testing framework
- **black** (>=21.6b0) - Code formatting (dev)
- **sphinx** (>=4.1.2) - Documentation generation (dev)

## Code Quality Tools
- **Ruff** - Fast Python linter and formatter (replaces flake8, black)
- **Pre-commit** - Git hooks for code quality
- **Target**: Python 3.8 compatibility
- **Line Length**: 88 characters (Black standard)
- **Quote Style**: Double quotes

## Common Commands

### Development Setup
```bash
# Install in development mode
pip install -e .

# Install with optional dependencies
pip install -e .[siesta,test,dev]
```

### Testing
```bash
# Run tests
pytest

# Run specific test modules
python examples/SIESTA/read_Fe.py
```

### Code Quality
```bash
# Run pre-commit hooks
pre-commit run --all-files

# Manual linting and formatting
ruff check --fix --extend-select I .
ruff format .
```

### Building and Distribution
```bash
# Build package
python -m build

# Upload to PyPI (via upload_to_pip.sh)
./upload_to_pip.sh
```

## Architecture Patterns
- Abstract base classes for extensibility (`BaseHamiltonian`, `Hamiltonian`)
- Factory pattern for different DFT code readers
- Dataclass usage for structured data
- NumPy array-centric design for performance
- Lazy evaluation where possible for large matrices