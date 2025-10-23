# FEM Implementation for Passive Diastolic Filling
A modular FEniCS-based implementation for solving the system of PDEs governing the passive diastolic filling of the left ventricle using the Holzapfel-Ogden constitutive model.
## Project Structure

```bash
.
├── config.py              # Configuration and parameters
├── mesh_utils.py          # Mesh loading and geometric utilities
├── material_model.py      # Holzapfel-Ogden material constitutive model
├── solver.py              # Finite element solver (Taylor-Hood formulation)
├── simulation.py          # Main simulation driver
├── postprocessing.py      # Post-processing and visualization utilities
├── main.py                # Main execution script
├── diagnostics.py         # Diagnostic and verification tools
└── README.md              
```

## Module Overview
```bash 
config.py

# Contains all configuration classes:

# a) MeshConfig: File paths for mesh, boundaries, and fiber directions
# b) MaterialParameters: Holzapfel-Ogden material parameters
# c) SolverParameters: Newton and Krylov solver settings
# d) SimulationParameters: Simulation control (pressure, steps, etc.)
```

```bash
mesh_utils.py

# Utilities for mesh and geometry:

load_ellipsoid_data(): # Load mesh, boundaries, and fiber fields
compute_cavity_volume(): # Calculate cavity volume via divergence theorem
compute_reference_volume(): # Get initial undeformed volume
```


```bash
material_model.py -> CardiacMaterial

# Implements Holzapfel-Ogden hyperelastic model:

# e) Passive stress with fiber, sheet, and cross-fiber contributions
# f) Active stress generation (fiber + transverse)
# g) Automatic kinematic computations
```

```bash
solver.py -> CardiacSolver

# Implements the solver based on the mixed Taylor-Hood formulation:

# h) Handles incompressibility constraint
# i) Implements boundary conditions
# j) Manages nonlinear solution process
```

```bash
simulation.py -> CardiacSimulation

# Defines and controls the overall simulation process
```

```bash
postprocessing.py

# Defines analysis and visualization functions:

compute_strain_stress(): # Green-Lagrange strain and PK2 stress
project_cauchy_stress(): # Project stress to DG0 space
plot_pv_curve(): # Generate pressure-volume curve
check_equilibrium_residual(): # Verify equilibrium
check_energy_balance(): # Verify work balance
check_normal_stress(): # Validate boundary conditions
```

```bash
main.py

# Execution file:

# k) Loads mesh and fiber data
# l) Initializes material and solver
# m) Runs pressure loading simulation
# n) Generates outputs and plots
```

```bash
diagnostics.py

Verification and debugging script that:

# o) Checks equilibrium residual
# p) Verifies energy balance
# q) Validates normal stress on boundaries
```
## Usage
### Basic Simulation
```bash 
python main.py 
```
This will:

Load the mesh from ```./HV2/```;
Run a 20-step pressure loading simulation;
Save displacement, strain, and stress fields;
Generate a pressure-volume curve to validate physiological behaviour.

### Run Diagnostics
```python diagnostics.py```
Performs verification checks on the solution.

### Custom Configuration
Modify parameters:
```from config import MaterialParameters, SimulationParameters```

## Change material stiffness
material_params = MaterialParameters(alpha=1.5)

## Change simulation parameters
sim_params = SimulationParameters()
sim_params.total_pressure = 15 * 133.3 / 1000.0  # 15 mmHg
sim_params.steps = 30
Output Files
After running main.py, the following files are generated:

HV2/results/u.pvd - Displacement field (ParaView)
HV2/results/strain.pvd - Green-Lagrange strain tensor
HV2/results/stress.pvd - 2nd Piola-Kirchhoff stress
HV2/results/output.xdmf - Cauchy stress (XDMF format)
pv_curve.pdf - Pressure-volume relationship plot

Requirements

FEniCS (tested with 2019.1.0)
NumPy
Matplotlib
Python 3.6+

Material Model
The implementation uses the Holzapfel-Ogden orthotropic model:
Passive stress:
σ = a·exp(b(I₁-3))·B 
    + 2·aₓ·(I₄ₓ-1)·exp(bₓ(I₄ₓ-1)²)·(f⊗f + s⊗s + ...)
    + aₓₛ·I₈ₓₛ·exp(bₓₛ·I₈ₓₛ²)·(f⊗s + s⊗f)
    - p·I

Active stress:
σₐ = Tₐ·(f⊗f + η·s⊗s + η·n⊗n)
where f, s, n are fiber, sheet, and normal directions.
Solver Configuration
The solver uses:

Mixed formulation: Taylor-Hood elements (P2-P1)
Nonlinear solver: Newton-Raphson
Linear solver: MUMPS (direct)
Convergence: Incremental criterion with tolerances of 1e-3

Boundary Conditions

Base (marker 10): Fixed (zero displacement)
Endocardium (marker 30): Pressure loading
Epicardium (marker 40): Traction-free

References

Holzapfel GA, Ogden RW (2009). Constitutive modelling of passive myocardium: a structurally based framework for material characterization. Phil Trans R Soc A 367:3445-3475.
Land S, et al. (2015). Verification of cardiac mechanics software: benchmark problems and solutions for testing active and passive material behaviour. Proc R Soc A 471:20150641.

License
MIT License - See LICENSE file for details.
Author
Converted from Jupyter notebook to modular structure.