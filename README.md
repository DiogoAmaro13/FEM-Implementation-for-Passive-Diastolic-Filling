# FEM Implementation for Passive Diastolic Filling

A modular FEniCS-based implementation for solving the system of PDEs governing the passive diastolic filling of the left ventricle using the Holzapfel-Ogden constitutive model. This system is given by:

$$
\begin{equation}
\begin{cases}
\nabla \cdot \boldsymbol{\sigma} + \mathbf{b} = 0 & \text{in } \Omega(t), \\
\boldsymbol{\sigma} \cdot \mathbf{n} = \mathbf{t} & \text{on } \Gamma^N, \\
\mathbf{u} = \mathbf{u}_0 & \text{on}\ \Gamma^D.
\end{cases}.
\end{equation}
$$

Here, $\mathbf{b}$ is the body force density per unit volume, $\mathbf{n}$ is the normal direction of $\partial \Omega$, $\mathbf{t}$ is the traction force, $\Gamma^N$ and $\Gamma^D$ are the Neumann and Dirichlet boundaries.

For further theoretical context, we refer to the paper accompanying this code, which can be accessed here.

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
## Usage
### Basic Simulation
```bash 
python main.py 
```
This will:

Load the mesh from ```./HV2/```; <br/>
Run a 20-step pressure loading simulation; <br/>
Save displacement, strain, and stress fields; <br/>
Generate a pressure-volume curve to validate physiological behaviour.

### Custom Configuration
To modify parameters:
```bash 
from config import MaterialParameters, SimulationParameters

# Change material stiffness
material_params = MaterialParameters(alpha=1.5)

# Change simulation parameters
sim_params = SimulationParameters()
sim_params.total_pressure = 15 * 133.3 / 1000.0  # 15 mmHg
sim_params.steps = 30
```
## Output Files
After running ```main.py```, the following files are generated:

```HV2/results/u.pvd``` - Displacement field <br/>
```HV2/results/strain.pvd``` - Green-Lagrange strain tensor <br/>
```HV2/results/stress.pvd``` - 2nd Piola-Kirchhoff stress <br/>
```HV2/results/output.xdmf``` - Cauchy stress <br/>
```pv_curve.pdf``` - Pressure-volume relationship plot

Each one of them must be visualized using an appropriate software. In this case, ParaView was used. It should be noted that each of the three ```.pvd``` files compile 20 different ```.vtu``` files, so they must be in the same directory when opening the ```.pvd``` files in ParaView. The same goes for the ```.xdmf``` file, which uses the auxiliary ```.h5``` file.

## Requirements

FEniCS (tested with 2019.1.0) <br/>
NumPy <br/>
Matplotlib <br/>
Python 3.6+

However, we strongly recommend activating the conda environment which already has all necessary dependencies. Just write
```bash
conda activate jax-fem-env
```
in the terminal.

## License
This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.