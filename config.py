"""
Configuration module for cardiac mechanics simulation.
Contains paths, material parameters, and solver settings.
"""
import os
from dolfin import Constant, parameters

# FEniCS compiler parameters
parameters["form_compiler"]["quadrature_degree"] = 2
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"


class MeshConfig:
    """Mesh and file paths configuration."""
    def __init__(self, mesh_dir='./HV2'):
        self.mesh_dir = mesh_dir
        self.mesh_name = os.path.join(mesh_dir, 'lvgmsh_mesh.xml')
        self.mf_name = os.path.join(mesh_dir, 'lvgmsh_fmaker.xml')
        self.fibre_name = os.path.join(mesh_dir, 'HV2_fibre_dir.xml')
        self.sheet_name = os.path.join(mesh_dir, 'HV2_sheet_dir.xml')
        self.normal_name = os.path.join(mesh_dir, 'HV2_normal_dir.xml')
        self.result_dir = os.path.join(mesh_dir, 'results/')
        
        # Create results directory if it doesn't exist
        os.makedirs(self.result_dir, exist_ok=True)


class MaterialParameters:
    """Material parameters for the cardiac tissue model."""
    def __init__(self, alpha=1.0):
        self.alpha = alpha
        self.eta = 0.1  # weighting of transversal fiber tension
        self.rho = Constant(1.0)  # kg
        
        # Holzapfel-Ogden model parameters
        self.a = 0.280 * alpha       # kPa
        self.a_f = 1.1685 * alpha    # kPa
        self.a_s = 0.116 * alpha     # kPa
        self.a_fs = 1.0 * alpha      # kPa
        self.b = 7.780               # dimensionless
        self.b_f = 11.83425          # dimensionless
        self.b_s = 8                 # dimensionless
        self.b_fs = 4                # dimensionless


class SolverParameters:
    """Solver configuration for nonlinear variational problems."""
    def __init__(self):
        self.ffc_options = {
            "optimize": True,
            "eliminate_zeros": True,
            "precompute_basis_const": True,
            "precompute_ip_const": True
        }
        
        self.form_compiler_parameters = {
            "keep_diagonal": True
        }
        
        self.newton_solver = {
            'absolute_tolerance': 1E-3,
            'relative_tolerance': 1E-3,
            'maximum_iterations': 50,
            'report': True,
            'error_on_nonconvergence': False,
            'linear_solver': 'mumps',
            'preconditioner': 'petsc_amg'
        }
        
        self.krylov_solver = {
            'absolute_tolerance': 1E-3,
            'relative_tolerance': 1E-3,
            'nonzero_initial_guess': True,
            'error_on_nonconvergence': True,
            'monitor_convergence': True,
            'report': True,
            'divergence_limit': 1E+12,
            'maximum_iterations': 10000
        }
    
    def apply_to_solver(self, solver):
        """Apply these parameters to a NonlinearVariationalSolver."""
        prm = solver.parameters
        prm['newton_solver']['convergence_criterion'] = 'incremental'
        prm['newton_solver']['absolute_tolerance'] = self.newton_solver['absolute_tolerance']
        prm['newton_solver']['relative_tolerance'] = self.newton_solver['relative_tolerance']
        prm['newton_solver']['maximum_iterations'] = self.newton_solver['maximum_iterations']
        prm['newton_solver']['report'] = self.newton_solver['report']
        prm['newton_solver']['error_on_nonconvergence'] = self.newton_solver['error_on_nonconvergence']
        prm['newton_solver']['linear_solver'] = self.newton_solver['linear_solver']
        
        prm['newton_solver']['krylov_solver']['absolute_tolerance'] = self.krylov_solver['absolute_tolerance']
        prm['newton_solver']['krylov_solver']['maximum_iterations'] = self.krylov_solver['maximum_iterations']
        prm['newton_solver']['krylov_solver']['relative_tolerance'] = self.krylov_solver['relative_tolerance']
        prm['newton_solver']['krylov_solver']['nonzero_initial_guess'] = self.krylov_solver['nonzero_initial_guess']
        prm['newton_solver']['krylov_solver']['error_on_nonconvergence'] = self.krylov_solver['error_on_nonconvergence']
        prm['newton_solver']['krylov_solver']['monitor_convergence'] = self.krylov_solver['monitor_convergence']
        prm['newton_solver']['krylov_solver']['report'] = self.krylov_solver['report']
        prm['newton_solver']['krylov_solver']['divergence_limit'] = self.krylov_solver['divergence_limit']


class SimulationParameters:
    """Simulation-specific parameters."""
    def __init__(self):
        self.total_pressure = 10 * 133.3 / 1000.0  # in kPa
        self.steps = 20
        self.mesh_scale = 0.1  # cm