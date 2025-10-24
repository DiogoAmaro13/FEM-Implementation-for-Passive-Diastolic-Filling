"""
Solver module for cardiac mechanics problems.
"""
from dolfin import (
    VectorElement, FiniteElement, FunctionSpace, Function, TestFunctions,
    DirichletBC, Constant, Measure, Identity, grad, det, inv, inner, dx, dot,
    derivative, NonlinearVariationalProblem, NonlinearVariationalSolver,
    split, FacetNormal
)

class CardiacSolver:
    """
    Solver for cardiac mechanics using mixed formulation (Taylor-Hood elements).
    """
    
    def __init__(self, mesh, boundary_markers, numbering, material_model, solver_params):
        """
        Initialize the cardiac mechanics solver.
        
        Args:
            mesh: FEniCS Mesh object
            boundary_markers: MeshFunction with boundary markers
            numbering: Dictionary mapping boundary names to markers
            material_model: CardiacMaterial instance
            solver_params: SolverParameters instance
        """
        self.mesh = mesh
        self.boundary_markers = boundary_markers
        self.numbering = numbering
        self.material = material_model
        self.solver_params = solver_params
        
        # Create function space
        self._setup_function_space()
        
        # Create boundary conditions
        self._setup_boundary_conditions()
        
        # Create pressure and active tension parameters
        self.p0 = Constant(0)
        self.T_a = Constant(0)
        
        # Build weak form
        self._setup_weak_form()
    
    def _setup_function_space(self):
        """Set up Taylor-Hood mixed element (quadratic displacement + linear pressure)."""
        V = VectorElement("CG", self.mesh.ufl_cell(), 2)
        P = FiniteElement("CG", self.mesh.ufl_cell(), 1)
        TH = V * P
        self.W = FunctionSpace(self.mesh, TH)
        
        # Define functions
        self.w = Function(self.W)
        self.u, self.p = split(self.w)
        self.v, self.q = TestFunctions(self.W)
    
    def _setup_boundary_conditions(self):
        """Set up boundary conditions (fixed base)."""
        zero_displacement = Constant(("0.0", "0.0", "0.0"))
        bcr = DirichletBC(self.W.sub(0), zero_displacement, 
                         self.boundary_markers, self.numbering['BASE'])
        self.bcs = [bcr]
    
    def _setup_weak_form(self):
        """Build the weak form of the problem."""
        # Compute kinematics
        kinematics = self.material.compute_kinematics(self.u)
        F = kinematics['F']
        J = kinematics['J']
        
        # Compute stresses
        passive_stress = self.material.passive_cauchy_stress(kinematics, self.p)
        active_stress = self.material.active_cauchy_stress(kinematics, self.T_a)
        total_cauchy = passive_stress + active_stress
        
        # First Piola-Kirchhoff stress
        P_p = self.material.first_piola_kirchhoff(passive_stress, kinematics)
        P_a = self.material.first_piola_kirchhoff(active_stress, kinematics)
        P = P_p + P_a
        
        # Define measure and normal
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundary_markers)
        n_mesh = FacetNormal(self.mesh)
        
        # Weak form: equilibrium + incompressibility + pressure boundary condition
        self.eq = (inner(P, grad(self.v)) * dx 
                   + inner(J - 1, self.q) * dx 
                   + dot(J * inv(F).T * n_mesh * self.p0, self.v) * ds(self.numbering['ENDO']))
        
        # Jacobian
        self.Jac = derivative(self.eq, self.w)
        
        # Store for postprocessing
        self.kinematics = kinematics
        self.passive_stress = passive_stress
        self.total_cauchy = total_cauchy
        self.P_p = P_p
    
    def solve(self, pressure_value, active_tension=0.0):
        """
        Solve the nonlinear problem for given pressure and active tension.
        
        Args:
            pressure_value: Endocardial pressure (kPa)
            active_tension: Active fiber tension (default: 0.0)
        
        Returns:
            tuple: (converged, num_iterations)
        """
        # Set parameters
        self.p0.assign(pressure_value)
        self.T_a.assign(active_tension)
        
        # Create problem and solver
        problem = NonlinearVariationalProblem(self.eq, self.w, self.bcs, J=self.Jac)
        solver = NonlinearVariationalSolver(problem)
        
        # Apply solver parameters
        self.solver_params.apply_to_solver(solver)
        
        # Solve
        (num_iterations, converged) = solver.solve()
        
        return converged, num_iterations
    
    def get_displacement(self):
        """Return the displacement field."""
        return self.w.split()[0]
    
    def get_pressure(self):
        """Return the pressure field."""
        return self.w.split()[1]

