"""
Main simulation driver for cardiac mechanics.
"""
import numpy as np
from mesh_utils import compute_cavity_volume
from postprocessing import compute_strain_stress, ResultsWriter


class CardiacSimulation:
    """
    Main simulation class coordinating the cardiac mechanics analysis.
    """
    
    def __init__(self, mesh, boundary_markers, numbering, solver, sim_params):
        """
        Initialize simulation.
        
        Args:
            mesh: FEniCS Mesh object
            boundary_markers: MeshFunction with boundary markers
            numbering: Dictionary mapping boundary names to markers
            solver: CardiacSolver instance
            sim_params: SimulationParameters instance
        """
        self.mesh = mesh
        self.boundary_markers = boundary_markers
        self.numbering = numbering
        self.solver = solver
        self.params = sim_params
        
        # Storage for results
        self.volume_history = None
        self.pressure_history = None
    
    def run_pressure_loading(self, result_dir):
        """
        Run pressure loading simulation.
        
        Args:
            result_dir: Directory to save results
        
        Returns:
            tuple: (volume_array, pressure_array)
        """
        steps = self.params.steps
        total_pressure = self.params.total_pressure
        
        # Initialize storage
        volume = np.zeros(steps + 1)
        pressure = np.zeros(steps + 1)
        
        # Initial state
        volume[0] = compute_cavity_volume(
            self.mesh, self.boundary_markers, self.numbering, 
            self.solver.get_displacement()
        )
        pressure[0] = 0.0
        
        # Initialize results writer
        writer = ResultsWriter(result_dir)
        
        # Incremental loading loop
        for i in range(steps):
            step_num = i + 1
            current_pressure = total_pressure / steps * step_num
            
            print(f'{step_num}/{steps}, pressure = {current_pressure:.6f} kPa')
            
            # Solve with current pressure
            converged, num_iter = self.solver.solve(current_pressure, active_tension=0.0)
            
            if not converged:
                print(f"Warning: Step {step_num} did not converge!")
            
            # Extract displacement
            u = self.solver.get_displacement()
            
            # Compute strain and stress
            eps, PK2 = compute_strain_stress(self.mesh, u)
            
            # Compute volume
            volume[step_num] = compute_cavity_volume(
                self.mesh, self.boundary_markers, self.numbering, u
            )
            pressure[step_num] = current_pressure
            
            # Write results
            writer.write_step(u, eps, PK2, step_num)
        
        print("\n✅ Simulation complete. Results saved in:")
        print(f"  {result_dir}/u.pvd")
        print(f"  {result_dir}/strain.pvd")
        print(f"  {result_dir}/stress.pvd")
        
        # Store results
        self.volume_history = volume
        self.pressure_history = pressure
        
        return volume, pressure
    
    def get_results(self):
        """
        Get simulation results.
        
        Returns:
            tuple: (volume_array, pressure_array)
        """
        if self.volume_history is None:
            raise RuntimeError("Simulation has not been run yet.")
        return self.volume_history, self.pressure_history