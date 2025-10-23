"""
Main execution script for cardiac mechanics simulation.
"""
from config import MeshConfig, MaterialParameters, SolverParameters, SimulationParameters
from mesh_utils import load_ellipsoid_data, compute_reference_volume
from material_model import CardiacMaterial
from solver import CardiacSolver
from simulation import CardiacSimulation
from postprocessing import plot_pv_curve, project_cauchy_stress, write_xdmf_results


def main():
    """Main execution function."""
    
    # ========== Configuration ==========
    print("=" * 60)
    print("Cardiac Mechanics Simulation")
    print("=" * 60)
    
    mesh_config = MeshConfig(mesh_dir='./HV2')
    material_params = MaterialParameters(alpha=1.0)
    solver_params = SolverParameters()
    sim_params = SimulationParameters()
    
    # ========== Load Mesh and Fiber Data ==========
    print("\n[1/6] Loading mesh and fiber data...")
    mesh, boundary_markers, numbering, fibers = load_ellipsoid_data(
        mesh_config.mesh_name,
        mesh_config.mf_name,
        mesh_config.fibre_name,
        mesh_config.sheet_name,
        mesh_config.normal_name,
        scale=sim_params.mesh_scale
    )
    print(f"  ✓ Mesh loaded: {mesh.num_vertices()} vertices, {mesh.num_cells()} cells")
    
    # ========== Compute Reference Volume ==========
    print("\n[2/6] Computing reference cavity volume...")
    vol_ref = compute_reference_volume(mesh, boundary_markers, numbering)
    print(f"  ✓ Reference volume: {vol_ref:.4f} mL")
    
    # ========== Initialize Material Model ==========
    print("\n[3/6] Initializing material model...")
    material = CardiacMaterial(material_params, fibers)
    print("  ✓ Holzapfel-Ogden material model initialized")
    
    # ========== Setup Solver ==========
    print("\n[4/6] Setting up solver...")
    solver = CardiacSolver(mesh, boundary_markers, numbering, material, solver_params)
    print("  ✓ Solver configured (Taylor-Hood elements)")
    
    # ========== Run Simulation ==========
    print("\n[5/6] Running pressure loading simulation...")
    print(f"  Target pressure: {sim_params.total_pressure:.4f} kPa")
    print(f"  Number of steps: {sim_params.steps}")
    
    simulation = CardiacSimulation(mesh, boundary_markers, numbering, solver, sim_params)
    volume, pressure = simulation.run_pressure_loading(mesh_config.result_dir)
    
    # ========== Post-processing ==========
    print("\n[6/6] Post-processing...")
    
    # Plot pressure-volume curve
    print("  → Plotting P-V curve...")
    plot_pv_curve(volume, pressure, output_file="pv_curve.pdf")
    
    # Project Cauchy stress
    print("  → Projecting Cauchy stress...")
    cauchy_stress = project_cauchy_stress(mesh, solver.w, material)
    
    # Write XDMF results
    print("  → Writing XDMF output...")
    write_xdmf_results(mesh, mesh_config.result_dir, solver.total_cauchy)
    
    print("\n" + "=" * 60)
    print("✅ All tasks completed successfully!")
    print("=" * 60)
    print(f"\nResults summary:")
    print(f"  Initial volume:  {volume[0]:.4f} mL")
    print(f"  Final volume:    {volume[-1]:.4f} mL")
    print(f"  Volume change:   {volume[-1] - volume[0]:.4f} mL")
    print(f"  Final pressure:  {pressure[-1]:.4f} kPa ({pressure[-1]*1000/133.3:.2f} mmHg)")
    print(f"\nOutput files:")
    print(f"  • {mesh_config.result_dir}/u.pvd")
    print(f"  • {mesh_config.result_dir}/strain.pvd")
    print(f"  • {mesh_config.result_dir}/stress.pvd")
    print(f"  • {mesh_config.result_dir}/output.xdmf")
    print(f"  • pv_curve.pdf")
    print()


if __name__ == "__main__":
    main()