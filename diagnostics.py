"""
Diagnostic and verification utilities for cardiac mechanics simulation.
Run this after main simulation to check equilibrium and energy balance.
"""
from dolfin import *
from config import MeshConfig, MaterialParameters, SolverParameters, SimulationParameters
from mesh_utils import load_ellipsoid_data
from material_model import CardiacMaterial
from solver import CardiacSolver
from postprocessing import (
    check_equilibrium_residual, 
    check_energy_balance, 
    check_normal_stress
)


def run_diagnostics():
    """Run diagnostic checks on the simulation results."""
    
    print("=" * 60)
    print("Cardiac Mechanics Diagnostics")
    print("=" * 60)
    
    # ========== Load Configuration ==========
    mesh_config = MeshConfig(mesh_dir='./HV2')
    material_params = MaterialParameters(alpha=1.0)
    solver_params = SolverParameters()
    sim_params = SimulationParameters()
    
    # ========== Load Mesh and Setup ==========
    print("\n[1/4] Loading mesh and initializing solver...")
    mesh, boundary_markers, numbering, fibers = load_ellipsoid_data(
        mesh_config.mesh_name,
        mesh_config.mf_name,
        mesh_config.fibre_name,
        mesh_config.sheet_name,
        mesh_config.normal_name,
        scale=sim_params.mesh_scale
    )
    
    material = CardiacMaterial(material_params, fibers)
    solver = CardiacSolver(mesh, boundary_markers, numbering, material, solver_params)
    
    # Run one solve to get a solution
    print("\n[2/4] Solving for final pressure state...")
    final_pressure = sim_params.total_pressure
    solver.solve(final_pressure, active_tension=0.0)
    print(f"  ✓ Solved at pressure = {final_pressure:.4f} kPa")
    
    # ========== Check Equilibrium Residual ==========
    print("\n[3/4] Checking equilibrium residual...")
    try:
        residual = check_equilibrium_residual(
            mesh,
            solver.P_p,
            solver.P_p,  # Using passive stress only
            solver.v,
            solver.q,
            solver.kinematics['J'],
            boundary_markers,
            numbering,
            solver.p0
        )
        print(f"  Equilibrium residual: {residual:.6e}")
        if residual < 1e-6:
            print("  ✓ Equilibrium satisfied (residual < 1e-6)")
        else:
            print("  ⚠ Warning: Equilibrium residual is high")
    except Exception as e:
        print(f"  ⚠ Could not compute equilibrium residual: {e}")
    
    # ========== Check Energy Balance ==========
    print("\n[4/4] Checking energy balance...")
    try:
        u = solver.get_displacement()
        p = solver.get_pressure()
        
        total_eq, internal_eq, external_eq = check_energy_balance(
            mesh,
            solver.P_p,
            u,
            p,
            solver.kinematics['J'],
            boundary_markers,
            numbering,
            solver.p0
        )
        
        print(f"  Total energy:    {total_eq:.6e}")
        print(f"  Internal work:   {internal_eq:.6e}")
        print(f"  External work:   {external_eq:.6e}")
        
        if abs(total_eq) < 1e-6:
            print("  ✓ Energy balance satisfied")
        else:
            print("  ⚠ Warning: Energy imbalance detected")
    except Exception as e:
        print(f"  ⚠ Could not compute energy balance: {e}")
    
    # ========== Check Normal Stress ==========
    print("\n[5/5] Checking normal stress on endocardium...")
    try:
        # Create deformed mesh
        V = VectorFunctionSpace(mesh, "CG", 1)
        u_int = project(u, V)
        deformed_mesh = Mesh(mesh)
        ALE.move(deformed_mesh, u_int)
        
        avg_stress, avg_pressure, endo_area = check_normal_stress(
            mesh,
            deformed_mesh,
            solver.total_cauchy,
            u,
            p,
            solver.p0,
            boundary_markers,
            numbering
        )
        
        print(f"  Endocardial area: {endo_area:.4f} cm²")
        print(f"  Average normal stress:  {avg_stress:.6f} kPa")
        print(f"  Average pressure load:  {avg_pressure:.6f} kPa")
        print(f"  Applied pressure:       {float(solver.p0):.6f} kPa")
        
        stress_error = abs(avg_stress - float(solver.p0))
        if stress_error < 0.01:
            print("  ✓ Normal stress matches applied pressure")
        else:
            print(f"  ⚠ Stress-pressure mismatch: {stress_error:.6f} kPa")
    except Exception as e:
        print(f"  ⚠ Could not compute normal stress: {e}")
    
    print("\n" + "=" * 60)
    print("✅ Diagnostics complete")
    print("=" * 60)


if __name__ == "__main__":
    run_diagnostics()