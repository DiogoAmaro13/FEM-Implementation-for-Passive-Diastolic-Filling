"""
Post-processing utilities for cardiac mechanics simulations.
"""
import numpy as np
#from dolfin import (
#    File, TensorFunctionSpace, VectorFunctionSpace, Function, project,
#    Identity, grad, inner, Constant, XDMFFile, FiniteElement, VectorElement,
#    TensorElement, FunctionSpace, Mesh, ALE, FacetNormal, Measure, assemble, dot
#)
from dolfin import *
import matplotlib.pyplot as plt
from matplotlib import rcParams


class ResultsWriter:
    """Handle writing results to files."""
    
    def __init__(self, result_dir):
        """
        Initialize results writer.
        
        Args:
            result_dir: Directory to save results
        """
        self.result_dir = result_dir
        self.save_u = File(f"{result_dir}/u.pvd")
        self.save_e = File(f"{result_dir}/strain.pvd")
        self.save_S = File(f"{result_dir}/stress.pvd")
    
    def write_step(self, u, eps, PK2, step):
        """
        Write fields for a single time step.
        
        Args:
            u: Displacement field
            eps: Strain tensor
            PK2: Second Piola-Kirchhoff stress
            step: Step number
        """
        u.rename('displacement', '')
        eps.rename('strain', '')
        PK2.rename('stress', '')
        
        self.save_u << (u, float(step))
        self.save_e << (eps, float(step))
        self.save_S << (PK2, float(step))


def compute_strain_stress(mesh, u):
    """
    Compute Green-Lagrange strain and 2nd Piola-Kirchhoff stress.
    
    Args:
        mesh: FEniCS Mesh
        u: Displacement field
    
    Returns:
        tuple: (strain_function, stress_function)
    """
    T = TensorFunctionSpace(mesh, 'P', 1)
    
    # Compute strain (Green-Lagrange)
    F = Identity(3) + grad(u)
    C = F.T * F
    E = 0.5 * (C - Identity(3))
    eps = Function(T, name="strain")
    eps.assign(project(E, T))
    
    # Compute 2nd Piola-Kirchhoff stress (placeholder Neo-Hookean)
    mu = Constant(10.0)
    lmbda = Constant(10.0)
    S = lmbda * inner(E, Identity(3)) * Identity(3) + 2 * mu * E
    PK2 = Function(T, name="PK2_stress")
    PK2.assign(project(S, T))
    
    return eps, PK2


def project_cauchy_stress(mesh, w, material_model):
    """
    Project Cauchy stress onto DG0 space.
    
    Args:
        mesh: FEniCS Mesh
        w: Solution function (mixed displacement-pressure)
        material_model: CardiacMaterial instance
    
    Returns:
        Function: Cauchy stress in DG0 space
    """
    Vsig = TensorFunctionSpace(mesh, "DG", degree=0)
    sig_num = Function(Vsig, name="Stress Numeric")
    
    u = w.split()[0]
    p = w.split()[1]
    
    # Recompute kinematics
    kinematics = material_model.compute_kinematics(u)
    
    # Recompute total Cauchy stress
    passive_stress = material_model.passive_cauchy_stress(kinematics, p)
    active_stress = material_model.active_cauchy_stress(kinematics, Constant(0))
    cauchy_stress = passive_stress + active_stress
    
    sig_num.assign(project(cauchy_stress, Vsig))
    return sig_num


def write_xdmf_results(mesh, result_dir, cauchy_stress):
    """
    Write results in XDMF format.
    
    Args:
        mesh: FEniCS Mesh
        result_dir: Directory to save results
        cauchy_stress: Cauchy stress tensor (UFL expression)
    """
    dFE = FiniteElement("DG", mesh.ufl_cell(), 0)
    tFE = TensorElement(dFE)
    
    fileResults = XDMFFile(result_dir + "output.xdmf")
    fileResults.parameters["flush_output"] = True
    fileResults.parameters["functions_share_mesh"] = True
    
    W = FunctionSpace(mesh, tFE)
    stress = Function(W, name='Stress')
    stress.assign(project(cauchy_stress, W))
    fileResults.write(stress, 0.)


def plot_pv_curve(volume, pressure, output_file="pv_curve.pdf"):
    """
    Plot pressure-volume curve.
    
    Args:
        volume: Array of volumes (mL)
        pressure: Array of pressures (kPa)
        output_file: Output filename
    """
    # Configure matplotlib for LaTeX rendering
    rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "axes.labelsize": 10,
        "axes.titlesize": 12,
        "legend.fontsize": 8,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9
    })
    
    # Convert pressure to mmHg
    pressure_mmhg = pressure * 1000.0 / 133.3
    
    plt.figure()
    plt.plot(volume, pressure_mmhg)
    plt.xlabel('Volume (mL)')
    plt.ylabel('Pressure (mmHg)')
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.show()


def check_equilibrium_residual(mesh, P_p, P, v, q, J, boundary_markers, numbering, p0):
    """
    Check equilibrium residual.
    
    Args:
        mesh: FEniCS Mesh
        P_p: Passive first Piola-Kirchhoff stress
        P: Total first Piola-Kirchhoff stress
        v, q: Test functions
        J: Jacobian determinant
        boundary_markers: MeshFunction
        numbering: Boundary numbering dict
        p0: Applied pressure
    
    Returns:
        float: Equilibrium residual energy
    """
    import dolfin
    from dolfin import dx, assemble, inner, grad, FacetNormal, Constant
    

    n_mesh = FacetNormal(mesh)
    ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    
    Div_P = dolfin.div(P_p)
    Psi_V = inner(Div_P, Div_P)
    
    N = FacetNormal(mesh)
    Jump_P_N = dolfin.jump(P, N)
    cell_h = Constant(mesh.hmin())
    Psi_F = inner(Jump_P_N, Jump_P_N) / cell_h
    dF = dolfin.Measure('dS', mesh)
    
    Psi_F_energy = assemble(Psi_F * dF)
    DivP_energy = assemble(Psi_V * dx)
    
    return DivP_energy + Psi_F_energy


def check_energy_balance(mesh, P_p, u, p, J, boundary_markers, numbering, p0):
    """
    Check energy balance.
    
    Args:
        mesh: FEniCS Mesh
        P_p: Passive first Piola-Kirchhoff stress
        u: Displacement field
        p: Pressure field
        J: Jacobian determinant
        boundary_markers: MeshFunction
        numbering: Boundary numbering dict
        p0: Applied pressure
    
    Returns:
        tuple: (total_energy, internal_energy, external_energy)
    """
    from dolfin import assemble, inner, grad, inv, FacetNormal, dot
    
    n_mesh = FacetNormal(mesh)
    ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    
    eq_internal = inner(P_p, grad(u)) * dx + inner(J - 1, p) * dx
    eq_external = dot(J * inv(Identity(3) + grad(u)).T * n_mesh * p0, u) * ds(numbering['ENDO'])
    
    total_eq_internal = assemble(eq_internal)
    total_eq_external = assemble(eq_external)
    total_eq = total_eq_internal + total_eq_external
    
    return total_eq, total_eq_internal, total_eq_external


def check_normal_stress(mesh, deformed_mesh, cauchy_stress, u, p, p0, boundary_markers, numbering):
    """
    Verify average normal stress on endocardium.
    
    Args:
        mesh: Original mesh
        deformed_mesh: Deformed mesh
        cauchy_stress: Cauchy stress tensor
        u: Displacement field
        p: Pressure field
        p0: Applied pressure
        boundary_markers: MeshFunction
        numbering: Boundary numbering dict
    
    Returns:
        tuple: (avg_stress_normal, avg_pressure_normal, endo_area)
    """
    from dolfin import assemble, inner, dot, inv, det, grad, Identity
    
    F = Identity(3) + grad(u)
    J = det(F)
    
    n_mesh = FacetNormal(mesh)
    ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    
    endo_area = assemble(1 * ds(numbering['ENDO']))
    
    P_p = J * cauchy_stress * inv(F).T
    stress_nor = dot(P_p * n_mesh, n_mesh) * ds(numbering['ENDO'])
    total_stress_nor = assemble(stress_nor)
    
    p_nor = dot(-J * p0 * inv(F).T * n_mesh, n_mesh) * ds(numbering['ENDO'])
    total_p_nor = assemble(p_nor)
    
    return total_stress_nor / endo_area, total_p_nor / endo_area, endo_area