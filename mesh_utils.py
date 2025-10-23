"""
Mesh loading and geometric utilities.
"""
from dolfin import *


def load_ellipsoid_data(mesh_name, mf_name, fibre_name, sheet_name, normal_name, scale=0.1):
    """
    Load mesh, boundary markers and fiber fields.
    
    Args:
        mesh_name: Path to mesh XML file
        mf_name: Path to mesh function (boundary markers) XML file
        fibre_name: Path to fiber direction XML file
        sheet_name: Path to sheet direction XML file
        normal_name: Path to normal direction XML file
        scale: Scaling factor for mesh coordinates (default: 0.1 for cm)
    
    Returns:
        tuple: (mesh, mf, numbering, [fiber, sheet, cross_sheet])
            - mesh: FEniCS Mesh object
            - mf: MeshFunction with boundary markers
            - numbering: Dictionary mapping boundary names to markers
            - list of fiber/sheet/normal Function objects
    """
    # Load mesh and mesh function
    mesh = Mesh(mesh_name)
    if scale is not None:
        mesh.scale(scale)

    cell = mesh.ufl_cell()
                
    # Load boundary markers (expects an XML MeshFunction file)
    mf = MeshFunction("size_t", mesh, mf_name)

    numbering = {
        "BASE": 10,
        "ENDO": 30,
        "EPI": 40
    }

    # Vector function space for fiber directions
    fiber_element = VectorElement("Lagrange", cell, degree=1)
    fiber_space = FunctionSpace(mesh, fiber_element)

    # Load fiber / sheet / normal as Functions (file must exist)
    fiber = Function(fiber_space, fibre_name)
    sheet = Function(fiber_space, sheet_name)
    cross_sheet = Function(fiber_space, normal_name)

    return mesh, mf, numbering, [fiber, sheet, cross_sheet]


def compute_cavity_volume(mesh, mf, numbering, u=None):
    """
    Compute cavity volume via surface integral using divergence theorem.
    
    For small displacements (linearized), use V = (1/3) ∫_{S} (X + u) · N dS (approx).
    If u is None, compute reference volume V0 = (1/3) ∫_S X·N dS.
    
    Args:
        mesh: FEniCS Mesh object
        mf: MeshFunction with boundary markers
        numbering: Dictionary with boundary marker numbers
        u: Displacement field (Function). If None, computes reference volume
    
    Returns:
        float: Cavity volume
    
    Note:
        Sign depends on orientation of N; check with a small test.
    """
    X = SpatialCoordinate(mesh) 
    N = FacetNormal(mesh)

    if u is not None:
        I = Identity(3)  # the identity matrix
        F = I + grad(u)  # the deformation gradient
        J = det(F)
        vol_form = (-1.0/3.0) * dot(X + u, J * inv(F).T * N)
    else:
        vol_form = (-1.0/3.0) * dot(X, N)
    
    ds = Measure('ds', domain=mesh, subdomain_data=mf)
    return assemble(vol_form * ds(numbering["ENDO"]))


def compute_reference_volume(mesh, boundary_markers, numbering):
    """
    Compute initial cavity volume of endocardial surface.
    
    Args:
        mesh: FEniCS Mesh object
        boundary_markers: MeshFunction with boundary markers
        numbering: Dictionary with boundary marker numbers
    
    Returns:
        float: Reference volume (absolute value)
    """
    n_mesh = FacetNormal(mesh)
    ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    X = SpatialCoordinate(mesh)
    vol = abs(assemble((-1.0/3.0) * dot(X, n_mesh) * ds(numbering['ENDO'])))
    return vol