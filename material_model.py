"""
Material constitutive models for cardiac tissue.
Implements Holzapfel-Ogden orthotropic hyperelastic model.
"""
# from dolfin import Identity, grad, det, inv, inner, tr, outer, exp, pow, conditional, ge
from dolfin import *
from ufl import variable


def subplus(x):
    """
    Ramp function: returns x if x >= 1, otherwise returns 1.
    Used for anisotropic invariants in constitutive model.
    """
    return conditional(ge(x, 1), x, 1.0)


class CardiacMaterial:
    """
    Holzapfel-Ogden orthotropic hyperelastic material model for cardiac tissue.
    """
    
    def __init__(self, material_params, fibers):
        """
        Initialize material model.
        
        Args:
            material_params: MaterialParameters object from config
            fibers: List of [fiber, sheet, normal] direction Functions
        """
        self.params = material_params
        self.f_0 = fibers[0]  # Fiber direction
        self.s_0 = fibers[1]  # Sheet direction
        self.n_0 = fibers[2]  # Normal direction
    
    def compute_kinematics(self, u):
        """
        Compute kinematic quantities from displacement field.
        
        Args:
            u: Displacement field (Function)
        
        Returns:
            dict: Dictionary containing F, C, B, J, I_1, and deformed fiber directions
        """
        d = u.geometric_dimension()
        I = Identity(d)
        F = I + grad(u)
        F = variable(F)
        C = F.T * F
        B = F * F.T
        
        I_1 = tr(C)
        J = det(F)
        
        # Deformed fiber directions
        f = F * self.f_0
        s = F * self.s_0
        n = F * self.n_0
        
        # Anisotropic invariants
        I_4f = inner(C * self.f_0, self.f_0)
        I_4s = inner(C * self.s_0, self.s_0)
        I_4n = inner(C * self.n_0, self.n_0)
        I8fs = inner(C * self.f_0, self.s_0)
        
        return {
            'F': F, 'C': C, 'B': B, 'J': J, 'I_1': I_1,
            'f': f, 's': s, 'n': n,
            'I_4f': I_4f, 'I_4s': I_4s, 'I_4n': I_4n, 'I8fs': I8fs
        }
    
    def passive_cauchy_stress(self, kinematics, p):
        """
        Compute passive Cauchy stress tensor.
        
        Args:
            kinematics: Dictionary from compute_kinematics
            p: Hydrostatic pressure (for incompressibility)
        
        Returns:
            Cauchy stress tensor (UFL expression)
        """
        B = kinematics['B']
        I_1 = kinematics['I_1']
        f = kinematics['f']
        s = kinematics['s']
        I_4f = kinematics['I_4f']
        I_4s = kinematics['I_4s']
        I8fs = kinematics['I8fs']
        I = Identity(3)
        
        a = self.params.a
        a_f = self.params.a_f
        a_s = self.params.a_s
        a_fs = self.params.a_fs
        b = self.params.b
        b_f = self.params.b_f
        b_s = self.params.b_s
        b_fs = self.params.b_fs
        
        sigma = (a * exp(b * (I_1 - 3)) * B
                 + 2 * a_f * (subplus(I_4f) - 1) * exp(b_f * pow(subplus(I_4f) - 1, 2)) * outer(f, f)
                 + 2 * a_s * (subplus(I_4s) - 1) * exp(b_s * pow(subplus(I_4s) - 1, 2)) * outer(s, s)
                 + a_fs * I8fs * exp(b_fs * pow(I8fs, 2)) * (outer(f, s) + outer(s, f))
                 - p * I)
        
        return sigma
    
    def active_cauchy_stress(self, kinematics, T_a):
        """
        Compute active Cauchy stress tensor.
        
        Args:
            kinematics: Dictionary from compute_kinematics
            T_a: Active tension (Constant or Expression)
        
        Returns:
            Active Cauchy stress tensor (UFL expression)
        """
        f = kinematics['f']
        s = kinematics['s']
        n = kinematics['n']
        eta = self.params.eta
        
        sigma_active = T_a * (outer(f, f) + eta * outer(s, s) + eta * outer(n, n))
        return sigma_active
    
    def first_piola_kirchhoff(self, cauchy_stress, kinematics):
        """
        Convert Cauchy stress to First Piola-Kirchhoff stress.
        
        Args:
            cauchy_stress: Cauchy stress tensor
            kinematics: Dictionary from compute_kinematics
        
        Returns:
            First Piola-Kirchhoff stress tensor
        """
        J = kinematics['J']
        F = kinematics['F']
        return J * cauchy_stress * inv(F).T