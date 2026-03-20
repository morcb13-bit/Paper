"""
b13phase — B13 Fractal Phase Library (BASE=3120)
=================================================

BASE = 3120 = 2⁴ × 3 × 5 × 13
            = 60 × 52
            = truncated icosahedron vertices × active states

Core modules:
  constants              — BASE, angular step constants, Z₁₃ coset sets
  level0_table           — COS_TABLE, SIN_TABLE (3120 entries)
  phase_digits           — exact radix-BASE digit arithmetic (arbitrary precision)
  phase_packed_u64       — packed 5-digit u64 arithmetic
  evaluator_proto        — fractal_cos_sin_proto (multi-digit phase evaluator)
  z13_structure          — Z₁₃ coset algebra, Fibonacci/Pell mod 13
  truncated_icosahedron  — vertex coordinates, central angles, Proposition D
  great_circle           — great circle analysis, φ-arithmetic, Observation R

Quick start:
  from b13phase import fractal_cos_sin_proto, digits_from_int
  digits = digits_from_int(312, 3)   # 36° in 3-digit BASE3120
  cx, cy = fractal_cos_sin_proto(digits)
"""

from .constants import (BASE, INT64_MAX,
                         STEP_Z5, STEP_Z12, STEP_H, STEP_Z13, STEP_MIN,
                         PHI_BASE, PHI_HALF_BASE,
                         Z13_ORDER, H_ORDER, Z5_ORDER,
                         VERTICES, PENTAGONS, HEXAGONS, ACTIVE,
                         F7, P7,
                         H_SET, H2_SET, H4_SET, Z13_STAR)

from .level0_table import COS_TABLE, SIN_TABLE, coset_of

from .phase_digits import (digits_zero, digits_from_int, digits_to_int,
                            digits_add, digits_inc, digits_sub,
                            digits_mul_scalar, digits_phase,
                            digits_negate, digits_half)

from .phase_packed_u64 import (pack_u64_digits, unpack_u64_digits,
                                add_u64_packed, sub_u64_packed,
                                negate_u64_packed,
                                BITS_PER_DIGIT, MASK, U64_DIGITS)

from .evaluator_proto import fractal_cos_sin_proto, norm_sq, norm_ratio

from .z13_structure import (coset_of, coset_index, is_forbidden,
                             z13_star_index, pentagon_index,
                             coset_pattern, fib_mod13, pell_mod13,
                             verify_4h_avoidance,
                             phi_residue_orbit, delta_s_orbit)

from .truncated_icosahedron import (PHI, R, R_SQ,
                                     vertex_coords_float,
                                     vertex_phase_index,
                                     vertex_cos_sin,
                                     all_vertices,
                                     inner_product,
                                     central_angle,
                                     classify_pair,
                                     observation_n_verify,
                                     opposite_face_angles)

from .great_circle import (great_circle_traversal,
                            theta_values_degrees,
                            phi_arithmetic_verify,
                            arc_length, central_angle_from_chord,
                            cos_theta_a_phi, cos_theta_5_phi, cos_theta_6_phi,
                            DENOM, H5, H6, H5_SQ, H6_SQ)

__version__ = "0.5.1"
__all__ = [
    # constants
    "BASE", "INT64_MAX",
    "STEP_Z5", "STEP_Z12", "STEP_H", "STEP_Z13", "STEP_MIN",
    "PHI_BASE", "PHI_HALF_BASE",
    "Z13_ORDER", "H_ORDER", "Z5_ORDER",
    "VERTICES", "PENTAGONS", "HEXAGONS", "ACTIVE",
    "F7", "P7", "H_SET", "H2_SET", "H4_SET", "Z13_STAR",
    # tables
    "COS_TABLE", "SIN_TABLE",
    # digit ops
    "digits_zero", "digits_from_int", "digits_to_int",
    "digits_add", "digits_inc", "digits_sub",
    "digits_mul_scalar", "digits_phase",
    "digits_negate", "digits_half",
    # packed u64
    "pack_u64_digits", "unpack_u64_digits",
    "add_u64_packed", "sub_u64_packed", "negate_u64_packed",
    "BITS_PER_DIGIT", "MASK", "U64_DIGITS",
    # evaluator
    "fractal_cos_sin_proto", "norm_sq", "norm_ratio",
    # Z₁₃
    "coset_of", "coset_index", "is_forbidden",
    "z13_star_index", "pentagon_index",
    "coset_pattern", "fib_mod13", "pell_mod13",
    "verify_4h_avoidance", "phi_residue_orbit", "delta_s_orbit",
    # truncated icosahedron
    "PHI", "R", "R_SQ",
    "vertex_coords_float", "vertex_phase_index", "vertex_cos_sin",
    "all_vertices", "inner_product", "central_angle",
    "classify_pair", "observation_n_verify", "opposite_face_angles",
    # great circle
    "great_circle_traversal", "theta_values_degrees",
    "phi_arithmetic_verify", "arc_length", "central_angle_from_chord",
    "cos_theta_a_phi", "cos_theta_5_phi", "cos_theta_6_phi",
    "DENOM", "H5", "H6", "H5_SQ", "H6_SQ",
]
