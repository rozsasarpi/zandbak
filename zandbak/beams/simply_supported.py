"""Simply supported beams.
"""
from typing import Iterable, Tuple, Union

import numpy as np

from zandbak import REAL_SCALAR, REAL_VECT


def hinged_hinged_beam_under_point_force(
    position_from_left_end: REAL_VECT,
    load_intensity: REAL_SCALAR,
    load_position_from_left_end: REAL_SCALAR,
    span_length: REAL_SCALAR,
    flexural_stiffness: REAL_VECT,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Deflection, rotation, and bending moment of a hinged-hinged single span beam with a
    single concentrated force.

    Reference: Figure 8;
    https://www.awc.org/pdf/codes-standards/publications/design-aids/AWC-DA6-BeamFormulas-0710.pdf

    Args:
        position_from_left_end: position(s) from the hinged end where the internal
            forces and displacements are calculated, can be a scalar and vector as
            well (`m` elements).
        load_intensity: point load intensity.
        load_position_from_left_end: point load position from the hinged end.
        span_length: span length.
        flexural_stiffness: elastic modulus * second moment of area, can be a scalar and
            vector as well (`n` elements).

    Returns:
        vertical translation, 2d numpy array of `n` x `m`.

        rotation, 2d numpy array of `n` x `m`.

        bending moment, 2d numpy array of `n` x `m`.

    Structural responses at `position_from_left_end`.
    """
    # ----------------------------------------------------
    # PRE-PROCESS
    # ----------------------------------------------------
    x = np.array(position_from_left_end).reshape((1, -1))
    p = load_intensity
    span = span_length
    a = load_position_from_left_end
    b = span - a
    ei = np.array(flexural_stiffness).reshape((-1, 1))

    m = x.size
    n = ei.size

    # ----------------------------------------------------
    # STRUCTURAL ANALYSIS
    # ----------------------------------------------------
    # bending moment
    m_a = p * a * b / span
    m_vec = np.array([0, m_a, 0])
    m_pos_vec = np.array([0, a, span])
    m_x_1 = np.interp(x=x.ravel(), xp=m_pos_vec, fp=m_vec)
    m_x = np.tile(m_x_1, (n, 1))

    # vertical translation (deflection)
    idx = x <= a
    idx_raveled = idx.ravel()
    d_x = np.empty((n, m))
    multiplier = 1 / (6 * ei * span)
    x_a = x[idx].reshape((1, -1))
    d_x[:, idx_raveled] = (
        p * b * np.matmul(multiplier, x_a * (span ** 2 - b ** 2 - x_a ** 2))
    )

    x_b = x[~idx].reshape((1, -1))
    d_x[:, ~idx_raveled] = (
        p
        * a
        * np.matmul(
            multiplier,
            (span - x_b) * (2 * span * x_b - x_b ** 2 - a ** 2),
        )
    )

    # rotation
    phi_x = np.empty((n, m))
    phi_x[:, idx_raveled] = (
        p * b * np.matmul(multiplier, span ** 2 - b ** 2 - 3 * x_a ** 2)
    )

    phi_x[:, ~idx_raveled] = (
        p
        * a
        * np.matmul(
            multiplier,
            2 * span ** 2 - 2 * span * x_b - 4 * span * x_b + 3 * x_b ** 2 + a ** 2,
        )
    )

    return d_x, phi_x, m_x


def hinged_hinged_beam_under_end_moment(
    position_from_left_end: Union[int, float, Iterable, np.ndarray],
    loaded_end: str,
    load_intensity: Union[int, float],
    span_length: Union[int, float],
    flexural_stiffness: Union[int, float, Iterable, np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """

    Reference: https://mechanicalc.com/reference/beam-deflection-tables

    Args:
        position_from_left_end: position(s) from the hinged end where the internal
            forces and displacements are calculated, can be a scalar and vector as
            well (`m` elements).
        load_intensity: point load (moment) intensity.
        loaded_end: "left" or "right".
        span_length: span length.
        flexural_stiffness: elastic modulus * second moment of area, can be a scalar and
            vector as well (`n` elements).

    Returns:
        vertical translation, 2d numpy array of `n` x `m`.

        rotation, 2d numpy array of `n` x `m`.

        bending moment, 2d numpy array of `n` x `m`.

    Structural responses at `position_from_left_end`.
    """
    # ----------------------------------------------------
    # PRE-PROCESS
    # ----------------------------------------------------
    if loaded_end not in ["left", "right"]:
        raise ValueError("`loaded_end` must be either 'left' or 'right'")

    x = np.array(position_from_left_end).reshape((1, -1))
    # phi_0 = load_intensity
    span = span_length
    ei = np.array(flexural_stiffness).reshape((-1, 1))

    if loaded_end == "right":
        x = np.sort(np.abs(x - span))

    n = ei.size
    # ----------------------------------------------------
    # STRUCTURAL ANALYSIS
    # ----------------------------------------------------

    # bending moment
    # m_0 = 3 * ei * phi_0 / l
    m_0 = load_intensity
    m_x_1 = m_0 * (span - x) / span
    m_x = np.tile(m_x_1, (n, 1))

    # vertical translation (deflection)
    multiplier = m_0 / (6 * ei * span)
    d_x = np.matmul(multiplier, x * (2 * span ** 2 - 3 * span * x + x ** 2))

    # rotation
    phi_x = np.matmul(multiplier, 2 * span ** 2 - 6 * span * x + 3 * x ** 2)

    if loaded_end == "right":
        d_x = np.fliplr(d_x)
        phi_x = np.fliplr(phi_x)
        m_x = np.fliplr(m_x)

    return d_x, phi_x, m_x


def hinged_clamped_beam_under_point_force(
    position_from_left_end: Union[int, float, Iterable, np.ndarray],
    load_intensity: Union[int, float],
    load_position_from_hinged_end: Union[int, float],
    span_length: Union[int, float],
    flexural_stiffness: Union[int, float, Iterable, np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """

    Reference: Figure 17;
    https://www.awc.org/pdf/codes-standards/publications/design-aids/AWC-DA6-BeamFormulas-0710.pdf

    Args:
        position_from_left_end: position(s) from the hinged end where the internal
            forces and displacements are calculated, can be a scalar and vector as
            well (`m` elements).
        load_intensity: point load intensity.
        load_position_from_hinged_end: point load position from the hinged end.
        span_length: span length.
        flexural_stiffness: elastic modulus * second moment of area, can be a scalar and
            vector as well (`n` elements).

    Returns:
        vertical translation, 2d numpy array of `n` x `m`.

        rotation, 2d numpy array of `n` x `m`.

        bending moment, 2d numpy array of `n` x `m`.

    Structural responses at `position_from_left_end`.
    """
    # ----------------------------------------------------
    # PRE-PROCESS
    # ----------------------------------------------------
    x = np.array(position_from_left_end).reshape((1, -1))
    p = load_intensity
    span = span_length
    a = load_position_from_hinged_end
    b = span - a
    ei = np.array(flexural_stiffness).reshape((-1, 1))

    m = x.size
    n = ei.size

    # ----------------------------------------------------
    # STRUCTURAL ANALYSIS
    # ----------------------------------------------------
    # vertical reaction force at the hinged support
    r_1 = p * b ** 2 / (2 * span ** 3) * (a + 2 * span)

    # bending moment
    m_1 = r_1 * a
    m_2 = p * a * b / (2 * span ** 2) * (a + span)
    m_vec = np.array([0, m_1, -m_2])
    m_pos_vec = np.array([0, a, span])
    m_x_1 = np.interp(x=x.ravel(), xp=m_pos_vec, fp=m_vec)
    m_x = np.tile(m_x_1, (n, 1))

    # vertical translation (deflection)
    idx = x <= a
    idx_raveled = idx.ravel()
    d_x = np.empty((n, m))
    multiplier = 1 / (12 * ei * span ** 3)
    x_a = x[idx].reshape((1, -1))
    d_x[:, idx_raveled] = (
        p
        * b ** 2
        * np.matmul(
            multiplier, x_a * (3 * a * span ** 2 - 2 * span * x_a ** 2 - a * x_a ** 2)
        )
    )

    x_b = x[~idx].reshape((1, -1))
    d_x[:, ~idx_raveled] = (
        p
        * a
        * np.matmul(
            multiplier,
            (span - x_b) ** 2
            * (3 * span ** 2 * x_b - a ** 2 * x_b - 2 * a ** 2 * span),
        )
    )

    # rotation
    phi_x = np.empty((n, m))
    phi_x[:, idx_raveled] = (
        p
        * b ** 2
        * np.matmul(
            multiplier, 3 * a * span ** 2 - 6 * span * x_a ** 2 - 3 * a * x_a ** 2
        )
    )

    phi_x[:, ~idx_raveled] = (
        p
        * a
        * np.matmul(
            multiplier,
            3 * (span - x_b) * (span ** 2 * (span - 3 * x_b) + a ** 2 * (span + x_b)),
        )
    )

    return d_x, phi_x, m_x


def clamped_clamped_beam_under_point_force(
    position_from_left_end: Union[int, float, Iterable, np.ndarray],
    load_intensity: Union[int, float],
    load_position_from_left_end: Union[int, float],
    span_length: Union[int, float],
    flexural_stiffness: Union[int, float, Iterable, np.ndarray],
):
    """
    Reference: Figure 25;
    https://www.awc.org/pdf/codes-standards/publications/design-aids/AWC-DA6-BeamFormulas-0710.pdf

    Args:
        position_from_left_end: position(s) from the hinged end where the internal
            forces and displacements are calculated, can be a scalar and vector as
            well (`m` elements).
        load_intensity: point load intensity.
        load_position_from_left_end: point load position from the hinged end.
        span_length: span length.
        flexural_stiffness: elastic modulus * second moment of area, can be a scalar and
            vector as well (`n` elements).

    Returns:
        vertical translation, 2d numpy array of `n` x `m`.

        rotation, 2d numpy array of `n` x `m`.

        bending moment, 2d numpy array of `n` x `m`.

    Structural responses at `position_from_left_end`.
    """
    # ----------------------------------------------------
    # PRE-PROCESS
    # ----------------------------------------------------
    x = np.array(position_from_left_end).reshape((1, -1))
    p = load_intensity
    span = span_length
    a = load_position_from_left_end
    b = span - a
    ei = np.array(flexural_stiffness).reshape((-1, 1))

    m = x.size
    n = ei.size

    # ----------------------------------------------------
    # STRUCTURAL ANALYSIS
    # ----------------------------------------------------
    # bending moment
    m_1 = p * a * b ** 2 / (span ** 2)
    m_2 = p * a ** 2 * b / (span ** 2)
    m_a = 2 * p * a ** 2 * b ** 2 / (span ** 3)
    m_vec = np.array([-m_1, m_a, -m_2])
    m_pos_vec = np.array([0, a, span])
    m_x_1 = np.interp(x=x.ravel(), xp=m_pos_vec, fp=m_vec)
    m_x = np.tile(m_x_1, (n, 1))

    # vertical translation (deflection) & rotation
    def deflection_rotation_left_to_load(x, a, b):
        idx = x <= a
        if any(~idx):
            print(f"{np.sum(~idx)} `x` values greater than `a` are ignored!")
            x = x[idx]
        multiplier = 1 / (6 * ei * span ** 3)
        x_a = x.reshape((1, -1))
        d_x_a = (
            p
            * b ** 2
            * np.matmul(multiplier, x_a ** 2 * (3 * a * span - 3 * a * x_a - b * x_a))
        )
        phi_x_a = (
            p
            * b ** 2
            * np.matmul(
                multiplier, 6 * a * span * x_a - 9 * a * x_a ** 2 - 3 * b * x_a ** 2
            )
        )
        return d_x_a, phi_x_a

    idx = x <= a
    idx_raveled = idx.ravel()
    d_x = np.empty((n, m))
    phi_x = np.empty((n, m))

    # left to the load
    d_x_a, phi_x_a = deflection_rotation_left_to_load(x=x[idx], a=a, b=b)
    d_x[:, idx_raveled] = d_x_a
    phi_x[:, idx_raveled] = phi_x_a

    # right to the load
    x_b_tr = np.sort(np.abs(x[~idx] - span))
    d_x_b_tr, phi_x_b_tr = deflection_rotation_left_to_load(x=x_b_tr, a=b, b=a)
    d_x_b = np.fliplr(d_x_b_tr)
    phi_x_b = -np.fliplr(phi_x_b_tr)
    d_x[:, ~idx_raveled] = d_x_b
    phi_x[:, ~idx_raveled] = phi_x_b

    return d_x, phi_x, m_x
