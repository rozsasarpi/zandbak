import numpy as np

from zandbak.beams.simply_supported import (
    clamped_clamped_beam_under_point_force,
    hinged_clamped_beam_under_point_force,
    hinged_hinged_beam_under_point_force,
)


def test_hinged_hinged_beam_under_point_force():
    span_length = 10
    flexural_stiffness = 3
    load_intensity = 2
    load_position_from_left_end = span_length / 2
    position_from_left_end = load_position_from_left_end

    d_x_expected = 1 / 48 * load_intensity * span_length ** 3 / flexural_stiffness
    m_x_expected = 1 / 4 * load_intensity * span_length

    d_x, _, m_x = hinged_hinged_beam_under_point_force(
        span_length=span_length,
        flexural_stiffness=flexural_stiffness,
        load_intensity=load_intensity,
        load_position_from_left_end=load_position_from_left_end,
        position_from_left_end=position_from_left_end,
    )

    np.testing.assert_almost_equal(d_x, d_x_expected)
    np.testing.assert_almost_equal(m_x, m_x_expected)


def test_hinged_clamped_beam_under_point_force():
    span_length = 10
    flexural_stiffness = 3
    load_intensity = 2
    load_position_from_left_end = span_length / 2
    position_from_left_end = load_position_from_left_end

    d_x_expected = 7 / 768 * load_intensity * span_length ** 3 / flexural_stiffness
    m_x_expected = 5 / 32 * load_intensity * span_length

    d_x, _, m_x = hinged_clamped_beam_under_point_force(
        span_length=span_length,
        flexural_stiffness=flexural_stiffness,
        load_intensity=load_intensity,
        load_position_from_hinged_end=load_position_from_left_end,
        position_from_left_end=position_from_left_end,
    )

    np.testing.assert_almost_equal(d_x, d_x_expected)
    np.testing.assert_almost_equal(m_x, m_x_expected)


def test_clamped_clamped_beam_under_point_force():
    span_length = 10
    flexural_stiffness = 3
    load_intensity = 2
    load_position_from_left_end = span_length / 2
    position_from_left_end = load_position_from_left_end

    d_x_expected = 1 / 192 * load_intensity * span_length ** 3 / flexural_stiffness
    m_x_expected = 1 / 8 * load_intensity * span_length

    d_x, _, m_x = clamped_clamped_beam_under_point_force(
        span_length=span_length,
        flexural_stiffness=flexural_stiffness,
        load_intensity=load_intensity,
        load_position_from_left_end=load_position_from_left_end,
        position_from_left_end=position_from_left_end,
    )

    np.testing.assert_almost_equal(d_x, d_x_expected)
    np.testing.assert_almost_equal(m_x, m_x_expected)
