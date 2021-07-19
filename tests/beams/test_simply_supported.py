from zandbak.beams.simply_supported import hinged_hinged_beam_under_point_force


def test_hinged_hinged_beam_under_point_force():
    span_length = 10
    flexural_stiffness = 3
    load_intensity = 2
    load_position_from_left_end = span_length / 2
    position_from_left_end = load_position_from_left_end

    d_x_expected = load_intensity * span_length ** 3 / (48 * flexural_stiffness)
    m_x_expected = load_intensity * position_from_left_end / 2

    d_x, _, m_x = hinged_hinged_beam_under_point_force(
        span_length=span_length,
        flexural_stiffness=flexural_stiffness,
        load_intensity=load_intensity,
        load_position_from_left_end=load_position_from_left_end,
        position_from_left_end=position_from_left_end,
    )

    assert d_x == d_x_expected
    assert m_x == m_x_expected
