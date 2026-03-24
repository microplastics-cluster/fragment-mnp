from fragmentmnp.regime import classify_release_regime


def test_regime_classifier_returns_expected_keys():
    out = classify_release_regime(
        D_p=1e-16,
        K_pw=1e4,
        k_frag=0.01,
        radius_m=1e-4,
    )
    assert "tau_diff" in out
    assert "tau_frag" in out
    assert "ratio_diff_to_frag" in out
    assert "label" in out