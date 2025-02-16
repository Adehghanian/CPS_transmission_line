import pytest
import numpy as np

from scattering_calculations import calculate_Z0, calculate_S_param, calculate_matrix, calculate_scattering_parameters

# Test data
frequencies = np.array([0.5e12, 1e12, 1.5e12])  # 0.5 THz, 1 THz, 1.5 THz
S_value = np.array([70e-6, 20e-6, 500e-6])  # Spacing between conductors
W_value = np.array([45e-6, 5e-6, 200e-6])  # Width of the conductors
L = 10e-6  # Length of transmission line in meters
h = 1e-6  # Thickness of the substrate


# Test different values for S_value, W_value, and h
@pytest.mark.parametrize(
    "S_value, W_value, h, expected_Z0_min, expected_Z0_max",
    [
        ([70e-6, 20e-6, 70e-6], [45e-6, 5e-6, 45e-6], 1e-6, 0, 500),
        ([70e-6, 20e-6, 200e-6], [45e-6, 5e-6, 100e-6], 100e-6, 0, 500),
        ([70e-6, 10e-6], [45e-6, 5e-6], 1e-6, 0, 1000),
        ([70e-5, 10e-5], [10e-4, 5e-6], 2e-6, 0, 1000),
        ([1e-5], [5e-4], 5e-4, 0, 5000),
    ]
)

def test_calculate_Z0_pass(S_value, W_value, h, expected_Z0_min, expected_Z0_max):
    eps_eff, Z0 = calculate_Z0(S_value, W_value, h)
    
    # Assertions
    assert eps_eff is not None, "Effective permittivity should not be None"
    assert Z0 is not None, "Characteristic impedance should not be None"
    assert isinstance(eps_eff, np.ndarray), "eps_eff should be an array"
    assert np.all(Z0 >= expected_Z0_min), f"Z0 values should be equal or greater than {expected_Z0_min}"
    assert np.all(Z0 < expected_Z0_max), f"Z0 values should be less than {expected_Z0_max}"
    #Check that eps_eff is between expected range (0 to 10)
    assert np.all(eps_eff >= 0), f"eps_eff values should be greater than or equal to 1"
    assert np.all(eps_eff <= 10), f"eps_eff values should be less than or equal to 10"

@pytest.mark.xfail(reason="Only positive Z0 values are valid")
@pytest.mark.parametrize(
    "S_value, W_value, h, expected_Z0_min, expected_Z0_max",
    [
        (np.array([2e-6, 3e-6, 4e-6]), np.array([1e-6, 1.5e-6, 2e-6]), 1e-6, -100, 0),  # Invalid negative Z0 case
        (np.array([1e-6, 2e-6]), np.array([0.5e-6, 1e-6]), 1e-6, 0, 0),             # Invalid negative Z0 case
    ]
)
def test_calculate_Z0_fail(S_value, W_value, h, expected_Z0_min, expected_Z0_max):
    eps_eff, Z0 = calculate_Z0(S_value, W_value, h)
    
    # Assertions
    assert eps_eff is not None, "Effective permittivity should not be None"
    assert Z0 is not None, "Characteristic impedance should not be None"
    assert isinstance(eps_eff, np.ndarray), "eps_eff should be an array"
    assert isinstance(Z0, np.ndarray), "Z0 should be an array"
    assert np.all(Z0 > expected_Z0_min), f"Z0 values should be greater than {expected_Z0_min}"
    assert np.all(Z0 < expected_Z0_max), f"Z0 values should be less than {expected_Z0_max}"



# Test calculate_S_param
def test_calculate_S_param():
    S11_dB, S21_dB, S11_phases, S21_phases = calculate_S_param(frequencies, S_value, W_value, L, h)
    assert S11_dB is not None, "S11_dB should not be None"
    assert S21_dB is not None, "S21_dB should not be None"
    assert S11_phases is not None, "S11_phases should not be None"
    assert S21_phases is not None, "S21_phases should not be None"
    assert len(S11_dB) == len(frequencies), "S11_dB length should match frequencies"
    assert len(S21_dB) == len(frequencies), "S21_dB length should match frequencies"
    assert len(S11_phases) == len(frequencies), "S11_phases length should match frequencies"
    assert len(S21_phases) == len(frequencies), "S21_phases length should match frequencies"

# Test different values for beta and Z0
@pytest.mark.parametrize(
    "beta, Z0, expected_shape, expected_iscomplex",
    [
        (2000, 2000, (2, 2), True),  # Regular test case
        (10000, 750, (2, 2), True),  # Different beta and Z0
        (10000, 100, (2, 2), True), # Larger Z0
        (2000, 250, (2, 2), True),  # Smaller Z0
    ]
)
def test_calculate_matrix(beta, Z0, expected_shape, expected_iscomplex):
    matrix = calculate_matrix(beta, L, Z0)
    
    # Assertions
    assert matrix is not None, "Matrix should not be None"
    assert matrix.shape == expected_shape, f"Matrix shape should be {expected_shape}"
    assert np.iscomplexobj(matrix) == expected_iscomplex, "Matrix should contain complex values"


# Test case for valid input and calculation (no error expected)
@pytest.mark.parametrize(
    "TL, Z0, expected_S11_S21_range",
    [
        (np.array([[1+1j, 0], [0, 1+1j]], dtype=complex), 200, (0, 1)),  # Normal case
        (np.array([[1+1j, 2], [3, 4+1j]], dtype=complex), 50, (0, 1)),  # Another valid case
        (np.array([[1+1j, 0], [1j, 1+1j]], dtype=complex), 1000, (0, 1))
    ]
)
def test_calculate_scattering_parameters_pass(TL, Z0, expected_S11_S21_range):
    S11, S21 = calculate_scattering_parameters(TL, Z0)
    
    # Assertions for valid cases (when no error is expected)
    assert S11 is not None, "S11 should not be None"
    assert S21 is not None, "S21 should not be None"
    assert isinstance(S11, complex), "S11 should be a complex number"
    assert isinstance(S21, complex), "S21 should be a complex number"
    assert np.abs(S11) <= expected_S11_S21_range[1], "S11 magnitude should be <= 1"
    assert np.abs(S21) <= expected_S11_S21_range[1], "S21 magnitude should be <= 1"
    assert np.abs(S11) >= expected_S11_S21_range[0], "S11 magnitude should be >= 0"
    assert np.abs(S21) >= expected_S11_S21_range[0], "S21 magnitude should be >= 0"

@pytest.mark.xfail
def test_calculate_scattering_parameters_fail():
    TL = np.array([[1, 2], [3, 4]], dtype=complex)  # Sample 2x2 transmission line matrix
    Z0 = -50  # Invalid negative Z0
    S11, S21 = calculate_scattering_parameters(TL, Z0)
    assert S11 is not None, "S11 should not be None"
    assert S21 is not None, "S21 should not be None"
    assert isinstance(S11, complex), "S11 should be a complex number"
    assert isinstance(S21, complex), "S21 should be a complex number"
    assert np.abs(S11) <= 1, "S11 magnitude should be <= 1"
    assert np.abs(S21) <= 1, "S21 magnitude should be <= 1"


# Test logging for errors
def test_logging_for_errors():
    # Test for calculate_Z0 with invalid inputs
    S_value_invalid = [None, 2e-6]  # Invalid S_value
    W_value_invalid = [None, 1e-6]  # Invalid W_value
    eps_eff, Z0 = calculate_Z0(S_value_invalid, W_value_invalid, h)
    assert eps_eff is None, "eps_eff should be None for invalid inputs"
    assert Z0 is None, "Z0 should be None for invalid inputs"

if __name__ == "__main__":
    pytest.main()
