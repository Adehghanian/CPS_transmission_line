
""" ............................
    .     --------------------        .
    .     --------------------        .
/                    S
    .     --------------------        .
    .     --------------------        .
    ............................
"""
import logging
import numpy as np
import scattering_calculations as SC

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def cosine_apodization(i, N_cells, apodization_factor):
    """
    Apply cosine apodization to the values.

    Args:
    i (int): The current cell index.
    N_cells (int): The total number of cells.
    apodization_factor (float): Scaling factor for apodization.

    Returns:
    float: Apodized value for the current cell.
    """
    return (1 - np.cos(np.pi * i / N_cells)) * apodization_factor


def gaussian_apodization(i, N_cells, apodization_factor):
    """
    Apply Gaussian apodization to the values.

    Args:
    i (int): The current cell index.
    N_cells (int): The total number of cells.
    apodization_factor (float): Scaling factor for apodization.

    Returns:
    float: Apodized value for the current cell.
    """
    center = N_cells / 2
    sigma = N_cells / 6  # Standard deviation for the Gaussian function
    return np.exp(-0.5 * ((i - center) / sigma) ** 2) * apodization_factor


def hyperbolic_apodization(i, N_cells, apodization_factor):
    """
    Apply Hyperbolic apodization to the values.

    Args:
    i (int): The current cell index.
    N_cells (int): The total number of cells.
    apodization_factor (float): Scaling factor for apodization.

    Returns:
    float: Apodized value for the current cell.
    """
    return (1 - np.tanh(i / N_cells)) * apodization_factor


def bragg_grating_apodized(N_cells: int, S_section1: float, S_section2: float, 
                            W_section1: float, W_section2: float, apodization_type: str, 
                            apodization_factor: float) -> tuple:
    """
    Generates the S and W value lists for a Bragg Grating structure with selected apodization.

    The apodization is applied by adjusting the S and W values based on the selected apodization function.

    Args:
    N_cells (int): Number of cells in the Bragg Grating.
    S_section1 (float): S value for section 1.
    S_section2 (float): S value for section 2.
    W_section1 (float): W value for section 1.
    W_section2 (float): W value for section 2.
    apodization_type (str): Type of apodization ('cosine', 'gaussian', 'hyperbolic').
    apodization_factor (float): Factor to scale apodization. Affects the S and W values across cells.

    Returns:
    tuple: Two lists - S values and W values for the apodized Bragg Grating structure.
    """
    S_value_list = []
    W_value_list = []

    # Select apodization function based on user choice
    if apodization_type == 'cosine':
        apodization_func = cosine_apodization
    elif apodization_type == 'gaussian':
        apodization_func = gaussian_apodization
    elif apodization_type == 'hyperbolic':
        apodization_func = hyperbolic_apodization
    else:
        raise ValueError("Invalid apodization_type. Choose from 'cosine', 'gaussian', or 'hyperbolic'.")

    # Generate apodized values using the selected apodization function
    for i in range(N_cells):
        S_value_list.append(S_section1 + (S_section2 - S_section1) * apodization_func(i, N_cells, apodization_factor))
        W_value_list.append(W_section1 + (W_section2 - W_section1) * apodization_func(i, N_cells, apodization_factor))

    logging.info(f"\n Generated apodized S_value_list: {S_value_list}")
    logging.info(f"\n Generated apodized W_value_list: {W_value_list}")

    return S_value_list, W_value_list


def main():
    """
    Main function to calculate and plot the S parameters of the Bragg Grating.

    This function generates the necessary parameters for the Bragg Grating,
    calculates the scattering parameters using the `CS` module, and plots the results.
    It allows the choice of apodization function.
    """
    # Define constants and parameters
    N_cells = 10
    S_value = [70e-6, 80e-6]  # Section 1 and Section 2 spacing values in meters
    W_value = [40e-6, 35e-6]  # Section 1 and Section 2 width values in meters
    L = 100e-6  # Length values for each section in meters
    h = 1e-6  # Substrate thickness in meters
    frequencies = np.arange(1e11, 1.4e12, 1e9)  # Frequency range from 100 GHz to 1.3 THz

    if N_cells < 1:
        raise ValueError("The total number of cells must be 1 or more")
    if any(val <= 0 for val in [S_value[0], S_value[1], W_value[0], W_value[1], L, h]):
        raise ValueError("The input values must be positive float numbers")
    # Option to choose the apodization type
    apodization_type = 'gaussian'  # Choose from 'cosine', 'gaussian', 'hyperbolic'
    apodization_factor = 1.0  # Scaling factor for the apodization
    if apodization_type not in ['cosine', 'gaussian', 'hyperbolic']:
        raise ValueError("apodization_type must be one of cosine, gaussian, or hyperbolic")

    logging.info("Starting Bragg Grating calculation...")

    # Generate the apodized S and W value lists
    S_value_list, W_value_list = bragg_grating_apodized(
        N_cells, S_value[0], S_value[1], W_value[0], W_value[1], apodization_type, apodization_factor
    )

    logging.info(f"\n Bragg Grating generated with {N_cells} cells and {apodization_type} apodization.")

    # Calculate the S parameters using the SC module
    S11_values_dB, S21_values_dB, S11_phases, S21_phases = SC.calculate_S_param(
        frequencies, S_value_list, W_value_list, L, h
    )

    # Plot the calculated S parameters
    SC.plot_S_parameters(frequencies, S11_values_dB, S21_values_dB, S11_phases, S21_phases)

    logging.info("\n Bragg Grating S parameters calculated and plotted.")

    TL = np.array([[1+1j, 0], [0, 1+1j]])
    Z0 = 200
    S11, S21 = SC.calculate_scattering_parameters(TL, Z0)
    print(S11, S21)


if __name__ == '__main__':
    main()
