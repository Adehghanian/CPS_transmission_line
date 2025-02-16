import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk
import logging
import matplotlib.colors as mcolors

# Configure logging for debugging purposes
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def calculate_Z0(S_value, W_value, h=1e-6):
    """
    Calculate the characteristic impedance (Z0) based on CPS parameters.
    
    Args:
    - S_value (float or array): Spacing between the conductors.
    - W_value (float or array): Width of the conductors.
    - h (float, optional): Thickness of the substrate. Default is 1e-6.

    Returns:
    - tuple: Effective permittivity (eps_eff) and calculated characteristic impedance (Z0).
    """
    try:
        S_value = np.array(S_value)
        W_value = np.array(W_value)
        
        K_value = S_value / (S_value + 2 * W_value)
        KP_value = np.sqrt(1 - K_value**2)
        K1_value = np.sinh(np.pi * S_value / (4 * h)) / np.sinh(np.pi * (S_value + 2 * W_value) / (4 * h))
        K1p_value = np.sqrt(1 - K1_value**2)
        eps_eff = 1 + ((1.3 - 1) / 2) * (ellipk(KP_value) / ellipk(K_value)) * (ellipk(K1_value) / ellipk(K1p_value))
        Z0_value = 120 * np.pi * ellipk(K_value) / (np.sqrt(eps_eff) * ellipk(KP_value))
        
        logging.debug(f"Calculated Z0: {Z0_value}, eps_eff: {eps_eff}")
        return eps_eff, Z0_value
    except Exception as e:
        logging.error(f"Error calculating Z0: {e}")
        return None, None


def calculate_S_param(frequencies, S_value, W_value, L, h):
    """
    Calculate S21 and S11 parameters in dB and their phases for given frequencies.
    
    Args:
    - frequencies (array): Array of frequencies (in Hz) at which to calculate the S parameters.
    - S_value (array): Spacing between the conductors.
    - W_value (array): Width of the conductors.
    - L (float): Length of the transmission line.
    - h (float): Thickness of the substrate.

    Returns:
    - tuple: Arrays of S11 and S21 magnitudes (dB), and their respective phases.
    """
    try:
        # Initialize lists to store results
        S21_values_dB, S11_values_dB, S21_phases, S11_phases = [], [], [], []
        
        # Log the frequency range and parameters
        logging.debug(f"Calculating S parameters for frequencies: {frequencies}, S_value: {S_value}, W_value: {W_value}, L: {L}, h: {h}")
        
        c = 3e8  # Speed of light in m/s
        eps_eff, Z0_values = calculate_Z0(S_value, W_value, h)
        
        if eps_eff is None:
            logging.error("Failed to calculate Z0, returning None for S parameters.")
            return None, None, None, None

        # Loop through frequencies to calculate the S parameters
        for F in frequencies:
            beta = (2 * np.pi * F * np.sqrt(eps_eff)) / c  # Phase constant
            TL = np.eye(2, dtype=complex)  # Transmission line matrix
            logging.debug(f"Frequency: {F}, Phase constant (beta): {beta}")
            a = len(W_value)
            # Loop through the conductor widths and calculate matrix for each section
            for i in range(len(W_value)):
                Z0 = Z0_values[i]
                beta_section = beta[i]
                section_length = L
                T1 = calculate_matrix(beta_section, section_length, Z0)
                TL = np.dot(T1, TL)

            # Calculate scattering parameters
            S11, S21 = calculate_scattering_parameters(TL, Z0 if Z0 != 0 else 1)

            # Convert to dB and phase
            S21_dB = 20 * np.log10(np.abs(S21))
            S11_dB = 20 * np.log10(np.abs(S11))
            S21_phase = np.angle(S21, deg=True)
            S11_phase = np.angle(S11, deg=True)

            # Append results
            S21_values_dB.append(S21_dB)
            S11_values_dB.append(S11_dB)
            S21_phases.append(S21_phase)
            S11_phases.append(S11_phase)

            logging.debug(f"At frequency {F} Hz: S11_dB = {S11_dB}, S21_dB = {S21_dB}, S11_phase = {S11_phase}, S21_phase = {S21_phase}")

        # Return results as numpy arrays
        return np.array(S11_values_dB), np.array(S21_values_dB), np.array(S11_phases), np.array(S21_phases)

    except Exception as e:
        logging.error(f"Error calculating S parameters: {e}")
        return None, None, None, None



def calculate_matrix(beta, L, Z0):
    """
    Calculate the transmission line matrix for a given section.

    Args:
    - beta (float): Phase constant.
    - L (float): Length of the section.
    - Z0 (float): Characteristic impedance.

    Returns:
    - np.ndarray: The 2x2 transmission line matrix for the section.
    """
    try:
        # Log the input values
        logging.debug(f"Calculating transmission line matrix with beta={beta}, L={L}, Z0={Z0}")

        # Calculate matrix components
        A = np.cos(beta * L)
        B = Z0 * 1j * np.sin(beta * L)
        C = 1j * np.sin(beta * L) / Z0 if Z0 != 0 else 0
        D = np.cos(beta * L)

        # Log the calculated values
        logging.debug(f"Calculated components: A={A}, B={B}, C={C}, D={D}")

        # Return the matrix
        matrix = np.array([[A, B], [C, D]]).reshape(2, 2)
        logging.debug(f"Transmission line matrix: {matrix}")
        return matrix

    except Exception as e:
        logging.error(f"Error calculating matrix: {e}")
        return None



import logging
import numpy as np

# Setup logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_scattering_parameters(TL, Z0):
    """
    Calculate S11 and S21 scattering parameters from transmission matrix.

    Args:
    - TL (np.ndarray): 2x2 transmission line matrix.
    - Z0 (float): Characteristic impedance.

    Returns:
    - tuple: S11 and S21 scattering parameters.
    """
    try:
        # Check if Z0 is zero to avoid division by zero
        if Z0 == 0:
            logging.debug(f"Z0 is zero, raising ZeroDivisionError")
            raise ZeroDivisionError("Z0 cannot be zero.")
        elif Z0 < 0:
            raise ValueError("negative Z0")
        # Log the input values
        logging.info(f"Calculating scattering parameters with Z0={Z0} and TL={TL}")

        # Calculate numerator and denominator for S11 and S21
        numerator_S11 = TL[0, 0] + TL[0, 1] / Z0 - TL[1, 0] * Z0 - TL[1, 1]
        denominator_S11_S21 = TL[0, 0] + TL[0, 1] / Z0 + TL[1, 0] * Z0 + TL[1, 1]


        # Calculate S11 and S21
        S11 = numerator_S11 / denominator_S11_S21
        S21 = 2 / denominator_S11_S21

        # Log the final results
        logging.info(f"S11: {S11}, S21: {S21}")

        return S11, S21

    except ZeroDivisionError as e:
        logging.error(f"Error calculating scattering parameters: {e}")
        return None, None
    except ValueError as e:
        logging.error(f"Error calculating scattering parameters: {e}")
        return None, None
    except Exception as e:
        logging.error(f"Unexpected error calculating scattering parameters: {e}")
        return None, None




def plot_S_parameters(frequencies, S11_values_dB, S21_values_dB, S11_phases, S21_phases):
    """
    Plot S11 and S21 parameters in dB and their phases.

    Args:
    - frequencies (array): Array of frequencies at which S parameters are calculated.
    - S11_values_dB (array): Array of S11 values in dB.
    - S21_values_dB (array): Array of S21 values in dB.
    - S11_phases (array): Array of S11 phase values.
    - S21_phases (array): Array of S21 phase values.
    """
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    axs[0, 0].plot(frequencies, S11_values_dB, label='S11 (dB)')
    axs[0, 0].set_title('S11 (dB)')
    axs[0, 0].set_xlabel('Frequency (Hz)')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].grid(True)
    axs[0, 0].legend()

    axs[0, 1].plot(frequencies, S21_values_dB, label='S21 (dB)')
    axs[0, 1].set_title('S21 (dB)')
    axs[0, 1].set_xlabel('Frequency (Hz)')
    axs[0, 1].set_ylabel('Magnitude (dB)')
    axs[0, 1].grid(True)
    axs[0, 1].legend()

    axs[1, 0].plot(frequencies, S11_phases, label='S11 Phase')
    axs[1, 0].set_title('S11 Phase')
    axs[1, 0].set_xlabel('Frequency (Hz)')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].grid(True)
    axs[1, 0].legend()

    axs[1, 1].plot(frequencies, S21_phases, label='S21 Phase')
    axs[1, 1].set_title('S21 Phase')
    axs[1, 1].set_xlabel('Frequency (Hz)')
    axs[1, 1].set_ylabel('Phase (Degrees)')
    axs[1, 1].grid(True)
    axs[1, 1].legend()

    plt.tight_layout()
    plt.show()
