# Apodized Bragg Grating

Bragg Grating
A Bragg Grating is a periodic structure embedded in a waveguide, typically used in optical and photonic systems to selectively reflect certain wavelengths of light while transmitting others. The concept is based on periodic modulation of the refractive index of a material, which causes constructive interference of light at specific wavelengths, effectively filtering and reflecting those wavelengths.

To calculate the scattering parameters of a Bragg Grating implemented on a Coplanar Strip Transmission Line (CPS TL), the first step is to compute the characteristic impedance for each section of the line. The method we use for this calculation is outlined below:

### Coplanar Strip Transmission Line (CPS TL) Analysis

The analysis of **Coplanar Strip Transmission Line (CPS TL)** can be conducted using a quasi-static approximation. This approximation allows for the estimation of the capacitance C  and characteristic impedance Z_CPS of the CPS TL, taking into account factors such as the thickness of the metal t , the finite thickness of the substrate h, and the relative permittivity epsilon_r.

#### Effective Permittivity Calculation

To determine the effective permittivity of the CPS, the calculation assumes that the total capacitance per unit length of the CPS is the combined sum of the partial capacitance arising from free space (without a dielectric substrate) and the partial capacitance resulting from the substrate with a thickness of H and relative permittivity epsilon_r - 1:

$$
C = C_a + C_s \tag{1}
$$

Where:
- \( C_a \) denotes the capacitance per unit length of the TL resulting from the presence of air without a dielectric substrate.
- \( C_s \) represents the capacitance per unit length of the TL attributed to the dielectric substrate with a relative permittivity of epsilon_r - 1.

#### Capacitance Calculation

To calculate \( C_a \) and \( C_s \), a series of conformal mappings specific to a CPS with a substrate of finite thickness is required. In this analysis, an imaginary plane positioned between the striplines is treated as an electric wall, while the interface between the substrate and air is regarded as a magnetic wall.

##### Capacitance from Air (\( C_a \))

The capacitance \( C_a \) between the striplines in the CPS, resulting from the presence of air, can be viewed as two identical capacitors connected in series between the striplines and the electric wall, and is given by:

$$
C_a = \varepsilon_0 \frac{K(k_1')}{K(k_1)} \tag{2}
$$

Where:

$$
k_1 = \frac{S}{S + 2W} \tag{3}
$$

And \( k_1' = \sqrt{1 - k_1^2} \).

##### Capacitance from Substrate (\( C_s \))

The capacitance resulting from the dielectric with a relative permittivity of epsilon_r - 1 is expressed as:

$$
C_s = \varepsilon_0 \frac{\varepsilon_r - 1}{2} \frac{K(k_2')}{K(k_2)} \tag{4}
$$

Where:

$$
k_2 = \frac{\tanh\left(\frac{\pi S}{4H}\right)}{\tanh\left(\frac{\pi (S + W)}{4H}\right)} \tag{5}
$$

And:

$$
k_2' = \sqrt{1 - k_2^2} \tag{6}
$$

#### Effective Permittivity

Using the equation for the total capacitance \( C = C_a + C_s \), we can derive the equation for the effective permittivity epsilon_eff:

$$
\varepsilon_{eff} = 1 + \frac{\varepsilon_r - 1}{2} \frac{K(k_1)}{K(k_1')} \frac{K(k_2')}{K(k_2)} \tag{7}
$$

#### Characteristic Impedance

The characteristic impedance of the CPS TL is given by:

$$
Z_{CPS} = \frac{120 \pi}{\sqrt{\varepsilon_{eff}}} \frac{K(k_1)}{K(k_1')} \tag{8}
$$

For more details and equations, please refer to the sources cited in the text: [Gevorgian et al. (2001)], [Gevorgian et al. (2003)].

Additionally, to compute the S-parameters, we employ the ABCD matrix method. As a result, both the magnitude and phase of the S21 and S11 parameters are calculated and plotted.
