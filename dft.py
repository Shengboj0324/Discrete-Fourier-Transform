import math
import cmath
from typing import List, Union, Sequence

def discrete_fourier_transform(signal: Sequence[Union[complex, float, int]]) -> List[complex]:
    """
    Transform time domain signal to frequency domain.
    Formula: X[k] = Σ(n=0 to N-1) x[n] * e^(-2πi * k * n / N)
    """
    # Check if input is a sequence (list, tuple, etc.)
    if not hasattr(signal, '__iter__') or isinstance(signal, (str, bytes)):
        raise TypeError("Input signal must be a sequence of numbers")

    if len(signal) == 0:
        raise ValueError("Input signal cannot be empty")

    # Convert everything to complex numbers so we can handle them uniformly
    try:
        complex_signal = [complex(x) for x in signal]
    except (TypeError, ValueError) as e:
        raise TypeError(f"All elements in signal must be convertible to complex numbers: {e}")

    N = len(complex_signal)
    dft_result = []

    # Go through each frequency bin
    for k in range(N):
        # Start with zero for this frequency
        frequency_component = complex(0, 0)

        # Add up contributions from each time sample
        for n in range(N):
            # Calculate the rotation factor
            angle = -2.0 * math.pi * k * n / N
            twiddle_factor = cmath.exp(1j * angle)

            # Multiply signal by rotation and add to sum
            frequency_component += complex_signal[n] * twiddle_factor

        # Store this frequency component
        dft_result.append(frequency_component)

    return dft_result


def inverse_discrete_fourier_transform(frequency_domain: Sequence[Union[complex, float, int]]) -> List[complex]:
    """
    Transform frequency domain back to time domain.
    Formula: x[n] = (1/N) * Σ(k=0 to N-1) X[k] * e^(2πi * k * n / N)
    """
    # Check if input is a sequence
    if not hasattr(frequency_domain, '__iter__') or isinstance(frequency_domain, (str, bytes)):
        raise TypeError("Input frequency_domain must be a sequence of numbers")

    # Can't transform empty data
    if len(frequency_domain) == 0:
        raise ValueError("Input frequency_domain cannot be empty")

    # Convert to complex numbers
    try:
        complex_freq = [complex(x) for x in frequency_domain]
    except (TypeError, ValueError) as e:
        raise TypeError(f"All elements must be convertible to complex numbers: {e}")

    N = len(complex_freq)
    idft_result = []

    # Go through each time sample
    for n in range(N):
        # Start with zero for this time point
        time_sample = complex(0, 0)

        # Add up contributions from each frequency bin
        for k in range(N):
            # Calculate rotation factor (positive angle for inverse)
            angle = 2.0 * math.pi * k * n / N
            twiddle_factor = cmath.exp(1j * angle)

            # Multiply frequency by rotation and add to sum
            time_sample += complex_freq[k] * twiddle_factor

        # Divide by N to get the right scale
        time_sample /= N
        idft_result.append(time_sample)

    return idft_result


def get_magnitude_spectrum(dft_result: List[complex]) -> List[float]:
    """Get the magnitude of each frequency component."""
    return [abs(component) for component in dft_result]


def get_phase_spectrum(dft_result: List[complex]) -> List[float]:
    """Get the phase of each frequency component in radians."""
    return [cmath.phase(component) for component in dft_result]


if __name__ == "__main__":
    # Quick demo
    print("Discrete Fourier Transform Demo")
    print("=" * 40)

    # Test with impulse signal
    impulse = [1, 0, 0, 0]
    print(f"Input signal: {impulse}")
    dft_impulse = discrete_fourier_transform(impulse)
    print(f"DFT result: {[f'{x:.3f}' for x in dft_impulse]}")
    print()

    # Check if we can get back the original
    reconstructed = inverse_discrete_fourier_transform(dft_impulse)
    print(f"IDFT result: {[f'{x:.3f}' for x in reconstructed]}")
    print("(Should match original input)")
