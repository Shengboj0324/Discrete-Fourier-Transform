"""
Discrete Fourier Transform Implementation

This module provides a clean, human-readable implementation of the Discrete Fourier Transform (DFT)
for educational and practical purposes. The implementation prioritizes clarity and correctness
over performance optimization.

Author: Generated for DFT Project
"""

import math
import cmath
from typing import List, Union


def discrete_fourier_transform(signal: List[Union[complex, float, int]]) -> List[complex]:
    """
    Compute the Discrete Fourier Transform of a sequence of complex numbers.
    
    The DFT transforms a sequence of N complex numbers into another sequence of N complex numbers
    according to the formula:
    
    X[k] = Σ(n=0 to N-1) x[n] * e^(-2πi * k * n / N)
    
    where:
    - x[n] is the input sequence
    - X[k] is the output sequence (frequency domain)
    - N is the length of the sequence
    - i is the imaginary unit
    
    Args:
        signal (List[Union[complex, float, int]]): Input sequence of numbers.
                Can contain complex numbers, floats, or integers.
    
    Returns:
        List[complex]: The DFT of the input sequence. Each element represents
                      the amplitude and phase of a frequency component.
    
    Raises:
        ValueError: If the input signal is empty.
        TypeError: If the input is not a list or contains invalid types.
    
    Example:
        >>> # Simple example with real numbers
        >>> signal = [1, 0, 0, 0]
        >>> result = discrete_fourier_transform(signal)
        >>> # Result should be [1+0j, 1+0j, 1+0j, 1+0j]
        
        >>> # Example with complex numbers
        >>> signal = [1+0j, 0+1j, -1+0j, 0-1j]
        >>> result = discrete_fourier_transform(signal)
    """
    # Input validation
    if not isinstance(signal, list):
        raise TypeError("Input signal must be a list")
    
    if len(signal) == 0:
        raise ValueError("Input signal cannot be empty")
    
    # Convert all input values to complex numbers for consistent processing
    try:
        complex_signal = [complex(x) for x in signal]
    except (TypeError, ValueError) as e:
        raise TypeError(f"All elements in signal must be convertible to complex numbers: {e}")
    
    N = len(complex_signal)
    dft_result = []
    
    # Compute DFT for each frequency bin k
    for k in range(N):
        # Initialize the sum for this frequency bin
        frequency_component = complex(0, 0)
        
        # Sum over all time samples n
        for n in range(N):
            # Calculate the complex exponential: e^(-2πi * k * n / N)
            angle = -2.0 * math.pi * k * n / N
            twiddle_factor = cmath.exp(1j * angle)
            
            # Multiply input sample by twiddle factor and add to sum
            frequency_component += complex_signal[n] * twiddle_factor
        
        dft_result.append(frequency_component)
    
    return dft_result


def inverse_discrete_fourier_transform(frequency_domain: List[Union[complex, float, int]]) -> List[complex]:
    """
    Compute the Inverse Discrete Fourier Transform (IDFT) of a frequency domain sequence.
    
    The IDFT transforms a sequence of N complex frequency components back to the time domain
    according to the formula:
    
    x[n] = (1/N) * Σ(k=0 to N-1) X[k] * e^(2πi * k * n / N)
    
    Args:
        frequency_domain (List[Union[complex, float, int]]): Frequency domain sequence.
    
    Returns:
        List[complex]: The time domain sequence reconstructed from the frequency components.
    
    Raises:
        ValueError: If the input sequence is empty.
        TypeError: If the input is not a list or contains invalid types.
    """
    # Input validation
    if not isinstance(frequency_domain, list):
        raise TypeError("Input frequency_domain must be a list")
    
    if len(frequency_domain) == 0:
        raise ValueError("Input frequency_domain cannot be empty")
    
    # Convert all input values to complex numbers
    try:
        complex_freq = [complex(x) for x in frequency_domain]
    except (TypeError, ValueError) as e:
        raise TypeError(f"All elements must be convertible to complex numbers: {e}")
    
    N = len(complex_freq)
    idft_result = []
    
    # Compute IDFT for each time sample n
    for n in range(N):
        # Initialize the sum for this time sample
        time_sample = complex(0, 0)
        
        # Sum over all frequency bins k
        for k in range(N):
            # Calculate the complex exponential: e^(2πi * k * n / N)
            # Note: positive sign for IDFT (opposite of DFT)
            angle = 2.0 * math.pi * k * n / N
            twiddle_factor = cmath.exp(1j * angle)
            
            # Multiply frequency component by twiddle factor and add to sum
            time_sample += complex_freq[k] * twiddle_factor
        
        # Scale by 1/N for proper normalization
        time_sample /= N
        idft_result.append(time_sample)
    
    return idft_result


def get_magnitude_spectrum(dft_result: List[complex]) -> List[float]:
    """
    Extract the magnitude spectrum from DFT results.
    
    Args:
        dft_result (List[complex]): Output from discrete_fourier_transform()
    
    Returns:
        List[float]: Magnitude of each frequency component
    """
    return [abs(component) for component in dft_result]


def get_phase_spectrum(dft_result: List[complex]) -> List[float]:
    """
    Extract the phase spectrum from DFT results.
    
    Args:
        dft_result (List[complex]): Output from discrete_fourier_transform()
    
    Returns:
        List[float]: Phase (in radians) of each frequency component
    """
    return [cmath.phase(component) for component in dft_result]


if __name__ == "__main__":
    # Simple demonstration
    print("Discrete Fourier Transform Demo")
    print("=" * 40)
    
    # Example 1: Simple impulse signal
    impulse = [1, 0, 0, 0]
    print(f"Input signal: {impulse}")
    dft_impulse = discrete_fourier_transform(impulse)
    print(f"DFT result: {[f'{x:.3f}' for x in dft_impulse]}")
    print()
    
    # Example 2: Verify with IDFT
    reconstructed = inverse_discrete_fourier_transform(dft_impulse)
    print(f"IDFT result: {[f'{x:.3f}' for x in reconstructed]}")
    print("(Should match original input)")
