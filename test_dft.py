"""
Test suite for the Discrete Fourier Transform implementation.

This module contains comprehensive tests to verify the correctness of the DFT implementation,
including edge cases, mathematical properties, and comparison with known results.
"""

import unittest
import math
import cmath
from dft import (
    discrete_fourier_transform, 
    inverse_discrete_fourier_transform,
    get_magnitude_spectrum,
    get_phase_spectrum
)


class TestDiscreteFourierTransform(unittest.TestCase):
    """Test cases for the DFT implementation."""
    
    def setUp(self):
        """Set up test fixtures with known signals and their expected DFT results."""
        # Tolerance for floating point comparisons
        self.tolerance = 1e-10
    
    def assertComplexListAlmostEqual(self, list1, list2, tolerance=None):
        """Helper method to compare lists of complex numbers with tolerance."""
        if tolerance is None:
            tolerance = self.tolerance
        
        self.assertEqual(len(list1), len(list2), "Lists must have the same length")
        
        for i, (a, b) in enumerate(zip(list1, list2)):
            self.assertAlmostEqual(a.real, b.real, delta=tolerance, 
                                 msg=f"Real part mismatch at index {i}")
            self.assertAlmostEqual(a.imag, b.imag, delta=tolerance, 
                                 msg=f"Imaginary part mismatch at index {i}")
    
    def test_impulse_signal(self):
        """Test DFT of an impulse signal (1, 0, 0, 0)."""
        signal = [1, 0, 0, 0]
        expected = [1+0j, 1+0j, 1+0j, 1+0j]
        
        result = discrete_fourier_transform(signal)
        self.assertComplexListAlmostEqual(result, expected)
    
    def test_constant_signal(self):
        """Test DFT of a constant signal."""
        signal = [1, 1, 1, 1]
        # For a constant signal, only DC component (k=0) should be non-zero
        expected = [4+0j, 0+0j, 0+0j, 0+0j]
        
        result = discrete_fourier_transform(signal)
        self.assertComplexListAlmostEqual(result, expected)
    
    def test_single_frequency_cosine(self):
        """Test DFT of a single frequency cosine wave."""
        # cos(2πn/4) for n = 0, 1, 2, 3
        signal = [1, 0, -1, 0]  # cos(0), cos(π/2), cos(π), cos(3π/2)
        
        result = discrete_fourier_transform(signal)
        
        # For a cosine at frequency k=1, we expect peaks at k=1 and k=3 (N-1)
        # with magnitude 2 each (due to symmetry)
        expected_magnitudes = [0, 2, 0, 2]
        actual_magnitudes = [abs(x) for x in result]
        
        for i, (expected, actual) in enumerate(zip(expected_magnitudes, actual_magnitudes)):
            self.assertAlmostEqual(actual, expected, delta=self.tolerance,
                                 msg=f"Magnitude mismatch at frequency bin {i}")
    
    def test_complex_exponential(self):
        """Test DFT of a complex exponential signal."""
        # e^(2πi*n/4) for n = 0, 1, 2, 3
        signal = [cmath.exp(2j * math.pi * n / 4) for n in range(4)]
        
        result = discrete_fourier_transform(signal)
        
        # Should have a single peak at k=1 with magnitude 4
        expected = [0+0j, 4+0j, 0+0j, 0+0j]
        self.assertComplexListAlmostEqual(result, expected)
    
    def test_linearity_property(self):
        """Test that DFT is linear: DFT(a*x + b*y) = a*DFT(x) + b*DFT(y)."""
        x = [1, 2, 3, 4]
        y = [0, 1, 0, 1]
        a, b = 2, 3
        
        # Compute DFT of linear combination
        linear_combination = [a*x[i] + b*y[i] for i in range(len(x))]
        dft_combination = discrete_fourier_transform(linear_combination)
        
        # Compute linear combination of DFTs
        dft_x = discrete_fourier_transform(x)
        dft_y = discrete_fourier_transform(y)
        expected = [a*dft_x[i] + b*dft_y[i] for i in range(len(dft_x))]
        
        self.assertComplexListAlmostEqual(dft_combination, expected)
    
    def test_parseval_theorem(self):
        """Test Parseval's theorem: energy in time domain equals energy in frequency domain."""
        signal = [1, 2, 3, 4]
        dft_result = discrete_fourier_transform(signal)
        
        # Energy in time domain
        time_energy = sum(abs(x)**2 for x in signal)
        
        # Energy in frequency domain (scaled by 1/N)
        freq_energy = sum(abs(X)**2 for X in dft_result) / len(signal)
        
        self.assertAlmostEqual(time_energy, freq_energy, delta=self.tolerance)
    
    def test_dft_idft_roundtrip(self):
        """Test that IDFT(DFT(x)) = x for various signals."""
        test_signals = [
            [1, 0, 0, 0],
            [1, 1, 1, 1],
            [1, 2, 3, 4],
            [1+1j, 2-1j, 3+2j, 4-2j],
            [0.5, -0.3, 1.7, -2.1]
        ]
        
        for signal in test_signals:
            with self.subTest(signal=signal):
                dft_result = discrete_fourier_transform(signal)
                reconstructed = inverse_discrete_fourier_transform(dft_result)
                
                # Convert original signal to complex for comparison
                original_complex = [complex(x) for x in signal]
                self.assertComplexListAlmostEqual(reconstructed, original_complex)
    
    def test_input_validation(self):
        """Test input validation and error handling."""
        # Test empty list
        with self.assertRaises(ValueError):
            discrete_fourier_transform([])
        
        # Test non-list input
        with self.assertRaises(TypeError):
            discrete_fourier_transform("not a list")
        
        # Test invalid element types
        with self.assertRaises(TypeError):
            discrete_fourier_transform([1, 2, "invalid", 4])
    
    def test_different_input_types(self):
        """Test that the function handles different numeric types correctly."""
        # Mix of int, float, and complex
        signal = [1, 2.5, 3+1j, 4.0]
        
        # Should not raise an exception
        result = discrete_fourier_transform(signal)
        self.assertEqual(len(result), 4)
        
        # All results should be complex numbers
        for x in result:
            self.assertIsInstance(x, complex)
    
    def test_magnitude_and_phase_spectrum(self):
        """Test magnitude and phase spectrum extraction functions."""
        signal = [1, 0, -1, 0]
        dft_result = discrete_fourier_transform(signal)
        
        magnitudes = get_magnitude_spectrum(dft_result)
        phases = get_phase_spectrum(dft_result)
        
        # Check that we get the right number of values
        self.assertEqual(len(magnitudes), len(dft_result))
        self.assertEqual(len(phases), len(dft_result))
        
        # Check that magnitudes are non-negative
        for mag in magnitudes:
            self.assertGreaterEqual(mag, 0)
        
        # Check that phases are in the correct range
        for phase in phases:
            self.assertGreaterEqual(phase, -math.pi)
            self.assertLessEqual(phase, math.pi)
    
    def test_symmetry_property_real_input(self):
        """Test that DFT of real signal has conjugate symmetry."""
        # Real signal
        signal = [1, 2, 3, 2]
        dft_result = discrete_fourier_transform(signal)
        
        N = len(dft_result)
        
        # For real input, X[k] = X*[N-k] for k = 1, 2, ..., N-1
        for k in range(1, N):
            conjugate_index = N - k
            expected_conjugate = dft_result[conjugate_index].conjugate()
            self.assertAlmostEqual(dft_result[k].real, expected_conjugate.real, delta=self.tolerance)
            self.assertAlmostEqual(dft_result[k].imag, expected_conjugate.imag, delta=self.tolerance)


def run_performance_test():
    """Simple performance test to ensure the implementation is reasonable."""
    import time
    
    print("\nPerformance Test")
    print("-" * 20)
    
    # Test with different signal lengths
    for N in [8, 16, 32, 64]:
        signal = [math.sin(2 * math.pi * i / N) for i in range(N)]
        
        start_time = time.time()
        result = discrete_fourier_transform(signal)
        end_time = time.time()
        
        print(f"N={N:2d}: {(end_time - start_time)*1000:.2f} ms")


if __name__ == "__main__":
    # Run the unit tests
    print("Running DFT Unit Tests")
    print("=" * 40)
    unittest.main(verbosity=2, exit=False)
    
    # Run performance test
    run_performance_test()
    
    print("\nAll tests completed!")
