
import unittest
from dft import discrete_fourier_transform, inverse_discrete_fourier_transform


class TestDFT(unittest.TestCase):

    def test_core_logic(self):
        tol = 1e-10

        # impulse
        result = discrete_fourier_transform([1, 0, 0, 0])
        expected = [1+0j, 1+0j, 1+0j, 1+0j]
        for i in range(4):
            self.assertAlmostEqual(result[i].real, expected[i].real, delta=tol)
            self.assertAlmostEqual(result[i].imag, expected[i].imag, delta=tol)

        #  constant
        result = discrete_fourier_transform([1, 1, 1, 1])
        self.assertAlmostEqual(result[0].real, 4, delta=tol)
        self.assertAlmostEqual(result[0].imag, 0, delta=tol)
        for i in range(1, 4):
            self.assertAlmostEqual(abs(result[i]), 0, delta=tol)

        # Linearity: DFT(a*x + b*y) = a*DFT(x) + b*DFT(y)
        x, y, a, b = [1, 2], [3, 4], 2, 3
        combo = [a*x[i] + b*y[i] for i in range(2)]
        dft_combo = discrete_fourier_transform(combo)
        dft_x, dft_y = discrete_fourier_transform(x), discrete_fourier_transform(y)
        expected = [a*dft_x[i] + b*dft_y[i] for i in range(2)]
        for i in range(2):
            self.assertAlmostEqual(dft_combo[i].real, expected[i].real, delta=tol)
            self.assertAlmostEqual(dft_combo[i].imag, expected[i].imag, delta=tol)

        # Round-check: IDFT(DFT(x)) = x
        signal = [1, 2, 3, 4]
        dft_result = discrete_fourier_transform(signal)
        reconstructed = inverse_discrete_fourier_transform(dft_result)
        for i in range(4):
            self.assertAlmostEqual(reconstructed[i].real, signal[i], delta=tol)
            self.assertAlmostEqual(reconstructed[i].imag, 0, delta=tol)


if __name__ == "__main__":
    unittest.main()
