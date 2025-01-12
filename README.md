# BDTIM-MLk SeisComP Plugin

A SeisComP plugin for calculating local magnitude (MLk) using KOERI's formulas.

## Description

This plugin implements KOERI's local magnitude calculation formulas for SeisComP, providing two different calculation methods based on epicentral distance. The plugin includes both amplitude measurement and magnitude calculation components.

### MLk Formulas

For distances ≤ 200 km:
```
MLk = Log₁₀A - (-1.118-0.0647Δ+0.00071Δ²-3.39x10⁻⁶Δ³+5.71x10⁻⁹Δ⁴)
```

For distances > 200 km:
```
MLk = Log₁₀A + 0.0082Δ - 5.9628x10⁻⁶Δ² + 2.1173
```

where:
- Δ is the epicentral distance in kilometers
- A is the zero-to-peak Wood-Anderson amplitude in millimeters

## Features

- Implements KOERI's local magnitude (MLk) formulas
- Provides distance-dependent magnitude calculation methods
- Converts peak-to-peak amplitudes to zero-to-peak values
- Configurable SNR threshold and maximum distance
- Optional waveform filtering
- Detailed debug logging for amplitude and magnitude calculations

## Prerequisites

- SeisComP 5.0 or later
- C++ compiler supporting C++11 or later
- CMake 3.0 or later

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/comoglu/bdtim-mlk.git
   cd bdtim-mlk
   ```

2. Create a build directory:
   ```bash
   mkdir build
   cd build
   ```

3. Configure and build:
   ```bash
   cmake ..
   make
   ```

4. Install:
   ```bash
   sudo make install
   ```

## Configuration

The plugin can be configured through SeisComP's configuration system. Available parameters:

### Magnitude Configuration

```ini
# In global.cfg or similar
magnitudes.MLk.minSNR = 2.0  # Minimum signal-to-noise ratio
```

### Amplitude Configuration

```ini
# In global.cfg or similar
amplitudes.MLk.maxDist = 8.0  # Maximum distance in degrees
amplitudes.MLk.filter = ""    # Optional waveform filter
```

## Usage

Once installed and configured, the plugin will automatically register itself with SeisComP. The MLk magnitude type will be available for:

- Real-time processing
- Manual analysis in scolv
- Amplitude reviews
- Magnitude calculations

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Author

Mustafa COMOGLU (comoglu@gmail.com)

## License

This project is licensed under the GPL-3.0 License - see the LICENSE file for details.

## Acknowledgments

- KOERI/BDTIM for the magnitude formulas
- SeisComP development team

## References

[To be included]