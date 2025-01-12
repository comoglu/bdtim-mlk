#define SEISCOMP_COMPONENT MLk

#include "mlk.h"

#include <seiscomp/logging/log.h>
#include <seiscomp/core/strings.h>
#include <seiscomp/math/geo.h>
#include <seiscomp/processing/magnitudeprocessor.h>

#include <math.h>

using namespace std;

namespace {

// Plugin author info
#define AUTHOR_NAME "Mustafa COMOGLU"
#define AUTHOR_EMAIL "comoglu@gmail.com"

}

// Register plugin to SeisComP
ADD_SC_PLUGIN(
    "MLk magnitude implementation by KOERI/BDTIM",
    AUTHOR_NAME " <" AUTHOR_EMAIL ">",
    1, 0, 0 // Version 1.0.0
)

// Register the amplitude processor
IMPLEMENT_SC_CLASS_DERIVED(Amplitude_MLk, AmplitudeProcessor, "Amplitude_MLk");
REGISTER_AMPLITUDEPROCESSOR(Amplitude_MLk, KOERI_ML_AMP_TYPE);

// Register the magnitude processor
REGISTER_MAGNITUDEPROCESSOR(Magnitude_MLk, KOERI_ML_MAG_TYPE);

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// MLk AMPLITUDE PROCESSOR IMPLEMENTATION

Amplitude_MLk::Amplitude_MLk(const std::string& type)
: Seiscomp::Processing::AmplitudeProcessor_MLv() {
    _type = type;
    setUnit("mm");
    setMinSNR(0);
}

bool Amplitude_MLk::setup(const Seiscomp::Processing::Settings &settings) {
    if (!AmplitudeProcessor_MLv::setup(settings)) {
        return false;
    }

    // Configure filter if specified
    std::string filterString;
    try {
        filterString = settings.getString("amplitudes." + _type + ".filter");
        SEISCOMP_DEBUG("MLk: Using configured filter: %s", filterString.c_str());
    }
    catch(...) {
        filterString = defaultFilter();
        SEISCOMP_DEBUG("MLk: Using default filter: %s", filterString.c_str());
    }

    if (!filterString.empty()) {
        _preFilter = filterString;
    }

    // Set maximum distance
    double maxDist;
    try {
        maxDist = settings.getDouble("amplitudes." + _type + ".maxDist");
        SEISCOMP_DEBUG("MLk: Using configured maxDist: %.2f", maxDist);
    }
    catch(...) {
        maxDist = 8.0;  // Default max distance in degrees
        SEISCOMP_DEBUG("MLk: Using default maxDist: %.2f", maxDist);
    }
    setMaxDist(maxDist);

    return true;
}

int Amplitude_MLk::capabilities() const {
    return NoCapability;
}

Seiscomp::Processing::AmplitudeProcessor::IDList
Amplitude_MLk::capabilityParameters(Capability cap) const {
    return IDList();
}

bool Amplitude_MLk::setParameter(Capability cap, const std::string &value) {
    return false;
}

bool Amplitude_MLk::computeAmplitude(
        const Seiscomp::DoubleArray &data,
        size_t i1, size_t i2,
        size_t si1, size_t si2,
        double offset,
        AmplitudeIndex *dt,
        AmplitudeValue *amplitude,
        double *period, double *snr) {

    SEISCOMP_DEBUG("MLk: Computing amplitude...");
    
    bool success = AmplitudeProcessor_MLv::computeAmplitude(
        data, i1, i2, si1, si2, offset,
        dt, amplitude, period, snr
    );

    if (success) {
        // Convert peak-to-peak to zero-to-peak
        double originalValue = amplitude->value;
        amplitude->value *= 0.5;

        SEISCOMP_DEBUG("MLk Amplitude Computation:");
        SEISCOMP_DEBUG("  Converting peak-to-peak to zero-to-peak");
        SEISCOMP_DEBUG("  Original value (peak-to-peak) = %.3f mm", originalValue);
        SEISCOMP_DEBUG("  Converted value (zero-to-peak) = %.3f mm", amplitude->value);
        if (period) SEISCOMP_DEBUG("  Period = %.3f s", *period);
        if (snr) SEISCOMP_DEBUG("  SNR = %.2f", *snr);
    }

    return success;
}

// END MLk AMPLITUDE PROCESSOR
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// MLk MAGNITUDE PROCESSOR IMPLEMENTATION

Magnitude_MLk::Magnitude_MLk(const std::string& type)
: Seiscomp::Processing::MagnitudeProcessor(type)
, _minSNR(2.0) {}

bool Magnitude_MLk::setup(const Seiscomp::Processing::Settings &settings) {
    if (!MagnitudeProcessor::setup(settings))
        return false;

    try {
        _minSNR = settings.getDouble("magnitudes." + _type + ".minSNR");
        SEISCOMP_DEBUG("MLk: Using configured minSNR: %.2f", _minSNR);
    }
    catch (...) {
        SEISCOMP_DEBUG("MLk: Using default minSNR: %.2f", _minSNR);
    }

    return true;
}

std::string Magnitude_MLk::amplitudeType() const {
    return KOERI_ML_AMP_TYPE;
}

double Magnitude_MLk::distance(double delta, double depth) {
    return Seiscomp::Math::Geo::deg2km(delta);
}

Seiscomp::Processing::MagnitudeProcessor::Status
Magnitude_MLk::computeMagShortDistance(double amplitude, double delta, double depth, double &value) {
    double r = distance(delta, depth);
    
    SEISCOMP_DEBUG("MLk Short Distance Formula (Δ ≤ 200km):");
    SEISCOMP_DEBUG("MLk = Log10A - (-1.118-0.0647Δ+0.00071Δ²-3.39x10⁻⁶Δ³+5.71x10⁻⁹Δ⁴)");
    SEISCOMP_DEBUG("where:");
    SEISCOMP_DEBUG("  A = %.4f mm (amplitude)", amplitude);
    SEISCOMP_DEBUG("  Δ = %.2f km (distance)", r);
    
    double correction = -1.118 - 
                       (0.0647 * r) + 
                       (0.00071 * pow(r, 2)) - 
                       (3.39e-6 * pow(r, 3)) + 
                       (5.71e-9 * pow(r, 4));
    
    value = log10(amplitude) - correction;
    
    SEISCOMP_DEBUG("Calculation steps:");
    SEISCOMP_DEBUG("  Log10(A) = %.4f", log10(amplitude));
    SEISCOMP_DEBUG("  Term1 = -1.118");
    SEISCOMP_DEBUG("  Term2 = %.4f (-0.0647Δ)", -0.0647 * r);
    SEISCOMP_DEBUG("  Term3 = %.4f (0.00071Δ²)", 0.00071 * pow(r, 2));
    SEISCOMP_DEBUG("  Term4 = %.4f (-3.39x10⁻⁶Δ³)", -3.39e-6 * pow(r, 3));
    SEISCOMP_DEBUG("  Term5 = %.4f (5.71x10⁻⁹Δ⁴)", 5.71e-9 * pow(r, 4));
    SEISCOMP_DEBUG("  Correction = %.4f", correction);
    SEISCOMP_DEBUG("  Final MLk = %.2f", value);
    
    return OK;
}

Seiscomp::Processing::MagnitudeProcessor::Status
Magnitude_MLk::computeMagLongDistance(double amplitude, double delta, double depth, double &value) {
    double r = distance(delta, depth);
    
    SEISCOMP_DEBUG("MLk Long Distance Formula (Δ > 200km):");
    SEISCOMP_DEBUG("MLk = Log10A + 0.0082Δ - 5.9628x10⁻⁶Δ² + 2.1173");
    SEISCOMP_DEBUG("where:");
    SEISCOMP_DEBUG("  A = %.4f mm (amplitude)", amplitude);
    SEISCOMP_DEBUG("  Δ = %.2f km (distance)", r);
    
    value = log10(amplitude) + 
            (0.0082 * r) - 
            (5.9628e-6 * pow(r, 2)) + 
            2.1173;
    
    SEISCOMP_DEBUG("Calculation steps:");
    SEISCOMP_DEBUG("  Log10(A) = %.4f", log10(amplitude));
    SEISCOMP_DEBUG("  Term1 = %.4f (0.0082Δ)", 0.0082 * r);
    SEISCOMP_DEBUG("  Term2 = %.4f (-5.9628x10⁻⁶Δ²)", -5.9628e-6 * pow(r, 2));
    SEISCOMP_DEBUG("  Term3 = 2.1173 (constant)");
    SEISCOMP_DEBUG("  Final MLk = %.2f", value);
    
    return OK;
}

Seiscomp::Processing::MagnitudeProcessor::Status
Magnitude_MLk::computeMagnitude(
        double amplitudeValue,
        const std::string &unit,
        double period,
        double snr,
        double delta,
        double depth,
        const Seiscomp::DataModel::Origin *hypocenter,
        const Seiscomp::DataModel::SensorLocation *receiver,
        const Seiscomp::DataModel::Amplitude *amplitude,
        const Seiscomp::Processing::MagnitudeProcessor::Locale *locale,
        double &value) {

    _validValue = false;

    if (amplitudeValue <= 0) {
        SEISCOMP_DEBUG("MLk: Invalid amplitude value (%.2f)", amplitudeValue);
        return AmplitudeOutOfRange;
    }

    double r = distance(delta, depth);
    SEISCOMP_DEBUG("MLk: Computing magnitude for distance %.2f km", r);
    SEISCOMP_DEBUG("MLk Input parameters:");
    SEISCOMP_DEBUG("  Amplitude = %.4f mm", amplitudeValue);
    SEISCOMP_DEBUG("  Period = %.2f s", period);
    SEISCOMP_DEBUG("  SNR = %.2f", snr);
    SEISCOMP_DEBUG("  Delta = %.2f degrees", delta);
    SEISCOMP_DEBUG("  Depth = %.2f km", depth);
    
    Status status;
    if (r <= 200) {
        SEISCOMP_DEBUG("MLk: Using short distance formula (Δ ≤ 200km)");
        status = computeMagShortDistance(amplitudeValue, delta, depth, value);
    } else {
        SEISCOMP_DEBUG("MLk: Using long distance formula (Δ > 200km)");
        status = computeMagLongDistance(amplitudeValue, delta, depth, value);
    }

    if (snr < _minSNR) {
        status = SNROutOfRange;
        _validValue = true;
        SEISCOMP_DEBUG("MLk: SNR %.2f below threshold %.2f", snr, _minSNR);
    }

    return status;
}

bool Magnitude_MLk::treatAsValidMagnitude() const {
    return _validValue;
}

// END MLk MAGNITUDE PROCESSOR
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<