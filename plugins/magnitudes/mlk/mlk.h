/*
 * File:   mlk.h
 * KOERI Local Magnitude Implementation
 */

#ifndef __MLK_PLUGIN_H__
#define __MLK_PLUGIN_H__

#define KOERI_ML_AMP_TYPE "MLk"
#define KOERI_ML_MAG_TYPE "MLk"

#include <seiscomp/processing/amplitudes/MLv.h>
#include <seiscomp/processing/magnitudeprocessor.h>
#include <seiscomp/core/plugin.h>
#include <seiscomp/geo/featureset.h>

#include <string>

/*
Calculates the MLk amplitude. This amplitude value is used by the MLk magnitude
processor. The amplitude value calculated is similar to MLv but returns
zero-to-peak value instead of peak-to-peak.
*/
class Amplitude_MLk : public Seiscomp::Processing::AmplitudeProcessor_MLv {
    DECLARE_SC_CLASS(Amplitude_MLk);
public:
    explicit Amplitude_MLk(const std::string& type=KOERI_ML_AMP_TYPE);
    virtual int capabilities() const;
    bool setup(const Seiscomp::Processing::Settings &settings);
    virtual IDList capabilityParameters(Capability cap) const;
    virtual bool setParameter(Capability cap, const std::string &value);

protected:
    virtual std::string defaultFilter() const { return ""; };
    bool computeAmplitude(const Seiscomp::DoubleArray &data,
            size_t i1, size_t i2,
            size_t si1, size_t si2,
            double offset,
            AmplitudeIndex *dt, AmplitudeValue *amplitude,
            double *period, double *snr);
};

/*
Calculates the MLk magnitude using KOERI's formulas.
For distances ≤ 200 km:
MLk = Log10A - (-1.118-0.0647Δ+0.00071Δ²-3.39x10⁻⁶Δ³+5.71x10⁻⁹Δ⁴)
For distances > 200 km:
MLk = Log10A + 0.0082Δ - 5.9628x10⁻⁶Δ² + 2.1173
*/
class Magnitude_MLk : public Seiscomp::Processing::MagnitudeProcessor {
    public:
        explicit Magnitude_MLk(const std::string& type = KOERI_ML_MAG_TYPE);
        bool setup(const Seiscomp::Processing::Settings &settings);
        std::string amplitudeType() const;
        Seiscomp::Processing::MagnitudeProcessor::Status computeMagnitude(
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
              double &value);

        bool treatAsValidMagnitude() const override;

    private:
        static double distance(double delta, double depth);
        bool _validValue;
        double _minSNR;

        // Helper functions for different distance ranges
        Seiscomp::Processing::MagnitudeProcessor::Status computeMagShortDistance(
              double amplitude,
              double delta,
              double depth,
              double &value);

        Seiscomp::Processing::MagnitudeProcessor::Status computeMagLongDistance(
              double amplitude,
              double delta,
              double depth,
              double &value);
};

#endif /* __MLK_PLUGIN_H__ */