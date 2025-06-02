/*
  ==============================================================================

    HRTFProcessor.h
    Created: 17 May 2025 4:08:04am
    Author:  admin

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include <BasicSOFA.hpp>
#include <vector>


class HRTFPanProcessor
{
public:
    HRTFPanProcessor();
    void processBlock(juce::AudioBuffer<float>& buffer);
    void setHRIRFromSOFA(const BasicSOFA::BasicSOFA& sofa);
    
    void convolve(juce::AudioBuffer<float>& buffer,int sampleNum,int i);

    void setHRTFPosition(float azimurhDeg, float elevationDeg, float radius);
    void setSampleRate(double newSampleRate);
private:
    float currentAzimuth = 0.0f;
    float currentElevation = 0.0f;
    float currentRadius = 0.0f;
    
    

    juce::dsp::Convolution convolverL;
    juce::dsp::Convolution convolverR;

    static constexpr int fftOrder = 11;
    static constexpr int fftSize = 1 << fftOrder;
    juce::dsp::FFT fft;

    std::vector<std::complex<float>> hrirL_Freq;
    std::vector<std::complex<float>> hrirR_Freq;

    std::vector<std::complex<float>> hrirL_Freq_new;
    std::vector<std::complex<float>> hrirR_Freq_new;
    float crossfadeAlpha = 1.0f;
    bool isCrossfading = false;

    int HRIRsize = 0;

    std::vector<float> saveOutL;
    std::vector<float> saveOutR;

    std::vector<float> windowBuffer;

    std::vector<float> fftInputBufferL;
    std::vector<float> fftInputBufferR;

    std::vector<float> outputRemainderBufferL;
    std::vector<float> outputRemainderBufferR;

    int fftInputBufferWritePos = 0;
    int outputRemainderReadPos = 0;

    static constexpr int overlapSize = fftSize / 2;
    static constexpr int hopSize = fftSize - overlapSize;

    void performFFT(std::vector<std::complex<float>>& out, const float* input);
    void performIFFT(float*out, const std::vector <std::complex<float>>& freqData);

    const BasicSOFA::BasicSOFA* currentSOFA = nullptr;

    double roundedThetaKey(double theta) { return std::round(theta); }
    double roundedPhiKey(double phi) { return std::round(phi); }
    double roundedRadiusKey(double radius) { return std::round(radius * 100.0) / 100.0; }

    double currentSampleRate = 44100.0;
};