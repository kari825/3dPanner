/*
  ==============================================================================

    HRTFProcessor.cpp
    Created: 17 May 2025 4:08:04am
    Author:  admin

  ==============================================================================
*/

#include <algorithm>
#include <array>
#include "HRTFProcessor.h"
#include "BasicSOFA.hpp"

HRTFPanProcessor::HRTFPanProcessor()
    :fft(fftOrder)
{
    //リサイズや初期化
    hrirL_Freq.resize(fftSize);
    hrirR_Freq.resize(fftSize);
    hrirL_Freq_new.resize(fftSize);
    hrirR_Freq_new.resize(fftSize);

    fftInputBufferL.resize(fftSize, 0.0f);
    fftInputBufferR.resize(fftSize, 0.0f);
    outputRemainderBufferL.resize(fftSize, 0.0f);
    outputRemainderBufferR.resize(fftSize, 0.0f);

    fftInputBufferWritePos = 0;
    outputRemainderReadPos = 0;

    juce::dsp::ProcessSpec spec;
    spec.sampleRate = currentSampleRate;         
    spec.maximumBlockSize = 512;       
    spec.numChannels = 1;              

    convolverL.prepare(spec);
    convolverR.prepare(spec);
}

//SOFAファイルからHRIRのデータを設定するためのもの
void HRTFPanProcessor::setHRIRFromSOFA(const BasicSOFA::BasicSOFA& sofa)
{
    currentSOFA = &sofa;

    //方位角、仰角、距離(外部スクリプトではそれぞれtheta,phi,radius)から左右のHRIRデータ取得
    const double theta = currentAzimuth;
    const double phi = currentElevation;
    const double radius = currentRadius;

    const double* hrirL = sofa.getHRIR(0, theta,phi, radius);
    const double* hrirR = sofa.getHRIR(1, theta, phi, radius);

    if (hrirL == nullptr || hrirR == nullptr)
    {
        juce::Logger::writeToLog("HRIR data is nullptr.");
        return;
    }

    sofa.getN();
    HRIRsize = sofa.getN();

    saveOutL.resize(HRIRsize - 1);
    saveOutL.clear();
    saveOutR.resize(HRIRsize - 1);
    saveOutR.clear();

    std::vector<float> timeL(fftSize, 0.0f);
    std::vector<float> timeR(fftSize, 0.0f);

    size_t hrirLength = static_cast<size_t>(sofa.getN());
    size_t numSamples = std::min(hrirLength, static_cast<size_t>(fftSize));

    for (size_t i = 0; i < numSamples; ++i)
    {
        timeL[i] = static_cast<float>(hrirL[i]);
        timeR[i] = static_cast<float>(hrirR[i]);
    }

    hrirL_Freq_new.resize(fftSize);
    hrirR_Freq_new.resize(fftSize);

    //フーリエ変換をして周波数領域のHRIRを保持
    performFFT(hrirL_Freq_new, timeL.data());
    performFFT(hrirR_Freq_new, timeR.data());

    hrirL_Freq = hrirL_Freq_new;
    hrirR_Freq = hrirR_Freq_new;


    isCrossfading = true;
    crossfadeAlpha = 1.0f;
}

void HRTFPanProcessor::setSampleRate(double newSampleRate)
{
    currentSampleRate = newSampleRate;
}

//オーディオブロックごとの処理
void HRTFPanProcessor::processBlock(juce::AudioBuffer<float>& buffer)
{
    const int numChannels = buffer.getNumChannels();
    int numSamples = buffer.getNumSamples();
    if (buffer.getNumChannels() < 2)
    {
        return;
    }

    float* const outL = buffer.getWritePointer(0);
    float* const outR = buffer.getWritePointer(1);
    const float* inL = buffer.getReadPointer(0);
    const float* inR = buffer.getReadPointer(1);
    
    isCrossfading = false;
    fftInputBufferWritePos = 0;

    for (int i = 0; i < numSamples; ++i)
    {
        fftInputBufferL[fftInputBufferWritePos] = inL[i];
        fftInputBufferR[fftInputBufferWritePos] = inR[i];
        fftInputBufferWritePos++;

        // FFTSizeのかたまりごとにFFTをかける
        if (fftInputBufferWritePos == fftSize)
        {
            convolve(buffer, fftSize, i);
        }
        // 最後のサンプルだった場合もたたみ込む
        else if (i == numSamples - 1)
        {
            convolve(buffer, fftInputBufferWritePos, i);
        }
    }

    return;

}

void HRTFPanProcessor::convolve(juce::AudioBuffer<float>& buffer,int sampleNum,int i)
{
    // インデックス番号のリセット
    fftInputBufferWritePos = 0;

    // FFTをかける
    std::vector<std::complex<float>> freqInputL(fftSize);
    std::vector<std::complex<float>> freqInputR(fftSize);
    performFFT(freqInputL, fftInputBufferL.data());
    performFFT(freqInputR, fftInputBufferR.data());

    // HRTFをたたみ込む
    std::vector<std::complex<float>> freqL_convolved(fftSize);
    std::vector<std::complex<float>> freqR_convolved(fftSize);
    for (size_t j = 0; j < fftSize; ++j)
    {
        freqL_convolved[j] = freqInputL[j] * hrirL_Freq[j];
        freqR_convolved[j] = freqInputR[j] * hrirR_Freq[j];
    }

    // IFFTをかける
    std::vector<float> tempOutL(fftSize);
    std::vector<float> tempOutR(fftSize);
    performIFFT(tempOutL.data(), freqL_convolved);
    performIFFT(tempOutR.data(), freqR_convolved);

    // たたみ込んだ結果を出力に代入する
    // overlapAdd法に基づく
    float* out0 = buffer.getWritePointer(0);
    float* out1 = buffer.getWritePointer(1);
    for (int k = 0; k < sampleNum; ++k)
    {
        out0[i - sampleNum + 1 + k] = tempOutL[k] * fftSize;
        out1[i - sampleNum + 1 + k] = tempOutR[k] * fftSize;
    }
    for (int k = 0; k < HRIRsize - 1;k++)
    {
        out0[i - sampleNum + 1+k] += saveOutL[k];
        out1[i - sampleNum + 1 + k] += saveOutR[k];
    }


    saveOutL.clear();
    saveOutR.clear();
    for (int k = 0; k < HRIRsize - 1; ++k)
    {
        saveOutL[k] = tempOutL[sampleNum + k] * fftSize;
        saveOutR[k] = tempOutR[sampleNum + k] * fftSize;
    }

    // FFT用に読み込んだサンプルをリセットする
    fftInputBufferL.clear();
    fftInputBufferR.clear();
}

//JUCEの機能のFFTを用いて時間の信号を周波数の信号に
void HRTFPanProcessor::performFFT(std::vector<std::complex<float>>& out, const float* input)
{
    std::array<float, fftSize * 2> temp;
    std::fill(temp.begin(), temp.end(), 0.0f);

    for (int i = 0;i < fftSize;++i)
        temp[i] = input[i];

    
    fft.performRealOnlyForwardTransform(temp.data());

    for (int i = 0; i < fftSize; ++i)
    {
        out[i] = std::complex<float>(temp[2 * i], temp[2 * i + 1]) ;
    }
}

//JUCEの機能を使った逆高速フーリエ変換
void HRTFPanProcessor::performIFFT(float* out, const std::vector<std::complex<float>>& freqData)
{

    std::vector<float> temp(fftSize * 2, 0.0f);

    for (int i = 0; i < fftSize; ++i)
    {
        temp[2 * i] = freqData[i].real();
        temp[2 * i + 1] = freqData[i].imag();
    }

    
    fft.performRealOnlyInverseTransform(temp.data());


    for (int i = 0; i < fftSize; ++i)
        out[i] = temp[i] / fftSize;

    float maxVal = *std::max_element(out, out + fftSize, [](float a, float b) {
        return std::abs(a) < std::abs(b);
        });
}



void HRTFPanProcessor::setHRTFPosition(float azimuthDeg, float elevationDeg, float radius)
{   
    // 音源の位置を設定。方位角、仰角、距離を更新
    currentAzimuth = azimuthDeg;
    currentElevation = elevationDeg;
    currentRadius = radius;

    if (currentSOFA != nullptr)
    {
        setHRIRFromSOFA(*currentSOFA);
        
        isCrossfading = true;
        crossfadeAlpha = 1.0f;
    }
}
