/*
  ==============================================================================

    AudioProcessor.h
    Created: 14 May 2025 2:26:31pm
    Author:  admin

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include <BasicSOFA.hpp>


class AppAudioProcessor
{
public:
    AppAudioProcessor();

    void setInputBuffer(const juce::AudioBuffer<float>& buffer);
    void process();
    const juce::AudioBuffer<float>& getOutputBuffer() const;

private:

    std::vector<float> windowTable;

    juce::AudioBuffer<float> inputBuffer;
    juce::AudioBuffer<float> outputBuffer;

   
    
};