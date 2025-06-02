/*
  ==============================================================================

    AudioProcessor.cpp
    Created: 14 May 2025 2:26:31pm

  ==============================================================================
*/

#include "AudioProcessor.h"
#include <JuceHeader.h>
#include <BasicSOFA.hpp>

//もともと窓関数を掛けるなどの処理をしていたスクリプトだったがその処理が必要なくなったため最低限の処理のスクリプトになってしまっている
AppAudioProcessor::AppAudioProcessor()
{

}

void AppAudioProcessor::setInputBuffer(const juce::AudioBuffer<float>& buffer)
{
    inputBuffer.makeCopyOf(buffer);
}



void AppAudioProcessor::process()
{
    if (inputBuffer.getNumSamples() == 0)
        return;

    juce::AudioBuffer<float> AAPBuffer;
    AAPBuffer.makeCopyOf(inputBuffer);

    //コピーされたバッファのチャンネル数が1(モノラル)だった場合ステレオにする
    if (AAPBuffer.getNumChannels() == 1)
    {
        AAPBuffer.setSize(2, AAPBuffer.getNumSamples(), true, true, true);
        AAPBuffer.copyFrom(1, 0, AAPBuffer, 0, 0, AAPBuffer.getNumSamples());
    }
    outputBuffer.makeCopyOf(AAPBuffer);
}

const juce::AudioBuffer<float>& AppAudioProcessor::getOutputBuffer() const
{
    return outputBuffer;
}
