#pragma once

#include <JuceHeader.h>
#include <BasicSOFA.hpp>
#include "HRTFProcessor.h"
#include "AudioProcessor.h"


//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent : public juce::AudioAppComponent,
    public juce::Button::Listener,
    public juce::Slider::Listener
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent() override;

    //==============================================================================
    void prepareToPlay(int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock(const juce::AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint(juce::Graphics& g) override;
    void resized() override;

    void buttonClicked(juce::Button* button) override;
    void mouseDown(const juce::MouseEvent& e) override;
    void mouseDrag(const juce::MouseEvent& e) override;

    void setElevation(float newElevation);
    void onSliderChanged();
    void sliderValueChanged(juce::Slider* slider) override;

    std::function<void(float azimuth, float elevation, float radiusNorm)> onPositionChanged;
private:

    //==============================================================================
    // Your private member variables go here...

    enum AudioState
    {
        Stopped,
        Starting,
        Playing,
        Stopping
    };
    void changeState(AudioState newState);
    void clearAudioData();

    juce::Point<float> padCenter;
    float padRadius = 150.0f;
    float currentAzimuth = 0.0f;
    float currentElevation = 0.0f;
    float currentRadiusNorm = 0.0f;
    double currentRadius = 0.0;
    float MaxRadius = 150.0f;
    void updateFromPad(const juce::Point<float>& pos);
    bool SOFAreaded = false;
    juce::Slider elevationSlider;

    BasicSOFA::BasicSOFA sofa;

    HRTFPanProcessor hrtfPanProcessor;
    juce::AudioBuffer<float> processingBuffer;

    juce::TextButton openFileButton;
    juce::TextButton playButton;
    juce::TextButton stopButton;

    juce::TextButton SOFAopenButton;

    std::unique_ptr<juce::FileChooser> chooser;
    std::unique_ptr <juce::FileChooser> SOFAchooser;

    juce::AudioBuffer<float> audioBuffer;
    AppAudioProcessor AAProcessor;
    juce::AudioBuffer<float> processedBuffer;
    int playHead = 0;

    std::vector<double>                 thetaList;
    std::vector<double>                 phiList;
    std::vector<double>                 radiusList;

    std::vector<float> hrirL_Freq;
    std::vector<float> hrirR_Freq;

    std::unique_ptr<juce::AudioFormatReader> audioReader;
    juce::File selectedFile;
    juce::AudioFormatManager formatManager;
    std::unique_ptr<juce::AudioFormatReaderSource> readerSource;
    AudioState state;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MainComponent)
};