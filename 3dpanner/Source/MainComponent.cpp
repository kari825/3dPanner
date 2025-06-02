#include <JuceHeader.h>
#include "MainComponent.h"
#include <BasicSOFA.hpp>


//==============================================================================
MainComponent::MainComponent():state(Stopped)
{
    // UIのボタンやスライダーの作成とリスナーを設定
    addAndMakeVisible(&openFileButton);
    openFileButton.setButtonText("Open WAV File");
    openFileButton.addListener(this);

    addAndMakeVisible(&playButton);
    playButton.setButtonText("Play");
    playButton.addListener(this);
    playButton.setColour(juce::TextButton::buttonColourId, juce::Colours::green);
    playButton.setEnabled(false);

    addAndMakeVisible(&stopButton);
    stopButton.setButtonText("Stop");
    stopButton.addListener(this);
    stopButton.setColour(juce::TextButton::buttonColourId, juce::Colours::red);
    stopButton.setEnabled(false);


    SOFAopenButton.setButtonText("Open SOFA File");
    SOFAopenButton.addListener(this);
    addAndMakeVisible(&SOFAopenButton);

    elevationSlider.setSliderStyle(juce::Slider::LinearVertical);
    elevationSlider.setRange(-90.0, 90.0, 1.0);
    elevationSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 60, 20);
    elevationSlider.setValue(0.0);
    addAndMakeVisible(elevationSlider);
    elevationSlider.addListener(this);

    formatManager.registerBasicFormats();
    // Make sure you set the size of the component after
    // you add any child components.
    setSize (800, 400);

    // Some platforms require permissions to open input channels so request that here
    if (juce::RuntimePermissions::isRequired (juce::RuntimePermissions::recordAudio)
        && ! juce::RuntimePermissions::isGranted (juce::RuntimePermissions::recordAudio))
    {
        juce::RuntimePermissions::request (juce::RuntimePermissions::recordAudio,
                                           [&] (bool granted) { setAudioChannels (granted ? 2 : 0, 2); });
    }
    else
    {
        // Specify the number of input and output channels that we want to open
        setAudioChannels (0, 2);
    }
    
    // UIを変更した際にsetHRTFPositionを呼び出す
    onPositionChanged = [this](float azimuth, float elevation, float currentRadiusNorm)
    {
        hrtfPanProcessor.setHRTFPosition(360 - azimuth, elevation, currentRadiusNorm);        
    };
}

MainComponent::~MainComponent()
{
    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    // This function will be called when the audio device is started, or when
    // its settings (i.e. sample rate, block size, etc) are changed.

    // You can use this function to initialise any resources you might need,
    // but be careful - it will be called on the audio thread, not the GUI thread.

    // For more details, see the help for AudioProcessor::prepareToPlay()
    hrtfPanProcessor.setSampleRate(sampleRate);
}

void MainComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    // オーディオのメイン処理
    bufferToFill.clearActiveBufferRegion();
    
    if (state == Playing && audioBuffer.getNumSamples() > 0)
    {
        const int numSampleToProcess = bufferToFill.numSamples;
        int remainingSamplesInFile = audioBuffer.getNumSamples() - playHead;
        if (remainingSamplesInFile <= 0)
        {
            changeState(Stopped);
            return;
        }
        const int samplesThisBlock = std::min(numSampleToProcess, remainingSamplesInFile);
        juce::AudioBuffer<float> inputBlock(audioBuffer.getNumChannels(), samplesThisBlock);
        for (int channel = 0; channel < audioBuffer.getNumChannels();++channel)
        {
            inputBlock.copyFrom(channel, 0, audioBuffer.getReadPointer(channel) + playHead, samplesThisBlock);
        }
        // 受け取ったオーディオデータがモノラルならモノラルからステレオに
        AAProcessor.setInputBuffer(inputBlock);
        AAProcessor.process();
        processedBuffer = AAProcessor.getOutputBuffer();

        if (processedBuffer.getNumSamples() == 0 || processedBuffer.getNumChannels() == 0)
        {
            return;
        }
    
        juce::AudioBuffer<float> hrtfInputBuffer;
        hrtfInputBuffer.makeCopyOf(AAProcessor.getOutputBuffer());
        // HRTFの畳み込み処理を行う
        hrtfPanProcessor.processBlock(hrtfInputBuffer);
        
        // 処理が終わったものを出力バッファにコピー
        for (int channel = 0;channel < bufferToFill.buffer->getNumChannels(); ++channel)
        {
            bufferToFill.buffer->copyFrom(channel, 0, hrtfInputBuffer, channel, 0, samplesThisBlock);
        }
        // 再生位置の更新
        playHead += samplesThisBlock;
    }
    else
    {
        bufferToFill.clearActiveBufferRegion();
    }

}

void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
}

//==============================================================================
void MainComponent::paint (juce::Graphics& g)
{
    // UIのメイン処理
    // 円形パッドや動かす黄色い円の描画
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    g.setFont(juce::Font(16.0f));
    
    g.setColour(juce::Colours::lightblue);
    g.fillEllipse(padCenter.x - padRadius, padCenter.y - padRadius, padRadius * 2, padRadius * 2);
    g.setColour(juce::Colours::white);
    float lineThickness = 3.0f;
    g.drawLine(padCenter.x, 19 * getHeight() / 20, padCenter.x, getHeight() / 20,lineThickness);
    g.drawLine(padCenter.x-9*getHeight()/20, padCenter.y, padCenter.x+9*getHeight() / 20, padCenter.y,lineThickness);
    float angleRad = juce::degreesToRadians(currentAzimuth);
    juce::Point<float> point = padCenter + juce::Point<float>(std::sin(angleRad), -std::cos(angleRad)) * (currentRadiusNorm * padRadius);
    float r = currentRadiusNorm * padRadius;
    juce::Point<float> offset(std::sin(angleRad), -std::cos(angleRad));

    float ballsize = 20;

    g.setColour(juce::Colours::white);
    g.fillEllipse(padCenter.x-ballsize/2-2, padCenter.y-ballsize/2-2, ballsize+4, ballsize+4);

    g.setColour(juce::Colours::yellow);
    g.fillEllipse(point.x-ballsize/2 , point.y-ballsize/2 , ballsize, ballsize);

}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.


    openFileButton.setBounds(4*getWidth() / 5, getHeight() / 6, 100, 40);
    playButton.setBounds(4 * getWidth() / 5, 2 * getHeight() / 6, 100, 40);
    stopButton.setBounds(4 * getWidth() / 5, 3 * getHeight() / 6, 100, 40);
    SOFAopenButton.setBounds(4 * getWidth() / 5, 4 * getHeight() / 6, 100, 40);
    elevationSlider.setBounds(3*getWidth() / 5, 10,30, getHeight() - 20);
    padCenter = juce::Point<float>(getWidth() / 3.0f, getHeight() / 2.0f);
    MaxRadius = std::min(getWidth() / 3.0f, getHeight() / 2.0f) - 10.0f;
}

// UIでボタンを押したときの処理
void MainComponent::buttonClicked(juce::Button* button)
{
    //wavファイルを開いたとき
    if (button == &openFileButton)
    {
        clearAudioData();
        // ファイルチューザーでwavファイルを選ばせる
        chooser = std::make_unique<juce::FileChooser>("Select a Wave file to play...",
            juce::File{},
            "*.wav");
        auto chooserFlags = juce::FileBrowserComponent::openMode
            | juce::FileBrowserComponent::canSelectFiles;

        chooser->launchAsync(chooserFlags, [this](const juce::FileChooser& fc)
        {
            auto file = fc.getResult();

            if (file != juce::File{})
            {
                juce::String filePath = file.getFullPathName();

                if (!file.exists()) {
                    juce::AlertWindow::showMessageBoxAsync(juce::AlertWindow::WarningIcon, "Error", "File does not exist.");
                    return;
                }

                std::unique_ptr<juce::AudioFormatReader> reader(formatManager.createReaderFor(file));
                if (reader != nullptr)
                {
                    // 選択したファイルを読み込み、audioBufferに格納
                    audioBuffer.setSize((int)reader->numChannels, (int)reader->lengthInSamples);
                    reader->read(&audioBuffer, 0, (int)reader->lengthInSamples, 0, true, true);
                    AAProcessor.setInputBuffer(audioBuffer);
                    AAProcessor.process();
                    processedBuffer = AAProcessor.getOutputBuffer();
                    playButton.setEnabled(true);                
                }
                else
                {
                    juce::AlertWindow::showMessageBoxAsync(juce::AlertWindow::WarningIcon, "Error", "Could not read the file.");
                }
            }
        });
    }
    else if (button == &playButton)
    {
        // playボタンを押した時
        playHead = 0;
        changeState(Starting);
    }
    else if (button == &stopButton)
    {
        // stopボタンを押した時
        changeState(Stopping);
    }
    else if (button == &SOFAopenButton)
    {
        // SOFAファイルをファイルチューザーで選ばせる
        SOFAchooser = std::make_unique<juce::FileChooser>("Select a SOFA file...",
            juce::File{},
            "*.sofa");
        
        SOFAchooser->launchAsync(
        juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles,
        [this](const juce::FileChooser& fc)
        {
            auto file = fc.getResult();
       
            if (file.existsAsFile())
            {
                std::string path = file.getFullPathName().toStdString();
                if (sofa.readSOFAFile(path))
                {
                    // SOFAファイルのHRIRデータを設定
                    hrtfPanProcessor.setHRIRFromSOFA(sofa);
                }
            }
        });
    }
}

// オーディオの再生状態の変更とUIの有効/無効を切り替え
void MainComponent::changeState(AudioState newState)
{
    if (state != newState)
    {
        state = newState;

        switch (state)
        {
        case Stopped:
            stopButton.setEnabled(false);
            playButton.setEnabled(true);
            break;

        case Starting:
            playButton.setEnabled(false);
            stopButton.setEnabled(true);
            if (audioBuffer.getNumSamples()> 0) { 
                changeState(Playing);
            }
            else {
                juce::AlertWindow::showMessageBoxAsync(juce::AlertWindow::WarningIcon, "Error", "No file loaded.");
                changeState(Stopped);
            }
            break;

        case Playing:
            stopButton.setEnabled(true);
            break;

        case Stopping:
            stopButton.setEnabled(false);
            playButton.setEnabled(true);
            break;
        }
    }
}

// マウスクリック、マウスドラッグを押したときに再描画
void MainComponent::mouseDown(const juce::MouseEvent& event)
{
    updateFromPad(event.position);
    repaint();
}

void MainComponent::mouseDrag(const juce::MouseEvent& e)
{
    updateFromPad(e.position);
    repaint();
}

// マウスカーソルの位置から音源の方位角と距離を計算
void MainComponent::updateFromPad(const juce::Point<float>& pos)
{
    auto delta = pos - padCenter;
    float angle = std::atan2(delta.x, -delta.y);
    float newAzimuth= juce::radiansToDegrees(angle);
    if (newAzimuth < 0.0f)
        newAzimuth += 360.0f;
    float distance = delta.getDistanceFromOrigin();
   
    currentRadiusNorm = juce::jlimit(0.0f, 1.0f, distance / MaxRadius);

    if (std::abs(newAzimuth - currentAzimuth) > 0.001f)
    {
        currentAzimuth = newAzimuth;
        if (onPositionChanged)
                onPositionChanged(currentAzimuth, currentElevation, currentRadiusNorm);
    }   
}

// 仰角のスライダーをもとに仰角を更新
void MainComponent::setElevation(float newElevation)
{
    if (std::abs(newElevation - currentElevation) > 0.001f) 
    {
        currentElevation = juce::jlimit(-90.0f,90.0f,newElevation);

        if (onPositionChanged)
            onPositionChanged(currentAzimuth, currentElevation, currentRadiusNorm);
    }
    repaint();
}


void MainComponent::onSliderChanged()
{
    float newElevation = elevationSlider.getValue();
    setElevation(newElevation);
}

void MainComponent::sliderValueChanged(juce::Slider* slider)
{
    if (slider == &elevationSlider)
    {
        onSliderChanged();
    }
}


// wavファイルを再び読み込んだ時にすでにあるデータをクリア
void MainComponent::clearAudioData()
{
    readerSource.reset();
    audioReader.reset();
    selectedFile = juce::File{};
    playButton.setEnabled(false);
    stopButton.setEnabled(false);
    state = Stopped;
}