function [] = runASPP()
    global AUDIO_DEVICE_INDEX ASPP_DIR_PATH
    global sliderUpdateMode getParamForSliderValueFunctions parametricASPPFunctions stftFunction inverseSTFTFunction windowSize codebook asppSimplexPoint inputUtteranceIDs inputNoiseTypes inputSNRs
        
    % User interface parameters
    sliderUpdateMode = 'continuous'; % 'onRelease'
    asppSimplexPoint = [1,0,0];
    showSaveButon = 0;
    
    % Input example parameters
    inputUtteranceIDs = {'dr2_fajw0_si633','dr2_mwsb0_si2256','dr3_fcmg0_si1242','dr3_makr0_si722','dr4_fcag0_si1503','dr5_mvlo0_si517','dr6_flag0_si2094','dr6_flag0_si834','dr6_meal0_si1547','dr6_meal0_si917','dr7_mbom0_si2274','dr7_mtml0_si1065'};
    inputNoiseTypes = {'Babble','Ventilation'};
    inputSNRs = {'5 dB','0 dB','-5 dB'};
    
    % STFT parameters
    windowSize = 1024;
    hopSize = 256;
    stftFunction = @(samples) getSTFT(samples, windowSize, hopSize);
    inverseSTFTFunction = @(complexSpectrogram) getISTFT(complexSpectrogram, windowSize, hopSize);
    
    % NMF parameters
    dictionarySize = 1024;
    % Load NMF dictionary
    dictionaryFilePath = sprintf( [ASPP_DIR_PATH filesep 'resources' filesep 'PrelearnedCodebooks' filesep 'PrelearnedNMFDictionary_%d.mat'], dictionarySize) ;
    matFileW = load(dictionaryFilePath);
    W = matFileW.W;
    codebook = W ./ sum(W,1); % normalize such that individual atoms to sum to 1
    
    % ASPP parameters
    ildSigmaRange = [0.1, 50.0];
    ipdKappaRange = [0.1, 100.0];
    icmBetaRange = [2, 500];
    % Define ASPP functions for different binaural cues
    getParamForSliderValue = @(sliderValue, minValue, maxValue) minValue * exp( log(maxValue/minValue) * sliderValue );
    getILDParamForSliderValue = @(sliderValue) getParamForSliderValue(sliderValue, ildSigmaRange(1), ildSigmaRange(2));
    getIPDParamForSliderValue = @(sliderValue) getParamForSliderValue(sliderValue, ipdKappaRange(1), ipdKappaRange(2));
    getICMParamForSliderValue = @(sliderValue) getParamForSliderValue(sliderValue, icmBetaRange(1), icmBetaRange(2));
    getParamForSliderValueFunctions = {getILDParamForSliderValue, getIPDParamForSliderValue, getICMParamForSliderValue};
    
    parametricILDASPPFunction = @(atomicILDDistances,sliderValue) getGaussianASPP(atomicILDDistances, getILDParamForSliderValue(sliderValue), 40);
    parametricIPDASPPFunction = @(atomicIPDDistances,sliderValue) getVonMisesASPP(atomicIPDDistances, getIPDParamForSliderValue(sliderValue));
    parametricICMASPPFunction = @(atomicICMDistances,sliderValue) getBetaASPP(atomicICMDistances, getICMParamForSliderValue(sliderValue));
    parametricASPPFunctions = {parametricILDASPPFunction, parametricIPDASPPFunction, parametricICMASPPFunction};
    
    setInputExample(inputUtteranceIDs{1}, 'Babble', '5 dB');
                       
    constructGrahpicalInterface(getParamForSliderValueFunctions, showSaveButon);
end

function [noisyWavFilePath] = getNoisyWavFilePath(utteranceID, noiseType, snr)
    global ASPP_DIR_PATH
    if strcmp(noiseType, 'Babble')
        noiseDir = 'HRIR_TIMIT_bab';
    else
        noiseDir = 'HRIR_TIMIT_fan';
    end
    snrDir = snr(find(~isspace(snr)));
        
    noisyWavFilePath = [ASPP_DIR_PATH filesep 'resources' filesep noiseDir filesep snrDir filesep [utteranceID '_mix.wav']];
end

function [] = setInputExample(utteranceID, noiseType, snr)
    global AUDIO_DEVICE_INDEX cleanSpeechAudioPlayer noisySpeechAudioPlayer createRecSpeechAudioPlayer cleanSpeechSamples noisySpeechSamples sampleRate
    global codebook stftFunction windowSize parametricASPPFunctions codebookMaskSliders
    
    clear global spectrograms codebookCoefficients atomicTargetDistances codebookMasksList 
    global spectrograms codebookCoefficients atomicTargetDistances codebookMasksList
    
    noisyWavFilePath = getNoisyWavFilePath(utteranceID, noiseType, snr);
    cleanWavFilePath = strrep(noisyWavFilePath, '_mix.wav', '_sim.wav');
    
    % Load speech signals
    [cleanSpeechSamples, sampleRate] = audioread(cleanWavFilePath);
    cleanSpeechAudioPlayer = audioplayer(cleanSpeechSamples, sampleRate, 16, AUDIO_DEVICE_INDEX);
    [noisySpeechSamples, sampleRate] = audioread(noisyWavFilePath);
    noisySpeechAudioPlayer = audioplayer(noisySpeechSamples, sampleRate, 16, AUDIO_DEVICE_INDEX);
    createRecSpeechAudioPlayer = @(samples) audioplayer(samples, sampleRate, 16, AUDIO_DEVICE_INDEX);
    
    [numSamples, numChannels] = size(noisySpeechSamples);
    for channelIndex = 1:numChannels
        spectrograms(channelIndex,:,:) = stftFunction(noisySpeechSamples(:,channelIndex));
    end
    spectrograms = spectrograms / max(abs(spectrograms(:)));
    [numChannels, numFreq, numTime] = size(spectrograms);
    
    for channelIndex = 1:numChannels
        codebookCoefficients(channelIndex,:,:) = inferCodebookCoefficients( abs(squeeze(spectrograms(channelIndex,:,:))), codebook, 100, 0, 1e-16, 1969);
    end
    
    % Load target binaural cues
    targetDOA = getTargetDOA(noisyWavFilePath);
    [targetIPDs, targetILDs] = getAnechoicITFs(targetDOA, windowSize);

    % Compute input binaural cues
    ilds = getILDs(spectrograms);
    ipds = getIPDs(spectrograms);
    icms = getICMs(spectrograms);

    % Compute distance between input and target binaural cues
    ildDifferences = abs(targetILDs - ilds);
    ipdDifferences = abs(targetIPDs - ipds);
    
    % Compute atomic target binaural cue distances
    binauralCueDifferences = {ildDifferences, ipdDifferences, icms};
    numBinauralCues = length(binauralCueDifferences);
    for index = 1:numBinauralCues
        atomicTargetDistances{index} = codebook.' * binauralCueDifferences{index};
    end
    
    % Compute codebook masks
    for index = 1:numBinauralCues
        try 
            sliderSetting = get(codebookMaskSliders{index}, 'Value');
        catch
            sliderSetting = 0.5;
        end
        
        parametricASPPFunction = parametricASPPFunctions{index};
        codebookMasksList{index} = parametricASPPFunction(atomicTargetDistances{index},sliderSetting);
    end
    
    recomputeMasksAndReconstruction();
end

function [ipds] = getIPDs(spectrogram)
    ipds = squeeze( angle( spectrogram(1,:,:) .* conj(spectrogram(2,:,:)) ./ abs(spectrogram(1,:,:)) ./ abs(spectrogram(2,:,:)) ) );
end

function [ilds] = getILDs(spectrogram)
    ilds = squeeze( 20 * log10( abs(spectrogram(1,:,:)) ./ abs(spectrogram(2,:,:)) ) );
end

function [icms] = getICMs(spectrogram)
    averagingWindowSize = [1,4]; % [numFreq,numTime] mean-filter size

    crossPSD = squeeze( spectrogram(1,:,:) .* conj(spectrogram(2,:,:)) );
    autoPSDL = squeeze( spectrogram(1,:,:) .* conj(spectrogram(1,:,:)) );
    autoPSDR = squeeze( spectrogram(2,:,:) .* conj(spectrogram(2,:,:)) );
    
    if numel(averagingWindowSize) == 1
        meanFilter = ones(1,averagingWindowSize) / averagingWindowSize;
    else
        meanFilter = ones(averagingWindowSize) / prod(averagingWindowSize);
    end
    crossPSD = conv2(crossPSD, meanFilter, 'same');
    autoPSDL = conv2(autoPSDL, meanFilter, 'same');
    autoPSDR = conv2(autoPSDR, meanFilter, 'same');
    
    icms = abs( crossPSD ./ sqrt(autoPSDL .* autoPSDR) );
end

function [y] = gaussian(x,mu,sigma)
    y = 1 / sqrt(2*pi*sigma^2) * exp( -(x-mu).^2 / (2*sigma^2) );
end

function [y] = vonMises(x,mu,kappa)
    y = exp( kappa * cos(x-mu) ) / (2*pi*besseli(0,kappa));
end

function [aspp] = getGaussianASPP(targetDistances,sigma,uniformDistribution)
    lambda = uniformDistribution * gaussian(targetDistances,0,sigma);
    aspp = lambda ./ (1+lambda);
end

function [aspp] = getVonMisesASPP(targetDistances,kappa)
    lambda = 2*pi*vonMises(targetDistances,0,kappa);
    aspp = lambda ./ (1+lambda);
end

function [aspp] = getBetaASPP(targetDistances,alpha)
    beta = 2;
    lambda = 1 * betapdf(targetDistances,alpha,beta);
    aspp = lambda ./ (1+lambda);
end

function [] = recomputeMasksAndReconstruction()
    global spectrograms codebook codebookCoefficients codebookMasksList asppSimplexPoint wienerFilters recSpectrograms
    
    numMasks = length(codebookMasksList);
    codebookMask = asppSimplexPoint(1) * codebookMasksList{1};
    for index = 2:numMasks
        codebookMask = codebookMask + asppSimplexPoint(index) * codebookMasksList{index};
    end
    
    [numChannels, numFreq, numTime] = size(spectrograms);
    wienerFilters = zeros(numChannels, numFreq, numTime);
    numChannels = size(spectrograms,1);
    for channelIndex = 1:numChannels
        channelCoefficients = squeeze(codebookCoefficients(channelIndex,:,:));
        wienerFilters(channelIndex,:,:) = (codebook * (channelCoefficients .* codebookMask)) ./ (codebook * channelCoefficients);
    end
    recSpectrograms = spectrograms .* wienerFilters;
    recSpectrograms = recSpectrograms / max(abs(recSpectrograms(:)));
end

function [] = constructGrahpicalInterface(getParamForSliderValueFunctions, showSaveButton)
    global spectrograms wienerFilters recSpectrograms sliderUpdateMode codebookMasksList parametricASPPFunctions
    global inputUtteranceIDs inputNoiseTypes inputSNRs codebookMaskSliders paramStringTemplates
    
    binauralCueLabels = {'ILD-ASPP', 'IPD-ASPP','ICM-ASPP'};
    paramStringTemplates = {'\\sigma=%.4f', '\\kappa=%.4f', '\\alpha=%.4f, \\beta=1', '\\kappa=%.4f'};
    
    window = figure(1);
    set(gcf, 'Name', 'ASPP Demo');
    set(gcf, 'MenuBar', 'none');
    set(gcf, 'Toolbar', 'none');
    set(gcf, 'NumberTitle', 'off');
    
    colormap('jet')
    
    windowPadding = 2;
    horizontalSpacing = 20;
    ioWidthSpacing = 10;
    ioHeightSpacing = 10;
    ioLabelSizeProportion = 25;
    headerSizeProportion = 40;
    asppLabelSizeProportion = 16;
    spectrogramLimits = [-60, 0];
    wienerLimits = [-30, 0];
    
    % Main GUI structure
    mainBox = uix.HBox('Parent',window,'Padding',windowPadding,'Spacing',horizontalSpacing);
    asppBox = uix.VBox('Parent',mainBox,'Padding',0,'Spacing',0);
    simplexControlBox = uix.VBox('Parent',mainBox,'Padding',0,'Spacing',0);
    ioBox = uix.VBox('Parent',mainBox,'Padding',0,'Spacing',0);
    set(mainBox,'Widths',[-1,-1,-2]);
    
    % Input / output Visualization Box
    ioBoxHeader = uix.HBox('Parent',ioBox,'Padding',0,'Spacing',0);
    ioBoxMain = uix.VBox('Parent',ioBox,'Padding',0,'Spacing',ioHeightSpacing);
    ioBoxInput = uix.HBox('Parent',ioBoxMain,'Padding',0,'Spacing',0);
    ioBoxWiener = uix.HBox('Parent',ioBoxMain,'Padding',0,'Spacing',0);
    ioBoxOutput = uix.HBox('Parent',ioBoxMain,'Padding',0,'Spacing',0);
    set(ioBox,'Heights',[-1,-headerSizeProportion]);
    
    axes(ioBoxHeader,'Position',[0,0,1,1]);
    axis('off');
    axes(ioBoxHeader,'Position',[0,0,1,1]);
    text(0.5,0.5,'Left','HorizontalAlignment','center','VerticalAlignment','middle');
    axis('off');
    axes(ioBoxHeader,'Position',[0,0,1,1]);
    text(0.5,0.5,'Right','HorizontalAlignment','center','VerticalAlignment','middle');
    axis('off');
    set(ioBoxHeader,'Widths',[-1,-ioLabelSizeProportion/2, -ioLabelSizeProportion/2]);
    
    axes(ioBoxInput,'Position',[0,0,1,1]);
    text(0.5,0.5,'Noisy Speech','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
    axis('off');
    ioBoxSpectrograms = uix.HBox('Parent',ioBoxInput,'Padding',0,'Spacing',ioWidthSpacing);
    ax = axes('Parent',ioBoxSpectrograms,'Position',[0,0,1,1]);
    spectrogramLImage = imshow(ax, 20*log10(abs(squeeze(spectrograms(1,:,:)))), spectrogramLimits);
    ax = axes('Parent',ioBoxSpectrograms,'Position',[0,0,1,1]);
    spectrogramRImage = imshow(ax, 20*log10(abs(squeeze(spectrograms(2,:,:)))), spectrogramLimits);
    set(ioBoxInput,'Widths',[-1,-ioLabelSizeProportion]);
    
    axes(ioBoxWiener,'Position',[0,0,1,1]);
    text(0.5,0.5,'Wiener Filter','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
    axis('off');
    ioBoxWienerFilters = uix.HBox('Parent',ioBoxWiener,'Padding',0,'Spacing',ioWidthSpacing);
    ax = axes('Parent',ioBoxWienerFilters,'Position',[0,0,1,1]);
    wienerFilterLImage = imshow(ax, 20*log10(squeeze(wienerFilters(1,:,:))), wienerLimits);
    ax = axes('Parent',ioBoxWienerFilters,'Position',[0,0,1,1]);
    wienerFilterRImage = imshow(ax, 20*log10(squeeze(wienerFilters(2,:,:))), wienerLimits);
    set(ioBoxWiener,'Widths',[-1,-ioLabelSizeProportion]);
    
    axes(ioBoxOutput,'Position',[0,0,1,1]);
    text(0.5,0.5,'Speech Estimate','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
    axis('off');
    ioBoxOutputSpectrograms = uix.HBox('Parent',ioBoxOutput,'Padding',0,'Spacing',ioWidthSpacing);
    ax = axes('Parent',ioBoxOutputSpectrograms,'Position',[0,0,1,1]);
    recSpectrogramLImage = imshow(ax, 20*log10(abs(squeeze(recSpectrograms(1,:,:)))), spectrogramLimits);
    ax = axes('Parent',ioBoxOutputSpectrograms,'Position',[0,0,1,1]);
    recSpectrogramRImage = imshow(ax, 20*log10(abs(squeeze(recSpectrograms(2,:,:)))), spectrogramLimits);
    set(ioBoxOutput,'Widths',[-1,-ioLabelSizeProportion]);
    
    % Simplex Controls Box
    axes(simplexControlBox);
    axis('off');
    
    global simplexSelectionAx
    simplexSelectionAx = axes('Parent',simplexControlBox,'Position',[0,0,1,1]);
    setASPPSelection(0,0)
    axes(simplexControlBox);
    axis('off');
    
    global playCleanButton playNoisyButton playEnhancedButton
    playButtonsBox = uix.HBox('Parent',simplexControlBox,'Padding',0,'Spacing',10);
    playCleanButton = uicontrol('Parent',playButtonsBox,'Style','togglebutton','String','Play Clean','Position',[0,0,1,1]);
    playNoisyButton = uicontrol('Parent',playButtonsBox,'Style','togglebutton','String','Play Noisy','Position',[0,0,1,1]);
    playEnhancedButton = uicontrol('Parent',playButtonsBox,'Style','togglebutton','String','Play Enhanced','Position',[0,0,1,1]);
    if showSaveButton
        saveButton = uicontrol('Parent',playButtonsBox,'Style','togglebutton','String','Save','Position',[0,0,1,1]);
    end
    
    axes(simplexControlBox);
    axis('off');
    
    inputParametersControlsBox = uix.VBox('Parent',simplexControlBox,'Padding',0,'Spacing',0);
    inputParametersLabelsBox = uix.HBox('Parent',inputParametersControlsBox,'Padding',0,'Spacing',0);
    uicontrol('Parent',inputParametersLabelsBox,'Style','text','String','Utterance:');%,'Position',[0,0,1,1]);
    uicontrol('Parent',inputParametersLabelsBox,'Style','text','String','Noise type:','Position',[0,0,1,1]);
    uicontrol('Parent',inputParametersLabelsBox,'Style','text','String','SNR:','Position',[0,0,1,1]);
    inputParametersDropDownsBox = uix.HBox('Parent',inputParametersControlsBox,'Padding',0,'Spacing',10);
    global utteranceDropDown noiseTypeDropDown snrDropDown
    utteranceDropDown = uicontrol('Parent',inputParametersDropDownsBox,'Style','popupmenu','String',inputUtteranceIDs,'Value',1,'Position',[0,0,1,1]);
    noiseTypeDropDown = uicontrol('Parent',inputParametersDropDownsBox,'Style','popupmenu','String',inputNoiseTypes,'Value',1,'Position',[0,0,1,1]);
    snrDropDown = uicontrol('Parent',inputParametersDropDownsBox,'Style','popupmenu','String',inputSNRs,'Value',1,'Position',[0,0,1,1]);
    
    axes(simplexControlBox);
    axis('off');
    
    set(simplexControlBox,'Heights',[-2,-8,-2,-1,-1,-1,-4]);
    
    % ASPP Visualization Box
    asppHeaderBox = uix.HBox('Parent',asppBox,'Padding',0,'Spacing',0);
    axes(asppHeaderBox,'Position',[0,0,1,1]);
    axis('off')
    axes(asppHeaderBox,'Position',[0,0,1,1]);
    text(0.5,0.5,'Atomic Speech Presence Probability (ASPP) Estimates','HorizontalAlignment','center','VerticalAlignment','middle');
    axis('off');
    set(asppHeaderBox,'Widths',[-1,-asppLabelSizeProportion]);
    
    asppBoxMain = uix.VBox('Parent',asppBox,'Padding',0,'Spacing',0);
    set(asppBox,'Heights',[-1,-headerSizeProportion]);
    
    codebookMaskImages = {};
    codebookMaskSliders = {};
    numCols = length(codebookMasksList);
    for index = 1:numCols
        currentASPPBox = uix.HBox('Parent',asppBoxMain,'Padding',0,'Spacing',0);
        
        binauralCueLabelBox = uix.VBox('Parent',currentASPPBox,'Padding',0,'Spacing',0);
        axes(binauralCueLabelBox,'Position',[0,0,1,1]);
        text(0.5,0.5,binauralCueLabels{index},'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
        axis('off');
        axes(binauralCueLabelBox,'Position',[0,0,1,1]);
        axis('off');
        axes(binauralCueLabelBox,'Position',[0,0,1,1]);
        axis('off');
        set(binauralCueLabelBox,'Heights',[-10,-1,-1]);
        
        currentASPPMainBox = uix.VBox('Parent',currentASPPBox,'Padding',0,'Spacing',0);
        set(currentASPPBox,'Widths',[-1,-asppLabelSizeProportion]);
        
        ax = axes('Parent',currentASPPMainBox,'Position',[0,0,1,1]);
        codebookMaskImage = imshow(ax, codebookMasksList{index}, [0,1]);
        codebookMaskImages{index} = codebookMaskImage;
        
        slider = uicontrol('Parent',currentASPPMainBox,'Style','slider','Position',[0,0,1,1],'value',0.5,'min',0,'max',1);
        codebookMaskSliders{index} = slider;
        
        paramStringTemplate = paramStringTemplates{index};
        paramValueFunction = getParamForSliderValueFunctions{index};
        paramString = sprintf(paramStringTemplate, paramValueFunction(0.5));
        axes(currentASPPMainBox,'Position',[0,0,1,1]);
        textObject = text(0.5,0.5,paramString,'HorizontalAlignment','center','VerticalAlignment','middle');
        axis('off');
        
        set(currentASPPMainBox,'Heights',[-10,-1,-1]);
        
        if strcmp(sliderUpdateMode,'continuous')
            addlistener(slider,'ContinuousValueChange',@(es, ed) updateSlider(es, codebookMaskImage, wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage, textObject, ...
                                                                              index, parametricASPPFunctions{index}, getParamForSliderValueFunctions{index}, paramStringTemplates{index}) );
        else
            slider.Callback = @(es, ed) updateSlider(es, codebookMaskImage, wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage, textObject, ...
                                                     index, parametricASPPFunctions{index}, getParamForSliderValueFunctions{index}, paramStringTemplates{index});
        end
    end
    
    playCleanButton.Callback = @(src,event) playCleanButtonPushed();
    playNoisyButton.Callback = @(src,event) playNoisyButtonPushed();
    playEnhancedButton.Callback = @(src,event) playEnhancedButtonPushed();
    if showSaveButton
        saveButton.Callback = @(src,event) saveExampleButtonPushed();
    end
    
    utteranceDropDown.Callback = @(es, ed) updateInputDropDowns(codebookMaskImages, spectrogramLImage, spectrogramRImage, wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage);
    noiseTypeDropDown.Callback = @(es, ed) updateInputDropDowns(codebookMaskImages, spectrogramLImage, spectrogramRImage, wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage);
    snrDropDown.Callback = @(es, ed) updateInputDropDowns(codebookMaskImages, spectrogramLImage, spectrogramRImage, wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage);
    
    global MOUSE_BUTTON_DOWN
    MOUSE_BUTTON_DOWN = 0;
    set(gcf, 'WindowButtonMotionFcn', @(object, eventdata) mouseMove(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage));
    set(gcf, 'WindowButtonDownFcn', @(object, eventdata) mouseDown(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage));
    set(gcf, 'WindowButtonUpFcn', @(object, eventdata) mouseUp(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage));
end

function [] = updateInputDropDowns(codebookMaskImages, spectrogramLImage, spectrogramRImage, wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage);
	global utteranceDropDown noiseTypeDropDown snrDropDown cleanSpeechAudioPlayer noisySpeechAudioPlayer recSpeechAudioPlayer playCleanButton playNoisyButton playEnhancedButton
    
    set(gcf, 'pointer', 'watch');
    drawnow;
    
    % stop any playing audio
    if cleanSpeechAudioPlayer.isplaying
        set(playCleanButton,'Value',0)
        stop(cleanSpeechAudioPlayer)
    end
    if noisySpeechAudioPlayer.isplaying
        set(playNoisyButton,'Value',0)
        stop(noisySpeechAudioPlayer)
    end
    if ~isempty(recSpeechAudioPlayer) && recSpeechAudioPlayer.isplaying
        set(playEnhancedButton,'Value',0)
        stop(recSpeechAudioPlayer)
    end
    
    utteranceIDs = get(utteranceDropDown,'String');
    utteranceID = utteranceIDs{ get(utteranceDropDown,'Value') };
    noiseTypes = get(noiseTypeDropDown,'String');
    noiseType = noiseTypes{ get(noiseTypeDropDown, 'Value') };
    snrs = get(snrDropDown,'String');
    snr = snrs{ get(snrDropDown, 'Value') };

    setInputExample(utteranceID, noiseType, snr);
    updateAllImages(codebookMaskImages, spectrogramLImage, spectrogramRImage, wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage);
    
    set(gcf, 'pointer', 'arrow');
end

function [] = playCleanButtonPushed()
    global cleanSpeechAudioPlayer noisySpeechAudioPlayer recSpeechAudioPlayer playCleanButton playNoisyButton playEnhancedButton
    
    if ~isempty(recSpeechAudioPlayer) && recSpeechAudioPlayer.isplaying
        set(playEnhancedButton,'Value',0)
        stop(recSpeechAudioPlayer)
    end
    
    if noisySpeechAudioPlayer.isplaying
        set(playNoisyButton,'Value',0)
        stop(noisySpeechAudioPlayer)
    end
    
    if cleanSpeechAudioPlayer.isplaying
        set(playCleanButton,'Value',0)
        stop(cleanSpeechAudioPlayer)
    else
        set(playCleanButton,'Value',1)
        play(cleanSpeechAudioPlayer)
    end
end

function [] = playNoisyButtonPushed()
    global cleanSpeechAudioPlayer noisySpeechAudioPlayer recSpeechAudioPlayer playCleanButton playNoisyButton playEnhancedButton
    
    if ~isempty(recSpeechAudioPlayer) && recSpeechAudioPlayer.isplaying
        set(playEnhancedButton,'Value',0)
        stop(recSpeechAudioPlayer)
    end
    
    if cleanSpeechAudioPlayer.isplaying
        set(playCleanButton,'Value',0)
        stop(cleanSpeechAudioPlayer)
    end
    
    if noisySpeechAudioPlayer.isplaying
        set(playNoisyButton,'Value',0)
        stop(noisySpeechAudioPlayer)
    else
        set(playNoisyButton,'Value',1)
        play(noisySpeechAudioPlayer)
    end
end

function [] = playEnhancedButtonPushed()
    global recSamples cleanSpeechAudioPlayer noisySpeechAudioPlayer playCleanButton playNoisyButton playEnhancedButton recSpeechAudioPlayer createRecSpeechAudioPlayer
    
    if cleanSpeechAudioPlayer.isplaying
        set(playCleanButton,'Value',0)
        stop(cleanSpeechAudioPlayer)
    end
    
    if noisySpeechAudioPlayer.isplaying
        set(playNoisyButton,'Value',0)
        stop(noisySpeechAudioPlayer)
    end

    if ~isempty(recSpeechAudioPlayer) && recSpeechAudioPlayer.isplaying
        set(playEnhancedButton,'Value',0)
        stop(recSpeechAudioPlayer)
    else
        set(playEnhancedButton,'Value',1)
        updateRecSamples()
        recSpeechAudioPlayer = createRecSpeechAudioPlayer(recSamples);
        play(recSpeechAudioPlayer)
    end
end

function [] = updateRecSamples()
    clear global recSamples
    global recSamples recSpectrograms inverseSTFTFunction
    for channelIndex = 1:size(recSpectrograms,1)
        recSamples(channelIndex,:) = inverseSTFTFunction(squeeze(recSpectrograms(channelIndex,:,:)));
    end
    recSamples  = recSamples / max(abs(recSamples(:))) * 0.9;
end

function [] = saveExampleButtonPushed()
    updateRecSamples();
    global ASPP_DIR_PATH cleanSpeechSamples noisySpeechSamples recSamples sampleRate utteranceDropDown noiseTypeDropDown snrDropDown asppSimplexPoint
    global codebookMaskSliders getParamForSliderValueFunctions paramStringTemplates
    
    outputDir = [ASPP_DIR_PATH filesep 'output'];
    status = mkdir(outputDir);
    
    utteranceIDs = get(utteranceDropDown,'String');
    utteranceID = utteranceIDs{ get(utteranceDropDown,'Value') };
    noiseTypes = get(noiseTypeDropDown,'String');
    noiseType = noiseTypes{ get(noiseTypeDropDown, 'Value') };
    snrs = get(snrDropDown,'String');
    snr = snrs{ get(snrDropDown, 'Value') };
    
    paramValueStrings = {};
    numParams = length(codebookMaskSliders);
    for index = 1:numParams
        codebookMaskSlider = codebookMaskSliders{index};
        getParamForSliderValueFunction = getParamForSliderValueFunctions{index};
        paramStringTemplate = paramStringTemplates{index};
        sliderValue = codebookMaskSlider.Value;
        paramValue = getParamForSliderValueFunction(sliderValue);
        paramValueString = sprintf(paramStringTemplate, paramValue);
        paramValueString = strrep(paramValueString, '\', '');
        paramValueString = strrep(paramValueString, ' ', '');
        paramValueStrings{index} = paramValueString;
    end
    paramValuesString = string(join(paramValueStrings, '_'));
    paramValuesString = paramValuesString{1};
    
    asppString = sprintf('[%.4f,%.4f,%.4f]', asppSimplexPoint(1), asppSimplexPoint(2), asppSimplexPoint(3));
    
    fileNameBase = [utteranceID '_' noiseType '_' snr '_' asppString '_' paramValuesString];
    
    cleanFilePath = [outputDir filesep [fileNameBase '_Clean.wav']];
    noisyFilePath = [outputDir filesep [fileNameBase '_Noisy.wav']];
    enhancedFilePath = [outputDir filesep [fileNameBase '_Enhanced.wav']];
    
    audiowrite(cleanFilePath, cleanSpeechSamples, sampleRate);
    audiowrite(noisyFilePath, noisySpeechSamples, sampleRate);
    audiowrite(enhancedFilePath, recSamples.', sampleRate);
end

function [] = updateSlider(es, codebookMaskImage, wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage, paramValueText, ...
                           codebookMaskIndex, asppFunction, getParamForSliderValueFunction, paramStringTemplate)
    global codebookMasksList atomicTargetDistances
    
    sliderValue = es.Value;
    paramValue = getParamForSliderValueFunction(sliderValue);
    set(paramValueText, 'String', sprintf(paramStringTemplate, paramValue));
    
    codebookMasksList{codebookMaskIndex} = asppFunction(atomicTargetDistances{codebookMaskIndex}, sliderValue);
    set(codebookMaskImage, 'CData', codebookMasksList{codebookMaskIndex});
    
    updateReconstructionImages(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage)
end

function [] = updateReconstructionImages(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage)
    global wienerFilters recSpectrograms
    recomputeMasksAndReconstruction();
    
    set(wienerFilterLImage, 'CData', 20*log10(squeeze(wienerFilters(1,:,:))));
    set(wienerFilterRImage, 'CData', 20*log10(squeeze(wienerFilters(2,:,:))));
    set(recSpectrogramLImage, 'CData', 20*log10( abs(squeeze(recSpectrograms(1,:,:)))) );
    set(recSpectrogramRImage, 'CData', 20*log10( abs(squeeze(recSpectrograms(2,:,:)))) );
end

function [] = updateAllImages(codebookMaskImages, spectrogramLImage, spectrogramRImage, wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage);
    global spectrograms wienerFilters recSpectrograms codebookMasksList
    
    numFrames = size(spectrograms,3);
    numCodebookMasks = length(codebookMasksList);
    for maskIndex = 1:numCodebookMasks
        set(codebookMaskImages{maskIndex}, 'CData', codebookMasksList{maskIndex});
        set(codebookMaskImages{maskIndex}.Parent, 'XLim', [1,numFrames]);
    end
    
    set(spectrogramLImage, 'CData', 20*log10(squeeze(abs(spectrograms(1,:,:)))));
    set(spectrogramLImage.Parent, 'XLim', [1,numFrames]);
    set(spectrogramRImage, 'CData', 20*log10(squeeze(abs(spectrograms(2,:,:)))));
    set(spectrogramRImage.Parent, 'XLim', [1,numFrames]);
    set(wienerFilterLImage, 'CData', 20*log10(squeeze(wienerFilters(1,:,:))));
    set(wienerFilterLImage.Parent, 'XLim', [1,numFrames]);
    set(wienerFilterRImage, 'CData', 20*log10(squeeze(wienerFilters(2,:,:))));
    set(wienerFilterRImage.Parent, 'XLim', [1,numFrames]);
    set(recSpectrogramLImage, 'CData', 20*log10( abs(squeeze(recSpectrograms(1,:,:)))) );
    set(recSpectrogramLImage.Parent, 'XLim', [1,numFrames]);
    set(recSpectrogramRImage, 'CData', 20*log10( abs(squeeze(recSpectrograms(2,:,:)))) );
    set(recSpectrogramRImage.Parent, 'XLim', [1,numFrames]);
end

function [image] = imshow(ax, x, clim)
    image = imagesc(ax, x);
    if nargin == 3
        caxis(clim)
    end
    set(ax,'xtick',[])
    set(ax,'ytick',[])
    set(ax,'YDir','normal')
end

function [] = mouseDown(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage)
    global MOUSE_BUTTON_DOWN
    MOUSE_BUTTON_DOWN = 1;
    [x,y] = getCurrentSimplexSelection();
    setASPPSelection(x,y);
    updateReconstructionImages(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage);
end

function [] = mouseUp(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage)
    global MOUSE_BUTTON_DOWN
    MOUSE_BUTTON_DOWN = 0;
    [x,y] = getCurrentSimplexSelection();
    setASPPSelection(x,y);
    updateReconstructionImages(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage);
end

function [] = mouseMove(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage)
    global MOUSE_BUTTON_DOWN
    if MOUSE_BUTTON_DOWN == 0
        return
    end
    [x,y] = getCurrentSimplexSelection();
    setASPPSelection(x,y);
    updateReconstructionImages(wienerFilterLImage, wienerFilterRImage, recSpectrogramLImage, recSpectrogramRImage);
end

function [x,y] = getCurrentSimplexSelection()
    C = get(gca, 'CurrentPoint');
    x = C(1,1);
    y = C(1,2);
end

function [] = setASPPSelection(x,y)
    [a1,a2,a3] = getSimplexPoint(x,y);
    invalidPoint = a1<0 || a2 <0 || a3 <0;
    if invalidPoint
        return
    end
    
    global asppSimplexPoint simplexSelectionAx
    asppSimplexPoint = [a1,a2,a3];
    
    simplex = polyshape([0,0.5,1,0],[0,sqrt(3)/2,0,0]);
    plot(simplexSelectionAx, simplex, 'FaceColor','white')
    hold('on')
    scatter(simplexSelectionAx, x, y, 40, 'b', 'filled', 'Marker', 'o')
    hold('off')

    ylim( [-0.05, sqrt(3)/2+0.05] )
    xlim( [-0.05, 1.05] )
    titleString = sprintf('ILD, IPD, ICM:\n(%.3f, %.3f, %.3f)', a1, a2, a3);
    xlabel(gca, titleString, 'FontSize', 12, 'FontWeight','bold');
    text(0,-0.05,'ILD','HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize',12);
    text(1,-0.05,'IPD','HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize',12);
    text(0.5, sqrt(3)/2+0.04, 'ICM','HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize',12);
    set(gca,'visible','off')
    set(findall(gca, 'type', 'text'), 'visible', 'on')
end

function [a1,a2,a3] = getSimplexPoint(x,y)
    p.x = x;
    p.y = y;
    p0.x = 0;
    p0.y = 0;
    p1.x = 1;
    p1.y = 0;
    p2.x = 0.5;
    p2.y = sqrt(3)/2;
    
    Area = 0.5 *(-p1.y*p2.x + p0.y*(-p1.x + p2.x) + p0.x*(p1.y - p2.y) + p1.x*p2.y);
    a3 = 1/(2*Area)*(p0.x*p1.y - p0.y*p1.x + (p0.y - p1.y)*p.x + (p1.x - p0.x)*p.y);
    a2 = 1/(2*Area)*(p0.y*p2.x - p0.x*p2.y + (p2.y - p0.y)*p.x + (p0.x - p2.x)*p.y);
    a1 = 1 - a3 - a2;
end