function [hrirs] = getAnechoicHRIRs(azimuth, sampleRate)
    if nargin < 2
        sampleRate = 16000;
    end

    global ASPP_DIR_PATH
    ANECHOIC_HRIR_DIR = [ASPP_DIR_PATH filesep 'resources' filesep 'HRIR' filesep 'HRIR_database_wav' filesep 'hrir' filesep 'anechoic'];
    ANECHOIC_HRIR_FILE_NAME_TEMPLATE = 'anechoic_distcm_80_el_0_az_%d.wav';
    HRIR_FRONT_CHANNEL_INDEXES = [3,4];

    hrirFilePath = [ANECHOIC_HRIR_DIR filesep sprintf(ANECHOIC_HRIR_FILE_NAME_TEMPLATE, azimuth)];
    [inputHRIRSamples, hrirSampleRate] = audioread(hrirFilePath);
    sampleRateRatio = hrirSampleRate / sampleRate;

    for hrirIndex = 1:length(HRIR_FRONT_CHANNEL_INDEXES)
        channelIndex = HRIR_FRONT_CHANNEL_INDEXES(hrirIndex);
        hrirs(hrirIndex,:) = decimate(inputHRIRSamples(:,channelIndex), sampleRateRatio);
    end
end