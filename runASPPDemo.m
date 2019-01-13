function [] = runASPPDemo(audioDeviceIndex)
    global ASPP_DIR_PATH AUDIO_DEVICE_INDEX

    [ASPP_DIR_PATH,~,~] = fileparts(mfilename('fullpath'));
    addpath([ASPP_DIR_PATH filesep 'src']);

    if nargin >= 1
        AUDIO_DEVICE_INDEX = audioDeviceIndex;
    elseif ~exist('AUDIO_DEVICE_INDEX', 'var') || isempty(AUDIO_DEVICE_INDEX)
        AUDIO_DEVICE_INDEX = -1; % Default audio device
    end
    printAudioDevices(AUDIO_DEVICE_INDEX)
    
    runASPP();
end