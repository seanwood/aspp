function [] = printAudioDevices(audioDeviceIndex)

audioDeviceInfo = audiodevinfo;
fprintf('Available audio devices:\n');
for deviceInfo = audioDeviceInfo.output
    if deviceInfo.ID == audioDeviceIndex
        fprintf('* ');
    end
    fprintf('Device index %d: %s\n', deviceInfo.ID, deviceInfo.Name);
end