function [targetDOA] = getTargetDOA(wavFilePath)

global ASPP_DIR_PATH
HRIR_SISEC_SOURCE_DIR = [ASPP_DIR_PATH filesep 'resources' filesep 'HRIR_TIMIT_Source'];
CAFETERIA_SOURCE_POSITION_ANGLE_DICT.A = 0;
CAFETERIA_SOURCE_POSITION_ANGLE_DICT.B = -30;
CAFETERIA_SOURCE_POSITION_ANGLE_DICT.C = -90;
CAFETERIA_SOURCE_POSITION_ANGLE_DICT.D = 90;

[~, wavFileBasename, ~] = fileparts(wavFilePath);
targetPositionFilePath = [HRIR_SISEC_SOURCE_DIR filesep strrep(wavFileBasename, 'mix', 'POS') '.txt'];
fileID = fopen(targetPositionFilePath, 'r');
targetPositionString = textscan(fileID, '%s');
targetPositionString = char(targetPositionString{1});
targetDOA = CAFETERIA_SOURCE_POSITION_ANGLE_DICT.(targetPositionString);
fclose(fileID);

end

