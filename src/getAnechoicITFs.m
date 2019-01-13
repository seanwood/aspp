function [ipds, ilds] = getAnechoicITFs(angleOfArrival, windowSize)

hrirs = getAnechoicHRIRs(angleOfArrival);
hrirFFTs(1,:) = rfft(hrirs(1,:), windowSize);
hrirFFTs(2,:) = rfft(hrirs(2,:), windowSize);

ipds = squeeze( angle(  hrirFFTs(1,:) .* conj(hrirFFTs(2,:)) ./ abs(hrirFFTs(1,:)) ./ abs(hrirFFTs(2,:)) ) ).';
ilds = squeeze( 20 * log10( abs(hrirFFTs(1,:)) ./ abs(hrirFFTs(2,:)) ) ).';

end