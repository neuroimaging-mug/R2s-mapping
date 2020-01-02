function [dDataOut, sName, sUnitFormat] = fSNR(dData)

sName = 'SNR';
sUnitFormat = '';
dDataOut = mean(dData)./std(dData);
% =========================================================================
% *** END OF FUNCTION fEvalROIMean
% =========================================================================