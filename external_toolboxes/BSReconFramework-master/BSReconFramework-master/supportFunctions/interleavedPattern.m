function [P, ranges] = interleavedPattern(data, AF, ACL, symmetric, PhaseEncDir)


if nargin < 4
    symmetric = false;
    if nargin < 3
        ACL = 0;    % autocalibration lines
        if nargin < 2
            AF = 2;     % acceleration factor
        end
    end
end


sizeData = size(data);
switch PhaseEncDir
    case 'COL'
        dimPhase = sizeData(1);
    case 'ROW'
        dimPhase = sizeData(2);
    otherwise
        error('Phase encoding direction error')
end
dimSlice = sizeData(3);
noEchoes = sizeData(4);

if length(sizeData) > 5
    error('interleavedPattern:argChk','dataset dimension is not supported yet');
end

P = zeros(size(data));
ranges = cell(noEchoes,1);
start = 1;

% derive range of autocalibration lines
refRange = (dimPhase/2 - ACL/2) + 1 : (dimPhase/2 + ACL/2);
minSlice = dimSlice/2-ACL/2 + 1;
maxSlice = dimSlice/2+ACL/2; 

for n = 1 : noEchoes
    startSlice = start;
    for iSlice = 1:dimSlice
        
        % derive range of undersampling pattern
        if symmetric
            AFRange1 = startSlice:AF:dimPhase/2;
            AFRange2 = dimPhase - (startSlice - 1) : -AF : dimPhase/2;
            if iSlice >= minSlice && iSlice <= maxSlice
                range = [refRange, AFRange1, AFRange2];
            else
                range = [AFRange1, AFRange2];
            end
        else
            AFRange = startSlice:AF:dimPhase;
            if iSlice >= minSlice && iSlice <= maxSlice
                range = [refRange, AFRange];
            else
                range = AFRange;
            end
        end
        
        range = removeDuplicates(range);
        
        switch PhaseEncDir
            case 'COL'
                P(range, :, iSlice, n, :) = data(range, :, iSlice, n, :);
            case 'ROW'
                P(:, range, iSlice, n, :) = data(:, range, iSlice, n, :);
        end
        if startSlice >= AF
            startSlice = 1;
        else
            startSlice = startSlice + 1;
        end
    end
    %     P(range, :, n, :) = data(range, :, n, :);
    if start >= AF
        start = 1;
    else
        start = start + 1;
    end
    ranges{n} = range;
end


    function range = removeDuplicates(range)
        range = sort(range);
        i = 2;
        while i <= length(range)
            if range(i) == range(i-1)
                range(i-1) = [];
            end
            i = i + 1;
        end
    end
end