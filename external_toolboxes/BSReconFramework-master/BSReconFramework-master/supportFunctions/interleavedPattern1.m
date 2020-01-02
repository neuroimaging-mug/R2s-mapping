function [P, ranges] = interleavedPattern1(data, AF, ACL, symmetric)


if nargin < 4
    symmetric = false;
    if nargin < 3
        ACL = 0;    % autocalibration lines
        if nargin < 2
            AF = 2;     % acceleration factor
        end
    end
end


si = size(data);
dimY = si(1);
noEchoes = si(3);

if length(si) > 4
    error('interleavedPattern:argChk','dataset dimension is not supported yet');
end

P = zeros(size(data));
ranges = cell(noEchoes,1);
start = 1;

% derive range of autocalibration lines

for n = 1 : noEchoes
    refRange = (dimY/2 - ACL/2) + 1 : (dimY/2 + ACL/2);
    
    % derive range of undersampling pattern
    if symmetric
        %if mod(n,2) == 1
            %AFRange1 = start:AF:(dimY/2);
            %AFRange2 = dimY - (start - 1) : -AF : (dimY/2);
            range = refRange;%, AFRange1, AFRange2];
       % else
           % AFRange1 = start:AF:dimY/2;
            %AFRange2 = dimY - (start - 1) : -AF : dimY/2;
        %    range = [AFRange1, AFRange2, dimY/2-1, dimY/2, dimY/2+1];
       % end
    else
        if mod(n,2) == 1
            AFRange = start:AF:dimY;
            range = [refRange, AFRange];
        else
            range = start:AF:dimY;
        end
    end
    
    range = removeDuplicates(range);
    
    P(range, :, n, :) = data(range, :, n, :);
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