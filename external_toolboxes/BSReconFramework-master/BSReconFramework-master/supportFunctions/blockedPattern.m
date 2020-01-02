function [uData, sampleIndices] = blockedPattern(data, varargin)

si = size(data);
if length(si) >4 || length(si) < 3
    error('blockedPattern:argChk','only 3D or 4D data array supported yet');
end

[dimY, dimX, NI, NC] = size(data);


if nargin == 2
    % if we have only 2 arguments,
    % interpret second argument as acceleration factor
    AF = varargin{1};
else
    AF = NI;
end


%---------------------------------
% evaluate varargin
%---------------------------------
deactPe     = false;
sortMethod  = 'linear';
startAt     = 'top';
nTE         = 0;
if nargin > 2
    for i = 1 : floor(length(varargin)/2)
        option       = varargin{i*2-1};
        option_value = varargin{i*2};
        
        switch lower(option)
            case {'af'}   % acceleration factor
                if option_value <= NI
                    AF = option_value;
                else
                    error('blockedPattern:argChk','rare factor cannot be greater than number of images');
                end
            case 'deactpe'
                deactPe = option_value;
            case {'order','ordering','sortmethod', 'method'}
                % can be:
                %   1  linear,
                %   2  centric,
                sortMethod = option_value;
            case 'startat'
                startAt = option_value; % this is relevant for linear ordering
                % options are
                % 'top', 'center', 'centerblock'
                
            otherwise
                error('blockedPattern:argChk','unknown option [%s]!\n',option);
        end
    end
end
if nTE > 0
    sortMethod = 'nte';
end
if mod(dimY, AF)
    warning('blockedPattern:argChk','Modulus(rareFactor, dimY) should be zero');
end

%---------------------------------
% simulate deactivation of phase encoding gradients
%---------------------------------
if deactPe
    yCenter = ceil(dimY/2);
    data = repmat(data(yCenter,:,:,:),[dimY,1,1,1]);
end



%---------------------------------
% sort echo train to form one single RARE image
%---------------------------------
sampleIndices = cell([NI,1]);
uData = zeros([dimY, dimX, NI, NC]);

switch sortMethod
    
    case {'linear', 'overl'}
        % simple 'stairway'
        blockSize = ceil(dimY / AF);
        pf = false;
        
        switch startAt
            case 'top'
                start = 1;
            case 'center'
                start = ceil(dimY/2) + 1;
            case 'centerm1b'
                start = ceil(dimY/2) + 1 - blockSize;
            case 'centerm1pf' % one block above center
                center = ceil(dimY/2) + 1;
                start = center;
                pf = true;
            case 'centerpf'
                center = ceil(dimY/2) + 1;
                start = center;
                pf = true;
            case 'centerblock'
                start = ceil(dimY/2) + 1 - ceil(blockSize / 2);
            case 'centerblockm1'
                start = ceil(dimY/2) + 1 - ceil(3*blockSize / 2);
                if start < 1
                    start = start + dimY;
                end
            case 'centerblockm2'
                start = ceil(dimY/2) + 1 - ceil(5*blockSize / 2);
                if start < 1
                    start = start + dimY;
                end
                
            case 'centerblockp1'
                start = ceil(dimY/2) + 1 + ceil(blockSize / 2);
                if start < 1
                    start = start + dimY;
                end                
            case 'centerblockpf'
                center = ceil(dimY/2) + 1 - ceil(blockSize / 2);
                start = center;
                pf = true;
                
        end
        
        
        if strcmp(sortMethod, 'overl')
            divi = 2;
        else
            divi = 1;
        end
        
            
        for n = 1 : NI
            stop = start + blockSize - 1;
            if stop > dimY
                % wrap around case
                if pf
                    stop = stop - dimY + center - 1;
                    topRange = center : stop ;
                else
                    stop = stop - dimY;
                    topRange = 1 : stop ;
                end
                botRange = start : dimY;

                range = [topRange, botRange];

                % next start
                start = stop + 1;
            else
                % normal case
                range = start : stop;

                % next start
                if stop >= dimY
                    if pf
                        start = center;
                    else
                        start = 1;
                    end
                else
                    start = start + floor(blockSize/divi);
                end

            end
            uData(range,:,n,:) =  data(range,:,n,:);
            sampleIndices{n} = range;                
        end

                        
            
        
        
    case 'centric'
        % put first echo to center (spin density weighted)
        blockSize     = ceil(dimY / AF);
        halfBlockSize = ceil(blockSize/2);
        
        top = ceil(dimY/2);
        bot = ceil(dimY/2) + 1;
        for n = 1 : NI
            topRange =  top - halfBlockSize + 1 : top;
            botRange =  bot : bot + halfBlockSize -1;
            
            topRange(topRange < 1) = topRange(topRange < 1) + ceil(dimY/2);
            botRange(botRange > dimY) = botRange(botRange > dimY) - (ceil(dimY/2));
            
            range = [topRange, botRange];
            
            uData(range,:,n,:) =  data(range,:,n,:);
            
            top = top - halfBlockSize;
            bot = bot + halfBlockSize;
            if top < 1
                % restart
                top = top + ceil(dimY/2);
                bot = bot - (ceil(dimY/2));
            end
            sampleIndices{n} = range;
        end
        

        
        
        
    otherwise
        error('blockedPattern:method','unknown sort method');
        
end
end