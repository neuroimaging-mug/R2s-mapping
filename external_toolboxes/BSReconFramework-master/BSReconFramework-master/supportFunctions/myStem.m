function myStem(ah, varargin)
    if ishandle(ah)
        stem(ah,varargin{:});
    end
end