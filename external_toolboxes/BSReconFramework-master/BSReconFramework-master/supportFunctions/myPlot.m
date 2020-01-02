function myPlot(ah,varargin)
    if ishandle(ah) 
        plot(ah,varargin{:});
    end
end