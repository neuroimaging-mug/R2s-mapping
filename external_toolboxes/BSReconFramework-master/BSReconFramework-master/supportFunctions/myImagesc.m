function myImagesc(ah,varargin)
if ishandle(ah)
    imagesc(varargin{:},'parent',ah);
%     set(ah,'DataAspectRatio',[1,1,1]);
    axis(ah,'off');
    title(ah,inputname(2));
end
end
