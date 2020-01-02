function str = clockStr(dateOrTime)
if nargin < 1
    dateOrTime = 'both';
end

c = fix(clock);

d = sprintf('%d.%d.%d',c(3),c(2),c(1));
t = sprintf('%d:%d:%d',c(4),c(5),c(6));

switch lower(dateOrTime)
    case 'date'
        str = d;
    case 'time'
        str = t;
    case 'both'
        str = [d, ' ' , t];
end


end