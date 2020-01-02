function [pattern, R_True1] = getVdrandomPattern(dims, AF, p, numSlice)
    [pdf1, ~] = genPDF([dims(1),numSlice],p,1/AF,2,0,0);
    [pattern, R_True1] = gen_random_pattern(pdf1,AF,100000000);
    pattern = reshape(pattern, [dims(1),1,numSlice]);
    pattern = repmat(pattern, [1, dims(2), 1]);
    sliceDiff = (dims(3) - numSlice)/2;
    pattern = cat(3,zeros(dims(1),dims(2),sliceDiff), cat(3, pattern, zeros(dims(1),dims(2),sliceDiff)));
    pattern = repmat(pattern, [1,1,1,dims(4),dims(5)]);