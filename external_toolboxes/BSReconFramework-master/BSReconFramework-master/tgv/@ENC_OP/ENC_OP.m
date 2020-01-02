function  res = ENC_OP( csm, k, w, type )

dims = size(csm);

if length(dims) == 4
    res.is3d =1;
else 
    res.is3d = 0;
end

%w = w*prod([size(csm,1),size(csm,2)]);
if type 
    %res.nuFT = NUFFT(k, w, 1, 0, [size(csm,1),size(csm,2)], 2);
    res.nuFT = NUFFT(k, 1, 1, 0, [size(csm,1),size(csm,2)], 2);

end

res.type = type;
res.k = k;
res.w = w;
res.csm = csm;
res.adjoint = 0;

res = class(res,'ENC_OP');

