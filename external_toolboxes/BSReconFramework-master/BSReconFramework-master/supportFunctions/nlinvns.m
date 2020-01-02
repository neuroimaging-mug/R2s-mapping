% Written and invented
% by Martin Uecker <muecker@gwdg.de> in 2008-09-22
%
% Modifications by Tilman Sumpf 2012 <tsumpf@gwdg.de>:
%	- removed fftshift during reconstruction (ns = "no shift")
%	- added switch to return coil profiles
%	- added switch to force the image estimate to be real
%	- use of vectorized operation rather than "for" loops
%
% Version 0.1
%
% Biomedizinische NMR Forschungs GmbH am
% Max-Planck-Institut fuer biophysikalische Chemie




function R = nlinvns(Y, n, returnProfiles, realConstr)
% usage: R = nlinvns(Y, n, returnProfiles, realConstr)
%

if nargin < 4
    realConstr = false;
    if (nargin < 3)
        returnProfiles = 0;
    end
end


disp('Start...')

alpha = 1.;

[x, y, c] = size(Y);

if returnProfiles
    R = zeros(x, y, n, c+2);
else
    R = zeros(x, y, n,2);
end


% initialization x-vector

X0 = zeros(x, y, c + 1);
X0(:,:,1) = 1.;	%object part

% initialize mask and weights
P = ones(size(Y(:,:,1)));
P(Y(:,:,1) ==  0) = 0;

W = weights(x, y);

P = fftshift2(P);
W = fftshift2(W);
Y = fftshift2(Y);

% normalize data vector
yscale = 100. / sqrt(scal(Y,Y));
YS = Y * yscale;

XT = zeros(x, y, c + 1);
XN = X0;
fprintf('Residuum ');

tic;
for i = 1:n,  % n is the number of iterations
    
    % the application of the weights matrix to XN
    % is moved out of the operator and the derivative
    
    XT(:, :, 1) = XN(:, :, 1);
    
    XT(:,:,2:end) = apweightsns( W,XN(:,:,2:end));
    
    RES = YS - opns(P, XT);
    
    fprintf(' %2.1f,  ',(sqrt(scal(RES, RES))));
    
    
    % calculate rhs
    
    r = derHns(P, W, XT, RES, realConstr);
    r = r + alpha * (X0 - XN);
    
    
    
    % begin CG
    
    z = zeros(size(r));
    
    d = r;
    dnew = scal(r, r);
    dnot = dnew;
    
    
    for j= 1 : 500
        
        % regularized normal equations
        
        q = derHns(P, W, XT, derns(P, W, XT, d), realConstr) + alpha * d;
        
        
        a = dnew / real(scal(d, q));
        
        z = z + a * d;
        r = r - a * q;
        
        dold = dnew;
        dnew = scal(r, r);
        
        d = d * (dnew / dold) + r;
        
        
        if (sqrt(dnew) < 1.e-2 * dnot)
            break
        end
        
    end
    fprintf('(%d)',j);
    % end CG
    
    XN = XN + z;
    
    alpha = alpha / 3.;
    
    % end
    
    % postprocessing...
    
    CR = apweightsns(W,XN(:,:,2:end));
    if returnProfiles
        R(:,:,i,3:end) = CR / yscale;
    end
    C = sum(conj(CR) .* CR,3);
    R(:, :, i, 1) = XN(:, : , 1) .* sqrt(C) / yscale;
    R(:, :, i, 2) = XN(:, : , 1);
    
    
end
R = fftshift2(R);
fprintf('\ndone in %s\n',humanize.seconds(toc));



end

function v = scal(a, b)
v = sum(sum(sum(conj(a) .* b))); % Bienchen sum herum
end


function C = apweightsns(W, CT)
C = nsIfft(ftimes(W, CT));
end

function C = apweightsnsH(W, CT)
C = ftimes(conj(W), nsFft(CT));
end




function K = opns(P, X)

K = ftimes(P, nsFft(ftimes(  X(:,:,1), X(:,:,2:end) ) ) );

end




function K = derns(P, W, X0, DX)

K = ftimes(X0(:,:,1),apweightsns(W,DX(:,:,2:end)));
K = K + ftimes(DX(:,:,1),X0(:,:,2:end));
K = ftimes(P,nsFft(K));


end

function DX = derHns(P, W, X0, DK, realConstr)

K = nsIfft(ftimes(P,DK));

if realConstr
    DXrho = sum( real( K .* conj(X0(:,:,2:end)))  ,3);
else
    DXrho = sum( K .* conj(X0(:,:,2:end))  ,3);
end
DXc = apweightsnsH(W,ftimes(K,conj(X0(:,:,1))));

DX = cat(3,DXrho, DXc);


end



function K = nsFft(M)
si = size(M);
a = 1 / (sqrt(si(1)) * sqrt(si(2)));
K = fft2(M) * a;

end

function K = nsIfft(M)
si = size(M);
a = sqrt(si(1)) * sqrt(si(2));
K = ifft2(M) * a;

end


function W = weights(x,y)
	W = zeros(x,y);	
	for i = 1:x
		for j = 1:y
			d = ((i - 1) / x - 0.5)^2 + ((j - 1) / y - 0.5)^2;
			W(i, j) = 1. / (1. + 220. * d)^16;
		end
    end
end


