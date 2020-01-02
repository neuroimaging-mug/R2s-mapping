function [u,e] = tgv2solve_pd_3d(data,encop,u0, alpha0, alpha1, maxits, reduction, inner_iter,imgdims)
% Primal dual TGV2 algorithm, as described in the TGV paper
%
% (c) 31.8.2010
% Florian Knoll (florian.knoll@tugraz.at)
% Kristian Bredies (kristian.bredies@uni-graz.at
% Thomas Pock (pock@icg.tugraz.at)
% Rudolf Stollberger (rudolf.stollberger@tugraz.at)
%
% If you consider this code to be useful for your research, please cite:
% Knoll, F.; Bredies, K.; Pock, T.; Stollberger, R.: Second Order Total
% Generalized Variation (TGV) for MRI: Magnetic Resonance in Medicine,
% to appear (2010)
%
% last change 01.10.2015: Matthias Schloegl
% -------------------------------------------------------------------------------------------------


% Initialize########################################

K = @(x) encop*x;
Kh = @(x) encop'*x;

% normalize data
% dscale = 100/norm(abs(data(:)));
% data = data * dscale;

alpha00 = alpha0;
alpha10 = alpha1;
alpha01 = alpha0*reduction;
alpha11 = alpha1*reduction;

% Initailize
%===============================================
M = imgdims(1); N = imgdims(2); L = imgdims(3); P = imgdims(4);

% primal variable
x = zeros(M,N,L,1+3);

% intial guess
%----------------------------------------------------------
if isempty(u0)
    x(:,:,:,1) = reshape( Kh(data(:)) , [M,N,L]);
    %x(:,:,:,1) = max(0, real( x(:,:,:,1)));
else
    x(:,:,:,1) = reshape(u0 , [M,N,L]); clear u0;
end

% extragradient
ext = x;

% dual variable
y = zeros(M,N,L,3+6);
z = zeros(size(data,1),1);

e = [];

% step-sizes
%----------------------------------------------------------
% estimate operator norm using power iteration
x1  = rand(M,N,L); 
y1 = Kh(K(x1));
for i=1:10
    if norm(y1(:))~=0
        x1 = y1./norm(y1(:));
    else
        x1 = y1;
    end
    [y1] = Kh(K(x1));
    l1 = y1(:)'*x1(:);
end
Lip = max(abs(l1));                                     % Lipschitz constant estimate

sigma = 1/sqrt(1 + sqrt(12 + Lip ) + 12);   % dual step size
tau    = sigma;                                             % primal step size    


% Algorithmic########################################

for k=0:maxits
    
    % reduce alpha0, alpha1 --> more focus on data-fidelity
    alpha0 = exp(k/maxits*log(alpha01) + (maxits-k)/maxits*log(alpha00));
    alpha1 = exp(k/maxits*log(alpha11) + (maxits-k)/maxits*log(alpha10));
    
    
    % Dual ascent step
    %===============================================
    
    % gradient
    %--------------------------------------------------------
    y(:,:,:,1:3) = y(:,:,:,1:3) + sigma*( cat(4, dxp3( ext(:,:,:,1) ), dyp3( ext(:,:,:,1) ), dzp3( ext(:,:,:,1) ) ) ...
        - ext(:,:,:,2:4) );
    
    % projection
    denom = max( 1,sqrt( ...
        abs(y(:,:,:,1)).^2 + abs(y(:,:,:,2)).^2 + abs(y(:,:,:,3)).^2 ...
                                        )/alpha1 );
    y(:,:,:,1) = y(:,:,:,1)./denom;  y(:,:,:,2) = y(:,:,:,2)./denom;  y(:,:,:,3) = y(:,:,:,3)./denom;
    
    % symmetrized gradient
    %--------------------------------------------------------
    y(:,:,:,4) = y(:,:,:,4) + sigma*dxm3( ext(:,:,:,2) ); % xx
    y(:,:,:,5) = y(:,:,:,5) + sigma*dym3( ext(:,:,:,3) ); % yy
    y(:,:,:,6) = y(:,:,:,6) + sigma*dzm3( ext(:,:,:,4) ); % zz
    
    y(:,:,:,7) = y(:,:,:,7) + sigma*( (dym3( ext(:,:,:,2) ) + dxm3( ext(:,:,:,3) ) )/2 );  % xy: (dy(x) + dx(y) )/2
    y(:,:,:,8) = y(:,:,:,8) + sigma*( (dym3( ext(:,:,:,4) ) + dzm3( ext(:,:,:,3) ) )/2 );  % yz: (dy(z) + dz(y) )/2
    y(:,:,:,9) = y(:,:,:,9) + sigma*( (dzm3( ext(:,:,:,2) ) + dxm3( ext(:,:,:,4) ) )/2 );  % xz: (dz(x) + dx(z) )/2
    
    % projection
    denom = max( 1,  sqrt(  ...
        abs( y(:,:,:,4) ).^2 + abs( y(:,:,:,5) ).^2  + abs( y(:,:,:,6) ).^2 ....
        + 2*abs( y(:,:,:,7) ).^2 + 2*abs( y(:,:,:,8) ).^2 + 2*abs( y(:,:,:,9) ).^2 ...
                                        ) /alpha0 );
    y(:,:,:,1)=y(:,:,:,1)./denom; y(:,:,:,2)=y(:,:,:,2)./denom; y(:,:,:,3)=y(:,:,:,3)./denom;
    y(:,:,:,4)=y(:,:,:,4)./denom; y(:,:,:,5)=y(:,:,:,5)./denom; y(:,:,:,6)=y(:,:,:,6)./denom;
    
    % operator
    %--------------------------------------------------------
    z = ( z + sigma*( K( ext(:,:,:,1) ) - data(:) ) ) ...
                                 / (1+sigma);
    
    % Primal descent step
    %===============================================
    
    % divergence
    %-------------------------------------------------------
    ext(:,:,:,1) = x(:,:,:,1) - tau*( - ( dxp3_( y(:,:,:,1) ) + dyp3_( y(:,:,:,2) ) + dzp3_( y(:,:,:,3) ) ) ...
                                                + reshape( Kh(z), [M,N,L]) );
    
    % symmetrized divergence
    %-------------------------------------------------------
    ext(:,:,:,2) = x(:,:,:,2) - tau*( - ( dxm3_( y(:,:,:,4) ) + dym3_( y(:,:,:,7) ) + dzm3_( y(:,:,:,9) ) ) ...
                                                - y(:,:,:,1) );
    ext(:,:,:,3) = x(:,:,:,3) - tau*( - ( dxm3_( y(:,:,:,7) ) + dym3_( y(:,:,:,5) ) + dzm3_( y(:,:,:,8) ) ) ...
                                                - y(:,:,:,2) );
    ext(:,:,:,4) = x(:,:,:,4) - tau*( - ( dxm3_( y(:,:,:,9) ) + dym3_( y(:,:,:,8) ) + dzm3_( y(:,:,:,6) ) ) ...
                                                - y(:,:,:,3) );

    
    %Set extragradient
    %===============================================
    x=2*ext - x;
    
    %Swap extragradient and primal variable
    [x,ext] = deal(ext,x);
    
    % TODO: adapt stepsize
    
    if mod(k,10) == 0

%         primal_dual_energy = 0;
%         dual_energy = sum(sum(abs(x(:,:,:,1)).^2))/2.0;
%         e = [e dual_energy];
%         temp = x(:,:,floor(end/2),1);
%         imwrite( abs(x(:,:,floor(end/2),1))./max(abs(temp(:))) , ['./tgvcheck/tgv_iter',num2str(k),'.png'] );
        fprintf('TGV2-L2-2D-PD: it = %4d, alpha0 = %f, rmse = %f\n', k, alpha0, 0);
    end
    
end

u = x(:,:,:,1);
