function [u, v, sig_out, tau_out,tvt,gap] = tgv2solve_pd_3d_new(data,encop,u0, lambda, alpha0, alpha1, maxits, reduction,imgdims,dx,dy,dz)
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

factor=20;


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

% set zero output
sig_out=0; tau_out=0;tvt=0;gap=0;

% step-sizes
%----------------------------------------------------------
% estimate operator norm using power iteration
x1  = rand(M,N,L); 
% y1 = Kh(K(x1));
% for i=1:10
%     if norm(y1(:))~=0
%         x1 = y1./norm(y1(:));
%     else
%         x1 = y1;
%     end
%     [y1] = Kh(K(x1));
%     l1 = y1(:)'*x1(:);
% end
% Lip = max(abs(l1));                                     % Lipschitz constant estimate

sigma  = 1/3;%1/sqrt(1 + sqrt(12 + Lip ) + 12);   % dual step size
tau    = 1/3;%sigma;                                             % primal step size    

sig_out(1)=sigma;
tau_out(1)=tau;

% Algorithmic########################################

for k=0:maxits
    
    % reduce alpha0, alpha1 --> more focus on data-fidelity
    %alpha0 = exp(k/maxits*log(alpha01) + (maxits-k)/maxits*log(alpha00));
    %alpha1 = exp(k/maxits*log(alpha11) + (maxits-k)/maxits*log(alpha10));
    
    
    % Dual ascent step
    %===============================================
    
    % gradient
    %--------------------------------------------------------
    y(:,:,:,1:3) = y(:,:,:,1:3) + sigma*( fgrad_3(ext(:,:,:,1),dx,dy,dz) - ext(:,:,:,2:4) );
    
    % projection
    denom = max( 1,sqrt( ...
        abs(y(:,:,:,1)).^2 + abs(y(:,:,:,2)).^2 + abs(y(:,:,:,3)).^2 ...
                                        )/alpha1 );
    y(:,:,:,1) = y(:,:,:,1)./denom;  y(:,:,:,2) = y(:,:,:,2)./denom;  y(:,:,:,3) = y(:,:,:,3)./denom;
    
    % symmetrized gradient
    %--------------------------------------------------------
    y(:,:,:,4:9) = y(:,:,:,4:9) + sigma*(sym_bgrad_3(ext(:,:,:,2:4),dx,dy,dz));
    
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
                                 / (1+sigma*lambda);
    
    
    % Primal descent step
    %===============================================
    
    % divergence
    %-------------------------------------------------------
    ext(:,:,:,1) = x(:,:,:,1) - tau*( -  bdiv_3( y(:,:,:,1:3),dx,dy,dz) + reshape( Kh(z), [M,N,L]) );
    
    % symmetrized divergence
    %-------------------------------------------------------
    ext(:,:,:,2:4) = x(:,:,:,2:4) - tau*( - fdiv_3( y(:,:,:,4:9),dx,dy,dz ) - y(:,:,:,1:3) );

    
    %Set extragradient
    %===============================================
    x=2*ext - x;
    
    %Swap extragradient and primal variable
    [x,ext] = deal(ext,x);
    
    %  adapt stepsize
    if (k<10) || (rem(k,50) == 0)
        [sigma,tau] = steps_tgv2(x,sigma,tau,K,dx,dy,dz);
        sig_out(2+k)=sigma;
        tau_out(k+2)=tau;
        %display(['sig=',num2str(sigma),' | tau=',num2str(tau)]);
    end
    
    % calculate pd-gap
     if rem(k,factor) == 0
				    tvt(1 + k/factor) = get_tgv2(x,alpha0,alpha1,dx,dy,dz);
				    gap(1 + k/factor) = abs( tvt(1 + k/factor) +  gstar_tgv2(x,y,z,K, Kh, data, dx,dy,dz,1) );
				    gap(1 + k/factor) = gap(1 + k/factor)./(N*M*L);
	 end
        
    if mod(k,50) == 0
        %imwrite( abs(x(:,:,floor(end/2),1))./max(abs(x(:))) , ['./tgvcheck/tgv_iter',num2str(k),'.png'] );
        fprintf('TGV2-L2-2D-PD: it = %4d, alpha0 = %f, rmse = %f\n', k, alpha0, 0);
    end
    
end

u = x(:,:,:,1);
v = x(:,:,:,2:4);
