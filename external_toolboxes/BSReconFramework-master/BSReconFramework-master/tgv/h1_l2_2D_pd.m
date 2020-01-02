function [u e] = h1_l2_2D_pd(g, f, K, w, alpha0, maxits, reduction)
% H1 regularized iterative coil sensitivity estimation
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

check_it = 1;

alpha00 = alpha0;
alpha01 = alpha0*reduction;

[M N] = size(g);

v = zeros(size(f));
xi = zeros(M,N,2);

u = zeros(M,N);
u_ = u;

e = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiindices of the spatial derivatives

% derivatives
% | uxx uxy |
% | uyx uyy |

% multiindices
% | 1 3 |
% | 3 2 |

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = sqrt(8);
tau_p = 1/8;
tau_d = 1/8;

for k=0:maxits
    
  alpha0 = exp(k/maxits*log(alpha01) + (maxits-k)/maxits*log(alpha00));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SAVE VARIABLES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  uold = u;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DUAL UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Ku0 = K*u_;
  Ku = sqrt(w).*Ku0;  
  r =  Ku - f;  
  
  v = (v + tau_d*r)/(1+ tau_d);
  
  ux = dxp(u_);
  uy = dyp(u_);
  
  xi(:,:,1) = (xi(:,:,1) - tau_d*ux)/(1+tau_d/alpha0);
  xi(:,:,2) = (xi(:,:,2) - tau_d*uy)/(1+tau_d/alpha0);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PRIMAL UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ww = K'*(sqrt(w).*v);
  div = dxm(xi(:,:,1)) + dym(xi(:,:,2));
  
  u = u - tau_p*(ww + div);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % AUXILIARY UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  u_ = 2*u - uold;
  
  if mod(k,check_it) == 0
    div = dxm(xi(:,:,1)) + dym(xi(:,:,2));
    primal_dual_energy = -sum(sum(real(u.*div))) + 0.5*sum(sum(abs(Ku-f).^2));
    dual_energy = sum(sum(abs(u).^2));
    e = [e dual_energy];
    rmse = sqrt(mean(mean(abs(u-g).^2)));
    fprintf('H1-L2-2D-PD: it = %4d,  e = %f, rmse = %f\n', k, e(end), rmse);    
    imshow(abs(u),[]);
    drawnow;
  end
  
end
