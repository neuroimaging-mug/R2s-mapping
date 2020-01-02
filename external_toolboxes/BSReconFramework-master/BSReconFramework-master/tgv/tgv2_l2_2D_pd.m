function [u e] = tgv2_l2_2D_pd(sens, f, K, w, alpha0, alpha1, maxits, reduction, inner_iter)
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

check_it = 1;

alpha00 = alpha0;
alpha10 = alpha1;
alpha01 = alpha0*reduction;
alpha11 = alpha1*reduction;

[M N P] = size(sens);

v = zeros(size(f));
p = zeros(M,N,2);
q = zeros(M,N,3);

u = zeros(M,N);
for ii=1:size(f,3)
      u = u + conj(sens(:,:,ii)).*(K'*(sqrt(w(:,:,ii)).*f(:,:,ii)));
end
u = max(0, real(u));

xi = zeros(M,N,2);
u_ = u;
xi_ = xi;

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

L = sqrt(64);
tau_p = 1/16;
tau_d = 1/8;

for k=0:maxits
  alpha0 = exp(k/maxits*log(alpha01) + (maxits-k)/maxits*log(alpha00));
  alpha1 = exp(k/maxits*log(alpha11) + (maxits-k)/maxits*log(alpha10));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SAVE VARIABLES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  uold = u; xiold = xi;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DUAL UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % operator
  Ku0 = zeros(size(f));
  for ii=1:P
     Ku0(:,:,ii) = K*(u_.*sens(:,:,ii));
  end
  Ku_ = sqrt(w).*Ku0;  
  r =  Ku_ - f;
  
  v = (v + tau_d*r)/(1+tau_d);
  
  % gradient
  ux = dxp(u_);
  uy = dyp(u_);
  
  p(:,:,1) = p(:,:,1) - tau_d*(ux + xi_(:,:,1));
  p(:,:,2) = p(:,:,2) - tau_d*(uy + xi_(:,:,2));
  
  % projection
  absp = sqrt(abs(p(:,:,1)).^2 + abs(p(:,:,2)).^2);
  denom = max(1,absp/alpha1);
  p(:,:,1) = p(:,:,1)./denom;
  p(:,:,2) = p(:,:,2)./denom;  
  
  % symmetrized gradient
  gradxi1 = dxm(xi_(:,:,1));
  gradxi2 = dym(xi_(:,:,2));
  gradxi3 = (dym(xi_(:,:,1)) + dxm(xi_(:,:,2)))/2;
  
  q(:,:,1) = q(:,:,1) - tau_d*gradxi1;
  q(:,:,2) = q(:,:,2) - tau_d*gradxi2;
  q(:,:,3) = q(:,:,3) - tau_d*gradxi3;
  
  % projection
  absq = sqrt(abs(q(:,:,1)).^2 + abs(q(:,:,2)).^2 + 2*abs(q(:,:,3)).^2);
  denom = max(1,absq/alpha0);
  q(:,:,1) = q(:,:,1)./denom;
  q(:,:,2) = q(:,:,2)./denom;
  q(:,:,3) = q(:,:,3)./denom;  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PRIMAL UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  % dual operator
  ww = zeros(size(u));
  for ii=1:size(f,3)
      ww = ww + conj(sens(:,:,ii)).*(K'*(sqrt(w(:,:,ii)).*v(:,:,ii)));
  end
  
  % divergence
  divp = dxm(p(:,:,1)) + dym(p(:,:,2));
  
  u = u - tau_p*(ww + divp);
  
  % divergence
  divq1 = dxp(q(:,:,1)) + dyp(q(:,:,3));
  divq2 = dxp(q(:,:,3)) + dyp(q(:,:,2));
  
  xi(:,:,1) = xi(:,:,1) - tau_p*(divq1 - p(:,:,1));
  xi(:,:,2) = xi(:,:,2) - tau_p*(divq2 - p(:,:,2));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % AUXILIARY UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  u_ = 2*u - uold;
  xi_ = 2*xi - xiold;
  
  if mod(k,check_it) == 0
    primal_dual_energy = 0;
    dual_energy = sum(sum(abs(u).^2))/2.0;
    e = [e dual_energy];
    fprintf('TGV2-L2-2D-PD: it = %4d, alpha0 = %f, rmse = %f\n', k, alpha0, 0);    
    imshow(abs(u),[]);
    drawnow;
  end
  
end
