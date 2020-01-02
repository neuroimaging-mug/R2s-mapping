index = [1:5];
for i = 4:5

switch par.traj
    case 'cart'
        encop2 = ENC_OP(repmat(TGV2_recon(:,:,:,index(i)),[1 1 1 ncoils]).*smaps,pattern(:,:,:,2),1,0); 
    case 'radial'
        encop2 = ENC_OP(repmat(TGV2_recon,[1 1 ncoils]).*smaps,pattern,w,1);
end

% BSrecon = h1_solve_cg(data,encop,mu,maxit,tol)
datain = (squeeze(par.y(:,:,:,2,:)));
datain = datain(repmat(pattern(:,:,:,2),[1 1 1 ncoils]));

mu = 1e-4;
phaseImage1(:,:,:,i) = h1_solve_cg(datain,encop2,mu,500,1e-6,[ny,nx,nz,ncoils], par.dz);

phaseImage(:,:,:,i) = angle(phaseImage1(:,:,:,i));
B1Map(:,:,:,i) = sqrt(phaseImage(:,:,:,i) / par.Kbs);
flipAngleMap(:,:,:,i) = B1Map(:,:,:,i) * par.gamma*sum(par.pulse)*par.deltat/max(par.pulse) *180/pi / par.alphaBS *100;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp =12:16;
TGV2_recon = zeros(ny,nx,nz,length(temp));
for i = 1:length(temp)
    alpha_temp = alpha*2^(temp(i));
    TGV2_recon(:,:,:,i) = tgv2solve_pd_3d_new(datain,encop, lowres(:), sqrt(3)*alpha_temp, alpha_temp, maxit, reduction, innerIter,[ny,nx,nz,ncoils], 1,1,par.dz);
    disp(num2str(i))
end
