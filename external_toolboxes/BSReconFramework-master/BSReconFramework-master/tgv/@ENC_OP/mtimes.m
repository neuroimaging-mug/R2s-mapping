function ress = mtimes(a,bb)

   %% =========================================
    if a.adjoint  % Backward (adjoint) operation:  kspace to image
    % =========================================

        %---------------------------------------------------------------------------------------
        if a.is3d % 3d data
        %---------------------------------------------------------------------------------------

            if ~a.type   % cartesian data
            %---------------------------------------------------------------------------------------

                scale = numel(a.k)/sum(a.k(:) > 0);
                ress = zeros(size(a.csm));
                ress( repmat( a.k, [1 1 1 size(a.csm,4)] ) ==1 ) = bb(:);
                ress = ifft3c( ress );%.*sqrt(scale);
                ress = sum( conj(a.csm).*ress,4 );

            else    % non-cartesian data
             %---------------------------------------------------------------------------------------

                error('not implemented')

            end % cartesian / non-cartesian

        %---------------------------------------------------------------------------------------
        else    % 2d data
        %---------------------------------------------------------------------------------------

            if ~a.type     % cartesian data
            %---------------------------------------------------------------------------------------

                scale = numel(a.k)/sum(a.k(:) > 0);
                ress = zeros(size(a.csm));
                ress(repmat(a.k,[1 1 size(a.csm,3)] )==1) = bb(:);
                ress = transform_kspace_to_image(ress,[1,2]);%.*sqrt(scale);
                ress = sum(conj(a.csm).*ress,3);

            else  % non-cartesian data
            %---------------------------------------------------------------------------------------

                ress = zeros(size(a.csm,1),size(a.csm,2));
                bb = reshape(bb,[size(a.k,1),size(a.k,2),size(a.csm,3)]);
                for i=1:size(a.csm,3)
                    ress = ress + conj(a.csm(:,:,i)).*(a.nuFT'*(bb(:,:,i).*sqrt(a.w)));
                end

            end % cartesian / non-cartesian

        end % 3d / 2d


        ress = ress(:);

   %% =========================================
    else % Forward operation: image to kspace 
    % =========================================

        %---------------------------------------------------------------------------------------
        if a.is3d   % 3d data
        %---------------------------------------------------------------------------------------

            if ~a.type   % cartesian data
            %---------------------------------------------------------------------------------------

                scale = numel(a.k)/sum(a.k(:) > 0);            
                ress = repmat( reshape(bb,size(a.csm,1),size(a.csm,2),size(a.csm,3) ), [1 1 1 size(a.csm,4)] ) .* a.csm;
                ress = fft3c( ress );% .*sqrt(scale);
                ress = ress( repmat(a.k,[1 1 1 size(a.csm,4)]) == 1 );


            else % non-cartesian data
             %---------------------------------------------------------------------------------------

               error('not implemented')

            end % cartesian / non-cartesian

        %---------------------------------------------------------------------------------------
        else % 2d data
        %---------------------------------------------------------------------------------------

            if ~a.type   % cartesian data
            %---------------------------------------------------------------------------------------

                scale = numel(a.k)/sum(a.k(:) > 0);        
                ress = repmat( reshape(bb,size(a.csm,1),size(a.csm,2)),[1 1 size(a.csm,3)] ) .* a.csm;
                ress = transform_image_to_kspace(ress, [1,2]);%*sqrt(scale);
                ress = ress(repmat(a.k,[1 1 size(a.csm,3)]) == 1);


            else % non-cartesian data
             %---------------------------------------------------------------------------------------

                bb = reshape(bb,[size(a.csm,1),size(a.csm,2)]);
                ress = a.nuFT*(repmat(bb,[1 1 size(a.csm,3)]).*a.csm);
                ress = ress.*repmat(sqrt(a.w),[1 1 size(a.csm,3)]);

            end % cartesian / non-cartesian

        end % 3d / 2d

        ress = ress(:);

    end

end



