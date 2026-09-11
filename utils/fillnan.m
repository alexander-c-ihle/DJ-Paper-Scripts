function A = fillnan(M)
    
    A = 0*M;
    [ny,nx,nz] = size(M);
    for k = 1:nz
        tmp = inpaint_nans(M(:,:,k),5);
        A(:,:,k) =tmp;
    end
    