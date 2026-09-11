%% smooth a gridded variable in x-y

function vout = smoothxy(v,wx,wy)

[ny nx nz] = size(v);

% loop over z
for k = 1:nz;
    
    vtmp = v(:,:,k);
    
    % pad the edges of v
    vpad=cat(2,vtmp(:,(end-wx+1):end,:),vtmp,vtmp(:,1:wx,:)); 
    
    % do smoothing
    [vsm,jnk]=ndnanfilter(vpad,'hamming',[wy wx],[],[],[],1);
    
    % remove padded edges
    vsm=vsm(:,(wx+1):(end-wx),:); 
    vout(:,:,k) = vsm;
    
end
    