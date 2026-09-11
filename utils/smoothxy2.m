%% smooth a gridded variable in x-y

function vout = smoothxy2(v,wx,wy)

[ny nx nz] = size(v);

% pad the edges of v
vpad=cat(2,v(:,(end-wx+1):end),v,v(:,1:wx));
vpad=cat(1,vpad((end-wy+1):end,:),vpad,vpad(1:wy,:));

% do smoothing
[vsm,jnk]=ndnanfilter(vpad,'hamming',[wy wx],[],[],[],1);

% remove padded edges
vout=vsm((wy+1):(end-wy),(wx+1):(end-wx),:);

