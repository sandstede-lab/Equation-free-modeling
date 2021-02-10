%% restriction operator based on diffusion map
% newData   - point to get restriction values for (column vector)
% evals     - square matrix of eigenvalues
% evecs     - matrix with eigenvectors in columns
% origData  - matrix of data in columns used to compute eigenvalues/vectors
% eps       - length scale for kernel
%
% returns diffusion map embedding for new data using the Nystrom extension
% in a column vector
function pnew = diffMapRestrict(newData,evals,evecs,origData,eps)

dist = pdist2(newData',origData')';     % calculate the pairwise distances between newData and origData
w = basicKernel(dist);
k = (1./sum(w)).*w;
pnew = (evecs' * k)./diag(evals);

    function af = basicKernel(s)
        af = exp(-s.^2/eps^2);
    end

end