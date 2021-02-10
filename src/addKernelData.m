function addKernelData(origData, eps)
%ADDKERNELDATA Stores self-distance kernel data to avoid having
%to compute it every time
selfDist = pdist2(origData',origData')';
selfKernel = basicKernel(selfDist, eps);
sumsSelfKernel = sum(selfKernel);
save('sumsSelfKernel.mat', 'sumsSelfKernel');
end

function af = basicKernel(s, eps)
    af = exp(-s.^2/eps^2);
end