%% two-dimensional lifting operator
% newVal    - 1x2 vector representing the desired Leslie and Ann values
% evec      - eigenvectors from diffusion map in columns
% eval      - eigenvalues from diffusion map
% eps       - epsilon from diffusion map
% v0        - v0 parameter value
% oldData   - data from which the diffusion map was constructed
%
% returns a 120 x 1 vector of car positions and velocities that restricts
% to close to newVal
function [lifted, idxMin] = diffMapLift(newVal, evec, eval, eps,v0, oldData, h)
lsqOptions = optimset('Display','Off', 'TolX', 1e-12);
lsqOptions.FunctionTolerance = 1e-8;
lsqOptions.StepTolerance = 1e-9;
lsqOptions.OptimalityTolerance = 1e-8;

%% find closest datapoints
validPoints = 10;
newDist = pdist2(newVal',evec);
[~, index] = sort(newDist, 2, 'ascend');
idxMin = index(1:validPoints);
closestPoints = oldData(:, idxMin); % use closest profiles for linear combination
idxMin = sort(idxMin);

% solve for best coefficients for linear combination
liftedCoeffs = lsqnonlin(@liftingEquations, 1/validPoints*ones(validPoints, 1), zeros(validPoints, 1), ones(validPoints,1), lsqOptions);
liftedHeadway = closestPoints * liftedCoeffs; % create new lifted profile
lifted = [hwayToPos(liftedHeadway) ; optimalVelocity(h, liftedHeadway, v0)];

%% sets up and solves a system of equations to find the best linear combination

    function out = liftingEquations(coeffs)
        liftGuess = closestPoints * coeffs;
        restrictGuess = diffMapRestrict(liftGuess,eval,evec,oldData,eps);
        out = [restrictGuess - newVal  ; sum(coeffs) - 1; zeros(validPoints - 1, 1)]; 
    end

end