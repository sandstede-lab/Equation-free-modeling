%% aligns the maximum point of a quadratic fit to data to the 
% position specified (middle is the default position)
% input:
%   data - the data points to align
%   position - the position to align the max data point to (optional)
% Returns:
%   aligned - the data aligned with the max of the interpolation at position
function aligned = alignMax(data, position)
    n = length(data);                   

    aligned = shiftMax(data, position);     % shift the data initally to avoid boundary issues
    [max3indx, max3val] = getTop3(aligned); 	% find the maximum 3 data points and indices
    maxfun = polyfit(max3indx,max3val,2);          	% fit a quadratic to the max 3 data points
    actualidx = -maxfun(2)/(2*maxfun(1));         	% find the index of the max of the wave (-b/2a)
    
    fnCars = csape(1:1:n, aligned, 'periodic');            % interpolate the wave profile
    newpts = mod(linspace(1,n,n) + actualidx,n);    % find the shift needed to move the max to the end
    aligned = fnval(fnCars,newpts)';                    % move the max to the end
    
    % if specified, align to position
    if(nargin > 1)
        aligned = shiftMax(aligned, position);
    end
    
    %% get the values and indices surrounding the maximum of each profile
    % hways - the profile
    %
    % returns:
    %       maxinds - the 3 indices around/including the maximum
    %       maxvals - the 3 values around/including the maximum
    function [maxinds, maxvals] = getTop3(hways)
        [maxhw,maxind] = max(hways,[],1);                   % find index and value of maximum
        ind1 = mod(maxind-2, n)+1;                     % find point to the left
        ind3 = mod(maxind, n)+1;                       % point to the right
        maxinds = [ind1 ; maxind ; ind3];                   % all of the indices together
        maxvals = [hways(ind1) ; maxhw ; hways(ind3)];      % values at the indices
    end
    
    %% function to circshift max to beginning
    function c = shiftMax(hways, center)
        if(nargin>1)
            shift = center + 1;
        else
            shift = 1;
        end
        [~, maxH] = max(hways,[],1);  % locate the max headway for each data point
        c = zeros(size(hways));
        
        % align all of the headways with the max in the front
        for iCar = 1:length(maxH)
            c(:,iCar) = circshift(hways(:,iCar), [-maxH(iCar)+shift,0]);
        end
    end

 end
 