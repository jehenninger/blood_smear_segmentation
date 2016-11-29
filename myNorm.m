function [ normChannel ] = myNorm(channel,max_threshold,min_threshold)

%   Normalizes mean intensity measurements by subtracting the minimum and
%   dividing by the maximum. Min_threshold is optional argument specifying what
%   value to normalize to
if size(channel,2)>1
    if ~exist('max_threshold','var')
        max_threshold = max(channel);
    end
    
    if ~exist('min_threshold','var')
        min_threshold = min(channel);
    end
    
    normChannel = channel - min_threshold;
    normChannel = normChannel./(max_threshold-min_threshold);
    
else
    
    if ~exist('max_threshold','var')
        max_threshold = max(channel(:));
    end
    
    if ~exist('min_threshold','var')
        min_threshold = min(channel(:));
    end
    
    
    normChannel = channel - min_threshold;
    normChannel = normChannel./(max_threshold-min_threshold);
    
end



end

