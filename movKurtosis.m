function output = movKurtosis(data, k, extensionValue)
%MOVKURTOSIS computes kurtosis on a 2d matrix over a sliding window of dimension k
% input:
% - data: the data of which we want to calculage the kurtosis
% - k: the window size
% - extensionValue: corresponding to the value to specify in the
% 'Endpoints' option for other sliding window functions like movmean

sz=size(data);
output=zeros(sz);

%extension of the original data
data=padarray(data, [k,k], extensionValue);

for y=k+1:sz(1)+k
    for x=k+1:sz(2)+k
        landmark=data(y-k:y+k,x-k:x+k); %cause it's what you see through the window
        landmark=landmark(:);
        if min(landmark)==max(landmark)
            output(y-k,x-k)=0;
        else
            output(y-k,x-k)=kurtosis(landmark(:)); %compute excess kurtosis of the landmark
        end
    end
end
end

