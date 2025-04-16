function y = invertedGaussian(x,params)
%INVERTEDGAUSSIAN used for kurtosis' mf
    y=1-gaussmf(x, params);
end

