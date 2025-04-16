function [xf,yf] = wholeArithmeticRecombination(x,y,alpha)
    xf=alpha*x+(1-alpha)*y;
    yf=alpha*y+(1-alpha)*x;
end

