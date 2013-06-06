function [slope,intercept,slerr,interr] = WeightedLSQFit(x,y,w)
%WeightedLinearLeastSquaresFit(x,y,w)
%   Take in arrays representing (x,y) values 
%   for a set of linearly varying data and 
%   an array of weights (w). Perform a 
%   weighted linear least squares regression. 
%   Return the resulting slope and intercept 
%   parameters of the best fit line with 
%   their uncertainties.
slope = (sum(w)*sum(w.*x.*y)-sum(w.*x)*sum(w.*y))/(sum(w)*sum(w.*x.^2)-(sum(w.*x))^2);
intercept = (sum(w.*x.^2)*sum(w.*y)-sum(w.*x)*sum(w.*x.*y))/(sum(w)*sum(w.*x.^2)-(sum(w.*x))^2);
slerr = sqrt(sum(w)/(sum(w)*sum(w.*x.^2)-(sum(w.*x)).^2));
interr = sqrt(sum(w.*x.^2)/(sum(w)*sum(w.*x.^2)-(sum(w.*x))^2));

end

