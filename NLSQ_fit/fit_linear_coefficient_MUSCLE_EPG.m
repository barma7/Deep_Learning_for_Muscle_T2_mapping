function [c,Rsq, SNR] = fit_linear_coefficient_MUSCLE_EPG(p, ydata, opt)

t2w = p(1);
b1 = p(2);

[Sw, Sf] = cpmg_muscle_FSEsig2(t2w, opt.T2fat, b1, opt.optw, opt.optf);

A = cat(2, Sw, Sf);

c = A\ydata;

% calculate R squared
yEst = A*c;

SStot = sum((ydata - mean(ydata)).^2);            % Total Sum-Of-Squares
SSres = sum((ydata - yEst).^2);                   % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;  
SNR = (sum(abs(ydata).^2)/length(ydata))/std(ydata - yEst);

