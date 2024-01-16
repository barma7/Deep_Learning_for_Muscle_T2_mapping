function [c,Rsq, SNR] = fit_linear_coefficient_FAT_EPG(p, ydata, opt)

t2w = p(1);
t2f = p(2);
b1 = p(3);

signal = cpmg_muscle_FSEsig(t2w, t2f, 0.9, b1, opt.optw, opt.optf);

A = signal;

c = A\ydata;

% calculate R squared
yEst = A*c;

SStot = sum((ydata - mean(ydata)).^2);            % Total Sum-Of-Squares
SSres = sum((ydata - yEst).^2);                   % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;  
SNR = (sum(abs(ydata).^2)/length(ydata))/std(ydata - yEst);

