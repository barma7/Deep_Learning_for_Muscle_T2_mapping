function yEst = fit_FAT_EPG(p, ydata, opt)

t2w = p(1);
t2f = p(2);
b1 = p(3);

signal = cpmg_muscle_FSEsig(t2w, t2f, 0.9, b1, opt.optw, opt.optf);

A = signal;

c = A\ydata;

yEst = A*c;

%figure(2)
%plot(yEst, '--', "LineWidth", 1.5); hold on;

