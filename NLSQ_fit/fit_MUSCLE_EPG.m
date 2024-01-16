function yEst = fit_MUSCLE_EPG(p, ydata, opt)

t2w = p(1);
b1 = p(2);

[Sw, Sf] = cpmg_muscle_FSEsig2(t2w, opt.T2fat, b1, opt.optw, opt.optf);

A = cat(2, Sw, Sf);

c = A\ydata;

yEst = A*c;

%figure(2)
%plot(yEst, '--', "LineWidth", 1.5); hold on;

