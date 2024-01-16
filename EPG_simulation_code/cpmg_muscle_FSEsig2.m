function [sig_w, sig_f] = cpmg_muscle_FSEsig2(T2w, T2f, B1, optw, optf)

sig_f = FSEsig(T2f,B1,1,optf);
sig_w = FSEsig(T2w, B1,1,optw);

