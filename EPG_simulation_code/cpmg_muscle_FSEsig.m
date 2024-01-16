function signal = cpmg_muscle_FSEsig(T2w, T2f, FF, B1, optw, optf)

sig_f = FSEsig(T2f,B1,1,optf);
sig_w = FSEsig(T2w, B1,1,optw);

signal = (1 - FF).*sig_w + FF.*sig_f;