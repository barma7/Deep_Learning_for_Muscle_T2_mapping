function [dictionary, LUT] = create_fat_dictionary(optw, optf, T2ws, T2fs, B1s, FF)

Nb_exmpls = length(T2ws)*length(T2fs)*length(B1s);

LUT = zeros(Nb_exmpls,3);
count = 1;
for i=1:length(T2fs)
    for k=1:length(T2ws)
        for l = 1: length(B1s)
            LUT(count,:) = cat(2,T2fs(i), T2ws(k), B1s(l));
            count = count+1;
        end
    end
end

dictionary = zeros(Nb_exmpls,optw.etl);

T2ws = LUT(:,2);
T2fs = LUT(:,1);
B1s = LUT(:,3);

tic
parfor i = 1:Nb_exmpls
    dictionary(i,:) = cpmg_muscle_FSEsig(T2ws(i), T2fs(i), FF, B1s(i), optw, optf);
end
toc