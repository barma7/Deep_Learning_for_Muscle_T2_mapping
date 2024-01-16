function [dictionary, LUT_combined] = create_muscle_dictionary(optw, optf, T2ws, T2fs, B1s, FFs)

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

dictionary_water = zeros(Nb_exmpls,optw.etl);
dictionary_fat = zeros(Nb_exmpls,optw.etl);

T2ws = LUT(:,2);
T2fs = LUT(:,1);
B1s = LUT(:,3);

tic
parfor i = 1:Nb_exmpls
    [dictionary_water(i,:), dictionary_fat(i, :)] = cpmg_muscle_FSEsig2(T2ws(i), T2fs(i), B1s(i), optw, optf);
end
toc 

%% Combine dictionary for FFs
dictionary = zeros(Nb_exmpls*length(FFs),optw.etl);
LUT_combined = zeros(Nb_exmpls*length(FFs),4);

count = 1;
tic
for i = 1:length(FFs)
    for j = 1:Nb_exmpls
        dictionary(count,:) = (1 - FFs(i)).*dictionary_water(j,:) + FFs(i).*dictionary_fat(j,:);
        LUT_combined(count,:) = cat(2, FFs(i), LUT(j,:));
        count = count+1;
    end
end
toc