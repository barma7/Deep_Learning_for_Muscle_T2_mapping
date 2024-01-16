function [FFfit, T2fit, B1fit, tElapsed] = dictionary_matching_dot(LUT, Dict, Data)

% Function to estimate FF, T2w and B1 from MESE data with dictonary
% matching with EPG simulamted data
% ----------------  INPUT PARAMETERS -------------------------------------
% LUT........Array with Look-Up-Table (FF, T2w, B1) [k x 3]
% Dict.......Array with noramlized simulated signals [k x n]
% Data...... Array with noramlized signals [m x n]

%------------------ OUTPUT PARAMETERS -------------------------------------
% FFfit.........Array with fitted FF [m x 1]
% T2fit.........Array with fitted T2 [m x 1]
% B1fit.........Array with fitted B1 [m x 1]

% k = # dictionary entries, m = #test_examples (pixels), n = EchoTrain Length


%% PART I: USING THE COSINE DISTANCE TO PERFORM MATCHING
Nb_tests = size(Data,1);

Divisore = 250;
% break the Nb_tests in block of 10000 test data
if (Nb_tests/Divisore - round(Nb_tests/Divisore)) < 0.5 && (Nb_tests/Divisore - round(Nb_tests/Divisore))>0
    Nb_blocks = round(Nb_tests/Divisore)+1;
else
    Nb_blocks = round(Nb_tests/Divisore);
end
%Nb_blocks

FFfit = zeros(Nb_tests, 1);
T2fit = zeros(Nb_tests, 1);
B1fit = zeros(Nb_tests, 1);

tstart = tic;

indx = 1;
stop_indx = Divisore;
start_indx = 1;

for b = 1:Nb_blocks
    if (b == Nb_blocks && stop_indx>Nb_tests)
        stop_indx = Nb_tests-(stop_indx-Divisore);
        dotMatrix = Data(start_indx:start_indx+stop_indx-1,:)*Dict';    
        [maxDot_array, Mtchd_indx_array] = max(abs(dotMatrix),[],2);

        for j=1:length(maxDot_array) 
            Mtch_indx = Mtchd_indx_array(j);
            FFfit(indx) = LUT(Mtch_indx,1);
            T2fit(indx) = LUT(Mtch_indx,2);
            B1fit(indx) = LUT(Mtch_indx,3);

            indx = indx+1;
        end
    else
        dotMatrix = Data(start_indx:stop_indx,:)*Dict';    
        [maxDot_array, Mtchd_indx_array] = max(abs(dotMatrix),[],2);

        for j=1:length(maxDot_array) 
            Mtch_indx = Mtchd_indx_array(j);
            FFfit(indx) = LUT(Mtch_indx,1);
            T2fit(indx) = LUT(Mtch_indx,2);
            B1fit(indx) = LUT(Mtch_indx,3);

            indx = indx+1;
        end
        start_indx = start_indx+Divisore;
        stop_indx = stop_indx+Divisore;
    end
end
tElapsed = toc(tstart);
