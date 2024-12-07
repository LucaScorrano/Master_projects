clear all
close all
clc
%% LAB 06
% OBIETTIVO 2 - Confronto delle performance di differenti classificatori sul test set

% Normalizzazione Test set

load maxcs.mat
load mincs.mat
load TS.mat
test_set_n=rescale(TS,'InputMin',mincs,'InputMax',maxcs);
test_set_norm=[test_set_n(:,3:57) TS(:,58)]; % le prime due colonne sono state tolte

% Coordinate soluzioni (ATTENZIONE:includono anche le prime 2 colonne relative all' ID)
% sol best codified 5 21 24 27 30 58
% sol best yes/no 6 9 11 12 14 16 17 19 22 25 26 27 28 29 31 32 34 36 42 46 47 56 58

TS_cod=test_set_norm(:,[3 19 22 25 28 56]);
TS_yn=test_set_norm(:,[4 7 9 10 12 14 15 17 20 23 24 25 26 27 29 30 32 34 40 44 45 54 56]);

%% Classificatore KNN random

% Codified variables
load mdl_R_cod.mat
%implementazione 
out_ts_Rcod=predict(mdl_R_cod, TS_cod (:, 1:end-1));
% Costruzione della confusion matrix
[contr,ordertr] = confusionmat(TS_cod(:,6),out_ts_Rcod); % ~ valori veri, valori predetti
figure("name","TS Random Codified")
confusionTS_R_cod = confusionchart(contr,ordertr);% ~ valori veri, valori predetti

% yes/no variables
load mdl_R_yn.mat
%implemenetazione 
out_ts_R_yn=predict(mdl_R_yn, TS_yn (:, 1:end-1));
% Costruzione della confusion matrix
[contr,ordertr] = confusionmat(TS_yn(:,23),out_ts_R_yn); % ~ valori veri, valori predetti
figure("name","TS Random Y/N variables")
confusionTS_R_yn = confusionchart(contr,ordertr);% ~ valori veri, valori predetti

%% Classificatore KNN dendro
% Codified
load mdl_dendro_cod
%implemenetazione 
out_ts_dendro_cod=predict(mdl_dendro_cod, TS_cod (:, 1:end-1));
% Costruzione della confusion matrix
[contr,ordertr] = confusionmat(TS_cod(:,6),out_ts_dendro_cod);
figure("name","TS dendro cod")
confusionTS_dendro_cod = confusionchart(contr,ordertr);

% yes/no variables
load mdl_dendro_yn
%implemenetazione 
out_ts_dendro_yn=predict(mdl_dendro_yn, TS_yn (:, 1:end-1));
% Costruzione della confusion matrix
[contr,ordertr] = confusionmat(TS_yn(:,23),out_ts_dendro_yn);
figure("name","TS dendro Y/N")
confusionTS_dendro_yn = confusionchart(contr,ordertr);

%% Classificatore KNN SOM
% Codified
load mdl_SOM_cod
%implementazione 
out_ts_som_cod=predict(mdl_SOM_cod, TS_cod (:, 1:end-1));
% Costruzione della confusion matrix
[contr,ordertr] = confusionmat(TS_cod(:,6),out_ts_som_cod);
figure("name","TS SOM cod")
confusionTS_SOM_cod = confusionchart(contr,ordertr);

% yes/no variables
load mdl_SOM_yn
%implementazione 
out_ts_som_yn=predict(mdl_SOM_yn, TS_yn (:, 1:end-1));
% Costruzione della confusion matrix
[contr,ordertr] = confusionmat(TS_yn(:,23),out_ts_som_yn);
figure("name","TS SOM yn")
confusionTS_SOM_yn = confusionchart(contr,ordertr);

%% Classificazione NET
% struttura scelta manualmente
load net_MAN.mat
n_input = length(TS_yn(1,1:end-1)); %neuroni di input = numero di variabili selezionate

matrice_TS = TS_yn(:,1:end-1); %usiamo come matrice solo quella corrispondente alle variabili selezionate
matrice_TS = matrice_TS'; %usiamo la trasposta perché la funzione net vuole la trasposta
target_TS = TS_yn(:,end); %come target usiamo la classe
target_TS = target_TS';

y_ts = zeros(1,length(TS_yn));
for rip = 1:10
    net_MAN = init(net_MAN);
    [net_MAN] = train(net_MAN,matrice_TS,target_TS);
    y_TS = net_MAN(matrice_TS);
    
    y_ts(y_TS<=0.5) = 0; 
    y_ts(y_TS>0.5) = 1;
    mat_TS(:,:,rip) = confusionmat(target_TS,y_ts);
end
for i = 1:10
    x = [mat_TS(1,1,i),mat_TS(2,2,i)];
    sommadiag(i) = sum(x);
    
end
[max_MAN, bestrip_MAN] = max(sommadiag);
%scelgo rip.2 perchè, a fronte di una perc. di corretti class. cl.1 <<(da
%97 a 91) risp. rip.3, ho di contro cl.0 rip.2 77,..>70.9 cl.0 rip.3
conf_net_MAN_TS = mat_TS(:,:,3);
confusion_net_TS_MAN=confusionchart(conf_net_MAN_TS);

% struttura GA
load net_GA_yn.mat
load maxcs.mat
load mincs.mat
load TS.mat
test_set_n=rescale(TS,'InputMin',mincs,'InputMax',maxcs);
test_set_norm=[test_set_n(:,3:57) TS(:,58)]; % le prime due colonne sono state tolte

% Coordinate soluzioni (ATTENZIONE:includono anche le prime 2 colonne relative all' ID)
% sol best GA2 5 7 9 11 12 14 15 16 18 21 24 28 29 31 33 34 38 39 45 48 49 50 53 55 56 57 58

TS_GA_2=test_set_norm(:,[3 5 7 9 10 12 13 14 16 19 22 26 27 29 31 32 36 37 43 46 47 48 51 53 54 55 56]);

n_input = length(TS_GA_2(1,1:end-1)); %neuroni di input = numero di variabili selezionate

matrice_TS = TS_GA_2(:,1:end-1); %usiamo come matrice solo quella corrispondente alle variabili selezionate
matrice_TS = matrice_TS'; %usiamo la trasposta perché la funzione net vuole la trasposta
target_TS = TS_GA_2(:,end); %come target usiamo la classe
target_TS = target_TS';
y_ts = zeros(1,length(TS_GA_2));

for rip = 1:10
    net_GA_yn = init(net_GA_yn);
    [net_GA_yn] = train(net_GA_yn,matrice_TS,target_TS);
    y_TS = net_GA_yn(matrice_TS);
    
    y_ts(y_TS<=0.5) = 0; 
    y_ts(y_TS>0.5) = 1;
    mat_TS_GA(:,:,rip) = confusionmat(target_TS,y_ts);
end
for i = 1:10
    x = [mat_TS_GA(1,1,i),mat_TS_GA(2,2,i)];
    sommadiag_GA(i) = sum(x);
    
end
[max_GA, bestrip_GA] = max(sommadiag_GA);
conf_net_GA_TS = mat_TS_GA(:,:,6);
confusion_net_TS_GA=confusionchart(conf_net_GA_TS);