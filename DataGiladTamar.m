clear;
close;
clc;
%function [mat] = DataGiladTamar(DataGilad,DataTamar,TstartGilad,TstartTamar)
%%
newfilename = "test2";
mkdir (["../data/"+newfilename])
filename = 1; % 1 for choosing the file interactively
folder = 'C:\Users\tamar\OneDrive - Technion\Documents\technion\semester 6\mabada\21April2021\data';
if filename == 0
    filenameG = '../data/02-5-2021_10-43-16-137.txt';
    filenameT = '../data/11May2021 (1).csv';
else
    filenameG = fullfile(folder, uigetfile('../data/*.txt'));
    filenameT = fullfile(folder, uigetfile('../data/*.csv'));
%     filenameG = uigetfile('../data/*.txt');
%     filenameT = uigetfile('../data/*.csv');
    % adding a path: [filename, path] = uigetfile('../data/*.txt'); 
end

DataGilad = readmatrix(filenameG); %28-4-2021_16-08-16-471.txt
% TstartGilad = duration('16:08:16.471','InputFormat',infmt,'Format',outfmt);
DataTamar = table2array(readtable(filenameT)); %Take 2021-04-28 04.08.15 PM.csv
DataTamar = DataTamar(:,(~isnan(DataTamar(1,:))));   % for nan - columns
% TstartTamar = duration('16:08:15.381','InputFormat',infmt,'Format',outfmt);

    %% episode identification
    tolT = 0.002;
    tolG = 0.7;
    % interp1
    DataGilad(:,1) = DataGilad(:,1) - DataGilad(1,1);
    DataGiladInterp = interp1(DataGilad(:,1),DataGilad(:,2:end),0:0.002:round(DataGilad(end-1,1),2));
    DataTamarInterp = interp1(DataTamar(:,2),DataTamar(:,3:end),0:0.002:round(DataGilad(end-1,1),2));
    DataGilad = [(0:0.002:round(DataGilad(end,1),2))',DataGiladInterp];
    DataTamar = [(0:0.002:round(DataGilad(end,1),2))',DataTamarInterp];
    % finding the index of the episode in the data
    startGidx = find(abs(DataGilad(:,2)) >= tolG,1,'first');
    startTidx = find(abs(DataTamar(:,2)-DataTamar(1,2)) >= tolT,1,'first');
    % the start time in [s] of the episode
    TstartG = DataGilad(startGidx,1);
    TstartT = DataTamar(startTidx,1);
    dt = TstartT - TstartG;
    
    if dt > 0
        [d, id] = min(abs(DataTamar(:,1) - dt));
        DataTamar = DataTamar(id:end,:);
        DataTamar(:,1) = DataTamar(:,1) - DataTamar(1,1);
    else
        [d, id] = min(abs(DataGilad(:,1) - dt));
        DataGilad = DataGilad(id:end,:);
        DataGilad(:,1) = DataGilad(:,1) - DataGilad(1,1);
    end
        
    [EIrowG,EIcolG] = size(DataGilad);
    [EIrowT,EIcolT] = size(DataTamar);
    markers = (EIcolT - 1)/3;
    if EIrowG > EIrowT
        r = EIrowG;
    else
        r = EIrowT;
    end
    
    TstartG_new = DataGilad(find(abs(DataGilad(:,2)) >= tolG,1,'first'),1);
    TstartT_new = DataTamar(find(abs(DataTamar(:,4) - DataTamar(1,4)) >= tolT,1,'first'),1);
    
    Dt = abs(TstartG_new -TstartT_new)
    mat = nan(r,EIcolG+EIcolT-1);
    mat(1:EIrowG,1:EIcolG) = DataGilad;
    mat(1:EIrowT,EIcolG+1:EIcolT+EIcolG-1) =  DataTamar(:,2:end);

%% figures
figure;
subplot(321); plot(mat(:,1),mat(:,2)); xlabel('time [s]'); ylabel('a_{x}');
subplot(322); plot(mat(:,1),mat(:,3)); xlabel('time [s]'); ylabel('a_{y}');
subplot(323); plot(mat(:,1),mat(:,4)); xlabel('time [s]'); ylabel('a_{z}');
subplot(324); plot(mat(:,1),mat(:,5)); xlabel('time [s]'); ylabel('g_{x}');
subplot(325); plot(mat(:,1),mat(:,6)); xlabel('time [s]'); ylabel('g_{y}');
subplot(326); plot(mat(:,1),mat(:,7)); xlabel('time [s]'); ylabel('g_{z}');
figure;
plot(mat(:,1),mat(:,8:3:end));
    % fft
if isempty(TstartT_new)
    TstartT_new = 1;
end
% TstartT_new = TstartT_new*(~isempty(TstartT_new)) + 1;        
figure;
for i = [1:10]
    y(:,i) = mat(TstartT_new:end,5+3*i);
    y1_norm = (y(:,i)-y(1,i))./y(1,i);
    y1_norm(isnan(y1_norm)) = [];
    x = mat(:,8);
    Fs = 500;
    Y = fft(y1_norm);
    T = 1/Fs;
    L = length(y1_norm);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    semilogy(f,P1)
    hold on
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    xlim([0,40]);
end

%% data to Tomer
for i = [1:markers]
    row1(2*i - 1:2*i) = ["Curve "+num2str(i),""];
    row2(2*i - 1:2*i) = ["",""];
    row3(2*i - 1:2*i) = ["Standard\Function\Function class","Time"];
    row4(2*i - 1:2*i) = ["Standard\Function\Point id","input:Z"+num2str(i)];
    row5(2*i - 1:2*i) = ["Standard\Function\Sample frequency","240 Hz"];
    row6(2*i - 1:2*i) = ["",""];
    row7(2*i - 1:2*i) = ["Linear","Real(...)"];
    mat_tomer(2:length(mat)+1,2*i-1) = mat(:,1);
    mat_tomer(2:length(mat)+1,2*i) = mat(:,3*i + 5);
%     VariableNames(2*i-1:2*i) =["Time"+num2str(i),"z"+num2str(i)] ;
    Names(2*i - 1:2*i) =["Time","z"+num2str(i)];
end
% mat_tomer(:,21:27) = mat(:,1:7);
% for i = [1:length(mat_lomer)]
TomersTable = array2table([row1; row2; row3; row4; row5; row6; row7; Names; mat_tomer]);
%% saving 
save("../data/"+newfilename+"/mat.mat",'mat');
writetable(TomersTable,"../data/"+newfilename +"/TomersTable.xlsx",'WriteVariableNames',false); 

