clear;
close;
clc;
%function [mat] = DataGiladTamar(DataGilad,DataTamar,TstartGilad,TstartTamar)
outfmt = 'hh:mm:ss.SSS';
infmt = 'hh:mm:ss.SSS';

filename = 0; % 1 for choosing the file interactively
if filename == 0
    filenameG = 'data/28-4-2021_17-05-45-652.txt';
    filenameT = 'data/Take 2021-05-05 04.36.50 PM_001.csv';
else
    filenameG = uigetfile('data/*.txt');
    filenameT = uigetfile('data/*.csv');
    % adding a path: [filename, path] = uigetfile('../data/*.txt'); 
end

DataGilad = readmatrix(filenameG); %28-4-2021_16-08-16-471.txt
TstartGilad = duration('16:08:16.471','InputFormat',infmt,'Format',outfmt);
DataTamar = table2array(readtable(filenameT)); %Take 2021-04-28 04.08.15 PM.csv
DataTamar = DataTamar(:,(~isnan(DataTamar(1,:))));   % for nan - columns
TstartTamar = duration('16:08:15.381','InputFormat',infmt,'Format',outfmt);

i = 0; %0 for ei 1 for delta t
if i == 0
    %% episode identification
    % interp1
    DataGilad(:,1) = DataGilad(:,1) - DataGilad(1,1);
    DataGiladInterp = interp1(DataGilad(:,1),DataGilad(:,2:end),0:0.002:round(DataGilad(end-1,1),2));
    DataTamarInterp = interp1(DataTamar(:,2),DataTamar(:,3:end),0:0.002:round(DataGilad(end-1,1),2));
    DataGilad = [(0:0.002:round(DataGilad(end,1),2))',DataGiladInterp];
    DataTamar = [(0:0.002:round(DataGilad(end,1),2))',DataTamarInterp];
    % finding the index of the episode in the data
    startGidx = find(abs(DataGilad(:,2)) >= 0.7,1,'first');
    startTidx = find(abs(DataTamar(:,4)-DataTamar(1,4)) >= 0.002,1,'first');
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
    if EIrowG > EIrowT
        r = EIrowG;
    else
        r = EIrowT;
    end
    
    TstartG_new = DataGilad(find(abs(DataGilad(:,2)) >= 0.7,1,'first'),1);
    TstartT_new = DataTamar(find(abs(DataTamar(:,4) - DataTamar(1,4)) >= 0.002,1,'first'),1);
    
    Dt = abs(TstartG_new -TstartT_new)
    mat = nan(r,EIcolG+EIcolT-1);
    mat(1:EIrowG,1:EIcolG) = DataGilad;
    mat(1:EIrowT,EIcolG+1:EIcolT+EIcolG-1) =  DataTamar(:,2:end);
    %tests
%     hold on;
%     plot(DataTamar(:,2),DataTamar(:,5)); plot(DataGilad(:,1),DataGilad(:,2));
else
%% delta t option - 0.4[s] inaccuracy :(
    deltaT = TstartTamar - (TstartGilad + seconds(DataGilad(1:1)));
        if deltaT > 0
            StartIdxG = find(DataGilad(:,1) >= seconds(abs(deltaT)),1,'first');
            DataGilad_new = DataGilad(StartIdxG:end,:);
            DataGilad_new(:,1) = DataGilad_new(:,1) - DataGilad_new(1,1);
            DataTamar_new = DataTamar;
            [rG, cG] = size(DataGilad_new);
            [rT, cT] = size(DataTamar_new);
        else
            StartIdxT = find(DataTamar(:,2) >= seconds(abs(deltaT)),1,'first');
            DataTamar(280,3) = nan;
            DataTamar_new = DataTamar(StartIdxT:end,:);
            DataTamar_new(:,2) = DataTamar_new(:,2) - DataTamar_new(1,2);
            DataGilad_new = DataGilad;
            DataGilad_new(:,1) = DataGilad(:,1) - DataGilad(1,1);
            [rT, cT] = size(DataTamar_new);
            [rG, cG] = size(DataGilad_new);
        end

        if [rG > rT]
            r = rG;
        else
            r = rT;
        end

        mat = nan(r,cG+cT-1);
        mat(1:rG,1:cG) = DataGilad_new;
        mat(1:rT,cG+1:cT+cG-1) =  DataTamar_new(:,2:end);
        
    %    figure;
    %    plot(mat(:,8),mat(:,11:3:end)-mat(1,11:3:end)); xlabel('time [s]'); ylabel('Z displacement [mm]');
    %    %     end
end
%% saving mat
save('data/mat.mat','mat');
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