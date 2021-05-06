clear;
close;
clc;
%function [mat] = DataGiladTamar(DataGilad,DataTamar,TstartGilad,TstartTamar)
outfmt = 'hh:mm:ss.SSS';
infmt = 'hh:mm:ss.SSS';
DataGilad = readmatrix('28-4-2021_17-05-45-652.txt'); %28-4-2021_16-08-16-471.txt
TstartGilad = duration('16:08:16.471','InputFormat',infmt,'Format',outfmt);
DataTamar = table2array(readtable('Take 2021-05-05 04.36.50 PM_001.csv')); %Take 2021-04-28 04.08.15 PM.csv
DataTamar = DataTamar(:,(~isnan(DataTamar(1,:))));   % for nan - columns
TstartTamar = duration('16:08:15.381','InputFormat',infmt,'Format',outfmt);

i = 0; %0 for ei 1 for delta t
if i == 0
    %% episode identification
    % finding the index of the episode in the data
    startGidx = find(abs(DataGilad(:,2)) >= 0.7,1,'first');
    startTidx = find(abs(DataTamar(:,5)-DataTamar(1,5)) >= 0.002,1,'first');
    % the start time in [s] of the episode
    TstartG = DataGilad(startGidx,1);
    TstartT = DataTamar(startTidx,2);
    dt = TstartT - TstartG;
    
    if dt > 0
        [d, id] = min(abs(DataTamar(:,2) - dt));
        DataTamar = DataTamar(id:end,:);
        DataTamar(:,2) = DataTamar(:,2) - DataTamar(1,2);
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
%     DataTamarnew(:,1) = resample(DataTamar(:,2),3,1);
%     DataTamarnew(:,2:64) = resample(DataTamar(:,3:end),3,1);
    TstartT_new = DataTamar(find(abs(DataTamar(:,5) - DataTamar(1,5)) >= 0.002,1,'first'),2);
    
    Dt = abs(TstartG_new -TstartT_new)
    mat = nan(r,EIcolG+EIcolT-1);
    mat(1:EIrowG,1:EIcolG) = DataGilad;
    mat(1:EIrowT,EIcolG+1:EIcolT+EIcolG-1) =  DataTamar(:,2:end);
    %tests
    hold on;
    plot(DataTamar(:,2),DataTamar(:,5)); plot(DataGilad(:,1),DataGilad(:,2));
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

       
    %    
    %    figure;
    %    plot(mat(:,8),mat(:,11:3:end)-mat(1,11:3:end)); xlabel('time [s]'); ylabel('Z displacement [mm]');
    %    %     end
end
 
figure;
subplot(321); plot(mat(:,1),mat(:,2)); xlabel('time [s]'); ylabel('a_{x}');
subplot(322); plot(mat(:,1),mat(:,3)); xlabel('time [s]'); ylabel('a_{y}');
subplot(323); plot(mat(:,1),mat(:,4)); xlabel('time [s]'); ylabel('a_{z}');
subplot(324); plot(mat(:,1),mat(:,5)); xlabel('time [s]'); ylabel('g_{x}');
subplot(325); plot(mat(:,1),mat(:,6)); xlabel('time [s]'); ylabel('g_{y}');
subplot(326); plot(mat(:,1),mat(:,7)); xlabel('time [s]'); ylabel('g_{z}');
figure;
plot(mat(:,8),mat(:,11:3:end));