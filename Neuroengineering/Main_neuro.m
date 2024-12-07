%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                NEUROENGINEERING                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

global dataset
global epoch

p = 40; % number of patients

dataset.data = cell(1,p);
dataset.signal = cell(1,p);
dataset.sign_filt = cell(1,p);

epoch.std = cell(1,p);
epoch.odd = cell(1,p);
epoch.clean_odd = cell(1,p);
epoch.clean_std = cell(1,p);

clean_odd_mean = zeros(30,358,p);
clean_std_mean = zeros(30,358,p);

latency_std = zeros(p,1);
latency_odd = zeros(p,1);
pkodd = zeros(p,1);
pkstd = zeros(p,1);
i_pkodd = zeros(p,1);
i_pkstd= zeros(p,1);

odd = zeros(3,30,p);
std = zeros(3,30,p);
mass = zeros(1,p);

labels = {'FP1','F3','F7','FC3','C3','C5','P3','P7','P9','PO7','PO3',...
           'O1','Oz','Pz','CPz','FP2','Fz','F4','F8','FC4','FCz','Cz',...
           'C4','C6','P4','P8','P10','PO8','PO4','O2'};

for i=1:p

    % i = 23;       % In the case, change patients as desired

    %% UPLOADING DATA
    disp(['Valutazione soggetto numero: ' ,int2str(i)]);
    % path = uigetdir('C:/');
    % path = strcat(path,'/S',num2str(i));
    path = strcat("C:\Users\Adriano\Desktop\TUTTO\Scuola\INGEGNERIA\MAGISTRALE\2_ANNO\NEUROENGINEERING\LABs\LAB04 - Project 2 - Event-related potentials\Data\",'\S',num2str(i));
    
    load(strcat(path,"/MNT.mat"));
    load(strcat(path,"/V_event.mat"));
    load(strcat(path,"/V_EEG.mat"));
    
    P9 = V_EEG.trial{1, 1}(9,:);
    P10 = V_EEG.trial{1, 1}(27,:);
    
    dataset.data{i} = V_EEG.trial{1, 1};
    dataset.data{i} = dataset.data{i}-mean([P9;P10]);
    fsample = V_EEG.fsample;
    asset = V_EEG.time{1, 1}; 

    % Let's separate eye movements
    dataset.eog.HEOG_l = V_EEG.trial{1,1}(31,:);
    dataset.eog.HEOG_r = V_EEG.trial{1,1}(32,:);
    dataset.eog.HEOG_low = V_EEG.trial{1,1}(33,:);


    %% BLINK REMOVAL
%     dataset.signal{i} = zeros(30,size(dataset.data{i},2));
%     ind = eyeblink(dataset.eog.HEOG_low,fsample);
%     data_supp = dataset.data{i};
%     data_supp(:,ind) = NaN;
%     for g = 1:30
%         dataset.signal{i}(g,:) = fillmissing(data_supp(g,:),'linear');
%     end

    dataset.signal{i} = dataset.data{i};      % To be commented on if you want to remove ocular artifact


    %% FILTERING
    % first high pass filter and then low pass filter
    fNy = fsample/2;  
    
    wp = 0.5/fNy;
    ws = 0.3/fNy;
    rp = 0.5;
    rs = 20;
    [n,wp]= cheb1ord(wp,ws,rp,rs);
    [b_high,a_high]= cheby1(n,rp,wp,'high');
    % figure,freqz(b_high,a_high,4098,fsample)
    dataset.sign_filt{i} = (filtfilt(b_high,a_high,(dataset.signal{i})'))';
   
    wp = 40/fNy;
    ws = 45/fNy;
    [n,wp] = cheb1ord(wp,ws,rp,rs);
    [b_low,a_low] = cheby1(n,rp,wp,'low');
    % figure,freqz(b_low,a_low,4098,fsample)
    dataset.sign_filt{i} = (filtfilt(b_low,a_low,(dataset.sign_filt{i})'))';


    %% DIVISION INTO EPOCHS ACCORDING TO THE TYPE OF STIMULUS
    value = [V_event.value]; % automatically removes the first line that is blank
    
    % I saved the indices of the normal and oddball stimuli in index, adding 1 because then we can compensate for the first line
    % I EVALUATE FROM THE TRIGGER ON
   
    index_std = find(value==2);
    index_odd = find(value==1);
    
    % except for the onset of events in onset and pre and post
    % I define the intervals in which I want to observe (from -200 ms to 1200 ms)
    onset = round([V_event.onset]*fsample);    % stimulus start (in samples)
    onset_st = onset(index_std);
    onset_odd = onset(index_odd);
    pre_st = onset_st-round(0.2*fsample);       % beginning of epoch
    post_st = onset_st+round(1.2*fsample)-1;    % end of epoch (1 s)
    pre_odd = onset_odd-round(0.2*fsample);
    post_odd = onset_odd+round(1.2*fsample)-1;

    camp = 358;     % number of samples = 1.4*fsamp
    
    % except all indices of the standard stimuli
    ind_st = zeros(1,length(onset_st)*camp);
    for j=0:length(onset(index_std))-1
        ind_st(j*camp+1:(j+1)*camp) = pre_st(j+1):post_st(j+1);
    end
    
    % except all indices of the oddball stimuli 
    ind_odd = zeros(1,length(onset_odd)*camp);
    for j=0:length(onset(index_odd))-1
        ind_odd(j*camp+1:(j+1)*camp) = pre_odd(j+1):post_odd(j+1);
    end
    
    % save the filtered signals and divided into epochs in a single matrix in the format
    % (epochs, 256 samples, channels)
    epoch.std{i} = zeros(length(onset_st),camp,30);
    epoch.odd{i} = zeros(length(onset_odd),camp,30);
    for k=1:30
        for j=0:length(onset_st)-1
            epoch.std{i}(j+1,:,k) = dataset.sign_filt{i}(k,ind_st(j*camp+1:(j+1)*camp));       % standard stimuli epochs  
        end
        for jj=0:length(onset_odd)-1
            epoch.odd{i}(jj+1,:,k) = dataset.sign_filt{i}(k,ind_odd(jj*camp+1:(jj+1)*camp));  % oddball stimuli epochs  
        end
    end

    
    %% DISPLAY ONE ODDBALL AND ONE STANDARD ON ALL CHANNELS
    
    t = linspace(-200,1200,camp);

%     ep = 30;    % oddball epochs
%     figure(1),clf
%     for ch=1:30     % channel
%         subplot(5,6,ch),plot(t,epoch.odd{i}(ep,:,ch),'r'),xlim([min(t) max(t)]),hold on,xline(0,'--k')
%         % title(['ODD - channel ',int2str(ch),' - ',cell2mat(V_EEG.label(ch))])
%         title([cell2mat(V_EEG.label(ch))])
%     end
%      
%     ep = 30;    % standard epochs
%     figure(2),clf
%     for ch=1:30     % channel
%         subplot(8,4,ch),plot(t,epoch.std{i}(ep,:,ch),'b'),xlim([min(t) max(t)]),hold on,xline(0,'--k')
%         title(['STD - channel ',int2str(ch),' - ',cell2mat(V_EEG.label(ch))])
%     end


    %%  DISPLAY ALL THE ODDBALL SIGNALS
    
    % channel Fp1
%     ch = 1;
%     figure(1)
%     clf
%     for ep=1:size(epoch.odd{i},1)
%         subplot(ceil(size(epoch.odd{i},1)/4),4,ep),plot(t,epoch.odd{i}(ep,:,ch),'r'),hold on,xline(0,'--k')%,hold on
%         title(['ODD - epoch ',int2str(ep)])
%     end
% 
%     % channel Fz
%     ch = 17;
%     figure(1)
%     for ep=1:size(epoch.odd{i},1)
%         subplot(ceil(size(epoch.odd{i},1)/4),4,ep),plot(t,epoch.odd{i}(ep,:,ch),'g'),xlim([min(t) max(t)]),hold off
%     end

    % If I wanted to see the standards as well...
%     ch = 17;
%     figure(2)
%     clf
%     for ep=1:size(epoch.std{i},1)
%         subplot(ceil(size(epoch.std{i},1)/10),10,ep),plot(t,epoch.std{i}(ep,:,ch),'b'),hold on,xline(0,'--k')
%         % title(['STD - epoch ',int2str(ep)])
%     end


    %% CORRELATION BETWEEN EPOCHS
    % we want to discard epochs that are very different from the average
    ch = 1;        % 1 = Fp1; 17 = Fz

    pk = max(epoch.odd{i}(:,:,ch),[],2);
    indpk = pk>400;       % I take out epochs with peaks over 400, difficult for them to really be visual stimuli...
    epoch.odd{i}(indpk,:,:) = [];
    
    odd_mean = mean(epoch.odd{i}(:,:,ch)); 

    % I see which epochs are most similar to the average by correlation
    corr = zeros(1,size(epoch.odd{i}(:,:,ch),1));
    for j=1:size(epoch.odd{i}(:,:,ch),1)
        c = xcorr(odd_mean,epoch.odd{i}(j,:,ch));
        corr(j) = max(c);
    end
    corr = corr/max(corr);  %  normalize against the one with maximum correlation

    % figure(3),clf,heatmap(corr)         % represent correlation values with heatmap for easy visualization

    % I take only epochs above threshold, identified with heatmap
    threshold = 0.2;
    index = find(corr>threshold);
    epoch.clean_odd{i} = epoch.odd{i}(index,:,:);

    % I display the epochs left after cleaning
%     ch = 1;
%     figure(4)
%     clf
%     for j=1:size(epoch.clean_odd{i},1)
%         subplot(ceil(size(epoch.clean_odd{i},1)/4),4,j),plot(t,epoch.clean_odd{i}(j,:,ch),'r'),hold on,xline(0,'--k')
%         title(['ODD - epoch ',int2str(index(j))])
%     end

    % I check that something has actually changed
    % figure(5),clf,plot(t,odd_mean,'--r'),hold on,plot(t,mean(epoch.clean_odd{i}(:,:,ch)),'r'),legend('Before cleaning','After cleaning')


    % I do the same thing for standard stimulus
    pk = max(epoch.std{i}(:,:,ch),[],2);
    indpk = find(pk>400);
    epoch.std{i}(indpk,:,:) = []; % only the amplitude control
    
%     std_mean = mean(epoch.std{i}(:,:,ch)); 
%     corr = zeros(1,size(epoch.std{i}(:,:,ch),1));
%     for j=1:size(epoch.std{i}(:,:,ch),1)
%         c = xcorr(std_mean,epoch.std{i}(j,:,ch));
%         corr(j) = max(c);
%     end
%     corr = corr/max(corr);
%     
%     threshold = 0.2;      
%     index = find(corr>threshold);
%     clean_st{i} = epoch.std{i}(index,:,:);
% 
%     % display the epochs left after cleaning
%     ch = 1;
%     figure(6)
%     clf
%     for j=1:size(clean_st{i},1)
%         subplot(ceil(size(clean_st{i},1)/6),6,j),plot(t,clean_st{i}(j,:,ch),'b'),hold on,xline(0,'--k')
%         %title(['STD - epoch ',int2str(index(j))])
%     end
% 
%     
%     figure(7),clf,plot(t,std_mean,'--r'),hold on,plot(t,mean(clean_st{i}(:,:,ch)),'r'),legend('Before cleaning','After cleaning')
    
    epoch.clean_std{i} = epoch.std{i};      


    %% LATENCY
    aux = mean(epoch.clean_odd{i});
    clean_odd_mean(:,:,i) = reshape(aux,camp,30)';  % 3D matrix, I bring it back to 2D
    aux = mean(epoch.clean_std{i});
    clean_std_mean(:,:,i) = reshape(aux,camp,30)';

    ch = 1;
    [pkodd(i), i_pkodd(i)] = max(clean_odd_mean(ch,91:307,i),[],2);  % I take as a window between 150 and 1000 ms, where I expect there to be the peak
    [pkstd(i), i_pkstd(i)] = max(clean_std_mean(ch,91:307,i),[],2);
    
    latency_odd(i) = (i_pkodd(i)+38)/fsample; % I add 38 to start from 0 seconds
    latency_std(i) = (i_pkstd(i)+38)/fsample;


    %% TOPOGRAFY 
    
    i_pkodd(i) = i_pkodd(i)+90;   % I'm adding the 90 taken out earlier
    %sign = (clean_odd_mean(:,i_pkodd(i)+90,i))';    
    % figure(8),clf,plot_topography(labels,sign);

    % I consider the peak odd for both to compare
    % I convert the values to logarithmic scale to see better how it is
    % distributed amplitude across channels and not just a big yellow spot in front
    odd(:,:,i) = abs(log10((clean_odd_mean(:,[i_pkodd(i)-38 i_pkodd(i) i_pkodd(i)+38],i))))';     % I take 150 ms before peak odd, peak instant and 150 ms after
    st(:,:,i) = abs(log10((clean_std_mean(:,[i_pkodd(i)-38 i_pkodd(i) i_pkodd(i)+38],i))))';     % I consider the peak odd for both
    
    % I seek the maximum to optimize scale topography
    if odd(2,1,i)>=st(2,1,i)
        mass(i) = odd(2,1,i);
    else
        mass(i) = st(2,1,i);
    end

end

%% SAVE VARIABLES FOR LATER, I AVOID RUN

% save variables fsample epoch labels odd st mass latency_std latency_odd pkstd pkodd t clean_std_mean clean_odd_mean V_EEG.label
% clear all

%load variables

%% LATENCY

% Plot the latencies for each patient, both standard and oddball
p = 40;
subjects = 1:p;
figure(9),bar(subjects,[latency_std*1000 latency_odd*1000]),title('Latency'),xlabel('# subject'),ylabel('Mean latency (ms)')
legend('Standard','Oddball')
% figure(10),bar([pkstd pkodd]),title('Peak amplitude'),legend('Standard','Oddball')


% find possible outliers
figure(11),boxplot([latency_std*1000 latency_odd*1000],{'standard latency','oddball latency'}),ylabel('time (ms)')
title('Boxplot - searching for outliers')
% figure(12),boxplot([pkstd pkodd],{'standard peak amplitude','amplitude peaks oddball'}),title('Search for possible ouliers')


%% GRAND AVERAGE
% average among all patients

grandavg_odd = mean(clean_odd_mean,3);
grandavg_std = mean(clean_std_mean,3);

ch = 1;
figure(13),plot(t,grandavg_std(ch,:),'b'),hold on,plot(t,grandavg_odd(ch,:),'r'),xline(0,'--k'),hold off
title(['Average across all subjects - channel ',cell2mat(V_EEG.label(ch))]),xlabel('time (s)'),ylabel('amplitude (uV)')
legend('standard','oddball')


%% TOPOGRAPHY
    
i = 8;      % patient
figure(14),clf
subplot(2,3,1),plot_topography(labels,odd(1,:,i));title('topography oddball averaged - 150 ms before'),colorbar,caxis([0 mass(i)])
subplot(2,3,2),plot_topography(labels,odd(2,:,i));title('topography oddball averaged - odd peak instant'),colorbar,caxis([0 mass(i)])
subplot(2,3,3),plot_topography(labels,odd(3,:,i));title('topography oddball averaged - 150 ms after'),colorbar,caxis([0 mass(i)])

subplot(2,3,4),plot_topography(labels,st(1,:,i));title('topography standard averaged - 150 ms before'),colorbar,caxis([0 mass(i)])
subplot(2,3,5),plot_topography(labels,st(2,:,i));title('topography standard averaged - picco odd instant'),colorbar,caxis([0 mass(i)])
subplot(2,3,6),plot_topography(labels,st(3,:,i));title('topography standard averaged - 150 ms after'),colorbar,caxis([0 mass(i)])


%% I DISPLAY MEAN AND STANDARD DEVIATION OF EPOCHS OF A PATIENT

i = 36;
ch = 17;
deviaz_std = std(epoch.std{i}(:,:,ch));
deviaz_odd = std(epoch.odd{i}(:,:,ch));
mean_std = mean(epoch.std{i}(:,:,ch));
mean_odd = mean(epoch.odd{i}(:,:,ch));

figure(15),clf

subplot(1,3,1)
patch([t fliplr(t)],[mean_std-deviaz_std fliplr(mean_std+deviaz_std)],'blue','FaceAlpha',.2),hold on
plot(t,mean_std,'k',LineWidth=1),hold on
plot(t,mean_std+deviaz_std,'b'),hold on
plot(t,mean_std-deviaz_std,'b'),hold off
xlim([min(t) max(t)])
title('Averaged standard averaged +- standard deviation'),xlabel('Time (s)'),ylabel('Amplitude (uV)')

subplot(1,3,2)
patch([t fliplr(t)],[mean_odd-deviaz_odd fliplr(mean_odd+deviaz_odd)],'red','FaceAlpha',.2),hold on
plot(t,mean_odd,'k',LineWidth=1),hold on
plot(t,mean_odd+deviaz_odd,'r'),hold on
plot(t,mean_odd-deviaz_odd,'r'),hold off
xlim([min(t) max(t)])
title('Averaged oddball +- standard deviation'),xlabel('Time (s)'),ylabel('Amplitude (uV)')

subplot(1,3,3)
plot(t,mean_odd-mean_std)
xlim([min(t) max(t)])
title('Averaged oddball - ageraved standard'),xlabel('Time (s)'),ylabel('Amplitude (uV)')


 %% CONTINUOUS WAVELET TRANSFORM
    
% I now consider all epochs, even those removed earlier

i = 4;     % random patient
ch = 17;   % Fz
sign_odd = epoch.odd{i}(30,:,ch);        % random epochs
sign_std = epoch.std{i}(69,:,ch);
[cfs_odd,f_odd] = cwt(sign_odd,fsample); % CWT calculation: continuous wavelet transform
[cfs_std,f_std] = cwt(sign_std,fsample);

figure(15),clf
subplot(3,1,1),plot(t,sign_odd,'r'),hold on,plot(t,sign_std,'b'),axis tight
title("Signals"),xlabel("Time (ms)"),ylabel("Amplitude (uV)")
legend('random oddball','random standard')

subplot(3,1,2),surface(t,f_odd,abs(cfs_odd)),axis tight
subtitle("CWT - ODDBALL"),xlabel("Time (ms)"),ylabel("Frequency (Hz)")
set(gca,"yscale","log"),shading flat
colorbar,caxis([0 max(max(abs(cfs_odd)))])

subplot(3,1,3),surface(t,f_std,abs(cfs_std)),axis tight
subtitle("CWT - STANDARD"),xlabel("Time (ms)"),ylabel("Frequency (Hz)")
set(gca,"yscale","log"),shading flat
colorbar,caxis([0 max(max(abs(cfs_odd)))])

