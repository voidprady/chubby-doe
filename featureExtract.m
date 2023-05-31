%matlab file for feature extraction from raw PPG signals.
%ENGR 295B : Master's project - Spring 2023.

addpath(genpath('lib/')); %for signal preprocessing 
SIGNAL_PATH = 'Data File/0_subject/'; %import data

% Patient Info for PPG Signals.

%opts = detectImportOptions('normalPatientInfo.csv');
%opts = detectImportOptions('diabetesPatientInfo.csv');
opts = detectImportOptions('HypertensionPatientInfo.csv');
opts.DataLines = 2;
%subject = readmatrix('normalPatientInfo.csv', opts);
subject  = readmatrix('HypertensionPatientInfo.csv',opts);
%subject = readmatrix('diabetesPatientInfo.csv', opts);

[m,n] = size(subject);
feature_vector = [];
disp(['total number of records : ', int2str(m)])

% Use a threshold to select the other points. Threshold dependent on the dataset.
% 10 - Hypertension; 5 - otherwise
THRESHOLD = 10;

% Loop over the Patients. 
for patient = 1:1:m

    % Load the Signal
    subject_ID  = int2str(subject(patient,1));
    selected    = int2str(subject(patient,15));
    tempSubject = strcat(SIGNAL_PATH,subject_ID,'_',selected,'.txt');  
    raw_PPG       = load(tempSubject);
    
    % Preprocessing - 4th Order Chebyshev Filter 
    new_PPGsignal = preprocessSignal(raw_PPG);
    PPGsignal     = new_PPGsignal;

    figure
    subplot(4,1,1)
    plot(raw_PPG),grid;
    title(strcat('Subject : ',subject_ID,' RAW PPG Signal'));
    
    % 1st Derivative
    x = 1:1:length(PPGsignal);
    dy = diff(PPGsignal)./diff(x);
    dy = preprocessSignal(dy); 
    subplot(4,1,2)
    plot(dy),grid;
    title('VPG WAVE');
    
    % 2nd derivative
    x2 = 1:1:length(dy);
    dy2 = diff(dy)./diff(x2);
    % Preprocessing the 2nd derivate - 4th Order Chebyshev Filter
    dy2 = preprocessSignal(dy2); 
    
    offset = min(dy2);
    new_dy2 = dy2 + abs(offset)*0 ; %offset -> 0 
    
    %figure
    subplot(4,1,3);
    new_dy2 = smooth(new_dy2);
    plot(new_dy2), grid on;
    title('Accelerated Photoplethysmography Signal (APG)');
    xlabel('Time');
    ylabel('Amplitude');
    hold on;

    % Detecting the Peaks of the APG Waveform. 
    % peak_loc = [];
    % peak_val = [];
    % trough_loc = [];
    % trough_val = [];
    
    [peak_val,peak_loc] = findpeaks(new_dy2);
    % Invert the signal to identify the Troughs.
    negAPGsignal = -1 * new_dy2;
    [trough_val, trough_loc] = findpeaks(negAPGsignal);

    % for i=2:1:length(new_dy2)-1
    %     cur = new_dy2(i);
    %     prev = new_dy2(i-1);
    %     next = new_dy2(i+1);
    %     if prev < cur && cur > next
    %         peak_loc = [peak_loc i];
    %         peak_val = [peak_val new_dy2(i)];
    %         %fprintf('a: [%s]\n', join(string(peak_loc), ','));
    %     end
    %     if prev > cur && cur < next
    %         trough_loc = [trough_loc i];
    %         trough_val = [trough_val new_dy2(i)];
    %         %fprintf('a: [%s]\n', join(string(trough_loc), ','));
    %     end
    % end

    % Find the first three maximum points in desceding order. 
    % Roughly there will be two peaks in 2.1s APG waveform.
    % Select the two maximum peaks and then need to find the first occurence of the alpha wave.
    % max3 used for the second peak detection.

    max_finder = sort(peak_val,'descend');
    d = 1; e = 1; f = 1;

    for i=1:1:length(max_finder)-2
        max1 = max_finder(i)
        max2 = max_finder(i+1)
        max3 = max_finder(i+2)
        
        %[max1, ind1] = max(max_finder);
        %max_finder(ind1) = -100000;
        %[max2, ind2] = max(max_finder);
        %max_finder(ind2) = -100000;
        %[max3, ind3] = max(max_finder);

        max1_idx = find(peak_val == max1);
        max2_idx = find(peak_val == max2);
        max3_idx = find(peak_val == max3);

        max1_time = peak_loc(max1_idx);
        max2_time = peak_loc(max2_idx);
        max3_time = peak_loc(max3_idx);
    
        %peak points identification

        a = min(max1_time, max2_time); %first peak point occured must be the a-peak index.

        for i=1:1:length(peak_loc)
            if peak_loc(i) > a
                c = peak_loc(i); % after a-peak, the next peak occured must be the c-peak index.
                for j=i+1:1:length(peak_loc)
                    if (peak_loc(j) - c )> THRESHOLD  %threshold
                        e = peak_loc(j); %if gt threshold, it is labelled as the e-peak
                    end
                    break;
                end
                if e~= 1
                    break;
                else
                    
                end    
            end
        end

        if e == 1
            continue;
        end
    
    %crest points identification
        for i=1:1:length(trough_loc)
            if trough_loc(i) > a 
                b = trough_loc(i); %first crest after the a-peak index must be the b-crest point
                for j=i+1:1:length(trough_loc)
                    % we can observe that the d-point occurs after the b-crest
                    % and c-peak. 
                    if ((trough_loc(j) - c) > 0 && (trough_loc(j) - b )> THRESHOLD)    %threshold
                        d = trough_loc(j);
                        for k=j+1:1:length(trough_loc)
                            if (trough_loc(k) - d) > THRESHOLD %threshold
                                f = trough_loc(k);
                            end
                            break;
                        end
                    end
                    break;
                end
                break;
            end
        end

        if d == 1 || f == 1
            continue;
        end

    %Identifying the Second Peak of the APG waveform.
    % Required when calculating the features ex: Peak to Peak
        a1 = max(max1_time, max2_time);
        % if ind3 occured before a1 then ind3 is a1
        if a < max3_time && max3_time < a1 && max3_time > f && max3 > new_dy2(e)
            a1 = max3_time;
        end

        if a1 == c
            continue;
        end

        if d ~= 1 && e ~= 1 && f ~= 1
            break;
        end
    end
    
    if d == 1 || e == 1 || f == 1
        disp("d,e,f are not found!");
    end
    %plot the identified points
    plot(a, new_dy2(a),'r.','MarkerSize',12);
    hold on;
    plot(b, new_dy2(b),'b.','MarkerSize',12);
    hold on;
    plot(c, new_dy2(c),'k.','MarkerSize',12);
    hold on;
    plot(d, new_dy2(d),'g.','MarkerSize',12);
    hold on;
    plot(e, new_dy2(e),'c.','MarkerSize',12);
    hold on;
    plot(f, new_dy2(f),'m.','MarkerSize',12);
    hold on;
    plot(a1, new_dy2(a1),'y.','MarkerSize',14);
    legend({'Signal','a','b','c','d','e','f', 'a2'},'Location','eastoutside','Orientation','vertical');


    % Feature extraction from the main PPG signal. some points can be infered from the APG waveform. 
    dicro_notch = e;
    dicro_peak = f;
    
    % find the peaks(1st and second) of the PPG signal.   
    PPGsignal = smooth(PPGsignal);
    [ppg_peak_val,ppg_peaks_loc] = findpeaks(PPGsignal);
    % Invert the signal to identify the Troughs.
    negPPGsignal = -1 * PPGsignal;
    [ppg_trough_val, ppg_trough_loc] = findpeaks(negPPGsignal);
    
    amp = 1;
    for i=1:1:length(ppg_peaks_loc)
        if ppg_peaks_loc(i) > a
            amp = ppg_peaks_loc(i);
            break;
        end
    end

    amp1 = 1;
    for i=1:1:length(ppg_peaks_loc)
        if ppg_peaks_loc(i) > a1
            amp1 = ppg_peaks_loc(i);
            break;
        end
    end

    low1 = -1; low2 = -1;
    for i=1:1:length(ppg_trough_loc)
        if ppg_trough_loc(i) < a
            low1 = ppg_trough_loc(i);
        end
        if ppg_trough_loc(i) < a1
            low2 = ppg_trough_loc(i);
        end
    end
    
    % Detrending the PPG. 
    temp_PPGsignal = detrend(PPGsignal,'linear');
    PPGsignal = temp_PPGsignal;

    % unknow issue with low becoming negative.
    if low2 < 0
        low2 = 1;
    end
    
    subplot(4,1,4)
    plot(PPGsignal),grid;
    title('Photoplethysmography Signal (PPG)');
    xlabel('Time');
    ylabel('Amplitude');
    hold on;
    
    plot(dicro_notch, PPGsignal(dicro_notch),'b.','MarkerSize',14);
    plot(dicro_peak, PPGsignal(dicro_peak),'k.','MarkerSize',14);
    plot(amp, PPGsignal(amp),'r.','MarkerSize',14);

    if low1 ~= -1
        plot(low1, PPGsignal(low1),'g.','MarkerSize',14);
    end
    if amp1 ~= -1
        plot(amp1, PPGsignal(amp1),'r.','MarkerSize',14);
    end
    if low2 ~= -1  
        plot(low2, PPGsignal(low2),'g.','MarkerSize',14);
    end
    % delete
    %legend(h,'Dicrotic Notch','Diastolic Peak','PPG Signal Start','Systolic Peak')
    
    %%% Morphological Feature Calculation. %%%
    % Use the extracted points if the two waveforms to claculate features. 
    
    %Heights of the a,b,c,d,e waves
    A = new_dy2(a);
    B = new_dy2(b);
    C = new_dy2(c);
    D = new_dy2(d);
    E = new_dy2(e);
    
    % PPG Features
    if low1 ~= -1
        SysAmp = PPGsignal(amp) - low1;
    elseif amp1 ~= -1
        SysAmp = PPGsignal(amp1) - low1;
    else 
        SysAmp = -1;
    end
    DysAmp = PPGsignal(dicro_peak);
    height = subject(patient,4);
   
    %calculating PI, If PI cant be calculated then the p2p is taken
    if low1 ~= -1
        PI = low2 - low1;
    elseif amp1 ~= -1
        PI = amp1 - amp;
    else
        PI = -1;
    end
    PI_Sys = PI /SysAmp ;
    AI = DysAmp / SysAmp;
    adj_AI = (SysAmp - DysAmp) / SysAmp ;
    deltaT = dicro_peak - amp;
    ArtStiff = height / deltaT;
    
    %Calculate the RT
    if low1 ~= -1
        RT = amp - low1;
    elseif amp1 ~= -1
        RT = amp1 - low2;
    else
        RT = -1;
    end
    
    % Pulse Area Related Features
    pulse2 = PPGsignal(amp:dicro_notch);
    area2 = trapz(pulse2);
    if low2 > dicro_notch
        pulse3 = PPGsignal(dicro_notch:low2);
        area3 = trapz(pulse3);
    else
        pulse3 = PPGsignal(low2:dicro_notch);
        area3 = trapz(pulse3);
    end    
    
    if low1 ~= -1 && low1 < amp
        pulse1 = PPGsignal(low1:amp);
        area1 = trapz(pulse1);
        TotArea = area1 + area2 + area3;
        AreaRatio = (area1+area2) / area3;
    elseif amp1 ~= -1 && low2 ~= -1 && low2 < amp1
        pulse1 = PPGsignal(low2:amp1);
        area1 = trapz(pulse1);
        TotArea = area1 + area2 + area3;
        AreaRatio = (area1+area2) / area3;
    else
        TotArea = -1;
        AreaRatio = -1;
    end
    
    %pulse width related features => can be incoporated in future. 
    %halfAmp = SysAmp / 2;
    %index = find(PPGsignal == halfAmp);
    %disp(index)
   
    extractedFeatures = [B/A C/A D/A E/A (B-C-D-E)/A (B-E)/A (B-C-D)/A (C+D-B)/A (a1-a) SysAmp TotArea AreaRatio PI PI_Sys AI adj_AI ArtStiff RT];    
    temp = [subject(patient,:) extractedFeatures ];
    feature_vector = [feature_vector ; temp];
    
end

% Write to CSV
feat_titles = {'subject_ID' 'Sex(M/F)' 'Age(year)' 'Height(cm)' 'Weight(kg)' 'Heart Rate(b/m)' 'BMI(kg/m^2)' 'Systolic Blood Pressure(mmHg)' 'Diastolic Blood Pressure(mmHg)' 'Hypertension' 'Diabetes' 'x1' 'x2' 'x3' 'Selected' 'B/A' 'C/A' 'D/A' 'E/A' '(B-C-D-E)/A' '(B-E)/A' '(B-C-D)/A' '(C+D-B)/A' '(a1-a)' 'SysAmp' 'TotArea' 'AreaRatio' 'PI' 'PI_Sys' 'AI' 'adj_AI' 'ArtStiff' 'RT'};
csv = [feat_titles; num2cell(feature_vector)];
writecell(csv, 'hypertensionAPGFeaturesV2.csv');
