%%
load Subject1.mat
load Subject2.mat

subject_1 = table2array(subject1)';
subject_1 = subject_1(1:19,:);

subject_2 = table2array(subject2)';
subject_2 = subject_2(1:19,:);


%%

fz_aft_filter1 = ALLEEG(4).data(5,:);
stem( fftshift( fft( fz_aft_filter1 ) ) , ".")

%%

fz_aft_filter2 = ALLEEG(10).data(5,:);
stem( fftshift( fft( fz_aft_filter2 ) ) , ".")
%%

epoch1 = create_epoch(ALLEEG(5).data(:,:));

%%

epoch2 = create_epoch(ALLEEG(11).data(:,:));

%%

final_data1 = ALLEEG(6).data([1 5 10 15],:,:);
save final1.mat final_data1

%%

final_data2 = ALLEEG(12).data([1 5 10 15],:,:);
save final2.mat final_data2

%%
load AD.mat
load Normal.mat

[AD_PLV_rare, AD_PLV_freq] = PLV(AD,2,3);
[Normal_PLV_rare, Normal_PLV_freq] = PLV(normal,2,3);

bxplt = [AD_PLV_freq; AD_PLV_rare; Normal_PLV_freq; Normal_PLV_rare];
G = [repmat({'AD freq'},14,1); repmat({'AD rare'},14,1); repmat({'normal freq'},14,1); repmat({'normal rare'},14,1)];
boxplot(bxplt,G)

%%
ADr_fit = fitdist(AD_PLV_rare,'Normal');
ADf_fit = fitdist(AD_PLV_freq,'Normal');
Normalr_fit = fitdist(Normal_PLV_rare,'Normal');
Normalf_fit = fitdist(Normal_PLV_freq,'Normal');
figure
subplot(2,2,1)
plot(ADr_fit)
title('AD rare')
subplot(2,2,2)
plot(ADf_fit)
title('AD freq')
subplot(2,2,3)
plot(Normalr_fit)
title('Normal rare')
subplot(2,2,4)
plot(Normalf_fit)
title('Normal freq')

%%
[~,p_rare] = ttest2(AD_PLV_rare, Normal_PLV_rare);
[~,p_freq] = ttest2(AD_PLV_freq, Normal_PLV_freq);


disp(['p-value for rare: ',num2str(p_rare),'  p-value for frequent:',num2str(p_freq)])
%%
figure
subplot(1,2,1)
polarhistogram( histogram(AD) )
title('AD')
subplot(1,2,2)
polarhistogram( histogram(normal) )
title('Normal')

figure
subplot(1,2,1)
polarhistogram( phaseMean(AD) );
title('AD')
subplot(1,2,2)
polarhistogram( phaseMean(normal) );
title('Normal')


%%

ADr_heatmap = [];
Normalr_heatmap = [];
ADf_heatmap = [];
Normalf_heatmap = [];

for i = 1:4
    for j = 1:i-1
        [ADr_hm, ADf_hm] = PLV(AD,i,j);
        [Normalr_hm, Normalf_hm] = PLV(normal,i,j);


        [~,p_rare_hm(i,j)] = ttest2(ADr_hm, Normalr_hm);
        [~,p_freq_hm(i,j)] = ttest2(ADf_hm, Normalf_hm);

        ADr_heatmap = [ADr_heatmap mean(ADr_hm)];
        Normalr_heatmap = [Normalr_heatmap mean(Normalr_hm)];
        ADf_heatmap = [ADf_heatmap mean(ADf_hm)];
        Normalf_heatmap = [Normalf_heatmap mean(Normalf_hm)];
    end
end
%%
xvalues = {'Fp1 & Fz','Fp1 & Cz','Fz & Cz','Fp1 & Pz','Fz & Pz','Cz & Pz'};
yvalues = {'AD','Normal'};
figure
heatmap(xvalues,yvalues,[ADf_heatmap;Normalf_heatmap])
figure
heatmap(xvalues,yvalues,[ADr_heatmap;Normalr_heatmap])


%%
load MCI.mat

[MCI_PLV_rare, MCI_PLV_freq] = PLV(MCI,2,3);

MCIr_heatmap = [];
MCIf_heatmap = [];

for i = 1:4
    for j = 1:i-1
        [MCIr_hm, MCIf_hm] = PLV(MCI,i,j);

        [~,p_rare_hm_MCI_AD(i,j)] = ttest2(MCIr_hm, ADr_hm);
        [~,p_freq_hm_MCI_AD(i,j)] = ttest2(MCIf_hm, ADf_hm);
        [~,p_rare_hm_MCI_N(i,j)] = ttest2(MCIr_hm, Normalr_hm);
        [~,p_freq_hm_MCI_N(i,j)] = ttest2(MCIf_hm, Normalf_hm);

        MCIr_heatmap = [MCIr_heatmap mean(MCIr_hm)];
        MCIf_heatmap = [MCIf_heatmap mean(MCIf_hm)];
    end
end

%%
xvalues1 = {'Fp1 & Fz','Fp1 & Cz','Fz & Cz','Fp1 & Pz','Fz & Pz','Cz & Pz'};
yvalues1 = {'AD','Normal','MCI'};
figure
heatmap(xvalues1,yvalues1,[ADf_heatmap;Normalf_heatmap;MCIf_heatmap])
figure
heatmap(xvalues1,yvalues1,[ADr_heatmap;Normalr_heatmap;MCIr_heatmap])

%%
figure
subplot(1,2,1)
polarhistogram( histogram(MCI) )

subplot(1,2,2)
polarhistogram( phaseMean(MCI) );

figure
bxplt1 = [MCI_PLV_freq; MCI_PLV_rare];
G1 = [repmat({'MCI freq'},7,1); repmat({'MCI rare'},7,1)];
boxplot(bxplt1,G1)


%%

[mif_AD , mir_AD] = MI(AD,2,3);
[mif_Normal , mir_Normal] = MI(normal,2,3);

bxplt2 = [mif_AD; mir_AD; mif_Normal; mir_Normal];
G2 = [repmat({'AD freq'},14,1); repmat({'AD rare'},14,1); repmat({'normal freq'},14,1); repmat({'normal rare'},14,1)];
boxplot(bxplt2,G2)

%%

ADr_heatmap_MI = [];
Normalr_heatmap_MI = [];
ADf_heatmap_MI = [];
Normalf_heatmap_MI = [];

for i = 1:4
    for j = 1:i-1
        [ADr_hm_MI, ADf_hm_MI] = MI(AD,i,j);
        [Normalr_hm_MI, Normalf_hm_MI] = MI(normal,i,j);


        [~,p_rare_hm(i,j)] = ttest2(ADr_hm_MI, Normalr_hm_MI);
        [~,p_freq_hm(i,j)] = ttest2(ADf_hm_MI, Normalf_hm_MI);

        ADr_heatmap_MI = [ADr_heatmap_MI mean(ADr_hm_MI)];
        Normalr_heatmap_MI = [Normalr_heatmap_MI mean(Normalr_hm_MI)];
        ADf_heatmap_MI = [ADf_heatmap_MI mean(ADf_hm_MI)];
        Normalf_heatmap_MI = [Normalf_heatmap_MI mean(Normalf_hm_MI)];
    end
end

xvalues2 = {'Fp1 & Fz','Fp1 & Cz','Fz & Cz','Fp1 & Pz','Fz & Pz','Cz & Pz'};
yvalues2 = {'AD','Normal'};
figure
heatmap(xvalues2,yvalues2,[ADf_heatmap_MI;Normalf_heatmap_MI])
figure
heatmap(xvalues2,yvalues2,[ADr_heatmap_MI;Normalr_heatmap_MI])

%%

function epoch = create_epoch(data)

    epch = zeros(19,600,120);
    for i = 1:120
        epch(:,1:600,i)=data(:,2000*(i-1)+1401:2000*i);
    end
    epoch = epch;
end



function [PLV_rare,PLV_freq] = PLV(struct,channel1,channel2)

sz = size(struct); 
PLVrare = zeros(sz(2),1);
PLVfreq = zeros(sz(2),1);

for i = 1:sz(2)
    epoch = struct(i).epoch;
    odor = struct(i).odor;
    noisy = struct(i).noisy;
    a = epoch(channel1,:,:);
    b = epoch(channel2,:,:);
    rare =[];
    frequent = [];

    for j = 1:size(odor)
        if(~ismember(j,noisy))
            PLV1 = PLVcalculator(unwrap(angle(hilbert(a(:,:,j)))), unwrap(angle(hilbert(b(:,:,j)))));
            if(odor(j)==1)
                rare = [rare,PLV1];

            else
                frequent = [frequent,PLV1];
            end
        end
    end

    PLVfreq(i) = mean(frequent);
    PLVrare(i) = mean(rare);
end
    PLV_rare = PLVrare;
    PLV_freq = PLVfreq;
end


function [ plv ] = PLVcalculator( phase_sig1, phase_sig2 )

[~, Ntrials] = size(phase_sig1);

e = exp(1i*(phase_sig1 - phase_sig2));
plv = abs(sum(e)) / Ntrials;
end


function hist = histogram(struct)
sz = size(struct); 
sample = randsample(sz(2),1);

epoch = struct(sample).epoch;
odor = struct(sample).odor;
noisy = struct(sample).noisy;
fz = epoch(2,:,:);
cz = epoch(3,:,:);
theta = [];

for j = 1:size(odor)
    if(~ismember(j,noisy))
        if(odor(j)==0)
            theta = [theta unwrap(angle(hilbert(fz(:,:,j))))-unwrap(angle(hilbert(cz(:,:,j))))];
        end
    end
end

hist = theta;

end


function mean_t = phaseMean(struct)
sz = size(struct); 
meanTheta = zeros(sz(2),1);

for i = 1:sz(2)
    epoch = struct(i).epoch;
    odor = struct(i).odor;
    noisy = struct(i).noisy;
    fz = epoch(2,:,:);
    cz = epoch(3,:,:);
    theta = [];

    for j = 1:size(odor)
        if(~ismember(j,noisy))
            if(odor(j)==0)
                theta = [theta unwrap(angle(hilbert(fz(:,:,j))))-unwrap(angle(hilbert(cz(:,:,j))))];
            end
        end
    end

    meanTheta(i) = mean(theta);
end
    mean_t = meanTheta;
end




function mi = MIcalculator(P,Q)

dist =0 ;
if size(Q,1)==1
    Q = Q ./sum(Q);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp,2);
    
    
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./Q);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    dist = abs(sum(temp,2));
end

[~, Ntrials] = size(P);
mi = phase(dist)/log10(Ntrials);

end


function [mif , mir] = MI(struct ,channel1, channel2)

sz = size(struct); 
MIrare = zeros(sz(2),1);
MIfreq = zeros(sz(2),1);

for i = 1:sz(2)
    epoch = struct(i).epoch;
    odor = struct(i).odor;
    noisy = struct(i).noisy;
    a = epoch(channel1,:,:);
    b = epoch(channel2,:,:);
    rare =[];
    frequent = [];

    for j = 1:size(odor)
        if(~ismember(j,noisy))
            MI1 = MIcalculator(a(:,:,j) , b(:,:,j));
            if(odor(j)==1)
                rare = [rare,MI1];

            else
                frequent = [frequent,MI1];
            end
        end
    end

    MIfreq(i) = mean(frequent);
    MIrare(i) = mean(rare);
end
    mir = MIrare;
    mif = MIfreq;
end