 fps=20;
rawF =  ms.RawTraces';
filtF = ms.FiltTraces';
%rawF = filtF;
neuron_num = size(rawF,1);
%spike = deconv_ca;

[name, filepath] = uigetfile('*.csv');
fid = importdata([filepath, name]);

%time_F = fid.data(1:82000,2);
time_F = fid.data(1:end,2);
raw_F = rawF(:,1:end);
xt_F = time_F(1):(1000/fps):time_F(end);
align_T = length(interp1(time_F,raw_F(1,:),xt_F,'nearest'));
alignF = nan(size(raw_F,1),align_T);
%alignspike = nan(size(raw_F,1),align_T);

for i =1:size(raw_F,1)
alignF(i,:) = interp1(time_F,raw_F(i,:),xt_F,'nearest');
%alignspike(i,:) = interp1(time_F,spike(i,:),xt_F,'nearest');
end
 %F=zscore(alignF,1);
%save('align_ca.mat','alignF','alignspike')
 save('align_ca.mat','alignF')
%  F=alignF;
 labels_resize = imresize(labels',[1 size(alignF,2)],'nearest');

%% average activity
labels_resize = imresize(labels',[1 size(alignF,2)],'nearest');
F = zscore(alignF);
neuron_num = size(F,1);
Ronset = find(diff(labels_resize)==-2)+1;
SR=1017.25;
epochLength = 5;

[td_ratio] = calculate_tdratio(EEG,SR,epochLength);
[power_hightheta_toall,power_lowtheta_toall,power_theta_toall,power_delta_toall] = calculate_theta(EEG,SR,epochLength);
resize_td = imresize(td_ratio',[1 size(F,2)],'nearest');
resize_delta = imresize(power_delta_toall',[1 size(F,2)],'nearest');
nr_d = resize_delta(find(labels_resize==3));
r_td = resize_td(find(labels_resize==1));
act_wake = intersect(find(resize_td>=mean(r_td)),find(labels_resize==2));
%qwake = intersect(find(nr_d>=0.5*mean(nr_d)),find(labels_resize==3));
qwake = setdiff(intersect(find(nr_d>=mean(nr_d)),find(labels_resize==3)),act_wake);
%qwake = find(labels_resize==2);

R_max = [];
NR_max = [];
W_max = [];
r_active = [];
r_inactive = [];
for i = 1:size(F,1)
   tmp = F(i,:);
   mean_r(i,1) =  mean(tmp(labels_resize==1));
   mean_nr(i,1) = mean(tmp(labels_resize==3));
   mean_w(i,1) = mean(tmp(labels_resize==2));
   mean_aw(i,1) = mean(tmp(act_wake));
   mean_qw(i,1) = mean(tmp(qwake));
   mean_all(i,1) = mean(tmp);
%    if max([mean_r(i,1) mean_nr(i,1) mean_w(i,1)])==mean_r(i,1)
%        R_max = [R_max i];
%    end
%    if max([mean_r(i,1) mean_nr(i,1) mean_w(i,1)])==mean_nr(i,1)
%        NR_max = [NR_max i];
%    end
%       if max([mean_r(i,1) mean_nr(i,1) mean_w(i,1)])==mean_w(i,1)
%        W_max = [W_max i];
%    end
ca_std(i,1) = std(F(i,:));
% if mean_r(i)-mean_nr(i)>=0.25*ca_std(i,1)&&mean_r(i)-mean_w(i)>=0.25*ca_std(i,1)
% r_active = [r_active;i];
% end
% if mean_r(i)-mean_nr(i)<-0.25*ca_std(i,1)&&mean_r(i)-mean_w(i)<-0.25*ca_std(i,1)
% r_inactive = [r_inactive;i];
% end
   clear tmp
end
r_active = union(find(mean_r-mean_nr>0.5),find(mean_r-mean_w>0.5));
r_inactive = union(find(mean_r-mean_nr<-0.5),find(mean_r-mean_w<-0.5));
othercells = setdiff(1:neuron_num',[r_active;r_inactive]);
figure
y1=-1.5;
y2=1.5;
line([y1 1.5],[y1 1.5],'Color','k')
hold on
line([y1 1.5],[0 0],'Color','k')
hold on
line([0 0],[y1 1.5],'Color','k')
% scatter(mean_r(r_active)-mean_nr(r_active),mean_r(r_active)-mean_w(r_active),[],[0.9843 0.5608 0.1686],'filled');%REM active
% hold on
% scatter(mean_r(r_inactive)-mean_nr(r_inactive),mean_r(r_inactive)-mean_w(r_inactive),[],[0.2706 0.4118 0.5647],'filled')
% hold on
% scatter(mean_r(othercells)-mean_nr(othercells),mean_r(othercells)-mean_w(othercells),[],[0.7 0.7 0.7],'filled')
scatter(mean_r(r_active)-mean_nr(r_active),mean_w(r_active)-mean_nr(r_active),[],[0.9843 0.5608 0.1686],'filled');%REM active
hold on
scatter(mean_r(r_inactive)-mean_nr(r_inactive),mean_w(r_inactive)-mean_nr(r_inactive),[],[0.2706 0.4118 0.5647],'filled')
hold on
scatter(mean_r(othercells)-mean_nr(othercells),mean_w(othercells)-mean_nr(othercells),[],[0.7 0.7 0.7],'filled')

axis tight
xlabel('REM-NREM')
ylabel('NREM-WAKE')
ylim([-1.5 1.5])
xlim([-1.5 1.5])
%states_avg = [mean_nr mean_r mean_w];
states_avg = [mean_r-mean_nr mean_r-mean_w];

figure
y1=-3;
y2=3;
line([y1 3],[0 0],'Color','k')
hold on
line([0 0],[y1 3],'Color','k')
scatter(mean_r(r_active)-mean_nr(r_active),mean_aw(r_active)-mean_nr(r_active),[],[0.9843 0.5608 0.1686],'filled');%REM active
hold on
scatter(mean_r(othercells)-mean_nr(othercells),mean_aw(othercells)-mean_nr(othercells),[],[0.7 0.7 0.7],'filled')
hold on
scatter(mean_r(r_inactive)-mean_nr(r_inactive),mean_aw(r_inactive)-mean_nr(r_inactive),[],[0.2706 0.4118 0.5647],'filled');%REM active
xlabel('REM - NREM')
ylabel('Active wake - NREM')

save('Ractive_R_inactive_idx.mat','r_active','r_inactive','othercells')
%% plot EEG brainstates with ca heatmap
labels_resize = imresize(labels',[1 size(alignF,2)],'nearest');
F=zscore(alignF);
fps=20;
%F=alignF;
%F = alignspike;
%F=processed_signal;
cell_num = size(F,1);
epochLength = 5; % brain state epoch length
SR = 1017.25; % EEG/EMG sampling rate
f1 = 1; % first frequency band of spectrogram to plot, in Hz
f2 = 30;
colors = [43 160 220;
            129 131 132
            255 194 61
            0 0 0] ./ 255;

T = length(labels_resize)/fps/3600;
x_t = T/size(F,2):T/size(F,2):T;
figure('Position', [200   200   1200   500]); %for 20 frames
a1 = axes('Position',[.05 .68 .9 .025]);
imagesc(x_t,1,labels_resize,[1 4])
xlabel('Time/h')
colormap(a1,colors);
axis off
box off
hold on
[s, t, f] = createSpectrogram(EEG, SR, epochLength);
[~,f1Idx] = min(abs(f - f1)); 
[~,f2Idx] = min(abs(f - f2)); 
s = s(1:end, f1Idx:f2Idx);
f = f(f1Idx:f2Idx);

% Plot EEG
t_eeg = T/size(s,1):T/size(s,1):T;
a_eeg = axes('Position',[.05 .74 .9 .18]);
imagesc(t_eeg,f,s'); axis xy
set(gca,'XTick',[],'YTick',[],'CLim',[0,max(max(s))/10]); 

%plot heatmap
a_cell = axes('Position',[.05 .05 .9 .6]);
t_ca = T/length(alignF):T/length(alignF):T;
% [sortedavge1, sortIndex1] =sortrows(mean_r,'descend');
% sort_trace = alignF(sortIndex1,:);
% imagesc(zscore(sort_trace,1))
% F=zscore(alignF')';
clusterF = F([r_active;setdiff(1:cell_num,[r_active;r_inactive])';r_inactive],:);
%clusterF = F([r_active;r_inactive],:);
%clusterF = F;
imagesc(t_ca,1:size(clusterF,1),clusterF);
%imagesc(t_ca,1:cell_num,F);
set(gca,'TickLength',[0,0.025]); 
%colormap(a_cell,hot)
caxis([-1 1.5])
linkaxes([a1 a_eeg a_cell],'x')
hold on
line([t_ca(1) t_ca(end)],[length(r_active)+0.5 length(r_active)+0.5],'Color','w','LineWidth',2)
line([t_ca(1) t_ca(end)],[length(r_active)+length(othercells)+0.5 length(r_active)+length(othercells)+0.5],'Color','w','LineWidth',2)
%line([t_ca(1) t_ca(end)],[length(r_active)+length(othercells)+length(r_inactive) length(r_active)+length(othercells)+length(r_inactive)],'Color','w','LineWidth',2)
%%
barplot = [mean_r mean_nr mean_w mean_aw,mean_all,mean_qw];

%% nrem to rem to wake
% cell_num = size(F,1);
% targ_diff=-2;
% after=1;
% before=3;
% slidin_window=30;
% window=30;
% wake_window=15;
% nrem_window=60;
% %draw_neurons = [7];
% cluster_color = [0.9843 0.5608 0.1686;
%     0.2706 0.4118 0.5647;
%     0.7 0.7 0.7];
% 
%     remonset = find(diff(labels_resize)==targ_diff)+1;
% 
%     for i = 1:length(remonset)
%        rem_dur(i) = find(labels_resize(remonset(i):end)~=after,1,'first');
%        if isempty(find(labels_resize(remonset(i)+rem_dur(i):end)~=2,1,'first'))
%            %postR_wake(i) = length(labels_resize) - remonset(i)+rem_dur(i)-2;
%            postR_wake(i) =size(F,2)-find(labels_resize(remonset(i):end)~=1,1,'first')-remonset(i);
% 
%        else
%        postR_wake(i) = find(labels_resize(remonset(i)+rem_dur(i):end)~=2,1,'first');
%        end
%        if isempty(find(labels_resize(1:remonset(i)-1)~=3,1,'last'))
%            beforeR_nrem(i)=1;
%        else
%        beforeR_nrem(i)=find(labels_resize(1:remonset(i)-1)~=3,1,'last');
%        end
%     end
%     
%     remain = find(rem_dur>=slidin_window*fps);
%     
%     rem=remonset(remain);
%     rem_duration = rem_dur(remain);
%     postR_wake = postR_wake(remain);
%     beforeR_nrem = beforeR_nrem(remain);
% 
%    
%     n_rem = length(rem);
%     norm_rem = nan(n_rem,slidin_window*fps);
%     norm_wake = nan(n_rem,wake_window*fps);
%     norm_nrem = nan(n_rem,nrem_window*fps);
% cell_n = length(draw_neurons);
% trace = nan(cell_n,(nrem_window+window+wake_window)*fps);
% for icell = 1:cell_n
%     for irem = 1:n_rem
%      tmp1= mean(F(icell,rem(irem):rem(irem)+rem_duration(irem)),1);
%      tmp2 = mean(F(icell,rem(irem)+rem_duration(irem):rem(irem)+rem_duration(irem)+postR_wake(irem)),1);
%      tmp3 = mean(F(icell,beforeR_nrem(irem):rem(irem)-1),1);
%      norm_rem(irem,:) = imresize(tmp1,[1 slidin_window*fps],'nearest');
%      if length(tmp2)<=size(norm_wake,2)
%          norm_wake(irem,1:length(tmp2)) = tmp2;
%      else
%          norm_wake(irem,:) = tmp2(1:wake_window*fps);
%      end
%      if length(tmp3)<=size(norm_nrem,2)
%          norm_nrem(irem,size(norm_nrem,2)+1-length(tmp3):end) = tmp3;
%      else
%          norm_nrem(irem,:) = tmp3(length(tmp3)-nrem_window*fps+1:end);
%      end
%     end
%     tmp_trace = cat(2,norm_nrem,norm_rem,norm_wake);
%     trace(icell,:) = tmp_trace;
% end
%     t=-nrem_window+1/fps:1/fps:window+wake_window;
%     avg = mean(trace,1);
%     sem = std(trace)/(sqrt(cell_n)-1);
% s_avg = smooth(avg,20)';
% s_sem =smooth(sem,20)';
%     %figure('Position',[50,50,400,350])
%     l = plot(t,s_avg,'Color',cluster_color(3,:),'LineWidth',2);hold on
%     h= fill([t fliplr(t)],[s_avg+s_sem fliplr(s_avg-s_sem)],cluster_color(3,:));
%     set(h,'edgealpha',0,'facealpha',0.2)
%     hold on
%     xlim([-30 50])
%     ylabel('zscore')
%     xlabel('Time from REM onset(s)')
%%
fps=20;
cell_num = size(F,1);
targ_diff=-2;
after=1;
before=3;
slidin_window=30;
window=30;
wakewindow=15;
nremwindow=60;
cluster_color = [0.9843 0.5608 0.1686;
    0.2706 0.4118 0.5647;
    0.7 0.7 0.7];

    remonset = find(diff(labels_resize)==targ_diff)+1;
    remoffset = intersect(find(diff(labels_resize)==1),find(labels_resize==1));

    remain=[];
    for i = 1:length(remonset)
        if all(labels_resize(remonset(i):remonset(i)+window*fps)==1)&&all(labels_resize(remonset(i)-nremwindow*fps+1:remonset(i)-1)==3)&&all(labels_resize(remoffset(i)+1:remoffset(i)+wakewindow*fps)==2)
            remain = [remain i];
        end
    end
    
    
    remonset=remonset(remain);
    remoffset = remoffset(remain);

   
    n_rem = length(remonset);
    cell_n = size(F,1);
    NtoR = nan(cell_n,(nremwindow+window)*fps);
    RtoW = nan(cell_n,(window+wakewindow)*fps);
    
for icell = 1:cell_n
    tmp_ntor = [];
    tmp_rtow = [];
    for irem = 1:n_rem
     tmp_r1= F(icell,remonset(irem):remonset(irem)+window*fps-1);
     tmp_nr = F(icell,remonset(irem)-nremwindow*fps:remonset(irem)-1);
     tmp_r2 = F(icell,remoffset(irem)-window*fps:remoffset(irem));
     tmp_w = F(icell,remoffset(irem)+1:remoffset(irem)+wakewindow*fps-1);
     tmp_ntor = [tmp_ntor;[tmp_nr tmp_r1]];
     tmp_rtow = [tmp_rtow;[tmp_r2 tmp_w]];
    end
    NtoR(icell,:)=mean(tmp_ntor,1);
    RtoW(icell,:)=mean(tmp_rtow,1);
end

save('transition.mat','NtoR','RtoW','nremwindow','window','wakewindow','fps')



%%
% 假设：
% F : N × T 的神经元活动矩阵
% labels_resize : 1 × T 的脑状态标签（1=REM, 2=wake, 3=NREM）
% fps = 20;

% 参数设置
fps = 20;
nrem_dur = 60;   % 单位：秒
wake_dur = 15;

n_nrem = nrem_dur * fps;   % 30秒对应帧数
n_wake = wake_dur * fps;   % 15秒对应帧数

N = size(F, 1);
T = size(F, 2);

% 初始化存储每次转换的片段
NtoW_segments = [];
WtoNR_segments = [];

% 遍历时间点
for t = 1:(T - n_nrem - n_wake)
    seg1 = labels_resize(t : t + n_nrem - 1);
    seg2 = labels_resize(t + n_nrem : t + n_nrem + n_wake - 1);

    % --- NREM → Wake ---
    if all(seg1 == 3) && all(seg2 == 2)
%         % 前后状态合法（不直接 REM→wake）
%         if t > 1 && labels_resize(t-1) == 2
%             continue;  % 跳过 wake→NREM→wake 快速回跳
%         end
        seg_data = F(:, t : t + n_nrem + n_wake - 1);
        NtoW_segments(:, :, end+1) = seg_data;  % N × L × M
    end

    % --- Wake → NREM ---
    if all(seg1 == 2) && all(seg2 == 3)
        if t > 1 && labels_resize(t-1) == 1
            continue;  % 跳过 REM→wake→NREM 不合规序列
        end
        seg_data = F(:, t : t + n_nrem + n_wake - 1);
        WtoNR_segments(:, :, end+1) = seg_data;  % N × L × M
    end
end

% 计算平均活动（在每个神经元上对每次 transition 平均）
% 结果是 N × (30+15)*fps 的平均轨迹
if ~isempty(NtoW_segments)
    NtoW = mean(NtoW_segments, 3);  % N × T_window
else
    NtoW = [];
end

if ~isempty(WtoNR_segments)
    WtoNR = mean(WtoNR_segments, 3);  % N × T_window
else
    WtoNR = [];
end
save('NRWake_transition.mat','NtoW','WtoNR','nrem_dur','wake_dur','fps')
