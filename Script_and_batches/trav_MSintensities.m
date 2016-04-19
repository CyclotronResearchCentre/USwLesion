
%% Play around with voxel intensities in GM/WM/lesion

% Refs:
% https://en.wikipedia.org/wiki/Skewness
% https://en.wikipedia.org/wiki/Kurtosis
% https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test

% Ideas:
% - compare the distribution of intensities with Two-sample Kolmogorov-
% Smirnov test. Function in Matlab is "kstest2"

load('D:\ccc_DATA\MS_ELommers\Data2_unwTPM\ExP_s7480-0005-00001-000176-00_MT_EPIB1.mat')
load('D:\ccc_DATA\MS_ELommers\Data1\ExP_kts7374-0005-00001-000176-00_MT_EPIB1.mat')

if isfield(res,'vMPMmsk')
    v_use = res.vMPMmsk;
else
    v_use = res.vMPM;
end

% Number of voxels per tissue class
Nvx = zeros(1,3);
for jj=1:3
    Nvx(jj) = size(v_use{jj},1);
end

% % Find min-max, mean, median, std, skewness ,kurtosis for each MPM
% mM = zeros(4,2); mM(:,1) = Inf; mM(:,2) = -Inf;
% meanVal = zeros(3,4); medVal = zeros(3,4); stdVal = zeros(3,4);
% skewVal = zeros(3,4); kurtVal = zeros(3,4);
% for ii=1:3 % tc
%     for jj=1:4 % mpm
%         % min/max total
%         tmp_m = min(v_use{ii}(:,jj));
%         if tmp_m<mM(jj,1), mM(jj,1) = tmp_m; end
%         tmp_M = max(v_use{ii}(:,jj));
%         if tmp_M>mM(jj,2), mM(jj,2) = tmp_M; end
%         % mean/meadian/std/skewness/kurtosis
%         meanVal(ii,jj) = mean(v_use{ii}(:,jj));
%         medVal(ii,jj) = median(v_use{ii}(:,jj));
%         stdVal(ii,jj) = std(v_use{ii}(:,jj));
%         skewVal(ii,jj) = skewness(v_use{ii}(:,jj));
%         kurtVal(ii,jj) = kurtosis(v_use{ii}(:,jj))-3;
%     end
% end
mM = res.vMPMstats.mM;
Nbins = 40;
bin = cell(1,4);
for ii=1:4
    bin{ii} = mM(ii,1):(mM(ii,2)-mM(ii,1))/Nbins:mM(ii,2);
end

lab_MPM = {'MT', 'A', 'R1', 'R2'};
lab_tc = {'GM', 'WM', 'Lesion'};

% Plot histograms from all 4 MPMs
N = zeros(Nbins+1,3,4);
figure,
for ii=1:4
    subplot(2,2,ii)
    for jj=1:3
        N(:,jj,ii) = histc(v_use{jj}(:,ii),bin{ii})/Nvx(jj);
    end
    bar(bin{ii},squeeze(N(:,:,ii)),'histc')
    legend(lab_tc)
    xlabel(lab_MPM{ii});
end
figure,
for ii=1:4
    subplot(2,2,ii)
    plot(bin{ii},cumsum(squeeze(N(:,:,ii))))
    ax = axis; ax([3 4]) = [0 1]; axis(ax);
    legend(lab_tc)
    xlabel(lab_MPM{ii});
end

% Plot joint distribution of MPMs for GM/WM/Lesion
figure,
f_ind = 0;
for ii=1:3
    for jj=(ii+1):4
        f_ind = f_ind+1;
        subplot(3,2,f_ind), hold on
        plot(v_use{1}(:,ii),v_use{1}(:,jj),'.b', ...
            v_use{2}(:,ii),v_use{2}(:,jj),'.g', ...
            v_use{3}(:,ii),v_use{3}(:,jj),'.r', ...
            'MarkerSize',.1)
        %         if f_ind==1, legend(lab_tc,'Location','EastOutside'); end
        xlabel(lab_MPM{ii});
        ylabel(lab_MPM{jj});
        axis square
        axis([mM(ii,:) mM(jj,:)])
    end
end
axes('position',[0.1 0.1 .8 .8])
plot(rand(2,1),rand(2,1),'.b',rand(2,1),rand(2,1),'.g',rand(2,1),rand(2,1),'.r','Visible','off')
set(gca,'Visible','off')
legend(lab_tc,'Location','North')

% Create joint histograms and plot them
% http://www.mathworks.com/matlabcentral/answers/101346-how-do-i-use-multiple-colormaps-in-a-single-figure
% MT vs A
descriptor = [mM(1,1), mM(1,2), 256 ; mM(2,1), mM(2,2), 256];
[jhist, descriptor]  = histogram2(v_use{1}(:,1)',v_use{1}(:,2)',descriptor);
jhist_GM = jhist/numel(v_use{1}(:,1));
x_axes = (0:descriptor(1,3)-1)/descriptor(1,3)*(descriptor(1,2)-descriptor(1,1))+descriptor(1,1);
y_axes = (0:descriptor(2,3)-1)/descriptor(2,3)*(descriptor(2,2)-descriptor(2,1))+descriptor(2,1);
% figure, imagesc(jhist), axis xy, axis image
figure, hold on
hs_GM = surf(x_axes,y_axes,jhist_GM);
XT = get(gca,'XTick'); YT = get(gca,'YTick');
set(hs_GM,'LineStyle','none')
colorbar
% set(gca,'XTickLabel',x_axes(XT))
% set(gca,'YTickLabel',y_axes(YT))

Nc = 64; ratio = 2;
hsv_b = zeros(Nc,3);
hsv_b(:,1) = .6;
hsv_b(:,2) = [(0:Nc/2-1)/(Nc/2) (Nc/2:-1:1)/(Nc/2)];
hsv_b(:,3) = (0:Nc-1)/(Nc-1)/ratio+(1-1/ratio);
cmap_b = hsv2rgb(hsv_b);
% colormap(cmap_b)

[jhist, descriptor]  = histogram2(v_use{2}(:,1)',v_use{2}(:,2)',descriptor);
jhist_WM = jhist/numel(v_use{2}(:,1));
hs_WM = surf(x_axes,y_axes,jhist_WM);
set(hs_WM,'LineStyle','none')
hsv_g = zeros(Nc,3);
hsv_g(:,1) = .3;
hsv_g(:,2) = [(0:Nc/2-1)/(Nc/2) (Nc/2:-1:1)/(Nc/2)];
hsv_g(:,3) = (0:Nc-1)/(Nc-1)/ratio+(1-1/ratio);
cmap_g = hsv2rgb(hsv_g);
% colormap(cmap_g)

[jhist, descriptor]  = histogram2(v_use{3}(:,1)',v_use{3}(:,2)',descriptor);
jhist_LE = jhist/numel(v_use{3}(:,1));
hs_LE = surf(x_axes,y_axes,jhist_LE);
set(hs_LE,'LineStyle','none')
hsv_r = zeros(Nc,3);
hsv_r(:,1) = .9;
hsv_r(:,2) = [(0:Nc/2-1)/(Nc/2) (Nc/2:-1:1)/(Nc/2)];
hsv_r(:,3) = (0:Nc-1)/(Nc-1)/ratio+(1-1/ratio);
cmap_r = hsv2rgb(hsv_r);

cmap_bgr = [cmap_b ; cmap_g ; cmap_r];
cd_GM = get(hs_GM,'CData');
cd_WM = get(hs_WM,'CData');
cd_LE = get(hs_LE,'CData');
min_cd = min([cd_GM(:) ; cd_WM(:) ; cd_LE(:)]);
max_cd = max([cd_GM(:) ; cd_WM(:) ; cd_LE(:)]);
cd_GM = min(Nc,round((Nc-1)*(cd_GM-min_cd)/(max_cd-min_cd))+1);
cd_WM = min(Nc,round((Nc-1)*(cd_WM-min_cd)/(max_cd-min_cd))+1)+Nc;
cd_LE = min(Nc,round((Nc-1)*(cd_LE-min_cd)/(max_cd-min_cd))+1)+2+Nc;
set(hs_GM,'CData',cd_GM)
set(hs_WM,'CData',cd_WM)
set(hs_LE,'CData',cd_LE)
caxis([0 3*Nc])
colormap(cmap_bgr)








%% Bits and pieces

% [jhist,xy_axes] = crc_jhist(v_use{1}(:,1),v_use{1}(:,2),mM([1 2],:));
% figure, imagesc(jhist), axis xy, axis image
% XT = get(gca,'XTick'); YT = get(gca,'YTick');
% set(gca,'XTickLabel',xy_axes(XT,1))
% set(gca,'YTickLabel',xy_axes(YT,2))




% % Plot in 3D for GM/WM/les
% jj=1;
% figure,
% plot3(v_use{1}(:,jj),v_use{2}(:,jj),v_use{3}(:,jj),'.')

% 
% 
% descriptor = [mM(1,1), mM(1,2), 256 ; mM(2,1), mM(2,2), 256];
% [jhist, descriptor]  = histogram2(v_use{1}(:,1)',v_use{1}(:,2)',descriptor);
% jhist_GM = jhist/numel(v_use{1}(:,1));
% x_axes = (0:descriptor(1,3)-1)/descriptor(1,3)*(descriptor(1,2)-descriptor(1,1))+descriptor(1,1);
% y_axes = (0:descriptor(2,3)-1)/descriptor(2,3)*(descriptor(2,2)-descriptor(2,1))+descriptor(2,1);
% % figure, imagesc(jhist), axis xy, axis image
% figure, hold on
% hs_GM = surf(x_axes,y_axes,jhist_GM);
% XT = get(gca,'XTick'); YT = get(gca,'YTick');
% set(hs_GM,'LineStyle','none')
% colorbar
% % set(gca,'XTickLabel',x_axes(XT))
% % set(gca,'YTickLabel',y_axes(YT))
% 
% Nc = 64; ratio = 2;
% hsv_b = zeros(Nc,3);
% hsv_b(:,1) = .6;
% hsv_b(:,2) = [(0:Nc/2-1)/(Nc/2) (Nc/2:-1:1)/(Nc/2)];
% hsv_b(:,3) = (0:Nc-1)/(Nc-1)/ratio+(1-1/ratio);
% cmap_b = hsv2rgb(hsv_b);
% % colormap(cmap_b)
% 
% [jhist, descriptor]  = histogram2(v_use{2}(:,1)',v_use{2}(:,2)',descriptor);
% jhist_WM = jhist/numel(v_use{2}(:,1));
% hs_WM = surf(x_axes,y_axes,jhist_WM);
% set(hs_WM,'LineStyle','none')
% hsv_g = zeros(Nc,3);
% hsv_g(:,1) = .3;
% hsv_g(:,2) = [(0:Nc/2-1)/(Nc/2) (Nc/2:-1:1)/(Nc/2)];
% hsv_g(:,3) = (0:Nc-1)/(Nc-1)/ratio+(1-1/ratio);
% cmap_g = hsv2rgb(hsv_g);
% % colormap(cmap_g)
% 
% cmap_bg = [cmap_b ; cmap_g];
% colormap(cmap_bg)
% cd_GM = get(hs_GM,'CData');
% cd_WM = get(hs_WM,'CData');
% min_cd = min([cd_GM(:) ; cd_WM(:)]);
% max_cd = max([cd_GM(:) ; cd_WM(:)]);
% cd_GM = min(Nc,round((Nc-1)*(cd_GM-min_cd)/(max_cd-min_cd))+1);
% cd_WM = min(Nc,round((Nc-1)*(cd_WM-min_cd)/(max_cd-min_cd))+1)+Nc;
% set(hs_GM,'CData',cd_GM)
% set(hs_WM,'CData',cd_WM)
% caxis([0 2*Nc])

