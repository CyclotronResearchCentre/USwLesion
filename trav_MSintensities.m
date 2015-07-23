
%% Play around with voxel intensities in GM/WM/lesion

% Refs:
% https://en.wikipedia.org/wiki/Skewness
% https://en.wikipedia.org/wiki/Kurtosis
% https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test

% Ideas:
% - compare the distribution of intensities with Two-sample Kolmogorov-
% Smirnov test. Function in Matlab is "kstest2"

load('D:\ccc_DATA\MS_ELommers\Data2_unwTPM\ExP_s7480-0005-00001-000176-00_MT_EPIB1.mat')

% Number of voxels per tissue class
Nvx = zeros(1,3);
for jj=1:3
    Nvx(jj) = size(res.vMPM{jj},1);
end

% Find min-max, mean, median, std, skewness ,kurtosis for each MPM
mM = zeros(4,2); mM(:,1) = Inf; mM(:,2) = -Inf;
meanVal = zeros(3,4); medVal = zeros(3,4); stdVal = zeros(3,4);
skewVal = zeros(3,4); kurtVal = zeros(3,4);
for ii=1:3 % tc
    for jj=1:4 % mpm
        % min/max total
        tmp_m = min(res.vMPM{ii}(:,jj));
        if tmp_m<mM(jj,1), mM(jj,1) = tmp_m; end
        tmp_M = max(res.vMPM{ii}(:,jj));
        if tmp_M>mM(jj,2), mM(jj,2) = tmp_M; end
        % mean/meadian/std/skewness/kurtosis
        meanVal(ii,jj) = mean(res.vMPM{ii}(:,jj));
        medVal(ii,jj) = median(res.vMPM{ii}(:,jj));
        stdVal(ii,jj) = std(res.vMPM{ii}(:,jj));
        skewVal(ii,jj) = skewness(res.vMPM{ii}(:,jj));
        kurtVal(ii,jj) = kurtosis(res.vMPM{ii}(:,jj))-3;
    end
end
Nbins = 40;
bin = cell(1,4);
for ii=1:4
    bin{ii} = mM(ii,1):(mM(ii,2)-mM(ii,1))/Nbins:mM(ii,2);
end

lab_MPM = {'MT', 'A', 'R1', 'R2'};
lab_tc = {'GM', 'WM', 'lesion'};

% Plot histograms from all 4 MPMs
N = zeros(Nbins+1,3,4);
figure,
for ii=1:4
    subplot(2,2,ii)
    for jj=1:3
        N(:,jj,ii) = histc(res.vMPM{jj}(:,ii),bin{ii})/Nvx(jj);
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

% % Plot in 3D for GM/WM/les
% jj=1;
% figure,
% plot3(res.vMPM{1}(:,jj),res.vMPM{2}(:,jj),res.vMPM{3}(:,jj),'.')




