
%% Play around with voxel intensities in GM/WM/lesion

load('D:\ccc_DATA\MS_ELommers\Data2_unwTPM\ExP_s7480-0005-00001-000176-00_MT_EPIB1.mat')

% Number of voxels per tissue class
Nvx = zeros(1,3);
for jj=1:3
    Nvx(jj) = size(res.vMPM{jj},1);
end

% Find min-max, mean, median, std for each MPM
mM = zeros(4,2); mM(:,1) = Inf; mM(:,2) = -Inf;
meanVal = zeros(3,4); medVal = zeros(3,4); stdVal = zeros(3,4);
for ii=1:3 % tc
    for jj=1:4 % mpm
        % min/max total
        tmp_m = min(res.vMPM{ii}(:,jj));
        if tmp_m<mM(jj,1), mM(jj,1) = tmp_m; end
        tmp_M = max(res.vMPM{ii}(:,jj));
        if tmp_M>mM(jj,2), mM(jj,2) = tmp_M; end
        % mean/meadian/std
        meanVal(ii,jj) = mean(res.vMPM{ii}(:,jj));
        medVal(ii,jj) = median(res.vMPM{ii}(:,jj));
        stdVal(ii,jj) = std(res.vMPM{ii}(:,jj));
    end
end
Nbins = 50;
bin = cell(1,4);
for ii=1:4
    bin{ii} = mM(ii,1):(mM(ii,2)-mM(ii,1))/Nbins:mM(ii,2);
end

% Plot histograms from all 4 MPMs
lab_MPM = {'MT', 'A', 'R1', 'R2'};
lab_tc = {'GM', 'WM', 'lesion'};
N = zeros(Nbins+1,4);
figure,
for ii=1:4
subplot(2,2,ii) % MT
    for jj=1:3
        N(:,jj) = histc(res.vMPM{jj}(:,ii),bin{ii})/Nvx(jj);
        bar(bin{ii},N,'histc')
        legend(lab_tc)
        xlabel(lab_MPM{ii});
    end
end


% % Plot in 3D for GM/WM/les
% jj=1;
% figure,
% plot3(res.vMPM{1}(:,jj),res.vMPM{2}(:,jj),res.vMPM{3}(:,jj),'.')




