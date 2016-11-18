% Script use to segment and analyze BRAT images
% All images have been converted to .nii format and are called VSD.nii
% the mask generated from the FLAIR is called voi.nii
%
% Segmentation
% For each subject, we load the SPM batch and run the segmentation on the
% T1, with the FLAIR added as second image. Since the BRAT image have no
% skull, the nb of gaussian is modelled as [3  2  3  2  1  0  1], i.e. 3
% gaussian for the tumour and oedema and 0 for skull. Finally the c3 is
% thresolded keeping the biggest cluster only (if more than one) and is
% used as a reference against the ground truth keeping voxels with higher
% prob than other tossie classes
%
% Data analysis
% The mask used for the segmentation and the resulting mask are both
% compared to the ground truth, this allows to estimate how much the
% segmentation works, and how much this is infludenced by the initial mask
%
% calls routines from https://github.com/CPernet/spmup to resize the template
%                from https://github.com/CPernet/Robust_Statistical_Toolbox for boxplots
%                from https://sourceforge.net/projects/robustcorrtool/ for Spearman corr


%% set the directory and defaults
% ------------------------------
current = pwd;
global defaults
defaults = spm('defaults','FMRI');
spm_jobman('initcfg')
BRAT_dir = uigetdir(pwd,'Select BRAT directory');
if BRAT_dir == 0
    return
else
    cd(BRAT_dir);
    try
        cd('HGG'); BRAT_dir1 = pwd;
        cd ..
        cd('LGG'); BRAT_dir2 = pwd;
    catch
        error('HGG directory not found');
    end
end

%% per patient, run the segmentation and compute similarity
% -----------------------------------------------------------

index = 1;
for tumour_type = 1:2
    % loop in High Grade Glioma or Low Grade Glioma
    % --------------------------------------------
    BRAT_dir = eval(['BRAT_dir' num2str(tumour_type)]);
    cd(BRAT_dir); folders = dir;
    if tumour_type == 1
        folders = folders(3:22);
    else
        folders = folders(3:12);
    end
    
    for patient = 1:size(folders,1)
        % loop for each patient
        % ---------------------
        patient_dir = [BRAT_dir filesep folders(patient).name filesep];
        cd(patient_dir); local = dir; today = date;
                
        % 1 segment/normalize image using the unified segmentation
        %---------------------------------------------------------
        clear matlabbatch; load([current filesep 'batch_standard.mat'])
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 0;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 1;
        clear vols
        for d=3:size(local,1)
            if ~isempty(findstr(local(d).name,'Flair.'))
                vols{2} = [patient_dir local(d).name filesep 'VSD.nii'];
            elseif ~isempty(findstr(local(d).name,'T1.'))
                vols{1} = [patient_dir local(d).name filesep 'VSD.nii'];
                normalized_standard{index} = {[patient_dir local(d).name filesep 'standard_segmentation' filesep 'wVSD.nii']};
            end
        end
        
        matlabbatch{1}.spm.spatial.preproc.channel.vols = vols';
        matlabbatch{2}.spm.spatial.normalise.write.subj.resample = {vols{1}};
        spm_jobman('run', matlabbatch);
        
        %cleanup
        cd(fileparts(vols{1}))
        local = dir; mkdir('standard_segmentation');
        for d=3:size(local,1)
            if strncmp(local(d).date,today,11) && local(d).isdir == 0
                movefile(local(d).name,[pwd filesep 'standard_segmentation'])
            end
        end
        
        % 2 segment/normalize image using the augmented unified segmentation
        %--------------------------------------------------------------------
        clear matlabbatch; load([current filesep 'batch_BT.mat'])
        % tissue class update
        matlabbatch{1}.spm.tools.USwLtools.uswl.options.NbGaussian = [3 2 3 2 1 0 1];
        cd(patient_dir); local = dir; today = date;
        for d=3:size(local,1)
            if ~isempty(findstr(local(d).name,'Flair.'))
                matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk = {[patient_dir local(d).name filesep 'voi.nii']};
                matlabbatch{1}.spm.tools.USwLtools.uswl.imgOth = {[patient_dir local(d).name filesep 'VSD.nii']};
                opt.mask = [patient_dir local(d).name filesep 'icv_kVSD.nii'];
            elseif ~isempty(findstr(local(d).name,'T1.'))
                matlabbatch{1}.spm.tools.USwLtools.uswl.imgRef = {[patient_dir local(d).name filesep 'VSD.nii']};
                tumour_mask = [patient_dir local(d).name filesep 'cleaned_c3VSD.nii'];
            elseif ~isempty(findstr(local(d).name,'more'))
                ground_truth = [patient_dir local(d).name filesep 'VSD.nii'];
            end
        end
        
        % patients with 2 tumours 
        if strcmp(folders(patient).name,'brats_2013_pat0005_1') || strcmp(folders(patient).name,'brats_2013_pat0026_1')
            matlabbatch{1}.spm.tools.USwLtools.uswl.options.thrLesion = 10;
        end
        spm_jobman('run', matlabbatch);
        
        
        % cleanup
        for folder = 1:2
            if folder == 1
                cd(fileparts(opt.mask)) % go into FLAIR dir
            else
                cd(fileparts(tumour_mask)) % go into T1 dir
                normalized_augmented{index} = {[pwd filesep 'augmented_segmentation' filesep 'wVSD.nii']};
            end
            
            local = dir;
            mkdir('augmented_segmentation')
            for d=3:size(local,1)
                if strncmp(local(d).date,today,11) && local(d).isdir == 0
                    movefile(local(d).name,[pwd filesep 'augmented_segmentation'])
                end
            end
        end
        
        % 3 analyze the masks from the augmented segmentation
        % ---------------------------------------------------
        % generate a lesion mask as (C3>C1)U(C3>C2)U(C3>C4) ie lesion prob > gray or white or csf
        class_path = fileparts(tumour_mask);
        C1 = [class_path filesep 'augmented_segmentation' filesep 'c1VSD.nii'];
        C2 = [class_path filesep 'augmented_segmentation' filesep 'c2VSD.nii'];
        C3 = spm_vol([class_path filesep 'augmented_segmentation' filesep 'cleaned_c3VSD.nii']);
        tmp = spm_read_vols(C3); 
        percentage_bigger_then_ninetynine_in_C3(index) = sum(tmp(:)>0.99) / sum(tmp(:)~=0);
        C4 = [class_path filesep 'augmented_segmentation' filesep 'c4VSD.nii'];
        
        % analyze
        opt.thr = 0; opt.mask = [opt.mask(1:end-12) 'augmented_segmentation' filesep 'icv_kVSD.nii'];
        target = spm_read_vols(spm_vol(ground_truth));
        
        % overlap rought initial mask / ground truth
        [mJ_before(index),mHD_before(index),overlap_before(index)] = image_overlap(cell2mat(matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk),target,opt);
        
        % overlap deciles of tissue class / ground truth
        parfor d=1:10
            source = tmp.*(tmp>(d/10)-0.1);
            [mJ(d),mHd(d),overlap(d)] = image_overlap(source,target,opt);
        end
        mJ_cont(index,:) = mJ; mHD_cont(index,:) = mHd;
        overlap_cont(index,:) = overlap.voxel;
        
        % overlap tissue class bigger than others / ground truth
         M = (spm_read_vols(C3) > spm_read_vols(spm_vol(C1))) .* (spm_read_vols(C3) > spm_read_vols(spm_vol(C2))) ...
            .* (spm_read_vols(C3) > spm_read_vols(spm_vol(C4)));
        C3.fname = [class_path filesep 'augmented_segmentation' filesep 'lesion_mask.nii'];
        C3.descrip = [C3.descrip 'and bigger than other prob'];
        lesion_mask = spm_write_vol(C3,M);      
        [mJ_after(index) ,mHD_after(index) ,overlap_after(index)] = image_overlap(lesion_mask.fname,target,opt);
        clear C1 C2 C3 C4 M source tmp target
        index= index+1; cd(BRAT_dir);     
                
    end
end
cd(current)
save tmp_results

%% check the results
% similarity to ground truth
for i=1:30
    Atp(i) = overlap_after(i).voxel.tp;
    Atn(i) = overlap_after(i).voxel.tn;
    Afp(i) = overlap_after(i).voxel.fp;
    Afn(i) = overlap_after(i).voxel.fn;
    Amcc(i) = overlap_after(i).voxel.mcc;
    ACK(i) = overlap_after(i).voxel.CK;
    Btp(i) = overlap_before(i).voxel.tp;
    Btn(i) = overlap_before(i).voxel.tn;
    Bfp(i) = overlap_before(i).voxel.fp;
    Bfn(i) = overlap_before(i).voxel.fn;
    Bmcc(i) = overlap_before(i).voxel.mcc;
    BCK(i) = overlap_before(i).voxel.CK;
end

[est_CM,HDI_CM]=rst_data_plot([Btp' Atp' Btn' Atn' Bfp' Afp' Bfn' Afn' ],'estimator','median');

[est_J,HDI_J]=rst_data_plot([mJ_before' mJ_after'],'estimator','median');
[est_H,HDI_H]=rst_data_plot([mHD_before' mHD_after'],'estimator','median');
[est_MCC,HDI_MCC]=rst_data_plot([Bmcc' Amcc'],'estimator','median');
[est_CK,HDI_CK]=rst_data_plot([BCK' ACK'],'estimator','median');
[httest,CIttest,pttest] = rst_pttest([mJ_after' mHD_after' Amcc' ACK'],[mJ_before' mHD_before' Bmcc' BCK'],'median',1,0.05,10000)
[r,t,h,outid,hboot,CI] = skipped_correlation([mJ_before' mJ_before' Bmcc' Bmcc'],[mJ_after' mJ_after'-mJ_before' Amcc' Amcc'-Bmcc']);


%% analysis of normalized volumes
% simply compute the mean difference to the template
distance = NaN(30,2);

template = [fileparts(which('spm')) filesep 'canonical' filesep 'avg152T1.nii'];
% resize
bb = [-78  -112   -70; 78    76    86]; vs = [1 1 1];
spmup_resize(template,bb,vs)
template = [fileparts(which('spm')) filesep 'canonical' filesep 'ravg152T1.nii'];
template = spm_read_vols(spm_vol(template));
template = (template.*100) ./ max(template(:));
for patient = 1:30
    img = spm_read_vols(spm_vol(cell2mat(normalized_standard{patient})));
    img = (img.*100) ./ max(img(:));
    dist = (template - img).^2;
    distance(patient,1) = sum(dist(~isnan(dist)))/numel(~isnan(dist));
end

template = [fileparts(which('spm')) filesep 'canonical' filesep 'avg152T1.nii'];
bb = [-78 -112 -70; 78 76 85]; vs = [1 1 1];
spmup_resize(template,bb,vs)
template = [fileparts(which('spm')) filesep 'canonical' filesep 'ravg152T1.nii'];
template = spm_read_vols(spm_vol(template));
template = (template.*100) ./ max(template(:));
for patient = 1:30
    img = spm_read_vols(spm_vol(cell2mat(normalized_augmented{patient})));
    img = (img.*100) ./ max(img(:));
    dist = (template - img).^2;
    distance(patient,2) = sum(dist(~isnan(dist)))/numel(~isnan(dist));
end

[med_dist, CI_distance]= rst_data_plot(distance,'estimator','median');
[hdist,CIdist,pdist] = rst_pttest(distance(:,1),distance(:,2),'median',1,0.05,10000)

% compute the difference between all pairs - take determinant
distance_matrix = zeros(30,30,2);

pairs = nchoosek([1:30],2);
for n=1:size(pairs,1)
    img1 = spm_read_vols(spm_vol(cell2mat(normalized_standard{pairs(n,1)})));
    img2 = spm_read_vols(spm_vol(cell2mat(normalized_standard{pairs(n,2)})));
    img1 = (img1.*100) ./ max(img1(:)); img2 = (img2.*100) ./ max(img2(:));
    dist = (img1 - img2).^2;
    distance_matrix(pairs(n,1),pairs(n,2),1) = sum(dist(~isnan(dist)))/numel(~isnan(dist));
    distance_matrix(pairs(n,2),pairs(n,1),1) = distance_matrix(pairs(n,1),pairs(n,2),1);
    
    img1 = spm_read_vols(spm_vol(cell2mat(normalized_augmented{pairs(n,1)})));
    img2 = spm_read_vols(spm_vol(cell2mat(normalized_augmented{pairs(n,2)})));
    img1 = (img1.*100) ./ max(img1(:)); img2 = (img2.*100) ./ max(img2(:));
    dist = (img1 - img2).^2;
    distance_matrix(pairs(n,1),pairs(n,2),2) = sum(dist(~isnan(dist)))/numel(~isnan(dist));
    distance_matrix(pairs(n,2),pairs(n,1),2) = distance_matrix(pairs(n,1),pairs(n,2),2);
end

A = triu(distance_matrix(:,:,1)); 
B = triu(distance_matrix(:,:,2));
[med_pdist, CI_pdistance] = rst_data_plot([A(find(A)),B(find(B))],'estimator','median');
[hpdist,CIpdist,ppdist] = rst_pttest(A(find(A)),B(find(B)),'median',1,0.05,10000)
[det(distance_matrix(:,:,1)) det(distance_matrix(:,:,2))]
figure; subplot(1,2,1); imagesc(A); 
subplot(1,2,2); imagesc(B); colormap(cubehelixmap('semi_continuous',64));


run Normalization_analysis
