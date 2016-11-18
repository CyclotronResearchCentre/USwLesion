% Normalization analysis
%
% This scripts is analyzing normalization procedures for patient data but
% using healthy participants
% 1. From the healthy brain, normalize == truth
% 2. reverse normalize a lesion/tumour already in standard space toward the
% healthy brain
% 3. create a 'patient brain' by averaging the healthy tissue (1/4) with
% the lesion (3/4) and renormalize using bith the standard approach and the
% augmented asegmentation approache == normalized patient
% 4 compare the 'normalized patient' to the 'truth' ' it will tell how much
% the normalization procedure varies once a lesion is present
% --> A. use structural similarity index (Wang et al 2004 IEEE Trans Image
%       Process 13, 600-612) for both the healthy tissue and for the lesion
% --> B. the standardized squared differences 
%
% here we use the BRAT data and tumours segmented and normalized using the
% USW_lesion toolbox and healthy brains from the 1000 connectome project
%
%% set the directory and defaults
% ------------------------------
current = pwd;
% global defaults
% defaults = spm('defaults','FMRI');
% spm_jobman('initcfg')
BRAT_dir = '/media/cpernet/SATA/BRAT/BRATS2015_Training';
% BRAT_dir = uigetdir(pwd,'Select BRAT directory');
% if BRAT_dir == 0
%     return
% else
%     cd(BRAT_dir);
%     try
%         cd('HGG'); BRAT_dir1 = pwd;
%         cd ..
%         cd('LGG'); BRAT_dir2 = pwd;
%     catch
%         error('HGG directory not found');
%     end
% end

Healthy_dir = '/media/cpernet/SATA/BRAT/F1000';
% Healthy_dir = uigetdir(pwd,'Select F1000 directory');

%% run the normalization on tumours
clear matlabbatch
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
index = 1;
for tumour_type = 1:2
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
        cd(patient_dir); local = dir;
         for d=3:size(local,1)
            if ~isempty(findstr(local(d).name,'T1.'))
                cd(local(d).name);
                matlabbatch{1}.spm.spatial.normalise.write.subj.def{1} = [pwd filesep 'augmented_segmentation' filesep 'y_VSD.nii'];
                matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1} = [pwd filesep 'augmented_segmentation' filesep 'lesion_mask.nii'];
                patient_T1{index} = [pwd filesep 'augmented_segmentation' filesep 'wVSD.nii'];
            end
         end
        lesion_masks{index} = spm_jobman('run', matlabbatch);        
        index = index+1; cd(BRAT_dir)
    end
end


%% run the normalization on healthy brains
cd(current); load('batch_standard');
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 2; % here we have skull 
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1]; % write forward and backward
matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
cd(Healthy_dir); local = dir; today = date;
for subjects= 1:30
    cd(local(subjects+2).name)
    matlabbatch{1}.spm.spatial.preproc.channel.vols = [];
    matlabbatch{1}.spm.spatial.preproc.channel.vols{1} = [pwd filesep 'anat.nii'];
    vols{1} = [pwd filesep 'anat.nii']; vols{2} = [pwd filesep 'c1anat.nii'];
    vols{3} = [pwd filesep 'c2anat.nii']; vols{4} = [pwd filesep 'c3anat.nii'];
    matlabbatch{2}.spm.spatial.normalise.write.subj.resample = vols';
    out = spm_jobman('run', matlabbatch);
    
    % cleanup
    check_data = dir; mkdir('standard_segmentation');
    try
        movefile(cell2mat(out{1}.param),[pwd filesep 'standard_segmentation'])
        movefile(cell2mat(out{1}.invdef),[pwd filesep 'standard_segmentation'])
        movefile(cell2mat(out{1}.fordef),[pwd filesep 'standard_segmentation'])
        for d=1:5
            movefile(cell2mat(out{1}.tiss(d).c),[pwd filesep 'standard_segmentation'])
        end
        for d=1:4
            movefile(cell2mat(out{2}.files(d)),[pwd filesep 'standard_segmentation'])
        end
    end
    cd ..
end

%% create new patients
cd(Healthy_dir); local = dir;
for subjects= 1:30
    cd(local(subjects+2).name) ; 
    % load tumour mask from a patient in standard space
    name = spm_vol(lesion_masks{subjects}{1}.files);
    tumour_mask = spm_read_vols(name{1});
    % use real T1 values in that mask
    T = spm_vol(patient_T1{subjects});
    tumour = spm_read_vols(T).*tumour_mask;
    T.fname = [pwd filesep 'tumour.nii'];
    spm_write_vol(T,tumour);
    % mix that tumour with the T1 of a healthy subject
    V = spm_vol([pwd filesep 'standard_segmentation' filesep 'wanat.nii']); 
    anat =  spm_read_vols(V);
    A = anat(find(tumour)); AM = max(A);
    B = tumour(find(tumour)); BM = max(B);
    if BM > AM
        scaling = BM/AM;
        anat(find(tumour)) = ((A.*0.25)+((B./scaling).*0.75)).*1.3;
   else
        scaling = AM/BM;
        anat(find(tumour)) = (((A./scaling)*0.25)+(B*0.75)).*1.3;
    end
    V.fname = [pwd filesep 'damaged_T1.nii'];
    spm_write_vol(V,anat); 
    
    % inverse normalize to be back in subject space
    load([current filesep 'normalize_write']);
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pwd filesep 'standard_segmentation' filesep 'iy_anat.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[pwd filesep 'damaged_T1.nii']};
    spm_jobman('run', matlabbatch);

    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[pwd filesep 'tumour.nii']};
    spm_jobman('run', matlabbatch); cd(Healthy_dir); 
end

%% run the normalization (twice) on the new patients
cd(Healthy_dir); local = dir; today = date;
for subjects= 1:30
    cd(local(subjects+2).name)
    load([current filesep 'batch_standard.mat']);
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 0;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.channel.vols = [];
    matlabbatch{1}.spm.spatial.preproc.channel.vols{1} = [pwd filesep 'wdamaged_T1.nii'];
    matlabbatch{2}.spm.spatial.normalise.write.subj.resample{1} = [pwd filesep 'wdamaged_T1.nii'];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    out = spm_jobman('run', matlabbatch);
    
    % cleanup
    mkdir('segmentation_damagedT1')
    movefile(cell2mat(out{1}.param),[pwd filesep 'segmentation_damagedT1'])
    movefile(cell2mat(out{1}.fordef),[pwd filesep 'segmentation_damagedT1'])
    for d=1:5        
            movefile(cell2mat(out{1}.tiss(d).c),[pwd filesep 'segmentation_damagedT1'])
    end
    movefile(cell2mat(out{2}.files),[pwd filesep 'segmentation_damagedT1'])
    
    % same using augmented segmentation
    clear matlabbatch; load([current filesep 'normalize_fake_patients.mat']); today = date;
    matlabbatch{1}.spm.tools.USwLtools.uswl.options.NbGaussian = [3 2 3 2 1 0 1];
    matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk = {[pwd filesep 'wtumour.nii']};
    matlabbatch{1}.spm.tools.USwLtools.uswl.imgRef = {[pwd filesep 'wdamaged_T1.nii']};
    spm_jobman('run', matlabbatch);
    cd ..
end


%% compare results
cd(Healthy_dir); local = dir;
for subjects= 1:30
    cd(local(subjects+2).name)
    % load the reference
    Ref = spm_read_vols(spm_vol([pwd filesep 'standard_segmentation' filesep 'wanat.nii']));
    Ref(isnan(Ref)) = 0;
    % load the standard normalization
    Standard = spm_read_vols(spm_vol([pwd filesep 'segmentation_damagedT1' filesep 'wwdamaged_T1.nii']));
    Standard(isnan(Standard)) = 0;
    % load the augmented seg/normalization
    Augmented = spm_read_vols(spm_vol([pwd filesep 'wwdamaged_T1.nii']));
    Augmented(isnan(Augmented)) = 0;
    
    % compute the mean square distance in the healthy tissue only
    % icv from initial anatomical and remove the tumour (all in MNI space)
    icv = spm_read_vols(spm_vol([pwd filesep 'standard_segmentation' filesep 'wc1anat.nii'])) + ...
        spm_read_vols(spm_vol([pwd filesep 'standard_segmentation' filesep 'wc2anat.nii'])) + ...
        spm_read_vols(spm_vol([pwd filesep 'standard_segmentation' filesep 'wc3anat.nii']));
    mask = smooth3((icv>0),'box',[3 3 3]); mask = mask>0; 
    
    name = spm_vol(lesion_masks{subjects}{1}.files);
    tumour_mask = logical(spm_read_vols(name{1}));  
    mask(tumour_mask) = 0; % this is brain mask minus tumour
    figure; for z=1:size(icv,3); imagesc(squeeze(mask(:,:,z)).*5+squeeze(tumour_mask(:,:,z))); pause(0.1); end; close
    
    % similarity 
    C = cov(Ref(mask),Standard(mask));
    sim(subjects,1) = ((2*nanmean(Ref(mask))*nanmean(Standard(mask))+100)*(2*C(1,2)+10)) / ...
        (nanmean(Ref(mask))^2+nanmean(Standard(mask))^2+100)*(nanvar(Ref(mask))+nanvar(Standard(mask))+10);
    
    C = cov(Ref(mask),Augmented(mask));    
    sim(subjects,2) = ((2*nanmean(Ref(mask))*nanmean(Augmented(mask))+100)*(2*C(1,2)+10)) / ...
        (nanmean(Ref(mask))^2+nanmean(Augmented(mask))^2+100)*(nanvar(Ref(mask))+nanvar(Augmented(mask))+10);
    
    
    C = cov(Ref(tumour_mask),Standard(tumour_mask));
    sim(subjects,3) = ((2*nanmean(Ref(tumour_mask))*nanmean(Standard(tumour_mask))+100)*(2*C(1,2)+10)) / ...
        (nanmean(Ref(tumour_mask))^2+nanmean(Standard(tumour_mask))^2+100)*(nanvar(Ref(tumour_mask))+nanvar(Standard(tumour_mask))+10);

    C = cov(Ref(tumour_mask),Augmented(tumour_mask));    
    sim(subjects,4) = ((2*nanmean(Ref(tumour_mask))*nanmean(Augmented(tumour_mask))+100)*(2*C(1,2)+10)) / ...
        (nanmean(Ref(tumour_mask))^2+nanmean(Augmented(tumour_mask))^2+100)*(nanvar(Ref(tumour_mask))+nanvar(Augmented(tumour_mask))+10);

    
     % intensity based differences
    Ref = (Ref.*100) ./ max(Ref(:));
    Standard = (Standard.*100) ./ max(Standard(:));
    Augmented = (Augmented.*100) ./ max(Augmented(:));
    distance(subjects,1) = nanmean((Ref(mask) - Standard(mask)).^2);
    distance(subjects,2) = nanmean((Ref(mask) - Augmented(mask)).^2);
    distance(subjects,3) = nanmean((Ref(tumour_mask) - Standard(tumour_mask)).^2);
    distance(subjects,4) = nanmean((Ref(tumour_mask) - Augmented(tumour_mask)).^2);
   
    
    cd ..
end

cd(current); 
save normalization_results

%% WARNING SUBJECT 19 EXCLUDED 
%% HAVEN'T FIGURED WHY BUT VALUES ARE OUTLYING MASSIVELY 
sim(19,:) = []; distance(19,:) = [];
[med_sim, CI_sim]= rst_data_plot(sim(:,[1 2]),'estimator','median');
[med_sim([3 4]), CI_sim(:,[3 4])]= rst_data_plot(sim(:,[3 4]),'estimator','median');
[hsim,CIsim,psim] = rst_pttest(sim(:,[2 4]),sim(:,[1 3]),'median',1,0.05,10000)

[med_dist, CI_distance]= rst_data_plot(distance(:,[1 2]),'estimator','median');
[med_dist([3 4]), CI_distance(:,[3 4])]= rst_data_plot(distance(:,[3 4]),'estimator','median');
[hdist,CIdist,pdist] = rst_pttest(distance(:,[2 4]),distance(:,[1 3]),'median',1,0.05,10000)
    
    
    
    
    

