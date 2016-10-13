function [Tumour_volume,T1_volume]=getBRAT_volumes(current, BRAT_dir)

% routine to extract volumes from BRAT data
% folders are organized as it (from Download) but VSD.mha were transformed
% to VSD.nii -- in addition the augmented segmentation was run and tissue
% classes 1 to 4 are further used to remove the remaining non brain tissue
%
% INPUTS current is the info dir used for analyses
%        BRAT_dir is the BRAT root directory
%        see data_analysis.m
%
% OUTPUTS are the volumes (in mm^3) 

cd(BRAT_dir);
cd('HGG'); BRAT_dir1 = pwd;
cd ..; cd('LGG'); BRAT_dir2 = pwd;
        
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


        % get the images
        % -------------------------------------------------------
         cd(patient_dir); local = dir; 
                
        clear vols
        for d=3:size(local,1)
            if ~isempty(findstr(local(d).name,'more.'))
                vols{6} = [patient_dir local(d).name filesep 'VSD.nii'];
            elseif ~isempty(findstr(local(d).name,'T1.'))
                vols{1} = [patient_dir local(d).name filesep 'VSD.nii'];
                vols{2} = [patient_dir local(d).name filesep 'augmented_segmentation' filesep 'c1VSD.nii'];
                vols{3} = [patient_dir local(d).name filesep 'augmented_segmentation' filesep 'c2VSD.nii'];
                vols{4} = [patient_dir local(d).name filesep 'augmented_segmentation' filesep 'c3VSD.nii'];
                vols{5} = [patient_dir local(d).name filesep 'augmented_segmentation' filesep 'c4VSD.nii'];
            end
        end
        
        % compute the mask, brain volume, tumour volume
        % -----------------------------------------------
        V = spm_vol(vols);
        mask = spm_read_vols(V{2})+spm_read_vols(V{3})+spm_read_vols(V{4})+spm_read_vols(V{5});
        T1 = spm_read_vols(V{1}).*mask;
        mm3 = diag(V{1}.mat); mm3 = prod(abs(mm3(1:3)));
        T1_volume(index) = sum(T1(:)~=0)*mm3;
        Tumour =  spm_read_vols(V{6}).*mask;
        mm3 = diag(V{6}.mat); mm3 = prod(abs(mm3(1:3)));
        Tumour_volume(index) = sum(Tumour(:)~=0)*mm3;
        index = index +1;
    end
end
       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
