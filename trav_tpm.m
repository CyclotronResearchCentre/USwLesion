
sz = size(tpm_orig)
pl = 80;

imat({log10(squeeze(tpm_orig(:,:,pl,1))) ...
    log10(squeeze(tpm_orig(:,:,pl,2))) ...
    log10(squeeze(tpm_orig(:,:,pl,3))) ...
    })


imat(squeeze(tpm_orig(:,:,pl,2))>=1e-3)

imat({squeeze(tpm_orig(:,:,pl,1)) ...
    squeeze(tpm_orig(:,:,pl,2)) ...
    squeeze(tpm_orig(:,:,pl,3)) ...
    squeeze(tpm_orig(:,:,pl,4)) ...
    squeeze(tpm_orig(:,:,pl,5)) ...
    })

imat({log10(squeeze(tpm_orig(:,:,pl,1))) ...
    log10(squeeze(tpm_orig(:,:,pl,2))) ...
    log10(squeeze(tpm_orig(:,:,pl,3))) ...
    log10(squeeze(tpm_orig(:,:,pl,4))) ...
    log10(squeeze(tpm_orig(:,:,pl,5))) ...
    })

imat({log10(squeeze(tpm_updt(:,:,pl,1))) ...
    log10(squeeze(tpm_updt(:,:,pl,2))) ...
    log10(squeeze(tpm_updt(:,:,pl,3))) ...
    log10(squeeze(tpm_updt(:,:,pl,4))) ...
    log10(squeeze(tpm_updt(:,:,pl,5))) ...
    log10(squeeze(tpm_updt(:,:,pl,6))) ...
    })

imat({ tpm_updt(:,:,pl,6) , log10(tpm_updt(:,:,pl,6)) })

imat({ tpm_ot(:,:,pl) , log10(tpm_ot(:,:,pl)) 
    stpm_ot(:,:,pl) , log10(stpm_ot(:,:,pl))})

for ii=1:6
%     tpm_ii = tpm_updt(:,:,:,ii);
    tpm_ii = tpm_orig(:,:,:,ii);
    nne(ii) = sum(tpm_ii(:)<0);
    nnz(ii) = sum(tpm_ii(:)==0);
end
[nne ; nnz]
tpm_ii = sum(tpm_orig,4);
tpm_ii = sum(tpm_updt,4);
nn1 = sum(tpm_ii(:)~=1)


% % Checkr #0's
% for ii=1:6
%     tpm_ii = squeeze(tpm_orig(:,:,:,ii));
%     sum(tpm_ii(:)<thr)
% end
% % Check tpm_orig, does NOT sum to 1 exactly everywhere
% S_orig = sum(tpm_orig,4);
% sum(S_orig(:)<1)
% sum(S_orig(:)>1)
% [min(S_orig(:)) max(S_orig(:))]

% % Deal with GM, WM, CSF -> set to thr outside head
% for ii = 1:3
%     tpm_ii = squeeze(tpm_orig(:,:,:,ii));
%     tpm_ii(tpm_ii<thr) = thr;
%     tpm_updt(:,:,:,ii) = tpm_ii;
% end
% % -> adjust 'other' outside the brain
% tpm_updt(:,:,:,6) = 1 - sum(tpm_updt(:,:,:,1:5),4);

% Deal with skull & Other > set to thr inside head/brain
for ii = [4 6]
    tpm_ii = squeeze(tpm_orig(:,:,:,ii));
    tpm_ii(tpm_ii<thr) = thr;
    tpm_updt(:,:,:,ii) = tpm_ii;
end
% Smooth a wee bit (logs of) 'other'
% tpm_ot = 1 - sum(tpm_updt(:,:,:,1:5),4); 
tpm_ot = tpm_updt(:,:,:,6) ;
stpm_ot = 10.^(smooth3(log10(tpm_ot),'gaussian'));
tpm_ot(tpm_ot<thr*100) = stpm_ot(tpm_ot<thr*100);
% tpm_updt(:,:,:,6) = stpm_ot;
tpm_updt(:,:,:,6) = tpm_ot;
% % -> adjsut GM/WM/CSF in brain
% EStpm = sum(tpm_updt,4) - 1;
% for ii=1:3
%     tpm_updt(:,:,:,ii) = tpm_updt(:,:,:,ii) - EStpm/3;
% end

% Deal with GM, WM, CSF -> set to thr outside head
for ii = 1:3
    tpm_ii = squeeze(tpm_orig(:,:,:,ii));
    tpm_ii(tpm_ii<thr) = thr;
    tpm_updt(:,:,:,ii) = tpm_ii;
end
% -> adjust 'other' outside the brain
% tpm_updt(:,:,:,6) = 1 - sum(tpm_updt(:,:,:,1:5),4);
tpm_updt(:,:,:,5) = 1 - sum(tpm_updt(:,:,:,[1:4 6]),4);


% Save updated TPMs
Vtpm_u = spm_vol( fullfile(spm('dir'),'tpm','TPM.nii'));
fn_TPM_upd = fullfile(spm('dir'),'tpm','unwTPM_sl2.nii');
for ii=1:6
    Vtpm_u(ii).fname = fn_TPM_upd;
    Vtpm_u(ii) = spm_create_vol(Vtpm_u(ii));
    Vtpm_u(ii) = spm_write_vol(Vtpm_u(ii),tpm_updt(:,:,:,ii));
end

%%%%%%%%%%%

vx = [113 57 7]
vx = [97 73 101]
val_orig = squeeze(tpm_orig(vx(1), vx(2),vx(3),:))'
sum(val_orig)-1
val_updt = squeeze(tpm_updt(vx(1), vx(2),vx(3),:))'
sum(val_updt)-1



for ii=1:5
    tpm_ii = squeeze(tpm_orig(:,:,:,ii));
    tpm_ii(tpm_ii<thr) = thr;
    tpm_updt(:,:,:,ii) = tpm_ii;
end
% Adjust 'other' and fix it
tpm_ot = 1 - sum(tpm_updt(:,:,:,1:5),4);
% stpm_ot = smooth3(tpm_ot,'gaussian',1);
% stpm_ot = 10.^(smooth3(log10(tpm_ot),'gaussian',1));
% stpm_ot = 10.^(smooth3(log10(tpm_ot)));
stpm_ot = smooth3(tpm_ot);
% tpm_ot(tpm_ot<thr*15) = thr;
% tpm_ot(tpm_ot<thr) = thr;
% tpm_ot(tpm_ot<7.3e-5) = thr;
% tpm_ot(tpm_ot<thr) = thr;
stpm_ot(tpm_ot<thr) = thr;
% tpm_updt(:,:,:,6) = tpm_ot;
tpm_updt(:,:,:,6) = stpm_ot;

% Adjust GM/WM/CSF for 'other', that was fixed mainly in the brain volume
EStpm = sum(tpm_updt,4) - 1;
for ii=1:5
    tpm_updt(:,:,:,ii) = tpm_updt(:,:,:,ii) - EStpm/5;
end

% Save updated TPMs
Vtpm_u = spm_vol( fullfile(spm('dir'),'tpm','TPM.nii'));
fn_TPM_upd = fullfile(spm('dir'),'tpm','unwTPM_sl2.nii');
for ii=1:6
    Vtpm_u(ii).fname = fn_TPM_upd;
    Vtpm_u(ii) = spm_create_vol(Vtpm_u(ii));
    Vtpm_u(ii) = spm_write_vol(Vtpm_u(ii),tpm_updt(:,:,:,ii));
end

%% Visual check
spm_check_registration(...
    fullfile(spm('dir'),'tpm','TPM.nii,1'),...
    fullfile(spm('dir'),'tpm','TPM.nii,2'),...
    fullfile(spm('dir'),'tpm','TPM.nii,3'),...
    fullfile(spm('dir'),'tpm','TPM.nii,4'),...
    fullfile(spm('dir'),'tpm','TPM.nii,5'),...
    fullfile(spm('dir'),'tpm','TPM.nii,6'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,1'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,2'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,3'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,4'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,5'),...
    fullfile(spm('dir'),'tpm','unwTPM_sl2.nii,6'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,1'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,2'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,3'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,4'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,5'),...
    fullfile(spm('dir'),'tpm','nwTPM_sl2.nii,6')...
    );
