function OUTEEG = preproc_and_do_amica(INEEG,stage);
%OUTEEG = preproc_and_do_amica(INEEG,stage);
%
%Stage:
% 1 - interpolate to the specified 'montage', and
% 2 - restore the original recording reference, and
% 3 - re-reference based on EEG/REF channel types, and
% 4 - do highpass filter, and
% 5 - compute amica, and
% 6 - get dipole locations
% 


if stage>0, 
 %Interpolate channels that are bad (25% of time); reject epochs that are bad (25% of channels)
 myconfig.type      =    'channel';
 myconfig.montage   = INEEG.etc.montage_file_scalp_only; 
 myconfig.rejthresh =           25; %should decrease this == strict epoch keep, but what if no trials left??
 myconfig.rejreject =            0;
 myconfig.rejpctbad =           50; %should increase this == relaxed channel keep
 myconfig.interpmaj =            1;
 tmpEEG             = proc_interp(INEEG,myconfig);
end

if stage>1,
 %Restore the original reference site, re-reference to the average potential of all channels, remove baseline
 curref = tmpEEG.ref(1);
 ref    = readlocs(INEEG.etc.montage_file_all_sites);
 if  ~isempty(strmatch('Ref',curref)),                       %added 12-23-18, to accommodate CZ              %<-klooge, what about other refs (e.g., M1)??
   tmpEEG.data(end+1,:,:) = zeros(size(tmpEEG.data(1,:,:))); %insert dummy-reference
   ref    = ref(strmatch(curref,{ref.labels}));
   fields = fieldnames(tmpEEG.chanlocs); 
   fdiff1 = setdiff(fields, fieldnames(ref));
   fdiff2 = setdiff(fieldnames(ref), fields);
   if ~isempty(fdiff1) || ~isempty(fdiff2), 
     for ff = 1:length(fdiff1), if isfield(tmpEEG.chanlocs,fdiff1(ff)), tmpEEG.chanlocs = rmfield(tmpEEG.chanlocs,fdiff1(ff)); end; end
     for ff = 1:length(fdiff2), if isfield(            ref,fdiff2(ff)),             ref = rmfield(ref,            fdiff2(ff)); end; end
   end
   tmpEEG.chanlocs = [tmpEEG.chanlocs ref];
   tmpEEG.reject.rejmanualE = [tmpEEG.reject.rejmanualE; ones(1,size(tmpEEG.reject.rejmanualE,2))];
   tmpEEG.reject.rejglobalE = [tmpEEG.reject.rejglobalE; ones(1,size(tmpEEG.reject.rejglobalE,2))];
   tmpEEG.nbchan = length(tmpEEG.chanlocs);
 end                                                         %added 12-23-18, to accommodate CZ
end

if stage>2,
 tmpEEG.etc.proc_reference = tmpEEG.ref;

 %Regardless of the above reference site, do average reference here
 tmpEEG = eeg_checkset(tmpEEG);  tmpEEG.ref = 'common';      %
 newrefs= find(ismember({tmpEEG.chanlocs.type},{'EEG','REF'}));   %%<-why include "REF" here??
 tmpEEG = pop_reref( tmpEEG,newrefs,'keepref','on'); tmpEEG.ref = 'averef';      %re-reference to average, keeping old reference
 tmpEEG = pop_rmbase(tmpEEG,[],1:tmpEEG.pnts);               %remove baseline
end

%Still, only considering the original channels for ICA computation purposes
tmpEEG.etc.clean_channel_mask = [sum(tmpEEG.reject.rejmanualE,2)<tmpEEG.trials]';

if stage>3,
 %Do high-pass filter
 tmpEEG = pop_firws(tmpEEG, 'fcutoff', 1, 'ftype', 'highpass', 'wtype', 'kaiser', 'warg', 7.85726, 'forder',1286);
end

if stage>4,
 propvar_to_retain = .998; %retain this proportion of variance (eigenvectors) using PCA, eventually
                           %redo, and change this to 999? 
 %Perform epoch-wise rank() on data; #ICs should never exceed data rank!
 allranks = []; 
 for tt = 1:tmpEEG.trials, 
     datM     = squeeze(tmpEEG.data(tmpEEG.etc.clean_channel_mask,:,tt));
     [U,S,V]  = svd(datM,'econ'); 
     allranks = [allranks; length(find( cumsum(diag(S)./sum(diag(S))) < propvar_to_retain))];  
 end;
 if isfield(tmpEEG.etc, 'clean_channel_mask'), ncomps = min([min(allranks) sum(tmpEEG.etc.clean_channel_mask)-1]); else, ncomps = min(allranks); end; %me thinks the "-1" is not necessary

 %Estimate AMICA model based on the above rank estimation 
 pts2use = 1:tmpEEG.pnts;
 tmpdata = tmpEEG.data(tmpEEG.etc.clean_channel_mask, pts2use, 1:tmpEEG.trials);
 tmpdata = reshape(permute(double(tmpdata), [2 3 1]), [size(tmpdata,2)*size(tmpdata,3) size(tmpdata,1)])';
 if isa(tmpEEG.data,'single'), single(tmpdata); end

 %Downsample the points going into AMICA
 target_k = 40;
 subsample= round((target_k / (size(tmpdata,2)/(size(tmpdata,1)^2))) * length(pts2use));
 if subsample<length(pts2use), tmpdata  = tmpdata(:, randsample(pts2use,subsample)); end; %randsample(s, pts2use, subsample)); end

 rng shuffle; %this was necessary for multiple/parallel computing, because MATLAB rng resets each new instance of MATLAB, conflicting file names...
 disp(['   ' mfilename '; (' tmpEEG.subject ') number of independent components to compute with rank estimation: ' num2str(ncomps) ' from ' num2str(size(tmpdata,1)) ' channels.']);
 if exist([tmpEEG.etc.amicaroot, '/', tmpEEG.subject, tmpEEG.etc.amicasfx]), 
    eval(['!rm -rf ' tmpEEG.etc.amicaroot, '/', tmpEEG.subject, tmpEEG.etc.amicasfx]);
 end
 [tmpEEG.icaweights,tmpEEG.icasphere,tmpEEG.etc.icamods] = runamica15(tmpdata, ...
     'num_chans', size(tmpdata,1), ... 
     'pcakeep',   ncomps,          ...
     'outdir',    [tmpEEG.etc.amicaroot, '/', tmpEEG.subject, tmpEEG.etc.amicasfx], ...
     'do_reject', 1,               ...
     'numrej',   15,               ...
     'rejsig',    3,               ...
     'rejint',    3,               ...
     'num_models',1); 
 tmpEEG.icasphere   = tmpEEG.icasphere(1:tmpEEG.etc.icamods.num_pcs,:);
 tmpEEG.icachansind = find(tmpEEG.etc.clean_channel_mask);
 tmpEEG = eeg_checkset(tmpEEG, 'ica');
end
 
if stage>5,
 %Dipole fitting - adapted from: https://sccn.ucsd.edu/wiki/Makoto's_useful_EEGLAB_code
 dipfitroot = fileparts(which('eegplugin_dipfit'));
 dipfitlocs = [dipfitroot '/standard_BEM/elec/standard_1005.elc'];
 [chanlocs_out,coord_transform] = coregister(tmpEEG.chanlocs(tmpEEG.etc.clean_channel_mask & ismember({tmpEEG.chanlocs.type},'EEG')), ...
                                             dipfitlocs,'warp', 'auto','manual','off'); 
 tmpEEG = pop_dipfit_settings( tmpEEG, ...
       'hdmfile', [dipfitroot '/standard_BEM/standard_vol.mat'], ...
       'coordformat', 'MNI',...
       'mrifile', [dipfitroot '/standard_BEM/standard_mri.mat'],...
       'chanfile', dipfitlocs, ...
       'coord_transform', coord_transform,...
       'chansel', tmpEEG.icachansind);
 tmpEEG = pop_multifit(tmpEEG, 1:length(tmpEEG.dipfit.chansel),'threshold', 100, 'dipplot','off','plotopt',{'normlen' 'on'});
 %tmpEEG = fitTwoDipoles(tmpEEG, 'LRR', 35); %in case it makes sense to also fit dual bilateral dipoles for comps

 dipfit = tmpEEG.dipfit;
 save([tmpEEG.etc.amicaroot, '/', tmpEEG.subject, tmpEEG.etc.amicasfx, '/dipfit'],'dipfit');

end
 
%Replace original EEG variable with tmpEEG.
OUTEEG = tmpEEG;

