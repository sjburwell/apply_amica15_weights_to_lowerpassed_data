
% Interpolate to full montage, recover reference channel, average re-reference, 1.0 Hz highpass, run AMICA, get dipoles
EEG.etc.amicaroot = cntdir; 
EEG.etc.amicasfx  = '_amicacnt';
if ~exist(     [EEG.etc.amicaroot, '/', EEG.subject, EEG.etc.amicasfx])  || ... 
   isempty(dir([EEG.etc.amicaroot, '/', EEG.subject, EEG.etc.amicasfx])) || ...
   overwriteamica==1,
   EEGica = preproc_and_do_amica(EEG,6);
end

% Insert pre-existing AMICA/DIPFIT information (regardless if overwriteamica==1 or overwriteamica==0)
EEG = preproc_and_do_amica(EEG,3);
icamod = loadmodout15([EEG.etc.amicaroot, '/', EEG.subject, EEG.etc.amicasfx]);
EEG.icaweights  = icamod.W; 
EEG.icasphere   = icamod.S(1:icamod.num_pcs,:);
EEG.icachansind = find(EEG.etc.clean_channel_mask);
EEG = eeg_do_icarms(EEG); %updated for 20200215
EEG = eeg_checkset(EEG,'ica');
if exist([EEG.etc.amicaroot, '/', EEG.subject, EEG.etc.amicasfx, '/dipfit.mat']),
   dipfit = load([EEG.etc.amicaroot, '/', EEG.subject, EEG.etc.amicasfx, '/dipfit.mat']);
   EEG.dipfit = dipfit.dipfit;
end
EEG = proc_select(EEG,'channel',{EEG.chanlocs(find(EEG.etc.clean_channel_mask)).labels});



