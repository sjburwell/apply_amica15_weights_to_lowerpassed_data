# apply_amica15_weights_to_lowerpassed_data
The below code will take an EEG structure containing continuous data that is filtered at some low frequency cutoff (e.g., 0.1 Hz), calculate AMICA on a higher-passed version of the data (e.g., 1Hz), and incorporate that (better) ICA decomposition into the lower-passed data. 

# Usage:
First, you need to have both EEGLAB, [eeg_commander](https://github.com/sjburwell/eeg_commander), and AMICA (vers. 15) in your path. Also, you should have the dipfit toolbox installed to EEGLAB -- this can be done using the EEGLAB Plug-in GUI.
```matlab
addpath /path/to/eeglab         %recommended to use version 14.1.1b or above
eeglab;                         %keeping EEGLAB GUI open is optional

addpath /path/to/eeg_commander  %obtained in bash by "git clone <repo>"
eeg_commander_startup;          %add necessary eeg_commander paths

addpath /path/to/amica15        %add necessary path to AMICA (version 15)
```

Next, once a continuous EEG structured variable is loaded into MATLAB memory, you need to provide the structure with pointers to some montage files which can be read using readlocs() (cf. [Importing Channel Locations](https://sccn.ucsd.edu/wiki/A03:_Importing_Channel_Locations)). One of these files should point to the most complete montage of channel locations that could exist for a given file, especially that which might contain recording reference sites. Another one contains only those locations thought to reflect "scalp channels" (e.g., Cz, Pz, etc.). 
```matlab
EEG.etc.montage_file_all_sites  = '/path/to/eeg_commander/misc/locs/montage10-10_sphrad1_n67.ced';  %contains all scalp chans, refs, EOGs, etc.
EEG.etc.montage_file_scalp_only = '/path/to/eeg_commander/misc/locs/montage10-10_sphrad1_n61.ced'; %contains only scalp channels
overwriteamica =           1 ;  %compute AMICA regardless if it's been generated before
cntdir         = EEG.filepath;  %directory containing the continuous EEG file (e.g., '.../subID12345'), will be the place where amica directory is output.
```

Finally, you can run the script run the script *run_AMICA_high2lowerHz* which will make a couple of calls to the function *preproc_and_do_amica()*. In that function, it will interpolate the data to a specified montage, restore the data to its original recording reference and then re-reference based on EEG/REF channel types, do a high-pass filter, compute AMICA and get corresponding dipoles. The final output will be a low-filter-cutoff (e.g., 0.1 Hz) EEG structure containing AMICA decomposition from a higher-filter-cutoff (e.g., 1.0 Hz) dataset. I.e.,  
```matlab
EEG % EEG data is 0.1 Hz highpass filtered, does or does not contain potentially lousy AMICA decomposition
run_AMICA_high2lowerHz;
EEG % EEG data is 0.1 Hz highpass filtered, contains probably *better* AMICA decomposition performed on 1 Hz highpass filtered data.
```

