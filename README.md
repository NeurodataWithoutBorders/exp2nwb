Introduction
------------
The Python framework stored in this folder implements next generation 
approach to practical conversion of experimental neurophysiology data from their 
original formats to the NWB format: EXP => NWB. The approach is expected to 
facilitate the conversion procedure from any original format, although will 
require partial re-implementation for new input data.

(Flow chart)

The distinctive features of this approach are:
- decomposition of the conversion project NWB file into a set of "elementary" 
  tasks, each of which can be performed independently. Each the task will produce 
  a partial NWB file, where only a certain group or subset of nested groups of 
  the target NWB file will be created and populated with datasets, while the rest 
  of the groups and datasets are either missing or empty;
- further decomposition of each the "elementary" task in two steps:
  1) (exp2dict): extract a relavant portion of experimental data from input file(s)
                 and store it in a Python dictionary. The keys in this dictionary 
                 will be the full paths to datasets or attributes in the target 
                 NWB file, and the values will be the contents of the datasets or 
                 attributes; and
  2) (dict2nwb): given a dictionary generated for a given task, produce the 
                 corresponding partial NWB file.
- each the partial NWB file thus produced is named accordingly with the 
  portion of data it contains, is stored in a project directory and  can 
  be subsequently reviewed using an HDF viewer. For example, the partial 
  NWB file where only data in the group "/general/subject" are populated 
  will be named "general.subject.h5". 
- after all the required partial NWB files have been produced, the final,
  full-size NWB file will be generated in one step via assembly of the partial
  files.

The exact set of the partial NWB files to be produced will depend on the type of 
the data being converted. Depending on the data, some of the partial files 
may or may not be produced. We discriminate between two major types of 
neurophysiology data:
1) electrophysiology ("ephys") data, where the brain activity is recorded using 
   either extracellular or intracellular electrode probes (or both); and
2) optophysiology ("ophys") data, where the brain activity is recorded using
   calcium fluirescence imaging.

The functions in the library dict2nwb are expected to be "standard",
i.e. not to change from one dataset to another, whreas the functions
in the library exp2dict should be re-implemented in conjunction
with domain scientists, depending on their data. As a sample, we provide 
implementation of exp2dict for two datasets, both generated at Svoboda lab
in Janelia and uploaded to the crcns.org website:
- ephys dataset (Li et al, Nature 2015), 
      https://portal.nersc.gov/project/crcns/download/alm-1     
      https://portal.nersc.gov/project/crcns/download/alm-2
  and
- ophys datsets (Peron et al, Neuron 2015) 
      https://portal.nersc.gov/project/crcns/download/ssc-1
According to the NWB API specification,

https://github.com/NeurodataWithoutBorders/specification/blob/master/version_1.0.5_beta/NWB_file_format_specification_1.0.5g_beta.pdf,

NWB data format is a HDF5 format with neurophysiology data-specific structure of
groups and datasets. In particular, at the top level of a NWB file, six groups 
with required/fixed names must be present ("acquisition", "analysis", "epochs", 
"general", "processing" and "stimulus"), together with several datasets, all of 
them also with required names. The top groups contain other (nested) groups and 
datasets, some of them also possessing required names, and others possessing 
user-specified names, depending on particular data.

Processing strings (p-strings)
------------------------------
converting the data of each of these types may be characterized by a set of 
dot-separated processing strings (p-strings), so that 
- each the p-string matches the name of the partial NWB file to be produced, and
- the dot-separated substrings in the p-string match the substrings in the
  slash-separated path to the data that will be stored in this file. 
The p-strings are actually used as inputs by our conversion approach, as will be 
described below, in section "CODE USAGE".

The p-strings that are expected to be common to both the ephys and ophys data are:
common_strings = ["top_datasets", 
                  "acquisition.timeseries.<behavioral_data_name>",
                  "analysis", 
                  "epochs", 
                  "general.devices", "general.subject", "general.top_datasets", 
                  "processing.<Module_name>.BehavioralEvents.<behavioral_data_name>",
                  "processing.<Module_name>.BehavioralTimeseries.<behavioral_data_name>",
                  "stimulus.presentation.<stimulus_name>"],
where 
- "top_datasets" represents a set of required datasets to be stored at 
  a given level in the NWB file, while all other strings correspond to groups.
- the user-specified <Module name> may be "Auditory", "Licks", "Pole", etc.,
  depending on the data
- the user-specified <behavioral_data_name> may be "reward_cue", "lick_left",
  "lick_right", "pole_accessible", etc., depending on the data
- the user-specified <stimulus_name> may be "photostimulus", "pole_in", "pole_out",
  etc., depending on the data

The expected ephys data-specific p-strings are:
ephys_strings = ["acquisition.timeseries.<ephys_type>cellular_traces",
                 "acquisition.timeseries.<behavioral_data_name>_trace",
                 "general.<ephys_type>cellular_ephys", 
                 "general.optogenetics",
                 "processing.<ephys_type>cellular_units.top_datasets",
                 "processing.<ephys_type>cellular_units.<ephys_data_name>"]
where
- the user-specified <ephys_type> may be "extra", "intra", etc., and
- the user-specified <ephys_data_name> may be one of ["EventWaveform", "UnitTimes"].

The expected ophys data-specific strings are:

ophys_strings = ["acquisition.images",
                 "acquisition.timeseries",
                 "general.optophysiology",
                 "processing.ROIs"],
where
- the partial NWB files acquisition.images.h5 and acquisition.timeseries.h5
  are supposed to store only the image-related, and not the behavioral data;
  the latter are supposed to be handled by other p-strings, as described above.

The code comprises a Python script exp2nwb.py, together with utility library 
files:
- libexp2dict.py defines the functions that create a dictionary for each 
  "elementary" task. Each the function is named according to the "elementary" 
  task it will be used by, or to the partial file that will eventually be 
  produced. For example, the function to create a dictionary that will be used 
  to produce the partial file .general.subject.h5. is named .general_subject. (i.e. with .dot. replaced by .underscore.);
- libdict2nwb.py, likewise, defines the functions that will produce partial 
  files, given their respective dictionaries. These functions have the same 
  names as the corresponding functions in the library file libexp2dict.py; and
While the definitions of the functions in the library libexp2dict.py are highly 
data-dependent and will have to be re-implemented for each particular 
experimental data source, the functions in library libdict2nwb.py are supposed 
to be .universal., i.e. applicable to conversion of potentially any data. 
The library
- libh5.py is only used for the purpose of illustration. It defines the utility 
functions that can only be used with the data produced by Svoboda lab at Janelia 
Research Campus. These functions have been used by our current implementation 
of the library libexp2dict.py. 

Code usage
----------
The code comprises a script exp2nwb.py together with three library files: 
libexp2dict.py, libdict2nwb.py and util.py. Given a processing string <pstring>, 
the script employs a function from the library libexp2dict.py to create 
a dictionary for each particular task and then, given the dictionary, a 
function from the library libdict2nwb.py to create a particual NWB file <pstring>.h5. 
The names of these functions are typically the same for the two libraries and match 
the processing string, but with "dots" replaced by "underscore" symbols. The 
library util.py have been employed by functions in libexp2dict.py in order to 
extract data from input file(s).

For the two datasets mentioned above, the full sets of supported p-strings are:
- for ephys datsets,
    ephis_pstrings = ["top_datasets", \
                      "analysis",\
                      "epochs", \
                      "general.devices", \
                      "general.subject", \
                      "general.top_datasets"]
                      "acquisition.timeseries.extracellular_traces",\
                      "acquisition.timeseries.lick_trace",\
                      "general.extracellular_ephys", \
                      "general.optogenetics",\
                      "processing.extracellular_units.top_datasets",\
                      "processing.extracellular_units.EventWaveform",\
                      "processing.extracellular_units.UnitTimes",\
                      "stimulus.presentation.auditory_cue",\
                      "stimulus.presentation.photostimulus",
                      "stimulus.presentation.pole_in",\
                      "stimulus.presentation.pole_out"]

- for ophys dataset,
    ophys_pstrings = ["top_datasets", \
                      "analysis",\
                      "epochs", \
                      "general.devices", \
                      "general.subject", \ 
                      "general.top_datasets", \
                      "acquisition.images",\
                      "acquisition.timeseries",\
                      "general.optophysiology",\
                      "processing.Auditory.BehavioralEvents.reward_cue",\
                      "processing.Licks.BehavioralEvents.lick_left",\
                      "processing.Licks.BehavioralEvents.lick_right",\
                      "processing.Pole.BehavioralEvents.pole_accessible",\
                      "processing.Reward.BehavioralEvents.water_left_reward", \
                      "processing.Reward.BehavioralEvents.water_right_reward",\
                      "processing.ROIs",\
                      "processing.Whisker.BehavioralEvents.pole_touch_protract", \
                      "processing.Whisker.BehavioralEvents.pole_touch_retract", \
                      "processing.Whisker.BehavioralTimeSeries.whisker_angle", \
                      "processing.Whisker.BehavioralTimeSeries.whisker_curvature", \
                      "stimulus.presentation.auditory_cue",\
                      "stimulus.presentation.pole_accessible", \
                      "stimulus.presentation.water_left", \
                      "stimulus.presentation.water_right",\
                      "stimulus.presentation.zaber_motor_pos"]

The current implementation of both exp2nwb.py and libexp2dict.py assumes the 
imput data files are in the HDF5 format. Those users who's data are actually
stored in MAT format, can use our sctipt 
    https://github.com/NeurodataWithoutBorders/mat2nwb/blob/master/mat2h5.py
in order to convet the MAT files into HDF5 files.

To see the usage of the script exp2nwb.py, type the name of the script:
    $ exp2nwb.py
Usage: 
    exp2nwb.py <data>.h5 [ <meta_data>.h5] [options (-h to list)]

All the available command line options can be viewed by typing
    $ exp2nwb.py -h

In particular, to produce a single partial NWB file <pstring>.h5 from the input 
data, type
    $ exp2nwb.py <data>.h5 [ <meta_data>.h5]  -s <pstring> -d <project_name>
The partial NWB file will be stored in the project folder <project_name>, 
which by default is named <data>.

To assemble all the partial NWB files in folder <project_name>, type
    $ exp2nwb.py <project_folder>.
This command will produce a full NWB file named <project_name>.nwb.
After the assembly, the project folder will be deleted by default, unless
the debugging mode has been used, as specified by the option "-D".


Requirements for the re-implementation of the library exp2dict.py
-----------------------------------------------------------------
While the library libdict2nwb.py may in some cases employ generic
functions to produce partial NWB files storing time series data
with a given name, the names of the functions in the library exp2dict.py 
always exactly match the names of the processing strings, but with the "dots"
replaced by the "underscore" sysmbols. Listed below are the dictionary keys
that are required to be created by each the function in the library exp2dict.py.

top_datasets keys:
    file_name
    start_time
    identifier
    description
    mode 

acquisition_images keys:
    acquisition.images.description
    acquisition.images.fov_###_green
    acquisition.images.fov_###_green.attrs
    acquisition.images.fov_###_red
    acquisition.images.fov_###_red.attrs
    acquisition.images.<dataset>                  
    acquisition.images.<dataset>.attrs           

acquisition_timeseries keys:
    acquisition.timeseries.description
    acquisition.timeseries.fov_###.dimension
    acquisition.timeseries.fov_###.external_file
    acquisition.timeseries.fov_###.external_file.attrs # data attributes
    acquisition.timeseries.fov_###.field_of_view
    acquisition.timeseries.fov_###.format
    acquisition.timeseries.fov_###.imaging_plane
    acquisition.timeseries.fov_###.scan_line_rate
    acquisition.timeseries.fov_###.timestamps

acquisition_timeseries_extracellular_traces keys:
    acquisition.timeseries.extracellular_traces

acquisition_timeseries_<series_name> keys:
    acquisition.timeseries.<series_name>.attrs         # group attributes
    acquisition.timeseries.<series_name>.data.attrs    # data attributes
    acquisition.timeseries.<series_name>.num_samples
    acquisition.timeseries.<series_name>.timestamps

analysis keys:
    analysis.description
    analysis.good_trials
    analysis.trial_start_times
    analysis.trial_type_mat
    analysis.trial_type_string

epochs keys:
    epochs.trial_###.description
    epochs.trial_###.start_time
    epochs.trial_###.stop_time
    epochs.trial_###.tags
    epochs.trial_###.units_present

general_divices keys:
    general.divices.ephys_acquisition
    general.divices.photostim_source

general_extracellular keys:
    general.extracellular.ADunit
    general.extracellular.data_types
    general.extracellular.ground_coordinates
    general.extracellular.impedance
    general.extracellular.penetration_num
    general.extracellular.recording_marker
    general.extracellular.recording_type

general_subject keys:
    general.subject.description
    general.subject.genotype
    general.subject.subject_id

general_optogenetics keys:
    general.optogenetics.site_1.device
    general.optogenetics.site_1.excitation_lambda
    general.optogenetics.site_1.location
    general.optogenetics.site_1.stimulation_method

general_optophysiology keys:
    general.optophysiology.description
    general.optophysiology.fov_###.channel_green.description
    general.optophysiology.fov_###.channel_green.emission_lambda
    general.optophysiology.fov_###.channel_red.description
    general.optophysiology.fov_###.channel_red.emission_lambda
    general.optophysiology.fov_###.description
    general.optophysiology.fov_###.device
    general.optophysiology.fov_###.excitation_lambda
    general.optophysiology.fov_###.imaging_rate
    general.optophysiology.fov_###.indicator
    general.optophysiology.fov_###.location
    general.optophysiology.fov_###.manifold
    general.optophysiology.fov_###.manifold.attrs
    general.optophysiology.fov_###.reference_frame

processing_extracellular_units_top_datasets keys:
    processing.extracellular_units.description
    processing.extracellular_units.identification_method
    processing.extracellular_units.spike_sorting
 
processing_extracellular_units_EventWaveform keys:
    processing.extracellular_units.EventWaveform.attrs                          # group attrs
    processing.extracellular_units.EventWaveform.unit_##.data
    processing.extracellular_units.EventWaveform.unit_##.data.attrs
    processing.extracellular_units.EventWaveform.unit_##.electrode_idx
    processing.extracellular_units.EventWaveform.unit_##.num_samples
    processing.extracellular_units.EventWaveform.unit_##.sample_length
    processing.extracellular_units.EventWaveform.unit_##.source
    processing.extracellular_units.EventWaveform.unit_##.timestamps

processing_extracellular_units_UnitTimes keys:
    processing.extracellular_units.UnitTimes.attrs                              # group attrs
    processing.extracellular_units.UnitTimes.cell_types
    processing.extracellular_units.UnitTimes.unit_##.source
    processing.extracellular_units.UnitTimes.unit_##.times
    processing.extracellular_units.UnitTimes.unit_##.trial_ids
    processing.extracellular_units.UnitTimes.unit_##.unit_description

processing_<Module_name>_BehavioralEvents_<behavioral_data_name> keys:
    processing.<Module_name>.BehavioralEvents.<behavioral_data_name>.attrs      # group attrs
    processing.<Module_name>.BehavioralEvents.<behavioral_data_name>.data
    processing.<Module_name>.BehavioralEvents.<behavioral_data_name>.data.attrs
    processing.<Module_name>.BehavioralEvents.<behavioral_data_name>.num_samples
    processing.<Module_name>.BehavioralEvents.<behavioral_data_name>.timestamps

processing_<Module_name>_BehavioralTimeSeries_<behavioral_data_name> keys:
    processing.<Module_name>.BehavioralTimeseries.<behavioral_data_name>.attrs  # group attrs
    processing.<Module_name>.BehavioralTimeseries.<behavioral_data_name>.data
    processing.<Module_name>.BehavioralTimeseries.<behavioral_data_name>.data.attrs
    processing.<Module_name>.BehavioralTimeseries.<behavioral_data_name>.num_samples
    processing.<Module_name>.BehavioralTimeseries.<behavioral_data_name>.timestamps

processing_ROIs keys:
    processing.ROIs.DfOverF.attrs                                               # group attrs
    processing.ROIs.DfOverF.fov_###.data
    processing.ROIs.DfOverF.fov_###.data.attrs                                  # data attrs
    processing.ROIs.DfOverF.fov_###.roi_names
    processing.ROIs.DfOverF.fov_###.timestamps
    processing.ROIs.DfOverF.fov_###.trial_ids
    processing.ROIs.ImageSegmentation.attrs                                     # group attrs
    processing.ROIs.ImageSegmentation.fov_###.<roi_id_##>.master1_shape
    processing.ROIs.ImageSegmentation.fov_###.<roi_id_##>.pixmap
    processing.ROIs.ImageSegmentation.fov_###.<roi_id_##>.weight
    processing.ROIs.ImageSegmentation.fov_###.description
    processing.ROIs.ImageSegmentation.fov_###.imaging_plane_name
    processing.ROIs.ImageSegmentation.fov_###.ref_image_green
    processing.ROIs.ImageSegmentation.fov_###.ref_image_red
    processing.ROIs.ImageSegmentation.fov_###.roi_ids
    processing.ROIs.description

stimulus_presentation_<stimulus_data_name> keys:
    stimulus.presentation.<stimulus_data_name>.attrs                            # group attrs
    stimulus.presentation.<stimulus_data_name>.data
    stimulus.presentation.<stimulus_data_name>.data.attrs                       # data attrs
    stimulus.presentation.<stimulus_data_name>.num_samples
    stimulus.presentation.<stimulus_data_name>.timestamps

