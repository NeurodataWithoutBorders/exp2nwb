Janelia Research Campus
-----------------------

![jrc_logo_180x40](https://cloud.githubusercontent.com/assets/1093770/24422906/2ba9caae-13c9-11e7-9177-c54f5c2c1f62.png)

Copyright (c) 2017, HHMI-Janelia Research Campus
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
      
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
      
    * Neither the name HHMI-Janelia Research Campus nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL HHMI-JANELIA RESEARCH CAMPUS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Contributor: Gennady Denisov

Acknowledgements: Karel Svoboda, Diego Gutnisky and Jeff Teeters 


Introduction
------------
According to the NWB API specification,

https://github.com/NeurodataWithoutBorders/specification/blob/master/version_1.0.5_beta/NWB_file_format_specification_1.0.5g_beta.pdf,

NWB data format is a HDF5 format with neurophysiology data-specific structure of
groups and datasets. In particular, at the top level of a NWB file, six groups 
with required/fixed names must be present ("acquisition", "analysis", "epochs", 
"general", "processing" and "stimulus"), together with several datasets, all of 
them also with required names. The top groups contain other (nested) groups and 
datasets, some of them also possessing required names, and others possessing 
user-specified names, depending on particular data.

The Python framework stored in this folder implements next generation 
approach to practical conversion of experimental neurophysiology data from their 
original formats to the NWB format: EXP => NWB. The approach is expected to 
facilitate the conversion procedure from any original format, although will 
require partial re-implementation for new input data.

![exp2nwb_flowchart_new](https://cloud.githubusercontent.com/assets/1093770/24525299/c3dc6738-1567-11e7-82c9-7bd2f4f761b1.png)

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
with domain scientists, depending on their data. 

As a sample, we provide implementation of exp2dict.py for two datasets, both generated at Svoboda lab
in Janelia and uploaded to the crcns.org website:
- ephys dataset (Li et al, Nature 2015, doi: 10.1038/nature14178), 

      https://portal.nersc.gov/project/crcns/download/alm-1     

      https://portal.nersc.gov/project/crcns/download/alm-2

 and

- ophys datsets (Peron et al, Neuron 2015, doi: 10.1016/j.neuron.2015.03.027) 

      https://portal.nersc.gov/project/crcns/download/ssc-1
      
A list of dictionary keys that are exoected to be defined when re-implementing the functions in library exp2dict is provided in doc/README_exp2dict_keys.

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

For the two datasets mentioned above, the full sets of supported p-strings is listed in doc/README_p-strings.

The current implementation of both exp2nwb.py and libexp2dict.py assumes the 
imput data files are in the HDF5 format. The experimental data at are actually
stored in MAT format, can be converted to HDF5 format using our script
    https://github.com/NeurodataWithoutBorders/mat2nwb/blob/master/mat2h5.py

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

