#!/usr/local/python-2.7.11/bin/python
#
# Copyright (c) 2017, HHMI-Janelia Research Campus All rights reserved. 
#
# exp2nwb.py 
#
# Conversion of experimental neurophysiology data from their
# original formats to the NWB format
#
import sys
import os

import nwb
from nwb import nwb_file
from nwb import nwb_utils
import h5py                    # used to extract data from the input h5 file(s)
import datetime
import unicodedata
import getpass
import numpy as np
from sets import Set
import re 
import optparse
import shutil

import util  
import libexp2dict
import libdict2nwb

# ------------------------------------------------------------------------------

def make_nwb_command_line_parser(parser):
    parser.add_option("-d", "--project_dir", dest="project_dir", help="project directory name", default="")
    parser.add_option("-D", "--debug",   action="store_true",  dest="debug", help="output debugging info; don't delete project dir", default=False)
    parser.add_option("-e", "--no_error_handling", action="store_false", dest="handle_errors", help="handle_errors", default=True)
    parser.add_option("-r", "--replace", action="store_true", dest="replace", help="if output file exists, replace it", default=False)
    parser.add_option("-s", "--pathstr", dest="path_str", help="dot-separated path string",metavar="path_str",default="")
    parser.add_option("-t", "--type",    dest="data_type", help="specify explicitly: 'ephys' or 'ophys'; otherwise will be determined automatically", default="")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="increase the verbosity level", default=False)

    return parser

# ------------------------------------------------------------------------------

def set_data_type(data_h5, meta_h5, options):
    if len(options.data_type) == 0:                 
        # Determine the data type automatically
        keys = [str(i) for i in util.get_key_list(data_h5["trialPropertiesHash"])]
        keys2 = util.get_key_list(data_h5["timeSeriesArrayHash"])
        if not "LickTime"   in util.get_key_list(data_h5["trialPropertiesHash"]) and \
           not "EphusVars"  in util.get_key_list(data_h5["timeSeriesArrayHash"]):
            options.data_type = "ophys"
        elif "EphusVars" in keys2:
            options.data_type = "ephys"
    if options.data_type not in ["ephys", "ophys"]:
        sys.exit("Could not determine the data type")

    if options.verbose:
        print("\nData type: " + options.data_type + "\n")
    return options

# ------------------------------------------------------------------------------

def produce_partial_nwb(data_path, metadata_path, options):
    if data_path.split(".")[-1] == "h5":
        data = h5py.File(data_path, "r")
        if len(metadata_path) > 0 and metadata_path.split(".")[-1] == "h5":
            metadata = h5py.File(metadata_path, "r")
        elif len(metadata_path) == 0:
            metadata = data
        else:
            sys.exit("\nInput metadata from files of this type is not supported")
    else:
        sys.exit("\nInput data from files of this type is not supported")
    options = set_data_type(data, metadata, options)

    # Extract (meta)data
    dict = libexp2dict.make_dict(data, metadata, options)

    # Populate (meta)data
    if len(dict.keys()) == 0:
        sys.exit("Computing the dictionary is not currently implemented")

    libdict2nwb.make_partial_nwb(options.path_str, dict, options.project_dir)

# ------------------------------------------------------------------------------

def extract_attrs(attr_pointer):
    attrs = {}
    for k in attr_pointer.keys():
        if not k in ["neurodata_type", "help", "interval", "ancestry"]:
            attrs[k] = attr_pointer[k]
    return attrs

# ------------------------------------------------------------------------------

#
# Recursively parse the group tree  in the H5 object
# and copy its contents to the NWB object
#
def update_nwb_object(last_item, curr_nwb_group, curr_h5_path, h5_object, nwb_object):  
    update_object = 1
    if hasattr(h5_object[curr_h5_path], '__dict__'):
        # Current path in h5_data points to a group    
        for k in h5_object[curr_h5_path].keys():
            curr_h5_path1 = curr_h5_path + "/" + str(k)
            my_attrs = extract_attrs(h5_object[curr_h5_path1].attrs)
#           print "curr_h5_path=", curr_h5_path, " k=", k, " my_attrs=", my_attrs
            num_keys = 0
            try:
                num_keys = len(h5_object[curr_h5_path1].keys())
            except:
                num_keys = 0
            if num_keys == 0:
                # Item is a dataset
#               print "    Item is a dataset"
                if k in ["nwb_core.py"]:
                    continue
                try:
                    if len(my_attrs.keys()) > 0:
                        curr_nwb_group.set_dataset(str(k), np.array(h5_object[curr_h5_path1]), attrs=my_attrs, abort=False)
                    else:
                        curr_nwb_group.set_dataset(str(k), np.array(h5_object[curr_h5_path1]), abort=False)
                except:
                    try:
                        curr_nwb_group.set_custom_dataset(str(k), np.array(h5_object[curr_h5_path1]))
                    except:
                        pass
            else:
#               print "    Item is a group"
                # Next-level item is a group
#               assert(len(h5_object[curr_h5_path1].keys()) > 0)
                try:
                    ancestry = "<" + h5_object[curr_h5_path1].attrs["ancestry"][-1] + ">"
                except:
                    ancestry = ""
#               print "ancestry=", ancestry, " k=", k, " curr_h5_path1=", curr_h5_path1, " attrs=", my_attrs
                if len(ancestry) > 0:
                   try:
                       curr_nwb_group1 = nwb_object(ancestry, str(k), \
                                             path=curr_h5_path, attrs=my_attrs, abort=False)
                   except:
                       try:
                           curr_nwb_group1 = curr_nwb_group.make_group(ancestry, str(k), \
                                                 attrs=my_attrs, abort=False)
                       except:
                           try:
                               curr_nwb_group1 = curr_nwb_group.make_group(ancestry, str(k), abort=False)
                           except:
                               try:
                                   curr_nwb_group1 = nwb_object.create_group(curr_h5_path)
#                                  h5_object.copy(curr_h5_path1, nwb_object.file_pointer[curr_h5_path])
                               except:
                                   # Group already exists
                                   try:
                                       curr_nwb_group1 = nwb_object.file_pointer[curr_h5_path1]
                                   except:
                                       curr_nwb_group1 = nwb_object.file_pointer[curr_h5_path]
                               try:
                                   h5_object.copy(curr_h5_path1, nwb_object.file_pointer[curr_h5_path])
                               except:
                                   pass
                else:  
                   try: 
                       curr_nwb_group1 = curr_nwb_group.make_group(str(k), \
                                         attrs=extract_attrs(h5_object[curr_h5_path].attrs), abort=False)
                   except:
                       try:
                           curr_nwb_group1 = nwb_object.create_group(curr_h5_path1)
                           h5_object.copy(curr_h5_path1, curr_nwb_group1)
                       except:
                           curr_nwb_group1 = nwb_object.file_pointer[curr_h5_path1]
                nwb_object = update_nwb_object(last_item, curr_nwb_group1, curr_h5_path1, h5_object, nwb_object)
    else:
        # Current path in h5_data points to a dataset
        path_items = curr_h5_path.split("/")
        upper_level_path = "/".join(path_items[0:-1])
        nwb_object.file_pointer[upper_level_path].create_dataset(k, data = np.array(h5_object[curr_h5_path], \
            attrs=extract_attrs(h5_object[curr_h5_path].attrs)))
#       print " dataset=", curr_h5_path, " attrs=", h5_object[curr_h5_path].attrs.keys()
    return nwb_object

# ------------------------------------------------------------------------------

def assemble_nwb_from_partial_files(project_dir, options):
    # Initialize NWB object
    vargs = {}
    try:
        first_infut_file = os.path.join(project_dir, "top_datasets.h5")
        top_h5           = h5py.File(first_infut_file, "r")
        vargs["file_name"]  = project_dir + ".nwb"
        vargs["identifier"] = nwb_utils.create_identifier(project_dir)
        vargs["mode"]       = "w"
        vargs["start_time"] = np.array(top_h5["session_start_time"]).tolist()[0]
        vargs["description"]= np.array(top_h5["session_description"]).tolist()[0]
        if options.verbose:
            print "\ndescription=", vargs["description"]
    except:
        sys.exit("Could not read info from file " + first_infut_file)
    top_h5.close()

    if os.path.exists(vargs["file_name"]):
        os.remove(vargs["file_name"])
    nwb_object = nwb_file.open(**vargs)

    for file in os.listdir(project_dir):
        if file == "top_datasets.h5":
            continue

        
        path_str = file[0:-3]
        if options.verbose:
            print "\nInput file=", file, " path_str=", path_str, "\n"

        file_path = os.path.join(project_dir, file)
        data_h5 = h5py.File(file_path, "r")

        # Create a pointer to the current location in H5 file
        groups = path_str.split(".")
        last_group = groups[-1]
        if last_group == "top_datasets":
            groups.pop(-1)
        h5_pointer = data_h5["/".join(groups)]

        # Sinchronize paths between the input H5 and output NWB files
        curr_h5_path = "/"      
        for i in range(len(groups)):
            if options.verbose:
                print "... Creating group ", groups[i], " curr_h5_path=", curr_h5_path

            if curr_h5_path in ["/stimulus/presentation", "/stimulus/templates", "/acquisition/timeseries"]: 
                # Extract the actual series_type from the partial .h5 file and use that in the nwb_object.make_group
                group_path = curr_h5_path + "/" + groups[i]
                ancestry = "<" + data_h5[group_path].attrs["ancestry"][-1] + ">"
                curr_nwb_group = nwb_object.make_group(ancestry, groups[i], path = curr_h5_path, abort=False)
            elif curr_h5_path == "/processing":
                curr_nwb_group = nwb_object.make_group("<Module>", groups[i], abort=False)     
            elif re.search("/processing/", curr_h5_path):
                if groups[i] in ["BehavioralEvents", "BehavioralTimeSeries",\
                                 "DfOverF", "ImageSegmentation", "EventWaveform", "UnitTimes"]:
                    curr_nwb_group = curr_nwb_group.make_group(groups[i], abort=False)
                else:
                    group_path = curr_h5_path + "/" + groups[i]
                    ancestry = "<" + data_h5[group_path].attrs["ancestry"][-1] + ">"
                    curr_nwb_group = curr_nwb_group.make_group(ancestry, groups[i]) 
            else:
                curr_nwb_group = nwb_object.make_group(groups[i], path = curr_h5_path, abort=False)
            if options.verbose:
                print "curr_nwb_group.name=", curr_nwb_group.name
            if i == 0:
                curr_h5_path += groups[i]
            else:
                curr_h5_path += "/" +  groups[i]

            if options.verbose:
                print "curr_h5_path=", curr_h5_path
           
        # Recoursively parse the groups tree in the H5 file 
        # and copy its contents to the NWB file
        nwb_object = update_nwb_object(last_group, curr_nwb_group, curr_h5_path, data_h5, nwb_object)            

    nwb_object.close()
    
# ------------------------------------------------------------------------------

def create_and_assemble_all_partial_files(data_path, metadata_path, options):

    data_h5 = h5py.File(data_path, "r")
    if len(metadata_path) > 0:
        meta_h5 = h5py.File(metadata_path, "r")
    else:
        meta_h5 = data_h5
    options = set_data_type(data_h5, meta_h5, options)

    common_pstrings = ["top_datasets", \
                       "analysis",\
                       "epochs", \
                       "general.devices", \
                       "general.subject", \
                       "general.top_datasets"]
    ephys_pstrings =  ["acquisition.timeseries.extracellular_traces",\
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
    ophys_pstrings =  ["acquisition.images",\
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
    pstrings = common_pstrings
    if options.data_type == "ephys":
        pstrings = pstrings + ephys_pstrings
    else:
        pstrings = pstrings + ophys_pstrings

    # Create partial NWB files
    pstrings = sorted(pstrings)
    for s in range(len(pstrings)):
        command = "/home/denisovg/work/Karel_Svoboda/exp2nwb.py " + data_path + " " + metadata_path + " " 
        command += " -s " + pstrings[s] + " -d " + options.project_dir
        if options.verbose:
            print "Running command: " + command
        os.system(command)
                    
    # Assemble all
    if options.debug:
        sys.stderr.write("\nStarting assembly...\n")
    os.system("python exp2nwb.py " + options.project_dir)
    if not options.debug:
        shutil.rmtree(options.project_dir)

# ------------------------------------------------------------------------------

def check_path_string(path_str):
    print "\npath_str=", path_str
    if len(path_str) == 0:
        sys.exit("\nPlease, specify a path string with option -s")
    if path_str not in ["acquisition.timeseries", "acquisition.timeseries.extracellular_traces",\
                        "analysis", "epochs", "general", \
                        "processing", "stimulus", "top_datasets", \
                        "acquisition.images", "acquisition.timeseries",\
                        "general.top_datasets", "general.devices",\
                        "general.extracellular_ephys", "general.optogenetics",\
                        "general.optophysiology", "general.subject"] \
       and not re.search("acquisition.timeseries", path_str):
        sys.exit("\nPlease, specify a correct path string with option -s")

# ------------------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog data_h5 [meta_data_h5] [options (-h to list)]"

    parser = optparse.OptionParser(usage=usage)
    parser = make_nwb_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.verbose:
        print("len(args)=", len(args))

    if len(args) == 1 and os.path.isdir(args[0]):
        assemble_nwb_from_partial_files(args[0], options)
    elif len(args) in [1, 2]:
        if options.verbose:
            check_path_string(options.path_str)
        if len(options.project_dir) == 0:
            options.project_dir = os.path.basename(args[0]).split(".")[0]
        if not os.path.isdir(options.project_dir):
            os.mkdir(options.project_dir)
        data_path     = args[0]
        if not re.search(".h5", data_path):
            sys.exit("\nScript make_mwb.py accepts as input .h5 data file")
        metadata_path = ""
        if len(args) == 2:
            metadata_path = args[1]
            if not re.search(".h5", metadata_path):
                sys.exit("\nScript make_mwb.py accepts as input .h5 metadata file")
        
        if not os.path.isdir(args[0]) and len(options.path_str) == 0:
            create_and_assemble_all_partial_files(data_path, metadata_path, options) 
        else:
            options.data_path = data_path
            data_basename = os.path.basename(data_path)
            produce_partial_nwb(data_path, metadata_path, options)
    else:
        parser.print_usage()
        sys.exit(2)

