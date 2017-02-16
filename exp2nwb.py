#!/usr/local/python-2.7.11/bin/python
#
# make_nwb.py
#
# Created by Claudia Friedsam on 2015-02-08.
# Redesigned by Gennady Denisov on 2016-03-28.

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
import libh5 
import libexp2dict
import libdict2nwb

# ------------------------------------------------------------------------------

def make_nwb_command_line_parser(parser):
    parser.add_option("-a", "--run_all", action="store_true",dest="run_all", help="create/assemble all project files", default=False)
    parser.add_option("-d", "--project_dir", dest="project_dir", help="project directory name", default="")
    parser.add_option("-D", "--debug",   action="store_true",  dest="debug", help="output debugging info", default=False)
    parser.add_option("-e", "--no_error_handling", action="store_false", dest="handle_errors", help="handle_errors", default=True)
    parser.add_option("-p", "--pathstr", dest="path_str", help="dot-separated path string",metavar="path_str",default="")
    parser.add_option("-r", "--replace", action="store_true", dest="replace", help="if output file exists, replace it", default=False)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="increase the verbosity level", default=False)

    return parser

# ------------------------------------------------------------------------------

def set_data_origin(data_h5, meta_h5, options):
    options.data_origin = "Unknown"
    keys = [str(i) for i in libh5.get_key_list(data_h5["trialPropertiesHash"])]
    keys2 = libh5.get_key_list(data_h5["timeSeriesArrayHash"])
    if 'GoodTrials' in keys and not 'LickTime' in keys and not 'EphusVars' in keys2:
        options.data_origin = "DG"
    elif not "LickTime"   in libh5.get_key_list(data_h5["trialPropertiesHash"]) and \
       not "EphusVars"  in libh5.get_key_list(data_h5["timeSeriesArrayHash"]):
        options.data_origin = "SP"
    elif "EphusVars" in keys2:
        # NL data
        options.data_origin = "NL"
    elif "intracellular" in meta_h5.keys():
        # JY data
        options.data_origin = "JY"
    else:
        sys.exit("Data origin not determined")

    if options.verbose:
        print("\nData origin: " + options.data_origin + "\n")
    return options

# ------------------------------------------------------------------------------

def produce_partial_nwb(data_path, metadata_path, options):
    data_h5 = h5py.File(data_path, "r")
    if len(metadata_path) > 0:
        meta_h5 = h5py.File(metadata_path, "r")
    else:
        meta_h5 = data_h5

    options = set_data_origin(data_h5, meta_h5, options)

    # Extract (meta)data
    if options.path_str == "top_datasets":
        dict = libexp2dict.top_datasets(data_h5, meta_h5, options)

    elif options.path_str.split(".")[0] == "acquisition":
        items_acquisition = ["acquisition.images", "acquisition.timeseries"]
        if not re.search("acquisition.timeseries", options.path_str) and\
           not re.search("acquisition.images", options.path_str):
            sys.exit("\nUnsupported path string: " + options.path_str)
        elif re.search("acquisition.timeseries", options.path_str):
            dict = libexp2dict.acquisition_timeseries(data_h5, options)
        elif re.search("acquisition.images", options.path_str):
            dict = libexp2dict.acquisition_images(data_h5, options)

    elif options.path_str == "analysis":
        dict = libexp2dict.collect_analysis_information(data_h5, options)

    elif options.path_str.split(".")[0] == "epochs":
        dict = libexp2dict.create_epochs(data_h5, options)

    elif options.path_str.split(".")[0] == "general":

        items_general = ["general.top_datasets", "general.devices", \
                         "general.subject", "general.extracellular_ephys", \
                         "general.optogenetics", "general.optophysiology"]
               
        if not options.path_str in items_general:
            sys.exit("\nUnsupported path string: " + options.path_str)

        # Extract metadata
        if options.path_str == "general.devices":
            dict = libexp2dict.general_devices_data()
        elif options.path_str == "general.extracellular_ephys":
            dict = libexp2dict.general_extracellular_ephys_data(meta_h5, options)
        elif options.path_str == "general.optophysiology":
            dict = libexp2dict.general_optophysiology_data(data_h5, options)
        else:
            dict = libexp2dict.general(meta_h5, options)

    elif options.path_str.split(".")[0] == "processing":
        if options.path_str.split(".")[1] == "extracellular_units":
            dict = libexp2dict.processing_extracellular_units_data(data_h5, meta_h5, options)
        elif options.path_str.split(".")[1] == "ROIs":
            dict = libexp2dict.processing_ROIs(data_h5, options)
        elif options.path_str.split(".")[2] == "BehavioralEvents":
            # e.g. processing.Auditory.BehavioralEvents.reward_cue
            #      processing.Licks.BehavioralEvents.lick_left, processing.Licks.BehavioralEvents.lick_right
            #      processing.Pole.BehavioralEvents.pole_accessible
            #      processing.Reward.BehavioralEvents.water_left_reward, processing.Reward.BehavioralEvents.water_right_raward
            #      processing.Whisker.BehavioralEvents.pole_touch_protract, processing.Whisker.BehavioralEvents.pole_touch_retract
            module_name = options.path_str.split(".")[1]
            series_name = options.path_str.split(".")[-1]
            dict = libexp2dict.processing_events(module_name, series_name, data_h5, options)
        elif options.path_str.split(".")[2] == "BehavioralTimeSeries":
            # e.g. processing.Whisker.BehavioralTimeseries.whisker_angle, processing.Whisker.BehavioralTimeseries.whisker_curve
            module_name = options.path_str.split(".")[1]
            series_name = options.path_str.split(".")[-1]
            dict = libexp2dict.processing_timeseries(module_name, series_name, data_h5, options)
        else:
            sys.exit("\nUnsupported path string " + options.path_str + ". Exit")

    elif options.path_str.split(".")[0] == "stimulus":
        if not len(options.path_str.split(".")) == 3 or not options.path_str.split(".")[1] == "presentation":
            # e.g. stimulus.presentation.auditory_cue
            #      stimulus.presentation.pole_accessible
            #      stimulus.presentation.water_left
            #      stimulus.presentation.water_right
            #      stimulus.presentation.zaber_motor_pos
            #      stimulus.presentation.pole_in
            #      stimulus.presentation.pole_out
            #      stimulus.presentation.photostimulus
            sys.exit("\nUnsupported path string " + options.path_str + ". Exit")
        else:
            series_name = options.path_str.split(".")[-1]
            dict = libexp2dict.stimulus_presentation_timeseries(series_name, data_h5, meta_h5, options)
    else:
        sys.exit("\nUnsupported path string " + options.path_str)

    # Populate (meta)data
    if len(dict.keys()) == 0:
        sys.exit("Computing the dictionary is not currently implemented")

    libdict2nwb.make_h5(options.path_str, dict, options.project_dir)

# ------------------------------------------------------------------------------
#
# Recursively parse the group tree  in the H5 object
# and copy its contents to the NWB object
#
def update_nwb_object(last_item, curr_nwb_group, curr_h5_path, h5_object, nwb_object):   
    if hasattr(h5_object[curr_h5_path], '__dict__'):
        # Current path in h5_data points to a group    
        for k in h5_object[curr_h5_path].keys():
            curr_h5_path1 = curr_h5_path + "/" + str(k)
            num_keys = 0
            try:
                num_keys = len(h5_object[curr_h5_path1].keys())
            except:
                num_keys = 0
            if num_keys == 0:
                # Item is a dataset
                try:
                    curr_nwb_group.set_dataset(str(k), np.array(h5_object[curr_h5_path1]))
                except:
                    curr_nwb_group.set_custom_dataset(str(k), np.array(h5_object[curr_h5_path1]))
            else:
                # Next-level item is a group
#               assert(len(h5_object[curr_h5_path1].keys()) > 0)
                if not last_item == "top_datasets":
                    if re.search("shank", k):
                        curr_nwb_group1 = curr_nwb_group.make_group("<electrode_group_X>", str(k), abort=False)
                    elif re.search("site_", k):
                        curr_nwb_group1 = curr_nwb_group.make_group("<site_X>", name=str(k), abort=False)
                    elif re.search("trial_", k):
                        curr_nwb_group1 = curr_nwb_group.make_group("<epoch_X>", name=str(k), abort=False)
                    elif re.search("unit", k) and re.search("EventWaveform", curr_h5_path):
                        curr_nwb_group1 = curr_nwb_group.make_group("<SpikeEventSeries>", name=str(k), abort=False)
                    elif re.search("unit", k) and re.search("UnitTimes", curr_h5_path):
                        curr_nwb_group1 = curr_nwb_group.make_group("<unit_N>", name=str(k), abort=False)
                    elif re.search("fov_", k) and re.search("timeseries", curr_h5_path):
                        curr_nwb_group1 = nwb_object.make_group("<TwoPhotonSeries>", name=str(k),\
                                           path='/acquisition/timeseries')
                    elif re.search("fov_", k) and re.search("DfOverF", curr_h5_path):
                        curr_nwb_group1 = curr_nwb_group.make_group("<RoiResponseSeries>", name=str(k), abort=False)
                    elif re.search("fov_", k):
                        curr_nwb_group1 = curr_nwb_group.make_group("<imaging_plane_X>", name=str(k), abort=False)
                    elif re.search("channel_", k):
                        curr_nwb_group1 = curr_nwb_group.make_group("<channel_X>", name=str(k), abort=False)
                    else:
                        curr_nwb_group1 = curr_nwb_group.make_group(str(k), path = curr_h5_path, abort=False)
                    nwb_object = update_nwb_object(last_item, curr_nwb_group1, curr_h5_path1, h5_object, nwb_object)
    else:
        # Current path in h5_data points to a dataset
        path_items = curr_h5_path.split("/")
        upper_level_path = "/".join(path_items[0:-1])
        nwb_object.file_pointer[upper_level_path].create_dataset(k, data = np.array(h5_object[curr_h5_path]))
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
                                 "DfOverF", "ImageSegmentation"]:
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
    options = set_data_origin(data_h5, meta_h5, options)

    task_strings = ["top_datasets", \
                    "analysis",\
                    "analysis", \
                    "general.devices", \
                    "general.subject", \
                    "general.top_datasets"]
    if options.data_origin == "NL":
        task_strings = task_strings + \
                   ["acquisition.timeseries.extracellular_traces",\
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
    elif  options.data_origin == "SP":
        task_strings = task_strings + \
                   ["acquisition.images",\
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

    # Create partial NWB files
    task_strings = sorted(task_strings)
    for s in range(len(task_strings)):
#       print "\nstring= ", task_strings[s]
        command = "/home/denisovg/work/Karel_Svoboda/exp2nwb.py " + data_path + " " + metadata_path + " " 
        command += " -p " + task_strings[s] + " -d " + options.project_dir
        if options.verbose:
            print "Running command: " + command
        os.system(command)
                    
    # Assemble all
    print "\nStaring assembly of the partial files...\n"
    os.system("python exp2nwb.py " + options.project_dir)

# ------------------------------------------------------------------------------

def check_path_string(path_str):
    print "\npath_str=", path_str
    if len(path_str) == 0:
        sys.exit("\nPlease, specify a path string with option -p")
    if path_str not in ["acquisition.timeseries", "acquisition.timeseries.extracellular_traces",\
                        "analysis", "epochs", "general", \
                        "processing", "stimulus", "top_datasets", \
                        "acquisition.images", "acquisition.timeseries",\
                        "general.top_datasets", "general.devices",\
                        "general.extracellular_ephys", "general.optogenetics",\
                        "general.optophysiology", "general.subject"] \
       and not re.search("acquisition.timeseries", path_str):
        sys.exit("\nPlease, specify a correct path string with option -p")

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
        
        if options.run_all:
            create_and_assemble_all_partial_files(data_path, metadata_path, options) 
        else:
            options.data_path = data_path
            data_basename = os.path.basename(data_path)
            produce_partial_nwb(data_path, metadata_path, options)
    else:
        parser.print_usage()
        sys.exit(2)

