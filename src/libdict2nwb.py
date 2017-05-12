# Copyright (c) 2017, HHMI-Janelia Research Campus All rights reserved.
#
# libdict2nwb.py
#
# Producing a partial NWB file given a dictionary 
#
import numpy as np
import re
import os, sys
import nwb
from nwb import nwb_file
from nwb import nwb_utils
import h5py
import datetime
import unicodedata
import getpass
from sets import Set
import util

verbose = 1    
debug   = 0

class Dict2NWB(object):

    def initialize_h5_object(self, h5_file, dict):
        vargs = {"start_time" : 0.0, "description" : "default", \
                 "file_name" : h5_file, "identifier" : h5_file[0:-4], \
                 "mode" : "w", "verbosity" : "none"}
        for k in ["start_time", "description", "identifier"]:
            if k in dict.keys():
                vargs[k] = dict[k]

        if  os.path.exists(h5_file):
            os.remove(h5_file)   
        h5_object = nwb_file.open(**vargs)
        return h5_object

# ------------------------------------------------------------------------------

    def top_datasets(self, string, dict, project_dir):
        dict["verbosity"] = "none"
        if debug:   
            util.print_dict_keys(dict)       
        h5_object = nwb_file.open(**dict)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def acquisition_images(self, string, dict, project_dir):
        # string = acquisition.images
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        ag  = h5_object.make_group("acquisition", abort=False)
        aig = ag.make_group("images", abort=False)
        for k in  dict.keys():
            if re.search("description", k):
                try:
                    aig.set_dataset("description", dict[k])
                except:
                    aig.set_custom_dataset("description", dict[k])
            elif not k.split(".")[-1] == "attrs":
                name = k.split(".")[-1]
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if re.search("fov_", k):
                    h5_object.set_dataset("<image_X>", dict[k], name=name,\
                                          dtype='uint8', attrs=attrs)
                else:
                    try:
                        aig.set_dataset(name, dict[k], attrs=attrs)
                    except:
                        aig.set_custom_dataset(name, dict[k], attrs=attrs)
                
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def acquisition_timeseries(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        if len(string.split(".")) == 3:
            # behavioral time series
            # string = acquisition.timeseries.<data_name>
            series_name = string.split(".")[-1]
            if string + ".attrs" in dict.keys() and \
               "series_type" in dict[string + ".attrs"].keys():
                series_type = dict[string + ".attrs"]["series_type"]
            else:   
                series_type = '<TimeSeries>'
            series_path = "/acquisition/timeseries"
            attrs = {"series_type" : "<TimeSeries>"}
            if string + ".attrs" in dict.keys():
                attrs = dict[string + ".attrs"]
            series_type = '<TimeSeries>'
            tsg = h5_object.make_group(series_type, series_name, path = series_path, attrs = attrs)
            for k in dict.keys():
                if len(k.split(".")) == len(string.split(".")) + 1:
                    name = k.split(".")[-1]
                    attrs = {}
                    if k + ".attrs" in dict.keys():
                        attrs = dict[k + ".attrs"]
                    if not name == "attrs":
                        try: 
                            tsg.set_dataset(name, dict[k], attrs = attrs)
                        except:
                            tsg.set_custom_dataset(name, dict[k], attrs = attrs)
        elif len(string.split(".")) == 2:
            # image time series - ophys data only     
            # string = acquisition.timeseries
            img_planes = []
            for k in dict.keys():
                if k.split(".")[-1] == "imaging_plane":
                    img_planes.append(k.split(".")[2])
            for img_plane in img_planes:
                group_attrs = {}
                if string + "." + img_plane + ".attrs" in dict.keys():
                    group_attrs = dict[string + ".attrs"]
                twop = h5_object.make_group("<TwoPhotonSeries>", img_plane, \
                           path='/acquisition/timeseries', attrs=group_attrs)
                for k in dict.keys():
                    if k.split(".")[2] == img_plane and not k.split(".")[-1] == "attrs":
                        name = k.split(".")[-1]
                        attrs = {}
                        if k + ".attrs" in dict.keys():
                            attrs = dict[k + ".attrs"]
                        try:
                            twop.set_dataset(name, dict[k], attrs = attrs)
                        except:
                            try:
                                twop.set_custom_dataset(name, dict[k], attrs = attrs)
                            except:
                                pass
        else:
            sys.exit("Unsupported path string " + string)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def acquisition_timeseries_extracellular_traces(self, string, dict, project_dir):
        # string = acquisition.timeseries.extracellular_traces
        if debug:
            util.print_dict_keys(dict)          

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        group_attrs = {} 
        if "acquisition.timeseries.extracellular_traces.attrs" in dict.keys():
            group_attrs = dict["acquisition.timeseries.extracellular_traces.attrs"]
        et = h5_object.make_custom_group("extracellular_traces", path = "/acquisition/timeseries", \
                                  attrs = group_attrs)
        for k in  dict.keys():
            if len(k.split(".")) == len(string.split(".")) + 1 and \
                   not k.split(".")[-1] == "attrs":
                name = k.split(".")[-1]
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if name == "ephys_raw_data":
                    et.set_custom_dataset(name, dict[k], attrs = attrs)
                else:
                    try:
                        et.set_dataset(name, dict[k], attrs = attrs)
                    except:
                        et.set_custom_dataset(name, dict[k], attrs = attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def analysis(self, string, dict, project_dir):
        # string = analysis
        if debug:
            util.print_dict_keys(dict)          
        analysis_datasets = ["description", "good_trials", "trial_start_times",\
                             "trial_type_mat", "trial_type_string"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        ag = h5_object.make_group("analysis", abort=False)
        path = "/" + ag.name
        for k in dict.keys():
            if not k.split(".")[-1] == "attrs":
                name = k.split(".")[-1]
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                try:
                    ag.set_dataset(name, dict[k], attrs = attrs)
                except:
                    ag.set_custom_dataset(name, dict[k], attrs = attrs)
                if not name in analysis_datasets:
                    print "    Warning: creating custom dataset analysis." + name
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def epochs(self, string, dict, project_dir):
        # string = epochs
        # no custom datasets are allowed
        if debug:
            util.print_dict_keys(dict)       
        epoch_datasets = ["start_time", "stop_time", "tags", "units_present", \
                          "description", "ROIs", "ROI_planes"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        eg = h5_object.make_group("epochs", abort=False)
        processed_trials = []
        tg = {}
        for k in dict.keys():
            if len(k.split(".")) == 2:
                name = k.split(".")[-1]
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if not re.search("rial", name):
                    print "    Warning: creating custom dataset epochs." + name
                    try:
                        eg.set_dataset(name, dict[k], attrs = attrs)
                    except:
                        eg.set_custom_dataset(name, dict[k], attrs = attrs)
            else:
                trial = k.split(".")[-2]
                if trial not in processed_trials:
                    processed_trials.append(trial)
                    tg[trial] = h5_object.make_group("<epoch_X>", trial, attrs = \
                            {"description" : "Data that belong to " + trial,\
                             "ancestry" : np.array(["epoch_X"])})
                if len(k.split(".")) == 3 and not k.split(".")[-1] == "attrs":
                    name = k.split(".")[-1]
                    attrs = {}
                    if k + ".attrs" in dict.keys():
                        attrs = dict[k + ".attrs"]
                    try:
                        tg[trial].set_dataset(name, dict[k], attrs = attrs)
                    except:
                        tg[trial].set_custom_dataset(name, dict[k], attrs = attrs)
                    if not name in epoch_datasets:
                        print "    Warning: creating custom dataset epochs." + trial + "." + name
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def general_top_datasets(self, string, dict, project_dir):
        # string = general.top_datasets
        if debug:
            util.print_dict_keys(dict)       
        top_datasets = ["data_collection", "experiment_description", "experimenter", \
                        "institution", "lab", "notes", "pharmacology", "protocol", \
                        "reference_atlas", "related_publications", "session_id", \
                        "slices", "source_script", "stimulus", "surgery", \
                        "task_keywords", "whisker_configuration", "virus"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        gg = h5_object.make_group("general", abort=False)
        for k in dict.keys():
            name = k.split(".")[len(k.split("."))-1]
            if not name in top_datasets:
                print "    Warning: creating custom dataset general." + name
            attrs = {}
            if k + ".attrs" in dict.keys():
                attrs = dict[k + ".attrs"]
            try:
                gg.set_dataset(name, dict[k], attrs = attrs)
            except:
                gg.set_custom_dataset(name, dict[k], attrs = attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def general_devices(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        for k in dict.keys():
            name = k.split(".")[-1]
            if not name == "attrs":
                attrs = attrs = {"ancestry" : np.array(["device_X"])}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                    attrs["ancestry"] = np.array(["device_X"])
                h5_object.set_dataset("<device_X>", dict[k], name=name, \
                                      attrs=attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def general_extracellular_ephys(self, string, dict, project_dir):
        # string = general.extracellular_ephys
        if debug:
            util.print_dict_keys(dict)       
        extracellular_ephys_datasets = ["ADunit", "electrode_group", "electrode_map", \
                                        "filtering"]
        shank_datasets = ["description", "device", "location"]
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)

        gg = h5_object.make_group("general", abort=False)
        eeg = gg.make_group("extracellular_ephys", abort=False)
        for k in dict.keys():
            if k.split(".")[-2] == "extracellular_ephys": 
                data_name = k.split(".")[-1]
                if not data_name == "attrs":
                    attrs = {}
                    if k + ".attrs" in dict.keys():
                        attrs = dict[k + ".attrs"]
                    if not data_name in extracellular_ephys_datasets:
                        print "    Warning: creating custom dataset general.extracellular_ephys." + data_name
                    try:
                        eeg.set_dataset(data_name, dict[k], attrs=attrs)
                    except:
                        eeg.set_custom_dataset(data_name, dict[k], attrs=attrs)
            elif re.search("shank", k.split(".")[-2]):
                sg = eeg.make_group("<electrode_group_X>", k.split(".")[-2], \
                                    attrs={"ancestry" : np.array(["electrode_group_X"])}, abort=False)
                data_name = k.split(".")[-1]
                if not data_name == "attrs":
                    attrs = {}
                    if k + ".attrs" in dict.keys():
                        attrs = dict[k + ".attrs"]
                    if not data_name in shank_datasets:
                        print "    Warning: creating custom dataset general.extracellular_ephys.shank." + data_name
                    try:
                        sg.set_dataset(data_name, dict[k], attrs=attrs)
                    except:
                        sg.set_custom_dataset(data_name, dict[k], attrs=attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def general_subject(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        subject_datasets = ["age", "description", "genotype", "sex", "species", \
                            "subject_id", "weight"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        gg = h5_object.make_group("general", abort=False)
        sg = gg.make_group("subject", abort=False)
        for k in dict.keys():
            data_name = k.split(".")[-1]
            if not data_name == "attrs":
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if not data_name in subject_datasets:
                    print "    Warning: creating custom dataset general.subject." + data_name
                try:
                    sg.set_dataset(data_name, dict[k], attrs=attrs)
                except:
                    sg.set_custom_dataset(data_name, dict[k], attrs=attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def general_optogenetics(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        optogenetics_datasets = ["device", "excitation_lambda", "location", "stimulation_method"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        gg = h5_object.make_group("general", abort=False)
        og = gg.make_group("optogenetics", abort=False)
        for k in dict.keys():
            data_name = k.split(".")[-1]
            if len(k.split(".")) == 3 and not data_name == "attrs":
                # Creating custom datasets in group general/optogenetics
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if not data_name in optogenetics_datasets:
                    print "    Warning: creating custom dataset general.optogenetics." + data_name
                try:
                    og.set_dataset(data_name, dict[k], attrs=attrs)
                except:
                    og.set_custom_dataset(data_name, dict[k], attrs=attrs)

            if k.split(".")[-3] == "optogenetics":
                osg = og.make_group("<site_X>", name=k.split(".")[-2], \
                                    attrs={"ancestry" : np.array(["site_X"])}, abort=False)
                if len(k.split(".")) == 4 and not data_name == "attrs":
                     # Creating datasets in group general/optogenetics/site_#
                    attrs = {}
                    if k + ".attrs" in dict.keys():
                        attrs = dict[k + ".attrs"]
                    if not data_name in optogenetics_datasets:
                        print "    Warning: creating custom dataset general.optogenetics." \
                              + k.split(".")[-2] + "." + data_name
                    try:
                        osg.set_dataset(data_name, dict[k], attrs=attrs)
                    except:
                        osg.set_custom_dataset(data_name, dict[k], attrs=attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def general_optophysiology(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        gg = h5_object.make_group("general", abort=False)
        op = gg.make_group("optophysiology", abort=False)
        for k in dict.keys():
            data_name = k.split(".")[-1]
            if not data_name == "attrs":
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if not re.search("fov_", k) and len(k.split(".")) == 3:
                    # Creating description and custom datasets in general/optophysiology
                    try:
                        op.set_dataset(data_name, dict[k], attrs = attrs)
                    except:
                        op.set_custom_dataset(data_name, dict[k], attrs = attrs)
                elif re.search("fov_", k):
                    group_name = k.split(".")[2]
                    img_plane_group = op.make_group("<imaging_plane_X>", group_name, \
                                                    attrs = {"ancestry" : np.array(["imaging_plane_X"])}, abort=False)
                    if not re.search("channel_", k):
                        try:
                            img_plane_group.set_dataset(data_name, dict[k], attrs = attrs)
                        except:
                            img_plane_group.set_custom_dataset(data_name, dict[k], attrs = attrs)
                    else:
                        name_channel = k.split(".")[3]
                        channel_group = img_plane_group.make_group("<channel_X>", name_channel, \
                                            attrs = {"ancestry" : np.array(["channel_X"])}, abort=False)
                        try:
                            channel_group.set_dataset(data_name, dict[k], attrs = attrs)
                        except:
                            channel_group.set_custom_dataset(data_name, dict[k], attrs = attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_extracellular_units(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        if re.search("top_datasets", string):
            self.processing_extracellular_units_top_datasets(string, dict, project_dir)
        elif re.search("EventWaveform", string):
            self.processing_extracellular_units_EventWaveform(string, dict, project_dir)
        elif re.search("UnitTimes", string):
            self.processing_extracellular_units_UnitTimes(string, dict, project_dir)

# ------------------------------------------------------------------------------

    def processing_extracellular_units_UnitTimes(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        module_name = "extracellular_units"
        module_attrs = {"ancestry" : "<Module>"}
        if "processing.extracellular_units.attrs" in dict.keys():
            module_attrs = dict["processing.extracellular_units.attrs"]
            module_attrs["ancestry"] = "<Module>"
        mod = h5_object.make_group("<Module>", module_name, attrs=module_attrs)
        group_attrs = {}
        if "processing.extracellular_units.UnitTimes.attrs" in dict.keys():
            group_attrs = dict["processing.extracellular_units.UnitTimes.attrs"]
        spk_times_iface = mod.make_group("UnitTimes", attrs=group_attrs)

        spk = {}
        for k in dict.keys():
            data_name = k.split(".")[-1]
            if not data_name == "attrs":
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if len(k.split(".")) == 4:
                    # Creating dataset cell_types, etc
                    try:
                        spk_times_iface.set_dataset(data_name, dict[k], attrs=attrs)
                    except:
                        spk_times_iface.set_custom_dataset(data_name, dict[k], attrs=attrs)
                elif len(k.split(".")) == 5:
                    unit = k.split(".")[-2]
                    attrs["ancestry"] = np.array(["<unit_N>"])
                    spk[unit] = spk_times_iface.make_group("<unit_N>", unit, \
                                                           attrs=attrs,abort = False)
                    try:
                        spk[unit].set_dataset(data_name, dict[k], attrs = attrs)
                    except:
                        spk[unit].set_custom_dataset(data_name, dict[k], attrs = attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_extracellular_units_EventWaveform(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        module_name = "extracellular_units"
        module_attrs={"ancestry" : np.array(["Module"])}
        if "processing.extracellular_units.attrs" in dict.keys():
            module_attrs = dict["processing.extracellular_units.attrs"]
            module_attrs["ancestry"] = "extracellular_units"
        mod = h5_object.make_group("<Module>", module_name, attrs=module_attrs)
        group_attrs = {}
        if "processing.extracellular_units.EventWaveform.attrs" in dict.keys():
            group_attrs = dict["processing.extracellular_units.EventWaveform.attrs"]
        spk_waves_iface = mod.make_group("EventWaveform")
        spk_waves_iface.set_attr("source", "EventWaveform in this module")

        spk = {}
        for k in dict.keys():
            data_name = k.split(".")[-1]
            if not data_name == "attrs":
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if len( k.split(".")) == 5:
                    unit = k.split(".")[-2]
                    group_attrs = {"source" : "---"}
                    spk[unit] = spk_waves_iface.make_group("<SpikeEventSeries>", unit, \
                                    attrs = group_attrs, abort = False)
                try:
                    spk[unit] = spk[unit].set_dataset(data_name, dict[k], attrs = attrs)
                except:
                    spk[unit] = spk[unit].set_custom_dataset(data_name, dict[k], attrs = attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_extracellular_units_top_datasets(self, string, dict, project_dir):
        print "dict.keys()=", dict.keys()
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        module_name = "extracellular_units"
        mod = h5_object.make_group("<Module>", module_name, attrs={"ancestry" : np.array(["Module>"])})

        for k in dict.keys():
            data_name = k.split(".")[-1]
            if not data_name == "attrs":
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if len( k.split(".")) == 3:
                    try:
                        mod.set_dataset(data_name, dict[k], attrsc=cattrs)
                    except:
                        mod.set_custom_dataset(data_name, dict[k], attrs=attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_timeseries(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)  
        module_name = string.split(".")[1]
        series_name = string.split(".")[-1]
        if string + ".attrs" in dict.keys() and \
           "series_type" in dict[string + ".attrs"].keys():
            series_type = dict[string + ".attrs"]["series_type"]
        else:
            series_type = "<TimeSeries>"
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        try:
            module_attrs = dict["processing." + module_name + ".attrs"]
        except:
            module_attrs = {}
        module_attrs["ancestry"] = np.array(["Module"])
        mod = h5_object.make_group("<Module>", module_name, attrs = module_attrs)

        iface = mod.make_group(string.split(".")[2], attrs=dict[string + ".attrs"])
        tsg   = iface.make_group(series_type, series_name, attrs = dict[string + ".attrs"])

        for k in dict.keys():
            data_name = k.split(".")[-1]
            if not data_name == "attrs" and len(k.split(".")) == 5:
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                try:
                    tsg.set_dataset(data_name, dict[k], attrs=attrs)
                except:
                    tsg.set_custom_dataset(data_name, dict[k], attrs=attrs)

        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_ROIs(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        h5_object = self.processing_ROIs_DfOverF(h5_object, string + ".DfOverF", dict, project_dir)
        dict2 = {}
        for k in dict.keys():
            if re.search("ImageSegmentation", k):
                dict2[k] = dict[k]
        h5_object = self.processing_ROIs_ImageSegmentation(h5_object, string + ".ImageSegmentation", \
                         dict2, project_dir)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_ROIs_DfOverF(self, h5_object, string, dict, project_dir):
        module_name = string.split(".")[1]        
        iface_name  = string.split(".")[2]

        mod = h5_object.make_group("<Module>", module_name, \
                                   attrs ={"ancestry" : np.array(["Module"])})
        img_plane_group = {}
        for k in dict.keys():
            data_name = k.split(".")[-1]
            if not data_name == "attrs":
                attrs = {}
                if k + ".attrs" in dict.keys():
                    attrs = dict[k + ".attrs"]
                if len(k.split(".")) == 3:
                    try:
                        mod.set_dataset(data_name, dict[k], attrs = attrs) # description etc.
                    except:
                        mod.set_custom_dataset(data_name, dict[k], attrs = attrs)
                elif len(k.split(".")) == 5 and re.search("fov_", k) and \
                     re.search("DfOverF", k):
                    group_attrs  = dict[string + ".attrs"]
                    dff_iface = mod.make_group("DfOverF", attrs=group_attrs, abort = False)
                    string2 = "processing.ROIs.ImageSegmentation"
                    group_attrs2 = dict[string2 + ".attrs"]
                    seg_iface = mod.make_group("ImageSegmentation", attrs=group_attrs2, abort = False)
                    img_plane = k.split(".")[3]
                    img_plane_group[img_plane] = dff_iface.make_group("<RoiResponseSeries>", img_plane, \
                                                    abort=False)
                    img_plane_group[img_plane].make_group("segmentation_interface", seg_iface, abort = False)
                    try:
                        img_plane_group[img_plane].set_dataset(data_name, dict[k], attrs = attrs)
                    except:
                        img_plane_group[img_plane].set_custom_dataset(data_name, dict[k], attrs = attrs)
        return h5_object

# ------------------------------------------------------------------------------

    def processing_ROIs_ImageSegmentation(self, h5_object, string, dict, project_dir):
        module_name = string.split(".")[1]
        iface_name  = string.split(".")[2]

        mod = h5_object.make_group("<Module>", module_name, attrs ={"ancestry" : "<Module>"}, abort=False)
        seg_iface = mod.make_group("ImageSegmentation", attrs=dict[string + ".attrs"], abort=False)
        processed_img_planes = []
        img_plane_group = {}
        for k in dict.keys():
            if re.search("fov_", k):
                img_plane = k.split(".")[3]
                img_plane_group[img_plane] = seg_iface.make_group("<image_plane>", img_plane, \
                                                    attrs = {"ancestry" : np.array(["image_plane"])}, abort=False)
                data_name = k.split(".")[-1]
                if not data_name == "attrs" and not re.search("ref_im", k):
                    attrs = {}
                    if k + ".attrs" in dict.keys():
                        attrs = dict[k + ".attrs"]
                    if len(k.split(".")) == 5:
                        try:
                            img_plane_group[img_plane].set_dataset(data_name, dict[k], attrs = attrs)
                        except:
                            img_plane_group[img_plane].set_custom_dataset(data_name, dict[k], attrs = attrs)

                if img_plane in processed_img_planes:
                    continue
                processed_img_planes.append(img_plane)
                path0 = string + "." + img_plane
                try:
                    nwb_utils.add_reference_image(seg_iface, img_plane, "%s_0001"%img_plane, \
                                                  dict[path0 + ".ref_image_green"], dict[path0 + ".ref_image_green.attrs"]["source"])
                except:
                    print("    Warning: cannot store red reference image")
                try:
                    nwb_utils.add_reference_image(seg_iface, img_plane, "%s_0002"%img_plane, \
                                                  dict[path0 + ".ref_image_red"], dict[path0 + ".ref_image_red.attrs"]["source"])
                except:
                    print("    Warning: cannot store green reference image")

                roi_ids = dict[path0 + ".roi_list"]
                for i in range(len(roi_ids)):
                    roi_id = roi_ids[i]
                    x      = dict[path0 + "." + str(roi_id) + ".x"]
                    weight = dict[path0 + "." + str(roi_id) + ".weight"]
                    pixmap = dict[path0 + "." + str(roi_id) + ".pixmap"]
                    master1_shape = dict[path0 + "." + str(roi_id) + ".master1_shape"]
#                   print "img_plane=", img_plane, " roi_id=", roi_id, " x=", x, " weight=", weight, " master1_shape=", master1_shape
                    nwb_utils.add_roi_mask_pixels(seg_iface, img_plane, "%d"%x, "ROI %d"%x, np.array(pixmap).astype('uint16'), \
                                                  weight, master1_shape[1], master1_shape[0])
        return h5_object

# ------------------------------------------------------------------------------

    def stimulus_presentation(self, string, dict, project_dir):
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        if debug:
            util.print_dict_keys(dict)       
        series_path = "/stimulus/presentation"
        for k in dict.keys():
            if len(k.split(".")) == 4:
                data_name = k.split(".")[-1]
                if not data_name == "attrs":
                    group_name = k.split(".")[-2]
                    # We assume that all groups have attributes
                    group_attrs = dict[".".join([k.split(".")[0], k.split(".")[1], group_name, "attrs"])]
                    series_type = group_attrs["series_type"]
                    tsg = h5_object.make_group(series_type, group_name, path = series_path, \
                                               attrs = group_attrs, abort=False)
                    data_attrs = {}
                    if k + ".attrs" in dict.keys():
                        data_attrs = dict[k + ".attrs"]
                    try:
                        tsg.set_dataset(data_name, dict[k], attrs = data_attrs)
                    except:
                        tsg.set_custom_dataset(data_name, dict[k], attrs = data_attrs)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

def make_partial_nwb(string, dict, project_dir):
    method_name = "unknown"
    if re.search("acquisition.timeseries", string) and \
       not string == "acquisition.timeseries.extracellular_traces":
       method_name = "acquisition_timeseries"
    elif re.search("processing", string) and \
         re.search("BehavioralTimeSeries", string) or re.search("BehavioralEvents", string):
            method_name = "processing_timeseries"
    elif re.search("stimulus.presentation", string):
       method_name = "stimulus_presentation"
    else:
        if len(string.split(".")) > 0:
            method_name = "_".join(string.split("."))
        else:
            method_name = string
    class_instance = Dict2NWB()
#   print "\nmethod_name=", method_name

    method = None
    try:
        method = getattr(class_instance, method_name)
    except:
        sys.exit("\nMethod libdict2nwb." + method_name + " has not been implemented")

    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
    method(string, dict, project_dir)

