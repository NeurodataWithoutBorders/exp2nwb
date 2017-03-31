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
        if debug:   
            util.print_dict_keys(dict)       
        h5_object = nwb_file.open(**dict)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def acquisition_images(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        acquisition_group     = h5_object.make_group("acquisition", abort=False)
        acquisition_img_group = acquisition_group.make_group("images", abort=False)
        for k in  dict.keys():
            if re.search("description", k):
                acquisition_img_group.set_custom_dataset("description", dict[k])
            elif not k.split(".")[-1] == "attrs":
                name = k.split(".")[-1]
                attrs = attrs=dict[k + ".attrs"]
                attrs["ancestry"] = np.array(["image_X"])
                h5_object.set_dataset("<image_X>", dict[k], name=name,\
                           dtype='uint8', attrs=attrs)
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
            series_name = string.split(".")[-1]
            series_type = dict[string + ".attrs"]["series_type"]
            series_path = "/acquisition/timeseries"
            tsg = h5_object.make_group(series_type, series_name, path = series_path, \
                              attrs = dict[string + ".attrs"])
            for k in dict.keys():
                if re.search("attr", k):
                    continue
                dataset_name = k.split(".")[-1]
                if dataset_name == "data":
                    tsg.set_dataset(dataset_name, dict[string + "." + dataset_name], attrs=dict[string + ".data_attrs"])
                elif dataset_name in ["timestamps", "num_samples"]:
                    tsg.set_dataset(dataset_name, dict[string + "." + dataset_name])
                else:
                    tsg.set_custom_dataset(dataset_name, dict[string + "." + dataset_name])
        elif len(string.split(".")) == 2:
            # image time series - ophys data only     
            path0 = "acquisition.timeseries."
            img_planes = []
            for k in dict.keys():
                if k.split(".")[-1] == "imaging_plane":
                    img_planes.append(k.split(".")[2])
            for img_plane in img_planes:
                group_attrs = dict[path0 + "attrs"]
                twop = h5_object.make_group("<TwoPhotonSeries>", img_plane, \
                           path='/acquisition/timeseries', attrs=group_attrs)
                for item in ["format", "external_file", "dimension", "scan_line_rate",\
                             "field_of_view", "imaging_plane", "timestamps"]:
                    if item == "external_file":
                        twop.set_dataset(item, \
                                  dict[path0 + img_plane +"." + item],\
                            attrs=dict[path0 + img_plane +"." + item + ".attrs"])
                    else:
                        twop.set_dataset(item, dict[path0 + img_plane +"." + item])
        else:
            sys.exit("Unsupported path string " + string)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def acquisition_timeseries_extracellular_traces(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)          
        datasets = ["ephys_raw_data"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        et = h5_object.make_group("<TimeSeries>", "extracellular_traces", path = "/acquisition/timeseries", \
             attrs = {"description" : "File containing the raw voltage data for this session"})
        for k in  dict.keys():
            if re.search("ephys_raw_data", k):
                et.set_dataset("ephys_raw_data", dict[k])
            elif not re.search("attr", k):
                name = k.split(".")[-1]
                et.set_custom_dataset(name, dict[k])
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def analysis(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)          
        analysis_datasets = ["description", "good_trials", "trial_start_times",\
                             "trial_type_mat", "trial_type_string"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        ag = h5_object.make_group("analysis", abort=False)
        for k in dict.keys():
            k1 = k.split(".")[len(k.split("."))-1]
            if k in analysis_datasets or k1 in analysis_datasets:
                try:
                    ag.set_dataset(k1, dict[k])
                except:
                    ag.set_custom_dataset(k1, dict[k])
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def epochs(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        epoch_datasets = ["start_time", "stop_time", "tags", "units_present", \
                          "description", "ROIs", "ROI_planes"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        eg = h5_object.make_group("epochs", abort=False)
        processed_trials = []
        for k in dict.keys():
            trial = k.split(".")[-2]
            if trial not in processed_trials:
                processed_trials.append(trial)
                start = dict["epochs." + trial + "." + "start_time"]
                stop  = dict["epochs." + trial + "." + "stop_time"]
#               epoch = nwb_utils.create_epoch(h5_object, trial, start, stop)
                epoch = h5_object.make_group("<epoch_X>", trial, attrs = \
                        {"description" : "Data that belong to " + trial,\
                         "ancestry" : np.array(["epoch_X"])})
                epoch.set_dataset("start_time", start)
                epoch.set_dataset("stop_time",  stop)
                for k1 in epoch_datasets:
                    k2 = "epochs." + trial + "." + k1
                    if k2 not in dict.keys():
                        continue
                    if k1 in ["tags", "description"]:
                        epoch.set_dataset(k1, dict[k2])
                    elif k1 in ["units_present", "ROIs", "ROI_planes"]:
                        epoch.set_custom_dataset(k1, dict[k2])
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def general_top_datasets(self, string, dict, project_dir):
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
            k1 = k.split(".")[len(k.split("."))-1]
            if k in top_datasets or k1 in top_datasets:
                try:
                    gg.set_dataset(k1, dict[k])
                except:
                    gg.set_custom_dataset(k1, dict[k])
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def general_devices(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        device_datasets = ["ephys_acquisition", "photostim_source"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        for k in dict.keys():
            h5_object.set_dataset("<device_X>", dict[k], name=k, \
                                  attrs={"ancestry" : np.array(["device_X"])})
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def general_extracellular_ephys(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        extracellular_ephys_datasets = ["ADunit", "electrode_group", "electrode_map", \
                                        "filtering"]
        shank_datasets = ["description", "device", "location"]
        
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)

        gg = h5_object.make_group("general", abort=False)
        eeg = gg.make_group("extracellular_ephys", abort=False)
        for k in dict.keys():
            if k.split(".")[-1] in extracellular_ephys_datasets or\
               k.split(".")[-2] == "extracellular_ephys":
               if verbose:
                   print "Creating dataset ", k.split(".")[-1], " value=", dict[k]
               eeg.set_dataset(k.split(".")[-1], dict[k])
                
            elif re.search("shank", k.split(".")[-2]):
                if verbose:
                    print "Creating group ", k.split(".")[-2]
                sg = eeg.make_group("<electrode_group_X>", k.split(".")[-2], \
                                    attrs={"ancestry" : np.array(["electrode_group_X"])}, abort=False)
                if k.split(".")[-1] in shank_datasets:
                    if verbose:
                        print "Creating dataset ", k.split(".")[-1], " value=", dict[k]
                    sg.set_dataset(k.split(".")[-1], dict[k])
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
            if k in subject_datasets or k.split(".")[-1] in subject_datasets:
                try:
                    sg.set_dataset(k.split(".")[-1], dict[k])
                except:
                    sg.set_custom_dataset(k.split(".")[-1], dict[k])
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
            if k.split(".")[-1] in optogenetics_datasets and k.split(".")[-3] == "optogenetics":
                osg = og.make_group("<site_X>", name=k.split(".")[-2], \
                                    attrs={"ancestry" : np.array(["site_X"])}, abort=False)
                try:
                    osg.set_dataset(k.split(".")[-1], dict[k])
                except:
                    osg.set_custom_dataset(k.split(".")[-1], dict[k])
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
            if not re.search("fov_", k) and len(k.split(".")) == 3 and re.search("description", k):
                op.set_custom_dataset("description", dict[k])
            elif re.search("fov_", k) and not re.search("channel_", k):
                name_group = k.split(".")[2]
                img_plane_group = op.make_group("<imaging_plane_X>", name_group, \
                                                attrs = {"ancestry" : np.array(["imaging_plane_X"])}, abort=False)
                name_data = k.split(".")[3]
                if not re.search("attrs", k):
                    if name_data == "manifold":
                        img_plane_group.set_dataset(name_data, dict[k], attrs = dict[k + ".attrs"])
                    else:
                        img_plane_group.set_dataset(name_data, dict[k])
            elif re.search("fov_", k) and re.search("channel_", k):
                name_group = k.split(".")[2]
                img_plane_group = op.make_group("<imaging_plane_X>", name_group, \
                                      attrs = {"ancestry" : np.array(["imaging_plane_X"])}, abort=False)
                name_channel = k.split(".")[3]
                channel_group = img_plane_group.make_group("<channel_X>", name_channel, \
                                    attrs = {"ancestry" : np.array(["channel_X"])}, abort=False)
                name_data = k.split(".")[4]
                channel_group.set_dataset(name_data, dict[k])
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
        mod = h5_object.make_group("<Module>", module_name, attrs={"ancestry" : "<Module>"})

        spk_times_iface = mod.make_group("UnitTimes")
        spk_times_iface.set_attr("source", dict[string + ".attrs"]["source"])

        processed_units = []

        for k in dict.keys():
            unit = k.split(".")[-2]
            if unit not in processed_units \
                and not re.search("attr", k)  \
                and not re.search("cell_types", k):             
                processed_units.append(unit)
                path = string + "." + unit
                spk = spk_times_iface.make_group("<unit_N>", unit, \
                                                 attrs={"ancestry" : np.array(["<unit_N>"])})
                spk.set_custom_dataset("times",            dict[path + ".times"])
                spk.set_custom_dataset("source",           dict[path + ".source"])
                spk.set_custom_dataset("trial_ids",        dict[path + ".trial_ids"])
                spk.set_custom_dataset("unit_description", dict[path + ".unit_description"])
        spk_times_iface.set_custom_dataset("cell_types",   dict[string + ".cell_types"])
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_extracellular_units_EventWaveform(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        module_name = "extracellular_units"
        mod = h5_object.make_group("<Module>", module_name, attrs={"ancestry" : np.array(["Module"])})

        spk_waves_iface = mod.make_group("EventWaveform")
        spk_waves_iface.set_attr("source", dict[string + ".attrs"]["source"])

        processed_units = []

        for k in dict.keys():
            unit = k.split(".")[-2]
            if not re.search("attr", k) and unit not in processed_units:
                processed_units.append(unit)
                path = string + "." + unit
                spk = spk_waves_iface.make_group("<SpikeEventSeries>", unit)
                spk.set_attr("source", "---")
                spk.set_custom_dataset("sample_length", dict[path + ".sample_length"])
                spk.set_custom_dataset("timestamps",    dict[path + ".timestamps"])
                spk.set_custom_dataset("data",          dict[path + ".data"])
                spk.set_custom_dataset("electrode_idx", dict[path + ".electrode_idx"])
                spk.set_custom_dataset("num_samples",   dict[path + ".num_samples"])

        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_extracellular_units_top_datasets(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        module_name = "extracellular_units"
        mod = h5_object.make_group("<Module>", module_name, attrs={"ancestry" : np.array(["Module>"])})

        path = ".".join(string.split(".")[0:-1])
        mod.set_custom_dataset('description',           dict[path + ".description"])
        mod.set_custom_dataset('spike_sorting',         dict[path + ".spike_sorting"])
        mod.set_custom_dataset('identification_method', dict[path + ".identification_method"])

        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_timeseries(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        module_name = string.split(".")[1]
        series_name = string.split(".")[-1]
        series_type = dict[string + ".attrs"]["series_type"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        try:
            module_attrs = dict["processing." + module_name + ".attrs"]
        except:
            module_attrs = {}
        module_attrs["ancestry"] = np.array(["Module"])
        mod = h5_object.make_group("<Module>", module_name, attrs = module_attrs)

        iface = mod.make_group(string.split(".")[2], attrs=dict[string + ".attrs"])
        tsg   = iface.make_group(series_type, series_name, attrs = dict[string + ".attrs"])

        tsg.set_dataset("data", dict[string + ".data"], attrs=dict[string + ".data.attrs"])
        tsg.set_dataset("timestamps", dict[string + ".timestamps"])
        tsg.set_dataset("num_samples", len(dict[string + ".timestamps"]))

        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_ROIs(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        h5_object = self.processing_ROIs_DfOverF(h5_object, string + ".DfOverF", dict, project_dir)
        h5_object = self.processing_ROIs_ImageSegmentation(h5_object, string + ".ImageSegmentation", dict, project_dir)
        if verbose:
            print "Creating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close() 

# ------------------------------------------------------------------------------

    def processing_ROIs_DfOverF(self, h5_object, string, dict, project_dir):
        module_name = string.split(".")[1]        
        iface_name  = string.split(".")[2]

        mod = h5_object.make_group("<Module>", module_name, \
                                   attrs ={"ancestry" : np.array(["Module"])})
        mod.set_custom_dataset("description", dict["processing.ROIs.description"])
        dff_iface = mod.make_group("DfOverF", attrs=dict[string + ".attrs"])
        string2   = "processing.ROIs.ImageSegmentation"
        seg_iface = mod.make_group("ImageSegmentation", 
                    attrs=dict[string2 + ".attrs"])
        processed_img_planes = []
        for k in dict.keys():
            if re.search("fov_", k):
                img_plane = k.split(".")[3]
                if img_plane in processed_img_planes:
                    continue
                processed_img_planes.append(img_plane)
                img_plane_group = dff_iface.make_group("<RoiResponseSeries>", img_plane, \
                                                    abort=False)
                path0 = string + "." + img_plane
                img_plane_group.set_dataset("data", dict[path0 + ".data"], attrs = \
                                         dict[path0 + ".data.attrs"])
                img_plane_group.set_dataset("timestamps", dict[path0 + ".timestamps"])
                img_plane_group.set_dataset("num_samples", len(dict[path0 + ".timestamps"]))
                img_plane_group.set_dataset("roi_names",  dict[path0 + ".roi_names"])
                img_plane_group.set_custom_dataset("trial_ids",  dict[path0 + ".trial_ids"])
                img_plane_group.make_group("segmentation_interface", seg_iface)
        return h5_object

# ------------------------------------------------------------------------------

    def processing_ROIs_ImageSegmentation(self, h5_object, string, dict, project_dir):
        module_name = string.split(".")[1]
        iface_name  = string.split(".")[2]

        mod = h5_object.make_group("<Module>", module_name, attrs ={"ancestry" : "<Module>"}, abort=False)
        seg_iface = mod.make_group("ImageSegmentation", attrs=dict[string + ".attrs"], abort=False)
        processed_img_planes = []
#       util.print_dict_keys(dict)       
        for k in dict.keys():
            if re.search("fov_", k):
                img_plane = k.split(".")[3]
                if img_plane in processed_img_planes:
                    continue
                processed_img_planes.append(img_plane)
                img_plane_group = seg_iface.make_group("<image_plane>", img_plane, \
                                                    attrs = {"ancestry" : np.array(["image_plane"])}, abort=False)
                path0 = string + "." + img_plane
                img_plane_group.set_dataset("description", dict[path0 + ".description"])
                img_plane_group.set_dataset("imaging_plane_name", dict[path0 + ".imaging_plane_name"])
                img_plane_group.set_dataset("roi_list",  dict[path0 + ".roi_ids"])
                try:
                    nwb_utils.add_reference_image(seg_iface, img_plane, "%s_0001"%img_plane, dict[path0 + ".ref_image_green"])
                except:
                    print("Warning: cannot store red reference image")
                try:
                    nwb_utils.add_reference_image(seg_iface, img_plane, "%s_0002"%img_plane, dict[path0 + ".ref_image_red"])
                except:
                    print("Warning: cannot store green reference image")
                
                roi_ids = dict[path0 + ".roi_ids"]
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

    def stimulus_presentation_timeseries(self, string, dict, project_dir):
        if debug:
            util.print_dict_keys(dict)       
        series_name = string.split(".")[-1]
        series_type = dict[string + ".attrs"]["series_type"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)

        series_path = "/stimulus/presentation"
        tsg = h5_object.make_group(series_type, series_name, path = series_path, \
                                  attrs = dict[string + ".attrs"])

        tsg.set_dataset("data", dict[string + ".data"], attrs=dict[string + ".data.attrs"])
        tsg.set_dataset("timestamps", dict[string + ".timestamps"])
        tsg.set_dataset("num_samples", len(dict[string + ".timestamps"]))
        for k in dict.keys():
            key = k.split(".")[-1]
            if key not in ["data", "timestamps", "num_samples"] and \
                   not re.search("attr", key):
                tsg.set_dataset(key, dict[string + "." + key], attrs=dict[string + ".data.attrs"])
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
       method_name = "stimulus_presentation_timeseries"
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

