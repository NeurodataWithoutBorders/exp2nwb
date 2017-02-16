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

verbose = 0    

class Exp2NWB(object):
    def initialize_h5_object(self, h5_file, dict):
        vargs = {"start_time" : 0.0, "description" : "default", "file_name" : h5_file, "identifier" : h5_file[0:-4], "mode" : "w"}
        for k in ["start_time", "description", "identifier"]:
            if k in dict.keys():
                vargs[k] = dict[k]

        if  os.path.exists(h5_file):
            os.remove(h5_file)   
        h5_object = nwb_file.open(**vargs)
        return h5_object

# ------------------------------------------------------------------------------

    def check_dict_keys(self, dict, target_datasets):
        for k in target_datasets:
            missing_key = 1
            for k1 in dict.keys():
                if k1 == k or k1.split(".")[-1] == k:
                    missing_key = 0
                    break
            if missing_key:
                if verbose:
                    print "   Warning: missing key: ", k

# ------------------------------------------------------------------------------

    def top_datasets(self, string, dict, project_dir):
        h5_object = nwb_file.open(**dict)
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def acquisition_timeseries(self, string, dict, project_dir):
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        if len(string.split(".")) == 3:
            series_name = string.split(".")[-1]
            series_type = dict[string + ".group_attrs"]["series_type"]
            series_path = "/acquisition/timeseries"
            tsg = h5_object.make_group(series_type, series_name, path = series_path, \
                              attrs = dict[string + ".group_attrs"])
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
            # Image time series
            path0 = "acquisition.timeseries."
            img_planes = []
            for k in dict.keys():
                if k.split(".")[-1] == "imaging_plane":
                    img_planes.append(k.split(".")[2])
            for img_plane in img_planes:
                twop = h5_object.make_group("<TwoPhotonSeries>", img_plane, \
                           path='/acquisition/timeseries', \
                           attrs=dict[path0 + "group_attrs"])
                for item in ["format", "external_file", "dimension", "scan_line_rate",\
                             "field_of_view", "imaging_plane", "timestamps"]:
                    if item == "external_file":
                        twop.set_dataset(item, dict[path0 + img_plane +"." + item],\
                                         attrs = dict[path0 + img_plane +".data_attrs"])
                    else:
                        twop.set_dataset(item, dict[path0 + img_plane +"." + item])
        else:
            sys.exit("Unsupported path string " + string)
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def acquisition_timeseries_extracellular_traces(self, string, dict, project_dir):
        datasets = ["ephys_raw_data"]

        self.check_dict_keys(dict, datasets)

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        et = h5_object.make_group("<TimeSeries>", "extracellular_traces", path = "/acquisition/timeseries", \
             attrs = {"description" : "File containing the raw voltage data for this session"})
        for k in  dict.keys():
            if re.search("ephys_raw_data", k):
                et.set_dataset("ephys_raw_data", dict[k])
            elif not re.search("attr", k):
                name = k.split(".")[-1]
                et.set_custom_dataset(name, dict[k])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def acquisition_images(self, string, dict, project_dir):

#       self.check_dict_keys(dict, datasets)

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        acquisition_group     = h5_object.make_group("acquisition", abort=False)
        acquisition_img_group = acquisition_group.make_group("images", abort=False)
        for k in  dict.keys():
            if re.search("description", k):
                acquisition_img_group.set_custom_dataset("description", dict[k])
            elif not k.split(".")[-1] == "attrs":
                name = k.split(".")[-1]
                h5_object.set_dataset("<image_X>", dict[k], name=name,\
                           dtype='uint8', attrs=dict[k + ".attrs"])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def analysis(self, string, dict, project_dir):
        analysis_datasets = ["description", "good_trials", "trial_start_times",\
                             "trial_type_mat", "trial_type_string"]
        self.check_dict_keys(dict, analysis_datasets)

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        ag = h5_object.make_group("analysis", abort=False)
        for k in dict.keys():
            k1 = k.split(".")[len(k.split("."))-1]
            if k in analysis_datasets or k1 in analysis_datasets:
                try:
                    ag.set_dataset(k1, dict[k])
                except:
                    ag.set_custom_dataset(k1, dict[k])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def epochs(self, string, dict, project_dir):
        epoch_datasets = ["start_time", "stop_time", "tags", "units_present", \
                          "description", "ROIs", "ROI_planes"]

        self.check_dict_keys(dict, epoch_datasets)

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        eg = h5_object.make_group("epochs", abort=False)
        processed_trials = []
        for k in dict.keys():
            trial = k.split(".")[-2]
            if trial not in processed_trials:
                processed_trials.append(trial)
                start = dict["epochs." + trial + "." + "start_time"]
                stop  = dict["epochs." + trial + "." + "stop_time"]
                epoch = nwb_utils.create_epoch(h5_object, trial, start, stop)
                for k1 in epoch_datasets:
                    k2 = "epochs." + trial + "." + k1
                    if k2 not in dict.keys():
                        continue
                    if k1 in ["tags", "description"]:
                        epoch.set_dataset(k1, dict[k2])
                    elif k1 in ["units_present", "ROIs", "ROI_planes"]:
                        epoch.set_custom_dataset(k1, dict[k2])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def general_top_datasets(self, string, dict, project_dir):
        top_datasets = ["data_collection", "experiment_description", "experimenter", \
                        "institution", "lab", "notes", "pharmacology", "protocol", \
                        "related_publications", "session_id", "slices", "source_script", \
                        "stimulus", "sugrery", "virus"]
        self.check_dict_keys(dict, top_datasets)

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        gg = h5_object.make_group("general", abort=False)
        for k in dict.keys():
            k1 = k.split(".")[len(k.split("."))-1]
            if k in top_datasets or k1 in top_datasets:
                try:
                    gg.set_dataset(k1, dict[k])
                except:
                    gg.set_custom_dataset(k1, dict[k])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def general_devices(self, string, dict, project_dir):
        device_datasets = ["ephys_acquisition", "photostim_source"]

        self.check_dict_keys(dict, device_datasets)

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        for k in dict.keys():
            h5_object.set_dataset("<device_X>", dict[k], name=k)
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def general_extracellular_ephys(self, string, dict, project_dir):
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
                   print "\nCreating dataset ", k.split(".")[-1], " value=", dict[k]
               eeg.set_dataset(k.split(".")[-1], dict[k])
                
            elif re.search("shank", k.split(".")[-2]):
                if verbose:
                    print "\nCreating group ", k.split(".")[-2]
                sg = eeg.make_group("<electrode_group_X>", k.split(".")[-2], abort=False)
                if k.split(".")[-1] in shank_datasets:
                    if verbose:
                        print "\nCreating dataset ", k.split(".")[-1], " value=", dict[k]
                    sg.set_dataset(k.split(".")[-1], dict[k])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def general_subject(self, string, dict, project_dir):
        subject_datasets = ["age", "description", "genotype", "sex", "species", \
                            "subject_id"]

        self.check_dict_keys(dict, subject_datasets)

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        gg = h5_object.make_group("general", abort=False)
        sg = gg.make_group("subject", abort=False)
        for k in dict.keys():
            if k in subject_datasets or k.split(".")[-1] in subject_datasets:
                try:
                    sg.set_dataset(k.split(".")[-1], dict[k])
                except:
                    sg.set_custom_dataset(k.split(".")[-1], dict[k])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def general_optogenetics(self, string, dict, project_dir):
        optogenetics_datasets = ["device", "excitation_lambda", "location", "stimulation_method"]

        self.check_dict_keys(dict, optogenetics_datasets)

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        gg = h5_object.make_group("general", abort=False)
        og = gg.make_group("optogenetics", abort=False)
        for k in dict.keys():
            if k.split(".")[-1] in optogenetics_datasets and k.split(".")[-3] == "optogenetics":
                osg = og.make_group("<site_X>", name=k.split(".")[-2], abort=False)
                try:
                    osg.set_dataset(k.split(".")[-1], dict[k])
                except:
                    osg.set_custom_dataset(k.split(".")[-1], dict[k])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def general_optophysiology(self, string, dict, project_dir):
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        gg = h5_object.make_group("general", abort=False)
        op = gg.make_group("optophysiology", abort=False)
        for k in dict.keys():
            if not re.search("fov_", k) and re.search("description", k):
                op.set_custom_dataset("description", dict[k])
            elif re.search("fov_", k) and not re.search("channel_", k):
                name_group = k.split(".")[2]
                img_plane_group = op.make_group("<imaging_plane_X>", name_group, abort=False)
                name_data = k.split(".")[3]
                if not re.search("attrs", k):
                    if name_data == "manifold":
                        img_plane_group.set_dataset(name_data, dict[k], attrs = dict[k + ".attrs"])
                    else:
                        img_plane_group.set_dataset(name_data, dict[k])
            elif re.search("fov_", k) and re.search("channel_", k):
                name_group = k.split(".")[2]
                img_plane_group = op.make_group("<imaging_plane_X>", name_group, abort=False)
                name_channel = k.split(".")[3]
                channel_group = img_plane_group.make_group("<channel_X>", name_channel, abort=False)
                name_data = k.split(".")[4]
                channel_group.set_dataset(name_data, dict[k])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def processing_extracellular_units(self, string, dict, project_dir):
        if re.search("top_datasets", string):
            self.processing_extracellular_units_top_datasets(string, dict, project_dir)
        elif re.search("EventWaveform", string):
            self.processing_extracellular_units_event_waveform(string, dict, project_dir)
        elif re.search("UnitTimes", string):
            self.processing_extracellular_units_unit_times(string, dict, project_dir)

# ------------------------------------------------------------------------------

    def processing_extracellular_units_unit_times(self, string, dict, project_dir):
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        module_name = "extracellular_units"
        mod = h5_object.make_group("<Module>", module_name)

        spk_times_iface = mod.make_group("UnitTimes")
        spk_times_iface.set_attr("source", dict[string + ".group_attrs"]["source"])

        processed_units = []

        for k in dict.keys():
            unit = k.split(".")[-2]
            if unit not in processed_units \
                and not re.search("attr", k)  \
                and not re.search("cell_types", k):             
                processed_units.append(unit)
                path = string + "." + unit
                spk = spk_times_iface.make_group("<unit_N>", unit)
                spk.set_custom_dataset("times",            dict[path + ".times"])
                spk.set_custom_dataset("source",           dict[path + ".source"])
                spk.set_custom_dataset("trial_ids",        dict[path + ".trial_ids"])
                spk.set_custom_dataset("unit_description", dict[path + ".unit_description"])

        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def processing_extracellular_units_event_waveform(self, string, dict, project_dir):
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        module_name = "extracellular_units"
        mod = h5_object.make_group("<Module>", module_name)

        spk_waves_iface = mod.make_group("EventWaveform")
        spk_waves_iface.set_attr("source", dict[string + ".group_attrs"]["source"])

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

        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def processing_extracellular_units_top_datasets(self, string, dict, project_dir):
        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        module_name = "extracellular_units"
        mod = h5_object.make_group("<Module>", module_name)

        path = ".".join(string.split(".")[0:-1])
        mod.set_custom_dataset('description',           dict[path + ".description"])
        mod.set_custom_dataset('spike_sorting',         dict[path + ".spike_sorting"])
        mod.set_custom_dataset('identification_method', dict[path + ".identification_method"])

        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def processing_timeseries(self, string, dict, project_dir):
        module_name = string.split(".")[1]
        series_name = string.split(".")[-1]
        if verbose:
            print "dict.keys()=", dict.keys()
        series_type = dict[string + ".group_attrs"]["series_type"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        try:
            mod = h5_object.make_group("<Module>", module_name, attrs = \
                                        dict[string + ".module_attrs"])
        except:
            mod = h5_object.make_group("<Module>", module_name)

        iface = mod.make_group(string.split(".")[2], attrs=dict[string + ".group_attrs"])
        tsg   = iface.make_group(series_type, series_name, attrs = dict[string + ".group_attrs"])

        tsg.set_dataset("data", dict[string + ".data"], attrs=dict[string + ".data_attrs"])
        tsg.set_dataset("timestamps", dict[string + ".timestamps"])
        tsg.set_dataset("num_samples", len(dict[string + ".timestamps"]))

        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def processing_ROIs(self, string, dict, project_dir):
        self.processing_ROIs_DfOverF(string + ".DfOverF", dict, project_dir)
        self.processing_ROIs_ImageSegmentation(string + ".ImageSegmentation", dict, project_dir)

# ------------------------------------------------------------------------------

    def processing_ROIs_DfOverF(self, string, dict, project_dir):
        module_name = string.split(".")[1]        
        iface_name  = string.split(".")[2]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        mod = h5_object.make_group("<Module>", module_name)
        mod.set_custom_dataset("description", dict["processing.ROIs.description"])
        dff_iface = mod.make_group("DfOverF", attrs=dict[string + ".group_attrs"])
        seg_iface = mod.make_group("ImageSegmentation", \
                    attrs=dict[string + ".group2_attrs"])
        processed_img_planes = []
        for k in dict.keys():
            if re.search("fov_", k):
                img_plane = k.split(".")[3]
                if img_plane in processed_img_planes:
                    continue
                processed_img_planes.append(img_plane)
                img_plane_gr = dff_iface.make_group("<RoiResponseSeries>", img_plane, abort=False)
                path0 = string + "." + img_plane
                img_plane_gr.set_dataset("data", dict[path0 + ".data"], attrs = \
                                         dict[path0 + ".data_attrs"])
                img_plane_gr.set_dataset("timestamps", dict[path0 + ".timestamps"])
                img_plane_gr.set_dataset("num_samples", len(dict[path0 + ".timestamps"]))
                img_plane_gr.set_dataset("roi_names",  dict[path0 + ".roi_names"])
                img_plane_gr.set_custom_dataset("trial_ids",  dict[path0 + ".trial_ids"])
                img_plane_gr.make_group("segmentation_interface", seg_iface)
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def processing_ROIs_ImageSegmentation(self, string, dict, project_dir):
        module_name = string.split(".")[1]
        iface_name  = string.split(".")[2]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)
        mod = h5_object.make_group("<Module>", module_name)

        seg_iface = mod.make_group("ImageSegmentation", attrs=dict[string + ".group_attrs"])
        processed_img_planes = []
        processed_roi_ids    = []
        for k in dict.keys():
            if re.search("fov_", k):
                img_plane = k.split(".")[3]
                if not img_plane in processed_img_planes:
                    continue
                processed_img_planes.append(img_plane)
                img_plane_gr = seg_iface.make_group("<image_plane>", img_plane, abort=False)
                path0 = string + "." + img_plane
                img_plane_gr.set_dataset("description", dict[path0 + ".description"])
                img_plane_gr.set_dataset("imaging_plane_name", dict[path0 + ".imaging_plane_name"])
                img_plane_gr.set_dataset("roi_list",  dict[path0 + ".roi_list"])
                roi_ids = dict[path0 + ".roi_list"]
                if len(k.split(".")) == 6 and (img_plane + "_" + roi_id) not in processed_roi_ids:
                    roi_id = k.split(".")[4]
                    processed_roi_ids.append(img_plane + "_" + roi_id)
                    nwb_utils.add_roi_mask_pixels(seg_iface, img_plane, "%d"%x, "ROI %d"%x, np.array(pixmap).astype('uint16'), \
            weight, master1_shape[1], master1_shape[0])
        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

    def stimulus_presentation_timeseries(self, string, dict, project_dir):
        series_name = ""
        for k in dict.keys():
            series_name = k.split(".")[-2]

        my_list = string.split(".")[0:-1] + [ series_name ]
        string = ".".join(my_list)
        series_type = dict[string + ".group_attrs"]["series_type"]

        h5_object = self.initialize_h5_object(os.path.join(project_dir, string+ ".h5"), dict)

        series_path = "/stimulus/presentation"
        tsg = h5_object.make_group(series_type, series_name, path = series_path, \
                                  attrs = dict[string + ".group_attrs"])

        tsg.set_dataset("data", dict[string + ".data"], attrs=dict[string + ".data_attrs"])
        tsg.set_dataset("timestamps", dict[string + ".timestamps"])
        tsg.set_dataset("num_samples", len(dict[string + ".timestamps"]))

        print "\nCreating partial NWB file: ", os.path.join(project_dir, string+ ".h5")
        h5_object.close2()

# ------------------------------------------------------------------------------

def make_h5(string, dict, project_dir):
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
    class_instance = Exp2NWB()
#   print "\nmethod_name=", method_name

    method = None
    try:
        method = getattr(class_instance, method_name)
    except:
        sys.exit("\nMethod " + method_name + "' has not been implemented")

    method(string, dict, project_dir)

