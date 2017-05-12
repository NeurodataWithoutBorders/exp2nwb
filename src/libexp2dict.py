# Copyright (c) 2017, HHMI-Janelia Research Campus All rights reserved.
#
# libexp2dict.py
#
# Extracting data dictionaries from experimental data structures
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
import util  

class Exp2Dict(object):

    def top_datasets(self, data, metadata, options):
        dict = {}
        session_id = options.project_dir
        if options.data_type == "ephys":
            dict["start_time"] = util.find_exp_time(metadata, options)
            dict["description"]= "Extracellular ephys recording of mouse doing discrimination " + \
                              "task (lick left/right), with optogenetic stimulation " +\
                              "plus pole and auditory stimulus"
        elif options.data_type == "ophys": 
            dict["start_time"] = util.find_exp_time(data, options)
            dict["description"]= " "
            print("    Warning: missing experiment description (to be added)")
        dict["file_name"]      = os.path.join(options.project_dir, "top_datasets.h5")
        dict["identifier"]     = nwb_utils.create_identifier(options.project_dir)
        dict["mode"] = "w"

        return dict
 
# ------------------------------------------------------------------------------

# For ophys data only

    def acquisition_images(self, data, metadata, options):
        dict = {}
        plane_map = util.create_plane_map(data, options)
        num_subareas = len(util.get_key_list(data['timeSeriesArrayHash'])) - 1
        path0 = "acquisition.images."
        description = "Mean images for each plane and color channel for the experimental session.  fov_VVPPP_color contains the image, in 8 bit integer format, for subvolume VV and plane PPP.  By convention, plane 001 (of four) is the flyback frame and therefore excluded.  See Peron et al 2015 Neuron for further details (filters, etc.)"
        dict[path0 + "description"] = description
        for subarea in range(num_subareas):
            plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
            if data[plane_path].keys()[0] == 'masterImage':
                num_planes = 1
            else:
                num_planes = len(data[plane_path].keys())
            for plane in range(num_planes):
                area_grp = data["timeSeriesArrayHash/descrHash"]["%d"%(subarea + 2)]
                if num_planes == 1:
                    plane_grp = area_grp["value/1"]
                else:
                    plane_grp = area_grp["value/1"]["%d"%(plane+1)]
                master = plane_grp["masterImage"]["masterImage"].value
                master_shape = np.array(master).shape
                oname = "area%d_plane%d" % (subarea+1, plane+1)
                image_plane = plane_map[oname]
                master_green = np.zeros([master_shape[0],master_shape[1]])
                master_red   = np.zeros([master_shape[0],master_shape[1]])
                if len(master_shape) == 3 or not options.handle_errors:
                    master_green = master[:,:,0]
                    master_red   = master[:,:,1]
                else:
                    master_green = master
                    master_red   = master
                # Store green and red referebce images
                desc_green = "Master image (green channel), in " + str(master_shape[0]) + \
                             "x" + str(master_shape[1]) + ", 8bit"
                dict[path0 + image_plane + "_green"] = master_green.astype('uint8')
                dict[path0 + image_plane + "_green.attrs"] = {"description":desc_green, "format": "raw"}
                desc_red = "Master image (red channel), in " + str(master_shape[0]) + \
                           "x" + str(master_shape[1]) + ", 8bit"
                dict[path0 + image_plane + "_red"] = master_red.astype('uint8')
                dict[path0 + image_plane + "_red.attrs"] = {"description":desc_red, "format": "raw"}
        return dict

# ------------------------------------------------------------------------------

    # For ophys data only

    def acquisition_timeseries(self, data, metadata, options):
        dict = {}
        path0 = "acquisition.timeseries."
        plane_map = util.create_plane_map(data, options)
        description = "Information about raw image data, which must be obtained from the original data files. 'external_file' provides a list of raw data files; 'external_file.attribute' starting_frame indicates, with 0 as the first value, the frame for this plane at which a given file starts.  Thus, if starting_frame is 0,52,..., then the first file contains frames 1-52, or 0-51. See processing/ROIs/DfOverF for fluorescence change for individual ROIs"
        dict[path0 + "description"] = description
        group_attrs = {"source": "Device 'two-photon microscope'",
                       "description": "2P image stack, one of many sequential stacks for this field of view"}
        dict[path0 + "attrs"] = group_attrs
        num_subareas = len(util.get_key_list(data['timeSeriesArrayHash'])) - 1
        for subarea in range(num_subareas):
            grp = data["timeSeriesArrayHash"]["value"]["%d"%(subarea+2)]
            t = 0.001 * grp["time"]["time"].value    # fetch time array
            plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
            if data[plane_path].keys()[0] == 'masterImage':
                num_planes = 1
            else:
                num_planes = len(data[plane_path].keys())
            for plane in range(num_planes):
                oname = "area%d_plane%d" % (subarea+1, plane+1)
                img_plane = plane_map[oname]
                try:
                    pgrp = grp["imagingPlane/%d"%(plane+1)] # move into imaging plane group to extract 2photon image stacks
                except:
                    pgrp = grp["imagingPlane"] # only one imaging plane is available (instead of 3)
                frame_idx = pgrp["sourceFileFrameIdx"]["sourceFileFrameIdx"].value
                lst = util.parse_h5_obj(pgrp["sourceFileList"])[0]
                cnt = 0
                srcfile = {}
                for k in range(len(lst)):
                    srcfile[str(k+1)] = lst[k]
                    cnt += 1
                assert len(t) == len(frame_idx[0])
                external_file, starting_frame, timestamps, nname = \
                    util.extract_frame_data(frame_idx, t, subarea, plane, plane_map, srcfile, cnt, options)
                dict[path0 + img_plane + ".format"]      = "external"
                # need to convert name to utf8, otherwise error generated:
                dict[path0 + img_plane + ".external_file"]  = external_file
                dict[path0 + img_plane + ".external_file.attrs"] = {"starting_frame": starting_frame}
                dict[path0 + img_plane + ".dimension"]      = [512, 512]
                dict[path0 + img_plane + ".scan_line_rate"] = 16000.0
                dict[path0 + img_plane + ".field_of_view"]  = [ 600e-6, 600e-6 ]
                dict[path0 + img_plane + ".imaging_plane"]  = img_plane
                dict[path0 + img_plane + ".timestamps"]     = timestamps
        return dict

# ------------------------------------------------------------------------------

    def acquisition_timeseries_lick_trace(self, data, metadata, options):
        dict = {}
        hash_group_pointer = data["timeSeriesArrayHash"]
        if options.data_type == "ephys":
            keyName     = "EphusVars"
            data        = hash_group_pointer["value/valueMatrix/valueMatrix"][:,0]
            timestamps  = hash_group_pointer["value/time/time"].value
            description = util.parse_h5_obj(hash_group_pointer["value/idStrDetailed/idStrDetailed"])[0][0]
        elif options.data_type == "ophys":
            keyName     = "lickVars"
            data        = hash_group_pointer["value/2/valueMatrix/valueMatrix"]
            timestamps  = hash_group_pointer["value/2/time/time"].value
            description = util.parse_h5_obj(hash_group_pointer["value/2/idStrDetailed/idStrDetailed"])[0][0]
        comment1 = keyName
        comment2 = util.get_description_by_key(hash_group_pointer, keyName)
        comments = comment1 + ": " + comment2
        data_attrs={"conversion":1.0, "unit":"unknown", "resolution":float('nan')}
        group_attrs={"description" : description, "comments" : comments, \
                    "source": "Times as reported in Nuo's data file", \
                    "series_type" : "<TimeSeries>"}

        series_path = "acquisition.timeseries.lick_trace"
        dict[series_path + ".attrs"]       = group_attrs
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".num_samples"] = len(timestamps)
        return dict

# ------------------------------------------------------------------------------

    def acquisition_timeseries_extracellular_traces(self, data, metadata, options):
        dict = {}
        fname = "voltage_filename" + os.path.basename(options.data_path)[13:-3] + ".mat"
        dict["acquisition.timeseries.extracellular_traces.ephys_raw_data"] = fname
        dict["acquisition.timeseries.extracellular_traces.attrs"] = {"source" : " ",\
             "description" : "File containing the raw voltage data for this session" }
        return dict

# ------------------------------------------------------------------------------

    def analysis(self, data, metadata, options):
        dict = {}
        if options.verbose:
            print("Collecting analysis information")

        trial_start_times = data["trialStartTimes/trialStartTimes"].value
        trial_types_all = []
        trial_type_strings = util.parse_h5_obj(data['trialTypeStr'])[0]
        # collect all trials (strings)
        for i in range(len(trial_type_strings)):
                trial_types_all.append(str(trial_type_strings[i]))
        trial_type_mat = data['trialTypeMat/trialTypeMat'].value

        if options.data_type == "ephys":
            good_trials = data['trialPropertiesHash/value/4/4'].value
        else:
            good_trials = [0] * len(trial_start_times)
            valid_whisker = util.get_valid_trials(data, "whisker", options)
            valid_Ca      = util.get_valid_trials(data, "Ca", options)
            for i in range(len(good_trials)):
                j = i + 1
                if j in valid_whisker and j in valid_Ca:
                    good_trials[i] = 1
        dict["analysis.trial_start_times"] = trial_start_times
        dict["analysis.trial_type_string"] = trial_types_all
        dict["analysis.trial_type_mat"]    = trial_type_mat
        dict["analysis.good_trials"]       = good_trials
        dict["analysis.description"]="Trial_type_string has six values, with three response possibilities and two correct responses.  The format is {response type}{correct response}, with response type taking the value Hit, Err, or NoLick, indicating that the animal got the trial right, wrong, or did nothing, respectively.  Corret response is either L or R, indicating that the animal should have licked left or right, respectively."

        return dict

# ------------------------------------------------------------------------------

    def epochs(self, data, metadata, options):
        dict = {}

        if options.verbose:
            print("    Getting trial types ...")
        util.get_trial_types(data, options)

        if options.data_type == "ephys":
            keys = data['eventSeriesHash/value'].keys()
            if options.verbose:
                print("    Getting trial units ...")
            keys = data['eventSeriesHash/value'].keys()
            unit_num = len(keys)
            util.get_trial_units(data, unit_num, options)

        if options.verbose:
            print("    Creating trials ...")
        dict = util.create_trials(dict, data, options)

        return dict

# ------------------------------------------------------------------------------

    def general_top_datasets(self, data, metadata, options):
        dict = {}
        dict_general = self.general(data, metadata, options)
        for k in dict_general.keys():
            if len(k.split(".")) == 2:
                dict[k] = dict_general[k]
        return dict

# ------------------------------------------------------------------------------

    def general_devices(self, data, metadata, options):
        dict = {}
        dict_general = self.general(data, metadata, options)
        for k in dict_general.keys():
            if len(k.split(".")) == 3 and k.split(".")[1] == "devices":
                dict[k] = dict_general[k]
        return dict

# ------------------------------------------------------------------------------

    def general_subject(self, data, metadata, options):
        dict = {}
        dict_general = self.general(data, metadata, options)
        for k in dict_general.keys():
            if len(k.split(".")) == 3 and k.split(".")[1] == "subject":
                dict[k] = dict_general[k]
        return dict

# ------------------------------------------------------------------------------

    def general_optogenetics(self, data, metadata, options):
        dict = {}
        dict_general = self.general(data, metadata, options)
        for k in dict_general.keys():
            if len(k.split(".")) == 4 and k.split(".")[1] == "optogenetics":
                dict[k] = dict_general[k]
        dict["general.optogenetics.source"] = " " 
        return dict

# ------------------------------------------------------------------------------

    def general(self, data, metadata, options):
        if options.verbose:
            print("Extracting metadata to a dictionary")
        dict = {}

        genotype = ""
        subject  = ""
        weight   = ""
        animalStrain = ""
        animalSource = ""
        animalGeneModification = ""
        value = ""
        DOB   = ""
        dict["general.institution"] = "Janelia Research Campus"
        dict["general.lab"]         =  "Svoboda lab"

        # Add citation info
        if options.data_type == "ephys":
            key_list = util.get_child_group_names(metadata)
            group_pointer = metadata
        else:
            key_list = util.get_key_list(metadata["metaDataHash"])
            group_pointer = metadata['/metaDataHash']

        for k in ["citation", "virus"]:
            if not k in key_list:
                key_list.append(k)

        if options.verbose:
            print("metadata key_list=", key_list)
        for key in key_list:
            if key in ["behavior", "virus", "fiber", "extracellular", "intracellular"]:
                # Note: these are two-level groups
                value = ""
                try:
                    key1_list = util.get_value_pointer_by_path_items(group_pointer, [key]).keys()
                    if options.verbose:
                        print("   key=" + key+ " key1_list="+ str(key1_list))
                    if re.search("tracellular", key) and not "impedance" in key1_list:
                        dict["general." + key + ".impedance"] = ["not recorded"]
                    for key1 in key1_list:
                        if options.verbose:
                            print("      key1="+ key1)
                        if key1 in ["siteLocations"]:
                            continue
                        elif key1 in ["ADunit", "penetrationN", "groundCoordinates", \
                                      "extracellularDataType", "recordingMarker", \
                                      "recordingType"]:
                            value1 = util.get_value_pointer_by_path_items(group_pointer, [key, key1, key1])
                            if key1 == "ADunit":
                                dict["general." + key + ".ADunit"] = str(value1[0])
                            elif key1 == "penetrationN":
                                dict["general." + key + ".penetration_num"] = str(value1[0])
                            elif key1 == "groundCoordinates":
                                dict["general." + key + ".ground_coordinates"]= value1
                            elif key1 == "extracellularDataType":
                                dict["general." + key + ".data_types"]= value1
                            elif key1 == "recordingMarker":
                                dict["general." + key + ".recording_marker"]= str(value1[0])
                            elif key1 == "recordingType":
                                dict["general." + key + ".recording_type"]= str(value1[0])
                            elif key1 == "impedance":
                                dict["general." + key + ".impedance"]= [str(value1[0])]
                        value1 = util.get_value_pointer_by_path_items(group_pointer, [key, key1, key1])
                        if options.verbose:
                            print("      key="+ key+ " key1="+ key1+ " value1="+ str(value1))
                        value2 = [str(v) for v in value1]
                        add_value = "       " + key1 + ": " + ",".join(value2) + "\n "
                        value += add_value
                        if options.verbose:
                            print("         key=", key, " key1=", key1, " add_value=", add_value, " now value=\n", value)
                except:
                    if key == "virus" and options.data_type == "ophys":
                        value = "virusID: AAV2/1 syn-GCaMP6s; for source, see the related_publications"
                        value += "\ninfectionCoordinates: 3.6 mm lateral (left), 1.5 mm posterior, 900 um AP extent, 1200 um ML extent"
                        value += "\ninfectionLocation: vibrissal S1"
                        value += "\ninjectionVolume: 20 nL per site, with 12 sites per anima"
                        value += "\nvirusLotNumber: Penn Vector Core #Av-1-PV2824"
                        value += "\nvirusSource: Penn vector core"
                        value += "\nvirusTiter: approx. 10^13"
                key1_list = []
    
            elif key == "photostim":
                # Exctract data from H5 file and count the number of sites
                key1_list = util.get_value_pointer_by_path_items(group_pointer, [key]).keys()
                value1 = {}
                if options.verbose:
                    print("key= photostim, key1_list=" + str(key1_list))
                num_sites = 1      
                photostim_loc1 = '' 
                for key1 in key1_list:
                    try:
                        value1[key1] = util.get_value_by_path_items(group_pointer, [key, key1])
                        if key1 == "photostimLocation":
                            if options.verbose:
                                print("photostimLocation=" + str(value1[key1]))
                            photostim_loc1 = value1[key1][0]
                            if len(value1[key1]) == 3:
                                num_sites = 2                   
                    except:
                        value1[key1] = []
                        for i in range(num_sites):
                            value1[key1].append("NA")
    
                # Create the site groups and populate the data fields in the NWB file
                for s in range(num_sites):
                    if num_sites == 1:
                        if photostim_loc1 == 'PONS':
                            group_name = "site_" + str(1)
                        else:
                            group_name = "site_" + str(2)
                for s in range(num_sites):
                    if options.data_type == "ephys":
                        dict["general.optogenetics.site_" + str(s+1) + ".device"] = "Ciel 250 mW 473nm Laser from Laser Quantum"
                        dict["general.optogenetics.site_" + str(s+1) + ".description"] = " "
                    for key1 in key1_list:
                        if key1 == "stimulationMethod":
                            dict["general.optogenetics.site_" + str(s+1) + ".stimulation_method"]=str(value1[key1][s])
                        elif key1 == "photostimCoordinates":
                            coord = str(value1[key1][s]) + " in mm"
                            if "photostimLocation" in key1_list:
                                coord += "\natlas location: " + str(value1["photostimLocation"][s]) 
                            dict["general.optogenetics.site_" + str(s+1) + ".location" ] = coord
                        elif key1 == "photostimWavelength":
                            dict["general.optogenetics.site_" + str(s+1) + ".excitation_lambda"] = str(value1[key1][0])
                        
            else:
                if "metaDataHash" in util.get_child_group_names(metadata) \
                    and not key in ["citation"]:
                    value = util.get_value_by_key(group_pointer,key)
                    if util.item_type(value) == "dataset":
                        value = np.array(value).tolist()
                    elif util.item_type(value) == "group":
                        value = value.name
                elif not "metaDataHash" in util.get_child_group_names(metadata):                          
                    value_list = np.array(util.get_value_pointer_by_path_items(group_pointer, [key, key])).tolist()
                    if len(value_list) == 1:
                        value = value_list[0]
                    elif len(value_list) > 0:
                        value_str_list = [str(v) for v in value_list]
                        value = ",".join(value_str_list)
                    else:
                        value = ""
            if options.verbose:
               print("            key=", key, " final value=", value)
    
            if key in ["animalGeneCopy", "animalGeneticBackground"] or \
               re.search("animalGeneModification", key) or \
               re.search("animalSource", key) or \
               re.search("animalStrain", key):
                genotype += key + ": " + str(value) + "\n"               
       
            elif re.search("animalStrain", key):
                animalStrain += key + ": " + value + "\n"

            elif re.search("animalSource", key):
                animalSource += key + ": " + value + "\n"

            elif key == "animalID":
                dict["general.subject.subject_id"] = value
 
            elif key == "dateOfBirth":
                subject  = key + ": " + str(value) + "\n"
                DOB = str(value)

            elif key == "dateOfExperiment":
                age = ""
                if options.verbose:
                    print "\noptions.data_type=", options.data_type     
                if options.data_type == "ephys":
                    DOE = str(value)
                    age = ""
                    if len(DOB) ==8 and len(DOE)== 8:                   
                        age = util.compute_age(DOB, DOE)
                    else:
                        age = "3 to 5 months"
                elif options.data_type == "ophys":
                    age = "6 to 8 weeks"
                if len(age) > 0:
                    dict["general.subject.age"]= age
    
            elif re.search("citation", key):
                # Add citation info
                if options.data_type == "ephys":
                    value = "doi: 10.1038/nature14178"
                else:
                    value = "doi: 10.1016/j.neuron.2015.03.027"
                dict["general.related_publications"]=value
    
            elif re.search("experimentType", key):
                dict["general.notes"]=value
    
            elif re.search("experimenters", key):
                dict["general.experimenter"]=value
    
            elif key in ["sex", "species", "age", "cell"]:
                dict["general.subject." + key]=value
    
            elif re.search("weight", key) and len(str(weight)) == 0:
                if len(str(value)) > 0:
                    weight = str(value)
                else:
                    weight = "not recorded"
                dict["general.subject.weight"]=weight
    
            elif re.search("referenceAtlas", key):
                dict["general.reference_atlas"]=value
    
            elif re.search("whiskerConfig", key):
                dict["general.whisker_configuration"] = value
    
            elif key in ["virus", "fiber"]:
                dict["general." + key]=value
    
            elif key == "behavior":
                task_kw = map(str,util.parse_h5_obj(metadata["behavior/task_keyword"])[0])
                dict["general.task_keywords"]=task_kw
    
        dict["general.subject.genotype"]    = genotype + "\n"
        dict["general.subject.description"] =  animalStrain + "\n" + animalSource + "  " + subject
        source_script="https://github.com/NeurodataWithoutBorders/exp2nwb"
        dict["general.source_script"] = source_script
        if options.data_type == "ephys":
            dict["general.surgery"]                = util.metadata_from_file("surgery.txt", options)
            dict["general.data_collection"]        = util.metadata_from_file("data_collection.txt", options)
            dict["general.experiment_description"] = util.metadata_from_file("experiment_description.txt", options)
        elif options.data_type == "ophys":
            dict["general.data_collection"]       ="doi: 10.1016/j.neuron.2015.03.027"
            dict["general.experiment_description"]="doi: 10.1016/j.neuron.2015.03.027"
            dict["general.surgery"]               ="doi: 10.1016/j.neuron.2015.03.027"
            dict["general.whisker_configuration"] = "see Table S1 in doi: 10.1016/j.neuron.2015.03.027"
            dict["general.subject.weight"]        = "not recorded"
        else:
            print("Missing  data_collection.txt, experiment_description.txt and surgery.txt")
    
        return dict

# ------------------------------------------------------------------------------

    def general_devices(self, data, metadata, options):
        dict = {}
        dict["ephys_acquisition"] = "32-electrode NeuroNexus silicon probes recorded on a PCI6133 National Instrimunts board. See 'general/experiment_description' for more information"
        dict["photostim_source"] = "Stimulating laser at 473 nm"
        return dict
    
# ------------------------------------------------------------------------------
    
    def general_extracellular_ephys(self, data, metadata, options):
        dict = {}
        M = util.read_probe_locations_matrix(metadata, options)
    
        # probe = M.tolist()
        probe = []
        sites = util.parse_h5_obj(util.check_entry(metadata, "extracellular/siteLocations"))
    #   assert len(sites) == 32, "Expected 32 electrode locations, found %d"%len(sites)
        for i in range(len(sites)):
            probe.append(sites[i])
            probe[-1] = probe[-1] * 1.0e-3
        probe = np.asarray(probe)
    
        num_shanks, shank_size, shank_coords = util.detect_shanks(M)
        shank = []
        for i in range(1, (num_shanks+1)):
            for j in range(shank_size):
                shank.append("shank" + str(i))
        dict["general.extracellular_ephys.impedance"] = ["not recorded"]
        dict["general.extracellular_ephys.electrode_map"] = probe
        dict["general.extracellular_ephys.electrode_group"] = shank
        dict["general.extracellular_ephys.filtering"]  = "Bandpass filtered 300-6K Hz"
    
        ephys_device_txt = "32-electrode NeuroNexus silicon probes recorded on a PCI6133 National Instrimunts board. See 'general/experiment_description' for more information"
        dict["general.devices.ephys_acquisition"] = ephys_device_txt
        dict["general.devices.photostim_source"]  = "Stimulating laser at 473 nm"
    
        # Creating the shank groups
        probe_type = util.get_value_pointer_by_path_items(metadata, \
                         ["extracellular", "probeType", "probeType"])
        rloc = util.get_value_pointer_by_path_items(metadata, \
                   ["extracellular", "recordingLocation", "recordingLocation"])
        description = util.get_description(metadata, options)
        device      = util.get_device(     metadata, options)
        if options.verbose:
            print("description=" + str(description))
            print("device="      + str(device))
            print("rloc="        + str(rloc))
        for i in range(num_shanks):
            loc = str(rloc[0])
            P = str(shank_coords[i][0])
            Lat = str(shank_coords[i][1])
            location = "loc: " + str(loc) + ", P: " + str(P) + ", Lat: " + str(Lat) + ", recordingLocation=" + str(rloc)
            dict["general.extracellular_ephys.shank_" + str(i) + ".location"]    = location
            dict["general.extracellular_ephys.shank_" + str(i) + ".description"] = description
            dict["general.extracellular_ephys.shank_" + str(i) + ".device"]      = device
        return dict
    
# ------------------------------------------------------------------------------
    
    def general_optophysiology(self, data, metadata, options):
#       dict = self.general(data, metadata, options)
        dict = {}
        plane_map = util.create_plane_map(data, options)
        num_subareas = len(util.get_key_list(data['timeSeriesArrayHash'])) - 1
        master_shape = util.extract_dict_master_shape(data, plane_map, options)
    
        dict["general.optophysiology.description"] = \
             "Imaging data were acquired using a custom two-photon microscope; " + \
             "see Peron et al., 2015 (DOI: 10.1016/j.neuron.2015.03.027) for details."
        for subarea in range(num_subareas):
            plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
            if data[plane_path].keys()[0] == 'masterImage':
                num_planes = 1
            else:
                num_planes = len(data[plane_path].keys())
            for plane in range(num_planes):
                oname = "area%d_plane%d" % (subarea+1, plane+1)
                manifold = util.define_manifold(master_shape[oname], plane_map[oname])
                plane_name = plane_map[oname]
                path0 = "general.optophysiology." + plane_name
                dict[path0 + ".description"] = "Imaging data were acquired using a custom two-photon microscope; see Peron et al., 2015 (DOI: 10.1016/j.neuron.2015.03.027) for details."
                dict[path0 + ".location"]    = "Vibrissal somatosensory cortex (approx. 1.6 mm posterior, 3.6 mm lateral Bregma), centered on spared whisker column"
                dict[path0 + ".indicator"] = "GCaMP6s"
                dict[path0 + ".excitation_lambda"] = "1000 nm"
                dict[path0 + ".imaging_rate"] = "16KHz (scan line rate); 7Hz per plane, 28Hz frame rate (4 planes, 1 of which is discarded flyback)"
                dict[path0 + ".device"] = "DOM3 (Customized Svoboda lab MIMMS microscope; see https://www.janelia.org/open-science/mimms)"
                dict[path0 + ".manifold"] = manifold
                dict[path0 + ".manifold.attrs"]={"note":"Manifold provided here is place-holder. Need to determine actual pixel location"}
                dict[path0 + ".channel_red.description"]     = "red channel description to be provided by Simon"
                dict[path0 + ".channel_red.emission_lambda"] = "675/70 after transmitted component of a Chroma 565DCXR dichroic"
                dict[path0 + ".channel_green.description"]   = "green channel description to be provided by Simon"
                dict[path0 + ".channel_green.emission_lambda"] = "BG22 after reflected component of a Chroma 565DCXR dichroic"
                dict[path0 + ".reference_frame"] = "3.6 mm lateral (left), 1.6 mm posterior Bregma, 0 mm depth"
        return dict

# ------------------------------------------------------------------------------
    
    def processing_extracellular_units(self, data, metadata, options):
        dict = {}
        module_name = "extracellular_units"
        data_item = options.path_str.split(".")[-1]
        if data_item == "top_datasets":
            dict = processing_extracellular_units_top_datasets(data, metadata, options)
        elif data_item == "EventWaveform":
            dict = processing_extracellular_units_EventWaveform(data, metadata, options)
        elif data_item == "UnitTimes":
            dict = processing_extracellular_units_UnitTimes(data, metadata, options)
        return dict
    
# ------------------------------------------------------------------------------
    
    def processing_ROIs(self, data, metadata, options):
        dict = {}
        dict1 = self.processing_ROIs_DfOverF(data, metadata, options)
        dict2 = self.processing_ROIs_ImageSegmentation(data, metadata, options)
        dict.update(dict1)
        dict.update(dict2)
        return dict
    
# ------------------------------------------------------------------------------
    
    def processing_ROIs_DfOverF(self, data, metadata, options):
        dict = {}
        plane_map = util.create_plane_map(data, options)
        num_subareas = len(util.get_key_list(data['timeSeriesArrayHash'])) - 1
        master_shape = util.extract_dict_master_shape(data, plane_map, \
                                                 options, num_planes = 3)
        dict["processing.ROIs.description"] = "Segmentation (pixel-lists) and dF/F (dffTSA) for all ROIs"
        path = "processing.ROIs.DfOverF."
        path2= "processing.ROIs.ImageSegmentation."
        dict[path  + "attrs"] = {"source" : "This module's DfOverF interface"}
        dict[path2 + "attrs"] = {"source" : "This module's ImageSegmentation interface"}
        for subarea in range(num_subareas):
            plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
            if data[plane_path].keys()[0] == 'masterImage':
                num_planes = 1
            else:
                num_planes = len(data[plane_path].keys())
            for plane in range(num_planes):
                tsah = data["timeSeriesArrayHash"]
                oname = "area%d_plane%d" % (subarea+1, plane+1)
                master1_shape = master_shape[oname]
                image_plane = plane_map[oname]
                area_grp = data["timeSeriesArrayHash/value"]["%d"%(subarea+2)]
                try:
                    plane_ids = area_grp["imagingPlane"][str(plane+1)]["ids/ids"].value
                except:
                    plane_ids = area_grp["imagingPlane"]["ids/ids"].value
                roi_names = []
                for i in range(len(plane_ids)):
                    roi_names.append("ROI%d"%plane_ids[i])
                dict[path + image_plane + ".roi_names"] = roi_names
                dict[path + image_plane + ".data"] = area_grp["valueMatrix/valueMatrix"].value
                dict[path + image_plane + ".data.attrs"] = \
                    {"unit":"dF/F", "conversion":1.0, "resolution":0.0}
                dict[path + image_plane + ".timestamps"] = area_grp["time/time"].value * 0.001
                dict[path + image_plane + ".trial_ids"]  = area_grp["trial/trial"].value
        return dict
    
# ------------------------------------------------------------------------------
    
    def processing_ROIs_ImageSegmentation(self, data, metadata, options):
        dict = {}
        plane_map = util.create_plane_map(data, options)
        num_subareas = len(util.get_key_list(data['timeSeriesArrayHash'])) - 1
        master_shape, ref_images_red, ref_images_green = \
            util.create_reference_images(plane_map, num_subareas, data, options)
        path = "processing.ROIs.ImageSegmentation."
        dict[path + "attrs"] = {"source" : "This module's ImageSegmentation interface"}
        for subarea in range(num_subareas):
            plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
            if data[plane_path].keys()[0] == 'masterImage':
                num_planes = 1
            else:
                num_planes = len(data[plane_path].keys())
            for plane in range(num_planes):
                tsah = data["timeSeriesArrayHash"]
                oname = "area%d_plane%d" % (subarea+1, plane+1)
                master1_shape = master_shape[oname]
                image_plane = plane_map[oname]
                ref_image_red   = ref_images_red[image_plane]
                ref_image_green = ref_images_green[image_plane]
                dict[path + image_plane + ".imaging_plane_name"] = image_plane
                dict[path + image_plane + ".description"]        = image_plane
                dict[path + image_plane + ".ref_image_red"]      = ref_image_red
                dict[path + image_plane + ".ref_image_green"]    = ref_image_green
                dict[path + image_plane + ".ref_image_red.attrs"]   = {"source" : "%s - red"   % oname}
                dict[path + image_plane + ".ref_image_green.attrs"] = {"source" : "%s - green" % oname}
                try:
                    ids = tsah["value"]["%d"%(subarea+2)]["imagingPlane"][str(plane+1)]["ids"]
                except:
                    ids = tsah["value"]["%d"%(subarea+2)]["imagingPlane"]["ids"]
                roi_ids = ids["ids"].value
                dict[path + image_plane + ".roi_list"] = np.array(roi_ids)
                for i in range(len(roi_ids)):
                    rid = roi_ids[i]
                    if num_planes == 1:
                        rois = tsah["descrHash"]["%d"%(subarea+2)]["value"]["1"]
                    else:
                        rois = tsah["descrHash"]["%d"%(subarea+2)]["value"]["1"]["%d"%(plane+1)]
                    try:
                        record = rois["rois"]["%s"%(1+i)]
                        x = int(util.parse_h5_obj(record["id"])[0])
                        assert x == int(roi_ids[i]) # make sure the ROI id is correct
                    except:
                        print("Missing ROI for area=" +str(subarea)+ " plane="+ str(plane) + " id=" +str(i))
                        continue
                    pix = util.parse_h5_obj(record["indicesWithinImage"])[0]
                    pixmap = []
                    for j in range(len(pix)):
                        v = pix[j]
                        px = int(v  / master1_shape[1])
                        py = int(v) % master1_shape[0]
                        pixmap.append([py,px])
                    dict[path + image_plane + "." + str(rid) + ".x"]      = x
                    dict[path + image_plane + "." + str(rid) + ".weight"] = np.zeros(len(pixmap)) + 1.0
                    dict[path + image_plane + "." + str(rid) + ".pixmap"] = np.array(pixmap).astype('uint16')
                    dict[path + image_plane + "." + str(rid) + ".master1_shape"] = master1_shape
        return dict
    
# ------------------------------------------------------------------------------
    
    def processing_extracellular_units_EventWaveform(self, data, metadata, options):
        dict = {}
        unit_descr = np.array(util.parse_h5_obj(data['eventSeriesHash/descr/descr'])[0]).tolist()
        unit_num = len(unit_descr)
    
        # Populating a dictionary
        series_path = "processing.extracellular_units.EventWaveform"
        group_attrs = {"source" :  "EventWaveform in this module"} 
        dict[series_path + ".attrs"] = group_attrs
        
        for i in range(unit_num):
            i = i+1
            if options.verbose:
                print "i=", i, " unit_num=", unit_num
            unit = "unit_%d%d" % (int(i/10), i%10)
            if unit_num > 1:
                grp_name = "eventSeriesHash/value/%d" % i
            grp_top_folder = data[grp_name]
            timestamps     = grp_top_folder["eventTimes/eventTimes"]
            trial_ids      = grp_top_folder["eventTrials/eventTrials"]
            waveforms = grp_top_folder["waveforms/waveforms"]
            sample_length = waveforms.shape[1]
            channel = grp_top_folder["channel/channel"].value
    
            cell_type = util.parse_h5_obj(grp_top_folder["cellType"])[0]
            if  'numpy' in str(type(cell_type)):
                cells_conc = ' and '.join(map(str, cell_type))
                cell_type = cells_conc
            else:
                cell_type = str(cell_type)
    
            data_attrs = {"description" : cell_type, "resolution":float('nan'), "unit":"Volts", "conversion":0.1}
    
            dict[series_path + "." + unit + ".data"]          = waveforms
            dict[series_path + "." + unit + ".data.attrs"]    = data_attrs
            dict[series_path + "." + unit + ".timestamps"]    = timestamps
            dict[series_path + "." + unit + ".sample_length"] = len(timestamps)
            dict[series_path + "." + unit + ".electrode_idx"] = [channel[0]]
            dict[series_path + "." + unit + ".num_samples"]   = len(timestamps)
    
        return dict
    
# ------------------------------------------------------------------------------
    
    def processing_extracellular_units_UnitTimes(self, data, metadata, options):
        dict = {}
    
        unit_descr = np.array(util.parse_h5_obj(data['eventSeriesHash/descr/descr'])[0]).tolist()
        unit_num = len(unit_descr)
        grp_name = "eventSeriesHash/value"
    
        cell_types = ['unclassified']*unit_num
        n = max(range(unit_num)) + 1
        electrode_depths = np.zeros(n)
    
        # Populating a dictionary
        series_path = "processing.extracellular_units.UnitTimes"
        group_attrs = {"source" :  "EventWaveform in this module"}
        dict[series_path + ".attrs"] = group_attrs
    
        series_path = "processing.extracellular_units.UnitTimes"
        for i in range(unit_num):
            i = i+1
            if options.verbose:
                print "i=", i, " unit_num=", unit_num
            unit = "unit_%d%d" % (int(i/10), i%10)
            if unit_num > 1:
                grp_name = "eventSeriesHash/value/%d" % i
            grp_top_folder = data[grp_name]
            timestamps     = grp_top_folder["eventTimes/eventTimes"]
            trial_ids      = grp_top_folder["eventTrials/eventTrials"]
    
            cell_type = util.parse_h5_obj(grp_top_folder["cellType"])[0]
            if  'numpy' in str(type(cell_type)):
                cells_conc = ' and '.join(map(str, cell_type))
                cell_type = cells_conc
            else:
                cell_type = str(cell_type)
            cell_types[i-1] = unit + " - " + cell_type
    
            dict[series_path + "." + unit + ".times"]  = timestamps
            dict[series_path + "." + unit + ".source"] = "Data from processed matlab file"
            dict[series_path + "." + unit + ".trial_ids"]        = trial_ids
            dict[series_path + "." + unit + ".unit_description"] = cell_type
    
        dict[series_path + ".cell_types"]            = cell_types
    
        return dict
    
# ------------------------------------------------------------------------------
    
    def processing_extracellular_units_top_datasets(self, data, metadata, options):
        dict = {}
        description = 'Spike times and waveforms'
        try:
            spike_sorting = util.get_value_pointer_by_path_items(metadata, ["extracellular", \
                               "spikeSorting", "spikeSorting"])
        except:
            spike_sorting = "method not recorded"
    
        ident_meth1 = str(np.array(util.get_value_pointer_by_path_items(metadata, ["extracellular", \
                                   "identificationMethod", "identificationMethod"])).tolist()[0])
        ident_meth2 = str(np.array(util.get_value_pointer_by_path_items(metadata, ["extracellular", \
                                   "identificationMethod", "identificationMethod"])).tolist()[2])
        ident_meth  = ident_meth1 + "\n" + ident_meth2
    
        # Populating a dictionary
        series_path = "processing.extracellular_units"
        dict[series_path + ".description"]                = description
        dict[series_path + ".spike_sorting"]              = spike_sorting
        dict[series_path + ".identification_method"]      = ident_meth
    
        return dict
    
# ------------------------------------------------------------------------------
    
    def processing_Whisker_BehavioralTimeSeries_whisker_angle(self, data, metadata, options):
        dict = {}
    
        if options.data_type == "ophys":
            keyName = 'whiskerVars'
            hash_group_pointer  = data['timeSeriesArrayHash']
            grp        = util.get_value_by_key(hash_group_pointer , keyName)
            timestamps = grp["time/time"].value * 0.001
    
            var   = grp["valueMatrix/valueMatrix"].value
            if options.handle_errors:
                try:
                    grp = data["timeSeriesArrayHash/value/1"]
                except:
                    grp = data["timeSeriesArrayHash/value"]
            else:
                grp = data["timeSeriesArrayHash/value/1"]
            descr = util.parse_h5_obj(grp["idStrs"])[0]
    
            data = var[0]
            if options.verbose:
                print "descr[0]=", descr[0]
            group_attrs={"description": "Angle of whiskers", \
                         "source": "Whisker angle as reported in Simon's data file",\
                         "series_type" : "<TimeSeries>"}
            data_attrs ={"unit": "degrees", "conversion":1.0, "resolution":0.001,
                         "keyName" : keyName}
    
       # Populating a dictionary
        series_path = "processing.Whisker.BehavioralTimeSeries.whisker_angle"
        dict[series_path + ".attrs"]       = group_attrs
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
    
        return dict
    
# ------------------------------------------------------------------------------
   
    def processing_Whisker_BehavioralTimeSeries_whisker_curvature(self, data, metadata, options):
        dict = {}

        if options.data_type == "ophys":
            keyName = 'whiskerVars'
            hash_group_pointer  = data['timeSeriesArrayHash']
            grp        = util.get_value_by_key(hash_group_pointer , keyName)
            timestamps = grp["time/time"].value * 0.001

            var   = grp["valueMatrix/valueMatrix"].value
            if options.handle_errors:
                try:
                    grp = data["timeSeriesArrayHash/value/1"]
                except:
                    grp = data["timeSeriesArrayHash/value"]
            else:
                grp = data["timeSeriesArrayHash/value/1"]
            descr = util.parse_h5_obj(grp["idStrs"])[0]

            data = var[1]
            if options.verbose:
                print "descr[1]=", descr[1]
            group_attrs={"description": descr[1],
                         "source": "Curvature (relative) of whiskers as reported in Simon's data file", \
                         "series_type" : "<TimeSeries>"}
            data_attrs={"unit":"Unknown", "conversion": 1.0, "resolution": 1.0,
                        "keyName" : keyName}

       # Populating a dictionary
        series_path = "processing.Whisker.BehavioralTimeSeries.whisker_curvature"
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs

        return dict

# ------------------------------------------------------------------------------
 
    def processing_Whisker_BehavioralEvents_pole_touch_protract(self, data, metadata, options):
        dict = {}
        keyName = "touches"
        hash_group_pointer  =  data["eventSeriesArrayHash"]
        grp = util.get_value_by_key(hash_group_pointer , keyName)
    
        timestamps = grp["eventTimes/1/1"].value * 0.001
        data = util.ones_with_alternating_sign(timestamps)
        group_attrs = {"description" : "Intervals that whisker touches pole (protract)",\
                       "source" : "Intervals are as reported in Simon's data file",
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}
   
        # Populating a dictionary
        series_path = "processing.Whisker.BehavioralEvents.pole_touch_protract"
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs
    
        return dict
    
# ------------------------------------------------------------------------------

    def processing_Whisker_BehavioralEvents_pole_touch_retract(self, data, metadata, options):
        dict = {}
        keyName = "touches"
        hash_group_pointer  =  data["eventSeriesArrayHash"]
        grp = util.get_value_by_key(hash_group_pointer , keyName)

        timestamps = grp["eventTimes/2/2"].value * 0.001
        data = util.ones_with_alternating_sign(timestamps)
        group_attrs = {"description" : "Intervals that whisker touches pole (retract)",\
                       "source" : "Intervals are as reported in Simon's data file",
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}

        # Populating a dictionary
        series_path = "processing.Whisker.BehavioralEvents.pole_touch_retract"
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs

        return dict

# ------------------------------------------------------------------------------
    
    def processing_Reward_BehavioralEvents_water_left_reward(self, data, metadata, options):
        dict = {}
        hash_group_pointer  = data['eventSeriesArrayHash']
        source = "Intervals are as reported in somatosensory cortex data file"
    
        keyName = "leftReward"
        description = util.get_description_by_key(hash_group_pointer, keyName)
        grp = util.get_value_by_key(hash_group_pointer , keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = util.ones_with_alternating_sign(timestamps)
        group_attrs = {"source" : source, "description" : description, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}
    
        # Populating a dictionary
        series_path = "processing.Reward.BehavioralEvents.water_left_reward"
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs
    
        return dict
    
# ------------------------------------------------------------------------------
   
    def processing_Reward_BehavioralEvents_water_right_reward(self, data, metadata, options):
        dict = {}
        hash_group_pointer  = data['eventSeriesArrayHash']
        source = "Intervals are as reported in somatosensory cortex data file"

        keyName = "rightReward"
        description = util.get_description_by_key(hash_group_pointer , keyName)
        grp = util.get_value_by_key(hash_group_pointer , keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = util.ones_with_alternating_sign(timestamps)
        group_attrs = {"source" : source, "description" : description, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}

        # Populating a dictionary
        series_path = "processing.Reward.BehavioralEvents.water_right_reward"
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs

        return dict

# ------------------------------------------------------------------------------
 
    def processing_Pole_BehavioralEvents_pole_accessible(self, data, metadata, options):
        dict = {}
        if options.data_type == "ophys":
            keyName = "poleInReach"
            hash_group_pointer = data['eventSeriesArrayHash']
            grp = util.get_value_pointer_by_key(hash_group_pointer , keyName, \
                                                  options.debug)
            timestamps = grp["eventTimes/eventTimes"].value * 0.001
            data = [1] * len(timestamps)
            description = util.get_description_by_key(hash_group_pointer , keyName)
            group_attrs = {"source" : "Intervals are as reported in Simon's data file", \
                           "description" : description, "series_type" : "<IntervalSeries>"}
            data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                           "keyName" : keyName}
    
            # Populating a dictionary
            series_path = "processing.Pole.BehavioralEvents.pole_accessible"
            dict[series_path + ".timestamps"]  = timestamps
            dict[series_path + ".data"]        = data
            dict[series_path + ".num_samples"] = len(timestamps)
            dict[series_path + ".data.attrs"]  = data_attrs
            dict[series_path + ".attrs"]       = group_attrs
    
        return dict
    
# ------------------------------------------------------------------------------
    
    def processing_Licks_BehavioralEvents_lick_left(self, data, metadata, options):
        dict = {}
        hash_group_pointer  = data['eventSeriesArrayHash']
    
        keyName = 'leftLicks'
        description = util.get_description_by_key(hash_group_pointer, keyName)
        grp = util.get_value_by_key(hash_group_pointer, keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = [1] * len(timestamps)
        data_attrs = {"description": description, \
                      "source": "Times as reported in Simon's data file", \
                      "unit":"licks", "conversion": 1.0, "resolution": 1.0, \
                      "keyName" : keyName}
        group_attrs = {"description": "Left lickport contact times (beam breaks left)",\
                       "source": "Times as reported in Simon's data file",\
                       "comments": "Timestamp array stores lick times", \
                       "series_type" : "<TimeSeries>"}
    
        # Populating a dictionary
        series_path = "processing.Licks.BehavioralEvents.lick_left"
        dict["processing.Licks"]           = {"description" : "Lick port contact times"}
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs
        return dict
    
# ------------------------------------------------------------------------------

    def processing_Licks_BehavioralEvents_lick_right(self, data, metadata, options):
        dict = {}
        hash_group_pointer  = data['eventSeriesArrayHash']

        keyName = 'rightLicks'
        description = util.get_description_by_key(hash_group_pointer, keyName)
        grp = util.get_value_by_key(hash_group_pointer, keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = [1] * len(timestamps)
        data_attrs = {"description": description, \
                      "source": "Times as reported in Simon's data file",
                      "unit":"licks", "conversion": 1.0, "resolution": 1.0, \
                      "keyName" : keyName}
        group_attrs = {"description": "Right lickport contact times (beam breaks right)",
                       "source": "Times as reported in Simon's data file",
                       "comments": "Timestamp array stores lick times", \
                       "series_type" : "<TimeSeries>"}

        # Populating a dictionary
        series_path = "processing.Licks.BehavioralEvents.lick_right"
        dict["processing.Licks"]           = {"description" : "Lick port contact times"}
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs
        return dict

# ------------------------------------------------------------------------------
    
    def processing_Auditory_BehavioralEvents_reward_cue(self, data, imetadata, options):
        dict = {}
        if options.data_type == "ophys":
            hash_group_pointer  = data['eventSeriesArrayHash']
            keyName = "rewardCue"
            description = util.get_description_by_key(hash_group_pointer, keyName)
            grp = util.get_value_by_key(hash_group_pointer, keyName)
            timestamps = grp["eventTimes/eventTimes"].value * 0.001
            data = [1] *len(timestamps)
            group_attrs = {"description" : "Intervals when auditory cue presented",\
                           "source": "Intervals are as reported in Simon's data file", \
                           "series_type" : "<IntervalSeries>"}
            data_attrs = {"unit": "None", "conversion": 1.0, "resolution": 0.0, \
                          "keyName" : keyName}
    
            # Populating a dictionary
            series_path = "processing.Auditory.BehavioralEvents.reward_cue"
            dict[series_path + ".timestamps"]  = timestamps
            dict[series_path + ".data"]        = data
            dict[series_path + ".num_samples"] = len(timestamps)
            dict[series_path + ".data.attrs"]  = data_attrs
            dict[series_path + ".attrs"]       = group_attrs
        return dict
    
# ------------------------------------------------------------------------------
    
    def stimulus_presentation_zaber_motor_pos(self, data, metadata, options):
        dict = {}
    
        keyName = "StimulusPosition"
        timestamps = data["trialStartTimes/trialStartTimes"].value * 0.001
        stim_pos = util.get_value_by_key(data['trialPropertiesHash'], keyName)
        description = util.get_description_by_key(data["trialPropertiesHash"], keyName)
        group_attrs = {"source" : "Simon's somatosensory cortex data file", "description" : description, \
                       "series_type" : "<TimeSeries>"}
        data_attrs  = {"unit":"unknown", "conversion": 1.0, "resolution":1.0, \
                       "keyName" : keyName}
    
        # Populating a dictionary
        series_path = "stimulus.presentation.zaber_motor_pos"
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = stim_pos
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs
    
        return dict
    
# ------------------------------------------------------------------------------
    
    def stimulus_presentation_water_left(self, data, metadata, options):
        dict = {}
        hash_group_pointer  = data['eventSeriesArrayHash']
        source = "Intervals are as reported in somatosensory cortex data file"
    
        keyName = "leftReward"
        description = util.get_description_by_key(hash_group_pointer, keyName)
        grp = util.get_value_by_key(hash_group_pointer , keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = util.ones_with_alternating_sign(timestamps)
        group_attrs = {"source" : source, "description" : description, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}
    
        # Populating a dictionary
        series_path = "stimulus.presentation.water_left"           
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs
    
        return dict
    
# ------------------------------------------------------------------------------

    def stimulus_presentation_water_right(self, data, metadata, options):
        dict = {}
        hash_group_pointer  = data['eventSeriesArrayHash']
        source = "Intervals are as reported in somatosensory cortex data file"

        keyName = "rightReward"
        description2 = util.get_description_by_key(hash_group_pointer , keyName)
        grp = util.get_value_by_key(hash_group_pointer , keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = util.ones_with_alternating_sign(timestamps)
        group_attrs = {"source" : source, "description" : description2, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0}

        # Populating a dictionary
        series_path = "stimulus.presentation.water_right"     
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs

        return dict

# ------------------------------------------------------------------------------

    def stimulus_presentation_pole_in(self, data, metadata, options):
        dict = {}

        hash_group_pointer = data["trialPropertiesHash"]
        source = "Times as reported in motor cortex data file, but relative to session start"
        trial_start_times = data["trialStartTimes/trialStartTimes"].value
        grp = data["trialPropertiesHash/value/"]
        data_attrs  = {"resolution":0.1,"conversion":1.0}

        keyName = "PoleInTime"
        time = grp["1/1"].value
        timestamps = time + trial_start_times
        data = [1] * len(timestamps)
        description = util.get_description_by_key(hash_group_pointer, keyName)
        data_attrs  = {"keyName" : keyName, "unit" : "sec"}
        group_attrs = {"source": source, "description" : description, \
                           "series_type" : "<TimeSeries>"}

        # Populating a dictionary
        series_path = "stimulus.presentation.pole_in"               
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs

        return dict

# ------------------------------------------------------------------------------

    def stimulus_presentation_pole_accessible(self, data, metadata, options):
        dict = {}

        # create time series for stimulus/presentation group
        keyName = "poleInReach"
        hash_group_pointer = data['eventSeriesArrayHash']
        grp = util.get_value_pointer_by_key(hash_group_pointer , keyName, \
                                              options.debug)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = util.ones_with_alternating_sign(timestamps)
        description = util.get_description_by_key(hash_group_pointer , keyName)
        group_attrs = {"source" : "Intervals are as reported in Simon's data file", \
                       "description" : description, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, "keyName" : keyName}

        # Populating a dictionary
        series_path = "stimulus.presentation.pole_accessible"
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs

        return dict

# ------------------------------------------------------------------------------
    
    def stimulus_presentation_pole_out(self, data, metadata, options):
        dict = {}
    
        hash_group_pointer = data["trialPropertiesHash"]
        source = "Times as reported in motor cortex data file, but relative to session start"
        trial_start_times = data["trialStartTimes/trialStartTimes"].value
        grp = data["trialPropertiesHash/value/"]
        data_attrs  = {"resolution":0.1,"conversion":1.0}
   
        keyName = "PoleOutTime"
        time = grp["2/2"].value
        timestamps = time + trial_start_times
        data = [1] * len(timestamps)
        description = util.get_description_by_key(hash_group_pointer, keyName)
        data_attrs  = {"keyName" : keyName, "unit" : "sec"}
        group_attrs = {"source": source, "description" : description, \
                       "series_type" : "<TimeSeries>"}
    
        # Populating a dictionary
        series_path = "stimulus.presentation.pole_out" 
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs
    
        return dict
    
# ------------------------------------------------------------------------------
    
    def stimulus_subtrace(orig_trace, time, trial_t, types, types_to_remove, options):
        new_trace = orig_trace
        if options.verbose:
            print("    num_stimulus_types=" + str(len(types)) + " types=" + str(types))
    
        for i in range(len(types)):
            if str(types[i]) in types_to_remove:
                print("    Setting to zero values for type=" + str(types[i]))
                start = trial_t[i]
                stop  = trial_t[i+1]
                new_trace[(time >= start) & (time < stop)] = 0.
                print("start=" + str(start) + " stop=" + str(stop) +  "  " + str(((time >= start) & (time < stop)).sum()) + " values has been reset")
        return new_trace
    
# ------------------------------------------------------------------------------
    
    def stimulus_presentation_photostimulus(self, data, metadata, options):
        dict = {}
        grp_name = "timeSeriesArrayHash/value/time/time"
        timestamps = np.array(data[grp_name].value)
        # calculate sampling rate
        rate = (timestamps[-1] - timestamps[0])/(len(timestamps)-1)
        # get descriptions
        comment1 = util.parse_h5_obj(data["timeSeriesArrayHash/keyNames"])[0][0]
        comment2 = util.parse_h5_obj(data["timeSeriesArrayHash/descr"])[0][0]
        comments = comment1 + ": " + comment2
        grp_name = "timeSeriesArrayHash/value/idStrDetailed"
        description = util.parse_h5_obj(data[grp_name])[0]
    
        # laser data
        keyName1 = "EphusVars"
        keyName2 = "laser_power"
        hash_group_pointer = data["timeSeriesArrayHash"]
        laser_power = util.get_value2_by_key2(hash_group_pointer, keyName1, "",keyName2)
    
        num_traces, trace_types, traces = \
            util.stimulus_subtraces_by_type(data, laser_power, timestamps, options)
        group_attrs = {"description" : description[2], "comments" : comments, \
                       "source" : "Nuo's data file", \
                       "series_type" : "<OptogeneticSeries>"}
        data_attrs = {"resolution":float('nan'), "unit":"mW", "conversion":1000.0, \
                      "keyName" : keyName2, "source" : "Nuo's data file"}
    
        if options.verbose:
            print("    num_traces="  + str(num_traces))
            print("    trace_types=" + str(trace_types))
        if num_traces > 0:
            for s in range(num_traces):
                try:
                    # For site_###, copy the data from /general/optogenetics/site_###/location
                    coord = str(util.get_value_by_path_items(metadata, ["photostim", "photostimCoordinates"]))
                    loc   = str(util.get_value_by_path_items(metadata, ["photostim", "photostimLocation"]))
                    series_path = "stimulus.presentation.photostimulus_" + str(trace_types[s])
                    dict[series_path + ".timestamps"]  = timestamps
                    dict[series_path + ".data"]        = traces[s]
                    dict[series_path + ".num_samples"] = len(timestamps)
                    dict[series_path + ".site"] = str(coord + " in mm\natlas location: " + loc)
                    dict[series_path + ".data.attrs"]  = data_attrs
                    dict[series_path + ".attrs"]       = group_attrs
                except:
                    print "No photostimulus info found in the metadata file"
        return dict
    
# ------------------------------------------------------------------------------
    
    def stimulus_presentation_auditory_cue(self, data, metadata, options):
        dict = {}
    
        if options.data_type == "ephys":
            keyName = "CueTime"
            hash_group_pointer  = data['trialPropertiesHash']
            trial_start_times = data["trialStartTimes/trialStartTimes"].value
            timestamps = hash_group_pointer["value/3/3"].value + trial_start_times
            data = util.ones_with_alternating_sign(timestamps)
            description = hash_group_pointer["descr/descr"][2]
            description = util.get_description_by_key(hash_group_pointer, keyName)
            group_attrs = {"comments" : description,\
                           "description" : keyName, \
                           "source" : "Times are as reported in Nuo's data file, but relative to session time", \
                           "series_type" : "<TimeSeries>"}
            data_attrs = {"unit": "None", "conversion": 1.0, "resolution": 0.0, "duration" : 0.1, \
                          "keyName" : keyName}
        elif options.data_type == "ophys":
            keyName = "rewardCue"
            hash_group_pointer  = data['eventSeriesArrayHash']
            description = util.get_description_by_key(hash_group_pointer, keyName)
            grp = util.get_value_by_key(hash_group_pointer, keyName)
            timestamps = grp["eventTimes/eventTimes"].value * 0.001
            data = [1] * len(timestamps)
            group_attrs = {"description" : "Intervals when auditory cue presented",\
                           "source": "Intervals are as reported in Simon's data file", \
                           "series_type" : "<IntervalSeries>"}
            data_attrs = {"unit": "None", "conversion": 1.0, "resolution": 0.0, \
                          "keyName" : keyName}
    
        # Populating a dictionary
        series_path = "stimulus.presentation.auditory_cue"
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data.attrs"]  = data_attrs
        dict[series_path + ".attrs"]       = group_attrs
    
        return dict
    
# ------------------------------------------------------------------------------
    
def make_dict(data, metadata, options):
    method_name = "unknown"
    if len(options.path_str.split(".")) > 0:
        method_name = "_".join(options.path_str.split("."))
    else:
        method_name = options.path_str
    class_instance = Exp2Dict()

    method = None
    try:
        method = getattr(class_instance, method_name)
    except:
        sys.exit("\nMethod libexp2dict." + method_name + "' has not been implemented")

    dict = method(data, metadata, options)

    return dict
