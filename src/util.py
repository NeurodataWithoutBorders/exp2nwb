# Copyright (c) 2017, HHMI-Janelia Research Campus All rights reserved.
#
# util.py
#
# Utility functions for extracting data from HDF5 files 
#
import numpy as np
import re
import os, sys
import nwb
from nwb import nwb_utils
import datetime
import h5py
from sets import Set

verbose = 0    

# ------------------------------------------------------------------------------
def check_entry(file_name,obj):
    try:
        return file_name[obj]
    except KeyError:
        print(str(obj) +" does not exist")
        return []

# ------------------------------------------------------------------------------

def print_dict_keys(dict):
    keys = sorted(dict.keys())
    print "keys=\n"
    for i in range(len(keys)):
        print "   ", keys[i]

# ------------------------------------------------------------------------------

def parse_h5_obj(obj, level = 0, output = [], verbose = 0):
    if level == 0:
        output = []
    try:
        if isinstance(obj, h5py.highlevel.Dataset):
            level = level+1
            if obj.value.any():
                output.append(obj.value)
            else:
                full_name = obj.name.split("/")
                output.append(full_name[-1])
        elif isinstance(obj, h5py.highlevel.Group):
            level = level+1
            if not obj.keys():
                output.append([])
            else:
                for key in obj.keys():
                    parse_h5_obj(obj[key], level, output, verbose)
        else:
            output.append([])
    except KeyError:
        print("Can't find" + str(obj))
        output.append([])
    return output

# ------------------------------------------------------------------------------

def ones_with_alternating_sign(timestamps):
    on_off = np.int_(np.zeros(len(timestamps)))
    on_off += -1
    on_off[::2] *= -1
    data = on_off
    return data

# ------------------------------------------------------------------------------

# Extract the start time/date of experiment
def find_exp_time(input_h5, options):
    child_groups = get_child_group_names(input_h5)
    if options.data_type == "ophys":
        d = get_value_by_key(input_h5['/metaDataHash'], \
                         "dateOfExperiment")
        t = get_value_by_key(input_h5['/metaDataHash'], \
                         "timeOfExperiment")
        dt=datetime.datetime.strptime(d+t, "%Y%m%d%H%M%S")

    elif options.data_type == "ephys":
        d = np.array(get_value_pointer_by_path_items(input_h5, \
                     ["dateOfExperiment", "dateOfExperiment"])).tolist()[0]
        t = np.array(get_value_pointer_by_path_items(input_h5, \
                     ["timeOfExperiment", "timeOfExperiment"])).tolist()[0]
        try:
            if re.search("not recorded", t):
                dt = datetime.datetime.strptime(d, "%Y%m%d")
            else:
                dt = datetime.datetime.strptime(d+t, "%Y%m%d%H%M%S")
        except:
            d = re.sub('\x00','', d)
            t = re.sub('\x00','', t)
            dt = datetime.datetime.strptime(d+t, "%Y%m%d")
    else:
        sys.exit("Data origin is unknown")

    return dt.strftime("%a %b %d %Y %H:%M:%S")

# ------------------------------------------------------------------------------

def compute_age(DOB, DOE):
    years  = 0
    months = 0
    days   = 0
    DOB_year  = int(DOB[:4])
    DOE_year  = int(DOE[:4])
    DOB_month = int(DOB[4:6])
    DOE_month = int(DOE[4:6])
    DOB_day   = int(DOB[6:])
    DOE_day   = int(DOE[6:])
    if DOE_day >= DOB_day:
       days = DOE_day - DOB_day
       if DOE_month >= DOB_month:
           months = DOE_month - DOB_month
           years  = DOE_year  - DOB_year
       else:
           months = DOE_month - DOB_month + 12
           years  = DOE_year  - DOB_year  - 1
    else:
       days = DOE_month - DOB_month + 30
       if DOE_month >= DOB_month + 1:
           months = DOE_month - DOB_month - 1
           years  = DOE_year  - DOB_year
       else:
           months = DOE_month - DOB_month - 1 + 12
           years  = DOE_year  - DOB_year  - 1
    age = ""
    if years > 0:
        age += str(years)  + " years "
    if months > 0:
        age += str(months) + " months "
    if days > 0:
        age += str(days)   + " days "
    return age

# ------------------------------------------------------------------------------

def metadata_from_file(file_name, options):
    if os.environ['NWB_DATA'] is not None:
        file_path = os.path.join(os.environ['NWB_DATA'],file_name)
    else:
        data_dir = os.path.join(",".join(options.data_path.split("/")[0:-1]))
        file_path = os.path.join(data_dir, file_name)
    description = " "
    if os.path.isfile(file_path):
        description = nwb_utils.load_file(file_path)
    else:
        print "Warning: missing file " + file_name
    return description

# ------------------------------------------------------------------------------

def get_description(metadata, options):
    extra_grp = get_value_pointer_by_path_items(metadata, \
                                ["extracellular"])
    value = ""
    for key in extra_grp.keys():
        if key in ["recordingType", "penetrationN", "groundCoordinates","referenceCoordinates", \
                   "extracellularDataType", "cellType", "identificationMethod","spikeSorting"]:
            try:
                # getting a dictionary
                key1_list = get_value_pointer_by_path_items(extra_grp, [key, key]).keys()
                for key1 in key1_list:
                    value1 = get_value_pointer_by_path_items(extra_grp, [key, key1])[:]
                    value2 = [str(v) for v in value1]
                    value += "       " + key1 + ": " + ",".join(value2) + "\n "
            except:
                # getting a dataset
                value = get_value_pointer_by_path_items(extra_grp, [key, key])
    return value

# ------------------------------------------------------------------------------

def get_device(metadata, options):
    extra_grp = get_value_pointer_by_path_items(metadata, \
                                ["extracellular"])
    value = ""
    for key in extra_grp.keys():
        if key in ["probeSource", "probeType", "ADunit", "amplifierRolloff"]:
            key1_list = get_value_pointer_by_path_items(extra_grp, [key]).keys()
            for key1 in key1_list:
                value1 = get_value_pointer_by_path_items(extra_grp, [key, key1])[:]
                value2 = [str(v) for v in value1]
                value += "       " + key1 + ": " + ",".join(value2) + "\n "
    return value

# ------------------------------------------------------------------------------

def read_probe_locations_matrix(metadata, options):
    num_locs = len(get_value_pointer_by_path_items(metadata, \
                   ["extracellular", "siteLocations"]).keys())
    M = np.zeros([num_locs, 3])
    for i in range(num_locs):
        probe_id = i + 1
        coords = np.array(get_value_pointer_by_path_items(metadata, \
                   ["extracellular", "siteLocations", \
                    str(probe_id), str(probe_id)])).tolist()
        M[i][:] = coords[:]
    return M

# ------------------------------------------------------------------------------

def detect_shanks(M):
    P = np.zeros([M.shape[0], 1]) # vector of already processed row ids
    x1 = M[0,0]
    y1 = M[0,1]
    cols1 = np.where(M[:, 0] == x1)
    cols2 = np.where(M[:, 1] == y1)
    cols = np.intersect1d(cols1, cols2)
    shank_size = np.prod(cols.shape)
    curr_shank_id = 0
    num_shanks = 0
    shank_coords = []
    for i in range(M.shape[0]):
        if P[i] == 1:
            continue
        shank_coords.append([M[i, 0], M[i, 1]])
        cols1 = np.where(M[:,0] == M[i, 0])
        cols2 = np.where(M[:,1] == M[i, 1])
        cols = np.intersect1d(cols1, cols2)
        P[cols] = 1
        num_shanks = num_shanks + 1
    return (num_shanks, shank_size, shank_coords)

# ------------------------------------------------------------------------------

def get_description(metadata, options):
    extra_grp = get_value_pointer_by_path_items(metadata, \
                                ["extracellular"])
    value = ""
    for key in extra_grp.keys():
        if key in ["recordingType", "penetrationN", "groundCoordinates","referenceCoordinates", \
                   "extracellularDataType", "cellType", "identificationMethod","spikeSorting"]:
            try:
                # getting a dictionary
                key1_list = get_value_pointer_by_path_items(extra_grp, [key, key]).keys()
                for key1 in key1_list:
                    value1 = get_value_pointer_by_path_items(extra_grp, [key, key1])[:]
                    value2 = [str(v) for v in value1]
                    value += "       " + key1 + ": " + ",".join(value2) + "\n "
            except:
                # getting a dataset
                value = get_value_pointer_by_path_items(extra_grp, [key, key])
    return value

# ------------------------------------------------------------------------------

def get_device(metadata, options):
    extra_grp = get_value_pointer_by_path_items(metadata, \
                                ["extracellular"])
    value = ""
    for key in extra_grp.keys():
        if key in ["probeSource", "probeType", "ADunit", "amplifierRolloff"]:
            key1_list = get_value_pointer_by_path_items(extra_grp, [key]).keys()
            for key1 in key1_list:
                value1 = get_value_pointer_by_path_items(extra_grp, [key, key1])[:]
                value2 = [str(v) for v in value1]
                value += "       " + key1 + ": " + ",".join(value2) + "\n "
    return value

# ------------------------------------------------------------------------------

def define_manifold(master_shape, plane_map):
    manifold = np.zeros((master_shape[0], master_shape[1], 3), dtype=np.float32)
    VV  = int(plane_map[4:6])
    PPP = int(plane_map[6:9])
    for i in range(master_shape[0]):
        for j in range(master_shape[0]):
            manifold[i][j][0] = 1.0 * np.float(j) * (600./512.)
            manifold[i][j][1] = 1.0 * np.float(i) * (600./512.)
            if VV < 10:
                manifold[i][j][2] = 45*(VV- 1) + 15*(PPP-2)
            else:
                manifold[i][j][2] = 45*(VV-10) + 15*(PPP-2)
    return manifold

# ------------------------------------------------------------------------------

def get_valid_trials(data, var, options):
    ts_path = "timeSeriesArrayHash/descrHash/"
    val = []
    if options.handle_errors:
        num_subareas = 1
        if '1' in data['timeSeriesArrayHash/descrHash'].keys():
            num_subareas = len(data['timeSeriesArrayHash/descrHash'].keys()) - 1
    else:
        num_subareas = len(data['timeSeriesArrayHash/descrHash'].keys()) - 1

    if var == "whisker":
        if options.handle_errors:
            try:
                ids = parse_h5_obj(data[ts_path + '1/value'])[0]
            except:
                ids = parse_h5_obj(data[ts_path + 'value'])[0]
        else:
            ids = parse_h5_obj(data[ts_path + '1/value'])[0]
        val = val + list(ids)
        val = list(Set(val))
    if var == "Ca":
        for i in range(2,num_subareas+1):
            ids_path = ts_path + "%d/value/2" % i
            ids = parse_h5_obj(data[ids_path])[0]
            val = val + list(ids)
        val = list(Set(val))
    return val

# ------------------------------------------------------------------------------

# collect unit information for a given trial
epoch_units = {}
def get_trial_units(data_h5, unit_num, options):
    for i in range(unit_num):
        i = i+1
        unit = "unit_%d%d" % (int(i/10), i%10)
        if unit_num > 0:
            grp_name = "eventSeriesHash/value/%d" % i
        else:
            grp_name = "eventSeriesHash/value"
        grp_top_folder = data_h5[grp_name]
        trial_ids = grp_top_folder["eventTrials/eventTrials"].value
        trial_ids = Set(trial_ids)
        for trial_num in trial_ids:
            tid = trial_num
            trial_name = "trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
            if trial_name not in epoch_units.keys():
                epoch_units[trial_name] = []
            epoch_units[trial_name].append(unit)

 # ------------------------------------------------------------------------------


# add trial types to epoch for indexing
epoch_tags = {}
epoch_trial_types = {}
def get_trial_types(data_h5, options):
    trial_id = data_h5["trialIds/trialIds"].value
    trial_types_all = []
    trial_type_strings = parse_h5_obj(data_h5['trialTypeStr'])[0]
    num_trial_types    = len(trial_type_strings)
    if options.data_type == "ephys":
        photostim_types = get_value_by_key(data_h5["trialPropertiesHash"], "PhotostimulationType")
    elif options.data_type == "ophys":
        valid_whisker = get_valid_trials(data_h5, "whisker", options)
        valid_Ca      = get_valid_trials(data_h5, "Ca", options)

    # collect all trials (strings)
    for i in range(num_trial_types):
        trial_types_all.append(str(trial_type_strings[i]))
    # write specific entries for the given trial
    for i in range(len(trial_id)):
        tid = trial_id[i]
        trial_name = "trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
        epoch_tags[trial_name] = []
        trial_types = []
        trial_type_mat = parse_h5_obj(data_h5['trialTypeMat'])[0]
        for j in range(num_trial_types):
            if trial_type_mat[j,i] == 1:
                epoch_tags[trial_name].append(trial_types_all[j])
                trial_types.append(trial_types_all[j])
        if options.data_type == "ephys":
            ps_type_value = photostim_types[i]
            if ps_type_value == 0:
                photostim_type = "non-stimulation trial"
            elif ps_type_value == 1:
                photostim_type = "PT axonal stimulation"
            elif ps_type_value == 2:
                photostim_type = "IT axonal stimulation"
            else:
                photostim_type = "discard"
            epoch_tags[trial_name].append(photostim_type)
        elif options.data_type == "ophys":
            if i in valid_whisker:
                trial_types.append("Valid whisker data")
            else:
                trial_types.append("Invalid whisker data")
            if i in valid_Ca:
                trial_types.append("Valid Ca data")
            else:
                trial_types.append("Invalid Ca data")
            epoch_trial_types[trial_name] = trial_types

# ------------------------------------------------------------------------------

def create_trial_roi_map(orig_h5, plane_map, options):
    epoch_roi_list   = {}
    epoch_roi_planes = {}
    num_subareas = len(orig_h5['timeSeriesArrayHash/descrHash'].keys()) - 1
    for i in range(num_subareas):
        area = i + 1
        grp_name = "timeSeriesArrayHash/value/%d" % (area + 1)
        block = orig_h5[grp_name]
        trials = block["trial/trial"].value
        ids = block["ids/ids"].value
        # create way to map ROI onto plane
        planemap = {}
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(area + 1)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        if num_planes > 1:
            print("  Warning: in create_trial_roi_map num_planes == 1")

        for j in range(num_planes):
            plane = j + 1
            try:
                grp_name = "imagingPlane/%d/ids/ids" % plane
                plane_id = block[grp_name].value
            except:
                print("Warning: only one imaging plane detected for area " + str(area))
                grp_name = "imagingPlane/ids/ids"
                plane_id = block[grp_name].value
            for k in range(len(plane_id)):
                planemap["%d"%plane_id[k]] = "%d" % plane
        trial_list = {}
        for j in range(len(trials)):
            tid = trials[j]
            name = "trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
            if name not in trial_list:
                trial_list[name] = j
        for trial_name in trial_list.keys():
            roi_list = []
            plane_list = []
            for k in range(len(ids)):
                roi = "%d" % ids[k]
                plane = int(planemap[roi])
                # convert from file-specific area/plane mapping to
                #   inter-session naming convention
                #imaging_plane = "area%d_plane%d" % (area, plane)
                oname = "area%d_plane%d" % (area, plane)
                imaging_plane = plane_map[oname]
                s = "processing/ROIs/DfOverF/%s/%s" % (imaging_plane, roi)
                roi_list.append(ids[k])
                plane_list.append(imaging_plane)
            epoch_roi_list[trial_name] = roi_list
            epoch_roi_planes[trial_name] = plane_list
    return (epoch_roi_list, epoch_roi_planes)

# ------------------------------------------------------------------------------

def create_trials(dict, data_h5, options):
    trial_id = data_h5["trialIds/trialIds"].value
    if options.verbose:
        print("\nCreating trials with ids: " + str(trial_id))
    trial_t  = data_h5["trialStartTimes/trialStartTimes"].value
    # trial stop isn't stored. assume that it's twice the duration of other
    #   trials -- padding on the high side shouldn't matter
    ival     = (trial_t[-1] - trial_t[0]) / (len(trial_t) - 1)
    trial_t  = np.append(trial_t, trial_t[-1] + 2*ival)

    if options.data_type == "ephys":
        good_trials = get_value_by_key(data_h5["/trialPropertiesHash"],  "GoodTrials")
        ignore_ivals_start = [time for (time, good_trial) in zip(trial_t,good_trials) if good_trial == 0]
        # trial stop isn't stored. assume that it's twice the duration of other
        #   trials -- padding on the high side shouldn't matter
        ignore_ivals_stop = [time for (time, good_trial) in zip(trial_t[1:],good_trials) if good_trial == 0]
        ignore_intervals  = [ignore_ivals_start, ignore_ivals_stop]

        keyName3 = "PhotostimulationType"
        hash_group_pointer2 = data_h5["/trialPropertiesHash"]
        stimulus_types = np.array(get_value_by_key(hash_group_pointer2, keyName3)).tolist()
        count_1 = stimulus_types.count(1)
        count_2 = stimulus_types.count(2)
    elif options.data_type == "ophys":
        plane_map = create_plane_map(data_h5, options)
        epoch_roi_list,epoch_roi_planes = create_trial_roi_map(data_h5, plane_map, options)
    for i in range(len(trial_id)):
        tid = trial_id[i]
        trial = "trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
        dict["epochs." + trial + ".description"] = "Data that belong to " + trial
        if options.data_type == "ephys":
            start = trial_t[i]
            stop  = trial_t[i+1]
            dict["epochs." + trial + ".start_time"] = start
            dict["epochs." + trial + ".stop_time"]  = stop
            tags = []
            if good_trials[i] == 1:
                tags.append("Good trial")
            else:
                tags.append("Non-performing")
            for j in range(len(epoch_tags[trial])):
                tags.append(epoch_tags[trial][j])
            try:
                dict["epochs." + trial + ".tags"] = tags
            except:
                sys.err("   Unable to create dataset 'tag' containing " + str(tags))
            # keep with tradition and create a units field, even if it's empty
            if trial not in epoch_units:
                units = ["NA"]
            else:
                units = epoch_units[trial]
            try:
                dict["epochs." + trial + ".units_present"] = units
            except:
                print "   Unable to create dataset 'units_present' containing ", units

            raw_path = "descrHash/value/%d" % (trial_id[i])
            raw_file = parse_h5_obj(data_h5[raw_path])[0]
            if len(raw_file) == 1:
                raw_file = 'na'
            else:
                raw_file = str(raw_file)

        elif options.data_type == "ophys":
            start = trial_t[i]
            stop  = trial_t[i+1]

            dict["epochs." + trial + ".start_time"] = start
            dict["epochs." + trial + ".stop_time"]  = stop

            if trial in epoch_roi_list.keys():
                dict["epochs." + trial + ".ROIs"]       = epoch_roi_list[trial]
                dict["epochs." + trial + ".ROI_planes"] = epoch_roi_planes[trial]
            tags = []
            if trial in epoch_trial_types:
                for j in range(len(epoch_trial_types[trial])):
                    tags.append(epoch_trial_types[trial][j])
            dict["epochs." + trial + ".tags"] = tags

    return dict


# ------------------------------------------------------------------------------

def extract_dict_master_shape(data_h5, plane_map, options, num_planes = 3):
    master_shape = {}
    num_subareas = len(get_key_list(data_h5['timeSeriesArrayHash'])) - 1
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if data_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(data_h5[plane_path].keys())
        area_grp = data_h5["timeSeriesArrayHash/descrHash"]["%d"%(2+subarea)]
        for plane in range(num_planes):
            if num_planes == 1:
                plane_grp = area_grp["value/1"]
            else:
                plane_grp = area_grp["value/1"]["%d"%(plane+1)]
            master = plane_grp["masterImage"]["masterImage"].value
            master1_shape = np.array(master).shape
            oname = "area%d_plane%d" % (subarea+1, plane+1)
            master_shape[oname] = master1_shape

    return master_shape

# ------------------------------------------------------------------------------

def extract_stimulus_subtrace(orig_trace, time, trial_t, types, types_to_remove, options):
    new_trace = orig_trace
    if options.verbose:
        print("    num_stimulus_types=" + str(len(types)) + " types=" + str(types))

    for i in range(len(types)):
        if str(types[i]) in types_to_remove:
            start = trial_t[i]
            stop  = trial_t[i+1]
            new_trace[(time >= start) & (time < stop)] = 0.
            if options.verbose:
                print("start=" + str(start) + " stop=" + str(stop) +  \
                      "  " + str(((time >= start) & (time < stop)).sum()) + \
                      " values has been reset")
    return new_trace

# ------------------------------------------------------------------------------


def stimulus_subtraces_by_type(data, laser_power, time, options):
    keyName = "PhotostimulationType"
    hash_group = data["/trialPropertiesHash"]
    types = get_value_pointer_by_key(hash_group, keyName, options.debug)
    if options.verbose:
        print("    types="+ str(types))

    trial_t  = data["trialStartTimes/trialStartTimes"].value
    # trial stop isn't stored. assume that it's twice the duration of other
    #   trials -- padding on the high side shouldn't matter
    ival     = (trial_t[-1] - trial_t[0]) / (len(trial_t) - 1)
    trial_t  = np.append(trial_t, trial_t[-1] + 2*ival)

    count_1   = types.count(1)
    count_2   = types.count(2)
    count_NaN = types.count('NaN')
    sub_traces = []
    num_types = 0
    trace_types = []
    if count_1 == 0 and count_2 > 0:
        num_types = 1
        trace_types = [2]
        sub_traces.append(extract_stimulus_subtrace(laser_power, time, trial_t, types, ['NaN'], options))
    elif count_2 == 0 and count_1 > 0:
        num_types = 1
        trace_types = [1]
        sub_traces.append(extract_stimulus_subtrace(laser_power, time, trial_t, types, ['NaN'], options))
    elif count_1 > 0 and count_2 > 0:
        num_types = 2
        trace_types = [1, 2]
        sub_traces.append(extract_stimulus_subtrace(laser_power, time, trial_t, types, ['2','NaN'], options))
        sub_traces.append(extract_stimulus_subtrace(laser_power, time, trial_t, types, ['1','NaN'], options))
    return (num_types, trace_types, sub_traces)

# ------------------------------------------------------------------------------


# pull out masterImage arrays and create entries for each in
#   /acquisition/images
# masterImages are store in:
#        tsah::descrHash::[2-7]::value::1::[1-3]::masterImage
# each image has 2 color channels, green and red

def create_reference_images(plane_map, num_subareas, data_h5, options):
    master_shapes   = {}
    ref_images_red   = {}
    ref_images_green = {}
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if data_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(data_h5[plane_path].keys())
        for plane in range(num_planes):
            area_grp = data_h5["timeSeriesArrayHash/descrHash"]["%d"%(2+subarea)]
            if num_planes == 1:
                plane_grp = area_grp["value/1"]
            else:
                plane_grp = area_grp["value/1"]["%d"%(plane+1)]
            master = plane_grp["masterImage"]["masterImage"].value
            master_shape = np.array(master).shape

            if len(master_shape) == 3 or not options.handle_errors:
                green = np.squeeze(master[:,:,0])
                red   = np.squeeze(master[:,:,1])
            else:
                green = master
                red   = master

            # convert from file-specific area/plane mapping to
            #   inter-session naming convention
            oname = "area%d_plane%d" % (subarea+1, plane+1)
            master_shapes[oname] = master_shape
            image_plane = plane_map[oname]
            ref_images_green[image_plane] = green
            ref_images_red[image_plane] = red
    return (master_shapes, ref_images_red, ref_images_green)

# ------------------------------------------------------------------------------

# each of nwb_object's hdf5 files have imaging planes and subareas
#   labels consistent within the file, but inconsistent between
#   files. create a map between the h5 plane name and the
#   identifier used between files
# plane_map = {}
def add_plane_map_entry(plane_map, h5_plane_name, filename, options):
    toks = filename.split("fov_")
    if len(toks) != 2:
        print("Error parsing %s for imaging plane name" % filename)
        sys.exit(1)
    univ_name = "fov_" + toks[1][:5]
    if univ_name not in plane_map:
        plane_map[h5_plane_name] = univ_name
    return univ_name

# ------------------------------------------------------------------------------

def create_plane_map(data_h5, options):
    plane_map = {}

    if options.handle_errors:
        num_subareas = 1
        if '1' in data_h5['timeSeriesArrayHash/descrHash'].keys():
            num_subareas = len(data_h5['timeSeriesArrayHash/descrHash'].keys()) - 1
    else:
        num_subareas = len(data_h5['timeSeriesArrayHash/descrHash'].keys()) - 1

    for subarea in range(num_subareas):
        # fetch time array
        if options.handle_errors:
            try:
                grp = data_h5['timeSeriesArrayHash/value/%d/imagingPlane' %(subarea + 2)]
                grp2 = data_h5['timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)]
            except:
                try:
                    grp = data_h5['timeSeriesArrayHash/value/imagingPlane']
                    grp2 = data_h5['timeSeriesArrayHash/descrHash/value']
                except:
                    print("Cannot create plane map")
                    break
        else:
            grp = data_h5['timeSeriesArrayHash/value/%d/imagingPlane' %(subarea + 2)]
            grp2 = data_h5['timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)]

        if grp2.keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(grp2.keys())
        for plane in range(num_planes):
            if options.handle_errors:
                try:
                    pgrp = grp["%d"%(plane+1)]
                except:
                    # Warning: only one imaging plane available (instead of 3)
                    pgrp = grp
            else:
                pgrp = grp["%d"%(plane+1)]

            old_name = "area%d_plane%d" % (subarea+1, plane+1)
            frame_idx = pgrp["sourceFileFrameIdx"]["sourceFileFrameIdx"].value
            lst = parse_h5_obj(pgrp["sourceFileList"])[0]
            for k in lst:
                # srcfile = str(lst[k][k].value)
                srcfile = str(k)
                add_plane_map_entry(plane_map, old_name, srcfile, options)
                break
    return plane_map

# ------------------------------------------------------------------------------

def save_2p_frames(external_file, starting_frame, timestamps, fname, stack_t):
    starting_frame.append(len(timestamps))
    timestamps.extend(stack_t)
    external_file.append(fname.encode('utf8'))

# ------------------------------------------------------------------------------

def extract_frame_data(frame_idx, t, subarea, plane, plane_map, srcfile, cnt, options):
    lastfile =-1
    lastframe = 1
    stack_t = []
    nname = None
    fname = None
    zero = np.zeros(1)
    assert len(t) == len(frame_idx[0])
    # following arrays used to make external_file as an array, reducing number of image_series
    external_file = []
    starting_frame = []
    timestamps = []
    for i in range(len(frame_idx[0])):
        filenum = frame_idx[0][i]
        if lastfile < 0:
            lastfile = filenum
        framenum = frame_idx[1][i]
        stack_t.append(t[i])
        # check for embedded NaNs
        if np.isnan(framenum):
            continue
        # use fname as a flag. if it's not None then there's data
        #   to write
        if fname is None:
            # convert from file-specific area/plane mapping to
            #   inter-session naming convention
            oname = "area%d_plane%d" % (subarea+1, plane+1)
            nname = plane_map[oname]
            name = "%s_%d" % (nname, filenum)
            fname = os.path.basename(srcfile["%d"%filenum])
        # make sure frames and file numbers are sequential
        if not (lastfile == filenum and framenum == lastframe+1) and \
           not                          framenum == 1:
            # Warning: framenum or filenum does not start from 1
            continue
#       assert (lastfile == filenum and framenum == lastframe+1) or framenum == 1
        if lastfile != filenum:
            if i>0:
                if not np.isnan(frame_idx[0][i-1] ) and \
                   not np.isnan(frame_idx[1][i-1]):
                    if options.handle_errors:
                        try:
                             save_2p_frames(external_file, starting_frame, \
                                            timestamps, fname, stack_t)
                        except:
                            print("Warning: unable to create_2p_ts for name=" + str(name))
                    else:
                        save_2p_frames(external_file, starting_frame, \
                                       timestamps, fname, stack_t)
                    stack_t = []
                    fname = None
        lastframe = framenum
        lastfile = filenum
    # make sure we write out the last entry
    if fname is not None:
        if options.handle_errors:
            try:
                save_2p_frames(external_file, starting_frame, \
                               timestamps, fname, stack_t)
            except:
                print("Warning: unable to create_2p_ts for name=" + str(name))
        else:
            save_2p_frames(external_file, starting_frame, \
                           timestamps, fname, stack_t)
    return (external_file, starting_frame, timestamps, nname)

# ------------------------------------------------------------------------------

def get_child_group_names(parent_group_pointer):
    return parent_group_pointer.keys()

# ------------------------------------------------------------------------------

def get_group_keys(group_pointer):
    key_list = []
    for k in group_pointer.keys():
        if not re.search('#', k):
            key_list.append(k)
    return sorted(key_list)

# ------------------------------------------------------------------------------

def get_data_from_refs_dataset(orig_h5, dataset_pointer):
    path = dataset_pointer.name
    if dataset_pointer.dtype.kind == 'O':
        ref = orig_h5[path][0][0]
        obj = orig_h5[ref]
        value = ''.join(i for i in obj[:])
    else:
        obj = orig_h5[path]
        value = ''.join(i.encode(ascii) for i in obj[:])
    return value

# ------------------------------------------------------------------------------

def get_key_list(hash_group_pointer):
    try:
        keyNames_dataset = hash_group_pointer['keyNames/keyNames']
        key_list = np.array(keyNames_dataset).tolist()
    except:
        key_list = []
    return key_list                                         

# ------------------------------------------------------------------------------

def get_value_by_key(hash_group_pointer, key):
    if verbose:
        print "In get_value_by_key: group_pointer.name=", hash_group_pointer.name
        print "                                    key=", key
    key_list = get_key_list(hash_group_pointer)
    if verbose:
        print "In get_value_by_key: partial_key=", key, " key_list=", key_list
    key_name = ""
    for k in key_list:
        if re.search(key, k):
            key_name = k
    if verbose:
        print "In get_value_by_key: key_name=", key_name, "\n"
    ind = key_list.index(key_name)
    if verbose:
        print "ind=", ind, "\n"
    value_group = hash_group_pointer['value']
    if verbose:
        print "\nIn get_value_by_key: value_group.name=",  value_group.name
        print "\n                   value_group.keys()=", value_group.keys()
        print "\nlen(value_group.keys())=", len(value_group.keys()), " ind=", ind  
    if len(value_group.keys()) == 1 and "value" in value_group.keys():
        value_pointer = value_group['value']
        if item_type(value_pointer) == "dataset":
            if verbose:
                print "Case1: returning an element of a dataset"
            value = value_pointer[ind]
        else:
            if verbose:
                print "Case2: value/value type is ", item_type(value_pointer)
            sys.exit("Unsapported case in get_value_by_key()")
    else:
        if verbose:
            print "level2 key_list=", value_group.keys(), ";  ind+1=", ind+1
        if str(ind+1) in value_group.keys():
            value_pointer = value_group[str(ind+1)]
            if item_type(value_pointer) == "dataset":
                if verbose:
                    print "Case3: value/" + str(ind+1) + " is a dataset"
                value = np.array(value_pointer).tolist()
            else:
                if str(ind+1) in value_pointer.keys() and \
                   item_type(value_pointer[str(ind+1)]) == "dataset":
                    if verbose:
                        print "Case4: value/" + str(ind+1) + "/" + str(ind+1) + " is ", item_type(value_pointer)
                    value = np.array(value_pointer[str(ind+1)]).tolist()
                else:
                    if verbose:
                        print "Case5: value/" + str(ind+1) + " is ", item_type(value_pointer)
                    value =  value_pointer
        else:
            value = value_group
    return value

# ------------------------------------------------------------------------------

def get_value2_by_key2(hash_group1_pointer, key1, hash_group2, key2):
    if verbose:
        print "\n\nEntering get_value2_by_key2 ..."
    key_list = get_key_list(hash_group1_pointer)
    key_name = ""
    for k in key_list:
        if re.search(key1, k):
            key_name = k
    if verbose:
        print "In get_value_by_key: key_name=", key_name, "\n"
    ind = key_list.index(key_name)
    if verbose:
        print "ind=", ind, "\n"
    value_group = hash_group1_pointer['value']
    if verbose:
        print "\nIn get_value_by_key: value_group.name=",  value_group.name
        print "\n                   value_group.keys()=", value_group.keys()
        print "\nlen(value_group.keys())=", len(value_group.keys()), " ind=", ind
        print "In get_value2_by_key2: level2 key_list=", value_group.keys()
    if str(ind+1) in value_group.keys():
        value_pointer = value_group[str(ind+1)]
    else:
        value_pointer = value_group
    if verbose:
        print "value_pointer.keys()=", value_pointer.keys(), "; hash_group2=", hash_group2
    if 'valueMatrix' in value_pointer.keys() and \
       'idStr'       in value_pointer.keys():
       # ephys data
       if verbose:
           print " ephys data"
       key_list2 = np.array(value_pointer['idStr/idStr']).tolist()
       ind2 = key_list2.index(key2)
       if verbose:
           print "ind2=", ind2, "\n"
       value2 = value_pointer['valueMatrix/valueMatrix'][:,ind2]
    else: 
        try:
            value_pointer2 = value_pointer[hash_group2]
            value2 = get_value_by_key(value_pointer2, key2)
        except:
            sys.exit("\nCannot determine value2 in get_value2_by_key2")
    return value2

# ------------------------------------------------------------------------------

def get_value_pointer_by_path_items(orig_h5, path_items):  
    path = orig_h5.name
    if verbose:
        print "Enter get_value_pointer_by_path_items: orig_h5.name=", orig_h5.name, " path_items=", path_items
    if len(path_items) > 0 and not path_items == ['']:
        iter_len = len(path_items)
        i = 0
        while i < iter_len:
            if path_items[i] in orig_h5[path].keys():
                if verbose:
                    print "Case 1, item=", path_items[i]
                path += '/' + str(path_items[i])
                if verbose:
                    print "Case1, path=", path
                value_pointer = orig_h5[path]
                if verbose:
                    print "item_type(value_pointer)=", item_type(value_pointer)
                if item_type(value_pointer) == "dataset":
                        if verbose:
                            print "dataset_shape=", value_pointer.shape
                        if value_pointer.shape[0] > 1:
                            value = np.transpose(value_pointer)
                        else:
                            value = np.array(value_pointer)
                        if verbose:
                            print "data=", value, " data_len=", len(value), " dataset_shape=", value.shape, \
                                  " i=", i, " iter_len=", iter_len
                        if i == iter_len - 2:
                            return value[int(path_items[i+1])]
                        elif i == iter_len - 3:
                            return value[int(path_items[i+1])][int(path_items[i+2])]
                        return value
            else:
                iter_len = len(path_items) + 1
                value = []
                if verbose:
                    print "Case2, i=", i, " iter_len=", iter_len
                for k in orig_h5[path].keys():
                    path2 = path + '/' + str(k)
                    value_pointer = orig_h5[path2]
                    if verbose:
                        print "item_type(value_pointer)=", item_type(value_pointer)
                    if verbose and item_type(orig_h5[path2]) == "dataset":
                        print "dataset_shape=", orig_h5[path2].shape
                    elif item_type(orig_h5[path2]) == "group":
                        for k1 in orig_h5[path2].keys():
                            path3 = path2 + '/' + str(k1)
                            if verbose:
                                print "    \npath2=", path2, " path3=", path3
                            value_pointer = orig_h5[path3]
                            if verbose:
                                print "item_type(value_pointer1)=", item_type(value_pointer)
                            if item_type(orig_h5[path3]) == "dataset":
                                if verbose:
                                    print "dataset_shape=", orig_h5[path3].shape
                                if value_pointer.shape[0] > 1:
                                    value1 = np.transpose(np.array(value_pointer)).tolist()
                                else:
                                    value1 = np.array(value_pointer).tolist()
                                    if verbose:
                                        print "data=", value
                                value.append(value1)
                    if verbose:
                        print "Case2, path=", path, " value=", value
                return value
            i += 1
        if verbose:
            print "path=", path
#       print "value_pointer.name=", value_pointer.name
    else:
        value_pointer = orig_h5
    if verbose:
        print "    In get_value_pointer_by_path_items: path=", path
        if item_type(value_pointer) == "group":
            print "                                        key_list=", value_pointer.keys()
    return value_pointer

# ------------------------------------------------------------------------------

def get_key_index(key_list, partial_key):
    key_name = ""
    for k in key_list:
        if re.search(partial_key, k):
            key_name = k
    ind = key_list.index(key_name)
    return ind

# ------------------------------------------------------------------------------
# Returns a pointer to the value 
def get_value_pointer_by_key(hash_group_pointer, partial_key, verbose):
    key_list = get_key_list(hash_group_pointer)            
    if len(key_list) > 1: 
        ind = get_key_index(key_list, partial_key) + 1
    else:
        ind = 1
    value_group_pointer = hash_group_pointer['value']
    if str(ind) in hash_group_pointer['value/'].keys():
        value_group_pointer = hash_group_pointer['value/' + str(ind)]
    if verbose:
        print "   In get_value_pointer_by_key: value_group_items=", value_group_pointer.keys(), " ind=", ind
        print "                                item_type(value_group_pointer)=", item_type(value_group_pointer)
    if str(ind) in value_group_pointer.keys():
        if verbose:
            print "    In get_value_pointer_by_key: case 1"
        value_pointer1 = value_group_pointer[str(ind)]     
        if item_type(value_pointer1) == "group" and str(ind) in value_pointer1.keys():
            if verbose:
                print "    In get_value_pointer_by_key: case 11"
            value_pointer = value_pointer1[str(ind)]
        else:
            if verbose:
                print "    In get_value_pointer_by_key: case 12"
                print "                                 item_type(value_pointer1)=", item_type(value_pointer1)
                print "                                 value_pointer1.name=", value_pointer1.name
            value_pointer = value_pointer1 
    else:
        if verbose:
            print "    In get_value_pointer_by_key: case 2"
        value_pointer = value_group_pointer
    if item_type(value_pointer) == "dataset":
        value_pointer = np.array(value_pointer).tolist()              
    return value_pointer  

# ------------------------------------------------------------------------------

def get_value_by_path_items(orig_h5, path_items):
    value_pointer = get_value_pointer_by_path_items(orig_h5, path_items)
    if verbose:
        print "    In get_value_by_path_items: value_pointer.name=", value_pointer.name
    if hasattr(value_pointer , '__dict__') and len(value_pointer.keys()) > 1:
        if verbose:
            print "    Case 1"
        data = []
        for k in value_pointer.keys():
            data1 = np.array(value_pointer[k + "/" + k]).tolist()
            data.append(data1)
    else:
        if verbose:
            print "    Case 2"
        data = np.array(value_pointer[path_items[len(path_items)-1]]).tolist()
    return data

# ------------------------------------------------------------------------------

def get_all_keys(orig_h5, metadata):
    all_keys = []
    # Extract data keys
    top_groups = get_child_group_names(orig_h5)
    for group in top_groups:
        if not re.search("Hash", group):
            continue
        group_keys = orig_h5[group].keys()
        if len(group_keys) < 2:
            continue
        path_items = [group, "keyNames", "keyNames"]
        print "path_items=", path_items
        key_list = get_value_pointer_by_path_items(orig_h5, path_items)
        for k in key_list:
            if not k in all_keys:
                all_keys.append(k)
    # Extract metadata keys
    if len(metadata) > 0:
        top_groups = get_child_group_names(metadata)
        for k in top_groups:
            if not k in key_list:
                all_keys.append(k)
    return all_keys

# ------------------------------------------------------------------------------

def get_description_by_key(hash_group_pointer, partial_key):
    key_list = get_key_list(hash_group_pointer)
    ind = get_key_index(key_list, partial_key) + 1
    descr_data = np.array(hash_group_pointer['descr/descr']).tolist()
    if verbose:
        print "\nIn get_description_by_key: len(descr_data)=", len(descr_data), \
              "\n     descr_data_items=", descr_data, " ind=", ind
        print "\nFunction get_description_by_key returns: ", descr_data[ind-1]
    return descr_data[ind-1]

# ------------------------------------------------------------------------------

def get_file_list(data, match_string):
    # Compile a file list
    file_list = []
    if os.path.isfile(data) and (re.search(".h5", data) or re.search(".nwb", data) or re.search(".borg", data)):
        file_list = [ data ]
    elif os.path.isdir(data):
        file_list = []
        if verbose:
            print "num_files=", len(os.listdir(data))
        for file in os.listdir(data):
            if len(match_string) > 0 and not re.search(match_string, file):
                continue
            if (re.search(".h5", file) or re.search(".nwb", data)) \
                and os.path.isfile(os.path.join(data, file)):
                file_path = os.path.join(data, file)
                file_list.append(file_path)
    else:
        sys.exit("Cannot process " + data)
    return file_list
   
# ------------------------------------------------------------------------------

def item_type(item_pointer):
    item_type = ""
    try:
        keys = item_pointer.keys()
        item_type = "group"
    except:
        try:
            keys = np.array(item_pointer).tolist() 
            item_type = "dataset"
        except:
            item_type = "data_item"
    return item_type


# ------------------------------------------------------------------------------

def get_data_by_key(nwb_root, item_path, verbose):
        item_pointer = h5_root[item_path]
        try:
            # item is group
            keys = item_pointer.keys()
            if verbose:
#               print "group: path=", item_path, " members=" + str(keys)
                print "group: path=", item_path, " num_members=", len(keys)
            for k in keys:
                if len(item_path) == 1:                # item_path == '/'
                    item_path1 = item_path + k
                else:
                    item_path1 = item_path + '/' + k
                parse_h5_item(h5_root, item_path1, verbose)
        except:
            # item is dataset
            try:
                data = np.array(item_pointer)
                if verbose:
                    print "dataset: path=", item_path, " , shape=", data.shape, \
                          " , dtype=", data.dtype, " data=", data.tolist()
            except:
                sys.exit("Path " + path + " is not valid")

