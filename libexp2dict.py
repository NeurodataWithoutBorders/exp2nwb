#!/usr/local/python-2.7.11/bin/python
#
# libexp2dict.py
#
# Extracting data dictionaries from experimental data structures

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
import libexp2nwb 

# ------------------------------------------------------------------------------
def check_entry(file_name,obj):
    try:
        return file_name[obj]
    except KeyError:
        print(str(obj) +" does not exist")
        return []

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
    child_groups = libh5.get_child_group_names(input_h5)
    if options.data_origin == "SP":
        d = libh5.get_value_by_key(input_h5['/metaDataHash'], \
                         "dateOfExperiment")
        t = libh5.get_value_by_key(input_h5['/metaDataHash'], \
                         "timeOfExperiment")
        dt=datetime.datetime.strptime(d+t, "%Y%m%d%H%M%S")

    elif options.data_origin in ["NL", "DG"]:
        d = np.array(libh5.get_value_pointer_by_path_items(input_h5, \
                     ["dateOfExperiment", "dateOfExperiment"])).tolist()[0]
        t = np.array(libh5.get_value_pointer_by_path_items(input_h5, \
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
    elif options.data_origin == "JY":
        d = np.array(libh5.get_value_pointer_by_path_items(input_h5, \
                         ["dateOfExperiment", "dateOfExperiment"])).tolist()[0]
        print("date=", "20"+d)
        dt = datetime.datetime.strptime("20"+d, "%Y%m%d")
    else:
        sys.exit("Data origin is unknown")
       
    return dt.strftime("%a %b %d %Y %H:%M:%S")

# ------------------------------------------------------------------------------

def set_metadata(group, keyname, value):
    if keyname in ["extracellular", "intracellular"]:
        group.set_dataset("description", value)
    else:
        group.set_custom_dataset(keyname, value)   

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

def top_datasets(data_h5, meta_h5, options):
    vargs = {}
    session_id = options.project_dir
    if options.data_origin in ["NL", "JY", "DG"]:
        vargs["start_time"] = find_exp_time(meta_h5, options)
        vargs["description"]= "Extracellular ephys recording of mouse doing discrimination " + \
                              "task (lick left/right), with optogenetic stimulation " +\
                              "plus pole and auditory stimulus"
    else:
        vargs["start_time"] = find_exp_time(data_h5, options)
        vargs["description"]= " "
        print("Warning: missing experiment description (to be added)")
    vargs["file_name"]      = os.path.join(options.project_dir, "top_datasets.h5")
    vargs["identifier"]     = nwb_utils.create_identifier(options.project_dir)
    vargs["mode"] = "w"

    return vargs
 
# ------------------------------------------------------------------------------

def general(input_h5, options):
    if options.verbose:
        print("Extracting metadata to a dictionary")
    dict_metadata = {}

    genotype = ""
    subject  = ""
    weight   = ""
    animalStrain = ""
    animalSource = ""
    animalGeneModification = ""
    value = ""
    DOB   = ""
    dict_metadata["institution"] = "Janelia Research Campus"
    dict_metadata["lab"]         =  "Svoboda lab"

    # Add citation info
    if options.data_origin in ["NL", "JY", "DG"]:
        key_list = libh5.get_child_group_names(input_h5)
        group_pointer = input_h5
    elif options.data_origin == "SP":
        key_list = libh5.get_key_list(input_h5["metaDataHash"])
        group_pointer = input_h5['/metaDataHash']

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
                key1_list = libh5.get_value_pointer_by_path_items(group_pointer, [key]).keys()
                if options.verbose:
                    print("   key=" + key+ " key1_list="+ str(key1_list))
                if re.search("tracellular", key) and not "impedance" in key1_list:
                    dict[key + ".impedance"] = "not recorded"
                for key1 in key1_list:
                    if options.verbose:
                        print("      key1="+ key1)
                    if key1 in ["siteLocations"]:
                        continue
                    elif key1 in ["ADunit", "penetrationN", "groundCoordinates", \
                                  "extracellularDataType", "recordingMarker", \
                                  "recordingType", "impedance"]:
                        value1 = libh5.get_value_pointer_by_path_items(group_pointer, [key, key1, key1])
                        if key1 == "ADunit":
                            dict[key + ".ADunit"] = str(value1[0])
                        elif key1 == "penetrationN":
                            dict[key + ".penetration_num"] = str(value1[0])
                        elif key1 == "groundCoordinates":
                            dict[key + ".ground_coordinates"]= value1
                        elif key1 == "extracellularDataType":
                            dict[key + ".data_types"]= value1
                        elif key1 == "recordingMarker":
                            dict[key + ".recording_marker"]= str(value1[0])
                        elif key1 == "recordingType":
                            dict[key + ".recording_type"]= str(value1[0])
                        elif key1 == "impedance":
                            dict[key + ".impedance"]= str(value1[0])
                    value1 = libh5.get_value_pointer_by_path_items(group_pointer, [key, key1, key1])
                    if options.verbose:
                        print("      key="+ key+ " key1="+ key1+ " value1="+ str(value1))
                    value2 = [str(v) for v in value1]
                    add_value = "       " + key1 + ": " + ",".join(value2) + "\n "
                    value += add_value
                    if options.verbose:
                        print("         key=", key, " key1=", key1, " add_value=", add_value, " now value=\n", value)
            except:
                if key == "virus" and options.data_origin == "SP":
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
            key1_list = libh5.get_value_pointer_by_path_items(group_pointer, [key]).keys()
            value1 = {}
            if options.verbose:
                print("key= photostim, key1_list=" + str(key1_list))
            num_sites = 1      
            photostim_loc1 = '' 
            for key1 in key1_list:
                try:
                    value1[key1] = libh5.get_value_by_path_items(group_pointer, [key, key1])
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
                if options.data_origin == "NL":
                    dict_metadata["general.optogenetics.site_" + str(s+1) + ".device"] = "Ciel 250 mW 473nm Laser from Laser Quantum"
                for key1 in key1_list:
                    if key1 == "stimulationMethod":
                        dict_metadata["general.optogenetics.site_" + str(s+1) + ".stimulation_method"]=str(value1[key1][s])
                    elif key1 == "photostimCoordinates":
                        coord = str(value1[key1][s]) + " in mm"
                        if "photostimLocation" in key1_list:
                            coord += "\natlas location: " + str(value1["photostimLocation"][s]) 
                        dict_metadata["general.optogenetics.site_" + str(s+1) + ".location" ] = coord
                    elif key1 == "photostimWavelength":
                        dict_metadata["general.optogenetics.site_" + str(s+1) + ".excitation_lambda"] = str(value1[key1][0])
                    
        else:
            if "metaDataHash" in libh5.get_child_group_names(input_h5) \
                and not key in ["citation"]:
                value = libh5.get_value_by_key(group_pointer,key)
                if libh5.item_type(value) == "dataset":
                    value = np.array(value).tolist()
                elif libh5.item_type(value) == "group":
                    value = value.name
            elif not "metaDataHash" in libh5.get_child_group_names(input_h5):                          
                value_list = np.array(libh5.get_value_pointer_by_path_items(group_pointer, [key, key])).tolist()
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
            dict_metadata["subject_id"] = value
 
        elif key == "dateOfBirth":
            subject  = key + ": " + str(value) + "\n"
            DOB = str(value)

        elif key == "dateOfExperiment":
            age = ""
            if options.verbose:
                print "\noptions.data_origin=", options.data_origin
            if options.data_origin == "NL":
                DOE = str(value)
                age = ""
                if len(DOB) ==8 and len(DOE)== 8:                   
                    age = compute_age(DOB, DOE)
                else:
                    age = "3 to 5 months"
            elif options.data_origin == "SP":
                age = "6 to 8 weeks"
            if len(age) > 0:
                dict_metadata["age"]= age

        elif re.search("citation", key):
            # Add citation info
            if options.data_origin in ["NL", "SP"]:
                if options.data_origin == "NL":
                    value = "doi: 10.1038/nature14178"
                elif options.data_origin == "SP":
                    value = "doi: 10.1016/j.neuron.2015.03.027"
                dict_metadata["related_publications"]=value

        elif re.search("experimentType", key):
            dict_metadata["notes"]=value

        elif re.search("experimenters", key):
            dict_metadata["experimenter"]=value

        elif key in ["sex", "species", "age", "cell"]:
            dict_metadata[key]=value

        elif re.search("weight", key) and len(str(weight)) == 0 and len(str(value)) > 0:
            weight = str(value)
            dict_metadata["weight"]=weight

        elif re.search("referenceAtlas", key):
            dict_metadata["reference_atlas"]=value

        elif re.search("whiskerConfig", key):
            dict_metadata["whisker_configuration"] = value

        elif key in ["virus", "fiber"]:
            dict_metadata[key]=value

        elif key == "behavior":
            task_kw = map(str,parse_h5_obj(input_h5["behavior/task_keyword"])[0])
            dict_metadata["task_keywords"]=task_kw

#       elif key in ["extracellular", "intracellular"]:
#           ephys_group.set_custom_dataset(key, value)
    dict_metadata["subjectgenotype"]    = genotype + "\n"
    dict_metadata["subjectdescription"] =  animalStrain + "\n" + animalSource + "  " + subject
    source_script="https://github.com/NeurodataWithoutBorders/exp2nwb"
    dict_metadata["source_script"] = source_script
    if options.data_origin == "NL":
        dict_metadata["surgery"]                = metadata_from_file("surgery.txt", options)
        dict_metadata["data_collection"]        = metadata_from_file("data_collection.txt", options)
        dict_metadata["experiment_description"] = metadata_from_file("experiment_description.txt", options)
    elif options.data_origin == "SP":
        dict_metadata["data_collection"]       ="doi: 10.1016/j.neuron.2015.03.027"
        dict_metadata["experiment_description"]="doi: 10.1016/j.neuron.2015.03.027"
        dict_metadata["surgery"]               ="doi: 10.1016/j.neuron.2015.03.027"
        dict_metadata["whisker_configuration"] = "see Table S1 in doi: 10.1016/j.neuron.2015.03.027"
        dict_metadata["weight"]                = "not recorded"
    elif options.data_origin == "DG":
        dict_metadata["data_collection"]        = "doi:10.1038/nn.4412"
        dict_metadata["experiment_description"] = "doi:10.1038/nn.4412"
        dict_metadata["surgery"]                = "doi:10.1038/nn.4412"

    else:
        print("Missing  data_collection.txt, experiment_description.txt and surgery.txt")
    
    return dict_metadata

# ------------------------------------------------------------------------------

def read_probe_locations_matrix(meta_h5, options):
    num_locs = len(libh5.get_value_pointer_by_path_items(meta_h5, \
                   ["extracellular", "siteLocations"]).keys())
    M = np.zeros([num_locs, 3])
    for i in range(num_locs):
        probe_id = i + 1
        coords = np.array(libh5.get_value_pointer_by_path_items(meta_h5, \
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

def get_description(meta_h5, options):
    extra_grp = libh5.get_value_pointer_by_path_items(meta_h5, \
                                ["extracellular"])
    value = ""
    for key in extra_grp.keys():
        if key in ["recordingType", "penetrationN", "groundCoordinates","referenceCoordinates", \
                   "extracellularDataType", "cellType", "identificationMethod","spikeSorting"]:
            try:
                # getting a dictionary
                key1_list = libh5.get_value_pointer_by_path_items(extra_grp, [key, key]).keys()
                for key1 in key1_list:
                    value1 = libh5.get_value_pointer_by_path_items(extra_grp, [key, key1])[:]
                    value2 = [str(v) for v in value1]
                    value += "       " + key1 + ": " + ",".join(value2) + "\n "
            except:
                # getting a dataset
                value = libh5.get_value_pointer_by_path_items(extra_grp, [key, key])
    return value

# ------------------------------------------------------------------------------

def get_device(meta_h5, options):
    extra_grp = libh5.get_value_pointer_by_path_items(meta_h5, \
                                ["extracellular"])
    value = ""
    for key in extra_grp.keys():
        if key in ["probeSource", "probeType", "ADunit", "amplifierRolloff"]:
            key1_list = libh5.get_value_pointer_by_path_items(extra_grp, [key]).keys()
            for key1 in key1_list:
                value1 = libh5.get_value_pointer_by_path_items(extra_grp, [key, key1])[:]
                value2 = [str(v) for v in value1]
                value += "       " + key1 + ": " + ",".join(value2) + "\n "
    return value

# ------------------------------------------------------------------------------

def general_devices_data():
    dict_devices = {}
    dict_devices["ephys_acquisition"] = "32-electrode NeuroNexus silicon probes recorded on a PCI6133 National Instrimunts board. See 'general/experiment_description' for more information"
    dict_devices["photostim_source"] = "Stimulating laser at 473 nm"
    return dict_devices

# ------------------------------------------------------------------------------

def general_extracellular_ephys_data(meta_h5, options):
    dict = {}
    M = read_probe_locations_matrix(meta_h5, options)

    # probe = M.tolist()
    probe = []
    sites = parse_h5_obj(check_entry(meta_h5, "extracellular/siteLocations"))
#   assert len(sites) == 32, "Expected 32 electrode locations, found %d"%len(sites)
    for i in range(len(sites)):
        probe.append(sites[i])
        probe[-1] = probe[-1] * 1.0e-3
    probe = np.asarray(probe)

    num_shanks, shank_size, shank_coords = detect_shanks(M)
    shank = []
    for i in range(1, (num_shanks+1)):
        for j in range(shank_size):
            shank.append("shank" + str(i))
    dict["general.extracellular_ephys.electrode_map"] = probe
    dict["general.extracellular_ephys.electrode_group"] = shank
    dict["general.extracellular_ephys.filtering"]  = "Bandpass filtered 300-6K Hz"

    ephys_device_txt = "32-electrode NeuroNexus silicon probes recorded on a PCI6133 National Instrimunts board. See 'general/experiment_description' for more information"
    dict["general.devices.ephys_acquisition"] = ephys_device_txt
    dict["general.devices.photostim_source"]  = "Stimulating laser at 473 nm"

    # Creating the shank groups 
    probe_type = libh5.get_value_pointer_by_path_items(meta_h5, \
                     ["extracellular", "probeType", "probeType"])
    rloc = libh5.get_value_pointer_by_path_items(meta_h5, \
               ["extracellular", "recordingLocation", "recordingLocation"])
    description = get_description(meta_h5, options)
    device      = get_device(     meta_h5, options)
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

def general_optophysiology_data(data_h5, options):
    dict = {}
    plane_map = create_plane_map(data_h5, options)
    num_subareas = len(libh5.get_key_list(data_h5['timeSeriesArrayHash'])) - 1
    master_shape = extract_dict_master_shape(data_h5, plane_map, options)

    dict["general.optophysiology.description"] = \
         "Imaging data were acquired using a custom two-photon microscope; " + \
         "see Peron et al., 2015 (DOI: 10.1016/j.neuron.2015.03.027) for details."
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if data_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(data_h5[plane_path].keys())
        for plane in range(num_planes):
            oname = "area%d_plane%d" % (subarea+1, plane+1)
            manifold = define_manifold(master_shape[oname], plane_map[oname])
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

def get_valid_trials(data_h5, data, options):
    ts_path = "timeSeriesArrayHash/descrHash/"
    val = []
    if options.handle_errors:
        num_subareas = 1
        if '1' in data_h5['timeSeriesArrayHash/descrHash'].keys():
            num_subareas = len(data_h5['timeSeriesArrayHash/descrHash'].keys()) - 1
    else:
        num_subareas = len(data_h5['timeSeriesArrayHash/descrHash'].keys()) - 1

    if data == "whisker":
        if options.handle_errors:
            try:
                ids = parse_h5_obj(data_h5[ts_path + '1/value'])[0]
            except:
                ids = parse_h5_obj(data_h5[ts_path + 'value'])[0]
        else:
            ids = parse_h5_obj(data_h5[ts_path + '1/value'])[0]
        val = val + list(ids)
        val = list(Set(val))
    if data == "Ca":
        for i in range(2,num_subareas+1):
            ids_path = ts_path + "%d/value/2" % i
            ids = parse_h5_obj(data_h5[ids_path])[0]
            # ids = list(data_h5[ids_path].value)
            val = val + list(ids)
        val = list(Set(val))
    return val

# ------------------------------------------------------------------------------

def collect_analysis_information(data_h5, options):
    dict = {}
    if options.verbose:
        print("Collecting analysis information")

    trial_start_times = data_h5["trialStartTimes/trialStartTimes"].value
    trial_types_all = []
    trial_type_strings = parse_h5_obj(data_h5['trialTypeStr'])[0]
    # collect all trials (strings)
    for i in range(len(trial_type_strings)):
            trial_types_all.append(str(trial_type_strings[i]))
    trial_type_mat = data_h5['trialTypeMat/trialTypeMat'].value

    if options.data_origin == "NL":
        good_trials = data_h5['trialPropertiesHash/value/4/4'].value
    elif options.data_origin == "SP":
        good_trials = [0] * len(trial_start_times)
        valid_whisker = get_valid_trials(data_h5, "whisker", options)
        valid_Ca = get_valid_trials(data_h5, "Ca", options)
        for i in range(len(good_trials)):
            j = i + 1
            if j in valid_whisker and j in valid_Ca:
                good_trials[i] = 1
    elif options.data_origin == "SP":
        good_trials = data_h5['trialPropertiesHash/value/6/6'].value
    elif options.data_origin == "DG":
        good_trials = data_h5['trialPropertiesHash/value/5/5'].value
    dict["trial_start_times"] = trial_start_times
    dict["trial_type_string"] = trial_types_all
    dict["trial_type_mat"]    = trial_type_mat
    dict["good_trials"]       = good_trials
    if not options.data_origin == "DG":
        dict["description"]="Trial_type_string has six values, with three response possibilities and two correct responses.  The format is {response type}{correct response}, with response type taking the value Hit, Err, or NoLick, indicating that the animal got the trial right, wrong, or did nothing, respectively.  Corret response is either L or R, indicating that the animal should have licked left or right, respectively."
    else:
        dict["description"]="Trial_type_string has six values, with three response possibilities and two correct responses.  The format is {response type}{correct response}, with response type taking the value Hit, Err, or NoLick, indicating that the animal got the trial right, wrong, or did nothing, respectively.  Corret response is either L or R, indicating that the animal should have licked left or right, respectively."

    return dict

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

def create_epochs(data_h5, options):
    dict = {}
    
    if options.verbose:
        print("    Getting trial types ...")
    get_trial_types(data_h5, options)

    if options.data_origin in ["NL", "DG"]:
        keys = data_h5['eventSeriesHash/value'].keys()
        if options.verbose:
            print("    Getting trial units ...")
        keys = data_h5['eventSeriesHash/value'].keys()
        if options.data_origin == "NL":
            unit_num = len(keys)
        else:
            unit_num = 1
        get_trial_units(data_h5, unit_num, options)

    if options.verbose:
        print("    Creating trials ...")
    dict = create_trials(dict, data_h5, options)

    return dict

# ------------------------------------------------------------------------------

# add trial types to epoch for indexing
epoch_tags = {}
epoch_trial_types = {}
def get_trial_types(data_h5, options):
    trial_id = data_h5["trialIds/trialIds"].value
    trial_types_all = []
    trial_type_strings = parse_h5_obj(data_h5['trialTypeStr'])[0]
    num_trial_types    = len(trial_type_strings)
    if options.data_origin == "NL":
        photostim_types = libh5.get_value_by_key(data_h5["trialPropertiesHash"], "PhotostimulationType")
    elif options.data_origin == "SP":
        # SP data
        valid_whisker = get_valid_trials(data_h5, "whisker", options)
        valid_Ca      = get_valid_trials(data_h5, "Ca", options)
#       print "valid_whisker=", valid_whisker
#       print "valid_Ca=", valid_Ca
    elif options.data_origin == "JY":
        # JY data
        trial_types = libh5.get_value_by_key(data_h5["trialPropertiesHash"], 'TrialName')

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
        if options.data_origin == "NL":
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
        elif options.data_origin == "SP":
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
#           valid_whisker = get_valid_trials(orig_h5, "whisker", options)
#           valid_Ca = get_valid_trials(orig_h5, "Ca", options)
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

    if options.data_origin == "NL":
        good_trials = libh5.get_value_by_key(data_h5["/trialPropertiesHash"],  "GoodTrials")
        ignore_ivals_start = [time for (time, good_trial) in zip(trial_t,good_trials) if good_trial == 0]
        # trial stop isn't stored. assume that it's twice the duration of other
        #   trials -- padding on the high side shouldn't matter
        ignore_ivals_stop = [time for (time, good_trial) in zip(trial_t[1:],good_trials) if good_trial == 0]
        ignore_intervals  = [ignore_ivals_start, ignore_ivals_stop]

        if options.data_origin == "NL":
            keyName3 = "PhotostimulationType"
            hash_group_pointer2 = data_h5["/trialPropertiesHash"]
            stimulus_types = np.array(libh5.get_value_by_key(hash_group_pointer2, keyName3)).tolist()
            count_1 = stimulus_types.count(1)
            count_2 = stimulus_types.count(2)
    elif options.data_origin == "SP":
        plane_map = create_plane_map(data_h5, options)
        epoch_roi_list,epoch_roi_planes = create_trial_roi_map(data_h5, plane_map, options)
    elif options.data_origin == "DG":
        good_trials = data_h5['trialPropertiesHash/value/5/5'].value
    elif options.data_origin in ["JY"]:
        ephys_value_pointer = libh5.get_value_by_key(data_h5["/timeSeriesArrayHash"], "Ephys")
        time = np.array(libh5.get_value_pointer_by_path_items(ephys_value_pointer, ["time", "time"]))
    for i in range(len(trial_id)):
        tid = trial_id[i]
        trial = "trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
        dict["epochs." + trial + ".description"] = "Data that belong to " + trial
        if options.data_origin in ["NL", "DG"]:
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
#           print "trial=", trial, " epoch_units=", epoch_units
            if trial not in epoch_units:
                units = ["NA"]
            else:
                units = epoch_units[trial]
            try:
                dict["epochs." + trial + ".units_present"] = units
            except:
                print "   Unable to create dataset 'units_present' containing ", units

            if not options.data_origin == "DG":
                raw_path = "descrHash/value/%d" % (trial_id[i])
            else:
                raw_path = "descrHash/value/value"
            raw_file = parse_h5_obj(data_h5[raw_path])[0]
            if len(raw_file) == 1:
                raw_file = 'na'
            else:
                raw_file = str(raw_file)

        elif options.data_origin in ["JY"]:
            start = min(time[:,i].tolist())
            stop  = max(time[:,i].tolist())
            dict["epochs." + trial + ".start_time"] = start
            dict["epochs." + trial + ".stop_time"]  = stop
            # loop through whiskerVars
            keyName = 'whiskerVars'
            hash_group_pointer  = data_h5['timeSeriesArrayHash']
            grp   = libh5.get_value_by_key(hash_group_pointer , keyName)

            valueMatrix = np.array(grp["valueMatrix/valueMatrix"])
            idStr         = np.array(grp["idStr/idStr"])
            try:
                idStrDetailed = np.array(grp["idStrDetailed/idStrDetailed"])
            except:
                idStrDetailed = idStr
            num_vars = len(np.array(grp["id/id"]).tolist())

        elif options.data_origin == "SP":
            # SP data
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

def set_data_origin(data_h5, meta_h5, options):
    options.data_origin = "Unknown"
    keys = [str(i) for i in libh5.get_key_list(data_h5["trialPropertiesHash"])]
    if 'GoodTrials' in keys and not 'LickTime' in keys and not 'EphusVars' in keys:
        options.data_origin = "DG"
    elif not "LickTime"   in libh5.get_key_list(data_h5["trialPropertiesHash"]) and \
       not "EphusVars"  in libh5.get_key_list(data_h5["timeSeriesArrayHash"]):
        options.data_origin = "SP"
    elif "EphusVars" in libh5.get_key_list(data_h5["timeSeriesArrayHash"]):
        options.data_origin = "NL"
    elif "intracellular" in meta_h5.keys():
        # JY data
        options.data_origin = "JY"
    else:
        sys.exit("Data origin not determined")

    return options

# ------------------------------------------------------------------------------

def extract_dict_master_shape(data_h5, plane_map, options, num_planes = 3):
    master_shape = {}
    num_subareas = len(libh5.get_key_list(data_h5['timeSeriesArrayHash'])) - 1
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

def update_dict_master_shape(data_h5, master_shape, plane_map, area, \
                           plane, options, num_plane = 3):
    area_grp = data_h5["timeSeriesArrayHash/descrHash"]["%d"%(1+area)]
    if num_plane == 1:
        plane_grp = area_grp["value/1"]
    else:
        plane_grp = area_grp["value/1"]["%d"%(plane)]
    master = plane_grp["masterImage"]["masterImage"].value
    master1_shape = np.array(master).shape
    green = np.zeros([master1_shape[0],master1_shape[1]])
    red   = np.zeros([master1_shape[0],master1_shape[1]])
    for i in range(master1_shape[0]):
        for j in range(master1_shape[1]):
            if len(master1_shape) == 3 or not options.handle_errors:
                green[i][j] = master[i][j][0]
                red[i][j]   = master[i][j][1]
            else:
                # Warning: only one master image is available, so the green and red signals will be identical
                green[i][j] = master[i][j]
                red[i][j]   = master[i][j]
    # convert from file-specific area/plane mapping to
    #   inter-session naming convention
    oname = "area%d_plane%d" % (area, plane)
    master_shape[oname] = master1_shape

    return master_shape

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

            sys.stdout.write('.')
            sys.stdout.flush()
    return (master_shapes, ref_images_red, ref_images_green)

# ------------------------------------------------------------------------------

def extract_master_shapes(data_h5, options):
    master_shapes   = {}
    num_subareas = len(libh5.get_key_list(data_h5['timeSeriesArrayHash'])) - 1
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if data_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(data_h5[plane_path].keys())
        for plane in range(num_planes):
            area_grp = data_h5["timeSeriesArrayHash/descrHash"]["%d"%(subarea+2)]
            if num_planes == 1:
                plane_grp = area_grp["value/1"]
            else:
                plane_grp = area_grp["value/1"]["%d"%(plane+1)]
            master = plane_grp["masterImage"]["masterImage"].value
            master_shape = np.array(master).shape

            # convert from file-specific area/plane mapping to
            #   inter-session naming convention
            oname = "area%d_plane%d" % (subarea+1, plane+1)
            master_shapes[oname] = master_shape

            sys.stdout.write('.')
            sys.stdout.flush()
    return master_shapes

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
        #print filename + " -> " + univ_name
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
            # if num_planes == 1:
            #     pgrp = grp
            # else:
            #     pgrp = grp["%d"%(plane+1)]
#           print("\nsubarea=", subarea, " plane=", plane
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
#   print "\nfname=", fname
#   print "external_file=", external_file, "\n"
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

def acquisition_timeseries_images(data_h5, options):
    dict = {}
    path0 = "acquisition.timeseries."
    plane_map = create_plane_map(data_h5, options)
    description = "Information about raw image data, which must be obtained from the original data files. 'external_file' provides a list of raw data files; 'external_file.attribute' starting_frame indicates, with 0 as the first value, the frame for this plane at which a given file starts.  Thus, if starting_frame is 0,52,..., then the first file contains frames 1-52, or 0-51. See processing/ROIs/DfOverF for fluorescence change for individual ROIs"
    dict[path0 + "description"] = description
    group_attrs = {"source": "Device 'two-photon microscope'",
                   "description": "2P image stack, one of many sequential stacks for this field of view"}
    dict[path0 + "group_attrs"] = group_attrs
    num_subareas = len(libh5.get_key_list(data_h5['timeSeriesArrayHash'])) - 1
    for subarea in range(num_subareas):
        grp = data_h5["timeSeriesArrayHash"]["value"]["%d"%(subarea+2)]
        t = 0.001 * grp["time"]["time"].value    # fetch time array
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if data_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(data_h5[plane_path].keys())
        for plane in range(num_planes):
            oname = "area%d_plane%d" % (subarea+1, plane+1)
            img_plane = plane_map[oname]
            try:
                pgrp = grp["imagingPlane/%d"%(plane+1)] # move into imaging plane group to extract 2photon image stacks
            except:
                pgrp = grp["imagingPlane"] # only one imaging plane is available (instead of 3)
            frame_idx = pgrp["sourceFileFrameIdx"]["sourceFileFrameIdx"].value
            lst = parse_h5_obj(pgrp["sourceFileList"])[0]
            cnt = 0
            srcfile = {}
            for k in range(len(lst)):
                srcfile[str(k+1)] = lst[k]
                cnt += 1
            assert len(t) == len(frame_idx[0])
            external_file, starting_frame, timestamps, nname = \
                extract_frame_data(frame_idx, t, subarea, plane, plane_map, srcfile, cnt, options)
            dict[path0 + img_plane + ".format"]      = "external"
            # need to convert name to utf8, otherwise error generated:
            dict[path0 + img_plane + ".external_file"]  = external_file
            dict[path0 + img_plane + ".data_attrs"]     = {"starting_frame": starting_frame}
#           dict[path0 + ition/timeseries/" + img_plane + "/dimension"]      = [master_shape[0], master_shape[1]])
            dict[path0 + img_plane + ".dimension"]      = [512, 512]
            dict[path0 + img_plane + ".scan_line_rate"] = 16000.0
            dict[path0 + img_plane + ".field_of_view"]  = [ 600e-6, 600e-6 ]
            dict[path0 + img_plane + ".imaging_plane"]  = img_plane
            dict[path0 + img_plane + ".timestamps"]     = timestamps
            sys.stdout.write('.')
            sys.stdout.flush()
    return dict

# ------------------------------------------------------------------------------

# Extract master images for green and red channels
def acquisition_images(data_h5, options):
    dict = {}
    plane_map = create_plane_map(data_h5, options)
    num_subareas = len(libh5.get_key_list(data_h5['timeSeriesArrayHash'])) - 1
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if data_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(data_h5[plane_path].keys())
        for plane in range(num_planes):
            area_grp = data_h5["timeSeriesArrayHash/descrHash"]["%d"%(subarea + 2)]
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
            path0 = "acquisition.images."
            desc_green = "Master image (green channel), in " + str(master_shape[0]) + \
                         "x" + str(master_shape[1]) + ", 8bit"
            dict[path0 + image_plane + "_green"] = master_green.astype('uint8')
            dict[path0 + image_plane + "_green.attrs"] = {"description":desc_green, "format": "raw"}
            desc_red = "Master image (red channel), in " + str(master_shape[0]) + \
                       "x" + str(master_shape[1]) + ", 8bit"
            dict[path0 + image_plane + "_red"] = master_red.astype('uint8')
            dict[path0 + image_plane + "_red.attrs"] = {"description":desc_red, "format": "raw"}
            sys.stdout.write('.')
            sys.stdout.flush()
    return dict

# ------------------------------------------------------------------------------

def acquisition_timeseries(data_h5, options):
    dict = {}
    if len(options.path_str.split(".")) == 3:
        series_name = options.path_str.split(".")[-1]
        if series_name == "lick_trace":
            # NL or DG data
            dict = acquisition_timeseries_lick_trace(data_h5, options)
        elif series_name == "extracellular_traces":
            dict = acquisition_timeseries_extracellular_traces(data_h5, options)
    elif len(options.path_str.split(".")) == 2:
        # SP data
        dict = acquisition_timeseries_images(data_h5, options)
    else:
        sys.exit("\nUnsupported acquisition data type: " + data_name)
    return dict

# ------------------------------------------------------------------------------

def acquisition_timeseries_lick_trace(data_h5, options):
    dict = {}
    series_name = "lick_trace"
    hash_group_pointer = data_h5["timeSeriesArrayHash"]
    if options.data_origin == "NL":
        keyName     = "EphusVars"
        data        = hash_group_pointer["value/valueMatrix/valueMatrix"][:,0]
        timestamps  = hash_group_pointer["value/time/time"].value
        description = parse_h5_obj(hash_group_pointer["value/idStrDetailed/idStrDetailed"])[0][0]
    else:
        keyName     = "lickVars"
        data        = hash_group_pointer["value/2/valueMatrix/valueMatrix"]
        timestamps  = hash_group_pointer["value/2/time/time"].value
        description = parse_h5_obj(hash_group_pointer["value/2/idStrDetailed/idStrDetailed"])[0][0]
    comment1 = keyName
    comment2 = libh5.get_description_by_key(hash_group_pointer, keyName)
    comments = comment1 + ": " + comment2
    data_attrs={"conversion":1.0, "unit":"unknown", "resolution":float('nan')}
    group_attrs={"description" : description, "comments" : comments, \
                "source": "Times as reported in Nuo's data file", \
                "series_type" : "<TimeSeries>"}

    series_path = "acquisition.timeseries." + series_name + "."
    dict[series_path + "timestamps"]  = timestamps
    dict[series_path + "data"]        = data             
    dict[series_path + "num_samples"] = len(timestamps)
    dict[series_path + "data_attrs"]  = data_attrs
    dict[series_path + "group_attrs"] = group_attrs
    return dict

# ------------------------------------------------------------------------------

def acquisition_timeseries_extracellular_traces(data_h5, options):
    dict = {}
    if options.data_origin == "NL":
        fname = "voltage_filename" + os.path.basename(options.data_path)[13:-3] + ".mat"
    elif options.data_origin == "DG":
        fname = "voltageTraces_" + os.path.basename(options.data_path)[0:-3] + ".mat"
    dict["acquisition.timeseries.extracellular_traces"] = fname
    return dict

    return dict

# ------------------------------------------------------------------------------

def processing_extracellular_units_data(data_h5, meta_h5, options):
    dict = {}
    module_name = "extracellular_units"
    data_item = options.path_str.split(".")[-1]
    if data_item == "top_datasets":
        dict = processing_extracellular_units_top_datasets(data_h5, meta_h5, options)
    elif data_item == "EventWaveform":
        dict = processing_extracellular_units_event_waveform(data_h5, meta_h5, options)
    elif data_item == "UnitTimes":
        dict = processing_extracellular_units_unit_times(data_h5, meta_h5, options)
    return dict
    
# ------------------------------------------------------------------------------

def processing_ROIs(data_h5, options):
    dict = {}
    dict1 = processing_ROIs_DfOverF(data_h5, options)
    dict2 = processing_ROIs_ImageSegmentation(data_h5, options)
    dict.update(dict1)
    dict.update(dict2)
    return dict

# ------------------------------------------------------------------------------

def processing_ROIs_DfOverF(data_h5, options):
    dict = {}
    plane_map = create_plane_map(data_h5, options)
    num_subareas = len(libh5.get_key_list(data_h5['timeSeriesArrayHash'])) - 1
    master_shape = extract_dict_master_shape(data_h5, plane_map, \
                                             options, num_planes = 3)
    dict["processing.ROIs.description"] = "Segmentation (pixel-lists) and dF/F (dffTSA) for all ROIs"
    path = "processing.ROIs.DfOverF."
    dict[path + "group_attrs"]  = {"source" : "This module's DfOverF interface"}
    dict[path + "group2_attrs"] = {"source" : "This module's ImageSegmentation interface"}
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if data_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(data_h5[plane_path].keys())
        for plane in range(num_planes):
            tsah = data_h5["timeSeriesArrayHash"]
            oname = "area%d_plane%d" % (subarea+1, plane+1)
            master1_shape = master_shape[oname]
            image_plane = plane_map[oname]
            area_grp = data_h5["timeSeriesArrayHash/value"]["%d"%(subarea+2)]
            try:
                plane_ids = area_grp["imagingPlane"][str(plane+1)]["ids/ids"].value
            except:
                plane_ids = area_grp["imagingPlane"]["ids/ids"].value
            roi_names = []
            for i in range(len(plane_ids)):
                roi_names.append("ROI%d"%plane_ids[i])
            dict[path + image_plane + ".roi_names"] = roi_names  
            dict[path + image_plane + ".data"] = area_grp["valueMatrix/valueMatrix"].value
            dict[path + image_plane + ".data_attrs"] = \
                {"unit":"dF/F", "conversion":1.0, "resolution":0.0}
            dict[path + image_plane + ".timestamps"] = area_grp["time/time"].value * 0.001
            dict[path + image_plane + ".trial_ids"]  = area_grp["trial/trial"].value
    return dict 

# ------------------------------------------------------------------------------

def processing_ROIs_ImageSegmentation(data_h5, options):
    dict = {}
    plane_map = create_plane_map(data_h5, options)
    num_subareas = len(libh5.get_key_list(data_h5['timeSeriesArrayHash'])) - 1
    master_shape, ref_images_red, ref_images_green = \
        create_reference_images(plane_map, num_subareas, data_h5, options)
    path = "processing.ROIs.ImageSegmentation."
    dict[path + "group_attrs"] = {"source" : "This module's ImageSegmentation interface"}
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if data_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(data_h5[plane_path].keys())
        for plane in range(num_planes):
            tsah = data_h5["timeSeriesArrayHash"]
            oname = "area%d_plane%d" % (subarea+1, plane+1)
            master1_shape = master_shape[oname]
            image_plane = plane_map[oname]
            ref_image_red   = ref_images_red[image_plane]
            ref_image_green = ref_images_green[image_plane]
            dict[path + image_plane + ".imaging_plane_name"] = image_plane
            dict[path + image_plane + ".description"]        = image_plane
            dict[path + image_plane + ".ref_image_red"]      = ref_image_red
            dict[path + image_plane + ".ref_image_green"]    = ref_image_green
            try:
                ids = tsah["value"]["%d"%(subarea+2)]["imagingPlane"][str(plane+1)]["ids"]
            except:
                ids = tsah["value"]["%d"%(subarea+2)]["imagingPlane"]["ids"]
            roi_ids = ids["ids"].value
            dict[path + image_plane + ".roi_ids"] = np.array(roi_ids)
            for i in range(len(roi_ids)):
                rid = roi_ids[i]
                if num_planes == 1:
                    rois = tsah["descrHash"]["%d"%(subarea+2)]["value"]["1"]
                else:
                    rois = tsah["descrHash"]["%d"%(subarea+2)]["value"]["1"]["%d"%(plane+1)]
                try:
                    record = rois["rois"]["%s"%(1+i)]
                    x = int(parse_h5_obj(record["id"])[0])
                    assert x == int(roi_ids[i]) # make sure the ROI id is correct
                except:
                    print("Missing ROI for area=" +str(subarea)+ " plane="+ str(plane) + " id=" +str(i))
                    continue
                pix = parse_h5_obj(record["indicesWithinImage"])[0]
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

def processing_extracellular_units_event_waveform(data_h5, meta_h5, options):
    dict = {}
    group_attrs = {"source" :  "Data as reported in Nuo's file"}
    unit_descr = np.array(parse_h5_obj(data_h5['eventSeriesHash/descr/descr'])[0]).tolist()
    unit_num = len(unit_descr)

    # Populating a dictionary
    series_path = "processing.extracellular_units.EventWaveform"
    dict[series_path + ".group_attrs"] = group_attrs

    for i in range(unit_num):
        i = i+1 
        if options.verbose:
            print "i=", i, " unit_num=", unit_num
        unit = "unit_%d%d" % (int(i/10), i%10)
        if unit_num > 1:
            grp_name = "eventSeriesHash/value/%d" % i
        grp_top_folder = data_h5[grp_name]
        timestamps     = grp_top_folder["eventTimes/eventTimes"]
        trial_ids      = grp_top_folder["eventTrials/eventTrials"]
        waveforms = grp_top_folder["waveforms/waveforms"]
        sample_length = waveforms.shape[1]
        channel = grp_top_folder["channel/channel"].value

        cell_type = parse_h5_obj(grp_top_folder["cellType"])[0]
        if  'numpy' in str(type(cell_type)):
            cells_conc = ' and '.join(map(str, cell_type))
            cell_type = cells_conc
        else:
            cell_type = str(cell_type)

        data_attrs = {"description" : cell_type, "resolution":float('nan'), "unit":"Volts", "conversion":0.1}

        dict[series_path + "." + unit + ".data"]          = waveforms
        dict[series_path + "." + unit + ".timestamps"]    = timestamps
        dict[series_path + "." + unit + ".sample_length"] = len(timestamps)
        dict[series_path + "." + unit + ".data_attrs"]    = data_attrs
        dict[series_path + "." + unit + ".source"]        = "Data from processed matlab file"
        dict[series_path + "." + unit + ".electrode_idx"] = [channel[0]]
        dict[series_path + "." + unit + ".num_samples"]   = len(timestamps)

    return dict

# ------------------------------------------------------------------------------

def processing_extracellular_units_unit_times(data_h5, meta_h5, options):
    dict = {}
    group_attrs = {"source" :  "EventWaveform in this module"}
                           
    unit_descr = np.array(parse_h5_obj(data_h5['eventSeriesHash/descr/descr'])[0]).tolist()
    unit_num = len(unit_descr)
    grp_name = "eventSeriesHash/value"

    cell_types = ['unclassified']*unit_num
    n = max(range(unit_num)) + 1
    electrode_depths = np.zeros(n)

    # Populating a dictionary
    series_path = "processing.extracellular_units.UnitTimes"
    dict[series_path + ".group_attrs"] = group_attrs

    series_path = "processing.extracellular_units.UnitTimes"
    for i in range(unit_num):
        i = i+1
        if options.verbose:
            print "i=", i, " unit_num=", unit_num
        unit = "unit_%d%d" % (int(i/10), i%10)
        if unit_num > 1:
            grp_name = "eventSeriesHash/value/%d" % i
        grp_top_folder = data_h5[grp_name]
        timestamps     = grp_top_folder["eventTimes/eventTimes"]
        trial_ids      = grp_top_folder["eventTrials/eventTrials"]

        cell_type = parse_h5_obj(grp_top_folder["cellType"])[0]
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

def processing_extracellular_units_top_datasets(data_h5, meta_h5, options):
    dict = {}
    description = 'Spike times and waveforms'
    try:
        spike_sorting = libh5.get_value_pointer_by_path_items(meta_h5, ["extracellular", \
                           "spikeSorting", "spikeSorting"])
    except:
        spike_sorting = "method not recorded"

    if options.data_origin == "NL":
        ident_meth1 = str(np.array(libh5.get_value_pointer_by_path_items(meta_h5, ["extracellular", \
                                   "identificationMethod", "identificationMethod"])).tolist()[0])
        ident_meth2 = str(np.array(libh5.get_value_pointer_by_path_items(meta_h5, ["extracellular", \
                                   "identificationMethod", "identificationMethod"])).tolist()[2])
        ident_meth  = ident_meth1 + "\n" + ident_meth2
    else:
        ident_meth  =  str(np.array(libh5.get_value_pointer_by_path_items(meta_h5, ["extracellular", \
                                   "identificationMethod", "identificationMethod"])).tolist()[0])

    # Populating a dictionary
    series_path = "processing.extracellular_units"
    dict[series_path + ".description"]            = description
    dict[series_path + ".spike_sorting"]          = spike_sorting
    dict[series_path + ".identification_method"]  = ident_meth   

    return dict

# ------------------------------------------------------------------------------

def processing_events(module_name, series_name, data_h5, options):
    dict = {}
    if module_name == "Auditory" and series_name == "reward_cue":
        dict = processing_events_auditory_reward_cue(data_h5, options)
    elif module_name == "Licks":
        dict = processing_events_licks(series_name, data_h5, options)
    elif module_name == "Pole":
        dict = processing_events_pole_accassible(data_h5, options)
    elif module_name == "Reward":
        dict = processing_events_reward(series_name, data_h5, options)
    elif module_name == "Whisker":
        dict = processing_events_whisker(series_name, data_h5, options)
    else:
        sys.exit("Extracting processing BehavioralEvents " + series_name + " not implemented")
    return dict

# ------------------------------------------------------------------------------

def processing_timeseries(module_name, series_name, data_h5, options):
    dict = {}
    if module_name == "Whisker":
        dict = processing_timeseries_whisker(series_name, data_h5, options)
    else:
        sys.exit("Extracting processing BehavioralTimeSeries " + series_name + " not implemented")
    return dict

# ------------------------------------------------------------------------------

def processing_timeseries_whisker(series_name, data_h5, options):
    dict = {}

    if options.data_origin == "SP":
        keyName = 'whiskerVars'
        hash_group_pointer  = data_h5['timeSeriesArrayHash']
        grp        = libh5.get_value_by_key(hash_group_pointer , keyName)
        timestamps = grp["time/time"].value * 0.001

        var   = grp["valueMatrix/valueMatrix"].value
        if options.handle_errors:
            try:
                grp = data_h5["timeSeriesArrayHash/value/1"]
            except:
                grp = data_h5["timeSeriesArrayHash/value"]
        else:
            grp = data_h5["timeSeriesArrayHash/value/1"]
        descr = parse_h5_obj(grp["idStrs"])[0]

        if series_name == "whisker_angle":
            data = var[0]
            if options.verbose:
                print "descr[0]=", descr[0]
            group_attrs={"description": "Angle of whiskers", \
                         "source": "Whisker angle as reported in Simon's data file",\
                         "series_type" : "<TimeSeries>"}
            data_attrs ={"unit": "degrees", "conversion":1.0, "resolution":0.001,
                        "keyName" : keyName}

        elif series_name == "whisker_curvature":
            data = var[1]
            if options.verbose:
                print "descr[1]=", descr[1]
            group_attrs={"description": descr[1],
                         "source": "Curvature (relative) of whiskers as reported in Simon's data file", \
                         "series_type" : "<TimeSeries>"}
            data_attrs={"unit":"Unknown", "conversion": 1.0, "resolution": 1.0,
                        "keyName" : keyName}

   # Populating a dictionary
    series_path = "processing.Whisker.BehavioralTimeSeries." + series_name
    dict[series_path + ".timestamps"]  = timestamps
    dict[series_path + ".data"]        = data
    dict[series_path + ".num_samples"] = len(timestamps)
    dict[series_path + ".data_attrs"]  = data_attrs
    dict[series_path + ".group_attrs"] = group_attrs

    return dict

# ------------------------------------------------------------------------------

def processing_events_whisker(series_name, data_h5, options):
    dict = {}
    keyName = "touches"
    hash_group_pointer  =  data_h5["eventSeriesArrayHash"]
    grp = libh5.get_value_by_key(hash_group_pointer , keyName)

    if series_name == "pole_touch_protract":
        timestamps = grp["eventTimes/1/1"].value * 0.001
        data = ones_with_alternating_sign(timestamps)
        group_attrs = {"description" : "Intervals that whisker touches pole (protract)",\
                       "source" : "Intervals are as reported in Simon's data file",
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}
    elif series_name == "pole_touch_retract":
        timestamps = grp["eventTimes/2/2"].value * 0.001
        data = ones_with_alternating_sign(timestamps)
        group_attrs = {"description" : "Intervals that whisker touches pole (retract)",\
                       "source" : "Intervals are as reported in Simon's data file",
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}

    # Populating a dictionary
    series_path = "processing.Whisker.BehavioralEvents." + series_name
    dict[series_path + ".timestamps"]  = timestamps
    dict[series_path + ".data"]        = data
    dict[series_path + ".num_samples"] = len(timestamps)
    dict[series_path + ".data_attrs"]  = data_attrs
    dict[series_path + ".group_attrs"] = group_attrs

    return dict

# ------------------------------------------------------------------------------

def processing_events_reward(series_name, data_h5, options):
    dict = {}
    hash_group_pointer  = data_h5['eventSeriesArrayHash']
    source = "Intervals are as reported in somatosensory cortex data file"

    if options.data_origin == "SP" and series_name == "water_left_reward":
        keyName = "leftReward"
        description = libh5.get_description_by_key(hash_group_pointer, keyName)
        grp = libh5.get_value_by_key(hash_group_pointer , keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = ones_with_alternating_sign(timestamps)
        group_attrs = {"source" : source, "description" : description, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}
    elif options.data_origin == "SP" and series_name == "water_right_reward":
        keyName = "rightReward"
        description = libh5.get_description_by_key(hash_group_pointer , keyName)
        grp = libh5.get_value_by_key(hash_group_pointer , keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = ones_with_alternating_sign(timestamps)
        group_attrs = {"source" : source, "description" : description, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}
    else:
        sys.exit("Function processing_events_reward not implemented for series_name= " + series_name)

    # Populating a dictionary
    series_path = "processing.Reward.BehavioralEvents." + series_name
    dict[series_path + ".timestamps"]  = timestamps
    dict[series_path + ".data"]        = data
    dict[series_path + ".num_samples"] = len(timestamps)
    dict[series_path + ".data_attrs"]  = data_attrs
    dict[series_path + ".group_attrs"] = group_attrs

    return dict

# ------------------------------------------------------------------------------

def processing_events_pole_accassible(data_h5, options):
    dict = {}
    if options.data_origin == "SP":
        keyName = "poleInReach"
        hash_group_pointer = data_h5['eventSeriesArrayHash']
        grp = libh5.get_value_pointer_by_key(hash_group_pointer , keyName, \
                                              options.debug)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = [1] * len(timestamps)
        description = libh5.get_description_by_key(hash_group_pointer , keyName)
        group_attrs = {"source" : "Intervals are as reported in Simon's data file", \
                       "description" : description, "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}

        # Populating a dictionary
        series_path = "processing.Pole.BehavioralEvents.pole_accessible"
        dict[series_path + ".timestamps"]  = timestamps
        dict[series_path + ".data"]        = data
        dict[series_path + ".num_samples"] = len(timestamps)
        dict[series_path + ".data_attrs"]  = data_attrs
        dict[series_path + ".group_attrs"] = group_attrs

    return dict

# ------------------------------------------------------------------------------

def processing_events_licks(series_name, data_h5, options):
    dict = {}
    hash_group_pointer  = data_h5['eventSeriesArrayHash']

    if options.data_origin == "SP" and series_name == "lick_left":
        keyName = 'leftLicks'
        description = libh5.get_description_by_key(hash_group_pointer, keyName)
        grp = libh5.get_value_by_key(hash_group_pointer, keyName)
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

    elif options.data_origin == "SP" and series_name == "lick_right":
        keyName = 'rightLicks'
        description = libh5.get_description_by_key(hash_group_pointer, keyName)
        grp = libh5.get_value_by_key(hash_group_pointer, keyName)
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
    else:
        sys.exit("Function processing_events_licks not implemented for series_name= " + series_name)

    # Populating a dictionary
    series_path = "processing.Licks.BehavioralEvents." + series_name
    dict[series_path + ".timestamps"]  = timestamps
    dict[series_path + ".data"]        = data
    dict[series_path + ".num_samples"] = len(timestamps)
    dict[series_path + ".data_attrs"]  = data_attrs
    dict[series_path + ".group_attrs"] = group_attrs
    dict[series_path + ".module_attrs"]= {"description" : "Lick port contact times"}
    return dict

# ------------------------------------------------------------------------------

def processing_events_auditory_reward_cue(data_h5, options):
    dict = {}
    if options.data_origin == "SP":
        hash_group_pointer  = data_h5['eventSeriesArrayHash']
        keyName = "rewardCue"
        description = libh5.get_description_by_key(hash_group_pointer, keyName)
        grp = libh5.get_value_by_key(hash_group_pointer, keyName)
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
        dict[series_path + ".data_attrs"]  = data_attrs
        dict[series_path + ".group_attrs"] = group_attrs
    return dict

# ------------------------------------------------------------------------------

def stimulus_presentation_timeseries(series_name, data_h5, meta_h5, options):
    dict = {}

    if series_name == "auditory_cue":
        dict = stimulus_presentation_auditory_cue_timeseries(data_h5, options)
    elif re.search("photostimulus", series_name):
        dict = stimulus_presentation_photostimulus_timeseries(data_h5, meta_h5, options)
    elif re.search("pole", series_name):
        dict = stimulus_presentation_pole_timeseries(series_name, data_h5, options)
    elif re.search("water", series_name):
        dict = stimulus_presentation_water_timeseries(series_name, data_h5, options)
    elif series_name == "zaber_motor_pos":
        dict = stimulus_presentation_zaber_motor_pos_timeseries(series_name, data_h5, options)
    else:
        sys.exit("Extracting stimulus presentation " + series_name + " not implemented")

    return dict

# ------------------------------------------------------------------------------

def stimulus_presentation_zaber_motor_pos_timeseries(series_name, data_h5, options):
    dict = {}

    keyName = "StimulusPosition"
    timestamps = data_h5["trialStartTimes/trialStartTimes"].value * 0.001
    data = libh5.get_value_by_key(data_h5['trialPropertiesHash'], keyName)
    description = libh5.get_description_by_key(data_h5["trialPropertiesHash"], keyName)
    group_attrs = {"source" : "Simon's somatosensory cortex data file", "description" : description, \
                   "series_type" : "<TimeSeries>"}
    data_attrs  = {"unit":"unknown", "conversion": 1.0, "resolution":1.0, \
                   "keyName" : keyName}

    # Populating a dictionary
    series_path = "stimulus.presentation.zaber_motor_pos"
    dict[series_path + ".timestamps"]  = timestamps
    dict[series_path + ".data"]        = data
    dict[series_path + ".num_samples"] = len(timestamps)
    dict[series_path + ".data_attrs"]  = data_attrs
    dict[series_path + ".group_attrs"] = group_attrs

    return dict

# ------------------------------------------------------------------------------

def stimulus_presentation_water_timeseries(series_name, data_h5, options):
    dict = {}
    hash_group_pointer  = data_h5['eventSeriesArrayHash']
    source = "Intervals are as reported in somatosensory cortex data file"

    if options.data_origin == "SP" and series_name == "water_left":
        keyName = "leftReward"
        description = libh5.get_description_by_key(hash_group_pointer, keyName)
        grp = libh5.get_value_by_key(hash_group_pointer , keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = ones_with_alternating_sign(timestamps)
        group_attrs = {"source" : source, "description" : description, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "keyName" : keyName}
    elif options.data_origin == "SP" and series_name == "water_right":
        keyName = "rightReward"
        description2 = libh5.get_description_by_key(hash_group_pointer , keyName)
        grp = libh5.get_value_by_key(hash_group_pointer , keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = ones_with_alternating_sign(timestamps)
        group_attrs = {"source" : source, "description" : description2, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0}
    else:
        sys.exit("Function stimulus_presentation_water_timeseries not implemented for series_name= " + series_name)

    # Populating a dictionary
    series_path = "stimulus.presentation." + series_name
    dict[series_path + ".timestamps"]  = timestamps
    dict[series_path + ".data"]        = data
    dict[series_path + ".num_samples"] = len(timestamps)
    dict[series_path + ".data_attrs"]  = data_attrs
    dict[series_path + ".group_attrs"] = group_attrs

    return dict

# ------------------------------------------------------------------------------

def stimulus_presentation_pole_timeseries(series_name, data_h5, options):
    dict = {}

    if options.data_origin == "SP" and series_name == "pole_accessible":
        # create time series for stimulus/presentation group
        keyName = "poleInReach"
        hash_group_pointer = data_h5['eventSeriesArrayHash']
        grp = libh5.get_value_pointer_by_key(hash_group_pointer , keyName, \
                                              options.debug)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = ones_with_alternating_sign(timestamps)
        description = libh5.get_description_by_key(hash_group_pointer , keyName)
        group_attrs = {"source" : "Intervals are as reported in Simon's data file", \
                       "description" : description, \
                       "series_type" : "<IntervalSeries>"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, "keyName" : keyName}
    elif options.data_origin in ["NL", "JY", "DG"]:
        hash_group_pointer = data_h5["trialPropertiesHash"]
        source = "Times as reported in motor cortex data file, but relative to session start"
        trial_start_times = data_h5["trialStartTimes/trialStartTimes"].value
        grp = data_h5["trialPropertiesHash/value/"]
        data_attrs  = {"resolution":0.1,"conversion":1.0}

        if series_name == "pole_in":            
            keyName = "PoleInTime"
            if options.data_origin in ["NL", "JY"]:
                time = grp["1/1"].value
            else:
                time = grp["2/2"].value
            timestamps = time + trial_start_times
            data = [1] * len(timestamps)
            description = libh5.get_description_by_key(hash_group_pointer, keyName)
            data_attrs  = {"keyName" : keyName, "unit" : "sec"}
            group_attrs = {"source": source, "description" : description, \
                           "series_type" : "<TimeSeries>"}
        elif series_name == "pole_out":
            keyName = "PoleOutTime"
            if options.data_origin in ["NL", "JY"]:
                time = grp["2/2"].value
            else:
                time = grp["3/3"].value
            timestamps = time + trial_start_times
            data = [1] * len(timestamps)
            description = libh5.get_description_by_key(hash_group_pointer, keyName)
            data_attrs  = {"keyName" : keyName, "unit" : "sec"}
            group_attrs = {"source": source, "description" : description, \
                           "series_type" : "<TimeSeries>"}
    else:
        sys.exit("Function stimulus_presentation_pole_timeseries not implemented for series_name= " + series_name)

    # Populating a dictionary
    series_path = "stimulus.presentation." + series_name
    dict[series_path + ".timestamps"]  = timestamps
    dict[series_path + ".data"]        = data             
    dict[series_path + ".num_samples"] = len(timestamps)
    dict[series_path + ".data_attrs"]  = data_attrs
    dict[series_path + ".group_attrs"] = group_attrs

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

def stimulus_subtraces_by_type(data_h5, laser_power, time, options):
    keyName = "PhotostimulationType"
    hash_group = data_h5["/trialPropertiesHash"]
    types = libh5.get_value_pointer_by_key(hash_group, keyName, options.debug)
    if options.verbose:
        print("    types="+ str(types))

    trial_t  = data_h5["trialStartTimes/trialStartTimes"].value
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
        sub_traces.append(stimulus_subtrace(laser_power, time, trial_t, types, ['NaN'], options))
    elif count_2 == 0 and count_1 > 0:
        num_types = 1
        trace_types = [1]
        sub_traces.append(stimulus_subtrace(laser_power, time, trial_t, types, ['NaN'], options))
    elif count_1 > 0 and count_2 > 0:
        num_types = 2
        trace_types = [1, 2]
        sub_traces.append(stimulus_subtrace(laser_power, time, trial_t, types, ['2','NaN'], options))
        sub_traces.append(stimulus_subtrace(laser_power, time, trial_t, types, ['1','NaN'], options))
    return (num_types, trace_types, sub_traces)

# ------------------------------------------------------------------------------

def stimulus_presentation_photostimulus_timeseries(data_h5, meta_h5, options):
    dict = {}
    grp_name = "timeSeriesArrayHash/value/time/time"
    timestamps = np.array(data_h5[grp_name].value)
    # calculate sampling rate
    rate = (timestamps[-1] - timestamps[0])/(len(timestamps)-1)
    # get descriptions
    comment1 = parse_h5_obj(data_h5["timeSeriesArrayHash/keyNames"])[0][0]
    comment2 = parse_h5_obj(data_h5["timeSeriesArrayHash/descr"])[0][0]
    comments = comment1 + ": " + comment2
    grp_name = "timeSeriesArrayHash/value/idStrDetailed"
    description = parse_h5_obj(data_h5[grp_name])[0]

    # laser data
    keyName1 = "EphusVars"
    keyName2 = "laser_power"
    hash_group_pointer = data_h5["timeSeriesArrayHash"]
    laser_power = libh5.get_value2_by_key2(hash_group_pointer, keyName1, "",keyName2)

    num_traces, trace_types, traces = \
        stimulus_subtraces_by_type(data_h5, laser_power, timestamps, options)
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
                coord = str(libh5.get_value_by_path_items(meta_h5, ["photostim", "photostimCoordinates"]))
                loc   = str(libh5.get_value_by_path_items(meta_h5, ["photostim", "photostimLocation"]))
                series_path = "stimulus.presentation.photostimulus_" + str(trace_types[s])
                dict[series_path + ".timestamps"]  = timestamps
                dict[series_path + ".data"]        = traces[s]  
                dict[series_path + ".num_samples"] = len(timestamps)
                dict[series_path + ".site"] = str(coord + " in mm\natlas location: " + loc)
                dict[series_path + ".data_attrs"]   = data_attrs
                dict[series_path + ".group_attrs"]  = group_attrs
            except:
                print "No photostimulus info found in the metadata file"
    return dict

# ------------------------------------------------------------------------------

def stimulus_presentation_auditory_cue_timeseries(data_h5, options):
    dict = {}
    series_name = "auditory_cue"

    if options.data_origin == "NL":
        keyName = "CueTime"
        hash_group_pointer  = data_h5['trialPropertiesHash']
        trial_start_times = data_h5["trialStartTimes/trialStartTimes"].value
        timestamps = hash_group_pointer["value/3/3"].value + trial_start_times
        data = ones_with_alternating_sign(timestamps)
        description = hash_group_pointer["descr/descr"][2]
        description = libh5.get_description_by_key(hash_group_pointer, keyName)
        group_attrs = {"comments" : description,\
                       "description" : keyName, \
                       "source" : "Times are as reported in Nuo's data file, but relative to session time", \
                       "series_type" : "<TimeSeries>"}
        data_attrs = {"unit": "None", "conversion": 1.0, "resolution": 0.0, "duration" : 0.1, \
                      "keyName" : keyName}
    elif options.data_origin == "SP":
        keyName = "rewardCue"
        hash_group_pointer  = data_h5['eventSeriesArrayHash']
        description = libh5.get_description_by_key(hash_group_pointer, keyName)
        grp = libh5.get_value_by_key(hash_group_pointer, keyName)
        timestamps = grp["eventTimes/eventTimes"].value * 0.001
        data = [1] * len(timestamps)
        group_attrs = {"description" : "Intervals when auditory cue presented",\
                       "source": "Intervals are as reported in Simon's data file", \
                       "series_type" : "<IntervalSeries>"}
        data_attrs = {"unit": "None", "conversion": 1.0, "resolution": 0.0, \
                      "keyName" : keyName}
    else:
        sys.exit("Function stimulus_presentation_auditory_cue_timeseries not implemented for series_name= " + series_name)

    # Populating a dictionary
    series_path = "stimulus.presentation." + series_name
    dict[series_path + ".timestamps"]  = timestamps
    dict[series_path + ".data"]        = data
    dict[series_path + ".num_samples"] = len(timestamps)
    dict[series_path + ".data_attrs"]  = data_attrs
    dict[series_path + ".group_attrs"] = group_attrs

    return dict



