import ROOT
import glob

RDF = ROOT.ROOT.RDataFrame

# Due to fact that not all files from 'prefix' are used in original analysis 
# the desired ones are specified in the input_files_DS1.txt
source_file       = "input_files_DS1.txt"
input_ntuple_name = "TotemNtuple"
prefix            = "/eos/totem/data/cmstotem/2015/90m/Totem/Ntuple/version2/4495/"
input_files       = [prefix + line.rstrip('\n') for line in open(source_file)]

# Convert to PyROOT vector
vec_input_files = ROOT.vector('string')()
[vec_input_files.push_back(f) for f in input_files]
print("done")

# Returns filtered rdf
def filter_rdf(rdf, L_1_F, L_2_N, L_2_F, R_1_F, R_2_N, R_2_F):

    # Filter and define output branches
    r = rdf.Filter("({}.valid   + {}.valid   + {}.valid ) >= 2 && ".format(L_1_F, L_2_N, L_2_F) \
             + "({}.valid + {}.valid  + {}.valid ) >= 2".format(R_1_F, R_2_N, R_2_F))  \
        .Define("v_L_1_F", "{}.valid".format(L_1_F)) \
        .Define("v_L_2_N", "{}.valid".format(L_2_N)) \
        .Define("v_L_2_F", "{}.valid".format(L_2_F)) \
        .Define("v_R_1_F", "{}.valid".format(R_1_F)) \
        .Define("v_R_2_N", "{}.valid".format(R_2_N)) \
        .Define("v_R_2_F", "{}.valid".format(R_2_F)) \
        .Define("x_L_1_F", "{}.x".format(L_1_F)) \
        .Define("x_L_2_N", "{}.x".format(L_2_N)) \
        .Define("x_L_2_F", "{}.x".format(L_2_F)) \
        .Define("x_R_1_F", "{}.x".format(R_1_F)) \
        .Define("x_R_2_N", "{}.x".format(R_2_N)) \
        .Define("x_R_2_F", "{}.x".format(R_2_F)) \
        .Define("y_L_1_F", "{}.y".format(L_1_F)) \
        .Define("y_L_2_N", "{}.y".format(L_2_N)) \
        .Define("y_L_2_F", "{}.y".format(L_2_F)) \
        .Define("y_R_1_F", "{}.y".format(R_1_F)) \
        .Define("y_R_2_N", "{}.y".format(R_2_N)) \
        .Define("y_R_2_F", "{}.y".format(R_2_F)) \
        .Define("timestamp", "event_info.timestamp - 1444860000")  \
        .Define("run_num", "event_info.run_no")                    \
        .Define("bunch_num", "trigger_data.bunch_num")             \
        .Define("event_num", "trigger_data.event_num")             \
        .Define("trigger_num", "trigger_data.trigger_num")         \
        .Define("trigger_bits", "trigger_data.input_status_bits")
        
    print("Distilled events: %s" % r.Count().GetValue())
    return r

def save_rdf(filtered_rdf, diagonal):
    # Output tree, file and branches
    outTreeName = "distilled"
    outFileName = "distill_DS1_" + diagonal + "_new.root"
    branchList  = ["v_L_1_F", "v_L_2_N", "v_L_2_F", "v_R_1_F", "v_R_2_N", "v_R_2_F",
                   "x_L_1_F", "x_L_2_N", "x_L_2_F", "x_R_1_F", "x_R_2_N", "x_R_2_F",
                   "y_L_1_F", "y_L_2_N", "y_L_2_F", "y_R_1_F", "y_R_2_N", "y_R_2_F",
                   "timestamp", "run_num", "bunch_num", "event_num", "trigger_num", "trigger_bits"]
    
    # Convert to PyROOT vector
    vec_outbranchlist = ROOT.vector('string')()
    [vec_outbranchlist.push_back(b) for b in branchList]
    
    filtered_rdf.Snapshot(outTreeName, outFileName, vec_outbranchlist)
    

# 45b_56t
diagonal = "d45b_56t"

L_1_F = "track_rp_5"
L_2_N = "track_rp_21"
L_2_F = "track_rp_25"
R_1_F = "track_rp_104"
R_2_N = "track_rp_120"
R_2_F = "track_rp_124"

treename= "TotemNtuple"
rdf = RDF(treename, vec_input_files)

r = filter_rdf(rdf, L_1_F, L_2_N, L_2_F, R_1_F, R_2_N, R_2_F)
#save_rdf(r, diagonal)

# 45t_56b
diagonal = "d45t_56b"

L_1_F = "track_rp_4"
L_2_N = "track_rp_20"
L_2_F = "track_rp_24"
R_1_F = "track_rp_105"
R_2_N = "track_rp_121"
R_2_F = "track_rp_125"

treename= "TotemNtuple"
rdf = RDF(treename, vec_input_files)

r = filter_rdf(rdf, L_1_F, L_2_N, L_2_F, R_1_F, R_2_N, R_2_F)
#save_rdf(r, diagonal)

# 45b_56b
diagonal = "ad45b_56b"

L_1_F = "track_rp_5"
L_2_N = "track_rp_21"
L_2_F = "track_rp_25"
R_1_F = "track_rp_105"
R_2_N = "track_rp_121"
R_2_F = "track_rp_125"

treename= "TotemNtuple"
rdf = RDF(treename, vec_input_files)

r = filter_rdf(rdf, L_1_F, L_2_N, L_2_F, R_1_F, R_2_N, R_2_F)
#save_rdf(r, diagonal)

# 45t_56t
diagonal = "ad45t_56t"

L_1_F = "track_rp_4"
L_2_N = "track_rp_20"
L_2_F = "track_rp_24"
R_1_F = "track_rp_104"
R_2_N = "track_rp_120"
R_2_F = "track_rp_124"

treename= "TotemNtuple"
rdf = RDF(treename, vec_input_files)

r = filter_rdf(rdf, L_1_F, L_2_N, L_2_F, R_1_F, R_2_N, R_2_F)
#save_rdf(r, diagonal)