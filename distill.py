import ROOT
import glob

RDF = ROOT.ROOT.RDataFrame

# Extracted from: DS1/block1/input_files.h
input_ntuple_name = "TotemNtuple"
prefix            = "/eos/totem/data/cmstotem/2015/90m/Totem/Ntuple/version2/4495/"
input_files       = [prefix + line.rstrip('\n') for line in open(source_file)]

# Input branches
# Considering: verticals in 45 and 56 (d45b_56t)

# No need to pass them to RDF
# branches = ["event_info.*", "trigger_data.*",
#             "track_rp_5.*", "track_rp_21.*", "track_rp_25.*",
#             "track_rp_105.*", "track_rp_121.*", "track_rp_125.*",
#            ]

treename= "TotemNtuple"
rdf = RDF(treename, input_files)

# Output tree, file and branches
outTreeName = "distilled"
outFileName = "distill_DS1_new.root"
branchList  = ["v_L_1_F", "v_L_2_N", "v_L_2_F", "v_R_1_F", "v_R_2_N", "v_R_2_F",
               "x_L_1_F", "x_L_2_N", "x_L_2_F", "x_R_1_F", "x_R_2_N", "x_R_2_F",
               "y_L_1_F", "y_L_2_N", "y_L_2_F", "y_R_1_F", "y_R_2_N", "y_R_2_F",
               "run_num",
               "bunch_num",
               "event_num",
               "trigger_num",
               "trigger_bits"
              ]

# Convert to PyROOT vector
vec_outbranchlist = ROOT.vector('string')()
[vec_outbranchlist.push_back(b) for b in branchList]

# Filter and define output branches
# Diagonal d45b_56t
r = rdf.Filter("(track_rp_5.valid   + track_rp_21.valid   + track_rp_25.valid ) >= 2 && " \
             + "(track_rp_104.valid + track_rp_120.valid  + track_rp_124.valid ) >= 2")  \
       .Define("v_L_1_F",      "track_rp_5.valid")   \
       .Define("v_L_2_N",      "track_rp_21.valid")  \
       .Define("v_L_2_F",      "track_rp_25.valid")  \
       .Define("v_R_1_F",      "track_rp_104.valid") \
       .Define("v_R_2_N",      "track_rp_120.valid") \
       .Define("v_R_2_F",      "track_rp_124.valid") \
       .Define("x_L_1_F",      "track_rp_5.x")       \
       .Define("x_L_2_N",      "track_rp_21.x")      \
       .Define("x_L_2_F",      "track_rp_25.x")      \
       .Define("x_R_1_F",      "track_rp_104.x")     \
       .Define("x_R_2_N",      "track_rp_120.x")     \
       .Define("x_R_2_F",      "track_rp_124.x")     \
       .Define("y_L_1_F",      "track_rp_5.y")       \
       .Define("y_L_2_N",      "track_rp_21.y")      \
       .Define("y_L_2_F",      "track_rp_25.y")      \
       .Define("y_R_1_F",      "track_rp_104.y")     \
       .Define("y_R_2_N",      "track_rp_120.y")     \
       .Define("y_R_2_F",      "track_rp_124.y")     \
       .Define("run_num",      "event_info.run_no")             \
       .Define("bunch_num",    "trigger_data.bunch_num")        \
       .Define("event_num",    "trigger_data.event_num")        \
       .Define("trigger_num",  "trigger_data.trigger_num")      \
       .Define("trigger_bits", "trigger_data.input_status_bits")

# All above actions are not executed at the moment they are called,
# but they are lazy, i.e. delayed until the moment one of their results
# is accessed (in this case by .GetValue() )
print("Distilled events: %s" % r.Count().GetValue())

# Save output tree
r.Snapshot(outTreeName, outFileName, vec_outbranchlist)