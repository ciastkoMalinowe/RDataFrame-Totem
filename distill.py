import ROOT
import glob

RDF = ROOT.ROOT.RDataFrame

# Extracted from: DS1/block1/input_files.h
input_ntuple_name = "TotemNtuple"
prefix            = "/eos/totem/data/cmstotem/2015/90m/Totem/Ntuple/version2/4495/"
input_files       = prefix + "Totem*" # Expands to http://cern.ch/go/RT7v

# Input branches
# Considering: verticals in 45 and 56

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
               "y_L_1_F", "y_L_2_N", "y_L_2_F", "y_R_1_F", "y_R_2_N", "y_R_2_F"]

# Convert to PyROOT vector
vec_outbranchlist = ROOT.vector('string')()
[vec_outbranchlist.push_back(b) for b in branchList]

# Filter and define output branches
r = rdf.Filter("(track_rp_5.valid   + track_rp_21.valid   + track_rp_25.valid ) >= 2 && " \
             + "(track_rp_105.valid + track_rp_121.valid  + track_rp_125.valid ) >= 2")  \
       .Define("v_L_1_F", "track_rp_5.valid")   \
       .Define("v_L_2_N", "track_rp_21.valid")  \
       .Define("v_L_2_F", "track_rp_25.valid")  \
       .Define("v_R_1_F", "track_rp_105.valid") \
       .Define("v_R_2_N", "track_rp_121.valid") \
       .Define("v_R_2_F", "track_rp_125.valid") \
       .Define("x_L_1_F", "track_rp_5.x")       \
       .Define("x_L_2_N", "track_rp_21.x")      \
       .Define("x_L_2_F", "track_rp_25.x")      \
       .Define("x_R_1_F", "track_rp_105.x")     \
       .Define("x_R_2_N", "track_rp_121.x")     \
       .Define("x_R_2_F", "track_rp_125.x")     \
       .Define("y_L_1_F", "track_rp_5.y")       \
       .Define("y_L_2_N", "track_rp_21.y")      \
       .Define("y_L_2_F", "track_rp_25.y")      \
       .Define("y_R_1_F", "track_rp_105.y")     \
       .Define("y_R_2_N", "track_rp_121.y")     \
       .Define("y_R_2_F", "track_rp_125.y")

print("Distilled events: %s" % r.Count().GetValue())

# Save output tree
r.Snapshot(outTreeName, outFileName, vec_outbranchlist)
