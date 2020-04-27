"""

GISHydroNXT is a python version of GISHydro2000 legacy software. Legacy version was built using Avenue language, which is obsolete and no longer supported by ESRI for future development.
Current ESRI standard is python language. All tools and menu items here, as needed, are developed using python to replicate functionalities of GISHydro legacy version. For GUI, pythonaddins
(as part of default v10.1 installation) and wxPython (2.8.12.1 (msw-unicode)) are used. wxPython version can be obtained from:

http://sourceforge.net/projects/wxpython/files/wxPython/2.8.12.1/wxPython2.8-win32-unicode-2.8.12.1-py27.exe/download?use_mirror=iweb

and we  use "wxPython2.8-win64-unicode-py27"  for 64-bit Python 2.7 on windows 7. Please run setup and select "Full Installation" option during installation. Also configure and
install .pyc files as prompted. Also make sure to install Service Pack 1 for ArcGIS 10.1 (Detailed installation instructions are provided separately as well).

"""
#
#
# *******************************************************************************************************
# Purpose:      Master script to implement tools/buttons and call other modules/functions, as needed.
# Calls:        hydro, tasker
# By:           Javier Mardones & Ibraheem Kahn
# *******************************************************************************************************
#
Modifieddt = "February 13, 2020"       ### CHANGE UPDATE DATE HERE
ArcGIS_Version = "10.5.1"
Directory = r"C:\GISHydroNXT\umdgism/"        ### CHANGE TO MAKE HELP FUNCTIONS WORK
thomasversion = "2020"
#
# *******************************************************************************************************
# import all standard required modules and functions
# *******************************************************************************************************

import pythonaddins
import arcpy

try:
    from operator import mul
    import wx
    import glob
    import os
    import sys
    import numpy as np
    import time
    from datetime import datetime
    from itertools import groupby
    from operator import itemgetter
    import wx.lib.plot as plot
    import wx.grid as gridlib
    import shutil
    import shelve
    import csv
    import getpass
except ImportError as e:
    pythonaddins.MessageBox(e.message, "Python Libraries Missing")
app = wx.PySimpleApp()

local_path = os.path.dirname(__file__)
sys.path.insert(0, local_path)

# *******************************************************************************************************
# import hydro, hydro.BasinStat, hydro.OuputFRresults, tasker
# *******************************************************************************************************
try:
    import hydro
except ImportError as e:
    pythonaddins.MessageBox(e.message, "ImportError")
try:
    import tasker
except ImportError as e:
    pythonaddins.MessageBox(e.message, "ImportError")


# *******************************************************************************************************
# Global variables of selected input file types
# *******************************************************************************************************
soil = ""
landuse = ""
hyd = ""

# *******************************************************************************************************
# Global variables
# *******************************************************************************************************
optfolder = ""
scratchfolder = ""
proj = ""
demtot = ""
xoutletstring = ""
youtletstring = ""
FC = ""
LI = ""
areami2 = ""
pctAsoil = ""
pctBsoil = ""
pctCsoil = ""
pctDsoil = ""
basinrelief = ""
theslope = ""
landslope = ""
regionlist = ""
regionarea = ""
maprec = ""
extent = ""
provstring = ""
xoutlet = ""
youtlet = ""
avgCN = ""
tc = ""
lagtime = ""
p2yr = ""
lenlst = ""
slopeval = ""
no_subwatersheds = ""
cn_list = ""
tc_list = ""
outletelev = ""
gagelist = ""
trans_distance = ""
trans_dem = ""
TWL = ""
TWL_min = ""
TWL_max = ""
areami2_usda = ""
reachslope = ""
Wbf = ""
Dbf = ""
gageid = ""
Tc_method = ""
Tc_ns = ""
Tc_P = ""
Tc_L = ""
Tc_paved = ""
Tc_unpaved = ""
Tc_NHD = ""
Tc_infStreams = ""
Tc_sa = ""
Tc_nc = ""
Tc_cwCoef = ""
Tc_cwExp = ""
Tc_cdCoef = ""
Tc_cdExp = ""
Tc_caCoef = ""
Tc_caExp = ""
Tc_allsub = ""
Tc_selsub = ""
year = ""
critdur = ""
critavg = ""
cb_selected = ""
year_uSpecified = ""
prec_uSpecified = ""
overall_tc = ""
overland_tc = ""
swale_tc = ""
channel_tc = ""
over_seg = ""
swale_seg = ""
channel_seg = ""
firstover = ""
lastover = ""
firstswale = ""
lastswale = ""
firstchannel = ""
lastchannel = ""
seglist = ""
subValue = ""
arcid_global = ""
traceback = ""
attlist = ""
landedit = ""
struc_name = ""
reservoir_elev = ""
precip_str = ""


# *******************************************************************************************************
# Menu option classes with sequence of processing
# *******************************************************************************************************

#                                           ++++++++++++++++
#                                           +              +
#                                           +  MAIN LOGIC  +
#                                           +              +
#                                           ++++++++++++++++



class AddOutlets(object):
    """Implementation for GISHydroNXT_addin.button7 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        # *******************************************************************************************************
        # process "AddOutlets" raster to handle NoData values -- as it has to be added to "outlets" later
        # *******************************************************************************************************
        arcpy.PointToRaster_conversion(optfolder + "/AddasOutlets.shp", "FID", optfolder + "/Outlets_temp",
                                       "MOST_FREQUENT", "NONE", 30)
        outlets_adj = arcpy.sa.Plus(optfolder + "/Outlets_temp", 1)
        outlets_adj.save(optfolder + "/AddOutlets")
        outlets_custom = arcpy.sa.Con(arcpy.sa.IsNull(optfolder + "/AddOutlets"), 0, optfolder + "/AddOutlets")
        outlets_custom = arcpy.sa.SetNull(outlets_custom, outlets_custom, "VALUE = 0")
        outlets_custom.save(optfolder + "/outlets_user")

        #pythonaddins.MessageBox("Outlets added successfully","Add Outlets")
        print("test")
        # *******************************************************************************************************
        # Turn Add streams OFF and Delineate Subwatersheds ON
        # *******************************************************************************************************
        tool5.enabled = False
        button6.enabled = False
        button7.enabled = False
        save(optfolder)

class AddStreams(object):
    """Implementation for GISHydroNXT_addin.button6 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        # *******************************************************************************************************
        # add traced streams to view -- addOutputsToMap
        # *******************************************************************************************************
        arcpy.env.addOutputsToMap = True
        arcpy.env.extent = "MAXOF"
        arcpy.Merge_management(optfolder + "/AddasStreams.shp",
                               optfolder + "/StrmMerge.shp")  # no apparent benefit of using merge
        arcpy.env.snapRaster = optfolder + "/flowacc"
        arcpy.FeatureToRaster_conversion(optfolder + "/StrmMerge.shp", "Id", optfolder + "/ModStr", "#")
        ModStr = arcpy.sa.Times(optfolder + "/ModStr", optfolder + "/basingrid")
        Streams = arcpy.sa.Con(ModStr == 0, 1)
        traced = arcpy.sa.Times(Streams, optfolder + "/basingrid")
        arcpy.env.addOutputsToMap = False
        traced.save(optfolder + "/ModStreams")  # add it to TOC
        modstreams = optfolder + "/ModStreams"

        # *******************************************************************************************************
        # turn layers ON/OFF in current data frame
        # *******************************************************************************************************
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        modlayer = arcpy.mapping.Layer(modstreams)  # creating new layer
        arcpy.mapping.AddLayer(df, modlayer, "TOP")
        for lyr in layers:
            if lyr.name == "StrmMerge":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "flowdir_dem":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "line":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "ModStr":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "ModStreams":
                #arcpy.mapping.RemoveLayer(df, lyr)
                arcpy.ApplySymbologyFromLayer_management(lyr, r"" + Directory + "/data/mdfiles/legends/ModStreams.lyr")
            if lyr.name == "AddasStreams":
                lyr.visible = False

        #pythonaddins.MessageBox("Streams added successfully","Add Streams")

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()

        # *******************************************************************************************************
        # Turn Add streams OFF and Delineate Subwatersheds ON
        # *******************************************************************************************************

        tool4.enabled = False
        tool5.enabled = True
        button5.enabled = True
        button6.enabled = False
        button8.enabled = True
        save(optfolder)

class AddSubwatershedOutlets(object):
    """Implementation for GISHydroNXT_addin.tool5 (Tool)"""

    def __init__(self):
        self.enabled = False
        self.shape = 3

    def onMouseDownMap(self, x, y, button, shift):
        # *******************************************************************************************************
        # a) check if "AddOutlets" raster exist then delete it if re-running tool
        # b) add xy to "AddasOutlets.shp" recursively and convert it to a raster
        # *******************************************************************************************************
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        xy = (x, y)

        outletxy = arcpy.sa.ExtractByPoints(optfolder + "/InfStreams", [arcpy.Point(x, y)], "INSIDE")
        outletxy.save(optfolder + "/outletxy")

        aux_max = arcpy.sa.Raster(optfolder + "/outletxy").maximum
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        for lyr in arcpy.mapping.ListLayers(mxd, "", df):
            if lyr.name == "outletxy":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "AddasOutlets":
                arcpy.mapping.RemoveLayer(df, lyr)

        arcpy.Delete_management(optfolder + "/outletxy", "")
        del outletxy

        fullPath = optfolder + "/AddasOutlets.shp"
        addLayer = arcpy.mapping.Layer(fullPath)
        arcpy.mapping.AddLayer(df, addLayer)

        if aux_max <= 0:
            pythonaddins.MessageBox("Please select a valid Stream", "Outlet Error")
            return

        cursor = arcpy.da.InsertCursor(optfolder + "/AddasOutlets.shp", ("SHAPE@XY"))
        cursor.insertRow([xy])

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()

        # *******************************************************************************************************
        # Turn Add streams OFF and Delineate Subwatersheds ON
        # *******************************************************************************************************
        button7.enabled = True
        save(optfolder)

class AreaOfInterest(object):
    """Implementation for GISHydroNXT_addin.tool1 (Tool)"""

    def __init__(self):
        self.enabled = True
        self.checked = False
        self.cursor = 5
        self.shape = "Rectangle"

    def onRectangle(self, rectangle_geometry):
        """Initiated, when the rectangle is drawn and the mouse button is released.
        The rectangle is an extent object. County quads will be visible in the startup mxd file as a reference, to draw conservative rectangle."""
        # optfolder will be created once user select data of interest (moved under "Data of Interest" class)
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        global extent
        extent = rectangle_geometry
        DataSelectionFrame()


class AttributeTable(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Segment Attributes - Sub-area " + str(arcid_global + 1), size=(1240, 570))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)
        self.index = 0
        self.SetPosition((530, 420))

        global upval, thename_lst, thetype, downval, avgarea, dsarea, upelev, downelev, theslope, thechanwidth, thechandepth, thechanarea, inclength, downlength, thevel, thetinc, thettot
        # define attribute table frame and loop over attributes
        self.list_ctrl = wx.ListCtrl(panel, size=(1220, 520),
                                     style=wx.LC_REPORT | wx.LC_HRULES | wx.LC_VRULES | wx.BORDER_SUNKEN)
        self.list_ctrl.InsertColumn(0, "FID", width=70)
        #        self.list_ctrl.InsertColumn(1, "Shape", width=80)
        self.list_ctrl.InsertColumn(1, "UpPixel", width=70)
        self.list_ctrl.InsertColumn(2, "SegName", width=70)
        self.list_ctrl.InsertColumn(3, "Type", width=70)
        self.list_ctrl.InsertColumn(4, "DownPixel", width=70)
        self.list_ctrl.InsertColumn(5, "Avg. Area", width=70)
        #        self.list_ctrl.InsertColumn(7, "DS Area", width=80)
        self.list_ctrl.InsertColumn(6, "UpElev", width=70)
        self.list_ctrl.InsertColumn(7, "DownElev", width=70)
        self.list_ctrl.InsertColumn(8, "Slope", width=70)
        self.list_ctrl.InsertColumn(9, "Width", width=70)
        self.list_ctrl.InsertColumn(10, "Depth", width=70)
        self.list_ctrl.InsertColumn(11, "Xarea", width=70)
        self.list_ctrl.InsertColumn(12, "I_Length", width=70)
        self.list_ctrl.InsertColumn(13, "Tot_Length", width=70)
        self.list_ctrl.InsertColumn(14, "Vel.", width=70)
        self.list_ctrl.InsertColumn(15, "I_Time", width=70)
        self.list_ctrl.InsertColumn(16, "Tot_Time", width=70)

        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder

        shape = ["polyline"] * len(
            upval)  # since all shape types are "polyline" so simply multiplied it with length of records
        rowin = np.arange(0, len(upval), dtype=int)

        for a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q in reversed(
                zip(rowin, upval, thename_lst, thetype, downval, avgarea, upelev, downelev, theslope, thechanwidth,
                    thechandepth, thechanarea, inclength, downlength, thevel, thetinc, thettot)):
            self.list_ctrl.InsertStringItem(self.index, str(a))
            #            self.list_ctrl.SetStringItem(self.index, 1, str(b))
            self.list_ctrl.SetStringItem(self.index, 1, str(b))
            self.list_ctrl.SetStringItem(self.index, 2, str(c))
            self.list_ctrl.SetStringItem(self.index, 3, str(d))
            self.list_ctrl.SetStringItem(self.index, 4, str(e))
            self.list_ctrl.SetStringItem(self.index, 5, str(f))
            #            self.list_ctrl.SetStringItem(self.index, 7, str(h))
            self.list_ctrl.SetStringItem(self.index, 6, str(g))
            self.list_ctrl.SetStringItem(self.index, 7, str(h))
            self.list_ctrl.SetStringItem(self.index, 8, str(i))
            self.list_ctrl.SetStringItem(self.index, 9, str(j))
            self.list_ctrl.SetStringItem(self.index, 10, str(k))
            self.list_ctrl.SetStringItem(self.index, 11, str(l))
            self.list_ctrl.SetStringItem(self.index, 12, str(m))
            self.list_ctrl.SetStringItem(self.index, 13, str(n))
            self.list_ctrl.SetStringItem(self.index, 14, str(o))
            self.list_ctrl.SetStringItem(self.index, 15, str(p))
            self.list_ctrl.SetStringItem(self.index, 16, str(q))

        self.Show(True)

        global attlist
        if not os.path.exists(optfolder + "/vel_meth/attribute_tables"):
            os.makedirs(optfolder + "/vel_meth/attribute_tables")
        with open(optfolder + "/vel_meth/attribute_tables/" + "Segment_Attributes_Sub-area_" + str(
                arcid_global + 1) + ".csv", "wb") as f:
            w = csv.writer(f)
            attlist = list(
                zip(rowin, shape, upval, thename_lst, thetype, downval, avgarea, dsarea, upelev, downelev, theslope,
                    thechanwidth,
                    thechandepth, thechanarea, inclength, downlength, thevel, thetinc, thettot))
            w.writerow(["FID", "Shape", "UpPixel", "SegName", "Type", "DownPixel", "Avg. Area", "DS Area", "UpElev",
                        "DownElev", "Slope",
                        "Width", "Depth", "Xarea", "I_Length", "Tot_Length", "Vel.", "I_Time", "Tot_Time"])
            w.writerows(attlist)

    def OnClose(self, event):
        self.Show(False)


class BasinComposition(object):
    """Implementation for GISHydroNXT_addin.button2 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):

        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            ## input rasters
            lu = optfolder + "/landuse"
            soil = optfolder + "/Soils"
            wshed = optfolder + "/watershed.shp"

            ## create a folder for basin composition files
            basincomp = optfolder + "/basincomp"
            if not os.path.exists(basincomp):
                os.makedirs(basincomp)

            ## extract by mask to extent of watershed
            lu_ext = optfolder + "/basincomp/lu_ext"
            soil_ext = optfolder + "/basincomp/soil_ext"
            lu_out = arcpy.sa.ExtractByMask(lu, wshed)
            lu_out.save(lu_ext)
            soil_out = arcpy.sa.ExtractByMask(soil, wshed)
            soil_out.save(soil_ext)

            ## convert all clipped rasters to polygons
            lu_poly = optfolder + "/basincomp/lu_poly.shp"
            soil_poly = optfolder + "/basincomp/soil_poly.shp"
            arcpy.env.addOutputsToMap = False
            arcpy.RasterToPolygon_conversion(lu_ext, lu_poly, "NO_SIMPLIFY", "VALUE")
            arcpy.env.addOutputsToMap = False
            arcpy.RasterToPolygon_conversion(soil_ext, soil_poly, "NO_SIMPLIFY", "VALUE")

            ## intersect land use and soil to prepare two polygons: "lu_soil" and "lu_cn"
            arcpy.env.addOutputsToMap = False
            lu_soil = optfolder + "/basincomp/lu_soil.shp"
            arcpy.Intersect_analysis([lu_poly, soil_poly], lu_soil, "ALL", "#", "INPUT")

            ## dissolve above intersected polygon
            arcpy.env.addOutputsToMap = False
            lu_soil_diss = optfolder + "/basincomp/lu_soil_diss.shp"
            arcpy.Dissolve_management(lu_soil, lu_soil_diss, "GRIDCODE;GRIDCODE_1", "#", "MULTI_PART", "DISSOLVE_LINES")

            ## add filed to both of above dissolved polygons and compute area in acres
            if not len(arcpy.ListFields(lu_soil_diss, "area")) > 0:
                arcpy.AddField_management(lu_soil_diss, "area", "FLOAT", 15, 4)
            arcpy.CalculateField_management(lu_soil_diss, "area", "!shape.area@acres!", "PYTHON")

            # prepre a list of lu codes to feed into lu_description function in order to obtain matching descriptions list
            lu_match = []
            sc = arcpy.SearchCursor(lu_ext, "", "", "VALUE", "")
            for i in sc:
                v = i.getValue("VALUE")
                lu_match.append(v)

            # create list of lists with zeroes
            soil_acre_lists = [[0, 0, 0, 0] for i in range(len(lu_match))]

            # preapre a list of soil acreage using lu_match list
            lc_soil_diss = []
            soil_lc_diss = []
            lc_soil_aa = []
            sr = arcpy.SearchCursor(lu_soil_diss, "", "", "GRIDCODE;GRIDCODE_1;area", "")
            for s in sr:
                lc = s.getValue("GRIDCODE")
                lc_soil_diss.append(lc)
                sc = s.getValue("GRIDCODE_1")
                soil_lc_diss.append(sc)
                aa = s.getValue("area")
                lc_soil_aa.append(round(aa, 2))

            for idx, lu in enumerate(lu_match):
                for l, s, a in zip(lc_soil_diss, soil_lc_diss, lc_soil_aa):
                    if l == lu:
                        soil_acre_lists[idx][int(s) - 1] = a

            # prepare matching list of lu description using lu codes from lu raster of watershed
            if landuse == "NLCD 2011":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookuppoor.txt"

            if landuse == "NLCD 2006":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookuppoor.txt"

            if landuse == "NLCD 2001":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookuppoor.txt"

            if landuse == "1997 MOP":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookuppoor.txt"

            if landuse == "2002 MOP":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookuppoor.txt"

            if landuse == "2010 MOP":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookuppoor.txt"

            if landuse == "2002 MD/DE":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/mddelookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/mddelookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/mddelookuppoor.txt"

            if landuse == "Ultimate":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/zoninglookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/zoninglookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/zoninglookuppoor.txt"

            if landuse == "MRLC":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/mrlclookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/mrlclookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/mrlclookuppoor.txt"

            if landuse == "1970s USGS":
                if hyd == "Fair":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/usgslookupfair.txt"
                elif hyd == "Good":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/usgslookupgood.txt"
                elif hyd == "Poor":
                    lut_file = r"" + Directory + "/data/mdfiles/lookup/usgslookuppoor.txt"

            # run hydro function to obtain land use description of categories present in watershed
            lu_desc = hydro.lu_description(lut_file, lu_match)

            #change land use edited nubmers to custom names
            global landedit
            for i,lu in enumerate(lu_desc):
                for row in landedit:
                    if lu == row[0]:
                        lu_desc[i] = row[1]

            # Perform two tasks:
            # a) take "lu_desc" and "soil_acre_lists" and concatenate them
            # b) format according to basin composition file in legacy version
            land_soil_area = ""
            width = 30
            for ld, sg in zip(lu_desc, soil_acre_lists):
                land_soil_area = land_soil_area + "{: <{}}".format(ld, width) + str(sg[0]).rjust(10) + str(sg[1]).rjust(
                    10) + str(sg[2]).rjust(10) + str(sg[3]).rjust(10)
                land_soil_area = land_soil_area + "" "\n"

            # *******************************************************************************************************
            # Text file string variables
            # *******************************************************************************************************
            now = datetime.now()
            month = now.strftime("%B")
            day = now.strftime("%d")
            year = now.strftime("%Y")

            datastring = ""
            datastring = datastring + "GISHydro Release Version Date:    %s" "\n" % (Modifieddt)
            datastring = datastring + "Project Name:                     %s" % (proj)
            datastring = datastring + "" "\n"
            datastring = datastring + "Analysis Date:                    %s %s, %s " "\n" % (month, day, year)
            datastring = datastring + "" "\n"
            datastring = datastring + "" "\n"
            datastring = datastring + "Landuse and Soil Distributions for:" "\n"
            datastring = datastring + "" "\n"
            datastring = datastring + "Distribution of Landuse by Soil Group" "\n"
            datastring = datastring + "" "\n"
            datastring = datastring + "Acres on Indicated Soil Group".rjust(66)
            datastring = datastring + "" "\n"
            datastring = datastring + "Land Use".rjust(8) + "A-Soil".rjust(32) + "B-Soil".rjust(10) + "C-Soil".rjust(
                10) + "D-Soil".rjust(10)
            datastring = datastring + "" "\n"
            datastring = datastring + "" "\n"
            datastring = datastring + land_soil_area
            # datastring = datastring + "{: <{}}".format(ld, width) + str(sg[0]).rjust(10) + str(sg[1]).rjust(10) + str(sg[2]).rjust(10) + str(sg[3]).rjust(10)

            # sum list of lists separately and cat at the end of lu description
            total_area = [sum(i) for i in zip(*soil_acre_lists)]
            datastring = datastring + "{: <{}}".format("Total Area:", width) + str(total_area[0]).rjust(10) + str(
                total_area[1]).rjust(10) + str(total_area[2]).rjust(10) + str(total_area[3]).rjust(10)

            datastring = datastring + "" "\n"
            datastring = datastring + "" "\n"
            datastring = datastring + "Distribution of Land Use and Curve Numbers Used" "\n"
            datastring = datastring + "" "\n"
            datastring = datastring + "Land Use".rjust(8) + "Acres".rjust(34) + "Percent".rjust(10) + "A".rjust(
                4) + "B".rjust(4) + "C".rjust(4) + "D".rjust(4)
            datastring = datastring + "" "\n"
            datastring = datastring + "" "\n"

            # loop over land use, related total acreage, percent of land covered by this lu category, and A-B-C-D curve numbers
            curve_num = []
            for l in lu_match:
                with open(lut_file, "r") as f:
                    next(f)
                    for line in f:
                        luc = line.split("\t")[0]
                        if int(l) == int(luc):
                            temp = []
                            A = line.split("\t")[2]  # CN A
                            temp.append(A)
                            B = line.split("\t")[3]  # CN B
                            temp.append(B)
                            C = line.split("\t")[4]  # CN C
                            temp.append(C)
                            D = line.split("\t")[5]  # CN D
                            temp.append(D)
                            curve_num.append(temp)

            # sum areas for each sub-list individually
            acres = [sum(i) for i in soil_acre_lists]
            total_all = sum(total_area)
            try:
                percent = [round(float(x / total_all) * 100, 2) for x in acres]
            except:
                pythonaddins.MessageBox("Not enough land use cover","Basin Composition Error")
                return

            # add curve number to land edits
            for row in landedit:
                curve_num.append([row[4], row[5], row[6], row[7]])

            for l, a, p, cn in zip(lu_desc, acres, percent, curve_num):
                datastring = datastring + "{: <{}}".format(l, width) + str(a).rjust(12) + str(p).rjust(10) + str(
                    cn[0]).rjust(4) + str(cn[1]).rjust(4) + str(cn[2]).rjust(4) + str(cn[3]).rjust(4)
                datastring = datastring + "" "\n"

            # *******************************************************************************************************
            # turn layers ON/OFF in current data frame
            # *******************************************************************************************************
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "lu_ext":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "soil_ext":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "lu_poly":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "soil_poly":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "lu_soil":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "lu_soil_diss":
                    arcpy.mapping.RemoveLayer(df, lyr)

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()

            # *******************************************************************************************************
            # write strings to basin stat text file.
            # Message box containing datastring as message
            # *******************************************************************************************************
            defFN = optfolder + "/basincomp/basincomp.txt"
            compfile = open(defFN, "w")
            compfile.write(datastring)
            compfile.close()

            # *******************************************************************************************************
            # open "basincomp" file in text editor
            # *******************************************************************************************************
            hydro.openbrowser(defFN)

            # *******************************************************************************************************
            # turn peak discharge OFF
            # *******************************************************************************************************
            button2.enabled = False
            tool3.enabled = False
            button3.enabled = True
            arcpy.env.addOutputsToMap = True
            save(optfolder)

class BasinStatistics(object):
    """Implementation for GISHydroNXT_addin.button3 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            arcpy.env.scratchWorkspace = scratchfolder
            arcpy.env.workspace = optfolder
            button1.enabled = False
            button3.enabled = False
            opath = optfolder

            # *******************************************************************************************************
            # Warning messages
            # *******************************************************************************************************
            Impwarntext = """
                ******************************************************
                    IMPERVIOUS AREA IN WATERSHED EXCEEDS 10%.
                    Calculated discharges from USGS Regression
                    Equations may not be appropriate.
                ******************************************************
                             """
            provwarntext = """
                ******************************************************
                    Watershed is within 5km of physiographic
                    province boundary.  You should consider
                    sensitivity of discharges to region location.
                ******************************************************
                             """
            limewarntext = """
                ******************************************************
                    Watershed is within 1km of underlying limestone
                    geology.  You should consider sensitivity
                    of discharges to percent limestone calculated.
                ******************************************************
                             """

            # *******************************************************************************************************
            # Get outlet coordinates and prepare masked grids for calculations
            # *******************************************************************************************************
            outletcell = optfolder + "/outletcell"
            rast = arcpy.Raster(outletcell)
            cellsize = rast.meanCellWidth
            cellsq = cellsize * cellsize
            global xoutlet
            xoutlet = float(xoutletstring)

            global youtlet
            youtlet = float(youtletstring)

            basingrid = optfolder + "/basingrid"
            dirgrid = arcpy.sa.Times(optfolder + "/flowdir_dem", basingrid)
            elevgrid = arcpy.sa.Times(optfolder + "/dem", basingrid)
            lantype = arcpy.sa.Times(optfolder + "/landuse", basingrid)

            # Get basingrid count [number of pixels]
            shedtab = arcpy.SearchCursor(basingrid, "", "", "Count", "")
            for row in shedtab:
                basinarea = row.getValue("Count")
            del row

            # Set values for land edits (URB,IMP,FC,ST)
            FC_edit = 0
            ST_edit = 0
            IMP_edit = 0
            URB_edit = 0
            global landedit
            if os.path.exists(optfolder + "/aux_folder/landuse"):
                land1 = arcpy.sa.Times(optfolder + "/aux_folder/landuse", basingrid)
                land2 = arcpy.sa.Times(optfolder + "/landuse", basingrid)

                landminus = arcpy.sa.Minus(land2, land1)
                landedited = arcpy.sa.SetNull(landminus, land2, "VALUE = 0")
                lantype = arcpy.sa.SetNull(landminus, land1, "VALUE > 0")
                landedited.save(optfolder + "/aux_folder/landuse_edit")

                landsearch = arcpy.SearchCursor(optfolder + "aux_folder/landuse_edit", "","","VALUE;COUNT","")
                for editrast in landsearch:
                    for row in landedit:
                        if editrast.getValue("VALUE") == row[0]:
                            IMP_edit += float(row[3]*editrast.getValue("COUNT")/basinarea)
                            if row[2] == "FC":
                                FC_edit += float(editrast.getValue("COUNT")/basinarea*100)
                            elif row[2] == "ST":
                                ST_edit += float(editrast.getValue("COUNT")/basinarea*100)
                            elif row[2] == "UrbPct":
                                URB_edit += float(editrast.getValue("COUNT")/basinarea*100)

            # *******************************************************************************************************
            # Compute channel and land slope
            # *******************************************************************************************************
            theslope = hydro.channelslope(dirgrid, elevgrid, opath)  # already converted into feet/mile
            theslope_feet = float(theslope / 5280.0)

            global thechannelslope  # 1/18/2018: global variable created for use in discharge comparison
            thechannelslope = theslope_feet
            maxlength = float(hydro.maxlength)  # already converted into miles for use in channel slope function

            ### SLOPE GRID:
            dlgrid_temp1 = arcpy.sa.Log2(dirgrid)
            dlgrid_temp2 = dlgrid_temp1 % 2
            dlgrid_temp3 = arcpy.sa.Con(dlgrid_temp2 > 0, pow(2, 0.5), 1)
            dlgrid_temp4 = arcpy.sa.Times(dlgrid_temp3, cellsize)
            dlgrid = arcpy.sa.Times(dlgrid_temp4, basingrid)
            dir_shift_1 = arcpy.Shift_management(elevgrid, optfolder + "/aux_folder/dir_shift_1", -cellsize, 0)
            dir_shift_2 = arcpy.Shift_management(elevgrid, optfolder + "/aux_folder/dir_shift_2", -cellsize, cellsize)
            dir_shift_4 = arcpy.Shift_management(elevgrid, optfolder + "/aux_folder/dir_shift_4", 0, cellsize)
            dir_shift_8 = arcpy.Shift_management(elevgrid, optfolder + "/aux_folder/dir_shift_8", cellsize, cellsize)
            dir_shift_16 = arcpy.Shift_management(elevgrid, optfolder + "/aux_folder/dir_shift_16", cellsize, 0)
            dir_shift_32 = arcpy.Shift_management(elevgrid, optfolder + "/aux_folder/dir_shift_32", cellsize, -cellsize)
            dir_shift_64 = arcpy.Shift_management(elevgrid, optfolder + "/aux_folder/dir_shift_64", 0, -cellsize)
            dir_shift_128 = arcpy.Shift_management(elevgrid, optfolder + "/aux_folder/dir_shift_128", -cellsize, -cellsize)

            shift_temp1 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_1), dlgrid)
            shift_temp2 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_2), dlgrid)
            shift_temp3 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_4), dlgrid)
            shift_temp4 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_8), dlgrid)
            shift_temp5 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_16), dlgrid)
            shift_temp6 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_32), dlgrid)
            shift_temp7 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_64), dlgrid)
            shift_temp8 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_128), dlgrid)

            slope = arcpy.sa.Con(dirgrid == 1, shift_temp1,
                                 arcpy.sa.Con(dirgrid == 2, shift_temp2,
                                              arcpy.sa.Con(dirgrid == 4, shift_temp3,
                                                           arcpy.sa.Con(dirgrid == 8, shift_temp4,
                                                                        arcpy.sa.Con(dirgrid == 16, shift_temp5,
                                                                                     arcpy.sa.Con(dirgrid == 32,
                                                                                                  shift_temp6,
                                                                                                  arcpy.sa.Con(
                                                                                                      dirgrid == 64,
                                                                                                      shift_temp7,
                                                                                                      arcpy.sa.Con(
                                                                                                          dirgrid == 128,
                                                                                                          shift_temp8,
                                                                                                          0))))))))
            path = opath + "/slope_calc"
            if not os.path.exists(path):
                os.mkdir(path, 0o755)

            slope.save(optfolder + "/slope_calc/landslope")
            # Slope value was going below 0 which isn't correct. Condition is added to at least have 0.01
            slopegrid = arcpy.sa.Con(optfolder + "/slope_calc/landslope" > 0,
                                     optfolder + "/slope_calc/landslope", 0.01)

            landsloperesult = arcpy.GetRasterProperties_management(slopegrid, "MEAN")
            landslopevalue = float(landsloperesult.getOutput(0))
            global landslope
            landslope = float(landslopevalue) / 3.28084  # modified: 09/06/2018 (divided by 3.28084)

            # *******************************************************************************************************
            # Determine province of outlet location
            # *******************************************************************************************************
            outletcell = arcpy.sa.Int(optfolder + "/outletcell")
            outletcell.save(optfolder + "/outcell")
            arcpy.RasterToPolygon_conversion(optfolder + "/outcell", optfolder + "/outletpoly", "NO_SIMPLIFY", "COUNT")
            prov = r"" + Directory + "/data/maryland/Mdprov.shp"
            outletpoly = optfolder + "/outletpoly.shp"
            outpoint = optfolder + "/outletpoint.shp"
            arcpy.Clip_analysis(prov, outletpoly, outpoint)

            global provstring
            ProvTab = arcpy.SearchCursor(outpoint, "", "", "PROVINCE", "")
            for row in ProvTab:
                if row.getValue("PROVINCE") == "A":
                    provstring = "Appalachian Plateaus and Allegheny Ridges"
                elif row.getValue("PROVINCE") == "B":
                    provstring = "Blue Ridge and Great Valley"
                elif row.getValue("PROVINCE") == "P":
                    provstring = "Piedmont"
                elif row.getValue("PROVINCE") == "W":
                    provstring = "Western Coastal Plain"
                elif row.getValue("PROVINCE") == "E":
                    provstring = "Eastern Coastal Plain"
                else:
                    row.getValue("PROVINCE") == "Unknown"
                    provstring == "Unknown"

            # *******************************************************************************************************
            # Get percent soil types from soil dataset
            # *******************************************************************************************************
            ssurgo = arcpy.sa.Times(optfolder + "/soils", basingrid)
            pct = hydro.SSURGOPct(basinarea, ssurgo)
            pctSoil = map(float, pct)
            pctAsoil = float(pctSoil[0])
            pctBsoil = float(pctSoil[1])
            pctCsoil = float(pctSoil[2])
            pctDsoil = float(pctSoil[3])
            pctWsoil = float(pctSoil[4])

            # *******************************************************************************************************
            # Get LU count for Urban, Nil, Forest, and Storage -- it will be different for
            # MOP, MD/DE, MRLC, and USGS
            # *******************************************************************************************************

            if landuse == "NLCD 2011":
                count = hydro.GetLUCountNLCD(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])
            elif landuse == "NLCD 2006":
                count = hydro.GetLUCountNLCD(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])
            elif landuse == "NLCD 2001":
                count = hydro.GetLUCountNLCD(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])
            elif landuse == "2010 MOP":
                count = hydro.GetLUCountAnderson(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])
            elif landuse == "2002 MOP":
                count = hydro.GetLUCountAnderson(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])
            elif landuse == "1997 USGS":
                count = hydro.GetLUCountAnderson(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])
            elif landuse == "2002 MD/DE":
                count = hydro.GetLUCountMDDE(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])
            elif landuse == "Ultimate":
                count = hydro.GetLUCountUltimate(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])
            elif landuse == "MRLC":
                count = hydro.GetLUCountMRLC(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])
            else:  # landuse == "1970s USGS":
                count = hydro.GetLUCountUSGS(basinarea, lantype)
                LUcount = map(float, count)
                UrbPct = float(LUcount[0])
                FC = float(LUcount[2])
                ST = float(LUcount[3])

            # *******************************************************************************************************
            # Get Impervious count -- count will vary depending upon choice of input Landuse and Hyd condition
            # *******************************************************************************************************
            # For Landuse: "NLCD 2011 Landuse"
            if landuse == "NLCD 2011":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountNLCDFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountNLCDGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountNLCDPoor(lantype)

            # For Landuse: "NLCD 2006 Landuse"
            if landuse == "NLCD 2006":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountNLCDFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountNLCDGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountNLCDPoor(lantype)

            # For Landuse: "NLCD 2001 Landuse"
            if landuse == "NLCD 2001":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountNLCDFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountNLCDGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountNLCDPoor(lantype)

            # For Landuse: "2010 MOP Landuse", "2002 MOP Landuse", and "1997 USGS Landuse"
            if landuse == "2010 MOP":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountAndersonFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountAndersonGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountAndersonPoor(lantype)

            if landuse == "2002 MOP":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountAndersonFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountAndersonGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountAndersonPoor(lantype)

            if landuse == "1997 MOP":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountAndersonFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountAndersonGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountAndersonPoor(lantype)

            # For Landuse: "2002 MD/DE Landuse"
            if landuse == "2002 MD/DE":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountMDDEFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountMDDEGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountMDDEPoor(lantype)

            # For Landuse: "Ultiamte Landuse"
            if landuse == "Ultimate":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountUltimateFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountUltimateGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountUltimatePoor(lantype)

            # For Landuse: "MRLC Landuse"
            if landuse == "MRLC":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountMRLCFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountMRLCGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountMRLCPoor(lantype)

            # For Landuse: "1970s USGS Landuse"
            if landuse == "1970s USGS":
                if hyd == "Fair":
                    Impcount = hydro.GetImpCountUSGSFair(lantype)
                elif hyd == "Good":
                    Impcount = hydro.GetImpCountUSGSGood(lantype)
                elif hyd == "Poor":
                    Impcount = hydro.GetImpCountUSGSPoor(lantype)

            global IA
            IA = float((Impcount / basinarea) * 100)

            IA += IMP_edit
            ST += ST_edit
            FC += FC_edit
            UrbPct += URB_edit

            # *******************************************************************************************************
            # Get Limestone percent count
            # *******************************************************************************************************
            arcpy.env.extent = "MAXOF"
            limestonem = r"" + Directory + "/data/maryland/limestonem.shp"

            LIcnt = 0
            limegrid = arcpy.sa.ExtractByMask("basingrid", limestonem)
            limegrid.save(optfolder + "/limegrid")
            arcpy.BuildRasterAttributeTable_management(optfolder + "/limegrid", "Overwrite")
            with arcpy.da.SearchCursor(optfolder + "/limegrid", "Count") as rows:
                for row in rows:
                    LIcnt += row[0] or 0
                    global LI
                    LI = float((float(LIcnt) / basinarea) * 100)  # 10-23-2013
            LI = float((float(LIcnt) / basinarea) * 100)  # 10-23-2013

            global areami2
            areami2 = float((basinarea * cellsq) / 2588881)  # conversion into sq miles

            # *******************************************************************************************************
            # Get basi relief [it is difference of mean elevation and outlet elevation]
            # *******************************************************************************************************
            elev1 = arcpy.GetRasterProperties_management(elevgrid, "MEAN")
            mean_elev = float(elev1.getOutput(0))

            global basinrelief
            basinrelief = float(mean_elev - outletelev)  # Assuming it is already converted into feets

            # *******************************************************************************************************
            # Average CN number using above outgrid depending upon user choice of HydCon
            # *******************************************************************************************************
            cnGrid = arcpy.sa.Times(optfolder + "/curveNumber", basingrid)
            avgCN_val = arcpy.GetRasterProperties_management(cnGrid, "MEAN")  # figure out outgrid source
            global avgCN
            avgCN = float(avgCN_val.getOutput(0))

            # *******************************************************************************************************
            # percent soil types except SSURGO
            # *******************************************************************************************************
            if soil == "STATSGO":
                tempsoilgrid = arcpy.sa.Times(basingrid, optfolder + "/soilG_a")
                pctAR = arcpy.GetRasterProperties_management(tempsoilgrid, "MEAN")
                pctAsoilR = float(pctAR.getOutput(0))
                tempsoilgrid = arcpy.sa.Times(basingrid, optfolder + "/soilG_b")
                pctBR = arcpy.GetRasterProperties_management(tempsoilgrid, "MEAN")
                pctBsoilR = float(pctBR.getOutput(0))
                tempsoilgrid = arcpy.sa.Times(basingrid, optfolder + "/soilG_c")
                pctCR = arcpy.GetRasterProperties_management(tempsoilgrid, "MEAN")
                pctCsoilR = float(pctCR.getOutput(0))
                tempsoilgrid = arcpy.sa.Times(basingrid, optfolder + "/soilG_d")
                pctDR = arcpy.GetRasterProperties_management(tempsoilgrid, "MEAN")
                pctDsoilR = float(pctDR.getOutput(0))
            else:
                # cell count from SSURGO or RAGAN soil dataset
                # *****************************************************************************

                ssurgan = arcpy.sa.Times(optfolder + "/soils", basingrid)
                pctR = hydro.SoilPct(basinarea, ssurgan)
                ##            pctSoil = map(int,pctR)
                pctSoil = map(float, pctR)
                pctAR = float(pctSoil[0])
                pctAsoilR = float((pctAR / basinarea) * 100)
                pctBR = float(pctSoil[1])
                pctBsoilR = float((pctBR / basinarea) * 100)
                pctCR = float(pctSoil[2])
                pctCsoilR = float((pctCR / basinarea) * 100)
                pctDR = float(pctSoil[3])
                pctDsoilR = float((pctDR / basinarea) * 100)

            """
            The following code calculates the Time of Concentration. If multiple provinces
            are involved, tc is weighted average of area of watershed in each province.
            More correct would be to perform weighted average based on length of channel
            in each province.  This modification will be performed at a later time.  (GEM - 12/01/99)
            """

            # *******************************************************************************************************
            #  Create Zonal stats table based on shedpoly and Mdprov"s Province field
            #  and add fields
            # *******************************************************************************************************
            Mdprov = r"" + Directory + "/data/maryland/Mdprov.shp"
            theVTab = optfolder + "/theVTab.dbf"

            # *******************************************************************************************************
            # don"t add "theVTab" to TOC -- it will change list by drawing order to list
            # by source which will prohibit addition of new layers
            # *******************************************************************************************************
            arcpy.sa.ZonalStatisticsAsTable(Mdprov, "PROVINCE", "basingrid", "theVTab", "DATA", "ALL")
            arcpy.DeleteField_management("theVTab", "ZONE_CODE;MIN;MAX;RANGE;MEAN;STD;SUM;VARIETY;MAJORITY;MINORITY;MEDIAN")
            addFieldNameList = ["Q1.25", "Q1.50", "Q1.75", "Q2", "Q5", "Q10", "Q25", "Q50", "Q100", "Q200", "Q500"]
            for each in addFieldNameList:
                if not len(arcpy.ListFields("theVTab", each)) > 0:
                    arcpy.AddField_management("theVTab", each, "FLOAT", 10, 3)

            # *******************************************************************************************************
            # Create regionlist and regionarea and declare them as global for use in
            #  Thomas Discharge script (is it used in Thomas Discharge script?)
            # *******************************************************************************************************
            sumarea = 0
            theVTab = arcpy.SearchCursor("theVTab", "", "", "Count", "")
            for each in theVTab:
                count = each.getValue("Count")
                sumarea = sumarea + count
            sumArea = sumarea
            del each

            regionlist = []
            regionarea = []
            breakstring = ""
            theVTab = arcpy.SearchCursor("theVTab", "", "", "Province;Count", "")
            for row in theVTab:
                theProv = row.getValue("Province")
                theArea = float(row.getValue("Count"))
                areapercent = float((theArea / sumArea) * 100)
                areapercent = "{0:.2f}".format(areapercent)
                regionlist.append(theProv)
                regionarea.append(areapercent)
                if row.getValue("Province") == "A":
                    breakstring = breakstring + "       -Appalachian Plateaus and Allegheny Ridges %s percent of area" "\n" % (
                        areapercent)
                elif row.getValue("Province") == "B":
                    breakstring = breakstring + "       -Blue Ridge and Great Valley %s percent of area" "\n" % (
                        areapercent)
                elif row.getValue("Province") == "P":
                    breakstring = breakstring + "       -Piedmont %s percent of area" "\n" % (areapercent)
                elif row.getValue("Province") == "W":
                    breakstring = breakstring + "       -Western Coastal Plain %s percent of area" "\n" % (areapercent)
                elif row.getValue("Province") == "E":
                    breakstring = breakstring + "       -Eastern Coastal Plain %s percent of area" "\n" % (areapercent)
            del row

            # *******************************************************************************************************
            # Compute Time of Concentration:
            #                               1]  W.O. Thomas, Jr. Equation   [tc]
            #                               2]  SCS Lag equation * 1.67     [lagtime]
            # *******************************************************************************************************
            sumtc = 0
            theVTab = arcpy.SearchCursor("theVTab", "", "", "Province;Count", "")
            for row in theVTab:
                theProv = row.getValue("Province")
                theArea = row.getValue("Count")

                if row.getValue("Province") == "A":
                    temptc = 0.133 * ((maxlength) ** (0.475)) * ((theslope) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                            (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154)) * ((10) ** (0.194))
                elif row.getValue("Province") == "W":
                    temptc = 0.133 * ((maxlength) ** (0.475)) * ((theslope) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                            (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154)) * ((10) ** (0.366))
                elif row.getValue("Province") == "E":
                    temptc = 0.133 * ((maxlength) ** (0.475)) * ((theslope) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                            (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154)) * ((10) ** (0.366))
                else:
                    temptc = 0.133 * ((maxlength) ** (0.475)) * ((theslope) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                            (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154))
                sumtc = sumtc + (temptc * theArea)
            del row

            global tc
            tc = (sumtc / basinarea)

            # *******************************************************************************************************
            # Calculate lagtime
            # *******************************************************************************************************

            global lagtime
            lagtime = ((np.float64(100 * ((maxlength * 5280) ** (0.8)) * (((1000 / avgCN) - 9) ** (0.7))) / (
                    1900 * ((abs(landslope) * 100) ** (0.5)))) / 60)

            # *******************************************************************************************************
            # Calculate Mean Annual Precipitation
            # *******************************************************************************************************
            arcpy.env.scratchWorkspace = scratchfolder
            arcpy.env.workspace = optfolder
            mapstpm = r"" + Directory + "/data/maryland/map/mapstpm"
            prec_grid = r"" + Directory + "/data/prec/p2-24m"
            arcpy.env.cellSize = str(cellsize)
            basingrid_p = arcpy.Raster(basingrid)
            analysis_extent = basingrid_p.extent
            arcpy.env.extent = analysis_extent
            maprecbasin = arcpy.sa.Times(mapstpm,
                                         basingrid_p)  # Make sure basingrid has value 1 otherwise all precip will be 0
            theprec = arcpy.sa.Times(prec_grid, basingrid_p)
            precavg = arcpy.GetRasterProperties_management(maprecbasin, "MEAN")
            precavg = float(precavg.getOutput(0))
            avgprec = arcpy.GetRasterProperties_management(theprec, "MEAN")
            avgprec = float(avgprec.getOutput(0))

            global maprec
            maprec = float(precavg / (1000 * 2.54))

            global p2yr
            p2yr = float(avgprec / 1000)

            # *******************************************************************************************************
            # format precision before text file string settings
            # *******************************************************************************************************
            maxlength = "{0:.2f}".format(maxlength)
            theslope = "{0:.8f}".format(theslope)
            theslope_feet = "{0:.8f}".format(theslope_feet)
            landslope = "{0:.8f}".format(landslope)
            pctAsoil = "{0:.1f}".format(pctAsoil)
            pctBsoil = "{0:.1f}".format(pctBsoil)
            pctCsoil = "{0:.1f}".format(pctCsoil)
            pctDsoil = "{0:.1f}".format(pctDsoil)
            pctWsoil = "{0:.1f}".format(pctWsoil)
            UrbPct = "{0:.1f}".format(UrbPct)
            FC = "{0:.1f}".format(FC)
            ST = "{0:.1f}".format(ST)
            IA = "{0:.1f}".format(IA)
            LI = "{0:.1f}".format(LI)
            areami2 = "{0:.2f}".format(areami2)
            basinrelief = "{0:.2f}".format(basinrelief)
            avgCN = "{0:.1f}".format(avgCN)
            pctAsoilR = "{0:.1f}".format(pctAsoilR)
            pctBsoilR = "{0:.1f}".format(pctBsoilR)
            pctCsoilR = "{0:.1f}".format(pctCsoilR)
            pctDsoilR = "{0:.1f}".format(pctDsoilR)
            tc = "{0:.2f}".format(tc)
            lagtime = "{0:.2f}".format(lagtime)
            maprec = "{0:.2f}".format(maprec)
            p2yr = "{0:.2f}".format(p2yr)

            # *******************************************************************************************************
            # Basin statistics analysis date
            # *******************************************************************************************************
            now = datetime.now()
            month = now.strftime("%B")
            day = now.strftime("%d")
            year = now.strftime("%Y")

            # *******************************************************************************************************
            # Text file string variables
            # *******************************************************************************************************
            datastring = ""
            datastring = datastring + "GISHydro Release Version Date:    %s" "\n" % (Modifieddt)
            datastring = datastring + "Project Name:                     %s" % (proj)
            datastring = datastring + "" "\n"
            datastring = datastring + "Analysis Date:                    %s %s, %s " "\n" % (month, day, year)
            datastring = datastring + "Data Selected:" "\n"
            datastring = datastring + "     DEM Coverage:                %s" "\n" % (demtot)
            datastring = datastring + "     Land Use Coverage:           %s" "\n" % (landuse)
            datastring = datastring + "     Soil Coverage:               %s" "\n" % (soil)
            datastring = datastring + "     Hydrologic Condition:        %s" "\n" % (hyd)
            datastring = datastring + "     Impose NHD stream Locations: Yes" "\n"
            datastring = datastring + "     Outlet Easting:              %s m (MD Stateplane, NAD 1983)" "\n" % (
                int(xoutlet))
            datastring = datastring + "     Outlet Northing:             %s m (MD Stateplane, NAD 1983)" "\n" % (
                int(youtlet))
            datastring = datastring + "Findings: " "\n"
            datastring = datastring + "     Outlet Location:             %s" "\n" % (provstring)
            datastring = datastring + "     Outlet State:                Maryland" "\n"
            datastring = datastring + "     Drainage Area                %s square miles" "\n" % (areami2)
            datastring = datastring + breakstring
            datastring = datastring + "" "\n"
            datastring = datastring + "     Channel Slope:               %s feet/mile (%s feet/feet)" "\n" % (
                theslope, theslope_feet)
            datastring = datastring + "     Land Slope:                  %s feet/feet" "\n" % (landslope)
            datastring = datastring + "     Urban Area (percent):        %s " "\n" % (UrbPct)
            datastring = datastring + "     Impervious Area (percent):   %s " "\n" % (IA)

            # *******************************************************************************************************
            # Print out Impervious area warning message
            # *** warning message is included in for loop despite the fact that technically it could be printed
            #     twice. Since both "Appalachian Plateau" and "Eastern Coastal Plain" are far apart so it is
            #     impossible to have that big watershed while doing analysis with GISHydroNXT
            # *******************************************************************************************************
            theVTab = arcpy.SearchCursor("theVTab", "", "", "Province", "")
            for row in theVTab:
                if (row.getValue("Province") == "A") or (row.getValue("Province") == "E"):
                    if float(IA) >= 10:
                        pythonaddins.MessageBox("IMPERVIOUS AREA IN WATERSHED EXCEEDS 10%." "\n"
                                                "Calculated discharges from USGS Regression" "\n"
                                                "Equations may not be appropriate.", "Warning", 0)
                        datastring = datastring + Impwarntext
                    else:
                        datastring = datastring + ""
                datastring = datastring + ""

            # *******************************************************************************************************
            # Close to boundary condition for provinces -- Near tool isn"t available with
            # basic level license therefore a more crude method was emplyed here. It can
            # be improved in future by using "arcpy.Geometry()" tool to get distance
            # *******************************************************************************************************
            arcpy.Buffer_analysis(optfolder + "/watershed.shp", optfolder + "/wats_prov.shp", "5000", "#", "#", "ALL",
                                  "FID")
            province = r"" + Directory + "/data/maryland/provlines.shp"
            wats_prov = optfolder + "/wats_prov.shp"
            prov_int = optfolder + "/prov_int"
            arcpy.Intersect_analysis([province, wats_prov], prov_int, "ALL", "#", "INPUT")
            prov_cursor = arcpy.SearchCursor(optfolder + "/prov_int.shp", "", "", "FID", "")
            prov = prov_cursor.next()
            if prov != None:
                datastring = datastring + provwarntext

            # *******************************************************************************************************
            # Close to boundary condition for limestone -- Near tool isn"t available with
            # basic level license therefore a more crude method was emplyed here. It can
            # be improved in future by using "arcpy.Geometry()" tool to get distance
            # *******************************************************************************************************
            arcpy.Buffer_analysis(optfolder + "/watershed.shp", optfolder + "/wats_lime.shp", "1000", "#", "#", "ALL",
                                  "FID")
            wats_lime = optfolder + "/wats_lime.shp"
            lime_int = optfolder + "/lime_int"
            arcpy.Intersect_analysis([limestonem, wats_lime], lime_int, "ALL", "#", "INPUT")
            lime_cursor = arcpy.SearchCursor(optfolder + "/lime_int.shp", "", "", "FID", "")
            lime = lime_cursor.next()
            if lime != None:
                datastring = datastring + limewarntext

            # *******************************************************************************************************
            # Continued -- Text file string variables
            # *******************************************************************************************************
            datastring = datastring + "" "\n"
            datastring = datastring + "     Time of Concentration:       %s hours [W.O. Thomas, Jr. Equation]" "\n" % (tc)
            datastring = datastring + "     Time of Concentration:       %s hours  [From SCS Lag Equation * 1.67]" "\n" % (
                lagtime)
            datastring = datastring + "     Longest Flow Path:           %s miles" "\n" % (maxlength)
            datastring = datastring + "     Basin Relief:                %s feet" "\n" % (basinrelief)
            datastring = datastring + "     Average CN:                  %s" "\n" % (avgCN)
            datastring = datastring + "     Forest Cover (percent):      %s" "\n" % (FC)
            datastring = datastring + "     Storage (percent):           %s" "\n" % (ST)
            datastring = datastring + "     Limestone (percent):         %s" "\n" % (LI)
            datastring = datastring + "     Selected Soils Data Statistics Percent:" "\n"
            datastring = datastring + "        A Soils:                  %s" "\n" % (pctAsoilR)
            datastring = datastring + "        B Soils:                  %s" "\n" % (pctBsoilR)
            datastring = datastring + "        C Soils:                  %s" "\n" % (pctCsoilR)
            datastring = datastring + "        D Soils:                  %s" "\n" % (pctDsoilR)
            datastring = datastring + "     SSURGO Soils Data Statistics Percent (used in Regression Equations):" "\n"
            datastring = datastring + "        A Soils:                  %s" "\n" % (pctAsoil)
            datastring = datastring + "        B Soils:                  %s" "\n" % (pctBsoil)
            datastring = datastring + "        C Soils:                  %s" "\n" % (pctCsoil)
            datastring = datastring + "        D Soils:                  %s" "\n" % (pctDsoil)
            datastring = datastring + "     2-Year,24-hour Prec.:        %s inches" "\n" % (p2yr)
            datastring = datastring + "     Mean Annual Prec.:           %s inches" "\n" % (maprec)

            # *******************************************************************************************************
            # turn layers ON/OFF in current data frame
            # *******************************************************************************************************
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "landslope":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "ssurgo":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "limegrid":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "wats_prov":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "prov_int":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "wats_lime":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "lime_int":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "landuse_edit":
                    arcpy.mapping.RemoveLayer(df, lyr)
            if os.path.exists(optfolder + "/aux_folder/landuse_edit"):
                arcpy.Delete_management(optfolder + "/aux_folder/landuse_edit")

            # Deleting shift files for slope
            for i in range(8):
                if os.path.exists(optfolder + "/aux_folder/dir_shift_" + str(2 ** i)):
                    arcpy.Delete_management(optfolder + "/aux_folder/dir_shift_" + str(2 ** i))

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()
            arcpy.env.extent = "MAXOF"

            # *******************************************************************************************************
            # write strings to basin stat text file.
            # Message box containing datastring as message
            # *******************************************************************************************************
            defFN = optfolder + "/basinstat.txt"
            statfile = open(defFN, "w")
            statfile.write(datastring)
            statfile.close()

            # *******************************************************************************************************
            # open "basinstat" file in text editor
            # *******************************************************************************************************
            hydro.openbrowser(defFN)

            # *******************************************************************************************************
            # turn peak discharge OFF
            # *******************************************************************************************************
            button1.enabled = True
            tool7.enabled = False
            button4.enabled = True
            save(optfolder)

class CalculateAttributes(object):
    """Implementation for GISHydroNXT_addin.button10 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            # For some reason the workspace must be setted again
            arcpy.env.scratchWorkspace = scratchfolder
            arcpy.env.workspace = optfolder
            arcpy.env.snapRaster = optfolder + "/dem"
            # *******************************************************************************************************
            # input data needed for "subrivers.shp" and "subsheds.shp"
            # *******************************************************************************************************
            dem = optfolder + "/dem"
            subrivers = optfolder + "/subrivers.shp"
            slope_stats = optfolder + "/slope_stats.dbf"
            polyras = optfolder + "/polyras"
            elevzones = optfolder + "/elevzones.shp"
            elevmerge = optfolder + "/elevmerge.shp"
            flowdir = optfolder + "/flowdir"
            fc = optfolder + "/subshed.shp"  # shapefile
            subsheds = optfolder + "/subsheds"  # raster data
            outlets = optfolder + "/outlets"
            longfp = optfolder + "/longfp.dbf"
            slope_sheds = optfolder + "/slope_sheds.dbf"
            slopegrid = optfolder + "/slope_calc/landslope"
            cntable = optfolder + "/cntable.dbf"
            lu = optfolder + "/landuse"

            # *******************************************************************************************************
            # add length, and slope fields. Update "slope" field
            # "subrivers.shp" attributes calculation
            # *******************************************************************************************************
            if not len(arcpy.ListFields(subrivers, "Length")) > 0:
                arcpy.AddField_management(subrivers, "Length", "FLOAT", 15, 4)
            if not len(arcpy.ListFields(subrivers, "Slope")) > 0:
                arcpy.AddField_management(subrivers, "Slope", "Double", 15,
                                          9)  # 08/05/2013: "Float" to "Double" and field scale from 4 to 9
            arcpy.CalculateField_management(subrivers, "Length", "!shape.length@feet!",
                                            "PYTHON")  # length converted into feet
            arcpy.PolylineToRaster_conversion(subrivers, "ARCID", polyras, "MAXIMUM_LENGTH", "NONE", 30)
            arcpy.RasterToPolygon_conversion(polyras, elevzones, "NO_SIMPLIFY", "Value")
            arcpy.Dissolve_management(elevzones, elevmerge, "GRIDCODE")  # 08/05/2013: added to handle multiple polygons
            arcpy.sa.ZonalStatisticsAsTable(elevmerge, "FID", dem, slope_stats, "DATA",
                                            "ALL")  # NODATA is changed to DATA # 08/05/2013: elevzones changed to elevmerge

            zonalcur = arcpy.SearchCursor(slope_stats, "", "", "MIN;MAX", "")
            lencursor = arcpy.arcpy.SearchCursor(subrivers, "", "", "Length", "")

            minlst = []
            maxlst = []
            global lenlst
            lenlst = []

            for z in zonalcur:
                min_elev = (z.getValue("MIN"))  # removed multiplication factor of 3.2808 because
                max_elev = (z.getValue("MAX"))  # elevation units are already in feets
                minlst.append(min_elev)
                maxlst.append(max_elev)

            for l in lencursor:
                len_cur = l.getValue("Length")
                lenlst.append(len_cur)

            diff = [maxlstd - minlstd for maxlstd, minlstd in zip(maxlst, minlst)]
            global slopeval
            slopeval = [diffd / lenlstd for diffd, lenlstd in zip(diff, lenlst)]

            global no_subwatersheds
            no_subwatersheds = len(slopeval)

            # slope values should at least be 0.0001
            for index, item in enumerate(slopeval):
                if item == "0":
                    slopeval[index] = "0.0001"

            count = 0
            fc = optfolder + "/subrivers.shp"
            slopecur = arcpy.UpdateCursor(fc)
            for s in slopecur:
                val = slopeval[count]
                s.Slope = val
                slopecur.updateRow(s)
                count = count + 1

            # update "GRID_CODE" using "FROM_NODE" attribute values
            nodes = arcpy.UpdateCursor(subrivers)
            for node in nodes:
                node.setValue("GRID_CODE", node.getValue("FROM_NODE"))
                nodes.updateRow(node)

            del node
            del nodes

            # *******************************************************************************************************
            # Fields added:
            #               1] length
            #               2] slope
            #               3] Tc method
            #               4] Area in square miles
            #               5] CN
            #
            # Update "slope" field "subsheds.shp" attributes calculation
            # *******************************************************************************************************
            arcpy.env.scratchWorkspace = scratchfolder
            arcpy.env.workspace = optfolder

            flds = arcpy.sa.FlowLength(flowdir, "DOWNSTREAM", "#")
            minflds = arcpy.sa.ZonalStatistics(subsheds, "VALUE", flds, "MINIMUM", "DATA")
            fldswo = arcpy.sa.Minus(flds, minflds)
            NotOutlet = arcpy.sa.IsNull(outlets)
            fdrnooutlet = arcpy.sa.Divide(flowdir, NotOutlet)
            fluswb = arcpy.sa.FlowLength(fdrnooutlet, "UPSTREAM", "#")
            flplus = arcpy.sa.Plus(fldswo, fluswb)
            lengthlfp = arcpy.sa.ZonalStatistics(subsheds, "VALUE", flplus, "MAXIMUM", "DATA")
            length_lfp = lengthlfp * 3.28084

            # zonal statistics for "longfp" and "slope_sheds"
            slopegrid = arcpy.sa.Divide(slopegrid,
                                        3.28084)  # 09-07-2018: Landslope value should be factored into meters unit
            arcpy.sa.ZonalStatisticsAsTable(subsheds, "VALUE", length_lfp, longfp, "DATA",
                                            "ALL")  # 07-25-2013: "NODATA" was changed to "DATA"
            arcpy.sa.ZonalStatisticsAsTable(subsheds, "VALUE", slopegrid, slope_sheds, "DATA",
                                            "ALL")  # 07-25-2013: "NODATA" was changed to "DATA"

            # calculate area in square miles
            subshed_poly = optfolder + "/subshed.shp"
            if not len(arcpy.ListFields(subshed_poly, "AreaMi2")) > 0:
                arcpy.AddField_management(subshed_poly, "AreaMi2", "FLOAT", 15, 4)
            arcpy.CalculateField_management(subshed_poly, "AreaMi2", "!shape.area@squaremiles!", "PYTHON")

            # update TcMethod
            subshed_poly = optfolder + "/subshed.shp"
            if not len(arcpy.ListFields(subshed_poly, "TcMethod")) > 0:
                arcpy.AddField_management(subshed_poly, "TcMethod", "TEXT", 15, 4)
            subpoly = arcpy.UpdateCursor(subshed_poly)
            tc_name_lst = []
            tc_name_lst.append(Tc_method)
            tc_records = len(slopeval)  # 09-30-2013: Use slopeval list to get length of records to insert in Tc Method list
            tc_n_list = tc_name_lst * tc_records
            tc_index = 0
            for tc in subpoly:
                Method = tc_n_list[tc_index]
                tc.TcMethod = Method
                subpoly.updateRow(tc)
                tc_index = tc_index + 1

            # update CN
            cngrid = optfolder + "/curvenumber"
            arcpy.sa.ZonalStatisticsAsTable(subsheds, "VALUE", cngrid, cntable, "NODATA", "ALL")

            cntable = optfolder + "/cntable.dbf"
            cn_mean = arcpy.SearchCursor(cntable, "", "", "MEAN", "")

            global cn_list
            cn_list = []
            for cn in cn_mean:
                mean_cn = cn.getValue("MEAN")
                cn_list.append(mean_cn)

            cn_lst = 0
            subshed_poly = optfolder + "/subshed.shp"
            if not len(arcpy.ListFields(subshed_poly, "CurveNum")) > 0:
                arcpy.AddField_management(subshed_poly, "CurveNum", "FLOAT", 15, 4)
            cn_update = arcpy.UpdateCursor(subshed_poly)
            for cnum in cn_update:
                sheds_cn = cn_list[cn_lst]
                cnum.CurveNum = sheds_cn
                cn_update.updateRow(cnum)
                cn_lst = cn_lst + 1

            # update slope -- get mean values from "slope_sheds", add field in "subshed_poly", and update it
            slope_sheds = optfolder + "/slope_sheds.dbf"
            slope_mean = arcpy.SearchCursor(slope_sheds, "", "", "MEAN", "")

            slope_list = []
            for slp in slope_mean:
                mean_slope = slp.getValue("MEAN")
                slope_list.append(mean_slope)

            # assume a minimum slope condition [MDSHA.TimeOfConcentration.CalculateTc -- line 44]
            slope_list = [0.001 if i == 0 else i for i in slope_list]

            # adjust units to "feet/feet"
            slope_list = [i / 3.28084 for i in slope_list]

            # create a list holding subshed (raster) count fields
            subarea = []
            rows = arcpy.SearchCursor(subsheds, "", "", "Count", "")
            for row in rows:
                thecount = row.getValue("Count")
                subarea.append(thecount)

            # *******************************************************************************************************
            # calculate Tc at the end as we need Slope, CN, and LngFlwPth (L) for its computation
            # *******************************************************************************************************
            # SCS method
            if Tc_method == "SCS Lag Formula":
                # update LFP
                arcpy.env.addOutputsToMap = False
                longfp = optfolder + "/longfp.dbf"
                lfp_max = arcpy.SearchCursor(longfp, "", "", "MAX", "")  # already in feets

                lfp_list = []
                for lfp in lfp_max:
                    max_lfp = lfp.getValue("MAX")
                    lfp_list.append(max_lfp)

                lfp_lst = 0
                subshed_poly = optfolder + "/subshed.shp"
                if not len(arcpy.ListFields(subshed_poly, "LngFlwPth")) > 0:
                    arcpy.AddField_management(subshed_poly, "LngFlwPth", "FLOAT", 15, 4)
                lfp_update = arcpy.UpdateCursor(subshed_poly)

                for lonfp in lfp_update:
                    sheds_lfp = lfp_list[lfp_lst]
                    lonfp.LngFlwPth = sheds_lfp
                    lfp_update.updateRow(lonfp)
                    lfp_lst = lfp_lst + 1

                slp_lst = 0
                subshed_poly = optfolder + "/subshed.shp"
                if not len(arcpy.ListFields(subshed_poly, "Slope")) > 0:
                    arcpy.AddField_management(subshed_poly, "Slope", "FLOAT", "#", "#")
                slope_update = arcpy.UpdateCursor(subshed_poly)
                for slope in slope_update:
                    sheds_slope = slope_list[slp_lst]
                    slope.Slope = sheds_slope
                    slope_update.updateRow(slope)
                    slp_lst = slp_lst + 1

                global tc_list
                tc_list = []
                for a, b, c in zip(lfp_list, cn_list, slope_list):
                    tc_val = (np.float64(100 * ((a) ** 0.8) * (((1000 / b) - 9) ** 0.7)) / (
                            19000 * ((c) ** 0.5))) / 60  # convert float to np.float to handle zero division error
                    tc_list.append(tc_val)

                tc_lst = 0
                subshed_poly = optfolder + "/subshed.shp"
                if not len(arcpy.ListFields(subshed_poly, "Tc")) > 0:
                    arcpy.AddField_management(subshed_poly, "Tc", "FLOAT", 15, 4)
                tc_update = arcpy.UpdateCursor(subshed_poly)
                for tc in tc_update:
                    sheds_tc = tc_list[tc_lst]
                    tc.Tc = sheds_tc
                    tc_update.updateRow(tc)
                    tc_lst = tc_lst + 1
                arcpy.env.addOutputsToMap = True
            # Hydrology Panel Method
            elif Tc_method == "Hydrology Panel Tc Method":
                arcpy.env.addOutputsToMap = False
                # get "gc_lst" before extent definition to avoid intersection for only sub-basin
                Mdprov = r"" + Directory + "/data/maryland/Mdprov.shp"
                subshed = optfolder + "/subshed.shp"
                subshed_prov = optfolder + "/subshed_prov.shp"
                arcpy.Intersect_analysis([Mdprov, subshed], subshed_prov, "ALL", "#", "INPUT")
                if not len(arcpy.ListFields(subshed_poly, "Area")) > 0:
                    arcpy.AddField_management(subshed_prov, "Area", "FLOAT", 15, 4)
                arcpy.CalculateField_management(subshed_prov, "Area", "!shape.area@squaremeters!", "PYTHON")

                # clip and save sub-sheds 10% bigger than the original extent
                fc = optfolder + "/subshed.shp"
                desc = arcpy.Describe(fc)
                rows = arcpy.SearchCursor(fc)
                FC = []
                ST = []
                IA = []
                lst = 0
                for idx, row in enumerate(rows):
                    aPoly = row.getValue(desc.shapefieldname)
                    setExtent = aPoly.extent
                    expPercent = 0.1
                    bigExtent = hydro.expandExtent(setExtent, expPercent)
                    arcpy.env.extent = arcpy.Extent(int(bigExtent[0]), int(bigExtent[1]), int(bigExtent[2]),
                                                    int(bigExtent[3]))
                    facc = optfolder + "/flowacc"
                    faccdesc = arcpy.Describe(facc)
                    cellSize = faccdesc.meanCellHeight
                    arcpy.env.cellSize = cellSize

                    # create mask raster for each sub-basin [04-18-2014: Legacy does same masking and only allow calculations for sub-basin extent]
                    basingrid = optfolder + "/basingrid"
                    arcpy.MakeFeatureLayer_management(fc, "layer" + str(idx), ' "FID" = ' + str(idx))
                    mask = arcpy.Clip_management(basingrid, "#", optfolder + "/shedMask" + str(idx), "layer" + str(idx),
                                                 "0", "ClippingGeometry")

                    # clip and save sub-sheds 10% bigger than the original extent
                    arcpy.Clip_management(lu, "%f %f %f %f" % (bigExtent[0], bigExtent[1], bigExtent[2], bigExtent[3]),
                                          optfolder + "/Tc_temp" + str(idx), "#", "#", "NONE")
                    Tc_sub = arcpy.sa.Times(optfolder + "/Tc_temp" + str(idx), mask)
                    Tc_sub.save(optfolder + "/Tc_subshed" + str(idx))

                    # Determine FC, ST, Impervious counts, and save them to separate lists
                    FCcnt = hydro.TcFC(optfolder + "/Tc_subshed" + str(idx))

                    if landuse == "MRLC":
                        STcnt = hydro.TcST_MRLC(optfolder + "/Tc_subshed" + str(idx))
                    else:
                        STcnt = hydro.TcST(optfolder + "/Tc_subshed" + str(idx))

                    if landuse == "MRLC":
                        IMPcnt = hydro.TcImpMRLC(optfolder + "/Tc_subshed" + str(idx))
                    elif landuse == "Ultimate":
                        IMPcnt = hydro.TcImpUltimate(optfolder + "/Tc_subshed" + str(idx))
                    elif landuse == "1997 USGS" or landuse == "1970s USGS":
                        IMPcnt = hydro.TcImpUSGS(optfolder + "/Tc_subshed" + str(idx))
                    else:
                        IMPcnt = hydro.TcImpAnderson(optfolder + "/Tc_subshed" + str(idx))

                    FC.append(FCcnt)
                    ST.append(STcnt)
                    IA.append(IMPcnt)
                    lst = lst + 1

                # compute percent area of FC, ST, and IA in each subwatershed
                FC = [(float(a) / b) * 100 for a, b in zip(FC, subarea)]
                ST = [(float(a) / b) * 100 for a, b in zip(ST, subarea)]
                IA = [(float(a) / b) * 100 for a, b in zip(IA, subarea)]

                # create lists to store gridcode (for duplication of FC, ST, IA, maxlength, and theslope), province, and area
                gc_lst = []
                prov_lst = []
                area_lst = []
                sub_prov = arcpy.SearchCursor(subshed_prov, "", "", "GRIDCODE;PROVINCE;Area", "")
                for sub in sub_prov:
                    gc = sub.getValue("GRIDCODE")
                    gc_lst.append(gc)
                    prov = sub.getValue("PROVINCE")
                    prov_lst.append(prov)
                    area = sub.getValue("Area")
                    area_lst.append(area / 900)  # list with subshed area percent in different provinces

                del sub_prov

                # update length (maxlength) and slope (theslope) fields
                longfp = optfolder + "/longfp.dbf"
                lfp_max = arcpy.SearchCursor(longfp, "", "", "MAX", "")  # already in feets

                lfp_list = []
                for lfp in lfp_max:
                    max_lfp = lfp.getValue("MAX")
                    lfp_list.append(max_lfp)

                lfp_lst = 0
                subshed_poly = optfolder + "/subshed.shp"
                if not len(arcpy.ListFields(subshed_poly, "LngFlwPth")) > 0:
                    arcpy.AddField_management(subshed_poly, "LngFlwPth", "FLOAT", 15, 4)
                lfp_update = arcpy.UpdateCursor(subshed_poly)
                for lonfp in lfp_update:
                    sheds_lfp = lfp_list[lfp_lst]
                    lonfp.LngFlwPth = sheds_lfp
                    lfp_update.updateRow(lonfp)
                    lfp_lst = lfp_lst + 1

                slp_lst = 0
                subshed_poly = optfolder + "/subshed.shp"
                if not len(arcpy.ListFields(subshed_poly, "Slope")) > 0:
                    arcpy.AddField_management(subshed_poly, "Slope", "FLOAT", "#", "#")
                slope_update = arcpy.UpdateCursor(subshed_poly)
                for slope in slope_update:
                    sheds_slope = slope_list[slp_lst]
                    slope.Slope = sheds_slope
                    slope_update.updateRow(slope)
                    slp_lst = slp_lst + 1

                # adjustment to longest flow path and slope for Tc calculation
                maxlength = [float(x) / 5280 for x in lfp_list]  # maxlength = [float(x)/5280 for x in length]
                theslope = [x * 5280 for x in slope_list]  # theslope = [x*5280 for x in slope]

                # update lists to match intersect subshed with poly (for average weighting of Tc)
                FC_updated = []
                ST_updated = []
                IA_updated = []
                maxlength_updated = []
                theslope_updated = []
                tmp_tc = []
                for i, g in enumerate(groupby(gc_lst)):
                    FC_updated += [FC[i]] * len(list(g[1]))

                for i, g in enumerate(groupby(gc_lst)):
                    ST_updated += [ST[i]] * len(list(g[1]))

                for i, g in enumerate(groupby(gc_lst)):
                    IA_updated += [IA[i]] * len(list(g[1]))

                for i, g in enumerate(groupby(gc_lst)):
                    maxlength_updated += [maxlength[i]] * len(list(g[1]))

                for i, g in enumerate(groupby(gc_lst)):
                    theslope_updated += [theslope[i]] * len(list(g[1]))

                # looping over lists to compute "temptc" -- make sure that all lists are of equal length
                for a, b, c, d, e, f, g in zip(prov_lst, maxlength_updated, theslope_updated, FC_updated, IA_updated,
                                               ST_updated, area_lst):
                    if a == "A":
                        temptc = (0.133 * (b ** (0.475)) * (c ** (-0.187)) * ((101 - d) ** (-0.144)) * (
                                (101 - e) ** (0.861)) * ((f + 1) ** (0.154)) * ((10) ** (0.194)))
                    elif a == "W":
                        temptc = (0.133 * (b ** (0.475)) * (c ** (-0.187)) * ((101 - d) ** (-0.144)) * (
                                (101 - e) ** (0.861)) * ((f + 1) ** (0.154)) * ((10) ** (0.366)))
                    elif a == "E":
                        temptc = (0.133 * (b ** (0.475)) * (c ** (-0.187)) * ((101 - d) ** (-0.144)) * (
                                (101 - e) ** (0.861)) * ((f + 1) ** (0.154)) * ((10) ** (0.366)))
                    else:
                        temptc = (0.133 * (b ** (0.475)) * (c ** (-0.187)) * ((101 - d) ** (-0.144)) * (
                                (101 - e) ** (0.861)) * ((f + 1) ** (0.154)))
                    tmp_tc.append(g * temptc)

                tc = []
                for num, grp in groupby(enumerate(gc_lst), itemgetter(1)):
                    tmp_list = [tmp_tc[idx] for idx, _ in grp]
                    tc.append(sum(tmp_list))

                tc_list = [a / b for a, b in zip(tc, subarea)]

                tc_lst = 0
                subshed_poly = optfolder + "/subshed.shp"
                if not len(arcpy.ListFields(subshed_poly, "Tc")) > 0:
                    arcpy.AddField_management(subshed_poly, "Tc", "FLOAT", 15, 4)
                tc_update = arcpy.UpdateCursor(subshed_poly)
                for tc in tc_update:
                    sheds_tc = tc_list[tc_lst]
                    tc.Tc = sheds_tc
                    tc_update.updateRow(tc)
                    tc_lst = tc_lst + 1
                arcpy.env.addOutputsToMap = True
            # Velocity Method
            else:
                # LongPathSub"s schema is locked once added to data frame and "Reset" button fails due to this reason
                # If we don"t add LongPathSub to current data frame then we can remove "vel_meth" folder without any error
                arcpy.env.addOutputsToMap = False
                Tc_method == "Velocity Method Tc Calculation"
                # *******************************************************************************************************
                # Add fields to "subshed.shp" for subsequent processing during segment merge
                # *******************************************************************************************************
                subshed_poly = optfolder + "/subshed.shp"
                arcpy.AddField_management(subshed_poly, "sheet_n", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "sheet_P", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "sheet_L", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "shal_Paved", "TEXT", 15, 4)
                arcpy.AddField_management(subshed_poly, "channel_n", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "ChanDef", "TEXT", 15, 4)
                arcpy.AddField_management(subshed_poly, "ChanSA", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "WidthCoef", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "WidthExp", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "DepthCoef", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "DepthExp", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "XAreaCoef", "FLOAT", 15, 4)
                arcpy.AddField_management(subshed_poly, "XAreaExp", "FLOAT", 15, 4)

                # add "n_sheet"
                sheet_n_lst = []
                sheet_n_lst.append(Tc_ns)
                sheet_n_lst = sheet_n_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for n_sheet in subpoly:
                    relay = sheet_n_lst[index]
                    n_sheet.sheet_n = relay
                    subpoly.updateRow(n_sheet)
                    index = index + 1

                # add "P_sheet"
                sheet_P_lst = []
                sheet_P_lst.append(Tc_P)
                sheet_P_lst = sheet_P_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for P_sheet in subpoly:
                    relay = sheet_P_lst[index]
                    P_sheet.sheet_P = relay
                    subpoly.updateRow(P_sheet)
                    index = index + 1

                # add "L_sheet"
                sheet_L_lst = []
                sheet_L_lst.append(Tc_L)
                sheet_L_lst = sheet_L_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for L_sheet in subpoly:
                    relay = sheet_L_lst[index]
                    L_sheet.sheet_L = relay
                    subpoly.updateRow(L_sheet)
                    index = index + 1

                # add "Paved or Unpaved"
                shallow_Paved_lst = []
                if Tc_paved == True:
                    shallow_Paved_lst.append("Paved")
                    shallow_Paved_lst = shallow_Paved_lst * no_subwatersheds
                else:
                    shallow_Paved_lst.append("Unpaved")
                    shallow_Paved_lst = shallow_Paved_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for Paved in subpoly:
                    relay = shallow_Paved_lst[index]
                    Paved.shal_Paved = relay
                    subpoly.updateRow(Paved)
                    index = index + 1

                # use "NHD" or modified streams
                ChanDef_lst = []
                if Tc_NHD == True:
                    ChanDef_lst.append("NHD")
                    ChanDef_lst = ChanDef_lst * no_subwatersheds
                else:
                    ChanDef_lst.append("SourceArea")
                    ChanDef_lst = ChanDef_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for Def_Chan in subpoly:
                    relay = ChanDef_lst[index]
                    Def_Chan.ChanDef = relay
                    subpoly.updateRow(Def_Chan)
                    index = index + 1

                # add "n_channel"
                channel_n_lst = []
                channel_n_lst.append(Tc_nc)
                channel_n_lst = channel_n_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for n_channel in subpoly:
                    relay = channel_n_lst[index]
                    n_channel.channel_n = relay
                    subpoly.updateRow(n_channel)
                    index = index + 1

                # add "SA_sheet"
                ChanSA_lst = []
                ChanSA_lst.append(Tc_sa)
                ChanSA_lst = ChanSA_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for SA_Chan in subpoly:
                    relay = ChanSA_lst[index]
                    SA_Chan.ChanSA = relay
                    subpoly.updateRow(SA_Chan)
                    index = index + 1

                # add "cw_coef"
                WidthCoef_lst = []
                WidthCoef_lst.append(Tc_cwCoef)
                WidthCoef_lst = WidthCoef_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for cw_coef in subpoly:
                    relay = WidthCoef_lst[index]
                    cw_coef.WidthCoef = relay
                    subpoly.updateRow(cw_coef)
                    index = index + 1

                # add "cw_exp"
                WidthExp_lst = []
                WidthExp_lst.append(Tc_cwExp)
                WidthExp_lst = WidthExp_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for cw_exp in subpoly:
                    relay = WidthExp_lst[index]
                    cw_exp.WidthExp = relay
                    subpoly.updateRow(cw_exp)
                    index = index + 1

                # add "cd_coef"
                DepthCoef_lst = []
                DepthCoef_lst.append(Tc_cdCoef)
                DepthCoef_lst = DepthCoef_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for cd_coef in subpoly:
                    relay = DepthCoef_lst[index]
                    cd_coef.DepthCoef = relay
                    subpoly.updateRow(cd_coef)
                    index = index + 1

                # add "cd_exp"
                DepthExp_lst = []
                DepthExp_lst.append(Tc_cdExp)
                DepthExp_lst = DepthExp_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for cd_exp in subpoly:
                    relay = DepthExp_lst[index]
                    cd_exp.DepthExp = relay
                    subpoly.updateRow(cd_exp)
                    index = index + 1

                # add "ca_coef"
                XAreaCoef_lst = []
                XAreaCoef_lst.append(Tc_caCoef)
                XAreaCoef_lst = XAreaCoef_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for ca_coef in subpoly:
                    relay = XAreaCoef_lst[index]
                    ca_coef.XAreaCoef = relay
                    subpoly.updateRow(ca_coef)
                    index = index + 1

                # add "ca_exp"
                XAreaExp_lst = []
                XAreaExp_lst.append(Tc_caExp)
                XAreaExp_lst = XAreaExp_lst * no_subwatersheds
                subshed_poly = optfolder + "/subshed.shp"
                subpoly = arcpy.UpdateCursor(subshed_poly)
                index = 0
                for ca_exp in subpoly:
                    relay = XAreaExp_lst[index]
                    ca_exp.XAreaExp = relay
                    subpoly.updateRow(ca_exp)
                    index = index + 1

                # Add fields to "subshed.shp"
                longfp = optfolder + "/longfp.dbf"
                lfp_max = arcpy.SearchCursor(longfp, "", "", "MAX", "")  # already in feets
                lfp_list = []
                for lfp in lfp_max:
                    max_lfp = lfp.getValue("MAX")
                    lfp_list.append(max_lfp)

                lfp_lst = 0
                subshed_poly = optfolder + "/subshed.shp"
                arcpy.AddField_management(subshed_poly, "LngFlwPth", "FLOAT", 15, 4)
                lfp_update = arcpy.UpdateCursor(subshed_poly)
                for lonfp in lfp_update:
                    sheds_lfp = lfp_list[lfp_lst]
                    lonfp.LngFlwPth = sheds_lfp
                    lfp_update.updateRow(lonfp)
                    lfp_lst = lfp_lst + 1

                slope_sheds = optfolder + "/slope_sheds.dbf"
                slope_mean = arcpy.SearchCursor(slope_sheds, "", "", "MEAN", "")
                slope_list = []
                for slp in slope_mean:
                    mean_slope = slp.getValue("MEAN")
                    slope_list.append(mean_slope)

                slp_lst = 0
                subshed_poly = optfolder + "/subshed.shp"
                arcpy.AddField_management(subshed_poly, "Slope", "FLOAT", "#", "#")
                slope_update = arcpy.UpdateCursor(subshed_poly)
                for slope in slope_update:
                    sheds_slope = slope_list[slp_lst]
                    slope.Slope = sheds_slope
                    slope_update.updateRow(slope)
                    slp_lst = slp_lst + 1

                # create folder "vel_meth"
                path = optfolder + "/vel_meth"
                os.mkdir(path, 0o755)

                # read sheet global variables [some variables are re-defines just to keep consistency with Avenue code]
                sheet_n = float(Tc_ns)
                sheet_P = float(Tc_P)
                sheet_L = float(Tc_L)

                # read channel global variables
                chan_n = float(Tc_nc)
                thechanSA = float(Tc_sa)
                widthcoef = float(Tc_cwCoef)
                widthexp = float(Tc_cwExp)
                depthcoef = float(Tc_cdCoef)
                depthexp = float(Tc_cdExp)
                xareacoef = float(Tc_caCoef)
                xareaexp = float(Tc_caExp)

                # read data
                dem = optfolder + "/dem"
                dirgrid = optfolder + "/flowdir"
                fc = optfolder + "/subshed.shp"
                desc = arcpy.Describe(fc)
                rows = arcpy.SearchCursor(fc)

                # begin indexing and process data for each sub-basin
                tc_list = []
                for idx, row in enumerate(rows):

                    aPoly = row.getValue(desc.shapefieldname)
                    setExtent = aPoly.extent
                    expPercent = 0.1
                    bigExtent = hydro.expandExtent(setExtent, expPercent)
                    arcpy.env.extent = arcpy.Extent(int(bigExtent[0]), int(bigExtent[1]), int(bigExtent[2]),
                                                    int(bigExtent[3]))
                    facc = optfolder + "/flowacc"
                    faccdesc = arcpy.Describe(facc)
                    cellSize = faccdesc.meanCellHeight
                    arcpy.env.cellSize = cellSize

                    # create mask raster for each sub-basin
                    basingrid = optfolder + "/basingrid"
                    arcpy.MakeFeatureLayer_management(fc, "cliplayer", "FID = " + str(idx))
                    mask = arcpy.Clip_management(basingrid, "#", optfolder + "/vel_meth/mask" + str(idx), "cliplayer", "",
                                                 "ClippingGeometry")
                    # clip and save sub-sheds 10% bigger than the original extent

                    facc = arcpy.sa.Plus(facc, 1)

                    if Tc_NHD == True:
                        # 05-09-2017: temporarily disabling clipping of NHD streams
                        arcpy.PolylineToRaster_conversion(Directory + "/data/maryland/nhd_streamsm.shp", "FID",
                                                          optfolder + "/vel_meth/streamgrid" + str(idx), "MAXIMUM_LENGTH",
                                                          "NONE", "30")

                    if Tc_infStreams == True:
                        srcpixel = (thechanSA * pow(5280, 2)) / pow((3.28084 * cellSize), 2)
                        streamgrid_con = arcpy.sa.Con(facc >= srcpixel, 1, 0)
                        InfStr_null = arcpy.sa.SetNull(streamgrid_con, streamgrid_con, "Value = 0")
                        InfStr_min = arcpy.sa.Minus(InfStr_null, 1)
                        streamgrid_masked = arcpy.sa.Times(InfStr_min, mask)
                        streamgrid_masked.save(optfolder + "/vel_meth/streamgrid" + str(idx))

                    # get upstream and log2 base direction raster for sub-extents
                    dirgrid = optfolder + "/flowdir"
                    upgrid = arcpy.sa.FlowLength(dirgrid, "UPSTREAM", "")
                    upgrid = arcpy.sa.Times(upgrid, 3.28084)
                    upgrid = arcpy.sa.Times(upgrid, mask)
                    dlgrid_temp1 = arcpy.sa.Log2(dirgrid)
                    dlgrid_temp2 = dlgrid_temp1 % 2
                    dlgrid_temp3 = arcpy.sa.Con(dlgrid_temp2 > 0, pow(2, 0.5), 1)
                    dlgrid_temp4 = arcpy.sa.Times(dlgrid_temp3, cellSize)
                    dlgrid = arcpy.sa.Times(dlgrid_temp4, 3.28084)
                    dlgrid = arcpy.sa.Times(dlgrid, mask)
                    upgridp = upgrid + dlgrid

                    # calculate indicator grid ("indic = 3" means everything is swale)
                    indic = arcpy.sa.Times(basingrid, 3)  # create swale raster
                    indic = arcpy.sa.Con(upgridp <= sheet_L, 1, indic)  # substitute sheet flow
                    indic = arcpy.sa.Con((upgridp > sheet_L) & (upgrid < sheet_L), 2,
                                         indic)  # substitute pixels that are part sheet, part swale
                    channel_con = arcpy.sa.IsNull(optfolder + "/vel_meth/streamgrid" + str(idx))
                    indic = arcpy.sa.Con(channel_con == 0, 4, indic)  # substitute channel pixels
                    indic = arcpy.sa.Times(indic, mask)

                    # calculate Swale part of travel time
                    if Tc_paved:
                        swale_coef = 73181.52  # revised from 73440 based on 3600 times the value in Appendix F of TR-55 Manual
                    if Tc_unpaved:
                        swale_coef = 58084.20  # revised from 57600 based on 3600 times the value in Appendix F of TR-55 Manual

                    ### VEL METH *FORTRAN* CODE (2018/8/15)
                    flowdir_sub = arcpy.sa.Times(dirgrid, mask)
                    usflowlength = arcpy.sa.FlowLength(flowdir_sub, "UPSTREAM", "")
                    dsflowlength = arcpy.sa.FlowLength(flowdir_sub, "DOWNSTREAM", "")
                    longestpath = arcpy.sa.Plus(usflowlength, dsflowlength)
                    MaxValR = arcpy.GetRasterProperties_management(longestpath, "MAXIMUM")
                    MaxVal = float(MaxValR.getOutput(0))
                    longestpath = arcpy.sa.SetNull(longestpath, longestpath, "Value < " + str(MaxVal - cellSize / 10))
                    longestpath = arcpy.sa.Divide(longestpath, longestpath)
                    strlnk = arcpy.sa.StreamLink(longestpath, dirgrid)
                    streams_no = arcpy.GetCount_management(strlnk)
                    streams_no = int(streams_no.getOutput(0))
                    while streams_no > 2:
                        attributes = arcpy.SearchCursor(strlnk)
                        values = []
                        valsort = []
                        for i in attributes:
                            values.append(i.getValue("Count"))
                            valsort.append(i.getValue("Count"))
                        valsort.sort()
                        strrepeat = 0
                        for i in range(len(valsort) - 1):
                            if valsort[i] == valsort[i + 1]:
                                strrepeat = valsort[i]
                        if strrepeat > 0:
                            index = values.index(strrepeat)
                        else:
                            index = values.index(min(values))
                        longestpath = arcpy.sa.SetNull(strlnk, 1, "Value = " + str(index + 1))
                        strlnk = arcpy.sa.StreamLink(longestpath, dirgrid)
                        arcpy.BuildRasterAttributeTable_management(strlnk)
                        streams_no = arcpy.GetCount_management(strlnk)
                        streams_no = int(streams_no.getOutput(0))
                    longestpath = arcpy.sa.Times(longestpath, indic)
                    arcpy.CopyRaster_management(longestpath, optfolder + "/vel_meth/LongPathSub" + str(idx), "", "", "", "",
                                                "", "8_BIT_SIGNED")
                    arcpy.BuildRasterAttributeTable_management(optfolder + "/vel_meth/LongPathSub" + str(idx))
                    fcp = optfolder + "/vel_meth/pntvel.shp"
                    arcpy.RasterToPoint_conversion(longestpath, fcp, "VALUE")

                    # Must do the following to solve for floating points math issue (multiply by 1e6 and then divide by 1e6)
                    elevgrid1 = arcpy.sa.Int(dem)
                    elevgrid2 = arcpy.sa.Minus(dem, elevgrid1)
                    elevgrid2 = arcpy.sa.Times(elevgrid2, 1e6)

                    landslope = optfolder + "/slope_calc/landslope"
                    arcpy.sa.ExtractMultiValuesToPoints(fcp,
                                                        [[landslope, "Slope"], [elevgrid1, "Elev1"], [elevgrid2, "Elev2"],
                                                         [facc, "Acc"], [indic, "Type"], [dlgrid, "dL"]], "NONE")
                    rawlist = []
                    fields = ['Slope', 'Elev1', 'Elev2', 'Acc', 'Type', 'dL']
                    with arcpy.da.SearchCursor(fcp, fields) as cursor:
                        for rw in cursor:
                            rawlist.append(rw)

                    sorted_list = sorted(rawlist, key=itemgetter(3))
                    Pixel = list(range(1, len(sorted_list) + 1))
                    Slope = [item[0] / 3.28084 for item in sorted_list]
                    Elev1 = [item[1] for item in sorted_list]
                    Elev2 = [item[2] / 1e6 for item in sorted_list]
                    Elev = [x + y for x, y in zip(Elev1, Elev2)]  # concatenate two from the elev list
                    Drain_Area = [item[3] for item in sorted_list]
                    StrType = [item[4] for item in sorted_list]
                    I_Length = [item[5] for item in sorted_list]

                    ### Fix slope function
                    i = 1
                    while i < len(Elev):
                        i -= 1
                        aux_elev = [Elev[i]]
                        aux_dist = [I_Length[i]]
                        j = i + 1
                        while j < len(Elev) and Elev[i] <= Elev[j] + 0.1:
                            aux_elev.append(Elev[j])
                            aux_dist.append(I_Length[j])
                            j += 1
                        if len(aux_elev) > 1 and j < len(Elev):
                            aux_elev.append(Elev[j])
                        while i <= j:
                            if len(aux_elev) > 1 and i < len(Slope) and i < j:
                                slope_aux = (aux_elev[0] - aux_elev[-1]) / float(sum(aux_dist))
                                if slope_aux <= 0:
                                    slope_aux = 0.1 / float(sum(aux_dist))
                                Slope[i] = slope_aux
                            elif i == len(Elev) - 1:
                                Slope[i] = 0.1 / float(aux_dist[-1])
                            i += 1
                    Slope = [0.000001 if x <= 0 else x for x in Slope]

                    ### Time of concentration calc
                    velmeth_prop = zip(Slope, Elev, StrType, I_Length, Drain_Area)

                    Type = []
                    Mixed = []
                    Tot_Length = []
                    Vel = []
                    I_Time = []
                    Tot_Time = []
                    Width = []
                    Depth = []
                    Xarea = []
                    AvgArea = []
                    totlength = 0
                    tottime = 0
                    const = (pow(cellSize, 2) / 27878400) * (pow(3.28084, 2))
                    for vel_prop in velmeth_prop:
                        slope = vel_prop[0]
                        str_type = vel_prop[2]
                        dlg = vel_prop[3]
                        acc = vel_prop[4]
                        avgarea = acc * const
                        width = -1
                        depth = -1
                        xarea = -1
                        totlength = totlength + dlg
                        oldtotlength = totlength - dlg
                        themixed = ' No'
                        flowtype = 'overland'
                        if str_type == 1:
                            dt = 0.007 * (sheet_n * dlg) ** 0.8 / (sheet_P ** 0.5 * slope ** 0.4)
                            thevel = dlg / (dt * 3600)
                        elif str_type == 2:
                            overdist = sheet_L - oldtotlength
                            if flowtype != 'swale' and overdist > 0:
                                swaledist = totlength - sheet_L
                                overdt = 0.007 * (sheet_n * overdist) ** 0.8 / (sheet_P ** 0.5 * slope ** 0.4)
                                swalevel = (swale_coef * slope ** 0.5) / 3600
                                swaledt = swaledist / (swalevel * 3600)
                                dt = overdt + swaledt
                                thevel = dlg / (dt * 3600.0)
                                themixed = 'Yes'
                            else:
                                flowtype = 'swale'
                                thevel = (swale_coef * slope ** 0.5) / 3600
                                dt = dlg / thevel / 3600
                        elif str_type == 3 and flowtype != 'channel':
                            flowtype = 'swale'
                            thevel = (swale_coef * slope ** 0.5) / 3600
                            dt = dlg / thevel / 3600
                        else:
                            flowtype = 'channel'
                            width = widthcoef * (avgarea) ** widthexp
                            depth = depthcoef * (avgarea) ** depthexp
                            xarea = xareacoef * (avgarea) ** xareaexp
                            hr = xarea / (2.0 * depth + width)
                            thevel = 1.49 / chan_n * hr ** 0.66667 * slope ** 0.5
                            dt = dlg / thevel / 3600
                        tottime = tottime + dt
                        Type.append(flowtype)
                        Mixed.append(themixed)
                        Tot_Length.append(totlength)
                        Vel.append(thevel)
                        I_Time.append(dt)
                        Tot_Time.append(tottime)
                        Width.append(width)
                        Depth.append(depth)
                        Xarea.append(xarea)
                        AvgArea.append(avgarea)

                    Pixel = ['%.0f' % elem for elem in Pixel]
                    Drain_Area = ['%.0f' % elem for elem in Drain_Area]
                    Elev = ['%.1f' % elem for elem in Elev]
                    Slope = ['%.6f' % elem for elem in Slope]
                    Width = ['%.2f' % elem for elem in Width]
                    Depth = ['%.2f' % elem for elem in Depth]
                    Xarea = ['%.2f' % elem for elem in Xarea]
                    I_Length = ['%.1f' % elem for elem in I_Length]
                    Tot_Length = ['%.1f' % elem for elem in Tot_Length]
                    Vel = ['%.2f' % elem for elem in Vel]
                    I_Time = ['%.3f' % elem for elem in I_Time]
                    Tot_Time = ['%.3f' % elem for elem in Tot_Time]

                    velmeth_list = zip(Pixel, Type, Mixed, Drain_Area, Elev, Slope, AvgArea, Width,
                                       Depth, Xarea, I_Length, Tot_Length, Vel, I_Time, Tot_Time)

                    with open(optfolder + "/vel_meth/velmethtable" + str(idx) + ".csv", "wb") as f:
                        w = csv.writer(f)
                        w.writerow(["Pixel", "Type", "Mixed", "Drain_Area", "Elev", "Slope", "Avg_Area", "Width",
                                    "Depth", "Xarea", "I_Length", "Tot_Length", "Vel.", "I_Time", "Tot_Time"])
                        w.writerows(velmeth_list)

                    tc_list.append(Tot_Time[-1])
                del idx
                del row

                # add "tc_list" to "subshed.shp"
                subshed_poly = optfolder + "/subshed.shp"
                arcpy.AddField_management(subshed_poly, "Tc", "FLOAT", 15, 4)
                tc_update = arcpy.UpdateCursor(subshed_poly)
                for tc_lst,tc in enumerate(tc_update):
                    sheds_tc_list = tc_list[tc_lst]
                    sheds_tc = float(sheds_tc_list)  # 5/23/2017: type is adjusted to match above Field type
                    tc.Tc = sheds_tc
                    tc_update.updateRow(tc)
                    tc_lst = tc_lst + 1

                global arcid_list
                global arcid
                subshed = optfolder + "/subshed.shp"
                arcid_list = []
                shedtab = arcpy.SearchCursor(subshed, "", "", "ARCID", "")
                for s in shedtab:
                    arcid = s.getValue("ARCID")
                    arcid_list.append(str(int(arcid)))

                arcpy.env.addOutputsToMap = False
                arcpy.env.extent = "MAXOF"

                # re-reading because above variables lead to shapefile lock and prohibit file deletion
                fc = optfolder + "/subshed.shp"
                subs = arcpy.SearchCursor(fc)
                arcpy.env.qualifiedFieldNames = False  # to have attribute tables with original names
                for s, sub in enumerate(subs):
                    arcpy.env.addOutputsToMap = True
                    arcpy.RasterToPolyline_conversion(optfolder + "/vel_meth/LongPathSub" + str(s),
                                                      optfolder + "/aux_folder/Longest_Path_Sub_" + str(s + 1) + ".shp",
                                                      "ZERO", "0", "NO_SIMPLIFY", "Value")

            # *******************************************************************************************************
            # turn layers ON/OFF in current data frame
            # *******************************************************************************************************
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "subshed_prov":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "polyras":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "ModStreams":
                    arcpy.mapping.RemoveLayer(df, lyr)
                for i in range(0, no_subwatersheds, 1):
                    if lyr.name == "Longest_Path_Sub_" + str(i + 1) + "":
                        lyr.visible = True
                        arcpy.ApplySymbologyFromLayer_management(lyr,
                                                                 r"" + Directory + "/data/mdfiles/legends/longpathsub.lyr")
                if lyr.name == "elevmerge":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "elevzones":
                    arcpy.mapping.RemoveLayer(df, lyr)

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()

            # *******************************************************************************************************
            # clean "vel_meth" directory after indexed looping [delete everything except "velmethtables" and *.exe]
            # *******************************************************************************************************
            if os.path.exists(optfolder + "/vel_meth/pntvel.shp"):
                arcpy.Delete_management(optfolder + "/vel_meth/pntvel.shp")
            for i in range(0, no_subwatersheds, 1):
                arcpy.Delete_management(optfolder + "/vel_meth/streamgrid" + str(i))

            #pythonaddins.MessageBox("Attributes Calculation Complete!", "Calculate Attributes", 0)

            # *******************************************************************************************************
            # turn subwatershed delineation OFF and Stream cross-section or Precipitation Depths tool ON
            # *******************************************************************************************************
            button10.enabled = False

            FromNode = []
            ToNode = []
            subriver = optfolder + "/subrivers.shp"
            sr = arcpy.SearchCursor(subriver, "", "", "GRID_CODE;FROM_NODE;TO_NODE", "")
            for node in sr:
                fn = node.getValue("FROM_NODE")
                tn = node.getValue("TO_NODE")
                FromNode.append(fn)
                ToNode.append(tn)

            # "reach_check_lst" is to check if there any reaches.
            # If there are no reaches then watershed is treated as with single basin.
            reach_check_lst = [x for x in FromNode if x in ToNode]
            if not reach_check_lst:
                tool6.enabled = False  # changed from tool.5 to tool.8
                tool8.enabled = False
                button13.enabled = True
            else:
                tool6.enabled = True  # changed from tool.5 to tool.8
                tool8.enabled = True
                button13.enabled = False
                button11.enabled = True
                button12.enabled = True

            # If "Velocity method" is used then we need combine longest flow path,
            # otherwise grey out button
            if Tc_method == "Velocity Method Tc Calculation":
                button11.enabled = True
            else:
                button11.enabled = False
            arcpy.env.extent = "MAXOF"

        pythonaddins.MessageBox("Attributes Calculation Complete!", "GISHydroNXT")
        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()
        save(optfolder)

class ChooseStormDepths(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "IDF Precipitation Values", size=(490, 450))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)
        Depths_font = wx.StaticText(panel, -1, "All depths are in Inches", (160, 5))
        font_size = wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        Depths_font.SetFont(font_size)

        wx.StaticText(panel, -1, "1-yr", (25, 57))
        wx.StaticText(panel, -1, "6-hour", (80, 30))
        wx.StaticText(panel, -1, "12-hour", (178, 30))
        wx.StaticText(panel, -1, "24-hour", (278, 30))
        wx.StaticText(panel, -1, "48-hour", (377, 30))

        wx.StaticText(panel, -1, "2-yr", (25, 87))
        wx.StaticText(panel, -1, "5-yr", (25, 117))
        wx.StaticText(panel, -1, "10-yr", (25, 147))
        wx.StaticText(panel, -1, "25-yr", (25, 177))
        wx.StaticText(panel, -1, "50-yr", (25, 207))
        wx.StaticText(panel, -1, "100-yr", (25, 237))
        wx.StaticText(panel, -1, "200-yr", (25, 267))

        # create lists to hold default (-999), index (locations where to replace avg preci), and average precip depths
        default = []
        for i in range(0, 36):
            val = "-"
            default.append(val)
        precip_avg = ["%0.2f" %x for x in critavg]
        a, b, c = [default, cb_selected, precip_avg]
        d = dict(zip(b, c))
        depths = [d.get(i, j) for i, j in enumerate(a)]

        depthList = [(90, 55), (190, 55), (290, 55), (390, 55),
                     (90, 85), (190, 85), (290, 85), (390, 85),
                     (90, 115), (190, 115), (290, 115), (390, 115),
                     (90, 145), (190, 145), (290, 145), (390, 145),
                     (90, 175), (190, 175), (290, 175), (390, 175),
                     (90, 205), (190, 205), (290, 205), (390, 205),
                     (90, 235), (190, 235), (290, 235), (390, 235),
                     (90, 265), (190, 265), (290, 265), (390, 265),
                     (90, 295), (190, 295), (290, 295), (390, 295)]

        self.controls = []
        for pos, precip in zip(depthList, depths):
            control = wx.TextCtrl(panel, -1, value=str(precip), pos=pos, size=(60, 25))
            if precip == "-":
                control.Enable(False)
            self.controls.append(control)
        cbXYList = [(70, 60), (170, 60), (270, 60), (370, 60),
                    (70, 90), (170, 90), (270, 90), (370, 90),
                    (70, 120), (170, 120), (270, 120), (370, 120),
                    (70, 150), (170, 150), (270, 150), (370, 150),
                    (70, 180), (170, 180), (270, 180), (370, 180),
                    (70, 210), (170, 210), (270, 210), (370, 210),
                    (70, 240), (170, 240), (270, 240), (370, 240),
                    (70, 270), (170, 270), (270, 270), (370, 270),
                    (70, 300), (170, 300), (270, 300), (370, 300)]

        self.cb_list = []
        for pos, precip in zip(cbXYList, depths):
            cb = wx.CheckBox(panel, -1, "", pos)
            cb.SetValue(False)
            if precip == "-":
                cb.Enable(False)
            else:
                cb.SetValue(True)
            self.cb_list.append(cb)

        self.btnOK = wx.Button(panel, label="OK", pos=(190, 355))
        self.Bind(wx.EVT_BUTTON, self.OnOK, id=self.btnOK.GetId())

        self.Show(True)
        self.Centre(True)
        style = self.GetWindowStyle()
        self.SetWindowStyle(style | wx.STAY_ON_TOP)

    def OnOK(self, event):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        # ******************************************************************************************************
        # Use specified duration and depth for index adjustment
        # ******************************************************************************************************
        all_values = []
        for control in self.controls:
            all_values.append(control.GetValue())

        IDF = [str(item) for item in all_values]

        # user specified (or edited) duration and precipitation depths
        py = hydro.PrcpandYear(self.cb_list, IDF)
        global year_uSpecified
        year_uSpecified = py[0]

        global prec_uSpecified
        prec_uSpecified = py[1]

        self.Show(False)

    def OnClose(self, event):
        self.Show(False)


class CombineLongestFlowPathSegments(object):
    """Implementation for GISHydroNXT_addin.button11 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        SegMergeDialog()


class CompareDischarges(object):
    """Implementation for GISHydroNXT_addin.button5_1 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        # Input variables needed for discharge comparison
        # **************************************************
        BR = basinrelief
        FC = float(hydro.FC)  # LI is already declared as global variable
        ST = float(hydro.ST)
        DA = float(areami2)
        HA = float(hydro.pctAsoil)
        HC = float(hydro.pctCsoil)
        HD = float(hydro.pctDsoil)
        HCD = float(HC + HD)
        SLL = float(landslope)
        SLC = float(thechannelslope)
        pr = float(p2yr)

        # lime area
        # basin area
        watershed = optfolder + "/basingrid"
        shedtab = arcpy.SearchCursor(watershed, "", "", "Count", "")
        for row in shedtab:
            basinarea = row.getValue("Count")

        LIcnt = 0
        with arcpy.da.SearchCursor(optfolder + "/limegrid", "Count") as rows:
            for row in rows:
                LIcnt += row[0] or 0
                lime = float((float(LIcnt) / basinarea) * 100)  # 10-23-2013
        lime = float((float(LIcnt) / basinarea) * 100)  # 10-23-2013

        # IA
        imp_count = hydro.Impcount
        IA = float((imp_count / basinarea) * 100)

        LIME = float(lime)
        RCN = float(avgCN)

        # Create regionlist and regionarea
        sumarea = 0
        theVTab = optfolder + "/theVTab.dbf"
        theVTab = arcpy.SearchCursor("theVTab", "", "", "Count", "")
        for each in theVTab:
            count = each.getValue("Count")
            sumarea = sumarea + count
        sumArea = sumarea

        del each

        regionlist = []
        regionfrac = []
        theVTab = arcpy.SearchCursor("theVTab", "", "", "Province;Count", "")
        for row in theVTab:
            theProv = row.getValue("Province")
            theArea = float(row.getValue("Count"))
            areapercent = float(theArea / sumArea)
            regionlist.append(theProv)
            regionfrac.append(areapercent)

        del row

        # Get Limestone percent count
        arcpy.env.extent = "MAXOF"
        limestonem_old = r"" + Directory + "/data/maryland/oldlimestonem.shp"
        basingrid = optfolder + "/basingrid"

        LIcnt = 0
        limegrid = arcpy.sa.ExtractByMask(basingrid, limestonem_old)
        limegrid.save(optfolder + "/limegrid_old")
        arcpy.BuildRasterAttributeTable_management(optfolder + "/limegrid_old", "Overwrite")

        with arcpy.da.SearchCursor(optfolder + "/limegrid_old", "Count") as rows:
            for row in rows:
                LIcnt += row[0] or 0
        LIME_old = float((float(LIcnt) / basinarea) * 100)

        # *******************************************************************************************************
        # Discharge comparison text file
        # *******************************************************************************************************
        now = datetime.now()
        month = now.strftime("%B")
        day = now.strftime("%d")
        year = now.strftime("%Y")

        datastring = ""
        datastring = datastring + "GISHydro Release Version Date:    %s" "\n" % (Modifieddt)
        datastring = datastring + "Project Name:                     %s" % (proj)
        datastring = datastring + "" "\n"
        datastring = datastring + "Analysis Date:                    %s %s, %s " "\n" % (month, day, year)
        datastring = datastring + "" "\n"
        datastring = datastring + "" "\n"
        datastring = datastring + "Discharge Comparison for:" "\n"
        datastring = datastring + "" "\n"
        w = 16  # width of columns
        head = ["Return_Period", "Carpenter", "Carpenter+1SE", "Dillow", "Dillow+1SE", "Thomas2006", "Thomas+1SE_2006",
                "Thomas_2010", "Thomas+1SE_2010"]

        # year list
        year = ["1.25 year", "1.50 year", "1.75 year", "2 year", "5 year", "10 year", "25 year", "50 year", "100 year",
                "200 year", "500 year"]

        # Concatenate Header in a for loop to account for number of provinces
        for i in range(0, len(regionlist), 1):
            carpenter = hydro.CarpenterQ(regionlist[i], regionfrac[i], DA, FC, HA, HD, pr, SLC, ST)
            carp1 = carpenter[0]
            carp2 = carpenter[1]

            dillow = hydro.DillowQ(regionlist[i], regionfrac[i], DA, FC, BR, LIME, RCN, ST)
            dill1 = dillow[0]
            dill2 = dillow[1]

            fixreg2006 = hydro.FixedRegionQ2006(regionlist[i], regionfrac[i], DA, SLL, LIME_old, FC, IA, HD, BR, HA)
            fr2006_1 = fixreg2006[0]
            fr2006_2 = fixreg2006[1]

            fixreg2010 = hydro.FixedRegionQ2010(regionlist[i], regionfrac[i], DA, LIME, FC, IA, HCD, HA, SLL)
            fr2010_1 = fixreg2010[0]
            fr2010_2 = fixreg2010[1]

            datastring = datastring + "{: <{}}".format(head[0], w) + "{: <{}}".format(head[1], w) + "{: <{}}".format(
                head[2], w) + "{: <{}}".format(head[3], w) \
                         + "{: <{}}".format(head[4], w) + "{: <{}}".format(head[5], w) + "{: <{}}".format(head[6],
                                                                                                          w) + "{: <{}}".format(
                head[7], w) + "{: <{}}".format(head[8], w)
            datastring = datastring + "" "\n"
            for q in range(0, 11, 1):
                datastring = datastring + "{: <{}}".format(year[q], w) + "{: <{}}".format(carp1[q],
                                                                                          w) + "{: <{}}".format(
                    carp2[q], w) + "{: <{}}".format(dill1[q], w) \
                             + "{: <{}}".format(dill2[q], w) + "{: <{}}".format(fr2006_1[q], w) + "{: <{}}".format(
                    fr2006_2[q], w) + "{: <{}}".format(fr2010_1[q], w) \
                             + "{: <{}}".format(fr2010_2[q], w)
                datastring = datastring + "" "\n"

        # *******************************************************************************************************
        # turn layers ON/OFF in current data frame
        # *******************************************************************************************************
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "limegrid_old":
                arcpy.mapping.RemoveLayer(df, lyr)

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()

        # *******************************************************************************************************
        # write strings to discharge comparison text file.
        # *******************************************************************************************************
        defFN = optfolder + "/dischargecomparison.txt"
        statfile = open(defFN, "w")
        statfile.write(datastring)
        statfile.close()

        # *******************************************************************************************************
        # open "dischargecomparison" file in text editor
        # *******************************************************************************************************
        hydro.openbrowser(defFN)

        # *******************************************************************************************************
        # turn subwatershed delineation OFF
        # *******************************************************************************************************
        button5_1.enabled = False
        arcpy.env.extent = "MAXOF"
        save(optfolder)

class ControlBoundaries(object):
    """Implementation for GISHydroNXT_addin.tool3 (tool)"""

    def __init__(self):
        self.enabled = False
        self.cursor = 3
        self.shape = "Line"

    def onLine(self, line_geometry):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        controlbound_dir = optfolder + "/control_boundaries.shp"
        watershed_dir = optfolder + "/watershed.shp"
        watershed2_dir = optfolder + "/watershed2.shp"
        watershed_aux_dir = optfolder + "/watershed_aux.shp"
        controlbound_dir = optfolder + "/control_boundaries.shp"
        watershed_clip_dir = optfolder + "/watershed_clip.shp"
        outlet_ws_dir = optfolder + "/outlet_ws.shp"
        basin_dir = optfolder + "/basingrid"
        flowdir_dir = optfolder + "/flowdir"
        flowdir_cb_dir = optfolder + "/flowdir_cb"
        flowdir_dem_dir = optfolder + "/flowdir_dem"
        arcpy.CopyFeatures_management(line_geometry, controlbound_dir)

        user_prompt = pythonaddins.MessageBox("This tool will cut the watershed, are you sure?", "Control Boundaries", 4)

        if user_prompt == "Yes":
            arcpy.Rename_management(watershed_dir, watershed_aux_dir)

            arcpy.FeatureToPolygon_management([watershed_aux_dir, controlbound_dir], watershed_clip_dir, "",
                                              "NO_ATTRIBUTES", "")
            arcpy.MakeFeatureLayer_management(watershed_clip_dir, "clip_lyr")
            arcpy.SelectLayerByLocation_management("clip_lyr", "intersect", outlet_ws_dir)
            arcpy.CopyFeatures_management("clip_lyr", watershed2_dir)

            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            for lyr in arcpy.mapping.ListLayers(mxd, "", df):
                if lyr.name == "watershed_clip":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "watershed_aux":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "clip_lyr":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "flowdir":
                    arcpy.mapping.RemoveLayer(df, lyr)

            arcpy.Delete_management(watershed_clip_dir, "")
            arcpy.Delete_management(watershed_aux_dir, "")
            arcpy.Delete_management(basin_dir, "")
            arcpy.Delete_management(flowdir_dir, "")

            arcpy.Clip_management(flowdir_dem_dir, "#", flowdir_cb_dir, watershed2_dir, "", "ClippingGeometry")
            wshed = arcpy.sa.Watershed(flowdir_cb_dir, optfolder + "/outletcell", "VALUE")
            shed = arcpy.sa.Con(wshed >= 0, 1, arcpy.sa.IsNull(wshed))
            shed.save(basin_dir)
            arcpy.RasterToPolygon_conversion(basin_dir, watershed_dir, "NO_SIMPLIFY", "VALUE")
            arcpy.Clip_management(flowdir_dem_dir, "#", flowdir_dir, watershed_dir, "", "ClippingGeometry")

            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            for lyr in arcpy.mapping.ListLayers(mxd, "", df):
                if lyr.name == "watershed2":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "flowdir_cb":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "flowdir_cb2":
                    arcpy.mapping.RemoveLayer(df, lyr)
            arcpy.ApplySymbologyFromLayer_management(arcpy.mapping.Layer("watershed"),
                                                     r"" + Directory + "/data/mdfiles/legends/watershed.lyr")
            arcpy.Delete_management(watershed2_dir, "")
            arcpy.Delete_management(flowdir_cb_dir, "")

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "control_boundaries":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "basingrid":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "flowdir":
                lyr.visible = False
        arcpy.Delete_management(controlbound_dir, "")
        save(optfolder)

class ControlPanel(object):
    """Implementation for GISHydroNXT_addin.button14 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        TR20ControlPanel()


class CreateContours(object):
    """Implementation for GISHydroNXT_addin.button16 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        # Create contour dialog for user input on contour interval
        SurfaceContours()


# *******************************************************************************************************
# Interactive and selection based processing using Frames (wxPython)
# *******************************************************************************************************
#                                           +++++++++++++++++++
#                                           +                 +
#                                           + Dialogue Frames +
#                                           +                 +
#                                           +++++++++++++++++++

class DataSelectionFrame(wx.Frame):
    def __init__(self):

        # ******************************************************************************************************
        # Create Directory path and remove unwanted files in mxd location
        # ******************************************************************************************************
        global Directory
        mxdloc = arcpy.mapping.MapDocument("CURRENT")
        mxdloc = mxdloc.filePath
        mxdloc = mxdloc.replace("\\", "/")
        Directory = r"" + mxdloc.replace("mdinterface/GISHydroNXT.mxd", "")

        for filename in glob.glob(r"" + Directory + "/mdinterface/t*"):
            try:
                os.remove(filename)
            except:
                pass

        if not os.path.exists(r"" + Directory + "/data"):
            pythonaddins.MessageBox("Error: GISHydroNXT data folder is missing", "Load Error", 0)
            return

        # *******************************************************************************************************
        # checkout ArcGIS 10.1 extensions (Spatial and 3D Analyst)
        # *******************************************************************************************************
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            pythonaddins.MessageBox("Error: Couldn't get Spatial Analyst extension", "Load Error", 0)
            return
        if arcpy.CheckExtension("3D") == "Available":
            arcpy.CheckOutExtension("3D")
        else:
            pythonaddins.MessageBox("Error: Couldn't get 3D Analyst extension", "Load Error", 0)
            return

        ### REMOVE OLD FOLDERS

        #now = time.time()
        #old1 = now - 3600*24*30*6
        #old2 = now - 3600*24*30

        #dir_aux1 = Directory.replace("umdgism","temp")
        #dir_aux2 = Directory.replace("umdgism","intermediate")
        #for root, dirs, files in os.walk(dir_aux1, topdown=False):
        #    for _dir in dirs:
        #        if time.ctime(os.path.getmtime(_dir)) < old1:
        #            try:
        #                shutil.rmtree(_dir)
        #            except:
        #                pass
        #        if all[time.ctime(os.path.getmtime(_dir)) < old2, "_My_Project" in _dir[-15:]]:
        #            try:
        #                shutil.rmtree(_dir)
        #            except:
        #                pass

        ########## DATA INPUT FRAME

        wx.Frame.__init__(self, None, -1, "Input Data Selection", size=(390, 375))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)

        wx.StaticBox(panel, -1, "Project Name", (15, 15), size=(335, 60))
        self.Proj = wx.TextCtrl(panel, -1, value="My Project", pos=(40, 40), size=(285, 20))

        wx.StaticBox(panel, -1, "Data Type Selection", (15, 85), size=(335, 150))
        wx.StaticText(panel, -1, "Soil Data:", (40, 110))
        wx.StaticText(panel, -1, "Landuse Data:", (40, 140))
        wx.StaticText(panel, -1, "Hydrologic Condition:", (40, 170))
        wx.StaticText(panel, -1, "Digital Elevation Model:", (40, 200))
        wx.StaticBox(panel, -1, "Flow Properties", (15, 240), size=(240, 80))
        wx.StaticText(panel, -1, "Accumulation Threshold:", (40, 265))
        wx.StaticText(panel, -1, "Burn Streams:", (40, 295))

        SoilList = ["Ragan", "SSURGO 2000's", "SSURGO 201805"]  # JULY 2019 STATSGO REMOVED (HYDROLOGY PANEL MEETING)
        self.Soil = wx.ComboBox(panel, -1, "SSURGO 201805", (190, 105), (135, -1), SoilList, wx.CB_DROPDOWN)
        LanduseList = ["NLCD 2011", "NLCD 2006", "NLCD 2001", "2010 MOP",
                       "2002 MOP",
                       "2002 MD/DE", "Ultimate", "MRLC", "1997 MOP",
                       "1970s USGS"]
        self.Landuse = wx.ComboBox(panel, -1, "2010 MOP", (190, 135), (135, -1), LanduseList, wx.CB_DROPDOWN)
        HydList = [ "Good","Fair","Poor"]
        self.Hyd = wx.ComboBox(panel, -1, "Good", (190, 165), (135, -1), HydList, wx.CB_DROPDOWN)
        DEMList = ["NED DEM 201805"]  # HERE IT IS POSSIBLE TO ADD MORE DEM'S TO THE LIST
        self.DEMTOT = wx.ComboBox(panel, -1, "NED DEM 201805", (190, 195), (135, -1), DEMList,
                                  wx.CB_DROPDOWN)
        self.FA_Threshold = wx.TextCtrl(panel, -1, value="250", pos=(195, 260), size=(40, -1), style=wx.TE_RIGHT)
        self.burn = wx.CheckBox(panel, -1, "", pos=(215, 295))
        self.burn.SetValue(True)

        self.btnApply = wx.Button(panel, label="Apply", pos=(270, 247), size=(80, 73))
        self.Bind(wx.EVT_BUTTON, self.OnSet, id=self.btnApply.GetId())

        self.Show(True)
        self.Centre(True)
        style = self.GetWindowStyle()
        self.SetWindowStyle(style | wx.STAY_ON_TOP)

    def OnClose(self, event):
        self.Show(False)

    def OnSet(self, event):

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            self.Show(False)
            # *******************************************************************************************************
            # Generate warning message if no data selected or project name specified
            # *******************************************************************************************************
            if len(str(self.Proj.GetValue())) > 50:  # NXT does not handle large strings well)
                pythonaddins.MessageBox("Project Name Exceeds Maximum Number of Characters (50)", "Project Name Warning", 0)
                DataSelectionFrame()
                return

            elif (str(self.Soil.GetValue()) != "Ragan" and
                  str(self.Soil.GetValue()) != "SSURGO 2000's" and
                  str(self.Soil.GetValue()) != "SSURGO 201805" and
                  str(self.Soil.GetValue()) != "STATSGO"):
                pythonaddins.MessageBox("Please Select Soil Data", "Data Selection Warning", 0)
                DataSelectionFrame()
                return

            elif (str(self.Landuse.GetValue()) != "NLCD 2011" and
                  str(self.Landuse.GetValue()) != "NLCD 2006" and
                  str(self.Landuse.GetValue()) != "NLCD 2001" and
                  str(self.Landuse.GetValue()) != "2010 MOP" and
                  str(self.Landuse.GetValue()) != "2002 MOP" and
                  str(self.Landuse.GetValue()) != "2002 MD/DE" and
                  str(self.Landuse.GetValue()) != "Ultimate" and
                  str(self.Landuse.GetValue()) != "MRLC" and
                  str(self.Landuse.GetValue()) != "1997 MOP" and
                  str(self.Landuse.GetValue()) != "1970s USGS"):
                pythonaddins.MessageBox("Please Select Landuse Data", "Data Selection Warning", 0)
                DataSelectionFrame()
                return
            elif (str(self.Hyd.GetValue()) != "Fair" and
                  str(self.Hyd.GetValue()) != "Good" and
                  str(self.Hyd.GetValue()) != "Poor"):
                pythonaddins.MessageBox("Please Select Hydrological Condition", "Data Selection Warning", 0)
                DataSelectionFrame()
                return
            elif (int(self.FA_Threshold.GetValue()) <= 15):
                pythonaddins.MessageBox("Flow Accumulation Threshold Too Low", "Data Selection Warning", 0)
                DataSelectionFrame()
                return
            elif (str(self.Proj.GetValue()) == ""):
                pythonaddins.MessageBox("Please Specify Project Name", "Project Name Warning", 0)
                DataSelectionFrame()
                return
            elif (str(self.DEMTOT.GetValue()) != "NED DEM 201805"):
                pythonaddins.MessageBox("Please Select Elevation Dataset", "DEM Warning", 0)
                DataSelectionFrame()
                return
            else:

                # *******************************************************************************************************
                # Generate optfolder and clip data based on extent rectangle and choice of data selection
                # *******************************************************************************************************

                global soil
                soil = str(self.Soil.GetValue())
                global landuse
                landuse = str(self.Landuse.GetValue())
                global hyd
                hyd = str(self.Hyd.GetValue())
                global proj
                proj = str(self.Proj.GetValue())
                global demtot
                demtot = self.DEMTOT.GetValue()
                global fa_thres
                fa_thres = int(self.FA_Threshold.GetValue())
                global landedit
                landedit = []

                global burn
                burn = 0
                if self.burn.GetValue():
                    burn = 100

                # *******************************************************************************************************
                # Create temp and intermediate folders (previous to project naming)
                # *******************************************************************************************************
                global optfolder
                global scratchfolder
                global timestr
                proj_ = proj.replace(" ", "_")
                timestr = time.strftime("%Y%m%d_%H%M%S")
                data = Directory.replace("/umdgism/", "/") + "/temp/"
                optfolder = r"" + data + timestr + "_" + proj_ + "/"
                data = Directory.replace("/umdgism/", "/") + "/intermediate/"
                scratchfolder = r"" + data + timestr + "_" + proj_ + "/"

                # *******************************************************************************************************
                # environmental settings
                # *******************************************************************************************************

                os.makedirs(optfolder)  # Placed here instead of AoI class to avoid error if tool was closed without selection.
                os.makedirs(scratchfolder)  # Scratch folder will be created once data selection is made.
                os.mkdir(optfolder + "/aux_folder")

                try:
                    if os.path.exists(optfolder + "/aux_folder/user.txt"):
                        os.remove(optfolder + "/aux_folder/user.txt")
                    saveas = open(optfolder + "/aux_folder/user.txt", "w")
                    saveas.write(str(getpass.getuser()))
                    saveas.close()
                except:
                    pass

                arcpy.env.overwriteOutput = True
                arcpy.env.scratchWorkspace = scratchfolder
                arcpy.env.workspace = optfolder

                ######## AREA OF INTEREST
                if demtot == "NED DEM 201805":
                    demname = "neddem"
                # elif demtot == "INSERT DEM NAME"       ## ADD HERE OTHER DEM's

                arcpy.env.snapRaster = r"" + Directory + "/data/dems/neddem"
                array = arcpy.Array()
                array.add(arcpy.Point(extent.XMin, extent.YMin))
                array.add(arcpy.Point(extent.XMin, extent.YMax))
                array.add(arcpy.Point(extent.XMax, extent.YMax))
                array.add(arcpy.Point(extent.XMax, extent.YMin))
                array.add(arcpy.Point(extent.XMin, extent.YMin))
                mask = arcpy.Polygon(array)
                arcpy.CopyFeatures_management(mask, optfolder + "/mask_aux.shp")

                dem = arcpy.Raster(r"" + Directory + "/data/dems/neddem")
                arcpy.PolygonToRaster_conversion(optfolder + "/mask_aux.shp", "FID", optfolder + "/mask_dem", "#", "#", dem)
                arcpy.RasterToPolygon_conversion(optfolder + "/mask_dem", optfolder + "/mask.shp", "NO_SIMPLIFY")

                arcpy.Clip_management(r"" + Directory + "/data/dems/" + demname, "#", optfolder + "/dem",
                                      optfolder + "/mask.shp", "#", "NONE")

                scs = r"" + Directory + "/data/ssurgo/ssurgo_old"
                arcpy.Clip_management(scs, "#", optfolder + "/ssurgo", optfolder + "/mask.shp", "#", "ClippingGeometry","MAINTAIN_EXTENT")

                MOP = """

                The LandUse/Land Cover data set you have selected has been
                provided courtesy of the Maryland Department of Planning.
                Any use of that data set outside of this application without
                the permission of the Department of Planning is prohibited.
                The 2010 data are based on superior imagery and a refined
                classification system. The 2002 and earlier Land Use/Land
                Cover datasets are not reconciled with these improvements.
                For more information on Department of Planning data, please
                visit the MDP web site http://www.mdp.state.md.us or
                call (410) 767-4500.

                """

                global landuse_global
                if landuse == "2010 MOP":
                    pythonaddins.MessageBox(MOP, "Maryland Department of Planning Data Usage Agreement", 0)
                if landuse == "2002 MOP":
                    pythonaddins.MessageBox(MOP, "Maryland Department of Planning Data Usage Agreement", 0)
                if landuse == "2002 MOP":
                    lan = r"" + Directory + "/data/Landuse/mdplu2002"
                elif landuse == "2002 MD/DE":
                    lan = r"" + Directory + "/data/Landuse/mdde2002"
                elif landuse == "2010 MOP":
                    lan = r"" + Directory + "/data/Landuse/lu2010"
                elif landuse == "Ultimate":
                    lan = r"" + Directory + "/data/Landuse/luult"
                elif landuse == "MRLC":
                    lan = r"" + Directory + "/data/Landuse/mrlc"
                elif landuse == "1970s USGS":
                    lan = r"" + Directory + "/data/Landuse/lu70"
                elif landuse == "1997 MOP":
                    lan = r"" + Directory + "/data/Landuse/lu97m"
                elif landuse == "NLCD 2011":
                    lan = r"" + Directory + "/data/Landuse/nlcd2011"
                elif landuse == "NLCD 2006":
                    lan = r"" + Directory + "/data/Landuse/nlcd2006"
                elif landuse == "NLCD 2001":
                    lan = r"" + Directory + "/data/Landuse/nlcd2001"

                arcpy.Clip_management(lan, "#", optfolder + "/landuse", optfolder + "/mask.shp", "#", "NONE","MAINTAIN_EXTENT")
                landuse_global = landuse
                landuse_null = arcpy.sa.IsNull(optfolder + "/landuse")
                rows = arcpy.SearchCursor(landuse_null, "", "", "Value;Count", "")
                NoDataCount = 0
                for row in rows:
                    if int(row.getValue("Value")) == 1:
                        NoDataCount = int(row.getValue("Count"))
                columncount = arcpy.GetRasterProperties_management(landuse_null, "COLUMNCOUNT")
                rowcount = arcpy.GetRasterProperties_management(landuse_null, "ROWCOUNT")
                nodataperc = NoDataCount / float(int(columncount.getOutput(0)) * int(rowcount.getOutput(0))) * 100
                if nodataperc > 5:
                    pythonaddins.MessageBox("Land use does not cover all of the AOI extent.", "Warning")


                # *******************************************************************************************************
                # Select, clip, and generate CN using soil and landuse rasters
                # *******************************************************************************************************
                # Ragan data
                if soil == "Ragan":
                    scs = r"" + Directory + "/data/ragan/ragan"
                    arcpy.Clip_management(scs, "#", optfolder + "/Soils", optfolder + "/mask.shp", "#", "ClippingGeometry","MAINTAIN_EXTENT")

                    # Landuse and hydrologic condition for CN calculation
                    if (landuse == "NLCD 2011" and hyd == "Fair"):
                        outGrid = hydro.nlcdlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2011" and hyd == "Good"):
                        outGrid = hydro.nlcdlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2011" and hyd == "Poor"):
                        outGrid = hydro.nlcdlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "NLCD 2006" and hyd == "Fair"):
                        outGrid = hydro.nlcdlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2006" and hyd == "Good"):
                        outGrid = hydro.nlcdlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2006" and hyd == "Poor"):
                        outGrid = hydro.nlcdlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "NLCD 2001" and hyd == "Fair"):
                        outGrid = hydro.nlcdlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2001" and hyd == "Good"):
                        outGrid = hydro.nlcdlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2001" and hyd == "Poor"):
                        outGrid = hydro.nlcdlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "2010 MOP" and hyd == "Fair"):
                        outGrid = hydro.andlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2010 MOP" and hyd == "Good"):
                        outGrid = hydro.andlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2010 MOP" and hyd == "Poor"):
                        outGrid = hydro.andlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "2002 MOP" and hyd == "Fair"):
                        outGrid = hydro.andlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MOP" and hyd == "Good"):
                        outGrid = hydro.andlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MOP" and hyd == "Poor"):
                        outGrid = hydro.andlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "1997 MOP" and hyd == "Fair"):
                        outGrid = hydro.andlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1997 MOP" and hyd == "Good"):
                        outGrid = hydro.andlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1997 MOP" and hyd == "Poor"):
                        outGrid = hydro.andlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "2002 MD/DE" and hyd == "Fair"):
                        outGrid = hydro.mddelookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MD/DE" and hyd == "Good"):
                        outGrid = hydro.mddelookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MD/DE" and hyd == "Poor"):
                        outGrid = hydro.mddelookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "Ultimate" and hyd == "Fair"):
                        outGrid = hydro.zoninglookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "Ultimate" and hyd == "Good"):
                        outGrid = hydro.zoninglookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "Ultimate" and hyd == "Poor"):
                        outGrid = hydro.zoninglookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "MRLC" and hyd == "Fair"):
                        outGrid = hydro.mrlclookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "MRLC" and hyd == "Good"):
                        outGrid = hydro.mrlclookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "MRLC" and hyd == "Poor"):
                        outGrid = hydro.mrlclookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "1970s USGS" and hyd == "Fair"):
                        outGrid = hydro.usgslookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1970s USGS" and hyd == "Good"):
                        outGrid = hydro.usgslookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1970s USGS" and hyd == "Poor"):
                        outGrid = hydro.usgslookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                # SSURGO data
                elif soil == "SSURGO 2000's" or soil == "SSURGO 201805":
                    scs = r"" + Directory + "/data/ssurgo/ssurgo_old"

                    if soil == "SSURGO 201805":
                        scs = r"" + Directory + "/data/ssurgo/ssurgo_2018"


                    arcpy.Clip_management(scs, "#", optfolder + "/Soils", optfolder + "/mask.shp", "#", "ClippingGeometry","MAINTAIN_EXTENT")

                    # Landuse and hydrologic condition for CN calculation
                    if (landuse == "NLCD 2011" and hyd == "Fair"):
                        outGrid = hydro.nlcdlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2011" and hyd == "Good"):
                        outGrid = hydro.nlcdlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2011" and hyd == "Poor"):
                        outGrid = hydro.nlcdlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "NLCD 2006" and hyd == "Fair"):
                        outGrid = hydro.nlcdlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2006" and hyd == "Good"):
                        outGrid = hydro.nlcdlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2006" and hyd == "Poor"):
                        outGrid = hydro.nlcdlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "NLCD 2001" and hyd == "Fair"):
                        outGrid = hydro.nlcdlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2001" and hyd == "Good"):
                        outGrid = hydro.nlcdlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2001" and hyd == "Poor"):
                        outGrid = hydro.nlcdlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "2010 MOP" and hyd == "Fair"):
                        outGrid = hydro.andlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2010 MOP" and hyd == "Good"):
                        outGrid = hydro.andlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2010 MOP" and hyd == "Poor"):
                        outGrid = hydro.andlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "2002 MOP" and hyd == "Fair"):
                        outGrid = hydro.andlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MOP" and hyd == "Good"):
                        outGrid = hydro.andlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MOP" and hyd == "Poor"):
                        outGrid = hydro.andlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "1997 MOP" and hyd == "Fair"):
                        outGrid = hydro.andlookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1997 MOP" and hyd == "Good"):
                        outGrid = hydro.andlookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1997 MOP" and hyd == "Poor"):
                        outGrid = hydro.andlookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "2002 MD/DE" and hyd == "Fair"):
                        outGrid = hydro.mddelookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MD/DE" and hyd == "Good"):
                        outGrid = hydro.mddelookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MD/DE" and hyd == "Poor"):
                        outGrid = hydro.mddelookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "Ultimate" and hyd == "Fair"):
                        outGrid = hydro.zoninglookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "Ultimate" and hyd == "Good"):
                        outGrid = hydro.zoninglookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "Ultimate" and hyd == "Poor"):
                        outGrid = hydro.zoninglookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "MRLC" and hyd == "Fair"):
                        outGrid = hydro.mrlclookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "MRLC" and hyd == "Good"):
                        outGrid = hydro.mrlclookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "MRLC" and hyd == "Poor"):
                        outGrid = hydro.mrlclookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "1970s USGS" and hyd == "Fair"):
                        outGrid = hydro.usgslookupfair(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1970s USGS" and hyd == "Good"):
                        outGrid = hydro.usgslookupgood(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1970s USGS" and hyd == "Poor"):
                        outGrid = hydro.usgslookuppoor(optfolder + "/landuse", optfolder + "/Soils")
                        outGrid.save(optfolder + "/curveNumber")

                # STATSGO data
                else:
                    scs = r"" + Directory + "/data/maryland/statsgo_allm.shp"
                    output = optfolder + "/Soils.shp"
                    arcpy.Clip_analysis(scs, optfolder + "/mask.shp", output)
                    if (landuse == "NLCD 2011" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgonlcdlookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2011" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgonlcdlookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2011" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgonlcdlookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "NLCD 2006" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgonlcdlookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2006" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgonlcdlookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2006" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgonlcdlookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "NLCD 2001" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgonlcdlookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2001" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgonlcdlookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "NLCD 2001" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgonlcdlookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "2010 MOP" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgoandlookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2010 MOP" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgoandlookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2010 MOP" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgoandlookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "2002 MOP" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgoandlookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MOP" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgoandlookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MOP" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgoandlookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "1997 MOP" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgoandlookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1997 MOP" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgoandlookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1997 MOP" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgoandlookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "2002 MD/DE" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgomddelookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MD/DE" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgomddelookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "2002 MD/DE" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgomddelookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "Ultimate" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgozoninglookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "Ultimate" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgozoninglookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "Ultimate" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgozoninglookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "MRLC" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgomrlclookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "MRLC" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgomrlclookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "MRLC" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgomrlclookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

                    if (landuse == "1970s USGS" and hyd == "Fair"):
                        opath = optfolder
                        outGrid = hydro.statsgousgslookupfair(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1970s USGS" and hyd == "Good"):
                        opath = optfolder
                        outGrid = hydro.statsgousgslookupgood(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")
                    elif (landuse == "1970s USGS" and hyd == "Poor"):
                        opath = optfolder
                        outGrid = hydro.statsgousgslookuppoor(optfolder + "/landuse", optfolder + "/Soils.shp", opath)
                        outGrid.save(optfolder + "/curveNumber")

            # **************************************************************************
            # Flow accumulation and direction OTF
            # **************************************************************************
            arcpy.env.snapRaster = optfolder + "/dem"
            arcpy.env.scratchWorkspace = scratchfolder
            arcpy.env.workspace = optfolder

            arcpy.Clip_analysis(r"" + Directory + "/data/maryland/nhd_streamsm.shp", optfolder + "/mask.shp",
                                optfolder + "/nhdstreams.shp")
            arcpy.PolylineToRaster_conversion(optfolder + "/nhdstreams.shp", "FID", optfolder + "/nhd_rast", "", "",
                                              optfolder + "/dem")

            calc1 = arcpy.sa.IsNull(optfolder + "/nhd_rast")
            calc2 = arcpy.sa.Minus(calc1, 1)
            calc1 = arcpy.sa.Times(calc2, -1 * burn)
            calc2 = arcpy.sa.Minus(optfolder + "/dem", calc1)
            burned = arcpy.sa.Plus(calc2, burn)
            burned.save(optfolder + "/burned")

            fill = arcpy.sa.Fill(optfolder + "/burned")
            fill.save(optfolder + "/fill")
            flowdir = arcpy.sa.FlowDirection(optfolder + "/fill")
            flowdir.save(optfolder + "/flowdir")
            flowacc = arcpy.sa.FlowAccumulation(optfolder + "/flowdir")
            flowacc.save(optfolder + "/flowacc")
            infstr = arcpy.sa.Con(flowacc >= fa_thres, flowacc)
            infstr.save(optfolder + "/InfStreams")

            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            infstr = arcpy.mapping.Layer(optfolder + "/InfStreams")
            arcpy.mapping.AddLayer(df, infstr, "AUTO_ARRANGE")

            # Flowdir backup for control boundaries
            arcpy.CopyRaster_management(optfolder + "/flowdir", optfolder + "/flowdir_dem")

            # *******************************************************************************************************
            # create "AddasStreams.shp" for to be used in subwatershed delineation and
            # create "AddasOutlets.shp" to add user specified outlets
            # *******************************************************************************************************
            spatial_reference = arcpy.Describe(optfolder + "/dem").spatialReference
            arcpy.CreateFeatureclass_management(optfolder, "AddasStreams.shp", "POLYLINE", "", "ENABLED", "DISABLED",
                                                spatial_reference)
            arcpy.CreateFeatureclass_management(optfolder, "AddasOutlets.shp", "POINT", "", "ENABLED", "DISABLED",
                                                spatial_reference)
            arcpy.CreateFeatureclass_management(optfolder, "AddasReservoir.shp", "POINT", "", "ENABLED", "DISABLED",
                                                spatial_reference)

            # *******************************************************************************************************
            # turn layers ON/OFF in current data frame
            # *******************************************************************************************************
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "Soils":
                    lyr.visible = False
                if lyr.name == "md quads":
                    lyr.visible = False
                if lyr.name == "mask":
                    arcpy.ApplySymbologyFromLayer_management(lyr, r"" + Directory + "/data/mdfiles/legends/extent.lyr")
                if lyr.name == "InfStreams":
                    lyr.name = "Inferred Streams"
                    arcpy.ApplySymbologyFromLayer_management(lyr, r"" + Directory + "/data/mdfiles/legends/infstreams.lyr")
                if lyr.name == "landuse":
                    lyr.visible = False
                if lyr.name == "md dem":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "ssurgo":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "dem":
                    arcpy.ApplySymbologyFromLayer_management(lyr, r"" + Directory + "/data/mdfiles/legends/dem.lyr")
                if lyr.name == "soilG_a":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "soilG_b":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "soilG_c":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "soilG_d":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "fill":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "nhd_rast":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "burned":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "nhdstreams":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "flowdir_dem":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "mask_dem":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "mask_aux":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "md quads":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "water_grid":
                    arcpy.mapping.RemoveLayer(df, lyr)

            arcpy.Delete_management(optfolder + "/nhdstreams.shp", "")
            arcpy.Delete_management(optfolder + "/burned", "")
            arcpy.Delete_management(optfolder + "/nhd_rast", "")
            arcpy.Delete_management(optfolder + "/fill", "")
            arcpy.Delete_management(optfolder + "/mask_dem", "")
            arcpy.Delete_management(optfolder + "/mask_aux.shp", "")
            arcpy.Delete_management(optfolder + "/water_grid.shp", "")
            del flowacc, fill, burned, flowdir, infstr

            lyr = arcpy.mapping.ListLayers(mxd, "dem", df)[0]
            df.extent = lyr.getSelectedExtent()

            tables = arcpy.mapping.ListTableViews(mxd, "", df)
            for tbl in tables:
                if tbl.name == "theVTab":
                    arcpy.mapping.RemoveTableView(df, tbl)

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()

            # ******************************************************************************************************
            # turn AoI OFF and watershed delineation ON
            # ******************************************************************************************************
            tool1.enabled = False
            tool2.enabled = True
            button16.enabled = True
            button24.enabled = True
            tool7.enabled = True
            button17.enabled = True
            button18.enabled = True
            button19.enabled = True
            button20.enabled = True
            button21.enabled = True


class DelineateSubwatersheds(object):
    """Implementation for GISHydroNXT_addin.button8 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        # extent is smaller when subwatersheds are delineated again after Reset
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder


        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            # *******************************************************************************************************
            # add subwatersheds to view -- addOutputsToMap
            # *******************************************************************************************************
            arcpy.env.addOutputsToMap = False
            arcpy.env.extent = "MAXOF"
            basingrid = optfolder + "/basingrid"
            flowacc = optfolder + "/flowacc"
            flwdir = arcpy.sa.Times(optfolder + "/flowdir", basingrid)
            strlnk = arcpy.sa.StreamLink(optfolder + "/ModStreams", flwdir)
            zonemax = arcpy.sa.ZonalStatistics(strlnk, "Value", flowacc, "MAXIMUM", "NODATA")
            faccmax = arcpy.sa.Con(flowacc == zonemax, zonemax, arcpy.sa.IsNull(zonemax))
            outlets = arcpy.sa.Con(faccmax > 0, faccmax)
            outlets.save(optfolder + "/outlets")

            # if outlets are added by user then add mod stream and user specified outlets
            outlets_user   = optfolder + "/outlets_user"
            if os.path.exists(outlets_user):
                aux_out1 = arcpy.RasterToPoint_conversion(outlets)
                aux_out2 = arcpy.RasterToPoint_conversion(outlets_user)
                aux_out3 = arcpy.Merge_management([aux_out1,aux_out2])
                outlets = arcpy.PointToRaster_conversion(aux_out3, "FID", optfolder + "/outlets2", "#", "#", basingrid)

            arcpy.env.extent = "MAXOF"
            subwshed = arcpy.sa.Watershed(flwdir, outlets, "VALUE")
            subwshed = arcpy.sa.Times(subwshed, basingrid)
            arcpy.env.addOutputsToMap = True
            arcpy.RasterToPolygon_conversion(subwshed,optfolder + "/tmpsubwshd","NO_SIMPLIFY","VALUE")

            arcpy.env.extent = "MAXOF"
            tmpsub = optfolder + "/tmpsubwshd.shp"
            subshed = optfolder + "/subshed.shp"
            arcpy.Dissolve_management(tmpsub,subshed,"GRIDCODE","#","MULTI_PART","DISSOLVE_LINES")

            subrivers = optfolder + "/subrivers.shp"
            # new stream link raster to handle added outlets
            if os.path.exists(outlets_user):
                StrmFDRgrid = arcpy.sa.Times(optfolder + "/ModStreams",flwdir)
                newLnkGrid = arcpy.sa.Watershed(StrmFDRgrid, outlets, "VALUE")
                newLnkGrid = arcpy.sa.Times(optfolder + "/ModStreams",newLnkGrid)
                arcpy.sa.StreamToFeature(newLnkGrid, flwdir, subrivers, "NO_SIMPLIFY")
                del aux_out1,aux_out2,aux_out3
            else:
                arcpy.sa.StreamToFeature(optfolder + "/ModStreams", flwdir, subrivers, "NO_SIMPLIFY")

            # *******************************************************************************************************
            # correct sub-basin indexing -- order of sub-basins has to be the same as
            # appearing in "subriver.shp"
            # *******************************************************************************************************
            target1 = optfolder + "/subrivers.shp"
            joint1 = optfolder + "/subshed.shp"
            spatial_join1 = optfolder + "/subriver_subshed.shp"

            # add table using spatial join (target = subriver, join = subshed)
            fieldmappings = arcpy.FieldMappings()
            fieldmappings.addTable(target1)
            fieldmappings.addTable(joint1)
            arcpy.SpatialJoin_analysis(target1, joint1, spatial_join1, "#", "#", fieldmappings, "HAVE_THEIR_CENTER_IN", "#", "#")

            # get subriver FID and centroid (or midpoint) XY values
            sub_cur = arcpy.SearchCursor(spatial_join1, "", "", "Shape;ARCID;GRIDCODE", "")
            x = []
            y = []
            for row in sub_cur:
                XMidPoint = row.shape.positionAlongLine(0.50, True).firstPoint.X
                x.append(XMidPoint)
                YMidPoint = row.shape.positionAlongLine(0.50, True).firstPoint.Y
                y.append(YMidPoint)

            del sub_cur

            xy = list(zip(x, y))
            # create a point shapefile, add fields "ARCID" & "GRIDCODE", and add midpoint XY values to field "SHAPE@"
            spatial_reference = arcpy.Describe(optfolder + "/subriver_subshed.shp").spatialReference
            arcpy.CreateFeatureclass_management(optfolder, "subriver_xy.shp", "POINT", "", "ENABLED", "DISABLED",
                                                spatial_reference)
            subriver_xy = optfolder + "/subriver_xy.shp"

            # added on 09-30-2014 to debug field type error -- changed field type to "Double" with precision "10"
            if not len(arcpy.ListFields(subriver_xy, "ARCID")) > 0:
                arcpy.AddField_management(subriver_xy, "ARCID", "DOUBLE", 10, "")
                arcpy.AddField_management(subriver_xy, "GRIDCODE", "DOUBLE", 10, "")
            sr_xy = arcpy.da.InsertCursor(subriver_xy, ("SHAPE@XY"))
            for i in xy:
                sr_xy.insertRow([i])

            # update "ARCID" and "GRIDCODE" row values using "subriver_subshed" shapefile
            sub_cur = arcpy.SearchCursor(spatial_join1, "", "", "ARCID;GRIDCODE", "")
            ic = arcpy.UpdateCursor(subriver_xy, "", "", "ARCID;GRIDCODE", "")

            row1 = sub_cur.next()
            row2 = ic.next()

            while row1:
                row2.setValue("ARCID", row1.getValue("ARCID"))
                ic.updateRow(row2)
                row2.setValue("GRIDCODE", row1.getValue("GRIDCODE"))
                ic.updateRow(row2)
                row1 = sub_cur.next()
                row2 = ic.next()
            del sub_cur, ic

            # perform spatial join (target = subshed, join = subriver_xy) and lists extraction
            target2 = optfolder + "/subshed.shp"
            joint2 = optfolder + "/subriver_xy.shp"
            spatial_join2 = optfolder + "/subshed_subriver.shp"

            fieldmappings = arcpy.FieldMappings()
            fieldmappings.addTable(target2)
            fieldmappings.addTable(joint2)
            arcpy.SpatialJoin_analysis(target2, joint2, spatial_join2, "#", "#", fieldmappings, "INTERSECT", "#", "#")

            # create a polygon which will contain subriver FID, ARCID, and GRIDCODE (obtained from spatially joined shapefile)
            spatial_reference = arcpy.Describe(optfolder + "/subshed_subriver.shp").spatialReference
            arcpy.CreateFeatureclass_management(optfolder, "subshed.shp", "POLYGON", "", "DISABLED", "DISABLED",
                                                spatial_reference)

            sc = arcpy.da.SearchCursor(spatial_join2, ("ARCID", "GRIDCODE", "SHAPE@"))

            storage = []
            for row in sc:  # changed sF to sc [could be a mistake to have "sF"]
                storage.append(row)

            sortedlist = sorted(storage, key=lambda x: x[0])  # list sort is based on "ARCID"

            poly = optfolder + "/subshed.shp"

            # added on 09-30-2014 to debug field type error -- changed field type to "Double" with precision "10"
            if not len(arcpy.ListFields(poly, "ARCID")) > 0:
                arcpy.AddField_management(poly, "ARCID", "DOUBLE", 10, "")
                arcpy.AddField_management(poly, "GRIDCODE", "DOUBLE", 10, "")

            arcpy.DeleteField_management(poly, "Id")
            fcout = arcpy.da.InsertCursor(poly, ("ARCID", "GRIDCODE", "SHAPE@"))

            for i in sortedlist:
                fcout.insertRow(i)
            del fcout, sc

            # convert newly created and indexed sub-watershed shapefile to raster
            arcpy.PolygonToRaster_conversion(subshed, "ARCID", optfolder + "/subsheds", "", "ARCID",
                                             optfolder + "/dem")  # priority field is "ARCID"

            # *******************************************************************************************************
            # turn layers ON/OFF in current data frame
            # *******************************************************************************************************
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "watershed":
                    lyr.visible = False
                if lyr.name == "AddasOutlets":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "AddasRerservoir":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "AddasStreams":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "ModStreams":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "Outlets_temp":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "tmpsubwshd":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subsheds":
                    df.extent = lyr.getSelectedExtent()
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subshed_subriver":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subriver_subshed":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subriver_xy":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subshed":
                    arcpy.ApplySymbologyFromLayer_management(lyr, r"" + Directory + "/data/mdfiles/legends/watershed.lyr")
                if lyr.name == "subrivers":
                    arcpy.ApplySymbologyFromLayer_management(lyr, r"" + Directory + "/data/mdfiles/legends/subrivers.lyr")
            arcpy.env.extent = "MAXOF"
            if os.path.exists(optfolder + "/outlets2"):
                arcpy.Delete_management(optfolder + "/outlets2", "")

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()

            # *******************************************************************************************************
            # turn subwatershed delineation OFF
            # *******************************************************************************************************

            button7.enabled = False
            button8.enabled = False
            button9.enabled = True
            tool5.enabled = False
            save(optfolder)

class ExecuteTR20(object):
    """Implementation for GISHydroNXT_addin.button15 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):

        if os.path.exists(optfolder + "/WinTR20/TR20out.txt"):
            os.remove(optfolder + "/WinTR20/TR20out.txt")
        if os.path.exists(optfolder + "/WinTR20/TR20.inp"):
            os.remove(optfolder + "/WinTR20/TR20.inp")
        if os.path.exists(optfolder + "/WinTR20/WinTR20_V32.exe"):
            os.remove(optfolder + "/WinTR20/WinTR20_V32.exe")
        if os.path.exists(optfolder + "/WinTR20/TR20err.txt"):
            os.remove(optfolder + "/WinTR20/TR20err.txt")
        if os.path.exists(optfolder + "/WinTR20/TR20.dbg"):
            os.remove(optfolder + "/WinTR20/TR20.dbg")

        shutil.copy2(optfolder + "/WinTR20/TR20in.txt", optfolder + "/WinTR20/TR20.inp")
        shutil.copy2(Directory + "/data/mdfiles/WinTR20/WinTR20_V32.exe", optfolder + "/WinTR20/")
        os.chdir(optfolder + "/WinTR20")
        os.system(optfolder + "/WinTR20/WinTR20_V32.exe")
        if os.path.exists(optfolder + "/WinTR20/TR20.err"):
            os.rename(optfolder + "/WinTR20/TR20.err", optfolder + "/WinTR20/TR20err.txt")

        # *******************************************************************************************************
        # If TR20 error file size is greater than "0" then display a warning message to check error file
        # *******************************************************************************************************
        if os.path.exists(optfolder + "/WinTR20/TR20err.txt"):
            filestats = os.path.getsize(optfolder + "/WinTR20/TR20err.txt")
            if filestats > 0:
                TR20_prompt = pythonaddins.MessageBox("Please check file " + optfolder + "/WinTR20/TR20err.txt" +
                                                      " for possible errors." + "\n" +
                                                      "Press OK to open error text file, or Cancel to continue",
                                                      "WinTR20 Input file error", 1)
                if TR20_prompt == "OK":
                    hydro.openbrowser(optfolder + "/WinTR20/TR20err.txt")

        # *******************************************************************************************************
        # Change extension of "TR20in.out" file from ".out" to ".txt" and open it in text editor
        # *******************************************************************************************************
        if os.path.exists(optfolder + "/WinTR20/TR20.out"):
            os.rename(optfolder + "/WinTR20/TR20.out", optfolder + "/WinTR20/TR20out.txt")
            hydro.openbrowser(optfolder + "/WinTR20/TR20out.txt")


class FlowPaths(object):
    """Implementation for GISHydroNXT_addin.tool4 (Tool)"""

    def __init__(self):
        self.enabled = False
        self.cursor = 3

    def onMouseDownMap(self, x, y, button, shift):
        # *******************************************************************************************************
        # add line to view -- addOutputsToMap
        # *******************************************************************************************************
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        arcpy.env.addOutputsToMap = True
        point = arcpy.Point(x, y)
        ptGeometry = arcpy.PointGeometry(point)
        theLine = arcpy.sa.CostPath(ptGeometry, optfolder + "/dem", optfolder + "/flowdir_dem", "BEST_SINGLE")
        theLine_masked = arcpy.sa.Times(theLine, optfolder + "/basingrid")
        arcpy.sa.StreamToFeature(theLine_masked, optfolder + "/flowdir", "line.shp", "NO_SIMPLIFY")

        searchcursor = arcpy.da.SearchCursor("line", ("OID@", "SHAPE@"))
        try:
            firstrow = searchcursor.next()
            insertcursor = arcpy.da.InsertCursor(optfolder + "/AddasStreams.shp", ("OID@", "SHAPE@"))
            insertcursor.insertRow(firstrow)
            del firstrow, insertcursor, searchcursor
        except:
            pass


        # *******************************************************************************************************
        # change stream flow path symbology
        # *******************************************************************************************************
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "line":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "AddasStreams":
                arcpy.ApplySymbologyFromLayer_management(lyr, r"" + Directory + "/data/mdfiles/legends/streams.lyr")

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()

        button1.enabled = False
        button5_1.enabled = False
        button6.enabled = True
        save(optfolder)

class GageListThomas(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Please select gage for Tasker analysis adjustment", size=(300, 150))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)

        wx.StaticText(panel, -1, "Watershed contains USGS Gage(s)", (15, 15))
        GageList = gagelist + ["Perform No Adjustment"]
        self.gage = wx.ComboBox(panel, -1, "'Select'", (30, 30), wx.DefaultSize, GageList, wx.CB_DROPDOWN)

        self.btnApply = wx.Button(panel, label="Apply", pos=(80, 65))
        self.Bind(wx.EVT_BUTTON, self.OnSet, id=self.btnApply.GetId())

        self.Show(True)
        self.Centre(True)
        style = self.GetWindowStyle()
        self.SetWindowStyle(style | wx.STAY_ON_TOP)

    def OnClose(self, event):
        self.Show(False)
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "outletbuffer":
                arcpy.Delete_management(lyr)
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "gagefield":
                arcpy.Delete_management(lyr)
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "mdlayer":
                arcpy.Delete_management(lyr)
                arcpy.mapping.RemoveLayer(df, lyr)

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()

    def OnSet(self, event):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        global gageid
        gageid = str(self.gage.GetValue())
        self.Show(False)

        # ******************************************************************************************************
        # read global variables (defined in basin stat) into variable names for use
        # in Thomas peak discharge analysis
        # ******************************************************************************************************
        FC = float(hydro.FC)  # LI is already declared as global variable
        FC = "{0:.2f}".format(FC)
        DA = areami2
        HA = float(hydro.pctAsoil)
        HC = float(hydro.pctCsoil)
        HD = float(hydro.pctDsoil)
        SLL = float(landslope)
        HCD = float(HC + HD)
        ImpA = float(IA)

        # ******************************************************************************************************
        # define initial values and lists. add and index those lists as part of dictionary
        # ******************************************************************************************************
        sumarea = 0
        provstring = ""
        outtaskerstring = ""
        sQ1p25 = 0
        sQ1p50 = 0
        sQ1p75 = 0
        sQ2 = 0
        sQ5 = 0
        sQ10 = 0
        sQ25 = 0
        sQ50 = 0
        sQ100 = 0
        sQ200 = 0
        sQ500 = 0
        Q1p25list = [0, 0, 0, 0, 0, 0, 0, 0]
        Q1p50list = [0, 0, 0, 0, 0, 0, 0, 0]
        Q2list = [0, 0, 0, 0, 0, 0, 0, 0]
        Q5list = [0, 0, 0, 0, 0, 0, 0, 0]
        Q10list = [0, 0, 0, 0, 0, 0, 0, 0]
        Q25list = [0, 0, 0, 0, 0, 0, 0, 0]
        Q50list = [0, 0, 0, 0, 0, 0, 0, 0]
        Q100list = [0, 0, 0, 0, 0, 0, 0, 0]
        Q200list = [0, 0, 0, 0, 0, 0, 0, 0]
        Q500list = [0, 0, 0, 0, 0, 0, 0, 0]

        qlist = {"1": [], "2": [], "3": [], "4": [], "5": [], "6": [], "7": [], "8": [], "9": [], "10": []}
        qlist["1"].extend(Q1p25list)
        qlist["2"].extend(Q1p50list)
        qlist["3"].extend(Q2list)
        qlist["4"].extend(Q5list)
        qlist["5"].extend(Q10list)
        qlist["6"].extend(Q25list)
        qlist["7"].extend(Q50list)
        qlist["8"].extend(Q100list)
        qlist["9"].extend(Q200list)
        qlist["10"].extend(Q500list)

        # ******************************************************************************************************
        # copy root tasker files into optfolder (not needed after tasker was recreated in python), the folder still needs to be created
        # ******************************************************************************************************
        src = r"" + Directory + "/data/tasker"
        dst = optfolder + "/tasker/"
        if not os.path.exists(dst):
            os.mkdir(dst, 0o755)
        for the_file in os.listdir(src):
            file_path = os.path.join(src, the_file)
            if os.path.isfile(file_path):
                shutil.copyfile(file_path, dst + the_file)

        # ******************************************************************************************************
        # read zonal stat table (theVTab), declare fields into variables, and
        # loop through theVTab fields to getValue
        # ******************************************************************************************************
        theVTab = optfolder + "/theVTab.dbf"
        theVTab = arcpy.SearchCursor("theVTab", "", "", "Count", "")

        # ******************************************************************************************************
        # loop to get total count of pixels of basingrid
        # ******************************************************************************************************
        for each in theVTab:
            count = each.getValue("Count")
            sumarea = sumarea + count
        sumArea = sumarea

        del each

        theVTab = optfolder + "/theVTab.dbf"
        theVTab = arcpy.SearchCursor("theVTab", "", "",
                                     "Province;Count;Q1.25;Q1.50;Q1.75;Q2;Q5;Q10;Q25;Q50;Q100;Q200;Q500", "")

        # loop to begin tasker handling and area weighted analysis
        discharge_values = []
        flood_intervals = []
        for row in theVTab:
            AreaField = float(row.getValue("Count"))
            areapercent = float((AreaField / sumArea) * 100)
            if row.getValue("Province") == "A":
                provstring = provstring + "       -Appalachian Plateaus and Allegheny Ridges %s percent of area" "\n" % (
                    "{0:.2f}".format(areapercent))
            elif row.getValue("Province") == "B":
                provstring = provstring + "       -Blue Ridge and Great Valley %s percent of area" "\n" % (
                    "{0:.2f}".format(areapercent))
            elif row.getValue("Province") == "P":
                provstring = provstring + "       -Piedmont %s percent of area" "\n" % ("{0:.2f}".format(areapercent))
            elif row.getValue("Province") == "W":
                provstring = provstring + "       -Western Coastal Plain %s percent of area" "\n" % (
                    "{0:.2f}".format(areapercent))
            elif row.getValue("Province") == "E":
                provstring = provstring + "       -Eastern Coastal Plain %s percent of area" "\n" % (
                    "{0:.2f}".format(areapercent))
            else:
                pythonaddins.MessageBox("No Province Selected", "Problem...", 0)

            intaskerstring = "thomasout.txt" "\n"
            intaskerstring = intaskerstring + "" "\n"
            if row.getValue("Province") == "A":
                intaskerstring = intaskerstring + "a"  "\n"
                intaskerstring = intaskerstring + "%s" "\n" % (DA)
                intaskerstring = intaskerstring + "%s" "\n" % (SLL)
            elif row.getValue("Province") == "B":
                intaskerstring = intaskerstring + "p"  "\n"
                intaskerstring = intaskerstring + "%s" "\n" % (DA)
                intaskerstring = intaskerstring + "%s" "\n" % (FC)
                intaskerstring = intaskerstring + "%s" "\n" % (LI)
                intaskerstring = intaskerstring + "%s" "\n" % (ImpA)
            elif row.getValue("Province") == "P":
                intaskerstring = intaskerstring + "p"  "\n"
                intaskerstring = intaskerstring + "%s" "\n" % (DA)
                intaskerstring = intaskerstring + "%s" "\n" % (FC)
                intaskerstring = intaskerstring + "%s" "\n" % (LI)
                intaskerstring = intaskerstring + "%s" "\n" % (ImpA)

            elif row.getValue("Province") == "W":
                intaskerstring = intaskerstring + "wc" "\n"
                intaskerstring = intaskerstring + "%s" "\n" % (DA)
                intaskerstring = intaskerstring + "%s" "\n" % (ImpA)
                intaskerstring = intaskerstring + "%s" "\n" % (HCD)
            elif row.getValue("Province") == "E":
                intaskerstring = intaskerstring + "ec" "\n"
                intaskerstring = intaskerstring + "%s" "\n" % (DA)
                intaskerstring = intaskerstring + "%s" "\n" % (SLL)
                intaskerstring = intaskerstring + "%s" "\n" % (HA)

            if gageid == "Perform No Adjustment":
                intaskerstring = intaskerstring + "N" "\n"
                intaskerstring = intaskerstring + "N" "\n"
            else:
                intaskerstring = intaskerstring + "Y" "\n"
                intaskerstring = intaskerstring + "%s" "\n" % (gageid)
                intaskerstring = intaskerstring + "N" "\n"

            thomas2020 = optfolder + "/tasker/thomas2020.txt"
            thomasin = open(thomas2020, "w")
            thomasin.write(intaskerstring)
            thomasin.close()


            tasker.RRE(optfolder,thomas2020)

            # REMOVED AFTER TASKER.PY
            #os.chdir(optfolder + "/tasker")
            #os.system(optfolder + "/tasker/thomas2016auto.exe")
            #time.sleep(4)

            # *****************************************************************************
            # open "thomasout.txt" file and read peak discharges
            # *****************************************************************************

            infilename = optfolder + "/tasker/thomasout.txt"
            infile = open(infilename, "r").readlines()

            for i in infile:
                outtaskerstring = outtaskerstring + i

            # *****************************************************************************
            # Discharge values for each province for "theVTab" are not possible to write
            # at this time because update cursor can"t be used within search cursor "for"
            # loop -- appended to list "discharge_values" and will be used at the end to
            # update discharge values in "theVTab"
            # *****************************************************************************
            theline = infile[11].split()
            Q1p25 = theline[1]
            discharge_values.append(Q1p25)
            theline = infile[12].split()
            Q1p50 = theline[1]
            discharge_values.append(Q1p50)
            Q1p75 = -999
            theline = infile[13].split()
            Q2 = theline[1]
            discharge_values.append(Q2)
            theline = infile[14].split()
            Q5 = theline[1]
            discharge_values.append(Q5)
            theline = infile[15].split()
            Q10 = theline[1]
            discharge_values.append(Q10)
            theline = infile[16].split()
            Q25 = theline[1]
            discharge_values.append(Q25)
            theline = infile[17].split()
            Q50 = theline[1]
            discharge_values.append(Q50)
            theline = infile[18].split()
            Q100 = theline[1]
            discharge_values.append(Q100)
            theline = infile[19].split()
            Q200 = theline[1]
            discharge_values.append(Q200)
            theline = infile[20].split()
            Q500 = theline[1]
            discharge_values.append(Q500)

            # *****************************************************************************
            # compute discharge and assign confidence intervals to qlist entries
            # *****************************************************************************
            Q1p25 = float(Q1p25)
            sQ1p25 = sQ1p25 + (Q1p25 * AreaField)
            Q1p50 = float(Q1p50)
            sQ1p50 = sQ1p50 + (Q1p50 * AreaField)
            Q1p75 = float(Q1p75)
            sQ1p75 = sQ1p75 + (Q1p75 * AreaField)
            Q2 = float(Q2)
            sQ2 = sQ2 + (Q2 * AreaField)
            Q5 = float(Q5)
            sQ5 = sQ5 + (Q5 * AreaField)
            Q10 = float(Q10)
            sQ10 = sQ10 + (Q10 * AreaField)
            Q25 = float(Q25)
            sQ25 = sQ25 + (Q25 * AreaField)
            Q50 = float(Q50)
            sQ50 = sQ50 + (Q50 * AreaField)
            Q100 = float(Q100)
            sQ100 = sQ100 + (Q100 * AreaField)
            Q200 = float(Q200)
            sQ200 = sQ200 + (Q200 * AreaField)
            Q500 = float(Q500)
            sQ500 = sQ500 + (Q500 * AreaField)

            for i, lines in enumerate(infile):
                if i > 24 and i < 35:
                    line = lines.split()
                    a1 = line[1]
                    a2 = line[2]
                    a3 = line[3]
                    a4 = line[4]
                    a5 = line[5]
                    a6 = line[6]
                    a7 = line[7]
                    a8 = line[8]
                    c1 = [(float(a1) * areapercent) / 100]
                    c2 = [(float(a2) * areapercent) / 100]
                    c3 = [(float(a3) * areapercent) / 100]
                    c4 = [(float(a4) * areapercent) / 100]
                    c5 = [(float(a5) * areapercent) / 100]
                    c6 = [(float(a6) * areapercent) / 100]
                    c7 = [(float(a7) * areapercent) / 100]
                    c8 = [(float(a8) * areapercent) / 100]
                    qlist = [c1, c2, c3, c4, c5, c6, c7, c8]
                    flood_intervals.append(qlist)

        lists = {i: [el[0] for el in v] for i, v in enumerate(flood_intervals, start=1)}

        del row

        # ******************************************************************************************************
        # index and count number of sub-lists -- prepare for TaskerString function
        # ******************************************************************************************************
        number = hydro.FloodIntervalLists(lists)
        if number > 10:
            q_list1 = [sum(i) for i in zip(lists[1], lists[11])]
            q_list2 = [sum(i) for i in zip(lists[2], lists[12])]
            q_list3 = [sum(i) for i in zip(lists[3], lists[13])]
            q_list4 = [sum(i) for i in zip(lists[4], lists[14])]
            q_list5 = [sum(i) for i in zip(lists[5], lists[15])]
            q_list6 = [sum(i) for i in zip(lists[6], lists[16])]
            q_list7 = [sum(i) for i in zip(lists[7], lists[17])]
            q_list8 = [sum(i) for i in zip(lists[8], lists[18])]
            q_list9 = [sum(i) for i in zip(lists[9], lists[19])]
            q_list10 = [sum(i) for i in zip(lists[10], lists[20])]
        else:
            q_list1 = lists[1]
            q_list2 = lists[2]
            q_list3 = lists[3]
            q_list4 = lists[4]
            q_list5 = lists[5]
            q_list6 = lists[6]
            q_list7 = lists[7]
            q_list8 = lists[8]
            q_list9 = lists[9]
            q_list10 = lists[10]

        # ******************************************************************************************************
        # discharge computation based on province
        # ******************************************************************************************************
        Q1p25 = int(sQ1p25 / sumArea)
        Q1p50 = int(sQ1p50 / sumArea)
        Q1p75 = int(sQ1p75 / sumArea)
        Q2 = int(sQ2 / sumArea)
        Q5 = int(sQ5 / sumArea)
        Q10 = int(sQ10 / sumArea)
        Q25 = int(sQ25 / sumArea)
        Q50 = int(sQ50 / sumArea)
        Q100 = int(sQ100 / sumArea)
        Q200 = int(sQ200 / sumArea)
        Q500 = int(sQ500 / sumArea)

        # ******************************************************************************************************
        # Basin statistics analysis date
        # ******************************************************************************************************
        now = datetime.now()
        month = now.strftime("%B")
        day = now.strftime("%d")
        year = now.strftime("%Y")

        # ******************************************************************************************************
        # Text file string variables
        # ******************************************************************************************************
        datastring = ""
        datastring = datastring + "GISHydro Release Version Date: %s" "\n" % (Modifieddt)
        datastring = datastring + "Project Name:                  %s" % (proj)
        datastring = datastring + "" "\n"
        datastring = datastring + "Analysis Date:                 %s %s, %s " "\n" % (month, day, year)
        datastring = datastring + "Thomas Version:                %s " "\n \n" % (thomasversion)
        datastring = datastring + "Geographic Province(s):" "\n"
        datastring = datastring + provstring
        datastring = datastring + "" "\n"
        datastring = datastring + "Q(1.25):   %s cfs" "\n" % (Q1p25)
        datastring = datastring + "Q(1.50):   %s cfs" "\n" % (Q1p50)
        datastring = datastring + "Q(2):      %s cfs" "\n" % (Q2)
        datastring = datastring + "Q(5):      %s cfs" "\n" % (Q5)
        datastring = datastring + "Q(10):     %s cfs" "\n" % (Q10)
        datastring = datastring + "Q(25):     %s cfs" "\n" % (Q25)
        datastring = datastring + "Q(50):     %s cfs" "\n" % (Q50)
        datastring = datastring + "Q(100):    %s cfs" "\n" % (Q100)
        datastring = datastring + "Q(200):    %s cfs" "\n" % (Q200)
        datastring = datastring + "Q(500):    %s cfs" "\n" % (Q500)
        datastring = datastring + "" "\n"
        datastring = datastring + "Area Weighted Prediction Intervals (from Tasker)" "\n"
        datastring = datastring + " Return     50 PERCENT        67 PERCENT        90 PERCENT        95 PERCENT" "\n"
        datastring = datastring + " Period  lower    upper    lower    upper    lower    upper    lower    upper" "\n"
        datastring = datastring + hydro.TaskerString(1.25, q_list1)
        datastring = datastring + "" "\n"
        datastring = datastring + hydro.TaskerString(1.50, q_list2)
        datastring = datastring + "" "\n"
        datastring = datastring + hydro.TaskerString(2, q_list3)
        datastring = datastring + "" "\n"
        datastring = datastring + hydro.TaskerString(5, q_list4)
        datastring = datastring + "" "\n"
        datastring = datastring + hydro.TaskerString(10, q_list5)
        datastring = datastring + "" "\n"
        datastring = datastring + hydro.TaskerString(25, q_list6)
        datastring = datastring + "" "\n"
        datastring = datastring + hydro.TaskerString(50, q_list7)
        datastring = datastring + "" "\n"
        datastring = datastring + hydro.TaskerString(100, q_list8)
        datastring = datastring + "" "\n"
        datastring = datastring + hydro.TaskerString(200, q_list9)
        datastring = datastring + "" "\n"
        datastring = datastring + hydro.TaskerString(500, q_list10)
        datastring = datastring + "" "\n"
        datastring = datastring + "" "\n"
        datastring = datastring + "" "\n"
        datastring = datastring + "Individual Province Tasker Analyses Follow: " "\n"
        datastring = datastring + ""
        datastring = datastring + outtaskerstring

        # ******************************************************************************************************
        # write strings to basin stat text file.
        # ******************************************************************************************************
        defFN = optfolder + "/frdischarges.txt"
        fr = open(defFN, "w")
        fr.write(datastring)
        fr.close()

        # ******************************************************************************************************
        # open "frdischarges" file in text editor
        # ******************************************************************************************************

        hydro.openbrowser(defFN)

        # ******************************************************************************************************
        # turn layers ON/OFF in current data frame
        # ******************************************************************************************************
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "mdlayer":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "gagefield":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "outletpoint":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "outletpoly":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "limepoly":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "outletpoint":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "outletbuffer":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "gagefield":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "outletpoly":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "mask_ints":  # added on 10-22-2017: this layer is output of intersect tool
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "gauge_outlet":  # added on 10-22-2017: this layer is output of join and intersect tool
                arcpy.mapping.RemoveLayer(df, lyr)

        # ******************************************************************************************************
        # turn "peak discharge" OFF, and "S", "Reset Sub-watershed & Add Streams ON
        # ******************************************************************************************************
        button4.enabled = False
        button5_1.enabled = True
        button6.enabled = True
        tool4.enabled = True

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()
        save(optfolder)


class HelpAbout(object):
    """Implementation for GISHydroNXT_addin.button21 (Button)"""

    def __init__(self):
        self.enabled = True

    def onClick(self):
        pythonaddins.MessageBox("This is GISHydroNXT for ArcGIS " + str(ArcGIS_Version) + ".\n"\
                                + "Version: " + Modifieddt, "About")


class HelpHTML(object):
    """Implementation for GISHydroNXT_addin.button17 (Button)"""

    def __init__(self):
        self.enabled = True

    def onClick(self):
        hydro.openbrowser(Directory + r"data/help/NXTHTML/index.html")


class HelpManualToolGuide(object):
    """Implementation for GISHydroNXT_addin.button20 (Button)"""

    def __init__(self):
        self.enabled = True

    def onClick(self):
        hydro.openbrowser(Directory + r"data/help/Users_Manual_Tool_Guide.pdf")


class HelpTechRef(object):
    """Implementation for GISHydroNXT_addin.button18 (Button)"""

    def __init__(self):
        self.enabled = True

    def onClick(self):
        hydro.openbrowser(Directory + r"data/help/GISHydroNXT_TechnicalReference.pdf")


class HelpTrainingManual(object):
    """Implementation for GISHydroNXT_addin.button19 (Button)"""

    def __init__(self):
        self.enabled = True

    def onClick(self):
        hydro.openbrowser(Directory + r"data/help/TrainingManual_June2019.pdf")


class LandUseEditor(object):
    """Implementation for GISHydroNXT_addin.tool7 (tool)"""

    def __init__(self):
        self.enabled = False
        self.cursor = 3
        self.shape = "Line"

    def onLine(self, line_geometry):
        arcpy.env.addOutputsToMap = False
        arcpy.env.overwriteOutput = True
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        arcpy.env.extent = "MAXOF"

        aux_points = optfolder + "/aux_folder/landuse_points.shp"
        aux_polyline = optfolder + "/aux_folder/landuse_polyline.shp"
        aux_polygon = optfolder + "/aux_folder/landuse_polygon.shp"
        arcpy.CopyFeatures_management(line_geometry, aux_polyline)
        arcpy.FeatureVerticesToPoints_management(aux_polyline, aux_points)
        arcpy.PointsToLine_management(aux_points, aux_polyline, "#", "#", "CLOSE")
        arcpy.env.addOutputsToMap = True
        arcpy.FeatureToPolygon_management(aux_polyline, aux_polygon)

        LandUsePoly()


class LandUsePoly(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Land Use Editor", size=(380, 290))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)

        wx.StaticText(panel, -1, "Land Use Category Name:", (20, 20))
        self.categoryname = wx.TextCtrl(panel, -1, value="My Land Use", pos=(190, 20), size=(155, 20))

        wx.StaticText(panel, -1, "Enter Imperviousness (%):", (20, 50))
        self.lu_imp = wx.TextCtrl(panel, -1, value="39", pos=(190, 50), size=(40, 20), style=wx.TE_RIGHT)
        wx.StaticText(panel, -1, "A soil CN:", (20, 90))
        self.lu_a = wx.TextCtrl(panel, -1, value="39", pos=(85, 90), size=(40, 20), style=wx.TE_RIGHT)
        wx.StaticText(panel, -1, "B soil CN:", (20, 115))
        self.lu_b = wx.TextCtrl(panel, -1, value="61", pos=(85, 115), size=(40, 20), style=wx.TE_RIGHT)
        wx.StaticText(panel, -1, "C soil CN:", (20, 140))
        self.lu_c = wx.TextCtrl(panel, -1, value="75", pos=(85, 140), size=(40, 20), style=wx.TE_RIGHT)
        wx.StaticText(panel, -1, "D soil CN:", (20, 165))
        self.lu_d = wx.TextCtrl(panel, -1, value="80", pos=(85, 165), size=(40, 20), style=wx.TE_RIGHT)

        wx.StaticBox(panel, -1, "Major Land Use Category", (170, 80), size=(175, 105))
        self.landusecat1 = wx.RadioButton(panel, -1, "None", (225, 100), style=wx.RB_GROUP)
        self.landusecat2 = wx.RadioButton(panel, -1, "Storage", (225, 120))
        self.landusecat3 = wx.RadioButton(panel, -1, "Forest", (225, 140))
        self.landusecat4 = wx.RadioButton(panel, -1, "Urban", (225, 160))

        self.lu_imp.Bind(wx.EVT_TEXT, self.ChangeSoil)
        self.landusecat1.Bind(wx.EVT_RADIOBUTTON, self.ChangeSoil)
        self.landusecat2.Bind(wx.EVT_RADIOBUTTON, self.ChangeSoil)
        self.landusecat3.Bind(wx.EVT_RADIOBUTTON, self.ChangeSoil)
        self.landusecat4.Bind(wx.EVT_RADIOBUTTON, self.ChangeSoil)

        self.btnEdit = wx.Button(panel, label="Revise Curve Number", pos=(35, 200))
        self.Bind(wx.EVT_BUTTON, self.OnRevise, id=self.btnEdit.GetId())
        self.btnCancel = wx.Button(panel, label="Cancel", pos=(230, 200))
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=self.btnCancel.GetId())

        self.Show(True)
        self.Centre(True)
        style = self.GetWindowStyle()
        self.SetWindowStyle(style | wx.STAY_ON_TOP)

    def ChangeSoil(self, event):
        if self.lu_imp.GetValue().isdigit() == True:
            if int(self.lu_imp.GetValue()) > 100:
                self.lu_imp.SetValue("100")
            elif int(self.lu_imp.GetValue()) < 0:
                self.lu_imp.SetValue("0")
            soila_lu = str(
                int(float(self.lu_imp.GetValue()) * 98 / 100 + (1 - float(self.lu_imp.GetValue()) / 100) * 39))
            soilb_lu = str(
                int(float(self.lu_imp.GetValue()) * 98 / 100 + (1 - float(self.lu_imp.GetValue()) / 100) * 61))
            soilc_lu = str(
                int(float(self.lu_imp.GetValue()) * 98 / 100 + (1 - float(self.lu_imp.GetValue()) / 100) * 75))
            soild_lu = str(
                int(float(self.lu_imp.GetValue()) * 98 / 100 + (1 - float(self.lu_imp.GetValue()) / 100) * 80))
            self.lu_a.SetValue(soila_lu)
            self.lu_b.SetValue(soilb_lu)
            self.lu_c.SetValue(soilc_lu)
            self.lu_d.SetValue(soild_lu)
        else:
            self.lu_imp.SetValue("")
        if self.landusecat2.GetValue() == True:
            self.lu_imp.SetValue("100")
            self.lu_imp.Enable(False)
        elif self.landusecat3.GetValue() == True:
            self.lu_imp.SetValue("0")
            self.lu_imp.Enable(False)
        else:
            self.lu_imp.Enable(True)

    def OnRevise(self, event):
        arcpy.env.addOutputsToMap = False
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        self.Show(False)

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            cna = int(self.lu_a.GetValue())
            cnb = int(self.lu_b.GetValue())
            cnc = int(self.lu_c.GetValue())
            cnd = int(self.lu_d.GetValue())

            catname = self.categoryname.GetValue()
            imp = int(self.lu_imp.GetValue())
            if self.landusecat1.GetValue():
                cat = "NIL"
            elif self.landusecat2.GetValue():
                cat = "ST"
            elif self.landusecat3.GetValue():
                cat = "FC"
            else:
                cat = "UrbPct"

            mluc = 0
            with arcpy.da.UpdateCursor(optfolder + "/landuse", ["VALUE", "CLASS_NAME"]) as cursor:
                for row in cursor:
                    if row[0] >= mluc:
                        mluc = row[0] + 1000
                    if row[1] == catname:
                        catname = catname + " copy"

            lu_update = arcpy.UpdateCursor(optfolder + "/aux_folder/landuse_polygon.shp")
            for lu in lu_update:
                lu.ID = mluc
                lu_update.updateRow(lu)

            arcpy.env.snapRaster = optfolder + "/landuse"
            arcpy.PolygonToRaster_conversion(optfolder + "/aux_folder/landuse_polygon.shp", "ID",
                                             optfolder + "/aux_folder/auxlu_raster", "#", "#", optfolder + "/dem")
            arcpy.CopyRaster_management(optfolder + "/landuse", optfolder + "/aux_folder/landuse_aux")

            if not os.path.exists(optfolder + "/aux_folder/landuse"):
                arcpy.CopyRaster_management(optfolder + "/landuse", optfolder + "/aux_folder/landuse")

            global landedit
            landedit.append([mluc, catname, cat, imp, cna, cnb, cnc, cnd])

            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "landuse":
                    arcpy.mapping.RemoveLayer(df, lyr)
            if os.path.exists(optfolder + "/landuse"):
                arcpy.Delete_management(optfolder + "/landuse", "")

            arcpy.MosaicToNewRaster_management(optfolder + "/aux_folder/landuse_aux;" + optfolder + "/aux_folder/auxlu_raster", optfolder, "landuse", "#", "32_BIT_SIGNED", "#", "1", "LAST", "FIRST")

            arcpy.JoinField_management(optfolder + "/landuse", "VALUE", optfolder + "/aux_folder/landuse_aux", "Value",
                                       "CLASS_NAME;RED;GREEN;BLUE;OPACITY")
            with arcpy.da.UpdateCursor(optfolder + "/landuse", ["VALUE", "CLASS_NAME", "RED", "GREEN", "BLUE", "OPACITY"]) as cursor:
                for row in cursor:
                    if row[0] == mluc:
                        if self.landusecat1.GetValue():
                            row[1] = catname
                            row[2] = 1  # RED
                            row[3] = 1  # GREEN
                            row[4] = 1  # BLUE
                            row[5] = 1  # OPACITY
                        elif self.landusecat2.GetValue():
                            row[1] = catname
                            row[2] = 0  # RED
                            row[3] = 0  # GREEN
                            row[4] = 1  # BLUE
                            row[5] = 1  # OPACITY
                        elif self.landusecat3.GetValue():
                            row[1] = catname
                            row[2] = 0  # RED
                            row[3] = 1  # GREEN
                            row[4] = 0  # BLUE
                            row[5] = 1  # OPACITY
                        elif self.landusecat4.GetValue():
                            row[1] = catname
                            row[2] = 1  # RED
                            row[3] = 0  # GREEN
                            row[4] = 0  # BLUE
                            row[5] = 1  # OPACITY
                    cursor.updateRow(row)

            auxr = r"" + optfolder + "/aux_folder/auxcn_raster"
            arcpy.CopyRaster_management(optfolder + "/curvenumber", optfolder + "/aux_folder/curvenum_aux")
            arcpy.Clip_management(optfolder + "/soils", "#", auxr, optfolder + "/aux_folder/landuse_polygon.shp", "#",
                                  "ClippingGeometry")
            outCon = arcpy.sa.Con(arcpy.Raster(auxr) == 1, cna, arcpy.sa.Con(arcpy.Raster(auxr) == 2, cnb,
                                               arcpy.sa.Con(arcpy.Raster(auxr) == 3, cnc,
                                                            arcpy.sa.Con(arcpy.Raster(auxr) == 4, cnd, 0))))
            outCon.save(optfolder + "/aux_folder/outCon")
            arcpy.MosaicToNewRaster_management(
                optfolder + "/aux_folder/curvenum_aux" + ";" + optfolder + "/aux_folder/outCon", optfolder, "/curvenumber",
                "#", "#", "#", "1", "LAST", "FIRST")

            arcpy.env.addOutputsToMap = True
            arcpy.MakeRasterLayer_management(optfolder + "/landuse", "Land Use")

            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "landuse_points":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "landuse_polyline":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "landuse_polygon":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "Land Use":
                    lyr.visible = False
            if os.path.exists(auxr):
                arcpy.Delete_management(auxr, "")
            if os.path.exists(optfolder + "/aux_folder/landuse_points.shp"):
                arcpy.Delete_management(optfolder + "/aux_folder/landuse_points.shp", "")
            if os.path.exists(optfolder + "/aux_folder/landuse_polyline.shp"):
                arcpy.Delete_management(optfolder + "/aux_folder/landuse_polyline.shp")
            if os.path.exists(optfolder + "/aux_folder/landuse_polygon.shp"):
                arcpy.Delete_management(optfolder + "/aux_folder/landuse_polygon.shp")
            if os.path.exists(optfolder + "/slope_calc"):
                arcpy.Delete_management(optfolder + "/slope_calc", "")

            button1.enabled = True
            if button3.enabled:
                button2.enabled = True
            button3.enabled = False
            button4.enabled = False

        pythonaddins.MessageBox("Land Use Correctly Modified", "Land Use Editor")
        save(optfolder)

    def OnClose(self, event):
        self.Show(False)
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "landuse_points":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "landuse_polyline":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "landuse_polygon":
                arcpy.mapping.RemoveLayer(df, lyr)
        if os.path.exists(optfolder + "/aux_folder/landuse_points.shp"):
            arcpy.Delete_management(optfolder + "/aux_folder/landuse_points.shp", "")
        if os.path.exists(optfolder + "/aux_folder/landuse_polyline.shp"):
            arcpy.Delete_management(optfolder + "/aux_folder/landuse_polyline.shp")
        if os.path.exists(optfolder + "/aux_folder/landuse_polygonlanduse_polygon.shp"):
            arcpy.Delete_management(optfolder + "/aux_folder/landuse_polygon.shp")


class LoadButton(object):
    """Implementation for GISHydroNXT_addin.button23(Button)"""

    def __init__(self):
        self.enabled = True
        self.checked = False

    def onClick(self):
        global scratchfolder
        global optfolder

        mxdloc = arcpy.mapping.MapDocument("CURRENT")
        mxdloc = mxdloc.filePath
        mxdloc = mxdloc.replace("\\", "/")
        data = r"" + mxdloc.replace("/umdgism/mdinterface/GISHydroNXT.mxd", "/") + "/temp/"
        dir_aux = r"" + mxdloc.replace("mdinterface/GISHydroNXT.mxd", "data")

        if not os.path.exists(dir_aux):
            pythonaddins.MessageBox("Error: GISHydroNXT data folder is missing", "Load Error", 0)
            return

        # *******************************************************************************************************
        # checkout ArcGIS 10.1 extensions (Spatial and 3D Analyst)
        # *******************************************************************************************************
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            pythonaddins.MessageBox("Error: Couldn't get Spatial Analyst extension", "Load Error", 0)
            return
        if arcpy.CheckExtension("3D") == "Available":
            arcpy.CheckOutExtension("3D")
        else:
            pythonaddins.MessageBox("Error: Couldn't get 3D Analyst extension", "Load Error", 0)
            return

        choices = []
        for root, dirs, files in os.walk(data):
            for file in files:
                if file.endswith(".mxd"):
                     choices.append(file.replace(".mxd",""))
        dialog = wx.SingleChoiceDialog(None, 'Choose a project', 'Saved Projects', choices[::-1])
        if dialog.ShowModal() == wx.ID_OK:
            pass
            dialog.Destroy()
        else:
            dialog.Destroy()
            return

        selected_proj = dialog.GetStringSelection()
        optfolder = data + selected_proj + "/"
        scratchfolder = data.replace("temp", "intermediate")

        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        load(optfolder)

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()


class MergedSegTable(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Segment Merge Attributes - Sub-area " + str(arcid_global + 1),
                          size=(1240, 570))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)
        self.index = 0

        self.SetPosition((530, 420))

        # define attribute table frame and loop over attributes
        self.list_ctrl = wx.ListCtrl(panel, size=(1220, 520),
                                     style=wx.LC_REPORT | wx.LC_HRULES | wx.LC_VRULES | wx.BORDER_SUNKEN)
        self.list_ctrl.InsertColumn(1, "FID", width=70)
        #        self.list_ctrl.InsertColumn(2, "Shape", width=80)
        self.list_ctrl.InsertColumn(2, "UpPixel", width=70)
        self.list_ctrl.InsertColumn(3, "SegName", width=70)
        self.list_ctrl.InsertColumn(4, "Type", width=70)
        self.list_ctrl.InsertColumn(5, "DownPixel", width=70)
        self.list_ctrl.InsertColumn(6, "Avg. Area", width=70)
        # self.list_ctrl.InsertColumn(8, "DS Area", width=80)
        self.list_ctrl.InsertColumn(7, "UpElev", width=70)
        self.list_ctrl.InsertColumn(8, "DownElev", width=70)
        self.list_ctrl.InsertColumn(9, "Slope", width=70)
        self.list_ctrl.InsertColumn(10, "Width", width=70)
        self.list_ctrl.InsertColumn(11, "Depth", width=70)
        self.list_ctrl.InsertColumn(12, "Xarea", width=70)
        self.list_ctrl.InsertColumn(13, "I_Length", width=70)
        self.list_ctrl.InsertColumn(14, "Tot_Length", width=70)
        self.list_ctrl.InsertColumn(15, "Vel.", width=70)
        self.list_ctrl.InsertColumn(16, "I_Time", width=70)
        self.list_ctrl.InsertColumn(17, "Tot_Time", width=70)

        global attlist
        attlist_rev = reversed(attlist)
        for item in attlist_rev:
            self.list_ctrl.InsertStringItem(self.index, str(int(item[0])))
            #            self.list_ctrl.SetStringItem(self.index, 1, str(item[1]))
            self.list_ctrl.SetStringItem(self.index, 1, str(int(item[2])))
            self.list_ctrl.SetStringItem(self.index, 2, str(item[3]))
            self.list_ctrl.SetStringItem(self.index, 3, str(item[4]))
            self.list_ctrl.SetStringItem(self.index, 4, str(int(item[5])))
            self.list_ctrl.SetStringItem(self.index, 5, str(item[6]))
            #            self.list_ctrl.SetStringItem(self.index, 7, str(item[7]))
            self.list_ctrl.SetStringItem(self.index, 6, str(item[8]))
            self.list_ctrl.SetStringItem(self.index, 7, str(item[9]))
            self.list_ctrl.SetStringItem(self.index, 8, str(item[10]))
            self.list_ctrl.SetStringItem(self.index, 9, str(item[11]))
            self.list_ctrl.SetStringItem(self.index, 10, str(item[12]))
            self.list_ctrl.SetStringItem(self.index, 11, str(item[13]))
            self.list_ctrl.SetStringItem(self.index, 12, str(item[14]))
            self.list_ctrl.SetStringItem(self.index, 13, str(item[15]))
            self.list_ctrl.SetStringItem(self.index, 14, str(item[16]))
            self.list_ctrl.SetStringItem(self.index, 15, str(item[17]))
            self.list_ctrl.SetStringItem(self.index, 16, str(item[18]))
        self.Show(True)

        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder

        global notable
        if not os.path.exists(optfolder + "/vel_meth/attribute_tables/merged_tables"):
            os.makedirs(optfolder + "/vel_meth/attribute_tables/merged_tables")
        with open(optfolder + "/vel_meth/attribute_tables/merged_tables/" + "Segment_Merge_Attributes_Sub-area_" + str(
                arcid_global + 1) + "_(" + str(int(notable[arcid_global])) + ").csv", "wb") as f:
            w = csv.writer(f)
            w.writerow(["FID", "Shape", "UpPixel", "SegName", "Type", "DownPixel", "Avg. Area", "DS Area", "UpElev",
                        "DownElev", "Slope",
                        "Width", "Depth", "Xarea", "I_Length", "Tot_Length", "Vel.", "I_Time", "Tot_Time"])
            w.writerows(attlist)
        notable[arcid_global] = notable[arcid_global] + 1

    def OnClose(self, event):
        self.Show(False)


class PrcpFrequency(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Precipitation Frequency & Duration Selector", size=(470, 550))
        self.Bind(wx.EVT_CLOSE, self.onClose)
        panel = wx.Panel(self, -1)

        wx.StaticText(panel, -1, "Check desired storms:", (15, 15))
        wx.StaticText(panel, -1, "1-year", (30, 75))
        wx.StaticText(panel, -1, "6-hour", (115, 50))
        wx.StaticText(panel, -1, "12-hour", (188, 50))
        wx.StaticText(panel, -1, "24-hour", (268, 50))
        wx.StaticText(panel, -1, "48-hour", (346, 50))

        wx.StaticText(panel, -1, "2-year", (30, 105))
        wx.StaticText(panel, -1, "5-year", (30, 135))
        wx.StaticText(panel, -1, "10-year", (30, 165))
        wx.StaticText(panel, -1, "25-year", (30, 195))
        wx.StaticText(panel, -1, "50-year", (30, 225))
        wx.StaticText(panel, -1, "100-year", (30, 255))
        wx.StaticText(panel, -1, "200-year", (30, 285))
        wx.StaticText(panel, -1, "500-year", (30, 315))

        cbXYList = [(120, 75), (200, 75), (278, 75), (356, 75),
                    (120, 105), (200, 105), (278, 105), (356, 105),
                    (120, 135), (200, 135), (278, 135), (356, 135),
                    (120, 165), (200, 165), (278, 165), (356, 165),
                    (120, 195), (200, 195), (278, 195), (356, 195),
                    (120, 225), (200, 225), (278, 225), (356, 225),
                    (120, 255), (200, 255), (278, 255), (356, 255),
                    (120, 285), (200, 285), (278, 285), (356, 285),
                    (120, 315), (200, 315), (278, 315), (356, 315)]

        self.cb_list = []
        for pos in cbXYList:
            cb = wx.CheckBox(panel, -1, "", pos)
            cb.SetValue(False)
            self.cb_list.append(cb)

        wx.StaticText(panel, -1, "Output Storm Depths to File", (70, 355))
        self.cb = wx.CheckBox(panel, -1, "", (50, 355))
        self.cb.SetValue(True)

        self.btnSelect = wx.Button(panel, label="Select All", pos=(45, 395))
        self.Bind(wx.EVT_BUTTON, self.OnSelectAll, id=self.btnSelect.GetId())

        self.btnUnSelectAll = wx.Button(panel, label="Unselect All*", pos=(173, 395))
        self.Bind(wx.EVT_BUTTON, self.OnUnSelectAll, id=self.btnUnSelectAll.GetId())

        self.btnApply = wx.Button(panel, label="Apply/Close", pos=(305, 395))
        self.Bind(wx.EVT_BUTTON, self.OnApply, id=self.btnApply.GetId())

        wx.StaticText(panel, -1, "* Note: 'Unselect All' button will not unselect storms that" "\n"
                                 "   have already determined", (20, 450))

        self.Show(True)
        self.Centre(True)
        style = self.GetWindowStyle()
        self.SetWindowStyle(style | wx.STAY_ON_TOP)

    def OnSelectAll(self, event):
        for cb in self.cb_list:
            cb.SetValue(True)

    def OnUnSelectAll(self, event):
        for cb in self.cb_list:
            cb.SetValue(False)

    def OnApply(self, event):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        self.Show(False)

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            # ******************************************************************************************************
            # Use selected duration and year to compute average precipitation
            # ******************************************************************************************************
            global year
            year = []
            global critdur
            critdur = []
            global critavg
            critavg = []
            global cb_selected
            cb_selected = []
            pr_string = ""
            for i, cb in enumerate(self.cb_list):
                if cb.GetValue():
                    # "1" is added to selected boxes to start indexing from 1
                    cb_selected.append(i)  # "selected" index list from 1
                    dem = optfolder + "/dem"
                    basingrid = arcpy.Raster(dem)
                    cellsize = basingrid.meanCellWidth
                    arcpy.env.cellSize = str(cellsize)
                    analysis_extent = basingrid.extent
                    arcpy.env.extent = analysis_extent

                    # following year and duration list is to create avg prec list
                    yearlist = ["1", "2", "5", "10", "25", "50", "100", "200", "500"]
                    durlist = ["05n", "10n", "15n", "30n", "01", "02", "03", "06", "12", "24", "48"]
                    theyear = yearlist[i // 4]  # "//" will floor the value to get respective indexed year from above
                    year.append(theyear)
                    thecritdur = durlist[(i) % 4 + 7]
                    critdur.append(thecritdur)
                    pr_path = Directory + "/data/prec/"
                    PrecipFreq = hydro.GetPrecipFrequency(pr_path, durlist, theyear, thecritdur, optfolder)
                    thecritavg = round(PrecipFreq[1],2)
                    critavg.append(thecritavg)
                    pr_string = pr_string + "     " + "%s-year, %s-hour: %s inches" "\n" % (
                        theyear, thecritdur, thecritavg)

            # ******************************************************************************************************
            # Text file string variables
            # ******************************************************************************************************
            datastring = ""
            datastring = datastring + "GISHydro Release Version Date:    %s" "\n" % (Modifieddt)
            datastring = datastring + "Project Name:                     %s" % (proj)
            datastring = datastring + "" "\n"
            datastring = datastring + "Data Selected:" "\n"
            datastring = datastring + "     Outlet Easting:              %s m (MD Stateplane, NAD 1983)" "\n" % (
                xoutletstring)
            datastring = datastring + "     Outlet Northing:             %s m (MD Stateplane, NAD 1983)" "\n" % (
                youtletstring)
            datastring = datastring + "Precipitation Frequency-Duration Depths:" "\n"
            datastring = datastring + pr_string

            # ******************************************************************************************************
            # write strings to precipitation depth text file.
            # ******************************************************************************************************
            defFN = optfolder + "/precstat.txt"
            pr = open(defFN, "w")
            pr.write(datastring)
            pr.close()

            # ******************************************************************************************************
            # open "precstat" file in text editor and close selection dialog box
            # ******************************************************************************************************
            time.sleep(3)  # 5-24-2017: sleep time of few seconds to finish all processes before text file can be opened
            hydro.openbrowser(defFN)

            # ******************************************************************************************************
            # turn watershed delineation OFF and basin stat ON
            # ******************************************************************************************************
            button14.enabled = True
            arcpy.env.extent = "MAXOF"
            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()
            save(optfolder)

    def onClose(self, event):
        self.Show(False)

class PrecipitationDepths(object):
    """Implementation for GISHydroNXT_addin.button13 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        PrcpFrequency()


class Reservoir(object):
    """Implementation for GISHydroNXT_addin.tool8 (Tool)"""

    def __init__(self):
        self.enabled = False
        self.shape = 3

    def onMouseDownMap(self, x, y, button, shift):

        ### ADD RESERVOIR (STRUCTURE) FOR WINTR20 INPUT
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        xy = (x, y)
        arcpy.DeleteFeatures_management(optfolder + "/AddasReservoir.shp")
        cursor = arcpy.da.InsertCursor(optfolder + "/AddasReservoir.shp", ("SHAPE@XY"))
        cursor.insertRow([xy])

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()

        global rating
        in_target = optfolder + "/AddasReservoir.shp"
        in_join = optfolder + "/subshed.shp"
        out_feature_class = optfolder + "/aux_folder/intersect_aux.shp"
        arcpy.SpatialJoin_analysis(in_target, in_join, out_feature_class)
        sr = arcpy.SearchCursor(out_feature_class, "", "", "ARCID", "")
        for node in sr:
            rating = int(node.getValue("ARCID"))

        ### TRAP ERROR IN CASE TRANSECT IS IN INCORRECT REACH
        arcidnode = []
        fromnode = []
        tonode = []
        subriver_prop = arcpy.SearchCursor(optfolder + "/subrivers.shp", "", "", "ARCID;From_Node;To_Node", "")
        for node in subriver_prop:
            arcidnode.append(int(node.getValue("ARCID")))
            fromnode.append(int(node.getValue("From_Node")))
            tonode.append(int(node.getValue("To_Node")))
        subreach_lst = [fn for fn in fromnode if fn in tonode]
        arcid_list = []
        for sub in subreach_lst:
            index = fromnode.index(sub)
            arcid_list.append(arcidnode[index])
        line_trun = [optfolder + "/AddasReservoir.shp", optfolder + "/subshed.shp"]
        arcpy.Intersect_analysis(line_trun, optfolder + "/line_trun.shp", "ALL", "#", "INPUT")
        transect_prop = arcpy.SearchCursor(optfolder + "/line_trun.shp", "", "", "ARCID", "")
        transnodes = []
        for tra in transect_prop:
            transnodes.append(int(tra.getValue("ARCID")))
        condition_in = [tn for tn in transnodes if tn in arcid_list]

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "line_trun":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "intersect_aux":
                arcpy.mapping.RemoveLayer(df, lyr)
        arcpy.Delete_management(optfolder + "/line_trun.shp", "")
        arcpy.Delete_management(optfolder + "/intersect_aux.shp", "")

        if len(condition_in) == 0:
            pythonaddins.MessageBox("Structure must be placed within a routing reach subwatershed", "Reservoir Location Error")
            return

        global reservoir_elev
        reselev = arcpy.GetCellValue_management(optfolder + "/dem", '%s %s' % (x,y))
        reservoir_elev = str(round(float(reselev.getOutput(0)), 2))

        #Default name for reservoir
        global struc_name
        struc_name = "Struct" + str(rating)

        ReservoirFrame()
        save(optfolder)

class ReservoirFrame(wx.Frame):
    def __init__(self):

        wx.Frame.__init__(self, parent=None, title="Structure Rating", size=(420, 450))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.panel = wx.Panel(self)

        vbox = wx.BoxSizer(wx.VERTICAL)

        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        st1 = wx.StaticText(self.panel, label='Structure Identifier: ' + struc_name)
        hbox1.Add(st1, flag=wx.RIGHT, border=8)
        vbox.Add(hbox1, flag=wx.EXPAND|wx.LEFT|wx.RIGHT|wx.TOP, border=10)

        vbox.Add((-1, 10))

        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        st1 = wx.StaticText(self.panel, label='Reference Elevation: ' + reservoir_elev + " FASL")
        hbox1.Add(st1, flag=wx.RIGHT, border=8)
        vbox.Add(hbox1, flag=wx.EXPAND|wx.LEFT|wx.RIGHT|wx.TOP, border=10)

        vbox.Add((-1, 15))

        line1 = wx.BoxSizer(wx.HORIZONTAL)
        line = wx.StaticLine(self.panel)
        line1.Add(line, proportion=1)
        vbox.Add(line1, flag=wx.EXPAND|wx.LEFT|wx.RIGHT|wx.TOP, border=10)

        vbox.Add((-1, 25))

        self.mygrid = gridlib.Grid(self.panel)
        self.mygrid.CreateGrid(50, 3)

        for i in range(3):
            self.mygrid.SetColSize(i, 120)

        self.mygrid.SetColLabelValue(0, "Stage [ft]")
        self.mygrid.SetColLabelValue(1, "Discharge [cfs]")
        self.mygrid.SetColLabelValue(2, "Storage [ac ft]")

        self.mygrid.SetRowLabelSize(0)

        self.mygrid.SetColLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.mygrid.SetDefaultCellAlignment( wx.ALIGN_CENTRE, wx.ALIGN_TOP )
        # ******************************* #

        hbox2 = wx.BoxSizer(wx.VERTICAL)
        hbox2.Add(self.mygrid, proportion=1, flag=wx.EXPAND)
        vbox.Add(hbox2, proportion=1, flag=wx.LEFT|wx.RIGHT|wx.EXPAND, border=10)

        vbox.Add((-1, 25))

        hbox4 = wx.BoxSizer(wx.HORIZONTAL)
        btn4 = wx.Button(self.panel, label='Copy', size=(70, 30))
        hbox4.Add(btn4, flag=wx.LEFT|wx.BOTTOM, border=5)
        btn5 = wx.Button(self.panel, label='Paste', size=(70, 30))
        hbox4.Add(btn5, flag=wx.LEFT|wx.BOTTOM, border=5)
        btn6 = wx.Button(self.panel, label='Delete All', size=(70, 30))
        hbox4.Add(btn6, flag=wx.LEFT|wx.BOTTOM, border=5)
        vbox.Add(hbox4, flag=wx.ALIGN_RIGHT|wx.RIGHT, border=10)
        self.panel.SetSizer(vbox)

        vbox.Add((-1, 10))

        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        btn1 = wx.Button(self.panel, label='Add', size=(70, 30))
        hbox3.Add(btn1, flag=wx.LEFT|wx.BOTTOM, border=5)
        btn2 = wx.Button(self.panel, label='Display', size=(70, 30))
        hbox3.Add(btn2, flag=wx.LEFT|wx.BOTTOM, border=5)
        btn3 = wx.Button(self.panel, label='Cancel', size=(70, 30))
        hbox3.Add(btn3, flag=wx.LEFT|wx.BOTTOM, border=5)
        vbox.Add(hbox3, flag=wx.ALIGN_RIGHT|wx.RIGHT, border=10)
        self.panel.SetSizer(vbox)

        self.Bind(gridlib.EVT_GRID_CELL_CHANGE , self.onselectcell, self.mygrid)
        self.Bind(wx.EVT_BUTTON, self.OnSet, id=btn1.GetId())
        self.Bind(wx.EVT_BUTTON, self.OnPlot, id=btn2.GetId())
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=btn3.GetId())

        self.Bind(wx.EVT_BUTTON, self.copy, id=btn4.GetId())
        self.Bind(wx.EVT_BUTTON, self.paste, id=btn5.GetId())
        self.Bind(wx.EVT_BUTTON, self.delete, id=btn6.GetId())

        self.Show(True)
        self.Centre(True)

    def onselectcell(self, event):
        row = event.GetRow()
        col = event.GetCol()
        if self.mygrid.GetCellValue(row, col).isdigit() == False:
            self.mygrid.SetCellValue(row, col, "")
            event.Skip()

    def copy(self, event):
        # Number of rows and cols
        topleft = self.mygrid.GetSelectionBlockTopLeft()
        if list(topleft) == []:
            topleft = []
        else:
            topleft = list(topleft[0])
        bottomright = self.mygrid.GetSelectionBlockBottomRight()
        if list(bottomright) == []:
            bottomright = []
        else:
            bottomright = list(bottomright[0])
        if list(self.mygrid.GetSelectionBlockTopLeft()) == []:
            rows = 1
            cols = 1
            iscell = True
        else:
            rows = bottomright[0] - topleft[0] + 1
            cols = bottomright[1] - topleft[1] + 1
            iscell = False
        # data variable contain text that must be set in the clipboard
        data = ''
        # For each cell in selected range append the cell value in the data variable
        # Tabs '    ' for cols and '\r' for rows
        for r in range(rows):
            for c in range(cols):
                if iscell:
                    data += str(self.mygrid.GetCellValue(self.GetGridCursorRow() + r, self.mygrid.GetGridCursorCol() + c))
                else:
                    data += str(self.mygrid.GetCellValue(topleft[0] + r, topleft[1] + c))
                if c < cols - 1:
                    data += '    '
            data += '\n'
        # Create text data object
        clipboard = wx.TextDataObject()
        # Set data object value
        clipboard.SetText(data)
        # Put the data in the clipboard
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()
        else:
            wx.MessageBox("Can't open the clipboard", "Error")

    def paste(self, event):
        topleft = list(self.mygrid.GetSelectionBlockTopLeft())

        clipboard = wx.TextDataObject()
        if wx.TheClipboard.Open():
            wx.TheClipboard.GetData(clipboard)
            wx.TheClipboard.Close()
        else:
            wx.MessageBox("Can't open the clipboard", "Error")
        data = clipboard.GetText()
        if topleft == []:
            rowstart = self.mygrid.GetGridCursorRow()
            colstart = self.mygrid.GetGridCursorCol()
        else:
            rowstart = topleft[0][0]
            colstart = topleft[0][1]

        text4undo = ''
        # Convert text in a array of lines
        for y, r in enumerate(data.splitlines()):
            # Convert c in a array of text separated by tab
            for x, c in enumerate(r.split('    ')):
                if y + rowstart < self.mygrid.NumberRows and x + colstart < self.mygrid.NumberCols :
                    text4undo += str(self.mygrid.GetCellValue(rowstart + y, colstart + x)) + '    '
                    self.mygrid.SetCellValue(rowstart + y, colstart + x, c)
            text4undo = text4undo[:-1] + '\n'


    def delete(self, event):
        for i in range(50):
             for j in range(3):
                self.mygrid.SetCellValue(i, j, "")

    def OnClose(self, event):
        try:
            self.frm.Show(False)
        except:
            pass
        self.Show(False)

    def OnSet(self, event):

        stage = []
        discharge = []
        storage = []
        for i in range(50):
                stage.append(self.mygrid.GetCellValue(i, 0))
                discharge.append(self.mygrid.GetCellValue(i, 1))
                storage.append(self.mygrid.GetCellValue(i, 2))

        stage = filter(None, stage)
        discharge = filter(None, discharge)
        storage = filter(None, storage)

        try:
            discharge[0] = 0
            STA = [float(elem) for elem in stage]
            DIS = [float(elem) for elem in discharge]
            STO = [float(elem) for elem in storage]
        except:
            pythonaddins.MessageBox("Invalid values","Rating Error")

        if not all([len(STA)>2,len(DIS)>2,len(STO)>2]):
            pythonaddins.MessageBox("Rating table too short","Rating Error")
            self.Show(False)
            self.Show(True)
            return
        elif not all([len(DIS) == len(STA), len(DIS) == len(STO)]):
            pythonaddins.MessageBox("Columns are not the same size","Rating Error")
            self.Show(False)
            self.Show(True)
            return
        if not all(i < j for i, j in zip(STA, STA[1:])):
            pythonaddins.MessageBox("Stage values are not increasing","Rating Error")
            self.Show(False)
            self.Show(True)
            return
        elif not all(i < j for i, j in zip(DIS, DIS[1:])):
            pythonaddins.MessageBox("Discharge values are not increasing","Rating Error")
            self.Show(False)
            self.Show(True)
            return
        elif not all(i < j for i, j in zip(STO, STO[1:])):
            pythonaddins.MessageBox("Storage values are not increasing","Rating Error")
            self.Show(False)
            self.Show(True)
            return

        discharge[0] = 0
        STA = ['%.1f' % float(elem) for elem in STA]
        DIS = ['%.1f' % float(elem) for elem in DIS]
        STO = ['%.1f' % float(elem) for elem in STO]

        if not os.path.exists(optfolder + "/rating_table"):
            path = optfolder + "/rating_table"
            os.mkdir(path, 0o755)

        defFN = optfolder + "/rating_table/rattabout_reach" + str(rating) + ".txt"

        datastring = struc_name + "\n"
        datastring = datastring + STA[0] + "\n"
        datastring = datastring + "Reach" + str(rating) + "\n"
        datastring = datastring + "\n"
        datastring = datastring + '{0[0]:<12}{0[1]:<12}{0[2]:<12}'.format(
                ["Stage", "Discharge", "Storage"]) + "\n"
        datastring = datastring + '{0[0]:<12}{0[1]:<12}{0[2]:<12}'.format(
                ["-----", "-----", "-----"]) + "\n"
        datastring = datastring + "\n"

        rows = zip(STA, DIS, STO)
        for row in rows:
            datastring = datastring + '{0[0]:<12}{0[1]:<12}{0[2]:<12}'.format(row) + "\n"

        ratfile = open(defFN, "w")
        ratfile.write(datastring)
        ratfile.close()

        fn_lst = []
        tn_lst = []
        subriver = optfolder + "/subrivers.shp"
        sr = arcpy.SearchCursor(subriver, "", "", "FROM_NODE;TO_NODE", "")
        for node in sr:
            fn = node.getValue("FROM_NODE")
            tn = node.getValue("TO_NODE")
            fn_lst.append(fn)
            tn_lst.append(tn)

        # use reach numbers to name "userout" and "rattabout" files
        reach_lst = [x for x in fn_lst if x in tn_lst]

        hydro.openbrowser(defFN)

        # ******************************************************************************************************
        # turn layers ON/OFF in current data frame
        # ******************************************************************************************************
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "AddasReservoir":
                lyr.visible = False
            if lyr.name == "intersect_aux":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "/xline_reach" + str(rating):
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "/xpoint_reach" + str(rating):
                arcpy.mapping.RemoveLayer(df, lyr)
        if os.path.exists(optfolder + "/aux_folder/intersect_aux.shp"):
            arcpy.Delete_management(optfolder + "/aux_folder/intersect_aux.shp")
        if os.path.exists(optfolder + "/xline_reach" + str(rating) + ".shp"):
            arcpy.Delete_management(optfolder + "/xline_reach" + str(rating) + ".shp")
        if os.path.exists(optfolder + "/xpoint_reach" + str(rating) + ".shp"):
            arcpy.Delete_management(optfolder + "/xpoint_reach" + str(rating) + ".shp")

        # ******************************************************************************************************
        # turn Xeditor and reservoirs OFF and Precipitation Depths ON
        # ******************************************************************************************************
        myPath = optfolder + "/rating_table/"
        txtCounter = len(glob.glob1(myPath, "*.txt"))
        reachCounter = len(reach_lst)

        if txtCounter == reachCounter:
            user_prompt = pythonaddins.MessageBox("Process finished, continue?", "Reach Tool", 4)
            if user_prompt == "Yes":
                tool6.enabled = False
                tool8.enabled = False
                button13.enabled = True
                button25.enabled = True
                pythonaddins.MessageBox("Process finished!" + "\n" +
                                        "Continue to Precipitation Depths selection.", "Reach Tool")

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()
        self.Show(False)
        self.frm.Show(False)


    def OnPlot(self, event):

        stage = []
        discharge = []
        storage = []
        for i in range(50):
                stage.append(self.mygrid.GetCellValue(i, 0))
                discharge.append(self.mygrid.GetCellValue(i, 1))
                storage.append(self.mygrid.GetCellValue(i, 2))

        stage = []
        discharge = []
        storage = []
        for i in range(50):
                stage.append(self.mygrid.GetCellValue(i, 0))
                discharge.append(self.mygrid.GetCellValue(i, 1))
                storage.append(self.mygrid.GetCellValue(i, 2))

        stage = filter(None, stage)
        discharge = filter(None, discharge)
        storage = filter(None, storage)

        try:
            discharge[0] = 0
            STA = [float(elem) for elem in stage]
            DIS = [float(elem) for elem in discharge]
            STO = [float(elem) for elem in storage]
        except:
            pythonaddins.MessageBox("Invalid values","Rating Error")

        if not all([len(STA)>2,len(DIS)>2,len(STO)>2]):
            pythonaddins.MessageBox("Rating table too short","Rating Error")
            self.Show(False)
            self.Show(True)
            return
        elif not all([len(DIS) == len(STA), len(DIS) == len(STO)]):
            pythonaddins.MessageBox("Columns are not the same size","Rating Error")
            self.Show(False)
            self.Show(True)
            return
        if not all(i < j for i, j in zip(STA, STA[1:])):
            pythonaddins.MessageBox("Stage values are not increasing","Rating Error")
            self.Show(False)
            self.Show(True)
            return
        elif not all(i < j for i, j in zip(DIS, DIS[1:])):
            pythonaddins.MessageBox("Discharge values are not increasing","Rating Error")
            self.Show(False)
            self.Show(True)
            return
        elif not all(i < j for i, j in zip(STO, STO[1:])):
            pythonaddins.MessageBox("Storage values are not increasing","Rating Error")
            self.Show(False)
            self.Show(True)
            return

        discharge[0] = 0
        STA = [float(i) for i in stage]#map(int, stage)
        DIS = [float(i) for i in discharge]#map(int, discharge)
        STO = [float(i) for i in storage]#map(int, storage)

        datadis = zip(STA, DIS)
        datasto = zip(STA, STO)
        self.frm = wx.Frame(self, -1, "Structure Rating Curve", size=(600, 450))
        client = plot.PlotCanvas(self.frm)
        line1 = plot.PolyLine(datadis, legend="Discharge [cfs]", colour="blue", width=2)
        line2 = plot.PolyLine(datasto, legend="Storage [ac ft]", colour="red", width=2)
        gc = plot.PlotGraphics([line1,line2], "Structure: " + struc_name, "Elevation", "Flow - Storage")
        client.Draw(gc, xAxis=(min(STA), max(STA)), yAxis=(0, max([max(DIS),max(STO)])))
        client.SetEnableLegend(True)
        client.SetEnableGrid(True)
        self.frm.Show(True)


class ResetGlobal(object):
    """Implementation for GISHydroNXT_addin.button24 (button)"""

    def __init__(self):
        self.enabled = True
        self.checked = False

    def onClick(self):

        arcpy.env.addOutputsToMap = True

        restart = pythonaddins.MessageBox("This will restart GISHydroNXT, are you sure?", "Restart", 4)

        if restart == "Yes":

            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            tables = arcpy.mapping.ListTableViews(mxd, "", df)
            for lyr in layers:
                arcpy.mapping.RemoveLayer(df, lyr)
            for tbl in tables:
                arcpy.mapping.RemoveTableView(df, tbl)

            mxdloc = mxd.filePath
            mxdloc = mxdloc.replace("\\", "/")
            Directory = r"" + mxdloc.replace("mdinterface/GISHydroNXT.mxd", "/")

            if not os.path.exists(Directory):
                pythonaddins.MessageBox("Error: GISHydroNXT data folder is missing", "Load Error", 0)
                return

            dem = r"" + Directory + "/data/dems/neddem"
            arcpy.MakeRasterLayer_management(dem, "md dem")
            quads = r"" + Directory + "/data/maryland/quads83v3m.shp"
            arcpy.MakeFeatureLayer_management(quads, "md quads")
            roads = r"" + Directory + "/data/maryland/mjr-rdsstpm.shp"
            arcpy.MakeFeatureLayer_management(roads, "md roads")
            nhd = r"" + Directory + "/data/maryland/nhd_streamsm.shp"
            arcpy.MakeFeatureLayer_management(nhd, "nhd streams")

            layer = arcpy.mapping.Layer("md dem")
            arcpy.ApplySymbologyFromLayer_management(layer, r"" + Directory + "/data/mdfiles/legends/DEM_global.lyr")
            layer = arcpy.mapping.Layer("md quads")
            arcpy.ApplySymbologyFromLayer_management(layer, r"" + Directory + "/data/mdfiles/legends/quads.lyr")
            layer = arcpy.mapping.Layer("md roads")
            arcpy.ApplySymbologyFromLayer_management(layer, r"" + Directory + "/data/mdfiles/legends/roads.lyr")
            layer = arcpy.mapping.Layer("nhd streams")
            arcpy.ApplySymbologyFromLayer_management(layer, r"" + Directory + "/data/mdfiles/legends/nhd.lyr")

            tool1.enabled = True
            tool2.enabled = False
            tool3.enabled = False
            tool4.enabled = False
            tool5.enabled = False
            tool6.enabled = False
            tool7.enabled = False
            tool8.enabled = False
            button1.enabled = False
            button2.enabled = False
            button3.enabled = False
            button4.enabled = False
            button5_1.enabled = False
            button5.enabled = False
            button6.enabled = False
            button7.enabled = False
            button8.enabled = False
            button9.enabled = False
            button10.enabled = False
            button11.enabled = False
            button12.enabled = False
            button13.enabled = False
            button14.enabled = False
            button15.enabled = False
            button16.enabled = False
            button17.enabled = True
            button18.enabled = True
            button19.enabled = True
            button20.enabled = True
            button21.enabled = True
            button23.enabled = True
            button24.enabled = True
            button25.enabled = False

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()

        else:
            return

class ResetSubWatershed(object):
    """Implementation for GISHydroNXT_addin.button5 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        # *******************************************************************************************************
        # Remove layers and tables from data frame
        # *******************************************************************************************************

        restart = pythonaddins.MessageBox("This will restart CRWR-PrePro tools (including flow paths and outlets), are you sure?", "Restart", 4)

        if restart == "Yes":

            result_addas = arcpy.GetCount_management(optfolder + "/AddasStreams.shp")
            result_outlets = arcpy.GetCount_management(optfolder + "/AddasOutlets.shp")
            count_subsheds = int(result_addas.getOutput(0)) + int(result_outlets.getOutput(0))

            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            tables = arcpy.mapping.ListTableViews(mxd, "", df)
            for lyr in layers:
                if lyr.name == "watershed":
                    lyr.visible = True
                if lyr.name == "ModStreams":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subrivers":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "line":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "elevmerge":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "elevzones":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subshed":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "AddasStreams":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "AddasOutlets":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "AddasReservoir":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "polyras":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "Outlets_temp":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "line":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "modstr":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "outlets":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "strlnk":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "StrmMerge.shp":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subriver_subshed.shp":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subriver_xy.shp":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subshed_subriver.shp":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subshed_temp":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "subsheds":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "line":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "tmpsubwshd.shp":
                    arcpy.mapping.RemoveLayer(df, lyr)
                for i in range(0, 2 * count_subsheds, 1):
                    if lyr.name == "Longest_Path_Sub_" + str(i + 1):
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "xline_reach" + str(i + 1):
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "xpoint_reach" + str(i + 1):
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if os.path.exists(optfolder + "/aux_folder/mask" + str(i)):
                        arcpy.Delete_management(optfolder + "/aux_folder/mask" + str(i), "")
                    if os.path.exists(optfolder + "/xline_reach" + str(i + 1) + ".shp"):
                        arcpy.Delete_management(optfolder + "/xline_reach" + str(i + 1) + ".shp")
                    if os.path.exists(optfolder + "/xpoint_reach" + str(i + 1) + ".shp"):
                        arcpy.Delete_management(optfolder + "/xpoint_reach" + str(i + 1) + ".shp")
                    if os.path.exists(optfolder + "/aux_folder/longest_path_sub_" + str(i + 1)):
                        arcpy.Delete_management(optfolder + "/aux_folder/longest_path_sub_" + str(i + 1), "")

            for tbl in tables:
                if tbl.name == "cntable":
                    arcpy.mapping.RemoveTableView(df, tbl)
                if tbl.name == "longfp":
                    arcpy.mapping.RemoveTableView(df, tbl)
                if tbl.name == "slope_sheds":
                    arcpy.mapping.RemoveTableView(df, tbl)
                if tbl.name == "slope_stats":
                    arcpy.mapping.RemoveTableView(df, tbl)

            # *******************************************************************************************************
            # Delete files from optfolder
            # *******************************************************************************************************
            # Add Outlets file names
            # ********************************************
            outlets_tmp = optfolder + "/outlets_temp"
            addoutlets = optfolder + "/AddOutlets"
            outlets_usr = optfolder + "/outlets_user"
            addasstreams = optfolder + "/AddasStreams.shp"
            addasoutlets = optfolder + "/AddasOutlets.shp"
            line = optfolder + "/line.shp"

            if os.path.exists(outlets_tmp):
                arcpy.Delete_management(outlets_tmp, "")
            if os.path.exists(addoutlets):
                arcpy.Delete_management(addoutlets, "")
            if os.path.exists(outlets_usr):
                arcpy.Delete_management(outlets_usr, "")
            if os.path.exists(addasstreams):
                arcpy.Delete_management(addasstreams, "")
            if os.path.exists(addasstreams):
                arcpy.Delete_management(addasoutlets, "")
            if os.path.exists(line):
                arcpy.Delete_management(addasoutlets, "")

            # Add Streams file names
            # ********************************************
            strm_merge = optfolder + "/StrmMerge.shp"
            modstr = optfolder + "/modstr"
            modstreams = optfolder + "/modstreams"

            if os.path.exists(strm_merge):
                arcpy.Delete_management(strm_merge, "")
            if os.path.exists(modstr):
                arcpy.Delete_management(modstr, "")
            if os.path.exists(modstreams):
                arcpy.Delete_management(modstreams, "")

            # Delineate Subwatersheds file names
            # ********************************************
            strlnk = optfolder + "/strlnk"
            outlets = optfolder + "/outlets"
            outlets_str = optfolder + "/outlets_str"
            added_outl = optfolder + "/added_outlets"
            subshed_tmp = optfolder + "/subshed_temp"
            tmpsubw = optfolder + "/tmpsubwshd.shp"
            sub_shape = optfolder + "/subshed.shp"
            new_strlnk = optfolder + "/newstrlnk"
            subrivers = optfolder + "/subrivers.shp"
            subriv_sub = optfolder + "/subriver_subshed.shp"
            subriv_xy = optfolder + "/subriver_xy.shp"
            sub_subriv = optfolder + "/subshed_subriver.shp"
            sub_ras = optfolder + "/subsheds"

            if os.path.exists(strlnk):
                arcpy.Delete_management(strlnk, "")
            if os.path.exists(outlets):
                arcpy.Delete_management(outlets, "")
            if os.path.exists(outlets_str):
                arcpy.Delete_management(outlets_str, "")
            if os.path.exists(added_outl):
                arcpy.Delete_management(added_outl, "")
            if os.path.exists(subshed_tmp):
                arcpy.Delete_management(subshed_tmp, "")
            if os.path.exists(tmpsubw):
                arcpy.Delete_management(tmpsubw, "")
            if os.path.exists(sub_shape):
                arcpy.Delete_management(sub_shape, "")
            if os.path.exists(new_strlnk):
                arcpy.Delete_management(new_strlnk, "")
            if os.path.exists(subrivers):
                arcpy.Delete_management(subrivers, "")
            if os.path.exists(subriv_sub):
                arcpy.Delete_management(subriv_sub, "")
            if os.path.exists(subriv_xy):
                arcpy.Delete_management(subriv_xy, "")
            if os.path.exists(sub_subriv):
                arcpy.Delete_management(sub_subriv, "")
            if os.path.exists(sub_ras):
                arcpy.Delete_management(sub_ras, "")

            # "Calculate Attributes" file names
            # ********************************************
            if os.path.exists(optfolder + "/elevmerge.shp"):
                for i in range(0, no_subwatersheds, 1):
                    if os.path.exists(optfolder + "/Tc_subshed" + str(i)):
                        arcpy.Delete_management(optfolder + "/Tc_subshed" + str(i))
                    if os.path.exists(optfolder + "/vel_meth/LongPathSub" + str(i)):
                        arcpy.Delete_management(optfolder + "/vel_meth/LongPathSub" + str(i))

            polyras = optfolder + "/polyras"
            elevmerge = optfolder + "/elevmerge.shp"
            elevzones = optfolder + "/elevzones.shp"
            sub_prov = optfolder + "/subshed_prov.shp"
            cntable = optfolder + "/cntable.dbf"
            longfp = optfolder + "/longfp.dbf"
            slope_shed = optfolder + "/slope_sheds.dbf"
            slope_stat = optfolder + "/slope_stats.dbf"
            line = optfolder + "/line.shp"

            if os.path.exists(polyras):
                arcpy.Delete_management(polyras, "")
            if os.path.exists(elevmerge):
                arcpy.Delete_management(elevmerge, "")
            if os.path.exists(elevzones):
                arcpy.Delete_management(elevzones, "")
            if os.path.exists(sub_prov):
                arcpy.Delete_management(sub_prov, "")
            if os.path.exists(cntable):
                arcpy.Delete_management(cntable, "")
            if os.path.exists(longfp):
                arcpy.Delete_management(longfp, "")
            if os.path.exists(slope_shed):
                arcpy.Delete_management(slope_shed, "")
            if os.path.exists(slope_stat):
                arcpy.Delete_management(slope_stat, "")
            if os.path.exists(line):
                arcpy.Delete_management(line, "")

            # Deleting "vel_meth" folder after deleting "LongPathSub" files above. Edited on 5/24/2017 to delete vel_meth folder
            # which is under schema lock

            # Delete "vel_meth" folder
            dst = optfolder + "/vel_meth"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            # Delete "sub_basincomp" folder
            # ********************************************
            dst = optfolder + "/sub_basincomp"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            # Transect layers
            linetrum = optfolder + "/line_trun.shp"
            if os.path.exists(linetrum):
                arcpy.Delete_management(linetrum, "")

            # Delete "rating_table" and "elev_stage_profile" folder
            # ******************************************************
            dst = optfolder + "/rating_table"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            dst = optfolder + "/elev_stage_profile"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            # Delete WinTR20 folder
            dst = optfolder + "/WinTR20"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            # Delete Design Storm folder
            dst = optfolder + "/design_storm"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass
            # *******************************************************************************************************
            # Re-create "AddasStreams.shp" for to be used in subwatershed delineation
            # (this shapefile was created earlier during data selection step but reset
            #  button will delete it. Here it is re-created to store stream lines on run)
            # *******************************************************************************************************
            arcpy.env.addOutputsToMap = True
            spatial_reference = arcpy.Describe(optfolder + "/flowdir").spatialReference
            arcpy.CreateFeatureclass_management(optfolder, "AddasStreams.shp", "POLYLINE", "", "ENABLED", "DISABLED",
                                                spatial_reference)

            spatial_reference = arcpy.Describe(optfolder + "/flowdir").spatialReference
            arcpy.CreateFeatureclass_management(optfolder, "AddasOutlets.shp", "POINT", "", "ENABLED", "DISABLED",
                                                spatial_reference)

            try:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                df = arcpy.mapping.ListDataFrames(mxd, "Layers")[0]
                lyr = arcpy.mapping.ListLayers(mxd, "watershed", df)[0]
                df.extent = lyr.getSelectedExtent()
            except:
                pass

            # Remove precstat.txt and TR20in.txt
            # ********************************************
            precstat = optfolder + "/precstat.txt"
            TR20in = optfolder + "/TR20in.txt"
            if os.path.exists(precstat):
                arcpy.Delete_management(precstat, "")
            if os.path.exists(TR20in):
                arcpy.Delete_management(TR20in, "")

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()

            # *******************************************************************************************************
            # turn watershed delineation OFF and basin stat ON
            # *******************************************************************************************************
            button1.enabled = True
            tool6.enabled = False
            tool8.enabled = False
            tool5.enabled = False
            tool4.enabled = True
            button6.enabled = False
            button5.enabled = False
            button8.enabled = False
            button9.enabled = False
            button10.enabled = False
            button11.enabled = False
            button13.enabled = False
            button14.enabled = False
            button15.enabled = False
            save(optfolder)
        else:
            return


class ResetControlPanel(object):
    """Implementation for GISHydroNXT_addin.button5 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        # *******************************************************************************************************
        # Remove layers and tables from data frame
        # *******************************************************************************************************
        restart = pythonaddins.MessageBox("This will restart the WinTR20 reaches (transects and structures), are you sure?", "Restart", 4)

        if restart == "Yes":

            result_addas = arcpy.GetCount_management(optfolder + "/AddasStreams.shp")
            result_outlets = arcpy.GetCount_management(optfolder + "/AddasOutlets.shp")
            count_subsheds = int(result_addas.getOutput(0)) + int(result_outlets.getOutput(0))

            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "line":
                    arcpy.mapping.RemoveLayer(df, lyr)
                for i in range(0, 2 * count_subsheds, 1):
                    if lyr.name == "xline_reach" + str(i + 1):
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "xpoint_reach" + str(i + 1):
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if os.path.exists(optfolder + "/xline_reach" + str(i + 1) + ".shp"):
                        arcpy.Delete_management(optfolder + "/xline_reach" + str(i + 1) + ".shp")
                    if os.path.exists(optfolder + "/xpoint_reach" + str(i + 1) + ".shp"):
                        arcpy.Delete_management(optfolder + "/xpoint_reach" + str(i + 1) + ".shp")

            # *******************************************************************************************************
            # Delete files from optfolder
            # *******************************************************************************************************

            line = optfolder + "/line.shp"
            if os.path.exists(line):
                arcpy.Delete_management(line, "")

            # Transect layers
            linetrum = optfolder + "/line_trun.shp"
            if os.path.exists(linetrum):
                arcpy.Delete_management(linetrum, "")

            # Delete "rating_table" and "elev_stage_profile" folder
            # ******************************************************
            dst = optfolder + "/rating_table"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            dst = optfolder + "/elev_stage_profile"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            # Delete WinTR20 folder
            dst = optfolder + "/WinTR20"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            # Delete Design Storm folder
            dst = optfolder + "/design_storm"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            # Remove precstat.txt and TR20in.txt
            # ********************************************
            precstat = optfolder + "/precstat.txt"
            if os.path.exists(precstat):
                arcpy.Delete_management(precstat, "")

            # *******************************************************************************************************
            # turn watershed delineation OFF and basin stat ON
            # *******************************************************************************************************
            tool6.enabled = True
            tool8.enabled = True
            button13.enabled = False
            button14.enabled = False
            button15.enabled = False
            button25.enabled = False
            save(optfolder)
        else:
            return

class ResetWatershed(object):
    """Implementation for GISHydroNXT_addin.button1 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):

        restart = pythonaddins.MessageBox("This will restart the Hydro tools (including the watershed), are you sure?", "Restart", 4)

        if restart == "Yes":
            arcpy.env.scratchWorkspace = scratchfolder
            arcpy.env.workspace = optfolder

            global landedit
            landedit = []
            # *******************************************************************************************************
            # Remove watershed layer from data frame and Delete all temp files pertaining to watershed delineation
            # *******************************************************************************************************
            if os.path.exists(optfolder + "/aux_folder/curvenum_aux"):
                arcpy.env.addOutputsToMap = False
                arcpy.Delete_management(optfolder + "curveNumber", "")
                arcpy.CopyRaster_management(optfolder + "/aux_folder/curvenum_aux", optfolder + "curveNumber")
                arcpy.Delete_management(optfolder + "/aux_folder/curvenum_aux", "")
            if os.path.exists(optfolder + "/aux_folder/landuse"):
                arcpy.env.addOutputsToMap = True
                arcpy.Delete_management(optfolder + "landuse", "")
                arcpy.CopyRaster_management(optfolder + "/aux_folder/landuse", optfolder + "landuse")
                arcpy.Delete_management(optfolder + "/aux_folder/landuse", "")
            arcpy.env.addOutputsToMap = True
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            tables = arcpy.mapping.ListTableViews(mxd, "", df)
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "watershed":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "landuse":
                    lyr.visible = False
                if lyr.name == "outletpoint":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "outletpoly":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "flowdir":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "outlet_ws":
                    arcpy.mapping.RemoveLayer(df, lyr)

            for tbl in tables:
                if tbl.name == "thevTab":
                    arcpy.mapping.RemoveTableView(df, tbl)

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()

            # *******************************************************************************************************
            # Delete files from optfolder
            # *******************************************************************************************************
            # Basin Stat file names
            # ********************************************
            basingrid = optfolder + "/basingrid"
            watershed = optfolder + "/watershed.shp"
            maxlength = optfolder + "/maxlength"
            outcell = optfolder + "/outcell"
            outletpoly = optfolder + "/outletpoly.shp"
            outletpoint = optfolder + "/outletpoint.shp"
            thevtab = optfolder + "/thevTab.dbf"
            limegrid = optfolder + "/limegrid"
            wats_prov = optfolder + "/wats_prov.shp"
            prov_int = optfolder + "/prov_int.shp"
            wats_lime = optfolder + "/wats_lime.shp"
            lime_int = optfolder + "/lime_int.shp"
            basinstat = optfolder + "/basinstat.txt"
            landslope = optfolder + "/slope_calc/landslope"  # can"t delete it with "slope_calc" folder so a separate delete here
            old_lime = optfolder + "/limegrid_old.shp"
            flowdir = optfolder + "/flowdir"
            outlet_ws = optfolder + "/outlet_ws.shp"
            outletcell = optfolder + "/outletcell"
            contours = optfolder + "/contours"

            if os.path.exists(contours):
                arcpy.Delete_management(contours, "")
            if os.path.exists(flowdir):
                arcpy.Delete_management(flowdir, "")
            if os.path.exists(outletcell):
                arcpy.Delete_management(outletcell, "")
            if os.path.exists(outlet_ws):
                arcpy.Delete_management(outlet_ws, "")
            if os.path.exists(basingrid):
                arcpy.Delete_management(basingrid, "")
            if os.path.exists(watershed):
                arcpy.Delete_management(watershed, "")
            if os.path.exists(maxlength):
                arcpy.Delete_management(maxlength, "")
            if os.path.exists(outcell):
                arcpy.Delete_management(outcell, "")
            if os.path.exists(outletpoly):
                arcpy.Delete_management(outletpoly, "")
            if os.path.exists(outletpoint):
                arcpy.Delete_management(outletpoint, "")
            if os.path.exists(thevtab):
                arcpy.Delete_management(thevtab, "")
            if os.path.exists(limegrid):
                arcpy.Delete_management(limegrid, "")
            if os.path.exists(wats_prov):
                arcpy.Delete_management(wats_prov, "")
            if os.path.exists(prov_int):
                arcpy.Delete_management(prov_int, "")
            if os.path.exists(wats_lime):
                arcpy.Delete_management(wats_lime, "")
            if os.path.exists(lime_int):
                arcpy.Delete_management(lime_int, "")
            if os.path.exists(basinstat):
                arcpy.Delete_management(basinstat, "")
            if os.path.exists(landslope):
                arcpy.Delete_management(landslope, "")
            if os.path.exists(landslope):
                arcpy.Delete_management(old_lime, "")
            if os.path.exists(thevtab):
                arcpy.Delete_management(thevtab, "")

            # Delete "slope_calc" folder
            dst = optfolder + "/slope_calc"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            # Delete "tasker" folder
            dst = optfolder + "/tasker"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            # Thomas Discharge file names
            # ********************************************
            outletbuff = optfolder + "/outletbuffer.shp"
            gagefield = optfolder + "/gagefield.shp"
            thomasdisc = optfolder + "/frdischarges.txt"
            mask_ints = optfolder + "/mask_ints.shp"
            gauge_ol = optfolder + "/gauge_outlet.shp"

            if os.path.exists(outletbuff):
                arcpy.Delete_management(outletbuff, "")
            if os.path.exists(gagefield):
                arcpy.Delete_management(gagefield, "")
            if os.path.exists(thomasdisc):
                arcpy.Delete_management(thomasdisc, "")
            if os.path.exists(mask_ints):
                arcpy.Delete_management(mask_ints, "")
            if os.path.exists(gauge_ol):
                arcpy.Delete_management(gauge_ol, "")

            # Delete "basincomp" folder
            dst = optfolder + "/basincomp"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            lyr = arcpy.mapping.ListLayers(mxd, "dem", df)[0]
            df.extent = lyr.getSelectedExtent()

            # *******************************************************************************************************
            # turn watershed delineation OFF and basin stat ON
            # *******************************************************************************************************
            tool2.enabled = True
            tool7.enabled = True
            tool3.enabled = False
            button1.enabled = False
            button2.enabled = False
            button3.enabled = False
            button4.enabled = False
            button5_1.enabled = False
            button6.enabled = False
            tool4.enabled = False

            save(optfolder)
        else:
            return

class SegMergeDialog(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Velocity Method Segment Generator", size=(580, 500))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)

        self.SetPosition((30, 30))

        text1 = wx.StaticText(panel, -1, "Select Sub-Area", (30, 25))
        text1.SetForegroundColour((255, 0, 0))  # set text color

        self.subID = wx.ComboBox(panel, -1, "", (55, 45), (50, -1), arcid_list, wx.CB_DROPDOWN)

        wx.StaticBox(panel, -1, "Create/Update Segment", (15, 90), size=(250, 340))
        wx.StaticBox(panel, -1, "Quick Merge", (30, 120), size=(220, 120))
        wx.StaticText(panel, -1, "Single Overland", (100, 150))
        self.cb1 = wx.CheckBox(panel, -1, "", pos=(70, 150))
        wx.StaticText(panel, -1, "Single Swale", (100, 180))
        self.cb2 = wx.CheckBox(panel, -1, "", pos=(70, 180))
        wx.StaticText(panel, -1, "Single Channel", (100, 210))
        self.cb3 = wx.CheckBox(panel, -1, "", pos=(70, 210))

        wx.StaticBox(panel, -1, "Merge Specific Segment", (30, 250), size=(220, 110))
        wx.StaticText(panel, -1, "Upstream Pixel #", (40, 280))
        self.uspixel = wx.TextCtrl(panel, -1, value="", pos=(170, 275), size=(65, 25))
        wx.StaticText(panel, -1, "Downstream Pixel #", (40, 315))
        self.dspixel = wx.TextCtrl(panel, -1, value="", pos=(170, 310), size=(65, 25))

        wx.StaticBox(panel, -1, "Velocity Method Statistics", (280, 40), size=(260, 340))
        wx.StaticText(panel, -1, "Sub-Area #", (290, 70))
        self.sanumber = wx.TextCtrl(panel, -1, value="", pos=(455, 65), size=(75, 25),
                                    style=wx.TE_READONLY)  # Sub-area number
        wx.StaticText(panel, -1, "Overall Travel Time (hrs):", (290, 110))
        self.Overall_Tc = wx.TextCtrl(panel, -1, value="", pos=(455, 105), size=(75, 25),
                                      style=wx.TE_READONLY)  # overall Tc hrs
        wx.StaticText(panel, -1, "Overland Travel Time (hrs):", (290, 150))
        self.over_Tc = wx.TextCtrl(panel, -1, value="", pos=(455, 145), size=(75, 25),
                                   style=wx.TE_READONLY)  # overland Tc hrs
        wx.StaticText(panel, -1, "Swale Travel Time (hrs):", (290, 190))
        self.swale_Tc = wx.TextCtrl(panel, -1, value="", pos=(455, 185), size=(75, 25),
                                    style=wx.TE_READONLY)  # swale Tc hrs
        wx.StaticText(panel, -1, "Channel Travel Time (hrs):", (290, 230))
        self.chan_Tc = wx.TextCtrl(panel, -1, value="", pos=(455, 225), size=(75, 25),
                                   style=wx.TE_READONLY)  # channel Tc hrs
        wx.StaticText(panel, -1, "# Overland Segments:", (290, 270))
        self.over_seg = wx.TextCtrl(panel, -1, value="", pos=(455, 265), size=(75, 25),
                                    style=wx.TE_READONLY)  # overland segments
        wx.StaticText(panel, -1, "# Swale Segments:", (290, 310))
        self.swale_seg = wx.TextCtrl(panel, -1, value="", pos=(455, 305), size=(75, 25),
                                     style=wx.TE_READONLY)  # swale segments
        wx.StaticText(panel, -1, "# Channel Segments:", (290, 350))
        self.chan_seg = wx.TextCtrl(panel, -1, value="", pos=(455, 345), size=(75, 25),
                                    style=wx.TE_READONLY)  # channel segments

        self.btnApply = wx.Button(panel, label="Recalculate Tc", pos=(90, 380))
        self.Bind(wx.EVT_BUTTON, self.RecalculateTc, id=self.btnApply.GetId())

        self.btnApply = wx.Button(panel, label="Close Dialog", pos=(450, 400))
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=self.btnApply.GetId())

        self.btnApply = wx.Button(panel, label="Apply", pos=(140, 45))
        self.Bind(wx.EVT_BUTTON, self.OnApply, id=self.btnApply.GetId())

        self.Show(True)
        self.Centre(True)
        style = self.GetWindowStyle()
        self.SetWindowStyle(style | wx.STAY_ON_TOP)

        global notable
        notable = np.ones(int(arcid))

    def RecalculateTc(self, event):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder

        global attlist
        single_overland = 0
        single_swale = 0
        single_channel = 0

        if self.cb1.IsChecked():
            single_overland = 1
        if self.cb2.IsChecked():
            single_swale = 1
        if self.cb3.IsChecked():
            single_channel = 1

        UPPX = self.uspixel.GetValue()
        DWNPX = self.dspixel.GetValue()
        if UPPX == "" and DWNPX == "":
            UPPX = -1
            DWNPX = -1
        else:
            UPPX = int(UPPX)
            DWNPX = int(DWNPX)

        channeln = float(Tc_nc)
        overland_n = float(Tc_ns)
        w_coef = float(Tc_cwCoef)
        w_exp = float(Tc_cwExp)
        d_coef = float(Tc_cdCoef)
        d_exp = float(Tc_cdExp)
        a_coef = float(Tc_caCoef)
        a_exp = float(Tc_caExp)
        P2 = float(Tc_P)
        l_tc = float(Tc_L)

        attlist_merge = [list(elem) for elem in attlist]
        attlist = attlist_merge
        attlist1 = attlist_merge
        attlist2 = attlist_merge
        attlist3 = attlist_merge
        overland_list = [i for i in attlist1 if i[4] == "overland"]
        swale_list = [i for i in attlist2 if i[4] == "swale"]
        channel_list = [i for i in attlist3 if i[4] == "channel"]

        if len(overland_list) != 0 and len(swale_list) != 0 and UPPX < overland_list[-1][5] and DWNPX > swale_list[0][
            2] and UPPX != -1:
            pythonaddins.MessageBox(
                "Up and Down pixel flow types are different. Please select up and down pixels of same type.",
                "Segment Merge Error", 0)
            return
        if len(channel_list) != 0 and len(swale_list) != 0 and UPPX < swale_list[-1][5] and DWNPX > channel_list[0][
            2] and UPPX != -1:
            pythonaddins.MessageBox(
                "Up and Down pixel flow types are different. Please select up and down pixels of same type.",
                "Segment Merge Error", 0)
            return

        if single_overland == 1 and len(overland_list) != 0:
            shape = overland_list[0][1]
            uppxl = float(overland_list[0][2])
            Type = overland_list[0][4]
            dwpxl = float(overland_list[-1][5])
            avgarea = 0
            for col in overland_list:
                avgarea = (float(col[5]) - float(col[2])) * float(col[6]) + avgarea
            avgarea = avgarea / (float(overland_list[-1][5]) - float(overland_list[0][2]))
            dsarea = np.sum(map(float, [col[7] for col in overland_list]))
            upelev = float(overland_list[0][8])
            dwnelev = float(overland_list[-1][9])
            ilength = np.sum(map(float, [col[14] for col in overland_list]))
            slope = (upelev - dwnelev) / ilength
            width = np.mean(map(float, [col[11] for col in overland_list]))
            depth = np.mean(map(float, [col[12] for col in overland_list]))
            xarea = np.mean(map(float, [col[13] for col in overland_list]))
            totlength = ilength
            if ilength > l_tc:
                iswale = ilength - l_tc
                vel_s = 16.1345 * (slope ** 0.5)
                if Tc_unpaved:
                    vel_s = 20.3282 * (slope ** 0.5)
                swaletinc = iswale / (vel_s * 3600)
                overtinc = 0.007 * ((overland_n * l_tc) ** 0.8) / ((P2 ** 0.5) * (slope ** 0.4))
                itime = overtinc + swaletinc
                vel = ilength / (itime * 3600)
            else:
                itime = 0.007 * ((overland_n * ilength) ** 0.8) / ((P2 ** 0.5) * (slope ** 0.4))
                vel = ilength / (itime * 3600)
            tottime = itime
            if np.isnan(vel) or vel == 0:
                vel = 0.001
            overland_list = [[[], shape, uppxl, [], Type, dwpxl, avgarea, dsarea, upelev, dwnelev,
                              slope, width, depth, xarea, ilength, totlength, vel, itime, tottime]]

        if single_swale == 1 and len(swale_list) != 0:
            shape = swale_list[0][1]
            uppxl = float(swale_list[0][2])
            Type = swale_list[0][4]
            dwpxl = float(swale_list[-1][5])
            avgarea = 0
            for col in swale_list:
                avgarea = (float(col[5]) - float(col[2])) * float(col[6]) + avgarea
            avgarea = avgarea / (float(swale_list[-1][5]) - float(swale_list[0][2]))
            dsarea = np.sum(map(float, [col[7] for col in swale_list]))
            upelev = float(swale_list[0][8])
            dwnelev = float(swale_list[-1][9])
            ilength = np.sum(map(float, [col[14] for col in swale_list]))
            slope = (upelev - dwnelev) / ilength
            width = np.mean(map(float, [col[11] for col in swale_list]))
            depth = np.mean(map(float, [col[12] for col in swale_list]))
            xarea = np.mean(map(float, [col[13] for col in swale_list]))
            if len(overland_list) != 0:
                totlength = ilength + float(overland_list[-1][15])
            else:
                totlength = ilength
            vel = 16.1345 * slope ** 0.5
            if Tc_unpaved:
                vel = 20.3282 * slope ** 0.5
            if np.isnan(vel) or vel == 0:
                vel = 0.001
            itime = ilength / vel / 3600
            if len(overland_list) != 0:
                tottime = itime + float(overland_list[-1][18])
            else:
                tottime = itime
            swale_list = [[[], shape, uppxl, [], Type, dwpxl, avgarea, dsarea, upelev, dwnelev,
                           slope, width, depth, xarea, ilength, totlength, vel, itime, tottime]]

        if single_swale == 0 and len(swale_list) != 0 and UPPX >= float(
                swale_list[0][2]) and DWNPX > UPPX and DWNPX <= float(swale_list[-1][5]):
            swale_list_1 = [i for i in swale_list if float(i[2]) < UPPX]
            swale_list_2 = [i for i in swale_list if float(i[2]) >= UPPX and float(i[5]) <= DWNPX]
            swale_list_3 = [i for i in swale_list if float(i[2]) >= DWNPX]

            shape = swale_list_2[0][1]
            uppxl = float(swale_list_2[0][2])
            Type = swale_list_2[0][4]
            dwpxl = float(swale_list_2[-1][5])
            avgarea = 0
            for col in swale_list_2:
                avgarea = (float(col[5]) - float(col[2])) * float(col[6]) + avgarea
            avgarea = avgarea / (float(swale_list_2[-1][5]) - float(swale_list_2[0][2]))
            dsarea = np.sum(map(float, [col[7] for col in swale_list_2]))
            upelev = float(swale_list_2[0][8])
            dwnelev = float(swale_list_2[-1][9])
            ilength = np.sum(map(float, [col[14] for col in swale_list_2]))
            slope = (upelev - dwnelev) / ilength
            width = np.mean(map(float, [col[11] for col in swale_list_2]))
            depth = np.mean(map(float, [col[12] for col in swale_list_2]))
            xarea = np.mean(map(float, [col[13] for col in swale_list_2]))
            if UPPX == float(swale_list[0][2]):
                if len(overland_list) != 0:
                    totlength = ilength + float(overland_list[-1][15])
                else:
                    totlength = ilength
            else:
                totlength = ilength + float(swale_list_1[-1][15])
            vel = 16.1345 * (slope ** 0.5)
            if Tc_unpaved:
                vel = 20.3282 * (slope ** 0.5)
            if np.isnan(vel) or vel == 0:
                vel = 0.001
            itime = ilength / vel / 3600
            if UPPX == float(swale_list[0][2]):
                if len(overland_list) != 0:
                    tottime = itime + float(overland_list[-1][18])
                else:
                    tottime = itime
            else:
                tottime = itime + float(swale_list_1[-1][18])

            swale_list_2 = [[[], shape, uppxl, [], Type, dwpxl, avgarea, dsarea, upelev, dwnelev,
                             slope, width, depth, xarea, ilength, totlength, vel, itime, tottime]]

            for i in range(len(swale_list_3)):
                shape = swale_list_3[i][1]
                uppxl = float(swale_list_3[i][2])
                Type = swale_list_3[i][4]
                dwpxl = swale_list_3[i][5]
                avgarea = swale_list_3[i][6]
                dsarea = swale_list_3[i][7]
                upelev = swale_list_3[i][8]
                dwnelev = swale_list_3[i][9]
                ilength = swale_list_3[i][14]
                slope = swale_list_3[i][10]
                width = swale_list_3[i][11]
                depth = swale_list_3[i][12]
                xarea = swale_list_3[i][13]
                if i == 0:
                    totlength = float(ilength) + float(swale_list_2[-1][15])
                else:
                    totlength = float(ilength) + float(swale_list_3[i - 1][15])
                vel = swale_list_3[i][16]
                if np.isnan(vel) or vel == 0:
                    vel = 0.001
                itime = swale_list_3[i][17]
                if i == 0:
                    tottime = float(itime) + float(swale_list_2[-1][18])
                else:
                    tottime = float(itime) + float(swale_list_3[i - 1][18])
                swale_list_3[i] = [[], shape, uppxl, [], Type, dwpxl, avgarea, dsarea, upelev, dwnelev,
                                   slope, width, depth, xarea, ilength, totlength, vel, itime, tottime]

            swale_list = swale_list_1 + swale_list_2 + swale_list_3

        if single_channel == 1 and len(channel_list) != 0:
            shape = channel_list[0][1]
            uppxl = float(channel_list[0][2])
            Type = channel_list[0][4]
            dwpxl = float(channel_list[-1][5])
            avgarea = 0
            for col in channel_list:
                avgarea = (float(col[5]) - float(col[2])) * float(col[6]) + avgarea
            avgarea = avgarea / (float(channel_list[-1][5]) - float(channel_list[0][2]))
            dsarea = np.sum(map(float, [col[7] for col in channel_list]))
            upelev = float(channel_list[0][8])
            dwnelev = float(channel_list[-1][9])
            ilength = np.sum(map(float, [col[14] for col in channel_list]))
            slope = (upelev - dwnelev) / ilength
            width = w_coef * (avgarea ** w_exp)
            depth = d_coef * (avgarea ** d_exp)
            xarea = a_coef * (avgarea ** a_exp)
            if len(swale_list) != 0:
                totlength = ilength + float(swale_list[-1][15])
            elif len(overland_list) != 0:
                totlength = ilength + float(overland_list[-1][15])
            else:
                totlength = ilength
            vel = 1.49 / channeln * (xarea / (width + 2 * depth)) ** (0.66667) * (slope ** 0.5)
            if np.isnan(vel) or vel == 0:
                vel = 0.001
            itime = ilength / vel / 3600
            if len(swale_list) != 0:
                tottime = itime + float(swale_list[-1][18])
            elif len(overland_list) != 0:
                tottime = itime + float(overland_list[-1][18])
            else:
                tottime = itime

            channel_list = [[[], shape, uppxl, [], Type, dwpxl, avgarea, dsarea, upelev, dwnelev,
                             slope, width, depth, xarea, ilength, totlength, vel, itime, tottime]]

        if single_channel == 0 and len(channel_list) != 0 and UPPX >= float(
                channel_list[0][2]) and DWNPX > UPPX and DWNPX <= float(channel_list[-1][5]):
            channel_list_1 = [i for i in channel_list if float(i[2]) < UPPX]
            channel_list_2 = [i for i in channel_list if float(i[2]) >= UPPX and float(i[5]) <= DWNPX]
            channel_list_3 = [i for i in channel_list if float(i[2]) >= DWNPX]

            shape = channel_list_2[0][1]
            uppxl = float(channel_list_2[0][2])
            Type = channel_list_2[0][4]
            dwpxl = float(channel_list_2[-1][5])
            avgarea = 0
            for col in channel_list_2:
                avgarea = (float(col[5]) - float(col[2])) * float(col[6]) + avgarea
            avgarea = avgarea / (float(channel_list_2[-1][5]) - float(channel_list_2[0][2]))
            dsarea = np.sum(map(float, [col[7] for col in channel_list_2]))
            upelev = float(channel_list_2[0][8])
            dwnelev = float(channel_list_2[-1][9])
            ilength = np.sum(map(float, [col[14] for col in channel_list_2]))
            slope = (upelev - dwnelev) / ilength
            width = w_coef * (avgarea ** w_exp)
            depth = d_coef * (avgarea ** d_exp)
            xarea = a_coef * (avgarea ** a_exp)
            if UPPX == float(channel_list[0][2]):
                if len(swale_list) != 0:
                    totlength = ilength + float(swale_list[-1][15])
                elif len(overland_list) != 0:
                    totlength = ilength + float(overland_list[-1][15])
                else:
                    totlength = ilength
            else:
                totlength = ilength + float(channel_list_1[-1][15])
            vel = 1.49 / channeln * (xarea / (width + 2 * depth)) ** (0.66667) * (slope ** 0.5)
            if np.isnan(vel) or vel == 0:
                vel = 0.001
            itime = ilength / vel / 3600
            if UPPX == float(channel_list[0][2]):
                if len(swale_list) != 0:
                    tottime = itime + float(swale_list[-1][18])
                elif len(overland_list) != 0:
                    tottime = itime + float(overland_list[-1][18])
                else:
                    tottime = itime
            else:
                tottime = itime + float(channel_list_1[-1][18])

            channel_list_2 = [[[], shape, uppxl, [], Type, dwpxl, avgarea, dsarea, upelev, dwnelev,
                               slope, width, depth, xarea, ilength, totlength, vel, itime, tottime]]

            for i in range(len(channel_list_3)):
                shape = channel_list_3[i][1]
                uppxl = float(channel_list_3[i][2])
                Type = channel_list_3[i][4]
                dwpxl = channel_list_3[i][5]
                avgarea = channel_list_3[i][6]
                dsarea = channel_list_3[i][7]
                upelev = channel_list_3[i][8]
                dwnelev = channel_list_3[i][9]
                ilength = channel_list_3[i][14]
                slope = channel_list_3[i][10]
                width = channel_list_3[i][11]
                depth = channel_list_3[i][12]
                xarea = channel_list_3[i][13]
                if i == 0:
                    totlength = float(ilength) + float(channel_list_2[-1][15])
                else:
                    totlength = float(ilength) + float(channel_list_3[i - 1][15])
                vel = channel_list_3[i][16]
                itime = channel_list_3[i][17]
                if i == 0:
                    tottime = float(itime) + float(channel_list_2[-1][18])
                else:
                    tottime = float(itime) + float(channel_list_3[i - 1][18])
                channel_list_3[i] = [[], shape, uppxl, [], Type, dwpxl, avgarea, dsarea, upelev, dwnelev,
                                     slope, width, depth, xarea, ilength, totlength, vel, itime, tottime]

            channel_list = channel_list_1 + channel_list_2 + channel_list_3

        if len(overland_list) != 0 and len(swale_list) != 0:
            swale_list[0][15] = overland_list[-1][15] + swale_list[0][14]
            swale_list[0][18] = overland_list[-1][18] + swale_list[0][17]
            for i in range(len(swale_list) - 1):
                swale_list[i + 1][15] = swale_list[i + 1][14] + swale_list[i][15]
                swale_list[i + 1][18] = swale_list[i + 1][17] + swale_list[i][18]

        if len(channel_list) != 0 and len(swale_list) != 0:
            channel_list[0][15] = swale_list[-1][15] + channel_list[0][14]
            channel_list[0][18] = swale_list[-1][18] + channel_list[0][17]
            for i in range(len(channel_list) - 1):
                channel_list[i + 1][15] = channel_list[i + 1][14] + channel_list[i][15]
                channel_list[i + 1][18] = channel_list[i + 1][17] + channel_list[i][18]
        elif len(channel_list) != 0 and len(overland_list) != 0:
            channel_list[0][15] = overland_list[-1][15] + channel_list[0][14]
            channel_list[0][18] = overland_list[-1][18] + channel_list[0][17]
            for i in range(len(channel_list) - 1):
                channel_list[i + 1][15] = channel_list[i + 1][14] + channel_list[i][15]
                channel_list[i + 1][18] = channel_list[i + 1][17] + channel_list[i][18]

        for i in range(len(overland_list)):
            overland_list[i][3] = "M" + str(i + 1)
        for i in range(len(swale_list)):
            swale_list[i][3] = "S" + str(i + 1)
        for i in range(len(channel_list)):
            channel_list[i][3] = "C" + str(i + 1)

        attlist = overland_list + swale_list + channel_list
        for i in range(len(attlist)):
            attlist[i][0] = i

        # obtain variables to update segment merge dialog after merge
        nover = len(overland_list)
        nswale = len(swale_list)
        nchannel = len(channel_list)
        total_tc = "{0:.2f}".format(float(attlist[-1][-1]))
        ov_tc = 0
        sw_tc = 0
        ch_tc = 0

        if len(overland_list) != 0:
            ov_tc = "{0:.2f}".format(np.sum(map(float, [col[17] for col in overland_list])))
        if len(swale_list) != 0:
            sw_tc = "{0:.2f}".format(np.sum(map(float, [col[17] for col in swale_list])))
        if len(channel_list) != 0:
            ch_tc = "{0:.2f}".format(np.sum(map(float, [col[17] for col in channel_list])))

        # populate segment merge dialog box with sub-area, Tc and number of segments for each of overland, swale, and channel
        self.sanumber.SetValue(str(arcid_global + 1))
        self.over_seg.SetValue(str(nover))
        self.swale_seg.SetValue(str(nswale))
        self.chan_seg.SetValue(str(nchannel))
        self.Overall_Tc.SetValue(str(total_tc))
        self.over_Tc.SetValue(str(ov_tc))
        self.swale_Tc.SetValue(str(sw_tc))
        self.chan_Tc.SetValue(str(ch_tc))

        # create tc list values using subshed Tc field
        tc_list = []
        subshed_poly = optfolder + "/subshed.shp"
        shedtab_tc = arcpy.SearchCursor(subshed_poly, "", "", "Tc", "")
        for t in shedtab_tc:
            get_tc = t.getValue("Tc")
            tc_list.append(get_tc)

        # update list value at index (arcid_global)
        tc_list[arcid_global] = float(attlist[-1][-1])

        # update subshed based on new Tc values
        # add "tc_list" to "subshed.shp"
        tc_lst = 0
        subshed_poly = optfolder + "/subshed.shp"
        tc_update = arcpy.UpdateCursor(subshed_poly)
        for tc in tc_update:
            sheds_tc = tc_list[tc_lst]
            tc.Tc = sheds_tc
            tc_update.updateRow(tc)
            tc_lst = tc_lst + 1

        MergedSegTable()

        # clear up and down pixels
        self.uspixel.Clear()
        self.dspixel.Clear()

    def OnApply(self, event):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        self.uspixel.Clear()
        self.dspixel.Clear()

        global subValue
        subValue = self.subID.GetValue()
        sub_num = int(subValue) - 1
        global arcid_global
        arcid_global = sub_num

        velmethexe = []
        with open(optfolder + "/vel_meth/velmethtable" + str(arcid_global) + ".csv", "rb") as inputfile:
            reader = csv.reader(inputfile)
            next(reader, None)
            for row in csv.reader(inputfile):
                velmethexe.append([x.strip(" ") for x in row])

        # create empty lists to store variable values for attribute table display

        global upval, thename_lst, thetype, downval, avgarea, dsarea, upelev, downelev, theslope, thechanwidth, thechandepth, thechanarea, inclength, downlength, thevel, thetinc, thettot

        upval = []
        thename_lst = []
        thetype = []
        downval = []
        avgarea = []
        dsarea = []
        upelev = []
        downelev = []
        theslope = []
        thechanwidth = []
        thechandepth = []
        thechanarea = []
        inclength = []
        downlength = []
        thevel = []
        thetinc = []
        thettot = []

        nopxo = 0
        nopxs = 0
        nopxc = 0
        tc_over = 0
        tc_swale = 0
        tc_channel = 0
        for i in range(len(velmethexe)):
            upval.append(int(velmethexe[i][0]) - 1)
            downval.append(int(velmethexe[i][0]))
            thetype.append(str(velmethexe[i][1]))
            if str(velmethexe[i][1]) == "overland":
                thename_lst.append("M" + str(nopxo + 1))
                nopxo = nopxo + 1
                tc_over = float(velmethexe[i][13]) + tc_over
            elif str(velmethexe[i][1]) == "swale":
                thename_lst.append("S" + str(nopxs + 1))
                nopxs = nopxs + 1
                tc_swale = float(velmethexe[i][13]) + tc_swale
            else:
                thename_lst.append("C" + str(nopxc + 1))
                nopxc = nopxc + 1
                tc_channel = float(velmethexe[i][13]) + tc_channel
            avgarea.append(float(velmethexe[i][6]))
            dsarea.append(float(velmethexe[i][3]))
            upelev.append(float(velmethexe[i][4]))
            try:
                downelev.append(float(velmethexe[i + 1][4]))
            except:
                downelev.append(float(velmethexe[i][4]) - 0.1)
            theslope.append(float(velmethexe[i][5]))
            thechanwidth.append(float(velmethexe[i][7]))
            thechandepth.append(float(velmethexe[i][8]))
            thechanarea.append(float(velmethexe[i][9]))
            inclength.append(float(velmethexe[i][10]))
            downlength.append(float(velmethexe[i][11]))
            thevel.append(float(velmethexe[i][12]))
            thetinc.append(float(velmethexe[i][13]))
            thettot.append(float(velmethexe[i][14]))

        # define global variables for number of segments of overland, swale, and channel
        global over_seg, swale_seg, channel_seg
        over_seg = nopxo
        swale_seg = nopxs
        channel_seg = nopxc

        # populate segment merge dialog box with sub-area, Tc and number of segments for each of overland, swale, and channel
        self.sanumber.SetValue(str(arcid_global + 1))
        self.over_seg.SetValue(str(over_seg))
        self.swale_seg.SetValue(str(swale_seg))
        self.chan_seg.SetValue(str(channel_seg))

        # conditions to gray out quick merge options in case of no segments
        if over_seg == 0:
            self.cb1.Disable()
        else:
            self.cb1.Enable()
        if swale_seg == 0:
            self.cb2.Disable()
        else:
            self.cb2.Enable()
        if channel_seg == 0:
            self.cb3.Disable()
        else:
            self.cb3.Enable()

        # global variables for various travel times
        global overall_tc, overland_tc, swale_tc, channel_tc
        overall_tc = "{0:.2f}".format(thettot[-1])
        overland_tc = "{0:.2f}".format(tc_over)
        swale_tc = "{0:.2f}".format(tc_swale)
        channel_tc = "{0:.2f}".format(tc_channel)

        # populate segment merge dialog box with travel time values
        self.Overall_Tc.SetValue(str(overall_tc))
        self.over_Tc.SetValue(str(overland_tc))
        self.swale_Tc.SetValue(str(swale_tc))
        self.chan_Tc.SetValue(str(channel_tc))

        # class to generate new frame for attribute table
        AttributeTable()

        # create tc list values using subshed Tc field
        tc_list = []
        subshed_poly = optfolder + "/subshed.shp"
        shedtab_tc = arcpy.SearchCursor(subshed_poly, "", "", "Tc", "")
        for t in shedtab_tc:
            get_tc = t.getValue("Tc")
            tc_list.append(get_tc)

        # update list value at index (arcid_global)
        tc_list[arcid_global] = float(attlist[-1][-1])

        # update subshed based on new Tc values
        # add "tc_list" to "subshed.shp"
        tc_lst = 0
        subshed_poly = optfolder + "/subshed.shp"
        tc_update = arcpy.UpdateCursor(subshed_poly)
        for tc in tc_update:
            sheds_tc = tc_list[tc_lst]
            tc.Tc = sheds_tc
            tc_update.updateRow(tc)
            tc_lst = tc_lst + 1

        # keep segments merge window on top
        self.Show(False)
        self.Show(True)

    def OnClose(self, event):
        self.Show(False)


class SetTcParameters(object):
    """Implementation for GISHydroNXT_addin.button9 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        # *******************************************************************************************************
        # Hydraulic coefficients for channel geometry -- default to outlet"s physiographic province
        # *******************************************************************************************************
        a = [13.87, 14.78, 10.3]
        b = [0.44, 0.39, 0.38]
        c = [0.95, 1.18, 1.01]
        d = [0.31, 0.34, 0.32]
        e = [13.17, 17.42, 10.34]
        f = [0.75, 0.73, 0.70]
        theVTab = optfolder + "/theVTab.dbf"
        theVTab = arcpy.SearchCursor("theVTab", "", "", "Province", "")
        for row in theVTab:
            if row.getValue("Province") == "A":
                Coef_W = a[0]
                Exp_W = b[0]
                Coef_D = c[0]
                Exp_D = d[0]
                Coef_A = e[0]
                Exp_A = f[0]
            elif row.getValue("Province") == "B":
                Coef_W = a[0]
                Exp_W = b[0]
                Coef_D = c[0]
                Exp_D = d[0]
                Coef_A = e[0]
                Exp_A = f[0]
            elif row.getValue("Province") == "P":
                Coef_W = a[1]
                Exp_W = b[1]
                Coef_D = c[1]
                Exp_D = d[1]
                Coef_A = e[1]
                Exp_A = f[1]
            elif row.getValue("Province") == "W":
                Coef_W = a[2]
                Exp_W = b[2]
                Coef_D = c[2]
                Exp_D = d[2]
                Coef_A = e[2]
                Exp_A = f[2]
            elif row.getValue("Province") == "E":
                Coef_W = a[2]
                Exp_W = b[2]
                Coef_D = c[2]
                Exp_D = d[2]
                Coef_A = e[2]
                Exp_A = f[2]

        del row

        # Hydraulic Coefficients -- declared globally to default in TcFrame
        global a_w
        a_w = Coef_W
        global b_w
        b_w = Exp_W
        global c_d
        c_d = Coef_D
        global d_d
        d_d = Exp_D
        global e_a
        e_a = Coef_A
        global f_a
        f_a = Exp_A


        # Extract precipitation values from tables (p2-24m)
        pr_path = Directory + "/data/prec/"
        thefilename = pr_path+ "p2-24m"
        basingrid = arcpy.sa.Times(optfolder + "/basingrid",thefilename)
        precavg = arcpy.GetRasterProperties_management(basingrid, "MEAN")
        precavg = float(precavg.getOutput(0))
        global precip_str
        precip_str = str(round(precavg/1000,2))

        # TcFrame dialog box
        TcFrame()


class SurfaceContours(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Contour Parameters", size=(400, 220))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)

        wx.StaticBox(panel, -1, "Enter parameters:", (15, 15), size=(350, 100))
        wx.StaticText(panel, -1, "Contour interval:", (25, 40))
        self.ContourInterval = wx.TextCtrl(panel, -1, value="10", pos=(130, 40), size=(130, 20))

        wx.StaticText(panel, -1, "Base contour:", (25, 70))
        self.BaseContour = wx.TextCtrl(panel, -1, value="0", pos=(130, 70), size=(130, 20))

        self.btnOK = wx.Button(panel, label="OK", pos=(50, 130))
        self.Bind(wx.EVT_BUTTON, self.OnOK, id=self.btnOK.GetId())

        self.btnCancel = wx.Button(panel, label="Cancel", pos=(225, 130))
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=self.btnCancel.GetId())

        self.Show(True)
        self.Centre(True)
        style = self.GetWindowStyle()
        self.SetWindowStyle(style | wx.STAY_ON_TOP)

    def OnOK(self, event):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        # contour and base interval
        contourInterval = float(self.ContourInterval.GetValue())
        baseContour = float(self.BaseContour.GetValue())
        # Execute Contour
        arcpy.env.addOutputsToMap = True
        demgrid = arcpy.Raster(optfolder + "/dem")
        analysis_extent = demgrid.extent
        arcpy.env.extent = analysis_extent
        arcpy.sa.Contour(optfolder + "/dem", optfolder + "/Contours", contourInterval, baseContour)
        self.Show(False)

        # change contour symbology
        layer = arcpy.mapping.Layer("Contours")
        layer.name = "Contours " + str(contourInterval) + " ft"
        arcpy.ApplySymbologyFromLayer_management(layer, r"" + Directory + "/data/mdfiles/legends/contour.lyr")

        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()

    def OnClose(self, event):
        self.Show(False)


class TcFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Time of Concentration Calculation", size=(550, 450))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)

        wx.StaticBox(panel, -1, "Select Tc Method:", (15, 15), size=(235, 120))
        TcList = ["SCS Lag Formula", "Hydrology Panel Tc Method", "Velocity Method Tc Calculation"]
        self.Tc = wx.ComboBox(panel, -1, "Velocity Method Tc Calculation", (25, 40), (205, -1), TcList, wx.CB_DROPDOWN)

        wx.StaticBox(panel, -1, "Sheet Flow", (15, 140), size=(115, 140))
        wx.StaticText(panel, -1, "ns:", (25, 175))
        self.nsValue = wx.TextCtrl(panel, -1, value="0.1", pos=(60, 170), size=(65, 25))
        wx.StaticText(panel, -1, "P[in]:", (25, 210))
        self.PValue = wx.TextCtrl(panel, -1, value=precip_str, pos=(60, 205), size=(65, 25))
        wx.StaticText(panel, -1, "L[ft]:", (25, 245))
        self.LValue = wx.TextCtrl(panel, -1, value="100", pos=(60, 240), size=(65, 25))

        wx.StaticBox(panel, -1, "Shallow Flow", (135, 140), size=(115, 140))
        self.unpaved = wx.RadioButton(panel, -1, "Unpaved", (150, 175), style=wx.RB_GROUP)
        self.paved = wx.RadioButton(panel, -1, "Paved", (150, 225))

        wx.StaticBox(panel, -1, "Channel Flow", (260, 15), size=(255, 265))
        self.NHD = wx.RadioButton(panel, -1, "Use NHD Streams", (265, 40), style=wx.RB_GROUP)
        self.infStreams = wx.RadioButton(panel, -1, "Use Inferred Streams", (265, 65))

        wx.StaticText(panel, -1, "Source Area", (425, 65))
        self.sa = wx.TextCtrl(panel, -1, value="0.0897", pos=(425, 90), size=(65, 25))

        wx.StaticText(panel, -1, "nc", (285, 95))
        self.nc = wx.TextCtrl(panel, -1, value="0.05", pos=(305, 90), size=(65, 25))

        wx.StaticText(panel, -1, "Channel Width:", (270, 120))
        wx.StaticText(panel, -1, "Coef.", (270, 143))
        self.cwCoef = wx.TextCtrl(panel, -1, value=str(a_w), pos=(305, 138), size=(65, 25))
        wx.StaticText(panel, -1, "Exp.", (395, 140))
        self.cwExp = wx.TextCtrl(panel, -1, value=str(b_w), pos=(425, 135), size=(65, 25))

        wx.StaticText(panel, -1, "Channel Depth:", (270, 170))
        wx.StaticText(panel, -1, "Coef.", (270, 193))
        self.cdCoef = wx.TextCtrl(panel, -1, value=str(c_d), pos=(305, 188), size=(65, 25))
        wx.StaticText(panel, -1, "Exp.", (390, 193))
        self.cdExp = wx.TextCtrl(panel, -1, value=str(d_d), pos=(425, 188), size=(65, 25))

        wx.StaticText(panel, -1, "Channel Area:", (270, 220))
        wx.StaticText(panel, -1, "Coef.", (270, 243))
        self.caCoef = wx.TextCtrl(panel, -1, value=str(e_a), pos=(305, 238), size=(65, 25))
        wx.StaticText(panel, -1, "Exp.", (390, 243))
        self.caExp = wx.TextCtrl(panel, -1, value=str(f_a), pos=(425, 238), size=(65, 25))

        wx.StaticBox(panel, -1, "Apply To:", (15, 295), size=(500, 55))
        self.allsub = wx.RadioButton(panel, -1, "ALL Sub-Areas", (60, 320), style=wx.RB_GROUP)
        self.selsub = wx.RadioButton(panel, -1, "ONLY Selected Sub-Areas", (300, 320))
        self.selsub.Disable()

        self.btnApply = wx.Button(panel, label="Cancel", pos=(120, 360))
        self.Bind(wx.EVT_BUTTON, self.OnSet, id=self.btnApply.GetId())

        self.btnApply = wx.Button(panel, label="Apply", pos=(300, 360))
        self.Bind(wx.EVT_BUTTON, self.OnSet, id=self.btnApply.GetId())

        self.Show(True)
        self.Centre(True)
        style = self.GetWindowStyle()
        self.SetWindowStyle(style | wx.STAY_ON_TOP)

    def OnClose(self, event):
        self.Show(False)

    def OnSet(self, event):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        # ******************************************************************************************************
        # Get Tc method selection and parameters for velocity method calculation
        # ******************************************************************************************************
        global Tc_method
        Tc_method = str(self.Tc.GetValue())
        if (Tc_method != "SCS Lag Formula" and
                Tc_method != "Hydrology Panel Tc Method" and
                Tc_method != "Velocity Method Tc Calculation"):
            pythonaddins.MessageBox("Please select Time of Concentration Method", "Tc Method Selection Warning", 0)
            self.Show(False)
            self.Show(True)

        else:
            self.Show(False)

            # define global variable for all Tc parameters
            global Tc_ns
            Tc_ns = self.nsValue.GetValue()
            global Tc_P
            Tc_P = self.PValue.GetValue()
            global Tc_L
            Tc_L = self.LValue.GetValue()
            global Tc_paved
            Tc_paved = self.paved.GetValue()
            global Tc_unpaved
            Tc_unpaved = self.unpaved.GetValue()
            global Tc_NHD
            Tc_NHD = self.NHD.GetValue()
            global Tc_infStreams
            Tc_infStreams = self.infStreams.GetValue()
            global Tc_sa
            Tc_sa = self.sa.GetValue()
            global Tc_nc
            Tc_nc = self.nc.GetValue()
            global Tc_cwCoef
            Tc_cwCoef = self.cwCoef.GetValue()
            global Tc_cwExp
            Tc_cwExp = self.cwExp.GetValue()
            global Tc_cdCoef
            Tc_cdCoef = self.cdCoef.GetValue()
            global Tc_cdExp
            Tc_cdExp = self.cdExp.GetValue()
            global Tc_caCoef
            Tc_caCoef = self.caCoef.GetValue()
            global Tc_caExp
            Tc_caExp = self.caExp.GetValue()
            global Tc_allsub
            Tc_allsub = self.allsub.GetValue()
            global Tc_selsub
            Tc_selsub = self.selsub.GetValue()

        # lable subwatersheds with "ARCID"
        mxd = arcpy.mapping.MapDocument("CURRENT")
        layer = arcpy.mapping.ListLayers(mxd, "subshed")[0]
        if layer.supports("LABELCLASSES"):
            for lblclass in layer.labelClasses:
                lblclass.showClassLabels = True

        lblclass.expression = "[ARCID]"
        layer.showLabels = True
        arcpy.RefreshActiveView()

        # ******************************************************************************************************
        # turn watershed delineation OFF and basin stat ON
        # ******************************************************************************************************
        button9.enabled = False
        button10.enabled = True
        save(optfolder)

class ThomasDischarge(object):
    """Implementation for GISHydroNXT_addin.button4 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        # *******************************************************************************************************
        # gage checking and selection
        # *******************************************************************************************************
        # added on 10-13-2017: for to include in intersect analysis
        usgsgages = r"" + Directory + "/data/maryland/usgsgagesm.shp"
        mdgages = r"" + Directory + "/data/maryland/mdgagedstreams2016.shp"
        outletpoint = optfolder + "/outletpoint.shp"
        mask = optfolder + "/mask.shp"
        gagefound = hydro.CheckGages(usgsgages, mdgages, outletpoint, mask, optfolder)

        global gagelist
        gagelist = gagefound

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            # *******************************************************************************************************
            # Check if gauge exist in the "gagefound" list and implement if/else for
            # gage dialog to pop up
            # *******************************************************************************************************
            if not gagelist:
                # *****************************************************************************
                # read global variables (defined in basin stat) into variable names for use
                # in Thomas peak discharge analysis
                FC = float(hydro.FC)  # LI is already declared as global variable
                FC = "{0:.2f}".format(FC)
                DA = areami2
                HA = float(hydro.pctAsoil)
                HC = float(hydro.pctCsoil)
                HD = float(hydro.pctDsoil)
                SLL = float(landslope)
                HCD = float(HC + HD)
                ImpA = float(IA)

                # *****************************************************************************
                # Define initial values and lists. add and index those lists as part of
                # dictionary
                # *****************************************************************************
                sumarea = 0
                provstring = ""
                outtaskerstring = ""
                sQ1p25 = 0
                sQ1p50 = 0
                sQ1p75 = 0
                sQ2 = 0
                sQ5 = 0
                sQ10 = 0
                sQ25 = 0
                sQ50 = 0
                sQ100 = 0
                sQ200 = 0
                sQ500 = 0
                Q1p25list = [0, 0, 0, 0, 0, 0, 0, 0]
                Q1p50list = [0, 0, 0, 0, 0, 0, 0, 0]
                Q2list = [0, 0, 0, 0, 0, 0, 0, 0]
                Q5list = [0, 0, 0, 0, 0, 0, 0, 0]
                Q10list = [0, 0, 0, 0, 0, 0, 0, 0]
                Q25list = [0, 0, 0, 0, 0, 0, 0, 0]
                Q50list = [0, 0, 0, 0, 0, 0, 0, 0]
                Q100list = [0, 0, 0, 0, 0, 0, 0, 0]
                Q200list = [0, 0, 0, 0, 0, 0, 0, 0]
                Q500list = [0, 0, 0, 0, 0, 0, 0, 0]

                qlist = {"1": [], "2": [], "3": [], "4": [], "5": [], "6": [], "7": [], "8": [], "9": [], "10": []}
                qlist["1"].extend(Q1p25list)
                qlist["2"].extend(Q1p50list)
                qlist["3"].extend(Q2list)
                qlist["4"].extend(Q5list)
                qlist["5"].extend(Q10list)
                qlist["6"].extend(Q25list)
                qlist["7"].extend(Q50list)
                qlist["8"].extend(Q100list)
                qlist["9"].extend(Q200list)
                qlist["10"].extend(Q500list)

                # *****************************************************************************
                # copy root tasker files into optfolder (not needed after tasker was recreated in python), the folder still needs to be created
                # *****************************************************************************
                src = r"" + Directory + "/data/tasker"
                dst = optfolder + "/tasker/"
                if not os.path.exists(dst):
                    os.mkdir(dst, 0o755)
                for the_file in os.listdir(src):
                    file_path = os.path.join(src, the_file)
                    if os.path.isfile(file_path):
                        shutil.copyfile(file_path, dst + the_file)

                # *****************************************************************************
                # read zonal stat table (theVTab), declare fields into variables, and
                # loop through theVTab fields to getValue
                # *****************************************************************************
                theVTab = optfolder + "/theVTab.dbf"
                theVTab = arcpy.SearchCursor("theVTab", "", "", "Count", "")

                # *****************************************************************************
                # loop to get total count of pixels of basingrid
                # *****************************************************************************
                for each in theVTab:
                    count = each.getValue("Count")
                    sumarea = sumarea + count
                sumArea = sumarea

                del each

                theVTab = optfolder + "/theVTab.dbf"
                theVTab = arcpy.SearchCursor("theVTab", "", "",
                                             "Province;Count;Q1.25;Q1.50;Q1.75;Q2;Q5;Q10;Q25;Q50;Q100;Q200;Q500", "")

                # loop to begin tasker handling and area weighted analysis
                discharge_values = []
                flood_intervals = []
                for row in theVTab:
                    AreaField = float(row.getValue("Count"))
                    areapercent = float((AreaField / sumArea) * 100)
                    if row.getValue("Province") == "A":
                        provstring = provstring + "       -Appalachian Plateaus and Allegheny Ridges %s percent of area" "\n" % (
                            "{0:.2f}".format(areapercent))
                    elif row.getValue("Province") == "B":
                        provstring = provstring + "       -Blue Ridge and Great Valley %s percent of area" "\n" % (
                            "{0:.2f}".format(areapercent))
                    elif row.getValue("Province") == "P":
                        provstring = provstring + "       -Piedmont %s percent of area" "\n" % (
                            "{0:.2f}".format(areapercent))
                    elif row.getValue("Province") == "W":
                        provstring = provstring + "       -Western Coastal Plain %s percent of area" "\n" % (
                            "{0:.2f}".format(areapercent))
                    elif row.getValue("Province") == "E":
                        provstring = provstring + "       -Eastern Coastal Plain %s percent of area" "\n" % (
                            "{0:.2f}".format(areapercent))

                    else:
                        pythonaddins.MessageBox("No Province Selected", "Problem...", 0)

                    intaskerstring = "thomasout.txt" "\n"
                    intaskerstring = intaskerstring + "" "\n"
                    if row.getValue("Province") == "A":
                        intaskerstring = intaskerstring + "a"  "\n"
                        intaskerstring = intaskerstring + "%s" "\n" % (DA)
                        intaskerstring = intaskerstring + "%s" "\n" % (SLL)
                    elif row.getValue("Province") == "B":
                        intaskerstring = intaskerstring + "p"  "\n"
                        intaskerstring = intaskerstring + "%s" "\n" % (DA)
                        intaskerstring = intaskerstring + "%s" "\n" % (FC)
                        intaskerstring = intaskerstring + "%s" "\n" % (LI)
                        intaskerstring = intaskerstring + "%s" "\n" % (ImpA)
                    elif row.getValue("Province") == "P":
                        intaskerstring = intaskerstring + "p"  "\n"
                        intaskerstring = intaskerstring + "%s" "\n" % (DA)
                        intaskerstring = intaskerstring + "%s" "\n" % (FC)
                        intaskerstring = intaskerstring + "%s" "\n" % (LI)
                        intaskerstring = intaskerstring + "%s" "\n" % (ImpA)

                    elif row.getValue("Province") == "W":
                        intaskerstring = intaskerstring + "wc" "\n"
                        intaskerstring = intaskerstring + "%s" "\n" % (DA)
                        intaskerstring = intaskerstring + "%s" "\n" % (ImpA)
                        intaskerstring = intaskerstring + "%s" "\n" % (HCD)
                    elif row.getValue("Province") == "E":
                        intaskerstring = intaskerstring + "ec" "\n"
                        intaskerstring = intaskerstring + "%s" "\n" % (DA)
                        intaskerstring = intaskerstring + "%s" "\n" % (SLL)
                        intaskerstring = intaskerstring + "%s" "\n" % (HA)

                    intaskerstring = intaskerstring + "N" "\n"
                    intaskerstring = intaskerstring + "N" "\n"
                    thomas2020 = optfolder + "/tasker/thomas2020.txt"
                    thomasin = open(thomas2020, "w")
                    thomasin.write(intaskerstring)
                    thomasin.close()

                    tasker.RRE(optfolder,thomas2020)        # Added in Feb 2020 replacing the

                    # REMOVED AFTER TASKER.PY
                    #os.chdir(optfolder + "/tasker")
                    #os.system(optfolder + "/tasker/thomas2016auto.exe")
                    #time.sleep(4)

                    # *****************************************************************************
                    # open "thomasout.txt" file and read peak discharges
                    # *****************************************************************************
                    infilename = optfolder + "/tasker/thomasout.txt"
                    infile = open(infilename, "r").readlines()

                    for i in infile:
                        outtaskerstring = outtaskerstring + i

                    # *****************************************************************************
                    # Discharge values for each province for "theVTab" are not possible to write
                    # at this time because update cursor can"t be used within search cursor "for"
                    # loop -- appended to list "discharge_values" and will be used at the end to
                    # update discharge values in "theVTab"
                    # *****************************************************************************
                    theline = infile[11].split()
                    Q1p25 = theline[1]
                    discharge_values.append(Q1p25)
                    theline = infile[12].split()
                    Q1p50 = theline[1]
                    discharge_values.append(Q1p50)
                    Q1p75 = -999
                    theline = infile[13].split()
                    Q2 = theline[1]
                    discharge_values.append(Q2)
                    theline = infile[14].split()
                    Q5 = theline[1]
                    discharge_values.append(Q5)
                    theline = infile[15].split()
                    Q10 = theline[1]
                    discharge_values.append(Q10)
                    theline = infile[16].split()
                    Q25 = theline[1]
                    discharge_values.append(Q25)
                    theline = infile[17].split()
                    Q50 = theline[1]
                    discharge_values.append(Q50)
                    theline = infile[18].split()
                    Q100 = theline[1]
                    discharge_values.append(Q100)
                    theline = infile[19].split()
                    Q200 = theline[1]
                    discharge_values.append(Q200)
                    theline = infile[20].split()
                    Q500 = theline[1]
                    discharge_values.append(Q500)

                    # *****************************************************************************
                    # compute discharge and assign confidence intervals to qlist entries
                    # *****************************************************************************
                    Q1p25 = float(Q1p25)
                    sQ1p25 = sQ1p25 + (Q1p25 * AreaField)
                    Q1p50 = float(Q1p50)
                    sQ1p50 = sQ1p50 + (Q1p50 * AreaField)
                    Q1p75 = float(Q1p75)
                    sQ1p75 = sQ1p75 + (Q1p75 * AreaField)
                    Q2 = float(Q2)
                    sQ2 = sQ2 + (Q2 * AreaField)
                    Q5 = float(Q5)
                    sQ5 = sQ5 + (Q5 * AreaField)
                    Q10 = float(Q10)
                    sQ10 = sQ10 + (Q10 * AreaField)
                    Q25 = float(Q25)
                    sQ25 = sQ25 + (Q25 * AreaField)
                    Q50 = float(Q50)
                    sQ50 = sQ50 + (Q50 * AreaField)
                    Q100 = float(Q100)
                    sQ100 = sQ100 + (Q100 * AreaField)
                    Q200 = float(Q200)
                    sQ200 = sQ200 + (Q200 * AreaField)
                    Q500 = float(Q500)
                    sQ500 = sQ500 + (Q500 * AreaField)

                    for i, lines in enumerate(infile):
                        if i > 24 and i < 35:
                            line = lines.split()
                            a1 = line[1]
                            a2 = line[2]
                            a3 = line[3]
                            a4 = line[4]
                            a5 = line[5]
                            a6 = line[6]
                            a7 = line[7]
                            a8 = line[8]
                            c1 = [(float(a1) * areapercent) / 100]
                            c2 = [(float(a2) * areapercent) / 100]
                            c3 = [(float(a3) * areapercent) / 100]
                            c4 = [(float(a4) * areapercent) / 100]
                            c5 = [(float(a5) * areapercent) / 100]
                            c6 = [(float(a6) * areapercent) / 100]
                            c7 = [(float(a7) * areapercent) / 100]
                            c8 = [(float(a8) * areapercent) / 100]
                            qlist = [c1, c2, c3, c4, c5, c6, c7, c8]
                            flood_intervals.append(qlist)

                lists = {i: [el[0] for el in v] for i, v in enumerate(flood_intervals, start=1)}

                del row

                # *****************************************************************************
                # index and count number of sub-lists -- prepare for TaskerString function
                # *****************************************************************************
                number = hydro.FloodIntervalLists(lists)
                if number > 10:
                    q_list1 = [sum(i) for i in zip(lists[1], lists[11])]
                    q_list2 = [sum(i) for i in zip(lists[2], lists[12])]
                    q_list3 = [sum(i) for i in zip(lists[3], lists[13])]
                    q_list4 = [sum(i) for i in zip(lists[4], lists[14])]
                    q_list5 = [sum(i) for i in zip(lists[5], lists[15])]
                    q_list6 = [sum(i) for i in zip(lists[6], lists[16])]
                    q_list7 = [sum(i) for i in zip(lists[7], lists[17])]
                    q_list8 = [sum(i) for i in zip(lists[8], lists[18])]
                    q_list9 = [sum(i) for i in zip(lists[9], lists[19])]
                    q_list10 = [sum(i) for i in zip(lists[10], lists[20])]
                else:
                    q_list1 = lists[1]
                    q_list2 = lists[2]
                    q_list3 = lists[3]
                    q_list4 = lists[4]
                    q_list5 = lists[5]
                    q_list6 = lists[6]
                    q_list7 = lists[7]
                    q_list8 = lists[8]
                    q_list9 = lists[9]
                    q_list10 = lists[10]

                # *****************************************************************************
                # discharge computation based on province
                # *****************************************************************************
                Q1p25 = int(sQ1p25 / sumArea)
                Q1p50 = int(sQ1p50 / sumArea)
                Q1p75 = int(sQ1p75 / sumArea)
                Q2 = int(sQ2 / sumArea)
                Q5 = int(sQ5 / sumArea)
                Q10 = int(sQ10 / sumArea)
                Q25 = int(sQ25 / sumArea)
                Q50 = int(sQ50 / sumArea)
                Q100 = int(sQ100 / sumArea)
                Q200 = int(sQ200 / sumArea)
                Q500 = int(sQ500 / sumArea)

                # *****************************************************************************
                # Basin statistics analysis date
                # *****************************************************************************
                now = datetime.now()
                month = now.strftime("%B")
                day = now.strftime("%d")
                year = now.strftime("%Y")

                # *****************************************************************************
                # Text file string variables
                # *****************************************************************************
                datastring = ""
                datastring = datastring + "GISHydro Release Version Date: %s" "\n" % (Modifieddt)
                datastring = datastring + "Project Name:                  %s" % (proj)
                datastring = datastring + "" "\n"
                datastring = datastring + "Analysis Date:                 %s %s, %s " "\n" % (month, day, year)
                datastring = datastring + "Thomas Version:                %s " "\n \n" % (thomasversion)
                datastring = datastring + "Geographic Province(s):" "\n"
                datastring = datastring + provstring
                datastring = datastring + "" "\n"
                datastring = datastring + "Q(1.25):   %s cfs" "\n" % (Q1p25)
                datastring = datastring + "Q(1.50):   %s cfs" "\n" % (Q1p50)
                datastring = datastring + "Q(2):      %s cfs" "\n" % (Q2)
                datastring = datastring + "Q(5):      %s cfs" "\n" % (Q5)
                datastring = datastring + "Q(10):     %s cfs" "\n" % (Q10)
                datastring = datastring + "Q(25):     %s cfs" "\n" % (Q25)
                datastring = datastring + "Q(50):     %s cfs" "\n" % (Q50)
                datastring = datastring + "Q(100):    %s cfs" "\n" % (Q100)
                datastring = datastring + "Q(200):    %s cfs" "\n" % (Q200)
                datastring = datastring + "Q(500):    %s cfs" "\n" % (Q500)

                datastring = datastring + "" "\n"
                datastring = datastring + "Area Weighted Prediction Intervals (from Tasker)" "\n"
                datastring = datastring + " Return     50 PERCENT        67 PERCENT        90 PERCENT        95 PERCENT" "\n"
                datastring = datastring + " Period  lower    upper    lower    upper    lower    upper    lower    upper" "\n"
                datastring = datastring + hydro.TaskerString(1.25, q_list1)
                datastring = datastring + "" "\n"
                datastring = datastring + hydro.TaskerString(1.50, q_list2)
                datastring = datastring + "" "\n"
                datastring = datastring + hydro.TaskerString(2, q_list3)
                datastring = datastring + "" "\n"
                datastring = datastring + hydro.TaskerString(5, q_list4)
                datastring = datastring + "" "\n"
                datastring = datastring + hydro.TaskerString(10, q_list5)
                datastring = datastring + "" "\n"
                datastring = datastring + hydro.TaskerString(25, q_list6)
                datastring = datastring + "" "\n"
                datastring = datastring + hydro.TaskerString(50, q_list7)
                datastring = datastring + "" "\n"
                datastring = datastring + hydro.TaskerString(100, q_list8)
                datastring = datastring + "" "\n"
                datastring = datastring + hydro.TaskerString(200, q_list9)
                datastring = datastring + "" "\n"
                datastring = datastring + hydro.TaskerString(500, q_list10)
                datastring = datastring + "" "\n"

                datastring = datastring + "" "\n"
                datastring = datastring + "" "\n"
                datastring = datastring + "Individual Province Tasker Analyses Follow: " "\n"
                datastring = datastring + ""
                datastring = datastring + outtaskerstring

                # *****************************************************************************
                # write strings to basin stat text file.
                # *****************************************************************************
                defFN = optfolder + "/frdischarges.txt"
                fr = open(defFN, "w")
                fr.write(datastring)
                fr.close()

                # *****************************************************************************
                # open "frdischarges" file in text editor
                # *****************************************************************************
                hydro.openbrowser(defFN)

                # *****************************************************************************
                # turn layers ON/OFF in current data frame
                # *****************************************************************************
                mxd = arcpy.mapping.MapDocument("CURRENT")
                df = arcpy.mapping.ListDataFrames(mxd)[0]
                layers = arcpy.mapping.ListLayers(mxd, "", df)
                for lyr in layers:
                    if lyr.name == "mdlayer":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "gagefield":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletpoint":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletpoly":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "limepoly":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletpoint":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletbuffer":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "gagefield":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletpoly":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "mask_ints":  # added on 10-22-2017: this layer is output of intersect tool
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "gauge_outlet":  # added on 10-22-2017: this layer is output of join and intersect tool
                        arcpy.mapping.RemoveLayer(df, lyr)

                arcpy.RefreshTOC()
                arcpy.RefreshActiveView()

                # *****************************************************************************
                # turn peak discharge OFF, and "S" & Add Streams ON
                # *****************************************************************************
                button4.enabled = False
                tool7.enabled = False
                button5_1.enabled = True
                button6.enabled = True
                tool4.enabled = True
            else:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                df = arcpy.mapping.ListDataFrames(mxd)[0]
                layers = arcpy.mapping.ListLayers(mxd, "", df)
                for lyr in layers:
                    if lyr.name == "mdlayer":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "gagefield":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletpoint":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletpoly":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "limepoly":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletpoint":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletbuffer":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "gagefield":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "outletpoly":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "mask_ints":  # added on 10-22-2017: this layer is output of intersect tool
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "gauge_outlet":  # added on 10-22-2017: this layer is output of join and intersect tool
                        arcpy.mapping.RemoveLayer(df, lyr)
                GageListThomas()


class TransectLine(object):
    """Implementation for GISHydroNXT_addin.tool6 (Tool)"""

    def __init__(self):
        self.enabled = False
        self.cursor = 1
        self.shape = "Line"

    def onLine(self, line_geometry):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        arcpy.env.extent = "MAXOF"
        # *******************************************************************************************************
        # Add DEM elevation and distance fields to transect line
        # *******************************************************************************************************
        arcpy.CopyFeatures_management(line_geometry, optfolder + "/trans_line.shp")  # saving transect line

        ### TRAP ERROR IN CASE TRANSECT IS IN INCORRECT REACH
        arcidnode = []
        fromnode = []
        tonode = []
        subriver_prop = arcpy.SearchCursor(optfolder + "/subrivers.shp", "", "", "ARCID;From_Node;To_Node", "")
        for node in subriver_prop:
            arcidnode.append(int(node.getValue("ARCID")))
            fromnode.append(int(node.getValue("From_Node")))
            tonode.append(int(node.getValue("To_Node")))
        subreach_lst = [x for x in fromnode if x in tonode]
        arcid_list = []
        for sub in subreach_lst:
            index = fromnode.index(sub)
            arcid_list.append(arcidnode[index])
        line_trun = [optfolder + "/trans_line.shp", optfolder + "/subshed.shp"]
        arcpy.Intersect_analysis(line_trun, optfolder + "/line_trun.shp", "ALL", "#", "INPUT")
        transect_prop = arcpy.SearchCursor(optfolder + "/line_trun.shp", "", "", "ARCID", "")
        transnodes = []
        for tra in transect_prop:
            transnodes.append(int(tra.getValue("ARCID")))
        condition_in = [x for x in transnodes if x in arcid_list]
        if len(condition_in) == 0:
            pythonaddins.MessageBox("Transect must be drawn into a routing reach", "Transect Digitize Error")
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "trans_line":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "line_trun":
                    arcpy.mapping.RemoveLayer(df, lyr)
            arcpy.Delete_management(optfolder + "/trans_line.shp", "")
            arcpy.Delete_management(optfolder + "/line_trun.shp", "")
            return

        global provstring
        arbitrary_length = 500
        if "Coastal" in provstring:
            arbitrary_length = 800

        with arcpy.da.SearchCursor(optfolder + "/trans_line.shp", "SHAPE@") as rows:
            for row in rows:
                if row[0].length > arbitrary_length:
                    pythonaddins.MessageBox("Transect length is " + str(round((row[0].length), 1)) + " meters. Maximum allowable length is " +
                                            str(arbitrary_length) + "m. Please draw it again.",
                                            "Transect Length Warning", 0)

                    mxd = arcpy.mapping.MapDocument("CURRENT")
                    df = arcpy.mapping.ListDataFrames(mxd)[0]
                    layers = arcpy.mapping.ListLayers(mxd, "", df)
                    for lyr in layers:
                        if lyr.name == "trans_line":
                            arcpy.mapping.RemoveLayer(df, lyr)
                    arcpy.Delete_management(optfolder + "/trans_line.shp", "")
                    return

                else:
                    arcpy.CreateRandomPoints_management(optfolder, "transect", optfolder + "/line_trun.shp", "#", "#",
                                                        6, "#",
                                                        "#")  # 02-15-2015: changed sampling interval from 3 to 6
                    if not len(arcpy.ListFields(optfolder + "/transect.shp", "Distance")) > 0:
                        arcpy.AddField_management(optfolder + "/transect.shp", "Distance", "FLOAT", 15, 4)

                    # add distance along transect
                    rows = arcpy.UpdateCursor(optfolder + "/transect.shp")
                    distance = 0
                    count = 0
                    for row in rows:
                        count = count + 1
                        distance = distance + 9.8424
                        row.Distance = distance
                        rows.updateRow(row)

                    # delete field "CID"
                    arcpy.DeleteField_management(optfolder + "/transect.shp", "CID")

                    arcpy.CopyFeatures_management(optfolder + "/transect.shp", optfolder + "/aux_folder/transect_aux.shp")

                    # extract multivalues to points using DEM
                    arcpy.sa.ExtractMultiValuesToPoints(optfolder + "/aux_folder/transect_aux.shp", optfolder + "/dem", "NONE")

                    arcpy.sa.ExtractMultiValuesToPoints(optfolder + "/aux_folder/transect_aux.shp", optfolder + "/flowacc", "NONE")

                    arcpy.CopyFeatures_management(optfolder + "/aux_folder/transect_aux.shp", optfolder + "/transect.shp")

                    # correct start to end distance and elevation values in feet
                    rows = arcpy.UpdateCursor(optfolder + "/transect.shp")
                    for row in rows:
                        row.setValue("Distance", row.getValue("Distance") - 9.8424)
                        row.setValue("dem", row.getValue(
                            "dem"))  # elevation values are already in feet -- no need for conversion
                        rows.updateRow(row)

                    dem_list = []
                    distance_list = []
                    flowacc_list = []
                    upstreamDA = optfolder + "/transect.shp"
                    usda = arcpy.SearchCursor(upstreamDA, "", "", "dem;Distance;flowacc", "")
                    for da in usda:
                        dem = da.getValue("dem")
                        dem_list.append(dem)
                        distance = da.getValue("Distance")
                        distance_list.append(distance)
                        flow = da.getValue("flowacc")
                        flowacc_list.append(flow)

                    global trans_distance
                    trans_distance = [round(elem, 1) for elem in distance_list]

                    global trans_dem
                    trans_dem = [round(elem, 1) for elem in dem_list]
                    # *******************************************************************************************************
                    # Transection line width, min, and max elevation
                    # TWL = Transection Width Length
                    # TWE_min = minimum Transection Width Elevation
                    # TWE_max = maximum Transection Width Elevation
                    # *******************************************************************************************************
                    global TWL
                    TWL_max = float(max(distance_list))
                    TWL = "{0:.2f}".format(TWL_max)

                    global TWE_min
                    TWE1 = float(min(dem_list))
                    TWE_min = "{0:.2f}".format(TWE1)

                    global TWE_max
                    TWE2 = float(max(dem_list))
                    TWE_max = "{0:.2f}".format(TWE2)

                    # *******************************************************************************************************
                    # Calculate upstream drainage area
                    # *******************************************************************************************************
                    upda_max = max(flowacc_list)
                    flowacc_cell = optfolder + "/flowacc"
                    rast = arcpy.Raster(flowacc_cell)
                    cellsize = rast.meanCellWidth
                    cellsq = cellsize * cellsize

                    global areami2_usda
                    areami2_usda = float((upda_max * cellsq) / 2588881)  # conversion into sq miles
                    areami2_usda = "{0:.2f}".format(areami2_usda)

                    # *******************************************************************************************************
                    # get reach slope using subsheds, transect, and subriver shapefiles
                    # *******************************************************************************************************
                    arcpy.env.qualifiedFieldNames = False  # to have attribute tables with original names
                    rSlope_intersect = [optfolder + "/line_trun.shp", optfolder + "/subshed.shp"]
                    arcpy.Intersect_analysis(rSlope_intersect, optfolder + "/reachslope.shp", "ALL", "#", "INPUT")
                    arcpy.DeleteField_management(optfolder + "/reachslope.shp",
                                                 "FID_trans_;Id;AreaMi2;TcMethod;CurveNum;LngFlwPth;Tc;Slope")
                    arcpy.MakeFeatureLayer_management(optfolder + "/reachslope.shp", "reachlayer")
                    FID = "ARCID"
                    Code = "ARCID"
                    arcpy.AddJoin_management("reachlayer", FID, optfolder + "/subrivers.shp", Code)
                    arcpy.FeatureClassToShapefile_conversion("reachlayer", optfolder)
                    reach = arcpy.SearchCursor(optfolder + "/reachlayer.shp", "", "", "Slope", "")
                    for row in reach:
                        reach_slope = row.getValue("Slope")
                    del row
                    del reach

                    global reachslope
                    reachslope = "{0:.5f}".format(reach_slope)

                    # *******************************************************************************************************
                    # get bankfull channel width and depth
                    # *******************************************************************************************************
                    sumarea = 0
                    theVTab = optfolder + "/theVTab.dbf"
                    theVTab = arcpy.SearchCursor("theVTab", "", "", "Count", "")
                    for each in theVTab:
                        count = each.getValue("Count")
                        sumarea = sumarea + count
                    sumArea = sumarea

                    del each

                    AParea = 0
                    PDarea = 0
                    CParea = 0
                    theVTab = optfolder + "/theVTab.dbf"
                    theVTab = arcpy.SearchCursor("theVTab", "", "", "Province;Count", "")
                    for row in theVTab:
                        theArea = float(row.getValue("Count"))
                        areapercent = float((theArea / sumArea) * 100)
                        if row.getValue("Province") == "A":
                            AParea = AParea + float(areapercent)
                        elif row.getValue("Province") == "B":
                            AParea = AParea + float(areapercent)
                        elif row.getValue("Province") == "P":
                            PDarea = PDarea + float(areapercent)
                        elif row.getValue("Province") == "W":
                            CParea = CParea + float(areapercent)
                        elif row.getValue("Province") == "E":
                            CParea = CParea + float(areapercent)

                    del row

                    regionarea = []
                    AParea = AParea / 100
                    PDarea = PDarea / 100
                    CParea = CParea / 100
                    regionarea.append(AParea)
                    regionarea.append(PDarea)
                    regionarea.append(CParea)

                    a = [13.87, 14.78, 10.3]
                    b = [0.44, 0.39, 0.38]
                    c = [0.95, 1.18, 1.01]
                    d = [0.31, 0.34, 0.32]

                    Coef_W = sum(map(mul, regionarea, a))
                    Exp_W = sum(map(mul, regionarea, b))

                    global Wbf
                    WidthBF = Coef_W * (float(areami2_usda) ** Exp_W)
                    Wbf = "{0:.2f}".format(WidthBF)

                    Coef_D = sum(map(mul, regionarea, c))
                    Exp_D = sum(map(mul, regionarea, d))

                    global Dbf
                    DepthBF = Coef_D * (float(areami2_usda) ** Exp_D)
                    Dbf = "{0:.2f}".format(DepthBF)

                    # *******************************************************************************************************
                    # turn layers ON/OFF in current data frame
                    # *******************************************************************************************************
                    mxd = arcpy.mapping.MapDocument("CURRENT")
                    df = arcpy.mapping.ListDataFrames(mxd)[0]
                    layers = arcpy.mapping.ListLayers(mxd, "", df)
                    for lyr in layers:
                        if lyr.name == "transect":
                            arcpy.mapping.RemoveLayer(df, lyr)
                        if lyr.name == "reachlayer":
                            arcpy.mapping.RemoveLayer(df, lyr)
                        if lyr.name == "reachslope":
                            arcpy.mapping.RemoveLayer(df, lyr)
                        if lyr.name == "elevzones":
                            arcpy.mapping.RemoveLayer(df, lyr)
                        if lyr.name == "elevmerge":
                            arcpy.mapping.RemoveLayer(df, lyr)
                        if lyr.name == "line_trun":
                            arcpy.mapping.RemoveLayer(df, lyr)
                        if os.path.exists(optfolder + "/reachlayer.shp"):
                            arcpy.Delete_management(optfolder + "/reachlayer.shp", "")
                        if os.path.exists(optfolder + "/aux_folder/transect_aux.shp"):
                            arcpy.Delete_management(optfolder + "/aux_folder/transect_aux.shp", "")

                    # *******************************************************************************************************
                    # call Cross-Section Editor dialog
                    # *******************************************************************************************************
                    Xeditor()
                    save(optfolder)


class TR20ControlPanel(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "GISHydroNXT - WinTR-20 Control Panel", size=(390, 540))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)
        self.data = [(1, 2), (2, 3), (3, 5)]

        wx.StaticBox(panel, -1, "WinTR-20 Input/Output File Locations", (15, 15), size=(340, 55))
        wx.StaticText(panel, -1, "Input File:", (40, 35))
        InputFile = optfolder + "WinTR20/TR20in.txt"
        self.InputFile = wx.TextCtrl(panel, -1, value=InputFile, pos=(110, 35), size=(220, 18))

        wx.StaticBox(panel, -1, "Input Options", (15, 80), size=(340, 55))
        wx.StaticText(panel, -1, "Hydrograph:", (40, 105))
        HydrographList = ["Standard PRF 484", "DelMarVa PRF 284"]
        self.DelMarVa = wx.ComboBox(panel, -1, "Standard PRF 484", (200, 105), (130, -1), HydrographList,
                                    wx.CB_DROPDOWN)

        wx.StaticBox(panel, -1, "Standard Control Output Options", (15, 145), size=(340, 80))
        wx.StaticText(panel, -1, "Hydrograph:", (40, 185))
        self.summary = wx.RadioButton(panel, -1, "Summary Output", (220, 170), style=wx.RB_GROUP)
        self.detail = wx.RadioButton(panel, -1, "Detailed Output", (220, 200))

        wx.StaticBox(panel, -1, "Executive Control Options", (15, 235), size=(340, 70))
        wx.StaticText(panel, -1, "Main time Increment:", (40, 265))
        self.Main0 = wx.TextCtrl(panel, -1, value="0.1", pos=(250, 265), size=(40, 20))
        wx.StaticText(panel, -1, "hrs.", (295, 265))

        wx.StaticBox(panel, -1, "Rainfall", (15, 315), size=(340, 100))
        wx.StaticText(panel, -1, "ARC:", (210, 355))
        ARCList = ["1", "2", "3"]
        self.ARC = wx.ComboBox(panel, -1, "2", (250, 350), (40, -1), ARCList, wx.CB_DROPDOWN)
        wx.StaticText(panel, -1, "Perform Areal Reduction", (130, 380))
        self.ARF = wx.CheckBox(panel, -1, "", pos=(110, 380))
        self.ARF.SetValue(True)

        self.btnEdit = wx.Button(panel, label="Choose Storm Depth(s)", pos=(33, 350))
        self.Bind(wx.EVT_BUTTON, self.OnEdit, id=self.btnEdit.GetId())

        self.btnOK = wx.Button(panel, label="OK", pos=(50, 440))
        self.Bind(wx.EVT_BUTTON, self.OnOK, id=self.btnOK.GetId())

        self.btnCancel = wx.Button(panel, label="Cancel", pos=(225, 440))
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=self.btnCancel.GetId())

        self.Show(True)
        self.Centre(True)
        style = self.GetWindowStyle()
        self.SetWindowStyle(style | wx.STAY_ON_TOP)

        global year_uSpecified
        year_uSpecified = cb_selected
        global prec_uSpecified
        prec_uSpecified = critavg

    def OnClose(self, event):
        self.Show(False)

    def OnEdit(self, event):
        ChooseStormDepths()

    def OnOK(self, event):
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        self.Show(False)
        Main_increment = self.Main0.GetValue()

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            # ******************************************************************************************************
            # TR20 input text files with and without subwatersheds will be generated
            # ******************************************************************************************************
            # create "FROM NODE", "TO NODE", and "GRID CODE" lists
            fn_lst = []
            ai_lst = []
            ln_lst = []
            tn_lst = []
            subriver = optfolder + "/subrivers.shp"
            subriver_len = arcpy.SearchCursor(subriver, "", "", "ARCID;FROM_NODE;TO_NODE;Length", "")
            for node in subriver_len:
                ai_lst.append(int(node.getValue("ARCID")))
                fn_lst.append(int(node.getValue("FROM_NODE")))
                tn_lst.append(int(node.getValue("TO_NODE")))
                ln_lst.append(int(node.getValue("Length")))
            subreach_lst = [x for x in fn_lst if x in tn_lst]
            reach_lst = []
            len_lst = []
            for sub in subreach_lst:
                index = fn_lst.index(sub)
                reach_lst.append(ai_lst[index])
                len_lst.append(ln_lst[index])

            reach_list = []
            for tn in tn_lst:
                aux = False
                for ind, fn in enumerate(fn_lst):
                    if tn == fn:
                        reach_list.append(ai_lst[ind])
                        aux = True
                if aux == False:
                    reach_list.append(0)

            lst_reach = []
            for rch in reach_list:
                if rch == 0:
                    lst_reach.append("Outlet")
                else:
                    lst_reach.append("Reach" + str(rch))

            # create "slope", "area (mi^2)", "CN", and "Tc" lists using sub watersheds shapefile
            sp_lst = []
            da_lst = []
            cn_lst = []
            tc_lst = []
            subshed = optfolder + "/subshed.shp"
            ss = arcpy.SearchCursor(subshed, "", "", "CurveNum;Slope;Tc;AreaMi2", "")
            for att in ss:
                sp_lst.append(att.getValue("Slope"))
                da_lst.append(att.getValue("AreaMi2"))
                cn_lst.append(att.getValue("CurveNum"))
                tc_lst.append(att.getValue("Tc"))

            # create "GAGE" and YY lists with matching length as gridcode in sub area
            GAGE = []
            NN_sub = []
            YY_sub = []
            for n in enumerate(fn_lst):
                GAGE.append("GAGE")
                NN_sub.append("NN")
                YY_sub.append("YY")

            # convert all integer/float lists into list of strings for justified text writing
            ai_lst = map(str, ai_lst)
            sp_lst = map(str, ["%0.2f" % x for x in sp_lst])
            da_lst = map(str, ["%0.2f" % x for x in da_lst])
            cn_lst = map(str, ["%0.1f" % x for x in cn_lst])
            tc_lst = map(str, ["%0.2f" % x for x in tc_lst])

            # write sub-area strings
            sub = ""
            if self.detail.GetValue() == True:
                for a, b, c, d, e, f, g in zip(ai_lst, lst_reach, GAGE, da_lst, cn_lst, tc_lst, YY_sub):
                    sub += "          " + a.ljust(10) + b.ljust(10) + c.ljust(11) + d.ljust(9) + e.ljust(10) + f.ljust(
                        10) + g + "\n"

            else:
                for a, b, c, d, e, f, g in zip(ai_lst, lst_reach, GAGE, da_lst, cn_lst, tc_lst, NN_sub):
                    sub += "          " + a.ljust(10) + b.ljust(10) + c.ljust(11) + d.ljust(9) + e.ljust(10) + f.ljust(
                        10) + g + "\n"

            # write reach strings
            NN_reach = []
            for n in enumerate(len_lst):
                NN_reach.append("NN")

            # changed "reachstring" to "reachstring_1" to create formatted string in two steps
            # (due to unequal length of lists: reach_blk and len_lst)
            reachstring = ""
            reach_str = map(str, reach_lst)
            reach_blk = ["Reach" + i for i in reach_str]
            reach_blk.insert(len(reach_blk) + 1, "Outlet")  # reach_blk = reach block

            len_lst = ["%0.1f" % x for x in len_lst]

            for h, k, l in zip(zip(reach_blk, reach_blk[1:]), reach_lst, len_lst):
                r1 = h[0]
                r2 = h[1]
                if os.path.exists(optfolder + "/xline_reach" + str(k) + ".shp"):
                    xs = "XS" + str(k)
                    reachstring += r1.rjust(16) + r2.rjust(10) + xs.rjust(7) + str(l).rjust(17 + len(str(l))) + "NN".rjust(
                        22 - len(str(l))) + "\n"  # changed log logic to length logic
                else:
                    xs = "Struct" + str(k)
                    reachstring += r1.rjust(16) + r2.rjust(10) + xs.rjust(21) + "\n"

            stormstring = ""

            ARC = str(self.ARC.GetValue())
            if self.ARF.GetValue() == True:
                RF_precip = hydro.ArealRF(year_uSpecified, prec_uSpecified, areami2)
            else:
                RF_precip = prec_uSpecified

            storm_number = len(prec_uSpecified)
            designstorm = ""
            designstorm = designstorm + "%s" "\n" % (storm_number)
            designstorm = designstorm + "%s" "\n" % (str(self.ARC.GetValue()))

            precdist = []
            prec_selected = [float(i) for i in prec_uSpecified]
            for i, p in zip(year_uSpecified, prec_selected):
                if i in cb_selected:
                    durlist = ["05n", "10n", "15n", "30n", "01", "02", "03", "06", "12", "24", "48"]
                    yearlist = ["1", "2", "5", "10", "25", "50", "100", "200", "500"]
                    year = yearlist[i // 4]
                    critdur = durlist[(i) % 4 + 7]
                    designstorm = designstorm + "p%s-%s" "\n" % (year, critdur)
                    designstorm = designstorm + "%s" "\n" % (critdur)

                    for d in durlist:
                        thefilename = Directory + "/data/prec/p" + str(year) + "-" + str(d) + "m"
                        precavg = arcpy.GetRasterProperties_management(thefilename, "MEAN")
                        precavg = float(precavg.getOutput(0))
                        theavg = precavg / 1000
                        precdist.append(theavg)
                        designstorm = designstorm + "%s\n" % (theavg)
                else:
                    durlist = ["05n", "10n", "15n", "30n", "01", "02", "03", "06", "12", "24", "48"]
                    yearlist = ["1", "2", "5", "10", "25", "50", "100", "200", "500"]
                    nodata_list = [-999] * 11
                    year = yearlist[i // 4]
                    critdur = durlist[(i) % 4 + 7]
                    if critdur == "06":
                        nodata_list[7] = p
                    if critdur == "12":
                        nodata_list[8] = p
                    if critdur == "24":
                        nodata_list[9] = p
                    if critdur == "48":
                        nodata_list[10] = p

                    designstorm = designstorm + "p%s-%s" "\n" % (year, critdur)
                    designstorm = designstorm + "%s" "\n" % (critdur)
                    for n in nodata_list:
                        designstorm = designstorm + "%s\n" % (n)

            if os.path.exists(optfolder + "/design_storm"):
                path = optfolder + "/design_storm"
                os.chdir(path)
                stormFN = optfolder + "/design_storm/stormdata.txt"
                stormfile = open(stormFN, "w")
                stormfile.write(designstorm)
                stormfile.close()
                os.system(path + "/design_storm_v2.exe")
                time.sleep(4)
            else:
                path = optfolder + "/design_storm"
                os.mkdir(path, 0o755)
                src = r"" + Directory + "/data/design_storm/design_storm_v2.exe"
                shutil.copy2(src, path)
                os.chdir(path)
                stormFN = optfolder + "/design_storm/stormdata.txt"
                stormfile = open(stormFN, "w")
                stormfile.write(designstorm)
                stormfile.close()
                os.system(path + "/design_storm_v2.exe")
                time.sleep(4)

            stormanalysis = optfolder + "/design_storm/stormanalysis.txt"
            storm_type = []
            with open(stormanalysis) as f:
                next(f)
                for lines in f:
                    line = lines.split()
                    storm = line[4]
                    storm_type.append(storm)

            year = []
            prcp = []
            duration = []
            for i, p, s in zip(year_uSpecified, RF_precip, storm_type):
                yearlist = ["1", "2", "5", "10", "25", "50", "100", "200", "500"]
                durlist = ["05n", "10n", "15n", "30n", "01", "02", "03", "06", "12", "24", "48"]
                theyear = yearlist[i // 4]
                year.append(theyear)
                thecritdur = durlist[(i) % 4 + 7]
                duration.append(thecritdur)
                precavg = str(p)
                prcp.append(precavg)

                if (theyear == "1") or (theyear == "2") or (theyear == "5"):
                    if float(precavg) <= 9.99:
                        if s == "Type":
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(9) + precavg.rjust(20) + "Type II".rjust(13) + ARC.rjust(4) + "\n"
                        else:
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(9) + precavg.rjust(20) + "rtp".rjust(9) + theyear.rjust(1) + "-".rjust(
                                1) + thecritdur.rjust(1) + ARC.rjust(4) + "\n"
                    else:
                        if s == "Type":
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(9) + precavg.rjust(21) + "Type II".rjust(12) + ARC.rjust(4) + "\n"
                        else:
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(9) + precavg.rjust(21) + "rtp".rjust(8) + theyear.rjust(1) + "-".rjust(
                                1) + thecritdur.rjust(1) + ARC.rjust(4) + "\n"

                elif (theyear == "10") or (theyear == "25") or (theyear == "50"):
                    if float(precavg) <= 9.99:
                        if s == "Type":
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(8) + precavg.rjust(20) + "Type II".rjust(13) + ARC.rjust(4) + "\n"
                        else:
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(8) + precavg.rjust(20) + "rtp".rjust(9) + theyear.rjust(1) + "-".rjust(
                                1) + thecritdur.rjust(1) + ARC.rjust(3) + "\n"
                    else:
                        if s == "Type":
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(8) + precavg.rjust(21) + "Type II".rjust(12) + ARC.rjust(4) + "\n"
                        else:
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(8) + precavg.rjust(21) + "rtp".rjust(8) + theyear.rjust(1) + "-".rjust(
                                1) + thecritdur.rjust(1) + ARC.rjust(3) + "\n"

                else:
                    if float(precavg) <= 9.99:
                        if s == "Type":
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(7) + precavg.rjust(20) + "Type II".rjust(13) + ARC.rjust(4) + "\n"
                        else:
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(7) + precavg.rjust(20) + "rtp".rjust(9) + theyear.rjust(1) + "-".rjust(
                                1) + thecritdur.rjust(1) + ARC.rjust(2) + "\n"
                    else:
                        if s == "Type":
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(7) + precavg.rjust(21) + "Type II".rjust(12) + ARC.rjust(4) + "\n"
                        else:
                            stormstring = stormstring + "p".rjust(11) + theyear.rjust(1) + "-".rjust(1) + thecritdur.rjust(
                                1) + "GAGE".rjust(7) + precavg.rjust(21) + "rtp".rjust(8) + theyear.rjust(1) + "-".rjust(
                                1) + thecritdur.rjust(1) + ARC.rjust(2) + "\n"

            # write rainfall distribution from designstorm text file
            rainfall_string = ""
            storm_infilename = optfolder + "/design_storm/designstorm.txt"
            with open(storm_infilename) as dist:
                next(dist)
                for lines in dist:
                    rainfall_string = rainfall_string + lines

            # write stream cross-section block from rating table folder
            rattabstring = ""
            rattabstring_res = ""
            for i in reach_lst:
                infilename = optfolder + "/rating_table/rattabout_reach" + str(i) + ".txt"
                infile = open(infilename, "r").readlines()

                header_line1 = infile[0].split()
                header_line2 = infile[1].split()
                header_line3 = infile[2].split()
                reach_drop = header_line1[0]
                reach_elev = header_line2[0]
                reach_no = header_line3[0]

                if header_line1[0][0] == "X":
                    rattab_header = "          " + reach_drop.ljust(10) + reach_elev.ljust(8) + reach_no.rjust(32)
                    rattabstring = rattabstring + rattab_header + "\n"

                    with open(infilename) as f:
                        lines_after = f.readlines()[7:]
                        for lines in lines_after:
                            line = lines.split()

                            # added to format new rattab.exe output
                            rattab = line[0].rjust(26) + "    " + line[1].ljust(10) + line[2].ljust(10) \
                                    + line[3].ljust(10) + line[4].ljust(10) + "\n"
                            rattabstring = rattabstring + rattab

                if header_line1[0][0] == "S":

                    rattab_header = "          " + reach_drop.ljust(12) + reach_elev
                    rattabstring_res = rattabstring_res + rattab_header + "\n"

                    with open(infilename) as f:
                        lines_after = f.readlines()[7:]
                        for lines in lines_after:
                            line = lines.split()

                            # added to format new rattab.exe output
                            rattab = line[0].rjust(26) + line[1].rjust(12) + line[2].rjust(12) + "\n"
                            rattabstring_res = rattabstring_res + rattab

            # Write TR20 Input file
            inputstring = ""
            inputstring = inputstring + "WinTR-20: Version 3.20" + "{:>18}{:>1}{:>9}{:>1}{:>9}{:>1}{:>7}{:>1}".format(
                    "", "0", "", "0", "", "1.0", "", "0")
            inputstring = inputstring + "" "\n"
            inputstring = inputstring + "GISHydroNXT - [folder: %sWinTR20]" "\n" % optfolder
            inputstring = inputstring + "" "\n"
            inputstring = inputstring + "SUB-AREA:" "\n"
            inputstring = inputstring + sub
            inputstring = inputstring + "" "\n"
            inputstring = inputstring + "" "\n"
            # DelMarVa Hydrograph
            if self.DelMarVa.GetValue() == "DelMarVa PRF 284":
                table_dmv = open(r"" + Directory + "/data/mdfiles/lookup/Table_DMV.txt")
                dmv_head = table_dmv.read()
                dmv_hyd = dmv_head[39:]
                inputstring = inputstring + dmv_hyd
                inputstring = inputstring + "" "\n"

            if os.path.exists(optfolder + "/rating_table"):
                inputstring = inputstring + "STREAM REACH:" "\n"
                inputstring = inputstring + reachstring
                inputstring = inputstring + "" "\n"
                inputstring = inputstring + "" "\n"

            inputstring = inputstring + "STORM ANALYSIS:" "\n"
            inputstring = inputstring + stormstring
            inputstring = inputstring + "" "\n"
            inputstring = inputstring + "" "\n"
            inputstring = inputstring + "RAINFALL DISTRIBUTION:" "\n"
            inputstring = inputstring + rainfall_string
            inputstring = inputstring + "" "\n"

            if not rattabstring == "":
                inputstring = inputstring + "STREAM CROSS SECTION:" "\n"
                inputstring = inputstring + rattabstring
                inputstring = inputstring + "" "\n"

            if not rattabstring_res == "":
                inputstring = inputstring + "STRUCTURE RATING:" "\n"
                inputstring = inputstring + rattabstring_res
                inputstring = inputstring + "" "\n"

            inputstring = inputstring + "" "\n"
            inputstring = inputstring + "" "\n"
            inputstring = inputstring + "GLOBAL OUTPUT:" "\n"
            inputstring = inputstring + "{:>20}{:>1}{:>8}{:>1}{:>7}{:>1}{:>5}{:>1}".format("", "1.", "",
                                                                                             float(Main_increment), "",
                                                                                             "YNNNN", "", "YNNNNN")
            dst = optfolder + "/WinTR20/"
            if os.path.exists(dst):
                for the_file in os.listdir(dst):
                    file_path = os.path.join(dst, the_file)
                    if os.path.isfile(file_path):
                        try:
                            os.remove(file_path)
                            shutil.rmtree(dst)
                        except:
                            pass

            path = optfolder + "/WinTR20"
            if not os.path.exists(path):
                os.mkdir(path, 0o755)
            defFN = optfolder + "/WinTR20/TR20in.txt"
            f = open(defFN, "w")
            f.write(inputstring)
            f.close()

            # open "WinTR20in" file in text editor
            # *****************************************************************************
            time.sleep(3)  # 5-24-2017: sleep time of few seconds to finish all processes before text file can be opened
            hydro.openbrowser(defFN)
            self.Show(False)

            # turn ON Win-TR20 button
            # *****************************************************************************
            button15.enabled = True
            save(optfolder)

class WatershedDelineation(object):
    """Implementation for GISHydroNXT_addin.tool2 (Tool)"""

    def __init__(self):
        self.enabled = False
        self.cursor = 3

    def onMouseDownMap(self, x, y, button, shift):
        # *******************************************************************************************************
        # Extract outlet points from point geometry
        # *******************************************************************************************************
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder
        arcpy.env.extent = "MAXOF"
        arcpy.env.addOutputsToMap = True
        outletxy = arcpy.sa.ExtractByPoints(optfolder + "/InfStreams", [arcpy.Point(x, y)], "INSIDE")
        outletxy.save(optfolder + "/outletcell")
        if arcpy.sa.Raster(optfolder + "/outletcell").maximum <= 0:
            pythonaddins.MessageBox("Please select a valid Stream", "Outlet Error")
            return

        arcpy.CopyFeatures_management(arcpy.PointGeometry(arcpy.Point(x, y)), optfolder + "/outlet_ws.shp")
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "outlet_ws":
                arcpy.mapping.RemoveLayer(df, lyr)

        # outlet from point geometry
        # *****************************************************************************
        pnt = arcpy.Point(x, y)
        global xoutletstring
        xoutletstring = str(pnt.X)

        global youtletstring
        youtletstring = str(pnt.Y)
        outletcell = arcpy.management.GetCellValue(optfolder + "/dem", "{} {}".format(x, y))

        global outletelev
        try:
            outletelev = float(outletcell.getOutput(0))
        except:
            pythonaddins.MessageBox("Please select a valid Stream", "Outlet Error")
            return
        flowdir_dir = optfolder + "/flowdir"
        flowdir_dem_dir = optfolder + "/flowdir_dem"
        watershed_dir = optfolder + "/watershed.shp"
        basin_dir = optfolder + "/basingrid"

        wshed = arcpy.sa.Watershed(flowdir_dem_dir, optfolder + "/outletcell", "VALUE")  # watershed tool execution
        shed = arcpy.sa.Con(wshed >= 0, 1, arcpy.sa.IsNull(wshed))
        shed.save(basin_dir)
        try:
            arcpy.RasterToPolygon_conversion(basin_dir, watershed_dir, "NO_SIMPLIFY", "VALUE")
        except:
            pythonaddins.MessageBox("Please select a valid Stream", "Outlet Error")
            return

        if not os.path.exists(flowdir_dir):
            arcpy.CopyRaster_management(flowdir_dem_dir, flowdir_dir)

        # *******************************************************************************************************
        # turn layers ON/OFF in current data frame
        # *******************************************************************************************************

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        try:
            wshExtent = df.extent
            wshExtent.XMin, wshExtent.YMin = extent.XMin, extent.YMin
            wshExtent.XMax, wshExtent.YMax = extent.XMax, extent.YMax
            df.extent = wshExtent
        except:
            pass
        for lyr in arcpy.mapping.ListLayers(mxd, "", df):
            if lyr.name == "flowdir":
                arcpy.mapping.RemoveLayer(df, lyr)
            if lyr.name == "watershed":
                arcpy.ApplySymbologyFromLayer_management(lyr, r"" + Directory + "/data/mdfiles/legends/watershed.lyr")

        # *******************************************************************************************************
        # turn watershed delineation OFF and Basin Composition ON
        # *******************************************************************************************************
        tool2.enabled = False
        tool3.enabled = True
        button1.enabled = True
        button2.enabled = True

        try:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd, "Layers")[0]
            lyr = arcpy.mapping.ListLayers(mxd, "watershed", df)[0]
            df.extent = lyr.getSelectedExtent()
        except:
            return
        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()
        save(optfolder)


class WriteSubAreaLandUseDistribution(object):
    """Implementation for GISHydroNXT_addin.button12 (Button)"""

    def __init__(self):
        self.enabled = False
        self.checked = False

    def onClick(self):

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            arcpy.env.scratchWorkspace = scratchfolder
            arcpy.env.workspace = optfolder
            # read datetime time stamp single time (so it should be here outside loop)
            now = datetime.now()
            month = now.strftime("%B")
            day = now.strftime("%d")
            year = now.strftime("%Y")

            # write initial lines of sub-basin composition outside for loop
            # (concatenate all sub-areas basin compositions)
            datastring = ""
            datastring = datastring + "GISHydro Release Version Date:    %s" "\n" % (Modifieddt)
            datastring = datastring + "Project Name:                     %s" % (proj)
            datastring = datastring + "" "\n"
            datastring = datastring + "Analysis Date:                    %s %s, %s " "\n" % (month, day, year)
            datastring = datastring + "" "\n"

            # create a folder for basin composition files
            sub_basincomp = optfolder + "/sub_basincomp"
            if not os.path.exists(sub_basincomp):
                os.makedirs(sub_basincomp)

            fc = optfolder + "/subshed.shp"
            desc = arcpy.Describe(fc)
            rows = arcpy.SearchCursor(fc)
            for idx, row in enumerate(rows):
                aPoly = row.getValue(desc.shapefieldname)
                setExtent = aPoly.extent
                arcpy.env.extent = setExtent

                # create mask raster for each sub-basin
                basingrid = optfolder + "/basingrid"
                arcpy.MakeFeatureLayer_management(fc, "layer" + str(idx), ' "FID" = ' + str(idx))
                mask = arcpy.Clip_management(basingrid, "#", optfolder + "/aux_folder/mask" + str(idx), "layer" + str(idx),
                                             "0", "ClippingGeometry")

                ## grids and shapefile used for sub-basin land distribution
                lu = optfolder + "/landuse"
                soil = optfolder + "/Soils"
                nhd = r"" + Directory + "/data/maryland/nhd_streamsm.shp"

                lu_ext = optfolder + "/sub_basincomp/lu_ext"
                soil_ext = optfolder + "/sub_basincomp/soil_ext"
                lu_out = arcpy.sa.Times(lu, mask)
                lu_out.save(lu_ext + str(idx))
                soil_out = arcpy.sa.Times(soil, mask)
                soil_out.save(soil_ext + str(idx))

                # 1] intersect sub-basins and NHD streams to obtain string of stream names in each sub-basin
                # 2] obtain name strings
                arcpy.env.addOutputsToMap = False
                nhd_sub = optfolder + "/sub_basincomp/nhd" + str(idx) + ".shp"
                arcpy.Intersect_analysis([nhd, "layer" + str(idx)], nhd_sub, "ALL", "#", "INPUT")
                nhd_sc = arcpy.SearchCursor(nhd_sub, "", "", "NAME", "")
                nhd_streams = []
                for s in nhd_sc:
                    nhd_names = s.getValue("NAME")
                    nhd_streams.append(nhd_names)

                # unicode to string conversion and then replace " "  in list to "Unidentified Reach"
                nhd_streams_tmp = [x.encode("UTF8") for x in nhd_streams]
                nhd_streams = ["Unidentified Reach" if i == " " else i for i in nhd_streams_tmp]  # replace "" in a list

                ## convert all clipped rasters to polygons
                lu_poly = optfolder + "/sub_basincomp/lu_poly" + str(idx) + ".shp"
                soil_poly = optfolder + "/sub_basincomp/soil_poly" + str(idx) + ".shp"
                arcpy.env.addOutputsToMap = False
                arcpy.RasterToPolygon_conversion(lu_ext + str(idx), lu_poly, "NO_SIMPLIFY", "VALUE")
                arcpy.env.addOutputsToMap = False
                arcpy.RasterToPolygon_conversion(soil_ext + str(idx), soil_poly, "NO_SIMPLIFY", "VALUE")

                ## intersect land use and soil to prepare two polygons: "lu_soil" and "lu_cn"
                arcpy.env.addOutputsToMap = False
                lu_soil = optfolder + "/sub_basincomp/lu_soil" + str(idx) + ".shp"
                arcpy.Intersect_analysis([lu_poly, soil_poly], lu_soil, "ALL", "#", "INPUT")

                ## dissolve above intersected polygon
                arcpy.env.addOutputsToMap = False
                lu_soil_diss = optfolder + "/sub_basincomp/lu_soil_diss" + str(idx) + ".shp"
                arcpy.Dissolve_management(lu_soil, lu_soil_diss, "GRIDCODE;GRIDCODE_1", "#", "MULTI_PART", "DISSOLVE_LINES")

                ## add filed to both of above dissolved polygons and compute area in acres
                if not len(arcpy.ListFields(lu_soil_diss, "area")) > 0:
                    arcpy.AddField_management(lu_soil_diss, "area", "FLOAT", 15, 4)
                arcpy.CalculateField_management(lu_soil_diss, "area", "!shape.area@acres!", "PYTHON")

                # prepre a list of lu codes to feed into lu_description function in order to obtain matching descriptions list
                lu_match = []
                sc = arcpy.SearchCursor(lu_ext + str(idx), "", "", "VALUE", "")
                for i in sc:
                    v = i.getValue("VALUE")
                    lu_match.append(v)

                # create list of lists with zeroes
                soil_acre_lists = [[0, 0, 0, 0] for i in range(len(lu_match))]

                # preapre a list of soil acreage using lu_match list
                lc_soil_diss = []
                soil_lc_diss = []
                lc_soil_aa = []
                sr = arcpy.SearchCursor(lu_soil_diss, "", "", "GRIDCODE;GRIDCODE_1;area", "")
                for s in sr:
                    lc = s.getValue("GRIDCODE")
                    lc_soil_diss.append(lc)
                    sc = s.getValue("GRIDCODE_1")
                    soil_lc_diss.append(sc)
                    aa = s.getValue("area")
                    lc_soil_aa.append(round(aa, 2))

                for n, lu in enumerate(lu_match):
                    for l, s, a in zip(lc_soil_diss, soil_lc_diss, lc_soil_aa):
                        if l == lu:
                            soil_acre_lists[n][int(s) - 1] = a

                # prepare matching list of lu description using lu codes from lu raster of watershed
                if landuse == "NLCD 2011":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookuppoor.txt"

                if landuse == "NLCD 2006":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookuppoor.txt"

                if landuse == "NLCD 2001":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/nlcdlookuppoor.txt"

                if landuse == "1997 MOP":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookuppoor.txt"

                if landuse == "2002 MOP":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookuppoor.txt"

                if landuse == "2010 MOP":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/andlookuppoor.txt"

                if landuse == "2002 MD/DE":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/mddelookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/mddelookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/mddelookuppoor.txt"

                if landuse == "Ultimate":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/zoninglookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/zoninglookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/zoninglookuppoor.txt"

                if landuse == "MRLC":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/mrlclookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/mrlclookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/mrlclookuppoor.txt"

                if landuse == "1970s USGS":
                    if hyd == "Fair":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/usgslookupfair.txt"
                    elif hyd == "Good":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/usgslookupgood.txt"
                    elif hyd == "Poor":
                        lut_file = r"" + Directory + "/data/mdfiles/lookup/usgslookuppoor.txt"

                # run hydro function to obtain land use description of categories present in watershed
                lu_desc = hydro.lu_description(lut_file, lu_match)

                # Perform two tasks:
                # a) take "lu_desc" and "soil_acre_lists" and concatenate them
                # b) format according to basin composition file in legacy version
                land_soil_area = ""
                width = 30
                for ld, sg in zip(lu_desc, soil_acre_lists):
                    land_soil_area = land_soil_area + "{: <{}}".format(ld, width) + str(sg[0]).rjust(10) + str(sg[1]).rjust(
                        10) + str(sg[2]).rjust(10) + str(sg[3]).rjust(10)
                    land_soil_area = land_soil_area + "" "\n"

                # *******************************************************************************************************
                # Text file string variables
                # *******************************************************************************************************
                datastring = datastring + "" "\n"
                datastring = datastring + "Landuse and Soil Distributions for: Sub-Area %s" "\n" % (
                    str(idx + 1))  # 1/2/2018: because we want sub-basin index to start from 1 so "idx+1"
                datastring = datastring + "" "\n"
                datastring = datastring + "Streams located in this sub-area:"
                datastring = datastring + "" "\n"

                # add stream names
                for sn in nhd_streams:
                    datastring = datastring + str(sn)
                    datastring = datastring + "" "\n"

                datastring = datastring + "" "\n"
                datastring = datastring + "Distribution of Landuse by Soil Group" "\n"
                datastring = datastring + "" "\n"
                datastring = datastring + "Acres on Indicated Soil Group".rjust(66)
                datastring = datastring + "" "\n"
                datastring = datastring + "Land Use".rjust(8) + "A-Soil".rjust(32) + "B-Soil".rjust(10) + "C-Soil".rjust(
                    10) + "D-Soil".rjust(10)
                datastring = datastring + "" "\n"
                datastring = datastring + "" "\n"
                datastring = datastring + land_soil_area

                # sum list of lists separately and cat at the end of lu description
                total_area = [sum(i) for i in zip(*soil_acre_lists)]
                datastring = datastring + "{: <{}}".format("Total Area:", width) + str(total_area[0]).rjust(10) + str(
                    total_area[1]).rjust(10) + str(total_area[2]).rjust(10) + str(total_area[3]).rjust(10)

                datastring = datastring + "" "\n"
                datastring = datastring + "" "\n"
                datastring = datastring + "Distribution of Land Use and Curve Numbers Used" "\n"
                datastring = datastring + "" "\n"
                datastring = datastring + "Land Use".rjust(8) + "Acres".rjust(34) + "Percent".rjust(10) + "A".rjust(
                    4) + "B".rjust(4) + "C".rjust(4) + "D".rjust(4)
                datastring = datastring + "" "\n"
                datastring = datastring + "" "\n"

                # loop over land use, related total acreage, percent of land covered by this lu category, and A-B-C-D curve numbers
                curve_num = []
                for l in lu_match:
                    with open(lut_file, "r") as f:
                        next(f)
                        for line in f:
                            luc = line.split("\t")[0]
                            if int(l) == int(luc):
                                temp = []
                                A = line.split("\t")[2]  # CN A
                                temp.append(A)
                                B = line.split("\t")[3]  # CN B
                                temp.append(B)
                                C = line.split("\t")[4]  # CN C
                                temp.append(C)
                                D = line.split("\t")[5]  # CN D
                                temp.append(D)
                                curve_num.append(temp)

                # sum areas for each sub-list individually
                acres = [sum(i) for i in soil_acre_lists]
                total_all = sum(total_area)
                percent = [round(float(x / total_all) * 100, 2) for x in acres]
                for l, a, p, cn in zip(lu_desc, acres, percent, curve_num):
                    datastring = datastring + "{: <{}}".format(l, width) + str(a).rjust(12) + str(p).rjust(10) + str(
                        cn[0]).rjust(4) + str(cn[1]).rjust(4) + str(cn[2]).rjust(4) + str(cn[3]).rjust(4)
                    datastring = datastring + "" "\n"

                # *******************************************************************************************************
                # turn layers ON/OFF in current data frame
                # *******************************************************************************************************
                mxd = arcpy.mapping.MapDocument("CURRENT")
                df = arcpy.mapping.ListDataFrames(mxd)[0]
                layers = arcpy.mapping.ListLayers(mxd, "", df)
                for lyr in layers:
                    if lyr.name == "lu_ext":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "soil_ext":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "lu_poly":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "soil_poly":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "lu_soil":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "lu_soil_diss":
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "layer" + str(idx):
                        arcpy.mapping.RemoveLayer(df, lyr)
                    if lyr.name == "mask" + str(idx):
                        arcpy.mapping.RemoveLayer(df, lyr)

                arcpy.RefreshTOC()
                arcpy.RefreshActiveView()

            # *******************************************************************************************************
            # write strings to basin stat text file.
            # Message box containing datastring as message
            # *******************************************************************************************************
            defFN = optfolder + "/sub_basincomp.txt"
            compfile = open(defFN, "w")
            compfile.write(datastring)
            compfile.close()

            # *******************************************************************************************************
            # open "basincomp" file in text editor
            # *******************************************************************************************************
            hydro.openbrowser(defFN)

            # *******************************************************************************************************
            # turn peak discharge OFF
            # *******************************************************************************************************
            button12.enabled = False
            arcpy.env.extent = "MAXOF"
            save(optfolder)

class Xeditor(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Cross Section Editor", size=(590, 390))
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        panel = wx.Panel(self, -1)
        wx.StaticBox(panel, -1, "Transect Line Geometry", (20, 15), size=(260, 120))
        wx.StaticText(panel, -1, "Transect Line Width:", (40, 35))
        wx.StaticText(panel, -1, str(TWL), (200, 35))
        wx.StaticText(panel, -1, "ft", (255, 35))
        wx.StaticText(panel, -1, "Maximum Elevation:", (40, 60))
        wx.StaticText(panel, -1, str(TWE_max), (200, 60))
        wx.StaticText(panel, -1, "ft", (255, 60))
        wx.StaticText(panel, -1, "Minimum Elevation:", (40, 85))
        wx.StaticText(panel, -1, str(TWE_min), (200, 85))
        wx.StaticText(panel, -1, "ft", (255, 85))
        wx.StaticText(panel, -1, "Upstream Discharge Area:", (40, 110))
        wx.StaticText(panel, -1, str(areami2_usda), (200, 110))
        wx.StaticText(panel, -1, "mi^2", (245, 110))

        wx.StaticBox(panel, -1, "Reach Characteristics", (20, 140), size=(260, 90))
        wx.StaticText(panel, -1, "Reach Slope", (40, 170))
        self.reachslope = wx.TextCtrl(panel, -1, value=str(reachslope), pos=(150, 165), size=(75, 25))
        wx.StaticText(panel, -1, "ft/ft", (245, 170))
        wx.StaticText(panel, -1, "Bankfull Elevation", (40, 200))
        self.Ebf = wx.TextCtrl(panel, -1, value=str(TWE_min), pos=(150, 195),
                               size=(75, 25))  # Ebf = Bank Full Elevation
        wx.StaticText(panel, -1, "ft", (245, 200))

        wx.StaticBox(panel, -1, "Channel Geometry", (20, 235), size=(260, 90))
        wx.StaticText(panel, -1, "Bankful Channel Width", (40, 260))
        self.Wbf = wx.TextCtrl(panel, -1, value=str(Wbf), pos=(180, 255),
                               size=(70, 25))  # Wbf = Bank Full Channel Width
        wx.StaticText(panel, -1, "ft", (260, 260))
        wx.StaticText(panel, -1, "Bankful Channel Depth", (40, 290))
        self.Dbf = wx.TextCtrl(panel, -1, value=str(Dbf), pos=(180, 285),
                               size=(70, 25))  # Dbf = Bank Full Channel Depth
        wx.StaticText(panel, -1, "ft", (260, 290))

        wx.StaticBox(panel, -1, "Roughness Characteristics", (290, 15), size=(260, 140))
        wx.StaticText(panel, -1, "Main Channel n Value", (310, 35))
        self.nMainValue = wx.TextCtrl(panel, -1, value="0.050", pos=(460, 30), size=(70, 25))  # Main Channel n Value
        wx.StaticText(panel, -1, "Left* Overbank n Value", (310, 65))
        self.nLeftValue = wx.TextCtrl(panel, -1, value="0.100", pos=(460, 60), size=(70, 25))  # Left Overbank n Value
        wx.StaticText(panel, -1, "Right* Overbank n Value", (310, 95))
        self.nRightValue = wx.TextCtrl(panel, -1, value="0.100", pos=(460, 90), size=(70, 25))  # Right Overbank n Value
        wx.StaticText(panel, -1, "* Facing Downstream", (400, 120))

        self.btnExport = wx.Button(panel, label="Export Cross Section", pos=(290, 240))
        self.Bind(wx.EVT_BUTTON, self.OnExport, id=self.btnExport.GetId())

        self.btnPlot = wx.Button(panel, label="Plot Cross Section", pos=(440, 240))
        self.Bind(wx.EVT_BUTTON, self.OnPlot, id=self.btnPlot.GetId())

        self.btnApply = wx.Button(panel, label="Apply", pos=(290, 290))
        self.Bind(wx.EVT_BUTTON, self.OnApply, id=self.btnApply.GetId())

        self.btnClose = wx.Button(panel, label="Close", pos=(440, 290))
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=self.btnClose.GetId())

        self.Show(True)
        self.Centre(True)
        self.style = self.GetWindowStyle()
        self.SetWindowStyle(self.style | wx.STAY_ON_TOP)

    def OnClose(self, event):
        # ******************************************************************************************************
        # delete layers to avoid duplication conflict
        # ******************************************************************************************************

        if os.path.exists(optfolder + "/trans_line.shp"):
            arcpy.Delete_management(optfolder + "/trans_line.shp", "")
        if os.path.exists(optfolder + "/transect.shp"):
            arcpy.Delete_management(optfolder + "/transect.shp", "")
        if os.path.exists(optfolder + "/reachslope.shp"):
            arcpy.Delete_management(optfolder + "/reachslope.shp", "")
        if os.path.exists(optfolder + "/rating_table/ratTabIn.txt"):
            os.remove(optfolder + "/rating_table/ratTabIn.txt")
        if os.path.exists(optfolder + "/rating_table/userout.txt"):
            os.remove(optfolder + "/rating_table/userout.txt")
        if os.path.exists(optfolder + "/rating_table/rattabout.txt"):
            os.remove(optfolder + "/rating_table/rattabout.txt")
        if os.path.exists(optfolder + "/Elevation_Stage_Profile"):
            arcpy.Delete_management(optfolder + "/Elevation_Stage_Profile", "")
        if os.path.exists(optfolder + "/aux_folder/intersect_aux.shp"):
            arcpy.Delete_management(optfolder + "/aux_folder/intersect_aux.shp", "")

        # ******************************************************************************************************
        # open "rating table output" file in text editor
        # [No. of Xsections = No. of Reaches]
        # ******************************************************************************************************
        fn_lst = []
        tn_lst = []
        subriver = optfolder + "/subrivers.shp"
        sr = arcpy.SearchCursor(subriver, "", "", "FROM_NODE;TO_NODE", "")
        for node in sr:
            fn = node.getValue("FROM_NODE")
            tn = node.getValue("TO_NODE")
            fn_lst.append(fn)
            tn_lst.append(tn)

        # use reach numbers to name "userout" and "rattabout" files
        reach_lst = [x for x in fn_lst if x in tn_lst]

        # ******************************************************************************************************
        # turn Xeditor OFF and Precipitation Depths ON
        # ******************************************************************************************************
        # To make sure we have all the rating table text files for WinTR20 input text file
        # Please make sure that text files pop up after transect tool have values NOT blank
        myPath = optfolder + "/rating_table/"
        txtCounter = len(glob.glob1(myPath, "*.txt"))
        reachCounter = len(reach_lst)
        if txtCounter == reachCounter:
            user_prompt = pythonaddins.MessageBox("Process finished, continue?", "Reach Tool", 4)
            if user_prompt == "Yes":
                tool6.enabled = False
                tool8.enabled = False
                button13.enabled = True
                button25.enabled = True
                pythonaddins.MessageBox("Process finished!" + "\n" +
                                        "Continue to Precipitation Depths selection.", "Reach Tool")
        try:
            self.frm.Show(False)
        except:
            pass
        self.Show(False)  # close dialog box after deleting files and activating precip depth button

    def OnPlot(self, event):

        ## get elevation and flow accumulation values along transect line
        length = []
        elev = []
        arcpy.StackProfile_3d(optfolder + '/trans_line.shp', optfolder + '/dem', optfolder + '/Elevation_Stage_Profile')
        cursor = arcpy.da.SearchCursor(optfolder + '/Elevation_Stage_Profile', ['FIRST_DIST', 'FIRST_Z'])
        for row in cursor:
            length.append(row[0] * 3.28084)
            elev.append(row[1])

        wbf_val = self.Wbf.GetValue()
        wbf = float(wbf_val)
        dbf_val = self.Dbf.GetValue()
        dbf = float(dbf_val)
        u = 0.65

        minindx = elev.index(min(elev))
        bfelev = min(elev)
        alpha = (wbf / 2) / (dbf ** u)
        xcentr = length[minindx]
        del elev[minindx]
        del length[minindx]

        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 4)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 3)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 2)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 1)
        elev.insert(minindx, bfelev - dbf)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 1)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 2)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 3)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 4)

        length.insert(minindx, xcentr + (alpha * (((dbf / 4) * 4) ** u)))
        length.insert(minindx, xcentr + (alpha * (((dbf / 4) * 3) ** u)))
        length.insert(minindx, xcentr + (alpha * (((dbf / 4) * 2) ** u)))
        length.insert(minindx, xcentr + (alpha * (((dbf / 4) * 1) ** u)))
        length.insert(minindx, xcentr)
        length.insert(minindx, xcentr - (alpha * (((dbf / 4) * 1) ** u)))
        length.insert(minindx, xcentr - (alpha * (((dbf / 4) * 2) ** u)))
        length.insert(minindx, xcentr - (alpha * (((dbf / 4) * 3) ** u)))
        length.insert(minindx, xcentr - (alpha * (((dbf / 4) * 4) ** u)))

        self.data = zip(length, elev)

        # "self.data" holds both x and y lists for plot
        self.frm = wx.Frame(self, -1, "Elevation stage profile", size=(600, 450))
        client = plot.PlotCanvas(self.frm)
        line = plot.PolyLine(self.data, legend="", colour="red", width=1)
        gc = plot.PlotGraphics([line], "Height Cross-section", "Station (ft)", "Elevation (ft)")
        client.Draw(gc, xAxis=(min(length) - 20, max(length) + 20), yAxis=(min(elev) - 10, max(elev) + 10))
        self.frm.Show(True)

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        tables = arcpy.mapping.ListTableViews(mxd, "", df)
        for tbl in tables:
            if tbl.name == "Elevation_Stage_Profile":
                arcpy.mapping.RemoveTableView(df, tbl)
                arcpy.Delete_management(tbl)

    def OnExport(self, event):

        length = []
        elev = []
        arcpy.StackProfile_3d(optfolder + "/trans_line.shp", optfolder + "/dem", optfolder + "/Elevation_Stage_Profile")
        cursor = arcpy.da.SearchCursor(optfolder + "/Elevation_Stage_Profile", ['FIRST_DIST', 'FIRST_Z'])
        for row in cursor:
            length.append(row[0] * 3.28084)
            elev.append(row[1])

        wbf_val = self.Wbf.GetValue()
        wbf = float(wbf_val)
        dbf_val = self.Dbf.GetValue()
        dbf = float(dbf_val)
        u = 0.65

        minindx = elev.index(min(elev))
        bfelev = min(elev)
        alpha = (wbf / 2) / (dbf ** u)
        xcentr = length[minindx]
        del elev[minindx]
        del length[minindx]

        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 4)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 3)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 2)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 1)
        elev.insert(minindx, bfelev - dbf)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 1)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 2)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 3)
        elev.insert(minindx, (bfelev - dbf) + (dbf / 4) * 4)

        length.insert(minindx, xcentr + (alpha * (((dbf / 4) * 4) ** u)))
        length.insert(minindx, xcentr + (alpha * (((dbf / 4) * 3) ** u)))
        length.insert(minindx, xcentr + (alpha * (((dbf / 4) * 2) ** u)))
        length.insert(minindx, xcentr + (alpha * (((dbf / 4) * 1) ** u)))
        length.insert(minindx, xcentr)
        length.insert(minindx, xcentr - (alpha * (((dbf / 4) * 1) ** u)))
        length.insert(minindx, xcentr - (alpha * (((dbf / 4) * 2) ** u)))
        length.insert(minindx, xcentr - (alpha * (((dbf / 4) * 3) ** u)))
        length.insert(minindx, xcentr - (alpha * (((dbf / 4) * 4) ** u)))

        xmod = ["%0.2f" % i for i in length]
        ymod = ["%0.2f" % i for i in elev]

        # use reach numbers to name "userout" and "rattabout" files
        in_target = optfolder + "/trans_line.shp"
        in_join = optfolder + "/subrivers.shp"
        out_feature_class = optfolder + "/aux_folder/intersect_aux.shp"
        arcpy.SpatialJoin_analysis(in_target, in_join, out_feature_class)
        sr = arcpy.SearchCursor(out_feature_class, "", "", "ARCID", "")
        for node in sr:
            rating = int(node.getValue("ARCID"))

        # ******************************************************************************************************
        # Write elevation stage profile and distance to text file
        # ******************************************************************************************************
        datastring = ""
        for a, b in zip(xmod, ymod):
            datastring = datastring + "{:<6}{:>7}".format(a, b)
            datastring = datastring + "\n"

        if os.path.exists(optfolder + "/elev_stage_profile"):
            elev_stage = optfolder + "/elev_stage_profile/elev_stage_reach" + str(rating) + ".txt"
            elev_stage_file = open(elev_stage, "w")
            elev_stage_file.write(datastring)
            elev_stage_file.close()
        else:
            path = optfolder + "/elev_stage_profile"
            os.mkdir(path, 0o755)
            elev_stage = optfolder + "/elev_stage_profile/elev_stage_reach" + str(rating) + ".txt"
            elev_stage_file = open(elev_stage, "w")
            elev_stage_file.write(datastring)
            elev_stage_file.close()

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        tables = arcpy.mapping.ListTableViews(mxd, "", df)
        for tbl in tables:
            if tbl.name == "Elevation_Stage_Profile":
                arcpy.mapping.RemoveTableView(df, tbl)
                arcpy.Delete_management(tbl)
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "intersect_aux":
                arcpy.mapping.RemoveLayer(df, lyr)
                arcpy.Delete_management(lyr)

        pythonaddins.MessageBox("Elevation stage profile is successfully exported to %s" % elev_stage,
                                "Elevation Stage Proile", 0)

    def OnApply(self, event):
        self.Show(False)
        arcpy.env.scratchWorkspace = scratchfolder
        arcpy.env.workspace = optfolder

        with pythonaddins.ProgressDialog as dialogprogress:
            dialogprogress.title = "Loading"
            dialogprogress.description = "GISHydroNXT is working, please wait..."
            dialogprogress.animation = "Spiral"

            # ******************************************************************************************************
            # get reach slope, Wbf, Dbf, Manning numbers, reach number, stream it drains to number,
            # distance, and elevation to include in rating table input text file
            # ******************************************************************************************************
            rs_val = float(self.reachslope.GetValue())
            wbf_val = float(self.Wbf.GetValue())
            dbf_val = float(self.Dbf.GetValue())
            nMain_val = float(self.nMainValue.GetValue())
            nLeft_val = float(self.nLeftValue.GetValue())
            nRight_val = float(self.nRightValue.GetValue())

            global rating
            in_target = optfolder + "/trans_line.shp"
            in_join = optfolder + "/subrivers.shp"
            out_feature_class = optfolder + "/aux_folder/intersect_aux.shp"
            arcpy.SpatialJoin_analysis(in_target, in_join, out_feature_class)
            sr = arcpy.SearchCursor(out_feature_class, "", "", "ARCID", "")
            for node in sr:
                rating = int(node.getValue("ARCID"))

            # ******************************************************************************************************
            # run rating table algorithm (replaced rattab.exe)
            # ******************************************************************************************************

            mean_aux = []
            elev = []
            length = []
            for i in range(len(trans_dem)-1):
                if trans_dem[i] == trans_dem[i+1]:
                    mean_aux.append(trans_distance[i])
                else:
                    elev.append(trans_dem[i])
                    mean_aux.append(trans_distance[i])
                    length.append(sum(mean_aux) / float(len(mean_aux)))
                    mean_aux = []

            u = 0.65
            minindx = elev.index(min(elev))
            bfelev = min(elev)
            alpha = (wbf_val / 2) / (dbf_val ** u)
            xcentr = length[minindx]
            del elev[minindx]
            del length[minindx]

            p  = (((dbf_val / 4) * 1 - (dbf_val / 4) * 0)**2 +
                  (wbf_val / 2 - alpha * (((dbf_val / 4) * 3) ** u) - (wbf_val / 2 - alpha * (((dbf_val / 4) * 4) ** u)))**2)**0.5
            p += (((dbf_val / 4) * 2 - (dbf_val / 4) * 1)**2 +
                  (wbf_val / 2 - alpha * (((dbf_val / 4) * 2) ** u) - (wbf_val / 2 - alpha * (((dbf_val / 4) * 3) ** u)))**2)**0.5
            p += (((dbf_val / 4) * 3 - (dbf_val / 4) * 2)**2 +
                  (wbf_val / 2 - alpha * (((dbf_val / 4) * 1) ** u) - (wbf_val / 2 - alpha * (((dbf_val / 4) * 2) ** u)))**2)**0.5
            p += (((dbf_val / 4) * 4 - (dbf_val / 4) * 3)**2 +
                  (wbf_val / 2 - alpha * (((dbf_val / 4) * 0) ** u) - (wbf_val / 2 - alpha * (((dbf_val / 4) * 1) ** u)))**2)**0.5
            pch = p*2

            alpha*dbf_val**u - alpha * (((dbf_val / 4) * 3) ** u)

            elev.insert(minindx, (bfelev - dbf_val) + (dbf_val / 4) * 4)
            elev.insert(minindx, (bfelev - dbf_val) + (dbf_val / 4) * 3)
            elev.insert(minindx, (bfelev - dbf_val) + (dbf_val / 4) * 2)
            elev.insert(minindx, (bfelev - dbf_val) + (dbf_val / 4) * 1)
            elev.insert(minindx, (bfelev - dbf_val) + (dbf_val / 4) * 0)
            elev.insert(minindx, (bfelev - dbf_val) + (dbf_val / 4) * 1)
            elev.insert(minindx, (bfelev - dbf_val) + (dbf_val / 4) * 2)
            elev.insert(minindx, (bfelev - dbf_val) + (dbf_val / 4) * 3)
            elev.insert(minindx, (bfelev - dbf_val) + (dbf_val / 4) * 4)

            length.insert(minindx, xcentr + (alpha * (((dbf_val / 4) * 4) ** u)))
            length.insert(minindx, xcentr + (alpha * (((dbf_val / 4) * 3) ** u)))
            length.insert(minindx, xcentr + (alpha * (((dbf_val / 4) * 2) ** u)))
            length.insert(minindx, xcentr + (alpha * (((dbf_val / 4) * 1) ** u)))
            length.insert(minindx, xcentr + (alpha * (((dbf_val / 4) * 0) ** u)))
            length.insert(minindx, xcentr - (alpha * (((dbf_val / 4) * 1) ** u)))
            length.insert(minindx, xcentr - (alpha * (((dbf_val / 4) * 2) ** u)))
            length.insert(minindx, xcentr - (alpha * (((dbf_val / 4) * 3) ** u)))
            length.insert(minindx, xcentr - (alpha * (((dbf_val / 4) * 4) ** u)))

            n = 30  #number of iterations
            maxleft = max(elev[0:elev.index(min(elev))])
            maxright = max(elev[elev.index(min(elev))+1:])
            dx = (min(maxleft, maxright) - min(elev) + dbf_val)/float(n)

            stage = [min(elev)]
            discharge = [0]
            endarea = [0]
            topwidth = [0]

            for i in range(n):
                areal = 0
                periml = 0
                arear = 0
                perimr = 0
                width = 0
                bar = min(elev) + dbf_val + dx*i
                for j in range(len(elev)-1):
                    xj = 0
                    ajl = 0
                    ajr = 0
                    pjl = 0
                    pjr = 0
                    diff = bar - elev[j+1]
                    if length[j] < length[elev.index(min(elev))]:
                        if diff > 0:
                            if elev[j] >= bar:
                                yj = diff
                                xj = (length[j+1]-length[j])*yj/(elev[j]-elev[j+1])
                                ajl = yj*xj*0.5
                                pjl = (yj**2 + xj**2)**0.5
                            else:
                                yj = (diff + bar - elev[j])/2
                                xj = length[j+1] - length[j]
                                ajl = yj*xj
                                pjl = ((elev[j]-elev[j+1])**2 + (xj)**2)**0.5
                        elif elev[j] < bar:
                                yj = bar - elev[j]
                                xj = (length[j+1]-length[j])*yj/(elev[j+1]-elev[j])
                                ajl = yj*xj*0.5
                                pjl = (yj**2 + xj**2)**0.5
                    else:
                        if diff > 0:
                            if elev[j] >= bar:
                                yj = diff
                                xj = (length[j+1]-length[j])*yj/(elev[j]-elev[j+1])
                                ajr = yj*xj*0.5
                                pjr = (yj**2 + xj**2)**0.5
                            else:
                                yj = (diff + bar - elev[j])/2
                                xj = length[j+1] - length[j]
                                ajr = yj*xj
                                pjr = ((elev[j]-elev[j+1])**2 + (xj)**2)**0.5
                        elif elev[j] < bar:
                                yj = bar - elev[j]
                                xj = (length[j+1]-length[j])*yj/(elev[j+1]-elev[j])
                                ajr = yj*xj*0.5
                                pjr = (yj**2 + xj**2)**0.5
                    width += xj
                    areal += ajl
                    periml += pjl
                    arear += ajr
                    perimr += pjr

                if width > 0 and bar < max(elev):
                    hr = (areal + arear)/(periml + perimr)
                    ncomp = ((periml * (nLeft_val)**(1.5) + pch * nMain_val**(1.5) + perimr * (nRight_val)**(1.5)) / (periml + perimr))**(0.66667) ## Equal velocity method

                    stage.append(bar)
                    discharge.append(1.49 * (areal + arear) * hr**(0.66667) * rs_val**0.5 / ncomp)
                    endarea.append(areal + arear)
                    topwidth.append(width)

                    if discharge[-1] < discharge[-2]:
                        discharge[-1] = discharge[-2] + 1

            if not os.path.exists(optfolder + "/rating_table"):
                path = optfolder + "/rating_table"
                os.mkdir(path, 0o755)

            defFN = optfolder + "/rating_table/rattabout_reach" + str(rating) + ".txt"

            discharge = [float(elem) for elem in discharge]
            discharge.sort()

            ST = ['%.1f' % elem for elem in stage]
            DI = ['%.1f' % elem for elem in discharge]
            EN = ['%.1f' % elem for elem in endarea]
            TO = ['%.1f' % elem for elem in topwidth]
            RS = ["%.4f" % rs_val]*len(ST)

            datastring = "XS" + str(rating) + "\n"
            datastring = datastring + str(min(stage) + dbf_val) + "\n"
            datastring = datastring + "Reach" + str(rating) + "\n"
            datastring = datastring + "\n"
            datastring = datastring + '{0[0]:<12}{0[1]:<12}{0[2]:<12}{0[3]:<12}{0[4]:<12}'.format(
                    ["Stage", "Discharge", "End-Area", "Topwidth", "Slope"]) + "\n"
            datastring = datastring + '{0[0]:<12}{0[1]:<12}{0[2]:<12}{0[3]:<12}{0[4]:<12}'.format(
                    ["-----", "-----", "-----", "-----", "-----"]) + "\n"
            datastring = datastring + "\n"

            rows = zip(ST, DI, EN, TO, RS)
            for row in rows:
                datastring = datastring + '{0[0]:<12}{0[1]:<12}{0[2]:<12}{0[3]:<12}{0[4]:<12}'.format(row) + "\n"

            ratfile = open(defFN, "w")
            ratfile.write(datastring)
            ratfile.close()

            # ******************************************************************************************************
            # open "rating table output" file in text editor
            # [No. of Xsections = No. of Reaches]
            # ******************************************************************************************************

            arcpy.CopyFeatures_management(optfolder + "/transect.shp", optfolder + "/xpoint_reach" + str(rating) + ".shp")
            arcpy.CopyFeatures_management(optfolder + "/line_trun.shp", optfolder + "/xline_reach" + str(rating) + ".shp")

            # ******************************************************************************************************
            # turn layers ON/OFF in current data frame
            # ******************************************************************************************************
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = arcpy.mapping.ListLayers(mxd, "", df)
            for lyr in layers:
                if lyr.name == "xpoint_reach" + str(rating):
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "intersect_aux":
                    arcpy.mapping.RemoveLayer(df, lyr)
                if lyr.name == "trans_line":
                    arcpy.mapping.RemoveLayer(df, lyr)
            tables = arcpy.mapping.ListTableViews(mxd, "", df)
            for tbl in tables:
                if tbl.name == "Elevation_Stage_Profile":
                    arcpy.mapping.RemoveTableView(df, tbl)
                    arcpy.Delete_management(tbl)
            if os.path.exists(optfolder + "/aux_folder/intersect_aux.shp"):
                arcpy.Delete_management(optfolder + "/aux_folder/intersect_aux.shp")

            arcpy.RefreshTOC()
            arcpy.RefreshActiveView()

        hydro.openbrowser(defFN)

        self.OnClose(event)


def save(optfolder):
    ### SAVE PROGRESS FUNCTION
    proj_ = proj.replace(" ", "_")
    mxd = arcpy.mapping.MapDocument("CURRENT")
    if os.path.exists(optfolder + "/" + timestr + "_" + proj_ + ".mxd"):
        os.remove(optfolder + "/" + timestr + "_" + proj_ + ".mxd")
    mxd.saveACopy(optfolder + "/" + timestr + "_" + proj_ + ".mxd")

    savelines = str(tool1.enabled) + "\n"
    savelines = savelines + str(tool2.enabled) + "\n"
    savelines = savelines + str(tool3.enabled) + "\n"
    savelines = savelines + str(tool4.enabled) + "\n"
    savelines = savelines + str(tool5.enabled) + "\n"
    savelines = savelines + str(tool6.enabled) + "\n"
    savelines = savelines + str(tool7.enabled) + "\n"
    savelines = savelines + str(tool8.enabled) + "\n"
    savelines = savelines + str(button1.enabled) + "\n"
    savelines = savelines + str(button2.enabled) + "\n"
    savelines = savelines + str(button3.enabled) + "\n"
    savelines = savelines + str(button4.enabled) + "\n"
    savelines = savelines + str(button5.enabled) + "\n"
    savelines = savelines + str(button5_1.enabled) + "\n"
    savelines = savelines + str(button6.enabled) + "\n"
    savelines = savelines + str(button7.enabled) + "\n"
    savelines = savelines + str(button8.enabled) + "\n"
    savelines = savelines + str(button9.enabled) + "\n"
    savelines = savelines + str(button10.enabled) + "\n"
    savelines = savelines + str(button11.enabled) + "\n"
    savelines = savelines + str(button12.enabled) + "\n"
    savelines = savelines + str(button13.enabled) + "\n"
    savelines = savelines + str(button14.enabled) + "\n"
    savelines = savelines + str(button15.enabled) + "\n"
    savelines = savelines + str(button16.enabled) + "\n"
    savelines = savelines + str(button17.enabled) + "\n"
    savelines = savelines + str(button18.enabled) + "\n"
    savelines = savelines + str(button19.enabled) + "\n"
    savelines = savelines + str(button20.enabled) + "\n"
    savelines = savelines + str(button21.enabled) + "\n"
    savelines = savelines + str(button23.enabled) + "\n"
    savelines = savelines + str(button24.enabled) + "\n"
    savelines = savelines + str(button25.enabled)

    #saving button/tool state
    savefile = optfolder + "/saveconfig1.txt"
    if os.path.exists(savefile):
        os.remove(savefile)
    saveas = open(savefile, "w")
    saveas.write(savelines)
    saveas.close()

    #saving variables
    filename = optfolder + '/savevar.out'
    globals_ = globals()
    my_shelf = shelve.open(filename, 'n')
    for key, value in globals_.items():
        if not key.startswith('__'):
            try:
                my_shelf[key] = value
            except:
                pass
    my_shelf.close()

    # save layers
    savelines2 = ""
    savelines3 = ""
    if not os.path.exists(optfolder + "/lyrs"):
        os.mkdir(optfolder + "/lyrs", 0o755)
    mxd = arcpy.mapping.MapDocument("CURRENT")
    layers = arcpy.mapping.ListLayers(mxd)
    for lyr in layers:
        try:
            if lyr.supports("dataSource"): # some layers might not support the property "dataSource"
                savelines2 = savelines2 + str(lyr.dataSource) + "\n"
            savelines3 = savelines3 + str(lyr.name) + "\n"
            lyr.saveACopy(optfolder + "/lyrs/" + lyr.name + ".lyr")()
        except:
            pass
    savelines2 = savelines2 + "end"
    savelines3 = savelines3 + "end"

    savefile = optfolder + "/saveconfig2.txt"
    if os.path.exists(savefile):
        os.remove(savefile)
    saveas = open(savefile, "w")
    saveas.write(savelines2)
    saveas.close()

    savefile = optfolder + "/saveconfig3.txt"
    if os.path.exists(savefile):
        os.remove(savefile)
    saveas = open(savefile, "w")
    saveas.write(savelines3)
    saveas.close()

def load(optfolder):

    with pythonaddins.ProgressDialog as dialogprogress:
        dialogprogress.title = "Loading"
        dialogprogress.description = "Loading project, please wait..."
        dialogprogress.animation = "Spiral"

        #load saved variables
        filename = r"" + optfolder + '/savevar.out'
        try:
            my_shelf = shelve.open(filename)
        except:
            pythonaddins.MessageBox("Incorrect User", "Error Loading Project")
            return
        for key in my_shelf:
            try:
                if isinstance(my_shelf[key], bool):
                    globals()[key]=my_shelf[key]
                if isinstance(my_shelf[key], str):
                    globals()[key]=my_shelf[key]
                if isinstance(my_shelf[key], int):
                    globals()[key]=my_shelf[key]
                if isinstance(my_shelf[key], float):
                    globals()[key]=my_shelf[key]
                if isinstance(my_shelf[key], list):
                    globals()[key]=my_shelf[key]
                if isinstance(my_shelf[key], bytes):
                    globals()[key]=my_shelf[key]
            except:
                pass
        my_shelf.close()

        # load toolbar status
        loadlist = []
        savefile = optfolder + "/saveconfig1.txt"
        loadlist_aux = open(savefile,'r').read().split('\n')
        for loadbool in loadlist_aux:
            if 'True' in loadbool:
                loadlist.append(True)
            else:
                loadlist.append(False)

        i = 0
        tool1.enabled = loadlist[i]; i+=1
        tool2.enabled = loadlist[i]; i+=1
        tool3.enabled = loadlist[i]; i+=1
        tool4.enabled = loadlist[i]; i+=1
        tool5.enabled = loadlist[i]; i+=1
        tool6.enabled = loadlist[i]; i+=1
        tool7.enabled = loadlist[i]; i+=1
        tool8.enabled = loadlist[i]; i+=1
        button1.enabled = loadlist[i]; i+=1
        button2.enabled = loadlist[i]; i+=1
        button3.enabled = loadlist[i]; i+=1
        button4.enabled = loadlist[i]; i+=1
        button5.enabled = loadlist[i]; i+=1
        button5_1.enabled = loadlist[i]; i+=1
        button6.enabled = loadlist[i]; i+=1
        button7.enabled = loadlist[i]; i+=1
        button8.enabled = loadlist[i]; i+=1
        button9.enabled = loadlist[i]; i+=1
        button10.enabled = loadlist[i]; i+=1
        button11.enabled = loadlist[i]; i+=1
        button12.enabled = loadlist[i]; i+=1
        button13.enabled = loadlist[i]; i+=1
        button14.enabled = loadlist[i]; i+=1
        button15.enabled = loadlist[i]; i+=1
        button16.enabled = loadlist[i]; i+=1
        button17.enabled = loadlist[i]; i+=1
        button18.enabled = loadlist[i]; i+=1
        button19.enabled = loadlist[i]; i+=1
        button20.enabled = loadlist[i]; i+=1
        button21.enabled = loadlist[i]; i+=1
        button23.enabled = loadlist[i]; i+=1
        button24.enabled = loadlist[i]; i+=1
        button25.enabled = loadlist[i]

        # remove all layers
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        tables = arcpy.mapping.ListTableViews(mxd, "", df)
        for lyr in layers:
            arcpy.mapping.RemoveLayer(df, lyr)
        for tbl in tables:
            arcpy.mapping.RemoveTableView(df, tbl)

        # add saved layers
        savefile = optfolder + "/saveconfig2.txt"
        loadlist2 = open(savefile,'r').read().split('\n')
        savefile = optfolder + "/saveconfig3.txt"
        loadlist3 = open(savefile,'r').read().split('\n')

        for a,b in zip(loadlist2,loadlist3):
            try:
                try:
                    arcpy.MakeFeatureLayer_management(a,"temp")
                except:
                    arcpy.MakeRasterLayer_management(a,"temp")
                layer = arcpy.mapping.Layer("temp")
                layer.name = b
                arcpy.ApplySymbologyFromLayer_management(layer, r"" + optfolder + "/lyrs/" + b + ".lyr")
            except:
                pass
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        layers = arcpy.mapping.ListLayers(mxd, "", df)
        for lyr in layers:
            if lyr.name == "Soils":
                lyr.visible = False
            if lyr.name == "soils":
                lyr.visible = False
            if lyr.name == "Landuse":
                lyr.visible = False
            if lyr.name == "landuse":
                lyr.visible = False
            if lyr.name == "AddasReservoir":
                lyr.visible = False

        try:
            lyr = arcpy.mapping.ListLayers(mxd, "dem", df)[0]
            df.extent = lyr.getSelectedExtent()
        except:
            pass

class wxPython(object):
    """Implementation for GISHydroNXT_addin.wxPython (Extension)"""

    def __init__(self):
        # For performance considerations, please remove all unused methods in this class.
        self.enabled = True

    def startup(self):
        try:
            from wx import PySimpleApp
            self._wxApp = PySimpleApp()
            self._wxApp.MainLoop()
        except:
            sMsg = "Error starting extension:\n" + traceback.format_exc()
            pythonaddins.MessageBox(sMsg, "TestAddIn")

# ******************************************************************************************************
# Initialize buttons to avoid "global name definition error"
# ******************************************************************************************************
tool1 = AreaOfInterest()
tool2 = WatershedDelineation()
tool3 = ControlBoundaries()
tool4 = FlowPaths()
tool5 = AddSubwatershedOutlets()
tool6 = TransectLine()
tool7 = LandUseEditor()
tool8 = Reservoir()
button1 = ResetWatershed()
button2 = BasinComposition()
button3 = BasinStatistics()
button4 = ThomasDischarge()
button5_1 = CompareDischarges()
button5 = ResetSubWatershed()
button6 = AddStreams()
button7 = AddOutlets()
button8 = DelineateSubwatersheds()
button9 = SetTcParameters()
button10 = CalculateAttributes()
button11 = CombineLongestFlowPathSegments()
button12 = WriteSubAreaLandUseDistribution()
button13 = PrecipitationDepths()
button14 = ControlPanel()
button15 = ExecuteTR20()
button16 = CreateContours()
button17 = HelpHTML()
button18 = HelpTechRef()
button19 = HelpTrainingManual()
button20 = HelpManualToolGuide()
button21 = HelpAbout()
button23 = LoadButton()
button24 = ResetGlobal()
button25 = ResetControlPanel()



