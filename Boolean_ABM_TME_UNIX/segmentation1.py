import matplotlib.pyplot as plt
import sys
from sklearn.datasets import make_circles
from sklearn.cluster import DBSCAN
import numpy as np
import pandas as pd
import os
import math
import alphashape
from PIL import Image, ImageDraw
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely.geometry import Point

'''
this script averages the vessel concentration and immune cell concentration in the tumor and stroma compartments 
for each day over ten replicates of a given simulation. the user inputs the folder containing the 10 sets
and the maximum day the simulation runs till
The script outputs a csv file containing the averages for each day for vessel and cell concentrations in the tumor
and stroma compartments as well as a png graph charting theve averages (4 lines total)
'''

# simple function to read in a csv to a pandas dataframe
def readData(folder:str, name:str):
    fileName = folder + "/" + name + ".csv"
    mydata = pd.read_csv(fileName)
    return mydata

# determines a label for clusters of cancer cells, all other cells labelled -1
def addsegmentation(mydata):
    cancerCells = mydata[mydata["cell_type"]==0].copy()
    locations = cancerCells[["xloc","yloc"]].to_numpy()
    dbscan = DBSCAN(eps=100, min_samples=50)
    clusters = dbscan.fit_predict(locations)
    cancerCells["labels"] = clusters

    otherCells = mydata[mydata["cell_type"] != 0].copy()
    otherCells["labels"] = -1

    # Combine both dataframes
    combined = pd.concat([cancerCells, otherCells], ignore_index=True)

    return combined


# creates tumor and stroma region as shape objects,
# tumor region is based on cancer cells (labelled not -1, if code is working)
# stroma region is a buffer around tumor region
def defineRegions(mydata):
    locations = mydata[["xloc", "yloc"]].to_numpy()
    unique_labels = mydata['labels'].unique().tolist()
    area = 0
    tumor_shapes = []
    buffer_shapes = []
    for label in unique_labels:
        if label != -1:
            clusterpoints = locations[mydata["labels"] == label]
            if len(clusterpoints) > 2:
                alpha = 0.009
                tumor_shape = alphashape.alphashape(clusterpoints, alpha).buffer(0)
                buffer_shape = tumor_shape.buffer(1000).buffer(0)
                #stroma_ring = buffer_shape.difference(tumor_shape)
                tumor_shapes.append(tumor_shape)
                buffer_shapes.append(buffer_shape)
                #print("Tumor valid?", tumor_shape.is_valid)
                #print("Buffer valid?", buffer_shape.is_valid)

    tumor_union = unary_union(tumor_shapes)
    buffer_union = unary_union([buffer_shapes])
    stroma_union = buffer_union.difference(tumor_union)
    #print("Ring valid?", stroma_union.is_valid)
    #print("Stroma intersects tumor?", stroma_union.intersects(tumor_union))
    return tumor_union, stroma_union


def addOutline(mydata, draw, tumor_shape, stroma_shape, width, height, scale):
    def draw_outline(shape, color):
        if isinstance(shape, Polygon):
            shapes = [shape]
        elif isinstance(shape, MultiPolygon):
            shapes = list(shape.geoms)
        else:
            print("error")
            return

        for poly in shapes:
            x, y = poly.exterior.xy
            outline_points = [(width / 2 + (px / scale), height / 2 - (py / scale)) for px, py in zip(x, y)]
            draw.line(outline_points + [outline_points[0]], fill=color, width=4)


    def fill_shape(shape, fill_color):
        if isinstance(shape, Polygon):
            shapes = [shape]
        elif isinstance(shape, MultiPolygon):
            shapes = list(shape.geoms)
        else:
            print("error")
            return
        for poly in shapes:
            x, y = poly.exterior.xy
            coords = [(width / 2 + (px / scale), height / 2 - (py / scale)) for px, py in zip(x, y)]
            draw.polygon(coords, fill=fill_color)
            for hole in poly.interiors:
                hole_coords = [(width / 2 + x / scale, height / 2 - y / scale) for x, y in hole.coords]
                draw.polygon(hole_coords, fill="black")

    #print("Stroma intersects tumor?", stroma_shape.intersects(tumor_shape))
    fill_shape(stroma_shape, (173, 216, 230, 100))
    draw_outline(tumor_shape, "white")
    draw_outline(stroma_shape, "red")



# function determines whether active vessels are in the tumor or stroma region or neither
# calculates the tumor area and stroma area
# and calculates the vessel concentraiton in each of those area
def calcActiveVesselConc(mydata, tumor_shape, stroma_shape):
    activeVessels = mydata[np.any([mydata["phenotype"] == -2, mydata["phenotype"] == str(-2)], axis=0)].copy()
    activeVessels.reset_index(inplace=True)


    tumorActiveVessels  = activeVessels[activeVessels.apply(lambda row: tumor_shape.contains(Point(row["xloc"], row["yloc"])), axis=1)]
    stromaActiveVessels = activeVessels[activeVessels.apply(lambda row: stroma_shape.contains(Point(row["xloc"], row["yloc"])), axis=1)]

    tumorarea = tumor_shape.area
    stromaarea = stroma_shape.area

    # concentration in vessels/mm^, so have to convert um^2 to mm^2 (hence the 1000000)
    tumorvesselconc = len(tumorActiveVessels) * (1000000) / tumorarea
    stromavesselconc = len(stromaActiveVessels)* (1000000) / stromaarea

    #print("Active Vessel Number in the stroma: " + str(len(stromaActiveVessels)))

    #print("Active Vessel Number in the tumor: " + str(len(tumorActiveVessels)))
    return tumorvesselconc, stromavesselconc


# function determines whether immune cells are in the tumor or stroma region or neither
# calculates the tumor area and stroma area
# and calculates the immune cell concentraiton in each of those area
def calcImmuneConc(mydata, tumorshape, stromashape):
    immuneCells = mydata[
        np.any([mydata["cell_type"] == 1, mydata["cell_type"] == 2, mydata["cell_type"] == 3], axis=0)].copy()
    immuneCells.reset_index(inplace=True)

    stromaimmune= immuneCells[immuneCells.apply(lambda row: stromashape.contains(Point(row["xloc"], row["yloc"])), axis=1)]
    tumorimmune = immuneCells[immuneCells.apply(lambda row: tumorshape.contains(Point(row["xloc"], row["yloc"])), axis=1)]
    # concentration in vessels/mm^, so have to convert um^2 to mm^2 (hence the 1000000)

    tumorarea = tumorshape.area
    stromaarea = stromashape.area

    tumorimmuneconc = len(tumorimmune) * (1000000) / tumorarea
    stromaimmuneconc = len(stromaimmune) * (1000000) / stromaarea

    return tumorimmuneconc, stromaimmuneconc
    #print("Immune Cell Number in the stroma: " + str(len(stromaimmune)))

    #print("Immune Cell Number in the tumor: " + str(len(tumorimmune)))


# function to create a jpg image of one day, including the outlines demarcating tumor and stroma regions
def makeImage2(mydata, folder:str, name:str, scale:int, tumor_shape, stroma_shape):
    maxX = max(abs(mydata["xloc"]))
    maxY = max(abs(mydata["yloc"]))
    width = int((maxX * 2 + 200) / scale)
    height = int((maxY * 2 + 200) / scale)
    # put data in order of cells to draw
    data_inorder = mydata[mydata["cell_type"]==0]
    data_inorder.reset_index(inplace=True)
    xloc = data_inorder["xloc"].mean()
    yloc = data_inorder["yloc"].mean()
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"]==-2]],ignore_index = True)
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"]==1]],ignore_index = True)
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"]==2]],ignore_index = True)
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"]==3]],ignore_index = True)
    #data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"]==4]],ignore_index = True)

    # Create a new image with the specified dimensions
    im = Image.new("RGBA", (int(width), int(height)), "black")
    image = ImageDraw.Draw(im, "RGBA")

    addOutline(mydata, image, tumor_shape, stroma_shape, width, height, scale)
    drawCells(data_inorder,image,width,height,scale)


    saveName = folder + "/Segmentation_" + name + ".png"
    im.save(saveName)

# function returns a color based on the state (aka phenotype) of a cell
'''
-3 = collapsed vessel//vessel that does not allow extravasation
-2 = active vessel
-1 = dead cell
0 = M0 macrophages
1 = M1 macrophages
2 = M2 macrophages
3 = Cancer
4 = CD4 T helper cell
5 = CD4 T regulatory cell
6 = naive CD8 T cell
7 = exhausted/suppressed CD8 T cell
8 = Promemory CD8 T cell 
9 = 
10 = 
'''
def setColor(cs:int):
    INDIGO = (51, 34, 136);
    CYAN = (136, 204, 238);
    TEAL = (64, 170, 153);
    GREEN = (17, 119, 51);
    OLIVE = (152, 153, 51);
    SAND = (221, 204, 119);
    ROSE = (204, 102, 119);
    WINE = (136, 34, 85);
    PURPLE = (170, 68, 153);
    PALE_GREY = (221, 221, 221);
    BLACK = (0, 0, 0);
    WHITE = (255, 255, 255);

    brightRed = (220, 56, 58);
    brownRed = (130, 50, 25);

    lightBlue = (166,206,227);
    darkBlue = (31,120,180);
    lightGreen = (178,223,138);
    darkGreen = (51,160,44);
    lightRed = (251,154,153);
    darkRed = (227, 26, 28);
    lightOrange = (253,191,111);
    darkOrange = (255,127,0);
    lightPurple = (202,178,214);
    darkPurple = (106,61,154);
    yellow = (255,255,153);

    match cs:
        case -3:
            return brownRed;
        case -2:
            return brightRed;
        case -1:
            #return WHITE;
            return BLACK;
        case 0:
            #return INDIGO
            return lightGreen
        case 1:
            #return CYAN
            return lightGreen
        case 2:
            #return TEAL
            return darkGreen
        case 3:
            #return GREEN
            return darkBlue
        case 4:
            #return ROSE
            return lightRed
        case 5:
            #return WINE
            return darkRed
        case 6: # naive Tcell
            #return PURPLE
            return lightOrange
        case 7: # exhausted/suppressed
            #return SAND
            return darkOrange
        case 8: # promemory
            #return PALE_GREY
            return lightPurple
        case 9:
            return darkPurple
        case 10:
            return yellow
        case _:
            return WHITE


# function which draws the all the cells in the given dataframe (image object already made)
def drawCells(data_inorder, image, width_um:float, height_um:float,scale:int):
    i = 0
    while i < len(data_inorder):
        xLoc = width_um / 2 + data_inorder.at[i, "xloc"] / scale
        #xLoc = data_inorder.at[i, "xloc"] / scale
        yLoc = height_um / 2 - data_inorder.at[i, "yloc"] / scale
        #yLoc =  data_inorder.at[i, "yloc"] / scale
        if data_inorder.at[i, "phenotype"] == "E":
            data_inorder.at[i, "phenotype"] = 7
        if  data_inorder.at[i, "phenotype"] == "N":
            data_inorder.at[i, "phenotype"] = 6
        if data_inorder.at[i, "phenotype"] == "M":
            data_inorder.at[i, "phenotype"] = 6
        #print(int(data_inorder.at[i, "phenotype"]))
        color = setColor(int(data_inorder.at[i, "phenotype"]))
        radius = data_inorder.at[i, "radius"] / scale
        image.circle((xLoc, yLoc), radius, fill=color, outline=color)

        i += 1

# function to get average vessel/immune cell concentration in stroma and tumor areas for each day over 10 replicate sets
def getAvgTSData(days,folder, reps:int):
    tumorvesselAvgs = [0] * (days+1)
    stromavesselAvgs = [0] * (days+1)
    tumorAvgs = [0] * (days+1)
    stromaAvgs = [0] * (days+1)
    # iterates over each set
    for i in range(1,reps+1):
        #print("Rep: " + str(i))
        folderpath = folder + "/set_" + str(i) + "/cellLists"
        # gets the vessel/immune cell concentration from each day from that set
        tumorvesselConcs, stromavesselConcs, tumorConcs, stromaConcs = getTSVesselData(folderpath, int(days))
        # add these new concentrations to the previously summed concentrations
        tumorvesselAvgs[:] = [x + y for x, y in zip(tumorvesselAvgs, tumorvesselConcs)]
        stromavesselAvgs[:] = [x + y for x, y in zip(stromavesselAvgs, stromavesselConcs)]
        tumorAvgs[:] = [x + y for x, y in zip(tumorAvgs, tumorConcs)]
        stromaAvgs[:] = [x + y for x, y in zip(stromaAvgs, stromaConcs)]
    #divided these sums by the number of sets
    stromaAvgs[:] = [x / reps for x in stromaAvgs]
    tumorAvgs[:] = [x / reps for x in tumorAvgs]
    stromavesselAvgs[:] = [x / reps for x in stromavesselAvgs]
    tumorvesselAvgs[:] = [x / reps for x in tumorvesselAvgs]

    return tumorvesselAvgs, stromavesselAvgs, tumorAvgs, stromaAvgs

# function to segment all the days of a given set, returning lists of day to day vessel/immune cell concentrations for tumor and
# stroma compartments
def getTSVesselData(folder, days:int):
    tumorvesselConcs = []
    stromavesselConcs = []
    tumorConcs = []
    stromaConcs = []
    cancerCellList = []
    # iterates over the days
    for i in range(days+1):
        print(i)
        name = "cells_day_"+str(i)
        # reads the data
        mydata = readData(folder, name)
        cancerCells = mydata[mydata["cell_type"]==0].copy()
        # segments based on cancer cells, adds a labels column to the data
        mydata = addsegmentation(mydata)

        unique_labels = mydata['labels'].unique().tolist()
        # tumor and stroma shape objects
        tumor_shape, stroma_shape = defineRegions(mydata)
        # draws the image, with the tumor and stroma shapes
        makeImage2(mydata, folder, name, 2, tumor_shape, stroma_shape)
        # calculates concentration of active vessels in tumor and stroma based on tumor and stroma shape
        tumorvesselconc, stromavesselconc = calcActiveVesselConc(mydata, tumor_shape, stroma_shape)
        # calculates concentration of immune cells in tumor and stroma based on tumor and stroma shape
        tumorconc, stromaconc = calcImmuneConc(mydata, tumor_shape, stroma_shape)
        # appends these values to a list, one list for each concentration, each value representing the next day
        tumorvesselConcs.append(tumorvesselconc)
        stromavesselConcs.append(stromavesselconc)
        tumorConcs.append(tumorconc)
        stromaConcs.append(stromaconc)
        # also saving number of cancer cells each day
        cancerCellList.append(len(cancerCells))

    # turning individual lists into one data frame to save to csv
    dataFrame = pd.DataFrame(
        {"Days": range(days + 1), "TumorVesselConc": tumorvesselConcs, "StromaVesselConc": stromavesselConcs})
    saveName = folder + "/vesselTS.csv"
    dataFrame.to_csv(saveName)
    dataFrame1 = pd.DataFrame(
        {"Days": range(days + 1), "TumorImmuneConc": tumorConcs, "StromaImmuneConc": stromaConcs,
         "Cancer": cancerCellList})
    saveName1 = folder + "/immuneTS.csv"
    dataFrame1.to_csv(saveName1)

    # returns individual lists
    return tumorvesselConcs, stromavesselConcs, tumorConcs, stromaConcs

# function to plot and save a png of the vessel/immune cell concentrations in the tumor or stroma region over time
def plotTS(days, folder, tumorvesselConcs, stromavesselConcs, tumorConcs, stromaConcs):
    dataFrame = pd.DataFrame(
        {"Days": range(int(days) + 1), "TumorVesselConc": tumorvesselConcs, "StromaVesselConc": stromavesselConcs})
    saveName =  folder + "/vesselTS" + str(days) + ".csv"
    dataFrame.to_csv(saveName)
    dataFrame1 = pd.DataFrame(
        {"Days": range(int(days) + 1), "TumorImmuneConc": tumorConcs, "StromaImmuneConc": stromaConcs})
    saveName1 = folder + "/immuneTS" + str(days) + ".csv"
    dataFrame1.to_csv(saveName1)

    plt.plot(tumorvesselConcs,label="Tumor Vessels")
    plt.plot(stromavesselConcs, label="Stroma Vessels")

    plt.legend()
    plt.title(f"Active Vessel Concentration Over Time for {folder}")
    plt.xlabel("Time (Days)")
    plt.ylabel("Vessel Concentration (vessel/mm^2)")
    vesselfigureName = folder+"/vesselTS" + str(days) + ".png"
    plt.savefig(vesselfigureName)
    plt.clf()

    plt.plot(tumorConcs, label="Tumor Immune Cells")
    plt.plot(stromaConcs, label="Stroma Immune Cells")

    plt.legend()
    plt.title(f"Immune Cell Concentration Over Time for {folder}")
    plt.xlabel("Time (Days)")
    plt.ylabel("Immune Concentration (cell/mm^2)")
    figureName = folder + "/immuneTS" + str(days) + ".png"
    plt.savefig(figureName)
    plt.clf()




# User input: folder name, maximum day of simulation
# assumes there are 10 sets of replicates for a given simulation
folder = sys.argv[1]
days = int(sys.argv[2])
reps = 10
scale = 2


tumorVesselAvgs, stromaVesselAvgs, tumorAvgs, stromaAvgs = getAvgTSData(days, folder, reps)
plotTS(days, folder, tumorVesselAvgs, stromaVesselAvgs, tumorAvgs, stromaAvgs)


