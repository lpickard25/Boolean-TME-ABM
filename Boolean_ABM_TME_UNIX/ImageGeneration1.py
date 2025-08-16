import pandas as pd
import numpy as np
import math
import sys
import os
import re
from PIL import Image, ImageDraw

# This script generates the .jpg images of the cancer cells from the last day of each set in a given folder
# as well as the .jpg images of all the cells w/ vessels from each day of set number 1,
# as well as the .jpg images of onlt the cancer cells w/ vessels from each day of set number 1,
# the user inputs the desired folder, the number of the last day, and the number of sets in that folder
#


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

    brightRed = (255, 34, 68);
    brownRed = (200, 90, 50);

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
            #return WINE;
            return brownRed;
        case -2:
            return brightRed;
        case -1:
            #return WHITE;
            return BLACK;
        case 0:
            return INDIGO
            #return lightGreen
        case 1:
            return CYAN
            #return lightGreen
        case 2:
            return TEAL
            #return darkGreen
        case 3: #cancer
            return GREEN
            #return darkBlue
        case 4: #cd4 t helper
            return ROSE
            #return lightRed
        case 5: #cd4 T reg
            return WINE
            #return darkRed
        case 6: # naive Tcell
            return PURPLE
            #return lightOrange
        case 7: # exhausted/suppressed
            return SAND
            #return darkOrange
        case 8: # promemory
            return PALE_GREY
            #return lightPurple
        case 9:
            return darkPurple
        case 10:
            return yellow
        case _:
            return WHITE



# function which draws the all the cells in the given dataframe (image object already made)
# inputs are a list of the cells, the image object, height and width of image, and scale
def drawCells(data_inorder, image, width_um:float, height_um:float,scale:int):
    i = 0
    # iterates over each cell row in the list
    while i < len(data_inorder):
        xLoc = width_um / 2 + data_inorder.at[i, "xloc"] / scale
        yLoc = height_um / 2 - data_inorder.at[i, "yloc"] / scale
        # changing CD8 T cell states to numerical values
        # E = exhausted, N = naive, M = promemory
        if data_inorder.at[i, "phenotype"] == "E":
            data_inorder.at[i, "phenotype"] = 7
        if  data_inorder.at[i, "phenotype"] == "N":
            data_inorder.at[i, "phenotype"] = 6
        if data_inorder.at[i, "phenotype"] == "M":
            data_inorder.at[i, "phenotype"] = 8

        color = setColor(int(data_inorder.at[i, "phenotype"]))
        radius = data_inorder.at[i, "radius"] / scale
        image.circle((xLoc, yLoc), radius, fill=color, outline=color)
        i += 1


## function to create, draw, and save one image, from given dataframe automatically deciding width and height
# based on maximum x/y locations of cells in given dataframe
def makeImage(mydata, folder:str, scale:int, save:str):

    # determine max x and y to set dimensions
    maxX = max(abs(mydata["xloc"]))
    maxY = max(abs(mydata["yloc"]))

    # Set the desired dimensions in micrometers (µm)
    # if scale is 1 then 1 pixel = 1 um
    # if scale is 2 then 1 pixel = 2 um
    width_um = int((maxX*2+200)/scale)
    height_um = int((maxY*2+200)/scale)

    # Create a new image with the specified dimensions
    im = Image.new("RGB", (width_um, height_um), "black")
    image = ImageDraw.Draw(im)

    drawCells(mydata,image,width_um,height_um,scale)
    saveImage(folder,save,im)

# function to save an image object to a jpg
def saveImage(folder:str,name:str,im):
    saveName =  folder + "/" + name + ".jpg"
    im.save(saveName)

# function to read in desired data, given the folder it is located in and the name of the file
# also sorts the cells into order of drawing, assuming all cells are drawn
# assumes columns = ["cell_type", "xloc", "yloc", "radius", "phenotype", "pdl1"]
def readData(folder:str, name:str):

    fileName =  folder + "/" + name + ".csv"
    
    mydata = pd.read_csv(fileName)

    #mydata = pd.read_csv(fileName, header=None)
    #mydata.columns = ["cell_type", "xloc", "yloc", "radius", "phenotype", "pdl1"]
    # put data in order of cells to draw
    data_inorder = mydata[mydata["cell_type"]==0].copy() # cancer cells
    data_inorder.reset_index(inplace=True)
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"] == -2]], ignore_index=True) # vessels
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"] == 1]], ignore_index=True) # macrophage
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"] == 5]], ignore_index=True)
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"] == 2]], ignore_index=True) # CD4 t cells
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"] == 3]], ignore_index=True) # CD8 t cells
    data_inorder = pd.concat([data_inorder, mydata[mydata["cell_type"] == 4]], ignore_index=True)
    return data_inorder



## function to make and save one image, with width and height previously determined
# necessary to ensure that images of different days from same set are the same size
# used in timeCourseImages()
def makeImage2(mydata, folder:str, name:str, scale:int, width:float, height:float):

    # Create a new image with the specified dimensions
    im = Image.new("RGB", (int(width), int(height)), "black")
    image = ImageDraw.Draw(im)

    drawCells(mydata,image,width,height,scale)

    saveImage(folder,name,im)
    return im

# function to read in data and extract only one cell type + vessels to draw
# used in timeCourseImagesOneCell()
def plotOneCell(folder:str, name:str, cellState:int, width_um, height_um, scale:int):
    fileName = folder + "/" + name + ".csv"
    mydata = pd.read_csv(fileName)

    # extract only the cells of the desired phenotype (aka cellState)
    cellData = mydata[np.any([mydata["phenotype"] == cellState, mydata["phenotype"] == str(cellState)], axis=0)].copy()
    cellData.reset_index(inplace=True)
    # add back the vessels
    cellData = pd.concat([cellData, mydata[mydata["cell_type"] == -2]], ignore_index=True)

    im = Image.new("RGB", (width_um, height_um), "black")
    image = ImageDraw.Draw(im)

    drawCells(cellData,image,width_um,height_um,scale)

    newName = name + "_cellState_" + str(cellState)
    newFolder = folder + "/cellState_" +str(cellState)
    saveName =  newFolder + "/" + newName + ".jpg"
    im.save(saveName)

# function to draw jpgs of every cell for every day of a given simulation folder
# determines width and height of all images based on final (largest) image
def timeCourseImages(folderName:str,dataRoot:str,scale:int,days:int,outFolder:str):
    fileName = folderName + "/" + dataRoot + str(days) +".csv"
    mydata = pd.read_csv(fileName)

    maxY = max(abs(mydata["yloc"]))
    maxX = max(abs(mydata["xloc"]))
    width_um = float((maxX * 2 + 200) / scale)
    height_um = float((maxY * 2 + 200) / scale)

    for i in range(0,int(days)+1):
        dataName = dataRoot + str(i)
        outputName = "day_"+str(i)+"_cells"
        mydata = readData(folderName, dataName)
        makeImage2(mydata, outFolder, outputName, scale, width_um, height_um)

# function to draw jpgs of one cell type for every day of a given simulation folder
# determines width and height of all images based on final (largest) image
def timeCourseImagesOneCell(folderName:str,dataRoot:str, days, cellState:int,scale:int):
    fileName = folderName + "/" + dataRoot + str(days) + ".csv"
    mydata = pd.read_csv(fileName)
    maxX = max(abs(mydata["xloc"]))
    maxY = max(abs(mydata["yloc"]))
    width_um = int((maxX * 2 + 200) / scale)
    height_um = int((maxY * 2 + 200) / scale)
 
    try:
        os.mkdir(folderName + "/cellState_" + str(cellState))
        print(f"Directory '{savefolder}' created successfully.")
    except FileExistsError:
        print(f"Directory '{savefolder}' already exists.")
    for i in range(0,int(days)+1):
        mydata1 = dataRoot + str(i)
        plotOneCell(folderName, mydata1, cellState, width_um, height_um, scale)



# User input: folder in which the sets are contained, lastday value, number of sets
fld = sys.argv[1] # like "./NewVesselRuns/ControlRun"
lastDay = re.sub(r'[\x00-\x1F\x7F]', '', str(sys.argv[2]))
sets = re.sub(r'[\x00-\x1F\x7F]', '', str(sys.argv[3]))

scale = 2
cellState = 3
# name of the save folder, created within the input folder
savefolder = fld + "/Day" + str(lastDay)+ "Tumors"

#create the save folder
try:
    os.mkdir(savefolder)
    print(f"Directory '{savefolder}' created successfully.")
except FileExistsError:
    print(f"Directory '{savefolder}' already exists.")

# iterate from 1 to the set # given in the input folder
for i in range(1,int(sets)+1):
    # get the final day of the cancer cells only
    # this naming system is based on current set up of the model I used
    folder = fld + "/set_" + str(i) + "/cellLists"
    name = "cells_day_" + str(lastDay)
    mydata = readData(folder, name)
    cellData = mydata[np.any([mydata["phenotype"] == cellState, mydata["phenotype"] == str(cellState)], axis=0)].copy()
    cellData.reset_index(inplace=True)
    newName = name + "_set_ " + str(i) + "_tumor"
    makeImage(cellData, savefolder, scale, newName)

# get representative time series of set 1
# both all the cells w/ vessels and only cancer cells w/ vessels
folder = fld +"/set_1"
outFolder = folder + "/images"
try:
    os.mkdir(outFolder)
    print(f"Directory '{outFolder}' created successfully.")
except FileExistsError:
    print(f"Directory '{outFolder}' already exists.")
    
basename = "cells_day_"
# only cancer cells w/ vessels
timeCourseImagesOneCell(folder,basename, lastDay, 3, scale)
# all cells w/ vessels
timeCourseImages(folder,basename,scale,lastDay,outFolder)

