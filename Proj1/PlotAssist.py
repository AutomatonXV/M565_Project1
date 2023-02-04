import os

import matplotlib.pyplot as plt 
import colorsys
import matplotlib.ticker as mtick
import matplotlib.font_manager
import numpy as np

#pip install matplotlib-label-lines
from labellines import labelLine, labelLines

class HigsPlot:
    def __init__(self,xscale = 'linear', yscale = 'linear', xnot = 'default', ynot = 'default',):
        #x or yscale can be set to 'log'
        #x or ynot sets the notation mode, 'scientific' (uses e^), 'default', 'decimal'
        self.Projected3D = False, 
    
        self.XLabel = None,
        self.YLabel = None,
        self.ZLabel = None,
        self.Title = None, #No need

        self.XTick = {'Start':None, 'End':None, 'Mid':None}
        self.YTick = {'Start':None, 'End':None, 'Mid':None}
        #Plot Limits
        self.SetLimits = False
        self.Left = None ; self.Right = None,
        self.Top = None ; self.Bottom = None,

        #notation
        self.XNot = xnot
        self.YNot = ynot

        

        #Install the font, print(matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf'))
        #find the file dir in which the font is located in the print, copy paste it into font_dir
        font_dir = ["C:\\Users\\SinhA\\AppData\\Local\\Microsoft\\Windows\\Fonts"]
        for font in matplotlib.font_manager.findSystemFonts(font_dir):
            matplotlib.font_manager.fontManager.addfont(font)
        #print(matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf'))

        # Set font family globally
        matplotlib.rcParams['mathtext.fontset'] = 'custom' 
        matplotlib.rcParams['mathtext.rm'] = 'XCharter' #Roman
        matplotlib.rcParams['mathtext.it'] = 'XCharter:italic' #italic
        matplotlib.rcParams['mathtext.bf'] = 'XCharter:bold' #bold
        matplotlib.rcParams['font.family'] = 'XCharter'
        matplotlib.rcParams['font.size'] = 12 

        self.plt = plt        
        self.fig = plt.figure()
        self.ax = plt.gca()
        self.ax.set_xscale(xscale)
        self.ax.set_yscale(yscale)
        self.ax.minorticks_on()

        
        #Experimental, havent played much with this.
        if self.XNot == 'decimal':
            print("entering dec")
            self.ax.get_xaxis().set_major_formatter(
            mtick.FormatStrFormatter('%.1f')) #change f to e for scientific
        elif self.XNot == 'scientific':
            self.ax.get_xaxis().set_major_formatter(
            mtick.FormatStrFormatter('%.1e')) #change f to e for scientific

        if self.YNot == 'decimal':
            self.ax.get_yaxis().set_major_formatter(
            mtick.FormatStrFormatter('%.1f')) #change f to e for scientific
        elif self.YNot == 'scientific':
            self.ax.get_yaxis().set_major_formatter(
            mtick.FormatStrFormatter('%.1e')) #change f to e for scientific
        #
        
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)

    def Plot(self,XY,Label = "", LineWidth = 1.5, Color = None, LineStyle = '-', Marker = 'None', MarkerSize = 1.5):
        if not Color: raise ValueError("Missing EZColor 'Color' property.")   #Color is a EZColor object
        X,Y = XY
        plt.plot(X,Y, label = Label, linewidth = LineWidth, color = Color.RGB, linestyle = LineStyle, marker = Marker, markersize = MarkerSize)

    def Mode3D(self, ):
        self.plt.axis('off')
        fig = plt.figure() 
        ax = fig.add_subplot(111,projection='3d') 
        #ax = self.plt.axes(projection = '3d')
        #self.ax = ax
        

    def Plot3D(self, XYZ, Label = "", LineWidth = 1.5, Color = None, LineStyle = '-', Marker = 'None', MarkerSize = 1.5):
        if not Color: raise ValueError("Missing EZColor 'Color' property.")   #Color is a EZColor object
        
        self.Projected3D = True
        X,Y,Z = XYZ
        #ax = self.plt.axes(projection = '3d')
        self.ax
        plt.plot(X,Y,Z, color = Color.RGB)
        #self.ax = ax

    def Scatter(self, XY, Label = "", Color = None, Marker = 'None', MarkerSize = 1.5):
        if not Color: raise ValueError("Missing EZColor 'Color' property.")   #Color is a EZColor object
        X,Y = XY
        plt.scatter(X,Y, label = Label, color = Color.RGB, marker = Marker, s = MarkerSize)

    def SetTicks(self,Axis, Start, End, Mid):
        #Add comas at thousands (ONLY GREATER THAN 10 000!)
        if Axis == 'X':
            self.XTick['Start'] = Start
            self.XTick['End'] = End
            self.XTick['Mid'] = Mid
        elif Axis == 'Y':
            self.YTick['Start'] = Start
            self.YTick['End'] = End
            self.YTick['Mid'] = Mid
        else: raise ValueError("Axis is either 'X' or' Y'!")

    def ActuateTicks(self,which = 'X'):
        Start,End,Mid = None, None, None
        if which == 'X':
            Start,End,Mid = self.XTick['Start'], self.XTick['End'], self.XTick['Mid']
        else:
            Start,End,Mid = self.YTick['Start'], self.YTick['End'], self.YTick['Mid']

        xticks = np.arange(Start, End, Mid)
        xlist = []
        for x in xticks:
            mystr = str(x)
            if x >= 10000:
                mystr = "{:,}".format(x)
            xlist.append(mystr)

        #This formatting style puts a , after 10 000's place. This is the 'default' style. 
        
        if which == 'X' and self.XNot == 'default':
            self.ax.set_xticklabels(xlist)
            plt.xticks(np.arange(Start, End, Mid))

        elif which == 'Y' and self.YNot == 'default':
            self.ax.set_yticklabels(xlist)
            plt.yticks(np.arange(Start, End, Mid))

    def SetLim(self,Left = None, Right = None, Top = None, Bottom = None, Front = None, Back = None):
        self.Left = Left
        self.Right = Right + Right/900
        self.Top = Top
        self.Bottom = Bottom
        self.Front = Front
        self.Back = Back
        self.SetLimits = True

    def ActuateLims(self,):
        if self.SetLimits == True:
            plt.xlim(left = self.Left, right = self.Right)
            plt.ylim(top = self.Top, bottom = self.Bottom)

    def AxLabels(self,X = "X",Y= "Y", Z = None, Title = None):
        self.XLabel = X
        self.YLabel = Y
        self.ZLabel = Z
        self.Title = Title

    def Finalize3D(self,):
        #print("3d finalizing..")
        self.ax.set_xticks(minor=True)
        self.ax.set_yticks(minor=True)
        self.ax.set_zticks(minor=True)
        self.ax.set_xticks([0,5,10,15,20,25,30,35])
        self.ax.set_yticks([0,5,10,15,20,25,30,35])
        self.ax.set_zticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

        self.ax.axes.set_xlim3d(left=self.Left, right=self.Right-self.Right/900) 
        self.ax.axes.set_ylim3d(bottom=self.Bottom, top=self.Top-self.Top/900) 
        self.ax.axes.set_zlim3d(bottom=self.Back, top=self.Front-self.Top/900)
        
        self.ax.axes.set_xlabel(self.XLabel)
        self.ax.axes.set_ylabel(self.YLabel)
        self.ax.set_zlabel(self.ZLabel)

        
        pass

    def Finalize(self,):
        print("finalizing"); 
        #if self.Projected3D == True: self.Finalize3D(); return 

        if self.XLabel and self.YLabel:
            plt.xlabel(self.XLabel); plt.ylabel(self.YLabel)

        if self.Title != None:
            plt.title(self.Title)

        if self.XTick['Start']:
            self.ActuateTicks('X')
        if self.YTick['Start']:
            self.ActuateTicks('Y')
        self.ActuateLims()

    def Show(self, SaveName = "ISP_EqvFrz", Format = "pdf"):
        #plt.savefig("Figs\{file}".format(file = SaveName + "." + Format))
        plt.show()