#I hate the way colorsys works

import colorsys

clrLabels = {
    'white'   :   (1,1,1),                'black'   :   (0,0,0),
    'red'     :   (1,0,0),                'green'   :   (0,1,0),
    'blue'    :   (0,0,1),                'cyan'    :   (0,1,1),
    'magenta' :   (1,0,1),                'yellow'  :   (1,1,0),
    'violet'  :   (148/255,0,211/255),    'indigo'  :   (75/255,0,130/255),
    'orange'  :   (1,0.5,0),              'magenta' :   (1, 51/255, 131/255)
}

class CustomColors:
    def __init__(self,colorLabel = None, RGB = None):
        if colorLabel == None and RGB == None:
            raise ValueError("Need to provide either the RGB or the colorlabel")
        
        if RGB == None:
            RGB = clrLabels[colorLabel]
        else:
            RGB = (RGB[0]/255, RGB[1]/255, RGB[2]/255)  #Normalize RGB

        self.RGB = RGB
        self.HSV = colorsys.rgb_to_hsv(self.RGB[0],self.RGB[1],self.RGB[2])
        pass

    def Shade(self,percent):
        H,S,V = self.HSV
        self.HSV = (H,S,percent)
        self.RGB = colorsys.hsv_to_rgb(self.HSV[0],self.HSV[1],self.HSV[2],)
        pass

    def HueShift(self,Val = None, Percent = None):
        if Val == None and Percent == None:
            raise ValueError("Need either Value or percent")
        if Val != None:
            #Normalize Val to 360
            Val = Val/360
        if Percent != None:
            Val = Percent
        H,S,V = self.HSV
        self.HSV = H+Val, S, V
        self.RGB = colorsys.hsv_to_rgb(self.HSV[0],self.HSV[1],self.HSV[2],)
        pass

    def Invert(self):
        self.HueShift(Percent = 50)
        pass
