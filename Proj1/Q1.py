import numpy as np
from EqFuncs import *
import PlotAssist, EZColors
import matplotlib.pyplot as plt


L1 = 5000; # pipe length, in meter
L2 = 15000 # pipe length, in  meter

D1 = 8 * 0.0254; #pipe dia, meter from in
D2 = 4 * 0.0254; #pipe dia, meter from in
A1 = np.pi * (1/4) * (D1)**2
A2 = np.pi * (1/4) * (D2)**2

eps = 4 * 10 ** (-6); #epsilon roughness, in meters
g0 = 9.81; # gravity, in m/s2
etol = 0.01; #convergeance criteria

dZ = np.linspace(66,88,100) #height, in meters
#dZ = np.array([66])
Trange = np.linspace(10,40,4)+273.15; #temperature, in kelvin
#Trange = np.array([10+273.15])


Hplot = PlotAssist.HigsPlot()
Hplot.AxLabels(X = "Elevation Difference (m)", Y = "Discharge Rate (m $^3$/s)")
Hplot.SetLim(Left = 65, Right = 90, Top = .25, Bottom = 0.1)


Tcounter = 0 #Plot hue incrementer
for T in Trange: #Loop over temperature contour lines
    Tcounter +=1

    #obtain properties at correct temperature
    mu = getVisc(T); #visocisty, in Pa.s
    rho = getDensity(T-273.15); #density, in kg/m3

    
    V2_List = [] #this gathers v's for each height z.
    for Z in dZ:

        #Guess values for 1st iter
        print("Running at", Z)
        fdarcy_1 = compFdarcy_guess(eps,D1); print("f_Darcy1 = ",fdarcy_1)
        fdarcy_2 = compFdarcy_guess(eps,D2); print("f_Darcy2 = ",fdarcy_2)
        Ksc = findKsc()

        V_2 = 0 #storage value 
    
        darcyDiff1 = 10000; darcyDiff2 = 10000 #convergeance criterion past value

        counter = 0
        while True:
            counter+=1
            print("Iterations: ",counter)
            if darcyDiff1 < etol and darcyDiff2 < etol: print("+"*20); break

            #Get head losses normalized by V^2
            Hentrance = 0.5 * (1/(2*g0))
            Hfric1 = fdarcy_1 * (L1/D1) * (1/(2*g0))
            Hfric2 = fdarcy_2 * (L2/D2) * (1/(2*g0))
            Hsc = Ksc * (1/(2*g0))
            Hdisch = (1/(2*g0))

            #the formula is dZ = V^2[A]
            A = Hentrance * (0.25**2) + Hfric1 * (0.25 **2) + Hsc + Hfric2 + Hdisch
            V_2 = np.sqrt(Z/A)
            V_1 = V_2 * A2/A1
            Re_1 = compRe(rho, mu, D1, V_1)
            Re_2 = compRe(rho,mu,D2,V_2) 
            
            print("Obtained V1", V_1)
            print("Obtained V2", V_2)
            print("Calculated Re1", Re_1)
            print("Calculated Re1", Re_2)

            Newfdarcy_1 = compFdarcy_Explicit(eps, D1, Re_1); 
            Newfdarcy_2 = compFdarcy_Explicit(eps, D2, Re_2); 
            

            darcyDiff1 = np.abs(Newfdarcy_1 - fdarcy_1) / Newfdarcy_1
            darcyDiff2 = np.abs(Newfdarcy_2 - fdarcy_2) / Newfdarcy_2
            
            print("darcy_diff", darcyDiff1, darcyDiff2)

            fdarcy_1 = Newfdarcy_1; fdarcy_2 = Newfdarcy_2; Ksc = findKsc(V_2); 
            print("f_Darcy1 = ",fdarcy_1)
            print("f_Darcy2 = ",fdarcy_2)
            print("New ksc", Ksc)
            print("-"*20)
        
        #add to list
        V2_List.append(V_2)
    
    #print("Here is V2", V2_List)
    Q2 = A2 * np.array(V2_List) / A1
    #print("Here is Q2", Q2)
    print("Here is", V2_List)
    Clr = EZColors.CustomColors(colorLabel= 'red')
    Clr.HueShift(Percent=.045*(5-Tcounter))
    Hplot.Plot((dZ, Q2), Color = Clr)



Hplot.plt.legend(["T = 10 $^o C$","T = 20 $^o C$","T = 30 $^o C$","T = 40 $^o C$"], frameon = False)
Hplot.Finalize()
Hplot.Show()


    