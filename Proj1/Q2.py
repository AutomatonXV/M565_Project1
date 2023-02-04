import numpy as np
from EqFuncs import *
import PlotAssist, EZColors

Z_A = 120; # meters
Z_B = 100; # meters
Z_C = 80;  # meters

L1 = 1000; # meters
L2 = 4000; # meters
L3 = 2000; # meters

D1 = 12 * 0.0254; # meters from in
D2 = 24 * 0.0254; # meters from in
D3 = 18 * 0.0254; # meters from in

A1 = np.pi * (1/4) * D1 ** 2; # m2 
A2 = np.pi * (1/4) * D2 ** 2; # m2 
A3 = np.pi * (1/4) * D3 ** 2; # m2 

Trange = np.linspace(10,40,4)+273.15 #Temperature, in Kelvin
#Trange = np.array([10+273.15]) #Temperature, in Kelvin

etol = 0.01; # tolerance. If etol  | Q1 - Q2 - Q3 | < etol, life is good
eps = 4 * 10 ** (-6); #epsilon roughness, in meters
g0 = 9.81; # gravity, in m/s2

Q1List = []
Q2List = []
Q3List = []
ZpList = []

for T in Trange:
    print("Running for Temperature:", T-273.15)
    #obtain properties at correct temperature
    mu = getVisc(T); #visocisty, in Pa.s
    rho = getDensity(T-273.15); #density, in kg/m3

    dZ = .33; # meters, decrement value. Ensure it can never reach a whole number

    #GUESS VALUES
    Zp = 102 ; # meters, guess value. Iterate this
    fd1 = compFdarcy_guess(eps,D1)
    fd2 = compFdarcy_guess(eps,D2)
    fd3 = compFdarcy_guess(eps,D3)
    Qerr = 10000; #dummy storage value
    
    Count = 0
    while True: 
        Count+=1
        
        
        #Need to iterate darcy friction
        def DarcyIterator(fd0, L,D,dZ):
            #print("Invoked for ", L, D, dZ)
            #dZ is Za - Zp for example
            #L is length, D is dia, fd is guess darcy
            #iterate darcy and return obtained V
            etol_darcy = 0.01
            fd = fd0 #set guess
            
            
            while True:
                
                #compute V from guess fd
                Hf = fd * (L/D) * (1/ (2*g0))  #normalied by V^2
                #print(Hf, dZ)
                V = np.sqrt(dZ/Hf)
                #compute Re (from V) and check
                Re = compRe(rho,mu,D,V); 
                if(Re < 4000): print("CRITICAL FAIL: REYNOLD NUMBER FAIL")
                #use this Re to compute a more accurate f_d
                newfd = compFdarcy_Explicit(eps,D,Re)
                diff = np.abs(newfd - fd) / fd
                if (diff < etol_darcy): return V, newfd, Re
                #if new f_d has not converged relative to last fd, repeat loop.
                fd = newfd
            
        #Flip the signs based on Zp.
        s1 = 1; s2 = 1; s3 = 1 #sign multipliers
        if Zp < Z_A: s1 = -1
        if Zp < Z_B: s2 = -1
        if Zp < Z_C: s3 = -1

        V1, fd1, Re1 = DarcyIterator(fd1, L1, D1, s1*(Zp - Z_A))
        V2, fd2, Re2 = DarcyIterator(fd2, L2, D2, s2*(Zp - Z_B))
        V3, fd3, Re3 = DarcyIterator(fd3, L3, D3, s3*(Zp - Z_C))
        
        Q1 = A1 * V1; Q2 = A2 * V2; Q3 = A3 * V3
        
        
        
        Q_err = s1*Q1 + s2*Q2 + s3*Q3
        #print("Qerror is", Q_err)
        #print("Discharge Q1 is", s1 * Q1)
        #print("Discharge Q2 is", s2 * Q2)
        #print("Discharge Q3 is", s3 * Q3)

        #Exit condition reached.
        if (np.abs(Q_err) < etol or Count > 10000): Q1List.append(s1*Q1); Q2List.append(s2*Q2); Q3List.append(s3*Q3); ZpList.append(Zp); break
            
        #Not reached? If Qerr is positive, then bring zp down, else up
        if Q_err > 0: 
            Zp -= dZ * 0.9
        else:
            Zp += dZ * 0.9

print("+"*20)
Q1List = np.array(Q1List) ; Q2List = np.array(Q2List) ; Q3List = np.array(Q3List)
QSum = Q1List + Q2List + Q3List
print("QList1   ", Q1List)
print("QList2   ", Q2List)
print("QList3   ", Q3List)
print("Net Sum", QSum)
print("ZList", ZpList)



Hplot = PlotAssist.HigsPlot()
Hplot.AxLabels(X = "Temperature ($^o C$)", Y = "Discharge Rate (m $^3$/s)")
Hplot.SetLim(Left = 10, Right = 40, Top = .5, Bottom = 0.1)
Clr1 = EZColors.CustomColors(colorLabel= 'red')
Clr2 = EZColors.CustomColors(colorLabel= 'orange')
Clr3 = EZColors.CustomColors(colorLabel= 'violet')
Hplot.Plot((Trange-273.15, np.abs(Q1List)), Color=Clr1)
Hplot.Plot((Trange-273.15, np.abs(Q2List)), Color=Clr2)
Hplot.Plot((Trange-273.15, np.abs(Q3List)), Color=Clr3)
Hplot.plt.legend(["(-) $Q_{1}$","(-) $Q_{2}$","(+) $Q_{3}$"], frameon = False)
Hplot.Finalize()
Hplot.Show()


Clr = EZColors.CustomColors(colorLabel='black')
Hplot2 = PlotAssist.HigsPlot()
Hplot2.AxLabels(X = "Temperature ($^o C$)", Y = "Net Sum Discharge Rate (m $^3$/s)")
Hplot2.SetLim(Left = 10, Right = 40, Top = .05, Bottom = -0.05)
Hplot2.Plot((Trange-273.15, QSum),Color = Clr)
Hplot2.Finalize()
Hplot2.Show()



