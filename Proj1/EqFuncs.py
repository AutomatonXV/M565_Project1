import numpy as np

#print(np.log(10), np.log10(10)), log() is ln, log10 is log base 10

#get viscosity 
def getVisc(T): #for water
    #T is temperature in Kelvin
    #from Baliga LH1, table 1-4:
    # a = 29.76; b = -5.24
    # return np.exp(a + b * np.log(T)) # IN CENTIPOISE!!! DIVIDE BY 1000 for Pa.s

    #T is temperature in Kelvin,  baliga lh1 table 1-12
    z = 273/T; mu0 = 1.788 * 10**(-3)
    return mu0 * np.exp(-1.704 - 5.306 * z + 7.003 * z**2)

def getDensity(T):
    #T is temperature in celcius, 
    #from baliga LH1, table 1-12, for WATER
    return 1000 - 0.0178 * np.abs(T - 4)**(1.7)

def compRe(rho,mu,D,V):
    return (rho * V * D) / mu

def compFdarcy_guess(e,d):
    A = (e/d)/(3.7)
    B = -2 * np.log10(A)
    return B**(-2)

def compFdarcy_Explicit(e,d,Re):
    A = (e/d)/(3.7)
    B = (5.02/Re) * np.log10(A + (14.49/Re))

    return (-2 * np.log10(A - B)) ** (-2)
    
def findKsc(V = 6):
    #V = 6 is a guess, returns Ksc = 0.33 (guess value)
    #THIS IS ONLY FOR D2/D1 = 2
    Vlist = np.array([.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    Klist = np.array([.38, .38, .37, .36, .35, .34, .33, .33, .32, .31, .3])

    for vi in range(0,Vlist.size-1):
        currV = Vlist[vi]
        nextV = Vlist[vi+1]
        if currV == V: return Klist[vi]
        if nextV == V: return Klist[vi+1]
        if currV < V and V < nextV:            
            currK = Klist[vi]
            nextK = Klist[vi+1]
            #interpolate:
            return currK + (V-currV) * (nextK - currK)/(nextV - currV)
        
    return 0.3 #all case failed, V is very big
            

        


