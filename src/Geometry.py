# Random point placement in the spherical shell : Coordinates(CoreDia, NPDia, Vff)
# Compute scattering amplitude of point scaterrer: RayleighAlpha(R, WL, n, K)
import numpy as np

def Coordinates(CoreDia, NPDia, Vff):

    LigLen = 0.9800
    ShellThickness = NPDia * 3
    R = (NPDia / 2 + LigLen) 

    SmallSphereRadius = R
    ShellOuterRad = CoreDia /2  + ShellThickness
    ShellInnerRadius = CoreDia/2

    #Vs = (4/3)* np.pi * (R ** 3)
    #Floor returns rounding to nearest highest value of the floating number
    ND = np.floor((ShellOuterRad ** 3 - ShellInnerRadius**3) / (R**3) * Vff)
    N = int(ND) # Number of nanoparticles in integre
    #This part of the code calcuates random distributions of the coordinates of the
    #Allocating memory for coordinates
    R_N = np.full((N,3), float('nan')) 
    # Allocating memory for specific nanoparticles coordinates
    Xcoord = np.full((N,1), float('nan'))
    Ycoord = np.full((N,1), float('nan'))
    Zcoord = np.full((N,1), float('nan'))
    for i in range(N):
        t = 2 * np.pi * np.random.rand() #Elevation angle
        p = np.arccos(2* np.random.rand()-1) # Azimuthal Angle
        dr = ShellInnerRadius +(np.random.rand())**(1/3) * (ShellOuterRad - ShellInnerRadius) # Distance along z
        Xcoord[i] = dr * np.cos(t) * np.sin(p)
        Ycoord[i] = dr * np.sin(t) * np.sin(p)
        Zcoord[i] = dr * np.cos(p)
        I = []
        # finding vacancies and checking for overlap 
        NoLap = 0 # Seting overlap flag
        #m = 0; 
        #Compute distances between the naoparticles 
        if i > 1 :
            
            dist = np.sqrt((Xcoord - Xcoord[i])**2 + (Ycoord - Ycoord[i])**2 + (Zcoord - Zcoord[i])**2)
           #looking for the overlap 
            idx = [idx for idx in range(i)]
            I = np.where(dist[idx] < SmallSphereRadius)
            I = np.asanyarray(I)
            if I.size == 0:
                NoLap = 0 #No overlap
            else:
                NoLap = 1 # Overlap Present
        if NoLap == 1:
            while NoLap == 1:
                
                t = 2 * np.pi * np.random.rand() #Elevation angle
                p = np.arccos(2* np.random.rand()-1) # Azimuthal Angle
                dr = ShellInnerRadius +(np.random.rand())**(1/3) * (ShellOuterRad - ShellInnerRadius )     
                Xcoord[i] = dr * np.cos(t) * np.sin(p)
                Ycoord[i] = dr * np.sin(t) * np.sin(p)
                Zcoord[i] = dr * np.cos(p)
                dist = np.sqrt((Xcoord - Xcoord[i])**2 + (Ycoord - Ycoord[i])**2 + (Zcoord - Zcoord[i])**2)
                idx = [idx for idx in range(i)]
                I = np.where(dist[idx] < SmallSphereRadius)
                I = np.asanyarray(I)
                if I.size == 0:
                    NoLap = 0
                else:
                    NoLap = 1

        X_coords = Xcoord
        Y_coords = Ycoord
        Z_coords = Zcoord 


    R_N[:,0:1] = X_coords
    R_N[:,1:2] = Y_coords
    R_N[:,2:3] = Z_coords
    return R_N

def RayleighAlpha(R, WL, n, K):
      
      nb = 1.0 # Refractive index of the Background medium
      N = np.zeros(1, dtype=complex)
      m = np.zeros(1,dtype=complex)
      K0 = np.zeros(1,dtype=complex)
      sigma_s = np.zeros(1,dtype=complex)
      sigma_abs = np.zeros(1,dtype=complex)
      sigma_t = np.zeros(1,dtype=complex)
      alpha_imag = np.zeros(1,dtype=complex)
      alpha_r = np.zeros(1,dtype=complex)
      alphaN = np.zeros(1,dtype=complex)
      # R Radius of the nanoparticle

      N = np.complex128(n,K)
      m = N / nb
      K0 = (2 * np.pi / WL)
      sigma_s = (((2*np.pi**5/3)*((2*R)**6)*(np.absolute((m**2-1)/(m**2+2)))**2) /(WL**4))
      sigma_abs = ((((8*np.pi**2)*(1/WL)*R**3*np.imag((m**2-1)/(m**2+2)))))
      sigma_t = sigma_s+ sigma_abs
      alpha_imag = (K0)*sigma_t
      alpha_r = np.absolute(np.sqrt((4*np.pi*sigma_s) - K0*alpha_imag**2))
      alphaN = (4*np.pi)*(alpha_r + 1j*alpha_imag)

      return alphaN


