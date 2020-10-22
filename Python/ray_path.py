import numpy as np
import matplotlib.pyplot as plt


I    = 2482      # Rows of file ModelS
J    = 6         # Columns of file ModelS
PI   = np.pi
RSun = 69.634e9  # cm


rad2, c, rho, p, Gamma, T = np.loadtxt("modelS.txt",unpack=True)
log_rho = np.log(rho)

rad = rad2*RSun
der_log_rho = np.zeros(I)
H = np.zeros(I)
wc = np.zeros(I)
der_H = np.zeros(I)

for i in range(0,I,1):
    if (i==0):
        der_log_rho[i] = (log_rho[i+1]-log_rho[i])/(rad[i+1]-rad[i])
    elif (i==I-1):
        der_log_rho[i] = (log_rho[i]-log_rho[i-1])/(rad[i]-rad[i-1])
    else:
        der_log_rho[i] = (log_rho[i+1]-log_rho[i-1])/(2*(rad[i+1]-rad[i]))
    H[i] = -1.0*(der_log_rho[i])**(-1)
    wc[i] = c[i]/(2*H[i])


for i in range(0,I,1):
    if (i==0):
        der_H[i] = (H[i+1]-H[i])/(rad[i+1]-rad[i])
    elif (i==I-1):
        der_H[i] = (H[i]-H[i-1])/(rad[i]-rad[i-1])
    else:
        der_H[i] = (H[i+1]-H[i-1])/(2*(rad[i+1]-rad[i]))
    #wc[i] = (c[i]/(2*H[i]))*(1-2*der_H[i])**0.5



# ************************************************************
# ************************************************************


def Lamb(l):
    L = (l*(l+1))**0.5
    return L



def tau(c, L, rad, w):
    res = 2.0*RSun/(c*( 1.0- (L*c/(rad*w))**2.0 )**0.5)
    return res





def Delta(c, L, rad, w, wc):
    #res = 2*(c*RSun/rad)/(pow(pow(w/L,2) - pow(c/rad,2),0.5));
    res = 2.0*RSun/(( (w*rad/(L*c))**2.0 - 1.0)**0.5)
    #res = 2*(c*RSun/rad)/( ((w/L)**2 - (wc/L)**2 - (c/rad)**2)**0.5);
    return res



# ************************************************************
# ************************************************************

w = 2*PI*5e-3


for l in range(200,4001,200):
    f = open("function_l_py/T_D_" +str(l)+".txt","w")

    for w in np.arange(2*PI*1e-3,2*PI*11e-3,2*PI*1e-4):
        
        tau_area = 0.0
        result_t = 0.0
        Del_area = 0.0
        result_D = 0.0

        for r in range(I-2,1-1,-1):

            res_t = tau(c[r], Lamb(l), rad[r], w)
            res_D = Delta(c[r], Lamb(l), rad[r], w, wc[r])

            if np.isnan(res_t) and np.isnan(res_D):
                tau_area += 0.0
                Del_area += 0.0 
            else: 
                b = (rad[r-1]-rad[r])/RSun
                if np.isnan(result_t):
                    result_t = 0
                    result_D = 0
                h_t = res_t + result_t
                h_D = res_D + result_D
                area_t = b*h_t/2
                area_D = b*h_D/2
                tau_area += area_t
                Del_area += area_D

            result_t = res_t
            result_D = res_D

        f.write(str(tau_area/60)+ " " + str(Del_area/1e8)+"\n")
    f.close()








