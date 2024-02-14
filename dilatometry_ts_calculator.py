import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
from scipy.ndimage import gaussian_filter1d

def mass2mole(mass):
    # define molecular weights of elements
    molweight = {'Fe':55.85,'C':12.01,'Si':28.09,'Mn':54.94,'Ni':58.69,'Cr':52,'Mo':95.94,'W':183.85,'Co':58.93,'V':50.94,'Nb':92.91,'Cu':63.55,'Al':26.98,'Ti':47.88,'O':16,'N':14.01,'B':10.81,'P':30.97,'S':32.06,'As':74.92}
    
    # check that all elements are considered, if not included set comp of element = 0
    if 'Fe' in mass.keys():
        pass
    else:
        mass['Fe'] = 100 - sum(mass.values())
    for e in molweight:
        if e in mass.keys():
            pass
        else:
            mass[e] = 0
            
    # convert weight percent to mole fraction
    a = []
    for elm in mass.keys():
        mole = mass[elm]/molweight[elm]
        a.append(mole)
    tot_moles = sum(a)
    moles = {}
    n = 0
    for mol in a:
        moles[list(mass.keys())[n]] = mol/tot_moles
        n += 1
    return moles

def perc2frac(wperc):
    # convert weight percent to weight fraction
    wegtf = {}
    for e in wperc:
        wegtf[e] = wperc[e]/100
    return wegtf

def a_f0(molef):
    # check that all elements are considered, if not included set comp of element = 0
    elements = ['C','Si','Mn','Ni','Cr','Mo','V']
    for e in elements:
        if e in molef.keys():
            pass
        else:
            molef[e] = 0
            
    # calculate ferrite lattice parameter in nm at 298K
    aFe = 0.28664 
    return aFe + ( ((aFe-(0.0279*molef['C']))**2)*(aFe+(0.2496*molef['C'])) - (aFe**3) )/(3*(aFe**2)) - (0.003*molef['Si']) + (0.006*molef['Mn']) + (0.007*molef['Ni']) + (0.031*molef['Mo']) + (0.005*molef['Cr']) + (0.0096*molef['V'])

def a_a0(wegtf):
    # check that all elements are considered, if not included set comp of element = 0
    elements = ['C','Mn','Ni','Cr','N','Al','Co','Cu','Mo','Nb','Ti','V','W']
    for e in elements:
        if e in wegtf.keys():
            pass
        else:
            wegtf[e] = 0
            
    # calculate austenite lattice parameter in nm at 298K
    return 0.3573 + (0.33*wegtf['C']) + (0.0095*wegtf['Mn']) - (0.002*wegtf['Ni']) + (0.006*wegtf['Cr']) + (0.22*wegtf['N']) + (0.056*wegtf['Al']) - (0.004*wegtf['Co']) + (0.015*wegtf['Cu']) + (0.031*wegtf['Mo']) + (0.051*wegtf['Nb']) + (0.039*wegtf['Ti']) + (0.018*wegtf['V']) + (0.018*wegtf['W'])

def latt_param(comp,p,T):
    # check if inputs are valid
    if p not in ['f','a','m']:
        return print('Invalid microstructure (p) defined. Valid inputs include: f - ferrite, a - austenite, m - martensite.')
    
    # convert comp (in wt.%) to mole fraction (molef) and weight fraction (wegtf)
    
    # austenite weight fraction (wegtf_a)
    wegtf_a = perc2frac(comp)

    # ferrite mole fraction (molef_f) assuming 0.03 wt.% C in solution
    wperc_f = copy.deepcopy(comp)
    wperc_f['C'] = 0.03
    molef_f = mass2mole(wperc_f)

    # martensite mole fraction (molef_m) assuming supersaturated w/ C
    wperc_m = copy.deepcopy(comp)
    molef_m = mass2mole(wperc_m)
    
    # define coefficients of thermal expansion (K-1) for ferrite, f, and austenite, a
    alpha = {'f':1.59*10**(-5),'a':2.03*10**(-5),'m':1.59*10**(-5)}
    
    # define initial lattice parameters @ RT (298K) 
    a0 = {'f':a_f0(molef_f),'a':a_a0(wegtf_a),'m':a_f0(molef_m)}
    
    # calculate lattice parameter at T (in Celcius)
    return ((alpha[p]*(T-25))+1)*a0[p]

def offset_Calc(comp,p,T,X):
    
    # define number of atoms in bcc and fcc unit cells
    N = {'f':2,'a':4,'m':2}
    
    # calculate the volume, V, of the new system (i.e., X amount of microstructure in austenite)
    V = ((1-X)*(1/N['a'])*(latt_param(comp,'a',T)**3)) + (X*(1/N[p])*(latt_param(comp,p,T)**3))
    
    # calculate the volume of the austenite, V0
    V0 = (1/N['a'])*(latt_param(comp,'a',T)**3)
    
    # calculate the strain difference of adding X fraction of new microstructure, p, into austenite at temperature, T
    return ((V**(1/3))/(V0**(1/3)))-1

def dilatometry_curve_plotter(file,L0,N):
    # read data in file
    data = pd.read_csv(file, delimiter = r"\s+", header = None, skiprows=4, engine='python',encoding= 'unicode_escape')
    
    # isolate first 2 columns (which should be 'Temperature' & 'Change in Length')
    for n in range(len(data.columns)):
        if n == 1 or n == 2:
            pass
        else:
            del data[n]
    data.columns = ['T', 'x']
    
    # calculate strain using initial length (L0)
    data['e'] = (data['x']*10**(-6))/L0
    # isolate cooling curve from data
    dataset = data[N:]
    
    # plot cooling curve as 'Temperature' vs. 'Strain'
    plt.figure(figsize=(15,7))
    plt.plot(dataset['T'],dataset['e'],color='mediumseagreen',linewidth=2,label='Dilatometry Curve')
    plt.xlabel('Temperature (\N{DEGREE SIGN}C)', fontsize=16)
    plt.ylabel('Strain (mm/mm)', fontsize=16)
    plt.tick_params(axis='x', labelsize=13)
    plt.tick_params(axis='y', labelsize=13)
    plt.legend(loc='upper left',fontsize=13)
    return

def offset_method(file,L0,Nc,Ti,dT,comp,p,X):
    # read data in file
    data = pd.read_csv(file, delimiter = r"\s+", header = None, skiprows=4, engine='python',encoding= 'unicode_escape')
    
    # isolate first 2 columns (which should be 'Temperature' & 'Change in Length')
    for n in range(len(data.columns)):
        if n == 1 or n == 2:
            pass
        else:
            del data[n]
    data.columns = ['T', 'x']
    # calculate strain using initial length (L0)
    data['e'] = (data['x']*10**(-6))/L0
    
    # isolate cooling curve from data
    dataset = data[Nc:]
    
    # isolate temperature and strain data as Python lists
    temp = list(dataset['T'])
    strain = list(dataset['e'])

    # define temperature ranges for offset lines using Ti and dT
    temp_range = [Ti,Ti+dT,Ti+(2*dT)]
    # create empty list for interception points to be appended to
    interceptions = []
    
    # create figure space to plot into
    fig, axs = plt.subplots(1, 3, figsize = (16,5))
    
    # set count start (i = 0)
    i = 0
    # loop through temperature ranges and calculate interception points, appending to pre-defined list
    for T1 in temp_range:
        # determine index values (n) for temperature range (allowing direct relation to strain list)
        n0, n1 = None, None
        j1, j2 = 0, 0
        while n1 == None:
            try:
                n1 = list(np.around(temp,0)).index(round(T1+j1,0))
            except ValueError:
                j1 += 1
        
        while n0 == None:
            try:
                n0 = list(np.around(temp,0)).index(round(T1+50+j2,0))
            except ValueError:
                j2 += 1
        
        # define equation for straight line (y = m.x + c) by calculating m and c using temperature range
        m = (strain[n0]-strain[n1])/(temp[n0]-temp[n1])
        c = strain[n1] - m*temp[n1]
        
        # create empty lists for the straight line data and offset line data
        lin_strain, off_strain = [], []
        
        # loop through temperature values to fill lists with data
        for T in temp:
            lin_strain.append(m*T + c)
            
        # determine start of transformation (i.e, when 0 < X < 0.01)
        idx1 = np.argwhere(np.gradient(np.sign(np.array(lin_strain) - np.array(strain)))).flatten() 
        To = temp[idx1[-1]]
        
        # calculate offset at To
        offset = offset_Calc(comp,p,To,X)
        
        # loop through temperature values to fill lists with data
        for T in temp:
            off_strain.append(m*T + c + offset)
        # determine intersection points between cooling curve ('strain') and offset line ('off_strain')
        idx2 = np.argwhere(np.gradient(np.sign(np.array(off_strain) - np.array(strain)))).flatten() 
        intercept = temp[idx2[-1]]
        
        # append temperature of interception to list
        interceptions.append(intercept)
        
        # plot analysis on graph for user visuals
        axs[i].plot(temp,strain,color='mediumseagreen',zorder=0,label='Dilatometry Curve')
        axs[i].plot(temp,lin_strain,color='k',linestyle=':',zorder=1,label='Gradient Line')
        axs[i].plot(temp,off_strain,color='steelblue',linestyle='-',zorder=2,label='Offset Line')
        axs[i].scatter(intercept,m*intercept + c + offset,color='tomato',marker='x',label='$T_\mathrm{s}$')
        axs[i].text(intercept+25, m*intercept + c + offset, str(round(intercept,1))+'\N{DEGREE SIGN}C', fontsize = 12, bbox={'facecolor': 'white'}, color = 'red')
        axs[i].set_xlabel('Temperature (\N{DEGREE SIGN}C)',fontsize=14)
        axs[0].set_ylabel('Strain (mm/mm)',fontsize=14)
        axs[i].legend(loc='best',fontsize=11.5)
        axs[i].set_title('Dilatometry Curve \n Gradient range: '+str(T1)+' - '+str(T1+50)+'\N{DEGREE SIGN}C',fontsize=15)
        minX, maxX = T1-200, T1+50
        axs[i].set_xlim(minX,maxX)
        axs[i].set_ylim(m*(minX)+c,m*(maxX)+c)
        # increase count
        i += 1
        
    # calculate and return the average interception temperature (i.e., Ts) and standard deviation as list -> [Ts,Std]
    return [round(np.mean(interceptions),1), round(np.std(interceptions),1)]

def second_deriv_method(file,L0,N,s,dT):

    # read data
    data = pd.read_csv(file, delimiter = r"\s+", header = None, skiprows=4, usecols = range(1,3), engine='python', encoding= 'unicode_escape')
    data.columns = ['T', 'x']
    data['e'] = (data['x']*10**(-6))/L0

    # isolate even spacing between temperature values (i.e., a 1C spacing) - this is for calculating the gradient later
    x, y = np.array(data['T'][N:]), np.array(data['e'][N:])
    yi = []
    n0 = 0
    xi = np.linspace(int(max(x)),int(min(x)),(int(max(x))-int(min(x)))+1)
    for T in xi:
        idx = closest_index(x,T)
        yi.append(y[idx])

    # smooth data using a Gaussian filter using sigma, s
    smooth = gaussian_filter1d(yi, s)

    # compute second derivative of smoothed data
    smooth_d2 = np.gradient(np.gradient(smooth))

    fig, axs = plt.subplots(1, 2, figsize = (11,5))

    # find position of start and finish temperatures in dT
    n2, n1 = list(xi).index(dT[0]), list(xi).index(dT[-1])

    # crop data around start and finish temperatures, dT
    xi_crop, yi_crop, smooth_crop, smooth_d2_crop = xi[n1:n2], yi[n1:n2], smooth[n1:n2], smooth_d2[n1:n2]

    # calculate the temperature/s when 2nd derivative = 0
    y0 = np.linspace(0,0,len(xi[n1:n2]))
    idx = np.argwhere(np.gradient(np.sign(y0 - smooth_d2_crop))).flatten() 
    intercepts = []

    # average temperatures when the 2nd derivative = 0, and find standard deviation
    for n in range(len(idx)):
        if smooth_d2_crop[idx[n]-1] < smooth_d2_crop[idx[n]]:
            intercepts.append(xi_crop[idx[n]])
    interception_point = np.mean(intercepts)
    std = np.std(intercepts)

    # plot dilatometry cooling curve, the Gaussian (smoothed) curve and the interception point
    n3, n0 = list(xi).index(int(interception_point-50)), list(xi).index(int(interception_point+50))
    axs[0].scatter(xi[n0:n3], yi[n0:n3], marker='.',color='mediumseagreen',s=7.5,zorder=1, label='Raw Data')
    axs[0].plot(xi[n0:n3], smooth[n0:n3],color='k',zorder=0, label='Gaussian Filtered (σ='+str(s)+')',linewidth=1)
    axs[0].errorbar(interception_point,yi_crop[list(xi_crop).index(int(interception_point))], marker = 'x', color = 'r', xerr = std, ecolor = 'k', capsize=5, linestyle="None",zorder=2,label='$T_\mathrm{s}$')    
    axs[0].legend(loc='best',fontsize=11.5)
    axs[0].set_ylabel('Strain, ε (mm/mm)',fontsize=14)
    axs[0].set_xlabel('Temperature (\N{DEGREE SIGN}C)',fontsize=14)
    axs[0].set_title('Dilatometry Curve',fontsize=15)

    # plot second derivative data and intercept point found using the Gaussian curve
    axs[1].plot(xi_crop, smooth_d2_crop,label='$\mathrm{d}^2ε/\mathrm{d}T^2$')
    axs[1].errorbar(interception_point,0, marker = 'x', color = 'r', xerr = std, ecolor = 'k', capsize=5, linestyle="None", zorder=2,label='$T_\mathrm{s}$')
    axs[1].text(interception_point-5, (min(smooth_d2_crop))*(0.35), str(round(interception_point,1))+'±'+str(round(std,1))+'\N{DEGREE SIGN}C', fontsize = 13, bbox={'facecolor': 'white'}, color = 'red')
    axs[1].set_ylabel('$\mathrm{d}^2ε/\mathrm{d}T^2$',fontsize=14)
    axs[1].set_xlabel('Temperature (\N{DEGREE SIGN}C)',fontsize=14)
    axs[1].set_title('Second Derivative Curve',fontsize=15)
    axs[1].legend(loc='upper right',fontsize=11.5)
    
    # return Ts and standard deviation
    return [round(interception_point,1), round(std,1)]

def closest_index(lst, K):
    lst = np.asarray(lst)
    idx = (np.abs(lst - K)).argmin()
    return idx
