'''
Viraj Pandya, Saurabh Jha, Curtis McCully, Yssavo Camacho, Brandon Patel

supersynova is a plotting script for use with SYN++.
 
Inputs:
- SYN++ YAML file (only one setup allowed)
- Optional: standardized/normalized rest-wavelength ASCII observed/comparison
        spectrum with 3 columns (in this order): wavelength, flux, flux_error
Output:
- Single image with full and single-ion synthetic spectra, and optional input spectrum.

December 20, 2013: Created. -VP
'''

import os, sys, shutil, subprocess, numpy as np, matplotlib.pyplot as plt
from pyraf import iraf
from iraf import onedspec
from time import strftime

'''Dictionary which maps ion SYN++ label numbers to actual ion names'''
def getIonName(labelNumber):
    ionLabels = ({'100':'H I','200':'He I','600':'C I','601':'C II','602':'C III','800':'O I','801':'O II','1100':'Na I','1200':'Mg I',
                 '1201':'Mg II','1400':'Si I','1401':'Si II','1402':'Si III','1600':'S I',
             '1601':'S II','2000':'Ca I','2001':'Ca II','2002':'Ca III','2101':'Sc II','2201':'Ti II',
             '2401':'Cr II','2600':'Fe I','2601':'Fe II','2602':'Fe III','2701':'Co II','2702':'Co III','2801':'Ni II'})
    return(ionLabels.get(labelNumber,'???'))

'''
The user's choice of x-axis units will define any horizontal offsets and the x-axis label.
'0' is Angstroms (default), and '1' is log(Angstroms). 
The elements are tuples: (xLeftOffset,xRightOffset,ionLabelOffset,xLabel)
'''
def getProps(wtype):
    waveProps = ({'linear':(-100,600,100,r'$\mathrm{Rest\,Wavelength} (\AA)$'),
                 'log':(-100,600,100,r'$\mathrm{Rest\,Wavelength} (\AA)$')})
    return(waveProps[wtype])
    
    

'''MAIN PLOTTING FUNCTION'''
def run(yaml,spectrum='',offset=1.0,units='linear',continuum=True,objectname='',plusoffset=0.0):
    timestamp = strftime("%Y%m%d-%H-%M-%S")
    os.mkdir(timestamp)
    
    plusoffset = float(plusoffset)
    if offset == 1.0 and continuum == True:
        offset = 1.5
    print("Copying original YAML file and removing any setup columns for inactive ions...")
    f = open(yaml,'r')
    activeLine = ''
    for line in f: # Find the line for the 'active' setup parameter.
        if 'active' in line:
            activeLine = line[:]
        else:
            continue
    f.close()
    if activeLine == '':
        print('ERROR: YAML file structure incorrect.')
        return
    activeStatus = (activeLine[activeLine.find('[')+1:activeLine.find(']')]).split(',')
    inactiveCols = [] # Add column numbers of inactive ions to this list.
    for c,s in enumerate(activeStatus):
        if (s.strip()).lower() == 'no':
            inactiveCols.append(c)
    f = open(yaml,'r')
    f2 = open(timestamp+'/supersyn.yaml','w') # Copy old YAML file but keep only active ion columns.
    for line in f:
        if 'log_tau_min' in line: 
            f2.write(line)
        elif ('ions' in line or 'active' in line or 'log_tau' in line or
              'v_min' in line or 'v_max' in line or 'aux' in line or 'temp' in line):
            prefix = line[0:line.find('[')+1]
            Cols2 = (line[line.find('[')+1:line.find(']')]).split(',')
            Cols = [b for a,b in enumerate(Cols2) if a not in inactiveCols]
            newRow = prefix + ','.join(Cols) + ']\n'
            f2.write(newRow)
            if 'active' in line: # Will need this to coordinate single-ion SYN++ runs later.
                activeCols = Cols[:]
        else: 
            f2.write(line)
    f.close()
    f2.close()
    '''
    Run supersyn.yaml file through SYN++ once with all ions active.
    '''
    plt.ioff() # Prevents figure from opening interactively when script is run/imported in ipython --pylab
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    print("Running SYN++ to create full synthetic spectrum...")
    if spectrum != '':
        r = subprocess.call("syn++ --wl-from="+spectrum+" "+timestamp+"/supersyn.yaml > "+timestamp+"/supersynfull_"+timestamp+".dat",shell=True)
        wobs,fobs,errobs = np.genfromtxt(spectrum,unpack=True)
        wsyn,fsyn,errsyn = np.genfromtxt(timestamp+'/supersynfull_'+timestamp+'.dat',unpack=True)
        ax.set_xlim([wsyn.min()+getProps(units)[0],wsyn.max()+getProps(units)[1]])
        ax.plot(wobs,fobs,c='black',alpha=1.0,linewidth=2)
        ax.plot(wsyn,fsyn,c='red',linewidth=1.5)
        ax.text(wsyn.max()+getProps(units)[2],0.35,'All',color='red',fontsize=12)
        if continuum == True:
            ax.text(6860.0,0.9*offset,'$\oplus$',color='black',fontsize=12) # Observed telluric absorption B-band
            ax.text(7590.0,0.9*offset,'$\oplus$',color='black',fontsize=12) # Observed telluric absorption A-band   
        else:
            ax.text(6860.0,1.1*offset,'$\oplus$',color='black',fontsize=12) # Observed telluric absorption B-band
            ax.text(7590.0,1.1*offset,'$\oplus$',color='black',fontsize=12) # Observed telluric absorption A-band       
    else:
        r = subprocess.call("syn++ "+timestamp+"/supersyn.yaml > "+timestamp+"/supersynfull_"+timestamp+".dat",shell=True)
        wsyn,fsyn,errsyn = np.genfromtxt(timestamp+'/supersynfull_'+timestamp+'.dat',unpack=True)
        ax.set_xlim([wsyn.min()+getProps(units)[0],wsyn.max()+getProps(units)[1]])
        ax.plot(wsyn,fsyn,c='red',linewidth=1.5)
        ax.text(wsyn.max()+getProps(units)[2],0.35,'All',color='red',fontsize=12)
    '''
    Now, turn every ion off except one, going from left to right, to create the
    single-ion synthetic spectra. Do all of this in a loop which reads and
    plots each synthetic spectrum immediately after it is created. 
    '''
    shutil.copyfile(timestamp+'/supersyn.yaml',timestamp+'/supersynfull_'+timestamp+'.yaml')
    voffset = offset
    for col,active in enumerate(activeCols):
        currentCols = activeCols[:] # Set all ions (columns) to inactive except column col
        for c,a in enumerate(currentCols):
            if c == col:
                currentCols[c] = 'Yes'
            else:
                currentCols[c] = 'No'
        f = open(timestamp+'/supersyn.yaml','r') # Replace active list with new single-ion-activated list 
        ftmp = open(timestamp+'/supersyntmp.yaml','w')
        for line in f:
            if 'active' in line:
                prefix = line[0:line.find('[')+1]
                newRow = prefix + ','.join(currentCols) + ']\n'
                ftmp.write(newRow)
            elif 'ions' in line: # Use ionLabels dictionary mapping to find ion name
                ionCols = (line[line.find('[')+1:line.find(']')]).split(',')
                ionName = getIonName(ionCols[col].strip())
                ftmp.write(line)
            else:
                ftmp.write(line)
        f.close()
        ftmp.close()
        os.remove(timestamp+'/supersyn.yaml')
        shutil.move(timestamp+'/supersyntmp.yaml',timestamp+'/supersyn.yaml')
        # Run SYN++ with this single-ion-activated YAML file and add to the existing plot (+ offset).
        # Also save this single-ion spectrum in a separate file for interactive access later.
        print("Running SYN++ to create "+ionName+"-only synthetic spectrum...")
        if spectrum != '':
            r = subprocess.call("syn++ --wl-from="+spectrum+" "+timestamp+"/supersyn.yaml > "+timestamp+"/supersyn"+ionName.replace(" ","")+"_"+timestamp+".dat",shell=True)
            wsyn,fsyn,errsyn = np.genfromtxt(timestamp+'/supersyn'+ionName.replace(" ","")+'_'+timestamp+'.dat',unpack=True)
            ax.plot(wsyn,fsyn+voffset,c='blue',linewidth=1.5)
            if continuum == True:
                ax.text(wsyn.max()+getProps(units)[2],voffset,ionName,color='blue',fontsize=12)
            else:
                ax.text(wsyn.max()+getProps(units)[2],voffset+offset*.75,ionName,color='blue',fontsize=12)
        else:
            r = subprocess.call("syn++ "+timestamp+"/supersyn.yaml > "+timestamp+"/supersyn"+ionName.replace(" ","")+"_"+timestamp+".dat",shell=True)
            wsyn,fsyn,errsyn = np.genfromtxt(timestamp+'/supersyn'+ionName.replace(" ","")+'_'+timestamp+'.dat',unpack=True)
            ax.plot(wsyn,fsyn+voffset,c='blue',linewidth=1.5)
            if continuum == True:
                ax.text(wsyn.max()+getProps(units)[2],voffset,ionName,color='blue',fontsize=12)
            else:
                ax.text(wsyn.max()+getProps(units)[2],voffset+offset*.75,ionName,color='blue',fontsize=12)
        voffset += offset
    '''
    if plusoffset != 0.0: # temporary last-minute quick fix (new run parameter) for 13gr
        ax.set_ylim([0.0,ax.get_ylim()[1]+plusoffset]) 
    elif continuum == True:
        ax.set_ylim([0.0,ax.get_ylim()[1]-offset]) # To remove excessive white space at the bottom and top      
    else:
        ax.set_ylim([0.0,ax.get_ylim()[1]+offset]) # To remove excessive white space at the bottom and top
    '''
    ax.set_xlabel(getProps(units)[3],fontsize=26)
    ax.set_ylabel(r'$\mathrm{Scaled \,Flux \,+\, Offset}$',fontsize=26)
    #ax.grid(True)
    #ax.set_title(r'$\mathrm{Decomposed \,Synthetic \,Spectrum}$',fontsize=24)
    #ax.text(0.60,0.95,objectname,transform=ax.transAxes,verticalalignment='top',fontsize=24)
    ax.set_title(objectname,fontsize=24)
    ax.tick_params(axis='y',labelleft='off')
    ax.tick_params(axis='x',labelsize=18)
    if units == 'log': # log scale
        print("Setting logarithmic scale...")
        ax.set_xscale('log')
        ax.set_xticks([3000,4000,5000,6000,7000,8000,9000,10000])
        ax.set_xticklabels( ['3000','4000','5000','6000','7000','8000','9000','10000'] )
        ax.set_autoscale_on(False)
        ax.set_xlim([wsyn.min()-100,wsyn.max()+650])    
    fig.tight_layout()
    fig_name = timestamp+'/supersynplot_'+strftime("%Y%m%d-%H-%M-%S")+'.png'
    plt.savefig(fig_name)
    print("FINISHED: created "+fig_name)
    plt.clf()
    
    # return name of figure in case
    return fig_name

'''
Supplements the driver below by getting four necessary inputs from the user.
'''
def getInputs():
    while True:
        yamlFile = raw_input("Please enter a SYN++ YAML filename: ")

        if os.path.exists(yamlFile):

            break
        else:
            print "ERROR: file does not exist."
    while True:
        spectrumFile = raw_input("Please enter a comparison spectrum filename (OPTIONAL): ")
        if os.path.exists(spectrumFile) == True or spectrumFile == '' or spectrumFile == ' ':
            break
        else:
            print "ERROR: file does not exist."
    while True:
        offsetInput = raw_input("Please enter a vertical offset for spectra (default=1.0): ")
        if offsetInput == '' or offsetInput == ' ':
            offsetFloat = 1.0
            break
        else:
            try: 
                offsetFloat = float(offsetInput)
                break
            except:
                print "ERROR: must enter a number or leave blank."
    while True:
        unitsInput = raw_input("Please enter 'linear' to use Angstroms or 'log' to use log(Angstroms) (default=Angstroms): ")
        if unitsInput == 'linear' or unitsInput == 'log':
            break
        elif unitsInput == '' or unitsInput == ' ':
            unitsInput = 'linear'
            break
        else:
            print "ERROR: must enter either 0, 1, or leave blank."
    return (yamlFile,spectrumFile,offsetFloat,unitsInput)

'''
DRIVER function.
'''
def driver():
    '''Prompt user for SYN++ YAML file and optional input spectrum (if script is run interactively).'''
    print("Welcome to the supersynova plotting script!")
    print("Here is the menu: ")
    while True:
        print("1. Run supersynova with SYN++ YAML file and optional spectrum as-is.")
        print("2. Run supersynova after automatic continuum subtraction for input spectrum.")
        print("3. Exit.")
        choice = raw_input("Please enter a number from the menu: ")
        if choice == '1':
            inputs = getInputs()
            run(yaml=inputs[0],spectrum=inputs[1],offset=inputs[2],units=inputs[3])
        elif choice == '2': 
            print("Under construction")
            '''
            inputs = getInputs()
            # Automatically run iraf.onedspec.continuum with order-2 spline3 division
            onedspec.continuum.unlearn()
            onedspec.continuum(
            run(yaml=inputs[0],spectrum=inputs[1],offset=inputs[2],units=inputs[3],continuum=False)
            '''
        elif choice == '3':
            sys.exit("Thanks for using supersynova!")
        else:
            print("ERROR: you must enter a number from the menu.")
            
    
    
    
'''User has the option of running this script interactively, so call driver() in that case.'''
if __name__ == "__main__":
    driver()
    
