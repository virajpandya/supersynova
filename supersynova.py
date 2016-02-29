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

import os, sys, shutil, subprocess,
import numpy as np
import matplotlib.pyplot as plt
from pyraf import iraf
from iraf import onedspec
from time import strftime

def getIonName(labelNumber):
    '''
    This dictionary maps the SYN++ label for a ion to that ion's actual name (for plotting).
    '''
    ionLabels = ({'100':'H I','200':'He I','600':'C I','601':'C II','602':'C III','800':'O I','801':'O II','1100':'Na I','1200':'Mg I',
                 '1201':'Mg II','1400':'Si I','1401':'Si II','1402':'Si III','1600':'S I',
             '1601':'S II','2000':'Ca I','2001':'Ca II','2002':'Ca III','2101':'Sc II','2201':'Ti II',
             '2401':'Cr II','2600':'Fe I','2601':'Fe II','2602':'Fe III','2701':'Co II','2702':'Co III','2801':'Ni II'})
    return(ionLabels.get(labelNumber,'???'))

def run(yaml,spectrum=''):
    '''
    This runs SYN++ with the input YAML file and plots the composite
    and single-ion spectra.
    '''

    # copy input YAML file and remove any inactive ion columns for a composite run
    print '-----> Removing inactive ion columns...'
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
    f2 = open('supersyn.yaml','w')
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

    # now run syn++ with all active ions turned on to create composite synthetic spectrum
    plt.ioff()
    fig,ax = plt.subplots(1,figsize=(8,8))
    print('-----> Running SYN++ for composite synthetic spectrum...')

    if spectrum != '': # comparison observed spectrum
        r = subprocess.call('syn++ --wl-from='+spectrum+' supersyn.yaml > supersynfull.dat',shell=True)
        wobs,fobs,errobs = np.genfromtxt(spectrum,unpack=True)
        wsyn,fsyn,errsyn = np.genfromtxt('supersynfull.dat',unpack=True)
        ax.plot(wobs,fobs,c='0.7')
        ax.plot(wsyn,fsyn,c='orange')
    else:
        r = subprocess.call('syn++ supersyn.yaml > supersynfull.dat',shell=True)
        wsyn,fsyn,errsyn = np.genfromtxt('supersynfull.dat',unpack=True)
        ax.plot(wsyn,fsyn,c='red',linewidth=1.5)

    # now create single ion spectra (one active ion turned on at a time)
    if os.path.exists('singleIon/') == False: os.makedirs('singleIon/')
    shutil.copyfile('supersyn.yaml','supersyn_original.yaml')

    for col,active in enumerate(activeCols):
        currentCols = activeCols[:]
        # Set all ion columns to inactive except for current ion
        for c,a in enumerate(currentCols):
            if c == col:
                currentCols[c] = 'Yes'
            else:
                currentCols[c] = 'No'
        f = open('supersyn.yaml','r') # Replace active list with new single-ion-activated list
        ftmp = open('supersyntmp.yaml','w')
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

        print('-----> Running SYN++ to create '+ionName+'-only synthetic spectrum...')
        if spectrum != '':
            r = subprocess.call('syn++ --wl-from='+spectrum+' supersyntmp.yaml > singleIon/supersyn'+ionName.replace(' ','')+'.dat',shell=True)
            wsyn,fsyn,errsyn = np.genfromtxt('singleIon/supersyn'+ionName.replace(' ','')+'.dat',unpack=True)
            ax.plot(wsyn,fsyn-col,c='blue')
            ax.text(wsyn.max()+100,np.median(fsyn)-col,ionName,color='blue',fontsize=12)
        else:
            r = subprocess.call('syn++ supersyntmp.yaml > singleIon/supersyn'+ionName.replace(' ','')+'.dat',shell=True)
            wsyn,fsyn,errsyn = np.genfromtxt('singleIon/supersyn'+ionName.replace(' ','')+'.dat',unpack=True)
            ax.plot(wsyn,fsyn+col,c='blue')
            ax.text(wsyn.max()+100,np.median(fsyn)-col,ionName,color='blue',fontsize=12)

    # set axis labels, ticks, etc.
    ax.set_xlabel('Wavelength [$\AA$]',fontsize=16)
    ax.set_ylabel(r'Scaled Flux $\plus$ Offset}',fontsize=16)
    #print('Setting logarithmic scale...')
    #ax.set_xscale('log')
    #ax.set_xticks([3000,4000,5000,6000,7000,8000,9000,10000])
    #ax.set_xticklabels( ['3000','4000','5000','6000','7000','8000','9000','10000'] )
    #ax.set_autoscale_on(False)
    #ax.set_xlim([wsyn.min()-100,wsyn.max()+650])
    plt.savefig('supersynova.png',bbox_inches='tight',dpi=300)
    plt.close('all')

if __name__ == '__main__':
    # parse arguments: first = yaml file, and second = observed spectrum ASCII file
    if len(sys.argv) != 1 or len(sys.argv) != 2:
        sys.exit('Need arguments: first=yaml file, second=observed spectrum ASCII file')
    if sys.argv[1][-5:] != 'yaml':
        sys.exit('First argument must end with .yaml')
    try:
        sys.argv[2]
    except:
        sys.argv[2] = ''

    # call the main run function
    run(yaml=sys.argv[1],spectrum=sys.argv[2])









##
