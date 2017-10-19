import LFPy
import numpy as np
import matplotlib.pylab as pl
import sys
import os
import random as random

#setting a random seed to get the same result all the time
np.random.seed(1988)
#what this program need fro running
#where is and will be the data
f2 = open('cellname', 'r')
celln = [line.strip() for line in f2]
cellname = celln[0]
f2.close()
os.chdir(cellname)
#morphology
#sigma
#electrode coordinates
#morphology
f3 = open('morphology.txt', 'r')
morph = [line.strip() for line in f3]
morpho = morph[0]
f3.close()
#active channel description
f4 = open('active.txt', 'r')
actw = [line.strip() for line in f4]
activewhere = actw[0]
f4.close()

#electrode coordinates
felec = open('elcoord_x_y_z', 'r')
elec = [line.strip() for line in felec]
felec.close()
elc=np.hstack((elec))
elcor=elc.reshape(3,-1)
#tissue properties
f1 = open('elprop', 'r')
elp = [line.strip() for line in f1]
sigma =float(elp[1])
f1.close()
######################x
#synaptic inputs

def stationary_poisson(nsyn,lambd,tstart,tstop):
    ''' Generates nsyn stationary possion processes with rate lambda between tstart and tstop'''
    interval_s = (tstop-tstart)*.001
    spiketimes = []
    for i in xrange(nsyn):
        spikecount = np.random.poisson(interval_s*lambd)
        spikevec = np.empty(spikecount)
        if spikecount==0:
            spiketimes.append(spikevec)
        else:
            spikevec = tstart + (tstop-tstart)*np.random.random(spikecount)
            spiketimes.append(np.sort(spikevec)) #sort them too!

    return spiketimes

###############
cell_parameters = {         
   
	'morphology' : morpho ,    
	'Ra': 123,
        'tstartms' : 0.,                 # start time of simulation, recorders start at t=0
        'tstopms' : 850.,                   # stop simulation at 200 ms. 
	'passive' : True,
    	'v_init' : -65,             # initial crossmembrane potential
    	'e_pas' : -65,              # reversal potential passive mechs
	'nsegs_method' :  'fixed_length',
	'max_nsegs_length':10, 
#	'lambda_f' : 1000,           # segments are isopotential at this frequency
    'custom_code'  : [activewhere], # will run this file
}

#electrode coordinates
x=elcor[0,]
x= x.astype('Float64')
y=elcor[1,]
y= y.astype('Float64')
z=elcor[2,]
z= z.astype('Float64')

#	y = pl.zeros(X.size)

#define parameters for extracellular recording electrode, using optional method
electrodeParameters = {
    'sigma' : sigma,              # extracellular conductivity
    'x' : x,        # x,y,z-coordinates of contact points
    'y' : y,
    'z' : z,
#     'method' : 'som_as_point',  #treat soma segment as sphere source
#     'method' : 'pointsource'
     'method' : 'linesource'
}
   





#pointprocess= {
#        'idx' : 1,
#        'record_current' : True,
#        'pptype' : 'IClamp',
#        'amp' : 0.05,
#        #'amp' : 0.2,
#        'dur' : 30,
 #       'delay' : 15
#    }



##create extracellular electrode object
electrode = LFPy.RecExtElectrode(**electrodeParameters)

simulationParameters = {
	'electrode' : electrode, 
	'rec_imem' : True,  # Record Membrane currents during simulation
	'rec_isyn' : True,  # Record synaptic currents
}




#Initialize cell instance, using the LFPy.Cell class
cell = LFPy.Cell(**cell_parameters)

#rotating the cell
#rotation = {'x' : np.pi/2, 'y' : 0, 'z' : 0}
#cell.set_rotation(**rotation)

#set the position of midpoint in soma to Origo (not needed, this is the default)
#Why I am doing this? #That might modifies the plots and the setups!!!!!!!!!!!!!!44
cell.set_pos(xpos = LFPy.cell.neuron.h.x3d(0) , ypos = LFPy.cell.neuron.h.y3d(0) , zpos = LFPy.cell.neuron.h.z3d(0))
#cell.set_pos(xpos = xpontok[1], ypos = ypontok[1], zpos = zpontok[1])



#Synaptic inputs


	    # Define synapse parameters
synapse_parameters = {
      'idx' : cell.get_closest_idx(x=100., y=0., z=500.), #100 0.00 500
      'e' : 0.,                   # reversal potential
      'syntype' : 'ExpSyn',       # synapse type
      #'tau' : 10.,                # syn. time constant
	'tau' : 2.,
#      'weight' : .001,            # syn. weight
      'weight' : .04,            # syn. weight 
      'record_current' : True,
}

synapse_parameters2 = {
  'idx' : cell.get_closest_idx(x=-100., y=0., z=400.), #100 0.00 500
      'e' : 0.,                   # reversal potential
      'syntype' : 'ExpSyn',       # synapse type
      #'tau' : 10.,                # syn. time constant
	'tau' : 2.,
#      'weight' : .001,            # syn. weight
      'weight' : .03,            # syn. weight 
      'record_current' : True,
}


synapse_parameters3 = {
      'idx' : cell.get_closest_idx(x=40., y=0., z=200.), #100 0.00 500
      'e' : 0.,                   # reversal potential
      'syntype' : 'ExpSyn',       # synapse type
      #'tau' : 10.,                # syn. time constant
	'tau' : 2.,
#      'weight' : .001,            # syn. weight
      'weight' : .04,            # syn. weight 
      'record_current' : True,
}



############################################x
# Define synapse parameters
synapse_parameters_random = {
    'idx' : 0, # to be set later
    'e' : 0.,                   # reversal potential
    'syntype' : 'ExpSyn',       # synapse type
    'tau' : 2.,                 # syn. time constant
    'weight' : .01,            # syn. weight
    'record_current' : True,
}

#synaptic spike times
n_pre_syn = 1000
pre_syn_sptimes = stationary_poisson(nsyn=n_pre_syn, lambd=2, tstart=0, tstop=400)#70)

#assign spike times to different units
n_synapses = 1000


# Create synapse and set time of synaptic input
pre_syn_pick = np.random.permutation(np.arange(n_pre_syn))[0:n_synapses]

for i_syn in xrange(n_synapses):
    syn_idx = int(cell.get_rand_idx_area_norm())
    synapse_parameters_random.update({'idx' : syn_idx})
    synapse = LFPy.Synapse(cell, **synapse_parameters_random)
    synapse.set_spike_times(pre_syn_sptimes[pre_syn_pick[i_syn]])
##############################################################


# Create synapse and set time of synaptic input
#synapse = LFPy.Synapse(cell, **synapse_parameters)
#synapse2 = LFPy.Synapse(cell, **synapse_parameters2)
#synapse3 = LFPy.Synapse(cell, **synapse_parameters3)
#insert_synapses(synapse_parameters_2, **insert_synapses_2)
#synapse.set_spike_times(np.array([5.,25., 44.]))
#synapse2.set_spike_times(np.array([15.,20., 40.]))
#synapse3.set_spike_times(np.array([10.,33.]))


TimesStim= np.arange(850)

for istim in xrange(850):
    pointprocess= {
        'idx' : 0,
#        'record_current' : True,
        'pptype' : 'IClamp',
        'amp' :  np.array(3.6*np.sin(2.*3.141*6.5*TimesStim/1000.))[istim], #3.6
#        #'amp' : 0.2,
        'dur' : 1.,
        'delay' : istim,
    }
    stimulus = LFPy.StimIntElectrode(cell, **pointprocess)





#stimulus = LFPy.StimIntElectrode(cell, **pointprocess)
	

#perform NEURON simulation, results saved as attributes in the cell instance
cell.simulate(**simulationParameters)



#np.savetxt( 'Istim',stimulus.i)
np.savetxt( 'membcurr',cell.imem)
np.savetxt( 'myLFP', electrode.LFP)

np.savetxt( 'somav.txt', cell.somav)

coords = np.hstack(
	(cell.xmid, cell.ymid, cell.zmid) 
)

np.savetxt( 'coordsmid_x_y_z',coords)


#coordinates of the segment's beginning
coordsstart = np.hstack(
	(cell.xstart, cell.ystart, cell.zstart) 
)

np.savetxt( 'coordsstart_x_y_z',coordsstart)

#coordinates of the segment's end
coordsend = np.hstack(
	(cell.xend, cell.yend, cell.zend) 
)

np.savetxt( 'coordsend_x_y_z',coordsend)

#sdiameter of the segments
segdiam = np.hstack(
	(cell.diam) 
)

np.savetxt( 'segdiam_x_y_z',segdiam)

##########x
#elec = np.hstack(
#	(electrode.x, electrode.y, electrode.z) 
#)

#np.savetxt(outname,' + 'elcoord_x_y_z',elec)

#length of segments
np.savetxt( 'seglength',cell.length)
#time in the simulation
np.savetxt( 'time',cell.tvec)

#lets write to file the simulation locations

np.savetxt( 'synapse_locations',pre_syn_pick)


#elprop=np.hstack((d,electrode.sigma))
#np.savetxt( 'elprop',electrode.sigma)
# Plotting of simulation results:


################################x

#LFPy.cell.neuron.h...
#h = LFPy.cell.neuron.h
#hossz=len(cell.allsecnames)
#f4 = open(outname,' + '/segcoords/branchnum', 'w')
#f4.write(str(hossz))
#f4.close()
#b=0
#for x in cell.allseclist:
#	b=b+1
#	xc=list()
#	yc=list()
#	zc=list()
#	for i in range(int(h.n3d())):
 #               #print h.x3d(i)
#		xc.append(h.x3d(i))
#		yc.append(h.y3d(i))
#		zc.append(h.z3d(i))	
#	np.savetxt(outname,' + '/segcoords/segcord'+str(b),np.hstack((xc,yc,zc)))

