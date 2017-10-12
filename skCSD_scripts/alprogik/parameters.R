#Parameters


#where is the data? LFP,electrode properties and the files needed in case we know the original membrane currents
data.location<-'/media/BA0ED4600ED416EB/agy/kCSD/progik/bs_futtat/branching_simple2014/morpho1_el144_valami2'
#where to save the output
output.location<-'/media/BA0ED4600ED416EB/agy/kCSD/progik/bs_futtat/branching_simple2014/output'
#where is the morphology
morpho.location<-'/media/BA0ED4600ED416EB/agy/kCSD/progik/bs_futtat/branching_simple2014/simulation/morphology/morpho1.swc'
#name of the simulation
today <- Sys.Date()
simulation.name<-paste('morpho1_',format(today, format="%Y_%m_%d"),sep='')
#Do we know the orginal membrane currents? y-yes, n-no
#in case yes, we need the membrane currents('membcurr'), 'coordsstart_x_y_z','coordsend_x_y_z','seglength'
memb.currents.known<-'y' #'n'

#parameters of basis functions
#widths
basis.width.min<-20
basis.width.max<-200
basis.width.step<-40
#number
basis.number.min<-100
basis.number.max<-400
basis.number.step<-50
