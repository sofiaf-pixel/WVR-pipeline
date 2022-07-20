

import numpy as np
import scipy.optimize
import pylab as pl
pl.ion()

def sim(x, p):
    a, b, c  = p
    return np.exp(-b * (x-a)**2) + c

def err(p, x, y):
    return sim(x, p) - y


# set up the data
data_x = np.linspace(0, 40, 50)
p1 = [20, .4, 0.5]       # parameters for the first trajectory
p2 = [25, .15, 0.2]       # parameters for the second trajectory, same b
data_y1 = sim(data_x, p1)
data_y2 = sim(data_x, p2)
ndata_y1 = data_y1 + np.random.normal(size=len(data_y1), scale=0.01)
ndata_y2 = data_y2 + np.random.normal(size=len(data_y2), scale=0.01)

pl.figure(1)
pl.clf()
pl.plot(data_x,ndata_y1,'s')
pl.plot(data_x,ndata_y2,'s')

# independent fitting of the two trajectories
print ("Independent fitting")
p_best, ier = scipy.optimize.leastsq(err, p1, args=(data_x, ndata_y1))
print ("Best fit parameter for first trajectory", p_best )

pl.plot(data_x,sim(data_x, p_best))

p_best, ier = scipy.optimize.leastsq(err, p2, args=(data_x, ndata_y2))
print( "Best fit parameter for second trajectory", p_best )
pl.plot(data_x,sim(data_x, p_best))
# global fit

# new err functions which takes a global fit
def err_global(p, x, y1, y2):
    # p is now a_1, b, c_1, a_2, c_2, with b shared between the two
    p1 = p[0], p[1], p[2]
    p2 = p[3], p[1], p[4]
    
    err1 = err(p1, x, y1)
    err2 = err(p2, x, y2)
    return np.concatenate((err1, err2))

p_global = [20, .3, 0.5, 25., 0.2]    # a_1, b, c_1, a_2, c_2
p_best, ier = scipy.optimize.leastsq(err_global, p_global, 
                                    args=(data_x, ndata_y1, ndata_y2))

p_best_1 = p_best[0], p_best[1], p_best[2]
p_best_2 = p_best[3], p_best[1], p_best[4]
print ( "Global fit results")
print ("Best fit parameters for first trajectory:", p_best_1 )
print ("Best fit parameters for second trajectory:", p_best_2 )

#pl.figure(2)
#pl.clf()
#pl.plot(data_x,ndata_y1,'rs')
pl.plot(data_x,sim(data_x,p_best_1),'b')
#pl.plot(data_x,ndata_y2,'bs')
pl.plot(data_x,sim(data_x,p_best_2),'b')
pl.text(-0.5,1.2,'Two Gaussians fit for  width,\n location and $\sigma$', fontsize=16, color='k')
pl.text(-0.5,.95,'The blue curves are constrained to \n have the same width, $\sigma$', fontsize=16, color='b')
pl.savefig('dual_fits.png')
