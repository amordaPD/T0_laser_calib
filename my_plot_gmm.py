import sys

import itertools

import numpy as np
from scipy import linalg
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd 
from pandas.plotting import scatter_matrix

from sklearn import mixture

color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                              'darkorange','darkgreen','red','brown','violet'])


def plot_results(X, Y_, means, covariances, index, title):
    splot = plt.subplot(2, 1, 1 + index)
    for i, (mean, covar, color) in enumerate(zip(
            means, covariances, color_iter)):
        v, w = linalg.eigh(covar)
        v = 2. * np.sqrt(2.) * np.sqrt(v)
        u = w[0] / linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color,label=(' %i category' % i))
        """
        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180. * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)
        """
    plt.xlim(30., 300)
    #plt.ylim(57, 59)
    plt.ylim(6.5, 8)
    plt.xlabel('amplitude')
    plt.ylabel('time [ns]')
    plt.legend(loc='best')
    #plt.xticks(())
    #plt.yticks(())
    plt.title(title)


def plot_1D_results(X, Y_,var,index_one_ph,title):
    if var==0 :
        range_var=(0,300)
        variable="amplitude (ADC counts)"
        for i in range(0,10):     
            if not np.any(Y_ == i):
                continue
            plt.hist(X[Y_ == i, 0],label=(' %i category' % i), range=range_var, bins=200,alpha=0.3)
    if var==1 :
        range_var=(6.5, 8)#(57,59)
        variable="time [ns]"
        n, bins, patches = plt.hist(X[Y_ == index_one_ph, 1],label=(' %i category' % index_one_ph), range=range_var, bins=200,alpha=0.3,normed=1)
        (mu, sigma) = norm.fit(X[Y_ == index_one_ph, var])
        # add a 'best fit' line
        y = mlab.normpdf( bins, mu, sigma)
        l = plt.plot(bins, y, 'r--', linewidth=2)
        #plot
        plt.xlabel('time')
        plt.ylabel('Probability')
        plt.title(r'$\mathrm{Histogram\ of\ time:}\ \mu_T=%.3f,\ \sigma_T=%.3f$' %(mu, sigma))
        plt.grid(True)

    plt.xlabel(variable)
    plt.ylabel('AU')
    plt.legend(loc='best')
    #plt.title(title)



#run_number = str(sys.argv[1])
    
# Number of samples per componentRead input data
print("reading input data")
listtemp = []
listtime = []
listamps = []
#with open('times_vs_amps.txt', 'r') as f:
#with open('D:\Padova\TOP\Time_calibration_codes\dati\stabilita-thr--20-T50-F2000-1.txt', 'r') as f:
with open('D:\Padova\TOP\Time_calibration_codes\dati\long_run_T50.txt', 'r') as f:
    content = f.readlines()
    content = content[3:-2]
    for x in content:
        row = x.split("*")
        if(float(row[3])<52.0 and float(row[3])>51.4):
            listtime.append(float(row[3]))
            listamps.append(float(row[4]))
            listtemp.append(float(row[5]))
print("input data read")
#print(listtime)
#print("producing scatter matrix")
#inputs = {'amplitude': listamps, ' time ' :  listtime }
#df_inputs = pd.DataFrame(data=inputs)
#scatter_matrix(df_inputs, alpha=0.2, figsize=(6, 6), diagonal='kde')
#plt.show()



print("initializing input matrix")
X=np.column_stack((listamps,listtime))
print(X)
print("input matrix initialized")
'''
print("--------------")
XY=np.column_stack((listamps,listtime,listtemp))
print(XY)
'''
print("--------------")
XX=np.column_stack((X,listtemp))
print(XX)

'''
questa definizione in due steps (prima X) e poi (XX) dell'insieme di dati è resa necessaria dal fatto che 
l'algoritmo utilizza tutte le variabili dell'insieme multivariato (X per internderci)
che gli do come input. Pertanto se io voglio tenere delle variabili spettatrici, 
devo dapprima definire l'insieme multivariato degli input dell'algoritmo, e quindi appendere
a questo la variabile spettatrice.

'''




'''
#Performing BIC
print("Performing model selection BIC")
lowest_bic = np.infty
bic = []
n_components_range = range(3, 8)
cv_types = ['spherical', 'tied', 'diag', 'full']
for cv_type in cv_types:
    for n_components in n_components_range:
        # Fit a Gaussian mixture with EM
        gmm = mixture.GaussianMixture(n_components=n_components,
                                      covariance_type=cv_type)
        gmm.fit(X)
        bic.append(gmm.bic(X))
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            best_gmm = gmm

bic = np.array(bic)
color_iter = itertools.cycle(['navy', 'turquoise', 'cornflowerblue',
                              'darkorange','green','red'])
clf = best_gmm
bars = []

# Plot the BIC scores
spl = plt.subplot(2, 1, 1)
for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
    xpos = np.array(n_components_range) + .2 * (i - 2)
    bars.append(plt.bar(xpos, bic[i * len(n_components_range):
                                  (i + 1) * len(n_components_range)],
                        width=.2, color=color))
plt.xticks(n_components_range)
plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
plt.title('BIC score per model')
xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
    .2 * np.floor(bic.argmin() / len(n_components_range))
plt.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
spl.set_xlabel('Number of components')
spl.legend([b[0] for b in bars], cv_types)

# Plot the winner
splot = plt.subplot(2, 1, 2)
Y_ = clf.predict(X)
for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_,
                                           color_iter)):
    v, w = linalg.eigh(cov)
    if not np.any(Y_ == i):
        continue
    plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)
    
    # Plot an ellipse to show the Gaussian component
    angle = np.arctan2(w[0][1], w[0][0])
    angle = 180. * angle / np.pi  # convert to degrees
    v = 2. * np.sqrt(2.) * np.sqrt(v)
    ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color)
    ell.set_clip_box(splot.bbox)
    ell.set_alpha(.5)
    splot.add_artist(ell)
    
plt.xlim(30., 200)
plt.ylim(57, 59)
plt.title('Selected GMM: full model, 2 components')
plt.subplots_adjust(hspace=.35, bottom=.02)
plt.show()


'''






# Fit a Gaussian mixture with EM using five components
ncomps=9

print("doing the fit with simple Gaussian Mixture")
gmm = mixture.GaussianMixture(n_components=ncomps, covariance_type='full').fit(X)
print("measured means")
print(gmm.means_)
print("measured covariances")
print(gmm.covariances_)





print("#######################")


# Fit a Dirichlet process Gaussian mixture using five components
print("doing the fit with Bayes Gaussian Mixture")
dpgmm = mixture.BayesianGaussianMixture(n_components=ncomps,covariance_type='full',init_params='kmeans').fit(X)
print("measured means")
print(dpgmm.means_)
print("measured covariances")
print(dpgmm.covariances_)

index_one_ph=-9
min_amp = 10000
for j in range(0,ncomps):
    if(dpgmm.means_[j,0])<min_amp:
        min_amp=dpgmm.means_[j,0]
        index_one_ph=j





if(sys.argv[1]==True):
    print("now plotting results")
    plot_results(X, gmm.predict(X), gmm.means_, gmm.covariances_, 0,'Gaussian Mixture')
    plot_results(X, dpgmm.predict(X), dpgmm.means_, dpgmm.covariances_, 1,
             'Bayesian Gaussian Mixture with a Dirichlet process prior')
    plt.show()

if(sys.argv[1]==True):
    plot_1D_results(X, dpgmm.predict(X),0,index_one_ph,'Bayesian Gaussian Mixture with a Dirichlet process prior')
    plt.show()

if(sys.argv[1]==True):
    plot_1D_results(X, dpgmm.predict(X),1,index_one_ph,'Bayesian Gaussian Mixture with a Dirichlet process prior')
    plt.show()

'''
output = open("file.txt", "w")
for iterator in range(0,ncomps):
    amps_cat=X[dpgmm.predict(X) == iterator, 0]
    times_cat=X[dpgmm.predict(X) == iterator, 1]
    for it in X[dpgmm.predict(X) == iterator, 0]:
        output.write(str(it)+' '+str(iterator)+ '\n')
output.close()
'''
'''
output = open("dati\long_run_T50_out.txt", "w")
for iterator in range(0,ncomps):
    for it in X[dpgmm.predict(X) == iterator]:
        output.write(str(it[0])+' '+str(it[1])+' '+str(iterator)+ '\n')
output.close()
'''
output = open("dati\long_run_T50_out.txt", "w")
for iterator in range(0,ncomps):
    for it in XX[dpgmm.predict(X) == iterator]:
        output.write(str(it[0])+' '+str(it[1])+' '+str(it[2])+' '+str(iterator)+ '\n')
output.close()





#times->SetScanField(0);times->Scan("times:amp"); > run_lungo.txt
