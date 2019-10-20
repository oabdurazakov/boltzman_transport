#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('pylab', 'inline')


# In[41]:


def fermi(x,T):    
    return 1./(exp(x/T) + 1.0)


# In[3]:


def bose(x,T):
    return 1./(exp(x/T) - 1.0)


# In[4]:


def electron_dos_linear(x):
    r = zeros(len(x), dtype = float)
    cond = logical_and(x >= band_bottom, x <= band_top)
    r[cond] = x[cond]
#    r[cond] = x[cond]**2f
#    r[cond] = 1.0
    return r


# In[5]:


def phonon_dos(x):
    r = zeros(len(x), dtype = float)
    c1 = logical_and(x>=0.0, x <= 0.004)
    c2 = logical_and(x>0.004, x <= 0.018)
    r[c1] = 6.25*1e4*x[c1]**2
    r[c2] = 1.0
    return r


# In[6]:


def phonon_dos_gauss(x):
    # center
    mu = 0.010
    # width
    sigma = 0.002
    return exp(-power(x-mu,2.) / (2 * power(sigma,2)))


# In[7]:


band_bottom = 0.0
band_top = 2.0
phonon_band_width = 0.020
margin = 0.050

tmin = 0.0
tmax = 2000

dt = 0.2
dx = 0.001
dw = dx
Nt = int(round(tmax/dt)) + 1
Nx = int(round((band_top - band_bottom + 2 * margin)/dx))+1
Nw = int(round((phonon_band_width - band_bottom)/dx))+1

# Grids
t = linspace(tmin,tmax,Nt)
x = linspace(band_bottom - margin, band_top + margin,Nx)
w = linspace(band_bottom, phonon_band_width, Nw)
    
print "Electron Band: [",band_bottom,"-", band_top, "] eV"
print "Maximum Time: ", tmax


# In[8]:


# Electron DOS
N = electron_dos_linear(x)
# Phonon DOS
F = phonon_dos(w)


# In[9]:


fig,ax = subplots(1,1)
ax.plot(x,N,'b-',linewidth=3)
#ax.set_xlim(0,0.4)
#ax.set_ylim(0,1)
ax.set_xlabel(r'$x\mathrm{[eV]}$',fontsize=28)
ax.set_ylabel(r'$N(x)$',fontsize=28)
ax.tick_params(axis='both',labelsize=14)
savefig('edos.png',bbox_inches='tight')


# In[10]:


fig,ax = subplots(1,1)
ax.plot(w*1000,F,'r-o',linewidth=3)
#ax.set_ylim(0,1.1)
ax.set_xlim(0,0.02*1000)
ax.set_xlabel(r'$\omega\,\mathrm{[meV]}$',fontsize=28)
ax.set_ylabel(r'$F\,(\omega)$',fontsize=28)
ax.tick_params(axis='both',labelsize=14)
savefig('pdos.png',bbox_inches='tight')


# In[11]:


def Icoll(f,Tph):
    ### Collision integral due to  the eletron-phonon scattering
    
    coll = zeros(len(x), dtype = float)

    for ix, xx in enumerate(x):
        
        if xx >= band_bottom and xx <= band_top:
            fa,fe,ee,ea = 0.0,0.0,0.0,0.0
            
            for iw,ww in enumerate(w):
                if ww!=0.0:
                    
                    fa +=      bose(ww,Tph)  * f[ix-iw] * (1 - occ*f[ix]   ) * N[ix-iw] * F[iw]
                    fe += (1 + bose(ww,Tph)) * f[ix+iw] * (1 - occ*f[ix]   ) * N[ix+iw] * F[iw]
                    ee += (1 + bose(ww,Tph)) * f[ix]    * (1 - occ*f[ix-iw]) * N[ix-iw] * F[iw]
                    ea +=      bose(ww,Tph)  * f[ix]    * (1 - occ*f[ix+iw]) * N[ix+iw] * F[iw]
            coll[ix] = lam * (fa + fe - ea - ee) * dw - G_escape * f[ix]
#             print coll[ix],fa,fe,ea,ee
    return coll


# In[12]:


def evolve(Tel,Tph):
        ### Time evolution of the electron density by the finite-difference method
        
        dens = zeros([len(t), len(x)], dtype = float)
        f = zeros(len(x), dtype = float)
        cond = logical_and(x >= band_bottom, x <= band_top)
        ## Initial electron distribution
        if fermi==True:
            f[cond] = fermi(x[cond],Tel) 
        else:            
            const = 0.1
            excit_height = 0.2
            f[cond] = const * exp(-x[cond]/excit_height)
        print sum(f*N)*dx
        for it,tt in enumerate(t):

                dens[it:,] = f
                f += Icoll(f,Tph)*dt
                if it%(int(round(Nt/10)))==0:
                        print sum(f*N)*dx, " Time: ", tt 
        return dens


# In[14]:


occ = True
fermi =False
print "Occupied Band? ", occ
lam = 70.0
G_escape = 0.000
print "Escape rate from the band:", G_escape
Tel = 0.2
#Tph = 0.004309 # 50 K
Tph = 0.008617 # 100 K
#Tph = 0.01723 # 200 K
#Tph = 0.02585 # 300 K
print "Particle number in the band:"
import time
t0=time.clock()
pop = evolve(Tel,Tph)
t1=time.clock()
print "Time spent for evolution:",t1-t0," s"


# In[28]:


# fig,ax = subplots(1,1,figsize=(8,6))
# i = 90
# print x[i]
# ax.plot(t,pop[:,i],linewidth=1,label='python')
# xx,yy = loadtxt('../../lex_sobota/p04.dat',unpack=True)
# ax.plot(xx,yy,'k--',linewidth=3,label='c++')
# ax.legend(loc=0)
# savefig('trace.png',bbox_inches='tight')


# In[29]:


def check():

    fig,ax=subplots(1,2,figsize=(16,6))
    #ax[0].set_xlim(amin(t),amax(t))
    ax[0].set_xlabel(r'$\mathrm{Delay\,time\,[arb.u]}$',fontsize=24)
    #ax[0].set_ylim(0,0.55)
    ax[0].set_ylabel(r'$f(\mathrm{t})\mathrm{\,[arb.\,u.]}$',fontsize=28)
    #ax[1].set_xlim(amin(x),amax(x))
    ax[1].set_xlim(0,1)

    ax[1].set_ylabel(r'$f(E-E_\mathrm{F})\mathrm{\,[arb.\,u.]}$',fontsize=28)
    ax[1].set_xlabel(r'$E-E_F\,\mathrm{[eV]}$',fontsize=24)
    j=0
    for i in range(Nx/10,Nx/5,4):
        ax[0].plot(t,pop[:,i],linewidth=1)
        j += 1
    j=0
    for i in range(0,Nt,Nt/100):
#         ax[1].plot(x,pop[i,:],linewidth=1)
        ax[1].plot(x,pop[i,:],linewidth=1)

        j +=1
#     ax[1].set_ylim(0.001,2)
    #ax[1].legend(loc=1,fontsize=16)
    #ax[0].legend(loc=1,fontsize=16)
    ax[0].tick_params(axis='both',labelsize=16)
    ax[1].tick_params(axis='both',labelsize=16)
#    ax[1].plot(x,fermi(x-0.068571,T=Tph),'k--',linewidth=3,label=r'$n_\mathrm{F}(T_\mathrm{ph})$')
    ax[1].legend(loc=0,fontsize=24)

#    savefig('evolve_fermi_noescape.png',bbox_inches='tight')


# In[30]:


check()


# In[31]:


import matplotlib.colors as mcolors

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

darkblue = (0,0,139/256.)
darkgreen = (0,128/256.,64/256.)
blue = (128/256.,128/256.,256/256.)
green = (0/256.,128/256., 64/256.)
yellow = (1,1,0)
white = (1,1,1)
sienna4 = (128/256.,64/256.,20/256.)
rvb = make_colormap(        [ white, darkblue, 0.075, darkblue, darkgreen, 0.3, darkgreen, yellow, 0.48, yellow, sienna4, 0.52, sienna4,         white, 0.9, white, white, 1.0, white] )
rvb1 = make_colormap(        [ white, darkblue, 0.2, darkblue, darkgreen, 0.2, darkgreen, yellow, 0.4, yellow, sienna4, 0.52, sienna4,         white, 0.9, white, white, 1.0, white] )

rvb2 = make_colormap([white, blue, 0.15/0.25-0.06, blue, green, 0.15/0.25+0.06, green, green, 0.75, green, white, 1., white])


# In[32]:


fig,ax = subplots(1,1,figsize=(8,6))
im = ax.imshow(pop.T, extent=(0,tmax,band_bottom-margin,band_top+margin), origin='lower',aspect='auto', cmap=rvb)
ax.set_ylim(0,0.1)
ax.set_xlabel(r'$\mathrm{Delay\,time\,[arb.u]}$',fontsize=24)
ax.set_ylabel(r'$E-E_F\,\mathrm{[eV]}$',fontsize=24)
ax.tick_params(axis='both',labelsize=14)
#savefig('falsecolor_fermi_noescape.png',bbox_inches='tight')


# In[42]:


def get_temp(pop,ff=False):   
    
    def boltzmann_fit(x,y,p0):
        
        from scipy.optimize import curve_fit   
        def fitfunc(x,T,a):
            return  a*exp(-x/T)
        
        if p0==None:
            p0 = [0.05,1.0]
        itmin = argmin(abs(x-0))
        itmax = argmin(abs(x-1.0))
        xx = copy(x[itmin:itmax])
        yy = y[itmin:itmax]
        (popt, pcov) = curve_fit(fitfunc,xx,yy,p0=p0)  
        p0 = popt
        
        return popt[0],sqrt(pcov[0,0]),popt[1]

    def fermi_fit(x,y,p0):
        
        from scipy.optimize import curve_fit   
        def fitfunc(x,T,mu,a):
            return  a*ffermi(x-mu,T)
        if p0==None:
            p0 = [0.1,0.1,1.1]
        itmin = argmin(abs(x-0))
        itmax = argmin(abs(x-0.9))
        xx = copy(x[itmin:itmax])
        yy = y[itmin:itmax]
        (popt, pcov) = curve_fit(fitfunc,xx,yy,p0=p0)    
        p0=popt
        return popt[0],sqrt(pcov[0,0]),popt[1]
    
    gamma = []
    err = []
    xl =[]
    p0=None
    for it,tt in enumerate(t):
        if it>0 and it%5==0:
            z = pop[it,:]
            #tt = copy(t[450:550])
            #zz = copy(z[450:550])
            #g,yerr=expo_fit(tt,zz)
            if ff:
                g,yerr,mu = fermi_fit(x,z,p0)
            else:
                g,yerr,mu = boltzmann_fit(x,z,p0)
            gamma.append(g)
            err.append(yerr)
            xl.append(tt)
            if it==2450:
                print mu
   
    gamma = array(gamma)
    err = array(err)
    xl = array(xl)
    return xl,gamma,err  


# In[43]:


x1,y1,z1 = get_temp(pop,True)


# In[44]:


fig,ax = subplots(1,1)
con = 11604 
ax.errorbar(x1,y1*con,yerr=z1,fmt='r-',color = 'r',markeredgewidth=0,linewidth=2)
ax.axhline(y=Tph*con,color='k',linestyle='--',linewidth=2,label=r'$T_\mathrm{ph}$')
ax.set_xlabel(r'$t\,\mathrm{[arb.u]}$',fontsize=28)
ax.set_ylabel(r'$T_\mathrm{el}(t)$',fontsize=28)
ax.tick_params(axis='both',labelsize=14)
ax.legend(loc=0,fontsize=20)
print 0.1*con
#ax.set_xlim(0,9)
ax.set_ylim(0,y1[0]*con+200)
#ax.annotate(r'$T_\mathrm{el}(0)=580\, \mathrm{K}$',xy=(0.6,0.6),xycoords='axes fraction',color='k',fontsize=20) 
#savefig('temp_fermi_noescape.png',bbox_inches='tight')


# In[122]:


def get_rates(pop):   

    def expo_fit(x,y,p0):
        
        from scipy.optimize import curve_fit   
        def fitfunc(x,a,b,c):
            return a*exp(-x*b)+c  
        def convolve_exp_norm(x,alpha, mu, sigma,a,b):
            co = alpha/2.0 * exp( alpha*mu+ alpha*alpha*sigma*sigma/2.0)
            x_erf = (mu + alpha*sigma*sigma - x)/(sqrt(2.0)*sigma)
            r = a*co * exp(-alpha*x) * (1.0 - scipy.special.erf(x_erf))+b
            return r
        
        if p0==None:
            p0 = [y[0],0.002,y[-1]]        
        (popt, pcov) = curve_fit(fitfunc,x,y,p0=p0)    
        p0=popt
        #ax.plot(x,fitfunc(x,popt[0],popt[1],popt[2]),'r--',label='Fit',linewidth=2)      
        return popt[1],sqrt(pcov[1,1]),p0
    gamma = []
    err = []
    xl =[]
    p0 = None
    for ixx,xx in enumerate(x):
        if xx>band_bottom and xx<=band_top and ixx%5==0:
            z = pop[:,ixx]
            z1 = 0.2*(amax(z)-z[-1])+z[-1]
            z2 = 0.002*(amax(z)-z[-1]) + z[-1]
            icen = argmin(abs(z-amax(z1)))
            imin = argmin(abs(z-z1))
            imax = argmin(abs(z-z2))
            tt = copy(t[imin:imax])
            zz = copy(z[imin:imax])
            g,yerr,p0=expo_fit(tt,zz,p0)
            print p0
#            g,yerr=expo_fit(t,z,p0)
            gamma.append(g)
            err.append(yerr)
            xl.append(xx)
   
    gamma = array(gamma)
    err = array(err)
    xl = array(xl)
    return xl,gamma,err  



x1,y1,z1 = get_rates(pop)

fig,ax=subplots(1,1,figsize=(7,6))
#scaling factor to make the ylabel look better
sc = 200
ax.errorbar(x1,y1*sc,yerr=z1,fmt='d-',color = 'r',markersize=8,markeredgewidth=0,linewidth=2,label=r'$100\mathrm{K}$')
ax.set_xlabel(r'$E-E_F\,\mathrm{[eV]}$',fontsize=28)
ax.set_ylabel(r'$\mathrm{Decay\: Rate\,[arb.\,u.]}$',fontsize=28)  
ax.tick_params('both', colors='k',labelsize=16)
ax.legend(loc=0,fontsize=24,frameon=False)
ax.set_ylim(0.0,0.03*sc)
ax.set_xlim(0,0.2)
ax.annotate(r'$\Gamma_\mathrm{escape}=%.1f$'%(G_escape*sc),xy=(0.6,0.1),xycoords='axes fraction',color='k',fontsize=28) 
#savefig('rates_temp_occ.png',bbox_inches='tight')










