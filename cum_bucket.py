import math
rho = 0.0018 #slugs/in^3
a = 0.0270
gamma = 1.2
cstar = 4890 # ft/s
R = 2000 #psi.ft^3.slugs^-1.R^-1
T0 = 4700 # R
g0 = 32.1741 #ft/s^2
n = 0.30


def simulate(num_grains, grain_length, grain_inner_diameter, throat_diameter):
    mass = rho*num_grains*grain_length*math.pi*(grain_inner_diameter/2)^2
    pass


def specificimpulse(impulse, propmass, gravconst):
    return ((impulse) / (propmass * gravconst))

def thrust(massflow, exhaustvelocity, exhaustarea, exhaustpressure, ambientpressure):
    return ((massflow * exhaustvelocity) + ((exhaustpressure - ambientpressure) * exhaustarea))

def machnumber(velocity, speedsound):
    return (velocity / speedsound)

def stagnationtemperature(adiabaticconst, machnumber):
    return (1+((adiabaticconst-1)/2)*(machnumber**2))

def stagnationpressure(adiabaticconst, machnumber):
    return ((1+((adiabaticconst-1)/2)*(machnumber**2))**(adiabaticconst/(adiabaticconst-1)))

def stagnationdensity(adiabaticconst, machnumber):
    return ((1+((adiabaticconst-1)/2)*(machnumber**2))**(1/(adiabaticconst-1)))

def areamachnumberrelation(adiabaticconst, machnumber)
    return ((1/(machnumber**2))* ((2/(adiabaticconst+1))*(1+((adiabaticconst-1)/2)*(machnumber**2)))**((adiabaticconst+1)/(adiabaticconst-1)))

def mass_flow(burn_area, burn_rate):
    return rho*burn_area*burn_rate

def burn_rate(chamber_pressure):
    return a*(chamber_pressure**n)

def chamber_pressure(burn_area, throat_area):
    return ((burn_area/throat_area)*cstar*a*rho)**(1/(1-n))
