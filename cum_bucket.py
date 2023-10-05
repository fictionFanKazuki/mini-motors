import math
import scipy.integrate as integrate
rho = 0.0018 #slugs/in^3
a = 0.0270
gamma = 1.2
cstar = 4890 # ft/s
R = 2000 #psi.ft^3.slugs^-1.R^-1
T0 = 4700 # R
g0 = 32.1741 #ft/s^2
n = 0.30
speed_of_sound = 1100 #ft/s


def simulate(num_grains, grain_length, grain_inner_diameter, throat_diameter):
    mass = rho*num_grains*grain_length*math.pi*((3.239/2)**2 - (grain_inner_diameter/2)**2)
    throat_area = math.pi*((throat_diameter/2)**2)
    burn_area = 2*math.pi*(((3.239/2)**2) + grain_length)
    while mass:
        meop = stagnationpressure(machnumber())
    pass


def specificimpulse(impulse, propmass, gravconst):
    return ((impulse) / (propmass * gravconst))

def thrust(massflow, exhaustvelocity, exhaustarea, exhaustpressure, ambientpressure):
    return ((massflow * exhaustvelocity) + ((exhaustpressure - ambientpressure) * exhaustarea))

def machnumber(flow_velocity, speedsound):
    return (flow_velocity / speedsound)

def stagnationtemperature(machnumber):
    return (1+((gamma-1)/2)*(machnumber**2))

def stagnationpressure(machnumber):
    return ((1+((gamma-1)/2)*(machnumber**2))**(gamma/(gamma-1)))

def stagnationdensity(machnumber):
    return ((1+((gamma-1)/2)*(machnumber**2))**(1/(gamma-1)))

def areamachnumberrelation(machnumber):
    return ((1/(machnumber**2))* ((2/(gamma+1))*(1+((gamma-1)/2)*(machnumber**2)))**((gamma+1)/(gamma-1)))

def mass_flow(burn_area, burn_rate):
    return rho*burn_area*burn_rate

def burn_rate(chamber_pressure):
    return a*(chamber_pressure**n)

def chamber_pressure(burn_area, throat_area):
    return ((burn_area/throat_area)*cstar*a*rho)**(1/(1-n))

def impulse(thrust, time):
    return integrate(thrust, 0, time)

def rocket_equation(exit_velocity, initial_mass, empty_mass):
    return exit_velocity*math.log(initial_mass/empty_mass)