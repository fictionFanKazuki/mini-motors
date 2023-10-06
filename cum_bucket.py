import math
import scipy.integrate as integrate
import scipy.optimize as optimize
TIME_INC = 0.001
rho = 0.0018 #slugs/in^3
a = 0.0270
gamma = 1.2
cstar = 4890 # ft/s
R = 2000 #psi.ft^3.slugs^-1.R^-1
T0 = 4700 # R
g0 = 32.1741 #ft/s^2
n = 0.30
speed_of_sound = math.sqrt(gamma*R*T0) #ft/s
ambient_pressure = 14.6959
exitarea = 0 # i had to add these variables for some fucking optimization stuff
throatarea = 0

def simulate(num_grains, grain_length, grain_inner_diameter, throat_diameter, exit_area):
    mass = rho*num_grains*grain_length*math.pi*((3.239/2)**2 - (grain_inner_diameter/2)**2)
    throat_area = math.pi*((throat_diameter/2)**2)
    burn_area = 2*math.pi*(((3.239/2)**2) + grain_length)*num_grains
    t = 0
    impulse = 0
    global exitarea 
    exitarea = exit_area
    global throatarea
    throatarea = throat_area

    optimized_mach = optimize.root_scalar(optimize_this)
    exit_mach_number = optimized_mach["root"]
    print(optimized_mach["converged"])
    while burn_area>0:
        t += TIME_INC
        mass = mass - mass_flow(burn_area, burn_rate(chamber_pressure(burn_area, throat_area)))*t
        burn_area = burn_area - burn_rate(chamber_pressure(burn_area, throat_area))*t
        current_thrust = thrust(mass_flow(burn_area, burn_rate(chamber_pressure(burn_area, throat_area))), exit_velocity(exit_mach_number), exit_area, stagnationpressure(exit_mach_number), ambient_pressure)
        impulse = current_thrust*TIME_INC
    specific_impulse = specificimpulse(impulse, mass)
    return 1/specific_impulse

        
def optimize_this(mach_number):
    return (exitarea/throatarea)**2 - areamach_numberrelation(mach_number)

def specificimpulse(impulse, propmass):
    return ((impulse) / (propmass * g0))

def thrust(massflow, exhaustvelocity, exhaustarea, exhaustpressure, ambientpressure):
    return ((massflow * exhaustvelocity) + ((exhaustpressure - ambientpressure) * exhaustarea))

def exit_velocity(mach_number):
    return mach_number*speed_of_sound

def stagnationtemperature(mach_number):
    return (1+((gamma-1)/2)*(mach_number**2))

def stagnationpressure(mach_number):
    return ((1+((gamma-1)/2)*(mach_number**2))**(gamma/(gamma-1)))

def stagnationdensity(mach_number):
    return ((1+((gamma-1)/2)*(mach_number**2))**(1/(gamma-1)))

def areamach_numberrelation(mach_number):
    return ((1/(mach_number**2))* ((2/(gamma+1))*(1+((gamma-1)/2)*(mach_number**2)))**((gamma+1)/(gamma-1)))

def mass_flow(burn_area, burn_rate):
    return rho*burn_area*burn_rate

def burn_rate(chamber_pressure):
    return a*(chamber_pressure**n)

def chamber_pressure(burn_area, throat_area):
    return ((burn_area/throat_area)*cstar*a*rho)**(1/(1-n))

def rocket_equation(exit_velocity, initial_mass, empty_mass):
    return exit_velocity*math.log(initial_mass/empty_mass)

print(optimize.differential_evolution(simulate)["x"])