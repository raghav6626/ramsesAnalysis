import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sc
from rapack import *
import rapack.constants as const
from dataclasses import dataclass, field


#Data class for conversion factors from namelist units to cgs units.
@dataclass
class units:
    density: float
    time: float  #(M,N): M = time, N = radius
    length: float
    velocity: float
    Bfield: float
    pressure: float
    mass: float
    temperature: float

'''
Do not confuse yourself with the word units. 
Look at the definition of the three variables: {units_density, units_time, units_length}. 
They are defined in terms of cgs units. 
Thus, the word "units" is just a multiplicative factor that converts the physical parameters in the namelist 
into c.g.s variables. 
You can perform the inverse operation by taking the inverse of this multiplicative factor. 
In short: 
parameter file variables --> c.g.s units; multiply by factors returned by the following program. 
c.g.s units --> parameter file variables; multiply by 1/factor returned by the following program.

Checked the conversion factors again on 11/11/2027
'''
def get_namelist_to_cgs_factor(units_density, units_time, units_length, mean_molecular_weight = 1.4):
    units_velocity = units_length/units_time
    units_Bfield = (np.sqrt(4*np.pi))*np.sqrt(units_density*units_velocity**2); #there is a factor of 4pi because of how RAMSES treats magnetic fields.
    units_pressure = units_density*units_velocity**2;
    units_mass = units_density*(units_length)**3
    units_temp = ((units_velocity**2)*mean_molecular_weight)/(const.boltzmann_cgs*const.avogadro_no); 
    
    # print(f"Units mass [msun]: {units_mass/const.msun:.3e}");
    # print(f"Units time in [Myr]: {units_time/const.one_myr:.3e}");
    # print(f"Units velocity [km/sec]:{units_velocity/1e5:.3e}");
    # print(f"Units temperature/mu [K]:{units_temp:.3e}");
    # print(f"Units B-field [G]:{units_Bfield:.3e}");
    ramses_units = units(density = units_density, time = units_time, length = units_length, 
                         velocity = units_velocity, Bfield = units_Bfield, pressure = units_pressure, 
                         mass = units_mass, temperature = units_temp);
    return ramses_units;


# =========================================
# Defining the base-units
'''
These should match the parameter.nml file in RAMSES simulations
'''
#==========================================
units_density = 1.573391528e-26			# g/cm^3 = 1 H / cm^3;	 
units_time = 3.0868724e16
units_length = 3.086568025e21 

mean_molecular_weight = (const.avogadro_no*const.atomic_mass_constant)/0.76 #equating it to the default value in RAMSES
print('Mean molecular weight = ', mean_molecular_weight);

namelist_to_cgs_factor = get_namelist_to_cgs_factor(units_density, units_time, units_length, mean_molecular_weight = mean_molecular_weight);



# =========================================
# The resolution
#==========================================

boxlen = 164*const.one_kpc
lmin = 9;
lmax = 12;

print('Resolution of the base-grid cell: ', boxlen/(2**lmin*const.one_pc), 'pc');
print('Resolution of the finest-grid cell: ', boxlen/(2**lmax*const.one_pc), 'pc');


# =========================================
# Memory requirement per process/thread
#==========================================

nmpi = 48;
numPart = 212500;

maxNumberOfCells = 8**(lmax);
maxOcts =  maxNumberOfCells/8

maxOctsPerMPI = maxOcts/nmpi;
npartPerMPI = numPart/nmpi;

mem_per_core = 0.8*( 2*(maxOctsPerMPI/10**6) + (npartPerMPI/10**7) );

print('Max octs per mpi process:', int(maxOctsPerMPI));
print('Max particles per mpi process:', int(npartPerMPI));
print('Memory per cpu (GB):', mem_per_core)

# =========================================
# Here we have some example conversions 
# quantity_namelist = quantity_cgs/namelist_to_cgs_factor
#==========================================;

d_cgm_ppcm3 = 10**(-6) #cgm density
d_cgs = d_cgm_ppcm3*mean_molecular_weight/const.avogadro_no;
d_cgm_namelist = d_cgs/namelist_to_cgs_factor.density;
print(rf'd_CGM_namelist = {d_cgm_namelist:.3e}');


t_cgm_cgs = 10**6 #K 
p_cgm_cgs = t_cgm_cgs*(d_cgs*const.boltzmann_cgs*const.avogadro_no)/mean_molecular_weight; 
p_cgm_namelist = p_cgm_cgs/namelist_to_cgs_factor.pressure
print(rf'pressure_CGM_namelist = {p_cgm_namelist:.3e}');

m_cgs = 8.593e4 * const.solar_mass/8; #grams From Kim 2016. We divide by 8 due to RAMSES parameter file conventions
m_namelist = m_cgs/namelist_to_cgs_factor.mass;
print('Refine mass scale in code units = ', m_namelist);
print('Mass scale in Msun = ', namelist_to_cgs_factor.mass/const.solar_mass);


mg_disk = 8.593*(10**9)*const.msun; #Msun 
m_namelist = mg_disk/namelist_to_cgs_factor.mass
print('Gaseous disc mass in parameter file: ', m_namelist);

B_cgs = 1e-10; #0.1 N Gauss
B_namelist = B_cgs/namelist_to_cgs_factor.Bfield;
print('B-field in namelist user-units for B_cgs = ', B_cgs, '; B_namelist = ', B_namelist);

