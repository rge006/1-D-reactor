# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 15:44:28 2020

@author: Robert
"""

import csv
import math
import numpy as np

import cantera as ct
from matplotlib import pyplot as plt

#######################Reactor details########################################
#Initial design case set at 823 K and 1500 kPa
#For testing against the University of Ghent packed bed test cases

#Packed bed geometry:
    # height: 0.50355 m
    # width: 0.4521 m
    # length: 0.20 m
#Catalyst pellet geometry
    # diameter: 0.005 m
    # coating thickness: 0.0004 m
    # core material: Al2O3
    # core material density: 3950 kg/m³
    # coating density: 6110 kg/m³
    # pellet density (coating+core): 3506.25 kg/m³
#Gas phase properties:
    # kinetics and thermo: Pd/CO3O4
    # density: ideal gas law
#Solid phase properties:
    # kinetics: Pd/CO3O4
    # density: 3506.25 kg/m³
    # diameter: 0.0005 m / 500 micron
    # specific heat capacity: 5 J/kgK #This is used to reach a fast thermo-equilibrium, only useful to shorten simulation time
    # Volume fraction (1-epsilon): 0.5
    # porosity: 0
#Boundary conditions:
    # inlet composition: CO2/O2/H2O/CH4 = 0.96624156/0.01500375/0.01500375/0.00375094
    # inlet temperature: 823K (gas phase only)
    # inlet velocity: 4.427545835 m/s
    # inlet pressure: 15 bar
    # adiabatic walls set at the beginning and end of the packed bed reactor

f = open ('results.csv', 'w')
wtr = csv.writer(f,lineterminator='\n')

#Test case temperature and pressure
T1 = 843 #[K]
T_in = 843 #[K] Inlet temperature
P = 1422810.977 #[Pa]

#Reactor information
r_length = 0.20 #[m] Packed bed length
# r_i_diameter = 0.005 #[m] Packed bed diameter
r_i_height = 0.50355 #[m]
r_i_width = 0.4521 #[m]
velocity = 5.082990544 #[m/s] space velocity (superficial)
porosity = 0.5 #[m³space/m³reactor] Packed bed porosity
r_i_area = r_i_height * r_i_width #[m²] Packed bed cross-sectional area

cat_porosity = 0 #Porosity of the catalyst, this is an assumption

#Importing kinetics and properties file
#Pre-exponential in kinetics file need to be modified each time to accomodate for catalyst mass
cti_file = 'pdc - ualberta.cti'

# PBR will be simulated by a chain of CSTRs
n_reactors = 50

#######################Reactor details########################################

#######################Catalyst details#######################################
#Using 500 micron spherical particles
p_diameter = 0.005 #[m] particle diameter

#This is not necessary for this simulation, there are no surface reactions
cat_area_per_vol = 5.9e6 #the catalytic surface area per solid volume, this property is irrelevant to the kinetics used in this simulation

#Packed bed surface area, cited from https://www.sciencedirect.com/topics/engineering/ergun-equation
#Using 5000 micron spherical particles
A_p = 4*(p_diameter/2)**2*math.pi #[m²] particle surface area
V_p = 4/3 * (p_diameter/2)**3 * math.pi #[m³] particle volume
S = (1-porosity)/V_p * A_p #[m²_particle/m³_bed] specific surface area of the packed bed

#Characteristic length of catalyst particles
d_char = V_p / A_p

#Surface-volume particle diameter
d_sv = 6* V_p / A_p #Sauter mean diameter

#Catalyst pellet properties
#From test case
cp_cat = 5 #[J/kgK] Pellet specific heat capacity #This is used to achieve a fast steady-state, not the actual heat capacity of the catalyst pellet
k_cat = 5 #[W/mK] Pellet thermal conductivity, taken at 5 W/mK from F B Lewis and N H Saunders 1973 J. Phys. C: Solid State Phys. 6 2525
rho_cat = 3506.25 #[kg_pellet/m³_pellet] Pellet density, spent catalyst bulk density

#Test case did not include radiation energy transfer
epsilon_r = 0.5 #[-] Pellet surface emissivity, half of all radiation is absorbed and half is refrlected

#Unused parameters in this simulation
# r_pore = 1.4E-9 #[m] Catalyst pore radius, converted from 1.4 nm
# tau_epsilon = 4.0 #[-] Tortuosity / porosity of catalytic layer
# t_cat = 0.15E-3 #[m] Thickness of catalyst layer, converted from 0.15 mm

#######################Catalyst details#######################################

#######################Diffusivity parameters#################################
#This section was not used for this simulation, as the diffusivity is calculated as D = mu_g/rho_g
#Molar mass
# m_ch4 = 16.04 #[g/mol], methane molar mass
# m_o2 = 16 #[g/mol], oxygen gas molar mass
# m_co2 = 44.01 #[g/mol] CO2 molar mass
# m_co = 28.01 #[g/mol] CO molar mass
# m_h2o = 18.02 #[g/mol] H2O molar mass
# m_n2 = 14 #[g/mol], nitrogen gas molar mass

# #Average collision diameter
# sigma_ch4 = 3.780 #[angstrom, 10^-10 m]
# sigma_o2 = 3.433 #[angstrom, 10^-10 m]
# sigma_co2 = 3.996 #[angstrom, 10^-10 m]
# sigma_co = 3.59 #[angstrom, 10^-10 m]
# sigma_h2o = 3.996 # this is the same as co2, I didn't find the value for h2o
# sigma_n2 = 3.667 #[angstrom, 10^-10 m]

# #For calculating collision integral
# ek_ch4 = 154
# ek_o2 = 113
# ek_co2 = 190
# ek_co = 110
# ek_h2o = 190 # this is the same as co2, I didn't find the value for h2o
# ek_n2 = 99.8
#######################Diffusivity parameters#################################

#######################Simulation set-up######################################
#Import gas model and set the initial conditions
bulk_gas = ct.Solution(cti_file, 'gas') #Phase object to represent the gas phase
surf_gas = ct.Solution(cti_file, 'gas') #Phase object to represent the solid "gas" phase

#Setting initial conditions for both phase objects
bulk_gas.TPY = T1, P, 'CH4:0.001881, O2:0.007523, H2O:0.019216, CO2:0.971381' #based on 2 vol% O2, adding equilmolar CH4 to consume half of total O2
surf_gas.TPY = T1, P, 'CH4:0.001881, O2:0.007523, H2O:0.019216, CO2:0.971381' #based on 2 vol% O2, adding equilmolar CH4 to consume half of total O2

# import the surface model
surf = ct.Interface(cti_file,'surface1', [surf_gas])
surf.TP = T1, P #Solid phase set to be same T & P as gas phase to shorten simulation time, ie assuming packed bed is preheated

#Simulated reactor information
dx = (r_length)/(n_reactors) #Length of one CSTR

#Volume of 1 finite element
r_volume = dx * r_i_area

#catalytic surface area of one reactor (finite volume)
cat_area = cat_area_per_vol * r_volume * (1-porosity) #This property is irrelevant to the kinetics used in this simulation

dt = 5e-6 #Time step size

#Header for csv output file
wtr.writerow(np.concatenate(('Time', ['T_g' for i in range (n_reactors)],['T_s' for i in range (n_reactors)],[np.concatenate((bulk_gas.species_names,surf_gas.species_names,surf.species_names)) for i in range (n_reactors)], 'out_P'),axis=None))
#######################Simulation set-up######################################

#######################Equation definition####################################

#Making list of molar weights for gas phase
MM_g = bulk_gas.molecular_weights

# MM = np.hstack((surf_gas.molecular_weights,surf.molecular_weights)) #array of molecular weights for all species

x_axis = np.arange(0, r_length+dx, dx) #x-axis variable used for plotting results

#To calculate gas phase heat capacity using cantera
def cp_eq(T, p, Y):
    # print(T)
    surf_gas.TPY = T, p, Y
    return surf_gas.cp_mass

def rho_eq_(T, p, Y):
    y_mm = np.divide(Y,MM_g)
    avg_MM = 1/np.sum(y_mm)
    rho = p * avg_MM / (8.314 * T * 1000)
    return rho

def mu_eq(T, p, Y):
    surf_gas.TPY = T, p, Y
    return surf_gas.viscosity

def k_eq(T, p, Y):
    surf_gas.TPY = T, p, Y
    return surf_gas.thermal_conductivity

def R_eq(T_g_r, T_s_r, p_r, Y_g_r, Y_s_r, Y_surf_r):
    
    #Setting up reaction phase conditions
    bulk_gas.TPY = T_g_r, p_r, Y_g_r #Setting up gas phase T, P and species fraction, T_g is used because the reaction temperature should match the gas phase temperature
    surf_gas.TPY = T_s_r, p_r, Y_surf_r #Setting up solid phase T, P and species fraction, T_s is used because the reaction temperature should match the solid phase temperature

    #Setting up solid reactor simulation
    rsurf.coverages = Y_s_r #Setting up occupied active site fraction by surface species in the solid phase
    # rsurf.area = cat_area

    r.volume = r_volume #Setting up reactor volume
    r.syncState() #Sync reactor state to that of the one set
    sim = ct.ReactorNet([r]) #Creating solid reactor simulation

    #Solid reactor info prior to reaction
    m_c = r.mass*r.Y #Array holding mass of each gas phase species before reaction
    coverage_b = rsurf.coverages #Occupied active site fraction by surface species in the solid phase before reaction
    e_b = r.thermo.cp_mass*r.thermo.T*r.mass #Total enthalpy of the reactor before reaction [J]

    # m_b = r.mass #Solid reactor mass before reaction [kg]

    #Integrating the solid reactor object by the timestep dt
    sim.atol = 1E-10 #Absolute tolerance of 1E-10 is used, this is smaller than the default tolerance to allow the integration to occur
    sim.advance(dt) #Integrating the Cantera batch reactor by one time step #Solid reactor

    #Solid reactor info post reaction
    coverage_p = rsurf.coverages #Occupied active site fraction of surface species in the solid phase after reaction
    e_p = r.thermo.cp_mass*r.thermo.T*r.mass #Total solid reactor enthalpy of the reactor after reaction [J]

    #Add in a cpT dm/dt term to the q term as a correction (applying multiplication rule for differentiation)

    #Change in solid reactor in energy, gas phase species, and surface coverage
    q = (e_p - e_b)/ (r_volume*dt)#Calculating enthalpy change in the solid reactor due to reaction per unit volume [J/m³s], this is the energy generation due to reaction
    rxn_g = (r.mass*r.Y-m_c)/(r.mass*dt) #Calculating change in gas phase species mass fraction due to reaction [-/s], mass used is from after reaction
    rxn_s = (coverage_p-coverage_b)/dt #Calculating change in occupied active site fraction by surface species in the solid phase due to reaction [-/s]

    #Setting up gas reactor simulation
    gas_r.volume = r_volume
    gas_r.syncState()
    gas_sim = ct.ReactorNet([gas_r])
    
    #Gas reactor info prior to reaction
    m_c_g = gas_r.mass*gas_r.Y #Array holding mass of each gas phase species before reaction
    e_b_g = gas_r.thermo.cp_mass*gas_r.thermo.T*gas_r.mass #Total enthalpy of the reactor before reaction [J]

    # m_b_g = gas_r.mass #Gas reactor mass before reaction [kg]

    #Integrating the gas reactor object by the timestep dt
    gas_sim.atol = 1E-10
    gas_sim.advance(dt) #Integrating the Cantera batch reactor by one time steop #Gas reactor

    #Gas reactor info post reaction
    e_p_g = gas_r.thermo.cp_mass*gas_r.thermo.T*gas_r.mass #Total gas reactor enthalpy of the reactor after reaction [J]
    
    #Change in gas reactor in energy, gas phase species
    gas_q = (e_p_g - e_b_g)/ (r_volume*dt) #Calculating enthalpy change in the gas reactor due to reaction per unit volume [J/m³s], this is the energy generation due to reaction
    rxn_gas = (gas_r.mass*gas_r.Y-m_c_g)/(gas_r.mass*dt) #Calculating change in gas phase species mass fraction due to reaction [-/s], mass used is from after reaction

    return np.hstack((q, rxn_g, rxn_s, gas_q, rxn_gas))

#Calculating for radiation heat transfer from solid phase to gas phase
def q_rad(surf_T, epsilon):
    """
    Parameters
    ----------
    surf_T : Float, units of K
        Solid phase temperature
    epsilon : float, units of [-]
        Packed bed porosity
    Returns
    -------
    q_rad : Float, units of W/m²
        Radiation energy transfer
    Equation
    ---------
    .. math::
        S_r = 1 + 1.84 * (1 - \epsilon) + 3.15 * (1 - \epsilon)^2 + 7.2 * (1 - \epsilon)^3 for \space \epsilon > 0.3\n
        \\beta = 1.5 * \epsilon_r * (1 - \epsilon) * S_r / d_{sv} \n
        q^{''}_{rad} = -16 * \sigma * T_s^3 / (3 * \\beta)
    >>> where:
        S_r is the scaling factor [m]
        epsilon is the porosity [-]
        epsilon_r is the emissivity [-]
        d_sv is the surface to volume diameter [m]
        beta is the extinction coefficient [-]
        sigma is the Stefan-Boltzmann constant [W/(m²·K^-4)]
        k_g is the gas phase thermal conductivity [W/mK]
        c_p_g is the gas phase constant pressure heat capacity [J/kg]
        
    Applicability
    -------------
    Used to calculate radiation heat transfer (without the dT/dx) [W/m²].
    """
    sigma_SB = 5.67E-8 #[W/(m^2*K^-4)] Stefan-Boltzmann constant
    #scaling factor for porosity > 0.3
    S_r = 1 + 1.84 * (1 - epsilon) + 3.15 * (1 - epsilon)**2 + 7.2 * (1 - epsilon)**3
    #extinction coefficient
    beta_e = 1.5 * epsilon_r * (1 - porosity) * S_r / d_sv
    q_rad = 16 * sigma_SB * surf_T**3 / (3 * beta_e)
    return q_rad

#Calculating diffusion coefficient of species 1 in species 2
# def D_eq(T, p, MM1, MM2, sigma, omega):
#     #(conversion from cm²/s to m²/s) #CH4 in N2, [m²/s]
#     D = 1.858E-3 * T**(3/2) * (1/MM1 + 1/MM2)**(1/2) / ( (p/101325) * sigma**2 * omega) * 1/10000
#     return D

#Calculating omega of CH4/O2 in N2 using a linear equation, this is to change the diffusivity as temperature increases
#Obtained by linearizing the collision integral table over T = 750-1500
# def omega_eq(kTe):
#     return - 0.0151 * kTe + 0.9009

#Calculating inter-phase heat transfer coefficient
def h_eq(rho, u, cp, mu, k):

    """
    Parameters
    ----------
    u : Float, units of m/s
        Superficial gas velocity
    Returns
    -------
    h_inter : Float, units of W/m²·K
        Interphase heat transfer coefficient
    Equation
    ---------
    .. math::
        Re = \\rho_g * u / (S * f * \mu_g) \n
        J_H = 0.91 * f * Re^{-0.51} for \space 0.01 < Re < 50 \n
        J_H = 0.61 * f * Re^{-0.41} for \space 50 < Re \n
        h_{inter} = J_H * \\rho_g * u * c_{p_{g}} * (\mu_g * c_{p_{g}} / k_g) ^ {-2/3}
    >>> where:
        Re is the Reynolds number
        J_H is the Colburn J-Factor
        f is the shape factor of the catalyst pellet, 0.91 for cylindrical pellets [-]
        u is the superficial gas velocity [m/s]
        S is the specific surface area of the packed bed [m²/m³]
        mu is the gas viscosity [Pa·s]
        rho_g is the gas phase density [kg/m³]
        k_g is the gas phase thermal conductivity [W/mK]
        c_p_g is the gas phase constant pressure heat capacity [J/kg]
        
    Applicability
    -------------
    For use to calculate interphase heat transfer coefficient.
    """
    #Calculating interphase heat transfer coefficient
    #Gunn analogy used for this simulation
    Re = rho * u * d_char / (mu*porosity)
    Pr = cp * mu / k
    Nu_s = (7 - 10*porosity + 5*porosity**2) * (1 + 0.7*Re**0.2*Pr**(1/3)) + (1.33 - 2.4*porosity + 1.2*porosity**2)*Re**0.7*Pr**(1/3)
    h_inter_ = Nu_s * k / d_char #[W/m²]
    return h_inter_

def u_eq(rho_g_):
    #Calculating the gas flow velocity using the conservation of mass [m/s]
    return mass_flow_rate / (r_i_area * rho_g_)

def p_eq(p, u_, mu_g):
    #Using the Ergun equation to calculate pressure drop [Pa]
    return p - (150 * mu_g/d_sv**2 * (1-porosity)**2/porosity**3*u_ + 1.75*mass_flow_rate/(d_sv*r_i_area) * (1-porosity)/porosity**3*u_)*dx

#For calculating mass conservation
def y_eq (y_g, y_g_prev, rho_g, s_mx, r_p, r_g, r_s):
    #Calculating advection term
    np.subtract(y_g, y_g_prev)/rho_g
    m_adv = mass_flow_rate/r_i_area * np.subtract(y_g, y_g_prev)/rho_g/dx
    
    # Calculating the accumulation term in the gas phase by adding up advection, reaction and interphase mass transfer terms
    dcdt_g = r_g - m_adv + s_mx / rho_g
    
    #Calculating the accumulation term in the pore phase by adding up reaction and interphase mass transfer term
    dcdt_p = r_p - s_mx/rho_g
    
    #Calculating the accumulation term in the solid phase by adding up reaction term
    dcdt_s = r_s
    return np.concatenate((dcdt_g, dcdt_p, dcdt_s))

def k_inter_eq(T, p, rho, u, cp, mu, k):
    #Calculating interphase mass transfer coefficient
    #Gunn analogy used here
    D = mu/rho #As stated in supporting information of the University of Ghent paper http://pubs.acs.org/doi/suppl/10.1021/acs.energyfuels.0c02824/suppl_file/ef0c02824_si_001.pdf
    Re = rho * u * d_char / (mu*porosity)
    Pr = cp * mu / k
    Sc = 1 #Set to 1 as stated in the supporting information
    Nu_s = (7 - 10*porosity + 5*porosity**2) * (1 + 0.7*Re**0.2*Pr**(1/3)) + (1.33 - 2.4*porosity + 1.2*porosity**2)*Re**0.7*Pr**(1/3) #Gunn heat transfer analogy
    
    Sh = Nu_s * (Sc/Pr)**(1/3)
    k_inter = Sh * D / d_char

    k_ = k_inter * np.ones(bulk_gas.n_species)
    return k_

def smx_eq(k_inter,rho_g, rho_s, y_g, y_s):
    #Calculating interphase mass transfer source term
    s_mx = k_inter * S * (rho_g*y_s - rho_g*y_g)
    return s_mx

#######################Equation definition####################################

#######################Variable declaration###################################

#Declare a cantera batch reactor object for solid phase, install a surface
r = ct.IdealGasConstPressureReactor(surf_gas, energy="on")
r.volume = r_volume
rsurf = ct.ReactorSurface(surf, r, A=cat_area)

#Declare a cantera batch reactor object for gas phase
gas_r = ct.IdealGasConstPressureReactor(bulk_gas, energy="on")
gas_r.volume = r_volume

#Inlet states
#Bulk gas properties is used as it's the bulk gas phase
Y_g_in = bulk_gas.Y #setting up inlet gas phase composition
cp_in = bulk_gas.cp_mass #set up inlet gas phase heat capacity [J/kgK]
H_in = bulk_gas.density_mass * cp_in * T_in #set up inlet gas phase enthalpy [kg/m³] * [J/kg]
rho_g_in = rho_eq_(T1, P, Y_g_in) #setting up inlet gas density

#Use these if the initial state is different than the inlet states
surf_gas.TPY = T1, P, 'CH4:0.001881, O2:0.007523, H2O:0.019216, CO2:0.971381' #based on 0.35 CH4 mole equivalence ratio to O2
# bulk_gas.TPY = T1, P, 'N2:1'

#Initial states arrays
Y_g_0 = bulk_gas.Y #initial gas phase composition
Y_surf = Y_g_0.copy() #Setting up array to track catalyst pore phase mass fractions
Y_s_0 = rsurf.coverages #mass fractions in the solid phase tracked as fraction of catalytic sites occupied by each species
Y0 = np.hstack((Y_g_0, Y_surf, Y_s_0)) #setting the initial gas, pore, and solid phase composition

#Setting up initial conditions across all cells
T_g_0 = T1 * np.ones(n_reactors) #setting up initial gas phase temperature of all cells
T_s_0 = T1 * np.ones(n_reactors) #setting up initial solid phase temperature of all cells
rho_g_0 = rho_eq_(T1, P, Y_g_0) #setting up initial gas phase density
cp_0 = bulk_gas.cp_mass #setting up initial gas phase heat capacity
k_0 = bulk_gas.thermal_conductivity #setting up initial gas phase thermal conductivity
h_0 = np.multiply(cp_0, T1) * np.ones(n_reactors) #setting up initial gas phase enthalpy
H_0 = np.multiply(rho_g_0, h_0) #setting up initial gas phase rho*H [kg/m³] * [J/kg]
mu_0 = mu_eq(T1, P, Y_g_0) #setting up initial gas phase viscosity

#Tracking states, global variables that are used to track the properties of the gas phase in each cell
u_ = velocity * np.ones(n_reactors) #setting up gas phase velocity of all cells (used for tracking)
p_ = P * np.ones(n_reactors) #setting up gas phase pressure of all cells (used for tracking)
rho_g_ = rho_g_0 * np.ones(n_reactors) #setting up gas phase density of all cells (used for tracking)
cp_g_ = cp_0 * np.ones(n_reactors) #setting up gas phase heat capacity of all cells (used for tracking)
mu_g = mu_0 * np.ones(n_reactors) #Setting up gas phase viscosity of all cells (used for tracking)
k_g = np.ones(n_reactors) #Setting up gas phase thermal conductivity of all cells (used for tracking)

rho_s_ = rho_g_0 * np.ones(n_reactors) #setting up gas phase density of all cells (used for tracking)

#Initial conditions used for solving PDE
Y = np.hstack([Y0[:] for i in range (n_reactors)]) #initial gas and solid phase composition for all reactors
y0 = np.hstack((H_0, T_s_0, Y)) #setting up array to pass into PDE

mass_flow_rate = velocity * rho_g_in * r_i_area #defining mass flowrate

rad_cont = np.zeros(n_reactors)
cond_cont = np.zeros(n_reactors)
conv_cont = np.zeros(n_reactors)
tot_cont = np.zeros(n_reactors)
q_tracker = np.zeros(n_reactors)

#Turning off gas phase reactions in the solid "gas" phase
surf_gas.set_multiplier(0)
#######################Variable declaration###################################

def consv_eqs(dt, y_):
    
# =============================================================================
#     Radiation was turned ON at PARTICLE EMISSIVITY = 1 for this simulation
# =============================================================================
    
    solution = np.zeros((n_reactors, 2)) #Solution array for holding rho_g*H and T_s
    c_solution = np.zeros((n_reactors, surf.n_total_species)) #Solution array for holding mass fraction (both gas and solid phase)

    H_ = y_[:n_reactors] #Gas phase rho_g * H_g
    
    T_s = y_[n_reactors:n_reactors*2] #Solid phase temperature

    # c_ = np.reshape(y_[n_reactors*2:], (-1, surf.n_total_species)) #Gas phase mass fraction and solid phase coverage fraction
    
    c_ = np.reshape(y_[n_reactors*2:], (-1, (bulk_gas.n_total_species+surf.n_total_species))) #Gas phase mass fraction, catalyst pore phase mass fraction, and solid phase coverage fraction
    
    c_g = c_.T[:][0:bulk_gas.n_total_species] #Gas phase mass fraction
    c_pore = c_.T[:][bulk_gas.n_total_species:(bulk_gas.n_total_species+surf_gas.n_total_species)] #Catalyst pore phase mass fraction
    c_s = c_.T[:][(bulk_gas.n_total_species+surf_gas.n_total_species):] #Solid phase coverage fraction

    # c_g = c_.T[:][0:surf_gas.n_total_species] #Gas phase mass fraction
    # c_s = c_.T[:][surf_gas.n_total_species:] #Solid phase coverage fraction

    global rho_g_, cp_g_, mu_g, k_g, u_, p_#, rad_cont

    T_g = np.divide(H_, np.multiply(rho_g_,cp_g_)) #Calculating gas phase temperature from rho_g*H_g

    Y_g = c_g.T.copy()
    Y_pore = c_pore.T.copy()
    Y_s = c_s.T.copy()

    #Calculating gas phase properties
    cp_g_ = list(map(cp_eq, T_g, p_, Y_g)) #Calculating all species gas phase heat capacities across all cells
    mu_g = list(map(mu_eq,T_g, p_, Y_g)) #Calculating all species viscosity across all cells
    k_g = list(map(k_eq,T_g, p_, Y_g)) #Calculating all species thermal conductivity across all cells
    h_inter = list(map(h_eq, rho_g_, u_, cp_g_, mu_g, k_g)) #Calculating convective heat transfer coefficient across all cells
    rho_g_ = list(map(rho_eq_, T_g, p_, Y_g)) #Calculating bulk gas phase density across all cells

    rho_s_ = list(map(rho_eq_, T_s, p_, Y_pore))

    #Calculating interphase mass transfer
    k_inter = list(map(k_inter_eq, T_g, p_, rho_g_, u_, cp_g_, mu_g, k_g)) #Calculating interphase mass transfer coefficient across all cells
    s_mx = np.array(list(map(smx_eq, k_inter, rho_g_, rho_s_, Y_g, Y_pore))) #Mass change from interphase exchange

    #Calculating terms related to reaction
    R = np.array(list(map(R_eq, T_g, T_s, p_, Y_g, Y_s, Y_pore))) #Calculating all species reaction rates and heat of reaction across all cells
    
    #Re-arranging array to separate into heat of reaction, gas phase and solid phase reaction rate
    q = R.T[:][0] #Rate of change for energy in solid reactor
    rxn_g = R.T[1:surf_gas.n_species+1][:].T #Variable used to track gas phase species rate of change in solid reactor
    rxn_s = R.T[surf_gas.n_species+1:surf.n_total_species+1][:].T  #Variable used to track solid phase coverage fraction rate of change in solid reactor
    gas_q = R.T[surf.n_total_species+1:surf.n_total_species+2][:].T #Rate of change for energy in gas reactor
    rxn_gas = R.T[surf.n_total_species+2:][:].T #Variable used to track gas phase species rate of change in gas reactor

    #Momentum conservation
    prev_P = np.hstack((P, p_[0:-1])) #Making new array that contains inlet pressure and all pressures except for last cell
    
    u_next = np.array(list(map(u_eq, rho_g_))) #Calculating gas velocity across all cells
    p_next = np.array(list(map(p_eq, prev_P, u_, mu_g))) #Calculating pressure across all cells
    
    #Setting velocity and pressure to the calculated values
    u_ = u_next
    p_ = p_next
    
    prev_y = np.vstack((Y_g_in, Y_g[0:-1])) #Making new array that contains inlet mass fraction and all mass fractions except for last cell
    
    #Mass conservation, d/dt (c_g_i[i])
    #Evaluating mass conservation equations for gas and solid phase
    c_solution = np.array(list(map(y_eq, Y_g, prev_y, rho_g_, s_mx, rxn_g, rxn_gas, rxn_s)))

    #cell 1##########
    #gas phase
    eg_adv = (H_[0]*u_[0]-velocity*H_in)/(dx) #Calculating for advection transfer
    inter_g = h_inter[0]*(T_s[0]-T_g[0])*S #Calculating for interphase heat transfer
    dH_gdt = (-eg_adv + inter_g + gas_q[0]) / (porosity) #Calculating for accumulation

    #Solid, (1-epsilon)*rho_s * d/dt(H_s[i])
    rad_s = q_rad(T_s[0],porosity) * (T_s[1]-T_s[0])/dx**2 #Calculating for radiation transfer, Adiabatic wall BC
    es_cond = (1-porosity)*(1-cat_porosity) * k_cat * (T_s[1]-T_s[0])/dx**2 #Calculating for conduction transfer, Adiabatic wall BC
    inter_s = h_inter[0]*(T_g[0]-T_s[0])*S #Calculating for interphase heat transfer
    dT_sdt = (es_cond + rad_s + inter_s + q[0]) / ((1-porosity)*(1-cat_porosity)*rho_cat*cp_cat) #Calculating for accumulation
    solution[0][0] = dH_gdt
    solution[0][1] = dT_sdt

    #cells 1 to n-1
    for i in range (1, n_reactors-1):
        #Energy conservation
        #Gas, epsilon * d/dt(rho_g[i]*H_g[i])
        eg_adv = (H_[i]*u_[i]-H_[i-1]*u_[i-1])/dx #Calculating for advection transfer
        inter_g = h_inter[i]*(T_s[i]-T_g[i])*S #Calculating for interphase heat transfer
        dH_gdt = (-eg_adv + inter_g + gas_q[i]) / (porosity) #Calculating for accumulation
        
        #Solid, (1-epsilon)*rho_s * d/dt(H_s[i])
        rad_s = q_rad(T_s[i],porosity) * (T_s[i+1]-2*T_s[i]+T_s[i-1])/dx**2 #Calculating for radiation transfer, Adiabatic wall BC
        es_cond = (1-porosity)*(1-cat_porosity) * k_cat * (T_s[i+1]-2*T_s[i]+T_s[i-1])/dx**2 #Calculating for conduction transfer, Adiabatic wall BC
        inter_s = h_inter[i]*(T_g[i]-T_s[i])*S #Calculating for interphase heat transfer
        dT_sdt = (es_cond + rad_s + inter_s + q[i]) / ((1-porosity)*(1-cat_porosity)*rho_cat*cp_cat) #Calculating for accumulation
    
        solution[i][0] = dH_gdt
        solution[i][1] = dT_sdt

    #cell n
    #Energy conservation
    #Gas, epsilon * d/dt(rho_g[i]*H_g[i])
    eg_adv = (H_[-1]*u_[-1] - H_[-2]*u_[-2]) / dx #Calculating for advection transfer
    inter_g = h_inter[-1]*(T_s[-1]-T_g[-1])*S #Calculating for interphase heat transfer
    dH_gdt = (-eg_adv + inter_g + gas_q[-1]) / (porosity) #Calculating for accumulation
    
    #Solid, (1-epsilon)*rho_s * d/dt(H_s[i])
    rad_s = q_rad(T_s[-1],porosity) * (T_s[n_reactors-2]-T_s[-1])/dx**2 #Adiabatic wall BC
    es_cond = (1-porosity)*(1-cat_porosity) * k_cat * (T_s[n_reactors-2]-T_s[-1])/dx**2 #Adiabatic wall BC
    inter_s = h_inter[-1]*(T_g[-1]-T_s[-1])*S #Calculating for interphase heat transfer
    dT_sdt = (es_cond + rad_s + inter_s + q[-1]) / ((1-porosity)*(1-cat_porosity)*rho_cat*cp_cat) #Calculating for accumulation
    
    solution[-1][0] = dH_gdt
    solution[-1][1] = dT_sdt

    states_ = np.ravel(solution, 'F')
    mass_c = np.ravel(c_solution, 'C')
    
    y_out = np.concatenate((states_,mass_c))
    return y_out

def euler_ (fun, t_span, y0, dt):
    t = 0
    dy_e = [] #explicit solution
    dy_i = [] #implicit solution
    new_y = []
    new_new_y = []
    
    t_arr = []
    y_arr = []
    
    solArr = []
    
    counter = 0
    plotcount = 0
    global rho_g_, cp_g_
    
    while t <= t_span:
    #Returning d/dt
        dy_e = fun(dt, y0)
        new_y = np.add(y0, np.multiply(dy_e,dt))
        for i in range (4):
            dy_i = fun(dt, new_y)
            new_new_y = np.add(y0, np.multiply(np.add(dy_e,dy_i)/2,dt))
            new_y = new_new_y
        t = t + dt
        # co_ = (u_*dt)/dx
        # print("courant number")
        # print(co_)
        t_arr.append(t)
        y_arr.append(np.concatenate((np.divide(y0[:n_reactors],np.multiply(cp_g_,rho_g_)),new_new_y[n_reactors:], p_[-1]), axis=None))
        y0 = new_new_y
        # print("loop " + str(counter))
        counter+=1
        
        print("time: " + str(t))

        if counter % 50 == 0:
            
            plotcount+=1
            
            H_ = y0[:n_reactors] #Gas phase rho_g * H_g
            T_s = y0[n_reactors:n_reactors*2] #Solid phase temperature
        
            # c_ = np.reshape(y0[n_reactors*2:], (-1, surf.n_total_species)) #Gas phase mass fraction and solid phase coverage fraction
            c_ = np.reshape(y0[n_reactors*2:], (-1, (bulk_gas.n_total_species+surf.n_total_species))) #Gas phase mass fraction, catalyst pore phase mass fraction, and solid phase coverage fraction
    
            c_g = c_.T[:][0:bulk_gas.n_total_species] #Gas phase mass fraction
            # c_pore = c_.T[:][bulk_gas.n_total_species:(bulk_gas.n_total_species+surf_gas.n_total_species)] #Catalyst pore phase mass fraction

            c_g = c_.T[:][0:surf_gas.n_total_species] #Gas phase mass fraction
            
            c_o2 = np.hstack((0.007523,c_g[:][3]))
            c_ch4 = np.hstack((0.001881,c_g[:][13]))
            # p_o2 = np.hstack((0.048746,c_pore[:][2]))
            # p_ch4 = np.hstack((0.0977562,c_pore[:][10]))
            # c_co2 = np.hstack((0.966242,c_g[:][15]))
            # c_h2 = np.hstack((0,c_g[:][0]))
            
            x_arr = np.arange(0,r_length+dx,dx)
            Tg_arr = np.hstack((T1,np.divide(H_,np.multiply(cp_g_,rho_g_))))
            Ts_arr = np.hstack((T_s[0],T_s))
            
            fig,ax = plt.subplots()
            ax.set_title('Axial temperature at t = ' + str(t))
            ax.set_xlabel('Axial distance [m]')
            ax.set_ylabel('Temperature [K]')
            ax.plot(x_arr, Tg_arr, c='b', label='gas phase T')
            ax.plot(x_arr, Ts_arr, c='r', label='solid phase T')
            ax.legend(loc='upper left')
            ax.set_xlim([0,r_length])
            ax.set_ylim([800,1200])
            ax2 = ax.twinx()
            ax2.set_ylabel('Mass fraction [-]')
            ax2.plot(x_arr, c_o2, c='y', label='O2 mass fraction')
            ax2.plot(x_arr, c_ch4, c='c', label='CH4 mass fraction')
            # ax2.plot(x_arr, c_co2, c='m', label='CO2 mass fraction')
            
            # ax2.plot(x_arr, c_co2, c='g', label='CO2 mass fraction')
            # ax2.plot(x_arr, c_co, c='m', label='CO mass fraction')
            # ax2.plot(x_arr, c_h2, c='k', label='H2 mass fraction')
            ax2.legend(loc='upper right')
            ax2.set_ylim([0,0.02])
            
            fileNameTemplate = r'new plots\Plot{0:02d}.png'
            
            fig.savefig(fileNameTemplate.format(plotcount), format='png')
            
            fig.clf()
            plt.close()
        
    solArr = np.column_stack((t_arr,y_arr))
    wtr.writerows(solArr)
    f.close()
    # print(co_)

    return solArr

# plt.show()

#0.15s
# y0 = 516377.7507,516642.1997,516912.8182,517192.1015,517479.8709,517775.9958,518080.5409,518393.6853,518715.6459,519046.5499,519385.992,519731.3299,520078.7073,520422.5252,520777.1547,521143.1978,521521.2564,521911.9736,522316.0465,522734.2301,523167.3396,523616.2588,524081.9414,524565.4145,525067.7896,525590.257,526134.0957,526700.6692,527291.4281,527907.9013,528551.6872,529224.4374,529927.8259,530663.5129,531433.0852,532237.9852,533079.4014,533958.1525,534874.5249,535828.092,536817.5234,537840.3776,538892.9192,539969.9563,541064.691,542168.4138,543269.4565,544349.279,545368.4055,546217.0863,977.0868919,978.7628461,980.6910849,982.6977304,984.7634178,986.8879793,989.0737449,991.3234233,993.6395602,996.0232749,998.4692177,1000.945524,1003.538837,1006.241116,1009.038235,1011.933295,1014.93197,1018.040708,1021.266586,1024.617285,1028.101191,1031.727438,1035.505982,1039.447666,1043.564298,1047.868727,1052.374899,1057.097925,1062.054129,1067.261044,1072.737391,1078.502982,1084.578524,1090.985311,1097.744729,1104.877553,1112.402968,1120.337267,1128.69221,1137.473055,1146.67626,1156.287097,1166.277219,1176.602363,1187.199737,1197.982708,1208.823851,1219.494524,1229.450691,1237.082064,5.33E-09,1.43E-10,0.048511153,7.57E-11,0.000153771,4.90E-10,2.58E-11,8.11E-09,7.57E-05,6.87E-05,0.09767708,5.42E-06,1.43E-13,1.12E-06,1.51E-09,1.68E-09,3.10E-14,5.94E-07,6.07E-12,3.89E-06,5.12E-10,4.92E-08,3.53E-14,2.55E-13,1.33E-10,8.90E-14,5.84E-13,1.17E-11,2.52E-18,0.853497829,7.19E-09,7.49E-11,0.048272468,6.06E-11,0.000310108,9.96E-10,2.48E-11,1.08E-09,0.000152776,0.000138594,0.097596712,1.16E-05,5.53E-13,2.26E-06,2.23E-09,7.57E-10,1.13E-13,1.05E-06,5.20E-11,7.17E-06,7.03E-10,9.80E-08,3.48E-13,5.10E-13,2.36E-10,3.22E-13,3.09E-12,2.21E-11,8.37E-18,0.853497829,0.128011221,0.21622275,0.654120547,5.62E-09,7.68E-08,0.001645316,1.05E-09,4.48E-10,8.24E-08,4.14E-16,3.69E-11,1.53E-08,3.15E-10,0.048270032,1.58E-10,0.000311071,7.78E-10,8.92E-11,2.12E-08,0.000150608,0.000143536,0.097598127,8.12E-06,2.77E-13,1.63E-06,4.07E-09,3.65E-09,6.41E-14,1.30E-06,1.21E-11,8.83E-06,1.25E-09,1.31E-07,8.21E-14,8.12E-13,3.99E-10,2.25E-13,1.39E-12,3.63E-11,1.06E-17,0.853497829,1.77E-08,1.04E-10,0.048025133,7.66E-11,0.000470852,1.10E-09,8.22E-11,2.54E-09,0.000226684,0.000219565,0.097518034,1.23E-05,6.21E-13,2.34E-06,4.90E-09,1.01E-09,1.31E-13,1.74E-06,8.58E-11,1.24E-05,1.20E-09,2.11E-07,7.53E-13,1.37E-12,5.74E-10,4.91E-13,6.67E-12,5.72E-11,2.33E-17,0.853497829,0.12735165,0.211830341,0.659140729,8.25E-09,7.89E-08,0.001677067,1.10E-09,4.96E-10,1.25E-07,8.86E-16,7.76E-11,2.83E-08,4.26E-10,0.048022752,2.11E-10,0.00047195,9.43E-10,1.84E-10,3.32E-08,0.000224349,0.000224158,0.097518392,9.51E-06,3.75E-13,1.84E-06,7.28E-09,4.85E-09,9.02E-14,2.04E-06,1.89E-11,1.45E-05,1.99E-09,2.49E-07,1.50E-13,1.81E-12,8.13E-10,3.76E-13,2.55E-12,7.82E-11,2.68E-17,0.853497829,3.10E-08,1.25E-10,0.047771777,8.81E-11,0.000635231,1.20E-09,1.68E-10,3.89E-09,0.000299226,0.000306016,0.097437604,1.28E-05,6.85E-13,2.35E-06,8.20E-09,1.16E-09,1.49E-13,2.45E-06,1.22E-10,1.82E-05,1.76E-09,3.65E-07,1.33E-12,2.80E-12,1.07E-09,6.64E-13,1.18E-11,1.13E-10,4.92E-17,0.853497829,0.127703648,0.208787433,0.661823199,1.09E-08,8.06E-08,0.001685458,1.14E-09,5.41E-10,1.68E-07,1.54E-15,1.32E-10,4.29E-08,4.84E-10,0.047769448,2.39E-10,0.000636345,1.04E-09,3.00E-10,4.21E-08,0.000296882,0.000310359,0.097437468,1.03E-05,4.44E-13,1.92E-06,1.08E-08,5.41E-09,1.09E-13,2.78E-06,2.59E-11,2.06E-05,2.68E-09,4.08E-07,2.43E-13,3.41E-12,1.38E-09,5.26E-13,4.11E-12,1.42E-10,5.35E-17,0.853497829,4.59E-08,1.38E-10,0.047512558,9.61E-11,0.00080306,1.28E-09,2.72E-10,4.89E-09,0.000370472,0.000397811,0.097355557,1.33E-05,7.47E-13,2.34E-06,1.18E-08,1.24E-09,1.67E-13,3.15E-06,1.59E-10,2.45E-05,2.37E-09,5.63E-07,2.12E-12,4.99E-12,1.72E-09,8.38E-13,1.85E-11,1.96E-10,9.00E-17,0.853497829,0.128528458,0.206425064,0.663363033,1.37E-08,8.20E-08,0.001683133,1.17E-09,5.83E-10,2.14E-07,2.38E-15,1.99E-10,5.85E-08,5.11E-10,0.047510019,2.52E-10,0.000804338,1.10E-09,4.25E-10,4.80E-08,0.000368286,0.000402011,0.097355102,1.08E-05,4.97E-13,1.93E-06,1.46E-08,5.63E-09,1.23E-13,3.50E-06,3.29E-11,2.71E-05,3.34E-09,6.11E-07,3.66E-13,5.79E-12,2.10E-09,6.74E-13,6.14E-12,2.35E-10,9.41E-17,0.853497829,6.17E-08,1.48E-10,0.047247132,1.02E-10,0.000974561,1.36E-09,3.85E-10,5.58E-09,0.000440672,0.000494917,0.097271807,1.37E-05,8.11E-13,2.32E-06,1.57E-08,1.27E-09,1.84E-13,3.83E-06,1.97E-10,3.12E-05,3.02E-09,8.09E-07,3.14E-12,8.16E-12,2.53E-09,1.01E-12,2.72E-11,3.12E-10,1.51E-16,0.853497829,0.129644745,0.204448773,0.664230115,1.65E-08,8.33E-08,0.001676005,1.20E-09,6.25E-10,2.61E-07,3.43E-15,2.81E-10,7.47E-08,5.24E-10,0.047244216,2.59E-10,0.000976105,1.15E-09,5.54E-10,5.17E-08,0.000438726,0.000499047,0.097271096,1.11E-05,5.43E-13,1.92E-06,1.85E-08,5.68E-09,1.36E-13,4.19E-06,4.01E-11,3.39E-05,3.99E-09,8.62E-07,5.24E-13,9.20E-12,2.98E-09,8.19E-13,8.70E-12,3.62E-10,1.53E-16,0.853497829,7.80E-08,1.56E-10,0.046975096,1.08E-10,0.00115001,1.44E-09,5.00E-10,6.03E-09,0.000510074,0.000597327,0.097186217,1.41E-05,8.79E-13,2.28E-06,1.96E-08,1.29E-09,2.02E-13,4.47E-06,2.35E-10,3.81E-05,3.72E-09,1.11E-06,4.45E-12,1.26E-11,3.49E-09,1.19E-12,3.82E-11,4.68E-10,2.38E-16,0.853497829,0.130980271,0.202724131,0.664628445,1.94E-08,8.46E-08,0.001666739,1.23E-09,6.65E-10,3.09E-07,4.72E-15,3.76E-10,9.13E-08,5.30E-10,0.0469717,2.63E-10,0.00115188,1.19E-09,6.82E-10,5.42E-08,0.0005084,0.000601442,0.097185274,1.15E-05,5.89E-13,1.89E-06,2.24E-08,5.65E-09,1.47E-13,4.84E-06,4.72E-11,4.10E-05,4.67E-09,1.16E-06,7.24E-13,1.39E-11,4.01E-09,9.64E-13,1.19E-11,5.30E-10,2.36E-16,0.853497829,9.48E-08,1.64E-10,0.046696024,1.13E-10,0.001329691,1.52E-09,6.14E-10,6.33E-09,0.000578912,0.000705058,0.097098635,1.45E-05,9.53E-13,2.25E-06,2.36E-08,1.29E-09,2.23E-13,5.09E-06,2.72E-10,4.53E-05,4.47E-09,1.46E-06,6.09E-12,1.85E-11,4.62E-09,1.38E-12,5.17E-11,6.72E-10,3.61E-16,0.853497829,0.132504647,0.201185126,0.664653242,2.24E-08,8.60E-08,0.001656516,1.26E-09,7.05E-10,3.59E-07,6.25E-15,4.86E-10,1.08E-07,5.33E-10,0.046692085,2.66E-10,0.001331925,1.24E-09,8.07E-10,5.58E-08,0.000577518,0.000709202,0.097097466,1.18E-05,6.36E-13,1.86E-06,2.64E-08,5.59E-09,1.59E-13,5.46E-06,5.44E-11,4.83E-05,5.38E-09,1.52E-06,9.74E-13,2.02E-11,5.21E-09,1.11E-12,1.58E-11,7.50E-10,3.50E-16,0.853497829,1.12E-07,1.72E-10,0.046409479,1.19E-10,0.001513899,1.61E-09,7.25E-10,6.55E-09,0.000647406,0.000818151,0.097008898,1.49E-05,1.03E-12,2.22E-06,2.77E-08,1.29E-09,2.45E-13,5.67E-06,3.10E-10,5.29E-05,5.28E-09,1.87E-06,8.14E-12,2.65E-11,5.91E-09,1.57E-12,6.82E-11,9.34E-10,5.29E-16,0.853497829,0.134205393,0.199795615,0.664352624,2.55E-08,8.74E-08,0.001645842,1.29E-09,7.45E-10,4.10E-07,8.06E-15,6.11E-10,1.26E-07,5.35E-10,0.046404946,2.68E-10,0.001516527,1.28E-09,9.28E-10,5.71E-08,0.000646291,0.000822363,0.097007498,1.21E-05,6.88E-13,1.83E-06,3.05E-08,5.52E-09,1.71E-13,6.05E-06,6.18E-11,5.60E-05,6.14E-09,1.94E-06,1.29E-12,2.86E-11,6.57E-09,1.26E-12,2.05E-11,1.03E-09,5.05E-16,0.853497829,1.30E-07,1.80E-10,0.046115003,1.25E-10,0.001702938,1.69E-09,8.32E-10,6.73E-09,0.000715767,0.000936668,0.096916833,1.54E-05,1.13E-12,2.18E-06,3.18E-08,1.29E-09,2.68E-13,6.22E-06,3.47E-10,6.08E-05,6.16E-09,2.35E-06,1.07E-11,3.68E-11,7.37E-09,1.77E-12,8.81E-11,1.27E-09,7.58E-16,0.853497829,0.136077909,0.198532816,0.663753768,2.87E-08,8.89E-08,0.001634924,1.33E-09,7.85E-10,4.63E-07,1.02E-14,7.51E-10,1.44E-07,5.37E-10,0.046109836,2.70E-10,0.001705983,1.32E-09,1.04E-09,5.81E-08,0.000714926,0.00094098,0.096915197,1.25E-05,7.44E-13,1.79E-06,3.47E-08,5.43E-09,1.84E-13,6.61E-06,6.92E-11,6.40E-05,6.96E-09,2.43E-06,1.67E-12,3.95E-11,8.11E-09,1.41E-12,2.61E-11,1.38E-09,7.12E-16,0.853497829,1.48E-07,1.89E-10,0.045812124,1.32E-10,0.001897119,1.78E-09,9.34E-10,6.88E-09,0.000784197,0.001060687,0.096822265,1.59E-05,1.23E-12,2.14E-06,3.60E-08,1.29E-09,2.94E-13,6.73E-06,3.84E-10,6.90E-05,7.12E-09,2.90E-06,1.37E-11,5.02E-11,9.01E-09,1.97E-12,1.12E-10,1.68E-09,1.06E-15,0.853497829,0.138119627,0.197379565,0.662876317,3.20E-08,9.04E-08,0.001623849,1.36E-09,8.26E-10,5.17E-07,1.26E-14,9.08E-10,1.62E-07,5.39E-10,0.045806318,2.72E-10,0.001900582,1.36E-09,1.16E-09,5.90E-08,0.000783612,0.001065124,0.096820395,1.28E-05,8.06E-13,1.76E-06,3.89E-08,5.35E-09,1.98E-13,7.13E-06,7.68E-11,7.24E-05,7.84E-09,2.99E-06,2.13E-12,5.35E-11,9.82E-09,1.57E-12,3.29E-11,1.82E-09,9.86E-16,0.853497829,1.66E-07,1.98E-10,0.045500419,1.39E-10,0.002096722,1.89E-09,1.03E-09,7.02E-09,0.00085287,0.001190288,0.096725032,1.64E-05,1.34E-12,2.11E-06,4.03E-08,1.29E-09,3.25E-13,7.21E-06,4.22E-10,7.75E-05,8.15E-09,3.53E-06,1.75E-11,6.72E-11,1.08E-08,2.19E-12,1.41E-10,2.19E-09,1.47E-15,0.853497829,0.140321455,0.196319104,0.661746085,3.55E-08,9.20E-08,0.001612653,1.40E-09,8.66E-10,5.72E-07,1.55E-14,1.08E-09,1.81E-07,5.41E-10,0.045494097,2.74E-10,0.002100513,1.41E-09,1.26E-09,6.00E-08,0.000852483,0.001194836,0.096722974,1.32E-05,8.75E-13,1.72E-06,4.32E-08,5.27E-09,2.13E-13,7.62E-06,8.45E-11,8.11E-05,8.79E-09,3.62E-06,2.70E-12,7.12E-11,1.17E-08,1.73E-12,4.09E-11,2.36E-09,1.35E-15,0.853497829,1.85E-07,2.08E-10,0.045179733,1.47E-10,0.002301833,1.98E-09,1.13E-09,7.17E-09,0.000921859,0.001325487,0.09662507,1.69E-05,1.47E-12,2.07E-06,4.47E-08,1.29E-09,3.54E-13,7.67E-06,4.59E-10,8.65E-05,9.26E-09,4.25E-06,2.21E-11,8.87E-11,1.28E-08,2.42E-12,1.74E-10,2.81E-09,2.00E-15,0.853497829,0.142639365,0.195324911,0.660433542,3.91E-08,9.36E-08,0.001601416,1.44E-09,9.07E-10,6.29E-07,1.87E-14,1.27E-09,2.00E-07,5.44E-10,0.045172453,2.77E-10,0.002306276,1.46E-09,1.37E-09,6.09E-08,0.000921805,0.001330279,0.096622668,1.36E-05,9.51E-13,1.68E-06,4.76E-08,5.18E-09,2.30E-13,8.08E-06,9.24E-11,9.03E-05,9.81E-09,4.35E-06,3.40E-12,9.37E-11,1.38E-08,1.90E-12,5.04E-11,3.02E-09,1.81E-15,0.853497829,2.05E-07,2.19E-10,0.04484893,1.55E-10,0.002513261,2.08E-09,1.22E-09,7.31E-09,0.000991566,0.001466566,0.096521963,1.75E-05,1.61E-12,2.03E-06,4.92E-08,1.28E-09,3.89E-13,8.09E-06,4.96E-10,9.59E-05,1.05E-08,5.06E-06,2.76E-11,1.16E-10,1.50E-08,2.66E-12,2.15E-10,3.57E-09,2.69E-15,0.853497829,0.145188956,0.194423754,0.658796566,4.29E-08,9.53E-08,0.001589895,1.48E-09,9.49E-10,6.88E-07,2.25E-14,1.48E-09,2.20E-07,5.46E-10,0.044840671,2.80E-10,0.002518345,1.52E-09,1.47E-09,6.19E-08,0.00099183,0.001471631,0.096519216,1.40E-05,1.04E-12,1.65E-06,5.20E-08,5.08E-09,2.48E-13,8.51E-06,1.01E-10,9.98E-05,1.09E-08,5.17E-06,4.27E-12,1.22E-10,1.61E-08,2.08E-12,6.17E-11,3.81E-09,2.42E-15,0.853497829,2.25E-07,2.31E-10,0.044507302,1.65E-10,0.002731451,2.19E-09,1.30E-09,7.47E-09,0.001062225,0.001613716,0.096415457,1.81E-05,1.77E-12,2.00E-06,5.37E-08,1.28E-09,4.29E-13,8.49E-06,5.34E-10,0.000105712,1.18E-08,5.97E-06,3.44E-11,1.50E-10,1.75E-08,2.91E-12,2.62E-10,4.49E-09,3.60E-15,0.853497829,0.147965365,0.19360465,0.656850947,4.68E-08,9.72E-08,0.001578142,1.53E-09,9.91E-10,7.48E-07,2.67E-14,1.71E-09,2.40E-07,5.50E-10,0.04449808,2.83E-10,0.002737169,1.57E-09,1.57E-09,6.30E-08,0.001062784,0.001619069,0.096412371,1.44E-05,1.13E-12,1.61E-06,5.66E-08,4.99E-09,2.68E-13,8.91E-06,1.09E-10,0.0001099,1.21E-08,6.10E-06,5.32E-12,1.57E-10,1.86E-08,2.26E-12,7.50E-11,4.78E-09,3.22E-15,0.853497829,2.46E-07,2.44E-10,0.044154224,1.76E-10,0.002956825,2.31E-09,1.39E-09,7.64E-09,0.001134036,0.001767112,0.096305321,1.87E-05,1.95E-12,1.96E-06,5.84E-08,1.27E-09,4.74E-13,8.85E-06,5.72E-10,0.000116048,1.33E-08,7.01E-06,4.25E-11,1.92E-10,2.01E-08,3.19E-12,3.19E-10,5.59E-09,4.78E-15,0.853497829,0.150950691,0.192848479,0.654633689,5.10E-08,9.91E-08,0.001566178,1.57E-09,1.03E-09,8.09E-07,3.16E-14,1.96E-09,2.61E-07,5.55E-10,0.044143973,2.87E-10,0.002963222,1.64E-09,1.67E-09,6.41E-08,0.001134897,0.001772789,0.096301873,1.49E-05,1.24E-12,1.57E-06,6.13E-08,4.90E-09,2.91E-13,9.29E-06,1.18E-10,0.000120465,1.34E-08,7.15E-06,6.60E-12,2.02E-10,2.13E-08,2.46E-12,9.06E-11,5.93E-09,4.24E-15,0.853497829,2.67E-07,2.58E-10,0.043788956,1.88E-10,0.003189882,2.44E-09,1.47E-09,7.82E-09,0.001207232,0.001926969,0.096191282,1.94E-05,2.17E-12,1.92E-06,6.31E-08,1.27E-09,5.25E-13,9.19E-06,6.11E-10,0.000126922,1.49E-08,8.17E-06,5.25E-11,2.45E-10,2.29E-08,3.48E-12,3.85E-10,6.91E-09,6.32E-15,0.853497829,0.154152649,0.192146513,0.652145833,5.53E-08,1.01E-07,0.001553971,1.63E-09,1.08E-09,8.72E-07,3.71E-14,2.24E-09,2.83E-07,5.60E-10,0.043777586,2.92E-10,0.003197018,1.70E-09,1.76E-09,6.54E-08,0.001208408,0.001933009,0.096187441,1.53E-05,1.37E-12,1.54E-06,6.61E-08,4.80E-09,3.16E-13,9.64E-06,1.27E-10,0.000131587,1.49E-08,8.34E-06,8.16E-12,2.57E-10,2.43E-08,2.66E-12,1.09E-10,7.32E-09,5.55E-15,0.853497829,2.90E-07,2.74E-10,0.043410682,2.01E-10,0.003431174,2.57E-09,1.55E-09,8.01E-09,0.001282062,0.002093524,0.09607304,2.01E-05,2.41E-12,1.88E-06,6.80E-08,1.26E-09,5.80E-13,9.51E-06,6.50E-10,0.000138379,1.66E-08,9.48E-06,6.45E-11,3.11E-10,2.60E-08,3.80E-12,4.63E-10,8.49E-09,8.31E-15,0.853497829,0.15758373,0.191491672,0.649381996,6.00E-08,1.03E-07,0.001541497,1.68E-09,1.13E-09,9.37E-07,4.34E-14,2.53E-09,3.06E-07,5.66E-10,0.043398095,2.97E-10,0.003439117,1.78E-09,1.86E-09,6.67E-08,0.001283571,0.00209997,0.09606877,1.59E-05,1.51E-12,1.50E-06,7.10E-08,4.71E-09,3.44E-13,9.96E-06,1.36E-10,0.000143315,1.64E-08,9.67E-06,1.01E-11,3.26E-10,2.75E-08,2.88E-12,1.30E-10,8.97E-09,7.25E-15,0.853497829,3.13E-07,2.92E-10,0.043018521,2.16E-10,0.0036813,2.72E-09,1.63E-09,8.21E-09,0.001358791,0.002267042,0.095950268,2.09E-05,2.68E-12,1.84E-06,7.31E-08,1.26E-09,6.45E-13,9.80E-06,6.90E-10,0.000150471,1.85E-08,1.10E-05,7.92E-11,3.92E-10,2.94E-08,4.14E-12,5.55E-10,1.04E-08,1.09E-14,0.853497829,0.161258751,0.190877424,0.646333929,6.48E-08,1.05E-07,0.001528716,1.74E-09,1.17E-09,1.00E-06,5.06E-14,2.85E-09,3.30E-07,5.74E-10,0.043004602,3.03E-10,0.003690128,1.85E-09,1.95E-09,6.82E-08,0.001360657,0.00227394,0.09594553,1.64E-05,1.67E-12,1.46E-06,7.61E-08,4.61E-09,3.75E-13,1.03E-05,1.46E-10,0.0001557,1.81E-08,1.12E-05,1.24E-11,4.11E-10,3.10E-08,3.12E-12,1.56E-10,1.09E-08,9.44E-15,0.853497829,3.38E-07,3.12E-10,0.04261151,2.33E-10,0.003940917,2.88E-09,1.70E-09,8.43E-09,0.001437706,0.002447813,0.09582261,2.17E-05,3.00E-12,1.80E-06,7.82E-08,1.25E-09,7.18E-13,1.01E-05,7.31E-10,0.000163254,2.06E-08,1.26E-05,9.71E-11,4.94E-10,3.30E-08,4.52E-12,6.63E-10,1.26E-08,1.42E-14,0.853497829,0.165194332,0.190297049,0.642991749,7.00E-08,1.08E-07,0.001515614,1.80E-09,1.22E-09,1.07E-06,5.87E-14,3.20E-09,3.55E-07,5.83E-10,0.042596137,3.09E-10,0.003950713,1.94E-09,2.05E-09,6.97E-08,0.001439953,0.002455212,0.095817357,1.70E-05,1.86E-12,1.42E-06,8.13E-08,4.52E-09,4.11E-13,1.05E-05,1.56E-10,0.000168804,2.00E-08,1.28E-05,1.53E-11,5.18E-10,3.47E-08,3.37E-12,1.85E-10,1.33E-08,1.23E-14,0.853497829,3.64E-07,3.34E-10,0.042188607,2.52E-10,0.004210737,3.05E-09,1.78E-09,8.67E-09,0.001519111,0.002636154,0.095689675,2.26E-05,3.37E-12,1.76E-06,8.36E-08,1.25E-09,8.02E-13,1.03E-05,7.72E-10,0.00017679,2.30E-08,1.45E-05,1.19E-10,6.21E-10,3.69E-08,4.93E-12,7.91E-10,1.53E-08,1.86E-14,0.853497829,0.169409547,0.189744217,0.639342761,7.55E-08,1.10E-07,0.001502142,1.86E-09,1.27E-09,1.14E-06,6.79E-14,3.57E-09,3.82E-07,5.93E-10,0.04217164,3.17E-10,0.004221599,2.03E-09,2.15E-09,7.14E-08,0.001521772,0.002644105,0.095683856,1.76E-05,2.07E-12,1.38E-06,8.67E-08,4.42E-09,4.51E-13,1.08E-05,1.67E-10,0.00018269,2.20E-08,1.47E-05,1.88E-11,6.50E-10,3.87E-08,3.63E-12,2.20E-10,1.61E-08,1.59E-14,0.853497829,3.91E-07,3.59E-10,0.041748676,2.74E-10,0.004491544,3.24E-09,1.86E-09,8.93E-09,0.001603341,0.002832409,0.095551037,2.36E-05,3.80E-12,1.71E-06,8.91E-08,1.24E-09,8.98E-13,1.05E-05,8.14E-10,0.000191149,2.55E-08,1.65E-05,1.46E-10,7.79E-10,4.11E-08,5.38E-12,9.41E-10,1.84E-08,2.42E-14,0.853497829,0.173925623,0.189212164,0.635372525,8.13E-08,1.13E-07,0.001488275,1.93E-09,1.32E-09,1.21E-06,7.85E-14,3.97E-09,4.10E-07,6.05E-10,0.041729959,3.26E-10,0.004503578,2.13E-09,2.25E-09,7.32E-08,0.001606451,0.002840968,0.095544593,1.83E-05,2.31E-12,1.34E-06,9.22E-08,4.32E-09,4.98E-13,1.10E-05,1.79E-10,0.000197433,2.42E-08,1.69E-05,2.32E-11,8.16E-10,4.31E-08,3.92E-12,2.60E-10,1.93E-08,2.06E-14,0.853497829,4.20E-07,3.88E-10,0.041290475,2.99E-10,0.004784194,3.44E-09,1.94E-09,9.20E-09,0.001690758,0.003036951,0.095406225,2.47E-05,4.30E-12,1.67E-06,9.49E-08,1.24E-09,1.01E-12,1.07E-05,8.58E-10,0.000206411,2.84E-08,1.89E-05,1.78E-10,9.76E-10,4.57E-08,5.89E-12,1.12E-09,2.21E-08,3.15E-14,0.853497829,0.178766506,0.188693984,0.631064041,8.75E-08,1.16E-07,0.001473973,2.01E-09,1.38E-09,1.28E-06,9.04E-14,4.41E-09,4.39E-07,6.20E-10,0.041269835,3.36E-10,0.00479752,2.24E-09,2.35E-09,7.51E-08,0.001694361,0.003046178,0.09539909,1.90E-05,2.60E-12,1.30E-06,9.80E-08,4.21E-09,5.50E-13,1.12E-05,1.91E-10,0.000213117,2.66E-08,1.92E-05,2.87E-11,1.02E-09,4.78E-08,4.23E-12,3.07E-10,2.32E-08,2.67E-14,0.853497829,4.51E-07,4.20E-10,0.040812652,3.27E-10,0.005089626,3.67E-09,2.02E-09,9.50E-09,0.00178176,0.003250183,0.095254721,2.58E-05,4.88E-12,1.63E-06,1.01E-07,1.23E-09,1.14E-12,1.09E-05,9.03E-10,0.000222663,3.16E-08,2.15E-05,2.19E-10,1.22E-09,5.05E-08,6.45E-12,1.32E-09,2.65E-08,4.11E-14,0.853497829,0.183958956,0.188182282,0.626397979,9.41E-08,1.19E-07,0.001459201,2.09E-09,1.43E-09,1.36E-06,1.04E-13,4.88E-09,4.71E-07,6.36E-10,0.040789895,3.47E-10,0.00510438,2.36E-09,2.46E-09,7.72E-08,0.001785906,0.003260143,0.095246821,1.98E-05,2.93E-12,1.26E-06,1.04E-07,4.11E-09,6.12E-13,1.14E-05,2.04E-10,0.000229832,2.93E-08,2.19E-05,3.55E-11,1.28E-09,5.28E-08,4.58E-12,3.62E-10,2.79E-08,3.46E-14,0.853497829,4.84E-07,4.57E-10,0.040313727,3.60E-10,0.005408872,3.92E-09,2.10E-09,9.82E-09,0.001876785,0.00347254,0.095095956,2.70E-05,5.57E-12,1.58E-06,1.07E-07,1.22E-09,1.30E-12,1.10E-05,9.49E-10,0.000240004,3.51E-08,2.45E-05,2.69E-10,1.54E-09,5.58E-08,7.08E-12,1.57E-09,3.18E-08,5.37E-14,0.853497829,0.189532929,0.187669125,0.621352361,1.01E-07,1.22E-07,0.001443918,2.18E-09,1.49E-09,1.44E-06,1.20E-13,5.38E-09,5.05E-07,6.56E-10,0.040288639,3.60E-10,0.005425204,2.49E-09,2.58E-09,7.94E-08,0.00188153,0.0034833,0.095087205,2.07E-05,3.31E-12,1.22E-06,1.10E-07,4.00E-09,6.82E-13,1.16E-05,2.18E-10,0.000247686,3.22E-08,2.50E-05,4.41E-11,1.61E-09,5.83E-08,4.95E-12,4.27E-10,3.34E-08,4.49E-14,0.853497829,5.19E-07,4.99E-10,0.039792087,3.97E-10,0.005743064,4.19E-09,2.19E-09,1.02E-08,0.001976314,0.003704486,0.094929298,2.83E-05,6.40E-12,1.54E-06,1.14E-07,1.21E-09,1.47E-12,1.12E-05,9.97E-10,0.000258547,3.91E-08,2.79E-05,3.31E-10,1.93E-09,6.15E-08,7.80E-12,1.86E-09,3.80E-08,7.02E-14,0.853497829,0.195521814,0.187145779,0.615902574,1.09E-07,1.25E-07,0.001428076,2.27E-09,1.56E-09,1.51E-06,1.38E-13,5.93E-09,5.42E-07,6.79E-10,0.039764428,3.76E-10,0.005761141,2.63E-09,2.70E-09,8.18E-08,0.001981726,0.00371612,0.094919601,2.16E-05,3.77E-12,1.18E-06,1.17E-07,3.89E-09,7.65E-13,1.17E-05,2.33E-10,0.000266795,3.54E-08,2.84E-05,5.48E-11,2.03E-09,6.41E-08,5.36E-12,5.02E-10,3.98E-08,5.84E-14,0.853497829,5.58E-07,5.48E-10,0.039245965,4.41E-10,0.006093445,4.50E-09,2.28E-09,1.05E-08,0.002080884,0.003946517,0.094754049,2.97E-05,7.38E-12,1.49E-06,1.21E-07,1.20E-09,1.69E-12,1.13E-05,1.05E-09,0.000278419,4.36E-08,3.17E-05,4.08E-10,2.43E-09,6.76E-08,8.61E-12,2.20E-09,4.53E-08,9.21E-14,0.853497829,0.201962795,0.18660252,0.610021201,1.17E-07,1.28E-07,0.001411634,2.37E-09,1.62E-09,1.59E-06,1.58E-13,6.51E-09,5.82E-07,7.07E-10,0.039215469,3.94E-10,0.006113457,2.79E-09,2.82E-09,8.44E-08,0.002087038,0.003959101,0.094743299,2.25E-05,4.31E-12,1.14E-06,1.24E-07,3.78E-09,8.62E-13,1.18E-05,2.49E-10,0.000287295,3.89E-08,3.23E-05,6.82E-11,2.56E-09,7.04E-08,5.81E-12,5.90E-10,4.76E-08,7.62E-14,0.853497829,6.00E-07,6.04E-10,0.038673435,4.91E-10,0.006461381,4.84E-09,2.37E-09,1.09E-08,0.002191085,0.004199157,0.094569438,3.13E-05,8.55E-12,1.44E-06,1.28E-07,1.19E-09,1.94E-12,1.14E-05,1.10E-09,0.000299762,4.86E-08,3.59E-05,5.04E-10,3.07E-09,7.42E-08,9.54E-12,2.60E-09,5.41E-08,1.21E-13,0.853497829,0.208897208,0.186028475,0.603677836,1.26E-07,1.32E-07,0.001394538,2.48E-09,1.69E-09,1.67E-06,1.82E-13,7.15E-09,6.27E-07,7.39E-10,0.038639807,4.14E-10,0.006483539,2.97E-09,2.96E-09,8.72E-08,0.002198069,0.004212772,0.094557511,2.36E-05,4.95E-12,1.09E-06,1.31E-07,3.68E-09,9.76E-13,1.19E-05,2.66E-10,0.000309336,4.28E-08,3.66E-05,8.51E-11,3.24E-09,7.72E-08,6.31E-12,6.92E-10,5.67E-08,9.96E-14,0.853497829,6.47E-07,6.70E-10,0.038072395,5.50E-10,0.006848372,5.22E-09,2.47E-09,1.14E-08,0.002307577,0.004462958,0.094374608,3.30E-05,9.97E-12,1.39E-06,1.36E-07,1.18E-09,2.25E-12,1.15E-05,1.16E-09,0.000322741,5.43E-08,4.08E-05,6.26E-10,3.89E-09,8.13E-08,1.06E-11,3.07E-09,6.45E-08,1.60E-13,0.853497829,0.216370832,0.185411305,0.596839093,1.35E-07,1.36E-07,0.001376731,2.60E-09,1.76E-09,1.76E-06,2.09E-13,7.84E-09,6.75E-07,7.78E-10,0.038035311,4.38E-10,0.00687291,3.17E-09,3.11E-09,9.01E-08,0.00231549,0.004477687,0.094361366,2.47E-05,5.72E-12,1.05E-06,1.39E-07,3.56E-09,1.11E-12,1.20E-05,2.85E-10,0.000333092,4.71E-08,4.16E-05,1.07E-10,4.11E-09,8.46E-08,6.87E-12,8.12E-10,6.76E-08,1.31E-13,0.853497829,6.99E-07,7.48E-10,0.03744056,6.20E-10,0.007256064,5.66E-09,2.58E-09,1.18E-08,0.002431092,0.004738492,0.09416861,3.48E-05,1.17E-11,1.34E-06,1.44E-07,1.17E-09,2.62E-12,1.15E-05,1.21E-09,0.00034754,6.07E-08,4.64E-05,7.80E-10,4.94E-09,8.90E-08,1.19E-11,3.62E-09,7.68E-08,2.12E-13,0.853497829,0.224434156,0.184736887,0.589468659,1.45E-07,1.40E-07,0.001358159,2.73E-09,1.84E-09,1.84E-06,2.40E-13,8.58E-09,7.30E-07,8.23E-10,0.037399662,4.66E-10,0.007283241,3.39E-09,3.27E-09,9.33E-08,0.002440048,0.00475442,0.094153898,2.60E-05,6.65E-12,1.00E-06,1.48E-07,3.45E-09,1.28E-12,1.21E-05,3.06E-10,0.00035876,5.19E-08,4.73E-05,1.35E-10,5.25E-09,9.25E-08,7.50E-12,9.51E-10,8.06E-08,1.73E-13,0.853497829,7.57E-07,8.39E-10,0.036775447,7.01E-10,0.007686254,6.15E-09,2.69E-09,1.24E-08,0.002562443,0.005026347,0.093950396,3.69E-05,1.38E-11,1.29E-06,1.53E-07,1.16E-09,3.08E-12,1.16E-05,1.28E-09,0.00037437,6.80E-08,5.27E-05,9.77E-10,6.32E-09,9.73E-08,1.33E-11,4.28E-09,9.16E-08,2.82E-13,0.853497829,0.233142692,0.183989131,0.581527181,1.56E-07,1.45E-07,0.001338756,2.87E-09,1.92E-09,1.93E-06,2.76E-13,9.38E-09,7.91E-07,8.77E-10,0.036730349,4.99E-10,0.007716357,3.64E-09,3.44E-09,9.67E-08,0.00257257,0.005043559,0.09393404,2.73E-05,7.77E-12,9.60E-07,1.57E-07,3.33E-09,1.47E-12,1.22E-05,3.29E-10,0.000386562,5.71E-08,5.38E-05,1.71E-10,6.73E-09,1.01E-07,8.21E-12,1.11E-09,9.61E-08,2.29E-13,0.853497829,8.22E-07,9.48E-10,0.036074381,7.98E-10,0.008140903,6.71E-09,2.81E-09,1.29E-08,0.002702535,0.005327116,0.093718806,3.91E-05,1.64E-11,1.24E-06,1.62E-07,1.15E-09,3.64E-12,1.16E-05,1.34E-09,0.000403471,7.63E-08,5.99E-05,1.23E-09,8.13E-09,1.06E-07,1.50E-11,5.05E-09,1.09E-07,3.78E-13,0.853497829,0.242557038,0.183149582,0.572972578,1.68E-07,1.49E-07,0.001318458,3.02E-09,2.00E-09,2.01E-06,3.18E-13,1.02E-08,8.61E-07,9.42E-10,0.036024663,5.37E-10,0.008174245,3.93E-09,3.62E-09,1.00E-07,0.002713976,0.005345692,0.093700611,2.88E-05,9.15E-12,9.14E-07,1.67E-07,3.22E-09,1.71E-12,1.22E-05,3.54E-10,0.000416752,6.30E-08,6.12E-05,2.17E-10,8.69E-09,1.10E-07,9.01E-12,1.30E-09,1.15E-07,3.05E-13,0.853497829,8.97E-07,1.08E-09,0.035334487,9.13E-10,0.00862214,7.36E-09,2.94E-09,1.35E-08,0.002852366,0.005641381,0.093472562,4.15E-05,1.97E-11,1.19E-06,1.73E-07,1.14E-09,4.34E-12,1.16E-05,1.42E-09,0.000435115,8.59E-08,6.82E-05,1.55E-09,1.05E-08,1.16E-07,1.70E-11,5.97E-09,1.31E-07,5.10E-13,0.853497829,0.252742766,0.182197081,0.563760508,1.81E-07,1.54E-07,0.001297194,3.18E-09,2.09E-09,2.10E-06,3.68E-13,1.12E-08,9.40E-07,1.02E-09,0.035279706,5.83E-10,0.008659058,4.24E-09,3.83E-09,1.04E-07,0.002865283,0.005661394,0.093452314,3.04E-05,1.08E-11,8.68E-07,1.77E-07,3.10E-09,2.01E-12,1.23E-05,3.82E-10,0.00044962,6.95E-08,6.97E-05,2.77E-10,1.13E-08,1.21E-07,9.93E-12,1.52E-09,1.37E-07,4.10E-13,0.853497829,9.84E-07,1.23E-09,0.034552711,1.05E-09,0.009132257,8.11E-09,3.07E-09,1.41E-08,0.003013039,0.005969694,0.093210261,4.42E-05,2.38E-11,1.13E-06,1.84E-07,1.13E-09,5.21E-12,1.16E-05,1.49E-09,0.000469612,9.68E-08,7.77E-05,1.97E-09,1.37E-08,1.27E-07,1.94E-11,7.04E-09,1.56E-07,6.93E-13,0.853497829,0.263770092,0.181107494,0.553844964,1.95E-07,1.60E-07,0.001274892,3.36E-09,2.18E-09,2.19E-06,4.25E-13,1.22E-08,1.03E-06,1.11E-09,0.034492405,6.36E-10,0.00917311,4.61E-09,4.05E-09,1.08E-07,0.003027605,0.005991205,0.093187725,3.21E-05,1.30E-11,8.22E-07,1.89E-07,2.98E-09,2.37E-12,1.23E-05,4.13E-10,0.000485492,7.67E-08,7.94E-05,3.57E-10,1.48E-08,1.31E-07,1.10E-11,1.77E-09,1.64E-07,5.54E-13,0.853497829,1.08E-06,1.42E-09,0.033725847,1.21E-09,0.0096737,8.97E-09,3.22E-09,1.48E-08,0.003185757,0.006312551,0.092930377,4.72E-05,2.89E-11,1.08E-06,1.96E-07,1.11E-09,6.32E-12,1.17E-05,1.58E-09,0.00050731,1.09E-07,8.86E-05,2.52E-09,1.81E-08,1.38E-07,2.23E-11,8.31E-09,1.87E-07,9.50E-13,0.853497829,0.275713013,0.179853373,0.54317947,2.10E-07,1.65E-07,0.001251477,3.56E-09,2.28E-09,2.27E-06,4.93E-13,1.33E-08,1.14E-06,1.22E-09,0.033659547,6.99E-10,0.009718861,5.02E-09,4.29E-09,1.12E-07,0.003202161,0.0063356,0.092905297,3.41E-05,1.56E-11,7.75E-07,2.02E-07,2.87E-09,2.82E-12,1.24E-05,4.49E-10,0.000524737,8.47E-08,9.07E-05,4.63E-10,1.95E-08,1.43E-07,1.22E-11,2.07E-09,1.96E-07,7.55E-13,0.853497829,1.20E-06,1.66E-09,0.032850579,1.41E-09,0.010249044,9.98E-09,3.38E-09,1.55E-08,0.003371825,0.006670355,0.092631256,5.05E-05,3.54E-11,1.02E-06,2.10E-07,1.10E-09,7.73E-12,1.17E-05,1.67E-09,0.000548604,1.24E-07,0.000101306,3.24E-09,2.40E-08,1.51E-07,2.58E-11,9.79E-09,2.24E-07,1.31E-12,0.853497829,0.288647968,0.17840383,0.531718547,2.27E-07,1.72E-07,0.001226877,3.78E-09,2.39E-09,2.36E-06,5.73E-13,1.45E-08,1.27E-06,1.36E-09,0.032777833,7.74E-10,0.01029889,5.50E-09,4.56E-09,1.17E-07,0.003390261,0.006694952,0.092603366,3.62E-05,1.89E-11,7.28E-07,2.16E-07,2.75E-09,3.39E-12,1.24E-05,4.90E-10,0.000567771,9.36E-08,0.00010382,6.04E-10,2.61E-08,1.56E-07,1.36E-11,2.41E-09,2.36E-07,1.04E-12,0.853497829,1.34E-06,1.94E-09,0.03192357,1.65E-09,0.010860953,1.12E-08,3.55E-09,1.63E-08,0.003572631,0.007043372,0.092311135,5.41E-05,4.38E-11,9.62E-07,2.25E-07,1.08E-09,9.56E-12,1.17E-05,1.78E-09,0.000593934,1.41E-07,0.000116021,4.20E-09,3.21E-08,1.64E-07,3.00E-11,1.15E-08,2.69E-07,1.83E-12,0.853497829,0.30265161,0.176724476,0.519420005,2.45E-07,1.78E-07,0.00120102,4.02E-09,2.50E-09,2.44E-06,6.67E-13,1.58E-08,1.42E-06,1.52E-09,0.031843963,8.62E-10,0.010915844,6.05E-09,4.86E-09,1.21E-07,0.003593298,0.007069491,0.092280158,3.85E-05,2.31E-11,6.81E-07,2.32E-07,2.63E-09,4.12E-12,1.25E-05,5.36E-10,0.000615053,1.04E-07,0.000119015,7.93E-10,3.51E-08,1.71E-07,1.52E-11,2.80E-09,2.84E-07,1.44E-12,0.853497829,1.51E-06,2.29E-09,0.030941564,1.94E-09,0.011512107,1.25E-08,3.74E-09,1.71E-08,0.003789628,0.00743168,0.091968165,5.81E-05,5.47E-11,9.03E-07,2.43E-07,1.06E-09,1.19E-11,1.17E-05,1.89E-09,0.000643786,1.61E-07,0.00013314,5.48E-09,4.35E-08,1.79E-07,3.51E-11,1.36E-08,3.25E-07,2.58E-12,0.853497829,0.317797652,0.17477765,0.506247853,2.64E-07,1.85E-07,0.001173843,4.28E-09,2.62E-09,2.53E-06,7.79E-13,1.73E-08,1.60E-06,1.72E-09,0.030854761,9.68E-10,0.011572366,6.69E-09,5.18E-09,1.26E-07,0.003812711,0.007459241,0.091933828,4.10E-05,2.85E-11,6.34E-07,2.49E-07,2.51E-09,5.05E-12,1.25E-05,5.90E-10,0.000667092,1.15E-07,0.000136714,1.05E-09,4.78E-08,1.86E-07,1.71E-11,3.25E-09,3.42E-07,2.02E-12,0.853497829,1.71E-06,2.72E-09,0.029901548,2.29E-09,0.012205108,1.41E-08,3.94E-09,1.80E-08,0.004024288,0.007835099,0.091600445,6.25E-05,6.88E-11,8.44E-07,2.62E-07,1.04E-09,1.50E-11,1.17E-05,2.02E-09,0.000698688,1.84E-07,0.000153108,7.18E-09,5.95E-08,1.96E-07,4.12E-11,1.59E-08,3.92E-07,3.67E-12,0.853497829,0.334152553,0.172523052,0.492175988,2.86E-07,1.92E-07,0.001145294,4.56E-09,2.74E-09,2.61E-06,9.13E-13,1.89E-08,1.82E-06,1.96E-09,0.029807339,1.09E-09,0.012270983,7.42E-09,5.53E-09,1.31E-07,0.004049944,0.007863959,0.091562495,4.38E-05,3.55E-11,5.87E-07,2.70E-07,2.39E-09,6.26E-12,1.25E-05,6.52E-10,0.000724433,1.27E-07,0.00015738,1.40E-09,6.58E-08,2.03E-07,1.93E-11,3.76E-09,4.14E-07,2.87E-12,0.853497829,1.96E-06,3.26E-09,0.028800965,2.71E-09,0.012942337,1.61E-08,4.15E-09,1.89E-08,0.004278039,0.008253123,0.091206091,6.74E-05,8.74E-11,7.84E-07,2.84E-07,1.02E-09,1.92E-11,1.17E-05,2.17E-09,0.000759207,2.11E-07,0.000176448,9.49E-09,8.24E-08,2.14E-07,4.87E-11,1.86E-08,4.75E-07,5.29E-12,0.853497829,0.351769884,0.169918824,0.477192731,3.09E-07,2.00E-07,0.001115339,4.87E-09,2.87E-09,2.69E-06,1.07E-12,2.06E-08,2.09E-06,2.25E-09,0.028699321,1.24E-09,0.013013966,8.28E-09,5.92E-09,1.36E-07,0.004306378,0.008283067,0.091164317,4.69E-05,4.46E-11,5.40E-07,2.93E-07,2.27E-09,7.83E-12,1.26E-05,7.25E-10,0.000787653,1.41E-07,0.000181555,1.88E-09,9.16E-08,2.22E-07,2.19E-11,4.35E-09,5.01E-07,4.12E-12,0.853497829,2.27E-06,3.94E-09,0.02763797,3.23E-09,0.013725779,1.83E-08,4.38E-09,1.99E-08,0.004552181,0.008684847,0.090783322,7.28E-05,1.12E-10,7.24E-07,3.10E-07,1.00E-09,2.46E-11,1.18E-05,2.33E-09,0.000825926,2.43E-07,0.000203771,1.26E-08,1.15E-07,2.34E-07,5.78E-11,2.17E-08,5.76E-07,7.70E-12,0.853497829,0.37068347,0.166923237,0.461305996,3.33E-07,2.08E-07,0.001083968,5.21E-09,3.00E-09,2.76E-06,1.26E-12,2.24E-08,2.43E-06,2.60E-09,0.027529114,1.41E-09,0.013803132,9.28E-09,6.35E-09,1.41E-07,0.00458323,0.008715574,0.090737578,5.02E-05,5.64E-11,4.94E-07,3.19E-07,2.15E-09,9.91E-12,1.27E-05,8.11E-10,0.000857339,1.57E-07,0.00020987,2.55E-09,1.29E-07,2.43E-07,2.50E-11,5.00E-09,6.09E-07,5.97E-12,0.853497829,2.65E-06,4.79E-09,0.026411762,3.84E-09,0.014556787,2.10E-08,4.63E-09,2.09E-08,0.004847761,0.009128889,0.090330571,7.87E-05,1.44E-10,6.65E-07,3.39E-07,9.81E-10,3.21E-11,1.18E-05,2.52E-09,0.000899428,2.80E-07,0.000235779,1.68E-08,1.63E-07,2.56E-07,6.89E-11,2.52E-08,7.01E-07,1.14E-11,0.853497829,0.390899703,0.163497175,0.444548487,3.60E-07,2.16E-07,0.001051204,5.58E-09,3.13E-09,2.82E-06,1.49E-12,2.45E-08,2.85E-06,3.03E-09,0.026296226,1.62E-09,0.014639621,1.04E-08,6.81E-09,1.46E-07,0.00488144,0.009160019,0.090280817,5.39E-05,7.19E-11,4.49E-07,3.49E-07,2.03E-09,1.26E-11,1.28E-05,9.11E-10,0.000934064,1.74E-07,0.000243045,3.46E-09,1.83E-07,2.65E-07,2.85E-11,5.74E-09,7.41E-07,8.76E-12,0.853497829,3.12E-06,5.86E-09,0.025122925,4.57E-09,0.015435824,2.41E-08,4.89E-09,2.20E-08,0.005165432,0.009583342,0.089846642,8.52E-05,1.87E-10,6.05E-07,3.72E-07,9.57E-10,4.19E-11,1.19E-05,2.74E-09,0.000980258,3.23E-07,0.000273262,2.25E-08,2.33E-07,2.79E-07,8.21E-11,2.90E-08,8.54E-07,1.69E-11,0.853497829,0.412389481,0.159607135,0.426982742,3.89E-07,2.24E-07,0.001017113,5.98E-09,3.27E-09,2.88E-06,1.77E-12,2.67E-08,3.37E-06,3.55E-09,0.025001621,1.85E-09,0.015523626,1.18E-08,7.30E-09,1.51E-07,0.005201516,0.009614418,0.089792972,5.78E-05,9.23E-11,4.05E-07,3.84E-07,1.92E-09,1.63E-11,1.29E-05,1.03E-09,0.001018342,1.93E-07,0.000281887,4.73E-09,2.63E-07,2.90E-07,3.27E-11,6.55E-09,9.03E-07,1.30E-11,0.853497829,3.72E-06,7.20E-09,0.023773822,5.44E-09,0.016362168,2.78E-08,5.16E-09,2.30E-08,0.005505282,0.010045744,0.089330875,9.22E-05,2.44E-10,5.47E-07,4.11E-07,9.32E-10,5.54E-11,1.20E-05,2.99E-09,0.001068876,3.73E-07,0.000317091,3.03E-08,3.36E-07,3.06E-07,9.78E-11,3.32E-08,1.04E-06,2.55E-11,0.853497829,0.435080862,0.155228855,0.408704858,4.19E-07,2.33E-07,0.000981803,6.40E-09,3.41E-09,2.93E-06,2.09E-12,2.91E-08,4.02E-06,4.18E-09,0.023648075,2.12E-09,0.016454113,1.33E-08,7.82E-09,1.56E-07,0.005543372,0.010076252,0.089273553,6.20E-05,1.19E-10,3.63E-07,4.24E-07,1.81E-09,2.12E-11,1.31E-05,1.17E-09,0.001110584,2.14E-07,0.00032727,6.47E-09,3.81E-07,3.17E-07,3.75E-11,7.43E-09,1.10E-06,1.95E-11,0.853497829,4.46E-06,8.88E-09,0.02236893,6.44E-09,0.01733362,3.20E-08,5.44E-09,2.41E-08,0.005866655,0.010513092,0.088783337,9.97E-05,3.21E-10,4.90E-07,4.55E-07,9.06E-10,7.36E-11,1.21E-05,3.28E-09,0.001165604,4.31E-07,0.000368189,4.07E-08,4.87E-07,3.35E-07,1.16E-10,3.77E-08,1.27E-06,3.88E-11,0.853497829,0.458853357,0.150351072,0.389846428,4.51E-07,2.42E-07,0.000945437,6.85E-09,3.55E-09,2.97E-06,2.48E-12,3.18E-08,4.84E-06,4.94E-09,0.022240501,2.43E-09,0.01742854,1.51E-08,8.36E-09,1.60E-07,0.005906157,0.010542494,0.088722838,6.65E-05,1.54E-10,3.22E-07,4.70E-07,1.71E-09,2.77E-11,1.33E-05,1.33E-09,0.001211035,2.38E-07,0.00038011,8.86E-09,5.56E-07,3.47E-07,4.30E-11,8.37E-09,1.34E-06,2.95E-11,0.853497829,5.40E-06,1.10E-08,0.020915138,7.58E-09,0.018346233,3.69E-08,5.73E-09,2.51E-08,0.006247976,0.010981912,0.088205031,0.00010761,4.21E-10,4.36E-07,5.06E-07,8.80E-10,9.81E-11,1.23E-05,3.60E-09,0.001270559,4.97E-07,0.000427489,5.46E-08,7.11E-07,3.66E-07,1.37E-10,4.23E-08,1.55E-06,5.96E-11,0.853497829,0.483533701,0.144979079,0.3705752,4.85E-07,2.50E-07,0.000908236,7.31E-09,3.69E-09,3.00E-06,2.94E-12,3.46E-08,5.87E-06,5.84E-09,0.020786232,2.76E-09,0.018442589,1.71E-08,8.92E-09,1.63E-07,0.006288075,0.011009693,0.088142072,7.12E-05,2.01E-10,2.83E-07,5.24E-07,1.61E-09,3.64E-11,1.35E-05,1.51E-09,0.0013197,2.64E-07,0.000441314,1.21E-08,8.12E-07,3.79E-07,4.92E-11,9.34E-09,1.63E-06,4.50E-11,0.853497829,6.58E-06,1.35E-08,0.019421977,8.84E-09,0.019394066,4.24E-08,6.01E-09,2.60E-08,0.006646579,0.011448393,0.087598103,0.000115927,5.54E-10,3.83E-07,5.65E-07,8.52E-10,1.31E-10,1.25E-05,3.97E-09,0.001383573,5.72E-07,0.00049587,7.29E-08,1.04E-06,4.00E-07,1.61E-10,4.68E-08,1.88E-06,9.19E-11,0.853497829,0.508889348,0.139138252,0.351098067,5.19E-07,2.58E-07,0.000870482,7.79E-09,3.82E-09,3.02E-06,3.49E-12,3.77E-08,7.15E-06,6.91E-09,0.0192953,3.12E-09,0.019489866,1.93E-08,9.46E-09,1.66E-07,0.006686178,0.011474123,0.087533741,7.61E-05,2.61E-10,2.47E-07,5.84E-07,1.51E-09,4.79E-11,1.37E-05,1.73E-09,0.001436252,2.91E-07,0.000511698,1.64E-08,1.19E-06,4.14E-07,5.61E-11,1.03E-08,1.98E-06,6.89E-11,0.853497829,8.04E-06,1.67E-08,0.01790188,1.02E-08,0.020468857,4.85E-08,6.28E-09,2.68E-08,0.007058465,0.011908604,0.086966171,0.000124445,7.25E-10,3.34E-07,6.31E-07,8.22E-10,1.75E-10,1.27E-05,4.38E-09,0.001504085,6.55E-07,0.000574051,9.67E-08,1.51E-06,4.36E-07,1.87E-10,5.10E-08,2.27E-06,1.42E-10,0.853497829,0.534602959,0.132879317,0.331681272,5.55E-07,2.66E-07,0.000832542,8.27E-09,3.94E-09,3.04E-06,4.12E-12,4.10E-08,8.74E-06,8.14E-09,0.017780991,3.49E-09,0.020561315,2.17E-08,9.96E-09,1.67E-07,0.007095959,0.01193206,0.0869021,8.11E-05,3.37E-10,2.14E-07,6.52E-07,1.41E-09,6.28E-11,1.39E-05,1.96E-09,0.001559846,3.21E-07,0.000591838,2.21E-08,1.72E-06,4.51E-07,6.34E-11,1.13E-08,2.38E-06,1.05E-10,0.853497829,9.85E-06,2.04E-08,0.016370978,1.15E-08,0.021559197,5.51E-08,6.52E-09,2.74E-08,0.007477729,0.01235889,0.086315073,0.000132877,9.43E-10,2.89E-07,7.06E-07,7.91E-10,2.31E-10,1.29E-05,4.81E-09,0.001630904,7.46E-07,0.000662373,1.27E-07,2.19E-06,4.75E-07,2.13E-10,5.47E-08,2.72E-06,2.19E-10,0.853497829,0.560170794,0.126293752,0.312736504,5.90E-07,2.73E-07,0.000794986,8.74E-09,4.06E-09,3.04E-06,4.85E-12,4.45E-08,1.07E-05,9.52E-09,0.016261631,3.84E-09,0.021643484,2.42E-08,1.04E-08,1.66E-07,0.007510165,0.01238041,0.086254679,8.58E-05,4.31E-10,1.84E-07,7.27E-07,1.32E-09,8.16E-11,1.41E-05,2.21E-09,0.001688634,3.51E-07,0.000681701,2.92E-08,2.48E-06,4.90E-07,7.09E-11,1.21E-08,2.84E-06,1.60E-10,0.853497829,1.20E-05,2.46E-08,0.014851953,1.28E-08,0.022647772,6.16E-08,6.72E-09,2.77E-08,0.007894658,0.012796827,0.085655256,0.000140569,1.20E-09,2.48E-07,7.87E-07,7.56E-10,3.01E-10,1.31E-05,5.25E-09,0.001761506,8.40E-07,0.00076026,1.63E-07,3.13E-06,5.15E-07,2.38E-10,5.75E-08,3.22E-06,3.32E-10,0.853497829,0.584534134,0.119568596,0.295134263,6.26E-07,2.78E-07,0.00075899,9.15E-09,4.17E-09,3.05E-06,5.67E-12,4.83E-08,1.29E-05,1.09E-08,0.01476665,4.10E-09,0.02271278,2.64E-08,1.07E-08,1.64E-07,0.007914878,0.012818501,0.08560718,8.96E-05,5.32E-10,1.58E-07,8.08E-07,1.22E-09,1.02E-10,1.43E-05,2.44E-09,0.001818207,3.80E-07,0.000779547,3.73E-08,3.48E-06,5.29E-07,7.76E-11,1.28E-08,3.34E-06,2.36E-10,0.853497829,1.45E-05,2.87E-08,0.013383852,1.37E-08,0.023702083,6.68E-08,6.88E-09,2.76E-08,0.008289374,0.013224149,0.085009666,0.000145707,1.45E-09,2.12E-07,8.71E-07,7.13E-10,3.74E-10,1.33E-05,5.61E-09,0.001889758,9.28E-07,0.00086451,2.02E-07,4.34E-06,5.55E-07,2.58E-10,5.94E-08,3.75E-06,4.82E-10,0.853497829,0.604789986,0.113178364,0.281299824,6.61E-07,2.80E-07,0.000727737,9.41E-09,4.28E-09,3.08E-06,6.53E-12,5.26E-08

sol = euler_(consv_eqs,0.1,y0,dt)

