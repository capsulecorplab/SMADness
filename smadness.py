#!/usr/bin/env python3

import math
from collections import namedtuple
#from units import Units

# Gravitational Constant (universal) (m^2/(kg*s^2))
G = 6.674e-11

# Radius of earth's Orbit or 1 AU (m)
AU = 149.6e9

"""
class sun:
    def __init__(self, m, r_body):
        # Mass of Sun (kg)
        self.m = m
        # Radius of sun (m)
        self.r_body = r_body

class moon:
    def __init__(self, m, r_body, r_orbit):
        # Mass of moon (kg)
        self.m = m
        # Radius of moon (m)
        self.r_body = r_body
        # Semi-major axis of moon's orbit about planet (m)
        self.r_orbit = r_orbit

class planet:
    def __init__(self, m, r_body, r_orbit, b, phi_solar):
        # Mass of planet (kg)
        self.m = m
        # Radius of planet (m)
        self.r_body = r_body
        # Semi-major axis of planet's orbit about sun (m)
        self.r_orbit = r_orbit
        # Magnetic field (T)
        self.b = b
        # Solar flux constant at ~1AU (W/m^2)
        self.phi_solar = phi_solar

# Sun
sun = sun(m=1.989e30, r_body=695700e3)

# Earth
earth = planet(m=5.972e24, r_body=6371e3, r_orbit=1*AU, b=3e-5, phi_solar=1.3608e3)
# Moon
moon = moon(m=7.34e22, r_body=1.737e6, r_orbit=384748e3)

# Pluto
pluto = planet(m=1.48e15, r_body=6.2e3, r_orbit=39.52*AU, phi_solar=0.9)

# Jupiter
jupiter = planet(m=1.8987e27, r_body=69911e3, r_orbit=5.204*AU)
"""

# Gravitational Constant (local)
def g(G,M,R):
    return G*M/R**2

# Gravity well
def phi_gw(G,M,g,R):
    return G*M/(g*R)

# Escape Velocity (m/s)
def v_esc(G,M,R):
	return math.sqrt(2*G*M/R)

# Distance from Earth to Spacecraft
def l_e(d_s_e,M_s,M_e):
	return d_s_e/(1 + math.sqrt(M_s/M_e))

# Semimajor axis
def a(hp, ha):
	return R + (hp + ha)/2

# Orbital period, T = 2pi * sqrt(a^3/mu)
def T(a,G,M):
	return 2*math.pi*math.sqrt(a**3/(G*M))

# Semimajor axis (given T) (m)
def a_T(G,M,T):
	return (G*M*(T/(2*math.pi))**2)**(1./3)

# Orbital Velocity (vis-viva equation) for elliptical trajectory
def v(G,M,r,a):
	return math.sqrt(2*G*M/r - G*M/a)

# Orbital Velocity (vis-viva equation) for hyperbolic trajectory
def v_hyperbolic(G,M,r,a):
	return math.sqrt(2*G*M/r + G*M/a)

# Gravitational potential
def V(G,M,R):
    return G*M/R

# Homann Transfer, deltav1
def deltav1_ht(G,M,r1,r2):
	mu = G*M
	return math.sqrt( (2*mu*r2) / (r1*(r1 + r2)) ) - math.sqrt(mu/r1)

# Homann Transfer, deltav2
def deltav2_ht(G,M,r1,r2):
	mu = G*M
	return - math.sqrt( (2*mu*r1) / (r2*(r1 + r2)) ) + math.sqrt(mu/r2)

# Deltav required for change of orbital plane at equator crossing
def deltav_op(Vi,alpha):
	return 2*Vi*math.sin(alpha/2)

# Nodal Regression rate for LEO (in deg per mean solar day)
def omegap_degperday(i,a,e):
	i = i*180/math.pi
	a = a*1e-3
	return -2.06474e14*math.cos(i)/(a**3.5*(1-e**2)**2)

# Nodal Regression rate for LEO (rad/s)
def omegap(R,a,e,J2,omega,i):
#	return -3./2. * R**2 / (a * (1. - e**2) )**2 * J2 * omega * math.cos(i)
	return - 3./2 * R**2 / a**2 * J2 * omega * math.cos(i)

# Hyperbolic excess velocity (given departure velocity)
def v_inf(G,M,rd,vd):
	return math.sqrt( vd**2 - 2*G*M/rd )

# Departure velocity (given hyperbolic excess velocity)
def v_d(G,M,rd,v_inf):
	return math.sqrt( 2*G*M/rd + v_inf**2)

# Heliocentric departure velocity (slingshot maneuver)
def V_Sd(Vp,v_dinf,delta):
	return math.sqrt( Vp**2 + v_dinf**2 - 2*Vp*v_dinf*math.cos(delta) )

# deltav braking
def deltav_brake(G,M,vainf,rp,ai):
	return math.sqrt( 2*G*M/rp + vainf**2 ) - math.sqrt( 2*G*M/rp - G*M/ai )

# Tsiolkovsky (Rocket) equation
def deltav_rocket(g,Isp,mi,mf):
	return g*Isp*math.log(mi/mf)

# Mass of propellant, given initial spacecraft mass (kg)
def m_p(mi,deltav,g0,Isp):
	return mi*(1 - math.exp(-deltav/(g0*Isp)))

# Attitude maneuver about an axis - burn time for req'd theta
def t_b_theta(theta_m,n,F,L,I):
	return math.sqrt(theta_m*I/(n*F*L))

# Attitude maneuver about an axis - burn time for req'd omega
def t_b_omega(omega_max,n,F,L,I):
	return omega_max*I/(n*F*L)

# Attitude maneuver about an axis - Propellant used
def m_p(n,F,t_b,g,Isp):
	return 2*n*F*t_b/(g*Isp)

# Induced voltage from tether (Faraday's law of induction) assuming perpendicular
def U_tether(V,B,L):
	return V*B*L 

# Probability of Proper Functioning
def R(Lambda,t):
	return math.exp(-Lambda*t)


