#!/usr/bin/env python

import math
from collections import namedtuple
#from units import Units

G = 6.674e-11

'''
# Demois
M = 1.48e15
R = 6.2e3

# Pluto
M = 1.3e22
R = 1200e3
'''

# Earth
M_e = 5.972e24	# Mass of Earth (kg)
R_e = 6371e3	# Radius of Earth (m)
AU = 149.6e9	# Radius of earth's Orbit or 1AU (m)
B_e = 3e-5	# Earth's magnetic field (T)

# Sun
M_s = 1.989e30	# Mass of Sun (kg)
R_s = 695700e3	# Radius of Sun (m)

# Moon
M_m = 7.34e22	# Mass of moon (kg)
R_m = 1.737e6	# Radius of moon (m)

# Jupiter
M_J = 1.8987e27	# Mass of Jupiter
R_J = 5.204*AU	# Radius of Jupiter's orbit
r_J = 69911e3	# Mean radius of Jupiter

# Gravitational Constant
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

#hp = 250e3
#ha = 800e3

# Semimajor axis
def a(hp, ha):
	return R + (hp + ha)/2

# Orbital period, T = 2pi * sqrt(a^3/mu)
def T(a,G,M):
	return 2*math.pi*math.sqrt(a**3/(G*M))

# Semimajor axis (given T) (m)
#T = 86164
def a_T(G,M,T):
	return (G*M*(T/(2*math.pi))**2)**(1./3)

#v = 7.76e3

# periapsis | apoapsis
#a = 228e6
#e = 0.09336
#ra = (e + 1)*a
#rp = (1 - e)*a

#e = .3
#ra = 20000e3 + R
#a = ra/(1 + e)
#rp = (1 - e)*a
#r = rp

# Orbital Velocity (vis-viva equation) for elliptical trajectory
def v(G,M,r,a):
	return math.sqrt(2*G*M/r - G*M/a)

# Orbital Velocity (vis-viva equation) for hyperbolic trajectory
def v_hyperbolic(G,M,r,a):
	return math.sqrt(2*G*M/r + G*M/a)

# Gravitational potential
#V = G*M/R

#R = 700
#g = 140.15/R**2
#g = v**2/R

#a = R + 400e3

#time since periapsis
#E = 68.08*math.pi/180
#n = math.sqrt(G*M/a**3)
#t = (E-e*math.sin(E))/n
#v = math.sqrt(2*G*M/r - G*M/a)

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

## Nodal Regression rate for LEO (in deg per mean solar day)
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

#print (1 - U_tether(v(G,M_e,296e3+R_e,296e3+R_e), B_e, 20e3)*math.cos(30.*math.pi/180)/U_tether(v(G,M_e,296e3+R_e,296e3+R_e), B_e, 20e3))*100
#print U_tether(v(G,M_e,296e3+R_e,296e3+R_e), B_e, 20e3)*math.cos(30.*math.pi/180)

# HW 5.4.3
#print R(20**-1,12)

# HW 5.3.6
#print v(G, M_e, 296e3+R_e, 296e3+R_e)
#print U_tether(v(G,M_e,296e3+R_e,296e3+R_e), B_e, 20e3)

# HW 5.2.8
#print 2*t_b(math.pi/2,2,3.5,3.6/2,500)
#print m_p(2,3.5,t_b(math.pi/2,2,3.5,3.6/2,500),9.81,190)*1000

#mp = m_p(4050,11.50e3-7.78e3,9.81,320)
#print mp/5

#print deltav_brake(G,M_J,5.64e3,100000e3,100000e3)

#print v_d(G,M_e,200e3,8.79e3)
#print v_d(G,M_e,200e3,8.79e3)
#print 0.5*v_d(G,M_e,200e3,38.58e3)**2
#print V_Sd()
#print v_d(G,M_e,0.924e9,17e3)

#print deltav_rocket(9.81,311,280+56+10.45,10.45+3+20+53) + deltav_rocket(9.81,342,13+53.45,13.45)
#print deltav_rocket(9.81,342,313+33.45,33.45)

#wrong:print deltav_rocket(9.81,311,280+56+10,10+56) + deltav_rocket(9.81,342,56+10,10)
#wrong:print deltav_rocket(9.81,342,23+313+10,10)

#print V_Sd(37.4,37.4,27*math.pi/180)

'''
vinf = v_inf(G,M_e,1000e3+R_e,11e3)
print vinf
print v_hyperbolic(G,M_e,0.924e9,G*M_e/vinf**2)
print v_d(G,M_e,1000e3+R_e,vinf)
'''

#r1 = 230e3 + R
#r2 = 20000e3 + R

#deltav1 = deltav1_ht(G,M,r1,r2)
#deltav2 = deltav2_ht(G,M,r1,r2)
#print deltav1
#print deltav2

'''
T = 24.*60.*60.
#print a_T(G,M,T) - R
i = math.pi/2
a = 350e3 + R
e = 0.
omega = 2*math.pi/T
#J2 = 2*e/3 - R**3*omega**2/(3*G*M)
J2 = -1.08262668e-3

print R
print i
print J2
print omega
print omegap(R,a,e,J2,omega,i)
print omegap(R,a,e,J2,omega,i)*T*180/math.pi
print omegap_degperday(i,a,e)
print 5.19671760308e-19*T*180/math.pi

T = 365*24*60*60
omega = 2.*math.pi/T
print omega
print r_e
print G*M_s
print G*M_e
print G*M_m
'''


