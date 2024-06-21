export cgrav, m_sun, m_earth, au, r_sun, km, day, year, km_per_s, degree

const cgrav = 6.67430e-8 # CODATA 2018

# GM product for Earth and Sun (cm^3 s^-2)
# IAU Resolution B3 (https://arxiv.org/abs/1510.07674)
const mu_sun = 1.3271244e26
const mu_earth = 3.986004e20
const m_sun = mu_sun/cgrav     # in g  
const m_earth = mu_earth/cgrav # in g 

const au = 1.49597870700e13 # in cm, IAU 2009 system of astronomical constants (Luzum et al. 2011)
const r_sun = 6.9566e10 # in cm, Schmutz & Kosovichev (2008)
const km = 1.0e5 # in cm

const day = 24.0*3600.0  # in s 
const year = 365.25*day  # in s  

const km_per_s = km # just for convenience

const degree = pi/180 #internally we use radians