export cgrav, m_sun, m_earth, au, km, day_in_sec, year_in_sec

const cgrav = 6.67430e-8 # CODATA 2018

# GM product for Earth and Sun (cm^3 s^-2)
# IAU Resolution B3 (https://arxiv.org/abs/1510.07674)
const mu_sun = 1.3271244e26
const mu_earth = 3.986004e20
const m_sun = mu_sun/cgrav
const m_earth = mu_earth/cgrav

const au = 1.49597870700e13 # in cm, IAU 2009 system of astronomical constants (Luzum et al. 2011)
const km = 1e5 # in cm

const day_in_sec = 24*3600
const year_in_sec = 365.25*day_in_sec

