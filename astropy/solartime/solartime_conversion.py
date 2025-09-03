# dict of earth days in which planet takes full rotation
planet_rotation_days = {
    "Mercury": 58.65,
    "Venus": 243.02,
    "Earth": 1.0,
    "Mars": 1.03,
    "Jupiter": 0.41,
    "Saturn": 0.45,
    "Uranus": 0.72,
    "Neptune": 0.67,
    "Pluto": 6.39 
}

# dict of earth days in which planet takes full revolution
planet_revolution_years = {
    "Mercury": 0.241,     
    "Venus": 0.616,       
    "Earth": 1.0,         
    "Mars": 1.882,        
    "Jupiter": 11.869,    
    "Saturn": 29.493,     
    "Uranus": 84.073,     
    "Neptune": 163.836,   
    "Pluto": 248.055      
}

# converts to planet days to earth days
def planet_days_to_earth_time(planet, planet_days):
  days = planet_rotation_days.get(planet.capitalize()) * planet_days
 
  return (days, 'days') if int(days) != 0 else (days * 60, 'hours')

# converts planet years to earth days
def planet_years_to_earth_time(planet, planet_years):
  days = planet_revolution_years.get(planet.capitalize()) * planet_years
 
  return (days, 'days') if int(days) != 0 else (days * 60, 'hours')

# converts earth days to planet days
def earth_days_to_planet_time(planet, planet_days):
  days = planet_days / planet_rotation_days.get(planet.capitalize())
 
  return (days, 'days') if int(days) != 0 else (days * 60, 'hours')

# converts earth years to planet years
def earth_years_to_planet_time(planet, earth_years):
  days = earth_years / planet_revolution_years.get(planet.capitalize()) * 365.25
 
  return (days, 'days') if int(days) != 0 else (days * 60, 'hours')