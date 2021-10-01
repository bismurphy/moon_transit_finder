from skyfield.api import load,EarthSatellite,wgs84, Topos
from skyfield.constants import AU_KM
from skyfield.positionlib import Geocentric
from numpy.linalg import norm
import load_tle

def normalize(vector):
    return vector / norm(vector)
def find_earth_intersect_point(origin,vector,earth_radius):
    #Step 1: Turn that vector into a unit vector. Still pointing same direction.
    unit_vector_toward_ground = normalize(vec_from_moon_to_sat)
    #Step 2: Do a binary search to find required vector length.
    #Step 2a: First, extend the vector over and over until you punch through the earth.
    #Note: This step is risky if the vector is near-tangent with the earth.
    #Only use when you know the vector is coming in roughly vertical to the surface.
    length_guess = 1
    while(norm(origin + unit_vector_toward_ground * length_guess) > RE):
        length_guess *= 2
    #Once we find that number, it's under the earth. That's an upper bound on length.
    length_upper_bound = length_guess
    #Meanwhile, the lower bound is the previous guess, which is just half as long.
    length_lower_bound = length_guess/2
    #Step 2b: Now that we have upper and lower bounds, binary search to find the
    #exact required length.
    error = 9999
    #Tolerance of 0.1 km should get us quite close.
    while(abs(error) > 0.1):
        #Standard binary search from here, go higher or lower based on where the guess winds up.
        new_guess = (length_lower_bound + length_upper_bound) / 2
        new_vector_end = norm(origin + unit_vector_toward_ground * new_guess)
        if new_vector_end > RE:
            length_lower_bound = new_guess
        else:
            length_upper_bound = new_guess
        error = new_vector_end - RE
    print(new_guess)
    intersection_point = origin + unit_vector_toward_ground * new_guess
    print(f"norm = {norm(intersection_point)}")
    print(type(intersection_point[0]))
    print(norm(intersection_point))
    geo = Geocentric(intersection_point / AU_KM,t=t)
    print(geo.distance().km)
    subpoint = wgs84.subpoint(geo)
    print(subpoint.elevation.km)

    return subpoint.latitude.degrees,subpoint.longitude.degrees

ts = load.timescale()

t = ts.utc(2021,9,23,4,31,59)

CAMBRIDGE = 42.371539,-71.098857, 0
ODDSITE = 40.33825356986421,-72.86106736111111,0
DURANGO = 37.273267,-107.871692, 2000


RE = norm(wgs84.latlon(*DURANGO).at(t).position.km)
planets = load('de421.bsp')
earth = planets['earth']
moon = planets['moon']

sat_tle = load_tle.get_tle(25544)
sat = EarthSatellite(*sat_tle)
moon = moon - earth

#Get vector from the moon to the satellite; extending this vector will make
#it run into the earth.
vec_from_moon_to_sat = (sat - moon).at(t).position.km
lat,long = find_earth_intersect_point(sat.at(t).position.km,vec_from_moon_to_sat,RE)
print(lat,long)

