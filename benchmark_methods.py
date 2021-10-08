#It occurs to me that there are two methods of finding a nearby transit location.
#1. Find where the ISS is close to the moon in the sky
#2. Find where the moon-ISS vector lands and compare to where we are

import load_tle
from skyfield.api import load, Star, wgs84, EarthSatellite, Topos
import numpy as np
import time
import math
from numpy.linalg import norm
from skyfield.positionlib import Geocentric
from skyfield.constants import AU_KM

CAMBRIDGE = 42.371539,-71.098857, 20
LATITUDE, LONGITUDE, ELEVATION = CAMBRIDGE
ts=load.timescale()
TIME = [2021, 10, 11, 0, 0,0] #Remember to use UTC!
DURATION = 5

planets = load('de421.bsp')
earth = planets['earth']
moon = planets['moon']
sun = planets['sun']
#First step: Generate passes.

tle = load_tle.get_tle(25544)

sat = EarthSatellite(*tle)
observer = Topos(LATITUDE,LONGITUDE, elevation_m = ELEVATION)
start = ts.utc(*TIME)
end = ts.tt_jd(start.tt + DURATION)
times_and_events = sat.find_events(observer, start, end)
passes = []
last_start = 0
for i in zip(*times_and_events):
    event_time = i[0]
    event_type = i[1]
    if event_type == 0:
        last_start = event_time
    #double check that we've had a start
    if event_type == 2 and last_start != 0:
        passes.append([last_start,event_time])

def angular_separation(alt1,alt2,az1,az2):
    alt1 = alt1.radians
    alt2 = alt2.radians
    az1 = az1.radians
    az2 = az2.radians
    #http://spiff.rit.edu/classes/phys373/lectures/radec/radec.html
    cosprod = np.cos(np.pi/2 - alt1) * np.cos(np.pi/2 - alt2)
    sinprod = np.sin(np.pi/2 - alt1) * np.sin(np.pi/2 - alt2)
    sinprod = sinprod * np.cos(az1 - az2)
    return np.arccos(cosprod + sinprod)

#Now: Method 1: Find time ISS is close to moon, for the chosen observer.
start_1 = time.time()
closest_moon_dist = 9999
closest_moon_time = 0
for p in passes:
    drawtime = p[0]
    while(drawtime.tt < p[1].tt):
        drawtime = ts.tt_jd(drawtime.tt + 1/86400)
        alt,az,_ = (sat - observer).at(drawtime).altaz()
        moon_alt,moon_az,_ = (moon - (earth+observer)).at(drawtime).altaz()
        if moon_alt.degrees < 5:
            break #break out if the moon is crazy low for this pass
        #this is not proper math but good enough heuristic for now
        moon_dist = angular_separation(alt,moon_alt,az,moon_az)
        if moon_dist < closest_moon_dist and alt.degrees > 10 and moon_alt.degrees > 10:
            closest_moon_dist = moon_dist
            closest_moon_time = drawtime
end_1 = time.time()
duration_1 = end_1 - start_1
print(f"Elapsed time for method 1: {duration_1} s")
print(f"Result: {closest_moon_time.utc}")


#And method 2: Find time the ISS-moon ground point is close to observer.
def find_earth_intersect_point(origin,vector,earth_radius,length_limit,time,length_guess = 1):
    #Step 1: Turn that vector into a unit vector. Still pointing same direction.
    unit_vector_toward_ground = normalize(vector)
    #Step 2: Do a binary search to find required vector length.
    #Step 2a: First, extend the vector over and over until you punch through the earth.
    #Note: This step is risky if the vector is near-tangent with the earth.
    #Only use when you know the vector is coming in roughly vertical to the surface.
    while(norm(origin + unit_vector_toward_ground * length_guess) > earth_radius):
        length_guess *= 2
        if length_guess > length_limit*2:
            return None    
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
        if new_vector_end > earth_radius:
            length_lower_bound = new_guess
        else:
            length_upper_bound = new_guess
        error = new_vector_end - earth_radius
    intersection_point = origin + unit_vector_toward_ground * new_guess

    geo = Geocentric(intersection_point/AU_KM,t=time)
    subpoint = wgs84.subpoint(geo)
    return subpoint.latitude.degrees,subpoint.longitude.degrees, new_guess
def find_length_limit(sat_height,horizon_angle):
    RE = 6371
    #Lotta geometry to do.
    horizon_angle += 90
    horizon_angle *= (math.pi/180)
    law_of_sines_ratio = math.sin(horizon_angle) / (RE + sat_height)
    #Now use that to find A, the angle formed by two vectors:
    #Vector 1 is vector from sat to center of earth, vector 2 is
    #vector from sat to an observer.
    #sin(A) / RE = law of sines ratio
    A = math.asin(law_of_sines_ratio * RE)
    #Now subtract that from the total triangle being pi.
    #that gives B, the angle between observer and sat, seen from center of earth.
    B = math.pi - (A + horizon_angle)
    #And again, law of sines. sin(B) / limit = ratio.
    limit = math.sin(B) / law_of_sines_ratio
    return limit
def normalize(vector):
    return vector / norm(vector)
def angular_separation_degrees(alt1,alt2,az1,az2):
    alt1 = alt1*math.pi/180
    alt2 = alt2*math.pi/180
    az1 = az1*math.pi/180
    az2 = az2*math.pi/180
    #http://spiff.rit.edu/classes/phys373/lectures/radec/radec.html
    cosprod = np.cos(np.pi/2 - alt1) * np.cos(np.pi/2 - alt2)
    sinprod = np.sin(np.pi/2 - alt1) * np.sin(np.pi/2 - alt2)
    sinprod = sinprod * np.cos(az1 - az2)
    return np.arccos(cosprod + sinprod)

start_2 = time.time()
moon = moon - earth
RE = 6371
closest_moon_dist = 9999
closest_moon_time = 0
length_limit = find_length_limit(430,10)
last_length = 1
for p in passes:
    drawtime = p[0]
    while(drawtime.tt < p[1].tt):
        drawtime = ts.tt_jd(drawtime.tt + 1/86400)
        moon_alt,moon_az,_ = (moon-observer).at(drawtime).altaz()
        if moon_alt.degrees < 5:
            break #break out if the moon is crazy low for this pass
        vec_from_moon_to_sat = (sat - moon).at(drawtime).position.km
        result = find_earth_intersect_point(sat.at(drawtime).position.km,vec_from_moon_to_sat,RE,length_limit,drawtime,length_guess = last_length)
        if result is not None:
            int_lat, int_lon, last_length = result
            travel_distance_for_observer = angular_separation_degrees(int_lon, LONGITUDE, int_lat, LATITUDE)
            if travel_distance_for_observer < closest_moon_dist:
                closest_moon_dist = travel_distance_for_observer
                closest_moon_time = drawtime
end_2 = time.time()
duration_2 = end_2 - start_2
print(f"Elapsed time for method 2: {duration_2} s")
print(f"Result: {closest_moon_time.utc}")
