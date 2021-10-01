from skyfield.api import load, wgs84, EarthSatellite, Topos

import matplotlib.pyplot as plt
import numpy as np

INIT_LAT,INIT_LON, ELEVATION = 40.33825356986421,-72.86106736111111,0

#Starting time to look for transits
START_TIME = [2021, 9, 20, 17, 0,0] #Remember to use UTC!

DURATION = 3 #days to search through

def plot_sat(TLE, timestamp):
    global LATITUDE,LONGITUDE
    sat = EarthSatellite(*TLE)
    ground = Topos(LATITUDE, LONGITUDE, elevation_m = ELEVATION)
    alt,az,dist = (sat - ground).at(timestamp).altaz()

    #Also get where it was a bit ago so we can draw a little line
    a_bit_ago = ts.tt_jd(timestamp.tt -10/86400)

    prev_alt,prev_az,prev_dist = (sat - ground).at(a_bit_ago).altaz()
    trail = axes.plot([prev_az.degrees,az.degrees],[prev_alt.degrees,alt.degrees],color='w')
    dot = plt.Circle((az.degrees,alt.degrees),0.1,color='r')
    axes.add_artist(dot)
def plot_moon(timestamp):
    global LATITUDE,LONGITUDE
    ground = earth + Topos(LATITUDE, LONGITUDE, elevation_m = ELEVATION)
    alt,az,dist = (moon - ground).at(timestamp).altaz()

    moon_angular_radius = np.arcsin(1737.4/dist.km)*180/np.pi #moon radius
    moon_disk = plt.Circle((az.degrees,alt.degrees),moon_angular_radius,color='w')
    axes.add_artist(moon_disk)
def find_closest_approach(tle):
    global LATITUDE,LONGITUDE
    sat = EarthSatellite(*tle)
    observer = Topos(LATITUDE,LONGITUDE, elevation_m = ELEVATION)
    start = ts.utc(*START_TIME)
    end = ts.tt_jd(start.tt + DURATION)
    times_and_events = sat.find_events(observer, start, end)
    passes = []
    last_start = 0
    for i in zip(*times_and_events):
        event_time = i[0]
        event_type = i[1]
        if event_type == 0:
            last_start = event_time
        if event_type == 2:
            passes.append([last_start,event_time])
    closest_moon_dist = 9999
    closest_moon_time = 0
    for p in passes:
        alts = []
        azes = []
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
    return closest_moon_time
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
#Take a skyfield time object with floating seconds and round to nearest int second
def round_seconds(skyfield_time):
    time_utc = list(skyfield_time.utc)
    time_utc[5] = round(time_utc[5])
    return ts.utc(*time_utc)
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
axes.set_aspect(1)

ts = load.timescale()
global PLOT_TIME, LATITUDE, LONGITUDE
LATITUDE, LONGITUDE = INIT_LAT, INIT_LON
planets = load('de421.bsp')
earth = planets['earth']
moon = planets['moon']
sun = planets['sun']

sat_tle = [
"1 25544U 98067A   21273.48726090  .00001212  00000-0  30298-4 0  9993",
"2 25544  51.6456 182.0705 0003968  45.0544 125.1915 15.48872009304915"]
raw_closest = find_closest_approach(sat_tle)
PLOT_TIME = round_seconds(raw_closest)
axes.set_facecolor("black")
axes.set_xlim(0,360)
axes.set_ylim(0,90)
plot_moon(PLOT_TIME)
plot_sat(sat_tle,PLOT_TIME)
print(PLOT_TIME.utc)
plt.show()
