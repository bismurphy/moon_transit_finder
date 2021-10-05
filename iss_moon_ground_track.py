from skyfield.api import load,EarthSatellite,wgs84, Topos
from skyfield.constants import AU_KM
from skyfield.positionlib import Geocentric
from numpy.linalg import norm
import load_tle
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import math

def normalize(vector):
    return vector / norm(vector)
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
def draw_plot(matplotlib_axis, satellite, moon, timerange, annotate_line = True):
    length_limit = find_length_limit(430,10)
    RE = 6371
    lats = []
    longs = []
    datapoints = []
    all_lines = []
    last_length = 1
    for timeval in timerange:
        vec_from_moon_to_sat = (satellite - moon).at(timeval).position.km
        result = find_earth_intersect_point(satellite.at(timeval).position.km,vec_from_moon_to_sat,RE,length_limit,timeval,length_guess = last_length)
        if result is not None:
            lat_result,lon_result, last_length = result
            datapoints.append([lat_result,lon_result,timeval])
        #Should we plot?
        if result is None or timeval == timerange[-1] or (len(longs) > 0 and lon_result - longs[-1] < -300):
            if len(datapoints) > 0:
                if annotate_line:
                    for dp in datapoints:
                        timetext = dp[2].utc_strftime("%H:%M:%S")
                        matplotlib_axis.annotate(timetext,
                                (dp[1],dp[0]),
                                textcoords="offset points",
                                xytext=(5,-5),
                                color="w",
                                transform=ccrs.PlateCarree())
                plotlats,plotlons = list(zip(*datapoints))[:2]
                plotline = matplotlib_axis.plot(plotlons,plotlats,color="w",transform=ccrs.PlateCarree())
                datapoints = []
        if result is not None:
            lats.append(lat_result)
            longs.append(lon_result)
    #for sake of knowing size to make the plot
    return lats,longs

if __name__ == "__main__":   
    ts = load.timescale()
    time_range = ts.utc(2021,10,11,21,36,range(0,60,10))

    planets = load('de421.bsp')
    earth = planets['earth']
    moon = planets['moon']

    #sat_tle = load_tle.get_tle(25544)

    sat_tle = load_tle.get_tle(25544)
    sat = EarthSatellite(*sat_tle)
    moon = moon - earth

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_facecolor("black")
    ax.add_feature(cartopy.feature.COASTLINE,edgecolor='lightgreen')
    states_provinces = cartopy.feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='110m',
        edgecolor='lightgreen',facecolor='none')
    ax.add_feature(cartopy.feature.BORDERS,edgecolor='lightgreen')
    ax.add_feature(states_provinces)

    draw_plot(ax, sat, moon, time_range)
    plt.show()
