from skyfield.api import load,EarthSatellite,wgs84, Topos
from skyfield.constants import AU_KM
from skyfield.positionlib import Geocentric
from numpy.linalg import norm
import load_tle
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

def normalize(vector):
    return vector / norm(vector)
def find_earth_intersect_point(origin,vector,earth_radius):
    #Step 1: Turn that vector into a unit vector. Still pointing same direction.
    unit_vector_toward_ground = normalize(vec_from_moon_to_sat)
    #Step 2: Do a binary search to find required vector length.
    #Step 2a: First, extend the vector over and over until you punch through the earth.
    #Note: This step is risky if the vector is near-tangent with the earth.
    #Only use when you know the vector is coming in roughly vertical to the surface.
    length_guess = 1 / AU_KM
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
    while(abs(error) > 0.1 / AU_KM):
        #Standard binary search from here, go higher or lower based on where the guess winds up.
        new_guess = (length_lower_bound + length_upper_bound) / 2
        new_vector_end = norm(origin + unit_vector_toward_ground * new_guess)
        if new_vector_end > RE:
            length_lower_bound = new_guess
        else:
            length_upper_bound = new_guess
        error = new_vector_end - RE
    intersection_point = origin + unit_vector_toward_ground * new_guess

    geo = Geocentric(intersection_point,t=t)
    subpoint = wgs84.subpoint(geo)
    return subpoint.latitude.degrees,subpoint.longitude.degrees

ts = load.timescale()

t = ts.utc(2021,10,3,22,51,range(30,151))

DURANGO = 37.273267,-107.871692, 2000
LOS_ANGELES = 34.0,-118.2, 100
BOULDER = 40.015, -105.270556,1655
SEATTLE = 47.609722, -122.333056, 100
CAMBRIDGE = 42.371539,-71.098857, 20
NYC = 40.712778, -74.006111,20


RE = norm(wgs84.latlon(*SEATTLE).at(t).position.au)
planets = load('de421.bsp')
earth = planets['earth']
moon = planets['moon']

#sat_tle = load_tle.get_tle(25544)

sat_tle = [
"1 25544U 98067A   21273.48726090  .00001212  00000-0  30298-4 0  9993",
"2 25544  51.6456 182.0705 0003968  45.0544 125.1915 15.48872009304915"]
sat = EarthSatellite(*sat_tle)
moon = moon - earth

#Get vector from the moon to the satellite; extending this vector will make
#it run into the earth.

vec_from_moon_to_sat = (sat - moon).at(t).position.au
lat,long = find_earth_intersect_point(sat.at(t).position.au,vec_from_moon_to_sat,RE)

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
ax.set_extent(
    [min(long) - 2,
     max(long) + 2,
     min(lat) - 2,
     max(lat) + 2],
    crs=ccrs.PlateCarree())
#Plot initial location
ax.plot(long,lat,transform=ccrs.PlateCarree())
#Now for the good stuff: Labeling.
long_labels = long[::10]
lat_labels = lat[::10]
label_times = t[::10]
labels = zip(label_times,long_labels,lat_labels)
for l in labels:
    timetext = l[0].utc_strftime("%H:%M:%S")
    ax.annotate(timetext,
                (l[1],l[2]),
                textcoords="offset points",
                xytext=(5,-5),
                color="w",
                transform=ccrs.PlateCarree())
ax.scatter(long_labels,lat_labels,marker=".",transform=ccrs.PlateCarree())

plt.show()
