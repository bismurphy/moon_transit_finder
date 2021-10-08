#Find a location to take a photo of the ISS transiting the moon!

#Generate a view of the stars (and any other desired objects) from a given
#observer point, facing a given direction.
#This will use a distortion-heavy map projection with horizon at y=0 and
#zenith at y=90, with x being azimuth

from skyfield.api import load, Star, wgs84, EarthSatellite, Topos
from skyfield.data import hipparcos
from skyfield import named_stars
from skyfield.framelib import ecliptic_frame
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button
import cartopy
import cartopy.crs as ccrs
import requests
import load_tle
import iss_moon_ground_track

STAR_MAG_LIMIT = 3
SAT_SIZE = 0.1 #size to draw sats on starmap
DURANGO = 37.273267,-107.871692, 2000
LOS_ANGELES = 34.0,-118.2, 100
BOULDER = 40.015, -105.270556,1655
SEATTLE = 47.609722, -122.333056, 100
CAMBRIDGE = 42.371539,-71.098857, 20
NYC = 40.712778, -74.006111,20

OBS1 =   41.600121, -70.901738  , 0
OBS2 = 41.591638, -70.901611, 0
OBS3 = 41.539012, -70.945050,0
INIT_LAT,INIT_LON, ELEVATION = CAMBRIDGE
sat_tle = load_tle.get_tle(25544)

LAT_RANGE = 2
LON_RANGE = 2
#seconds you can scroll before and after start time
TIME_SLIDER_RANGE = 60


TIME = [2021, 10, 12, 0, 0,0] #Remember to use UTC!

DURATION = 1 #days to search through

SAT_IMAGE = plt.imread('iss_white.png')
MOON_IMAGE = plt.imread('moon.png')
image_extent = 0.2

def setup_plot():
    axes[0].set_facecolor("black")
    axes[1].set_facecolor("black")


    plt.subplots_adjust(bottom=0.1,left=0.1,top=0.9,right=1,wspace=0)
    global time_slider
    slider_min = PLOT_TIME.tt - TIME_SLIDER_RANGE/86400
    slider_max = PLOT_TIME.tt + TIME_SLIDER_RANGE/86400
    time_slider = Slider(
        ax=plt.axes([0.2, 0.05, 0.6, 0.03],facecolor='k'),
        label="Time",
        valmin=slider_min,
        valmax=slider_max,
        valstep = np.linspace(slider_min,slider_max,TIME_SLIDER_RANGE*2+1),
        valinit = PLOT_TIME.tt,#Start at the current second in our minute
        )
    time_slider.valtext.set_text(PLOT_TIME.utc_strftime("%Y-%m-%d\n%H:%M:%S"))
    time_slider.on_changed(time_update)
    latitude_slider = Slider(
        ax = plt.axes([0.05,0.25,0.02,0.6],facecolor='k'),
        label="Latitude",
        valmin=INIT_LAT - LAT_RANGE,
        valmax=INIT_LAT + LAT_RANGE,
        valinit=INIT_LAT,
        orientation="vertical")
    latitude_slider.on_changed(latitude_update)
    longitude_slider = Slider(
        ax = plt.axes([0.2,0.95,0.6,0.03],facecolor='k'),
        label="Longitude",
        valmin=INIT_LON - LAT_RANGE,
        valmax=INIT_LON + LAT_RANGE,
        valinit=INIT_LON)
    longitude_slider.on_changed(longitude_update)
    extent = [
        INIT_LON-LON_RANGE,
        INIT_LON+LON_RANGE,
        INIT_LAT-LAT_RANGE,
        INIT_LAT+LAT_RANGE]
    axes[1].set_extent(extent,crs=ccrs.PlateCarree())
    states_provinces = cartopy.feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='110m',
        edgecolor='lightgreen',facecolor='none')
    
    #axes[1].add_feature(cartopy.feature.LAND)
    axes[1].add_feature(cartopy.feature.COASTLINE,edgecolor='lightgreen')
    axes[1].add_feature(states_provinces)
    #Plot initial location
    axes[1].plot([LONGITUDE],[LATITUDE],'yo',markersize=10,transform=ccrs.PlateCarree())
    
    return [time_slider,latitude_slider,longitude_slider]
def get_stars_with_names(mag_limit=100):
    with load.open(hipparcos.URL) as f:
        df = hipparcos.load_dataframe(f)
        df_filtered = df[df['magnitude'] <= mag_limit]
        hip_numbers = df_filtered.index.values.tolist()
        #Produce a dictionary that takes HIP numbers as keys and returns names
        name_dict = {v: k for k, v in named_stars.named_star_dict.items()}
        star_names = []
        for star in hip_numbers:
            if star in name_dict:
                star_names.append(name_dict[star])
            else:
                star_names.append("HIP " + str(star))
    return Star.from_dataframe(df_filtered), star_names, df_filtered['magnitude']
def plot_stamp(target_x, target_y, target_image, size,zorder=0.5):
    image_extent = [target_x - size, target_x + size,target_y - size, target_y + size]
    #Gets vector from sat to sun and converts to radec
    return axes[0].imshow(target_image,extent=image_extent,zorder=zorder)
def plot_sat(TLE, image, timestamp):
    global LATITUDE,LONGITUDE
    sat = EarthSatellite(*TLE)
    ground = Topos(LATITUDE, LONGITUDE, elevation_m = ELEVATION)
    alt,az,dist = (sat - ground).at(timestamp).altaz()

    #Also get where it was a bit ago so we can draw a little line
    a_bit_ago = ts.tt_jd(timestamp.tt -10/86400)

    prev_alt,prev_az,prev_dist = (sat - ground).at(a_bit_ago).altaz()
    trail = axes[0].plot([prev_az.degrees,az.degrees],[prev_alt.degrees,alt.degrees],color='w')
    stamp = plot_stamp(az.degrees, alt.degrees, image, SAT_SIZE,zorder=0.6)
    return [trail,stamp]
def plot_moon(image,timestamp):
    global LATITUDE,LONGITUDE
    ground = earth + Topos(LATITUDE, LONGITUDE, elevation_m = ELEVATION)
    alt,az,dist = (moon - ground).at(timestamp).altaz()

    moon_angular_radius = np.arcsin(1737.4/dist.km)*180/np.pi #moon radius
    moon_phase = get_moon_phase(ground, moon, sun,timestamp)
    moonmask = draw_moonmask(moon_phase,moon_angular_radius,az.degrees, alt.degrees)
    crosshair = plot_stamp(az.degrees,alt.degrees,plt.imread('crosshair.png'),5)
    print(az.degrees,alt.degrees)
    return plot_stamp(az.degrees, alt.degrees, image, moon_angular_radius), moonmask,crosshair
def get_moon_phase(observer,moon,sun,time):
    sun_longitude = (sun - observer).at(time).frame_latlon(ecliptic_frame)[1]
    moon_longitude = (moon - observer).at(time).frame_latlon(ecliptic_frame)[1]
    return (moon_longitude.degrees - sun_longitude.degrees) % 360
#Make a shadow over the proper part of the moon to show its phase
def draw_moonmask(phase,moon_size,center_x,center_y):
    #Start from the top and draw a circle going down
    y_vals_1 = np.linspace(1,-1,1000)
    x_vals_1 = center_x + (-moon_size if phase < 180 else moon_size) * np.sqrt(1-y_vals_1**2)
    y_vals_1 *= moon_size
    y_vals_1 += center_y
    #Now from bottom up, draw the curve.
    y_vals_2 = np.linspace(-1,1,1000)
    #How far displaced is the non-circular portion from center
    x_displacement = 1-(phase % 180 / 90)
    x_displacement *= moon_size
    x_vals_2 = center_x + x_displacement*np.sqrt(1-y_vals_2**2)
    y_vals_2 *= moon_size
    y_vals_2 += center_y
    points = list(zip(x_vals_1,y_vals_1)) + list(zip(x_vals_2,y_vals_2))
    return axes[0].add_patch(Polygon(points,color='k',alpha=0.7,linewidth=None,zorder=0.55))

def plot_starfield(timestamp):
    global LATITUDE,LONGITUDE
    bright_stars, star_names, mags = get_stars_with_names(STAR_MAG_LIMIT)

    observer = earth + wgs84.latlon(LATITUDE,LONGITUDE)

    alt, az,_ = observer.at(timestamp).observe(bright_stars).apparent().altaz()
    mystars = list(zip(az.degrees, alt.degrees, mags)) #swap to put x first
    visible_stars = [s for s in mystars if s[1] > 0]
    mags = np.array([s[2] for s in visible_stars])
    visible_stars = [s[:2] for s in visible_stars]
    dimmest = max(mags)
    mags = dimmest - mags
    #Mags are logarithmic. Undo the logarithm.
    mags = 2.512 ** mags
    #Now we have a linear scale. Now normalize it from 0 to 1.
    smallest = min(mags) - 1
    largest = max(mags)
    mags = 100*(mags - smallest) / (largest - smallest)
    return axes[0].scatter(*zip(*visible_stars), s=mags,color='w',marker="*")

def update_plot(timestamp):
    global plotted_objects

    for i in plotted_objects:
        if type(i) == type([]) and len(i) == 1:
            i = i[0]
        i.remove()
    starfield = plot_starfield(timestamp)
    iss_trail_and_stamp = plot_sat(sat_tle, SAT_IMAGE, timestamp)
    moon_stamp_and_mask = plot_moon(MOON_IMAGE,timestamp)
    map_marker = update_map()
    if len(plotted_objects) == 0:
        iss_extent = iss_trail_and_stamp[1].get_extent()
        moon_extent = moon_stamp_and_mask[0].get_extent()
        iss_xvals = iss_extent[0:2]
        iss_yvals = iss_extent[2:4]
        moon_xvals = moon_extent[0:2]
        moon_yvals = moon_extent[2:4]

        minx = min(*iss_xvals,*moon_xvals)-10
        maxx = max(*iss_xvals,*moon_xvals)+10
        miny = min(*iss_yvals,*moon_yvals)-10
        maxy = max(*iss_yvals,*moon_yvals)+10

        axes[0].set_xlim(minx,maxx)
        axes[0].set_ylim(miny,maxy)
    plotted_objects = [starfield,*iss_trail_and_stamp,*moon_stamp_and_mask,map_marker]
def update_map():
    global LATITUDE,LONGITUDE
    print(LATITUDE,LONGITUDE)
    return axes[1].plot([LONGITUDE],[LATITUDE],'rX',markersize=10,transform=ccrs.PlateCarree())

def time_update(slider_position):
    global PLOT_TIME
    print(slider_position*86400)
    PLOT_TIME = ts.tt_jd(slider_position)
    time_slider.valtext.set_text(PLOT_TIME.utc_strftime("%Y-%m-%d\n%H:%M:%S"))
    update_plot(PLOT_TIME)#And now replot everything.
def latitude_update(slider_position):
    global LATITUDE
    LATITUDE = slider_position
    update_plot(PLOT_TIME)#And now replot everything.
def longitude_update(slider_position):
    global LONGITUDE
    LONGITUDE = slider_position
    update_plot(PLOT_TIME)#And now replot everything.
def find_closest_approach(tle):
    global LATITUDE,LONGITUDE
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
axes = [None,None]
axes[0] = fig.add_subplot(1,2,1)
axes[1] = fig.add_subplot(1,2,2, projection=ccrs.PlateCarree())

ts = load.timescale()
global PLOT_TIME, LATITUDE, LONGITUDE
LATITUDE, LONGITUDE = INIT_LAT, INIT_LON
planets = load('de421.bsp')
earth = planets['earth']
moon = planets['moon']
sun = planets['sun']

raw_closest = find_closest_approach(sat_tle)
print(raw_closest)
PLOT_TIME = round_seconds(raw_closest)
sliders = setup_plot()#assign to variable to keep them alive
sat = EarthSatellite(*sat_tle)
timerange = np.linspace(PLOT_TIME.tt - 60/86400, PLOT_TIME.tt + 60/86400,25)
timerange = [ts.tt_jd(value) for value in timerange]
draw_line_on_map = iss_moon_ground_track.draw_plot(axes[1],sat,moon - earth, timerange)


global plotted_objects
plotted_objects = []
update_plot(PLOT_TIME)

plt.show()
