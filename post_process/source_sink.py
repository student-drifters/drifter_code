# source_sink.py
# Authors: JiM, hxu, Dylan Montagu
# Compiled by Dylan Montagu

"""
Given time range, geographic range and/or drifter ID, and whether geographic
range is source or sink, program plots drifter tracks on google maps.
Output is html file called "drifter_tracks.html"
Note from JiM: still need to run a sed command to add key after running this
""""

import datetime as dt
import pytz
import numpy as np
import pandas as pd
import gmplot
from dateutil.parser import parse


def hexcolors(n):
    """ Taken from JiM/hxu's drifter_functions.
        Compute a list of distinct colors, each of which is represented as an #RRGGBB value.
        Useful for less than 100 numbers."""
    if pow(n, float(1)/3) % 1 == 0.0:
        n += 1
    #make sure number we get is more than we need.
    rgbcolors=[]
    x = pow(n, float(1)/3)
    a = int(x)
    b = int(x)
    c = int(x)
    if a * b * c <= n:
        a += 1
    if a * b * c < n:
        b += 1
        c += 1

    for i in range(a):
        r = 254/a * i
        for j in range(b):
            s = 254/b * j
            for k in range(c):
                t = 254/c * k
                color = r, s, t
                rgbcolors.append(color)
    hexcolor=[]
    for i in rgbcolors:
        hexcolor.append("#{:02x}{:02x}{:02x}".format(int(i[0]), int(i[1]), int(i[2])))

    return hexcolor


def getobs_drift_byrange(gbox, input_time):
    """
    Taken from JiM/hxu's drifter_functions
    Grabs data from erddap url, return id, latitude,longitude, and times
    gbox includes 4 values, maxlon, minlon,maxlat,minlat, like:  [-69.0,-73.0,41.0,40.82]
    input_time can either contain two values: start_time & end_time OR one  value:interval_days
    and they should be timezone aware
    example: input_time=[dt(2012,1,1,0,0,0,0,pytz.UTC),dt(2012,2,1,0,0,0,0,pytz.UTC)] """

    lon_max=gbox[0]
    lon_min=gbox[1]
    lat_max=gbox[2]
    lat_min=gbox[3]

    mintime = input_time[0].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')  # change time format
    maxtime = input_time[1].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')

    # open url to get data
    url='http://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csvp?id,time,latitude,longitude&time>='\
    +str(mintime)+'&time<='+str(maxtime)+'&latitude>='\
    +str(lat_min)+'&latitude<='+str(lat_max)+'&longitude>='+str(lon_min)+'&longitude<='+str(lon_max)+'&orderBy("id,time")'

    df = pd.read_csv(url, skiprows=[1])
    for k in range(len(df)):
        df['time (UTC)'].values[k] = parse(df['time (UTC)'].values[k])

    return df['time (UTC)'].values, df.id.values, df['latitude (degrees_north)'].values, df['longitude (degrees_east)'].values


def getobs_drift_byid(id,input_time):
    """
    Taken from JiM/hxu's drifter_functions.
    Function written by Huanxin and used by getdrifter_erddap.py
    get data from url, return ids latitude,longitude, times
    input_time can either contain two values: start_time & end_time OR one value:interval_days
    and they should be timezone aware
    example: input_time=[dt(2012,1,1,0,0,0,0,pytz.UTC),dt(2012,2,1,0,0,0,0,pytz.UTC)] """

    mintime = input_time[0].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')  # change time format
    maxtime = input_time[1].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')

    # open url to get data
    url='http://comet.nefsc.noaa.gov/erddap/tabledap/drifters.csvp?id,time,latitude,longitude&time>='\
    +str(mintime)+'&time<='+str(maxtime)+'&id="'+str(id)+'"&orderBy("time")'

    df=pd.read_csv(url,skiprows=[1])
    for k in range(len(df)):
        df['time (UTC)'].values[k] = parse(df['time (UTC)'].values[k])
    df = df[df['longitude (degrees_east)'] <= -20]

    return df['time (UTC)'].values, df.id.values, df['latitude (degrees_north)'].values, df['longitude (degrees_east)'].values


def point_in_poly(x, y, poly):
    """ judge whether a site is in or out a polygon
        returns boolean value of whether point is in polygon or not """
    n = len(poly)
    inside = False
    p1x, p1y = poly[0]

    for i in range(n+1):
        p2x, p2y = poly[i % n]
        if ((y > min(p1y, p2y)) and (y <= max(p1y, p2y)) and (x <= max(p1x, p2x))):
            if p1y != p2y:
                xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
            if p1x == p2x or x <= xints:
                inside = not inside
        p1x, p1y = p2x, p2y

    return inside


def main(input_time, gbox, months, days, gbox_location):
    """ creates google map plot of drifter tracks that meet input criteria
        in:
            input_time = dt.datetime    # start_time, end_time
            gbox = float[4]             # maxlon, minlon,maxlat,minlat
            months = int, int           # months of interest
            days = dt.timedelta(days)   # maximum length of track to plot where, for example, days = 6*7 is six weeks
            gbox_location = strip       # "source" or "sink"
        example input:
            input_time = [dt.datetime(1980,1,1,0,0,0,0,pytz.UTC),dt.datetime(2018,10,15,0,0,0,0,pytz.UTC)]
            gbox = [-69.33,-69.75,41.5,41.]
            months = [8,9]
            days = dt.timedelta(days=6*7)
            gbox_location = "source"
    """
    input_time = input_time
    maxlon, minlon, maxlat, minlat = gbox[0], gbox[1], gbox[2], gbox[3]
    months = months
    days = days
    gbox_location = gbox_location

    polygon=[(maxlon,maxlat), (maxlon,minlat), (minlon,minlat), (minlon,maxlat)] #set polygon
    time, ids, lats, lons = getobs_drift_byrange(gbox,input_time)
    id = list(set(ids))
    colors = hexcolors(len(id))  #get hex colors,like '00FF00'

    gmap3 = gmplot.GoogleMapPlotter(np.mean([maxlat,minlat]), np.mean([maxlon,minlon]), 9)

    for k in range(len(id)):
        time, ids, lat, lon = getobs_drift_byid(id[k], input_time)  # get data by id
        for z in range(len(lat)-1):  # make plotting drifter start in gbox
         # ignore tracks beyond "days" long and months outside of "months"
            if (time[z] - time[0] < days) and (time[z].month >= months[0]) and (time[z].month <= months[-1]):
                inside = point_in_poly(lon[z],lat[z],polygon)
                if inside:
                    if gbox_location == "source":
                        lat = lat[z:]  # delete data which are before coming in the polygon
                        lon = lon[z:]
                        time = time[z:]
                    elif gbox_location == "sink":
                        lat = lat[:z]  # delete data which are before after in the polygon
                        lon = lon[:z]
                        time = time[:z]
                    gmap3.plot(lat, lon, colors[k], edge_width = 2.5)
                    break

    ranges=[(maxlat,maxlon), (minlat,maxlon), (minlat,minlon), (maxlat,minlon), (maxlat,maxlon)] #plot range you gave
    gmap3.plot(ranges[0], ranges[1])
    gmap3.draw('drifter_tracks.html')



if __name__ == "__main__":
    input_time = [dt.datetime(1980,1,1,0,0,0,0,pytz.UTC), dt.datetime(2018,10,15,0,0,0,0,pytz.UTC)] # start time and end time
    #gbox=[-70.035594,-70.597883,42.766619,42.093197] #  maxlon, minlon,maxlat,minlat for Stellwagen
    gbox = [-69.33, -69.75, 41.5, 41.0] #  maxlon, minlon,maxlat,minlat for Nant shoals for Cap Brady
    months = [8, 9] # months of interest
    days = dt.timedelta(days=6*7) # maximum length of track to plot where, for example, 6*7 is six weeks
    gbox_location = "source"

    main(input_time, gbox, months, days, gbox_location)















#gmap does not have a apikey function, so lets substitute it.
# sed '/<\/body>/i <script async defer src="https:\/\/maps.googleapis.com\/maps\/api\/js?v=3&client=gme-noaa&channel=NMFS.NEFSC.OCB.DRIFTERS&callback=initMap"><\/script>' test.html > sink_stellwagen.html
'''
tempHolder=''
oldLine='</body>'
newLine='</script><script async defer src="https://maps.googleapis.com/maps/api/js?v=3&client=gme-noaa&channel=NMFS.NEFSC.OCB.DRIFTERS&callback=initMap"></script></body>'
#Open the file created by gmplot, do a find and replace.
#My file is in flask, so it is in static/map.html, change path to your file
with open('/net/pubweb_html/drifter/test.html') as fh:
   for line in fh:
       tempHolder += line.replace(oldLine,newLine)
fh.close
#Now open the file again and overwrite with the edited text
fh=open('/net/pubweb_html/drifter/test.html', 'w')
fh.write(tempHolder)
fh.close
'''
