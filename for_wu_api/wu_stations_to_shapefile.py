import urllib2
import json
import shapefile
import time

used_ids = []
station_list = []

date_range = [20130702,
              20131009,
              20140111,
              20140213,
              20140415,
              20140425,
              20140710,
              20140818,
              20140908,
              20140909,
              20140913,
              20141126,
              20141224,
              20150414,
              20150602,
              20150624,
              20150807,
              20150820,
              20150930,
              20151002]
api_key = "8842d598ec26dc90"
state = "VA"
city = "Virginia_Beach"
# get station data from wunderground for each date
for date in date_range:
    print "getting data for {}".format(str(date))
    f = urllib2.urlopen(
        'http://api.wunderground.com/api/{}/geolookup/history_{}/q/{}/{}.json'.format(api_key, str(date), state, city))
    json_string = f.read()
    parsed_json = json.loads(json_string)
    f.close()

    # write the data to a shapefile with 'name' and 'id' as the attributes

    for station in parsed_json['location']['nearby_weather_stations']['pws']['station']:
        if station['id'] == "KVAVIRGI131":
            raw_input('jjj')
        if station['id'] not in used_ids:
            station_list.append(station)
            used_ids.append(station['id'])
            print "added station id {}".format(station['id']),
            print "lon: {}".format(station['lon']), "lat: {}".format(station['lat'])
    time.sleep(6)

w = shapefile.Writer(shapefile.POINT)
w.field('id', 'C', 40)
w.field('name', 'C', 40)
w.field('lon', 'C', 40)
w.field('lat', 'C', 40)
for point in station_list:
    w.point(point['lon'], point['lat'])
    w.record(point['id'], point['neighborhood'], point['lon'], point['lat'])

path = 'shapefile/test/'
file_name = 'stations{}'.format("consolidated2")
ext = '.shp'
w.save('{}{}{}'.format(path, file_name, ext))


# create the PRJ file
ext = ".prj"
prj = open("{}{}{}".format(path, file_name, ext), "w")
epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
prj.write(epsg)
prj.close()
