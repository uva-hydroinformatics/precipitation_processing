import csv
import urllib2
import json
import time


def write_list_to_csv(raw_data_list, data_type):
    with open("{}.csv".format(data_type), "a") as csvfile:
        writer = csv.writer(csvfile)
        csv_data_list = [
            data_type,
            station['id'],
            station['neighborhood'],
            str(date)[:4],
            str(date)[4:6],
            str(date)[6:8],
            raw_data_list['date']['hour'],
            raw_data_list['date']['min']
        ]
        if data_type == "observation":
            csv_data_list.append(raw_data_list['precip_totalm'])
        else:
            csv_data_list.append(raw_data_list['precipm'])
        writer.writerow(csv_data_list)


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
date_range = [20140908]
api_key = "8842d598ec26dc90"
state = "VA"
city = "Virginia_Beach"

# get station data from wunderground for given date
for date in date_range:
    print "getting data for date {}".format(str(date))
    f = urllib2.urlopen(
            "http://api.wunderground.com/api/{}/geolookup/history_{}/q/{}/{}.json".format(api_key,
                                                                                          str(date),
            state,
            city))
    json_string = f.read()
    parsed_json = json.loads(json_string)
    f.close()

    # get rainfall data for each station on that date
    stations = parsed_json['location']['nearby_weather_stations']['pws']['station']
    for station in stations:
        if station['id'] == "KVAVIRGI71":
            raw_input("we found it")
        print "getting data for {} station (id: {})".format(station['neighborhood'], station['id'])
        f = urllib2.urlopen(
                'http://api.wunderground.com/api/{}/history_{}/q/pws:{}.json'.format(api_key, str(date), station['id']))
        json_string = f.read()
        parsed_json = json.loads(json_string)
        f.close()
        daily_summary = parsed_json['history']['dailysummary']
        observations = parsed_json['history']['observations']
        write_list_to_csv(daily_summary[0], 'dailysummary')
        for observation in observations:
            write_list_to_csv(observation, 'observation')
        time.sleep(10)
