import bs4
import pandas as pd
import requests


def parse_data(soup, map, data_tag, in_tag, val_tag):
    """
    parses xml data into a pandas dataframe according to a map dict
    :param soup: the xml data as a BeautifulSoup object
    :param map: a dictionary mapping desired pandas df column names to tags in xml
    :param data_tag: the tag under which a single data point is contained
    :param in_tag: bool - whether the info is in the tag itself (as an attribute) or outside the tag
    :param val_tag: the tag that the data value is under (
            so that it can be made a numeric not a str)
    :return: pandas dataframe with parsed data
    """
    datapoints = soup.find_all(data_tag)
    df = pd.DataFrame(columns=map.keys())
    for d in datapoints:
        data_dict = dict()
        for key, val in map.iteritems():
            if in_tag:
                v = d[val]
            else:
                try:
                    v = getattr(d, val).string
                except TypeError:
                    print "unknown file structure"
            v = float(v) if val == val_tag else v
            data_dict[key] = v
        df = df.append(data_dict,
                       ignore_index=True)
    return df


year = '2014'
url = "http://tidesandcurrents.noaa.gov/api/datagetter?begin_date={0}0101&" \
      "end_date={0}1231&station=8638610&product=high_low&datum=NAVD&units=metric&" \
      "time_zone=lst&application=web_services&format=xml".format(year)

response = requests.get(url)
soup = bs4.BeautifulSoup(response.text, 'lxml')
data_tag = "hl"
map = {'date': 't', 'type': 'ty', 'value': 'v'}
d = get_data(soup, map, data_tag, True, 'v')
