import requests
import bs4

# get the days news site todo: get other tabs (2) (3)

date = "2014/09/08"
response = requests.get('http://wavy.com/{}/'.format(date))
soup = bs4.BeautifulSoup(response.text, "lxml")

# get headlines from the day and take only the ones with flood in them
hs = soup.select('h1.entry-title')
flood_hs = [h for h in hs if 'flood' in h.string.lower()]
# get the links to the flood stories
flood_as = [h.a.get('href') for h in flood_hs]
for a in flood_as:
    response = requests.get(a)
    soup = bs4.BeautifulSoup(response.text, 'lxml')