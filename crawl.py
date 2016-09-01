import requests
import bs4
import time

# get the days news site todo: get other tabs (2) (3)

date = "2014/09/08"
base_url = 'http://wavy.com'
response = requests.get('{}/{}/'.format(base_url, date))
soup = bs4.BeautifulSoup(response.text, "lxml")

# get headlines from the day and take only the ones with flood in them
hs = soup.select('h1.entry-title')
flood_hs = [h for h in hs if 'flood' in h.string.lower()]
# get the links to the flood stories
flood_as = [h.a.get('href') for h in [flood_hs[3]]]
photo_as = [a for a in flood_as if 'photomojo' in a]

response = requests.get(photo_as[0])
soup = bs4.BeautifulSoup(response.text, 'lxml')
p = soup.find(id='paginate')
num_photos = int(p.text.strip().split(" ")[-1])
descriptions = []

for i in range(num_photos):
    if i > 0:
        p = soup.find(id='paginate')
        a = p.find_all('a')[-1].get('href')
        response = requests.get('http://interactives.wavy.com{}'.format(a))
        soup = bs4.BeautifulSoup(response.text, "lxml")
        time.sleep(0.5)
    d = soup.find(id='photo_desc').text
    city = soup.find(id='photo_title')
    if d not in descriptions:
        descriptions.append((city.contents[0], d))
    print descriptions[i]
