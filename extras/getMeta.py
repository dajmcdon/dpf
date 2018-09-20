import requests
import csv
page = requests.get('http://mazurka.org.uk/auto/earis/mazurka68-3/')
from bs4 import BeautifulSoup
soup = BeautifulSoup(page.content, 'html.parser')
outfile = csv.writer(open('mazurka ids.csv', 'w'))
outfile.writerow(['performer', 'year', 'pid'])
for perf in soup.find_all('b'):
    text = perf.contents[0]
    name = text[:text.find(' ')].replace('-', '.')
    year = text[text.find('(')+1:text.find(')')]
    pid = perf.next_sibling.next_sibling.next_sibling[4:]
    if name == 'Flière':
        name = 'Fliere'
    if name == 'François':
        name = 'Francois'
    if name == 'Tomšič':
        name = 'Tomsic'
    outfile.writerow([name, year, pid])
