#-----------------------------------------------------------------------
# By Hunter Brooks
#-----------------------------------------------------------------------
from catalog_module.importmodule import *
file = 'tot_goodstuff.csv'

print('------------------------------------------------')
print('Stared WRAP TEST')
print('------------------------------------------------')
print('')

raw_file = pd.read_csv(file)
ra = raw_file['ra'].tolist()
dec = raw_file['dec'].tolist()
radius = 150

from catalog_module.aw_search import allwise_image

# for i in range(len(ra)):
#     print(ra[i])
#     print(dec[i])
#     print(nsc_image(ra[i], dec[i], radius))
# print(nsc_image(164.410559533, -28.9760034474, radius))
# print(ra[0])
# print(dec[0])
print(allwise_image(1, 1, radius))

print('------------------------------------------------')
print('Finished WRAP TEST')
print('------------------------------------------------')
print('')

