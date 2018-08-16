import os
from glob import glob

list = ['samples\n']
for file in glob('*.cel'):
    filename = os.path.basename(file) + '\n'
    list.append(filename)

with open('celfile_list.csv', 'w') as file:
    file.writelines(list)
    file.close()