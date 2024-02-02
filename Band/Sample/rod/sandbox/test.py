import csv
import pprint

with open('te.dat')as f:
    reader = csv.reader(f)
    l = [row for row in reader]

print(l)