import pandas as pd
import pdb

folder = "/project/UKB_moore/UKB_50978/"
people_who_quit1 = pd.read_csv(folder + "w50978_20200204.csv", header = None)
people_who_quit2 = pd.read_csv(folder + "w50978_20200820.csv", header = None)
people_who_quit3 = pd.read_csv(folder + "w50978_20210201.csv", header = None)
quitters = [people_who_quit1, people_who_quit2, people_who_quit3]
people_who_quit = pd.concat(quitters).drop_duplicates()
people_who_quit[1] = people_who_quit[0]
people_who_quit.to_csv("people_who_quit.txt", sep = "\t", header = False, index = False)
