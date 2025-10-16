# -*- coding: cp936 -*-
import re
import linecache
import numpy as np
import os
a = []
i = 0
h = []

# data_name = 'push_side_bad_2'
#data_name = 'Exp_2/simulation/sweep/sweep_5'
data_name = 'Exp_2'
shoulder = []
elbow = []
wrist = []
add_value = 1
file_name = 'data_2/1-original_data/' + data_name + '/prediction.txt'
with open(file_name, 'r') as f:
    data = f.readlines()  # txt中所有字符串读入data

    for line in data:
        odom = line.split()  # 将单个数据分隔开存好
        if ((i - 16) % 20 == 0):
            if odom[0] != "-nan":
                add_value = 1
                shoulder.append(np.asarray([odom], dtype=float))
            else:
                add_value = 0
        if ((i - 19) % 20 == 0) & (add_value == 1):
            elbow.append(180 - np.asarray([odom[2]], dtype=float))
        if ((i - 18) % 20 == 0) & (add_value == 1):
            wrist.append(np.asarray([odom], dtype=float))
        i = i + 1


lenth = len(elbow)
shoulder = np.array(shoulder).reshape(lenth, -1)
zeros1 = np.zeros((lenth, 1))
zeros2 = np.zeros((lenth, 2))
elbow = np.array(elbow).reshape(lenth, -1)
wrist = np.array(wrist).reshape(lenth, -1)

all = np.hstack([zeros1, shoulder, zeros1, elbow, zeros2, wrist])

import csv

# open the file in the write mode
# write_down_path = input("Enter a write path: ")
f = open('./data_2/2-processed_data/' + data_name + '/all_angles.csv', 'w')
# f = open(write_down_path, 'w')
# create the csv writer
writer = csv.writer(f)

writer.writerows(all)

elbow = []
shoulder = []
wrist = []
hand = []
add_value = 1
i = 0
with open(file_name, 'r') as f:
    data = f.readlines()  # txt中所有字符串读入data

    for line in data:
        odom = line.split()  # 将单个数据分隔开存好

        if ((i - 3) % 20 == 0):
            if odom[0] != "-nan":
                add_value = 1
                shoulder.append(np.asarray([odom[:]], dtype=float))
            else:
                add_value = 0
        if ((i - 4) % 20 == 0) & (add_value == 1):
            elbow.append(np.asarray([odom[:]], dtype=float))
        if ((i - 5) % 20 == 0) & (add_value == 1):
            wrist.append(np.asarray([odom[:]], dtype=float))
        if ((i - 6) % 20 == 0) & (add_value == 1):
            hand.append(np.asarray([odom[:]], dtype=float))
        i = i + 1


lenth = len(elbow)
shoulder = np.array(shoulder).reshape(lenth, -1)
elbow = np.array(elbow).reshape(lenth, -1)
wrist = np.array(wrist).reshape(lenth, -1)
hand = np.array(hand).reshape(lenth, -1)
all = np.hstack([shoulder, elbow, wrist, hand])


m = open('./data_2/2-processed_data/' + data_name + '/all_cartesian_translation.csv', 'w')
writer = csv.writer(m)

writer.writerows(all)


#BaseWrist = []
#add_value = 1
#i = 0
#with open(file_name, 'r') as f:
#    BW = f.readlines()  # txt中所有字符串读入data

#    for line in BW:
#        odom = line.split()  # 将单个数据分隔开存好

#       if ((i - 20) % 22 == 0):
#           if odom[0] != "-nan":
#               add_value = 1
#                BaseWrist.append(np.asarray([odom[:]], dtype=float))
#            else:
#                add_value = 0
#        i = i + 1

#lenth = len(BaseWrist)
#BaseWrist = np.array(BaseWrist).reshape(lenth, -1)



#m = open('./data_2/2-processed_data/' + data_name + '/Base_to_Wrist.csv', 'w')
#writer = csv.writer(m)
#writer.writerows(BaseWrist)
