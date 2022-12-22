import os,sys
import glob
path = sys.argv[1]
#path = "/Users/ankitverma/Documents/ENCODE/"
extension = "*.df.txt"

for i in glob.glob(''.join(path + extension), recursive=False):
    temp_id = i.strip().split(path)
    temp_id = temp_id[1]
    if os.path.getsize(i) == 0:
        #print(type(i))
        print(i, temp_id, os.path.getsize(i))
    else:
        #print(type(i))
        temp = open(i, "r")
        temp = temp.read().strip().split("\n")
        print(i, temp_id, len(temp))

