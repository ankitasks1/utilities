from pybedtools import BedTool
import os, sys
path = "./"
sample_list1 = open(''.join(path + sys.argv[1]), "r")
sample_list2 = open(''.join(path + sys.argv[2]), "r")
# sample_list1 = open('summit_sample_list.txt', "r")
# sample_list2 = open('summit_sample_list.txt', "r")
sample_list1 = list(sample_list1)
sample_list2 = list(sample_list2)
#print(sample_list)
if os.path.exists('One_to_one_encode_closest_summit_counts_out14.txt'):
    os.remove('One_to_one_encode_closest_summit_counts_out14.txt')

for list1 in sample_list1:
    for list2 in sample_list2:
        overlaps = {}
        counts = {}
        mybed1_peak = ''.join(path + list1).strip()
        mybed1_peak_id = mybed1_peak.split('/')[-1:][0].strip('_summits_sorted_merged.bed')
        #print(mybed1_peak_id)
        bedfile1 = BedTool(mybed1_peak)
        #counts[mybed1_peak_id] = len(bedfile1)
        #print(type(bedfile1))
        mybed2_peak = ''.join(path + list2).strip()
        mybed2_peak_id = mybed2_peak.split('/')[-1:][0].strip('_summits_sorted_merged.bed')
        #print(mybed2_peak_id)
        bedfile2 = BedTool(mybed2_peak)
        #counts[mybed2_peak_id] = len(bedfile2)
        # print(type(bedfile2))
        closest_file = bedfile1.closest(bedfile2, d=True)
        # print(closest_file)
        count = 0
        # predict no overlap between files
        if len(closest_file) <= 0:
            print(''.join(mybed1_peak_id + "\t" + mybed2_peak_id + "\t" + str(len(bedfile1)) + "\t" + str(len(bedfile2)) + "\t" + str(0) + "\n"))
            with open("One_to_one_encode_closest_summit_counts_out14.txt", "a") as myfinalout:
                #write the output
                myfinalout.write(''.join(mybed1_peak_id + "\t" + mybed2_peak_id + "\t" + str(len(bedfile1)) + "\t" + str(len(bedfile2)) + "\t" + str(0) + "\n"))
                myfinalout.close()
        else:
            for intcont in closest_file:
                # predict complete overlap (0) and less than equal to 150 bp distance of first file from second file
                if int(intcont[6]) >= int(0) and int(intcont[6]) <= int(150):
                    count += 1
                # predict greater than equal to 150 bp distance of first file from second file
                elif int(intcont[6]) > int(150):
                    count += 0

            print(''.join(mybed1_peak_id + "\t" + mybed2_peak_id + "\t" + str(len(bedfile1)) + "\t" + str(len(bedfile2)) + "\t" + str(count) + "\n"))
            with open("One_to_one_encode_closest_summit_counts_out14.txt", "a") as myfinalout:
                #write the output
                myfinalout.write(''.join(mybed1_peak_id + "\t" + mybed2_peak_id + "\t" + str(len(bedfile1)) + "\t" + str(len(bedfile2)) + "\t" + str(count) + "\n"))
                myfinalout.close()

print('\n')
print('Distance calculation finished ')
print('The results are in One_to_one_encode_closest_summit_counts_out14.txt')
print('\n')
