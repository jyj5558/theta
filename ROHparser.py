"""=================================================================================================
    Title: ThetaProject_ROHfile_parser

    This python script will parse the output ROH file from BCFtools/roh in step5 script
    and calculate ROH and F(ROH) of length >

    Jong Yoon Jeon     June 19 2022

================================================================================================="""
#Set working directory, species name using python variable?

#What is the species' reference genome length?
len_ref = 2391065193

"""-------------------------------------------------------------------------------------------------
Do not edit below this line
-------------------------------------------------------------------------------------------------"""

input = open("C:/Users/jyj55/Desktop/ROH_pta_test.txt", 'r')
output = open("C:/Users/jyj55/Desktop/ROH_pta_test_result.txt", 'w')

num_roh_100kb = 0
num_roh_1mb = 0
len_roh_100kb = 0
len_roh_1mb = 0

for line in input:
    if line.startswith("RG"):
        line = line.rstrip('\n')
        field = line.split("\t")
        if float(field[7]) >= 20:
            if float(field[5]) >= 100000 and float(field[5]) < 1000000:
                num_roh_100kb += 1
                len_roh_100kb += float(field[5])
            elif float(field[5]) >= 1000000:
                num_roh_1mb += 1
                len_roh_1mb += float(field[5])
print("num_roh_100kb: " + str(num_roh_100kb) + "\tlen_roh_100kb: " + str(len_roh_100kb) + "\nnum_roh_1mb: " + str(num_roh_1mb) + "\tlen_roh_1mb: " + str(len_roh_1mb))

output.write("Number of ROH >= 100kb: " + str(num_roh_100kb) + "\tLength of ROH >= 100kb: " + str(len_roh_100kb) \
             + "\nNumber of ROH >= 1mb: " + str(num_roh_1mb) + "\tLength of ROH >= 1mb: " + str(len_roh_1mb))
output.close()
input.close()
exit(0)