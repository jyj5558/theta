"""=================================================================================================
    Title: ThetaProject_ROHfile_parser

    This python script will parse the output ROH file from BCFtools/roh in step5 script
    and calculate ROH and F(ROH) of length >

    Jong Yoon Jeon     June 19 2022

================================================================================================="""
#Set working directory, species name using python variable?

#What is the species' reference genome length?
len_ref = 2,391,065,193

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
        field = line.split("\t")
        if field[8] >= 20:
            if field[6] >= 100000 and field[6] < 1000000:
                num_roh_100kb += 1
                len_roh_100kb += field[6]
            elif field[6] >= 1000000:
                num_roh_1mb += 1
                len_roh_1mb += field[6]
print("num_roh_100kb: " + num_roh_100kb + "\tlen_roh_100kb: " + len_roh_100kb "\nnum_roh_1mb: " + num_roh_1mb + "\tlen_roh_1mb: " + len_roh_1mb)

output.write()
output.close()
input.close()
exit(0)