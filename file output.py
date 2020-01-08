#Output binding effect data to tab-separated txt file for visualisation
output_file = open("output.txt", "w")
output_string = ""
for mz, i in binding_effect.items():
    output_string += str(mz) + "\t" + str(i) + "\n"
output_file.write(output_string)
output_file.close()