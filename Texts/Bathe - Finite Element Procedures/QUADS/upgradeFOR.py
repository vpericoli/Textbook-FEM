FOR_file = "QUADS.FOR"

lines = open(FOR_file, "r").readlines()
print lines
new_FOR_file = "QUADS.FOR"
new          = open(new_FOR_file, "w")
for line in lines:
    #print line[72:]
    #print line[:72]
    new.write(line[:72] + "\n")

new.close()
