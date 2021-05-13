import json

fname = "test.txt"
nsims = 74000

file = open(fname, "r")
replaced_content = ""
old_content = ""

for line in file.readlines():
    if 'nsims' in line:
        nsims_temp = line[-6:]
        print(nsims_temp)
        nsims_old = float(nsims_temp)
        line = line.strip()
        new_line = line.replace(line, "nsims = %s" % nsims)
        replaced_content = replaced_content + new_line + "\n"
    else:
        old_content = old_content + line

content = old_content + replaced_content

file.close()
write_file = open("test.txt", "w")
write_file.write(content)
write_file.close()
