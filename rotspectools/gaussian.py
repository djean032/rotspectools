import re

file = open('MP2_6-311+G(2d,p)_(3-Cyano)Methylenecyclopropane_anharmonic_freq.txt', 'r')

lines = file.readlines()

# regex to find line where it states "Sextic Centrifugal Distortion Constants"
pattern = re.compile(r'Sextic Centrifugal Distortion Constants')

for idx, line in enumerate(lines):
    if pattern.search(line):
        line_number_ref = idx
        break

print(line_number_ref)


rot_const_list = lines[line_number_ref - 5:line_number_ref - 2]

rot_const = []
for line in rot_const_list:
    # regex pattern to return everything from A00 to A0 in line
    pattern = re.compile(r'\D00\s*=\s*(\d+(?:\.\d+)?)\s*\D0')
    match = re.search(pattern, line)
    rot_const.append(match.group(1))



quart_s_red = lines[line_number_ref - 47:line_number_ref - 42]

s_quart_const = []
for line in quart_s_red:
    pattern = re.compile(r'-?\d+\.\d+(?:[dD]-?\d+)?')
    matches = pattern.findall(line)
    s_quart_const.append(matches[1])


sext_s_red = lines[line_number_ref + 46 :line_number_ref + 53]

s_sext_const = []
for line in sext_s_red:
    pattern = re.compile(r'-?\d+\.\d+(?:[dD]-?\d+)?')
    matches = pattern.findall(line)
    s_sext_const.append(matches[1])
print(s_sext_const)




