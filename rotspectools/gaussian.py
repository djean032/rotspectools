import re
from tkinter import filedialog as fd
from typing import Tuple

def get_file() -> Tuple[int, list[str]]:
    filename = fd.askopenfilename()
    print(filename)
    file = open(filename, "r")

    lines = file.readlines()

    # regex to find line where it states "Sextic Centrifugal Distortion Constants"
    pattern = re.compile(r"Sextic Centrifugal Distortion Constants")
    match_line = 0
    for idx, line in enumerate(lines):
        if pattern.search(line):
            match_line = idx
            break
    return match_line, lines


def get_rot_const(line_number_ref: int, lines: list[str]) -> list[str]:
    rot_const_list = lines[line_number_ref - 5 : line_number_ref - 2]

    rot_const = []
    for line in rot_const_list:
        # regex pattern to return everything from A00 to A0 in line
        pattern = re.compile(r"\D00\s*=\s*(\d+(?:\.\d+)?)\s*\D0")
        match = re.search(pattern, line)
        rot_const.append(match.group(1))
    return rot_const


def get_quart_s_const(line_number_ref: int, lines: list[str]) -> list[str]:
    quart_s_red = lines[line_number_ref - 47 : line_number_ref - 42]

    s_quart_const = []
    for line in quart_s_red:
        pattern = re.compile(r"-?\d+\.\d+(?:[dD]-?\d+)?")
        matches = pattern.findall(line)
        matches[1] = matches[1].replace("D", "E")
        s_quart_const.append(matches[1])
    return s_quart_const


def get_sext_s_const(line_number_ref: int, lines: list[str]) -> list[str]:
    sext_s_red = lines[line_number_ref + 46 : line_number_ref + 53]

    s_sext_const = []
    for line in sext_s_red:
        pattern = re.compile(r"-?\d+\.\d+(?:[dD]-?\d+)?")
        matches = pattern.findall(line)
        matches[1] = matches[1].replace("D", "E")
        s_sext_const.append(matches[1])
    return s_sext_const

def get_quart_a_const(line_number_ref: int, lines: list[str]) -> list[str]:
    quart_a_red = lines[line_number_ref - 74 : line_number_ref - 69]

    a_quart_const = []
    for line in quart_a_red:
        pattern = re.compile(r"-?\d+\.\d+(?:[dD]-?\d+)?")
        matches = pattern.findall(line)
        matches[1] = matches[1].replace("D", "E")
        a_quart_const.append(matches[1])
    return a_quart_const


def get_sext_a_const(line_number_ref: int, lines: list[str]) -> list[str]:
    sext_a_red = lines[line_number_ref + 30 : line_number_ref + 37]

    a_sext_const = []
    for line in sext_a_red:
        pattern = re.compile(r"-?\d+\.\d+(?:[dD]-?\d+)?")
        matches = pattern.findall(line)
        matches[1] = matches[1].replace("D", "E")
        a_sext_const.append(matches[1])
    return a_sext_const


ref_line, lines = get_file()
rot_con = get_rot_const(ref_line, lines)
s_quart_const = get_quart_s_const(ref_line, lines)
s_sext_const = get_sext_s_const(ref_line, lines)
a_quart_const = get_quart_a_const(ref_line, lines)
s_sext_const = get_sext_a_const(ref_line, lines)

rot_con = [float(i) for i in rot_con]
s_quart_const = [float(i) for i in s_quart_const]
s_sext_const = [float(i) for i in s_sext_const]
a_quart_const = [float(i) for i in a_quart_const]
a_sext_const = [float(i) for i in s_sext_const]

print("Rotational Constants: " + str(rot_con))
print("Quartic S Constants: " + str(s_quart_const))
print("Sextic S Constants: " + str(s_sext_const))
print("Quartic A Constants: " + str(a_quart_const))
print("Sextic A Constants: " + str(s_sext_const))
