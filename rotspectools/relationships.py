from uncertainties import ufloat
from uncertainties.umath import *
from scipy import constants as c

# constants

amu_const = c.physical_constants['atomic mass constant'][0]
amu_unc = c.physical_constants['atomic mass constant'][2]
amu_conv = ufloat(amu_const, amu_unc)
m_angstrom = ufloat(1e10, 0)
planck_const = c.physical_constants['Planck constant'][0]
planck_unc = c.physical_constants['Planck constant'][2]
planck = ufloat(planck_const, planck_unc)
inertia_moments_conv = (planck / (8 * c.pi**2)) * (m_angstrom**2 / amu_conv) / 1e6


def parse_rotational_constant(constants: str) -> list[float]:
    split_list = constants.split('(')
    split_list[1] = split_list[1].strip(')')
    decimals = len(split_list[0].split('.')[1])
    zero_pad = decimals - len(split_list[1]) - 1
    split_list[1] = "." + "0" * zero_pad + split_list[1]
    parsed_list = [float(i) for i in split_list]
    return parsed_list

def kappa(A, B, C):
    return (2 * B - A - C) / (A - C)

def inertia_moments(A, B, C):
    Ia = inertia_moments_conv / A
    Ib = inertia_moments_conv / B
    Ic = inertia_moments_conv / C
    return Ia, Ib, Ic

def inertial_defect(Ia, Ib, Ic):
    return Ic - Ib - Ia

def Paa(Ia, Ib, Ic):
    return (Ib + Ic - Ia) / 2

def Pbb(Ia, Ib, Ic):
    return (Ia + Ic - Ib) / 2

def Pcc(Ia, Ib, Ic):
    return (Ia + Ib - Ic) / 2

no_of_lines = 3
lines = ""
for i in range(no_of_lines):
    lines += input("Please enter the rotational constants, using the same format as your .res file: ") + "\n"

lines = lines.split('\n')
print(lines)
a_input = lines[0]
b_input = lines[1]
c_input = lines[2]

a_list = parse_rotational_constant(a_input)
b_list = parse_rotational_constant(b_input)
c_list = parse_rotational_constant(c_input)

A = ufloat(a_list[0], a_list[1])
B = ufloat(b_list[0], b_list[1])
C = ufloat(c_list[0], c_list[1])
print(A)
print(B)
print(C)



k = kappa(A, B, C)
print('Kappa: ')
print(k)

Ia, Ib, Ic = inertia_moments(A, B, C)

for i in [Ia, Ib, Ic]:
    print(i)

delta_i = inertial_defect(Ia, Ib, Ic)
print(delta_i)

paa = Paa(Ia, Ib, Ic)
pbb = Pbb(Ia, Ib, Ic)
pcc = Pcc(Ia, Ib, Ic)

print(paa)
print(pbb)
print(pcc)



# kappa = (2 * B - A - C) / (A - C)
