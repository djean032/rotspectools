a_input = input("Please enter the A rotational constant, using the same format as your .res file: ")
b_input = input("Please enter the B rotational constant, using the same format as your .res file: ")
c_input = input("Please enter the C rotational constant, using the same format as your .res file: ")

a_tmp = '7780.57901(39)'
b_tmp = '2566.81361(12)'
c_tmp = '2132.05216(13)'


a_list = a_tmp.split('(')
for i in a_list:
    print(i)


# kappa = (2 * B - A - C) / (A - C)
