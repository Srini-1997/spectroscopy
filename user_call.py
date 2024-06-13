path = input("Enter the path: ")
with open('user_input.txt', 'w') as f:
    f.write(path)

bias_list = input("Enter the bias files: ")
with open('user_input.txt', 'a') as f:
    f.write("\n"+bias_list)

flat_list = input("Enter the flat frames: ")
with open("user_input.txt", 'a') as f:
    f.write("\n"+flat_list)

src = input("Enter the name of source frame: ")
with open("user_input.txt", "a") as f:
    f.write("\n"+src)

std = input("Enter the name of standard star frame: ")
with open("user_input.txt", "a") as f:
    f.write("\n"+std)


while True:
    arc = input("Do you have separate arc lamp for source and standard star(y/n): ")
    if arc=='y':
        src_arc_lamp = input("Enter the arc lamp for source")
        std_arc_lamp = input("Enter the arc lamp for standard")
        with open("user_input.txt", "a") as f:
            f.write("\n"+src_arc_lamp)
            f.write("\n"+std_arc_lamp)
        break
    elif arc=='n':
        src_arc_lamp = std_arc_lamp = input("Enter the lamp")
        with open("user_input.txt", "a") as f:
            f.write("\n"+src_arc_lamp)
            f.write("\n"+std_arc_lamp)
        break
    else:
        print("Wrong Choice")
        continue

readnoise = input("Enter the read noise of the CCD")
with open("user_input.txt", 'a') as f:
    f.write("\n"+readnoise)

gain = input("Enter the gain of the CCD")
with open("user_input.txt", 'a') as f:
    f.write("\n"+gain)
