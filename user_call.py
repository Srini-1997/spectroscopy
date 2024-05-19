path = input("Enter the path: ")
with open('user_input.txt', 'w') as f:
    f.write(path)

bias_list = input("Enter the bias files: ")
with open('user_input.txt', 'a') as f:
    f.write("\n"+bias_list)

flat_list = input("Enter the flat frames: ")
with open("user_input.txt", 'a') as f:
    f.write("\n"+flat_list)


