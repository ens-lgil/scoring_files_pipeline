import os
import argparse

def main():

    argparser = argparse.ArgumentParser()
    argparser.add_argument("--loc_file", help='Path to the (new) variant location output file', required=True, metavar='VAR_FILE')
    argparser.add_argument("--coord_file", help='Path to the variants file with coordinates', required=True, metavar='VAR_COORD_FILE')

    args = argparser.parse_args()

    loc_file = args.loc_file
    coord_file = args.coord_file

    if not os.path.isfile(loc_file):
        print("File '"+loc_file+"' can't be found")
        exit(1)

    if not os.path.isfile(coord_file):
        print("File '"+coord_file+"' can't be found")
        exit(1)



    coord_file_handler = open(coord_file, 'a')
    loc_file_handler = open(loc_file, 'r')
    coord_file_handler.write(loc_file_handler.read())
    coord_file_handler.close()


if __name__ == '__main__':
    main()