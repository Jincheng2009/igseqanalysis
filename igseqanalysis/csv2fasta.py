# -*- coding: utf-8 -*-
#==============================================================================
#     Copyright (C) 2017  MedImmune, LLC
#     
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#==============================================================================

"""
Created on Thu Jan 05 14:58:59 2017

@author: wuji
"""
import sys
import getopt

def main():
    argv = sys.argv[1:]
    idx = 5
    readFromFile = False
    try:
        opts, args = getopt.getopt(argv,"hi:p:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "-i":
            csvfile = arg
        elif opt == "-p":
            idx = int(arg)

    if readFromFile:
        filein = open(csvfile)
    else:
        filein = sys.stdin

    for line in filein:
        line = line.rstrip()
        line = line.replace("\t",",")
        tokens = line.split(",")
        if len(tokens) < 6 or len(tokens[5]) == 0:
            continue
        fastaid = tokens[0]
        germline = tokens[1] + "+" + tokens[2]
        cdr3 = tokens[idx]
        cdr3 = cdr3.replace("-","")
        if len(cdr3) > 1:
            sys.stdout.write(">" + fastaid + "||" + germline + "\n")
            sys.stdout.write(cdr3 + "\n")

def usage():
    print 'cat cdr.csv | python csv2fasta.py -p 5 full.csv'    

if __name__ == "__main__":
    main()

    