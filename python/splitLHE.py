# from https://github.com/cjohnson-phys/Python_Scripts

#!/usr/bin/python

import sys, getopt, os

def printUsage():
    print 'Usage: splitLHE.py -i <input LHE file> -n <number of events per output file>\nPlease try again'

def main(argv):

    inputfile = ''
    nEvt = 10

    # Take user arguments for input file to split and number of events per file
    # Do some basic error checking to see if args are there
    try:
        opts, args = getopt.getopt(argv,"hi:n:",["ifile=","num="])
    except getopt.GetoptError:
        printUsage()
        sys.exit(2)

    if ( len(sys.argv) == 1 ):
        printUsage()
        sys.exit(2)
    
    for opt, arg in opts:
      if opt == '-h':
        printUsage()
        sys.exit()
      elif opt in ("-i", "--ifile"):
        inputfile = arg
      elif opt in ("-n", "--num"):
        nEvt = arg

    print 'Input file is', inputfile
    print 'Number of events per file:', nEvt

    # Root for output files, becomes baseout_fileNum.lhe
    baseout = inputfile.replace(".lhe","")

    # Now open the input LHE file (testing to see if it exists), proceed to split up
    try:
        with open(inputfile,"r") as infile: 
            eventNum = 0  # Count number of events
            fileNum = 0   # Index for output files. Numbered so baseout_0.lhe has all the init stuff in it (we'll get rid of it at the end)

            outFile = open(baseout+"_"+str(fileNum)+".lhe","w")  # Open first output file to stick all init stuff in
            commonBlock = '' # Store init info etc

            # Loop over each line
            while 1:
                line = infile.readline()
            
                if not line:
                    break
                else :

                    # Get common block
                    if ( ( eventNum == 0 ) and not( "<event>" in line ) ):
                        commonBlock += line
                    else: pass
                        
                    # Decide if this line represents start of a new event    
                    if( "<event>" in line ):
                        if ( (eventNum) % int(nEvt) == 0 ):                          # If we've gone through nEvt events, time for a new file. 
                            outFile.write("</LesHouchesEvents>")                     # End file correctly 
                            outFile.close()                                          # Close old output lhe file
                            fileNum += 1
                            fileNum_str = "%05d" % fileNum
                            outFile = open(baseout+"_"+fileNum_str+".lhe","w")      # Open new output lhe file with increased index. Should prob catch exceptions...
                            outFile.write(commonBlock)                              # Add common block to start of each file
                            print "File #:", fileNum
                        eventNum += 1

                    outFile.write(line)                    
            
            print eventNum, "events in the original LHE file"
            os.remove(baseout+"_0.lhe") # Leftover with baseout_0.lhe with just init information. Delete this! 
            
            infile.close()
            outFile.close()

    except IOError as e:
       print 'File does not exist!'
       sys.exit(2)

if __name__ == "__main__":
   main(sys.argv[1:])
