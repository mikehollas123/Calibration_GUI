
import os
input_folder = r"C:\Data\I2MS\500plus True Positives\Ping processed"

dirs = os.listdir(input_folder)

for filename in dirs:
    with open(os.path.join(input_folder, filename),"r") as openFile:
        line = openFile.readline()
        line = openFile.readline()
        with open(os.path.join(input_folder, f"Extracted_Stori_{filename.split('.')[0]}.csv"),"w") as writeFile:
            writeFile.write(f"Scan\tIonnumber\tSegnumber\tM/Z\tFrequency\tSlope\tSlope R Squared\tTime of Birth\tTime of Death\n")
            while line:



                values = line.split('\t')
                scan = values[0]
                MZ = values[1]
                Frequency = values[2]
                Slope = values[3]
                SlopeR2 = values[4]
                TOD = values[5]
                writeFile.write(f"{scan}\t0\t0\t{MZ}\t{Frequency}\t{Slope}\t{SlopeR2}\t0\t{TOD}")
                line = openFile.readline()