# CCT Predictor for Low Alloy Steels

## 1. Inputs

The code requires **5 inputs** in order to run. These are; **comp**, **G**, **rates**, **title**, and **savedata**.

### Input 1 - comp

An input for the steel alloy composition.

Alloy composition should be inputted in a Python dictonary format, with curly brackets ({}) used to define the dictionary space. Dictionaries are made up of key:value pairs, where keys and values are separated by a colon (:), and a comma (,) is used to separate key:value pairs. In this case, alloying element symbols are inputted as the keys and their associated concentration in wt.% inputted as the values. 

An example **comp** input would be as follows:

    comp = {'C':0.1,'Si':0.2,'Mn':0.3,'Ni':0.4,'Cr':0.5,'Mo':0.6}
    
for an Fe-0.1C-0.2Si-0.3Mn-0.4Ni-0.5Cr-0.6Mo alloy.

### Input 2 - G

An input for the austenite grain size. 

Austenite grain size should be inputted as the ASTM grain size number. Both a conversion table and equations for calculating this value from SI units can be found in ASTM E112: Standard Test Methods for Determining Average Grain Size.

An example input for **G** would be as follows:

    G = 10
    
for ASTM grain size 10.

### Input 3 - rates

An input for the cooling rates to be tested. 

Cooling rates are to be inputted in a Python list format, with square brackets ([]) used to define the list space. Lists are made up of values separated by a comma (,). In this case, the desired cooling rates should be inputted in degrees per second (°C/s).

An example input for **rates** would be as follows:

    rates = [0.1, 1, 10, 100]
    
for cooling rates 0.1, 1, 10, and 100°C/s.

### Input 4 - title

An input for the name/title of the alloy.

The **title** input is for labelling and file saving purposes only, where files will be saved using the inputted **title**. The **title** should be inputted in a Python string format, where the title is contained within apostrophes (').

An example input for **title** would be as follows:
    
    title = 'Steel 1'
    
for a steel titled Steel 1.

### Input 5 - savedata

An input for deciding whether the raw data and CCT image should be saved locally.

The **savedata** input should be set as either **True** or **False**, depending on whether the user wants to save the raw data and CCT image to their local files. If **True**, a directory will be created within the current working directory (cwd) called **title** CCT data. Raw data files for each cooling rate will be saved as .csv files and the CCT image as a .png.

Example inputs for **savedata** would be as follows:

    savedata = True
    
if the user wants to save their data locally, or:

    savedata = False
    
if the user does not want to save their data locally.

## 2. Running the Code

The code can be run using the command **CCT_Calculator(comp,G,rates,title,savedata)** - where inputs are defined as discussed.

An example of this would be:

    comp = {'C':0.1,'Si':0.2,'Mn':0.3,'Ni':0.4,'Cr':0.5,'Mo':0.6}
    G = 10
    rates = [0.1, 1, 10, 100]
    title = 'Steel 1'
    savedata = True
    
    CCT_Calculator(comp,G,rates)