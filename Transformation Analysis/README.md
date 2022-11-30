# Transformation Start, Ts Analysis for Dilatometry Curves

This Python code allows the user to both isolate the cooling curve from simply dilatometry quenching data and analyse the transformation start temperatures, Ts using the offset method. The offset method was first proposed by Yang and Bhadeshia for predicting the temperature at which 1% of a transformation has occurred. For more information regarding the offset method see H-S. Yang and H.K.D.H. Bhadeshia, Materials Science and Technology 23 (2007) 556-560. https://www.tandfonline.com/doi/abs/10.1179/174328407X176857. 

The code is split into 2 functions; **dilation_plotter** for plotting the cooling curve, and **dilation_analyser** for analysing the curves. The first function is ideal for observing the cooling curve before the analysis step. The second function determines the cooling transformation start temperature using a linear offset line. The code calculates 3 distinct lines, using 3 different temperature ranges, allowing a range of Ts values to be calculated. These values are then averaged, giving the user a standard deviation of their measurement.

## 1. Inputs

The code requires **7 inputs** in order to run. These are; **file**, **L0**, **r**, **N**, **Ti**, **dT**, and **offset**.

### Input 1 - file

An input for the dilatometry .asc file containing the temperature and change in length measurements.

The data file should be inputted in a Python string formation, where the file name (+ file path) is contained within apostrophes ('). To ensure the smoothness of the code, the file should **ONLY** contain the 'Temperature [°C]' and 'Change in length [μm]' columns - and in that order. Python has difficulty reading the column headers of an .asc file and relying on the header name to read data can cause problems. Instead, the code just reads the 2nd and 3rd columns (i.e., 'Temperature [°C]' and 'Change in length [μm]'). The first column is usually automatically filled with index values when saved as an .asc file and, as such, is ingnored.

An example input for **file** would be as follows:

	file = 'Steel 1_quenched_25.12.21.asc'

for a file called 'Steel 1_quenched_25.12.21.asc' in the same directory as the current working notebook.

### Input 2 - L0

An input for the initial length of the dilatometry sample - used for calculating strain.

The length should be inputted as a number in meters.

An example input for **L0** would be as follows:

	L0 = 0.01

for a sample of initial length 0.01 m (or 10 mm).

### Input 3 - r

An input for the cooling rate used in the quenching test.

The cooling rate should be inputted in degrees Celcius per second (°C/s) and in number format. The cooling rate is not used for calculations, but rather for plotting, and is used to help the user differentiate between samples. 

An example input for **r** would be as follows:
	
	r = 10

for a 10°C/s cooling rate.

### Input 4 - N

An input for the row number/index at which the cooling curve begins.

The row number/index should be inputted as a number and is used to help isolate the cooling curve. Users can look at the raw .asc file to determine this.

An example input for **N** would be as follows:

	N = 12000

for a dataset where the cooling begins on the 12000th measurement.

### Input 5 - Ti

An input for the temperature at which an offset will be measured from.

The **Ti** value should be inputted as a number in degrees Celcius (°C). It is recommended that the user choose a number at least 50°C away from the measured transformation start temperature, Ts - if possible. This might require a bit of trial and error before an appropriate temperature is found. The linear cooling behaviour will then be measured across a 50°C range (50°C above the inputted **Ti**) and a linear relationship determined. 2 more 50°C ranges are also tested. The results from each range are then averaged and outputted.

An example input for **Ti** would be as follows:

	Ti = 750

for a temperature of 750°C - where the linear range taken will be between 750 and 800°C.

### Input 6 - dT

An input for the difference between the 3 temperature ranges tested.

The **dT** value should be inputted as a number in degrees Celcius (°C). The code measures 3 temperature ranges for determing the linear cooling relationship. The first range starts at the **Ti** temperature. The position of the second and third temperature ranges are then determined by the **dT** value - where **dT** is the difference between the ranges. 

An example input for **dT** would be as follows:

	dT = 25

for ranges 25°C apart. For example, if **Ti** was set at 750°C, the next range would start at 775°C, and the third at 800°C.

### Input 7 - offset

An input for the offset strain associated with 1% transformation.

The offset should be inputted as a unitless number and can be calculated from the steel composition using the '**Offset_Calculator.ipynb**' notebook found in this repository.

An example input for **offset** would be as follows:

	offset = 0.00011672

for an offset strain of 0.00011672.


## 2. Running the Code

### dilation_plotter

The code can be run using the command **dilation_plotter(file,L0,r,N)** - where the inputs are defined as discussed.

An example of this would be:

	file = 'Steel 1_quenched_25.12.21.asc'
	L0 = 0.01
	r = 10
	N = 12000

	dilation_plotter(file,L0,r,N)

where the output will be a 'Temperature (°C)' vs 'Strain (mm/mm)' cooling curve.

### dilation_analyser

The code can be run using the command **dilation_analyser(file,L0,r,N,Ti,dT,offset)** - where inputs are defined as discussed.

An example of this would be:

	file = 'Steel 1_quenched_25.12.21.asc'
	L0 = 0.01
	r = 10
	N = 12000
	Ti = 750
	dT = 25
	offset = 0.00011672

	dilation_analyser(file,L0,r,N,Ti,dT,offset)

where the output will be 3 'Temperature (°C)' vs 'Strain (mm/mm)' cooling curves, showing each Ts calculation and their mean average and standard deviation given as a list in the format:
	
	output = [mean, std]