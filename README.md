# Steel Dilatometry Analysis Code

Python code for analysing steel dilatometer curves. Functions for the code are found in "dilatometry_ts_calculator.py" and include: 'dilatometry_curve_plotter' for plotting a temperature-strain curve, 'offset_method' for measuring transformation start temperatures, $T_\mathrm{s}$, using an offset method, and 'second_deriv_method' which is a second technique for measuring $T_\mathrm{s}$ by finding the inflection point in a curve.

## 1. Functions

### Dilatometry Curve Plotter

A function for plotting the measured dilatometry temperature vs. strain curve. Output is a temperature-strain graph of the dilatometry data.

    dilatometry_curve_plotter(file,L0,N)
  
**Inputs:**

    file, L0, N
    
**file**: raw dilatometry data in .asc file format. Raw data should be saved with "Temperature" in the first column and "Change in Length" in the second column (ignoring the index column). File name (+ datapath to file) should be inputted as a string. 

**L0**: intial length of the dilatometer sample in meters, m. Initial length should be inputted as a float.

**N**: row index number to start plotting the data from. Ideally used to isolate cooling data (i.e., the tail end of the dataset). Index number should be inputted as a integer.

**Example:**

Fig. 1 shows an example output for a bainitic steel dilatometry cooling curve. The original sample length, L0, was 0.01 m and the cooling curve data started at index value, N, 12000.

Example inputs:

    file = '210423_SteelDilatometryData_test1.asc'
    
    L0 = 0.01
    
    N = 12000
    
Example code:
    
    dilatometry_curve_plotter(file='210423_SteelDilatometryData_test1.asc', L0=0.01, N=12000)
    
Example output:

<figure>
  <img
  src="example figures/dilatometry_curve_plotter_EXAMPLE.png"
  alt="."
  width="80%" 
  height="80%">
  <figcaption>Fig. 1 An example dilatometry curve output by the 'dilatometry_curve_plotter' function. In this example, N is set to a value at the start of the cooling curve as to isolate this section of the data.</figcaption>
</figure>

&nbsp;
    
### Offset Method

A function for analysing the dilatometer cooling curve transformations using an offset method, similar to the technique described by Yang and Bhadeshia [1]. A total of three separate temperature ranges are chosen to measure three different $T_\mathrm{s}$ values. Output is three temperature-strain curves for each temperature range tested, plus a Python list giving an average $T_\mathrm{s}$ and standard deviation separated by a comma (i.e., [Ts, std]).

    offset_method(file,L0,Nc,Ti,dT,comp,p,X)

**Inputs:**

    file, L0, Nc, Ti, dT, comp, p, X
    
**file**: raw .asc dilatometer data (as detailed above).

**L0**: original length of dilatometer sample (as detailed above).

**Nc**: row index value for the start of the cooling curve. If data before cooling is used then this could disrupt the function and its ability to calculate $T_mathrm{s}$. N should be inputted as an integer.
    
**Ti**: the lower temperature value used for calculating the gradient of a linear section of the dilatometer curve. A gradient will be taken over a 50&deg;C range. It is recommended this value is at least 50&deg;C higher than the measured $T_\mathrm{s}$ value. Ti should be inputted as a float/integer.

**dT**: the temperature difference between the ranges of temperatures used to calculate gradients. Three ranges of temperatures are used to calculate three separate gradient values to then measured three separate $T_\mathrm{s}$ values. dT is the temperature difference between these equally spaced ranges. dT should be inputted as a float/integer.

**comp**: the composition of the steel alloy in wt.%. This input is used to calculate the exact offset induced by a transformation to either ferrite, 'f', or martensite, 'm', in austenite. Composition should be inpuuted in dictionary format.

**p**: the phase/constituent forming in austenite. This should be either ferrite, 'f', or martensite, 'm'. p should be inputted as a string.

**X**: the mole fraction of phase/constituent forming in austenite. X should be inputted as a float.
    
**Explanation of the Offset Method:**

The function calculates $T_\mathrm{s}$ using an offset method. This technique first calculates the gradient of the linear section of the dilatometry curve (at least 50&deg;C before the cooling transformation). A 50&deg;C range is selected for this. This gradient is then used to plot an offset line, where $T_\mathrm{s}$ is measured at its interception with the cooling curve. To reduce uncertainty, three 50&deg;C temperature ranges are used and an average $T_\mathrm{s}$ is found. The position of these temperature ranges is dictated by the **Ti** and **dT** inputs. The **Ti** input dictates where the first 50&deg;C range will start from, whereas **dT** dictates the temperature distance between these ranges. Fig. 2 shows a visualisation of this.

<figure>
  <img
  src="example figures/selecting_gradient_ranges.png"
  alt="."
  width="50%" 
  height="50%">
  <figcaption>Fig. 2 A schematic showing how the 'offset_method' function takes the inputs Ti and dT to select temperature ranges for gradient calculation.</figcaption>
</figure>

&nbsp;

Once a gradient has been measured, an offset line is defined. The amount this line is offset from the original curve is dependent on the temperature of transformation, the composition of the steel, **comp**, the phase/constituent that is forming in the austenite, **p**, and the mole fraction of this new phase/constituent, **X**. The code will use these inputs to calculate the exact amount of strain induced by **X** fraction of phase/constituent, **p**, within the austenite with composition, **comp**, at a specific temperature. [Note: the code does not require the input of temperature here as it will already be considered.] If the transformation is ferritic, pearlitic or bainitic, **p** = 'f'. If the transformation is martensitic, **p** = 'm'. 

**Example:**

Fig. 3 shows an example output for the 'offset_method' function. The technique was used on the same data plotted in Fig. 1. All inputs are true except the input for comp which has been altered for simplicity.

Example inputs:

        file = '210423_SteelDilatometryData_test1.asc'
        
        L0 = 0.01
        
        Nc = 12000
        
        Ti = 550
        
        dT = 50
        
        comp = {'C': 0.2, 'Si': 0.1, 'Mn': 0.3}
        
        p = 'f'
        
        X = 0.01
        
Example code:

        offset_method(file='210423_SteelDilatometryData_test1.asc', L0=0.01, Nc=12000, Ti=550, dT=50, comp={'C': 0.2, 'Si': 0.1, 'Mn': 0.3}, p='f', X=0.01)
        
Example output:

        [463.5, 1.7]

<figure>
  <img
  src="example figures/offset_method_EXAMPLE.png"
  alt="."
  width="100%" 
  height="100%">
  <figcaption>Fig. 3 An example output by the 'offset_method' function.</figcaption>
</figure>

&nbsp;

### Second Derivative Method

A function for analysing the dilatometer cooling curve transformations by finding the inflection point of a curve. This technique is ideal when the offset method cannot be used - typically during a double transformation where there is no linear section before the 2nd transition. The output is a 1d gaussian filtered temperature-strain curve showing the measured $T_\mathrm{s}$ position, a temperature-2nd derivative curve, and a Python list giving an average $T_\mathrm{s}$ and standard deviation separated by a comma (i.e., [Ts, std]).

    second_deriv_method(file,L0,Nc,s,Trange)

**Inputs:**

    file, L0, Nc, s, Trange
    
**file**: raw .asc dilatometer data (as detailed above).

**L0**: original length of dilatometer sample (as detailed above).

**Nc**: row index number for the start of the cooling curve (as detailed above).

**s**: sigma value for the 1d gaussian filter (i.e., curve smoothing factor). The raw dilatometry data needs to be smoothed in order to better isolate inflection points within the curve. Higher values of sigma will increase the intensity of data smoothing. Sigma should be inputted as a float/integer.

**Trange**: the temperature range in which to look for an inflection point. Dilatometry curves will contained a multitude of inflection points so users are asked to provide a temperature range in which to look for a specific inflection point (i.e., at the $T_\mathrm{s}$). The temperature range should be inputted as a Python list with a minimum temperature, T1, and maximum temperature, T2, separated by a comma - as so [T1, T2].

**Example:**

Fig. 4 shows an example output for the 'second_deriv_method' function.

Example inputs:

        file = '210423_SteelDilatometryData_test1.asc'
        
        L0 = 0.01
        
        Nc = 12000
        
        s = 3
        
        Trange = [620, 670]
        
Example code:

        second_deriv_method(file='210423_SteelDilatometryData_test1.asc', L0=0.01, Nc=12000, s=3, Trange=[620,670])
        
Example output:

        [645.5, 0.5]
        
<figure>
  <img
  src="example figures/second_deriv_method_EXAMPLE.png"
  alt="."
  width="75%" 
  height="75%">
  <figcaption>Fig. 4 An example output by the 'second_deriv_method' function.</figcaption>
</figure>

## 2. Bibliography

[1]:  H-S. Yang and H. K. D. H. Bhadeshia. Uncertainties in dilatometric determination of martesite start temperature. Materials Science and Technology, 23:556–560, 2007

## Citing the code:

[![DOI](https://zenodo.org/badge/356209778.svg)](https://zenodo.org/badge/latestdoi/356209778)
