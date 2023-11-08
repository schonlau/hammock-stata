# stata-hammock
Hammock Plot - Stata implementation

# Hammock plot


## Description

    The hammock plot draws a graph to visualize categorical or mixed categorical / continuous data.
    Variables are lined up parallel to the vertical axis. Categories within a variable are spread out along a
    vertical line. Categories of adjacent variables are connected by boxes. (The boxes are parallelograms; we
    use boxes for brevity). The "width" of a box is proportional to the number of observations that correspond
    to that box (i.e. have the same values/categories for the two variables). The "width" of a box refers to the
    distance between the longer set of parallel lines rather than the vertical distance.

    If the boxes degenerate to a single line, and no labels or missing values are used the hammock plot
    corresponds to a parallel coordinate plot. Boxes degenerate into a single line if barwidth is so small that
    the boxes for categorical variables appear to be a single line. For continuous variables boxes will usually
    appear to be a single line because each category typically only contains one observation.

    The order of variables in varlist determines the order of variables in the graph.  All variables in varlist
    must be numerical. String variables should be converted to numerical variables first, e.g. using encode or
    destring.

---

[Getting started](#Getting-started) | [Syntax](#Syntax) | [Examples](#Examples) | [Historical Context](#Historical-Context) | [References](#References)

---


## Getting started

Install hammock from either `ssc` or `GitHub` (may be more recent):

SSC (**v1.25**):

```shell
ssc install hammock, replace
```

GitHub (**v1.25**):

```
net install hammock, from("https://raw.githubusercontent.com/schonlau/hammock-stata/main/installation/") replace
```

If you have the hammock package installed, you can check for updates: `ado update, update`.


## Syntax 

The syntax is as follows:

```
hammock varlist [if] [in],
               [
                Missing BARwidth(real 1) MINBARfreq(int 1) /// 
                hivar(str) HIVALues(numlist  missingokay) SPAce(real 0.0) ///
                LABel labelopt(str) label_min_dist(real 3.0) ///
                SAMEscale(varlist) ///
                ASPECTratio(real 0.72727) COLorlist(str) shape(str) no_outline *
               ]
```

See the help file `help hammock` for details.

The most basic use is as follows:

```
hammock varlist
```

where `varlist` are is the list of variables to visualize in that order.


## Examples

You may want to change the scheme, the overall look of the graph. I prefer a plain white look:

```
set scheme s1color
```

Load the Shakespeare data set to your current directory:
```stata
net get hammock, from("https://raw.githubusercontent.com/schonlau/hammock-stata/main/installation/") replace
```

Each observation in this data set represents one of Shakespeare's plays. 
speaker1 and speaker2 refer to the social status of the first two speakers to speak at the beginning of the play. characters is the number of different persons in the play.
In this first example, we add labels and change the background color to a grey tone (gs5) to make the labels more readable:
```
hammock  type characters speaker1 speaker2 sex1 sex2, label color(gs5)
```
<img src="figures/grey.pdf" alt="Basic Hammock plot" width="600"/>

The ordering of the child-adolescent-adult variable is not in the desired order; adult should not be in the middle. We now specify a specific order, child-adolescent-adult. 

```python
var = ["hospitalizations","group","gender","comorbidities"]
group_dict= {1: "child", 2: "adolescent",3: "adult"}
value_order = {"group": group_dict}
hammock = hammock_plot.Hammock(data_df = df)
ax = hammock.plot(var=var, value_order=value_order )
```

<!--- to restrict image size, I am using a an html command, rather than the standard ![](image.png) --->
<!---    ![Hammock plot ](image/asthma1.png)   --->
<img src="image/asthma_value_order.png" alt="Hammock plot" width="600"/>

We highlight observations with comorbidities=0  in red:

```python
ax = hammock.plot(var=var, value_order=value_order ,hi_var="comorbidities", hi_value=[0], color=["red"])
```

<!---   ![Hammock plot with highlighting](image/asthma_highlighting.png)    --->
<img src="image/asthma_highlighting.png" alt="Hammock plot with highlighting" width="600"/>


### Example Satisfaction scales for the diabetes data

We import the diabetes dataset:

```python
import hammock_plot
import pandas as pd
df = pd.read_csv('../examples/diabetes_outlier/diabetes_for_python.csv')
```

The three variables represent different ordinal scales for satisfaction. We are checking for missing values: 
```python
var = ["sataces","satcomm","satrate"]
hammock = hammock_plot.Hammock(data_df = df)
ax = hammock.plot(var=var,  default_color="blue", missing=True) 
```

<img src="image/diabetes.png" alt="Hammock plot for the Diabetes Data" width="600"/>

The missing value category is shown at the bottom for each variable. We find missing values for all 3 variables, but fewest for the last one. We also see a phenomenon called "top coding", where 
satisfied respondents simply choose the highest value.




## Historical context

In 1898, Sankey diagrams were developed to visualize flows of energy and materials. 

In 1985, Inselberg popularized parallel coordinates to visualize continuous variables only. The central contribution is the use of parallel axes.

In 2003, Schonlau proposed the hammock plot. This was the first plot to visualize categorical data (or mixed categorical continuous data) on parallel axes. 

In 2010, Rosvall proposed alluvial plots to visualize network variables over time. Rather than using bars to connect axes, alluvial plots use rounded curves. Alluvial plots are now also used to visualize categorical data.

There are several additional variations that also visualize categorical data including Parallel Set plots (Bendix et al, 2005), Right Angle plots (Hofmann and Vendettuoli, 2013),
and generalized parallel coordinate plots (GPCPs) (popularized by VanderPlas et al., 2023). 

### References 
Bendix, F., Kosara, R., & Hauser, H. (2005). Parallel sets: visual analysis of categorical data. In IEEE Symposium on Information Visualization, 2005. INFOVIS 2005. 133-140. 

Hofmann, H., & Vendettuoli, M. (2013). Common angle plots as perception-true visualizations of categorical associations. IEEE transactions on visualization and computer graphics, 19(12), 2297-2305.

Inselberg, A., & Dimsdale, B. (2009). Parallel coordinates. Human-Machine Interactive Systems, 199-233.

Rosvall, Martin, & Bergstrom, C.T. (2010) "Mapping change in large networks." PloS one 5.1: e8694.

Sankey, H. (1898). Introductory note on the thermal efficiency of steam-engines. report of
the committee appointed on the 31st march, 1896, to consider and report to the council
upon the subject of the definition of a standard or standards of thermal efficiency for
steam-engines: With an introductory note. In Minutes of proceedings of the institution
of civil engineers, Volume 134, pp. 278–283.

Schonlau M. 
*[Visualizing Categorical Data Arising in the Health Sciences Using Hammock Plots.](http://www.schonlau.net/publication/03jsm_hammockplot.pdf)* 
In Proceedings of the Section on Statistical Graphics, American Statistical Association; 2003

VanderPlas, S., Ge, Y., Unwin, A., & Hofmann, H. (2023). 
Penguins Go Parallel: a grammar of graphics framework for generalized parallel coordinate plots. 
Journal of Computational and Graphical Statistics, 1-16. (online first)

### Other implementations of the hammock plot 
There is also a Python implementation `hammock_plot` and an R implementation as part of the package `ggparallel` (The R implementation does not currently allow for quantitative variables).



