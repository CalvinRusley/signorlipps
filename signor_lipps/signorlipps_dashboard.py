import numpy as np
from numpy import random
import pandas as pd
import holoviews as hv
import bokeh.io
import panel as pn
#pn.extension(sizing_mode="stretch_width", template="fast")
pn.extension(template="fast")
hv.extension('bokeh')

taphonomic_randomness = np.random.randn(30) #generate 30 normally distributed random numbers ahead of time so that you don't get a new set every time you refresh the plots
secondrowplot_width = 600

#use a heatmap to approximate a matplotlib sparsity plot that will behave better with panel
def hacky_sparsity(matrix, name):
    return hv.HeatMap(([str(i) for i in range(1, 31)], [str(i) for i in range(1, 51)], matrix)).opts(invert_yaxis=True, frame_width=250, frame_height=300).opts(cmap='blues').opts(xlabel = 'Taxa', ylabel = 'Beds', title = name)

###initialize each of the control mechanisms for the interactive plots
file_selector = pn.widgets.Select(name='Input file', options=['true_ranges.txt', 'true_ranges2.txt']) #file_selector allows you to toggle which data file you're using
mu_slider = pn.widgets.FloatSlider(name="µ", start=0, end=1, step=0.01, value=0.75) #mu_slider controls the mean likelihood of preservation for the 30 taxa
sigma_slider = pn.widgets.FloatSlider(name="σ", start=0, end=1, step=0.01, value=0.25) #sigma_slider controls the standard deviation in the likelihood of preservation for the 30 taxa
widgets = pn.Column(pn.Spacer(height=30), file_selector, pn.Spacer(height=30), mu_slider, pn.Spacer(height=15), sigma_slider, width=300) #group the controls together

@pn.depends(file_selector.param.value, mu_slider.param.value, sigma_slider.param.value) #declares that the following function will depend on the values of the widgets previously defined
def calc_fossil_ranges(file, mu, sigma):
    true_ranges = np.loadtxt(file, delimiter='\t') #load the chosen file as a matrix
    num_beds, num_taxa = true_ranges.shape #grab the numbers of taxa and beds for later use
    x = random.rand(num_beds,num_taxa) #generate a matrix of the same size as the input data filled with random numbers between 0 and 1
    y = x*true_ranges #multiply the random matrix by the input data, which eliminates the random values outside of the "true ranges" of each of the taxa
    Dtaph = mu + sigma*taphonomic_randomness #using the previously selected random numbers, shift and scale them according to the user-selected mu and sigma
    Dtaph_array = np.tile(Dtaph, (num_beds, 1))*true_ranges #clone the vector into an array equal in size to the input data
    fossil_occurences = np.greater(Dtaph_array,y) #whereever the likelyhood of preservation for a taxon exceeds the corresponding value in the "noise" matrix, a fossil is preserved
    df = pd.DataFrame(fossil_occurences, columns = ["Taxon_"+str(i) for i in range(1, 31)]) #turn the fossil occurences matrix into a dataframe for easier handling
    fossil_ranges = df.applymap(lambda x: 1 if x else np.nan).ffill(axis=0).fillna(0).to_numpy() #fill all deeper beds after the last occurence datum with "true" to reflect the ranges, not the isolated preservations
    return fossil_occurences, fossil_ranges, true_ranges, num_beds, num_taxa, Dtaph

@pn.depends(file_selector.param.value, mu_slider.param.value, sigma_slider.param.value)
def plot_fossil_range(file, mu, sigma):
    fossil_occurences, fossil_ranges, true_ranges, num_beds, num_taxa, Dtaph = calc_fossil_ranges(file, mu, sigma) #compute the stats you'll need
    plot_df = pd.DataFrame({'Fossil diversity' : fossil_ranges.sum(axis=1),'True diversity' : true_ranges.sum(axis=1),'Stratigraphic height (bed #)' : np.arange(num_beds)}) #calculate the diversities and make a dataframe for plotting ease
    return hv.Curve(data=plot_df,kdims=['Stratigraphic height (bed #)','Fossil diversity']).opts(frame_height=300,frame_width=secondrowplot_width,tools=['hover'],invert_xaxis=True) #create the plot

@pn.depends(file_selector.param.value, mu_slider.param.value, sigma_slider.param.value)
def plot_stem(file, mu, sigma):
    fossil_occurences, fossil_ranges, true_ranges, num_beds, num_taxa, Dtaph = calc_fossil_ranges(file, mu, sigma) #compute the stats you'll need
    fossil_ranges[np.where(fossil_ranges==0)]=1 #when fossil_occurences = 0, fossil_ranges = 0, so change the latter spots to 1s so that you don't get a divide-by-zero error
    preservation = fossil_occurences.sum(axis=0)/fossil_ranges.sum(axis=0) #calculate preservation
    taxa_vector = np.arange(num_taxa) #make a vector of a number for each taxon
    plot_df = pd.DataFrame({'Taxa' : np.arange(num_taxa),'Preservation' : preservation, 'Preservability' : Dtaph}) #make a dataframe to ease plotting
    return hv.Spikes(data=plot_df,kdims=['Taxa'])*hv.Points(data=plot_df,kdims=['Taxa','Preservation'], label='Preservation').opts(size=10).opts(ylabel = 'Preservation/Preservability', 
frame_height=300,frame_width=secondrowplot_width,tools=['hover'])*hv.Points(data=plot_df,kdims=['Taxa','Preservability'], label='Preservability').opts(size=10, tools=['hover']) #create a stem plot

@pn.depends(file_selector.param.value, mu_slider.param.value, sigma_slider.param.value)
def plot_fossil_occ(file, mu, sigma):
    fossil_occurences = calc_fossil_ranges(file, mu, sigma)[0] #compute fossil occurences
    return hacky_sparsity(fossil_occurences, 'Fossil Occurences') #make a sparsity plot of the matrix

@pn.depends(file_selector.param.value, mu_slider.param.value, sigma_slider.param.value)
def plot_fossil_ranges(file, mu, sigma):
    fossil_ranges = calc_fossil_ranges(file, mu, sigma)[1] #compute fossil ranges
    return hacky_sparsity(fossil_ranges, 'Fossil Ranges') #make a sparsity plot of the matrix

@pn.depends(file_selector.param.value, mu_slider.param.value, sigma_slider.param.value)
def plot_true_ranges(file, mu, sigma):
    true_ranges = calc_fossil_ranges(file, mu, sigma)[2] #compute true_ranges
    return hacky_sparsity(true_ranges, 'True Ranges') #make a sparsity plot of the matrix

#add the text of the page
first_text = pn.pane.Markdown("""

# Paleontological Evidence for a Bolide Impact:
## An exploration of the Signor–Lipps Effect
 
Extinctions are classified as one of two end-member types: gradual or catastrophic; gradual extinctions are generally attributed to a slow, steady forcing factor driven by Earth processes (e.g., climate change) whereas catastrophic 
extinctions necessarily invoke a single event, such as a bolide impact (Smit and Hertogen, 1980). For the Cretaceous-Paleogene extinction, extinction on both of these time scales (protracted and immediate) are required to depict the full 
range of biological crisis at this time. On the one hand, paleontological data clearly indicate gradual declines in diversity within various invertebrate groups over the last 1 to 20 million years prior to the K-P boundary. For these 
protracted extinctions, stratigraphic analysis suggests that the initiation of the terminal Cretaceous extinctions may have been connected to:

1. Widespread deterioration of global marine environments during large-scale eustatic fall and epicontinental regression.
2. Rapid changes in ocean chemistry and ocean circulation patterns.
3. Rapid temperature fluctuations which occurred during the middle and late Maastrichtian, i.e. last several million years of the Cretaceous (Smit and Hertogen, 1980).
 
On the other hand, there exists substantial paleontological evidence that some groups underwent sudden extinctions precisely at the iridium anomaly and exhibit no pre-extinction gradual demise. When studying the fossil record, we have to 
remember that it is very incomplete. The preservation of a taxon in the fossil record depends on many factors, including the abundance of the organisms, whether it had a skeleton or not, the environment in which it lived, the rate of 
sediment accumulation, the efficiency of scavenging organisms, the chemistry of the sediments and porewaters, and the diagenetic history of the sediments. For this portion of the lab we are going to concentrate on the following question: 
given a probability of preservation that is less than 100% and independent evidence for a mass extinction, what information can we obtain from an incomplete fossil record?

## Simulation Explanation
 
We begin with a known fossil record and ‘degrade’ it to simulate taphonomy and the uneven preservation associated with the fossil record. There are 30 taxa with 50 time points. After 30 time points we impose an extinction event on all 
taxa. Look at the plot of the "True Ranges". This is a 30x50 matrix where ones record a genus presence and zeros record absence. We are going to simulate their preservation as fossils in a sedimentary sequence that includes 50 sedimentary 
beds (increasing with depth, like in a bore hole) that capture a sample of the biota at a given time point with the event of sedimentation.

""")

second_text = pn.pane.Markdown("""

If you'd like, take a look at the code blocks on this page that produce these interactive graphics. Note that the code is commented so you can see the logic behind each step. All good code is clearly commented. Go through each line 
individually and make sure you understand what each step is accomplishing. In this exercise, the preservation of each taxon is the number of beds with fossils divided by the fossil range, and the diversity of each bed is the total number 
of taxa whose range includes that bed.
 
## Simulation Analysis
 
### Question 1:

Change the mean taxon specific taphonomy factor (µ) to 0.55 and the standard deviation (σ) to 0.15 and observe the resulting figures. How do they differ from the default values?

### Question 2: 

Which taxon has the most complete fossil record in your synthetic dataset? Which has the least? Do these correspond to taxa with the highest and lowest preservability? Why might they not? Also, give two examples of extant taxa that have 
good preservation potential and poor preservation potential.

### Question 3: 

It is time to weigh in on the accuracy of the fossil record. Does your plot of diversity over time (assuming constant sedimentation rate for the whole section) reveal a rapid mass extinction of all the taxa between beds 20 and 21? Explain 
why or why not. 

### Question 4: 

Run your script again with a more realistic preservation probability (µ and σ). Hint: recall from an earlier class lecture. What values did you choose? Compare these results with those from Question 11. Soft-bodied organisms are always 
characterized by low preservation potential. But, how might you learn about soft-bodied organisms with no body fossils to study?
 
### Question 5: 

Now switch the input file to true_ranges2.txt, and simulate the fossil record for a mass extinction that occurs in two pulses. What combination of taxon specific taphonomy factors (mu and std) are required to clearly resolve the two pulses 
(i.e. – what is the lower limit on preservation probability to see these features)? 
 
### Summary: 

The combination of an incomplete fossil record coupled with sudden or catastrophic events can lead to what paleontologists call the Signor-Lipps effect (named for its discoverers). When there is a sudden extinction, like the one that 
occurred between beds 20 and 21, the small but finite probability of preservation smears out the pattern, making a sudden extinction look gradual. You can probably imagine a preservational change making a gradual extinction look sudden as 
well. There is a similar effect with the origination of taxa, due to the small probability of preserving the very first organism of a group. This effect is affectionately known as the Sppil-Rongis effect.

""")

row1 = pn.Row(plot_true_ranges, plot_fossil_occ, plot_fossil_ranges, widgets)
row2 = pn.Row(plot_stem, plot_fossil_range)
dashboard = pn.Column(first_text, pn.Spacer(height=30), row1, pn.Spacer(height=35), row2, pn.Spacer(height = 25), second_text)
dashboard.servable()

