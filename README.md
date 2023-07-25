### Molviz_app - Molecular visualisation web application

#### **Bit of history**

Originally I had an idea of incorporating mols2grid library within Shiny for Python web app framework (after seeing an example of a similar app in Streamlit previously), so I worked on a few ideas, but obviously mols2grid was designed to work inside Jupyter Notebook/Lab and Shiny for Python was only out of alpha at that stage so things were still being developed, and unfortunately mols2grid wasn't directly compatible with the Shiny for Python at that time. I then went away to work on another project on molecular scaffolds and left this mini project aside. However, recently I had another idea of trying to build a Shiny for Python app from the scratch (my first ever app relating to cheminformatics or chemical information), so that users in relevant fields can view 2D images of small molecules in a web browser environment instead of only inside a Jupyter Notebook/Lab, with planned functions of viewing, saving and substructure highlighting of these molecules. So this was where it began. 

Note: I'm sure this is not the first ever app for this type of use as there are many great apps out there already, but I thought to try out the Shiny for Python framework - to put it to test and show that it should work in the chemistry world.


#### **Overall plan**

- Shiny for Python web application for viewing dataframe and 2D images of small molecules of interests - planned functions: viewing, saving and substructure highlighting of target molecules
- Blog post regarding to this work may also show another version of the interactive table (using ITables) in the Jupyter notebook environment (Shinylive app embedding method currently not working for this due to my personal preference of keeping the index column, but app deployment to other platforms should be working - not tested yet though)
- Molviz_app may be deployed via Shinyapps.io or HuggingFace - to be decided (not to be embedded inside Quarto doc as RDKit does not have pure wheel file since it's not written purely in Python, it's written using Python and C++, therefore Shinylive version will not work)


#### **Current status - end of July 2023**

- Showing individual PNG images of selected compounds
- Showing 4 PNG images of selected compounds in table grid format if running app_image.py
- Highlighting atoms and bonds via position number (or based on RDKit's MolToImage())
- Saving function of selected compounds with highlighted atoms & bonds for all 4 compounds in the table grid as PNG image in working directory

<br>

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.