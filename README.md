### Molviz_app - Molecular visualisation web application

#### **Bit of history**

Originally I had an idea of incorporating mols2grid library within Shiny for Python apps, so I worked on a few ideas, but obviously mols2grid was designed to work inside Jupyter Notebook/Lab, so it wasn't directly compatible with the Shiny for Python web application framework. Then I went away to work on another project and left this mini project aside for a while. However, recently I had another idea of trying to build a Shiny for Python app from the scratch instead, so that people can view 2D images of small molecules in a web browser environment instead of only inside a Jupyter Notebook/Lab. So this was where it began.


#### **Overall plan**

- Shiny for Python web application for viewing and searching in dataframe (ITables) and viewing 2D images of small molecules of interests
- The data table and image parts have been separated
- The interactive table will now be within the Jupyter notebook environment (Shinylive app embedding method currently not working, but app deployment to other platforms e.g. shinyapps.io should be working)
- Molviz_app will be deployed via Shinyapps.io (not as embedding inside Quarto doc as RDKit does not have pure wheel file since it's not written purely in Python, with also C++)


#### **Current status - end of July 2023**

- Showing individual PNG images of selected compounds side-by-side
- Showing 4 PNG images of selected compounds in table grid format if running app_image.py
- Saving function of selected compounds as PNG files in working directory for both individual selected compounds and all 4 together in a table grid (reactive effect after pressing confirm button, which will also lead to a direct appearance of compound PNG image in app)


<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.