#### Molviz_app - Molecular visualisation web application

##### **Bit of history**
Originally I had an idea of incorporating mols2grid library within Shiny in Python apps, so I worked on a few ideas, but obviously mols2grid was designed to work inside Jupyter Notebook/Lab, so it wasn't directly compatible with the Shiny in Python web app frameworks. Then I went on to work on another project and have left this mini project aside for a while, and then recently I had another idea of trying to build a Shiny in Python app from the scratch so that people can view 2D images of small molecules in a web browser, instead of only inside a Jupyter Notebook/Lab, with also a saving function available to save the compound images.

---

##### **Overall plan**
- Shiny for Python web application for viewing and searching in dataframe (ITables) and viewing 2D images of small molecules of interests
- The data table and image parts have been separated
- The interactive table will now be within the Jupyter notebook environment (app embedding method currently not working, but app deployment to other platforms should be working)
- Molviz app is still currently going to be embedded inside Quarto document (to be tested)

---

##### **Current status - end of July 2023**
- Showing at least two PNG images of small molecules in app if using app_image.py (plan to up to four if possible, so this may form the 4-grid view)
- Saving function of input compounds as PNG files in working directory (reactive effect after pressing confirm button, which will also lead to a direct appearance of compound PNG image in app)
- Open source licence pending