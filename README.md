### Molviz_app - Molecular visualisation web application

#### **Bit of history**

Originally I had an idea of incorporating mols2grid library within Shiny for Python apps, so I worked on a few ideas, but obviously mols2grid was designed to work inside Jupyter Notebook/Lab, so it wasn't directly compatible with the Shiny for Python web application framework. Then I went away to work on another project and left this mini project aside for a while. However, recently I had another idea of trying to build a Shiny for Python app from the scratch instead, so that people can view 2D images of small molecules in a web browser environment instead of only inside a Jupyter Notebook/Lab. So this was where it began.

---

#### **Overall plan**

- Shiny for Python web application for viewing and searching in dataframe (ITables) and viewing 2D images of small molecules of interests
- The data table and image parts have been separated
- The interactive table will now be within the Jupyter notebook environment (Shinylive app embedding method currently not working, but app deployment to other platforms e.g. shinyapps.io should be working)
- Molviz_app is still currently going to be embedded inside Quarto document (to be tested)

---

#### **Current status - end of July 2023**

- Showing at least two PNG images of small molecules in app if using app_image.py (plan to increase to four if possible, so this may form a 4-grid view)
- Saving function of input compounds as PNG files in working directory (reactive effect after pressing confirm button, which will also lead to a direct appearance of compound PNG image in app)
- Open source licence pending