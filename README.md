#### Molviz_app - Molecular visualisation web application

##### **Overall plan**
Shiny in Python web application for viewing and searching in dataframe (ITables) and viewing 2D images of small molecules of interests - embedding two apps inside a Quarto document for public access.

---

##### **Bit of history**
Originally I had an idea of incorporating mols2grid library within Shiny in Python apps, so I worked on a few ideas, but obviously mols2grid was designed to work inside Jupyter Notebook/Lab, so it wasn't directly compatible with the Shiny in Python web app frameworks. Then I went on to work on another project and have left this mini project aside for a while, and then recently I had another idea of trying to build a Shiny in Python app from the scratch so that people can view 2D images of small molecules in a web browser, instead of only inside a Jupyter Notebook etc. 

---

##### **Current status**
- Managed to show at least two PNG images of small molecules in app, with saving function as PNG files as automatic/reactive response from input
- open source license pending