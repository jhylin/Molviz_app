### Molviz_app - Molecular visualisation web application

##### **Overall plan**
Shiny in Python web application for viewing and searching in dataframe (ITables) and viewing 2D images of small molecules of interests - embedding two apps inside a Quarto document for public access.

---

##### **Bit of history**
Originally I had an idea of incorporating mols2grid library within Shiny in Python apps, so I worked on a few ideas, but obviously mols2grid was designed to work inside Jupyter Notebook/Lab, so it wasn't directly compatible with the Shiny in Python web app frameworks. Then I went on to work on another project and have left this mini project aside for a while, and then recently I had another idea of trying to build a Shiny in Python app from the scratch so that people can view 2D images of small molecules in a web browser, instead of only inside a Jupyter Notebook etc. 

---

##### **Current status**
- Working on app server/output part - testing different presentation options for images and data tables
- Data input currently not tested yet as there may be a few different choices, but my current idea sits more with this - the data input will be restricted to allow selections of a few different small-molecule datasets for now, once the app itself is more stable and I've made myself more familiar with the data/web security side, then I may look at allowing uploads of .csv files (with possibly also other file formats e.g. SDFs etc.)
- I'll look into which open source license is most suitable for this project hopefully soon (if it works out as planned) so that there may be more helpers for it