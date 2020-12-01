# engtools
A collection of engineering and database tools that include methods to calculate bioenergy primary conversion reactor yields.

## Chemical Engineering Tools
### SatSteam
A saturated steam state object that is used to deliver corresponsing temperature, pressure, enthalpy, etc.

### WaterViscosity
Viscosity of water at given temperature.

### WaterDensity
Density of water at given temperature.

### henry_constant
Calculate Henry constant for a given gas and temperature.

## Database Tools
### dbupdate_sqlite
This function is built with some presets for convenience, but it really provides a general process for creating and updating local sqlite databases from a growing directory of text files (.csv, .tsv) from automated data loggers, such as OptoDisplay. Instead of the time-consuming process of reading in all files and then filtering the results to a specific query, a database can be maintained and queried as needed much more efficiently.

### query_sqlite
The partner of dbupdate_sqlite, this function runs efficient queries on the created databases to return tags and date ranges as requested. Additional filters may also be passed as sql query language.

### read_data and batch_read
Tools for reading in .csv, .tsv, and excel files and parsing them appropriately. These are fairly simple, but provide batch capability and also are able to handle *null* characters in OptoDisplay files.

## Process Calcuations
### fis
Calculate slurry fraction of insoluble solids via no-wash method.

### Pretreater
Class for computing pretreatment reactor results. Pre-flash total solids can be estimated thermodynamically, and conversion yields can be calculated and plotted. Utilizes the `total solids conservation` method.

### Enzymer
Class for computing enzymatic hydrolysis results. Can select from full formal mass and mole closure model or the total solids conservation method, similar to the Pretreater, which reduces the need for analytical data.
