# EngTools

A collection of engineering, database, and miscellaneous dataframe tools that include methods to calculate bioenergy primary conversion reactor yields.

## Installation

Use the following to install.

    pip install git+https://github.com/dasievers/engtools

## General Data Tools

Miscellaneous conversion tools are available that do things like converting date-time indices to elapsed time, clipping a dataframe with date-time bookends, smoothing. Also available is a simple timing function to evaluate code snippets.

### dtrange

A function that takes a Pandas dataframe or series and returns only the portion of interest, specified by start/end date-times.

### df_smooth

A function that takes a Pandas dataframe or series and returns a smoothed version of the data based on user parameters. Smoothing is done using a moving window convolution, similar to a moving average.

### indexconvert

A function that takes a Pandas dataframe or series and returns a version with index converted from datetime64[ns] to float values. Useful for converting to e.g. runtime instead of absolute datetime. Various options are available including removing large gaps in time.

## Chemical Engineering Tools

### SatSteam & Water

Classes that enable looking up saturated steam properties and also density and viscosity of liquid water.

### henry_constant

Calculate Henry constant for a given gas and temperature.

### Hydrogen

A class that calculates the Z-factor of hydrogen gas given temperature and pressure.

### adiabatic_process

Function that models an adiabatic process of an ideal gas.

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
