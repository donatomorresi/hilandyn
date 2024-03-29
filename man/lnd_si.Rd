% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lnd_si.R
\docType{data}
\name{lnd_si}
\alias{lnd_si}
\title{Landsat time series of spectral indices.}
\format{
A \code{PackedSpatRaster} object with 108 layers, three per year.
}
\usage{
lnd_si
}
\description{
A raster dataset containing annual spectral indices from 1985 to 2020.
Spectral indices are:
the Moisture Stress Index \insertCite{@MSI; @hunt1989detection}{hilandyn};
the Tasseled cap Wetness \insertCite{@TCW; @crist1985tm}{hilandyn};
the Tasseled cap Angle \insertCite{@TCA; @powell2010quantification}{hilandyn}.
The spatial subset covers a wildfire occurred in the Aosta Valley (Italy) during the winter of 2005.
Landsat imagery was acquired by the Thematic Mapper (TM), Enhanced Thematic Mapper + (ETM+) and Operational Land Imager (OLI)
sensors. Surface reflectance was derived from Level 1 Tier 1 imagery from Landsat Collection 1 in the USGS archive using the
Framework for Operational Radiometric Correction for Environmental monitoring (FORCE) software \insertCite{frantz2019force}{hilandyn}.
Annual reflectance composites were produced using all the images acquired during the growing season (between June 1 and September 30).
The geometric median approach \insertCite{roberts2017high}{hilandyn} was employed to generate reflectance composites at the pixel level.
}
\references{
\insertRef{crist1985tm}{hilandyn}

\insertRef{hunt1989detection}{hilandyn}

\insertRef{powell2010quantification}{hilandyn}
}
\keyword{datasets}
