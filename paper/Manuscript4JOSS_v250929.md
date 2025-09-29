**cationCalc4EPMA: Efficient processing tool for EPMA dataset**

Kazuki Matsuyama<sup>1,2, \*</sup>, Yumiko Harigane<sup>3</sup>, and Yoshihiro Nakamura<sup>3</sup>

<sup>1</sup>Department of Earth and Planetary Sciences, Graduate School of Environmental Studies, Nagoya University, Nagoya, 464-8601, Japan

<sup>2</sup>Géosciences Montpellier, Université de Montpellier, Montpellier, 34095, France

<sup>3</sup>Research Institute of Geology and Geoinformation, Geological Survey of Japan, National Institute of Advanced Industrial Science and Technology, Tsukuba, 305-8567, Japan

\*Corresponding Author

Email: <matsuyama.kazuki.i8@s.mail.nagoya-u.ac.jp>

**Summary**

The cationCalc4EPMA is a freely available MATLAB® tool designed to process electron probe microanalysis (EPMA) data by converting oxide weight percentages into cation-based formulas. This conversion is crucial for petrological tasks such as determining structural formulas, estimating Fe³⁺/Fe²⁺ ratios, recognizing solid-solution series, and performing thermobarometric modeling. The tool features automatic mineral identification using empirical composition thresholds and includes built-in functions for calculating commonly used petrological indices, like the Mg# (Mg / Mg + Fe) for olivine. Additionally, it provides batch averaging based on sample codes to streamline data interpretation. By replacing manual or proprietary workflows, cationCalc4EPMA improves reproducibility, transparency, and accessibility in the geoscientific analysis of EPMA datasets.

**Statement of Need**

Electron probe microanalysis (EPMA) is a contemporary technique grounded in the physical principles of electron-induced X-ray emission and wavelength-dispersive spectroscopy (WDS), initially developed by Raymond Castaing during his doctoral research in Paris (Castaing, 1951). This method stands out as one of the most effective microanalytical approaches for conducting precise, non-destructive, and quantitative elemental analysis of solid substances, such as minerals, at the micrometer level (Yang et al., 2022). As a result, EPMA is extensively utilized and has profoundly influenced geoscientific studies (e.g., Sweatman and Long, 1969; Wiedenbeck et al., 2007; Kouketsu et al., 2014; Basch et al., 2025). The compositional information derived from EPMA is crucial for interpreting pressure-temperature-time (_P-T-t)_ histories (e.g., Ozawa, 2004), reconstructing fluid-rock interaction processes (e.g., Whattam et al., 2022), and assessing geodynamic conditions within the Earth's interior (e.g., Chardelin et al., 2024).

To extract meaningful petrological insights from EPMA data, it is crucial to convert weight percent values into cation-based formulas, which are typically standardized to a specific number of oxygen atoms or cations, depending on the mineral group. This transformation is essential for calculating structural formulas, evaluating Fe³⁺/Fe²⁺ ratios, identifying solid-solution series, and performing thermobarometric modeling. Ensuring the reproducibility of this conversion method is particularly crucial, as several minerals, including amphibole (Holland and Blundy, 1994; Leake et al., 1997; Putirka, 2016), mica (Rieder, 1998), and garnet (Grew et al., 2013), are named based on the calculated cation distribution. Despite its importance, many researchers continue to rely on outdated spreadsheet macros, manual calculations, or proprietary software, all of which are prone to human error and difficult to verify or replicate.

To address this issue, we have created cationCalc4EPMA, a cation calculator for EPMA datasets, designed for MATLAB®. This tool imports .csv files produced by EPMA and transforms weight percent data into cations, adhering to open science principles to enhance reproducibility and transparency. The calculation procedures in cationCalc4EPMA are primally based on those of Deer et al. (2013), except for the Fe<sup>3+</sup> calculation.

To estimate Fe<sup>3+</sup>, the charge balance method or site allocation method is employed. This approach is tailored for specific minerals: pyroxenes (Papike et al., 1974), spinel and inverse spinel (Droop, 1987), amphibole (Holland and Blundy, 1994), epidote (Masumoto et al., 2014), and garnet (Enami, 2012). Nonetheless, if verifying the exact amount of Fe³⁺ proves challenging or if its calculated presence might influence the end-member ratio, users have the option to assume Fe³⁺ = 0 by activating the commented-out code.

A key feature of this calculator is its ability to automatically identify the mineral phase from the measured data using empirical mass percentages and subsequently compute the cations. Additionally, the calculator is composed of one module and one .m file, enabling users to easily modify the code to suit their specific requirements, such as for minerals with an extensive list or for calculations tailored to particular data based on location. This design caters to a broader spectrum of users who rely on EPMA analysis.

**Methods**

cationCalc4EPMA processes EPMA datasets using the following steps:

1\. Importing Excel files:

Users choose a raw data file (.csv or .xlsx) from a list dialog. The module Cation_moduli.xlsx, which works alongside the main script, is loaded automatically. This module includes two sheets: one that holds the molar weights, cation numbers, oxygen atom counts, and names of oxides/anions; and another that contains stoichiometric details for target minerals. To add more minerals to the list, users need to update this module accordingly.

2\. Creating a mineral list:

Users select the mineral phases anticipated to be present in the raw dataset to create a list of a list of these mineral phases. The potential candidates for this list can be found in the Cation_moduli.xlsx module. This initial list comprises olivine, orthopyroxene (Opx), clinopyroxene (Cpx), spinel, plagioclase, ilmenite, magnetite, quartz, amphibole, chlorite, epidote, biotite, garnet, and apatite, which are typically found in metamorphic and ultramafic rocks. This list plays a crucial role later in the process for automatic mineral identification (Step 9), so it is essential to include the minerals that have been measured.

3\. Pre-processing Excel files:

The script is designed to work directly with the EPMA output files. Initial processing steps address any missing data and ready the information for further analysis.

4\. Selecting data for cation calculation:

Users have the option to choose which elements to factor into the calculation. The script presumes that the sum of the selected columns (wt%) will be close to 100%. If there are any repeated element names, the process stops, and an error is displayed.

5\. Filtering by Detection Limit (DL):

If the file contains DL values, any measurements that fall below the DL are automatically omitted.

6\. Renaming dataset variables:

Column (variable) names are renamed to improve readability.

7\. Handling missing elements:

For elements (e.g., P, V, S) that may not be included in every analysis, NaN-filled dummy columns are added to prevent errors during calculation.

8\. Creating a total wt% column:

The total of the selected wt% columns is calculated for each row.

9\. Identifying mineral phases:

Mineral phases are determined automatically based on empirical composition thresholds. Calculations are then performed for each identified phase. Users must ensure that this classification matches the actual mineral assemblage. Furthermore, since this function is constructed using the if statement, the sequence in which identification is executed is crucial. In the initial script, minerals with distinct compositions are given precedence, followed by the evaluation of silicates. When users introduce new mineral phases, it is advisable for them to determine the "identification" sequence by consulting the commented-out sections within the script.

10\. Splitting data by mineral phase:

The dataset is split into separate tables by assigned mineral names.

11\. Calculating elemental cation numbers:

The cation (or anion) numbers for each element are calculated by dividing the measured data (wt%) by the molar weight.

12\. Calculating oxygen molar amounts:

Oxygen molar quantities are determined using oxide weight percentages and module data. The results are compiled into a single Excel file, organized by mineral.

13\. Normalizing elemental cation numbers:

By utilizing the oxygen molar values and module data, the cation numbers for each element are determined according to the mineral's theoretical formula. Additionally, petrologically significant indices, such as Mg# for olivine and Almandine% for garnet, are computed.

14\. Outputting cation calculation results:

Results are exported to an Excel file, with separate sheets for each mineral phase.

15\. Averaging values by sample code:

EPMA data often uses sample identifiers like EPMA01_01, EPMA01_02, and so forth, which consist of the rock sample name (sample code) followed by a serial number. The script identifies common prefixes and calculates average values for each mineral associated with each sample code. This approach provides a useful estimate of typical compositions, though users should be aware that zoning or variations in composition within grains could affect interpretation.

**Package Summary**

cationCalc4EPMA enables:

1\. Processing of EPMA data:

Facilitates conversion of raw wt% data to cation/anion molar values.

2\. Automatic mineral identification:

Classifies mineral phases based on empirical mass% composition thresholds. These criteria can be modified by the user. The results are output in separate sheets per mineral.

3\. Calculation of common petrological indices:

Calculates commonly used indices, such as Mg# and almandine%, directly from cation data.

4\. Sample-based averaging:

Automatically averages values by mineral phase and sample code prefix.

**Acknowledgement**

This study was carried out using the electron probe microanalysis (EPMA; JXA-iPF200H) at GSJ-Lab in AIST, Japan. The authors wish to thank Hiroaki Koge, for recommending this submission to JOSS. We also thank Tomoaki Morishita, Yui Kouketsu, Katsuyoshi Michibayashi, and Benoit Ildefonse for their support during this study. This study was supported by grants from JST SPRING (JPMJSP2125), the Sasakawa Scientific Research Grant from The Japan Science Society (2024-6011), and the Japan Society for the Promotion of Science, Japan (JP25KJ1408).

**References**

Basch, V., M. Godard, A. Tommasi, and E. Rampone (2024), Melt/rock ratios and melt fluxes during reactive percolation: from matrix- to melt-controlled dynamics, Contrib Mineral Petr, 180(1), doi:10.1007/s00410-024-02194-1.

Castaing, R. (1951), Application of electron probes to local chemical and crystallographic analysis, University of Paris.

Chardelin, M., A. Tommasi, and J. A. Padron-Navarta (2024), Progressive Strain Localization and Fluid Focusing in Mantle Shear Zones during Rifting: Petrostructural Constraints from the Zabargad Peridotites, Red Sea, J Petrol, 65(8), doi:ARTN egae081

10.1093/petrology/egae081.

Deer, W. A., R. A. Howie, and J. Zussman (2013), An Introduction to the Rock-Forming Minerals, doi:10.1180/dhz.

Droop, G. T. R. (2018), A general equation for estimating Fe3+ concentrations in ferromagnesian silicates and oxides from microprobe analyses, using stoichiometric criteria, Mineralogical Magazine, 51(361), 431-435, doi:10.1180/minmag.1987.051.361.10.

Enami, M. (2012), Influence of garnet hosts on the Raman spectra of quartz inclusions, Journal of Mineralogical and Petrological Sciences, 107(4), 173-180, doi:10.2465/jmps.111216.

Grew, E. S., A. J. Locock, S. J. Mills, I. O. Galuskina, E. V. Galuskin, and U. Halenius (2013), Nomenclature of the garnet supergroup, American Mineralogist, 98(4), 785-811, doi:10.2138/am.2013.4201.

Holland, T., and J. Blundy (1994), Non-ideal interactions in calcic amphiboles and their bearing on amphibole-plagioclase thermometry, Contrib Mineral Petr, 116(4), 433-447, doi:10.1007/bf00310910.

Kouketsu, Y., M. Enami, T. Mouri, M. Okamura, and T. Sakurai (2014), Composite metamorphic history recorded in garnet porphyroblasts of Sambagawa metasediments in the Besshi region, central Shikoku, Southwest Japan, Island Arc, 23(4), 263-280, doi:10.1111/iar.12075.

Leake, B. E., et al. (2018), Nomenclature of Amphiboles; Report of the Subcommittee on Amphiboles of the International Mineralogical Association Commission on New Minerals and Mineral Names, Mineralogical Magazine, 61(405), 295-310, doi:10.1180/minmag.1997.061.405.13.

Masumoto, Y., M. Enami, M. Tsuboi, and M. Hong (2014), Magmatic zoisite and epidote in tonalite of the Ryoke belt, central Japan, European Journal of Mineralogy, 26(2), 279-291, doi:10.1127/0935-1221/2014/0026-2360.

Ozawa, K. (2004), Thermal History of the Horoman Peridotite Complex: A Record of Thermal Perturbation in the Lithospheric Mantle, J Petrol, 45(2), 253-273, doi:10.1093/petrology/egg110.

Papike, J. J., Cameron, K.L. and Baldwin, K. (1974), Amphiboles and Pyroxenes: Characterization of Other than Quadrilateral Components and Estimates of Ferric Iron from Microprobe Data, in Geological Society of America, edited, pp. 1053-1054.

Putirka, K. (2016), Amphibole thermometers and barometers for igneous systems and some implications for eruption mechanisms of felsic magmas at arc volcanoes, American Mineralogist, 101(4), 841-858, doi:10.2138/am-2016-5506.

Rieder, M., et al. (2024), Nomenclature of the Micas, Clays and Clay Minerals, 46(5), 586-595, doi:10.1346/ccmn.1998.0460513.

Sweatman, T. R., and J. V. P. Long (1969), Quantitative Electron-probe Microanalysis of Rock-forming Minerals, J Petrol, 10(2), 332-379, doi:10.1093/petrology/10.2.332.

Whattam, S. A., J. C. M. De Hoog, M. I. Leybourne, and M. Z. Khedr (2022), Link between melt-impregnation and metamorphism of Atlantis Massif peridotite (IODP Expedition 357), Contrib Mineral Petr, 177(11), doi:10.1007/s00410-022-01968-9.

Wiedenbeck, M., et al. (2007), Further Characterisation of the 91500 Zircon Crystal, Geostandards and Geoanalytical Research, 28(1), 9-39, doi:10.1111/j.1751-908X.2004.tb01041.x.

Yang, S.-Y. (2022), Electron Probe Microanalysis In Geosciences: Analytical Procedures And Recent Advances, Atomic Spectroscopy, 43(1), doi:10.46770/as.2021.912.