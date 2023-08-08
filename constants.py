# FRAGMENTER PARAMS
DEFAULT_PEPTIDE = '(164.0700)FDSFGDLSSASAIM(16)GNPK'
DEFAULT_MIN_CHARGE = 1
DEFAULT_MAX_CHARGE = 2
DEFAULT_MASS_TYPE = 'monoisotopic'
DEFAULT_FRAGMENT_TYPES = 'by'

MIN_PEPTIDE_CHARGE = 1
MAX_PEPTIDE_CHARGE = 5
MAX_PEPTIDE_AA_COUNT = 50
MAX_PEPTIDE_LENGTH = 1000

# SPECTRA PARAMS
DEFAULT_TOLERANCE_TH = 0.5
DEFAULT_TOLERANCE_PPM = 50.0
DEFAULT_MIN_INTENSITY = 50.0

TOLERANCE_OPTIONS = ['th', 'ppm']
DEFAULT_TOLERANCE_TYPE = 'th'
DEFAULT_TOLERANCE_TYPE_INDEX = TOLERANCE_OPTIONS.index(DEFAULT_TOLERANCE_TYPE)
TOLERANCE_STEP_PPM = 1.0
TOLERANCE_STEP_TH = 0.1

MAX_TOLERANCE_VALUE_TH = 1.0
MAX_TOLERANCE_VALUE_PPM = 1000.0
MIN_TOLERANCE_VALUE = 0.0


WIKI = """

# Peptide Fragmentation

## Introduction

Peptide fragmentation refers to the process by which peptides (short chains of amino acid monomers linked by peptide bonds) are broken into smaller fragments. This occurs during mass spectrometry analysis to provide useful information about the peptide's sequence, structure, and identity. The process uses methods such as collision-induced dissociation (CID), electron-transfer dissociation (ETD), or higher-energy collisional dissociation (HCD).

## Basic Concept

Before we delve into peptide fragmentation, let's understand a peptide molecule:

A peptide molecule consists of multiple amino acids linked together by peptide bonds. A peptide bond is formed between the carboxyl group (-COOH) of one amino acid and the amino group (-NH2) of another.

## Fragmentation Process

When peptides are subjected to fragmentation in mass spectrometry, the peptide bonds are targeted. The energy imparted to the molecule causes these bonds to break, creating different fragment ions. These fragments are typically classified into six types based on the location of the bond break: a, b, and c ions are the N-terminal fragments, and x, y, and z ions are the C-terminal fragments.


### Example abcxyz Fragment Ions

For example, consider the following peptide sequence: PEPTIDE. It can prodice the the following theoretical fragment ion sequences.
```
a/b/c ions (start from the front, just like the alphabet)
1 - P
2 - PE
3 - PEP
4 - PEPT
5 - PEPTI
6 - PEPTID
7 - PEPTIDE

x/y/z ions (start from the back, just like the alphabet)
1 - E
2 - DE
3 - IDE
4 - TIDE
5 - PTIDE
6 - EPTIDE
7 - PEPTIDE
```

## Collision-Induced Dissociation (CID)

CID is the most commonly used method of peptide fragmentation. In CID, peptides are accelerated in an electric field to give them kinetic energy, and then they collide with a neutral gas molecule, causing them to fragment.

The major fragment ions produced in CID are b and y ions. In general, CID is excellent for generating sequence information for peptides and works well for a broad range of peptides.

## Electron-Transfer Dissociation (ETD)

ETD is a fragmentation method that can preserve post-translational modifications, making it valuable for the study of protein modifications. In ETD, peptides are fragmented by transferring an electron to the peptide, which causes it to break apart.

ETD commonly results in c and z ions. This method is especially useful for longer and more charged peptides.

## Higher-Energy Collisional Dissociation (HCD)

HCD is a hybrid method that uses CID-type fragmentation, but with higher energy. This method generates a broader distribution of ion types and can lead to more complete sequence coverage.

Like CID, HCD commonly results in b and y ions, but with more secondary fragmentation.

## Analysis

The resulting fragment ions are then analyzed by the mass spectrometer. By examining the m/z (mass-to-charge) ratio of the fragment ions, researchers can deduce the amino acid sequence of the original peptide. This process is crucial in protein identification and characterization.



## Understanding Internal Fragment Ions

These ions originate when fragmentation transpires at two separate points along the peptide, yielding a fragment that is 'internal' to the peptide's sequence. This process entails the breaking of not one, but two peptide bonds, leading to the creation of a peptide fragment that isn't connected to either the N-terminus or C-terminus of the original peptide.

### Example Internal Fragmentation

Let's take the peptide sequence 'PEPTIDE' as an example. Here is its 'b' ion series, signifying that all fragment ions originate at the N-terminus:
```
b1 - P
b2 - PE
b3 - PEP
b4 - PEPT
b5 - PEPTI
b6 - PEPTID
b7 - PEPTIDE
```
In comparison, the internal fragment ions for the same peptide sequence 'PEPTIDE' are shown below. The internal fragment ions for each 'b' ion are denoted in parentheses. Notice the shift in the pattern of fragmentation, now requiring two fragmentations to occur to produce the internal fragment ions.
```
b1 - P ()
b2 - PE (E)
b3 - PEP (EP, E)
b4 - PEPT (EPT, PT, T)
b5 - PEPTI (EPTI, PTI, TI, I)
b6 - PEPTID (EPTID, PTID, TID, ID, D)
b7 - PEPTIDE (EPTIDE, PTIDE, TIDE, IDE, DE, E)
```

"""

HELP = """
# Help

The **Peptide Fragmenter** application is a tool for in silico fragmentation of peptide sequences. The app takes an amino acid sequence and calculates the fragment ions for a given charge range. 

## Input

- **Peptide Sequence**: This is the amino acid sequence to be fragmented. The sequence can include modification masses in parentheses. For example, (-3.14)PEP(123.456)TIDE contains a -3.14 N-Term modification and a 123.456 modification on the second Proline (P). The length of the peptide sequence cannot exceed 50 amino acids.

- **Charge Range**: Define the minimum and maximum charges for fragmentation. The minimum charge must be less than or equal to the maximum charge. 

- **Mass Type**: Choose between 'monoisotopic' or 'average' mass type to use for fragment calculation.

- **Fragment Types**: Select the types of fragments to calculate: 'a', 'b', 'c', 'x', 'y', 'z'. 

- **Internal Fragments**: Check this box if you want to include internal fragments in the calculation. 

## Outputs

- **Results Tab**: This tab presents the calculated fragment ions in a table, and a plot that shows the fragment segments in a sequence versus mass plot. You can download the fragment ion data as a CSV file.

- **Spectra Tab**: This tab presents an input form to add spectra to match fragment ions to. The format for spectra is: {m/z} {intensity}, one per line. There is also an option to adjust the offset and its type (Da or ppm), and the minimum intensity. If spectra data is provided, the tab displays a plot with matching ions marked and a table of the spectra data. This data can also be downloaded as a CSV file.

- **Wiki Tab**: This tab presents a wiki page with general information on peptide fragmentation.

- **Help Tab**: This tab presents a help page on how to use the application (You're currently here).

If you encounter any issues or have suggestions for improvement, please contact pgarrett@scripps.edu.

This is a work in progress and your feedback is greatly appreciated!
"""