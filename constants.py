import os


def get_env_int(var_name, default):
    return int(os.getenv(var_name, default))


def get_env_float(var_name, default):
    return float(os.getenv(var_name, default))


def get_env_str(var_name, default):
    return os.getenv(var_name, default)


# FRAGMENTER PARAMS
DEFAULT_PEPTIDE = '[164.0700]-FDSFGDLSSASAIM[16]GNPK'
DEFAULT_MIN_CHARGE = 1
DEFAULT_MAX_CHARGE = 2
DEFAULT_MASS_TYPE = 'monoisotopic'
DEFAULT_FRAGMENT_TYPES = 'a;b;c;x;y;z'

MIN_PEPTIDE_CHARGE = 0
MAX_PEPTIDE_CHARGE = get_env_int('MAX_PEPTIDE_CHARGE', 100)
MAX_PEPTIDE_AA_COUNT = get_env_int('MAX_PEPTIDE_AA_COUNT', 500)
MAX_PEPTIDE_LENGTH = get_env_int('MAX_PEPTIDE_LENGTH', 2000)

# SPECTRA PARAMS
DEFAULT_TOLERANCE_TH = 0.5
DEFAULT_TOLERANCE_PPM = 30.0
DEFAULT_MIN_INTENSITY = 0.0

TOLERANCE_OPTIONS = ['th', 'ppm']
DEFAULT_TOLERANCE_TYPE = 'ppm'
DEFAULT_TOLERANCE_TYPE_INDEX = TOLERANCE_OPTIONS.index(DEFAULT_TOLERANCE_TYPE)
TOLERANCE_STEP_PPM = 1.0
TOLERANCE_STEP_TH = 0.1

MAX_TOLERANCE_VALUE_TH = get_env_float('MAX_TOLERANCE_VALUE_TH', 1.0)
MAX_TOLERANCE_VALUE_PPM = get_env_float('MAX_TOLERANCE_VALUE_PPM', 1000.0)
MIN_TOLERANCE_VALUE = 0.0

BASE_URL = get_env_str('BASE_URL', 'http://localhost:8502')

WIKI = """
# Peptide Fragmentation  
  
### Introduction  
  
Peptide fragmentation refers to the process by which peptides are broken into smaller fragments. This process occurs within the mass spectrometer and provides useful information regarding the peptide's sequence and structure. Peptide fragmentation can be induced by a number of methods, such as: collision-induced dissociation (CID), electron-transfer dissociation (ETD), and/or higher-energy collisional dissociation (HCD). These techniques effectively add energy to peptides, inducing their dissociation into smaller fragments.

The fragmentation patterns observed are remarkably consistent, yet pose significant challenges for accurate modeling. Such patterns are influenced by multiple factors, including the applied collision energy, the peptide's charge state, the specific collision method employed, and the composition and sequence of the peptide's amino acids.

Despite the complexity inherent in predicting fragmentation patterns, advancements in deep learning have markedly improved these predictions, achieving near-perfect accuracy in some cases. This success highlights deep learning's potential in capturing the complex interdependencies and nuanced dynamics governing peptide fragmentation, significantly outperforming traditional algorithmic prediction methods.
  
### Peptides
  
A peptide is a chain of amino acids connected by peptide bonds, which are formed between the carboxyl group (-COOH) of one amino acid and the amino group (-NH2) of another. These bonds are established via a dehydration synthesis, a process that releases water molecules during bond formation.  Its also important to highlight that these peptide bonds are the weakest bonds within the peptide, and as such require the least amount of energy to break.
  
![image](https://biologydictionary.net/wp-content/uploads/2017/01/Peptide-Bond-Formation.jpg)  

## Fragmentation
  
When peptides undergo fragmentation in mass spectrometry, the process targets specific bonds within the molecule, varying by the chosen fragmentation method. The energy applied causes the peptide to break into smaller fragments. The type of ions produced from this fragmentation depends on where the break occurs within the peptide and which part of the molecule retains the charge. If the charge remains on the N-terminal side of the peptide, the resulting ions are classified as a, b, or c ions. Conversely, if the charge is on the C-terminal side, the ions are identified as x, y, or z ions.

### Terminal Fragment Ions

As you can see in the figure below, each fragmentation site produces two complementary ions, which when added together make up the composition of the original peptide sequence (This is only true for the neutral fragments seen below).
    
**The fragment sites for all terminal fragment ions:**

![image](https://www.matrixscience.com/images/cleavages.gif)  

The peptide bond is notably the weakest and, therefore, the most likely to break during fragmentation. This breakage predominantly results in the formation of b and y ions. These ions are especially prevalent in Collision-Induced Dissociation (CID) and Higher-energy Collisional Dissociation (HCD), which are classified as soft ionization techniques. In contrast, employing higher energy fragmentation methods, such as Electron Transfer Dissociation (ETD), facilitates the generation of additional types of ions, including c and z ions.
  
While the charge state of a original peptide ion is dictated by the presence of protons, the same cannot be said for the fragment ions, a notion often misunderstood. Specifically, the +1 charge state of a, b, x, and y fragment ions results from an electron loss instead of a proton gain. In contrast, for c and y ions, the +1 charge state is achieved by acquiring a hydrogen atom and a proton from their corresponding neutral fragment. For charge states greater than 1, the remaining charges after +1 can be attributed to additional protons. 
  
**The structure for all +1 terminal fragment ions:**

![image](https://www.matrixscience.com/images/abcxyz.gif)     
  
### Internal Fragment Ions  
  
Internal fragment ions originate when fragmentation occurs at two separate points along the peptide, yielding a fragment that is 'internal' to the peptide's sequence. This process entails the breaking multiple peptide bonds, leading to the creation of a peptide fragment that isn't connected to either the N-terminus or C-terminus. The number of terminal fragment ions grows linearly with the length of the peptides sequence, while the number internal fragment ions grow exponentially  For this reason, internal fragments are typically ignored since they often overcomplicate the interpretation of the spectra. 
  
**The structure for a +1 BY internal fragment ion:**

![image](https://www.matrixscience.com/images/internal.gif)  
  
### Immonium Ions  
  
Immonium ions are unique types of fragment ions that consist of just one amino acid. Immonium ions won't tell you how the amino acids are arranged, they can confirm the inclusion of certain amino acids in the spectra.

**The structure for a +1 immonium ion:**

![image](https://www.matrixscience.com/images/immonium.gif)  
  
## Fragmentation Methods  
  
### Collision-Induced Dissociation (CID)  
  
CID fragments peptides by accelerating the ions using an electrical field to increase their kinetic energy, then allowing them to collide with neutral molecules (typically helium, nitrogen, or argon). The collision converts kinetic energy to internal energy, causing bonds to break and the molecule to fragmentation into smaller pieces. As previously stated CID is a soft ionization technique and there is only enough energy provided to break the peptide bond. CID produces mainly b and y ions, and sometimes a ions. In general, CID is excellent for generating sequence information for peptides and works well for a broad range of peptides.  
  
### Higher-Energy Collisional Dissociation (HCD)  
  
In Higher-energy Collisional Dissociation (HCD), the term "higher energy" refers to an increased radiofrequency (RF) voltage applied, not the energy used for fragmentation. The higher RF energy is used to better retain fragment ions. HCD uses a similar collision energy to CID, and as a result also produces mainly b and y ions. 

### Electron-Transfer Dissociation (ETD)  
  
In ETD an electron is transferred to the positively-charged protein or peptide, leading to fragmentation along the peptide backbone. ETD induces more energy into the peptide that CID, and as a result, produces fragments by breaking other bonds in the peptides. ETD commonly results in c and z ions. This method is especially useful for longer and more charged peptides.  

### Resources:
https://www.matrixscience.com/help/fragmentation_help.html
"""

HELP = """
# Help

The **Peptide Fragmenter** application is a tool for in silico fragmentation of peptide sequences. The app takes an 
amino acid sequence and calculates the fragment ions for a given charge range. 

The input sequence should be a valid amino acid sequence, with no additional characters or spaces. The sequence can
contain the following amino acids: ARNDCEQGHILKMFPSTWYVJX. Ambiguous amino acids
like B, Z are not supported. Proforma notation also
supports ambiguity of the sequence and modifications. Since the app is designed for fragmentation, the sequence 
must not contain any ambiguity.

Other than a lack of ambiguity, full proforma 2.0 notation is supported.
"""