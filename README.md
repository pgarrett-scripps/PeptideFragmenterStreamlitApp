# PeptideFragmenter

PeptideFragmenter is a web application built with Streamlit that provides peptide fragmentation functionality. Users can input a peptide sequence and obtain calculated fragment ions for a specified charge range. The application also allows users to match fragment ions to input spectra.

https://peptidefragmenter.streamlit.app

## Features

- **Peptide Fragmentation**: Input an amino acid sequence and get the corresponding fragment ions for a given charge range. The app also supports modifications, which can be provided in parentheses with the mass difference in Daltons.
- **Spectra Matching**: Add spectra to match fragment ions. The matching is done based on user-defined tolerances, and the results are visualized using a bar chart.
- **Interactive Visualization**: The app offers interactive visualizations for fragment segments and spectra matching results.
- **Downloadable Data**: Users can download fragment ion data and spectra results in CSV format.
- **Wiki and Help Tabs**: The app provides a Wiki for general information and a Help section to guide users through the app's functionalities.

## Usage

1. **Input Parameters**: On the sidebar, provide the peptide sequence, charge range, mass type, and fragment types. You can also opt to include internal fragments.
2. **View Fragment Ions**: The `Results` tab displays fragment ions based on the input parameters. You can visualize fragment segments and download fragment ion data.
3. **Input Spectra**: In the `Spectra` tab, provide your spectra data and set the desired tolerances. The app will match the fragment ions to the input spectra and provide a visualization of the results. You can also view scores (e.g., Hyperscore, Binomial Score) and download spectra results.
4. **Additional Information**: The `Wiki` and `Help` tabs offer general information and guidance on using the app, respectively.

## Installation

This application requires Python and certain libraries. Ensure you have the following libraries installed:
- `pandas`
- `peptacular`
- `streamlit`
- `plotly`

You can install these libraries using pip:
```bash
pip install pandas peptacular streamlit plotly
```

## Running the App

To run the application locally:
```bash
streamlit run app.py
```

---

## URL Query Parameters

PeptideFragmenter supports custom configurations via URL query parameters. By adjusting these parameters in the URL, users can quickly access the application with specific settings. This section describes the supported query parameters:

- **`sequence`**: Represents the peptide sequence to be fragmented.
    - Example: `...?sequence=PEPTIDE`

- **`min_charge`**: The minimum charge to be considered for fragmentation.
    - Example: `...?min_charge=1`

- **`max_charge`**: The maximum charge to be considered for fragmentation.
    - Example: `...?max_charge=3`

- **`mass_type`**: Specifies the type of mass to use for fragment calculations. Supported values are `monoisotopic` and `average`.
    - Example: `...?mass_type=monoisotopic`

- **`fragment_types`**: A list of fragment ion types to consider. Possible values include `a`, `b`, `c`, `x`, `y`, and `z`.
    - Example: `...?fragment_types=abc`

### Using Multiple Parameters

You can combine multiple query parameters using the `&` character. For example:
```
https://peptidefragmenter.streamlit.app?sequence=PEPTIDE&min_charge=1&max_charge=3&mass_type=monoisotopic&fragment_types=abc
```

By using these query parameters, you can bookmark specific configurations, share links with colleagues, or embed the application with preset settings in other web pages.
