import pandas as pd
from pandas import DataFrame
from tkinter import Tk
from tkinter.filedialog import askopenfilename


def read_cat() -> DataFrame:
    Tk().withdraw()
    filename = askopenfilename()
    dataframe: DataFrame = pd.read_fwf(  # type: ignore
     filename,
     header=None,
     names=[
         "Frequency",
         "Error",
         "Integrated Intensity",
         "Deg of Freedom",
         "Lower State Energy (cm-1)",
         "Upper State Deg",
         "Tag",
         "QNFMT",
         "Upper N",
         "Upper Ka",
         "Upper Kc",
         "Upper V",
         "Upper J",
         "Upper F",
         "Lower N",
         "Lower Ka",
         "Lower Kc",
         "Lower V",
         "Lower J",
         "Lower F",
     ],
     widths=[13, 8, 8, 2, 10, 3, 7, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    )
    return dataframe


def map_strings_to_numeric(df: DataFrame):
    df = df.astype(str)
    df.drop(
     columns=["Upper J", "Upper F", "Lower J", "Lower F"],
     inplace=True,
    )
    quantum_list = [
     "Upper N",
     "Upper Ka",
     "Upper Kc",
     "Upper V",
     "Lower N",
     "Lower Ka",
     "Lower Kc",
     "Lower V",
    ]
    for column_name in quantum_list:
     column = df[column_name]
     numeric_values = []

     for value in column:
         if value.isdigit():
             numeric_values.append(value)
             continue

         if pd.isnull(value):  # Skip empty values
             numeric_values.append(value)
             continue

         value = str(value)
         letter = value[0]  # Extract the first character (letter)
         number = int(value[1:])  # Extract the remaining characters as the number
         letter_value = (
             ord(letter.upper()) - 65
         )  # Convert letter to value (A=0, B=1, ...)
         if letter.islower():  # Check if the letter is lowercase
             numeric_value = (
                 -10 * (letter_value + 1) - number
             )  # Calculate the final numeric value for lowercase letter
         else:
             numeric_value = (
                 100 + letter_value * 10 + number
             )  # Calculate the final numeric value for uppercase letter
         numeric_values.append(numeric_value)

     df[column_name] = numeric_values
    df = df.apply(pd.to_numeric)
    return df

