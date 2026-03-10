from ccdc import io
from rich import print


def extract_data(refcode):
    # Retrieve the entry from the CSD database
    reader = io.EntryReader('CSD')
    entry = reader.entry(refcode)


    # Extract the data
    formula = entry.formula
    density = entry.calculated_density
    chemical_name = entry.chemical_name_as_html
    color = entry.color
    habit = entry.habit
    melting_point = entry.melting_point
    smiles = entry.molecule.smiles


    crystal = entry.crystal
    a = crystal.cell_lengths.a
    b = crystal.cell_lengths.b
    c = crystal.cell_lengths.c
    alpha = crystal.cell_angles.alpha
    beta = crystal.cell_angles.beta
    gamma = crystal.cell_angles.gamma
    z_value = crystal.z_value
    z_prime = crystal.z_prime
    spacegroup_symbol = crystal.spacegroup_symbol

    publication = entry.publication.doi


    # Print the extracted data
    print(f"Chemical Name: {chemical_name}")
    print(f"RefCode: {refcode}")
    print(f"Formula: {formula}")
    print(f"smiles: {smiles}")
    print(f"Melting Point: {melting_point}")
    print(f"Density: {density}")
    print(f"Color: {color}")
    print(f"Habit: {habit}")
    print(f"spacegroup_symbol: {spacegroup_symbol}")
    print(f"z_value: {z_value}")
    print(f"z_prime: {z_prime}")
    print(f"a: {a}")
    print(f"b: {b}")
    print(f"c: {c}")
    print(f"alpha: {alpha}")
    print(f"beta: {beta}")
    print(f"glamma: {gamma}")
    print(f"Publication: {publication}")


def main():
    refcode = input("Enter a refcode: ").upper()
    extract_data(refcode)


if __name__ == "__main__":
    main()