import re
import requests
from rich import print


def search_chemical_info_by_name(chemical_name):
    base_url = "https://commonchemistry.cas.org/api"
    search_endpoint = "/search"
    search_url = base_url + search_endpoint


    params = {"q": chemical_name}


    try:
        print(f"Querying for chemical name: {chemical_name}")
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data.get("count", 0) > 0:
            result = data["results"][0]
            cas_number = result.get("rn", "CAS number not found")
            return cas_number
        else:
            return "CAS number not found"
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")
        print(f"Response text: {response.text}")
        return "Failed to fetch data from the Common Chemistry API."
    except ValueError:
        print("Error decoding JSON response")
        return "Failed to decode JSON response."


def search_chemical_info_by_cas(cas_number):
    base_url = "https://commonchemistry.cas.org/api"
    search_endpoint = "/search"
    search_url = base_url + search_endpoint


    params = {"q": cas_number}


    try:
        print(f"Querying for CAS number: {cas_number}")
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data.get("count", 0) > 0:
            result = data["results"][0]
            chemical_name = result.get("name", "Chemical name not found")
            return chemical_name
        else:
            return "Chemical name not found"
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")
        print(f"Response text: {response.text}")
        return "Failed to fetch data from the Common Chemistry API."
    except ValueError:
        print("Error decoding JSON response")
        return "Failed to decode JSON response."


def get_physical_properties(cas_number):
    url = "https://commonchemistry.cas.org/api/detail"
    payload = {'cas_rn': cas_number}
    response = requests.get(url, params=payload)
    data = response.json()
    return data


def is_cas_number(query):
    cas_pattern = r"^\d{2,7}-\d{2}-\d$"
    return re.match(cas_pattern, query) is not None

def print_chemical_info(properties):
    # print(f"CAS Number: {properties['rn']}")
    # print(f"Name: {properties['name']}")
    for prop in properties['experimentalProperties']:
        if prop['name'] == 'Melting Point':
            print(f"Melting Point: {prop['property']}")
        elif prop['name'] == 'Boiling Point':
            print(f"Boiling Point: {prop['property']}")
        elif prop['name'] == 'Density':
            print(f"Density: {prop['property']}")
    print(f"Molecular Mass: {properties['molecularMass']}")
    print(f"SMILES: {properties['smile']}")
    formula = properties['molecularFormula']
# Define the subscript characters
    subscripts = {
        '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
        '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'
    }
    # Loop through numbers 1-99 and replace subscript tags
    for i in range(1, 100):
        num_str = str(i)
        subscript_tag = f'<sub>{num_str}</sub>'
        subscript_replacement = ''.join(subscripts[digit] for digit in num_str)
        formula = formula.replace(subscript_tag, subscript_replacement)
    print(f"Chemical Formula: {formula}")
    for synonym in properties['synonyms']:
        print(f"- {synonym}")
    # Add more properties as needed
# Call the function with the properties dictionary


def main():
    query = input("Enter a chemical name or CAS number: ").strip()


    if is_cas_number(query):
        result = search_chemical_info_by_cas(query)
        print(f"Chemical Name: {result}")
        properties = get_physical_properties(query)
        print_chemical_info(properties)
    else:
        result = search_chemical_info_by_name(query)
        print(f"CAS Number: {result}")
        properties = get_physical_properties(result)
        print_chemical_info(properties)


if __name__ == "__main__":
    main()