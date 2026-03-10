import csv
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
def main():
    input_file = 'rc.csv'
    output_file = 'output.csv'
    with open(input_file, 'r') as file:
        reader = csv.reader(file, delimiter='|')
        rows = list(reader)
    # Assuming the chemical names are in the first column
    for i, row in enumerate(rows[1:], start=1):
        chemical_name = row[1].strip()
        cas_number = search_chemical_info_by_name(chemical_name)
        rows[i].append(cas_number)
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='|')
        # Add header for the new column
        rows[1].append('CAS Number')
        writer.writerows(rows)
if __name__ == "__main__":
    main()