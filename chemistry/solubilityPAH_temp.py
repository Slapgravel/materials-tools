import csv
import re

def temp_solvent(data):
    lines = data.split('\n')
    components = []
    experimental_values = []
    i = 0
    while i < len(lines):
        if 'COMPONENTS:' in lines[i]:
            i += 1
            while i < len(lines) and len(components) < 3:
                if '(1)' in lines[i] or '(2)' in lines[i] or '(3)' in lines[i]:
                    component = lines[i].strip().split(';')[0].split('(')[1].split(')')[1].strip()
                    components.append(component)
                i += 1
                # print(components)
        elif 'EXPERIMENTAL' in lines[i]:
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if not line or not any(char.isdigit() for char in line):  # Stop if line is empty or doesn't contain digits
                    break
                # print(line)
                values = re.split(r'\s+', line)
                if len(values) == 6:
                    experimental_values.append([values[0], values[1], values[2]])
                    experimental_values.append([values[3], values[4], values[5]])
                elif len(values) == 3:
                    experimental_values.append([values[0], values[1], values[2]])
                i += 1
        else:
            i += 1
    return [components] + experimental_values

# def solvent_only(data):
#     lines = data.split('\n')
#     components = []
#     variables = None
#     experimental_values = []

#     i = 0
#     while i < len(lines):
#         if 'COMPONENTS:' in lines[i]:
#             i += 1
#             while i < len(lines) and ('(1)' in lines[i] or '(2)' in lines[i] or '(3)' in lines[i]):
#                 component = lines[i].strip().split(';')[0].split('(')[1].split(')')[1].strip()
#                 components.append(component)
#                 i += 1
#         elif 'VARIABLES:' in lines[i]:
#             i += 1
#             while i < len(lines) and 'T/K =' in lines[i]:
#                 variables = lines[i].strip().split('=')[1].split(',')[0].strip()
#                 i += 1
#         elif 'EXPERIMENTAL VALUES' in lines[i]:
#             i += 1
#             while i < len(lines) and not ('@ x (s);' in lines[i] or '-----------------------------' in lines[i]):
#                 values = lines[i].strip().split()
#                 if len(values) == 3 and all(val.replace('.', '', 1).replace('-', '', 1).isdigit() for val in values):
#                     experimental_values.append(values)
#                 i += 1
#         else:
#             i += 1

#     output = []
#     for values in experimental_values:
#         row = ';'.join(components + values + [variables])
#         output.append(row.split(';'))

#     return output

def parse_file(filename):
    with open(filename, 'r') as f:
        content = f.read()
    datasets = [dataset for dataset in content.split('-----------------------------\n') if dataset.strip()]
    parsed_data = []
    for dataset in datasets:
        try:
            lines = dataset.split('\n')
            variables_line = None
            for line in lines:
                if 'VARIABLES' in line:
                    variables_index = lines.index(line)
                    for i in range(variables_index + 1, len(lines)):
                        if lines[i].strip():
                            variables_line = lines[i]
                            break
                    break
            if variables_line is not None:
                words = variables_line.lower().split()
                required_words = ["temperature", "solvent", "composition"]
                if all(word in words for word in required_words):
                    parsed_data.append(temp_solvent(dataset))
        except Exception as e:
            print(f"Error parsing dataset: {e}")
    return parsed_data
def write_to_csv(data, filename):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(['Component 1', 'Component 2', 'Component 3', 'x2(s)', 'x2b', 'x1', 'T/K'])
        for dataset in data:
            for row in dataset[1:]:
                writer.writerow(dataset[0] + [row[1]] + ['0'] + [row[2]] + [row[0]])

def read_toc(input_file):
    """Reads the Table of Contents from the input file."""
    toc = []
    with open(input_file, 'r') as f:
        for line in f:
            if '-----------------------------' in line:
                break
            elif '+' in line and '.' not in line:
                match = re.match(r'\s*([\w\s,.-]+)\s*\+\s*([\w\s,.-]+)', line)
                if match:
                    toc.append((match.group(1).strip(), match.group(2).strip()))
    return toc

def read_csv(csv_file):
    """Reads the CSV file and returns a set of (solvent1, solvent2) tuples."""
    csv_data = set()
    with open(csv_file, 'r') as f:
        reader = csv.reader(f, delimiter=';')
        next(reader)  # Skip header
        for row in reader:
            csv_data.add((row[1].lower(), row[2].lower()))
    return csv_data

def compare_toc_with_csv(toc, csv_data):
    """Compares the ToC with the CSV data and prints the results."""
    present = []
    absent = []
    for solvent1, solvent2 in toc:
        key = (solvent1.lower(), solvent2.lower())
        if key in csv_data:
            present.append((solvent1, solvent2))
        else:
            absent.append((solvent1, solvent2))
    return present, absent

name = 'placeholderNapthalene'
data = parse_file(name + '.txt')
write_to_csv(data, name + '.csv')

toc = read_toc(name + '.txt')
csv_data = read_csv(name + '.csv')
present, absent = compare_toc_with_csv(toc, csv_data)

print("Total:", len(present) + len(absent))

print("Present:", len(present))
for solvent1, solvent2 in present:
    print(f"{solvent1} + {solvent2}")

print("\nAbsent:", len(absent))
for solvent1, solvent2 in absent:
    print(f"{solvent1} + {solvent2}")