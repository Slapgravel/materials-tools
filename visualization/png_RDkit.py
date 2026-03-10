import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import os
from PIL import Image

# Setup output folder
script_dir = os.path.dirname(os.path.abspath(__file__))
nested_folder = os.path.join(script_dir, 'molecular_images')
os.makedirs(nested_folder, exist_ok=True)

# Load your CSV
df = pd.read_csv('Crystal_backend.csv')
first_col = df.columns[0]
smiles_col = 'SMILES' if 'SMILES' in df.columns else df.columns[1]

for idx, row in df.iterrows():
    smiles = row[smiles_col]
    filename = f"{str(row[first_col]).zfill(4)}.png"
    filepath = os.path.join(nested_folder, filename)
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Draw molecule (default is white background, RGB)
            img = Draw.MolToImage(mol, size=(600, 600), bgColor=None)
            # Convert to RGBA and make white pixels transparent
            img = img.convert("RGBA")
            datas = img.getdata()
            newData = []
            for item in datas:
                # If pixel is white, make it transparent
                if item[0] > 240 and item[1] > 240 and item[2] > 240:
                    newData.append((255, 255, 255, 0))
                else:
                    newData.append(item)
            img.putdata(newData)
            # Crop to non-transparent area
            alpha = img.split()[-1]
            bbox = alpha.getbbox()
            if bbox:
                cropped_img = img.crop(bbox)
                cropped_img.save(filepath)
            else:
                img.save(filepath)
        else:
            print(f"Invalid SMILES at row {idx}: {smiles}")
    except Exception as e:
        print(f"Error at row {idx}: {e}")

print(f"Images saved to {nested_folder}")