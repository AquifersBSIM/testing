import requests
from bs4 import BeautifulSoup
from typing import List, Dict, Optional
import re

class RCSBLigandParser:
    def __init__(self):
        self.ligands = {}
        self.structure_id = None
        self.base_url = "https://www.rcsb.org"

    def fetch_structure_page(self, pdb_id: str) -> bool:
        """
        Fetch and parse the structure page from RCSB
        """
        try:
            self.structure_id = pdb_id.lower()
            url = f"https://www.rcsb.org/structure/{self.structure_id}"
            response = requests.get(url)
            
            if response.status_code == 200:
                # Parse the HTML content
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Find all ligand rows by the pattern ligand_row_X where X is the ligand ID
                ligand_rows = soup.find_all('tr', id=re.compile(r'ligand_row_\w+'))
                
                if ligand_rows:
                    self.parse_ligand_rows(ligand_rows)
                    return True
                else:
                    print("No ligand information found on the page")
                    return False
                    
            else:
                print(f"Failed to fetch page. Status code: {response.status_code}")
                return False
                
        except Exception as e:
            print(f"Error fetching structure page: {e}")
            return False
            
    def parse_ligand_rows(self, ligand_rows) -> None:
        """
        Parse individual ligand rows and extract information,
        including chain (auth) data.
        """
        try:
            for row in ligand_rows:
                # Extract ligand ID from the row ID attribute
                row_id = row.get('id', '')
                if row_id:
                    ligand_id = row_id.replace('ligand_row_', '')

                    # Find the ligand link within the row
                    ligand_link = row.find('a', href=re.compile(r'/ligand/'))

                    # Find the strong tag containing the ligand name
                    strong_tag = row.find('strong')
                    ligand_name = strong_tag.text if strong_tag else None

                    if ligand_link:
                        ligand_url = self.base_url + ligand_link['href']
                        self.ligands[ligand_id] = {
                            'id': ligand_id,
                            'url': ligand_url,
                            'name': ligand_name
                        }
                        
                        # --- NEW: Parse chain info and store it ---
                        chain_info = self.parse_ligand_auths(row)
                        self.ligands[ligand_id]['chains'] = chain_info

        except Exception as e:
            print(f"Error parsing ligand rows: {e}")

    def parse_ligand_auths(self, row):
        """
        Parse the chain/residue labels (e.g. "H [auth A], M [auth B], R [auth C], ...")
        from the ligand row and return them as a list of dicts. If the same chain
        (auth ID) appears more than once, we only keep one entry.
        """
        chain_info = []
        seen_chain_ids = set()
        
        # Attempt to find the second <td>, which often contains the chain info
        tds = row.find_all('td')
        if len(tds) > 1:
            # Get the entire text, merging <br> with spaces
            chain_text = tds[1].get_text(separator=" ", strip=True)
            # Example chain_text might be: "H [auth A], M [auth B], R [auth C], W [auth D]"
            
            # Use a regex to capture the pattern:
            #   (\w) \[auth ([A-Za-z0-9]+)\]
            pattern = re.compile(r'(\w)\s*\[auth\s+([A-Za-z0-9]+)\]')
            matches = pattern.findall(chain_text)
            
            # Each match is a tuple (residue_label, chain_id)
            for residue_label, chain_id in matches:
                if chain_id not in seen_chain_ids:
                    seen_chain_ids.add(chain_id)
                    chain_info.append({
                        'residue_label': residue_label,
                        'chain_id': chain_id
                    })
        
        return chain_info

    def display_results(self) -> None:
        """
        Display the found ligands and their information in a two-column table format
        (Ligand | Chains) with Unicode box-drawing characters.
        """
        if not self.ligands:
            print("No ligands found")
            return

        # Prepare the rows for two columns: (Ligand Info, Chain Info)
        header = (F"Ligand(s) found in structure {self.structure_id.upper()}:", "CHAIN(S)")
        rows = []

        # We’ll collect a top “title” line separately
        title_line = f"LIGAND: {self.structure_id.upper()}"
        #subtitle_line = f"Ligands found in structure {self.structure_id.upper()}:"

        for index, (ligand_id, info) in enumerate(self.ligands.items(), 1):
            # Build the left column text
            if info.get('name'):
                ligand_text = f"{index}. {ligand_id} ({info['name']})"
            else:
                ligand_text = f"{index}. {ligand_id}"

            # Build the right column text (chain info)
            chain_data = info.get('chains', [])
            if chain_data:
                # e.g., "H [auth A], M [auth B], R [auth C], ..."
                chain_str_list = [
                    f"CHAIN {cd['chain_id']}"
                    for cd in chain_data
                ]
                chain_text = ", ".join(chain_str_list)
            else:
                chain_text = "N/A"

            rows.append((ligand_text, chain_text))

        # Determine column widths
        col1_width = max(len(r[0]) for r in rows + [header])
        col2_width = max(len(r[1]) for r in rows + [header])

        # Optional: pad them a bit if you like
        col1_width += 2
        col2_width += 2

        # Box drawing characters
        top_left = "┏"
        top_right = "┓"
        bottom_left = "┗"
        bottom_right = "┛"
        horizontal = "━"
        vertical = "┃"
        left_connector = "┣"
        right_connector = "┫"
        middle_connector = "┳"
        middle_connector2 = "╋"
        vertical_line = "┃"

        # A helper to print a full horizontal line across both columns
        def print_horizontal_line(left_char, mid_char, right_char, width1, width2):
            print(
                left_char
                + horizontal * width1
                + mid_char
                + horizontal * width2
                + right_char
            )

        # Print title lines in a single “box”
        title_max_width = len(title_line) + 2

        print(top_left + horizontal * title_max_width + top_right)
        print(vertical + title_line.center(title_max_width) + vertical)
        print(bottom_left + horizontal * title_max_width + bottom_right)

        # Print the table header
        print_horizontal_line(top_left, middle_connector, top_right, col1_width, col2_width)
        # Header row
        print(
            f"{vertical} {header[0].ljust(col1_width - 1)}"
            f"{vertical} {header[1].ljust(col2_width - 1)}{vertical}"
        )
        # Separator after header
        print_horizontal_line(left_connector, middle_connector2, right_connector, col1_width, col2_width)

        # Print each row
        for row_left, row_right in rows:
            print(
                f"{vertical} {row_left.ljust(col1_width - 1)}"
                f"{vertical} {row_right.ljust(col2_width - 1)}{vertical}"
            )

        # Print bottom line
        print_horizontal_line(bottom_left, middle_connector, bottom_right, col1_width, col2_width)

