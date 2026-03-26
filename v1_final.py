#======================================
#Title: Retrospective analysis of clinical and environmental genotyping results revealing persistence of Pseudomonas aeruginosa in the water system of a large tertiary children's hospital in England.
#Author: Yu Wan, Esha Sheth
#Date: March 2026
#Version: 1.0

#Description: This script extracts data from UKHSA pdf reports into a tsv file.

import re
import os
import sys
from PyPDF2 import PdfReader

def extract_info_from_pdf(pdf_path):
    try:
        reader = PdfReader(pdf_path)
        text = ""
        for page in reader.pages:
            text += page.extract_text()

        # Extraction fields
        phe_ref = re.search(r'PHE ref\.?\s*No\.?\s+([A-Z0-9 ]+)', text, re.IGNORECASE)
        date_received = re.search(r'Date received\s+(\d{2}\.\d{2}\.\d{4})', text, re.IGNORECASE)
        date_collected = re.search(r'Date of collection\s+(\d{2}\.\d{2}\.\d{4})', text,re.IGNORECASE)
        hospital_ID = re.search(r'Hospital No\.?\s+([A-Z0-9]+)', text,re.IGNORECASE)
        vntr_profile = re.search(r'VNTR Profile:?\s+([0-9,\- ]+)', text,re.IGNORECASE)
        organism = re.search(r'Opportunistic Pathogens Section\s*[\n\r]+1\.\s*([A-Za-z\- ]+)', text, re.IGNORECASE)
        senders_ref = re.search(r"Sender'?s?\s*ref\.?\s*No\.?\s*([A-Za-z0-9\-/]+)", text, re.IGNORECASE)
        isolation_site = re.search(r'Isolation site\s+([A-Za-z\- ]+)', text, re.IGNORECASE)
        patient_name = re.search(r'Name\s+([A-Za-z\-\' ]+,\s*[A-Za-z\-\' ]+)', text, re.IGNORECASE)

        return [
            senders_ref.group(1).strip() if senders_ref else "",
            phe_ref.group(1).replace(' ', '').strip() if phe_ref else "Not found",
            date_received.group(1).strip() if date_received else "",
            date_collected.group(1).strip() if date_collected else "",
            hospital_ID.group(1).strip() if hospital_ID else "",
            vntr_profile.group(1).replace(' ', '').strip() if vntr_profile else "",
            organism.group(1).strip() if organism else "",
            isolation_site.group(1).strip() if isolation_site else "",
            patient_name.group(1).strip() if patient_name else "",
            os.path.basename(pdf_path)
        ]
    except Exception as e:
        print(f"❌ Error processing {pdf_path}: {e}")
        return ["Error", "", "", "", "", "", "", "", "", os.path.basename(pdf_path)]

def process_all_pdfs(folder_path):
    output_path = os.path.join(folder_path, "summary.tsv")

    with open(output_path, "w", encoding="utf-8") as out:
        out.write("Senders_Ref\tPHE_Ref\tDate_Received\tDate_Collected\tHospital_ID\tVNTR_Profile\tOrganism\tIsolation_Site\tPatient_Name\tFilename\n")
        for root, _, files in os.walk(folder_path):
            for file in files:
                if file.lower().endswith(".pdf"):
                    pdf_path = os.path.join(root, file)
                    info = extract_info_from_pdf(pdf_path)
                    out.write("\t".join(info) + "\n")

    print(f"\n✅ Extraction complete! Summary saved to:\n{output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python v1.py <folder_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    process_all_pdfs(folder_path)